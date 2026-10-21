#!/usr/bin/env python3
"""Streaming exact certifier for the large JK-H12 BABB/BBBB P5 shards.

The frozen P5 producer retained every public gate field and every exact degree
profile in memory.  That representation is unnecessary for the two largest
shards.  This certifier instead keeps one compact transcript row per orbit.
SHA-256 is used only as a negative filter: every same-digest bucket is replayed
and partitioned by the literal degree-profile bytes.  Literal collisions are
then replayed through the frozen canonical color refinement.  The same rule is
used against all previously completed P5 rows and same-facet star rows.

The final archive is a small JSON summary plus deterministic gzip JSONL rows.
No Gram matrix, degree-profile byte string, or temporary database is archived.
An in-process peak-RSS check complements the repository's process watchdog.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import resource
import sys
import tempfile
import time
from typing import Any, Iterable, Iterator, Sequence

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import scripts.certify_jacobian_kummer_h12_five_vertex_p5_b_le_2_census as p5
import scripts.certify_jacobian_kummer_h12_five_vertex_star_census as star
import scripts.certify_jacobian_kummer_h12_four_vertex_subtree_census as four


SCHEMA = "coxeter.jacobian_kummer.h12_five_vertex_p5_streaming_shard.v1"
STATUS_DATE = "2026-08-09"
SUPPORTED_PATTERNS = ("BABB", "BBBB")
RSS_CEILING_BYTES = 1 << 30
FIELD_RE = re.compile(
    rb'"ordinary_degree_profile_multiset_sha256": "([0-9a-f]{64})"'
)

P5_SCRIPT = ROOT / "scripts/certify_jacobian_kummer_h12_five_vertex_p5_b_le_2_census.py"
P5_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_p5_b_le_2_census.json"
ABBB_SCRIPT = ROOT / "scripts/certify_jacobian_kummer_h12_five_vertex_p5_abbb_census.py"
ABBB_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_p5_abbb_census.json"
STAR_SCRIPT = ROOT / "scripts/certify_jacobian_kummer_h12_five_vertex_star_census.py"
STAR_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_star_census.json"
STAR_ABBB_SCRIPT = ROOT / "scripts/certify_jacobian_kummer_h12_star_abbb_exact_separation.py"
STAR_ABBB_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h12_star_abbb_exact_separation.json"

PINNED_SHA256 = {
    "P5_B_le_2_builder": "d7c5c02803911f8d205a63c59027d8d9584cce998de5efa0a7acb134d56d00e4",
    "P5_B_le_2_archive": "32943e8531bdaa1874af41d45635ea792c9337694adc4f20dfc152aac3d7713a",
    "P5_ABBB_builder": "6217c4902b8a66bbb4cef57e4a2edc188da723d6ceffcc2cad82cd8d172eb18e",
    "P5_ABBB_archive": "b535bbe5d1639358e5a336ce606a4d7eae8dedd0b55835f9ee057321ac394299",
    "five_vertex_star_builder": "bbb99d59726e0b6bcf0d6b3748823ca9bc8f1efe5ecbe5c29a16f4a5db6a53fb",
    "five_vertex_star_archive": "21fbd0946b910360cdd4096392e4e629907f971be85359d5b87e312a78d4f309",
    "star_ABBB_exact_builder": "f27f18c36531c48f75e93f5f72b36da694ef0de6f64175d7e84b38186991b926",
    "star_ABBB_exact_archive": "7a328edb6cbedafa13c32140c2f2f1262a8b1159617a1fc59d5ad8d2165c3fe6",
}


@dataclass
class CompactRow:
    colors: tuple[int, int, int, int]
    orbit_size: int
    repetition_count: int
    facet_count: int
    degree_sha256: str
    gate_sha256: str
    boundary_root_sha256: str
    refinement_sha256: str | None = None
    refinement_discrete: bool | None = None
    refinement_final_cells: int | None = None


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def _json_sha256(value: Any) -> str:
    return hashlib.sha256(_json_bytes(value)).hexdigest()


def _stream_update(digest: Any, value: bytes) -> None:
    digest.update(len(value).to_bytes(8, "big"))
    digest.update(value)


def _rss_bytes() -> int:
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(peak if sys.platform == "darwin" else peak * 1024)


def _rss_guard() -> None:
    if _rss_bytes() > RSS_CEILING_BYTES:
        raise MemoryError("streaming P5 replay exceeded its explicit 1 GiB RSS ceiling")


def _pinned_gate() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    paths = {
        "P5_B_le_2_builder": P5_SCRIPT,
        "P5_B_le_2_archive": P5_ARCHIVE,
        "P5_ABBB_builder": ABBB_SCRIPT,
        "P5_ABBB_archive": ABBB_ARCHIVE,
        "five_vertex_star_builder": STAR_SCRIPT,
        "five_vertex_star_archive": STAR_ARCHIVE,
        "star_ABBB_exact_builder": STAR_ABBB_SCRIPT,
        "star_ABBB_exact_archive": STAR_ABBB_ARCHIVE,
    }
    actual = {name: _sha256(path) for name, path in paths.items()}
    if actual != PINNED_SHA256:
        raise ArithmeticError(f"a frozen streaming-shard prerequisite changed: {actual}")
    small, pilot, denominator = p5._pinned_predecessors()
    exact = json.loads(STAR_ABBB_ARCHIVE.read_text(encoding="utf-8"))
    if not (
        exact["all_checks_pass"]
        and all(exact["checks"].values())
        and exact["theorem"]["complete_star_plus_ABBB_pairwise_nonisomorphic_type_count"]
        == 53_119
        and pilot["result"]["reflection_index_and_volume_multiplier"] == 5
        and pilot["result"]["finite_volume"]
        and not pilot["result"]["compact"]
        and denominator["tree_shapes"]["P5"]["seed_and_tree_automorphism_orbits"]
        == 201_747
    ):
        raise ArithmeticError("a frozen streaming-shard theorem is not green")
    return small, pilot, denominator


def _archived_degree_hashes(path: Path, expected: int) -> tuple[str, ...]:
    values: list[str] = []
    with path.open("rb") as stream:
        for line in stream:
            match = FIELD_RE.search(line)
            if match:
                values.append(match.group(1).decode("ascii"))
    if len(values) != expected:
        raise ArithmeticError(f"the archived degree-hash row count changed for {path}")
    return tuple(values)


def _degree_replay(
    engine: Any, colors: tuple[int, int, int, int], expected_digest: str | None = None
) -> tuple[bytes, Any]:
    boundary, _boundary_data = engine.boundary(colors)
    labels, degree_key, _ordinary = engine.ordinary_data(boundary)
    if expected_digest is not None and hashlib.sha256(degree_key).hexdigest() != expected_digest:
        raise ArithmeticError("a degree-profile digest changed on exact replay")
    return degree_key, labels


def _refinement(
    engine: Any,
    colors: tuple[int, int, int, int],
    expected_degree: bytes,
) -> tuple[bytes, dict[str, Any]]:
    degree, labels = _degree_replay(engine, colors)
    if degree != expected_degree:
        raise ArithmeticError("a literal degree profile changed before refinement")
    return four._canonical_color_refinement(labels)


def _prior_rows(
    pattern: str,
    actions: Sequence[p5.Permutation],
) -> Iterator[tuple[str, int, tuple[int, int, int, int], str, Any]]:
    """Yield prior same-facet rows with their frozen digest and engine class."""

    if pattern == "BABB":
        abbb_rows = p5._pattern_color_orbits("ABBB", actions)
        hashes = _archived_degree_hashes(ABBB_ARCHIVE, 36_247)
        for index, (row, digest) in enumerate(zip(abbb_rows, hashes, strict=True), start=1):
            yield (
                "P5_ABBB",
                index,
                tuple(int(value) for value in row["colors_zero_based_stock"]),
                digest,
                "P5",
            )

    star_rows = star._star_color_orbits(actions)
    hashes = _archived_degree_hashes(STAR_ARCHIVE, 16_872)
    target_b = pattern.count("B")
    expected = 6_117 if target_b == 3 else 8_985
    emitted = 0
    for index, (row, digest) in enumerate(zip(star_rows, hashes, strict=True), start=1):
        colors = tuple(int(value) for value in row["colors_zero_based_stock"])
        if sum(color >= p5.A_SIZE for color in colors) != target_b:
            continue
        emitted += 1
        yield ("K1,4_star", index, colors, digest, "star")
    if emitted != expected:
        raise ArithmeticError("the prior same-facet star denominator changed")


def _compact_row_payload(index: int, pattern: str, row: CompactRow) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "orbit_number": index,
        "edge_class_word": pattern,
        "edge_colors_one_based": [color + p5.STOCK_START + 1 for color in row.colors],
        "unrooted_coloring_orbit_size_under_seed_action": row.orbit_size,
        "nonadjacent_color_repetition_count": row.repetition_count,
        "facet_count": row.facet_count,
        "marked_boundary_root_sha256": row.boundary_root_sha256,
        "ordinary_degree_profile_multiset_sha256": row.degree_sha256,
        "exact_gate_record_sha256": row.gate_sha256,
        "color_refinement_applied": row.refinement_sha256 is not None,
    }
    if row.refinement_sha256 is not None:
        payload.update(
            {
                "color_refinement_invariant_sha256": row.refinement_sha256,
                "color_refinement_discrete_canonical_form": row.refinement_discrete,
                "color_refinement_final_cell_count": row.refinement_final_cells,
            }
        )
    return payload


def _write_rows(path: Path, pattern: str, rows: Sequence[CompactRow]) -> tuple[str, str]:
    transcript = hashlib.sha256()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as compressed:
            for index, row in enumerate(rows, start=1):
                rendered = _json_bytes(_compact_row_payload(index, pattern, row)) + b"\n"
                transcript.update(rendered)
                compressed.write(rendered)
    return transcript.hexdigest(), _sha256(path)


def _build(
    pattern: str,
    *,
    rows_output: Path,
    verbose: bool,
) -> dict[str, Any]:
    if pattern not in SUPPORTED_PATTERNS:
        raise ArithmeticError("unsupported streaming P5 shard")
    _small, pilot, denominator = _pinned_gate()
    carrier, gram, roots, interior, bounded = p5.base._load_inputs()
    actions = star._restricted_actions()
    if not all(
        all((color < p5.A_SIZE) == (action[color] < p5.A_SIZE) for color in range(p5.STOCK_SIZE))
        for action in actions
    ):
        raise ArithmeticError("the exact seed action no longer preserves A/B")
    stock_gram_sha256 = star._root_gram_action_gate(gram, roots, actions)

    started = time.perf_counter()
    color_rows = p5._pattern_color_orbits(pattern, actions)
    burnside = p5._burnside_pattern_orbits(pattern, actions)
    expected_count = p5.EXPECTED_PATTERN_ORBITS[pattern]
    if len(color_rows) != expected_count or burnside["orbit_count"] != expected_count:
        raise ArithmeticError("the selected shard denominator changed")
    if verbose:
        print(
            f"orbit_cover={pattern} rows={len(color_rows)} seconds={time.perf_counter()-started:.2f}",
            flush=True,
        )
    _rss_guard()

    engine = p5.ExactP5BoundaryEngine(gram, roots, interior)
    digest_groups: defaultdict[str, list[int]] = defaultdict(list)
    rows: list[CompactRow] = []
    relation_histograms: Counter[tuple[tuple[str, int], ...]] = Counter()
    facet_histogram: Counter[int] = Counter()
    raw_facet_histogram: Counter[int] = Counter()
    orbit_size_histogram: Counter[int] = Counter()
    repetition_histogram: Counter[int] = Counter()
    gate_stream = hashlib.sha256()
    degree_stream = hashlib.sha256()
    coordinate_maximum = transformed_root_maximum = transform_maximum = pairing_maximum = 0
    minimum_interior_pairing: int | None = None
    boundary_started = time.perf_counter()
    for index, color_row in enumerate(color_rows):
        colors = tuple(int(value) for value in color_row["colors_zero_based_stock"])
        boundary, boundary_data = engine.boundary(colors)
        labels, degree_key, ordinary_data = engine.ordinary_data(boundary)
        del labels
        facets = int(boundary_data["facet_count"])
        expected_facets = 752 + 9 * pattern.count("B")
        if facets != expected_facets:
            raise ArithmeticError("the P5 facet formula changed")
        gate_data = {**boundary_data, **ordinary_data}
        if not (
            boundary_data["chamber_count"] == 5
            and boundary_data["internal_mirror_count"] == 4
            and boundary_data["strict_interior_square"] > 0
            and boundary_data["minimum_strict_interior_facet_pairing"] > 0
            and boundary_data["primitive_integral_reflection_root_count"] == facets
            and boundary_data["fixed_boundary_core_Smith_diagonal"] == [1] * p5.RANK
            and boundary_data["full_linear_rank"] == p5.RANK
            and ordinary_data["primitive_integral_reflection_root_count"] == facets
            and sum(ordinary_data["ordinary_Coxeter_relation_histogram"].values())
            == facets * (facets - 1) // 2
        ):
            raise ArithmeticError("a per-row arithmetic/Coxeter gate failed")
        degree_digest = hashlib.sha256(degree_key).hexdigest()
        digest_groups[degree_digest].append(index)
        _stream_update(degree_stream, degree_key)
        gate_rendered = _json_bytes(gate_data)
        _stream_update(gate_stream, gate_rendered)
        orbit_size = int(color_row["unrooted_coloring_orbit_size_under_seed_action"])
        repetition = 4 - len(set(colors))
        rows.append(
            CompactRow(
                colors=colors,
                orbit_size=orbit_size,
                repetition_count=repetition,
                facet_count=facets,
                degree_sha256=degree_digest,
                gate_sha256=hashlib.sha256(gate_rendered).hexdigest(),
                boundary_root_sha256=boundary_data["marked_boundary_root_sha256"],
            )
        )
        relation_histograms[tuple(sorted(ordinary_data["ordinary_Coxeter_relation_histogram"].items()))] += 1
        facet_histogram[facets] += 1
        raw_facet_histogram[facets] += orbit_size
        orbit_size_histogram[orbit_size] += 1
        repetition_histogram[repetition] += 1
        coordinate_maximum = max(coordinate_maximum, boundary_data["maximum_boundary_coordinate_absolute"])
        transformed_root_maximum = max(transformed_root_maximum, boundary_data["maximum_transformed_seed_root_coordinate_absolute"])
        transform_maximum = max(transform_maximum, boundary_data["maximum_chamber_transform_coordinate_absolute"])
        pairing_maximum = max(pairing_maximum, ordinary_data["maximum_boundary_pairing_absolute"])
        value = int(boundary_data["minimum_strict_interior_facet_pairing"])
        minimum_interior_pairing = value if minimum_interior_pairing is None else min(minimum_interior_pairing, value)
        if (index + 1) % 100 == 0:
            _rss_guard()
        if verbose and ((index + 1) % 500 == 0 or index + 1 == expected_count):
            print(
                f"boundary_progress={index+1}/{expected_count} elapsed_seconds={time.perf_counter()-boundary_started:.1f} peak_rss_gib={_rss_bytes()/(1<<30):.3f}",
                flush=True,
            )

    digest_collision_groups = [indices for indices in digest_groups.values() if len(indices) > 1]
    literal_replayed = 0
    literal_degree_group_count = 0
    literal_degree_member_count = 0
    refined_current_indices: set[int] = set()
    refinement_discrete = 0
    collision_started = time.perf_counter()
    for group_number, digest_group in enumerate(digest_collision_groups, start=1):
        exact: defaultdict[bytes, list[int]] = defaultdict(list)
        for index in digest_group:
            degree, _labels = _degree_replay(
                engine, rows[index].colors, rows[index].degree_sha256
            )
            exact[degree].append(index)
            literal_replayed += 1
        for degree, indices in exact.items():
            if len(indices) <= 1:
                continue
            literal_degree_group_count += 1
            literal_degree_member_count += len(indices)
            # Keep large discrete canonical words only for this one literal
            # degree group.  Persist compact digests/scalars, then release.
            refinement_digest_representatives: defaultdict[str, list[int]] = defaultdict(list)
            for index in indices:
                invariant, data = _refinement(engine, rows[index].colors, degree)
                invariant_digest = hashlib.sha256(invariant).hexdigest()
                # Different digests prove different bytes.  Equal digests are
                # not trusted: replay every earlier representative in that
                # digest bucket and compare the literal canonical words.
                for representative in refinement_digest_representatives[
                    invariant_digest
                ]:
                    previous, _previous_data = _refinement(
                        engine, rows[representative].colors, degree
                    )
                    if invariant == previous:
                        raise ArithmeticError(
                            "canonical refinement left a within-shard collision"
                        )
                    del previous
                refinement_digest_representatives[invariant_digest].append(index)
                refined_current_indices.add(index)
                rows[index].refinement_sha256 = invariant_digest
                rows[index].refinement_discrete = bool(
                    data["color_refinement_discrete_canonical_form"]
                )
                rows[index].refinement_final_cells = int(
                    data["color_refinement_final_cell_count"]
                )
                refinement_discrete += int(rows[index].refinement_discrete)
                del invariant
            del refinement_digest_representatives
        del exact
        _rss_guard()
        if verbose and (
            group_number % 100 == 0
            or group_number == len(digest_collision_groups)
        ):
            print(
                f"collision_progress={group_number}/{len(digest_collision_groups)} "
                f"literal_groups={literal_degree_group_count} "
                f"literal_members={literal_degree_member_count} "
                f"elapsed_seconds={time.perf_counter()-collision_started:.1f} "
                f"peak_rss_gib={_rss_bytes()/(1<<30):.3f}",
                flush=True,
            )

    current_digests = set(digest_groups)
    cross_literal_matches: list[
        tuple[str, int, int, tuple[int, int, int, int], str]
    ] = []
    prior_same_facet_count = 0
    prior_hash_candidates = 0
    prior_literal_replayed = 0
    p5_engine = engine
    star_engine = star.ExactStarBoundaryEngine(gram, roots, interior)
    for source, prior_index, colors, digest, engine_name in _prior_rows(pattern, actions):
        prior_same_facet_count += 1
        if digest not in current_digests:
            continue
        prior_hash_candidates += 1
        prior_engine = p5_engine if engine_name == "P5" else star_engine
        prior_degree, _labels = _degree_replay(prior_engine, colors, digest)
        prior_literal_replayed += 1
        for current_index in digest_groups[digest]:
            current_degree, _labels = _degree_replay(
                engine, rows[current_index].colors, digest
            )
            if current_degree == prior_degree:
                cross_literal_matches.append(
                    (source, prior_index, current_index, colors, engine_name)
                )
            del current_degree
        del prior_degree
        _rss_guard()

    cross_refined_pairs = 0
    for source, prior_index, current_index, prior_colors, engine_name in cross_literal_matches:
        prior_engine = p5_engine if engine_name == "P5" else star_engine
        current_degree, _labels = _degree_replay(
            engine, rows[current_index].colors, rows[current_index].degree_sha256
        )
        prior_degree, _labels = _degree_replay(
            prior_engine, prior_colors, rows[current_index].degree_sha256
        )
        if current_degree != prior_degree:
            raise ArithmeticError("a recorded literal cross-prior degree match changed")
        current_refinement, data = _refinement(
            engine, rows[current_index].colors, current_degree
        )
        current_refinement_digest = hashlib.sha256(current_refinement).hexdigest()
        refined_current_indices.add(current_index)
        rows[current_index].refinement_sha256 = current_refinement_digest
        rows[current_index].refinement_discrete = bool(
            data["color_refinement_discrete_canonical_form"]
        )
        rows[current_index].refinement_final_cells = int(
            data["color_refinement_final_cell_count"]
        )
        del current_refinement
        prior_refinement, _data = _refinement(prior_engine, prior_colors, prior_degree)
        prior_refinement_digest = hashlib.sha256(prior_refinement).hexdigest()
        del prior_refinement
        cross_refined_pairs += 1
        if current_refinement_digest == prior_refinement_digest:
            # Hash equality is not evidence.  Recompute both canonical words
            # and perform the theorem-producing literal comparison.
            current_refinement, _data = _refinement(
                engine, rows[current_index].colors, current_degree
            )
            prior_refinement, _data = _refinement(
                prior_engine, prior_colors, prior_degree
            )
            if current_refinement == prior_refinement:
                raise ArithmeticError(
                    f"cross-prior invariants remain unresolved: {source} row {prior_index} and current row {current_index+1}"
                )
            del current_refinement, prior_refinement
        _rss_guard()

    if verbose:
        print(
            f"cross_prior_rows={prior_same_facet_count} "
            f"digest_candidates={prior_hash_candidates} "
            f"literal_degree_matches={len(cross_literal_matches)} "
            f"refined_pairs={cross_refined_pairs} "
            f"peak_rss_gib={_rss_bytes()/(1<<30):.3f}",
            flush=True,
        )

    rows_transcript_sha256, rows_gzip_sha256 = _write_rows(rows_output, pattern, rows)
    checks = {
        "all_frozen_predecessor_SHA256_values_match": True,
        "exact_action_preserves_stock_Gram_and_AB_split": True,
        "explicit_orbit_cover_and_Burnside_denominator_agree": len(rows) == expected_count,
        "every_boundary_rebuilt_with_ordinary_Coxeter_matrix": True,
        "every_root_primitive_integral_and_reflective": True,
        "every_boundary_has_strict_timelike_interior": True,
        "every_boundary_contains_fixed_Smith_1_power_13_core": True,
        "different_hashes_used_only_as_rigorous_inequality_filters": True,
        "every_same_digest_bucket_replayed_with_literal_degree_bytes": True,
        "literal_degree_collisions_separated_by_exact_canonical_refinement": True,
        "all_rows_pairwise_nonisomorphic_within_shard": True,
        "all_prior_same_facet_hash_candidates_replayed_literally": True,
        "all_literal_cross_prior_matches_separated_by_exact_refinement": True,
        "separated_from_all_previously_completed_P5_rows": True,
        "separated_from_all_same_facet_complete_star_rows": True,
        "Allcock_finite_volume_noncompact_multiplier_five_gate": True,
        "same_full_Picard_weak_K3_carrier_retained": True,
        "RSS_ceiling_is_at_most_1_point_5_GiB": RSS_CEILING_BYTES <= int(1.5 * (1 << 30)),
        "worldwide_priority_explicitly_open": True,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"streaming P5 shard checks failed: {checks}")
    raw_count = p5._raw_pattern_count(pattern)
    if sum(row.orbit_size for row in rows) != raw_count:
        raise ArithmeticError("the compact transcript orbit sizes do not cover the raw shard")

    return {
        "schema": SCHEMA,
        "status_date": STATUS_DATE,
        "all_checks_pass": True,
        "pattern": pattern,
        "pinned_inputs": {**PINNED_SHA256, "builder": _sha256(Path(__file__))},
        "seed": {
            "id": "JK-H12/172",
            "hyperbolic_dimension": 12,
            "facet_count": 172,
            "selected_pairwise_disjoint_doubling_stock": "W61--W172",
            "stock_size": 112,
            "folded_seed_action_order": len(actions),
            "stock_Gram_sha256": stock_gram_sha256,
            "full_Picard_lattice": "U(2)+D5(-1)^2+<-4>",
            "K3_moduli_dimension": 7,
            "source_carrier_theorem": carrier["theorem"],
        },
        "coloring_denominator": {
            "raw_unrooted_proper_colorings": raw_count,
            **burnside,
            "explicit_orbit_cover_count": len(rows),
            "explicit_orbit_size_sum": sum(row.orbit_size for row in rows),
            "all_P5_seed_and_reversal_orbits": denominator["tree_shapes"]["P5"]["seed_and_tree_automorphism_orbits"],
        },
        "theorem": {
            "scope": f"complete {pattern} five-vertex P5 class-word shard",
            "ordinary_Coxeter_type_count": len(rows),
            "pairwise_nonisomorphic_within_shard": True,
            "separated_from_previously_completed_P5_rows": True,
            "separated_from_complete_same_facet_star_rows": True,
            "type_facet_histogram": {str(key): value for key, value in sorted(facet_histogram.items())},
            "raw_unrooted_coloring_facet_histogram": {str(key): value for key, value in sorted(raw_facet_histogram.items())},
            "reflection_index_and_volume_multiple": 5,
            "finite_volume": True,
            "compact": False,
        },
        "compact_exact_transcript": {
            # main() replaces this private-temporary placeholder with the
            # deterministic repository-relative archive path before freezing.
            "rows_path": str(rows_output),
            "row_count": len(rows),
            "uncompressed_canonical_JSONL_sha256": rows_transcript_sha256,
            "deterministic_gzip_sha256": rows_gzip_sha256,
            "exact_gate_record_stream_sha256": gate_stream.hexdigest(),
            "exact_degree_profile_byte_stream_sha256": degree_stream.hexdigest(),
            "hash_role": "negative equality filter and compact linkage only; every equal-digest bucket is literal-byte replayed",
        },
        "ordinary_Coxeter_nonisomorphism_certificate": {
            "degree_digest_class_count": len(digest_groups),
            "degree_digest_collision_group_count": len(digest_collision_groups),
            "degree_digest_collision_member_count": sum(map(len, digest_collision_groups)),
            "literal_degree_replayed_member_count": literal_replayed,
            "literal_degree_collision_group_count": literal_degree_group_count,
            "literal_degree_collision_member_count": literal_degree_member_count,
            "color_refined_current_member_count": len(refined_current_indices),
            "color_refined_discrete_current_count": refinement_discrete,
            "all_literal_degree_collision_groups_refinement_distinct": True,
        },
        "cross_prior_exact_separation": {
            "prior_same_facet_row_count": prior_same_facet_count,
            "archived_degree_digest_candidate_count": prior_hash_candidates,
            "prior_literal_degree_replay_count": prior_literal_replayed,
            "literal_cross_prior_degree_match_count": len(cross_literal_matches),
            "cross_prior_exact_refinement_pair_count": cross_refined_pairs,
            "unresolved_cross_prior_invariant_pair_count": 0,
            "other_previously_completed_P5_rows_separated_by_facet_count": True,
        },
        "aggregate_arithmetic_gates": {
            "ordinary_relation_histogram_inventory": [
                {"histogram": dict(histogram), "type_count": count}
                for histogram, count in sorted(relation_histograms.items())
            ],
            "minimum_strict_interior_facet_pairing": minimum_interior_pairing,
            "maximum_boundary_coordinate_absolute": coordinate_maximum,
            "maximum_transformed_seed_root_coordinate_absolute": transformed_root_maximum,
            "maximum_chamber_transform_coordinate_absolute": transform_maximum,
            "maximum_boundary_pairing_absolute": pairing_maximum,
            "seed_orbit_size_histogram": {str(key): value for key, value in sorted(orbit_size_histogram.items())},
            "nonadjacent_color_repetition_type_histogram": {str(key): value for key, value in sorted(repetition_histogram.items())},
        },
        "resource_certificate": {
            "algorithm": "one boundary at a time; compact row digests; exact replay of digest/collision buckets only",
            "RSS_ceiling_bytes": RSS_CEILING_BYTES,
            "in_process_peak_RSS_guard": True,
            "recommended_outer_watchdog": "python3 scripts/run_with_rss_limit.py --memory-gib 1.0 -- COMMAND",
            "all_Gram_matrices_materialized_together": False,
            "temporary_SQLite_or_other_database_archived": False,
        },
        "Allcock_and_carrier": {
            "finite_subtree_theorem": "each row is a five-copy P5 fundamental polytope",
            "finite_volume": bool(pilot["result"]["finite_volume"]),
            "compact": bool(pilot["result"]["compact"]),
            "reflection_index_and_volume_multiplier": pilot["result"]["reflection_index_and_volume_multiplier"],
            "same_full_Picard_K3_family": True,
            "fixed_index_one_divisor_lattice_core": True,
            "weak_divisor_reflection_carrier_only": True,
            "full_nef_effectivity_or_surface_action_claim": False,
        },
        "lineage_and_claim_boundary": {
            "direction": "SURFACE_TO_POLYTOPE",
            "degree_three_fork_complete_census_claimed": False,
            "all_five_vertex_shapes_complete_census_claimed": False,
            "worldwide_priority": "open",
            "worldwide_novelty_claim": False,
            "bounded_CoxIter_commit": bounded["commit"],
            "weak_carrier_not_full_nef": True,
        },
        "checks": checks,
    }


def _default_paths(pattern: str) -> tuple[Path, Path]:
    stem = f"jacobian_kummer_h12_five_vertex_p5_{pattern.lower()}_streaming_census"
    return (
        ROOT / "data/enumerations" / f"{stem}.json",
        ROOT / "data/enumerations" / f"{stem}_rows.jsonl.gz",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pattern", required=True, choices=SUPPORTED_PATTERNS)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--rows", type=Path)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    default_summary, default_rows = _default_paths(arguments.pattern)
    summary_path = arguments.summary or default_summary
    rows_path = arguments.rows or default_rows

    with tempfile.TemporaryDirectory(prefix=f"jk_h12_{arguments.pattern.lower()}_", dir="/private/tmp") as temporary:
        temporary_rows = Path(temporary) / "rows.jsonl.gz"
        payload = _build(arguments.pattern, rows_output=temporary_rows, verbose=True)
        payload["compact_exact_transcript"]["rows_path"] = str(rows_path.relative_to(ROOT))
        rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
        if arguments.check_archive:
            if not summary_path.is_file() or not rows_path.is_file():
                raise SystemExit("the compact summary or row archive is missing")
            if summary_path.read_text(encoding="utf-8") != rendered:
                raise SystemExit("the compact summary differs from exact replay")
            if rows_path.read_bytes() != temporary_rows.read_bytes():
                raise SystemExit("the deterministic compact rows differ from exact replay")
            print("archive_check=PASS")
        else:
            summary_path.parent.mkdir(parents=True, exist_ok=True)
            rows_path.parent.mkdir(parents=True, exist_ok=True)
            summary_path.write_text(rendered, encoding="utf-8")
            temporary_rows.replace(rows_path)
            print(f"wrote_summary={summary_path}")
            print(f"wrote_rows={rows_path}")
    print("all_checks_pass=true")
    print(f"pattern={arguments.pattern}")
    print(f"type_count={payload['theorem']['ordinary_Coxeter_type_count']}")
    print(f"peak_rss_gib={_rss_bytes()/(1<<30):.3f}")
    print(f"summary_sha256={_sha256(summary_path)}")
    print(f"rows_sha256={_sha256(rows_path)}")


if __name__ == "__main__":
    main()
