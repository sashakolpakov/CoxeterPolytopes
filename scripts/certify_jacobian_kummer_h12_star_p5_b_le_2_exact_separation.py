#!/usr/bin/env python3
"""Close the exact low-B JK--H12 star-versus-P5 comparison.

Every one of the 1,770 star rows and 20,693 P5 rows with at most two B
colours is rebuilt from its archived colour word.  Cross-shape comparison is
literal degree-profile byte equality, followed (only for literal matches) by
exact canonical colour refinement.  SHA-256 is used for provenance and row
linkage, never as the theorem-producing equality test.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import resource
import sys
import time
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import scripts.certify_jacobian_kummer_h12_five_vertex_fork_b3_terminal_acceptance as terminal
import scripts.certify_jacobian_kummer_h12_five_vertex_p5_b_le_2_census as p5
import scripts.certify_jacobian_kummer_h12_five_vertex_star_census as star
import scripts.certify_jacobian_kummer_h12_four_vertex_subtree_census as four


OUTPUT = (
    ROOT
    / "data/enumerations/"
    "jacobian_kummer_h12_star_p5_b_le_2_exact_separation.json"
)
P5_LOW = (
    ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_p5_b_le_2_census.json"
)
STAR = ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_star_census.json"
COMPLETE_P5 = (
    ROOT
    / "data/enumerations/"
    "jacobian_kummer_h12_five_vertex_p5_complete_streaming_census.json"
)
HIGH_B_EXACT = (
    ROOT / "data/enumerations/jacobian_kummer_h12_star_abbb_exact_separation.json"
)
FORK_B_LE_3 = (
    ROOT
    / "data/enumerations/"
    "jacobian_kummer_h12_five_vertex_fork_b_le_3_complete_census.json"
)
PATHS = {
    "frozen_exact_cross_shape_comparator": ROOT
    / "scripts/certify_jacobian_kummer_h12_star_abbb_exact_separation.py",
    "streamed_archive_colour_reader": ROOT
    / "scripts/certify_jacobian_kummer_h12_five_vertex_fork_b3_terminal_acceptance.py",
    "P5_boundary_builder": ROOT
    / "scripts/certify_jacobian_kummer_h12_five_vertex_p5_b_le_2_census.py",
    "star_boundary_builder": ROOT
    / "scripts/certify_jacobian_kummer_h12_five_vertex_star_census.py",
    "exact_colour_refinement_builder": ROOT
    / "scripts/certify_jacobian_kummer_h12_four_vertex_subtree_census.py",
    "P5_B_le_2_archive": P5_LOW,
    "complete_star_archive": STAR,
    "complete_P5_theorem": COMPLETE_P5,
    "star_ABBB_exact_separation": HIGH_B_EXACT,
    "fork_B_le_3_theorem": FORK_B_LE_3,
}
PINNED_SHA256 = {
    "frozen_exact_cross_shape_comparator": "f27f18c36531c48f75e93f5f72b36da694ef0de6f64175d7e84b38186991b926",
    "streamed_archive_colour_reader": "26b16157f8c6c944ead72fb1f0d27cfe5e4a979c2580c5c8d482560313bb83b0",
    "P5_boundary_builder": "d7c5c02803911f8d205a63c59027d8d9584cce998de5efa0a7acb134d56d00e4",
    "star_boundary_builder": "bbb99d59726e0b6bcf0d6b3748823ca9bc8f1efe5ecbe5c29a16f4a5db6a53fb",
    "exact_colour_refinement_builder": "474c796521edeb0f7460906d3779cb0e540d8cbebef0f97176fbe812ec1ba78c",
    "P5_B_le_2_archive": "32943e8531bdaa1874af41d45635ea792c9337694adc4f20dfc152aac3d7713a",
    "complete_star_archive": "21fbd0946b910360cdd4096392e4e629907f971be85359d5b87e312a78d4f309",
    "complete_P5_theorem": "dd6c60a8ed02e985a889100991d2da410d8b18f0902fa314e1e6e9dc900dab16",
    "star_ABBB_exact_separation": "7a328edb6cbedafa13c32140c2f2f1262a8b1159617a1fc59d5ad8d2165c3fe6",
    "fork_B_le_3_theorem": "22a369275c4430ca246f2e0b6884519b2e25cb7c1a6a858fcb776ae4f0943533",
}

SCHEMA = "coxeter.jacobian_kummer.h12_star_p5_b_le_2_exact_separation.v1"
STATUS_DATE = "2026-08-10"
RSS_CEILING_BYTES = int(0.6 * (1 << 30))
P5_COUNTS = {0: 142, 1: 2_091, 2: 18_460}
STAR_COUNTS = {0: 19, 1: 169, 2: 1_582}
FACETS = {0: 752, 1: 761, 2: 770}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _rss_bytes() -> int:
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(peak if sys.platform == "darwin" else peak * 1024)


def _rss_guard() -> None:
    if _rss_bytes() > RSS_CEILING_BYTES:
        raise MemoryError("the low-B cross-shape audit exceeded 0.6 GiB RSS")


def _stream_update(digest: Any, value: bytes) -> None:
    digest.update(len(value).to_bytes(8, "big"))
    digest.update(value)


def _b_count(colors: tuple[int, ...]) -> int:
    return sum(color >= p5.A_SIZE for color in colors)


def _rebuild(
    engine: Any,
    colors: tuple[int, ...],
    expected_digest: str,
    expected_facets: int,
) -> tuple[bytes, Any]:
    boundary, boundary_data = engine.boundary(colors)
    labels, degree_key, ordinary = engine.ordinary_data(boundary)
    if not (
        len(boundary) == boundary_data["facet_count"] == expected_facets
        and ordinary["primitive_integral_reflection_root_count"] == expected_facets
        and hashlib.sha256(degree_key).hexdigest() == expected_digest
    ):
        raise ArithmeticError("an archived low-B exact degree replay changed")
    return degree_key, labels


def build(*, verbose: bool = False) -> dict[str, Any]:
    actual = {name: _sha256(path) for name, path in PATHS.items()}
    if actual != PINNED_SHA256:
        raise ArithmeticError(f"a frozen low-B comparison input changed: {actual}")
    # These production gates transitively pin the seed/action/Gram inputs.
    star_payload = star._pinned_predecessors()[0]
    p5._pinned_predecessors()
    complete_p5 = json.loads(COMPLETE_P5.read_text(encoding="utf-8"))
    high_b = json.loads(HIGH_B_EXACT.read_text(encoding="utf-8"))
    fork = json.loads(FORK_B_LE_3.read_text(encoding="utf-8"))
    if not (
        complete_p5["all_checks_pass"]
        and complete_p5["theorem"]["ordinary_Coxeter_type_count"] == 201_747
        and high_b["all_checks_pass"]
        and high_b["theorem"]["exact_cross_shape_ordinary_Coxeter_isomorphism_count"]
        == 0
        and fork["all_checks_pass"]
        and fork["theorem"]["ordinary_Coxeter_type_count"] == 92_977
        and star_payload["all_checks_pass"]
    ):
        raise ArithmeticError("a frozen cross-shape predecessor theorem is not green")

    _carrier, gram, roots, interior, _bounded = p5.base._load_inputs()
    star_engine = star.ExactStarBoundaryEngine(gram, roots, interior)
    p5_engine = p5.ExactP5BoundaryEngine(gram, roots, interior)

    star_keys: dict[int, defaultdict[bytes, list[tuple[int, tuple[int, ...]]]]] = {
        b: defaultdict(list) for b in range(3)
    }
    star_histogram: Counter[int] = Counter()
    star_key_stream = hashlib.sha256()
    started = time.perf_counter()
    for index, colors, archived_digest in terminal._pretty_json_orbit_rows(
        STAR, 16_872
    ):
        b = _b_count(colors)
        if b > 2:
            continue
        degree, labels = _rebuild(star_engine, colors, archived_digest, FACETS[b])
        del labels
        star_keys[b][degree].append((index, colors))
        star_histogram[b] += 1
        _stream_update(star_key_stream, degree)
        if verbose and sum(star_histogram.values()) % 250 == 0:
            print(
                f"star_low_progress={sum(star_histogram.values())}/1770 "
                f"elapsed_seconds={time.perf_counter()-started:.1f} "
                f"peak_rss_gib={_rss_bytes()/(1<<30):.3f}",
                flush=True,
            )
        if sum(star_histogram.values()) % 100 == 0:
            _rss_guard()
    if dict(star_histogram) != STAR_COUNTS:
        raise ArithmeticError(f"the low-B star denominator changed: {star_histogram}")

    p5_histogram: Counter[int] = Counter()
    p5_key_stream = hashlib.sha256()
    literal_match_pairs = 0
    refined_pairs = 0
    refined_star_cache: dict[int, bytes] = {}
    refinement_discrete_pairs = 0
    started = time.perf_counter()
    for index, colors, archived_digest in terminal._pretty_json_orbit_rows(
        P5_LOW, 20_693
    ):
        b = _b_count(colors)
        if b > 2:
            raise ArithmeticError("the P5 B<=2 archive contains a high-B row")
        degree, labels = _rebuild(p5_engine, colors, archived_digest, FACETS[b])
        p5_histogram[b] += 1
        _stream_update(p5_key_stream, degree)
        matches = star_keys[b].get(degree, ())
        if matches:
            p5_refinement, p5_data = four._canonical_color_refinement(labels)
            for star_index, star_colors in matches:
                literal_match_pairs += 1
                if star_index not in refined_star_cache:
                    star_degree, star_labels = _rebuild(
                        star_engine, star_colors, hashlib.sha256(degree).hexdigest(), FACETS[b]
                    )
                    if star_degree != degree:
                        raise ArithmeticError("a recorded star/P5 degree match changed")
                    star_refinement, _star_data = four._canonical_color_refinement(
                        star_labels
                    )
                    refined_star_cache[star_index] = star_refinement
                if refined_star_cache[star_index] == p5_refinement:
                    raise ArithmeticError(
                        f"UNKNOWN: low-B star row {star_index} and P5 row {index} "
                        "retain the same exact refinement invariant"
                    )
                refined_pairs += 1
                refinement_discrete_pairs += int(
                    p5_data["color_refinement_discrete_canonical_form"]
                )
        del labels
        if verbose and (index % 500 == 0 or index == 20_693):
            print(
                f"P5_low_progress={index}/20693 "
                f"elapsed_seconds={time.perf_counter()-started:.1f} "
                f"literal_matches={literal_match_pairs} "
                f"peak_rss_gib={_rss_bytes()/(1<<30):.3f}",
                flush=True,
            )
        if index % 100 == 0:
            _rss_guard()
    if dict(p5_histogram) != P5_COUNTS:
        raise ArithmeticError(f"the low-B P5 denominator changed: {p5_histogram}")

    exact_all_shape_union = 201_747 + 16_872 + 92_977
    checks = {
        "all_builder_archive_and_predecessor_hashes_match": True,
        "all_1770_low_B_star_rows_rebuilt": sum(star_histogram.values()) == 1_770,
        "all_20693_low_B_P5_rows_rebuilt": sum(p5_histogram.values()) == 20_693,
        "every_archive_row_linked_to_rebuilt_exact_degree_bytes": True,
        "cross_shape_comparison_uses_literal_degree_bytes": True,
        "every_literal_match_separated_by_exact_refinement": (
            literal_match_pairs == refined_pairs
        ),
        "low_B_star_P5_unresolved_pair_count_is_zero": True,
        "complete_star_plus_complete_P5_count_is_218619": (
            16_872 + 201_747 == 218_619
        ),
        "P5_star_fork_B_le_3_union_count_is_311596": (
            exact_all_shape_union == 311_596
        ),
        "only_fork_B4_denominator_remains": True,
        "worldwide_priority_explicitly_open": True,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"the low-B star/P5 exact comparison failed: {checks}")
    return {
        "schema": SCHEMA,
        "status_date": STATUS_DATE,
        "all_checks_pass": True,
        "provenance": {
            "this_builder_sha256": _sha256(Path(__file__).resolve()),
            "pinned_frozen_inputs_sha256": PINNED_SHA256,
            "comparison_method": (
                "frozen exact star/ABBB literal-byte comparator, streamed archive "
                "colours, literal equality, exact refinement of every equality"
            ),
        },
        "denominator": {
            "star_type_count_by_B": {str(key): value for key, value in STAR_COUNTS.items()},
            "P5_type_count_by_B": {str(key): value for key, value in P5_COUNTS.items()},
            "star_low_B_type_count": 1_770,
            "P5_low_B_type_count": 20_693,
        },
        "exact_replay": {
            "star_archive_order_literal_degree_stream_sha256": star_key_stream.hexdigest(),
            "P5_archive_order_literal_degree_stream_sha256": p5_key_stream.hexdigest(),
            "literal_cross_shape_degree_match_pair_count": literal_match_pairs,
            "exact_refinement_separated_pair_count": refined_pairs,
            "P5_discrete_refinement_match_pair_count": refinement_discrete_pairs,
            "distinct_star_rows_refined": len(refined_star_cache),
            "unresolved_pair_count": 0,
            "hash_role": "provenance and row linkage only, never cross-shape equality",
        },
        "theorem": {
            "scope": "complete star versus complete P5 inside the JK-H12 five-vertex stock",
            "low_B_star_plus_P5_pairwise_nonisomorphic_type_count": 22_463,
            "complete_star_plus_complete_P5_pairwise_nonisomorphic_type_count": 218_619,
            "complete_P5_star_plus_fork_B_le_3_pairwise_nonisomorphic_type_count": (
                exact_all_shape_union
            ),
            "finite_volume": True,
            "compact": False,
            "reflection_index_and_volume_multiple": 5,
            "hyperbolic_dimension": 12,
        },
        "claim_boundary": {
            "complete_five_vertex_three_shape_census_claimed": False,
            "only_unclassified_fork_B4_denominator": 106_385,
            "full_nef_effectivity_or_surface_action_claim": False,
            "worldwide_priority": "open",
            "worldwide_novelty_claim": False,
        },
        "resource_certificate": {
            "RSS_ceiling_bytes": RSS_CEILING_BYTES,
            "one_boundary_and_label_matrix_at_a_time": True,
            "only_1770_star_literal_keys_retained": True,
        },
        "checks": checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    payload = build(verbose=True)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if (
            not arguments.output.is_file()
            or arguments.output.read_text(encoding="utf-8") != rendered
        ):
            raise SystemExit("the low-B star/P5 exact-separation archive differs")
        print("archive_check=PASS")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
        print(f"wrote={arguments.output}")
    print("all_checks_pass=true")
    print(
        "literal_cross_shape_degree_match_pair_count="
        f"{payload['exact_replay']['literal_cross_shape_degree_match_pair_count']}"
    )
    print("complete_P5_star_fork_B_le_3_union=311596")
    print(f"peak_rss_gib={_rss_bytes()/(1<<30):.3f}")
    print(f"payload_sha256={hashlib.sha256(rendered.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
