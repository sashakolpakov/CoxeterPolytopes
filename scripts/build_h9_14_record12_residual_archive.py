#!/usr/bin/env python3
"""Freeze the corrected exhaustive record-12 abstract residual census.

The staged streams implement the sound order of operations forced by the
withdrawal of the provisional 344-orbit census: enumerate every locally
ridge-feasible branch, apply only the proved special-meeting blocker lemma,
couple both branches with both cut orientations, form the full blocker, and
only then quotient the surviving global complexes.  The output is an exact
abstract incidence census, not a PL, convex-polytopal, or Coxeter result.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record12_residual import (
    EXPECTED_BALANCED_BRANCH_TYPE_COUNT,
    EXPECTED_COLOURED_BRANCH_ORBIT_COUNT,
    EXPECTED_COLOURED_BRANCH_STREAM_SHA256,
    EXPECTED_COUPLED_PAIR_STREAM_SHA256,
    EXPECTED_FULL_BLOCKER_PASS_COUNT,
    EXPECTED_FULL_ORBIT_COUNT,
    EXPECTED_FULL_ORBIT_STREAM_SHA256,
    EXPECTED_LOCAL_PAIR_FEASIBLE_STATE_COUNT,
    EXPECTED_NORMALIZED_CUT_STATE_COUNT,
    EXPECTED_PROFILE_FACET_TYPE_PAIR_COUNT,
    EXPECTED_RELATIVE_EMBEDDING_COUNT,
    EXPECTED_SAFE_SPECIAL_BRANCH_STATE_COUNT,
    EXPECTED_SAFE_SPECIAL_STREAM_SHA256,
    EXPECTED_SCHEDULE_ORBIT_COUNTS,
    EXPECTED_SURVIVING_ORDERED_TYPE_COUNT,
    EXPECTED_UNFILTERED_LOCAL_STREAM_SHA256,
    EXPECTED_X_GRAPH_COUNT,
    TARGET_EDGE_MASKS,
    adjacency_from_graph6,
    branch_has_local_pair_extension,
    branch_triples,
    canonical_digest,
    incidence_audit,
    minimal_transversals,
    pair_mask_from_adjacency,
    pair_mask_to_matrix,
    switched_triangle_free,
    verify_archive,
)


ROOT = Path(__file__).resolve().parents[1]
X_ARCHIVE = ROOT / "data/enumerations/h9_14_record12_x_graphs.json"
OUTPUT = ROOT / "data/enumerations/h9_14_record12_residual.json"
WITHDRAWAL_REPORT = ROOT / "results/h9_14_record12_residual_hostile_audit.md"
HOSTILE_REPORT = ROOT / "results/h9_14_record12_corrected_371_hostile_audit.md"

DEFAULT_UNFILTERED_STREAM = Path("/tmp/record12_local_all.txt")
DEFAULT_SAFE_STREAM = Path("/tmp/record12_branch_safe_special_correct.txt")
DEFAULT_COLOURED_STREAM = Path(
    "/tmp/record12_branch_safe_special_correct_orbits.txt"
)
DEFAULT_PAIR_STREAM = Path("/tmp/record12_safe_pairs_all_sorted.txt")
DEFAULT_ORBIT_STREAM = Path("/tmp/record12_safe_full_orbits_v2.txt")

SPECIALS = (1 << 10) | (1 << 11)
FULL_LOCAL = (1 << 10) - 1


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _line_count(path: Path) -> int:
    with path.open("rb") as handle:
        return sum(1 for _ in handle)


def _check_stream(path: Path, count: int, digest: str) -> None:
    if _line_count(path) != count or _file_sha256(path) != digest:
        raise AssertionError(f"record-12 staged stream changed: {path}")


def _histogram(rows: Iterable[int]) -> dict[str, int]:
    return {
        str(size): count
        for size, count in sorted(Counter(row.bit_count() for row in rows).items())
    }


def _safe_special_rows(rows: Sequence[int]) -> bool:
    return all(not row & SPECIALS or row.bit_count() in (4, 5) for row in rows)


def _verify_safe_raw(path: Path) -> None:
    seen = set()
    for line in path.read_text(encoding="ascii").splitlines():
        word, raw_p, raw_y = line.split()
        p_mask, y_mask = int(raw_p), int(raw_y)
        key = (word, p_mask, y_mask)
        if key in seen:
            raise AssertionError("record-12 safe branch stream repeats a state")
        seen.add(key)
        x_mask = pair_mask_from_adjacency(adjacency_from_graph6(word))
        if (
            not switched_triangle_free(x_mask, 0)
            or not switched_triangle_free(x_mask, p_mask)
            or not branch_has_local_pair_extension(x_mask, p_mask, y_mask)
        ):
            raise AssertionError("record-12 safe branch failed a ridge condition")
        blocker = minimal_transversals(branch_triples(x_mask, p_mask, y_mask))
        if not _safe_special_rows(blocker):
            raise AssertionError("record-12 safe branch violates the proved filter")


def _branch_records(path: Path) -> list[dict[str, object]]:
    records = []
    seen = set()
    multiplicity_sum = 0
    for index, line in enumerate(path.read_text(encoding="ascii").splitlines()):
        word, raw_p, raw_y, raw_multiplicity = line.split()
        p_mask, y_mask = int(raw_p), int(raw_y)
        multiplicity = int(raw_multiplicity)
        colour_values = tuple(
            2 * int(bool(p_mask & (1 << vertex)))
            + int(bool(y_mask & (1 << vertex)))
            for vertex in range(10)
        )
        if colour_values != tuple(sorted(colour_values)):
            raise AssertionError("record-12 coloured representative is not normalized")
        key = (word, p_mask, y_mask)
        if key in seen or multiplicity <= 0:
            raise AssertionError("record-12 coloured branch row is malformed")
        seen.add(key)
        x_mask = pair_mask_from_adjacency(adjacency_from_graph6(word))
        if (
            not switched_triangle_free(x_mask, 0)
            or not switched_triangle_free(x_mask, p_mask)
            or not branch_has_local_pair_extension(x_mask, p_mask, y_mask)
        ):
            raise AssertionError("record-12 coloured branch failed a ridge condition")
        triples = branch_triples(x_mask, p_mask, y_mask)
        blocker = minimal_transversals(triples)
        if not _safe_special_rows(blocker):
            raise AssertionError("record-12 coloured branch violates the safe lemma")
        special_rows = tuple(row for row in blocker if row & SPECIALS)
        records.append(
            {
                "index": index,
                "x_graph6": word,
                "p_mask": p_mask,
                "y_mask": y_mask,
                "colour_partition": "".join(
                    chr(ord("a") + value) for value in colour_values
                ),
                "raw_representation_count": multiplicity,
                "branch_triple_count": len(triples),
                "branch_blocker_histogram": _histogram(blocker),
                "special_meeting_blocker_histogram": _histogram(special_rows),
            }
        )
        multiplicity_sum += multiplicity
    if (
        len(records) != EXPECTED_COLOURED_BRANCH_ORBIT_COUNT
        or multiplicity_sum != EXPECTED_SAFE_SPECIAL_BRANCH_STATE_COUNT
    ):
        raise AssertionError("record-12 coloured branch quotient changed")
    return records


def _minimalize(rows: Iterable[int]) -> tuple[int, ...]:
    result = []
    for row in sorted(set(rows), key=lambda value: (value.bit_count(), value)):
        if not any(prior & row == prior for prior in result):
            result.append(row)
    return tuple(result)


def _global_branch_blocker(
    local_rows: Sequence[int], *, p_branch: bool
) -> tuple[int, ...]:
    if p_branch:
        result = [1 << 0, 1 << 1]
        for row in local_rows:
            global_row = (row & FULL_LOCAL) << 2
            if row & (1 << 10):
                global_row |= 1 << 12
            if row & (1 << 11):
                global_row |= 1 << 13
            result.append(global_row)
    else:
        result = [1 << 12, 1 << 13]
        for row in local_rows:
            global_row = (row & FULL_LOCAL) << 2
            if row & (1 << 10):
                global_row |= 1 << 0
            if row & (1 << 11):
                global_row |= 1 << 1
            result.append(global_row)
    return tuple(result)


def _full_blocker(
    x_mask: int, z_mask: int, p_mask: int, y_mask: int,
    cache: dict[tuple[int, int, int], tuple[int, ...]],
) -> tuple[int, ...]:
    p_key = (x_mask, p_mask, y_mask)
    q_key = (z_mask, p_mask ^ y_mask, y_mask)
    if p_key not in cache:
        cache[p_key] = minimal_transversals(
            branch_triples(x_mask, p_mask, y_mask)
        )
    if q_key not in cache:
        cache[q_key] = minimal_transversals(
            branch_triples(z_mask, p_mask ^ y_mask, y_mask)
        )
    p_rows = _global_branch_blocker(cache[p_key], p_branch=True)
    q_rows = _global_branch_blocker(cache[q_key], p_branch=False)
    return _minimalize(first | second for first in p_rows for second in q_rows)


def _pair_summary(
    path: Path, branches: Sequence[Mapping[str, object]]
) -> dict[str, object]:
    blocker_cache: dict[tuple[int, int, int], tuple[int, ...]] = {}
    triples_cache: dict[tuple[int, int, int], int] = {}
    states = set()
    type_pairs = set()
    p_indices = set()
    q_indices = set()
    previous = None
    for line in path.read_text(encoding="ascii").splitlines():
        if previous is not None and line < previous:
            raise AssertionError("record-12 coupled-pair stream is not canonical-sorted")
        previous = line
        fields = line.split()
        if len(fields) != 10 or fields[0] != "PAIR":
            raise ValueError("malformed record-12 coupled-pair row")
        p_index, q_index, x_mask, z_mask, p_mask, y_mask, facets, a, b = map(
            int, fields[1:]
        )
        if not (0 <= p_index < len(branches) and 0 <= q_index < len(branches)):
            raise AssertionError("record-12 pair branch index is out of range")
        p_state, q_state = branches[p_index], branches[q_index]
        observed_x = pair_mask_from_adjacency(
            adjacency_from_graph6(str(p_state["x_graph6"]))
        )
        if (
            observed_x != x_mask
            or int(p_state["p_mask"]) != p_mask
            or int(p_state["y_mask"]) != y_mask
        ):
            raise AssertionError("record-12 pair P-state payload changed")
        # Every globally surviving pair is balanced.  Matrix conversion also
        # rejects any hidden within-side edge.
        pair_mask_to_matrix(x_mask)
        pair_mask_to_matrix(z_mask)
        if p_mask != 0b1111100000 or y_mask != 0:
            raise AssertionError("record-12 global survivor left the balanced branch")
        q_source_mask = pair_mask_from_adjacency(
            adjacency_from_graph6(str(q_state["x_graph6"]))
        )
        p_key = (x_mask, p_mask, y_mask)
        q_source_key = (
            q_source_mask,
            int(q_state["p_mask"]),
            int(q_state["y_mask"]),
        )
        if p_key not in triples_cache:
            triples_cache[p_key] = len(branch_triples(*p_key))
        if q_source_key not in triples_cache:
            triples_cache[q_source_key] = len(branch_triples(*q_source_key))
        if facets != triples_cache[p_key] + triples_cache[q_source_key] - y_mask.bit_count():
            raise AssertionError("record-12 pair facet count changed")
        blocker = _full_blocker(x_mask, z_mask, p_mask, y_mask, blocker_cache)
        histogram = Counter(row.bit_count() for row in blocker)
        twos = {row for row in blocker if row.bit_count() == 2}
        if (
            set(histogram) - {2, 4, 5}
            or twos != set(TARGET_EDGE_MASKS)
            or histogram[4] != a
            or histogram[5] != b
            or facets != 130 - 2 * a
        ):
            raise AssertionError("record-12 coupled full blocker changed")
        state_key = (x_mask, z_mask, p_mask, y_mask)
        if state_key in states:
            raise AssertionError("record-12 coupled-pair stream repeats an embedding")
        states.add(state_key)
        type_pairs.add((p_index, q_index))
        p_indices.add(p_index)
        q_indices.add(q_index)
    if (
        len(states) != EXPECTED_FULL_BLOCKER_PASS_COUNT
        or len(type_pairs) != EXPECTED_SURVIVING_ORDERED_TYPE_COUNT
        or len(p_indices | q_indices) != EXPECTED_BALANCED_BRANCH_TYPE_COUNT
    ):
        raise AssertionError("record-12 globally surviving pair census changed")
    return {
        "surviving_ordered_type_pairs": [list(pair) for pair in sorted(type_pairs)],
        "surviving_branch_type_indices": sorted(p_indices | q_indices),
        "surviving_branch_type_count": len(p_indices | q_indices),
        "full_blocker_passing_relative_embeddings": len(states),
    }


def _residual_orbits(path: Path) -> list[dict[str, object]]:
    records = []
    seen = set()
    multiplicity_sum = 0
    for line in path.read_text(encoding="ascii").splitlines():
        marker, raw_x, raw_z, raw_multiplicity = line.split()
        if marker != "ORBIT":
            raise ValueError("malformed record-12 full-orbit row")
        x_matrix, z_matrix = int(raw_x), int(raw_z)
        multiplicity = int(raw_multiplicity)
        key = (x_matrix, z_matrix)
        if key in seen or multiplicity <= 0:
            raise AssertionError("record-12 full-orbit representative repeats")
        seen.add(key)
        audit = incidence_audit(x_matrix, z_matrix)
        if not audit["verified"]:
            raise AssertionError("record-12 residual failed its incidence replay")
        records.append(
            {
                "x_matrix": x_matrix,
                "z_matrix": z_matrix,
                "relative_embedding_multiplicity": multiplicity,
                "missing_four_count": audit["missing_four_count"],
                "facet_count": audit["facet_count"],
                "minimal_nonface_histogram": {
                    str(key): value
                    for key, value in audit["minimal_nonface_histogram"].items()
                },
                "face_vector": list(audit["face_vector"]),
                "selected_complement_sha256": audit["selected_complement_sha256"],
            }
        )
        multiplicity_sum += multiplicity
    if len(records) != EXPECTED_FULL_ORBIT_COUNT or multiplicity_sum != EXPECTED_FULL_BLOCKER_PASS_COUNT:
        raise AssertionError("record-12 full automorphism quotient changed")
    return records


def build(
    unfiltered_stream: Path = DEFAULT_UNFILTERED_STREAM,
    safe_stream: Path = DEFAULT_SAFE_STREAM,
    coloured_stream: Path = DEFAULT_COLOURED_STREAM,
    pair_stream: Path = DEFAULT_PAIR_STREAM,
    orbit_stream: Path = DEFAULT_ORBIT_STREAM,
) -> dict[str, object]:
    _check_stream(
        unfiltered_stream,
        EXPECTED_LOCAL_PAIR_FEASIBLE_STATE_COUNT,
        EXPECTED_UNFILTERED_LOCAL_STREAM_SHA256,
    )
    _check_stream(
        safe_stream,
        EXPECTED_SAFE_SPECIAL_BRANCH_STATE_COUNT,
        EXPECTED_SAFE_SPECIAL_STREAM_SHA256,
    )
    _check_stream(
        coloured_stream,
        EXPECTED_COLOURED_BRANCH_ORBIT_COUNT,
        EXPECTED_COLOURED_BRANCH_STREAM_SHA256,
    )
    _check_stream(
        pair_stream,
        EXPECTED_FULL_BLOCKER_PASS_COUNT,
        EXPECTED_COUPLED_PAIR_STREAM_SHA256,
    )
    _check_stream(
        orbit_stream,
        EXPECTED_FULL_ORBIT_COUNT,
        EXPECTED_FULL_ORBIT_STREAM_SHA256,
    )
    _verify_safe_raw(safe_stream)
    branches = _branch_records(coloured_stream)
    pairing = _pair_summary(pair_stream, branches)
    residuals = _residual_orbits(orbit_stream)
    schedule = Counter(int(row["missing_four_count"]) for row in residuals)
    if dict(sorted(schedule.items())) != EXPECTED_SCHEDULE_ORBIT_COUNTS:
        raise AssertionError("record-12 corrected schedule histogram changed")

    x_document = json.loads(X_ARCHIVE.read_text(encoding="utf-8"))
    document: dict[str, object] = {
        "schema_version": "2.0.0",
        "kind": "h9_14_record12_exact_371_orbit_whole_shadow_residual",
        "record_index": 12,
        "graph6": "M???????????o?o??",
        "theorem_status": "certified_exact_371_orbit_abstract_residual",
        "claim": (
            "Under the published compact-H9/14 missing-face reductions, the "
            "entire record-12 graph shadow has exactly 371 abstract residual "
            "incidence-complex orbits. Their a=0,...,4 histogram is "
            "48,95,144,63,21. Every orbit passes the archived incidence, "
            "local g/Macaulay, Euler-link, and mod-two homology-link audits."
        ),
        "safe_special_filter_lemma": {
            "statement": (
                "If r is a minimal transversal of one branch triple clutter "
                "and r meets a special of that branch, then r is a minimal "
                "transversal of the full coupled facet-complement clutter."
            ),
            "proof": (
                "The row r hits its own branch by definition and the opposite "
                "branch because every opposite complement contains both branch "
                "specials. For each v in r, branch minimality supplies a branch "
                "triple disjoint from r without v; adjoining the opposite "
                "two-set remains disjoint because r contains none of those "
                "vertices. Thus every deletion of v fails globally."
            ),
            "consequence": (
                "Only special-meeting branch blocker rows are filtered to "
                "orders four and five. Isolate-only rows are never filtered "
                "before the two branches are coupled."
            ),
        },
        "x_graph_archive": {
            "path": str(X_ARCHIVE.relative_to(ROOT)),
            "certificate_sha256": x_document["certificate_sha256"],
            "raw_stream_sha256": x_document["raw_stream_sha256"],
            "graph_count": x_document["graph_count"],
        },
        "generation_streams": {
            "unfiltered_local_ridge_states": {
                "command": (
                    "geng -q -t 10 | /tmp/record12_branch_safe -1 -2 > "
                    "/tmp/record12_local_all.txt"
                ),
                "generator_source": "scripts/explore_h9_14_record12_branch.cpp",
                "line_count": EXPECTED_LOCAL_PAIR_FEASIBLE_STATE_COUNT,
                "raw_stream_sha256": EXPECTED_UNFILTERED_LOCAL_STREAM_SHA256,
            },
            "safe_special_branch_states": {
                "command": (
                    "geng -q -t 10 | /tmp/record12_branch_safe -1 -3 > "
                    "/tmp/record12_branch_safe_special_correct.txt"
                ),
                "generator_source": "scripts/explore_h9_14_record12_branch.cpp",
                "line_count": EXPECTED_SAFE_SPECIAL_BRANCH_STATE_COUNT,
                "raw_stream_sha256": EXPECTED_SAFE_SPECIAL_STREAM_SHA256,
            },
            "coloured_branch_orbits": {
                "command": (
                    "PYTHONPATH=src:. python3 -u "
                    "scripts/explore_h9_14_record12_unfiltered_orbits.py "
                    "/tmp/record12_branch_safe_special_correct.txt "
                    "/tmp/record12_branch_safe_special_correct_orbits.txt"
                ),
                "generator_source": "scripts/explore_h9_14_record12_unfiltered_orbits.py",
                "line_count": EXPECTED_COLOURED_BRANCH_ORBIT_COUNT,
                "raw_stream_sha256": EXPECTED_COLOURED_BRANCH_STREAM_SHA256,
            },
            "coupled_full_blocker_passes": {
                "command": (
                    "/tmp/record12_safe_pairs < "
                    "/tmp/record12_branch_safe_special_correct_orbits.txt | "
                    "LC_ALL=C sort > /tmp/record12_safe_pairs_all_sorted.txt"
                ),
                "generator_source": "scripts/explore_h9_14_record12_safe_pairs.cpp",
                "line_count": EXPECTED_FULL_BLOCKER_PASS_COUNT,
                "raw_stream_sha256": EXPECTED_COUPLED_PAIR_STREAM_SHA256,
            },
            "full_automorphism_orbits": {
                "command": (
                    "awk '{print \"BALANCED_PAIR\",$4,$5}' "
                    "/tmp/record12_safe_pairs_all_sorted.txt | "
                    "/tmp/record12_pair_orbits > "
                    "/tmp/record12_safe_full_orbits_v2.txt"
                ),
                "generator_source": "scripts/explore_h9_14_record12_pair_orbits.cpp",
                "line_count": EXPECTED_FULL_ORBIT_COUNT,
                "raw_stream_sha256": EXPECTED_FULL_ORBIT_STREAM_SHA256,
            },
        },
        "census": {
            "x_graph_count": EXPECTED_X_GRAPH_COUNT,
            "normalized_cut_states": EXPECTED_NORMALIZED_CUT_STATE_COUNT,
            "locally_pair_feasible_states": EXPECTED_LOCAL_PAIR_FEASIBLE_STATE_COUNT,
            "safe_special_branch_states": EXPECTED_SAFE_SPECIAL_BRANCH_STATE_COUNT,
            "coloured_branch_orbits": EXPECTED_COLOURED_BRANCH_ORBIT_COUNT,
            "profile_and_facet_compatible_ordered_type_pairs": (
                EXPECTED_PROFILE_FACET_TYPE_PAIR_COUNT
            ),
            "relative_embedded_graphs": EXPECTED_RELATIVE_EMBEDDING_COUNT,
            "full_blocker_passing_relative_embeddings": (
                EXPECTED_FULL_BLOCKER_PASS_COUNT
            ),
            "surviving_ordered_type_pairs": EXPECTED_SURVIVING_ORDERED_TYPE_COUNT,
            "surviving_balanced_branch_types": EXPECTED_BALANCED_BRANCH_TYPE_COUNT,
            "full_graph_shadow_orbits": EXPECTED_FULL_ORBIT_COUNT,
        },
        "pairing_audit": pairing,
        "cut_normalization_audit": {
            "both_q_and_complement_q_enumerated": True,
            "direct_profile_pairs": 13_516,
            "complement_profile_pairs": 9_598,
            "both_profile_pairs": 256,
            "direct_images": 1_381_203,
            "complement_images": 1_376_516,
            "image_overlap": 32_000,
            "union_relative_embeddings": EXPECTED_RELATIVE_EMBEDDING_COUNT,
            "correction": (
                "The complementary q-cut acts as a relative transpose of the "
                "second 5x5 branch matrix. Adding these orientations enlarges "
                "the provisional 344-orbit quotient to 371."
            ),
        },
        "automorphism_quotient": {
            "missing_edge_graph_automorphism_group": "D8 x S10",
            "balanced_cut_stabilizer": "(S5 wr C2) x D8",
            "input_relative_embeddings": EXPECTED_FULL_BLOCKER_PASS_COUNT,
            "full_orbits": EXPECTED_FULL_ORBIT_COUNT,
            "actions": [
                "row and column permutations of the 5x5 branch matrices",
                "simultaneous transpose induced by swapping the cut sides",
                "independent bipartite complementation of X and Z",
                "interchange of X and Z induced by swapping the K2,2 parts",
            ],
        },
        "schedule_orbit_counts": {
            str(key): value for key, value in sorted(schedule.items())
        },
        "coloured_branch_orbits": branches,
        "residual_orbits": residuals,
        "correction_history": {
            "withdrawn_restricted_count": 344,
            "corrected_whole_shadow_count": 371,
            "new_relative_transpose_orbits": 27,
            "new_orbit_schedule_counts": {"2": 15, "3": 9, "4": 3},
            "withdrawal_report": {
                "path": str(WITHDRAWAL_REPORT.relative_to(ROOT)),
                "sha256": _file_sha256(WITHDRAWAL_REPORT),
            },
            "corrected_hostile_report": {
                "path": str(HOSTILE_REPORT.relative_to(ROOT)),
                "sha256": _file_sha256(HOSTILE_REPORT),
            },
        },
        "relation_to_known_example": (
            "A separate exact rational certificate supplies one ordinary-"
            "polytopal a=2 member of the 371-orbit set. No other orbit is "
            "asserted polytopal here."
        ),
        "nonclaims": [
            "The 371 residuals are abstract incidence and mod-two homology-sphere candidates, not certified PL spheres.",
            "Except for the separately certified single rational witness, ordinary convex-polytopal realizability is not claimed.",
            "No Coxeter labels or Lorentzian Gram matrices are supplied.",
            "No compact or finite-volume hyperbolic Coxeter polytope is claimed.",
        ],
    }
    document["certificate_sha256"] = canonical_digest(document)
    verify_archive(document)
    return document


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--unfiltered-stream", type=Path, default=DEFAULT_UNFILTERED_STREAM)
    parser.add_argument("--safe-stream", type=Path, default=DEFAULT_SAFE_STREAM)
    parser.add_argument("--coloured-stream", type=Path, default=DEFAULT_COLOURED_STREAM)
    parser.add_argument("--pair-stream", type=Path, default=DEFAULT_PAIR_STREAM)
    parser.add_argument("--orbit-stream", type=Path, default=DEFAULT_ORBIT_STREAM)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args(argv)
    if args.verify_only:
        document = json.loads(OUTPUT.read_text(encoding="utf-8"))
        summary = verify_archive(document)
    else:
        document = build(
            args.unfiltered_stream,
            args.safe_stream,
            args.coloured_stream,
            args.pair_stream,
            args.orbit_stream,
        )
        OUTPUT.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        summary = verify_archive(document)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
