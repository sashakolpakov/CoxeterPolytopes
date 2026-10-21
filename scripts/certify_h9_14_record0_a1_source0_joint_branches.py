#!/usr/bin/env python3
"""Freeze and replay the nine-way compact label exhaustion for Record 0/a=1/source 0.

The expensive CP-SAT searches are represented by nine complete terminal
documents.  This verifier independently reconstructs the 108 classical
order-four Lanner placements, their full-S4 orbit partition, the labelled-A4
stabilizers, the joint-model counts, and the deletion-minimal nonfinite
pattern expansion.  Only all nine terminal INFEASIBLE documents imply the
fixed-source compact no-go.  The finite-volume noncompact problem is outside
the certificate.
"""

from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
import json
from math import factorial, isfinite
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from scripts.certify_h9_14_record0_a1_source1_joint_branches import (
    _check_minimal_nonfinite,
    _code_orbit,
    _digest,
    _file_sha256,
    _lanner4_placements,
    _pair,
    _permuted_code,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data/enumerations/h9_14_record0_a1_source0_joint_branch_certificate.json"
)

N = 14
GRAPH6 = "M?????????????o??"
MISSING_EDGES = (
    sum(1 << value for value in (0, 13)),
    sum(1 << value for value in (1, 13)),
)
FIXED_A4_VERTICES = (2, 3, 4, 5)
FIXED_A4 = sum(1 << value for value in FIXED_A4_VERTICES)
EXPECTED_BRANCH_INDICES = tuple(range(9))
EXPECTED_REPRESENTATIVES = (
    (2, 2, 3, 2, 3, 5),
    (2, 2, 3, 3, 2, 5),
    (2, 2, 4, 5, 2, 3),
    (2, 2, 5, 5, 2, 3),
    (2, 3, 3, 3, 4, 2),
    (2, 3, 3, 3, 5, 2),
    (2, 3, 4, 4, 3, 2),
    (2, 3, 4, 5, 3, 2),
    (2, 3, 5, 5, 3, 2),
)
EXPECTED_ORBIT_SIZES = (12, 12, 24, 12, 12, 12, 6, 12, 6)
EXPECTED_PARTITION_DIGEST = (
    "bc69a04c15d0d2126e687e28e5134f232fd2f65c1dc721285840a93c1c5830cc"
)
EXPECTED_PATTERN_DIGEST = (
    "4004924328df160222cfc792a9cd51fb305bae053178914db5d1dabedc72fba0"
)
EXPECTED_CUT_COUNT = 237_543
EXPECTED_CUT_DIGEST = (
    "033efa726adbea0daa6d7963fd1547653a5ab1c9c887fd064755fe774b8ce1ec"
)
EXPECTED_SEED_TRANSCRIPT = [
    {
        "expanded_global_face_cut_count": EXPECTED_CUT_COUNT,
        "expanded_global_face_cut_digest": EXPECTED_CUT_DIGEST,
        "new_abstract_pattern_count": 10,
        "seed_round": 0,
    }
]
EXPECTED_SEED_TRANSCRIPT_DIGEST = (
    "0d8ea6f3be259614a075f955caece5cc0bae295e1b6c74653685e264db59d57b"
)
EXPECTED_MODEL_COUNTS = {
    "candidate_facet_complements": 719,
    "conditional_minimal_A5_compact_lanner_encodings": 1607,
    "face_witness_constraints": 1324,
    "global_meeting_label_variables": 89,
    "independent_five_sets": 1617,
    "lanner_order_constraints": 6669,
    "ridge_constraints": 1456,
}
EXPECTED_CERTIFICATE_SHA256 = (
    "fc8117901ca897ee21bc1a7becd7fcc505b4218ceb8f6a4eae55005ab58c5226"
)

SOURCE_CERTIFICATES = (
    {
        "source_orbit_index": 1,
        "path": "data/enumerations/h9_14_record0_a1_source1_joint_branch_certificate.json",
        "file_sha256": "77b9c2808073fcd8b7f8d63724f0c70952f7317ca139700c7b9f7c8899b0fbd9",
        "boolean_path": ("coverage", "fixed_source_compact_coxeter_no_go"),
    },
    {
        "source_orbit_index": 2,
        "path": "data/enumerations/h9_14_record0_a1_source2_incidence_no_go.json",
        "file_sha256": "e9f6d4c5c607b2fecc6633e9c2c3aac14ac174aa467295348dfd04bf481626af",
        "boolean_path": ("theorem", "fixed_source_compact_incidence_no_go"),
    },
    {
        "source_orbit_index": 3,
        "path": "data/enumerations/h9_14_record0_a1_source3_terminal_cegar_audit.json",
        "file_sha256": "afe350d601d76707026475c328479c99c4080809e4665a30b585b696fa61f1c4",
        "boolean_path": ("decision", "fixed_source_compact_coxeter_no_go"),
    },
)

BRANCH_SCOPE = (
    "Complete joint compact incidence and finite-label formulation for one "
    "fixed a=1 source/A4-label branch. Only terminal INFEASIBLE closes the "
    "branch; all placement branches must close to infer a fixed-source theorem. "
    "An all-elliptic label model still requires exact Gram."
)
LOSSLESS_REASON = (
    "Move the A4 placement to its unique placement-orbit representative, then "
    "choose a lexicographically least incidence vector under the exact-code "
    "stabilizer."
)


def _mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def _contains_base(value: int) -> bool:
    return any(value & row == row for row in (*MISSING_EDGES, FIXED_A4))


@lru_cache(maxsize=None)
def placement_orbits() -> tuple[tuple[tuple[int, ...], ...], ...]:
    actions = tuple(permutations(range(4)))
    remaining = set(_lanner4_placements())
    answer = []
    while remaining:
        representative = min(remaining)
        orbit = tuple(
            sorted(
                {
                    _permuted_code(4, representative, action)
                    for action in actions
                }
            )
        )
        if not set(orbit) <= remaining:
            raise AssertionError("full-S4 placement orbits overlap")
        remaining.difference_update(orbit)
        answer.append(orbit)
    result = tuple(answer)
    if tuple(row[0] for row in result) != EXPECTED_REPRESENTATIVES:
        raise AssertionError("source-0 placement representatives changed")
    if tuple(map(len, result)) != EXPECTED_ORBIT_SIZES:
        raise AssertionError("source-0 placement orbit sizes changed")
    if sum(map(len, result)) != 108 or _digest(result) != EXPECTED_PARTITION_DIGEST:
        raise AssertionError("source-0 placement partition changed")
    return result


def partition_record() -> dict[str, object]:
    orbits = placement_orbits()
    return {
        "literal_placement_count": 108,
        "orbit_count": 9,
        "orbit_sizes": list(map(len, orbits)),
        "orbit_representatives": [list(row[0]) for row in orbits],
        "partition_digest": _digest(orbits),
        "coverage_is_disjoint_and_complete": (
            len({code for orbit in orbits for code in orbit}) == 108
            and {code for orbit in orbits for code in orbit}
            == set(_lanner4_placements())
        ),
    }


@lru_cache(maxsize=None)
def _candidate_facets() -> tuple[int, ...]:
    independent_fives = tuple(
        _mask(chosen)
        for chosen in combinations(range(N), 5)
        if all(_mask(chosen) & edge != edge for edge in MISSING_EDGES)
    )
    shadow = tuple(value for value in independent_fives if value & FIXED_A4 == FIXED_A4)
    base = (*MISSING_EDGES, FIXED_A4, *shadow)
    return tuple(
        value
        for value in (_mask(chosen) for chosen in combinations(range(N), 5))
        if all(value & row for row in base)
    )


@lru_cache(maxsize=None)
def independent_model_counts() -> dict[str, int]:
    independent_fives = tuple(
        value
        for value in (_mask(chosen) for chosen in combinations(range(N), 5))
        if all(value & edge != edge for edge in MISSING_EDGES)
    )
    facets = _candidate_facets()
    counts = {
        "candidate_facet_complements": len(facets),
        "conditional_minimal_A5_compact_lanner_encodings": sum(
            not _contains_base(value) for value in independent_fives
        ),
        "face_witness_constraints": sum(
            not _contains_base(_mask(chosen))
            for order in range(5)
            for chosen in combinations(range(N), order)
        ),
        "global_meeting_label_variables": sum(
            _mask(edge) not in MISSING_EDGES for edge in combinations(range(N), 2)
        ),
        "independent_five_sets": len(independent_fives),
        "lanner_order_constraints": sum(
            not _contains_base(_mask(chosen))
            for order in range(6, N + 1)
            for chosen in combinations(range(N), order)
        ),
        "ridge_constraints": sum(
            any(facet & six == facet for facet in facets)
            for six in (_mask(chosen) for chosen in combinations(range(N), 6))
        ),
    }
    if counts != EXPECTED_MODEL_COUNTS:
        raise AssertionError("independent source-0 joint-model counts changed")
    return counts


def _move_mask(value: int, action: Sequence[int]) -> int:
    return sum(1 << action[v] for v in range(N) if value & (1 << v))


def _move_code(code: Sequence[int], action: Sequence[int]) -> tuple[int, ...]:
    source = dict(zip(combinations(FIXED_A4_VERTICES, 2), code, strict=True))
    target = {
        _pair(action[left], action[right]): label
        for (left, right), label in source.items()
    }
    return tuple(target[edge] for edge in combinations(FIXED_A4_VERTICES, 2))


def _swap(left: int, right: int) -> tuple[int, ...]:
    action = list(range(N))
    action[left], action[right] = action[right], action[left]
    return tuple(action)


def stabilizer_record(branch: int) -> dict[str, object]:
    orbit = placement_orbits()[branch]
    code = orbit[0]
    outside = tuple(range(6, 13))
    actions = [_swap(outside[i], outside[i + 1]) for i in range(len(outside) - 1)]
    identity = tuple(range(N))
    internal_actions = []
    for image in permutations(FIXED_A4_VERTICES):
        base = list(range(N))
        for source, target in zip(FIXED_A4_VERTICES, image, strict=True):
            base[source] = target
        for swap_active in (False, True):
            action = base.copy()
            if swap_active:
                action[0], action[1] = 1, 0
            action_tuple = tuple(action)
            if action_tuple != identity and _move_code(code, action_tuple) == code:
                internal_actions.append(action_tuple)
    actions.extend(internal_actions)
    actions = list(dict.fromkeys(actions))
    facets = set(_candidate_facets())
    for action in actions:
        if _move_mask(FIXED_A4, action) != FIXED_A4:
            raise AssertionError("stabilizer action moved the fixed A4")
        if _move_code(code, action) != code:
            raise AssertionError("stabilizer action moved the fixed code")
        if {_move_mask(facet, action) for facet in facets} != facets:
            raise AssertionError("stabilizer action moved the facet ground")
    internal_order = factorial(4) // len(orbit)
    exact_order = factorial(7) * 2 * internal_order
    if len(internal_actions) != 2 * internal_order - 1:
        raise AssertionError("internal labelled stabilizer enumeration changed")
    return {
        "source_action_structure": "S4({2,3,4,5}) x S7({6,...,12}) x C2({0,1})",
        "source_action_order": factorial(4) * factorial(7) * 2,
        "exact_internal_S4_code_stabilizer_order": internal_order,
        "exact_labelled_A4_stabilizer_order": exact_order,
        "generator_actions": [list(action) for action in actions],
        "generator_actions_digest": _digest(actions),
        "lex_generator_step_counts": [
            sum(_move_mask(facet, action) != facet for facet in facets)
            for action in actions
        ],
        "orbit_stabilizer_identity": f"{len(orbit)} * {internal_order} = 24",
        "all_actions_preserve_source_code_and_facet_ground": True,
    }


def _pattern_from_document(document: Mapping[str, object]) -> tuple[int, tuple[int, ...]]:
    order = int(document["order"])
    code = tuple(int(value) for value in document["canonical_code"])
    pattern = (order, code)
    if code != min(_code_orbit(order, code)):
        raise ValueError("nonfinite pattern is not canonical")
    _check_minimal_nonfinite(pattern)
    expected_positive = [
        [left, right, label]
        for (left, right), label in zip(combinations(range(order), 2), code, strict=True)
        if label > 2
    ]
    if document.get("positive_labels") != expected_positive:
        raise ValueError("nonfinite pattern positive-label record changed")
    if document.get("labelled_embedding_orbit_size") != len(_code_orbit(order, code)):
        raise ValueError("nonfinite pattern orbit size changed")
    return pattern


def audit_patterns(documents: Sequence[Mapping[str, object]]) -> dict[str, object]:
    patterns = tuple(sorted({_pattern_from_document(row) for row in documents}))
    canonical_documents = [
        {
            "order": order,
            "canonical_code": list(code),
            "positive_labels": [
                [left, right, label]
                for (left, right), label in zip(
                    combinations(range(order), 2), code, strict=True
                )
                if label > 2
            ],
            "labelled_embedding_orbit_size": len(_code_orbit(order, code)),
        }
        for order, code in patterns
    ]
    if canonical_documents != list(documents) or _digest(canonical_documents) != EXPECTED_PATTERN_DIGEST:
        raise ValueError("abstract nonfinite pattern library changed")
    stream = sha256()
    count = 0
    for order, code in patterns:
        for target in combinations(range(N), order):
            restriction = _mask(target)
            if _contains_base(restriction):
                continue
            target_pairs = tuple(combinations(target, 2))
            for placement in _code_orbit(order, code):
                positive = [
                    [left, right, label]
                    for (left, right), label in zip(
                        target_pairs, placement, strict=True
                    )
                    if label > 2
                ]
                stream.update(
                    json.dumps([restriction, positive], separators=(",", ":")).encode()
                    + b"\n"
                )
                count += 1
    if count != EXPECTED_CUT_COUNT or stream.hexdigest() != EXPECTED_CUT_DIGEST:
        raise ValueError("independent monotone-cut expansion changed")
    return {
        "abstract_nonfinite_pattern_count": len(patterns),
        "order_histogram": {
            str(order): number
            for order, number in sorted(Counter(order for order, _ in patterns).items())
        },
        "all_patterns_deletion_minimal_nonfinite": True,
        "abstract_nonfinite_patterns": canonical_documents,
        "abstract_nonfinite_pattern_digest": EXPECTED_PATTERN_DIGEST,
        "expanded_global_monotone_face_cut_count": count,
        "expanded_global_monotone_face_cut_digest": stream.hexdigest(),
        "independent_streaming_replay": True,
    }


def _valid_statistics(document: Mapping[str, object], branch: int) -> dict[str, object]:
    wall = document.get("terminal_solver_wall_time_seconds")
    branches = document.get("terminal_solver_branch_count")
    conflicts = document.get("terminal_solver_conflict_count")
    if (
        isinstance(wall, bool)
        or not isinstance(wall, (int, float))
        or not isfinite(wall)
        or wall < 0
        or isinstance(branches, bool)
        or not isinstance(branches, int)
        or branches < 0
        or isinstance(conflicts, bool)
        or not isinstance(conflicts, int)
        or conflicts < 0
    ):
        raise ValueError(f"branch {branch}: invalid terminal statistics")
    return {
        "solver_wall_time_seconds": wall,
        "solver_branch_count": branches,
        "solver_conflict_count": conflicts,
    }


def audit_branch_document(
    document: Mapping[str, object],
    branch: int,
    pattern_audit: Mapping[str, object],
) -> dict[str, object]:
    common = {
        "kind": "h9_14_record0_a1_joint_incidence_label_cegar_exploratory",
        "theorem_certificate": False,
        "record_index": 0,
        "graph6": GRAPH6,
        "schedule_a": 1,
        "source_orbit_index": 0,
        "representative_A4s": [list(FIXED_A4_VERTICES)],
        "a4_placement_partition": partition_record(),
        "dual_model_counts": independent_model_counts(),
        "terminal_solver_status": "INFEASIBLE",
        "stopped_reason": "terminal_joint_incidence_label_infeasible",
        "compact_label_survivor": None,
        "a4_branch_compact_coxeter_no_go": True,
        "fixed_source_compact_coxeter_no_go": False,
        "finite_volume_noncompact_checked": False,
        "scope": BRANCH_SCOPE,
    }
    for key, expected in common.items():
        if document.get(key) != expected:
            raise ValueError(f"branch {branch}: {key} changed")
    stabilizer = stabilizer_record(branch)
    expected_branch = {
        "orbit_index": branch,
        "representative_code": list(EXPECTED_REPRESENTATIVES[branch]),
        "placement_orbit_size": EXPECTED_ORBIT_SIZES[branch],
        "incidence_symmetry_slice": "exact labelled-A4 stabilizer",
        "lex_generator_step_counts": stabilizer["lex_generator_step_counts"],
        "lossless_reason": LOSSLESS_REASON,
    }
    if document.get("a4_branch") != expected_branch:
        raise ValueError(f"branch {branch}: branch metadata changed")
    patterns = pattern_audit["abstract_nonfinite_patterns"]
    if branch == 0:
        rounds = document.get("rounds")
        if not isinstance(rounds, list) or len(rounds) != 1:
            raise ValueError("branch 0 must contain the generating cut round")
        row = rounds[0]
        if (
            row.get("round") != 0
            or row.get("new_abstract_nonfinite_patterns") != patterns
            or row.get("new_abstract_nonfinite_pattern_count") != 10
            or row.get("abstract_nonfinite_pattern_digest") != EXPECTED_PATTERN_DIGEST
            or row.get("expanded_global_face_cut_count") != EXPECTED_CUT_COUNT
            or row.get("expanded_global_face_cut_digest") != EXPECTED_CUT_DIGEST
            or row.get("new_global_monotone_face_cut_count") != EXPECTED_CUT_COUNT
        ):
            raise ValueError("branch 0 cut-generating round changed")
        if (
            document.get("completed_cut_round_count") != 1
            or document.get("seeded_abstract_nonfinite_pattern_count") != 0
            or document.get("seeded_abstract_nonfinite_patterns") != []
            or document.get("seeded_global_monotone_face_cut_count") != 0
            or document.get("seed_expansion_transcript") != []
        ):
            raise ValueError("branch 0 seed/history semantics changed")
    else:
        if (
            document.get("completed_cut_round_count") != 0
            or document.get("rounds") != []
            or document.get("seeded_abstract_nonfinite_pattern_count") != 10
            or document.get("seeded_abstract_nonfinite_patterns") != patterns
            or document.get("seeded_abstract_nonfinite_pattern_digest") != EXPECTED_PATTERN_DIGEST
            or document.get("seeded_global_monotone_face_cut_count") != EXPECTED_CUT_COUNT
            or document.get("seed_expansion_transcript") != EXPECTED_SEED_TRANSCRIPT
            or document.get("seed_expansion_transcript_digest") != EXPECTED_SEED_TRANSCRIPT_DIGEST
        ):
            raise ValueError(f"branch {branch}: seeded cut library changed")
    if (
        document.get("installed_abstract_nonfinite_pattern_count") != 10
        or document.get("installed_global_monotone_face_cut_count") != EXPECTED_CUT_COUNT
    ):
        raise ValueError(f"branch {branch}: installed cut totals changed")
    return {
        "branch_index": branch,
        "representative_code": list(EXPECTED_REPRESENTATIVES[branch]),
        "placement_orbit_size": EXPECTED_ORBIT_SIZES[branch],
        "stabilizer": stabilizer,
        "terminal_solver_status": "INFEASIBLE",
        "compact_label_survivor": None,
        "branch_compact_coxeter_no_go": True,
        "terminal_statistics": _valid_statistics(document, branch),
    }


def audit_prior_sources() -> list[dict[str, object]]:
    answer = []
    for row in SOURCE_CERTIFICATES:
        path = ROOT / str(row["path"])
        if _file_sha256(path) != row["file_sha256"]:
            raise ValueError(f"prior source certificate hash changed: {path}")
        document = json.loads(path.read_text(encoding="utf-8"))
        value: object = document
        for key in row["boolean_path"]:
            if not isinstance(value, Mapping):
                raise ValueError("prior source certificate path is malformed")
            value = value[key]
        if value is not True:
            raise ValueError("prior fixed-source certificate is not terminal")
        answer.append(
            {
                "source_orbit_index": row["source_orbit_index"],
                "path": row["path"],
                "file_sha256": row["file_sha256"],
                "fixed_source_compact_no_go": True,
            }
        )
    return answer


def build_record(input_directory: Path) -> dict[str, object]:
    raw = []
    for branch in EXPECTED_BRANCH_INDICES:
        path = input_directory / f"h9_record0_source0_branch{branch:02d}.json"
        if not path.exists():
            raise ValueError(f"missing terminal branch checkpoint: {path}")
        raw.append(
            {
                "branch_index": branch,
                "input_file_sha256": _file_sha256(path),
                "document": json.loads(path.read_text(encoding="utf-8")),
            }
        )
    round_zero = raw[0]["document"]["rounds"][0]
    pattern_audit = audit_patterns(round_zero["new_abstract_nonfinite_patterns"])
    evidence = [
        {
            **audit_branch_document(row["document"], row["branch_index"], pattern_audit),
            "input_file_sha256": row["input_file_sha256"],
            "terminal_document": row["document"],
        }
        for row in raw
    ]
    prior = audit_prior_sources()
    document = {
        "schema_version": 1,
        "kind": "h9_14_record0_a1_source0_joint_branch_aggregate_certificate",
        "claim": (
            "For Record 0, compact schedule a=1, fixed source orbit 0, no joint "
            "incidence/Coxeter labelling survives any of the nine exact labelled-A4 branches."
        ),
        "fixed_source": {
            "record_index": 0,
            "graph6": GRAPH6,
            "schedule_a": 1,
            "source_orbit_index": 0,
            "representative_A4": list(FIXED_A4_VERTICES),
            "source_action_structure": "S4({2,3,4,5}) x S7({6,...,12}) x C2({0,1})",
            "source_action_order": factorial(4) * factorial(7) * 2,
        },
        "placement_partition": partition_record(),
        "independent_model_counts": independent_model_counts(),
        "cut_library_audit": pattern_audit,
        "branch_evidence": evidence,
        "coverage": {
            "expected_branch_indices": list(EXPECTED_BRANCH_INDICES),
            "present_branch_indices": list(EXPECTED_BRANCH_INDICES),
            "terminally_closed_branch_indices": list(EXPECTED_BRANCH_INDICES),
            "placement_orbit_size_sum": sum(EXPECTED_ORBIT_SIZES),
            "all_nine_terminal_branches_present_and_closed": True,
            "fixed_source_compact_coxeter_no_go": True,
        },
        "record0_a1_corollary": {
            "source0_certificate": "this document",
            "prior_source_certificates": prior,
            "closed_source_orbit_indices": [0, 1, 2, 3],
            "all_four_source_orbits_closed": True,
            "whole_record0_a1_compact_schedule_no_go": True,
        },
        "broader_frontier": {
            "whole_record0_compact_shadow_no_go": False,
            "record0_compact_schedules_still_open_after_this_result": [0, *range(2, 17)],
            "record0_compact_schedules_already_closed_elsewhere": [17, 18],
            "global_wholly_excluded_shadow_count": 134,
            "global_open_shadow_count": 8,
            "global_open_shadow_records": [0, 1, 3, 4, 5, 7, 9, 31],
            "global_partition_changed_by_this_schedule_result": False,
        },
        "theorem": {
            "fixed_source_compact_coxeter_no_go": True,
            "whole_record0_a1_compact_schedule_no_go": True,
            "whole_record0_compact_shadow_no_go": False,
            "finite_volume_noncompact_checked": False,
            "finite_volume_noncompact_chamber_excluded": False,
            "compact_label_or_gram_survivor": False,
            "exact_gram_gate_reached": False,
            "reason_gram_gate_not_reached": "Every joint incidence/finite-label branch is infeasible.",
        },
        "nonclaims": [
            "This does not exclude Record 0 in compact schedules a=0 or a=2,...,16.",
            "This does not close another whole H9/14 graph shadow or change the 134+8 compact partition.",
            "This does not test finite-volume noncompact chambers, ideal vertices, or Vinberg cusp endpoints.",
            "This constructs no hyperbolic Coxeter polytope.",
        ],
        "verifier": {
            "script": "scripts/certify_h9_14_record0_a1_source0_joint_branches.py",
            "independent_replay": (
                "Reconstructs the classical Lanner placement partition, exact stabilizers, "
                "model counts, and every deletion-minimal nonfinite cut expansion; authenticates "
                "the archived no-timeout CP-SAT terminal results without rerunning them."
            ),
        },
    }
    document["certificate_sha256"] = _digest(document)
    return document


def verify_record(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record0_a1_source0_joint_branch_aggregate_certificate":
        raise ValueError("wrong source-0 aggregate certificate kind")
    observed = _digest({key: value for key, value in document.items() if key != "certificate_sha256"})
    if document.get("certificate_sha256") != observed:
        raise ValueError("source-0 aggregate certificate digest changed")
    if EXPECTED_CERTIFICATE_SHA256 and observed != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("source-0 aggregate certificate is not frozen")
    partition = partition_record()
    if document.get("placement_partition") != partition:
        raise ValueError("archived placement partition changed")
    if document.get("independent_model_counts") != independent_model_counts():
        raise ValueError("archived model counts changed")
    cut_library = document.get("cut_library_audit")
    if not isinstance(cut_library, Mapping):
        raise ValueError("missing cut library audit")
    pattern_audit = audit_patterns(cut_library["abstract_nonfinite_patterns"])
    if cut_library != pattern_audit:
        raise ValueError("archived cut library audit changed")
    evidence = document.get("branch_evidence")
    if not isinstance(evidence, list) or len(evidence) != 9:
        raise ValueError("aggregate needs exactly nine branch documents")
    indices = []
    for row in evidence:
        branch = int(row["branch_index"])
        indices.append(branch)
        replay = audit_branch_document(row["terminal_document"], branch, pattern_audit)
        for key, expected in replay.items():
            if row.get(key) != expected:
                raise ValueError(f"branch {branch}: aggregate evidence changed")
    if indices != list(EXPECTED_BRANCH_INDICES):
        raise ValueError("aggregate branch order/coverage changed")
    prior = audit_prior_sources()
    corollary = document.get("record0_a1_corollary")
    if not isinstance(corollary, Mapping) or corollary.get("prior_source_certificates") != prior:
        raise ValueError("Record-0/a=1 prior-source ledger changed")
    coverage = document.get("coverage")
    theorem = document.get("theorem")
    if (
        not isinstance(coverage, Mapping)
        or coverage.get("all_nine_terminal_branches_present_and_closed") is not True
        or coverage.get("fixed_source_compact_coxeter_no_go") is not True
        or not isinstance(theorem, Mapping)
        or theorem.get("whole_record0_a1_compact_schedule_no_go") is not True
        or theorem.get("whole_record0_compact_shadow_no_go") is not False
        or theorem.get("finite_volume_noncompact_checked") is not False
    ):
        raise ValueError("theorem scope changed")
    return {
        "certificate_sha256": observed,
        "closed_branch_count": 9,
        "covered_classical_A4_placement_count": 108,
        "fixed_source_compact_coxeter_no_go": True,
        "whole_record0_a1_compact_schedule_no_go": True,
        "whole_record0_compact_shadow_no_go": False,
        "finite_volume_noncompact_checked": False,
        "compact_label_or_gram_survivor": False,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-from", type=Path)
    parser.add_argument("--write", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args(argv)
    if args.build_from is not None:
        document = build_record(args.build_from)
        args.write.parent.mkdir(parents=True, exist_ok=True)
        args.write.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    else:
        document = json.loads(args.write.read_text(encoding="utf-8"))
    summary = verify_record(document)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
