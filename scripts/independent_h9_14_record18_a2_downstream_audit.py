#!/usr/bin/env python3
"""Independent topology and compact-Coxeter replay for record-18 ``a=2``."""

from __future__ import annotations

import argparse
from collections import Counter
import importlib.util
import json
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parents[1]
A1_SOURCE = ROOT / "scripts/independent_h9_14_record18_a1_downstream_audit.py"
A1_SPEC = importlib.util.spec_from_file_location(
    "record18_a1_independent_primitives", A1_SOURCE
)
if A1_SPEC is None or A1_SPEC.loader is None:
    raise RuntimeError("could not load independent downstream primitives")
a1 = importlib.util.module_from_spec(A1_SPEC)
A1_SPEC.loader.exec_module(a1)
base = a1.base

RAW_CERTIFIER_SOURCE = (
    ROOT / "scripts/certify_h9_14_record18_a2_add4_ridge_raw_families.py"
)
RAW_SPEC = importlib.util.spec_from_file_location(
    "record18_a2_add4_raw_certificate", RAW_CERTIFIER_SOURCE
)
if RAW_SPEC is None or RAW_SPEC.loader is None:
    raise RuntimeError("could not load the a2/add4 literal-family verifier")
raw_certificate = importlib.util.module_from_spec(RAW_SPEC)
RAW_SPEC.loader.exec_module(raw_certificate)

DEFAULT_INPUT = (
    ROOT
    / "data/enumerations/h9_14_record18_high_addition_independent_audit.json"
)
DEFAULT_OUTPUT = (
    ROOT / "data/enumerations/h9_14_record18_a2_downstream_independent_audit.json"
)
DEFAULT_RAW_INPUT = (
    ROOT / "data/enumerations/h9_14_record18_a2_add4_ridge_raw_families.json"
)
EXPECTED_A2_REPRESENTATIVES = {
    (53, 197),
    (53, 54),
    (53, 86),
    (53, 198),
    (23, 101),
    (53, 58),
    (53, 90),
    (53, 202),
    (23, 57),
    (23, 105),
    (15, 53),
    (23, 27),
    (23, 43),
    (23, 45),
}
EXPECTED_F_VECTOR = (14, 86, 310, 733, 1192, 1344, 1020, 477, 106)


def addition_mask(row: Sequence[int]) -> int:
    return base.mask(map(int, row))


def topology_audit(
    fours: tuple[int, int], additions: tuple[int, ...]
) -> tuple[dict[str, object], tuple[int, ...]]:
    # The forced A5 shadow consists of supersets of the A4 rows, so it is
    # redundant in the minimal-nonface clutter used for direct dualization.
    clutter = tuple(sorted((*base.EDGE_MASKS, *fours, *additions)))
    complements = base.minimal_transversals(clutter)
    blocker_sizes = Counter(value.bit_count() for value in complements)
    table = base.face_table(complements)
    vector = base.face_vector(table)
    recovered = base.minimal_nonfaces(table)
    ridges = base.ridge_histogram(complements)
    components = base.dual_component_count(complements)
    triples = all(table[value] for value in base.independent_masks(3))
    strict_failure = base.first_strict_local_g_failure(table)
    macaulay_failure = base.first_local_macaulay_failure(table)
    euler_failure = base.first_link_euler_failure(table)
    homology_failure = base.first_mod2_homology_failure(table)
    valid = (
        blocker_sizes == {5: 106}
        and vector == EXPECTED_F_VECTOR
        and set(recovered) == set(clutter)
        and ridges == {2: 477}
        and components == 1
        and triples
        and strict_failure is None
        and macaulay_failure is None
        and euler_failure is None
        and homology_failure is None
    )
    if not valid:
        raise AssertionError(
            "independent a2 topology failed: "
            f"blocker={blocker_sizes}, f={vector}, ridges={ridges}, "
            f"components={components}, strict={strict_failure}, "
            f"Macaulay={macaulay_failure}, Euler={euler_failure}, "
            f"homology={homology_failure}"
        )
    return (
        {
            "facet_complement_count": len(complements),
            "blocker_size_histogram": dict(sorted(blocker_sizes.items())),
            "face_vector": list(vector),
            "minimal_nonface_size_histogram": dict(
                sorted(Counter(value.bit_count() for value in recovered).items())
            ),
            "minimal_nonfaces_match_input_clutter": True,
            "ridge_degree_histogram": ridges,
            "dual_component_count": components,
            "covers_all_independent_triples": triples,
            "strict_local_g_constraints_hold": True,
            "local_macaulay_constraints_hold": True,
            "all_links_have_spherical_euler_characteristic": True,
            "all_links_are_mod2_homology_spheres": True,
            "passes_connected_homology_audit": True,
        },
        complements,
    )


def build_record(
    source: Path, raw_source: Path = DEFAULT_RAW_INPUT
) -> dict[str, object]:
    document = json.loads(source.read_text(encoding="utf-8"))
    a2_cases = [case for case in document["cases"] if case["a"] == 2]
    observed_representatives = {
        tuple(case["representative_A4_masks"]) for case in a2_cases
    }
    if observed_representatives != EXPECTED_A2_REPRESENTATIVES:
        raise AssertionError("the fourteen a2 representatives are not exact")

    raw_document = json.loads(raw_source.read_text(encoding="utf-8"))
    raw_certificate.verify_record(raw_document)
    raw_sets = {
        tuple(entry["representative_A4_masks"]): tuple(
            tuple(map(int, family)) for family in entry["target_families"]
        )
        for entry in raw_document["normalized_raw_target_sets"]
    }
    add4_representatives = {
        tuple(case["representative_A4_masks"])
        for case in a2_cases
        if case["additional_A5_count"] == 4
    }
    if set(raw_sets) != add4_representatives:
        raise AssertionError("raw-family archive does not cover every add4 case")
    add5_cases = [case for case in a2_cases if case["additional_A5_count"] == 5]
    if len(add5_cases) != 2 or any(
        case["canonical_target_count"] != 0 for case in add5_cases
    ):
        raise AssertionError("the two add5 cases are not exact incidence no-gos")

    audits = []
    for fours in raw_certificate.REPRESENTATIVES:
        for literal_index, additions in enumerate(raw_sets[fours]):
            topology, complements = topology_audit(fours, additions)
            supports = tuple(base.vertices(row) for row in (*fours, *additions))
            labels = a1.label_gate(supports, complements)
            if not labels["unrestricted_compact_coxeter_no_go"]:
                raise AssertionError("an a2 target class survived compact labels")
            audits.append(
                {
                    "A4_masks": list(fours),
                    "literal_family_index": literal_index,
                    "representative_addition_masks": list(additions),
                    "representative_additions": [
                        list(base.vertices(row)) for row in additions
                    ],
                    "topology": topology,
                    "coxeter_labels": labels,
                }
            )

    observed = sorted(
        (
            tuple(audit["A4_masks"]),
            tuple(
                stage["passing_state_count"]
                for stage in audit["coxeter_labels"]["stages"]
            ),
        )
        for audit in audits
    )
    expected = sorted(
        [((53, 198), (108, 84, 4, 28, 196, 0))] * 6
        + [((53, 202), (108, 84, 8, 16, 0))] * 2
    )
    if observed != expected:
        raise AssertionError(f"a2 compact label transcripts changed: {observed}")

    return {
        "kind": "independent_h9_14_record18_a2_downstream_audit",
        "source": str(source),
        "raw_family_source": str(raw_source),
        "a2_representative_count": len(a2_cases),
        "a2_add4_target_raw_weight": len(audits),
        "a2_add4_target_subgroup_class_count": sum(
            case["canonical_target_count"]
            for case in a2_cases
            if case["additional_A5_count"] == 4
        ),
        "a2_add4_literal_family_count_audited": len(audits),
        "all_target_classes_pass_connected_mod2_homology_audit": True,
        "all_target_classes_have_unrestricted_compact_coxeter_no_go": True,
        "classes": audits,
        "a2_compact_branch_excluded": True,
        "compact_only": True,
        "finite_volume_noncompact_checked": False,
        "whole_shadow_compact_coxeter_theorem": False,
        "scope": (
            "Independent downstream audit of all eight literal a=2 incidence "
            "targets after exact production/independent set equality. "
            "No ideal-vertex or finite-volume noncompact conclusion."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, nargs="?", default=DEFAULT_INPUT)
    parser.add_argument("--raw-input", type=Path, default=DEFAULT_RAW_INPUT)
    parser.add_argument("--write", type=Path)
    args = parser.parse_args()
    record = build_record(args.input, args.raw_input)
    rendered = json.dumps(record, indent=2, sort_keys=True) + "\n"
    if args.write is not None:
        args.write.parent.mkdir(parents=True, exist_ok=True)
        args.write.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
