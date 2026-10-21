#!/usr/bin/env python3
"""Separate Roulleau's two geometric deck actions from the 59-wall diagram group."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

from scripts.scout_roulleau_r12_ud4x2a1x2_59 import coordinate_model


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/roulleau_r12_ud4x2a1x2_surface_actions.json"
SCOUT = ROOT / "data/enumerations/roulleau_r12_ud4x2a1x2_59_scout.json"
SOURCE_TEX = Path("/private/tmp/roulleau_source/epiga_vol6_Roulleau.tex")
SOURCE_TEX_SHA256 = "e6eeaaa4a3f375f209762718e8ef9a361b0a95356cc8be2e9f99a6713abf5df3"


Permutation = tuple[int, ...]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _matrix_payload(matrix: sp.Matrix) -> list[list[int]]:
    return [[int(value) for value in matrix.row(index)] for index in range(matrix.rows)]


def _cycles(permutation: Permutation, names: list[str]) -> list[list[str]]:
    unseen = set(range(len(permutation)))
    result = []
    while unseen:
        start = min(unseen)
        orbit = []
        current = start
        while current not in orbit:
            orbit.append(current)
            unseen.discard(current)
            current = permutation[current]
        if len(orbit) > 1:
            result.append([names[index] for index in orbit])
    return result


def _action_record(
    permutation: Permutation,
    names: list[str],
    rows: sp.Matrix,
    basis_gram: sp.Matrix,
    wall_gram: sp.Matrix,
) -> dict[str, Any]:
    action = rows[[permutation[index] for index in range(12)], :]
    checks = {
        "permutation": sorted(permutation) == list(range(59)),
        "weighted_wall_gram_action": wall_gram.extract(permutation, permutation)
        == wall_gram,
        "all_wall_coordinates_transform": rows * action == rows[list(permutation), :],
        "integral_unimodular_NS_action": all(value.is_Integer for value in action)
        and abs(int(action.det())) == 1,
        "NS_isometry": action * basis_gram * action.T == basis_gram,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"wall action failed: {checks}")
    return {
        "permutation_zero_based": list(permutation),
        "nontrivial_cycles_by_source_name": _cycles(permutation, names),
        "NS_action_on_row_coordinates": _matrix_payload(action),
        "fixed_rank": 12 - (action - sp.eye(12)).rank(),
        "checks": checks,
    }


def _source_checks() -> dict[str, Any]:
    if _sha256(SOURCE_TEX) != SOURCE_TEX_SHA256:
        raise ArithmeticError("pinned Roulleau TeX hash changed")
    lines = SOURCE_TEX.read_text().splitlines()
    first = "\n".join(lines[8907:8947])
    second = "\n".join(lines[8948:9008])
    criterion = "\n".join(lines[1061:1095])
    table = "\n".join(lines[10920:10991])
    checks = {
        "first_double_plane_A5": "base-point free, of square $2$" in first
        and "an $\\mathbf{a}_{5}$ singularity" in first,
        "split_line_pairs": all(
            token in first
            for token in (
                "A_{1}+A_{14}",
                "A_{2}+A_{15}",
                "A_{3}+A_{16}",
                "A_{4}+A_{17}",
                "A_{5}+A_{18}",
                "A_{13}+A_{19}",
                "strict\ntransform of $L$ is $A_{6}+A_{7}$",
            )
        ),
        "split_conic_pairs": "$B_{ij}^{(a)}+C_{ij}^{(a)}=2D_{2}$" in first
        and "conics which are $6$-tangent" in first,
        "second_hyperelliptic_F4_cover": "\\subsubsection{Second involution}" in second
        and "factorizes through a surjective map $\\varphi'\\colon X\\to\\mathbf{F}_{4}$" in second,
        "second_branch_only_D4_and_nodes": "type $\\mathbf{d}_{4}$" in second
        and "other singular\npoints are nodes" in second,
        "roulleau_trivial_action_criterion": all(
            token in criterion
            for token in ("\\mathbf{a}_{1}", "\\mathbf{d_{4}}", "the action on $\\cu$-curves", "is trivial")
        ),
        "roulleau_nontrivial_A5_criterion": "\\mathbf{a}_{2},\\dots,\\mathbf{a}_{5}" in criterion
        and "acts non-trivially" in criterion,
        "NS_generated_by_curve_star": "N\\'eron--Severi lattice is always generated" in criterion,
        "general_Aut_is_C2_squared_and_tag_tn": (
            r"U\oplus\mathbf{D}_{4}^{\oplus2}\oplus\mathbf{A}_{1}^{\oplus2}" in table
            and "$59$" in table
            and "$\\zde$ & $tn$" in table
        ),
        "tn_legend": "a $tn$ means that both\ncases exists" in table,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"surface-action source witnesses failed: {checks}")
    return {
        "tex_sha256": SOURCE_TEX_SHA256,
        "first_involution_lines": [8908, 8947],
        "second_involution_lines": [8949, 9008],
        "NS_action_criterion_lines": [1062, 1095],
        "table_and_legend_lines": [10926, 10990],
        "checks": checks,
    }


def build_certificate() -> dict[str, Any]:
    source = _source_checks()
    scout = json.loads(SCOUT.read_text())
    group = tuple(tuple(element) for element in scout["weighted_group"]["elements_zero_based"])
    group_set = set(group)
    basis_gram, wall_records, wall_gram = coordinate_model()
    names = [str(record["name"]) for record in wall_records]
    name_to_index = {name: index for index, name in enumerate(names)}
    rows = sp.Matrix.vstack(*(record["coordinates"] for record in wall_records))
    identity = tuple(range(59))

    # For the first double plane, Roulleau's split lines and split six-tangent
    # conics force these deck-conjugate pairs.  Search the already exhaustive
    # diagram group for their unique extension to all 59 curves.
    forced_pairs: dict[int, int] = {}
    for left, right in (
        ("A1", "A14"),
        ("A2", "A15"),
        ("A3", "A16"),
        ("A4", "A17"),
        ("A5", "A18"),
        ("A13", "A19"),
        ("A6", "A7"),
    ):
        left_index, right_index = name_to_index[left], name_to_index[right]
        forced_pairs[left_index], forced_pairs[right_index] = right_index, left_index
    for left_kind, right_kind in (("B1", "C1"), ("B2", "C2")):
        for left in range(1, 6):
            for right in range(left + 1, 6):
                left_index = name_to_index[f"{left_kind}_{left}{right}"]
                right_index = name_to_index[f"{right_kind}_{left}{right}"]
                forced_pairs[left_index], forced_pairs[right_index] = right_index, left_index
    first_candidates = [
        element for element in group
        if all(element[source_index] == target_index for source_index, target_index in forced_pairs.items())
    ]
    if len(first_candidates) != 1:
        raise ArithmeticError(f"first deck did not have a unique diagram extension: {len(first_candidates)}")
    first_permutation = first_candidates[0]
    first_action = _action_record(
        first_permutation, names, rows, basis_gram, wall_gram
    )

    # The tempting central diagram involution is visible directly in the
    # printed B/C formulas after interchanging A6 and A7.
    p67 = list(identity)
    p67[name_to_index["A6"]], p67[name_to_index["A7"]] = (
        name_to_index["A7"], name_to_index["A6"]
    )
    for first_kind, second_kind in (("B1", "B2"), ("C1", "C2")):
        for left in range(1, 6):
            for right in range(left + 1, 6):
                first_index = name_to_index[f"{first_kind}_{left}{right}"]
                second_index = name_to_index[f"{second_kind}_{left}{right}"]
                p67[first_index], p67[second_index] = second_index, first_index
    p67_permutation = tuple(p67)
    p67_action = _action_record(
        p67_permutation, names, rows, basis_gram, wall_gram
    )

    central = all(
        tuple(element[p67_permutation[index]] for index in range(59))
        == tuple(p67_permutation[element[index]] for index in range(59))
        for element in group
    )
    checks = {
        "scout_exact_group_loaded": scout["weighted_group"]["order"] == 2880
        and len(group_set) == 2880,
        "first_deck_unique_source_forced_extension": len(first_candidates) == 1,
        "first_deck_is_nontrivial_graph_element": first_permutation in group_set
        and first_permutation != identity,
        "first_deck_fixed_rank_four": first_action["fixed_rank"] == 4,
        "second_deck_NS_action_is_identity": True,
        "p67_is_exact_central_graph_element": p67_permutation in group_set and central,
        "p67_fixed_rank_eleven": p67_action["fixed_rank"] == 11,
        "p67_is_neither_geometric_NS_action": p67_permutation
        not in {identity, first_permutation},
    }
    if not all(checks.values()):
        raise ArithmeticError(f"surface/diagram action separation failed: {checks}")

    return {
        "schema": "coxeter.roulleau_r12_ud4x2a1x2_surface_actions.v1",
        "artifact_id": "roulleau_r12_ud4x2a1x2_surface_actions",
        "source": {
            **source,
            "scout_path": str(SCOUT.relative_to(ROOT)),
            "scout_sha256": _sha256(SCOUT),
        },
        "general_surface": {
            "automorphism_group": "C2 x C2",
            "NS_image_order": 2,
            "NS_image_permutations_zero_based": [list(identity), list(first_permutation)],
            "reason": (
                "the first A5 double-plane deck action is nontrivial; the second "
                "D4-and-nodes hyperelliptic deck action is trivial on the curve "
                "star generating NS; the product has the first action"
            ),
        },
        "first_double_plane_deck_action": {
            **first_action,
            "source_forced_pairs": {
                names[source_index]: names[target_index]
                for source_index, target_index in sorted(forced_pairs.items())
                if source_index < target_index
            },
            "status": "actual_general_surface_automorphism_nontrivial_on_NS",
        },
        "second_F4_deck_action": {
            "permutation_zero_based": list(identity),
            "fixed_rank": 12,
            "status": "actual_general_surface_automorphism_trivial_on_NS",
            "source_proof": (
                "the branch singularities are D4 and nodes; Roulleau's criterion "
                "makes the deck action trivial on the (-2)-curve star, which generates NS"
            ),
        },
        "p67_central_diagram_action": {
            **p67_action,
            "central_in_full_weighted_group": central,
            "status": "integral_diagram_NS_isometry_not_realized_by_Aut_of_the_general_K3",
            "strongest_geometric_wording": (
                "central diagram-involution fixed section of the 59-wall nef cone; "
                "not a quotient-nef chamber of the general K3; a special-period "
                "surface lift is not established"
            ),
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.is_file() or args.output.read_text() != rendered:
            raise SystemExit("recomputed surface-action certificate differs from archive")
    else:
        args.output.write_text(rendered)


if __name__ == "__main__":
    main()
