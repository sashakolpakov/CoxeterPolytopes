"""Exact wall and Coxeter screen for Jacobian--Kummer class_0007.

The class is the GAP subgroup conjugacy class 7 frozen by
``search_jacobian_kummer_chamber_folds.py``.  It is an involution in the full
lattice/chamber stabilizer O(S16,D16), not necessarily a surface
automorphism.  This script restricts all 316 published D16 wall inequalities
to its rank-11 fixed space, proves the irredundant facet set exactly, and
checks every facet pair against the crystallographic Coxeter-angle list.

Vinberg's finite-volume test is intentionally gated on the pair screen: a
single forbidden product already proves that the restricted chamber is not a
Coxeter polytope, so no volume claim is then made.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
from typing import Any

import sympy as sp

from scripts.search_apery_fermi_chamber_folds import _exact_nonnegative_dependence
from scripts.search_jacobian_kummer_chamber_folds import (
    OUTPUT as SUBGROUP_CENSUS,
    ROOT,
    SOURCE,
    _matrix_payload,
    _parse_source,
    _rational_text,
    _wall_permutation,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia
from scripts.verify_k3_quotient_nef_slices import _coxeter_label


OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "jacobian_kummer_h10_class0007_screen.json"
)


def _matrix_from_payload(payload: list[list[str]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(entry) for entry in row] for row in payload])


def _class_record(path: Path = SUBGROUP_CENSUS) -> dict[str, Any]:
    census = json.loads(path.read_text(encoding="utf-8"))
    records = [
        record
        for record in census["high_dimension_candidates"]
        if record["id"] == "class_0007"
    ]
    if len(records) != 1:
        raise ArithmeticError("the frozen census does not contain unique class_0007")
    record = records[0]
    expected = {
        "gap_class_index": 7,
        "subgroup_order": 2,
        "conjugacy_class_size": 60,
        "normalizer_order": 384,
        "fixed_rank": 11,
        "hyperbolic_dimension": 10,
    }
    if any(record[key] != value for key, value in expected.items()):
        raise ArithmeticError("class_0007 metadata differs from the frozen census")
    if len(record["representative_generator_matrices"]) != 1:
        raise ArithmeticError("class_0007 does not have one frozen generator")
    return record


def _wall_orbits(permutation: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    unseen = set(range(len(permutation)))
    orbits: list[tuple[int, ...]] = []
    while unseen:
        start = min(unseen)
        orbit = tuple(sorted({start, permutation[start]}))
        if any(permutation[index] not in orbit for index in orbit):
            raise ArithmeticError("the selected wall permutation is not an involution")
        unseen.difference_update(orbit)
        orbits.append(orbit)
    return tuple(orbits)


def build_screen(
    source_path: Path = SOURCE,
    census_path: Path = SUBGROUP_CENSUS,
) -> dict[str, Any]:
    source = _parse_source(source_path)
    gram: sp.Matrix = source["gram"]
    ample: sp.Matrix = source["ample"]
    walls = tuple(
        sp.Matrix(1, gram.rows, vector) for vector in source["wall_vectors"]
    )
    record = _class_record(census_path)
    involution = _matrix_from_payload(
        record["representative_generator_matrices"][0]
    )
    if not (
        involution.rows == involution.cols == gram.rows == 17
        and involution**2 == sp.eye(17)
        and involution != sp.eye(17)
        and involution * gram * involution.T == gram
        and ample * involution == ample
        and involution.trace() == 5
    ):
        raise ArithmeticError("the frozen class_0007 matrix failed exact replay")

    permutation = _wall_permutation(
        involution, gram, ample, source["wall_vectors"]
    )
    orbits = _wall_orbits(permutation)
    orbit_normals = tuple(
        sum((walls[index] for index in orbit), sp.zeros(1, gram.rows))
        for orbit in orbits
    )
    if any(normal * involution != normal for normal in orbit_normals):
        raise ArithmeticError("an orbit-sum normal is not fixed")
    stacked = sp.Matrix.vstack(*orbit_normals)
    fixed_rank = 17 - (involution.T - sp.eye(17)).rank()
    if fixed_rank != 11 or stacked.rank() != fixed_rank:
        raise ArithmeticError("the orbit normals do not span the rank-11 fixed space")
    _, pivot_columns = stacked.rref()
    coordinates = stacked[:, list(pivot_columns)]
    orbit_gram = stacked * gram * stacked.T
    inertia = exact_inertia(orbit_gram)
    if inertia != (1, 10, len(orbits) - 11):
        raise ArithmeticError("the class_0007 fixed wall span is not Lorentzian")

    facets: list[int] = []
    redundant: list[dict[str, Any]] = []
    facet_witnesses: list[dict[str, Any]] = []
    for target, normal in enumerate(orbit_normals):
        dependence = _exact_nonnegative_dependence(coordinates, target)
        if dependence is not None:
            redundant.append(
                {
                    "orbit_index": target + 1,
                    "terms": [
                        {
                            "orbit_index": support + 1,
                            "coefficient": _rational_text(coefficient),
                        }
                        for support, coefficient in dependence
                    ],
                }
            )
            continue

        norm = sp.Rational(orbit_gram[target, target])
        if norm >= 0:
            raise ArithmeticError("a nonspacelike extreme wall remains")
        ample_pairing = sp.Rational((ample * gram * normal.T)[0])
        witness = ample - ample_pairing / norm * normal
        witness_square = sp.Rational((witness * gram * witness.T)[0])
        pairings = tuple(
            sp.Rational(value) for value in list(witness * gram * stacked.T)
        )
        if not (
            witness_square > 0
            and pairings[target] == 0
            and all(
                value > 0
                for index, value in enumerate(pairings)
                if index != target
            )
        ):
            raise ArithmeticError(
                f"extreme orbit wall {target + 1} lacks the strict projection witness"
            )
        facets.append(target)
        facet_witnesses.append(
            {
                "orbit_index": target + 1,
                "witness_square": _rational_text(witness_square),
                "minimum_strict_pairing": _rational_text(
                    min(
                        value
                        for index, value in enumerate(pairings)
                        if index != target
                    )
                ),
            }
        )

    if len(facets) + len(redundant) != len(orbits):
        raise ArithmeticError("the exact wall reduction is not exhaustive")

    good_labels: Counter[str] = Counter()
    bad_ratios: Counter[str] = Counter()
    bad_reasons: Counter[str] = Counter()
    bad_pairs: list[dict[str, Any]] = []
    for left_position, left in enumerate(facets):
        for right in facets[left_position + 1 :]:
            left_norm = sp.Rational(orbit_gram[left, left])
            right_norm = sp.Rational(orbit_gram[right, right])
            product = sp.Rational(orbit_gram[left, right])
            ratio = sp.factor(product**2 / (left_norm * right_norm))
            try:
                label = _coxeter_label(left_norm, right_norm, product)
            except ValueError as error:
                ratio_text = _rational_text(sp.Rational(ratio))
                reason = str(error)
                bad_ratios[ratio_text] += 1
                bad_reasons[reason] += 1
                bad_pairs.append(
                    {
                        "facet_positions": [left_position + 1, facets.index(right) + 1],
                        "facet_orbit_indices": [left + 1, right + 1],
                        "source_wall_orbits": [
                            [index + 1 for index in orbits[left]],
                            [index + 1 for index in orbits[right]],
                        ],
                        "left_square": _rational_text(left_norm),
                        "right_square": _rational_text(right_norm),
                        "inner_product": _rational_text(product),
                        "normalized_squared_product": ratio_text,
                        "reason": reason,
                    }
                )
            else:
                good_labels[label] += 1

    facet_gram = orbit_gram.extract(facets, facets)
    all_coxeter = not bad_pairs
    return {
        "schema": "jacobian_kummer_h10_class0007_exact_screen_v1",
        "status": {
            "exact_fixed_slice": True,
            "exact_irredundant_facets": True,
            "exact_full_coxeter_pair_gate": True,
            "all_facet_pairs_are_coxeter": all_coxeter,
            "vinberg_finite_volume_audit_applicable": all_coxeter,
            "finite_volume_claimed": False,
            "decision": (
                "coxeter_candidate_requires_vinberg_audit"
                if all_coxeter
                else "not_a_coxeter_polytope"
            ),
        },
        "source": {
            "surface": "general Jacobian Kummer K3 surface Y16",
            "paper_doi": "10.1007/s12215-019-00464-7",
            "paper_arxiv": "1908.05390",
            "published_chamber": "D16",
            "published_wall_assignment": "D16WallRecs",
            "source_file": str(source_path.relative_to(ROOT)),
            "subgroup_census": str(census_path.relative_to(ROOT)),
        },
        "group_class": {
            "ambient_group": "O(S16,D16) = (C2)^5 : S6",
            "ambient_group_order": 23040,
            "surface_automorphism_required": False,
            "gap_subgroup_class_index": 7,
            "subgroup_structure": "C2",
            "conjugacy_class_size": 60,
            "normalizer_order": 384,
            "representative_matrix": _matrix_payload(involution),
            "representative_trace": int(involution.trace()),
            "fixed_rank": fixed_rank,
            "hyperbolic_dimension": fixed_rank - 1,
        },
        "wall_orbits": [[index + 1 for index in orbit] for orbit in orbits],
        "wall_orbit_count": len(orbits),
        "wall_orbit_size_histogram": {
            str(size): count
            for size, count in sorted(Counter(map(len, orbits)).items())
        },
        "orbit_normal_rank": stacked.rank(),
        "orbit_gram_inertia": list(inertia),
        "facet_count": len(facets),
        "facet_orbit_indices": [index + 1 for index in facets],
        "facet_source_wall_orbits": [
            [source_index + 1 for source_index in orbits[index]] for index in facets
        ],
        "redundant_orbit_count": len(redundant),
        "redundancy_certificates": redundant,
        "facet_witnesses": facet_witnesses,
        "facet_normal_square_histogram": {
            _rational_text(square): count
            for square, count in sorted(
                Counter(
                    sp.Rational(facet_gram[index, index])
                    for index in range(len(facets))
                ).items()
            )
        },
        "facet_pair_count": len(facets) * (len(facets) - 1) // 2,
        "coxeter_edge_counts": dict(sorted(good_labels.items())),
        "non_coxeter_pair_count": len(bad_pairs),
        "non_coxeter_ratio_histogram": dict(sorted(bad_ratios.items())),
        "non_coxeter_reason_histogram": dict(sorted(bad_reasons.items())),
        "non_coxeter_pairs": bad_pairs,
        "facet_gram_sign_convention": "negative-square orbit-sum wall normals",
        "facet_gram": _matrix_payload(facet_gram),
        "positive_wall_gram": _matrix_payload(-facet_gram),
        "checks": {
            "representative_is_an_exact_nontrivial_involution": True,
            "representative_preserves_gram_and_ample_class": True,
            "representative_permutes_all_316_published_walls": True,
            "orbit_normals_span_fixed_lorentzian_space": True,
            "every_removed_wall_has_an_exact_nonnegative_dependence": True,
            "every_retained_wall_has_a_strict_timelike_relative_interior_witness": True,
            "every_irredundant_facet_pair_was_checked_exactly": True,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--census", type=Path, default=SUBGROUP_CENSUS)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    payload = build_screen(arguments.source, arguments.census)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file() or arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archived artifact differs: {arguments.output}")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
    print(
        json.dumps(
            {
                "class": "class_0007",
                "dimension": payload["group_class"]["hyperbolic_dimension"],
                "wall_orbit_count": payload["wall_orbit_count"],
                "facet_count": payload["facet_count"],
                "redundant_orbit_count": payload["redundant_orbit_count"],
                "non_coxeter_pair_count": payload["non_coxeter_pair_count"],
                "non_coxeter_ratio_histogram": payload["non_coxeter_ratio_histogram"],
                "decision": payload["status"]["decision"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
