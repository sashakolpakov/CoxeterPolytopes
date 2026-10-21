"""Exact wall and Coxeter screen for Jacobian--Kummer class_0011.

The input chamber is the published 316-wall chamber ``D16`` of a general
Jacobian--Kummer K3 surface.  ``class_0011`` is one of the two conjugacy
classes of subgroups of the full lattice/chamber stabilizer ``O(S16,D16)``
whose common fixed lattice has rank 11.  Thus its fixed slice is hyperbolic
10-space.  No surface-automorphism condition is imposed.

Every orbit-sum wall is treated exactly.  A removed inequality comes with a
positive rational conic dependence, while every retained inequality comes
with a strict timelike relative-interior witness.  Only this certified facet
set is sent through the normalized Coxeter-product gate.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
from typing import Any

import sympy as sp

from scripts.search_apery_fermi_chamber_folds import (
    _exact_nonnegative_dependence,
    _orbits,
)
from scripts.search_jacobian_kummer_chamber_folds import (
    ROOT,
    SOURCE,
    _matrix_payload,
    _parse_source,
    _rational_text,
    _wall_permutation,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia
from scripts.verify_k3_quotient_nef_slices import _coxeter_label


CENSUS = (
    ROOT
    / "data"
    / "enumerations"
    / "jacobian_kummer_d16_subgroup_folds.json"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "jacobian_kummer_h10_class0011_wall_certificate.json"
)


def _matrix(payload: list[list[str]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(entry) for entry in row] for row in payload])


def _representative(path: Path = CENSUS) -> tuple[dict[str, Any], sp.Matrix]:
    census = json.loads(path.read_text(encoding="utf-8"))
    records = [
        record
        for record in census["high_dimension_candidates"]
        if record["id"] == "class_0011"
    ]
    if len(records) != 1:
        raise ArithmeticError("the frozen census does not contain unique class_0011")
    record = records[0]
    matrices = record["representative_generator_matrices"]
    if not (
        record["gap_class_index"] == 11
        and record["subgroup_order"] == 2
        and record["fixed_rank"] == 11
        and record["hyperbolic_dimension"] == 10
        and len(matrices) == 1
    ):
        raise ArithmeticError("the frozen class_0011 metadata changed")
    return record, _matrix(matrices[0])


def build_certificate(
    source_path: Path = SOURCE,
    census_path: Path = CENSUS,
) -> dict[str, Any]:
    source = _parse_source(source_path)
    record, involution = _representative(census_path)
    gram: sp.Matrix = source["gram"]
    ample: sp.Matrix = source["ample"]
    walls = tuple(
        sp.Matrix(1, gram.rows, vector) for vector in source["wall_vectors"]
    )
    identity = sp.eye(gram.rows)
    if not (
        involution**2 == identity
        and involution * gram * involution.T == gram
        and ample * involution == ample
        and gram.rows - (involution.T - identity).rank() == 11
    ):
        raise ArithmeticError("class_0011 is not the advertised rank-11 involution")

    permutation = _wall_permutation(
        involution, gram, ample, source["wall_vectors"]
    )
    subgroup = frozenset((tuple(range(len(permutation))), permutation))
    wall_orbits = _orbits(subgroup, len(walls))
    orbit_normals = tuple(
        sum((walls[index] for index in orbit), sp.zeros(1, gram.rows))
        for orbit in wall_orbits
    )
    stacked = sp.Matrix.vstack(*orbit_normals)
    if stacked.rank() != 11:
        raise ArithmeticError("class_0011 orbit normals do not span its fixed space")
    _, pivot_columns = stacked.rref()
    coordinates = stacked[:, list(pivot_columns)]
    orbit_gram = stacked * gram * stacked.T
    inertia = exact_inertia(orbit_gram)
    if inertia != (1, 10, len(wall_orbits) - 11):
        raise ArithmeticError("class_0011 fixed wall span is not Lorentzian H10")

    facets: list[int] = []
    facet_witnesses: list[dict[str, Any]] = []
    redundancy_certificates: list[dict[str, Any]] = []
    for target, normal in enumerate(orbit_normals):
        dependence = _exact_nonnegative_dependence(coordinates, target)
        if dependence is not None:
            redundancy_certificates.append(
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
            raise ArithmeticError("a nonspacelike extreme inequality remains")
        ample_pairing = sp.Rational((ample * gram * normal.T)[0])
        witness = ample - ample_pairing / norm * normal
        witness_square = sp.Rational((witness * gram * witness.T)[0])
        pairings = tuple(
            sp.Rational((witness * gram * other.T)[0])
            for other in orbit_normals
        )
        if not (
            witness_square > 0
            and pairings[target] == 0
            and all(
                pairing > 0
                for index, pairing in enumerate(pairings)
                if index != target
            )
        ):
            raise ArithmeticError(
                "an extreme class_0011 inequality lacks a strict timelike witness"
            )
        facets.append(target)
        facet_witnesses.append(
            {
                "orbit_index": target + 1,
                "witness_square": _rational_text(witness_square),
                "minimum_strict_pairing": _rational_text(
                    min(
                        pairing
                        for index, pairing in enumerate(pairings)
                        if index != target
                    )
                ),
            }
        )

    if len(facets) + len(redundancy_certificates) != len(wall_orbits):
        raise ArithmeticError("class_0011 wall reduction is not exhaustive")

    labels: Counter[str] = Counter()
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
                bad_pairs.append(
                    {
                        "facet_orbit_indices": [left + 1, right + 1],
                        "source_wall_orbits": [
                            [index + 1 for index in wall_orbits[left]],
                            [index + 1 for index in wall_orbits[right]],
                        ],
                        "left_square": _rational_text(left_norm),
                        "right_square": _rational_text(right_norm),
                        "inner_product": _rational_text(product),
                        "normalized_squared_product": _rational_text(
                            sp.Rational(ratio)
                        ),
                        "reason": str(error),
                    }
                )
            else:
                labels[label] += 1

    facet_gram = orbit_gram.extract(facets, facets)
    all_coxeter = not bad_pairs
    return {
        "schema": "jacobian_kummer_h10_class0011_wall_certificate_v1",
        "status": {
            "exact_fixed_slice": True,
            "exact_wall_reduction": True,
            "strict_facet_witnesses": True,
            "exact_coxeter_pair_gate": True,
            "all_facet_pairs_are_coxeter": all_coxeter,
            "vinberg_volume_test_applicable": all_coxeter,
            "decision": "coxeter_candidate" if all_coxeter else "not_coxeter",
        },
        "source": {
            "surface": "general Jacobian--Kummer K3 surface Y16",
            "divisor_lattice": "S16 = NS(Y16)",
            "published_chamber": "D16 with 316 walls",
            "paper_doi": "10.1007/s12215-019-00464-7",
            "paper_arxiv": "1908.05390",
            "source_file": str(source_path.relative_to(ROOT)),
            "surface_automorphism_required": False,
        },
        "group_class": {
            "ambient_group": "O(S16,D16) = (C2)^5 : S6",
            "ambient_group_order": 23040,
            "gap_subgroup_class_index": record["gap_class_index"],
            "id": record["id"],
            "subgroup_structure": "C2",
            "conjugacy_class_size": record["conjugacy_class_size"],
            "normalizer_order": record["normalizer_order"],
            "representative_matrix": _matrix_payload(involution),
            "fixed_rank": 11,
            "hyperbolic_dimension": 10,
        },
        "wall_orbits": [[index + 1 for index in orbit] for orbit in wall_orbits],
        "wall_orbit_count": len(wall_orbits),
        "wall_orbit_size_histogram": {
            str(size): count
            for size, count in sorted(Counter(map(len, wall_orbits)).items())
        },
        "orbit_normal_rank": stacked.rank(),
        "orbit_gram_inertia": list(inertia),
        "facet_count": len(facets),
        "facet_orbit_indices": [index + 1 for index in facets],
        "facet_source_wall_orbits": [
            [source_index + 1 for source_index in wall_orbits[index]]
            for index in facets
        ],
        "redundant_orbit_count": len(redundancy_certificates),
        "redundancy_certificates": redundancy_certificates,
        "facet_witnesses": facet_witnesses,
        "facet_normal_square_histogram": {
            _rational_text(square): count
            for square, count in sorted(
                Counter(
                    sp.Rational(facet_gram[index, index])
                    for index in range(facet_gram.rows)
                ).items()
            )
        },
        "coxeter_edge_counts": dict(sorted(labels.items())),
        "non_coxeter_pair_count": len(bad_pairs),
        "non_coxeter_normalized_product_histogram": {
            ratio: count
            for ratio, count in sorted(
                Counter(
                    pair["normalized_squared_product"] for pair in bad_pairs
                ).items()
            )
        },
        "non_coxeter_pairs": bad_pairs,
        "facet_gram_sign_convention": "negative-square orbit-sum wall normals",
        "facet_gram": _matrix_payload(facet_gram),
        "positive_wall_gram": _matrix_payload(-facet_gram),
        "checks": {
            "representative_is_exact_rank_11_lattice_involution": True,
            "representative_preserves_gram_and_ample_class": True,
            "representative_permutes_all_316_published_walls": True,
            "orbit_normals_span_fixed_lorentzian_space": True,
            "every_removed_wall_has_exact_nonnegative_dependence": True,
            "every_retained_wall_has_strict_timelike_relative_interior_witness": True,
            "all_retained_facet_pairs_checked_exactly": True,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--census", type=Path, default=CENSUS)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    payload = build_certificate(arguments.source, arguments.census)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file() or arguments.output.read_text() != rendered:
            raise SystemExit(f"archived artifact differs: {arguments.output}")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
    print(
        json.dumps(
            {
                "decision": payload["status"]["decision"],
                "wall_orbit_count": payload["wall_orbit_count"],
                "facet_count": payload["facet_count"],
                "redundant_orbit_count": payload["redundant_orbit_count"],
                "coxeter_edge_counts": payload["coxeter_edge_counts"],
                "non_coxeter_pair_count": payload["non_coxeter_pair_count"],
                "non_coxeter_normalized_product_histogram": payload[
                    "non_coxeter_normalized_product_histogram"
                ],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
