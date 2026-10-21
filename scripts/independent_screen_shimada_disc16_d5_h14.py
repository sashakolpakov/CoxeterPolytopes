"""Independent exact screen of Shimada disc16 chamber D5 in H14.

The external input is Shimada's official ``threediscs/disc16.txt`` atlas.
Within the archived order-20 stabilizer of chamber D5, the representative is
selected intrinsically as the unique element whose fixed space on NS has
rank 15.  No proposed Gram matrix or wall list from another search is read.

The script reconstructs the action on all 165 primitive dual-basis walls,
forms the 100 involution-orbit sums, proves the exact irredundant inequality
set, and checks every retained pair against the crystallographic Coxeter
angle list.  Finite volume is deliberately left to a separate exact-cone
certificate.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

from scripts.search_apery_fermi_chamber_folds import _exact_nonnegative_dependence
from scripts.search_jacobian_kummer_chamber_folds import (
    ROOT,
    _matrix_payload,
    _rational_text,
)
from scripts.shimada_disc11_atlas import ShimadaDisc11Atlas
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia
from scripts.verify_k3_quotient_nef_slices import _coxeter_label


SOURCE = Path("/private/tmp/threediscs/disc16.txt")
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_d5_h14_independent_screen.json"
)
CHAMBER_INDEX = 5
AMBIENT_RANK = 20
FIXED_RANK = 15


def _row_product(row: tuple[int, ...], matrix: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    if len(row) != len(matrix):
        raise ValueError("incompatible row and matrix dimensions")
    return tuple(
        sum(row[index] * matrix[index][column] for index in range(len(row)))
        for column in range(len(matrix[0]))
    )


def _source_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _select_representative(atlas: ShimadaDisc11Atlas) -> tuple[tuple[int, ...], ...]:
    identity = sp.eye(AMBIENT_RANK)
    candidates = []
    for raw in atlas.stabilizer_group(CHAMBER_INDEX):
        matrix = sp.Matrix(raw)
        fixed_rank = AMBIENT_RANK - (matrix - identity).rank()
        if fixed_rank == FIXED_RANK:
            candidates.append(raw)
    if len(candidates) != 1:
        raise ArithmeticError(
            "D5 does not have a unique stabilizer element of fixed rank 15"
        )
    return candidates[0]


def _involution_orbits(permutation: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    unseen = set(range(len(permutation)))
    result: list[tuple[int, ...]] = []
    while unseen:
        start = min(unseen)
        orbit = tuple(sorted({start, permutation[start]}))
        if any(permutation[index] not in orbit for index in orbit):
            raise ArithmeticError("the wall action is not an involution")
        unseen.difference_update(orbit)
        result.append(orbit)
    return tuple(result)


def build_screen(path: Path = SOURCE) -> dict[str, Any]:
    atlas = ShimadaDisc11Atlas(path)
    chamber = atlas.chamber(CHAMBER_INDEX)
    gram = sp.Matrix(atlas.gram)
    gram_inverse = gram.inv()
    if not (
        gram.rows == gram.cols == AMBIENT_RANK
        and gram.det() == -16
        and exact_inertia(gram) == (1, 19, 0)
        and chamber.order_aut == 20
        and chamber.wall_count == 165
    ):
        raise ArithmeticError("the primary disc16/D5 data have changed")

    raw = _select_representative(atlas)
    raw_matrix = sp.Matrix(raw)
    dual = atlas.dual_action(raw)
    dual_matrix = sp.Matrix(dual)
    if not (
        raw_matrix**2 == sp.eye(AMBIENT_RANK)
        and raw_matrix * gram * raw_matrix.T == gram
        and raw_matrix.trace() == 10
        and AMBIENT_RANK - (raw_matrix - sp.eye(AMBIENT_RANK)).rank() == FIXED_RANK
        and dual_matrix == gram_inverse * raw_matrix * gram
        and dual_matrix**2 == sp.eye(AMBIENT_RANK)
        and dual_matrix * gram_inverse * dual_matrix.T == gram_inverse
    ):
        raise ArithmeticError("the intrinsically selected D5 element failed replay")

    walls = tuple(wall for orbit in chamber.wall_orbits for wall in orbit)
    if len(set(walls)) != len(walls) != 0:
        raise ArithmeticError("D5 contains duplicate primitive walls")
    wall_index = {wall: index for index, wall in enumerate(walls)}
    images = tuple(_row_product(wall, dual) for wall in walls)
    if any(image not in wall_index for image in images):
        raise ArithmeticError("the selected element does not preserve all D5 walls")
    permutation = tuple(wall_index[image] for image in images)
    if len(set(permutation)) != len(walls):
        raise ArithmeticError("the induced D5 wall action is not a permutation")
    wall_orbits = _involution_orbits(permutation)

    wall_rows = tuple(sp.Matrix(1, AMBIENT_RANK, wall) for wall in walls)
    orbit_normals = tuple(
        sum((wall_rows[index] for index in orbit), sp.zeros(1, AMBIENT_RANK))
        for orbit in wall_orbits
    )
    if any(normal * dual_matrix != normal for normal in orbit_normals):
        raise ArithmeticError("an orbit-sum wall is not fixed by the dual action")
    stacked = sp.Matrix.vstack(*orbit_normals)
    if stacked.rank() != FIXED_RANK:
        raise ArithmeticError("the orbit sums do not span the complete fixed space")
    orbit_gram = stacked * gram_inverse * stacked.T
    inertia = exact_inertia(orbit_gram)
    if inertia != (1, 14, len(wall_orbits) - FIXED_RANK):
        raise ArithmeticError("the fixed orbit-wall span is not H14 Lorentzian")

    ample_tuple = _row_product(chamber.weyl, atlas.projection_to_s_dual)
    ample = sp.Matrix(1, AMBIENT_RANK, ample_tuple)
    if not (
        ample * dual_matrix == ample
        and (ample * gram_inverse * ample.T)[0] > 0
        and all((ample * gram_inverse * normal.T)[0] > 0 for normal in orbit_normals)
    ):
        raise ArithmeticError("the published Weyl vector is not a strict fixed interior point")

    _, pivot_columns = stacked.rref()
    coordinates = stacked[:, list(pivot_columns)]
    facets: list[int] = []
    facet_witnesses: list[dict[str, Any]] = []
    unresolved: list[int] = []
    for target, normal in enumerate(orbit_normals):
        norm = sp.Rational(orbit_gram[target, target])
        if norm >= 0:
            unresolved.append(target)
            continue
        ample_pairing = sp.Rational((ample * gram_inverse * normal.T)[0])
        witness = ample - ample_pairing / norm * normal
        witness_square = sp.Rational((witness * gram_inverse * witness.T)[0])
        pairings = tuple(
            sp.Rational(value)
            for value in list(witness * gram_inverse * stacked.T)
        )
        if (
            witness_square > 0
            and pairings[target] == 0
            and all(
                value > 0
                for index, value in enumerate(pairings)
                if index != target
            )
        ):
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
        else:
            unresolved.append(target)

    redundancy_certificates: list[dict[str, Any]] = []
    for target in unresolved:
        dependence = _exact_nonnegative_dependence(coordinates, target)
        if dependence is None:
            raise ArithmeticError(
                f"orbit inequality {target + 1} is neither witnessed nor redundant"
            )
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
    if len(facets) + len(redundancy_certificates) != len(wall_orbits):
        raise ArithmeticError("the independent wall reduction is not exhaustive")

    good_labels: Counter[str] = Counter()
    bad_ratios: Counter[str] = Counter()
    bad_pairs: list[dict[str, Any]] = []
    for left_position, left in enumerate(facets):
        for right_position in range(left_position + 1, len(facets)):
            right = facets[right_position]
            left_norm = sp.Rational(orbit_gram[left, left])
            right_norm = sp.Rational(orbit_gram[right, right])
            product = sp.Rational(orbit_gram[left, right])
            ratio = sp.factor(product**2 / (left_norm * right_norm))
            try:
                label = _coxeter_label(left_norm, right_norm, product)
            except ValueError as error:
                ratio_text = _rational_text(sp.Rational(ratio))
                bad_ratios[ratio_text] += 1
                bad_pairs.append(
                    {
                        "facet_positions": [left_position + 1, right_position + 1],
                        "orbit_indices": [left + 1, right + 1],
                        "left_square": _rational_text(left_norm),
                        "right_square": _rational_text(right_norm),
                        "inner_product": _rational_text(product),
                        "normalized_squared_product": ratio_text,
                        "reason": str(error),
                    }
                )
            else:
                good_labels[label] += 1

    facet_gram = orbit_gram.extract(facets, facets)
    return {
        "schema": "shimada_disc16_d5_h14_independent_screen_v1",
        "status": {
            "primary_data_reconstructed": True,
            "representative_intrinsically_selected": True,
            "exact_irredundant_facets": True,
            "exact_complete_pair_gate": True,
            "all_facet_pairs_are_coxeter": not bad_pairs,
            "finite_volume_claimed": False,
            "decision": (
                "coxeter_candidate_requires_volume_audit"
                if not bad_pairs
                else "not_a_coxeter_polytope"
            ),
        },
        "source": {
            "archive": "threediscs.zip",
            "archive_url": "https://home.hiroshima-u.ac.jp/ichiro-shimada/preprints/AlgoAutK3/threediscs.zip",
            "member": "threediscs/disc16.txt",
            "local_path": str(path),
            "sha256": _source_sha256(path),
            "assignments": [
                "GramMatS",
                "prMatSdual",
                "chamberdata[5][weyl]",
                "chamberdata[5][orderaut]",
                "chamberdata[5][generatorsaut]",
                "chamberdata[5][wallorbit][*]",
            ],
        },
        "ambient_lattice": {
            "rank": AMBIENT_RANK,
            "determinant": int(gram.det()),
            "inertia": list(exact_inertia(gram)),
            "gram": _matrix_payload(gram),
            "dual_gram": _matrix_payload(gram_inverse),
        },
        "chamber": {
            "index": CHAMBER_INDEX,
            "published_name": "D5",
            "published_stabilizer_order": chamber.order_aut,
            "published_wall_orbit_sizes": [len(orbit) for orbit in chamber.wall_orbits],
            "published_wall_count": len(walls),
            "published_primitive_dual_walls": [list(wall) for wall in walls],
            "ample_dual_coordinates": list(ample_tuple),
            "ample_square": _rational_text(
                sp.Rational((ample * gram_inverse * ample.T)[0])
            ),
        },
        "representative": {
            "selection": "unique element of the full archived D5 stabilizer with fixed rank 15",
            "raw_primal_row_matrix": _matrix_payload(raw_matrix),
            "dual_row_matrix": _matrix_payload(dual_matrix),
            "order": 2,
            "trace": int(raw_matrix.trace()),
            "fixed_rank": FIXED_RANK,
            "hyperbolic_dimension": FIXED_RANK - 1,
        },
        "wall_action_permutation": [index + 1 for index in permutation],
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
            [source + 1 for source in wall_orbits[index]] for index in facets
        ],
        "redundant_orbit_count": len(redundancy_certificates),
        "redundancy_certificates": redundancy_certificates,
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
        "non_coxeter_pairs": bad_pairs,
        "facet_gram_sign_convention": "negative-square orbit-sum dual wall normals",
        "facet_gram": _matrix_payload(facet_gram),
        "positive_wall_gram": _matrix_payload(-facet_gram),
        "checks": {
            "full_archived_order_20_stabilizer_enumerated": True,
            "unique_rank_15_element_selected": True,
            "representative_is_exact_involution_and_isometry": True,
            "all_165_primitive_walls_reconstructed": True,
            "representative_permutes_all_165_walls": True,
            "orbit_sums_span_the_complete_fixed_space": True,
            "every_retained_wall_has_a_strict_timelike_witness": True,
            "every_removed_wall_has_an_exact_positive_conic_dependence": True,
            "every_irredundant_facet_pair_checked_exactly": True,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("disc16", type=Path, nargs="?", default=SOURCE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    result = build_screen(arguments.disc16)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file() or arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archived artifact differs: {arguments.output}")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
    print(
        json.dumps(
            {
                "dimension": result["representative"]["hyperbolic_dimension"],
                "wall_orbit_count": result["wall_orbit_count"],
                "facet_count": result["facet_count"],
                "redundant_orbit_count": result["redundant_orbit_count"],
                "coxeter_edge_counts": result["coxeter_edge_counts"],
                "non_coxeter_pair_count": result["non_coxeter_pair_count"],
                "decision": result["status"]["decision"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
