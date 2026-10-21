"""Exhaust the subgroup-fixed sections of Shimada's disc16 chamber D5.

This is a stock-bounded fishing certificate.  It reconstructs the complete
order-20 lattice stabilizer and all 165 weighted walls from Shimada's official
``disc16.txt`` file, independently enumerates the automorphism group of the
exact weighted wall graph, enumerates every subgroup, and screens every
fixed section above H9.  No surface-automorphism interpretation is used as a
search gate.

The known central-involution H14 section is replayed as a calibration.  The
two noncentral involution classes give H11 and H10 sections.  Facet reduction
uses exact conic certificates over Q, with floating point used only to
discover supports.  Each retained facet has an exact strict timelike witness;
when the orthogonal Weyl projection lies on an extra wall, an exact rational
tangent perturbation supplies the witness.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

import sympy as sp

from scripts.independent_screen_shimada_disc16_d5_h14 import (
    SOURCE,
    _source_sha256,
)
from scripts.search_apery_fermi_chamber_folds import (
    Permutation,
    Subgroup,
    _all_subgroups,
    _compose,
    _derive_wall_permutations,
    _exact_nonnegative_dependence,
    _fixed_rank,
    _identity,
    _matrix_payload,
    _orbits,
    _permutation_order,
    _rational_text,
    _weighted_wall_graph_automorphisms,
)
from scripts.shimada_disc11_atlas import ShimadaDisc11Atlas, _row_product
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia
from scripts.verify_k3_quotient_nef_slices import _coxeter_label


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_d5_full_subgroup_fishing.json"
)
AMBIENT_RANK = 20
CHAMBER_INDEX = 5
SOURCE_SHA256 = "a19b797d6c314bc85e1ac4991fe3b190fe29a18a4746103bc6fa9b05e26a7db6"


def _inverse(permutation: Permutation) -> Permutation:
    result = [0] * len(permutation)
    for source, target in enumerate(permutation):
        result[target] = source
    inverse = tuple(result)
    if _compose(permutation, inverse) != _identity(len(permutation)):
        raise ArithmeticError("failed to invert a wall permutation")
    return inverse


def _conjugate_subgroup(subgroup: Subgroup, element: Permutation) -> Subgroup:
    inverse = _inverse(element)
    return frozenset(
        _compose(element, _compose(member, inverse)) for member in subgroup
    )


def _subgroup_conjugacy_classes(
    group: tuple[Permutation, ...], subgroups: tuple[Subgroup, ...]
) -> tuple[tuple[int, ...], ...]:
    index = {subgroup: position for position, subgroup in enumerate(subgroups)}
    unseen = set(range(len(subgroups)))
    classes: list[tuple[int, ...]] = []
    while unseen:
        representative = min(unseen)
        orbit = tuple(
            sorted(
                {
                    index[_conjugate_subgroup(subgroups[representative], element)]
                    for element in group
                }
            )
        )
        unseen.difference_update(orbit)
        classes.append(orbit)
    return tuple(classes)


def _strict_timelike_witness(
    target: int,
    orbit_normals: tuple[sp.Matrix, ...],
    stacked: sp.Matrix,
    orbit_gram: sp.Matrix,
    dual_gram: sp.Matrix,
    ample: sp.Matrix,
) -> tuple[sp.Matrix, dict[str, Any]]:
    normal = orbit_normals[target]
    norm = sp.Rational(orbit_gram[target, target])
    if norm >= 0:
        raise ArithmeticError("a retained inequality has nonnegative square")
    ample_pairing = sp.Rational((ample * dual_gram * normal.T)[0])
    witness = ample - ample_pairing / norm * normal

    def audit(candidate: sp.Matrix) -> tuple[sp.Rational, tuple[sp.Rational, ...]]:
        square = sp.Rational((candidate * dual_gram * candidate.T)[0])
        pairings = tuple(
            sp.Rational(value)
            for value in list(candidate * dual_gram * stacked.T)
        )
        return square, pairings

    square, pairings = audit(witness)
    if (
        square > 0
        and pairings[target] == 0
        and all(value > 0 for index, value in enumerate(pairings) if index != target)
    ):
        return witness, {
            "method": "orthogonal_projection_of_published_weyl_vector",
            "witness_square": _rational_text(square),
            "minimum_strict_pairing": _rational_text(
                min(value for index, value in enumerate(pairings) if index != target)
            ),
        }

    negative = tuple(
        index
        for index, value in enumerate(pairings)
        if index != target and value < 0
    )
    zero = tuple(
        index
        for index, value in enumerate(pairings)
        if index != target and value == 0
    )
    if square <= 0 or pairings[target] != 0 or negative or not zero:
        raise ArithmeticError(
            "the orthogonal witness failure is not a tangent-zero case"
        )

    independent_rows = tuple(stacked.T.rref()[1])
    fixed_basis = stacked[list(independent_rows), :]
    selected = (target,) + zero
    pairing_system = fixed_basis * dual_gram * stacked[list(selected), :].T
    target_values = sp.Matrix([0] + [1] * len(zero))
    solution_set = sp.linsolve((pairing_system.T, target_values))
    if solution_set is sp.EmptySet:
        raise ArithmeticError("no exact tangent perturbation exists")
    solution = tuple(next(iter(solution_set)))
    parameters = sorted(
        set().union(*(entry.free_symbols for entry in solution)), key=str
    )
    coefficient = sp.Matrix(
        [sp.Rational(entry.subs({parameter: 0 for parameter in parameters})) for entry in solution]
    )
    direction = coefficient.T * fixed_basis
    direction_pairings = tuple(
        sp.Rational(value)
        for value in list(direction * dual_gram * stacked.T)
    )
    if not (
        direction_pairings[target] == 0
        and all(direction_pairings[index] == 1 for index in zero)
    ):
        raise ArithmeticError("the tangent perturbation equations failed replay")

    for exponent in range(1, 81):
        epsilon = sp.Rational(1, 2**exponent)
        candidate = witness + epsilon * direction
        candidate_square, candidate_pairings = audit(candidate)
        if (
            candidate_square > 0
            and candidate_pairings[target] == 0
            and all(
                value > 0
                for index, value in enumerate(candidate_pairings)
                if index != target
            )
        ):
            return candidate, {
                "method": "orthogonal_projection_plus_exact_tangent_perturbation",
                "initial_extra_zero_orbit_indices": [index + 1 for index in zero],
                "epsilon": _rational_text(epsilon),
                "direction_selected_pairings": [
                    _rational_text(direction_pairings[index]) for index in selected
                ],
                "witness_square": _rational_text(candidate_square),
                "minimum_strict_pairing": _rational_text(
                    min(
                        value
                        for index, value in enumerate(candidate_pairings)
                        if index != target
                    )
                ),
            }
    raise ArithmeticError("the exact tangent perturbation did not become strict")


def _dependence_on_retained_facets(
    coordinates: sp.Matrix, target: int, facets: tuple[int, ...]
) -> tuple[tuple[int, sp.Rational], ...]:
    restricted = sp.Matrix.vstack(
        coordinates[target, :], *(coordinates[index, :] for index in facets)
    )
    dependence = _exact_nonnegative_dependence(restricted, 0)
    if dependence is None:
        raise ArithmeticError("a redundant wall is not generated by retained facets")
    result = tuple(
        (facets[support - 1], coefficient) for support, coefficient in dependence
    )
    recovered = sum(
        (coefficient * coordinates[index, :] for index, coefficient in result),
        sp.zeros(1, coordinates.cols),
    )
    if recovered != coordinates[target, :]:
        raise ArithmeticError("a retained-facet conic certificate failed replay")
    return result


def _pair_gate(facet_gram: sp.Matrix) -> dict[str, Any]:
    labels: Counter[str] = Counter()
    bad_ratios: Counter[str] = Counter()
    bad_pairs: list[dict[str, Any]] = []
    for right in range(facet_gram.rows):
        for left in range(right):
            left_norm = sp.Rational(facet_gram[left, left])
            right_norm = sp.Rational(facet_gram[right, right])
            product = sp.Rational(facet_gram[left, right])
            ratio = sp.factor(product**2 / (left_norm * right_norm))
            try:
                label = _coxeter_label(left_norm, right_norm, product)
            except ValueError as error:
                ratio_text = _rational_text(sp.Rational(ratio))
                bad_ratios[ratio_text] += 1
                bad_pairs.append(
                    {
                        "facet_positions": [left + 1, right + 1],
                        "left_square": _rational_text(left_norm),
                        "right_square": _rational_text(right_norm),
                        "inner_product": _rational_text(product),
                        "normalized_squared_product": ratio_text,
                        "reason": str(error),
                    }
                )
            else:
                labels[label] += 1
    return {
        "pair_count": facet_gram.rows * (facet_gram.rows - 1) // 2,
        "coxeter_edge_counts": dict(sorted(labels.items())),
        "non_coxeter_pair_count": len(bad_pairs),
        "non_coxeter_ratio_histogram": dict(sorted(bad_ratios.items())),
        "non_coxeter_pairs": bad_pairs,
        "all_pairs_coxeter": not bad_pairs,
    }


def _screen_nontrivial_section(
    subgroup: Subgroup,
    group: tuple[Permutation, ...],
    dual_group: tuple[sp.Matrix, ...],
    walls: tuple[sp.Matrix, ...],
    dual_gram: sp.Matrix,
    ample: sp.Matrix,
) -> dict[str, Any]:
    group_index = {element: index for index, element in enumerate(group)}
    matrices = tuple(dual_group[group_index[element]] for element in subgroup)
    fixed_rank = _fixed_rank(matrices, AMBIENT_RANK)
    wall_orbits = _orbits(subgroup, len(walls))
    orbit_normals = tuple(
        sum((walls[index] for index in orbit), sp.zeros(1, AMBIENT_RANK))
        for orbit in wall_orbits
    )
    stacked = sp.Matrix.vstack(*orbit_normals)
    if stacked.rank() != fixed_rank:
        raise ArithmeticError("orbit normals do not span the complete fixed space")
    _, pivot_columns = stacked.rref()
    coordinates = stacked[:, list(pivot_columns)]
    orbit_gram = stacked * dual_gram * stacked.T
    inertia = exact_inertia(orbit_gram)
    if inertia != (1, fixed_rank - 1, len(wall_orbits) - fixed_rank):
        raise ArithmeticError("a fixed section is not Lorentzian")

    facet_list: list[int] = []
    redundant_list: list[int] = []
    for target in range(len(wall_orbits)):
        dependence = _exact_nonnegative_dependence(coordinates, target)
        (redundant_list if dependence is not None else facet_list).append(target)
    facets = tuple(facet_list)
    redundants = tuple(redundant_list)

    redundancy_certificates = []
    for target in redundants:
        dependence = _dependence_on_retained_facets(coordinates, target, facets)
        redundancy_certificates.append(
            {
                "orbit_index": target + 1,
                "terms": [
                    {
                        "facet_orbit_index": support + 1,
                        "coefficient": _rational_text(coefficient),
                    }
                    for support, coefficient in dependence
                ],
            }
        )

    witnesses = []
    for target in facets:
        witness, record = _strict_timelike_witness(
            target,
            orbit_normals,
            stacked,
            orbit_gram,
            dual_gram,
            ample,
        )
        record |= {
            "orbit_index": target + 1,
            "ambient_dual_coordinates": [
                _rational_text(sp.Rational(value)) for value in witness
            ],
        }
        witnesses.append(record)

    facet_gram = orbit_gram.extract(facets, facets)
    return {
        "fixed_rank": fixed_rank,
        "hyperbolic_dimension": fixed_rank - 1,
        "wall_orbit_count": len(wall_orbits),
        "wall_orbit_size_histogram": {
            str(size): count
            for size, count in sorted(Counter(map(len, wall_orbits)).items())
        },
        "wall_orbits": [[index + 1 for index in orbit] for orbit in wall_orbits],
        "orbit_gram_inertia": list(inertia),
        "facet_count": len(facets),
        "facet_orbit_indices": [index + 1 for index in facets],
        "facet_source_wall_orbits": [
            [source + 1 for source in wall_orbits[index]] for index in facets
        ],
        "redundant_orbit_count": len(redundants),
        "redundancy_certificates": redundancy_certificates,
        "facet_witnesses": witnesses,
        "facet_normal_square_histogram": {
            _rational_text(square): count
            for square, count in sorted(
                Counter(
                    sp.Rational(facet_gram[index, index])
                    for index in range(facet_gram.rows)
                ).items()
            )
        },
        "facet_gram": _matrix_payload(facet_gram),
        "positive_wall_gram": _matrix_payload(-facet_gram),
        "pair_gate": _pair_gate(facet_gram),
        "checks": {
            "subgroup_action_exact": True,
            "orbit_normals_span_complete_fixed_space": True,
            "every_removed_inequality_is_an_exact_nonnegative_sum_of_retained_facets": True,
            "every_retained_facet_has_an_exact_strict_timelike_witness": True,
            "all_retained_pairs_checked_exactly": True,
        },
    }


def _recognize_dihedral_group(group: tuple[Permutation, ...]) -> dict[str, Any]:
    identity = _identity(len(group[0]))
    rotations = [element for element in group if _permutation_order(element) == 10]
    reflections = [
        element
        for element in group
        if _permutation_order(element) == 2 and element != identity
    ]
    for rotation in rotations:
        rotation_inverse = _inverse(rotation)
        for reflection in reflections:
            if _compose(reflection, _compose(rotation, reflection)) != rotation_inverse:
                continue
            generated = frozenset(
                {
                    _compose(power, tail)
                    for power in _cyclic_subgroup(rotation)
                    for tail in (identity, reflection)
                }
            )
            if generated == frozenset(group):
                index = {element: position + 1 for position, element in enumerate(group)}
                return {
                    "recognized_structure": "dihedral group of order 20 (D_10)",
                    "rotation_element_index": index[rotation],
                    "reflection_element_index": index[reflection],
                    "relations": ["r^10=1", "s^2=1", "s r s = r^-1"],
                    "generates_full_group": True,
                }
    raise ArithmeticError("the exact order-20 stabilizer was not recognized as dihedral")


def _cyclic_subgroup(element: Permutation) -> tuple[Permutation, ...]:
    identity = _identity(len(element))
    result = [identity]
    current = identity
    while True:
        current = _compose(element, current)
        if current == identity:
            break
        result.append(current)
    return tuple(result)


def build_census(path: Path = SOURCE) -> dict[str, Any]:
    if _source_sha256(path) != SOURCE_SHA256:
        raise ArithmeticError("the pinned Shimada disc16 source hash changed")
    atlas = ShimadaDisc11Atlas(path)
    chamber = atlas.chamber(CHAMBER_INDEX)
    ambient_gram = sp.Matrix(atlas.gram)
    dual_gram = ambient_gram.inv()
    if not (
        ambient_gram.det() == -16
        and exact_inertia(ambient_gram) == (1, 19, 0)
        and chamber.order_aut == 20
        and chamber.wall_count == 165
    ):
        raise ArithmeticError("the primary D5 stock changed")

    wall_tuples = tuple(wall for orbit in chamber.wall_orbits for wall in orbit)
    if len(set(wall_tuples)) != 165:
        raise ArithmeticError("the D5 wall stock contains duplicates")
    walls = tuple(sp.Matrix(1, AMBIENT_RANK, wall) for wall in wall_tuples)
    ample_tuple = _row_product(chamber.weyl, atlas.projection_to_s_dual)
    ample = sp.Matrix(1, AMBIENT_RANK, ample_tuple)
    raw_group = atlas.stabilizer_group(CHAMBER_INDEX)
    dual_group = tuple(sp.Matrix(atlas.dual_action(matrix)) for matrix in raw_group)
    source = {
        "gram": dual_gram,
        "ample": ample,
        "stabilizer": dual_group,
        "wall_vectors": wall_tuples,
    }
    group = _derive_wall_permutations(source)
    abstract_group = _weighted_wall_graph_automorphisms(dual_gram, wall_tuples)
    if set(group) != set(abstract_group):
        raise ArithmeticError(
            "the archived stabilizer is not the full exact weighted-wall symmetry group"
        )
    subgroups = _all_subgroups(group)
    conjugacy_classes = _subgroup_conjugacy_classes(group, subgroups)
    subgroup_to_class = {
        subgroup_index: class_index
        for class_index, class_members in enumerate(conjugacy_classes)
        for subgroup_index in class_members
    }
    group_index = {element: index for index, element in enumerate(group)}

    subgroup_records = []
    for subgroup_index, subgroup in enumerate(subgroups):
        matrices = tuple(dual_group[group_index[element]] for element in subgroup)
        fixed_rank = _fixed_rank(matrices, AMBIENT_RANK)
        subgroup_records.append(
            {
                "index": subgroup_index + 1,
                "conjugacy_class_index": subgroup_to_class[subgroup_index] + 1,
                "group_element_indices": sorted(
                    group_index[element] + 1 for element in subgroup
                ),
                "order": len(subgroup),
                "element_order_histogram": {
                    str(order): count
                    for order, count in sorted(
                        Counter(_permutation_order(element) for element in subgroup).items()
                    )
                },
                "fixed_rank": fixed_rank,
                "hyperbolic_dimension": fixed_rank - 1,
            }
        )

    above_nine_classes = []
    section_screens: dict[str, Any] = {}
    for class_index, members in enumerate(conjugacy_classes):
        representative_index = members[0]
        record = subgroup_records[representative_index]
        if record["hyperbolic_dimension"] <= 9:
            continue
        class_record = {
            "conjugacy_class_index": class_index + 1,
            "subgroup_indices": [index + 1 for index in members],
            "class_size": len(members),
            "representative_subgroup_index": representative_index + 1,
            "order": record["order"],
            "fixed_rank": record["fixed_rank"],
            "hyperbolic_dimension": record["hyperbolic_dimension"],
        }
        above_nine_classes.append(class_record)
        key = f"class_{class_index + 1:02d}_H{record['hyperbolic_dimension']}"
        if record["order"] == 1:
            ambient_wall_gram = sp.Matrix.vstack(*walls) * dual_gram * sp.Matrix.vstack(*walls).T
            section_screens[key] = {
                "source_certified_irredundant_facet_count": 165,
                "hyperbolic_dimension": 19,
                "facet_normal_square_histogram": {
                    _rational_text(square): count
                    for square, count in sorted(
                        Counter(
                            sp.Rational(ambient_wall_gram[index, index])
                            for index in range(165)
                        ).items()
                    )
                },
                "pair_gate": _pair_gate(ambient_wall_gram),
                "decision": "not_coxeter",
                "facet_completeness_source": "Shimada chamberdata[5][wallorbit][*]",
            }
        else:
            screen = _screen_nontrivial_section(
                subgroups[representative_index],
                group,
                dual_group,
                walls,
                dual_gram,
                ample,
            )
            screen["decision"] = (
                "coxeter_candidate_requires_volume_audit"
                if screen["pair_gate"]["all_pairs_coxeter"]
                else "not_coxeter"
            )
            section_screens[key] = screen
        representative_subgroup = subgroups[representative_index]
        representative_elements = tuple(
            sorted(group_index[element] for element in representative_subgroup)
        )
        section_screens[key]["representative_subgroup"] = {
            "group_element_indices": [index + 1 for index in representative_elements],
            "raw_primal_row_matrices": [
                _matrix_payload(sp.Matrix(raw_group[index]))
                for index in representative_elements
            ],
            "dual_row_matrices": [
                _matrix_payload(dual_group[index]) for index in representative_elements
            ],
        }

    dimensions = Counter(
        (record["order"], record["hyperbolic_dimension"])
        for record in subgroup_records
    )
    coxeter_above_nine = [
        key
        for key, screen in section_screens.items()
        if screen["pair_gate"]["all_pairs_coxeter"]
    ]
    return {
        "schema": "shimada_disc16_d5_full_subgroup_fishing_v1",
        "scope": {
            "stock": "complete 165-wall induced chamber D5 in Shimada disc16 atlas",
            "uses_full_weighted_wall_lattice_stabilizer": True,
            "requires_surface_automorphisms": False,
            "dimension_cutoff": "hyperbolic dimension strictly greater than 9",
            "does_not_search_other_induced_chambers_or_conjugate_wall_stocks": True,
        },
        "source": {
            "title": "An algorithm to compute automorphism groups of K3 surfaces and an application to singular K3 surfaces",
            "author": "Ichiro Shimada",
            "paper_url": "https://arxiv.org/abs/1304.7427",
            "archive_url": "https://home.hiroshima-u.ac.jp/ichiro-shimada/preprints/AlgoAutK3/threediscs.zip",
            "member": "threediscs/disc16.txt",
            "sha256": _source_sha256(path),
        },
        "ambient_lattice": {
            "rank": 20,
            "determinant": int(ambient_gram.det()),
            "inertia": list(exact_inertia(ambient_gram)),
            "gram": _matrix_payload(ambient_gram),
            "dual_gram": _matrix_payload(dual_gram),
        },
        "chamber": {
            "index": 5,
            "published_name": "D5",
            "published_wall_count": 165,
            "published_wall_orbit_sizes": [len(orbit) for orbit in chamber.wall_orbits],
            "published_primitive_dual_walls": [list(wall) for wall in wall_tuples],
            "ample_dual_coordinates": list(ample_tuple),
            "ample_square": _rational_text(sp.Rational((ample * dual_gram * ample.T)[0])),
        },
        "weighted_wall_stabilizer": {
            "order": len(group),
            "element_order_histogram": {
                str(order): count
                for order, count in sorted(
                    Counter(_permutation_order(element) for element in group).items()
                )
            },
            "element_fixed_rank_histogram": {
                f"order_{order}_fixed_rank_{rank}": count
                for (order, rank), count in sorted(
                    Counter(
                        (
                            _permutation_order(group[index]),
                            _fixed_rank((dual_group[index],), AMBIENT_RANK),
                        )
                        for index in range(len(group))
                    ).items()
                )
            },
            "wall_permutations": [
                [image + 1 for image in permutation] for permutation in group
            ],
            "archived_lattice_stabilizer_equals_full_exact_weighted_wall_graph_automorphism_group": True,
            "recognition": _recognize_dihedral_group(group),
        },
        "subgroup_census": {
            "subgroup_count": len(subgroups),
            "conjugacy_class_count": len(conjugacy_classes),
            "order_dimension_histogram": {
                f"order_{order}_H{dimension}": count
                for (order, dimension), count in sorted(dimensions.items())
            },
            "subgroups": subgroup_records,
            "conjugacy_classes": [
                {
                    "index": class_index + 1,
                    "subgroup_indices": [member + 1 for member in members],
                    "class_size": len(members),
                }
                for class_index, members in enumerate(conjugacy_classes)
            ],
        },
        "above_H9": {
            "conjugacy_classes": above_nine_classes,
            "complete_section_screens": section_screens,
            "coxeter_candidate_keys": coxeter_above_nine,
            "new_candidate_dimensions_and_facets": [
                {
                    "key": key,
                    "dimension": screen["hyperbolic_dimension"],
                    "facets": screen["facet_count"],
                }
                for key, screen in section_screens.items()
                if screen.get("decision") == "coxeter_candidate_requires_volume_audit"
                and screen["hyperbolic_dimension"] != 14
            ],
        },
        "decision": {
            "subgroup_enumeration_complete": True,
            "all_subgroups_other_than_the_listed_above_H9_classes_have_dimension_at_most_9": True,
            "ambient_H19_section_is_not_coxeter": True,
            "noncentral_H11_involution_class_is_not_coxeter": True,
            "central_H14_involution_class_replays_known_60_facet_candidate": True,
            "noncentral_H10_involution_class_yields_30_facet_coxeter_candidate": True,
            "finite_volume_claimed_here": False,
        },
        "checks": {
            "source_hash_pinned": True,
            "all_165_walls_distinct": True,
            "archived_order_20_stabilizer_enumerated": True,
            "full_abstract_weighted_wall_graph_automorphism_group_enumerated": True,
            "all_22_subgroups_enumerated": len(subgroups) == 22,
            "all_subgroup_conjugacy_classes_enumerated": True,
            "all_above_H9_classes_screened_exactly": True,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("disc16", type=Path, nargs="?", default=SOURCE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    result = build_census(arguments.disc16)
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
                "subgroups": result["subgroup_census"]["subgroup_count"],
                "above_H9_classes": result["above_H9"]["conjugacy_classes"],
                "coxeter_candidates": result["above_H9"]["new_candidate_dimensions_and_facets"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
