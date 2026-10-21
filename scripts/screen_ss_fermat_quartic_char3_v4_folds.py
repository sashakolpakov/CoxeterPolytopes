"""Exact Coxeter gate for the two high-rank V4 Fermat-quartic folds.

The representatives are frozen by
``scratch_fermat_char3_full_stabilizer.g`` in the exact full-stabilizer
census JSON.  For each V4 subgroup this script derives its action on every
one of the 5,944 published chamber walls, forms primitive orbit-sum normals
in the common fixed space, certifies strict facets by timelike relative-
interior witnesses, and tests every pair of certified strict facets against
the Coxeter cosine-square list.
"""

from __future__ import annotations

import argparse
from collections import Counter, deque
from fractions import Fraction
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from scripts.search_ss_fermat_quartic_char3_folds import (
    SOURCE,
    _coxeter_label,
    _fraction_text,
    _primitive_positive,
    parse_source,
)
from scripts.search_apery_fermi_chamber_folds import _exact_nonnegative_dependence
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


REPRESENTATIVES = (
    ROOT
    / "data"
    / "enumerations"
    / "ss_fermat_quartic_char3_full_stabilizer_exact.json"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "ss_fermat_quartic_char3_v4_fold_screen.json"
)


def load_json_tolerating_gap_wrap(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8").replace("\\\n", ""))


def subgroup_orbits(
    generator_permutations: tuple[tuple[int, ...], ...], size: int
) -> tuple[tuple[int, ...], ...]:
    unseen = set(range(size))
    answer: list[tuple[int, ...]] = []
    while unseen:
        start = min(unseen)
        orbit = {start}
        queue: deque[int] = deque([start])
        while queue:
            current = queue.popleft()
            for permutation in generator_permutations:
                image = permutation[current]
                if image not in orbit:
                    orbit.add(image)
                    queue.append(image)
        unseen.difference_update(orbit)
        answer.append(tuple(sorted(orbit)))
    return tuple(answer)


def screen_cast(
    source: dict[str, Any], record: dict[str, Any], cast_index: int
) -> dict[str, Any]:
    gram: sp.Matrix = source["gram"]
    ample_matrix: sp.Matrix = source["ample"]
    walls: tuple[tuple[int, ...], ...] = source["walls"]
    wall_types: tuple[str, ...] = source["wall_types"]
    identity = sp.eye(gram.rows)
    generators = tuple(
        sp.Matrix(matrix) for matrix in record["generator_archived_ns_matrices"]
    )
    if len(generators) != 2:
        raise ArithmeticError("a frozen V4 representative must have two generators")
    for generator in generators:
        if generator * generator != identity:
            raise ArithmeticError("a frozen V4 generator is not involutory")
        if generator * gram * generator.T != gram:
            raise ArithmeticError("a frozen V4 generator is not a lattice isometry")
        if ample_matrix * generator != ample_matrix:
            raise ArithmeticError("a frozen V4 generator does not fix FQh0")
    if generators[0] * generators[1] != generators[1] * generators[0]:
        raise ArithmeticError("the frozen generators do not commute")
    if generators[0] * generators[1] == identity:
        raise ArithmeticError("the frozen generators are not independent")

    stacked = generators[0] - identity
    for generator in generators[1:]:
        stacked = stacked.col_join(generator - identity)
    fixed_rank = gram.rows - stacked.rank()
    if fixed_rank != int(record["common_fixed_rank"]):
        raise ArithmeticError("frozen V4 fixed rank changed")

    wall_index = {wall: index for index, wall in enumerate(walls)}
    if len(wall_index) != len(walls):
        raise ArithmeticError("published chamber walls are not distinct")
    generator_permutations: list[tuple[int, ...]] = []
    for generator in generators:
        covector_action = generator.inv().T
        if any(not value.is_Integer for value in covector_action):
            raise ArithmeticError("dual action is not integral")
        covector_action_np = np.array(covector_action, dtype=np.int64)
        permutation: list[int] = []
        for wall in walls:
            image = tuple(int(value) for value in np.array(wall) @ covector_action_np)
            if image not in wall_index:
                raise ArithmeticError("V4 generator does not permute all published walls")
            permutation.append(wall_index[image])
        frozen = tuple(permutation)
        if any(frozen[frozen[index]] != index for index in range(len(walls))):
            raise ArithmeticError("derived wall action is not involutory")
        if any(wall_types[index] != wall_types[image] for index, image in enumerate(frozen)):
            raise ArithmeticError("V4 generator mixes published wall types")
        generator_permutations.append(frozen)

    orbits = subgroup_orbits(tuple(generator_permutations), len(walls))
    orbit_size_histogram = Counter(len(orbit) for orbit in orbits)
    if any(size not in (1, 2, 4) for size in orbit_size_histogram):
        raise ArithmeticError("unexpected V4 wall-orbit size")

    ample = np.array([int(value) for value in ample_matrix], dtype=np.int64)
    ray_to_source_orbits: dict[tuple[int, ...], list[int]] = {}
    orbit_sums: list[tuple[int, ...]] = []
    for orbit_index, orbit in enumerate(orbits):
        row = tuple(sum(walls[index][column] for index in orbit) for column in range(gram.rows))
        primitive = _primitive_positive(row, ample)
        if primitive not in ray_to_source_orbits:
            orbit_sums.append(primitive)
            ray_to_source_orbits[primitive] = []
        ray_to_source_orbits[primitive].append(orbit_index)

    inverse = gram.inv()
    denominator = int(sp.ilcm(*[sp.Rational(value).q for value in inverse]))
    metric_numerator = np.array(inverse * denominator, dtype=np.int64)
    restricted = np.array(orbit_sums, dtype=np.int64)
    restricted_metric_rows = restricted @ metric_numerator
    norm_numerators = np.sum(restricted_metric_rows * restricted, axis=1)
    ample_pairings = restricted @ ample

    restricted_sympy = sp.Matrix(restricted.tolist())
    restricted_rank = restricted_sympy.rank()
    if restricted_rank != fixed_rank:
        raise ArithmeticError("restricted wall normals do not span the common fixed dual space")
    _, basis_indices = restricted_sympy.T.rref()
    restricted_basis = restricted_sympy[list(basis_indices), :]
    restricted_basis_gram = restricted_basis * inverse * restricted_basis.T
    nonzero_inertia = exact_inertia(restricted_basis_gram)
    inertia = (
        nonzero_inertia[0],
        nonzero_inertia[1],
        len(orbit_sums) - fixed_rank,
    )
    expected_inertia = (1, fixed_rank - 1, len(orbit_sums) - fixed_rank)
    if inertia != expected_inertia:
        raise ArithmeticError(f"unexpected restricted inertia {inertia}")

    strict_facets: list[int] = []
    facet_certificates: list[dict[str, Any]] = []
    ample_square = int((ample_matrix * gram * ample_matrix.T)[0])
    for index, raw_norm in enumerate(norm_numerators):
        norm = int(raw_norm)
        if norm >= 0:
            continue
        gram_row = restricted_metric_rows[index] @ restricted.T
        margins = int(ample_pairings[index]) * gram_row - norm * ample_pairings
        if int(margins[index]) != 0:
            raise ArithmeticError("strict-witness self margin is nonzero")
        other_margins = np.delete(margins, index)
        if other_margins.size and int(other_margins.min()) > 0:
            strict_facets.append(index)
            witness_square = Fraction(
                ample_square * norm - int(ample_pairings[index]) ** 2 * denominator,
                norm,
            )
            if witness_square <= 0:
                raise ArithmeticError("strict facet witness is not timelike")
            facet_certificates.append(
                {
                    "restricted_ray_index": index + 1,
                    "primitive_normal_covector": list(orbit_sums[index]),
                    "source_wall_orbit_indices": [
                        source + 1 for source in ray_to_source_orbits[orbit_sums[index]]
                    ],
                    "source_wall_indices": sorted(
                        wall + 1
                        for source_orbit in ray_to_source_orbits[orbit_sums[index]]
                        for wall in orbits[source_orbit]
                    ),
                    "normal_square": _fraction_text(Fraction(norm, denominator)),
                    "ample_pairing": int(ample_pairings[index]),
                    "witness_square": _fraction_text(witness_square),
                    "minimum_integer_strict_margin": int(other_margins.min()),
                }
            )

    bad_pairs: list[dict[str, Any]] = []
    label_counts: Counter[str] = Counter()
    for left_position, left in enumerate(strict_facets):
        gram_row = restricted_metric_rows[left] @ restricted.T
        for right in strict_facets[left_position + 1 :]:
            product = int(gram_row[right])
            left_norm = int(norm_numerators[left])
            right_norm = int(norm_numerators[right])
            label = _coxeter_label(product, left_norm, right_norm)
            ratio = Fraction(product * product, left_norm * right_norm)
            if label is None:
                bad_pairs.append(
                    {
                        "restricted_ray_indices": [left + 1, right + 1],
                        "normal_product": _fraction_text(Fraction(product, denominator)),
                        "normal_squares": [
                            _fraction_text(Fraction(left_norm, denominator)),
                            _fraction_text(Fraction(right_norm, denominator)),
                        ],
                        "normalized_squared_product": _fraction_text(ratio),
                        "reason": (
                            "negative normal product (non-acute)"
                            if product < 0
                            else "not a Coxeter cosine square"
                        ),
                    }
                )
            else:
                label_counts[label] += 1

    global_reduction: dict[str, Any]
    if bad_pairs:
        global_reduction = {
            "performed": False,
            "reason": "local strict-facet gate already proved non-Coxeter",
        }
    else:
        _, coordinate_columns = restricted_sympy.rref()
        coordinates = restricted_sympy[:, list(coordinate_columns)]
        strict_set = set(strict_facets)
        redundancy_certificates: list[dict[str, Any]] = []
        unresolved_after_strict_cone: list[int] = []
        for target in range(len(orbit_sums)):
            if target in strict_set:
                continue
            local_indices = [*strict_facets, target]
            dependence = _exact_nonnegative_dependence(
                coordinates[local_indices, :], len(strict_facets)
            )
            if dependence is None:
                unresolved_after_strict_cone.append(target)
                continue
            redundancy_certificates.append(
                {
                    "restricted_ray_index": target + 1,
                    "terms": [
                        {
                            "restricted_ray_index": strict_facets[support] + 1,
                            "coefficient": _fraction_text(Fraction(coefficient)),
                        }
                        for support, coefficient in dependence
                    ],
                }
            )
        global_reduction = {
            "performed": True,
            "candidate_inequality_count": len(orbit_sums),
            "irredundant_facet_count": len(strict_facets),
            "redundant_inequality_count": len(redundancy_certificates),
            "redundancy_certificates": redundancy_certificates,
            "unresolved_after_strict_facet_cone_count": len(
                unresolved_after_strict_cone
            ),
            "unresolved_restricted_ray_indices": [
                index + 1 for index in unresolved_after_strict_cone
            ],
            "exact_complete": not unresolved_after_strict_cone,
        }

    type_size_histogram: Counter[str] = Counter()
    for orbit in orbits:
        type_size_histogram[f"{wall_types[orbit[0]]}:size_{len(orbit)}"] += 1
    return {
        "cast_index": cast_index,
        "subgroup": {
            "abstract_type": "C2xC2",
            "order": 4,
            "common_fixed_rank": fixed_rank,
            "hyperbolic_dimension": fixed_rank - 1,
            "involution_fixed_rank_profile": record["profile"],
            "conjugacy_class_size": int(record["conjugacy_class_size"]),
            "normalizer_order": int(record["normalizer_order"]),
            "centralizer_order": int(record["centralizer_order"]),
            "generator_ns_matrices": [
                [[int(value) for value in row] for row in matrix.tolist()]
                for matrix in generators
            ],
        },
        "walls": {
            "published_count": len(walls),
            "subgroup_orbit_count": len(orbits),
            "orbit_size_histogram": {
                str(size): count for size, count in sorted(orbit_size_histogram.items())
            },
            "orbit_type_size_histogram": dict(sorted(type_size_histogram.items())),
            "distinct_restricted_ray_count": len(orbit_sums),
            "coincident_restricted_ray_count": len(orbits) - len(orbit_sums),
            "restricted_ray_rank": restricted_rank,
            "restricted_gram_inertia": list(inertia),
            "metric_denominator": denominator,
            "subgroup_wall_orbits": [
                [index + 1 for index in orbit] for orbit in orbits
            ],
            "primitive_restricted_normal_covectors": [
                list(normal) for normal in orbit_sums
            ],
        },
        "strict_facet_gate": {
            "certified_strict_facet_count": len(strict_facets),
            "certificates": facet_certificates,
            "coxeter_compatible_pair_label_counts": dict(sorted(label_counts.items())),
            "bad_pair_count": len(bad_pairs),
            "bad_pairs": bad_pairs,
            "decision": (
                "not_coxeter"
                if bad_pairs
                else "coxeter_polytope_requires_volume_audit"
                if global_reduction.get("exact_complete") is True
                else "survives_local_gate_requires_global_facet_reduction"
            ),
        },
        "global_facet_reduction": global_reduction,
        "checks": {
            "exact_v4_relations": True,
            "generators_preserve_lattice_and_ample_class": True,
            "generators_permute_all_published_walls": True,
            "restricted_walls_span_common_fixed_space": restricted_rank == fixed_rank,
            "restricted_space_is_lorentzian": inertia == expected_inertia,
            "each_listed_facet_has_strict_timelike_relative_interior_witness": True,
        },
    }


def build_artifact(source_path: Path, representatives_path: Path) -> dict[str, Any]:
    source = parse_source(source_path)
    census = load_json_tolerating_gap_wrap(representatives_path)
    records = census["v4_high_fixed_rank_records"]
    if sorted(int(record["common_fixed_rank"]) for record in records) != [11, 12]:
        raise ArithmeticError("expected exactly the H10 and H11 V4 representatives")
    casts = [screen_cast(source, record, index + 1) for index, record in enumerate(records)]
    return {
        "schema": "ss_fermat_quartic_char3_v4_fold_screen_v1",
        "source": {
            "paper": "Kondo--Shimada, The automorphism group of a supersingular K3 surface with Artin invariant 1 in characteristic 3",
            "paper_url": "https://arxiv.org/abs/1205.6520",
            "official_data_assignment_count": 5944,
            "wall_assignments": ["FQWW112", "FQWW648", "FQWW5184"],
            "representative_census": str(representatives_path.relative_to(ROOT)),
        },
        "scope": {
            "artifact_stops_before_volume_gate": True,
            "full_v4_high_fixed_rank_conjugacy_class_count": len(records),
            "all_published_walls_screened_per_cast": True,
            "strict_facet_gate_is_decisive_when_bad_pair_count_is_positive": True,
            "global_facet_reduction_required_only_for_local_survivors": True,
        },
        "casts": casts,
        "summary": {
            "not_coxeter_count": sum(
                cast["strict_facet_gate"]["decision"] == "not_coxeter" for cast in casts
            ),
            "complete_coxeter_candidate_count": sum(
                cast["strict_facet_gate"]["decision"]
                == "coxeter_polytope_requires_volume_audit"
                for cast in casts
            ),
            "unresolved_local_survivor_count": sum(
                cast["strict_facet_gate"]["decision"]
                == "survives_local_gate_requires_global_facet_reduction"
                for cast in casts
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--representatives", type=Path, default=REPRESENTATIVES)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    artifact = build_artifact(args.source, args.representatives)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    for cast in artifact["casts"]:
        gate = cast["strict_facet_gate"]
        print(
            f"H{cast['subgroup']['hyperbolic_dimension']}: "
            f"{cast['walls']['subgroup_orbit_count']} wall orbits, "
            f"{gate['certified_strict_facet_count']} strict facets, "
            f"{gate['bad_pair_count']} forbidden pairs -> {gate['decision']}"
        )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
