#!/usr/bin/env python3
"""Exact first reduction of a parked rank-ten Alexeev fixed section.

The all-subgroup denominator archived this class below its old ``H^10``
threshold.  It has fourteen raw wall orbits, but that does not imply fourteen
hyperbolic facets.  This script reconstructs the subgroup, saturated fixed
lattice, folded inequalities, exact redundancies, and Coxeter pair gate.  It
is exploratory: finite volume and a geometric cohomology lift are separate
later gates.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import shutil
import sympy as sp

from scripts.certify_alexeev_high_dim_reflection_examples import _extract_matrix
from scripts.certify_alexeev_20_2_1_15_9_1_exact_cones import (
    _analyze_case,
    _parse_cdd_rays,
    _replay_case,
    _run_cdd,
)
from scripts.inventory_alexeev_20_2_1_15_9_1_wall_symmetries import (
    SOURCE,
    _ambient_reconstruction,
    _generated_subgroup,
)
from scripts.screen_alexeev_20_2_1_15_9_1_fixed_sections import (
    _matrix,
    _screen_case,
)
from scripts.certify_shimada_disc16_d5_h10_30_vinberg_faces import (
    DiagramEnumerator,
)


ROOT = Path(__file__).resolve().parents[1]
DENOMINATOR = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_ten_small_zero_high_denominator.json"
)
DEFAULT_KEY = (14, 10, 1)
DEFAULT_CLASS_INDEX = 1


def _time_orientation(
    case: dict[str, object],
    cone: dict[str, object],
    rays: tuple[tuple[int, ...], ...],
) -> dict[str, str | bool]:
    """Check that every causal ray lies on the strict interior's time sheet."""

    gram = cone["fixed_gram"]
    interior = sp.Matrix(
        [
            sp.Rational(value)
            for value in case["strict_fixed_interior"]["fixed_lattice_coordinates"]
        ]
    )
    pairings = tuple(
        sp.Rational((interior.T * gram * sp.Matrix(ray))[0]) for ray in rays
    )
    if not pairings or any(value <= 0 for value in pairings):
        raise ArithmeticError("an extreme ray lies outside the chosen time cone")
    return {
        "all_extreme_rays_same_time_cone": True,
        "minimum_strict_interior_pairing": str(min(pairings)),
        "maximum_strict_interior_pairing": str(max(pairings)),
    }


def _literal_vinberg(
    case: dict[str, object],
    cone: dict[str, object],
    rays: tuple[tuple[int, ...], ...],
) -> dict[str, int | bool]:
    """Replay Vinberg's edge/end-point and affine-neighbourhood criteria."""

    dimension = int(cone["dimension"])
    inequality = cone["inequality"]
    fixed_gram = cone["fixed_gram"]
    facet_gram = _matrix(case["restricted_mirrors"]["facet_gram"])
    positive_gram = -facet_gram
    facet_count = facet_gram.rows
    endpoint_masks: list[int] = []
    endpoint_kinds: list[str] = []
    for ray in rays:
        vector = sp.Matrix(ray)
        slacks = inequality * vector
        active = tuple(index for index, value in enumerate(slacks) if value == 0)
        square = (vector.T * fixed_gram * vector)[0]
        endpoint_masks.append(sum(1 << index for index in active))
        endpoint_kinds.append("finite" if square > 0 else "ideal" if square == 0 else "hyperideal")
    if "hyperideal" in endpoint_kinds:
        raise ArithmeticError("literal Vinberg replay received a hyperideal endpoint")
    finite_masks = {
        mask
        for mask, kind in zip(endpoint_masks, endpoint_kinds, strict=True)
        if kind == "finite"
    }
    ideal_masks = {
        mask
        for mask, kind in zip(endpoint_masks, endpoint_kinds, strict=True)
        if kind == "ideal"
    }
    enumerator = DiagramEnumerator(positive_gram, dimension=dimension)
    finite_edges = enumerator.finite_subdiagrams(dimension - 1)
    finite_vertices = enumerator.finite_subdiagrams(dimension)
    if set(finite_vertices) != finite_masks:
        raise ArithmeticError("finite endpoint diagrams disagree with cone rays")
    edge_endpoints: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for right in range(len(endpoint_masks)):
        for left in range(right):
            intersection = endpoint_masks[left] & endpoint_masks[right]
            if intersection.bit_count() >= dimension - 1:
                edge_endpoints[intersection].append((left, right))
    if not (
        set(edge_endpoints) == set(finite_edges)
        and all(mask.bit_count() == dimension - 1 for mask in edge_endpoints)
        and all(len(pairs) == 1 for pairs in edge_endpoints.values())
    ):
        raise ArithmeticError("an elliptic edge lacks exactly two endpoints")
    affine = enumerator.connected_affine_subdiagrams()
    affine_neighborhoods: dict[int, int] = {}
    for component_mask in affine:
        component = tuple(
            index for index in range(facet_count) if component_mask & (1 << index)
        )
        zero_neighborhood = tuple(
            facet
            for facet in range(facet_count)
            if facet not in component
            and all(positive_gram[facet, node] == 0 for node in component)
        )
        neighborhood = component_mask | sum(1 << index for index in zero_neighborhood)
        if neighborhood not in ideal_masks:
            raise ArithmeticError("an affine zero-neighborhood is not an ideal endpoint")
        affine_neighborhoods[component_mask] = neighborhood
    for cusp_mask in ideal_masks:
        components = enumerator._components(cusp_mask, all_nonorthogonal=True)
        if not (
            sum(component.bit_count() - 1 for component in components)
            == dimension - 1
            and all(component in affine for component in components)
            and all(affine_neighborhoods[component] == cusp_mask for component in components)
        ):
            raise ArithmeticError("an ideal endpoint fails the full affine-neighbourhood test")
    return {
        "literal_vinberg_passes": True,
        "elliptic_edge_count": len(finite_edges),
        "finite_endpoint_count": len(finite_masks),
        "ideal_endpoint_count": len(ideal_masks),
        "connected_affine_component_count": len(affine),
        "bad_edge_count": 0,
        "bad_parabolic_neighborhood_count": 0,
    }


def build(
    key: tuple[int, int, int] = DEFAULT_KEY,
    class_index: int = DEFAULT_CLASS_INDEX,
    denominator_path: Path = DENOMINATOR,
) -> dict[str, object]:
    denominator = json.loads(denominator_path.read_text(encoding="utf-8"))
    examples = denominator.get("examples")
    if examples is None:
        examples = denominator["inventory"]["examples"]
    example = next(
        row for row in examples if tuple(row["key"]) == key
    )
    archived = next(
        row
        for row in example["subgroup_classes"]
        if int(row["class_index"]) == class_index
    )
    generators = tuple(
        tuple(map(int, permutation))
        for permutation in archived["representative_generators"]
    )
    subgroup = _generated_subgroup(generators)
    class_record = dict(archived)
    class_record["representative_elements"] = [
        list(permutation) for permutation in sorted(subgroup)
    ]
    source_text = SOURCE.read_text(encoding="utf-8")
    wall_gram = _extract_matrix(source_text, key)
    # Most archived Alexeev bases use the dual of the leading roots as a
    # strict chamber interior.  The (16,14,1) basis does not; the sum of all
    # chamber roots is instead an exact strict timelike interior.  Detect
    # this basis issue intrinsically so a valid fixed section is not dropped.
    ambient_gram, wall_rows = _ambient_reconstruction(
        wall_gram, int(example["ambient"]["rank"])
    )
    default_interior = sp.ones(1, ambient_gram.rows) * ambient_gram.inv()
    default_pairings = default_interior * ambient_gram * wall_rows.T
    source_interior = None
    if not (
        (default_interior * ambient_gram * default_interior.T)[0] > 0
        and all(value > 0 for value in default_pairings)
    ):
        root_sum = sp.ones(1, wall_rows.rows) * wall_rows
        root_sum_pairings = root_sum * ambient_gram * wall_rows.T
        if not (
            (root_sum * ambient_gram * root_sum.T)[0] > 0
            and all(value > 0 for value in root_sum_pairings)
        ):
            raise ArithmeticError("no certified strict source interior was found")
        source_interior = root_sum
    case = _screen_case(
        key=key,
        example=example,
        class_record=class_record,
        wall_gram=wall_gram,
        source_interior=source_interior,
        minimum_fixed_rank=10,
    )
    return case


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--key",
        default=",".join(map(str, DEFAULT_KEY)),
        help="Alexeev lattice key r,a,delta",
    )
    parser.add_argument("--class-index", type=int, default=DEFAULT_CLASS_INDEX)
    parser.add_argument("--denominator", type=Path, default=DENOMINATOR)
    parser.add_argument(
        "--volume",
        action="store_true",
        help="run exact GMP double description and classify all extreme rays",
    )
    parser.add_argument(
        "--vinberg",
        action="store_true",
        help="also replay the literal Vinberg edge/cusp conditions",
    )
    args = parser.parse_args()
    key = tuple(map(int, args.key.split(",")))
    if len(key) != 3:
        raise ValueError("--key must be r,a,delta")
    case = build(key, args.class_index, args.denominator)
    reduction = case["exact_reduction"]
    summary = {
                "id": case["id"],
                "fixed_rank": case["case"]["fixed_rank"],
                "hyperbolic_dimension": case["case"]["hyperbolic_dimension"],
                "raw_wall_orbit_count": case["case"]["wall_orbit_count"],
                "active_hyperbolic_facet_count": reduction[
                    "active_hyperbolic_facet_count"
                ],
                "conic_redundant_count": reduction["conic_redundant_count"],
                "causal_inactive_count": len(
                    reduction["causally_inactive_orbit_indices"]
                ),
                "unresolved_orbit_indices": reduction["unresolved_orbit_indices"],
                "classification": case["decision"]["classification"],
                "coxeter_obstruction_count": case["coxeter_gate"][
                    "non_coxeter_pair_count"
                ],
                "all_active_reflections_integral": case["restricted_mirrors"][
                    "all_active_reflections_integral_on_saturated_fixed_lattice"
                ],
                "finite_volume_checked": False,
                "geometric_action_lift_checked": False,
            }
    if args.volume or args.vinberg:
        cone = _replay_case(case)
        executable = shutil.which("cddexec_gmp") or "cddexec_gmp"
        output = _run_cdd(cone["inequality_text"], executable)
        rays = _parse_cdd_rays(output, int(cone["rank"]))
        volume, _ = _analyze_case(case, cone, rays)
        time_orientation = _time_orientation(case, cone, rays)
        exact = volume["exact_double_description"]
        summary.update(
            {
                "finite_volume_checked": True,
                "extreme_ray_count": exact["extreme_ray_count"],
                "finite_vertex_count": exact["finite_vertex_count"],
                "ideal_vertex_count": exact["ideal_vertex_count"],
                "hyperideal_vertex_count": exact["hyperideal_vertex_count"],
                "volume_classification": volume["decision"]["classification"],
            }
        )
        summary.update(time_orientation)
        if args.vinberg:
            summary.update(_literal_vinberg(case, cone, rays))
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
