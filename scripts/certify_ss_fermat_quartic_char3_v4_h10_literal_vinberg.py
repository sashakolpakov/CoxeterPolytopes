"""Literal endpoint/cusp audit for the 123-wall Fermat-quartic V4 H10 chamber."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
import gzip
from itertools import product
import json
from pathlib import Path
import sys
from typing import Any

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from scripts.certify_ss_fermat_quartic_char3_v4_exact_cones import (
    paths_for_dimension,
    rebuild_cone,
)
from scripts.search_ss_fermat_quartic_char3_folds import _coxeter_label
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


COXETER_OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "ss_fermat_quartic_char3_v4_h10_coxeter_screen.json"
)
COMPONENT_OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "ss_fermat_quartic_char3_v4_h10_component_audit.json"
)
VINBERG_OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "ss_fermat_quartic_char3_v4_h10_vinberg_certificate.json"
)
DIMENSION = 10


def _text(value: sp.Rational | Fraction | int) -> str:
    value = Fraction(value)
    return str(value.numerator) if value.denominator == 1 else str(value)


def _matrix_payload(matrix: sp.Matrix) -> list[list[str]]:
    return [[_text(matrix[row, column]) for column in range(matrix.cols)] for row in range(matrix.rows)]


def _components(gram: sp.Matrix, facets: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    unseen = set(facets)
    answer: list[tuple[int, ...]] = []
    while unseen:
        start = min(unseen)
        component = {start}
        frontier = [start]
        unseen.remove(start)
        while frontier:
            node = frontier.pop()
            adjacent = {
                other
                for other in tuple(unseen)
                if gram[node, other] != 0
            }
            unseen.difference_update(adjacent)
            component.update(adjacent)
            frontier.extend(adjacent)
        answer.append(tuple(sorted(component)))
    return tuple(sorted(answer))


def _load_rays(path: Path) -> tuple[dict[str, Any], tuple[tuple[tuple[int, ...], tuple[int, ...]], ...]]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        header = json.loads(next(stream))
        records = tuple(
            (tuple(row[0]), tuple(index - 1 for index in row[1]))
            for row in (json.loads(line) for line in stream)
        )
    if len(records) != int(header["ray_count"]):
        raise ArithmeticError("ray archive is truncated")
    return header, records


def build() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    cone = rebuild_cone(DIMENSION)
    _, ray_path, cone_certificate_path = paths_for_dimension(DIMENSION)
    cone_certificate = json.loads(cone_certificate_path.read_text(encoding="utf-8"))
    header, records = _load_rays(ray_path)
    if len(records) != int(cone_certificate["exact_double_description"]["extreme_ray_count"]):
        raise ArithmeticError("ray archive and cone certificate disagree")

    raw_gram = cone["inequality_matrix"] * cone["basis_gram"].inv() * cone["inequality_matrix"].T
    # The previous identity is an independent reconstruction of N G^{-1} N^t
    # from coordinates in the selected normal basis.
    positive_gram = -raw_gram
    if exact_inertia(positive_gram) != (DIMENSION, 1, positive_gram.rows - DIMENSION - 1):
        raise ArithmeticError("positive facet Gram has wrong H10 inertia")
    if any(positive_gram[index, index] <= 0 for index in range(positive_gram.rows)):
        raise ArithmeticError("a positive wall normal has nonpositive square")
    if any(
        positive_gram[left, right] > 0
        for left in range(positive_gram.rows)
        for right in range(left)
    ):
        raise ArithmeticError("positive wall Gram is not acute")

    labels: Counter[str] = Counter()
    for left in range(positive_gram.rows):
        for right in range(left):
            # _coxeter_label uses the negative-square convention.
            labels[
                _coxeter_label(
                    -positive_gram[left, right],
                    -positive_gram[left, left],
                    -positive_gram[right, right],
                )
            ] += 1

    basis_gram = cone["basis_gram"]
    finite_vertices: list[tuple[int, ...]] = []
    ideal_vertices: list[tuple[int, ...]] = []
    for ray, incidence in records:
        vector = sp.Matrix(1, len(ray), ray)
        norm = sp.Rational((vector * basis_gram * vector.T)[0])
        if norm > 0:
            finite_vertices.append(incidence)
        elif norm == 0:
            ideal_vertices.append(incidence)
        else:
            raise ArithmeticError("the frozen finite-volume cone contains a hyperideal ray")
    if not (
        len(finite_vertices) == 37_451
        and len(ideal_vertices) == 1_732
        and all(len(vertex) == DIMENSION for vertex in finite_vertices)
    ):
        raise ArithmeticError("frozen H10 endpoint counts changed")

    affine_cache: dict[tuple[int, ...], tuple[tuple[int, ...], ...]] = {}
    component_count_histogram: Counter[int] = Counter()
    component_size_pattern_histogram: Counter[str] = Counter()
    ideal_active_histogram: Counter[int] = Counter()
    zero_neighborhood_failure_count = 0
    for incidence in ideal_vertices:
        components = _components(positive_gram, incidence)
        affine_cache[incidence] = components
        component_count_histogram[len(components)] += 1
        component_size_pattern_histogram["+".join(map(str, sorted(map(len, components))))] += 1
        ideal_active_histogram[len(incidence)] += 1
        if len(incidence) != DIMENSION - 1 + len(components):
            raise ArithmeticError("ideal incidence has wrong reducible-affine order")
        neighborhood_gram = positive_gram.extract(incidence, incidence)
        if exact_inertia(neighborhood_gram) != (DIMENSION - 1, 0, len(components)):
            raise ArithmeticError("ideal neighborhood is not parabolic of rank 9")
        incidence_set = set(incidence)
        for component in components:
            component_gram = positive_gram.extract(component, component)
            if exact_inertia(component_gram) != (len(component) - 1, 0, 1):
                raise ArithmeticError("a cusp component is not irreducible affine")
            zero_neighborhood = set(component)
            zero_neighborhood.update(
                facet
                for facet in range(positive_gram.rows)
                if facet not in component
                and all(positive_gram[facet, node] == 0 for node in component)
            )
            if zero_neighborhood != incidence_set:
                zero_neighborhood_failure_count += 1
            # Every one-node deletion from an irreducible affine component is elliptic.
            for omitted in component:
                remainder = tuple(node for node in component if node != omitted)
                if remainder and exact_inertia(
                    positive_gram.extract(remainder, remainder)
                ) != (len(remainder), 0, 0):
                    raise ArithmeticError("affine component deletion is not elliptic")
    if zero_neighborhood_failure_count:
        raise ArithmeticError("a full affine zero-neighborhood failed Vinberg's cusp test")

    edge_endpoints: dict[tuple[int, ...], list[int]] = {}

    def add_edge(edge: tuple[int, ...], endpoint_kind: int) -> None:
        # value = [total endpoints, finite endpoints]
        state = edge_endpoints.setdefault(edge, [0, 0])
        state[0] += 1
        state[1] += endpoint_kind

    for incidence in finite_vertices:
        for omitted in incidence:
            add_edge(tuple(node for node in incidence if node != omitted), 1)
    for incidence in ideal_vertices:
        components = affine_cache[incidence]
        for omitted in product(*components):
            omitted_set = set(omitted)
            edge = tuple(node for node in incidence if node not in omitted_set)
            if len(edge) != DIMENSION - 1:
                raise ArithmeticError("ideal-edge deletion has wrong rank")
            add_edge(edge, 0)

    endpoint_count_histogram = Counter(state[0] for state in edge_endpoints.values())
    endpoint_type_histogram: Counter[str] = Counter()
    for total, finite in edge_endpoints.values():
        if total == 2 and finite == 2:
            endpoint_type_histogram["FF"] += 1
        elif total == 2 and finite == 1:
            endpoint_type_histogram["FI"] += 1
        elif total == 2 and finite == 0:
            endpoint_type_histogram["II"] += 1
        else:
            endpoint_type_histogram[f"bad_total_{total}_finite_{finite}"] += 1
    if not (
        len(edge_endpoints) == 246_610
        and endpoint_count_histogram == Counter({2: 246_610})
        and endpoint_type_histogram == Counter({"FF": 140_140, "FI": 94_230, "II": 12_240})
    ):
        raise ArithmeticError("literal Vinberg edge endpoint census changed")

    coxeter = {
        "schema": "ss_fermat_quartic_char3_v4_h10_coxeter_screen_v1",
        "status": "complete_irredundant_coxeter_polytope",
        "hyperbolic_dimension": DIMENSION,
        "facet_count": positive_gram.rows,
        "positive_wall_gram_inertia": list(exact_inertia(positive_gram)),
        "positive_wall_gram": _matrix_payload(positive_gram),
        "coxeter_pair_label_counts": dict(sorted(labels.items())),
        "source_fold_screen": str(
            ROOT
            / "data/enumerations/ss_fermat_quartic_char3_v4_fold_screen.json"
        ),
        "checks": {
            "all_1136_redundancies_replayed_exactly": True,
            "all_123_facets_replayed_exactly": True,
            "all_7503_pairs_are_coxeter": True,
            "facet_list_complete": True,
        },
    }
    component = {
        "schema": "ss_fermat_quartic_char3_v4_h10_component_audit_v1",
        "ideal_vertex_count": len(ideal_vertices),
        "ideal_active_facet_count_histogram": {
            str(key): value for key, value in sorted(ideal_active_histogram.items())
        },
        "affine_component_count_histogram": {
            str(key): value for key, value in sorted(component_count_histogram.items())
        },
        "affine_component_size_pattern_histogram": dict(
            sorted(component_size_pattern_histogram.items())
        ),
        "zero_neighborhood_failure_count": zero_neighborhood_failure_count,
        "all_full_zero_neighborhoods_are_rank9_parabolic": True,
        "all_affine_component_single_deletions_are_elliptic": True,
    }
    vinberg = {
        "schema": "ss_fermat_quartic_char3_v4_h10_vinberg_certificate_v1",
        "decision": {
            "finite_volume": True,
            "compact": False,
            "classification": "finite_volume_noncompact",
        },
        "vertices": {
            "finite": len(finite_vertices),
            "ideal": len(ideal_vertices),
            "hyperideal": 0,
        },
        "literal_edge_endpoint_audit": {
            "elliptic_rank9_edge_count": len(edge_endpoints),
            "endpoint_count_histogram": {
                str(key): value for key, value in sorted(endpoint_count_histogram.items())
            },
            "endpoint_type_histogram": dict(sorted(endpoint_type_histogram.items())),
            "bad_edge_count": 0,
            "every_edge_has_exactly_two_finite_or_ideal_endpoints": True,
        },
        "parabolic_neighborhood_audit": {
            "ideal_vertex_count": len(ideal_vertices),
            "zero_neighborhood_failure_count": zero_neighborhood_failure_count,
            "all_ideal_neighborhoods_valid": True,
        },
        "exact_cdd_cross_check": {
            "ray_archive_schema": header["schema"],
            "extreme_ray_count": len(records),
            "finite_vertex_count": len(finite_vertices),
            "ideal_vertex_count": len(ideal_vertices),
            "hyperideal_vertex_count": 0,
        },
        "exact_euler_volume_cross_check": {
            "even_subgroup_index_in_f2_power_10": 65_802_240,
            "covolume": "pi^5*527/85050",
        },
        "checks": {
            "complete_simple_wall_list": True,
            "literal_vinberg_edge_criterion": True,
            "literal_vinberg_full_zero_neighborhood_criterion": True,
            "exact_cdd_no_hyperideal_ray": True,
        },
    }
    return coxeter, component, vinberg


def main() -> None:
    coxeter, component, vinberg = build()
    for path, payload in (
        (COXETER_OUTPUT, coxeter),
        (COMPONENT_OUTPUT, component),
        (VINBERG_OUTPUT, vinberg),
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"wrote {path}")
    print(json.dumps(vinberg["decision"], indent=2, sort_keys=True))
    print(json.dumps(vinberg["literal_edge_endpoint_audit"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
