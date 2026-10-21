#!/usr/bin/env python3
"""Certify the maximal disjoint doubling stocks in the JK-H12 chamber.

This verifier deliberately imports no Jacobian--Kummer production or scout
module.  It reparses Dolgachev--Shimada's published wall data, reconstructs
the trace-nine fixed lattice and its 172 facet roots, and then performs its
own exact graph enumeration.  All graph relations use integer Gram data.

Two relations are kept separate:

* disjoint: parallel or ultraparallel hyperplanes;
* strictly ultraparallel: disjoint even after ideal closure.

The distinction is the one used in Allcock's doubling theorems.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any, Callable, Iterable, Sequence

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form


ROOT = Path(__file__).resolve().parents[1]
PUBLISHED = (
    ROOT / "data/external/dolgachev_shimada_15_nodal_2019/15nodalcompdata.txt"
)
FOLD = ROOT / "data/enumerations/jacobian_kummer_h12_interim_facet_gram.json"
VINBERG = (
    ROOT / "data/enumerations/jacobian_kummer_h12_vinberg_volume_certificate.json"
)
EXACT_CONE = (
    ROOT / "data/enumerations/jacobian_kummer_h12_exact_cone_certificate.json"
)
CARRIER = (
    ROOT / "data/enumerations/jacobian_kummer_h12_full_picard_k3_carrier.json"
)
OUTPUT = (
    ROOT
    / "data/enumerations/"
    "jacobian_kummer_h12_maximal_disjoint_doubling_family.json"
)

SCHEMA = "coxeter.jacobian_kummer.h12_maximal_disjoint_doubling_family.v1"
STATUS_DATE = "2026-08-08"
PINNED_SHA256 = {
    "published_Dolgachev_Shimada_data": (
        "413d7363faf4574595e9c28329737eb4ad8010e0aebae53afb76762ce9f8966a"
    ),
    "H12_fold_orbits": (
        "53c39642b81bbd7f5be172408bfecb857bd4ed57467c91bdbaed63f9e734fb92"
    ),
    "literal_Vinberg_seed_gate": (
        "c3474e33120cab541205225dd8df31beea7ae8006cd1667c4fd66408f9b931f1"
    ),
    "exact_GMP_seed_gate": (
        "84d9660b5990d89fb4a78b73180fa815843ee7af46d4d8c56be73f51e37a7c0f"
    ),
    "full_Picard_K3_carrier": (
        "f2c57b9fb0423560c9e238541b0ac1a81ca141917c864929501ebdfa7f294ff2"
    ),
}

Vector = tuple[int, ...]
Matrix = tuple[tuple[int, ...], ...]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_sha256(payload: Any) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def _extract_assignment(text: str, name: str) -> str:
    marker = f"{name}:="
    start = text.index(marker) + len(marker)
    while text[start].isspace():
        start += 1
    if text[start] != "[":
        raise ArithmeticError(f"published assignment {name} is not a list")
    depth = 0
    quote: str | None = None
    for end in range(start, len(text)):
        character = text[end]
        if quote is not None:
            if character == quote and text[end - 1] != "\\":
                quote = None
            continue
        if character in {'"', "'"}:
            quote = character
        elif character == "[":
            depth += 1
        elif character == "]":
            depth -= 1
            if depth == 0:
                return text[start : end + 1]
    raise ArithmeticError(f"unterminated published assignment {name}")


def _content(values: Iterable[int]) -> int:
    entries = tuple(abs(int(value)) for value in values if int(value))
    return math.gcd(*entries) if entries else 0


def _primitive_column(vector: sp.Matrix) -> sp.Matrix:
    rationals = tuple(sp.Rational(value) for value in vector)
    denominator = math.lcm(*(int(value.q) for value in rationals))
    integers = tuple(int(value * denominator) for value in rationals)
    content = _content(integers)
    if not content:
        raise ArithmeticError("zero vector cannot be primitively normalized")
    return sp.Matrix([value // content for value in integers])


def _pair(gram: Matrix, left: Sequence[int], right: Sequence[int]) -> int:
    return sum(
        int(left[row]) * gram[row][column] * int(right[column])
        for row in range(len(gram))
        for column in range(len(gram))
    )


def _exact_inertia(matrix: sp.Matrix) -> tuple[int, int, int]:
    work = sp.Matrix(matrix)
    positive = negative = 0
    while work.rows:
        pivot = next((i for i in range(work.rows) if work[i, i]), None)
        if pivot is None:
            pair = next(
                (
                    (i, j)
                    for i in range(work.rows)
                    for j in range(i + 1, work.rows)
                    if work[i, j]
                ),
                None,
            )
            if pair is None:
                return positive, negative, work.rows
            order = list(pair) + [i for i in range(work.rows) if i not in pair]
            work = work.extract(order, order)
            change = sp.eye(work.rows)
            change[1, 0] = 1
            work = change.T * work * change
            pivot = 0
        if pivot:
            order = [pivot] + [i for i in range(work.rows) if i != pivot]
            work = work.extract(order, order)
        diagonal = sp.Rational(work[0, 0])
        positive += int(bool(diagonal > 0))
        negative += int(bool(diagonal < 0))
        if work.rows == 1:
            return positive, negative, 0
        column = work[1:, 0]
        work = sp.simplify(work[1:, 1:] - column * column.T / diagonal)
    return positive, negative, 0


def _parse_published() -> dict[str, Any]:
    text = PUBLISHED.read_text(encoding="utf-8").replace("\\\n", "")
    gram = sp.Matrix(ast.literal_eval(_extract_assignment(text, "GramS16")))
    ample = sp.Matrix(1, 17, ast.literal_eval(_extract_assignment(text, "ample16")))
    wall_text = _extract_assignment(text, "D16WallRecs")
    wall_numbers = tuple(
        int(value) for value in re.findall(r"\bno\s*:=\s*(\d+)", wall_text)
    )
    vector_texts = re.findall(r"\bvector\s*:=\s*\[([^\]]+)\]", wall_text)
    walls = tuple(
        sp.Matrix(
            1,
            17,
            [sp.Rational(entry.strip()) for entry in vector.split(",")],
        )
        for vector in vector_texts
    )
    if not (
        gram.shape == (17, 17)
        and gram == gram.T
        and ample.shape == (1, 17)
        and len(walls) == 316
        and wall_numbers == tuple(range(1, 317))
    ):
        raise ArithmeticError("published Jacobian--Kummer source parse failed")
    return {"gram": gram, "ample": ample, "walls": walls}


def _reconstruct_seed() -> dict[str, Any]:
    paths = {
        "published_Dolgachev_Shimada_data": PUBLISHED,
        "H12_fold_orbits": FOLD,
        "literal_Vinberg_seed_gate": VINBERG,
        "exact_GMP_seed_gate": EXACT_CONE,
        "full_Picard_K3_carrier": CARRIER,
    }
    hashes = {name: _sha256(path) for name, path in paths.items()}
    if hashes != PINNED_SHA256:
        raise ArithmeticError(f"a pinned input changed: {hashes}")

    source = _parse_published()
    fold = json.loads(FOLD.read_text(encoding="utf-8"))
    ambient_gram: sp.Matrix = source["gram"]
    ample: sp.Matrix = source["ample"]
    walls: tuple[sp.Matrix, ...] = source["walls"]
    involution = sp.Matrix(fold["group_class"]["representative_matrix"])

    wall_index = {tuple(wall): index for index, wall in enumerate(walls)}
    involution_permutation = tuple(
        wall_index[tuple(wall * involution)] for wall in walls
    )
    basis = sp.Matrix(
        [
            list(_primitive_column(vector))
            for vector in (involution.T - sp.eye(17)).nullspace()
        ]
    )
    fixed_gram_sp = sp.simplify(basis * ambient_gram * basis.T)
    fixed_gram: Matrix = tuple(
        tuple(int(fixed_gram_sp[row, column]) for column in range(13))
        for row in range(13)
    )
    # This displayed minor is an independent saturation witness.
    unit_minor_columns = (0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 13, 14)
    unit_minor = int(basis[:, unit_minor_columns].det())
    interior_column = sp.simplify(
        fixed_gram_sp.inv() * basis * ambient_gram * ample.T
    )
    if any(sp.Rational(value).q != 1 for value in interior_column):
        raise ArithmeticError("published ample vector is nonintegral in fixed basis")
    interior = tuple(int(value) for value in interior_column)

    roots: list[Vector] = []
    raw_sums: list[sp.Matrix] = []
    for orbit in fold["facet_source_wall_orbits"]:
        raw = sum((walls[int(index) - 1] for index in orbit), sp.zeros(1, 17))
        coordinates = _primitive_column(
            fixed_gram_sp.inv() * basis * ambient_gram * raw.T
        )
        root = tuple(int(value) for value in coordinates)
        if _pair(fixed_gram, interior, root) < 0:
            root = tuple(-value for value in root)
        roots.append(root)
        raw_sums.append(raw)

    raw_gram = sp.Matrix.vstack(*raw_sums) * ambient_gram * sp.Matrix.vstack(
        *raw_sums
    ).T
    checks = {
        "involution_is_order_two": involution**2 == sp.eye(17),
        "involution_preserves_published_Gram": (
            involution * ambient_gram * involution.T == ambient_gram
        ),
        "involution_fixes_published_ample": ample * involution == ample,
        "involution_permutes_all_316_published_walls": (
            len(set(involution_permutation)) == 316
        ),
        "fixed_basis_is_saturated_rank_13": basis.shape == (13, 17)
        and abs(unit_minor) == 1,
        "fixed_signature_is_1_12": _exact_inertia(fixed_gram_sp) == (1, 12, 0),
        "fixed_determinant_is_256": int(fixed_gram_sp.det()) == 256,
        "all_172_facet_roots_reconstructed": len(roots) == 172,
        "raw_172_Gram_matches_frozen_fold": raw_gram
        == sp.Matrix(fold["facet_gram"]),
        "all_roots_are_primitive_negative_and_oriented": all(
            _content(root) == 1
            and _pair(fixed_gram, root, root) < 0
            and _pair(fixed_gram, interior, root) > 0
            for root in roots
        ),
        "published_ample_is_timelike": _pair(
            fixed_gram, interior, interior
        )
        > 0,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"raw seed reconstruction failed: {checks}")
    return {
        "hashes": hashes,
        "fixed_gram": fixed_gram,
        "interior": interior,
        "roots": tuple(roots),
        "checks": checks,
    }


def _classify_relation(pairing: int, left_square: int, right_square: int) -> str:
    if pairing < 0 or left_square >= 0 or right_square >= 0:
        return "invalid"
    square = pairing * pairing
    product = left_square * right_square
    if pairing == 0:
        return "m2"
    if 4 * square == product:
        return "m3"
    if 2 * square == product:
        return "m4"
    if 4 * square == 3 * product:
        return "m6"
    if square == product:
        return "parallel"
    if square > product:
        return "ultraparallel"
    return "invalid"


def _reflection_is_integral_isometry(gram: Matrix, root: Vector) -> bool:
    size = len(root)
    square = _pair(gram, root, root)
    covector = tuple(
        sum(root[row] * gram[row][column] for row in range(size))
        for column in range(size)
    )
    entries: list[list[int]] = []
    for row in range(size):
        matrix_row: list[int] = []
        for column in range(size):
            value = sp.Rational(int(row == column)) - sp.Rational(
                2 * root[row] * covector[column], square
            )
            if value.q != 1:
                return False
            matrix_row.append(int(value))
        entries.append(matrix_row)
    reflection = sp.Matrix(entries)
    gram_sp = sp.Matrix(gram)
    return bool(
        reflection**2 == sp.eye(size)
        and reflection.T * gram_sp * reflection == gram_sp
        and reflection * sp.Matrix(root) == -sp.Matrix(root)
    )


def _doubling_stock(seed: dict[str, Any]) -> dict[str, Any]:
    gram: Matrix = seed["fixed_gram"]
    roots: tuple[Vector, ...] = seed["roots"]
    squares = tuple(_pair(gram, root, root) for root in roots)
    relations: list[tuple[str, ...]] = []
    doubling: list[int] = []
    whole_histogram: Counter[str] = Counter()
    for left, root in enumerate(roots):
        labels: list[str] = []
        for right, other in enumerate(roots):
            if left == right:
                labels.append("m1")
                continue
            label = _classify_relation(
                _pair(gram, root, other), squares[left], squares[right]
            )
            if label == "invalid":
                raise ArithmeticError(f"invalid Coxeter relation W{left+1},W{right+1}")
            labels.append(label)
            if right < left:
                whole_histogram[label] += 1
        relations.append(tuple(labels))
        # These data have no finite odd label other than m=3.
        if "m3" not in labels:
            doubling.append(left)

    stock_histogram = Counter(
        relations[left][right]
        for position, left in enumerate(doubling)
        for right in doubling[:position]
    )
    expected_doubling = tuple(
        [30, 31, 32, 35, 38, 41, 47, 50, 55] + list(range(57, 172))
    )
    checks = {
        "all_14706_seed_relations_are_exact_Coxeter_relations": (
            sum(whole_histogram.values()) == 14_706
        ),
        "whole_seed_relation_histogram_matches_independent_certificates": (
            whole_histogram
            == Counter(
                {
                    "m2": 2820,
                    "m3": 144,
                    "m4": 800,
                    "parallel": 1750,
                    "ultraparallel": 9192,
                }
            )
        ),
        "exactly_124_doubling_walls": len(doubling) == 124,
        "doubling_wall_labels_are_exact": tuple(doubling) == expected_doubling,
        "doubling_criterion_is_no_incident_odd_finite_label": all(
            "m3" not in relations[wall] for wall in doubling
        ),
        "all_doubling_roots_have_primitive_square_minus_four": all(
            squares[wall] == -4 and _content(roots[wall]) == 1
            for wall in doubling
        ),
        "all_124_doubling_reflections_are_integral_lattice_isometries": all(
            _reflection_is_integral_isometry(gram, roots[wall])
            for wall in doubling
        ),
        "doubling_stock_relation_histogram_is_exact": stock_histogram
        == Counter({"m2": 264, "parallel": 1138, "ultraparallel": 6224}),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"doubling-stock reconstruction failed: {checks}")
    return {
        "doubling": tuple(doubling),
        "relations": tuple(relations),
        "whole_histogram": dict(sorted(whole_histogram.items())),
        "stock_histogram": dict(sorted(stock_histogram.items())),
        "checks": checks,
    }


def _maximal_cliques(
    vertices: Sequence[int], adjacency: dict[int, set[int]]
) -> tuple[tuple[int, ...], ...]:
    """Deterministic Bron--Kerbosch enumeration with pivoting."""

    result: list[tuple[int, ...]] = []

    def visit(chosen: set[int], candidates: set[int], excluded: set[int]) -> None:
        if not candidates and not excluded:
            result.append(tuple(sorted(chosen)))
            return
        pool = candidates | excluded
        pivot = max(
            pool,
            key=lambda vertex: (len(candidates & adjacency[vertex]), -vertex),
            default=None,
        )
        branch = candidates - (adjacency[pivot] if pivot is not None else set())
        for vertex in tuple(sorted(branch)):
            if vertex not in candidates:
                continue
            visit(
                chosen | {vertex},
                candidates & adjacency[vertex],
                excluded & adjacency[vertex],
            )
            candidates.remove(vertex)
            excluded.add(vertex)

    visit(set(), set(vertices), set())
    return tuple(sorted(set(result), key=lambda clique: (len(clique), clique)))


DISJOINT_CONFLICT_MATCHING_ONE_BASED = (
    (31, 162),
    (32, 169),
    (33, 170),
    (36, 130),
    (39, 141),
    (42, 145),
    (48, 161),
    (51, 159),
    (56, 164),
    (58, 172),
    (59, 171),
    (60, 168),
)
ULTRAPARALLEL_CONFLICT_MATCHING_ONE_BASED = (
    (31, 172),
    (32, 169),
    (33, 170),
    (36, 159),
    (39, 161),
    (42, 164),
    (48, 168),
    (51, 167),
    (56, 171),
    (58, 166),
    (59, 165),
    (60, 163),
    (61, 112),
    (62, 162),
    (63, 158),
    (64, 160),
    (65, 150),
    (66, 154),
    (67, 156),
    (68, 155),
    (69, 153),
    (70, 152),
    (71, 147),
    (72, 149),
    (73, 157),
    (74, 151),
    (75, 148),
    (76, 146),
)


def _graph_theorem(
    seed: dict[str, Any],
    stock: dict[str, Any],
    *,
    name: str,
    edge_relation: Callable[[str], bool],
    expected_maximum_one_based: tuple[int, ...],
    expected_histogram: Counter[int],
    expected_degree_histogram: Counter[int],
    conflict_matching_one_based: tuple[tuple[int, int], ...],
) -> dict[str, Any]:
    vertices: tuple[int, ...] = stock["doubling"]
    relations: tuple[tuple[str, ...], ...] = stock["relations"]
    adjacency = {vertex: set() for vertex in vertices}
    for position, left in enumerate(vertices):
        for right in vertices[:position]:
            if edge_relation(relations[left][right]):
                adjacency[left].add(right)
                adjacency[right].add(left)

    cliques = _maximal_cliques(vertices, adjacency)
    size_histogram = Counter(map(len, cliques))
    maximum = max(size_histogram)
    maximum_cliques = tuple(clique for clique in cliques if len(clique) == maximum)
    expected_maximum = tuple(wall - 1 for wall in expected_maximum_one_based)
    excluded = tuple(vertex for vertex in vertices if vertex not in expected_maximum)
    max_containing_excluded = {
        vertex: max(len(clique) for clique in cliques if vertex in clique)
        for vertex in excluded
    }

    matching = tuple((left - 1, right - 1) for left, right in conflict_matching_one_based)
    matching_vertices = tuple(vertex for edge in matching for vertex in edge)
    matching_is_conflicting = all(
        right not in adjacency[left] for left, right in matching
    )
    matching_upper_bound = len(vertices) - len(matching)
    relation_histogram = Counter(
        relations[left][right]
        for position, left in enumerate(expected_maximum)
        for right in expected_maximum[:position]
    )
    clique_payload = [[wall + 1 for wall in clique] for clique in cliques]
    checks = {
        "every_enumerated_set_is_a_clique": all(
            all(right in adjacency[left] for i, left in enumerate(clique) for right in clique[:i])
            for clique in cliques
        ),
        "every_enumerated_clique_is_maximal": all(
            all(
                any(other not in adjacency[wall] for wall in clique)
                for other in vertices
                if other not in clique
            )
            for clique in cliques
        ),
        "maximal_clique_size_histogram_is_exact": size_histogram
        == expected_histogram,
        "degree_histogram_is_exact": Counter(map(len, adjacency.values()))
        == expected_degree_histogram,
        "unique_maximum_is_exact_displayed_set": maximum_cliques
        == (expected_maximum,),
        "excluded_vertex_condition_rules_out_second_maximum": all(
            bound < maximum for bound in max_containing_excluded.values()
        ),
        "conflict_matching_has_distinct_vertices": len(set(matching_vertices))
        == 2 * len(matching),
        "conflict_matching_edges_are_exact_graph_nonedges": matching_is_conflicting,
        "matching_upper_bound_equals_constructed_maximum": matching_upper_bound
        == maximum,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"{name} graph theorem failed: {checks}")
    return {
        "relation": name,
        "vertex_count": len(vertices),
        "edge_count": sum(map(len, adjacency.values())) // 2,
        "degree_histogram": {
            str(key): value
            for key, value in sorted(Counter(map(len, adjacency.values())).items())
        },
        "maximal_clique_count": len(cliques),
        "maximal_clique_size_histogram": {
            str(key): value for key, value in sorted(size_histogram.items())
        },
        "maximum_size": maximum,
        "unique_maximum_one_based": list(expected_maximum_one_based),
        "unique_maximum_relation_histogram": dict(sorted(relation_histogram.items())),
        "excluded_vertices_one_based": [wall + 1 for wall in excluded],
        "maximum_clique_size_containing_each_excluded_vertex": {
            str(wall + 1): bound
            for wall, bound in sorted(max_containing_excluded.items())
        },
        "conflict_matching_one_based": [list(edge) for edge in conflict_matching_one_based],
        "conflict_matching_size": len(matching),
        "matching_upper_bound": matching_upper_bound,
        "all_maximal_cliques_one_based": clique_payload,
        "all_maximal_cliques_sha256": _json_sha256(clique_payload),
        "selected_root_coordinates_sha256": _json_sha256(
            [list(seed["roots"][wall]) for wall in expected_maximum]
        ),
        "checks": checks,
        "all_exact_checks_pass": True,
    }


def _smith_diagonal(columns: Sequence[Vector]) -> list[int]:
    matrix = sp.Matrix.hstack(*(sp.Matrix(column) for column in columns))
    smith = smith_normal_form(matrix)
    return [
        abs(int(smith[index, index]))
        for index in range(min(smith.rows, smith.cols))
        if smith[index, index]
    ]


def _carrier_and_volume(
    seed: dict[str, Any], disjoint: dict[str, Any], ultra: dict[str, Any]
) -> dict[str, Any]:
    vinberg = json.loads(VINBERG.read_text(encoding="utf-8"))
    cone = json.loads(EXACT_CONE.read_text(encoding="utf-8"))
    carrier = json.loads(CARRIER.read_text(encoding="utf-8"))
    roots: tuple[Vector, ...] = seed["roots"]
    disjoint_set = {wall - 1 for wall in disjoint["unique_maximum_one_based"]}
    ultra_set = {wall - 1 for wall in ultra["unique_maximum_one_based"]}
    disjoint_core_indices = tuple(i for i in range(172) if i not in disjoint_set)
    ultra_core_indices = tuple(i for i in range(172) if i not in ultra_set)
    disjoint_smith = _smith_diagonal([roots[i] for i in disjoint_core_indices])
    ultra_smith = _smith_diagonal([roots[i] for i in ultra_core_indices])
    checks = {
        "literal_Vinberg_seed_is_finite_volume_noncompact": (
            vinberg["vinberg_decision"] == "finite_volume_noncompact"
            and vinberg["finite_volume"]
            and not vinberg["compact"]
            and vinberg["ideal_vertex_count"] == 6214
            and vinberg["bad_edge_count"] == 0
            and vinberg["unrecognized_endpoint_edge_count"] == 0
            and vinberg["invalid_full_zero_neighborhood_count"] == 0
        ),
        "independent_exact_cone_has_no_hyperideal_ray": (
            cone["decision"]["classification"] == "finite_volume_noncompact"
            and cone["decision"]["finite_volume"]
            and not cone["decision"]["compact"]
            and cone["exact_double_description"]["hyperideal_vertex_count"] == 0
            and all(cone["checks"].values())
        ),
        "full_Picard_K3_carrier_certificate_is_green": carrier["all_checks_pass"]
        and all(carrier["checks"].values()),
        "all_172_seed_reflections_are_integral_on_full_carrier": carrier["checks"]
        ["all_172_reflections_integral_on_full_carrier"],
        "112_generator_complement_is_exact_W1_through_W60": disjoint_core_indices
        == tuple(range(60)),
        "96_generator_complement_is_exact_W1_through_W76": ultra_core_indices
        == tuple(range(76)),
        "60_wall_boundary_core_spans_full_Picard_lattice": disjoint_smith
        == [1] * 13,
        "76_wall_boundary_core_spans_full_Picard_lattice": ultra_smith
        == [1] * 13,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"carrier/volume gate failed: {checks}")
    carrier_data = carrier["chamber"]["carrier"]
    return {
        "seed_classification": "finite_volume_noncompact",
        "seed_facet_count": 172,
        "seed_finite_vertex_count": vinberg["finite_vertex_count"],
        "seed_cusp_count": vinberg["ideal_vertex_count"],
        "carrier_lattice": carrier_data["name"],
        "carrier_rank": carrier_data["rank"],
        "carrier_signature": carrier_data["signature"],
        "carrier_moduli_dimension": carrier_data["very_general_period_consequence"]
        ["moduli_dimension"],
        "disjoint_112_boundary_core": {
            "walls_one_based": [index + 1 for index in disjoint_core_indices],
            "wall_count": len(disjoint_core_indices),
            "rank": 13,
            "Smith_diagonal": disjoint_smith,
            "index_in_full_Picard_lattice": 1,
        },
        "ultraparallel_96_boundary_core": {
            "walls_one_based": [index + 1 for index in ultra_core_indices],
            "wall_count": len(ultra_core_indices),
            "rank": 13,
            "Smith_diagonal": ultra_smith,
            "index_in_full_Picard_lattice": 1,
        },
        "inheritance_reason": (
            "after translating a finite subtree to contain the identity chamber, "
            "every seed wall outside the selected generator stock remains an "
            "exterior wall; these displayed cores span the Picard lattice with "
            "index one.  Every other boundary root is an image of a primitive "
            "seed root under integral Picard isometries"
        ),
        "scope": (
            "weak full-Picard divisor-reflection carrier only; no boundary root "
            "is asserted effective, no descendant is asserted to be a nef cone, "
            "and no Picard reflection is asserted to be a surface automorphism"
        ),
        "checks": checks,
        "all_exact_checks_pass": True,
    }


def build() -> dict[str, Any]:
    seed = _reconstruct_seed()
    stock = _doubling_stock(seed)
    disjoint = _graph_theorem(
        seed,
        stock,
        name="parallel_or_ultraparallel",
        edge_relation=lambda relation: relation in {"parallel", "ultraparallel"},
        expected_maximum_one_based=tuple(range(61, 173)),
        expected_histogram=Counter({36: 9, 55: 36, 74: 42, 93: 12, 112: 1}),
        expected_degree_histogram=Counter({99: 12, 120: 16, 121: 96}),
        conflict_matching_one_based=DISJOINT_CONFLICT_MATCHING_ONE_BASED,
    )
    ultra = _graph_theorem(
        seed,
        stock,
        name="strictly_ultraparallel",
        edge_relation=lambda relation: relation == "ultraparallel",
        expected_maximum_one_based=tuple(range(77, 173)),
        expected_histogram=Counter({32: 32, 57: 12, 61: 16, 96: 1}),
        expected_degree_histogram=Counter({56: 12, 64: 16, 112: 96}),
        conflict_matching_one_based=ULTRAPARALLEL_CONFLICT_MATCHING_ONE_BASED,
    )
    carrier = _carrier_and_volume(seed, disjoint, ultra)
    checks = {
        "raw_published_seed_reconstruction_green": all(seed["checks"].values()),
        "exact_124_wall_doubling_stock_green": all(stock["checks"].values()),
        "unique_maximal_112_disjoint_stock_green": all(disjoint["checks"].values()),
        "unique_maximal_96_ultraparallel_stock_green": all(ultra["checks"].values()),
        "ultraparallel_maximum_is_subset_of_disjoint_maximum": set(
            ultra["unique_maximum_one_based"]
        )
        < set(disjoint["unique_maximum_one_based"]),
        "finite_volume_and_full_Picard_inheritance_green": all(
            carrier["checks"].values()
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"maximal-disjoint theorem failed: {checks}")

    theorem = {
        "hyperbolic_dimension": 12,
        "seed": "JK-H12-172",
        "unique_maximal_disjoint_doubling_stock_one_based": list(range(61, 173)),
        "maximal_disjoint_doubling_stock_size": 112,
        "disjoint_definition": (
            "the wall hyperplanes do not meet in H12; both parallel and "
            "ultraparallel pairs count, exactly as in Allcock's hypothesis"
        ),
        "unique_maximal_strictly_ultraparallel_stock_one_based": list(
            range(77, 173)
        ),
        "maximal_strictly_ultraparallel_stock_size": 96,
        "finite_subtree_statement": (
            "for k=112 (or k=96), let Gamma be the k-regular Cayley tree of "
            "the free product of k copies of C2 generated by the selected wall "
            "reflections.  Every finite subtree T gives a finite-sided, finite-"
            "volume, noncompact H12 Coxeter polytope Q_T"
        ),
        "exact_index_and_volume": (
            "[W(Q):W(Q_T)]=|V(T)| and Vol(Q_T)=|V(T)|*Vol(Q)"
        ),
        "allcock_4_1_applicability": {
            "applies": True,
            "hypotheses": (
                "Q is a Coxeter polytope and the selected walls are pairwise "
                "disjoint doubling walls"
            ),
            "conclusion": "every subtree union is a Coxeter polyhedron",
        },
        "allcock_4_3_applicability": {
            "applies": True,
            "hypotheses": (
                "Q has finite volume and at least three pairwise disjoint "
                "doubling walls; in fact 96 are already pairwise ultraparallel"
            ),
            "conclusion": (
                "the number of index-I reflection subgroups up to conjugacy by "
                "hyperbolic isometries, and hence the associated nonisometric "
                "finite-subtree polytopes, has an exponential lower bound in I"
            ),
        },
        "carrier": carrier,
    }
    return {
        "schema": SCHEMA,
        "status_date": STATUS_DATE,
        "all_checks_pass": all(checks.values()),
        "independence": {
            "published_source_reparsed": True,
            "fixed_lattice_and_172_roots_reconstructed": True,
            "other_JK_H12_scout_or_certifier_imported": False,
            "graph_library_used": False,
            "maximal_clique_algorithm": "independent deterministic Bron--Kerbosch",
            "separate_upper_bound": "explicit vertex-disjoint conflict matchings",
        },
        "pinned_inputs": seed["hashes"],
        "seed_reconstruction": {
            "dimension": 12,
            "facet_count": 172,
            "fixed_rank": 13,
            "fixed_signature": [1, 12, 0],
            "fixed_determinant": 256,
            "primitive_root_coordinates_sha256": _json_sha256(
                [list(root) for root in seed["roots"]]
            ),
            "checks": seed["checks"],
        },
        "doubling_stock": {
            "criterion": (
                "every finite incident Coxeter order is even; in this exact "
                "diagram this is equivalent to having no incident m=3 edge"
            ),
            "doubling_wall_count": len(stock["doubling"]),
            "doubling_walls_one_based": [wall + 1 for wall in stock["doubling"]],
            "whole_seed_relation_histogram": stock["whole_histogram"],
            "doubling_stock_relation_histogram": stock["stock_histogram"],
            "checks": stock["checks"],
        },
        "maximal_disjoint_graph_theorem": disjoint,
        "maximal_ultraparallel_graph_theorem": ultra,
        "theorem": theorem,
        "lineage_and_scope": {
            "allcock_mechanism": (
                "classical: Theorems 4.1 and 4.3 of Allcock, Infinitely many "
                "hyperbolic Coxeter groups through dimension 19, Geom. Topol. "
                "10 (2006), 737-758, DOI 10.2140/gt.2006.10.737"
            ),
            "historical_H12_boundary": (
                "Allcock already proved infinitely and exponentially many "
                "finite-covolume H12 Coxeter examples using a Conway-face seed"
            ),
            "project_gain": (
                "the exact unique 112-wall disjoint and 96-wall ultraparallel "
                "maxima inside this surface-derived 124-wall doubling stock, "
                "together with full-Picard inheritance for every finite subtree"
            ),
            "maximum_scope": (
                "maximality and uniqueness are only inside the exact 124 "
                "doubling walls of JK-H12-172; descendant walls, other seeds, "
                "and arbitrary reflection hyperplanes are not classified"
            ),
            "worldwide_priority": "open",
            "worldwide_novelty_claimed": False,
            "new_dimension_claimed": False,
            "new_arithmetic_commensurability_class_claimed": False,
            "full_nef_claimed": False,
            "surface_automorphism_claimed": False,
        },
        "checks": checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    payload = build()
    if arguments.check_archive:
        archived = json.loads(arguments.output.read_text(encoding="utf-8"))
        if payload != archived:
            raise SystemExit("archive differs from exact rebuild")
        print(f"archive replay passed: {arguments.output}")
        return
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(arguments.output)
    print(f"sha256={_sha256(arguments.output)}")


if __name__ == "__main__":
    main()
