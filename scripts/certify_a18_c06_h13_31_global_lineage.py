#!/usr/bin/env python3
"""Exact denominator and named-lineage audit for ``A18_4_0_C06``.

The target is the 31-facet finite-volume chamber in H^13 extracted from
Alexeev's (18,4,0) wall stock.  This certificate does four things that the
earlier bounded header comparison did not do:

* freezes its literal face vector, Coxeter graph and integral mirror data;
* identifies the 15-wall maximal denominator with Tumarkin's H^13 pyramid
  by an explicit integral lattice isometry and an independent Vinberg audit;
* exhausts pure doubling descendants through 31 facets of the named H^13
  Vinberg--McLeod--Allcock seeds encoded below; and
* exhausts every subtree-gluing construction coming from the pairwise
  disjoint doubling walls of Allcock's published H^13 chambers F and Q.

The last two bullets are named-family exclusions, not a worldwide priority
theorem.  In particular, arbitrary finite-index reflection subgroups of the
15-generator denominator and unrelated constructions remain outside scope.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict, deque
import hashlib
from itertools import combinations
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Any, Iterable

import networkx as nx
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.certify_borcherds_char2_high_fixed_sections_vinberg import (
    _parse_output,
)


ROOT = Path(__file__).resolve().parents[1]
CONSTRUCTION = (
    ROOT
    / "data/enumerations/alexeev_remaining_high_fixed_section_construction.json"
)
VINBERG = ROOT / "data/enumerations/alexeev_remaining_high_fixed_section_vinberg.json"
STABLE = (
    ROOT
    / "data/enumerations/alexeev_remaining_high_fixed_section_stable_overlattices.json"
)
EARLIER_LINEAGE = (
    ROOT / "data/enumerations/alexeev_remaining_high_fixed_section_lineage.json"
)
CPP_SOURCE = ROOT / "scripts/certify_borcherds_char2_high_fixed_sections_vinberg.cpp"
OUTPUT = ROOT / "data/enumerations/a18_c06_h13_31_global_lineage.json"
SCHEMA = "coxeter.a18_c06_h13_31.global_lineage.v1"

SQRT2 = sp.sqrt(2)
SQRT3 = sp.sqrt(3)
EDGE_MATCH = nx.algorithms.isomorphism.categorical_edge_match("label", None)


# The integral isometry found by splitting a hyperbolic plane in both
# lattices and using an exact isometry of their negative-definite complements.
# It is frozen here so replay needs neither PARI/GP nor a heuristic search.
# With Gs the maximal stable lattice and Gt the Tumarkin quotient below,
# M.T * Gs * M == Gt.
STABLE_TO_TUMARKIN = sp.Matrix(
    [
        [-12, -1, 6, -1, 2, 4, 7, -55, 5, 2, 0, 8, -36, 0],
        [22, 2, -10, 0, -3, -7, -13, 100, -9, -4, 1, -16, 66, 0],
        [46, 4, -21, 0, -6, -15, -27, 210, -20, -7, 2, -34, 138, 0],
        [10, 0, -4, 0, -1, -4, -5, 45, -4, -2, 1, -8, 30, 0],
        [26, 1, -11, 0, -3, -9, -14, 115, -11, -4, 1, -18, 75, 0],
        [4, 0, -2, 0, -1, -1, -3, 20, -2, -1, 0, -2, 12, 1],
        [-19, -1, 8, 0, 2, 6, 11, -85, 8, 3, -1, 14, -56, 1],
        [-35, -3, 16, 0, 4, 12, 21, -160, 15, 6, -1, 24, -104, 0],
        [5, 0, -2, 0, -1, -1, -2, 20, -2, -1, 1, -4, 13, 0],
        [-23, -2, 11, -1, 3, 8, 14, -105, 10, 4, -1, 16, -68, 0],
        [-15, -1, 7, 0, 1, 5, 8, -65, 6, 2, -1, 12, -44, 0],
        [18, 1, -8, 0, -2, -6, -10, 80, -7, -3, 1, -13, 52, 0],
        [-6, -1, 3, 0, 1, 2, 4, -29, 3, 1, 0, 4, -19, 0],
        [-42, -3, 19, 0, 5, 14, 24, -191, 18, 7, -2, 30, -124, 0],
    ]
)


# CoxIter commit b800be48d8240fbabd2d1c9542e61054584ad7fe is used only as
# a provenance marker.  The exact edge lists are embedded below, so replay is
# independent of a temporary checkout.
COXITER_TRANSCRIPTS = {
    "Tumarkin_13_01": {
        "sha256": "d1de070e5aa7cc90292f21424625ee8ef7dfbe718d004fb8434cc8beee0b35c8",
        "facets": 15,
        "edges": (
            [(i, i + 1, 3) for i in range(11)]
            + [(11, 12, 4), (12, 13, 3), (2, 14, 3)]
        ),
    },
    "Tumarkin_13_02": {
        "sha256": "871d05dd8911ee8f2ef3ca305753e3dcc1e7c7d1489f9a2fa005f45e6ce5b8d1",
        "facets": 15,
        "edges": (
            [(i, i + 1, 3) for i in range(11)]
            + [(11, 12, 4), (10, 13, 3), (2, 14, 3)]
        ),
    },
    "Tumarkin_13_03": {
        "sha256": "7d11d1dab27f3c8c7c1ff51dc32df2d70a98962049587ffa3e4f9ce562274c39",
        "facets": 15,
        "edges": (
            [(i, i + 1, 3) for i in range(11)]
            + [(10, 12, 3), (10, 13, 3), (2, 14, 3)]
        ),
    },
    "Vinberg_f1_P13": {
        "sha256": "20bb824f5594699ac0bdbcf561785a2f06bd59e3d5aadf2619e28e4b6d935643",
        "facets": 15,
        "edges": (
            [(i, i + 1, 3) for i in range(11)]
            + [(11, 12, 4), (2, 13, 3), (10, 14, 3)]
        ),
    },
    "Vinberg_f2_H13": {
        "sha256": "70e12ac824ad23f00ac0e12b7d3cdc83cb0e551e74493571e95a9f270ee36e5f",
        "facets": 19,
        "edges": [
            *[(i, i + 1, 3) for i in range(11)],
            (11, 12, 4),
            (3, 13, 3),
            (0, 14, 0),
            (9, 15, 3),
            (14, 15, 0),
            (0, 16, 0),
            (11, 16, 3),
            (1, 17, 4),
            (12, 17, 0),
            (14, 17, 0),
            (12, 18, 0),
            (13, 18, 0),
            (14, 18, 0),
            (16, 18, 0),
        ],
    },
    "McLeod_f3_H13": {
        "sha256": "20cadc50f35becfe2ed25fa54b37d96bd8ca1610ed31a4acef0d995f45861c5b",
        "facets": 22,
        "edges": [
            *[(i, i + 1, 3) for i in range(11)],
            (11, 12, 4),
            (0, 13, 6),
            (4, 14, 3),
            (0, 15, 3),
            (10, 15, 3),
            (12, 16, 0),
            (13, 16, 0),
            (14, 16, 4),
            (8, 17, 6),
            (13, 17, 0),
            (16, 17, 0),
            (1, 18, 0),
            (12, 18, 4),
            (16, 18, 4),
            (17, 18, 0),
            (2, 19, 6),
            (12, 19, 0),
            (13, 19, 0),
            (17, 19, 0),
            (7, 20, 0),
            (12, 20, 4),
            (13, 20, 0),
            (16, 20, 4),
            (18, 20, 0),
            (19, 20, 0),
            (6, 21, 6),
            (12, 21, 0),
            (13, 21, 0),
            (17, 21, 0),
            (18, 21, 0),
            (19, 21, 0),
        ],
    },
}


ALLCOCK_NODES = tuple("ABCDEFGHIJKLMNOPQRSTU")
ALLCOCK_BONDS_BY_NAME = {
    # Double bonds.
    ("A", "B"): 4,
    ("B", "C"): 4,
    ("B", "D"): 4,
    ("B", "E"): 4,
    ("H", "J"): 4,
    ("I", "J"): 4,
    ("M", "J"): 4,
    ("N", "J"): 4,
    # Heavy/dashed bonds: Coxeter exponent infinity.
    ("D", "E"): 0,
    ("D", "F"): 0,
    ("E", "G"): 0,
    ("H", "F"): 0,
    ("I", "G"): 0,
    ("D", "I"): 0,
    ("E", "H"): 0,
    ("H", "I"): 0,
    ("H", "A"): 0,
    ("I", "C"): 0,
    ("M", "D"): 0,
    ("N", "E"): 0,
    # Single bonds.
    ("A", "K"): 3,
    ("M", "K"): 3,
    ("C", "L"): 3,
    ("N", "L"): 3,
    ("G", "P"): 3,
    ("F", "T"): 3,
    ("K", "O"): 3,
    ("O", "P"): 3,
    ("P", "Q"): 3,
    ("Q", "R"): 3,
    ("R", "S"): 3,
    ("S", "T"): 3,
    ("T", "U"): 3,
    ("U", "L"): 3,
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _payload_sha256(payload: Any) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def _matrix(payload: list[list[str | int]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(value) for value in row] for row in payload])


def _case(payload: dict[str, Any], identifier: str = "A18_4_0_C06") -> dict[str, Any]:
    return next(row for row in payload["cases"] if row["id"] == identifier)


def _relation(value: sp.Expr) -> int:
    value = sp.expand(value)
    if value == 0:
        return 2
    if value == -sp.Rational(1, 2):
        return 3
    if value == -SQRT2 / 2:
        return 4
    if value == -SQRT3 / 2:
        return 6
    if bool(value <= -1):
        return 0
    raise ArithmeticError(f"unexpected normalized root product: {value}")


def _seed_gram(facets: int, edges: Iterable[tuple[int, int, int]]) -> sp.Matrix:
    gram = sp.eye(facets)
    values = {
        0: -sp.Integer(1),
        3: -sp.Rational(1, 2),
        4: -SQRT2 / 2,
        6: -SQRT3 / 2,
    }
    for left, right, exponent in edges:
        gram[left, right] = gram[right, left] = values[exponent]
    return gram


def _allcock_gram() -> sp.Matrix:
    indices = {name: position for position, name in enumerate(ALLCOCK_NODES)}
    return _seed_gram(
        len(ALLCOCK_NODES),
        [
            (indices[left], indices[right], exponent)
            for (left, right), exponent in ALLCOCK_BONDS_BY_NAME.items()
        ],
    )


def _unit_roots(size: int) -> tuple[sp.Matrix, ...]:
    identity = sp.eye(size)
    return tuple(identity.col(index) for index in range(size))


def _pairing(left: sp.Matrix, ambient: sp.Matrix, right: sp.Matrix) -> sp.Expr:
    return sp.expand((left.T * ambient * right)[0])


def _relation_matrix(
    roots: tuple[sp.Matrix, ...], ambient: sp.Matrix
) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(
            1 if left == right else _relation(_pairing(alpha, ambient, beta))
            for right, beta in enumerate(roots)
        )
        for left, alpha in enumerate(roots)
    )


def _relation_graph(relations: tuple[tuple[int, ...], ...]) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(range(len(relations)))
    for left in range(len(relations)):
        for right in range(left + 1, len(relations)):
            if relations[left][right] != 2:
                graph.add_edge(left, right, label=str(relations[left][right]))
    return graph


def _relation_histogram(relations: tuple[tuple[int, ...], ...]) -> dict[str, int]:
    return {
        str(key): value
        for key, value in sorted(
            Counter(
                relations[left][right]
                for left in range(len(relations))
                for right in range(left + 1, len(relations))
            ).items()
        )
    }


def _graph_hash(relations: tuple[tuple[int, ...], ...]) -> str:
    return nx.weisfeiler_lehman_graph_hash(
        _relation_graph(relations), edge_attr="label"
    )


def _doubling_walls(relations: tuple[tuple[int, ...], ...]) -> list[int]:
    return [
        index
        for index, row in enumerate(relations)
        if all(value in (0, 1, 2) or value % 2 == 0 for value in row)
    ]


def _double_across(
    roots: tuple[sp.Matrix, ...], ambient: sp.Matrix, wall: int
) -> tuple[sp.Matrix, ...]:
    alpha = roots[wall]
    output: list[sp.Matrix] = []
    for index, beta in enumerate(roots):
        if index == wall:
            continue
        product = _pairing(beta, ambient, alpha)
        output.append(beta)
        if product != 0:
            reflected = (beta - 2 * product * alpha).applyfunc(sp.expand)
            if all(reflected != prior for prior in output):
                output.append(reflected)
    return tuple(output)


def _target_relations(case: dict[str, Any]) -> tuple[tuple[int, ...], ...]:
    gram = _matrix(case["restricted_mirrors"]["facet_gram"])
    relations = []
    for left in range(gram.rows):
        row = []
        for right in range(gram.cols):
            if left == right:
                row.append(1)
                continue
            ratio = sp.factor(
                gram[left, right] ** 2 / (gram[left, left] * gram[right, right])
            )
            if ratio == 0:
                row.append(2)
            elif ratio == sp.Rational(1, 4):
                row.append(3)
            elif ratio == sp.Rational(1, 2):
                row.append(4)
            elif ratio >= 1:
                row.append(0)
            else:
                raise ArithmeticError((left, right, ratio))
        relations.append(tuple(row))
    return tuple(relations)


def _same_relation_graph(
    left: tuple[tuple[int, ...], ...], right: tuple[tuple[int, ...], ...]
) -> bool:
    return nx.is_isomorphic(
        _relation_graph(left), _relation_graph(right), edge_match=EDGE_MATCH
    )


def _pure_doubling_audit(
    name: str,
    ambient: sp.Matrix,
    roots: tuple[sp.Matrix, ...],
    target: tuple[tuple[int, ...], ...],
) -> dict[str, Any]:
    queue: deque[tuple[tuple[sp.Matrix, ...], tuple[int, ...]]] = deque([(roots, ())])
    seen: dict[int, list[tuple[tuple[int, ...], ...]]] = defaultdict(list)
    target_paths: list[list[int]] = []
    while queue:
        chamber_roots, path = queue.popleft()
        relations = _relation_matrix(chamber_roots, ambient)
        size = len(chamber_roots)
        if any(_same_relation_graph(relations, prior) for prior in seen[size]):
            continue
        seen[size].append(relations)
        if size == len(target) and _same_relation_graph(relations, target):
            target_paths.append(list(path))
        for wall in _doubling_walls(relations):
            child = _double_across(chamber_roots, ambient, wall)
            if len(child) <= len(target):
                queue.append((child, path + (wall,)))
    size_fingerprints = [
        {
            "relation_histogram": _relation_histogram(relations),
            "wl_hash": _graph_hash(relations),
            "target_isomorphic": _same_relation_graph(relations, target),
        }
        for relations in seen.get(len(target), [])
    ]
    return {
        "seed": name,
        "seed_facet_count": len(roots),
        "seed_doubling_walls_zero_based": _doubling_walls(
            _relation_matrix(roots, ambient)
        ),
        "isomorphism_classes_by_facet_count": {
            str(size): len(rows) for size, rows in sorted(seen.items())
        },
        "explored_isomorphism_class_count": sum(map(len, seen.values())),
        "facet_31_fingerprints": size_fingerprints,
        "target_match_paths": target_paths,
        "target_absent": not target_paths,
    }


def _target_invariants(
    construction_case: dict[str, Any],
    vinberg_case: dict[str, Any],
    stable_case: dict[str, Any],
) -> dict[str, Any]:
    relations = _target_relations(construction_case)
    graph = _relation_graph(relations)
    automorphisms = list(
        nx.algorithms.isomorphism.GraphMatcher(
            graph, graph, edge_match=EDGE_MATCH
        ).isomorphisms_iter()
    )
    remaining = set(graph.nodes)
    orbits = []
    while remaining:
        seed = min(remaining)
        orbit = sorted({mapping[seed] for mapping in automorphisms})
        orbits.append(orbit)
        remaining.difference_update(orbit)

    special = vinberg_case["finite_special_subset_counts_by_rank"]
    face_vector = [
        special[-1] + vinberg_case["ideal_vertex_count"],
        *reversed(special[:-1]),
        1,
    ]
    if face_vector != [
        3671,
        28766,
        96280,
        194317,
        269179,
        270518,
        201943,
        112606,
        46543,
        13964,
        2925,
        399,
        31,
        1,
    ]:
        raise ArithmeticError("the target face vector changed")

    maximal = stable_case["unique_maximal_reflection_stable_overlattice"]
    maximal_gram = _matrix(maximal["gram"])
    basis = _matrix(maximal["basis_columns_in_original_fixed_lattice"])
    primitive_rows = []
    for record in construction_case["restricted_mirrors"]["records"]:
        if record["hyperbolic_status"] != "active_facet":
            continue
        old = sp.Matrix(record["primitive_integral_same_direction_coordinates"])
        coordinates = basis.inv() * old
        if not all(value.q == 1 for value in coordinates):
            raise ArithmeticError("a target root left the stable overclosure")
        content = math.gcd(*(abs(int(value)) for value in coordinates))
        primitive = coordinates / content
        norm = int((primitive.T * maximal_gram * primitive)[0])
        divisibility = math.gcd(
            *(abs(int(value)) for value in maximal_gram * primitive)
        )
        primitive_rows.append((norm, divisibility))
    primitive_histogram = {
        f"norm_{norm}_div_{divisibility}": count
        for (norm, divisibility), count in sorted(Counter(primitive_rows).items())
    }
    if primitive_histogram != {"norm_-4_div_2": 5, "norm_-2_div_1": 26}:
        raise ArithmeticError("primitive maximal-lattice mirror types changed")

    positive_gram = construction_case["restricted_mirrors"]["positive_wall_gram"]
    checks = {
        "construction_archive_passes": construction_case["construction_checks"][
            "coxeter_pairs"
        ],
        "literal_vinberg_finite_volume": (
            vinberg_case["classification"]
            == "finite_volume_noncompact_Coxeter_polytope"
            and vinberg_case["bad_edge_count"] == 0
            and vinberg_case["bad_parabolic_neighborhood_count"] == 0
        ),
        "target_relation_histogram": _relation_histogram(relations)
        == {"0": 66, "2": 327, "3": 44, "4": 28},
        "automorphism_group_order_48": len(automorphisms) == 48,
        "automorphism_orbits_frozen": orbits
        == [
            [0, 26],
            [1, 9, 18, 19, 22, 23, 24, 25, 27, 28, 29, 30],
            [2, 8, 10, 15, 16, 17, 20, 21],
            [3, 7, 11, 14],
            [4, 6, 12, 13],
            [5],
        ],
        "unique_index_four_stable_overclosure": (
            maximal["unique"]
            and maximal["index_over_original_lattice"] == 4
            and maximal["determinant"] == -4
        ),
        "primitive_mirror_types_frozen": primitive_histogram
        == {"norm_-4_div_2": 5, "norm_-2_div_1": 26},
    }
    if not all(checks.values()):
        raise ArithmeticError(f"target invariant audit changed: {checks}")
    return {
        "id": "A18-C2c-H13-31",
        "internal_id": "A18_4_0_C06",
        "hyperbolic_dimension": 13,
        "facet_count": 31,
        "volume_type": "finite-volume noncompact",
        "finite_vertex_count": vinberg_case["finite_vertex_count"],
        "ideal_vertex_count": vinberg_case["ideal_vertex_count"],
        "elliptic_edge_count": vinberg_case["elliptic_edge_count"],
        "edge_endpoint_type_histogram": vinberg_case[
            "edge_endpoint_type_histogram"
        ],
        "finite_special_subset_counts_by_rank": special,
        "face_vector_f0_through_f13": face_vector,
        "coxeter_relation_histogram": _relation_histogram(relations),
        "coxeter_relation_wl_hash": _graph_hash(relations),
        "positive_wall_gram_sha256": _payload_sha256(positive_gram),
        "coxeter_graph_automorphism_group_order": len(automorphisms),
        "coxeter_graph_automorphism_orbits_zero_based": orbits,
        "fixed_lattice": {
            "rank": 14,
            "determinant": construction_case["saturated_fixed_lattice"][
                "determinant"
            ],
            "smith_invariants": construction_case["saturated_fixed_lattice"][
                "smith_diagonal_absolute"
            ],
            "inertia": construction_case["saturated_fixed_lattice"]["inertia"],
        },
        "maximal_reflection_stable_overclosure": {
            "unique": maximal["unique"],
            "index_over_fixed_lattice": maximal["index_over_original_lattice"],
            "determinant": maximal["determinant"],
            "smith_invariants": maximal["smith_invariants"],
            "parity": maximal["parity"],
            "primitive_target_mirror_type_histogram": primitive_histogram,
        },
        "checks": checks,
    }


def _tumarkin_integer_gram() -> sp.Matrix:
    transcript = COXITER_TRANSCRIPTS["Tumarkin_13_01"]
    norms = [2] * 15
    norms[12] = norms[13] = 4
    gram = sp.zeros(15)
    for index, norm in enumerate(norms):
        gram[index, index] = norm
    for left, right, exponent in transcript["edges"]:
        if exponent == 3:
            product = -sp.sqrt(norms[left] * norms[right]) / 2
        elif exponent == 4:
            product = -sp.sqrt(2 * norms[left] * norms[right]) / 2
        else:
            raise ArithmeticError("unexpected Tumarkin bond")
        if not product.is_Integer:
            raise ArithmeticError("Tumarkin symmetrization ceased to be integral")
        gram[left, right] = gram[right, left] = product
    return gram


def _run_tumarkin_vinberg(gram: sp.Matrix, compiler: str | None = None) -> dict[str, Any]:
    selected = compiler or shutil.which("clang++") or shutil.which("g++")
    if selected is None:
        raise RuntimeError("a C++17 compiler is required for the Tumarkin audit")
    with tempfile.TemporaryDirectory(prefix="a18-c06-tumarkin-") as raw:
        executable = Path(raw) / "vinberg"
        subprocess.run(
            [
                selected,
                "-O3",
                "-std=c++17",
                "-Wall",
                "-Wextra",
                "-Werror",
                str(CPP_SOURCE),
                "-o",
                str(executable),
            ],
            check=True,
        )
        stream = "13 15\n" + "\n".join(
            " ".join(str(int(value)) for value in gram.row(index))
            for index in range(gram.rows)
        ) + "\n"
        completed = subprocess.run(
            [str(executable)],
            input=stream,
            text=True,
            capture_output=True,
            check=True,
        )
    return _parse_output(completed.stdout)


def _tumarkin_denominator(
    construction_case: dict[str, Any],
    stable_case: dict[str, Any],
    compiler: str | None = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    root_gram = _tumarkin_integer_gram()
    nullspace = root_gram.nullspace()
    if len(nullspace) != 1:
        raise ArithmeticError("Tumarkin diagram did not have corank one")
    denominator = sp.ilcm(*[value.q for value in nullspace[0]])
    radical = [int(value * denominator) for value in nullspace[0]]
    content = math.gcd(*(abs(value) for value in radical))
    radical = [value // content for value in radical]
    if radical[7] < 0:
        radical = [-value for value in radical]
    expected_radical = [2, 4, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -2, -1, 3]
    keep = [index for index in range(15) if index != 7]
    quotient_gram = -root_gram.extract(keep, keep)
    maximal = stable_case["unique_maximal_reflection_stable_overlattice"]
    stable_gram = _matrix(maximal["gram"])
    isometry = STABLE_TO_TUMARKIN
    smith = smith_normal_form(root_gram, domain=sp.ZZ)
    smith_nonzero = [
        abs(int(smith[index, index]))
        for index in range(min(smith.shape))
        if smith[index, index]
    ]
    volume = _run_tumarkin_vinberg(root_gram, compiler)
    expected_special = [
        15,
        105,
        455,
        1365,
        3002,
        4995,
        6390,
        6315,
        4794,
        2745,
        1140,
        315,
        44,
    ]
    maximal_lineage = stable_case["maximal_reflection_lineage"]
    checks = {
        "primitive_radical_frozen": radical == expected_radical,
        "radical_coefficient_one_allows_integral_quotient_basis": radical[7] == 1,
        "root_gram_smith_1^12_2^2_0": smith_nonzero == [1] * 12 + [2, 2],
        "quotient_determinant_minus_four": quotient_gram.det() == -4,
        "isometry_unimodular": abs(isometry.det()) == 1,
        "explicit_integral_isometry": isometry.T * stable_gram * isometry
        == quotient_gram,
        "vinberg_finite_special_census": volume[
            "finite_special_subset_counts_by_rank"
        ]
        == expected_special,
        "vinberg_two_cusps_no_bad_edges_or_neighborhoods": (
            volume["ideal_vertex_count"] == 2
            and volume["bad_edge_count"] == 0
            and volume["bad_parabolic_neighborhood_count"] == 0
        ),
        "matches_SK_List13_entry4_census": (
            maximal_lineage["entry_one_based"] == 4
            and maximal_lineage["genus_symbol"] == "II_{13,1}(2^{-2})"
            and maximal_lineage["full_chamber_facets"] == 15
            and maximal_lineage["full_chamber_finite_vertices"] == 44
            and maximal_lineage["full_chamber_cusps"] == 2
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"Tumarkin denominator identification changed: {checks}")
    record = {
        "primary_diagram": "Tumarkin 2004 Table 13, 13-tum04_13_01",
        "coxiter_transcript_sha256": COXITER_TRANSCRIPTS["Tumarkin_13_01"][
            "sha256"
        ],
        "root_norms": [2] * 12 + [4, 4, 2],
        "integral_symmetrized_root_gram": [
            [int(value) for value in root_gram.row(index)]
            for index in range(root_gram.rows)
        ],
        "primitive_radical_one_based_node_coefficients": radical,
        "quotient_basis": "delete one-based node 8, whose radical coefficient is 1",
        "quotient_gram": [
            [int(value) for value in quotient_gram.row(index)]
            for index in range(quotient_gram.rows)
        ],
        "stable_to_tumarkin_basis_isometry_M": [
            [int(value) for value in isometry.row(index)]
            for index in range(isometry.rows)
        ],
        "isometry_identity": "M^T G_stable M = G_Tumarkin_quotient",
        "vinberg_audit": {
            "finite_special_subset_counts_by_rank": volume[
                "finite_special_subset_counts_by_rank"
            ],
            "finite_vertex_count": volume["finite_special_subset_counts_by_rank"][-1],
            "ideal_vertex_count": volume["ideal_vertex_count"],
            "edge_endpoint_type_histogram": volume[
                "edge_endpoint_type_histogram"
            ],
            "bad_edge_count": volume["bad_edge_count"],
            "bad_parabolic_neighborhood_count": volume[
                "bad_parabolic_neighborhood_count"
            ],
        },
        "scharlau_kirschmer_denominator": {
            "list": "List13",
            "entry_one_based": 4,
            "genus_symbol": "II_{13,1}(2^{-2})",
            "facets": 15,
            "finite_vertices": 44,
            "cusps": 2,
        },
        "checks": checks,
    }
    auxiliary = {
        "root_gram": root_gram,
        "quotient_gram": quotient_gram,
        "keep": keep,
        "radical": radical,
        "stable_gram": stable_gram,
        "stable_basis": _matrix(maximal["basis_columns_in_original_fixed_lattice"]),
    }
    return record, auxiliary


def _reflect(vector: sp.Matrix, root: sp.Matrix, gram: sp.Matrix) -> sp.Matrix:
    return sp.expand(
        vector
        - 2 * (vector.T * gram * root)[0] / (root.T * gram * root)[0] * root
    )


def _same_primitive_line(left: sp.Matrix, right: sp.Matrix) -> bool:
    return left == right or left == -right


def _tumarkin_simple_roots(auxiliary: dict[str, Any]) -> list[sp.Matrix]:
    keep = auxiliary["keep"]
    radical = auxiliary["radical"]
    quotient = auxiliary["quotient_gram"]
    roots = []
    for node in range(15):
        if node == 7:
            root = sp.Matrix([-radical[old] for old in keep])
        else:
            root = sp.zeros(14, 1)
            root[keep.index(node)] = 1
        roots.append(root)
    matrix = sp.Matrix.hstack(*roots)
    if matrix.T * quotient * matrix != -auxiliary["root_gram"]:
        raise ArithmeticError("Tumarkin quotient simple roots changed")
    return roots


def _subgroup_presentation_handoff(
    construction_case: dict[str, Any], auxiliary: dict[str, Any]
) -> dict[str, Any]:
    gram = auxiliary["quotient_gram"]
    simple = _tumarkin_simple_roots(auxiliary)
    stable_basis = auxiliary["stable_basis"]
    isometry = STABLE_TO_TUMARKIN

    interior_original = sp.Matrix(
        [
            sp.Rational(value)
            for value in construction_case["strict_fixed_interior"][
                "fixed_lattice_coordinates"
            ]
        ]
    )
    interior = isometry.inv() * stable_basis.inv() * interior_original
    chamber_word: list[int] = []
    for _ in range(100000):
        pairings = [(interior.T * gram * root)[0] for root in simple]
        bad = [index for index, pairing in enumerate(pairings) if pairing > 0]
        if not bad:
            break
        index = bad[0]
        interior = _reflect(interior, simple[index], gram)
        chamber_word.append(index)
    else:
        raise ArithmeticError("target interior did not reduce to the Tumarkin chamber")
    if len(chamber_word) != 247:
        raise ArithmeticError("target-to-denominator chamber word changed")

    target_roots = []
    for record in construction_case["restricted_mirrors"]["records"]:
        if record["hyperbolic_status"] != "active_facet":
            continue
        old = sp.Matrix(record["primitive_integral_same_direction_coordinates"])
        stable = stable_basis.inv() * old
        content = math.gcd(*(abs(int(value)) for value in stable))
        root = isometry.inv() * (stable / content)
        for index in chamber_word:
            root = _reflect(root, simple[index], gram)
        if not all(value.q == 1 for value in root):
            raise ArithmeticError("a conjugated target root became nonintegral")
        target_roots.append(root)

    # A strict interior point for the denominator chamber.  The chosen
    # negative pairing column is derived from the signs of the affine radical.
    column = [
        -9 if auxiliary["radical"][node] > 0 else -30
        if auxiliary["radical"][node] < 0
        else -1
        for node in auxiliary["keep"]
    ]
    denominator_interior = gram.inv() * sp.Matrix(column)
    if not (
        (denominator_interior.T * gram * denominator_interior)[0] > 0
        and all(
            (denominator_interior.T * gram * root)[0] < 0 for root in simple
        )
    ):
        raise ArithmeticError("the denominator interior witness changed")

    records = []
    for target_index, original_root in enumerate(target_roots):
        root = original_root
        word: list[int] = []
        for _ in range(10000):
            simple_index = next(
                (
                    index
                    for index, candidate in enumerate(simple)
                    if _same_primitive_line(root, candidate)
                ),
                None,
            )
            if simple_index is not None:
                break
            height = abs(int((denominator_interior.T * gram * root)[0]))
            candidates = []
            for index, candidate in enumerate(simple):
                reflected = _reflect(root, candidate, gram)
                reflected_height = abs(
                    int((denominator_interior.T * gram * reflected)[0])
                )
                if reflected_height < height:
                    candidates.append((reflected_height, index, reflected))
            if not candidates:
                raise ArithmeticError(f"root reduction stuck at target {target_index}")
            _, index, root = min(candidates, key=lambda row: (row[0], row[1]))
            word.append(index)
        else:
            raise ArithmeticError("target root reduction exceeded its step bound")

        reconstructed = simple[simple_index]
        for index in reversed(word):
            reconstructed = _reflect(reconstructed, simple[index], gram)
        if not _same_primitive_line(reconstructed, original_root):
            raise ArithmeticError("a target reflection word did not reconstruct")
        records.append(
            {
                "target_facet_one_based": target_index + 1,
                "simple_reflection_one_based": simple_index + 1,
                "conjugating_word_one_based": [index + 1 for index in word],
                "word_length": len(word),
            }
        )

    length_histogram = {
        str(length): count
        for length, count in sorted(
            Counter(record["word_length"] for record in records).items()
        )
    }
    expected_histogram = {
        "0": 7,
        "1": 1,
        "2": 6,
        "6": 1,
        "7": 1,
        "14": 3,
        "28": 1,
        "29": 1,
        "34": 1,
        "35": 2,
        "36": 1,
        "43": 1,
        "44": 2,
        "45": 1,
        "52": 1,
        "53": 1,
    }
    if length_histogram != expected_histogram:
        raise ArithmeticError("target reflection-word histogram changed")

    presentation = {
        "generator_count": 15,
        "coxeter_edges_zero_based": COXITER_TRANSCRIPTS["Tumarkin_13_01"][
            "edges"
        ],
        "default_unlisted_exponent": 2,
        "target_reflection_conjugates": records,
    }
    return {
        "ambient_presentation": (
            "15 involutory generators with the embedded Tumarkin Coxeter exponents"
        ),
        "target_chamber_conjugation_word_length": len(chamber_word),
        "target_reflection_conjugates": records,
        "conjugating_word_length_histogram": length_histogram,
        "presentation_sha256": _payload_sha256(presentation),
        "exact_word_checks_pass": True,
        "coset_enumeration_experiment": {
            "software": "GAP 4.15.1 with ACE",
            "strategy": "hard",
            "workspace_words": 200000000,
            "maximum_temporary_cosets": 13333332,
            "active_cosets_at_overflow": 12951331,
            "result": "OVERFLOW; no complete coset table",
            "logical_inference": (
                "none: temporary cosets can later collapse, so this is neither "
                "an index lower bound nor evidence of infinite index"
            ),
        },
        "index_status": (
            "finite by chamber containment; exact numerical value unresolved"
        ),
    }


def _maximal_disjoint_doubling_sets(
    relations: tuple[tuple[int, ...], ...]
) -> list[list[int]]:
    walls = _doubling_walls(relations)
    graph = nx.Graph()
    graph.add_nodes_from(walls)
    graph.add_edges_from(
        (left, right)
        for left, right in combinations(walls, 2)
        if relations[left][right] == 0
    )
    return sorted((sorted(clique) for clique in nx.find_cliques(graph)))


def _gluing_cost(relations: tuple[tuple[int, ...], ...], wall: int) -> int:
    # Two copies of the gluing wall disappear; every other facet orthogonal
    # to it is coplanar with its reflected copy and merges once.
    return 2 + sum(
        relations[wall][other] == 2
        for other in range(len(relations))
        if other != wall
    )


def _transform_word(
    word: tuple[int, ...],
    vector: sp.Matrix,
    selected_roots: tuple[sp.Matrix, ...],
    ambient: sp.Matrix,
) -> sp.Matrix:
    for letter in reversed(word):
        vector = _reflect(vector, selected_roots[letter], ambient)
    return vector


def _oriented_root_key(vector: sp.Matrix) -> tuple[sp.Expr, ...]:
    values = tuple(sp.expand(value) for value in vector)
    for value in values:
        if value != 0:
            return values if bool(value > 0) else tuple(-entry for entry in values)
    raise ArithmeticError("zero root")


def _tree_boundary(
    chamber_roots: tuple[sp.Matrix, ...],
    ambient: sp.Matrix,
    selected_walls: tuple[int, ...],
    subtree: frozenset[tuple[int, ...]],
) -> tuple[sp.Matrix, ...]:
    selected_roots = tuple(chamber_roots[index] for index in selected_walls)
    internal = set()
    for word in subtree:
        if not word:
            continue
        prefix, letter = word[:-1], word[-1]
        internal.add(
            _oriented_root_key(
                _transform_word(prefix, selected_roots[letter], selected_roots, ambient)
            )
        )
    boundary: dict[tuple[sp.Expr, ...], sp.Matrix] = {}
    for word in sorted(subtree, key=lambda value: (len(value), value)):
        for root in chamber_roots:
            image = _transform_word(word, root, selected_roots, ambient)
            key = _oriented_root_key(image)
            if key not in internal:
                boundary.setdefault(key, image)
    return tuple(boundary.values())


def _allcock_subtree_audit(
    target: tuple[tuple[int, ...], ...]
) -> tuple[dict[str, Any], dict[str, Any]]:
    ambient = _allcock_gram()
    f_roots = _unit_roots(21)
    f_relations = _relation_matrix(f_roots, ambient)
    q_roots = _double_across(f_roots, ambient, 3)
    q_relations = _relation_matrix(q_roots, ambient)
    f_sets = _maximal_disjoint_doubling_sets(f_relations)
    q_sets = _maximal_disjoint_doubling_sets(q_relations)
    f_costs = {
        str(wall): _gluing_cost(f_relations, wall)
        for wall in _doubling_walls(f_relations)
    }
    q_costs = {
        str(wall): _gluing_cost(q_relations, wall)
        for wall in _doubling_walls(q_relations)
    }

    if not (
        f_sets == [[1], [3, 4], [3, 8], [4, 7], [7, 8], [9]]
        and q_sets
        == [
            [1, 5],
            [1, 11],
            [2, 4],
            [2, 10],
            [4, 5, 9],
            [4, 9, 11],
            [5, 9, 10],
            [9, 10, 11],
            [12],
        ]
    ):
        raise ArithmeticError("Allcock disjoint-doubling-wall inventory changed")

    # Q has 25 facets.  The tree formula is
    #   facets(Q_T) = 25 I - sum_i c_i E_i.
    # For I>=3, c_i<=20 gives at least 5I+20 >=35 facets.  At I=2 the
    # only 31-facet possibilities are the two walls of cost 19.
    q_candidates = []
    for wall in _doubling_walls(q_relations):
        child = _double_across(q_roots, ambient, wall)
        if len(child) == 31:
            relations = _relation_matrix(child, ambient)
            q_candidates.append(
                {
                    "wall_zero_based": wall,
                    "relation_histogram": _relation_histogram(relations),
                    "wl_hash": _graph_hash(relations),
                    "target_isomorphic": _same_relation_graph(relations, target),
                }
            )

    # F has 21 facets and no three disjoint doubling walls.  Every genuine
    # two-colour Cayley tree is a path; all four possible colour pairs have
    # cost 17 on both colours.  Hence its facet count is
    #   21I-17(I-1)=4I+17,
    # which is never 31.  The two old four-copy records were artifacts of a
    # Figure 3 transcript that omitted the heavy H-F and I-G bonds.
    f_candidates: list[dict[str, Any]] = []

    named_selected = (4, 5, 9)
    checks = {
        "F_has_21_facets": len(f_roots) == 21,
        "Q_is_the_25_facet_double_of_F_across_14_bar_203": len(q_roots) == 25,
        "Allcock_named_Q_walls_are_pairwise_disjoint_doubling_walls": (
            all(wall in _doubling_walls(q_relations) for wall in named_selected)
            and all(
                q_relations[left][right] == 0
                for left, right in combinations(named_selected, 2)
            )
        ),
        "all_Q_three_wall_cliques_frozen": [
            row for row in q_sets if len(row) == 3
        ]
        == [[4, 5, 9], [4, 9, 11], [5, 9, 10], [9, 10, 11]],
        "Q_only_31_facet_tree_candidates_are_two_doubles": (
            [row["wall_zero_based"] for row in q_candidates] == [1, 2]
            and all(
                row["relation_histogram"]
                == {"0": 110, "2": 321, "3": 20, "4": 14}
                for row in q_candidates
            )
        ),
        "Q_31_facet_candidates_do_not_match": all(
            not row["target_isomorphic"] for row in q_candidates
        ),
        "F_facet_formula_has_no_31_solution": not f_candidates
        and all(
            _gluing_cost(f_relations, left)
            == _gluing_cost(f_relations, right)
            == 17
            for left, right in ((3, 4), (3, 8), (4, 7), (7, 8))
        ),
        "F_31_facet_candidates_do_not_match": all(
            not row["target_isomorphic"] for row in f_candidates
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"Allcock subtree audit changed: {checks}")
    record = {
        "source": "Allcock 2006, Figure 3 and section 'Variations on doubling'",
        "source_transcript_sha256": {
            "main_a_tex": "74c1fcadaf0b2229e566044761acd35e3405db7011c961b83753c956ca6926f7",
            "figure_3_pdf": "347f53c2693172a6f967322fd249d58be88e8c6d22ee891a6437308b4d088c1e",
        },
        "transcript_correction": (
            "heavy bonds H-F and I-G, visible in the pinned vector Figure 3, "
            "were omitted from the previous project transcript"
        ),
        "facet_count_formula": (
            "m*I - sum_i c_i*E_i, c_i=2+(number of other facets "
            "orthogonal to selected wall i)"
        ),
        "F": {
            "facets": 21,
            "doubling_walls_zero_based": _doubling_walls(f_relations),
            "maximal_pairwise_disjoint_sets_zero_based": f_sets,
            "gluing_costs_by_wall_zero_based": f_costs,
            "31_facet_tree_candidates": f_candidates,
            "two_colour_facet_formula": "21I-17(I-1)=4I+17, never 31",
        },
        "Q": {
            "construction": "double F across zero-based wall 3 = 14|203",
            "facets": 25,
            "doubling_walls_zero_based": _doubling_walls(q_relations),
            "maximal_pairwise_disjoint_sets_zero_based": q_sets,
            "gluing_costs_by_wall_zero_based": q_costs,
            "Allcock_named_selected_walls_zero_based": list(named_selected),
            "Allcock_named_selected_gluing_costs": [
                _gluing_cost(q_relations, wall) for wall in named_selected
            ],
            "31_facet_tree_candidates": q_candidates,
            "I_at_least_3_lower_bound": "25I-20(I-1)=5I+20 >= 35",
        },
        "conclusion": (
            "the target is not in any subtree-gluing family obtained from a "
            "pairwise-disjoint set of doubling walls of the published H13 F or Q"
        ),
        "checks": checks,
    }
    auxiliary = {"ambient": ambient, "F_roots": f_roots}
    return record, auxiliary


def _named_pure_doubling_audits(
    target: tuple[tuple[int, ...], ...], allcock_auxiliary: dict[str, Any]
) -> dict[str, Any]:
    records = [
        _pure_doubling_audit(
            "Allcock_D6D6_F",
            allcock_auxiliary["ambient"],
            allcock_auxiliary["F_roots"],
            target,
        )
    ]
    for name, transcript in COXITER_TRANSCRIPTS.items():
        ambient = _seed_gram(transcript["facets"], transcript["edges"])
        records.append(
            _pure_doubling_audit(
                name, ambient, _unit_roots(transcript["facets"]), target
            )
        )
    by_name = {record["seed"]: record for record in records}
    expected_classes = {
        "Allcock_D6D6_F": {
            "21": 1,
            "24": 1,
            "25": 1,
            "27": 1,
            "30": 2,
            "31": 1,
        },
        "Vinberg_f1_P13": {"15": 2},
        "Vinberg_f2_H13": {
            "19": 1,
            "21": 1,
            "22": 1,
            "25": 1,
            "26": 2,
            "27": 1,
            "28": 1,
        },
        "McLeod_f3_H13": {"22": 1, "27": 3},
    }
    checks = {
        "all_named_pure_double_closures_exclude_target": all(
            record["target_absent"] for record in records
        ),
        "Allcock_class_counts_frozen": by_name["Allcock_D6D6_F"][
            "isomorphism_classes_by_facet_count"
        ]
        == expected_classes["Allcock_D6D6_F"],
        "Vinberg_f1_class_counts_frozen": by_name["Vinberg_f1_P13"][
            "isomorphism_classes_by_facet_count"
        ]
        == expected_classes["Vinberg_f1_P13"],
        "Vinberg_f2_class_counts_frozen": by_name["Vinberg_f2_H13"][
            "isomorphism_classes_by_facet_count"
        ]
        == expected_classes["Vinberg_f2_H13"],
        "McLeod_f3_class_counts_frozen": by_name["McLeod_f3_H13"][
            "isomorphism_classes_by_facet_count"
        ]
        == expected_classes["McLeod_f3_H13"],
        "Tumarkin_three_n_plus_2_seeds_have_no_31_facet_descendant": all(
            by_name[name]["target_absent"]
            and "31" not in by_name[name]["isomorphism_classes_by_facet_count"]
            for name in ("Tumarkin_13_01", "Tumarkin_13_02", "Tumarkin_13_03")
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"named pure-doubling audit changed: {checks}")
    return {
        "scope": (
            "all pure iterated doubles with at most 31 facets of the embedded "
            "Allcock F, Vinberg f1/f2, McLeod f3, and three Tumarkin n+2 seeds"
        ),
        "records": records,
        "checks": checks,
    }


def _source_and_scope() -> dict[str, Any]:
    return {
        "primary_sources": [
            {
                "name": "Allcock, Infinitely many hyperbolic Coxeter groups through dimension 19",
                "publication": "Geometry & Topology 10 (2006), 737-758",
                "url": "https://doi.org/10.2140/gt.2006.10.737",
                "use": "H13 D6D6 face, its double Q, and arbitrary-subtree theorem",
            },
            {
                "name": "Tumarkin, Hyperbolic Coxeter n-polytopes with n+2 facets",
                "publication": "Mathematical Notes 75 (2004), 848-854",
                "url": "https://arxiv.org/abs/math/0301133",
                "use": "complete n+2 classification and the H13/15 denominator diagram",
            },
            {
                "name": "McLeod, Hyperbolic reflection groups associated to the quadratic forms -3x0^2+x1^2+...+xn^2",
                "url": "https://arxiv.org/abs/1007.2299",
                "use": "the H13 f3 seed",
            },
            {
                "name": "Scharlau--Kirschmer Lorentzian lattice tables",
                "url": "https://www.math.rwth-aachen.de/~Markus.Kirschmer/lorentz/index.html",
                "use": "List13 entry 4 maximal reflection chamber census",
            },
        ],
        "pinned_machine_transcript": {
            "repository": "CoxIter",
            "commit": "b800be48d8240fbabd2d1c9542e61054584ad7fe",
            "embedded_file_sha256": {
                name: record["sha256"] for name, record in COXITER_TRANSCRIPTS.items()
            },
        },
        "proved_scope": [
            "the target is a literal finite-volume noncompact H13/31 Coxeter polytope",
            "its unique reflection-stable maximal integral overclosure is explicitly isometric to the quotient lattice of Tumarkin 13-tum04_13_01",
            "the named pure-doubling closures through 31 facets do not contain its Coxeter graph",
            "Allcock's published H13 F/Q arbitrary-subtree families do not contain its Coxeter graph",
        ],
        "not_proved": [
            "worldwide priority of the literal 31-wall Coxeter chamber",
            "absence from arbitrary finite-index reflection subgroups of the 15-generator denominator",
            "a numerical subgroup index in the maximal reflection group",
            "lifting the selected diagram involution to an automorphism of a K3 surface",
            "a quotient divisor cone or any surface-automorphism interpretation",
        ],
    }


def build(
    construction_path: Path = CONSTRUCTION,
    vinberg_path: Path = VINBERG,
    stable_path: Path = STABLE,
    earlier_lineage_path: Path = EARLIER_LINEAGE,
    compiler: str | None = None,
) -> dict[str, Any]:
    construction = json.loads(construction_path.read_text(encoding="utf-8"))
    vinberg = json.loads(vinberg_path.read_text(encoding="utf-8"))
    stable = json.loads(stable_path.read_text(encoding="utf-8"))
    earlier = json.loads(earlier_lineage_path.read_text(encoding="utf-8"))
    if not all(
        payload["all_checks_pass"] for payload in (construction, vinberg, stable, earlier)
    ):
        raise ArithmeticError("an upstream A18 certificate failed")
    construction_case = _case(construction)
    vinberg_case = _case(vinberg)
    stable_case = _case(stable)
    target_relations = _target_relations(construction_case)

    target = _target_invariants(construction_case, vinberg_case, stable_case)
    denominator, auxiliary = _tumarkin_denominator(
        construction_case, stable_case, compiler
    )
    subgroup = _subgroup_presentation_handoff(construction_case, auxiliary)
    allcock, allcock_auxiliary = _allcock_subtree_audit(target_relations)
    pure = _named_pure_doubling_audits(target_relations, allcock_auxiliary)
    checks = {
        "target_exact_invariants_pass": all(target["checks"].values()),
        "Tumarkin_denominator_identified": all(denominator["checks"].values()),
        "exact_subgroup_word_handoff_complete": subgroup[
            "exact_word_checks_pass"
        ]
        and len(subgroup["target_reflection_conjugates"]) == 31,
        "Allcock_F_Q_subtree_families_excluded": all(allcock["checks"].values()),
        "named_pure_doubling_families_excluded": all(pure["checks"].values()),
        "no_worldwide_priority_claim": "worldwide priority"
        in _source_and_scope()["not_proved"][0],
        "no_surface_automorphism_claim": (
            construction_case["geometry_scope"][
                "geometric_automorphism_membership_certified"
            ]
            is False
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"global H13/31 lineage audit changed: {checks}")
    return {
        "schema": SCHEMA,
        "sources": {
            "construction_path": str(construction_path.relative_to(ROOT)),
            "construction_sha256": _sha256(construction_path),
            "vinberg_path": str(vinberg_path.relative_to(ROOT)),
            "vinberg_sha256": _sha256(vinberg_path),
            "stable_overlattice_path": str(stable_path.relative_to(ROOT)),
            "stable_overlattice_sha256": _sha256(stable_path),
            "earlier_bounded_lineage_path": str(
                earlier_lineage_path.relative_to(ROOT)
            ),
            "earlier_bounded_lineage_sha256": _sha256(earlier_lineage_path),
            "vinberg_cpp_path": str(CPP_SOURCE.relative_to(ROOT)),
            "vinberg_cpp_sha256": _sha256(CPP_SOURCE),
        },
        "target": target,
        "maximal_denominator": denominator,
        "subgroup_presentation_handoff": subgroup,
        "Allcock_H13_subtree_gluing": allcock,
        "named_pure_doubling_closures": pure,
        "priority_scope": _source_and_scope(),
        "candidate_assessment": {
            "status": "strongest current worldwide-open literal-chamber candidate in this project",
            "reason": (
                "H13, only 31 facets, exact volume/lattice/denominator data, and "
                "the principal published H13 doubling/subtree lineage now closed"
            ),
            "careful_claim": (
                "newly extracted in this project; no matching diagram in the stated "
                "bounded corpora and named families; worldwide literal priority open"
            ),
        },
        "checks": checks,
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--construction", type=Path, default=CONSTRUCTION)
    parser.add_argument("--vinberg", type=Path, default=VINBERG)
    parser.add_argument("--stable", type=Path, default=STABLE)
    parser.add_argument("--earlier-lineage", type=Path, default=EARLIER_LINEAGE)
    parser.add_argument("--compiler")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build(
        args.construction,
        args.vinberg,
        args.stable,
        args.earlier_lineage,
        args.compiler,
    )
    if args.check_archive:
        archived = json.loads(args.output.read_text(encoding="utf-8"))
        if payload != archived:
            raise ArithmeticError("rebuilt H13/31 global-lineage archive changed")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["checks"], sort_keys=True))
    print(args.output)


if __name__ == "__main__":
    main()
