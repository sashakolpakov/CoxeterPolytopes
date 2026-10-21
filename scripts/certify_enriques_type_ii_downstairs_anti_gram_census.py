#!/usr/bin/env python3
"""Exhaust the Type-II two-sheet anti-Grams from downstairs data alone.

For two lifted Enriques curves over roots with downstairs intersection ``m``,
nonnegative same/cross intersections force

``delta_i.delta_j in {-2m,-2m+4,...,2m}``.

Shimada's Type-II (``12B``) graph has fifteen ``m=1`` pairs.  This script
constructs all ``2^15`` assignments, quotients exactly by independent sheet
swaps, computes exact inertia on every orbit representative, and applies the
half-sum glue/rootlessness tests visible from the downstairs wall vectors.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
from dataclasses import asdict
import hashlib
import itertools
import json
import math
from pathlib import Path
import shutil
import subprocess
from typing import Any, Iterable

import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form

try:
    from scripts.certify_enriques_variable_lift_package import parse_gap_records
except ModuleNotFoundError:  # direct execution
    from certify_enriques_variable_lift_package import parse_gap_records
from coxeter_polytopes.exact.matrix import inertia


ROOT = Path(__file__).resolve().parents[1]
SHIMADA_DATA = (
    ROOT / "data" / "external" / "shimada_l10l26note" / "L10L26compdata2.txt"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "enriques_type_ii_downstairs_anti_gram_census.json"
)
SHIMADA_SHA256 = "218c8188e467bde29369cdeda3771d767c3f70105d8fb56def8d698b90e75a6c"
BRANDHORST_SHIMADA = (
    ROOT
    / "data"
    / "external"
    / "brandhorst_shimada_aut_enriques"
    / "AutEnrVol.pdf"
)
BRANDHORST_SHIMADA_SHA256 = (
    "e11f06d376996bbb2aa3221c6e67ac82a9fa99e9e56096ebc35279cf206e610b"
)
MARTIN_PDF = (
    ROOT / "data" / "external" / "martin_enriques_2019" / "Martin2019.pdf"
)
MARTIN_PDF_SHA256 = (
    "7831fbf14ed90f4a30eec68a989dce660724e6384bd1da3f108d373929b905a9"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _matrix_key(matrix: Iterable[Iterable[int]]) -> tuple[int, ...]:
    rows = [list(map(int, row)) for row in matrix]
    return tuple(
        rows[left][right]
        for left in range(len(rows))
        for right in range(left, len(rows))
    )


def _matrix_from_key(key: tuple[int, ...], size: int) -> list[list[int]]:
    matrix = [[0] * size for _ in range(size)]
    cursor = 0
    for left in range(size):
        for right in range(left, size):
            matrix[left][right] = matrix[right][left] = int(key[cursor])
            cursor += 1
    if cursor != len(key):
        raise ArithmeticError("bad upper-triangle key length")
    return matrix


def _switch_canonical(
    matrix: list[list[int]],
) -> tuple[tuple[int, ...], list[list[int]], list[int]]:
    """Force sorted-BFS tree edges negative in every support component."""

    size = len(matrix)
    adjacency = [
        [other for other in range(size) if other != vertex and matrix[vertex][other]]
        for vertex in range(size)
    ]
    signs: list[int | None] = [None] * size
    for root in range(size):
        if signs[root] is not None:
            continue
        signs[root] = 1
        queue = [root]
        while queue:
            vertex = queue.pop(0)
            for other in sorted(adjacency[vertex]):
                if signs[other] is None:
                    edge_sign = 1 if matrix[vertex][other] > 0 else -1
                    signs[other] = -int(signs[vertex]) * edge_sign
                    queue.append(other)
    eps = [int(value) for value in signs]
    canonical = [
        [eps[left] * matrix[left][right] * eps[right] for right in range(size)]
        for left in range(size)
    ]
    return _matrix_key(canonical), canonical, eps


def _switch_orbit(matrix: list[list[int]]) -> set[tuple[int, ...]]:
    size = len(matrix)
    orbit = set()
    for initial in itertools.product((-1, 1), repeat=size - 1):
        signs = initial + (1,)
        switched = [
            [signs[left] * matrix[left][right] * signs[right] for right in range(size)]
            for left in range(size)
        ]
        orbit.add(_matrix_key(switched))
    return orbit


def _key_hash(keys: Iterable[tuple[int, ...]]) -> str:
    rendered = json.dumps(sorted(keys), separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(rendered).hexdigest()


def _gf2_rank(rows: Iterable[Iterable[int]]) -> int:
    packed = []
    width = 0
    for raw in rows:
        row = [int(value) & 1 for value in raw]
        width = max(width, len(row))
        packed.append(sum(value << column for column, value in enumerate(row)))
    rank = 0
    for column in range(width):
        pivot = next((index for index in range(rank, len(packed)) if (packed[index] >> column) & 1), None)
        if pivot is None:
            continue
        packed[rank], packed[pivot] = packed[pivot], packed[rank]
        for index in range(len(packed)):
            if index != rank and ((packed[index] >> column) & 1):
                packed[index] ^= packed[rank]
        rank += 1
    return rank


def _gf2_kernel_basis_of_transpose(rows: list[list[int]]) -> list[list[int]]:
    """Basis of ``{z: z*rows=0}`` over F2, checked by exhaustive replay."""

    size = len(rows)
    kernel = []
    for bits in itertools.product((0, 1), repeat=size):
        if all(
            sum(bits[row] * rows[row][column] for row in range(size)) % 2 == 0
            for column in range(len(rows[0]))
        ):
            kernel.append(list(bits))
    basis: list[list[int]] = []
    rank = 0
    for vector in kernel:
        trial = basis + [vector]
        trial_rank = _gf2_rank(trial)
        if trial_rank > rank:
            basis.append(vector)
            rank = trial_rank
    generated = {
        tuple(
            sum(coefficients[index] * basis[index][column] for index in range(len(basis))) % 2
            for column in range(size)
        )
        for coefficients in itertools.product((0, 1), repeat=len(basis))
    }
    if generated != {tuple(vector) for vector in kernel}:
        raise ArithmeticError("GF(2) kernel basis replay failed")
    return basis


def _smith_diagonal(matrix: sp.Matrix) -> list[int]:
    smith = smith_normal_form(matrix, domain=sp.ZZ)
    return [
        abs(int(smith[index, index]))
        for index in range(min(smith.rows, smith.cols))
        if smith[index, index]
    ]


def _gp_matrix(matrix: sp.Matrix) -> str:
    return "[" + ";".join(
        ",".join(str(int(matrix[row, column])) for column in range(matrix.cols))
        for row in range(matrix.rows)
    ) + "]"


def _pari_short_vectors(form: sp.Matrix, bound: int) -> tuple[list[list[int]], dict[str, Any]]:
    """Enumerate one representative from each +/- pair with PARI qfminim."""

    gp = shutil.which("gp")
    if gp is None:
        raise RuntimeError("PARI/GP is required for the exact minimum test")
    if form != form.T or any(value.q != 1 for value in form):
        raise ArithmeticError("qfminim received a nonintegral nonsymmetric form")
    exact = inertia([[int(form[i, j]) for j in range(form.cols)] for i in range(form.rows)])
    if exact.positive != form.rows or exact.negative or exact.zero:
        raise ArithmeticError("qfminim form is not positive definite")
    program = (
        f"A={_gp_matrix(form)};"
        f"R=qfminim(A,{int(bound)});"
        "print(R[1]);print(R[2]);print(matsize(R[3]));"
        "V=R[3];for(j=1,matsize(V)[2],print(Vec(V[,j])));quit()\n"
    )
    completed = subprocess.run(
        [gp, "-fq"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    if len(lines) < 3:
        raise ArithmeticError(f"truncated qfminim output: {completed.stdout!r}")
    signed_count = int(lines[0])
    maximum = int(lines[1])
    dimensions = ast.literal_eval(lines[2])
    vectors = [list(map(int, ast.literal_eval(line))) for line in lines[3:]]
    allowed_dimensions = [[form.rows, len(vectors)]]
    if not vectors:
        allowed_dimensions.append([0, 0])
    if dimensions not in allowed_dimensions:
        raise ArithmeticError(f"qfminim dimensions changed: {completed.stdout!r}")
    if signed_count != 2 * len(vectors):
        raise ArithmeticError("PARI sign convention changed")
    norms = [
        int((sp.Matrix([vector]) * form * sp.Matrix([vector]).T)[0, 0])
        for vector in vectors
    ]
    if any(norm > bound for norm in norms) or (norms and max(norms) != maximum):
        raise ArithmeticError("qfminim exact replay failed")
    return vectors, {
        "engine": "PARI/GP qfminim exact positive-definite enumeration",
        "bound": bound,
        "signed_vector_count": signed_count,
        "vectors_modulo_sign_count": len(vectors),
        "maximum_returned_norm": maximum,
        "exact_replayed_norm_histogram_modulo_sign": {
            str(norm): count for norm, count in sorted(Counter(norms).items())
        },
    }


def _generated_anti_lattice(
    delta_gram: list[list[int]], walls: list[list[int]]
) -> dict[str, Any]:
    """Compute the complete anti lattice inside ``<M,b_i>``.

    Writing ``b_i=(u_i+delta_i)/2``, an anti vector has the form
    ``sum z_i b_i-(sum z_i u_i)/2`` and is integral precisely when
    ``z*W=0 mod 2``.  The calculation below constructs that congruence
    lattice, maps it to the rank-nine delta span, and takes its exact HNF.
    """

    delta = sp.Matrix(delta_gram)
    wall_matrix = sp.Matrix(walls)
    kernel_mod_2 = _gf2_kernel_basis_of_transpose(walls)
    generators = sp.Matrix(
        [[2 if row == column else 0 for column in range(12)] for row in range(12)]
        + kernel_mod_2
    )
    congruence_basis = hermite_normal_form(generators.T).T
    if congruence_basis.shape != (12, 12):
        raise ArithmeticError("congruence HNF changed shape")
    if abs(int(congruence_basis.det())) != 2 ** _gf2_rank(walls):
        raise ArithmeticError("congruence lattice index disagrees with GF(2) rank")
    if any(
        int((congruence_basis * wall_matrix)[row, column]) % 2
        for row in range(12)
        for column in range(10)
    ):
        raise ArithmeticError("congruence basis violates zW=0 mod 2")

    independent = next(
        subset
        for subset in itertools.combinations(range(12), 9)
        if delta.extract(subset, subset).det()
    )
    basis_gram = delta.extract(independent, independent)
    coordinates = sp.Matrix.vstack(
        *[
            delta.extract([row], independent) * basis_gram.inv()
            for row in range(12)
        ]
    )
    if coordinates.extract(independent, range(9)) != sp.eye(9):
        raise ArithmeticError("delta coordinate factorization changed")
    if coordinates * basis_gram * coordinates.T != delta:
        raise ArithmeticError("delta coordinate factorization failed")

    rational_generators = congruence_basis * coordinates / 2
    denominator = math.lcm(
        *(int(value.q) for value in rational_generators)
    )
    integral_generators = rational_generators * denominator
    if any(value.q != 1 for value in integral_generators):
        raise ArithmeticError("coordinate denominator clearing failed")
    integral_generators = integral_generators.applyfunc(int)
    image_hnf = hermite_normal_form(integral_generators.T).T
    if image_hnf.shape != (9, 9):
        raise ArithmeticError("anti-image HNF changed shape")
    anti_basis_coordinates = image_hnf / denominator
    anti_gram = anti_basis_coordinates * basis_gram * anti_basis_coordinates.T
    if any(value.q != 1 for value in anti_gram):
        raise ArithmeticError("generated anti lattice is not integral")
    positive_form = -anti_gram
    roots, minimum_test = _pari_short_vectors(positive_form, 2)
    minimum_test["rootless_through_norm_two"] = not roots
    scaled_classification: dict[str, Any]
    if roots:
        root_rows = sp.Matrix(roots)
        divided_form = positive_form / 2
        if any(value.q != 1 for value in divided_form):
            raise ArithmeticError("rootful candidate is not two times an integral form")
        root_pairing = root_rows * divided_form * root_rows.T
        if root_rows.shape != (9, 9) or abs(int(root_rows.det())) != 1 or root_pairing != sp.eye(9):
            raise ArithmeticError("rootful candidate did not identify as the odd cubic lattice")
        scaled_classification = {
            "isometry_class": "<-2>^9",
            "proof": (
                "the nine PARI sign representatives form a unimodular basis and "
                "their Gram for (-A)/2 is I9"
            ),
            "unimodular_basis_rows_in_archived_anti_basis": roots,
            "basis_determinant_absolute": abs(int(root_rows.det())),
            "divided_basis_gram": [
                [int(root_pairing[i, j]) for j in range(9)] for i in range(9)
            ],
        }
    else:
        divided_form = positive_form / 2
        if any(value.q != 1 for value in divided_form):
            raise ArithmeticError("rootless candidate is not two times an integral form")
        norm_four_representatives, norm_four_test = _pari_short_vectors(positive_form, 4)
        signed_roots = {
            tuple(vector) for vector in norm_four_representatives
        } | {
            tuple(-value for value in vector) for vector in norm_four_representatives
        }
        if len(signed_roots) != 144:
            raise ArithmeticError("rootless candidate does not have the D9 root census")
        functional_base = next(
            base
            for base in range(2, 100)
            if all(
                sum((base**index) * value for index, value in enumerate(root))
                for root in signed_roots
            )
        )
        positive_roots = {
            root
            for root in signed_roots
            if sum((functional_base**index) * value for index, value in enumerate(root)) > 0
        }
        simple_roots = sorted(
            root
            for root in positive_roots
            if not any(
                tuple(root[index] - summand[index] for index in range(9)) in positive_roots
                for summand in positive_roots
            )
        )
        simple_rows = sp.Matrix(simple_roots)
        simple_gram = simple_rows * divided_form * simple_rows.T
        if simple_rows.shape != (9, 9) or abs(int(simple_rows.det())) != 1:
            raise ArithmeticError("D9 simple roots do not form a unimodular basis")
        if any(simple_gram[i, i] != 2 for i in range(9)) or any(
            simple_gram[i, j] not in (0, -1)
            for i in range(9)
            for j in range(9)
            if i != j
        ):
            raise ArithmeticError("simple-root Gram is not simply laced")
        adjacency = [
            [j for j in range(9) if j != i and simple_gram[i, j] == -1]
            for i in range(9)
        ]
        branch = [i for i in range(9) if len(adjacency[i]) == 3]
        leaves = [i for i in range(9) if len(adjacency[i]) == 1]
        if len(branch) != 1 or len(leaves) != 3 or sum(len(row) for row in adjacency) != 16:
            raise ArithmeticError("simple-root graph is not a D9 tree")
        distances = []
        for leaf in leaves:
            distance = 0
            previous = None
            vertex = leaf
            while vertex != branch[0]:
                next_vertices = [item for item in adjacency[vertex] if item != previous]
                if len(next_vertices) != 1:
                    raise ArithmeticError("bad D9 arm")
                previous, vertex = vertex, next_vertices[0]
                distance += 1
            distances.append(distance)
        if sorted(distances) != [1, 1, 6]:
            raise ArithmeticError("simple-root graph has the wrong D9 arm lengths")
        scaled_classification = {
            "isometry_class": "D9(-2)",
            "proof": (
                "(-A)/2 has 144 roots; a generic positive system has nine simple "
                "roots forming a unimodular basis with Dynkin arms 1,1,6"
            ),
            "norm_four_enumeration_in_A": norm_four_test,
            "signed_norm_four_vector_count": len(signed_roots),
            "positive_system_functional_base": functional_base,
            "simple_root_basis_rows_in_archived_anti_basis": [
                list(root) for root in simple_roots
            ],
            "basis_determinant_absolute": abs(int(simple_rows.det())),
            "divided_simple_root_gram": [
                [int(simple_gram[i, j]) for j in range(9)] for i in range(9)
            ],
            "dynkin_arm_lengths": sorted(distances),
        }
    return {
        "description": "anti lattice in the half-sum-glued lattice generated by M and the twelve chosen lifts",
        "formula": "A_gen={sum z_i*delta_i/2 : z in Z^12 and z*W=0 mod 2}",
        "wall_kernel_mod_2_dimension": len(kernel_mod_2),
        "wall_kernel_mod_2_basis": kernel_mod_2,
        "congruence_lattice_index_in_Z12": abs(int(congruence_basis.det())),
        "delta_basis_indices_one_based": [index + 1 for index in independent],
        "coordinate_common_denominator_before_hnf": denominator,
        "gram": [[int(anti_gram[i, j]) for j in range(9)] for i in range(9)],
        "rank": 9,
        "absolute_determinant": abs(int(anti_gram.det())),
        "smith_diagonal": _smith_diagonal(anti_gram),
        "even": all(int(anti_gram[i, i]) % 2 == 0 for i in range(9)),
        "minimum_test": minimum_test,
        "minus_two_roots_modulo_sign": roots,
        "scaled_lattice_classification": scaled_classification,
    }


def _binary_half_sum_roots(
    gram: list[list[int]], walls: list[list[int]]
) -> list[dict[str, Any]]:
    """Return all binary glue vectors yielding anti-invariant (-2)-roots."""

    roots = []
    for bits in itertools.product((0, 1), repeat=len(walls)):
        if not any(bits):
            continue
        downstairs_sum = [
            sum(bits[row] * walls[row][column] for row in range(len(walls)))
            for column in range(len(walls[0]))
        ]
        if any(value % 2 for value in downstairs_sum):
            continue
        numerator = sum(
            bits[left] * gram[left][right] * bits[right]
            for left in range(len(walls))
            for right in range(len(walls))
        )
        if numerator == -8:
            roots.append(
                {
                    "support_one_based": [index + 1 for index, bit in enumerate(bits) if bit],
                    "downstairs_half_sum": [value // 2 for value in downstairs_sum],
                    "four_times_square": numerator,
                    "square": -2,
                }
            )
    return roots


def _ternary_half_sum_root_witnesses(
    gram: list[list[int]], walls: list[list[int]]
) -> dict[str, Any]:
    """Independent small-coefficient search for explicit saturation roots."""

    witnesses = []
    for coefficients in itertools.product((-1, 0, 1), repeat=len(walls)):
        first = next((value for value in coefficients if value), 0)
        if first != 1:  # remove zero and the global +/- duplication
            continue
        downstairs_sum = [
            sum(coefficients[row] * walls[row][column] for row in range(len(walls)))
            for column in range(len(walls[0]))
        ]
        if any(value % 2 for value in downstairs_sum):
            continue
        numerator = sum(
            coefficients[left] * gram[left][right] * coefficients[right]
            for left in range(len(walls))
            for right in range(len(walls))
        )
        if numerator == -8:
            witnesses.append(
                {
                    "coefficients": list(coefficients),
                    "support_one_based": [
                        index + 1 for index, value in enumerate(coefficients) if value
                    ],
                    "downstairs_sum": downstairs_sum,
                    "downstairs_half_sum": [value // 2 for value in downstairs_sum],
                    "four_times_square": numerator,
                    "square": -2,
                }
            )
    witnesses.sort(key=lambda item: (len(item["support_one_based"]), item["coefficients"]))
    witness_hash = hashlib.sha256(
        json.dumps(
            [item["coefficients"] for item in witnesses], separators=(",", ":")
        ).encode("utf-8")
    ).hexdigest()
    minimum_support = len(witnesses[0]["support_one_based"]) if witnesses else None
    minimum_witnesses = (
        [item for item in witnesses if len(item["support_one_based"]) == minimum_support]
        if witnesses
        else []
    )
    return {
        "coefficient_box": [-1, 1],
        "global_sign_quotiented": True,
        "witness_count_modulo_sign": len(witnesses),
        "all_witness_coefficients_sha256": witness_hash,
        "minimum_support_size": minimum_support,
        "minimum_support_witness_count_modulo_sign": len(minimum_witnesses),
        "first_three_minimum_support_witnesses": minimum_witnesses[:3],
    }


def analyze() -> dict[str, Any]:
    if _sha256(SHIMADA_DATA) != SHIMADA_SHA256:
        raise AssertionError("Shimada input changed")
    if _sha256(BRANDHORST_SHIMADA) != BRANDHORST_SHIMADA_SHA256:
        raise AssertionError("Brandhorst--Shimada primary source changed")
    if _sha256(MARTIN_PDF) != MARTIN_PDF_SHA256:
        raise AssertionError("Martin primary source changed")
    record = parse_gap_records(SHIMADA_DATA)["12B"]
    configuration = sp.Matrix(record["configuration"])
    walls_matrix = sp.Matrix(record["walls"])
    gram_l10 = sp.Matrix(record["GramL10"])
    walls = [[int(value) for value in row] for row in record["walls"]]
    if configuration.shape != (12, 12) or walls_matrix.shape != (12, 10):
        raise ArithmeticError("Type-II datum has changed shape")
    if walls_matrix * gram_l10 * walls_matrix.T != configuration:
        raise ArithmeticError("wall Gram does not rebuild the configuration")
    edges = [
        (left, right)
        for left in range(12)
        for right in range(left + 1, 12)
        if configuration[left, right] == 1
    ]
    if len(edges) != 15 or any(
        configuration[left, right] not in (0, 1)
        for left in range(12)
        for right in range(left + 1, 12)
    ):
        raise ArithmeticError("Type-II edge alphabet changed")
    support_seen = {0}
    while True:
        enlarged = support_seen | {
            right if left in support_seen else left
            for left, right in edges
            if left in support_seen or right in support_seen
        }
        if enlarged == support_seen:
            break
        support_seen = enlarged
    if support_seen != set(range(12)):
        raise ArithmeticError("Type-II support graph is no longer connected")

    orbit_multiplicity: Counter[tuple[int, ...]] = Counter()
    representatives: dict[tuple[int, ...], list[list[int]]] = {}
    for values in itertools.product((-2, 2), repeat=len(edges)):
        matrix = [[-4 if left == right else 0 for right in range(12)] for left in range(12)]
        for (left, right), value in zip(edges, values):
            matrix[left][right] = matrix[right][left] = value
        key, representative, _ = _switch_canonical(matrix)
        orbit_multiplicity[key] += 1
        representatives[key] = representative
    if sum(orbit_multiplicity.values()) != 2**15 or len(representatives) != 16:
        raise ArithmeticError("raw denominator changed")

    weighted_inertia: Counter[tuple[int, int, int]] = Counter()
    orbit_inertia: Counter[tuple[int, int, int]] = Counter()
    orbit_records = []
    for key in sorted(representatives):
        representative = representatives[key]
        exact = inertia(representative)
        signature = (exact.positive, exact.negative, exact.zero)
        direct = _switch_orbit(representative)
        multiplicity = orbit_multiplicity[key]
        if len(direct) != 2**11 or multiplicity != 2**11:
            raise ArithmeticError("Type-II switch orbit size changed")
        if {
            _switch_canonical(_matrix_from_key(item, 12))[0] for item in direct
        } != {key}:
            raise ArithmeticError("direct switch orbit escapes its canonical class")
        weighted_inertia[signature] += multiplicity
        orbit_inertia[signature] += 1
        binary_roots = _binary_half_sum_roots(representative, walls)
        generated_anti_lattice = (
            _generated_anti_lattice(representative, walls)
            if exact.positive == 0
            else None
        )
        ternary_witnesses = (
            _ternary_half_sum_root_witnesses(representative, walls)
            if exact.positive == 0
            else None
        )
        orbit_records.append(
            {
                "canonical_gram": representative,
                "canonical_key_sha256": _key_hash([key]),
                "raw_assignment_multiplicity": multiplicity,
                "inertia": asdict(exact),
                "rank": exact.positive + exact.negative,
                "binary_half_sum_minus_two_roots": binary_roots,
                "binary_half_sum_minus_two_root_count": len(binary_roots),
                "generated_anti_lattice": generated_anti_lattice,
                "ternary_half_sum_root_search": ternary_witnesses,
            }
        )

    def histogram(counter: Counter[tuple[int, int, int]]) -> dict[str, int]:
        return {
            f"p{positive}_n{negative}_z{zero}": count
            for (positive, negative, zero), count in sorted(counter.items())
        }

    expected_weighted = Counter(
        {
            (0, 9, 3): 4096,
            (1, 9, 2): 2048,
            (1, 10, 1): 12288,
            (2, 10, 0): 14336,
        }
    )
    expected_orbits = Counter(
        {
            (0, 9, 3): 2,
            (1, 9, 2): 1,
            (1, 10, 1): 6,
            (2, 10, 0): 7,
        }
    )
    if weighted_inertia != expected_weighted or orbit_inertia != expected_orbits:
        raise ArithmeticError("Type-II exact inertia census changed")

    nonpositive = [item for item in orbit_records if item["inertia"]["positive"] == 0]
    survivors = [
        item
        for item in nonpositive
        if item["generated_anti_lattice"]["minimum_test"]["rootless_through_norm_two"]
    ]
    eliminated = [item for item in nonpositive if item not in survivors]
    if not (
        len(nonpositive) == 2
        and all(item["rank"] == 9 for item in nonpositive)
        and len(survivors) == 1
        and len(eliminated) == 1
    ):
        raise ArithmeticError("Type-II negative-semidefinite residue changed")
    survivor = survivors[0]
    killed = eliminated[0]
    survivor_lattice = survivor["generated_anti_lattice"]
    killed_lattice = killed["generated_anti_lattice"]
    if not (
        survivor_lattice["scaled_lattice_classification"]["isometry_class"]
        == "D9(-2)"
        and survivor_lattice["absolute_determinant"] == 2048
        and survivor_lattice["smith_diagonal"] == [2] * 8 + [8]
        and killed_lattice["scaled_lattice_classification"]["isometry_class"]
        == "<-2>^9"
        and killed_lattice["absolute_determinant"] == 512
        and killed_lattice["smith_diagonal"] == [2] * 9
    ):
        raise ArithmeticError("Type-II anti-lattice classification changed")
    if not (
        survivor["binary_half_sum_minus_two_root_count"] == 0
        and killed["binary_half_sum_minus_two_root_count"] == 0
        and survivor["ternary_half_sum_root_search"]["witness_count_modulo_sign"] == 0
        and killed["ternary_half_sum_root_search"]["witness_count_modulo_sign"] == 36
        and killed["ternary_half_sum_root_search"]["minimum_support_size"] == 9
        and killed_lattice["minimum_test"]["signed_vector_count"] == 18
    ):
        raise ArithmeticError("Type-II rootlessness witness census changed")

    checks = {
        "all_32768_downstairs_anti_grams_generated": True,
        "all_assignments_partitioned_into_16_complete_sheet_switch_orbits": True,
        "exact_inertia_computed_on_all_16_representatives": True,
        "exactly_two_negative_semidefinite_orbits_both_rank_nine": True,
        "full_half_sum_congruence_lattice_computed_by_exact_hnf": True,
        "one_orbit_eliminated_by_saturated_minus_two_roots": True,
        "surviving_generated_anti_lattice_is_D9_minus_2": True,
        "survivor_rootlessness_certified_by_exact_short_vector_enumeration": True,
        "D9_identification_has_explicit_unimodular_simple_root_basis": True,
        "primary_source_calibrates_Type_II_as_D9_D9_generic_12B": True,
    }
    return {
        "schema": "enriques_type_ii_downstairs_anti_gram_census_v1",
        "scope": {
            "input": (
                "Shimada 12B downstairs walls, nonnegative intersections of the "
                "two lifted curve sheets, and Enriques half-sum glue/rootlessness"
            ),
            "result": (
                "exactly one sheet-switch orbit survives; its generated anti "
                "lattice is D9(-2), exactly as in the very general Type-II carrier"
            ),
            "geometric_status": (
                "a wallwise forced realization of the known Type-II 12B chamber, "
                "not a new abstract Coxeter polytope"
            ),
            "not_claimed": [
                "a new Enriques surface or a new hyperbolic Coxeter polytope",
                "an independent reconstruction of the period or ample cone",
                "a curve-by-curve lift labelling extracted from Brandhorst--Shimada",
            ],
        },
        "pinned_inputs": {
            "shimada_path": str(SHIMADA_DATA.relative_to(ROOT)),
            "shimada_sha256": SHIMADA_SHA256,
            "record": "12B",
            "brandhorst_shimada_path": str(BRANDHORST_SHIMADA.relative_to(ROOT)),
            "brandhorst_shimada_sha256": BRANDHORST_SHIMADA_SHA256,
            "martin_path": str(MARTIN_PDF.relative_to(ROOT)),
            "martin_sha256": MARTIN_PDF_SHA256,
        },
        "imported_geometric_calibration": {
            "brandhorst_shimada": {
                "primary_source": (
                    "Brandhorst--Shimada, Automorphism Groups of Certain "
                    "Enriques Surfaces, Definitions 1.1 and 1.4, Propositions "
                    "1.2--1.3, Table 1.1 row 184, Section 7.4"
                ),
                "imported_statement": (
                    "a very general Type-II Enriques surface is (D9,D9)-generic "
                    "No.184; SX is the half-sum lattice M_R, its anti complement "
                    "is D9(-2), and its nef chamber is 12B"
                ),
                "wallwise_inference": (
                    "the actual twelve lift differences must occur in this exhaustive "
                    "census; rootlessness eliminates one of the two signature orbits, "
                    "so their Gram is the unique survivor up to sheet swaps"
                ),
            },
            "martin": {
                "primary_source": (
                    "Martin 2019, Type-II section: explicit one-parameter K3-cover "
                    "Weierstrass family and fixed-point-free involution"
                ),
                "role": "independent existence/equation lineage only; no lift Gram imported",
            },
        },
        "downstairs_alphabet": {
            "formula": (
                "for a+c=m, delta_i.delta_j=2(a-c)=4a-2m; "
                "Type II has fifteen m=1 pairs and no larger intersections"
            ),
            "m1_edges_zero_based": [list(edge) for edge in edges],
            "m1_values": [-2, 2],
            "raw_assignment_count": 2**15,
            "wall_span_rank_over_Q": int(walls_matrix.rank()),
            "wall_matrix_rank_mod_2": _gf2_rank(walls),
        },
        "sheet_switch_quotient": {
            "action": "delta_i -> epsilon_i*delta_i, epsilon_i in {+1,-1}",
            "canonicalization": "force sorted-BFS spanning-tree edges negative",
            "connected_support_cycle_rank": len(edges) - 12 + 1,
            "predicted_orbit_count_from_cycle_signs": 2 ** (len(edges) - 12 + 1),
            "orbit_count": len(orbit_records),
            "orbit_size_histogram": {
                str(size): count
                for size, count in sorted(Counter(orbit_multiplicity.values()).items())
            },
            "all_canonical_keys_sha256": _key_hash(representatives),
            "rootlessness_invariance": {
                "coefficient_transport": "z_i -> epsilon_i*z_i",
                "parity_reason": "epsilon_i=1 mod 2, hence zW=0 mod 2 is preserved",
                "norm_reason": "for D'=EDE and z'=zE, z'D'z'^t=zDz^t",
            },
        },
        "exact_inertia_census": {
            "raw_assignment_histogram": histogram(weighted_inertia),
            "sheet_switch_orbit_histogram": histogram(orbit_inertia),
            "orbit_records": orbit_records,
        },
        "negative_semidefinite_orbit_count": len(nonpositive),
        "projection_theorem": {
            "surviving_orbit_count": 1,
            "survivor_canonical_key_sha256": survivor["canonical_key_sha256"],
            "survivor_canonical_gram": survivor["canonical_gram"],
            "survivor_generated_anti_lattice": survivor_lattice,
            "eliminated_orbit_canonical_key_sha256": killed["canonical_key_sha256"],
            "eliminated_orbit_canonical_gram": killed["canonical_gram"],
            "eliminated_generated_anti_lattice": killed_lattice,
            "eliminated_explicit_ternary_witness_search": (
                killed["ternary_half_sum_root_search"]
            ),
            "conclusion": (
                "negative semidefiniteness plus full half-sum saturated rootlessness "
                "forces the D9(-2) anti-Gram orbit up to independent sheet switches"
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
    rendered = json.dumps(analyze(), indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.is_file() or args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archived certificate differs: {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(hashlib.sha256(rendered.encode("utf-8")).hexdigest())


if __name__ == "__main__":
    main()
