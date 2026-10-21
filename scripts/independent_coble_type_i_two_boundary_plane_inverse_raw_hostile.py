#!/usr/bin/env python3
"""Helper-free hostile replay of the raw H11/26 plane inverse.

This audit starts from the odd-unimodular Gram and the 26 primitive wall
vectors.  It independently enumerates every sequential numerical P2
endpoint, derives the four line--quintic branches, and constructs the ten
positive-degree cusp-lift curves which the original root-seed handoff did
not construct.  No project Python module is imported.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import resource
import shutil
import subprocess
import sys
from typing import Any, Callable, Sequence

import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/coble_type_i_two_boundary_full_ns_h11.json"
PRIMARY = ROOT / "scripts/certify_coble_type_i_two_boundary_plane_inverse.py"
OUTPUT = ROOT / "data/enumerations/coble_type_i_two_boundary_plane_inverse_raw_hostile_audit.json"
SCHEMA = "coxeter.coble_type_i_two_boundary_plane_inverse_raw_hostile_audit.v1"

X, Y, Z = sp.symbols("x y z")

TARGET_SIGNATURE = (
    (1, 5),
    ((0, 7), (1, 3)),
    ((-2, 7), (-1, 4)),
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def peak_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / (1024**3 if sys.platform == "darwin" else 1024**2)


def memory_gate() -> None:
    if peak_gib() >= 0.5:
        raise MemoryError("the hostile replay crossed its 0.5 GiB ceiling")


def dot(gram: sp.Matrix, left: Sequence[Any], right: Sequence[Any]) -> sp.Expr:
    return (sp.Matrix(left).T * gram * sp.Matrix(right))[0]


def primitive(vector: sp.Matrix, orient_first: bool = False) -> sp.Matrix:
    denominator = 1
    for entry in vector:
        denominator = int(sp.ilcm(denominator, sp.denom(entry)))
    values = [int(entry * denominator) for entry in vector]
    divisor = math.gcd(*(abs(value) for value in values))
    if divisor == 0:
        raise ArithmeticError("cannot primitive-normalize zero")
    values = [value // divisor for value in values]
    if orient_first and next(value for value in values if value) < 0:
        values = [-value for value in values]
    return sp.Matrix(values)


def rational_json(vector: Sequence[sp.Expr]) -> list[int | str]:
    return [int(entry) if sp.denom(entry) == 1 else str(entry) for entry in vector]


def normalize_raw(source: dict[str, Any]) -> dict[str, Any]:
    gram = sp.Matrix(source["full_NS_lattice"]["gram"])
    walls = [sp.Matrix(vector) for vector in source["facets"]["NS_wall_vectors"]]
    wall_gram = sp.Matrix([[dot(gram, left, right) for right in walls] for left in walls])
    norms = [int(wall_gram[index, index]) for index in range(len(walls))]
    boundaries = [index for index, norm in enumerate(norms) if norm == -4]
    roots = [index for index, norm in enumerate(norms) if norm == -2]
    minus_one = [index for index, norm in enumerate(norms) if norm == -1]
    if boundaries != [0, 1] or roots != list(range(2, 12)) or minus_one != list(range(12, 26)):
        raise ArithmeticError("the intrinsic wall taxonomy changed")
    canonical = -(walls[0] + walls[1]) / 2
    if any(sp.denom(entry) != 1 for entry in canonical):
        raise ArithmeticError("the canonical class is not integral")
    return {
        "gram": gram,
        "walls": walls,
        "wall_gram": wall_gram,
        "norms": norms,
        "boundaries": boundaries,
        "roots": roots,
        "minus_one": minus_one,
        "canonical": canonical,
    }


def plane_class(raw: dict[str, Any], contracted: Sequence[int]) -> sp.Matrix:
    exceptional = sp.Matrix.hstack(*(raw["walls"][index] for index in contracted))
    kernel = (exceptional.T * raw["gram"]).nullspace()
    if len(kernel) != 1:
        raise ArithmeticError("a numerical P2 endpoint has non-line complement")
    plane = primitive(kernel[0])
    if dot(raw["gram"], raw["canonical"], plane) > 0:
        plane = -plane
    if dot(raw["gram"], plane, plane) != 1 or dot(raw["gram"], raw["canonical"], plane) != -3:
        raise ArithmeticError("an endpoint complement is not a P2 line class")
    return plane


def sequential_endpoints(raw: dict[str, Any]) -> dict[int, tuple[tuple[int, ...], tuple[tuple[sp.Expr, ...], ...]]]:
    matrix = raw["wall_gram"]
    allowed = [index for index in range(26) if index not in raw["boundaries"]]
    current: dict[int, tuple[tuple[int, ...], tuple[tuple[sp.Expr, ...], ...]]] = {
        0: (tuple(), tuple())
    }
    for _level in range(11):
        following: dict[int, tuple[tuple[int, ...], tuple[tuple[sp.Expr, ...], ...]]] = {}
        for mask, (contracted, inverse) in current.items():
            size = len(contracted)
            for candidate in allowed:
                if mask & (1 << candidate):
                    continue
                vector = [matrix[index, candidate] for index in contracted]
                product = [
                    sum(inverse[left][right] * vector[right] for right in range(size))
                    for left in range(size)
                ]
                quotient_square = matrix[candidate, candidate] - sum(
                    vector[index] * product[index] for index in range(size)
                )
                if quotient_square != -1:
                    continue
                enlarged = [
                    [inverse[left][right] - product[left] * product[right] for right in range(size)]
                    + [product[left]]
                    for left in range(size)
                ]
                enlarged.append([*product, sp.Integer(-1)])
                temporary = [*contracted, candidate]
                order = sorted(range(size + 1), key=temporary.__getitem__)
                new_contracted = tuple(temporary[index] for index in order)
                new_inverse = tuple(
                    tuple(enlarged[left][right] for right in order) for left in order
                )
                following.setdefault(mask | (1 << candidate), (new_contracted, new_inverse))
        current = following
        memory_gate()
    return current


def endpoint_signature(raw: dict[str, Any], contracted: Sequence[int]) -> tuple[Any, ...]:
    plane = plane_class(raw, contracted)
    boundary_degrees = tuple(
        sorted(int(dot(raw["gram"], raw["walls"][index], plane)) for index in raw["boundaries"])
    )
    root_degrees = tuple(
        sorted(Counter(int(dot(raw["gram"], raw["walls"][index], plane)) for index in raw["roots"]).items())
    )
    norms = tuple(sorted(Counter(raw["norms"][index] for index in contracted).items()))
    return boundary_degrees, root_degrees, norms


def contraction_order(raw: dict[str, Any], target: Sequence[int]) -> list[int]:
    target_set = set(target)
    matrix = raw["wall_gram"]

    def recurse(order: list[int]) -> list[int] | None:
        if len(order) == len(target_set):
            return order
        inverse = None
        if order:
            inverse = sp.Matrix([[matrix[left, right] for right in order] for left in order]).inv()
        for candidate in sorted(target_set - set(order)):
            square = matrix[candidate, candidate]
            if inverse is not None:
                vector = sp.Matrix([matrix[index, candidate] for index in order])
                square -= (vector.T * inverse * vector)[0]
            if square == -1:
                answer = recurse([*order, candidate])
                if answer is not None:
                    return answer
        return None

    answer = recurse([])
    if answer is None:
        raise ArithmeticError("an enumerated endpoint has no exact contraction order")
    return answer


def plane_basis(raw: dict[str, Any], target: Sequence[int]) -> tuple[list[int], sp.Matrix]:
    order = contraction_order(raw, target)
    strict: list[sp.Matrix] = []
    totals: list[sp.Matrix] = []
    for index in order:
        current = raw["walls"][index]
        if strict:
            gram = sp.Matrix([[dot(raw["gram"], left, right) for right in strict] for left in strict])
            vector = sp.Matrix([dot(raw["gram"], previous, current) for previous in strict])
            projection = gram.inv() * vector
            total = current - sum(
                (projection[position] * strict[position] for position in range(len(strict))),
                sp.zeros(12, 1),
            )
        else:
            total = current
        if any(sp.denom(entry) != 1 for entry in total):
            raise ArithmeticError("a total exceptional transform is not integral")
        strict.append(current)
        totals.append(total)
    basis = sp.Matrix.hstack(plane_class(raw, target), *totals)
    if basis.T * raw["gram"] * basis != sp.diag(1, *([-1] * 11)) or abs(int(basis.det())) != 1:
        raise ArithmeticError("the independently derived plane basis is not unimodular I_(1,11)")
    return order, basis


def proximity_chains(root_coordinates: Sequence[tuple[int, ...]]) -> list[list[int]]:
    parent: dict[int, int] = {}
    for coordinates in root_coordinates:
        if coordinates[0] != 0:
            continue
        negative = [index for index, value in enumerate(coordinates[1:]) if value == -1]
        positive = [index for index, value in enumerate(coordinates[1:]) if value == 1]
        if len(negative) != 1 or len(positive) != 1:
            raise ArithmeticError("a degree-zero root is not a proximity difference")
        parent[negative[0]] = positive[0]
    roots = sorted(set(range(11)) - set(parent))
    chains: list[list[int]] = []
    for root in roots:
        chain = [root]
        while True:
            children = [child for child, old in parent.items() if old == chain[-1]]
            if not children:
                break
            if len(children) != 1:
                raise ArithmeticError("the proximity forest branches")
            chain.append(children[0])
        chains.append(chain)
    return sorted(chains, key=lambda value: (len(value), value))


def plane_vector(degree: int, multiplicities: dict[int, int]) -> tuple[int, ...]:
    result = [degree, *([0] * 11)]
    for position, multiplicity in multiplicities.items():
        result[1 + position] = -multiplicity
    return tuple(result)


def derive_program(raw: dict[str, Any], target: Sequence[int], coordinates: list[tuple[int, ...]]) -> dict[str, Any]:
    roots = [(index, coordinates[index]) for index in raw["roots"]]
    chains = proximity_chains([value for _index, value in roots])
    line_index = next(index for index in raw["boundaries"] if coordinates[index][0] == 1)
    quintic_index = next(index for index in raw["boundaries"] if coordinates[index][0] == 5)
    line_support = {position for position, value in enumerate(coordinates[line_index][1:]) if value == -1}
    line_chains = [chain for chain in chains if set(chain).issubset(line_support)]
    if sorted(map(len, line_chains)) != [2, 3]:
        raise ArithmeticError("the line boundary does not select 3+2 chains")
    chain_a = next(chain for chain in line_chains if len(chain) == 3)
    chain_b = next(chain for chain in line_chains if len(chain) == 2)
    singular = [chain for chain in chains if chain not in line_chains]
    if len(singular) != 2 or any(len(chain) != 3 for chain in singular):
        raise ArithmeticError("the remaining chains are not two triples")

    positive_roots = []
    for index, value in roots:
        if value[0] == 1:
            support = frozenset(position for position, coefficient in enumerate(value[1:]) if coefficient == -1)
            positive_roots.append((index, value, support))
    proper_singular = {chain[0] for chain in singular}
    bcd = [record for record in positive_roots if record[2] == frozenset({chain_b[0], *proper_singular})]
    if len(bcd) != 1:
        raise ArithmeticError("the BCD line is not intrinsic")
    tangent: dict[tuple[int, ...], tuple[int, tuple[int, ...], frozenset[int]]] = {}
    for chain in singular:
        support = frozenset({chain_a[0], chain[0], chain[1]})
        matches = [record for record in positive_roots if record != bcd[0] and record[2] == support]
        if len(matches) != 1:
            raise ArithmeticError("a singular chain has no unique tangent root line")
        tangent[tuple(chain)] = matches[0]
    ordered = sorted((tuple(chain) for chain in singular), key=lambda chain: (tuple(sorted(tangent[chain][2])), tangent[chain][0]))
    chain_c, chain_d = map(list, ordered)
    labeled = {
        **{f"A{i}": value for i, value in enumerate(chain_a)},
        **{f"B{i}": value for i, value in enumerate(chain_b)},
        **{f"C{i}": value for i, value in enumerate(chain_c)},
        **{f"D{i}": value for i, value in enumerate(chain_d)},
    }
    ab = {labeled[name] for name in ("A0", "A1", "A2", "B0", "B1")}
    cd = {labeled[name] for name in ("C0", "C1", "C2", "D0", "D1", "D2")}
    expected_boundaries = {
        plane_vector(1, {position: 1 for position in ab}),
        plane_vector(5, {**{position: 1 for position in ab}, **{position: 2 for position in cd}}),
    }
    if {coordinates[index] for index in raw["boundaries"]} != expected_boundaries:
        raise ArithmeticError("the derived boundary classes do not equal the raw boundaries")

    expected_roots: set[tuple[int, ...]] = set()
    for chain in (chain_a, chain_b, chain_c, chain_d):
        for edge in range(len(chain) - 1):
            vector = [0, *([0] * 11)]
            vector[1 + chain[edge]] = 1
            vector[1 + chain[edge + 1]] = -1
            expected_roots.add(tuple(vector))
    expected_roots |= {
        plane_vector(1, {labeled["A0"]: 1, labeled["C0"]: 1, labeled["C1"]: 1}),
        plane_vector(1, {labeled["A0"]: 1, labeled["D0"]: 1, labeled["D1"]: 1}),
        plane_vector(1, {labeled["B0"]: 1, labeled["C0"]: 1, labeled["D0"]: 1}),
    }
    if {coordinates[index] for index in raw["roots"]} != expected_roots or len(expected_roots) != 10:
        raise ArithmeticError("the seven proximity and three line roots do not equal the raw roots")

    contracted_minus_one = [index for index in target if raw["norms"][index] == -1]
    expected_deep = {
        plane_vector(0, {labeled["A2"]: -1}),
        plane_vector(0, {labeled["B1"]: -1}),
        plane_vector(0, {labeled["C2"]: -1}),
        plane_vector(0, {labeled["D2"]: -1}),
    }
    if {coordinates[index] for index in contracted_minus_one} != expected_deep:
        raise ArithmeticError("the four deep exceptionals do not equal the raw target walls")

    order = ("A0", "A1", "A2", "B0", "B1", "C0", "C1", "C2", "D0", "D1", "D2")
    canonical = [
        tuple([vector[0], *(vector[1 + labeled[name]] for name in order)])
        for vector in coordinates
    ]
    return {
        "chains": {"A": chain_a, "B": chain_b, "C": chain_c, "D": chain_d},
        "label_positions": labeled,
        "canonical_coordinates": canonical,
        "canonical_wall_set_sha256": hashlib.sha256(json.dumps(sorted(canonical), separators=(",", ":")).encode()).hexdigest(),
    }


def raw_endpoint_audit(raw: dict[str, Any]) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    endpoints = sequential_endpoints(raw)
    signatures: Counter[str] = Counter()
    selected: list[tuple[int, ...]] = []
    for contracted, _inverse in endpoints.values():
        signature = endpoint_signature(raw, contracted)
        signatures[json.dumps(signature)] += 1
        if signature == TARGET_SIGNATURE:
            selected.append(contracted)
    if len(endpoints) != 52 or len(signatures) != 9 or len(selected) != 4:
        raise ArithmeticError("the independent endpoint denominator is not 52/9/4")
    branches: list[dict[str, Any]] = []
    for target in sorted(selected):
        order, basis = plane_basis(raw, target)
        inverse = basis.inv()
        coordinates = []
        for wall in raw["walls"]:
            value = inverse * wall
            if any(sp.denom(entry) != 1 for entry in value):
                raise ArithmeticError("a wall left a determinant-one plane basis")
            coordinates.append(tuple(int(entry) for entry in value))
        if inverse * raw["canonical"] != sp.Matrix([-3, *([1] * 11)]):
            raise ArithmeticError("K is not -3H+sum(Ei)")
        reconstruction = basis * sp.Matrix.hstack(*(sp.Matrix(value) for value in coordinates))
        if reconstruction != sp.Matrix.hstack(*raw["walls"]):
            raise ArithmeticError("the complete raw wall transport failed")
        program = derive_program(raw, target, coordinates)
        branches.append(
            {
                "target_raw_wall_indices_zero_based": list(target),
                "contraction_order_zero_based": order,
                "plane_basis_determinant": int(basis.det()),
                "chains_proper_to_deep": program["chains"],
                "canonical_label_positions": program["label_positions"],
                "canonical_wall_set_sha256": program["canonical_wall_set_sha256"],
                "canonical_coordinates": program["canonical_coordinates"],
                "all_26_classes_reconstruct_raw_vectors": True,
            }
        )
    hashes = {branch["canonical_wall_set_sha256"] for branch in branches}
    if len(hashes) != 1:
        raise ArithmeticError("the four branches do not have one canonical 26-wall set")
    return {
        "all_sequential_P2_endpoint_count": len(endpoints),
        "endpoint_degree_signature_count": len(signatures),
        "selected_line_quintic_endpoint_count": len(selected),
        "signature_histogram": dict(sorted(signatures.items())),
        "all_four_plane_bases_unimodular": all(abs(branch["plane_basis_determinant"]) == 1 for branch in branches),
        "all_four_have_one_canonical_26_wall_set": len(hashes) == 1,
    }, branches


def low_terms(
    polynomial: sp.Expr, first: sp.Symbol, second: sp.Symbol, multiplicity: int
) -> tuple[list[sp.Expr], sp.Expr]:
    parsed = sp.Poly(sp.expand(polynomial), first, second)
    equations: list[sp.Expr] = []
    high = sp.Integer(0)
    for (left, right), coefficient in parsed.terms():
        term = coefficient * first**left * second**right
        if left + right < multiplicity:
            equations.append(coefficient)
        else:
            high += term
    return equations, sp.expand(high)


def chain_conditions(
    polynomial: sp.Expr,
    x_map: Callable[[sp.Symbol, sp.Symbol], sp.Expr],
    y_map: Callable[[sp.Symbol, sp.Symbol], sp.Expr],
    multiplicities: Sequence[int],
    later_slopes: Sequence[int],
) -> list[sp.Expr]:
    u, v = sp.symbols("u v")
    current = sp.expand(polynomial.subs({X: x_map(u, v), Y: y_map(u, v)}))
    transverse = v
    answer: list[sp.Expr] = []
    for level, multiplicity in enumerate(multiplicities):
        equations, high = low_terms(current, u, transverse, multiplicity)
        answer.extend(equations)
        if level + 1 == len(multiplicities):
            break
        next_transverse = sp.Symbol(f"raw_audit_w{level}")
        slope = 0 if level == 0 else later_slopes[level - 1]
        current = sp.cancel(
            sp.expand(high.subs(transverse, u * (next_transverse + slope))) / u**multiplicity
        )
        transverse = next_transverse
    return [sp.expand(entry) for entry in answer]


CHAIN_CHARTS: dict[
    str,
    tuple[
        Callable[[sp.Symbol, sp.Symbol], sp.Expr],
        Callable[[sp.Symbol, sp.Symbol], sp.Expr],
        tuple[int, ...],
    ],
] = {
    "A": (lambda u, v: u, lambda u, v: v, (0,)),
    "B": (lambda u, v: 1 + u, lambda u, v: v, tuple()),
    "C": (lambda u, v: v, lambda u, v: 1 + u, (1,)),
    "D": (lambda u, v: -1 + u, lambda u, v: 2 - 2 * u + v, (1,)),
}


def actual_chain_multiplicities(
    polynomial: sp.Expr,
    x_map: Callable[[sp.Symbol, sp.Symbol], sp.Expr],
    y_map: Callable[[sp.Symbol, sp.Symbol], sp.Expr],
    length: int,
    later_slopes: Sequence[int],
) -> tuple[int, ...]:
    u, v = sp.symbols("u v")
    current = sp.expand(polynomial.subs({X: x_map(u, v), Y: y_map(u, v)}))
    transverse = v
    answer: list[int] = []
    for level in range(length):
        parsed = sp.Poly(sp.expand(current), u, transverse)
        order = min(sum(exponents) for exponents, coefficient in parsed.terms() if coefficient)
        answer.append(order)
        if level + 1 == length:
            break
        next_transverse = sp.Symbol(f"raw_audit_m{level}")
        slope = 0 if level == 0 else later_slopes[level - 1]
        current = sp.cancel(
            sp.expand(current.subs(transverse, u * (next_transverse + slope))) / u**order
        )
        transverse = next_transverse
    return tuple(answer)


def homogeneous(polynomial: sp.Expr, degree: int) -> sp.Expr:
    return sp.expand(polynomial.subs({X: X / Z, Y: Y / Z}) * Z**degree)


def poly_string(polynomial: sp.Expr) -> str:
    return str(sp.expand(polynomial)).replace("**", "^")


def singular_absolute_audit(polynomials: list[tuple[int, sp.Expr]]) -> dict[int, dict[str, int | bool]]:
    executable = shutil.which("Singular")
    if executable is None:
        raise RuntimeError("Singular is required for the hostile absolute-factor audit")
    program = 'ring raw_audit=0,(x,y,z),dp; LIB "absfact.lib";\n'
    for index, polynomial in polynomials:
        program += f"poly f{index}={poly_string(polynomial)};\n"
        program += f"def ar{index}=absFactorize(f{index}); setring ar{index};\n"
        program += f'print("RAW_AUDIT_{index}"); print(size(absolute_factors[1])); print(absolute_factors[2][2]); print(absolute_factors[4]);\n'
        program += "setring raw_audit;\n"
    program += "exit;\n"
    completed = subprocess.run(
        [executable, "-q", "--no-rc"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=60,
    )
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    answer: dict[int, dict[str, int | bool]] = {}
    for index, _polynomial in polynomials:
        marker = f"RAW_AUDIT_{index}"
        position = lines.index(marker)
        count, exponent, components = map(int, lines[position + 1 : position + 4])
        answer[index] = {
            "absolute_factor_count_with_unit": count,
            "absolute_nonunit_exponent": exponent,
            "absolute_component_count": components,
            "absolutely_irreducible": (count, exponent, components) == (2, 1, 1),
        }
    return answer


def missing_curve_audit(raw: dict[str, Any], branch: dict[str, Any]) -> dict[str, Any]:
    coordinates = [tuple(value) for value in branch["canonical_coordinates"]]
    positive_cusps = [index for index in raw["minus_one"] if coordinates[index][0] > 0]
    deep = [index for index in raw["minus_one"] if coordinates[index][0] == 0]
    if len(positive_cusps) != 10 or len(deep) != 4:
        raise ArithmeticError("the 14 (-1)-walls did not split as ten plane curves plus four deep exceptionals")
    classes: list[dict[str, Any]] = []
    singular_input: list[tuple[int, sp.Expr]] = []
    for index in positive_cusps:
        vector = coordinates[index]
        degree = vector[0]
        multiplicities = {
            "A": tuple(-value for value in vector[1:4]),
            "B": tuple(-value for value in vector[4:6]),
            "C": tuple(-value for value in vector[6:9]),
            "D": tuple(-value for value in vector[9:12]),
        }
        if any(value < 0 for values in multiplicities.values() for value in values):
            raise ArithmeticError("a positive-degree cusp class has a positive exceptional coefficient")
        monomials = tuple(X**i * Y**j for i in range(degree + 1) for j in range(degree + 1 - i))
        coefficients = sp.symbols(f"raw_c{index}_0:{len(monomials)}")
        generic = sum(coefficient * monomial for coefficient, monomial in zip(coefficients, monomials))
        conditions: list[sp.Expr] = []
        for name in ("A", "B", "C", "D"):
            expected = multiplicities[name]
            while expected and expected[-1] == 0:
                expected = expected[:-1]
            if not expected:
                continue
            x_map, y_map, slopes = CHAIN_CHARTS[name]
            conditions.extend(chain_conditions(generic, x_map, y_map, expected, slopes))
        matrix = sp.Matrix(
            [[condition.coeff(coefficient) for coefficient in coefficients] for condition in conditions]
        )
        kernel = matrix.nullspace()
        if matrix.rank() != len(monomials) - 1 or len(kernel) != 1:
            raise ArithmeticError(f"raw cusp wall {index} does not have a unique plane curve")
        values = primitive(kernel[0], orient_first=True)
        polynomial = sp.expand(sum(value * monomial for value, monomial in zip(values, monomials)))
        actual = {}
        for name in ("A", "B", "C", "D"):
            x_map, y_map, slopes = CHAIN_CHARTS[name]
            actual[name] = actual_chain_multiplicities(polynomial, x_map, y_map, len(multiplicities[name]), slopes)
        if actual != multiplicities:
            raise ArithmeticError(f"raw cusp wall {index} acquired wrong exact multiplicities: {actual}")
        gram = sp.diag(1, *([-1] * 11))
        canonical = sp.Matrix([-3, *([1] * 11)])
        class_vector = sp.Matrix(vector)
        square = int(dot(gram, class_vector, class_vector))
        k_dot = int(dot(gram, canonical, class_vector))
        arithmetic_genus = 1 + (square + k_dot) // 2
        if (square, k_dot, arithmetic_genus) != (-1, -1, 0):
            raise ArithmeticError("a missing curve class is not an adjunction (-1)-class")
        projective = homogeneous(polynomial, degree)
        singular_input.append((index, projective))
        classes.append(
            {
                "raw_index_zero_based_on_first_branch": index,
                "degree": degree,
                "multiplicities_A_B_C_D": {name: list(values) for name, values in multiplicities.items()},
                "condition_count": len(conditions),
                "condition_rank": int(matrix.rank()),
                "linear_system_dimension_projective": 0,
                "affine_equation": str(polynomial),
                "homogeneous_equation": str(projective),
                "exact_cluster_multiplicities": {name: list(values) for name, values in actual.items()},
                "square": square,
                "K_dot": k_dot,
                "arithmetic_genus_after_blowup": arithmetic_genus,
            }
        )
    absolute = singular_absolute_audit(singular_input)
    for record in classes:
        index = record["raw_index_zero_based_on_first_branch"]
        record["absolute_factor_audit"] = absolute[index]
        record["strict_transform_smooth_geometrically_integral_rational"] = (
            absolute[index]["absolutely_irreducible"] and record["arithmetic_genus_after_blowup"] == 0
        )
    if not all(record["strict_transform_smooth_geometrically_integral_rational"] for record in classes):
        raise ArithmeticError("a missing cusp lift is not a smooth integral rational strict transform")
    return {
        "four_deep_exceptional_minus_one_curves": deep,
        "ten_positive_degree_minus_one_curve_count": len(classes),
        "degree_multiset": sorted(record["degree"] for record in classes),
        "all_ten_unique_plane_curves": True,
        "all_ten_absolutely_irreducible": True,
        "all_ten_exact_strict_transforms_smooth_rational": True,
        "curves": classes,
    }


def cdd_rays(inequalities: sp.Matrix) -> list[sp.Matrix]:
    executable = shutil.which("cddexec_gmp")
    if executable is None:
        raise RuntimeError("cddexec_gmp is required for the hostile cone replay")
    lines = [
        "H-representation",
        "begin",
        f"{inequalities.rows} {inequalities.cols + 1} rational",
        *("0 " + " ".join(str(value) for value in inequalities.row(index)) for index in range(inequalities.rows)),
        "end",
    ]
    completed = subprocess.run(
        [executable, "--rep"],
        input="\n".join(lines) + "\n",
        text=True,
        capture_output=True,
        check=True,
        timeout=60,
    )
    output = completed.stdout.splitlines()
    begin = output.index("begin")
    count, width, kind = output[begin + 1].split()
    if kind != "rational" or int(width) != inequalities.cols + 1:
        raise ArithmeticError("cdd returned an unexpected generator format")
    rays: list[sp.Matrix] = []
    for line in output[begin + 2 : begin + 2 + int(count)]:
        tokens = [Fraction(value) for value in line.split()]
        if tokens[0] != 0:
            raise ArithmeticError("a homogeneous cone acquired an affine vertex row")
        denominator = math.lcm(*(value.denominator for value in tokens[1:]))
        entries = [value.numerator * (denominator // value.denominator) for value in tokens[1:]]
        divisor = math.gcd(*(abs(value) for value in entries))
        rays.append(sp.Matrix([value // divisor for value in entries]))
    return rays


def rational_lattice_basis(generators: sp.Matrix) -> sp.Matrix:
    denominator = 1
    for entry in generators:
        denominator = int(sp.ilcm(denominator, sp.denom(entry)))
    integral = (denominator * generators).applyfunc(int)
    return hermite_normal_form(integral) / denominator


def cm_root_and_cusp_audit(raw: dict[str, Any]) -> dict[str, Any]:
    boundaries = [raw["walls"][index] for index in raw["boundaries"]]
    crossing = [
        index
        for index in raw["minus_one"]
        if sorted(int(raw["wall_gram"][index, boundary]) for boundary in raw["boundaries"])
        == [1, 1]
    ]
    cusp_wall_indices = [
        index
        for index in raw["minus_one"]
        if sorted(int(raw["wall_gram"][index, boundary]) for boundary in raw["boundaries"])
        == [0, 2]
    ]
    if crossing != [12, 13] or cusp_wall_indices != list(range(14, 26)):
        raise ArithmeticError("the numerical (-1)-wall taxonomy changed")
    roots = [raw["walls"][index] for index in raw["roots"]]
    roots.extend(2 * raw["walls"][index] + boundaries[0] / 2 + boundaries[1] / 2 for index in crossing)
    if len(roots) != 12 or any(dot(raw["gram"], root, root) != -2 for root in roots):
        raise ArithmeticError("the reconstructed CM roots changed")
    if any(dot(raw["gram"], root, boundary) != 0 for root in roots for boundary in boundaries):
        raise ArithmeticError("a CM root is not boundary-orthogonal")

    independent: list[sp.Matrix] = []
    for root in roots:
        if sp.Matrix.hstack(*independent, root).rank() > len(independent):
            independent.append(root)
        if len(independent) == 10:
            break
    ambient_basis = sp.Matrix.hstack(*independent)
    pivot_rows = (ambient_basis.T).rref()[1]
    square = ambient_basis.extract(pivot_rows, range(10))
    root_coordinates = []
    for root in roots:
        coordinates = square.inv() * root.extract(pivot_rows, [0])
        if ambient_basis * coordinates != root:
            raise ArithmeticError("a root coordinate solve failed")
        root_coordinates.append(coordinates)
    root_lattice_coordinates = rational_lattice_basis(sp.Matrix.hstack(*root_coordinates))
    root_lattice_basis = ambient_basis * root_lattice_coordinates
    root_lattice_gram = root_lattice_basis.T * raw["gram"] * root_lattice_basis
    if root_lattice_gram.det() != -4:
        raise ArithmeticError("the twelve roots no longer span an index-two CM sublattice")

    glues: list[tuple[int, ...]] = []
    for bits in itertools.product((0, 1), repeat=10):
        if not any(bits):
            continue
        vector = sp.Matrix(bits)
        pairings = root_lattice_gram * vector / 2
        square_value = (vector.T * root_lattice_gram * vector)[0] / 4
        if (
            all(sp.denom(entry) == 1 for entry in pairings)
            and sp.denom(square_value) == 1
            and int(square_value) % 2 == 0
        ):
            glues.append(bits)
    if len(glues) != 1:
        raise ArithmeticError("the even index-two CM overlattice is not unique")
    glue = root_lattice_coordinates * sp.Matrix(glues[0]) / 2
    cm_coordinates = rational_lattice_basis(sp.Matrix.hstack(root_lattice_coordinates, glue))
    cm_basis = ambient_basis * cm_coordinates
    cm_gram = cm_basis.T * raw["gram"] * cm_basis
    if cm_gram.det() != -1 or any(int(cm_gram[index, index]) % 2 for index in range(10)):
        raise ArithmeticError("the reconstructed CM lattice is not even unimodular E10")

    inequalities = sp.Matrix([list(root.T * raw["gram"] * cm_basis) for root in roots])
    rays = cdd_rays(inequalities)
    counts = Counter()
    finite_squares: Counter[int] = Counter()
    ideal: list[sp.Matrix] = []
    for ray in rays:
        square_value = int((ray.T * cm_gram * ray)[0])
        kind = "finite" if square_value > 0 else "ideal" if square_value == 0 else "hyperideal"
        counts[kind] += 1
        if kind == "finite":
            finite_squares[square_value] += 1
        elif kind == "ideal":
            ideal.append(ray)
    face_ranks = []
    for facet in range(12):
        on_facet = [ray for ray in rays if (inequalities * ray)[facet] == 0]
        face_ranks.append(sp.Matrix.hstack(*on_facet).rank())
    if (len(rays), counts["finite"], counts["ideal"], counts["hyperideal"]) != (23, 14, 9, 0):
        raise ArithmeticError("the independent CM cone is not 14/9/0")
    if face_ranks != [9] * 12:
        raise ArithmeticError("a CM root facet is redundant")

    regenerated: list[sp.Matrix] = []
    ideal_integrality: list[list[bool]] = []
    cusp_records: list[dict[str, Any]] = []
    for ray_index, ray in enumerate(ideal):
        f_vector = cm_basis * ray
        integrality: list[bool] = []
        lifts: list[list[int]] = []
        for boundary in boundaries:
            candidate = f_vector - boundary / 2
            integral = all(sp.denom(entry) == 1 for entry in candidate)
            integrality.append(integral)
            if integral:
                if dot(raw["gram"], candidate, candidate) != -1:
                    raise ArithmeticError("an integral CM cusp lift is not square -1")
                regenerated.append(candidate)
                lifts.append([int(entry) for entry in candidate])
        ideal_integrality.append(integrality)
        cusp_records.append(
            {
                "ideal_ray_index_zero_based": ray_index,
                "primitive_CM_coordinates": [int(entry) for entry in ray],
                "boundary_lift_integrality": integrality,
                "integral_lifts_in_raw_NS": lifts,
            }
        )
    if Counter(tuple(row) for row in ideal_integrality) != Counter({(True, True): 6, (False, False): 3}):
        raise ArithmeticError("the nine CM cusps do not split 6+3")
    raw_cusps = {tuple(int(entry) for entry in raw["walls"][index]) for index in cusp_wall_indices}
    regenerated_set = {tuple(int(entry) for entry in vector) for vector in regenerated}
    if regenerated_set != raw_cusps or len(regenerated) != 12:
        raise ArithmeticError("complete primitive CM-cusp regeneration is not the literal raw wall set")
    return {
        "CM_root_count": len(roots),
        "CM_root_rank": sp.Matrix.hstack(*roots).rank(),
        "root_sublattice_determinant": int(root_lattice_gram.det()),
        "unique_even_index_two_CM_overlattice": True,
        "CM_lattice_determinant": int(cm_gram.det()),
        "CM_extreme_rays": len(rays),
        "CM_finite_ideal_hyperideal": [counts["finite"], counts["ideal"], counts["hyperideal"]],
        "CM_finite_square_histogram": {str(key): value for key, value in sorted(finite_squares.items())},
        "all_twelve_CM_facets_irredundant": True,
        "split_and_together_cusps": [6, 3],
        "all_and_only_integral_cusp_lifts_equal_twelve_raw_walls": True,
        "cusp_records": cusp_records,
    }


def full_cone_audit(raw: dict[str, Any]) -> dict[str, Any]:
    wall_matrix = sp.Matrix.hstack(*raw["walls"])
    pivots: list[int] = []
    for index in range(26):
        if wall_matrix[:, [*pivots, index]].rank() > len(pivots):
            pivots.append(index)
        if len(pivots) == 12:
            break
    basis = wall_matrix[:, pivots]
    basis_gram = basis.T * raw["gram"] * basis
    inequalities = wall_matrix.T * raw["gram"] * basis
    rays = cdd_rays(inequalities)
    counts = Counter()
    for ray in rays:
        square = (ray.T * basis_gram * ray)[0]
        counts["finite" if square > 0 else "ideal" if square == 0 else "hyperideal"] += 1
    face_ranks = []
    for facet in range(26):
        on_facet = [ray for ray in rays if (inequalities * ray)[facet] == 0]
        face_ranks.append(sp.Matrix.hstack(*on_facet).rank())
    if (len(rays), counts["finite"], counts["ideal"], counts["hyperideal"]) != (550, 508, 42, 0):
        raise ArithmeticError("the fresh full-cone replay is not 508/42/0")
    if face_ranks != [11] * 26:
        raise ArithmeticError("a raw wall is not an irredundant facet")
    off_diagonal = [
        int(raw["wall_gram"][left, right])
        for left in range(26)
        for right in range(left + 1, 26)
    ]
    if min(off_diagonal) < 0:
        raise ArithmeticError("the raw wall system is not acute")
    return {
        "coordinate_basis_facets_zero_based": pivots,
        "extreme_ray_count": len(rays),
        "finite_ideal_hyperideal": [counts["finite"], counts["ideal"], counts["hyperideal"]],
        "all_26_facets_irredundant": True,
        "all_facet_face_ranks": face_ranks,
        "off_diagonal_pairing_range": [min(off_diagonal), max(off_diagonal)],
        "finite_volume_noncompact": True,
    }


def build_certificate() -> dict[str, Any]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    raw = normalize_raw(source)
    endpoint, branches = raw_endpoint_audit(raw)
    memory_gate()
    missing = missing_curve_audit(raw, branches[0])
    memory_gate()
    cm = cm_root_and_cusp_audit(raw)
    full = full_cone_audit(raw)
    checks = {
        "raw_endpoint_denominator_52_9_4": [
            endpoint["all_sequential_P2_endpoint_count"],
            endpoint["endpoint_degree_signature_count"],
            endpoint["selected_line_quintic_endpoint_count"],
        ]
        == [52, 9, 4],
        "four_unimodular_coefficientwise_transports": endpoint["all_four_plane_bases_unimodular"]
        and endpoint["all_four_have_one_canonical_26_wall_set"],
        "ten_missing_lifts_explicitly_constructed": missing["ten_positive_degree_minus_one_curve_count"] == 10
        and missing["all_ten_unique_plane_curves"],
        "ten_missing_lifts_geometrically_valid": missing["all_ten_absolutely_irreducible"]
        and missing["all_ten_exact_strict_transforms_smooth_rational"],
        "finite_index_CM_cone_replayed": cm["CM_finite_ideal_hyperideal"] == [14, 9, 0]
        and cm["all_twelve_CM_facets_irredundant"],
        "complete_cusp_lifts_literal": cm["all_and_only_integral_cusp_lifts_equal_twelve_raw_walls"],
        "full_finite_volume_cone_replayed": full["finite_ideal_hyperideal"] == [508, 42, 0]
        and full["all_26_facets_irredundant"],
    }
    certificate = {
        "schema": SCHEMA,
        "status_date": "2026-08-09",
        "verdict": {
            "original_primary_certificate": "FAIL_AS_WRITTEN: its hidden-wall theorem gate was phrase matching and the root-seed model constructed only 16 of 26 negative curves",
            "repaired_theorem": "PASS: the hostile replay explicitly constructs the ten missing cusp lifts, after which root exhaustion and the complete primitive-cusp gate prove exact same marked P",
        },
        "input_policy": "only the raw odd-unimodular Gram and 26 primitive wall vectors are used in endpoint selection and plane transport; no project Python helper is imported",
        "raw_endpoint_replay": endpoint,
        "branches": branches,
        "explicit_missing_curve_repair": missing,
        "CM_root_and_cusp_replay": cm,
        "full_cone_replay": full,
        "mathematical_closure": {
            "root_exhaustion": "the ten actual (-2)-curves plus the two crossing roots are distinct effective irreducible roots; their fresh 14/9/0 CM chamber has finite covolume, so Kondo Proposition 3.7 excludes every hidden root",
            "same_boundary_hidden_curve_exclusion": "a hidden (-1)-curve meeting one boundary twice gives a primitive nef isotropic CM ray; the complete nine-cusp replay lists its integral class, already represented by one of the twelve constructed cusp curves",
            "negative_curve_trichotomy": "on a strong Coble surface every negative curve is a boundary, a (-2)-curve, or a (-1)-curve",
            "literal_return": "all 26 raw classes are now actual curves and every negative curve has one of those classes, hence the marked nef cone is literally the raw chamber",
            "published_dependency": "Kondo, Coble surfaces with finite automorphism group, Proposition 3.7 and the genus-one/primitive-isotropic results",
        },
        "historical_primary_sha256_before_repair": "af2c6e7fecabdd9c9190cb040c838247a5e51c71c9703840986fa21fbced824b",
        "sources": {
            "raw_input_sha256": sha256(SOURCE),
            "Kondo": "https://arxiv.org/abs/2102.00591",
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "resource_audit": {
            "ceiling_gib": 0.5,
            "observed_outer_watchdog_peak_gib": 0.086,
            "in_process_gate_checked": True,
            "below_ceiling": peak_gib() < 0.5,
        },
    }
    if not certificate["all_checks_pass"] or not certificate["resource_audit"]["below_ceiling"]:
        raise ArithmeticError("the raw hostile repair did not close")
    return certificate


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-output", type=Path)
    arguments = parser.parse_args()
    certificate = build_certificate()
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if arguments.write:
        OUTPUT.write_text(rendered, encoding="utf-8")
    if arguments.check_output is not None:
        if arguments.check_output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("the frozen raw hostile certificate differs from a fresh replay")
    if not arguments.write and arguments.check_output is None:
        print(rendered, end="")


if __name__ == "__main__":
    main()
