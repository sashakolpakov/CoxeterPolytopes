#!/usr/bin/env python3
"""Certify the A3/H3 plane-arrangement Coble factory search.

This file deliberately starts from the two remaining irreducible rank-three
real reflection arrangements.  It reconstructs their projective intersection
lattices, imposes exact anticanonical and bi-anticanonical jet conditions on
proper/first-near clusters, and audits the finite negative-curve carriers by
exact cddlib.

The H3 coefficient field is K=Q(sqrt(5)); no floating-point incidence or rank
decision is used.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
import functools
import hashlib
import itertools
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any, Iterable, Sequence

from flint import fmpz_mat
import networkx as nx
import sympy as sp
from sympy.polys.matrices import DomainMatrix


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

OUTPUT = ROOT / "data/enumerations/a3_h3_plane_arrangement_coble_factories.json"
SCHEMA = "coxeter.a3_h3_plane_arrangement_coble_factories.v1"

X, Y, Z = sp.symbols("x y z")
SQRT5 = sp.sqrt(5)
PHI = (1 + SQRT5) / 2
PHI_INV = (SQRT5 - 1) / 2


def _canon_q5(value: sp.Expr) -> sp.Expr:
    value = sp.sympify(value)
    if value == 0:
        return sp.Integer(0)
    if not value.has(SQRT5):
        return sp.cancel(value)
    return sp.collect(sp.expand(sp.radsimp(sp.cancel(value))), SQRT5)


def _projective(vector: Sequence[sp.Expr]) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    entries = [_canon_q5(value) for value in vector]
    pivot = next((value for value in entries if value != 0), None)
    if pivot is None:
        raise ArithmeticError("zero vector has no projective normalization")
    return tuple(_canon_q5(value / pivot) for value in entries)  # type: ignore[return-value]


def _cross(left: Sequence[sp.Expr], right: Sequence[sp.Expr]) -> tuple[sp.Expr, ...]:
    return tuple(_canon_q5(value) for value in sp.Matrix(left).cross(sp.Matrix(right)))


def _dot(left: Sequence[sp.Expr], right: Sequence[sp.Expr]) -> sp.Expr:
    return _canon_q5(sum(a * b for a, b in zip(left, right)))


def _a3_lines() -> list[tuple[sp.Expr, sp.Expr, sp.Expr]]:
    return [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1),
        (1, -1, 0),
        (0, 1, -1),
        (1, 0, -1),
    ]


def _h3_lines() -> list[tuple[sp.Expr, sp.Expr, sp.Expr]]:
    lines: set[tuple[sp.Expr, sp.Expr, sp.Expr]] = {
        _projective((1, 0, 0)),
        _projective((0, 1, 0)),
        _projective((0, 0, 1)),
    }
    for second in (PHI, -PHI):
        for third in (PHI_INV, -PHI_INV):
            base = (sp.Integer(1), second, third)
            for shift in range(3):
                lines.add(_projective(base[shift:] + base[:shift]))
    if len(lines) != 15:
        raise ArithmeticError(f"H3 mirror construction returned {len(lines)} lines")
    return sorted(lines, key=lambda row: tuple(str(value) for value in row))


@dataclass(frozen=True)
class Arrangement:
    name: str
    lines: tuple[tuple[sp.Expr, sp.Expr, sp.Expr], ...]
    points: tuple[tuple[sp.Expr, sp.Expr, sp.Expr], ...]
    incident_lines: tuple[tuple[int, ...], ...]

    @property
    def multiplicity_histogram(self) -> dict[int, int]:
        return dict(sorted(Counter(map(len, self.incident_lines)).items()))

    @property
    def pair_identity(self) -> int:
        return sum(math.comb(len(lines), 2) for lines in self.incident_lines)


def _arrangement(name: str, lines: Sequence[Sequence[sp.Expr]]) -> Arrangement:
    normalized_lines = tuple(_projective(line) for line in lines)
    points = {
        _projective(_cross(normalized_lines[left], normalized_lines[right]))
        for left in range(len(lines))
        for right in range(left)
    }
    ordered_points = tuple(sorted(points, key=lambda row: tuple(str(value) for value in row)))
    incidence = tuple(
        tuple(index for index, line in enumerate(normalized_lines) if _dot(line, point) == 0)
        for point in ordered_points
    )
    result = Arrangement(name, normalized_lines, ordered_points, incidence)
    if result.pair_identity != math.comb(len(lines), 2):
        raise ArithmeticError(f"{name} intersection lattice fails the pair identity")
    return result


def _incidence_graph(arrangement: Arrangement) -> nx.Graph:
    graph = nx.Graph()
    for line in range(len(arrangement.lines)):
        graph.add_node(("L", line), kind="line")
    for point, incident in enumerate(arrangement.incident_lines):
        graph.add_node(("P", point), kind=f"point-{len(incident)}")
        for line in incident:
            graph.add_edge(("P", point), ("L", line))
    return graph


@functools.lru_cache(maxsize=None)
def _line_automorphisms(arrangement: Arrangement) -> tuple[tuple[int, ...], ...]:
    if arrangement.name == "A3":
        graph = _incidence_graph(arrangement)
        matcher = nx.algorithms.isomorphism.GraphMatcher(
            graph,
            graph,
            node_match=lambda left, right: left == right,
        )
        return tuple(sorted(
            {
                tuple(mapping[("L", index)][1] for index in range(len(arrangement.lines)))
                for mapping in matcher.isomorphisms_iter()
            }
        ))
    lookup = {line: index for index, line in enumerate(arrangement.lines)}
    generators = []
    identity = sp.eye(3)
    for normal in arrangement.lines:
        column = sp.Matrix(normal)
        reflection = identity - sp.Rational(2) * column * column.T / (column.T * column)[0]
        permutation = tuple(
            lookup[_projective(reflection * sp.Matrix(line))]
            for line in arrangement.lines
        )
        generators.append(permutation)
    neutral = tuple(range(len(arrangement.lines)))
    group = {neutral}
    queue = [neutral]
    while queue:
        left = queue.pop()
        for right in generators:
            product = tuple(left[right[index]] for index in range(len(left)))
            if product not in group:
                group.add(product)
                queue.append(product)
    return tuple(sorted(group))


def _six_line_profile(arrangement: Arrangement, subset: Sequence[int]) -> dict[str, Any]:
    selected = set(subset)
    selected_multiplicities = [len(selected.intersection(lines)) for lines in arrangement.incident_lines]
    histogram = Counter(value for value in selected_multiplicities if value >= 2)
    triple_points = tuple(index for index, value in enumerate(selected_multiplicities) if value == 3)
    maximum = max(selected_multiplicities)
    pair_count = sum(math.comb(value, 2) for value in selected_multiplicities)
    if pair_count != 15:
        raise ArithmeticError("six selected lines lost a pairwise intersection")
    ordinary_boundary = maximum <= 3
    t = len(triple_points)
    return {
        "subset": tuple(subset),
        "intersection_multiplicity_histogram": dict(sorted(histogram.items())),
        "maximum_selected_multiplicity": maximum,
        "triple_points": triple_points,
        "triple_count": t,
        "ordinary_reduced_boundary_admissible": ordinary_boundary and t <= 4,
        "boundary_component_count": 6 + t,
        "blowup_count": 15 + t,
        "hyperbolic_dimension": 15 + t,
    }


def _six_line_orbits(arrangement: Arrangement) -> list[dict[str, Any]]:
    automorphisms = _line_automorphisms(arrangement)
    unseen = set(itertools.combinations(range(len(arrangement.lines)), 6))
    records = []
    while unseen:
        representative = min(unseen)
        orbit = {
            tuple(sorted(permutation[index] for index in representative))
            for permutation in automorphisms
        }
        unseen.difference_update(orbit)
        profile = _six_line_profile(arrangement, representative)
        profile["orbit_size"] = len(orbit)
        profile["orbit_profile_consistent"] = len(
            {
                (
                    _six_line_profile(arrangement, member)["maximum_selected_multiplicity"],
                    _six_line_profile(arrangement, member)["triple_count"],
                )
                for member in orbit
            }
        ) == 1
        records.append(profile)
    return sorted(records, key=lambda row: (row["triple_count"], row["subset"]))


def _monomial_exponents(degree: int) -> list[tuple[int, int, int]]:
    return [
        (left, middle, degree - left - middle)
        for left in range(degree, -1, -1)
        for middle in range(degree - left, -1, -1)
    ]


def _derivative_row(
    degree: int,
    point: Sequence[sp.Expr],
    directions: Sequence[Sequence[sp.Expr]],
) -> list[sp.Expr]:
    variables = (X, Y, Z)
    row = []
    for exponents in _monomial_exponents(degree):
        polynomial = sp.prod(variable**exponent for variable, exponent in zip(variables, exponents))
        for direction in directions:
            polynomial = sum(
                direction[index] * sp.diff(polynomial, variables[index])
                for index in range(3)
            )
        row.append(_canon_q5(polynomial.subs(dict(zip(variables, point)))))
    return row


def _tangent_frame(
    point: Sequence[sp.Expr], line: Sequence[sp.Expr]
) -> tuple[tuple[sp.Expr, ...], tuple[sp.Expr, ...]]:
    along = _cross(line, point)
    if _dot(line, along) != 0 or sp.Matrix.hstack(sp.Matrix(point), sp.Matrix(along)).rank() != 2:
        raise ArithmeticError("failed to construct an incident tangent direction")
    transverse = None
    for trial in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
        if sp.Matrix.hstack(sp.Matrix(point), sp.Matrix(along), sp.Matrix(trial)).det() != 0:
            transverse = trial
            break
    if transverse is None:
        raise ArithmeticError("failed to complete a tangent frame")
    return tuple(along), tuple(map(sp.Integer, transverse))


def _cluster_rows(
    arrangement: Arrangement,
) -> tuple[dict[int, list[list[sp.Expr]]], dict[tuple[int, int], list[list[sp.Expr]]]]:
    proper: dict[int, list[list[sp.Expr]]] = {}
    near: dict[tuple[int, int], list[list[sp.Expr]]] = {}
    for point_index, point in enumerate(arrangement.points):
        # Cubic evaluation followed by the three sextic double-point rows.
        frame_line = arrangement.lines[arrangement.incident_lines[point_index][0]]
        first, second = _tangent_frame(point, frame_line)
        proper[point_index] = [
            _derivative_row(3, point, ()),
            _derivative_row(6, point, ()),
            _derivative_row(6, point, (first,)),
            _derivative_row(6, point, (second,)),
        ]
        for line_index in arrangement.incident_lines[point_index]:
            along, transverse = _tangent_frame(point, arrangement.lines[line_index])
            near[(point_index, line_index)] = [
                _derivative_row(3, point, (along,)),
                _derivative_row(6, point, (along, along)),
                _derivative_row(6, point, (along, along, along)),
                _derivative_row(6, point, (along, transverse)),
            ]
    return proper, near


def _integer_rows(rows: Sequence[Sequence[sp.Expr]]) -> list[list[int]]:
    result = []
    for row in rows:
        if any(value.has(SQRT5) for value in row):
            raise ArithmeticError("integer rank path received a Q(sqrt(5)) row")
        denominators = [sp.denom(value) for value in row]
        common = int(sp.ilcm(*denominators)) if denominators else 1
        integers = [int(value * common) for value in row]
        divisor = math.gcd(*(abs(value) for value in integers))
        result.append([value // max(divisor, 1) for value in integers])
    return result


def _fmpz_rank(rows: Sequence[Sequence[int]]) -> int:
    if not rows:
        return 0
    return int(fmpz_mat(len(rows), len(rows[0]), [value for row in rows for value in row]).rank())


def _a3_parent_closed_clusters(
    arrangement: Arrangement,
    minimum_size: int = 10,
    maximum_size: int = 19,
) -> Iterable[tuple[tuple[int, ...], tuple[tuple[int, int], ...]]]:
    options = []
    for point, lines in enumerate(arrangement.incident_lines):
        local = [(False, tuple())]
        local.extend(
            (True, tuple((point, line) for line in chosen))
            for count in range(len(lines) + 1)
            for chosen in itertools.combinations(lines, count)
        )
        options.append(local)
    for choices in itertools.product(*options):
        proper = tuple(index for index, (selected, _) in enumerate(choices) if selected)
        near = tuple(flag for _, flags in choices for flag in flags)
        size = len(proper) + len(near)
        if minimum_size <= size <= maximum_size:
            yield proper, near


def _a3_hilbert_census(arrangement: Arrangement, all_sizes: bool) -> dict[str, Any]:
    proper_rows, near_rows = _cluster_rows(arrangement)
    proper_integer = {
        index: _integer_rows(rows) for index, rows in proper_rows.items()
    }
    near_integer = {
        flag: _integer_rows(rows) for flag, rows in near_rows.items()
    }
    sizes = range(10, 20) if all_sizes else (10,)
    census: dict[int, Counter[tuple[int, int]]] = {size: Counter() for size in sizes}
    survivors: dict[int, list[dict[str, Any]]] = {size: [] for size in sizes}
    for proper, near in _a3_parent_closed_clusters(arrangement, min(sizes), max(sizes)):
        size = len(proper) + len(near)
        if size not in census:
            continue
        cubic = (
            [proper_integer[index][0] for index in proper]
            + [near_integer[flag][0] for flag in near]
        )
        sextic = (
            [row for index in proper for row in proper_integer[index][1:]]
            + [row for flag in near for row in near_integer[flag][1:]]
        )
        cubic_rank = _fmpz_rank(cubic)
        sextic_rank = _fmpz_rank(sextic)
        census[size][(cubic_rank, sextic_rank)] += 1
        if cubic_rank == 10 and sextic_rank < 28:
            survivors[size].append(
                {
                    "proper_points": list(proper),
                    "near_flags": [list(flag) for flag in near],
                    "cubic_rank": cubic_rank,
                    "sextic_rank": sextic_rank,
                    "h0_minus_K": 10 - cubic_rank,
                    "h0_minus_2K": 28 - sextic_rank,
                }
            )
    orbit_summary = _a3_survivor_orbits(arrangement, survivors)
    unique_boundary_audit = _a3_unique_boundary_orbit_audit(
        arrangement,
        orbit_summary,
        proper_integer,
        near_integer,
    )
    return {
        "scope": f"all parent-closed A3 proper/first-near clusters of sizes {min(sizes)}..{max(sizes)}",
        "rank_census_by_size": {
            str(size): {
                f"cubic_{key[0]}_sextic_{key[1]}": value
                for key, value in sorted(census[size].items())
            }
            for size in sizes
        },
        "strong_coble_hilbert_survivor_count_by_size": {
            str(size): len(survivors[size]) for size in sizes
        },
        "survivors": {str(size): survivors[size] for size in sizes if survivors[size]},
        "survivor_orbits": orbit_summary,
        "unique_boundary_orbit_audit": unique_boundary_audit,
    }


def _a3_survivor_orbits(
    arrangement: Arrangement,
    survivors: dict[int, list[dict[str, Any]]],
) -> dict[str, Any]:
    line_group = _line_automorphisms(arrangement)
    point_lookup = {
        frozenset(lines): index for index, lines in enumerate(arrangement.incident_lines)
    }
    actions = []
    for permutation in line_group:
        point_permutation = tuple(
            point_lookup[frozenset(permutation[line] for line in incident)]
            for incident in arrangement.incident_lines
        )
        actions.append((permutation, point_permutation))

    records = {}
    for size, rows in survivors.items():
        row_lookup = {
            (
                tuple(row["proper_points"]),
                tuple(tuple(flag) for flag in row["near_flags"]),
            ): row
            for row in rows
        }
        unseen = set(row_lookup)
        orbits = []
        while unseen:
            representative = min(unseen)
            proper, near = representative
            orbit = set()
            for line_permutation, point_permutation in actions:
                transformed_proper = tuple(sorted(point_permutation[point] for point in proper))
                transformed_near = tuple(
                    sorted(
                        (point_permutation[point], line_permutation[line])
                        for point, line in near
                    )
                )
                candidate = (transformed_proper, transformed_near)
                if candidate in row_lookup:
                    orbit.add(candidate)
            unseen.difference_update(orbit)
            row = row_lookup[representative]
            orbits.append(
                {
                    "proper_points": list(proper),
                    "near_flags": [list(flag) for flag in near],
                    "cubic_rank": row["cubic_rank"],
                    "sextic_rank": row["sextic_rank"],
                    "h0_minus_2K": row["h0_minus_2K"],
                    "orbit_size": len(orbit),
                }
            )
        records[str(size)] = {
            "orbit_count": len(orbits),
            "orbit_count_by_sextic_rank": dict(
                sorted(Counter(row["sextic_rank"] for row in orbits).items())
            ),
            "marked_count_by_sextic_rank": dict(
                sorted(
                    Counter(
                        {
                            rank: sum(
                                row["orbit_size"] for row in orbits if row["sextic_rank"] == rank
                            )
                            for rank in sorted({row["sextic_rank"] for row in orbits})
                        }
                    ).items()
                )
            ),
            "orbits": orbits,
        }
    return records


def _a3_unique_boundary_orbit_audit(
    arrangement: Arrangement,
    orbit_summary: dict[str, Any],
    proper_rows: dict[int, list[list[int]]],
    near_rows: dict[tuple[int, int], list[list[int]]],
) -> dict[str, Any]:
    records = []
    variables = (X, Y, Z)
    monomials = [X**left * Y**middle * Z**right for left, middle, right in _monomial_exponents(6)]
    for size, size_record in orbit_summary.items():
        for orbit in size_record["orbits"]:
            if orbit["sextic_rank"] != 27:
                continue
            proper = tuple(orbit["proper_points"])
            near = tuple(tuple(flag) for flag in orbit["near_flags"])
            rows = (
                [row for index in proper for row in proper_rows[index][1:]]
                + [row for flag in near for row in near_rows[flag][1:]]
            )
            matrix = fmpz_mat(
                len(rows),
                len(rows[0]),
                [value for row in rows for value in row],
            )
            nullspace, nullity = matrix.nullspace()
            if nullity != 1:
                raise ArithmeticError("a rank-27 A3 sextic matrix has nonunit nullity")
            coefficients = _primitive([int(nullspace[index, 0]) for index in range(28)])
            if next(value for value in coefficients if value) < 0:
                coefficients = tuple(-value for value in coefficients)
            polynomial = sp.expand(sum(value * monomial for value, monomial in zip(coefficients, monomials)))
            scalar, factors = sp.factor_list(polynomial, *variables)
            degree_partition = sorted(
                [sp.Poly(factor, *variables).total_degree() for factor, exponent in factors for _ in range(exponent)]
            )
            reduced = all(exponent == 1 for _, exponent in factors)
            irreducible = len(factors) == 1 and factors[0][1] == 1 and degree_partition == [6]

            near_by_point: dict[int, list[int]] = defaultdict(list)
            for point, line in near:
                near_by_point[point].append(line)
            local_records = []
            local_pass = irreducible
            for point_index in proper if irreducible else ():
                point = arrangement.points[point_index]
                flags = sorted(near_by_point[point_index])
                if len(flags) > 1:
                    local_records.append(
                        {
                            "point": point_index,
                            "near_flags": flags,
                            "status": "reject: more than one prescribed first-near double point",
                        }
                    )
                    local_pass = False
                    continue
                line = flags[0] if flags else arrangement.incident_lines[point_index][0]
                along, transverse = _tangent_frame(point, arrangement.lines[line])

                def value(directions: Sequence[Sequence[sp.Expr]]) -> sp.Expr:
                    row = _derivative_row(6, point, directions)
                    return _canon_q5(sum(row[index] * coefficients[index] for index in range(28)))

                vv = value((along, along))
                vn = value((along, transverse))
                nn = value((transverse, transverse))
                if not flags:
                    discriminant = _canon_q5(vn**2 - vv * nn)
                    passed = discriminant != 0
                    local_records.append(
                        {
                            "point": point_index,
                            "kind": "ordinary proper node",
                            "tangent_discriminant": str(discriminant),
                            "passed": passed,
                        }
                    )
                else:
                    vvv = value((along, along, along))
                    vvvv = value((along, along, along, along))
                    vv_n = value((along, along, transverse))
                    # q_near=(Dvvvv/24)t^2+(Dvvn/2)ts+(Dnn/2)s^2.
                    discriminant = _canon_q5(
                        (vv_n / 2) ** 2 - 4 * (vvvv / 24) * (nn / 2)
                    )
                    passed = vv == 0 and vn == 0 and vvv == 0 and nn != 0 and discriminant != 0
                    local_records.append(
                        {
                            "point": point_index,
                            "kind": "proper double tangent plus ordinary first-near node",
                            "near_line": line,
                            "proper_Dvv": str(vv),
                            "proper_Dvn": str(vn),
                            "proper_Dnn": str(nn),
                            "near_Dvvv": str(vvv),
                            "near_tangent_discriminant": str(discriminant),
                            "passed": passed,
                        }
                    )
                local_pass = local_pass and passed
            strong_coble = irreducible and local_pass and int(size) == 10
            records.append(
                {
                    "proper_points": list(proper),
                    "near_flags": [list(flag) for flag in near],
                    "marked_orbit_size": orbit["orbit_size"],
                    "polynomial": str(polynomial),
                    "factor_scalar": str(scalar),
                    "factor_degree_partition_with_multiplicity": degree_partition,
                    "factor_count": len(factors),
                    "reduced": reduced,
                    "irreducible_over_Q": irreducible,
                    "local_resolution_audit": local_records,
                    "local_resolution_passed": local_pass,
                    "nondegenerate_one_boundary_strong_coble": strong_coble,
                }
            )
    return {
        "scope": "all A3 arrangement-group orbits with h0(-K)=0 and unique |-2K|",
        "orbit_count": len(records),
        "factor_partition_orbit_histogram": dict(
            sorted(
                Counter(
                    "+".join(map(str, row["factor_degree_partition_with_multiplicity"]))
                    for row in records
                ).items()
            )
        ),
        "irreducible_orbit_count": sum(row["irreducible_over_Q"] for row in records),
        "local_resolution_pass_orbit_count": sum(row["local_resolution_passed"] for row in records),
        "nondegenerate_one_boundary_strong_coble_orbit_count": sum(
            row["nondegenerate_one_boundary_strong_coble"] for row in records
        ),
        "records": records,
    }


def _intersection(left: Sequence[int], right: Sequence[int]) -> int:
    return left[0] * right[0] - sum(a * b for a, b in zip(left[1:], right[1:]))


def _primitive(values: Sequence[int]) -> tuple[int, ...]:
    divisor = math.gcd(*(abs(value) for value in values))
    if divisor == 0:
        raise ArithmeticError("zero vector is not primitive")
    return tuple(value // divisor for value in values)


def _cdd_rays(inequalities: sp.Matrix, executable: str) -> list[tuple[int, ...]]:
    """Return exact homogeneous rays, accepting cddexec's mislabeled real header.

    Homebrew cddlib 0.94m reports ``real`` in the V-header even though its own
    trailer says ``Number Type = rational`` and every emitted token is exact.
    We parse only Fraction-compatible tokens and therefore retain an exact
    certificate.
    """
    lines = [
        "H-representation",
        "begin",
        f"{inequalities.rows} {inequalities.cols + 1} rational",
    ]
    lines.extend(
        "0 " + " ".join(str(inequalities[i, j]) for j in range(inequalities.cols))
        for i in range(inequalities.rows)
    )
    lines.append("end")
    completed = subprocess.run(
        [executable, "--rep"],
        input="\n".join(lines) + "\n",
        text=True,
        capture_output=True,
        check=False,
        timeout=600,
    )
    if completed.returncode:
        raise RuntimeError(completed.stderr[-4000:] or completed.stdout[-4000:])
    output = completed.stdout.splitlines()
    begin = output.index("begin")
    count, width, kind = output[begin + 1].split()
    if kind not in {"rational", "real"} or int(width) != inequalities.cols + 1:
        raise ArithmeticError("cddlib returned an unexpected generator format")
    result = []
    for line in output[begin + 2 : begin + 2 + int(count)]:
        tokens = [Fraction(token) for token in line.split()]
        if tokens[0] != 0:
            raise ArithmeticError("the homogeneous cone acquired a vertex row")
        denominator = math.lcm(*(value.denominator for value in tokens[1:]))
        integers = [value.numerator * (denominator // value.denominator) for value in tokens[1:]]
        result.append(_primitive(integers))
    if len(result) != len(set(result)):
        raise ArithmeticError("cddlib returned duplicate oriented rays")
    return result


def _cone_audit(
    walls: Sequence[Sequence[int]], cdd_executable: str, *, store_rays: bool = True
) -> dict[str, Any]:
    rank = len(walls[0])
    wall_matrix = sp.Matrix(walls)
    lattice = sp.diag(1, *([-1] * (rank - 1)))
    gram = wall_matrix * lattice * wall_matrix.T
    pivots = list(gram.rref()[1])
    if len(pivots) != rank:
        raise ArithmeticError(f"walls have rank {len(pivots)}, expected {rank}")
    inequalities = gram[:, pivots]
    rays = _cdd_rays(inequalities, cdd_executable)
    basis_gram = gram.extract(pivots, pivots)
    counts: Counter[str] = Counter()
    ray_records = []
    facet_ranks = []
    slacks_by_ray = []
    for ray in rays:
        vector = sp.Matrix(ray)
        slacks = inequalities * vector
        square = int((vector.T * basis_gram * vector)[0])
        kind = "finite" if square > 0 else "ideal" if square == 0 else "hyperideal"
        counts[kind] += 1
        slacks_by_ray.append(slacks)
        ray_records.append({"coordinates": list(ray), "square": square, "kind": kind})
    for facet in range(len(walls)):
        active = [list(rays[index]) for index, slacks in enumerate(slacks_by_ray) if slacks[facet] == 0]
        facet_ranks.append(_fmpz_rank(active) if active else 0)
    sample_hyperideal = []
    basis_walls = wall_matrix.extract(pivots, range(rank))
    for index, row in enumerate(ray_records):
        if row["kind"] != "hyperideal" or len(sample_hyperideal) == 20:
            continue
        active_facets = [
            facet for facet, slack in enumerate(slacks_by_ray[index]) if slack == 0
        ]
        active_rows = [
            [int(inequalities[facet, column]) for column in range(rank)]
            for facet in active_facets
        ]
        physical = basis_walls.T * sp.Matrix(rays[index])
        sample_hyperideal.append(
            {
                **row,
                "physical_lattice_vector": [int(value) for value in physical],
                "active_facet_count": len(active_facets),
                "active_facet_rank": _fmpz_rank(active_rows),
                "all_wall_pairings_nonnegative": all(
                    slack >= 0 for slack in slacks_by_ray[index]
                ),
            }
        )
    payload = json.dumps(sorted(tuple(ray) for ray in rays), separators=(",", ":")).encode()
    return {
        "wall_count": len(walls),
        "wall_rank": wall_matrix.rank(),
        "coordinate_basis_facets_zero_based": pivots,
        "extreme_ray_count": len(rays),
        "finite_vertex_count": counts["finite"],
        "ideal_vertex_count": counts["ideal"],
        "hyperideal_vertex_count": counts["hyperideal"],
        "all_facet_face_ranks": facet_ranks,
        "all_facets_irredundant": all(value == rank - 1 for value in facet_ranks),
        "primitive_ray_list_sha256": hashlib.sha256(payload).hexdigest(),
        "rays": ray_records if store_rays else [],
        "rays_stored": store_rays,
        "sample_hyperideal_rays": sample_hyperideal,
    }


def _h3_triple_orbit(arrangement: Arrangement) -> tuple[int, ...]:
    return tuple(index for index, lines in enumerate(arrangement.incident_lines) if len(lines) == 3)


def _h3_petersen_edges(arrangement: Arrangement, centers: Sequence[int]) -> list[tuple[int, int]]:
    center_position = {point: index for index, point in enumerate(centers)}
    edges = []
    for line in range(len(arrangement.lines)):
        support = [center_position[point] for point in centers if line in arrangement.incident_lines[point]]
        if len(support) != 2:
            raise ArithmeticError("an H3 mirror does not contain two triple-orbit points")
        edges.append(tuple(sorted(support)))
    if len(set(edges)) != 15:
        raise ArithmeticError("H3 triple orbit did not produce fifteen Petersen edges")
    return sorted(edges)


def _rank_qsqrt5(rows: Sequence[Sequence[sp.Expr]]) -> int:
    field = sp.QQ.algebraic_field(SQRT5)
    matrix = DomainMatrix.from_Matrix(sp.Matrix(rows)).convert_to(field)
    return int(matrix.rank())


def _h3_winger_surface(arrangement: Arrangement) -> dict[str, Any]:
    centers = _h3_triple_orbit(arrangement)
    cubic_rows = []
    sextic_rows = []
    for center in centers:
        point = arrangement.points[center]
        along, transverse = _tangent_frame(
            point,
            arrangement.lines[arrangement.incident_lines[center][0]],
        )
        cubic_rows.append(_derivative_row(3, point, ()))
        sextic_rows.extend(
            [
                _derivative_row(6, point, ()),
                _derivative_row(6, point, (along,)),
                _derivative_row(6, point, (transverse,)),
            ]
        )
    cubic_rank = _rank_qsqrt5(cubic_rows)
    sextic_rank = _rank_qsqrt5(sextic_rows)

    polynomial = (
        2 * (X**6 + Y**6 + Z**6)
        + (33 + 27 * SQRT5) * (X**4 * Y**2 + X**2 * Z**4 + Y**4 * Z**2)
        + (33 - 27 * SQRT5) * (X**4 * Z**2 + X**2 * Y**4 + Y**2 * Z**4)
        - 204 * X**2 * Y**2 * Z**2
    )
    coefficient_vector = [
        sp.Poly(polynomial, X, Y, Z).coeff_monomial(X**left * Y**middle * Z**right)
        for left, middle, right in _monomial_exponents(6)
    ]
    if any(
        _canon_q5(sum(row[index] * coefficient_vector[index] for index in range(28))) != 0
        for row in sextic_rows
    ):
        raise ArithmeticError("the displayed Winger sextic misses a prescribed double point")

    collinear_triples = []
    for triple in itertools.combinations(centers, 3):
        determinant = _canon_q5(
            sp.Matrix([arrangement.points[index] for index in triple]).det()
        )
        if determinant == 0:
            collinear_triples.append(triple)

    node_discriminants = []
    for center in centers:
        point = arrangement.points[center]
        along, transverse = _tangent_frame(
            point,
            arrangement.lines[arrangement.incident_lines[center][0]],
        )
        vv = _derivative_row(6, point, (along, along))
        vn = _derivative_row(6, point, (along, transverse))
        nn = _derivative_row(6, point, (transverse, transverse))
        values = [
            _canon_q5(sum(row[index] * coefficient_vector[index] for index in range(28)))
            for row in (vv, vn, nn)
        ]
        # For q=(Dvv/2)u^2+(Dvn)uv+(Dnn/2)v^2, four times the
        # usual discriminant is 4*Dvn^2-Dvv*Dnn.
        discriminant = _canon_q5(4 * values[1] ** 2 - values[0] * values[2])
        node_discriminants.append(str(discriminant))
        if discriminant == 0:
            raise ArithmeticError("a Winger singularity is not an ordinary node")

    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required for the Winger singular-scheme audit")
    singular_program = r'''
LIB "primdec.lib";
ring r=(0,a),(x,y,z),dp;
minpoly=a2-5;
poly f=2*(x6+y6+z6)+(33+27a)*(x4y2+x2z4+y4z2)+(33-27a)*(x4z2+x2y4+y2z4)-204*x2y2z2;
list fac=factorize(f);
ideal j=jacob(f);
ideal rj=radical(j);
print("FACTOR_COUNT="+string(size(fac[1])-1));
print("JAC_DIM="+string(dim(std(j))));
degree(std(j));
degree(std(rj));
print("MINASS="+string(size(minAssGTZ(j))));
'''
    completed = subprocess.run(
        [singular, "-q"],
        input=singular_program,
        text=True,
        capture_output=True,
        check=False,
        timeout=120,
    )
    singular_output = completed.stdout + completed.stderr
    expected_fragments = ["FACTOR_COUNT=1", "JAC_DIM=1", "MINASS=10"]
    if completed.returncode or any(fragment not in singular_output for fragment in expected_fragments):
        raise ArithmeticError("Singular failed the Winger factor/Jacobian audit")
    if singular_output.count("degree (proj.)   = 10") != 2:
        raise ArithmeticError("the Winger Jacobian or its radical has degree other than ten")

    return {
        "field": "Q(sqrt(5)), characteristic zero",
        "center_indices": list(centers),
        "center_coordinates": [
            [str(value) for value in arrangement.points[index]] for index in centers
        ],
        "no_three_centers_collinear": not collinear_triples,
        "collinear_center_triples": [list(row) for row in collinear_triples],
        "cubic_condition_matrix_shape": [10, 10],
        "cubic_rank": cubic_rank,
        "h0_minus_K": 10 - cubic_rank,
        "sextic_double_condition_matrix_shape": [30, 28],
        "sextic_rank": sextic_rank,
        "h0_minus_2K": 28 - sextic_rank,
        "unique_sextic": str(sp.expand(polynomial)),
        "node_tangent_discriminants": node_discriminants,
        "all_ten_prescribed_singularities_ordinary_nodes": all(
            value != "0" for value in node_discriminants
        ),
        "singular_exact_audit": {
            "jacobian_affine_cone_dimension": 1,
            "projective_jacobian_degree": 10,
            "projective_radical_degree": 10,
            "minimal_associated_prime_count": 10,
            "nonunit_factor_count_over_Qsqrt5": 1,
            "raw_output": singular_output.strip().splitlines(),
        },
        "absolute_irreducibility_argument": {
            "K_irreducible": True,
            "arithmetic_genus": 10,
            "node_delta_sum": 10,
            "normalization_genus": 0,
            "equal_degree_geometric_component_counts_to_exclude": [2, 3, 6],
            "exclusion": (
                "K-irreducibility makes absolute components a transitive equal-degree Galois orbit. "
                "For r components the nodal genus formula forces sum(g_i)=r-1: r=2 contradicts "
                "equal conjugate genera, while r=3 conics and r=6 lines have genus sum zero."
            ),
            "absolutely_irreducible": True,
        },
        "surface_conclusion": {
            "blowup_count": 10,
            "picard_rank": 11,
            "hyperbolic_dimension": 10,
            "boundary_class": [6, *([-2] * 10)],
            "boundary_square": -4,
            "boundary_components": 1,
            "boundary_smooth_rational": True,
            "strong_coble": cubic_rank == 10 and sextic_rank == 27,
        },
        "lineage": {
            "classical_name": "Winger sextic / A5 Coble surface",
            "new_surface_claim": False,
            "note": (
                "Dolgachev identifies the unique A5-invariant rational sextic with the ten-node "
                "orbit and states that blowing up the nodes gives a Coble surface; the present "
                "contribution is an exact arrangement/Hilbert/cone replay, not a new surface."
            ),
        },
    }


def _a3_full_six_line_carrier(
    arrangement: Arrangement,
    cdd_executable: str,
) -> dict[str, Any]:
    triple_points = tuple(
        index for index, incident in enumerate(arrangement.incident_lines) if len(incident) == 3
    )
    double_points = tuple(
        index for index, incident in enumerate(arrangement.incident_lines) if len(incident) == 2
    )
    near_flags = tuple(
        (point, line)
        for point in triple_points
        for line in arrangement.incident_lines[point]
    )
    centers: list[tuple[str, int, int | None]] = [
        ("proper", point, None) for point in range(len(arrangement.points))
    ] + [("near", point, line) for point, line in near_flags]
    position = {center: index for index, center in enumerate(centers)}
    rank = 1 + len(centers)

    def divisor(degree: int, coefficients: dict[tuple[str, int, int | None], int]) -> tuple[int, ...]:
        values = [0] * len(centers)
        for center, value in coefficients.items():
            values[position[center]] = value
        return (degree, *values)

    line_boundaries = []
    for line in range(len(arrangement.lines)):
        coefficients: dict[tuple[str, int, int | None], int] = {}
        for point, incident in enumerate(arrangement.incident_lines):
            if line not in incident:
                continue
            coefficients[("proper", point, None)] = -1
            if point in triple_points:
                coefficients[("near", point, line)] = -1
        line_boundaries.append(divisor(1, coefficients))
    exceptional_boundaries = [
        divisor(
            0,
            {
                ("proper", point, None): 1,
                **{("near", point, line): -1 for line in arrangement.incident_lines[point]},
            },
        )
        for point in triple_points
    ]
    ordinary_exceptionals = [
        divisor(0, {("proper", point, None): 1}) for point in double_points
    ]
    near_exceptionals = [
        divisor(0, {("near", point, line): 1}) for point, line in near_flags
    ]
    boundaries = [*line_boundaries, *exceptional_boundaries]
    exceptionals = [*ordinary_exceptionals, *near_exceptionals]
    walls = [*boundaries, *exceptionals]
    boundary_sum = tuple(sum(wall[index] for wall in boundaries) for index in range(rank))
    expected_boundary = (6, *([-2] * len(centers)))
    if boundary_sum != expected_boundary:
        raise ArithmeticError("the resolved A3 mirror boundary is not -2K")
    if any(_intersection(wall, wall) != -4 for wall in boundaries):
        raise ArithmeticError("an A3 boundary component is not a (-4)-curve")
    if any(_intersection(wall, wall) != -1 for wall in exceptionals):
        raise ArithmeticError("an A3 exceptional is not a (-1)-curve")
    if any(
        _intersection(boundaries[left], boundaries[right]) != 0
        for left in range(len(boundaries))
        for right in range(left)
    ):
        raise ArithmeticError("resolved A3 boundary components are not disjoint")
    audit = _cone_audit(walls, cdd_executable, store_rays=False)
    return {
        "field": "Q, characteristic zero",
        "proper_center_count": len(arrangement.points),
        "first_near_center_count": len(near_flags),
        "blowup_count": len(centers),
        "picard_rank": rank,
        "hyperbolic_dimension": rank - 1,
        "boundary_component_count": len(boundaries),
        "boundary_class_sum": list(boundary_sum),
        "boundary_norms": [_intersection(wall, wall) for wall in boundaries],
        "boundary_components_pairwise_disjoint": True,
        "obvious_exceptional_count": len(exceptionals),
        "obvious_exceptional_norms": [_intersection(wall, wall) for wall in exceptionals],
        "strong_coble_argument": (
            "The ten reduced disjoint rational (-4)-curves sum to -2K.  Negative intersection "
            "forces their sum to be the unique bi-anticanonical divisor, and reducedness excludes |-K|."
        ),
        "strong_coble": True,
        "declared_wall_order": "six line boundaries; four triple-exceptional boundaries; three double exceptionals; twelve near exceptionals",
        "declared_wall_vectors": [list(wall) for wall in walls],
        "declared_wall_cone": audit,
        "scope_warning": (
            "This 25-curve stock is an outer approximation to the full nef cone; hyperideal rays "
            "prove only that the declared stock is incomplete/infinite-volume."
        ),
    }


def _h3_winger_carriers(
    arrangement: Arrangement,
    cdd_executable: str,
    *,
    include_degree_two: bool = False,
) -> dict[str, Any]:
    centers = _h3_triple_orbit(arrangement)
    edges = _h3_petersen_edges(arrangement, centers)
    rank = 11
    boundary = (3, *([-1] * 10))
    exceptionals = [(0, *[1 if index == chosen else 0 for index in range(10)]) for chosen in range(10)]
    mirror_lines = [
        (1, *[-1 if index in edge else 0 for index in range(10)])
        for edge in edges
    ]
    all_pair_lines = [
        (1, *[-1 if index in pair else 0 for index in range(10)])
        for pair in itertools.combinations(range(10), 2)
    ]
    five_point_conics = [
        (2, *[-1 if index in subset else 0 for index in range(10)])
        for subset in itertools.combinations(range(10), 5)
    ]
    carriers: dict[str, list[tuple[int, ...]]] = {
        "petersen": [*exceptionals, *mirror_lines],
        "petersen_plus_boundary": [boundary, *exceptionals, *mirror_lines],
        "all_pair_lines_plus_boundary": [boundary, *exceptionals, *all_pair_lines],
    }
    if include_degree_two:
        carriers["degree_two_plus_boundary"] = [
            boundary,
            *exceptionals,
            *all_pair_lines,
            *five_point_conics,
        ]
    audits = {}
    for name, walls in carriers.items():
        audits[name] = _cone_audit(walls, cdd_executable)
    graph = nx.Graph()
    graph.add_nodes_from(range(10))
    graph.add_edges_from(edges)
    independent_sets = [
        subset
        for subset in itertools.combinations(range(10), 4)
        if all(not graph.has_edge(*pair) for pair in itertools.combinations(subset, 2))
    ]
    witnesses = []
    for subset in independent_sets:
        ray = (1, *[-1 if index in subset else 0 for index in range(10)])
        active = [wall for wall in carriers["petersen"] if _intersection(wall, ray) == 0]
        witnesses.append(
            {
                "independent_set": list(subset),
                "ray": list(ray),
                "square": _intersection(ray, ray),
                "active_wall_count": len(active),
                "active_wall_rank": sp.Matrix(active).rank(),
            }
        )
    independent_triples = [
        subset
        for subset in itertools.combinations(range(10), 3)
        if all(not graph.has_edge(*pair) for pair in itertools.combinations(subset, 2))
    ]
    boundary_witnesses = []
    for subset in independent_triples:
        ray = (1, *[-1 if index in subset else 0 for index in range(10)])
        active = [
            wall
            for wall in carriers["petersen_plus_boundary"]
            if _intersection(wall, ray) == 0
        ]
        boundary_witnesses.append(
            {
                "independent_set": list(subset),
                "ray": list(ray),
                "square": _intersection(ray, ray),
                "boundary_pairing": _intersection(boundary, ray),
                "all_wall_pairings_nonnegative": all(
                    _intersection(wall, ray) >= 0
                    for wall in carriers["petersen_plus_boundary"]
                ),
                "active_wall_count": len(active),
                "active_wall_rank": sp.Matrix(active).rank(),
            }
        )
    degree_two_witnesses = []
    degree_two_walls = [boundary, *exceptionals, *all_pair_lines, *five_point_conics]
    for subset in itertools.combinations(range(10), 7):
        ray = (5, *[-2 if index in subset else 0 for index in range(10)])
        active = [wall for wall in degree_two_walls if _intersection(wall, ray) == 0]
        degree_two_witnesses.append(
            {
                "seven_subset": list(subset),
                "ray": list(ray),
                "square": _intersection(ray, ray),
                "primitive_boundary_pairing": _intersection(boundary, ray),
                "all_wall_pairings_nonnegative": all(
                    _intersection(wall, ray) >= 0 for wall in degree_two_walls
                ),
                "active_wall_count": len(active),
                "active_wall_rank": sp.Matrix(active).rank(),
            }
        )
    return {
        "triple_centers": list(centers),
        "petersen_edges": [list(edge) for edge in edges],
        "petersen_graph_isomorphic": nx.is_isomorphic(graph, nx.petersen_graph()),
        "size_four_independent_hyperideal_witnesses": witnesses,
        "boundary_compatible_independent_triple_hyperideal_witnesses": boundary_witnesses,
        "degree_two_truncation_seven_subset_hyperideal_witnesses": degree_two_witnesses,
        "carriers": audits,
    }


def _arrangement_payload(arrangement: Arrangement) -> dict[str, Any]:
    return {
        "name": arrangement.name,
        "field": "Q" if arrangement.name == "A3" else "Q(sqrt(5))",
        "line_count": len(arrangement.lines),
        "point_count": len(arrangement.points),
        "intersection_multiplicity_histogram": arrangement.multiplicity_histogram,
        "pair_identity": arrangement.pair_identity,
        "expected_pair_count": math.comb(len(arrangement.lines), 2),
        "lines": [[str(value) for value in line] for line in arrangement.lines],
        "points": [
            {
                "coordinates": [str(value) for value in point],
                "incident_lines": list(arrangement.incident_lines[index]),
            }
            for index, point in enumerate(arrangement.points)
        ],
        "projective_arrangement_automorphism_order": len(_line_automorphisms(arrangement)),
    }


def build_certificate(
    *,
    all_a3_sizes: bool = False,
    include_degree_two: bool = False,
    skip_a3_hilbert: bool = False,
) -> dict[str, Any]:
    cdd = shutil.which("cddexec")
    if cdd is None:
        raise RuntimeError("cddexec is required")
    a3 = _arrangement("A3", _a3_lines())
    h3 = _arrangement("H3", _h3_lines())
    a3_orbits = _six_line_orbits(a3)
    h3_orbits = _six_line_orbits(h3)
    return {
        "schema": SCHEMA,
        "arrangements": {
            "A3": _arrangement_payload(a3),
            "H3": _arrangement_payload(h3),
        },
        "six_mirror_boundary_orbits": {
            "A3": a3_orbits,
            "H3": h3_orbits,
            "H3_summary": {
                "orbit_count": len(h3_orbits),
                "subset_count": sum(row["orbit_size"] for row in h3_orbits),
                "admissible_orbit_count": sum(row["ordinary_reduced_boundary_admissible"] for row in h3_orbits),
                "admissible_subset_count": sum(
                    row["orbit_size"] for row in h3_orbits if row["ordinary_reduced_boundary_admissible"]
                ),
                "admissible_subsets_by_triple_count": dict(
                    sorted(
                        Counter(
                            {
                                t: sum(
                                    row["orbit_size"]
                                    for row in h3_orbits
                                    if row["ordinary_reduced_boundary_admissible"] and row["triple_count"] == t
                                )
                                for t in range(5)
                            }
                        ).items()
                    )
                ),
            },
        },
        "A3_hilbert_census": (
            {"status": "deferred by --skip-a3-hilbert"}
            if skip_a3_hilbert
            else _a3_hilbert_census(a3, all_a3_sizes)
        ),
        "A3_full_six_mirror_surface": _a3_full_six_line_carrier(a3, cdd),
        "H3_Winger_surface": _h3_winger_surface(h3),
        "H3_winger_ten_node_lane": _h3_winger_carriers(
            h3,
            cdd,
            include_degree_two=include_degree_two,
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--all-a3-sizes", action="store_true")
    parser.add_argument("--include-degree-two", action="store_true")
    parser.add_argument("--skip-a3-hilbert", action="store_true")
    args = parser.parse_args()
    if args.skip_a3_hilbert and args.output == OUTPUT:
        args.output = OUTPUT.with_name("a3_h3_plane_arrangement_coble_factories_h3_only.json")
    certificate = build_certificate(
        all_a3_sizes=args.all_a3_sizes,
        include_degree_two=args.include_degree_two,
        skip_a3_hilbert=args.skip_a3_hilbert,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(args.output),
        "a3_survivors": certificate["A3_hilbert_census"].get("strong_coble_hilbert_survivor_count_by_size"),
        "h3_six_line_summary": certificate["six_mirror_boundary_orbits"]["H3_summary"],
        "h3_carriers": {
            key: {
                "walls": value["wall_count"],
                "rays": value["extreme_ray_count"],
                "finite": value["finite_vertex_count"],
                "ideal": value["ideal_vertex_count"],
                "hyperideal": value["hyperideal_vertex_count"],
            }
            for key, value in certificate["H3_winger_ten_node_lane"]["carriers"].items()
        },
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
