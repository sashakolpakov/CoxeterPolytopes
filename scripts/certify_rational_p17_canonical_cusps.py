#!/usr/bin/env python3
"""Exact dual-cone/cusp certificate for the rational P17 obstruction.

This is independent of the finite surface-MMP state enumeration in
``certify_rational_p17_mmp.py``.  Write a characteristic rational-surface
canonical class as

    K=(a, 1-2*c_1, ..., 1-2*c_17),

in the standard I_(1,17) basis.  Adjunction for the 22 P17 facet curves
gives a rational polyhedron in ``(a,c_1,...,c_17)``.  For every primitive
extreme nef ray N, cddlib decides exactly whether ``K.N <= -1`` is
linearly feasible and whether that branch is bounded in ``a``.

There are 807 nef rays.  Only 72 negative branches are feasible: 64 have
the exact global bound a<=123, and the other eight are one orbit of the
isotropic cusp H-E_1 under the full order-eight weighted-diagram group.
On H-E_1, adjunction with the line facet and K^2=-8 force K to be the
standard canonical class.  Conjugation proves the same statement on all
eight unbounded cusp branches.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import gcd
from pathlib import Path
import shutil
import subprocess
from typing import Sequence

from sympy import Matrix

try:
    from scripts.certify_rational_i1n_inverse_design import (
        enumerate_simple_roots,
        exact_extreme_rays,
    )
    from scripts.certify_rational_p17_mmp import weighted_diagram_automorphisms
except ModuleNotFoundError:  # direct ``python scripts/...`` replay
    from certify_rational_i1n_inverse_design import (
        enumerate_simple_roots,
        exact_extreme_rays,
    )
    from certify_rational_p17_mmp import weighted_diagram_automorphisms


Vector = tuple[int, ...]
RationalVector = tuple[Fraction, ...]


def surface_pair(left: Sequence[int], right: Sequence[int]) -> int:
    return left[0] * right[0] - sum(
        a * b for a, b in zip(left[1:], right[1:])
    )


def primitive(vector: Sequence[int]) -> Vector:
    divisor = 0
    for entry in vector:
        divisor = gcd(divisor, abs(entry))
    if divisor == 0:
        raise ValueError("zero vector")
    result = tuple(entry // divisor for entry in vector)
    first = next(entry for entry in result if entry)
    return result if first > 0 else tuple(-entry for entry in result)


def nef_extreme_rays(roots: Sequence[Vector]) -> tuple[Vector, ...]:
    # A curve (d,e_i) pairs with a nef class (h,x_i) as
    # d*h-sum(e_i*x_i), so these are the H-rows in actual divisor
    # coordinates (h,x_i).
    rows = [(root[0], *(-entry for entry in root[1:])) for root in roots]
    return exact_extreme_rays(rows, 807)


def adjunction_polyhedron_rows() -> tuple[tuple[int, ...], ...]:
    """Affine H-rows [constant,a,c1,...,c17], all interpreted >=0."""

    rows: list[tuple[int, ...]] = []

    # a>=-3 (also follows from the line inequality and c_i>=0).
    rows.append((3, 1, *((0,) * 17)))

    # c_1>=...>=c_17>=0 from the B17 chain adjunction inequalities.
    for index in range(16):
        coefficients = [0] * 18
        coefficients[1 + index] = 1
        coefficients[2 + index] = -1
        rows.append((0, *coefficients))
    coefficients = [0] * 18
    coefficients[-1] = 1
    rows.append((0, *coefficients))

    # L, C, Q, A, S.  Repeated indices encode multiplicity two.
    positive_facets = (
        (3, 1, (0, 1, 2)),
        (11, 3, tuple(range(11))),
        (16, 4, (0, *range(15))),
        (18, 4, tuple(range(17))),
        (24, 6, (*range(7), *range(17))),
    )
    for constant, a_coefficient, indices in positive_facets:
        coefficients = [0] * 18
        coefficients[0] = a_coefficient
        for index in indices:
            coefficients[1 + index] -= 2
        rows.append((constant, *coefficients))
    return tuple(rows)


def negative_nef_row(nef: Vector) -> tuple[int, ...]:
    """The affine inequality K.N<=-1 in the (a,c_i) coordinates."""

    height = nef[0]
    multiplicities = tuple(-entry for entry in nef[1:])
    if any(value < 0 for value in multiplicities):
        raise AssertionError("nef multiplicities are not nonnegative")
    coefficients = [-height, *(2 * value for value in multiplicities)]
    return (-1 - sum(multiplicities), *coefficients)


def affine_v_representation(rows: Sequence[Sequence[int]]) -> tuple[RationalVector, ...]:
    executable = shutil.which("cddexec_gmp")
    if executable is None:
        raise RuntimeError("cddexec_gmp is required for the exact P17 cusp replay")
    column_count = len(rows[0])
    payload = "\n".join(
        (
            "H-representation",
            "begin",
            f"{len(rows)} {column_count} rational",
            *(" ".join(str(entry) for entry in row) for row in rows),
            "end",
            "",
        )
    )
    completed = subprocess.run(
        [executable, "--rep"],
        input=payload,
        check=True,
        capture_output=True,
        text=True,
    )
    lines = [line.strip() for line in completed.stdout.splitlines()]
    representation = lines.index("V-representation")
    begin = lines.index("begin", representation)
    count = int(lines[begin + 1].split()[0])
    return tuple(
        tuple(Fraction(entry) for entry in line.split())
        for line in lines[begin + 2 : begin + 2 + count]
    )


def branch_bound(
    base_rows: Sequence[Sequence[int]], nef: Vector
) -> Fraction | str | None:
    vertices_and_rays = affine_v_representation(
        (*base_rows, negative_nef_row(nef))
    )
    if not vertices_and_rays:
        return None
    rays = [row for row in vertices_and_rays if row[0] == 0]
    # row[1] is the a-coordinate.
    if any(row[1] > 0 for row in rays):
        return "unbounded"
    vertices = [row for row in vertices_and_rays if row[0] == 1]
    if not vertices:
        raise AssertionError("nonempty pointed branch has no vertex")
    return max(row[1] for row in vertices)


def diagram_lift(
    roots: Sequence[Vector], permutation: Sequence[int]
) -> Matrix:
    # R_1,...,R_17,L form a unimodular row basis.
    basis_indices = (*range(17), 17)
    basis = Matrix([roots[index] for index in basis_indices])
    image_basis = Matrix([roots[permutation[index]] for index in basis_indices])
    transformation = basis.inv() * image_basis
    if any(entry.q != 1 for entry in transformation):
        raise AssertionError("diagram lift is not integral")
    all_roots = Matrix(roots)
    if all_roots * transformation != Matrix(
        [roots[permutation[index]] for index in range(len(roots))]
    ):
        raise AssertionError("diagram lift fails on a facet root")
    return transformation


def cusp_orbit(
    roots: Sequence[Vector], automorphisms: Sequence[Sequence[int]]
) -> tuple[Vector, ...]:
    standard_cusp = (1, -1, *((0,) * 16))
    orbit = set()
    for permutation in automorphisms:
        lift = diagram_lift(roots, permutation)
        image = Matrix([standard_cusp]) * lift
        orbit.add(primitive(tuple(int(entry) for entry in image.row(0))))
    return tuple(sorted(orbit))


def canonical_digest(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def build_certificate() -> dict[str, object]:
    roots = tuple(tuple(root) for root in enumerate_simple_roots(17))
    nef_rays = nef_extreme_rays(roots)
    base_rows = adjunction_polyhedron_rows()

    infeasible: list[int] = []
    bounded: list[dict[str, object]] = []
    unbounded: list[tuple[int, Vector]] = []
    for index, nef in enumerate(nef_rays):
        bound = branch_bound(base_rows, nef)
        if bound is None:
            infeasible.append(index)
        elif bound == "unbounded":
            unbounded.append((index, nef))
        else:
            if not isinstance(bound, Fraction):
                raise AssertionError(bound)
            bounded.append(
                {
                    "nef_ray_index": index,
                    "nef_ray": list(nef),
                    "maximum_a": str(bound),
                }
            )

    automorphisms = weighted_diagram_automorphisms(roots)
    orbit = cusp_orbit(roots, automorphisms)
    unbounded_rays = tuple(sorted(nef for _, nef in unbounded))

    checks = {
        "nef_extreme_ray_count_is_807": len(nef_rays) == 807,
        "negative_branch_partition_is_735_plus_64_plus_8": (
            len(infeasible), len(bounded), len(unbounded)
        )
        == (735, 64, 8),
        "largest_exact_bounded_a_is_123": max(
            Fraction(record["maximum_a"]) for record in bounded
        )
        == 123,
        "diagram_automorphism_order_is_8": len(automorphisms) == 8,
        "unbounded_branches_equal_standard_cusp_orbit": unbounded_rays == orbit,
        "all_unbounded_nef_rays_are_isotropic": all(
            surface_pair(nef, nef) == 0 for nef in unbounded_rays
        ),
    }
    if not all(checks.values()):
        raise AssertionError(checks)

    certificate: dict[str, object] = {
        "schema": "rational_p17_canonical_cusps_v1",
        "lattice": "I_(1,17)",
        "facet_count": len(roots),
        "nef_extreme_ray_count": len(nef_rays),
        "nef_extreme_rays_sha256": canonical_digest(nef_rays),
        "adjunction_H_rows": [list(row) for row in base_rows],
        "adjunction_coordinate_change": "K=(a,1-2*c_1,...,1-2*c_17)",
        "canonical_norm_equation": "sum_i c_i(c_i-1)=(a^2-9)/4",
        "negative_nef_branch_census": {
            "infeasible": len(infeasible),
            "bounded": len(bounded),
            "unbounded": len(unbounded),
            "largest_exact_bounded_a": "123",
        },
        "bounded_branches": bounded,
        "unbounded_branches": [
            {"nef_ray_index": index, "nef_ray": list(nef)}
            for index, nef in unbounded
        ],
        "weighted_diagram_automorphism_order": len(automorphisms),
        "standard_cusp": [1, -1, *((0,) * 16)],
        "standard_cusp_orbit": [list(nef) for nef in orbit],
        "standard_cusp_lemma": {
            "negative_pairing": (
                "K.(H-E1)=a+1-2*c1 is even; K.(H-E1)<=-1 therefore "
                "c1>=(a+3)/2"
            ),
            "line_adjunction": (
                "c1+c2+c3<=(a+3)/2; monotonicity then forces "
                "c1=(a+3)/2 and c2=...=c17=0"
            ),
            "norm_substitution": (
                "c1(c1-1)=(a^2-9)/4 with c1=(a+3)/2 gives a=-3"
            ),
            "conclusion": "K=(-3,1^17)",
        },
        "checks": checks,
        "role_in_no_go": (
            "The cusp lemma closes all linearly unbounded characteristic-class "
            "branches.  The independent exact surface-MMP certificate classifies "
            "the remaining markings and proves that every canonical class is in "
            "the same eight-element diagram orbit."
        ),
    }
    certificate["certificate_payload_sha256"] = canonical_digest(certificate)
    return certificate


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(build_certificate(), indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.write_text(rendered)


if __name__ == "__main__":
    main()
