"""Vertex-targeted first packet for ``<-(theta-1)> + N(A24)``."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from fractions import Fraction

from coxeter_polytopes.exact.cubic7 import CUBIC7_ZERO, THETA7, Cubic7
from coxeter_polytopes.search.cubic7_a1_golay import (
    A_UNIT,
    COMPLETE_MIRROR_NORM_PROOF,
    UNIT_A1_GOLAY_SEED,
    enumerate_control_levels,
)
from coxeter_polytopes.search.cubic7_a24_plus import (
    a24_basis_audit,
    a24_niemeier_membership,
)


FIRST_CONTROL = 2 + 5 * THETA7 + 2 * THETA7**2
FIRST_SPATIAL = 2 * THETA7 + THETA7**2
FIRST_HEIGHT = 5 + 14 * THETA7 + 6 * THETA7**2


@dataclass(frozen=True, slots=True)
class A24UnitAlcoveRay:
    node: int
    inverse_cartan_diagonal: Fraction
    norm: Cubic7
    kind: str


@dataclass(frozen=True, slots=True)
class A24PlusUnitFirstPacket:
    complete_control_height_bound: int
    control_level_count: int
    first_control: Cubic7
    first_height: Cubic7
    root_norm: Cubic7
    spatial_coset: str | None
    reflection_integral: bool
    cap: Cubic7
    rays: tuple[A24UnitAlcoveRay, ...]
    finite_vertex_count: int
    cusp_count: int
    hyperideal_ray_count: int
    finite_volume: bool
    compact: bool


@lru_cache(maxsize=1)
def unit_a24_plus_first_packet() -> A24PlusUnitFirstPacket:
    if not a24_basis_audit():
        raise AssertionError("the A24 Niemeier glue basis lost even unimodularity")
    levels = enumerate_control_levels(UNIT_A1_GOLAY_SEED, 100)
    level = next(
        item
        for item in levels
        if item.norm == 1 and item.control == FIRST_CONTROL
    )
    spatial = (
        FIRST_SPATIAL,
        *(CUBIC7_ZERO for _ in range(23)),
        -FIRST_SPATIAL,
    )
    norm = -A_UNIT * FIRST_CONTROL**2 + sum(
        (value * value for value in spatial), start=CUBIC7_ZERO
    )
    coset = a24_niemeier_membership(spatial)
    if norm != 1 or coset != "A24":
        raise AssertionError("the affine candidate is not a norm-one A24 root")
    reflection_integral = (2 * A_UNIT * FIRST_CONTROL).is_integral and all(
        (2 * value).is_integral for value in spatial
    )
    cap = A_UNIT * FIRST_CONTROL / FIRST_SPATIAL
    rays = [A24UnitAlcoveRay(0, Fraction(0), -A_UNIT, "finite")]
    for node in range(1, 25):
        diagonal = Fraction(node * (25 - node), 25)
        ray_norm = -A_UNIT + cap * cap * Cubic7(diagonal)
        sign = ray_norm.sign(0)
        rays.append(
            A24UnitAlcoveRay(
                node,
                diagonal,
                ray_norm,
                "finite" if sign < 0 else "cusp" if sign == 0 else "hyperideal",
            )
        )
    result = tuple(rays)
    finite = sum(ray.kind == "finite" for ray in result)
    cusps = sum(ray.kind == "cusp" for ray in result)
    hyperideal = sum(ray.kind == "hyperideal" for ray in result)
    return A24PlusUnitFirstPacket(
        100,
        len(levels),
        FIRST_CONTROL,
        FIRST_HEIGHT,
        norm,
        coset,
        reflection_integral,
        cap,
        result,
        finite,
        cusps,
        hyperideal,
        hyperideal == 0,
        hyperideal == 0 and cusps == 0,
    )


__all__ = ["unit_a24_plus_first_packet"]
