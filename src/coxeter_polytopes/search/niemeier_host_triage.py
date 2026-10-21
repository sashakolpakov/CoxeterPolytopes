"""Exact affine-chamber triage for the unscreened rooted Niemeier hosts.

This is a prioritisation tool, not a Vinberg completion.  For every rooted
rank-24 Niemeier type not already screened in this project it computes:

* root count, root-lattice determinant, and Niemeier glue index;
* the exact projective-ray census after adding the universal factorwise
  affine-highest-root packet for the silver unit coefficient
  ``a=1+sqrt(2)``;
* whether the discriminant glue is forced to cross the two components in
  the four best two-component cases.

The affine packet has roots ``(2, epsilon*theta_R)`` of norm two and height
``2+2sqrt(2)``.  Its dominant ``t=1`` section is the product of the affine
Weyl alcoves.  A vertex on finite node ``i`` has exact spatial energy

``4 * (C_R^-1)[i,i] / mark_i^2``.

Glue walls can precede or refine this packet, so the census is a rigorous
geometric triage invariant, not the claimed chamber of an unscreened host.
The full tetracode treatment of E6^4 lives in :mod:`silver_e6_tetracode`.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import isqrt
from typing import TypeAlias

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.bugaenko_vinberg import _fraction_inverse


Component: TypeAlias = tuple[str, int]
RationalMatrix: TypeAlias = tuple[tuple[Fraction, ...], ...]


@dataclass(frozen=True, slots=True)
class NiemeierHost:
    name: str
    components: tuple[Component, ...]
    glue_assessment: str
    minimum_coupling_support: int | None
    glue_evidence: str


UNSEARCHED_ROOTED_NIEMEIER_HOSTS = (
    NiemeierHost("A2^12", (("A", 2),) * 12, "many-factor code", None, "backend required"),
    NiemeierHost("A3^8", (("A", 3),) * 8, "many-factor code", None, "backend required"),
    NiemeierHost("A4^6", (("A", 4),) * 6, "many-factor code", None, "backend required"),
    NiemeierHost(
        "A5^4 D4", (("A", 5),) * 4 + (("D", 4),), "five-factor mixed glue", None, "backend required"
    ),
    NiemeierHost("A6^4", (("A", 6),) * 4, "four-factor glue", None, "backend required"),
    NiemeierHost(
        "A7^2 D5^2", (("A", 7),) * 2 + (("D", 5),) * 2, "four-factor mixed glue", None, "backend required"
    ),
    NiemeierHost("A8^3", (("A", 8),) * 3, "three-factor glue", None, "backend required"),
    NiemeierHost(
        "A9^2 D6", (("A", 9),) * 2 + (("D", 6),), "three-factor mixed glue", None, "backend required"
    ),
    NiemeierHost(
        "A11 D7 E6", (("A", 11), ("D", 7), ("E", 6)), "three-factor mixed glue", None, "backend required"
    ),
    NiemeierHost(
        "A12^2",
        (("A", 12),) * 2,
        "forced cross-component discriminant glue",
        2,
        "det(A12)=13 is prime; a unimodular index-13 extension is the graph of an anti-isometry",
    ),
    NiemeierHost(
        "A15 D9",
        (("A", 15), ("D", 9)),
        "forced cross-component discriminant glue",
        2,
        "neither rank-15 nor rank-9 factor admits an even unimodular lattice; the index-8 glue must cross",
    ),
    NiemeierHost(
        "A17 E7",
        (("A", 17), ("E", 7)),
        "forced cross-component discriminant glue",
        2,
        "neither rank-17 nor rank-7 factor admits an even unimodular lattice; the index-6 glue must cross",
    ),
    NiemeierHost("D4^6", (("D", 4),) * 6, "many-factor code", None, "backend required"),
    NiemeierHost("D6^4", (("D", 6),) * 4, "four-factor glue", None, "backend required"),
    NiemeierHost("D8^3", (("D", 8),) * 3, "three-factor glue", None, "backend required"),
    NiemeierHost(
        "D10 E7^2", (("D", 10), ("E", 7), ("E", 7)), "three-factor binary glue", None, "backend required"
    ),
    NiemeierHost(
        "D12^2",
        (("D", 12),) * 2,
        "forced cross-component discriminant glue",
        2,
        "D12 has no proper even integral same-rank overlattice, so the index-4 glue must cross",
    ),
    NiemeierHost(
        "D16 E8",
        (("D", 16), ("E", 8)),
        "decomposable D16+ plus E8 glue",
        None,
        "the index-two spinor extension is internal to D16 and does not couple E8",
    ),
    NiemeierHost(
        "E6^4",
        (("E", 6),) * 4,
        "explicit ternary tetracode",
        3,
        "the implemented [4,2,3] tetracode has eight nonzero words, all of support three",
    ),
)


def cartan_gram(component: Component) -> RationalMatrix:
    family, rank = component
    gram = [[Fraction(2 if i == j else 0) for j in range(rank)] for i in range(rank)]
    if family == "A":
        edges = tuple((i, i + 1) for i in range(rank - 1))
    elif family == "D" and rank >= 4:
        edges = (*tuple((i, i + 1) for i in range(rank - 3)), (rank - 3, rank - 2), (rank - 3, rank - 1))
    elif family == "E" and rank in (6, 7, 8):
        # Arms of lengths 2, rank-4, and 1 at vertex two.
        edges = (*tuple((i, i + 1) for i in range(rank - 2)), (2, rank - 1))
    else:
        raise ValueError(f"unsupported simply-laced component {component}")
    for left, right in edges:
        gram[left][right] = gram[right][left] = Fraction(-1)
    return tuple(tuple(row) for row in gram)


@lru_cache(maxsize=None)
def highest_root_marks(component: Component) -> tuple[int, ...]:
    family, rank = component
    if family == "A":
        return (1,) * rank
    if family == "D":
        return (1, *((2,) * (rank - 3)), 1, 1)
    gram = cartan_gram(component)
    roots: set[tuple[int, ...]] = set()
    stack: list[tuple[int, ...]] = []
    for index in range(rank):
        root = tuple(1 if entry == index else 0 for entry in range(rank))
        for signed in (root, tuple(-value for value in root)):
            roots.add(signed)
            stack.append(signed)
    while stack:
        root = stack.pop()
        for index in range(rank):
            pairing = sum(gram[index][column] * root[column] for column in range(rank))
            reflected = list(root)
            reflected[index] -= int(pairing)
            reflected_tuple = tuple(reflected)
            if reflected_tuple not in roots:
                roots.add(reflected_tuple)
                stack.append(reflected_tuple)
    positive = tuple(root for root in roots if all(value >= 0 for value in root))
    return max(positive, key=sum)


@lru_cache(maxsize=None)
def affine_vertex_energies(component: Component) -> tuple[Fraction, ...]:
    inverse = _fraction_inverse(cartan_gram(component))
    marks = highest_root_marks(component)
    return (
        Fraction(0),
        *tuple(
            Fraction(4) * inverse[index][index] / (marks[index] * marks[index])
            for index in range(len(marks))
        ),
    )


def _component_invariants(component: Component) -> tuple[int, int, int]:
    family, rank = component
    if family == "A":
        return rank + 1, rank + 1, rank * (rank + 1)
    if family == "D":
        return 4, 2 * rank - 2, rank * (2 * rank - 2)
    values = {6: (3, 12), 7: (2, 18), 8: (1, 30)}
    determinant, coxeter = values[rank]
    return determinant, coxeter, rank * coxeter


@dataclass(frozen=True, slots=True)
class NiemeierTriageRow:
    rank: int
    host: str
    components: tuple[Component, ...]
    component_count: int
    root_count: int
    root_lattice_determinant: int
    glue_index: int
    glue_assessment: str
    minimum_coupling_support: int | None
    glue_evidence: str
    affine_product_ray_count: int
    affine_finite_count: int
    affine_cusp_count: int
    affine_hyperideal_count: int
    priority_reason: str


def _affine_counts(host: NiemeierHost) -> tuple[int, int, int, int]:
    distribution = {Fraction(0): 1}
    for component in host.components:
        updated: dict[Fraction, int] = {}
        for energy, multiplicity in distribution.items():
            for vertex_energy in affine_vertex_energies(component):
                total = energy + vertex_energy
                updated[total] = updated.get(total, 0) + multiplicity
        distribution = updated
    coefficient = Silver(1, 1)
    finite = sum(
        multiplicity
        for energy, multiplicity in distribution.items()
        if (Silver(energy) - coefficient).sign(0) < 0
    )
    cusps = sum(
        multiplicity
        for energy, multiplicity in distribution.items()
        if (Silver(energy) - coefficient).sign(0) == 0
    )
    total = sum(distribution.values())
    return total, finite, cusps, total - finite - cusps


@lru_cache(maxsize=1)
def niemeier_triage() -> tuple[NiemeierTriageRow, ...]:
    preliminary = []
    for host in UNSEARCHED_ROOTED_NIEMEIER_HOSTS:
        determinant = 1
        root_count = 0
        for component in host.components:
            component_determinant, _coxeter, component_roots = _component_invariants(component)
            determinant *= component_determinant
            root_count += component_roots
        glue_index = isqrt(determinant)
        if glue_index * glue_index != determinant:
            raise AssertionError(f"{host.name} root determinant is not a square")
        total, finite, cusps, hyperideal = _affine_counts(host)
        if host.glue_assessment.startswith("forced cross"):
            tier = 0
            reason = "two factors, forced cross-glue, and the smallest exact affine obstruction"
        elif host.name == "E6^4":
            tier = 2
            reason = "explicit cheap glue backend; now screened through its second coupled packet"
        elif host.name == "D16 E8":
            tier = 4
            reason = "small affine census but its glue is internal and the spatial lattice decomposes"
        elif len(host.components) == 3:
            tier = 1
            reason = "three factors and a moderate affine product; decode glue next"
        else:
            tier = 3
            reason = "four or more factors give a much larger product obstruction"
        preliminary.append((tier, hyperideal, -finite, host.name, host, root_count, determinant, glue_index, total, finite, cusps))
    preliminary.sort(key=lambda item: item[:4])
    rows = []
    for priority, item in enumerate(preliminary, start=1):
        _, hyperideal, _negative_finite, _name, host, root_count, determinant, glue_index, total, finite, cusps = item
        if host.glue_assessment.startswith("forced cross"):
            reason = "two factors, forced cross-glue, and the smallest exact affine obstruction"
        elif host.name == "E6^4":
            reason = "explicit cheap glue backend; screened through its second coupled packet"
        elif host.name == "D16 E8":
            reason = "small affine census but internal glue and a decomposable spatial lattice"
        elif len(host.components) == 3:
            reason = "three factors and a moderate affine product; exact glue decoder is the next gate"
        else:
            reason = "four or more factors give a substantially larger product obstruction"
        rows.append(
            NiemeierTriageRow(
                priority,
                host.name,
                host.components,
                len(host.components),
                root_count,
                determinant,
                glue_index,
                host.glue_assessment,
                host.minimum_coupling_support,
                host.glue_evidence,
                total,
                finite,
                cusps,
                hyperideal,
                reason,
            )
        )
    return tuple(rows)


__all__ = [
    "NiemeierHost",
    "NiemeierTriageRow",
    "UNSEARCHED_ROOTED_NIEMEIER_HOSTS",
    "affine_vertex_energies",
    "cartan_gram",
    "highest_root_marks",
    "niemeier_triage",
]
