"""Uniform global incidence rows and vertex statistics.

Import has no scientific side effects. Coefficients are rational polynomials
in n. The global builder binds P,Q,A,S to its full vertex occupation measures,
together with every old corner feature. No free unbound attachment is exported.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb, prod


PORTS = ("P", "Q", "A", "S")
VARIABLES = ("f1", "f2", "f3", "f4", "B", "E3", "K3", "K4", "B4", *PORTS,
             *(f"h{k}" for k in range(5, 10)), *(f"b{k}" for k in range(5, 10)),
             "h_large", "K_large", "b_large")


def finite_column(ranks):
    if not ranks or any(r < 1 for r in ranks):
        raise ValueError("Positive spherical component ranks required")
    n = sum(ranks)
    return {"P": comb(n, 3), "Q": 0,
            "A": sum(a*b*c for a, b, c in combinations(ranks, 3))}


def ideal_column(ranks, deletion_charges=None):
    if not ranks or any(r < 1 for r in ranks):
        raise ValueError("Positive affine factor ranks required")
    rank, vertices = sum(ranks), prod(r+1 for r in ranks)
    column = {"P": F(vertices*sum(r*(r-1) for r in ranks), 6),
              "Q": F(vertices*(rank*rank-sum(r*r for r in ranks)), 8), "A": 0}
    if deletion_charges is not None:
        if len(deletion_charges) != len(ranks) or any(h < 0 for h in deletion_charges):
            raise ValueError("One nonnegative sum over node deletions per affine factor")
        column["S"] = F(vertices, 2)*sum((F(h, r+1) for h, r in
                                         zip(deletion_charges, ranks)), F(0))
    return column


def distance_weight(distance):
    return F(1) if 1 <= distance <= 7 else F(1, 3) if 8 <= distance <= 15 else F(0)


def classical_histogram(kind, rank):
    """All ranks, without a truncation of the component rank parameter."""
    if kind == "A" and rank >= 1:
        return {d: rank-d for d in range(1, min(15, rank-1)+1)}
    if kind == "D" and rank >= 4:
        return {d: max(rank-2-d, 0)+2+int(d == 2)
                for d in range(1, min(15, rank-2)+1)}
    raise ValueError("This closed histogram helper covers A_r and D_r only")


def finite_classical_column(components):
    ranks = [r for kind, r in components]
    column = finite_column(ranks)
    n = sum(ranks)
    column["S"] = sum((F(n-d-1, 2)*distance_weight(d)*count
                       for kind, rank in components
                       for d, count in classical_histogram(kind, rank).items()), F(0))
    return column


def edge_column(f7, f15):
    if not 0 <= f7 <= f15:
        raise ValueError("Nested unordered-pair counts required")
    return {"S": F(2*f7+f15, 3)}


def tangent_cusp_filter(u, w, common_triple_cusp, original_coxeter=False):
    """m5 on the explicit g_ac=-1 touching-pair incidence stratum.

    u,w are minus the two intervening (normalized) Gram entries. No
    discrete margin is imposed for projected non-Coxeter three-faces.
    """
    u, w = F(u), F(w)
    if not (0 <= u < 1 and 0 <= w < 1):
        raise ValueError("Intervening facets must meet at positive acute angles")
    if common_triple_cusp:
        return True
    return u+w >= F(1, 2) if original_coxeter else u+w > 0


def polynomial_rows(parity):
    if parity not in ("odd", "even"):
        raise ValueError("Split the integer tail by parity")
    rows = []

    def add(name, sense, coefficients, source):
        rows.append(dict(name=name, sense=sense,
            coefficients={key: tuple(F(c) for c in (value if isinstance(value, tuple)
                                                   else (value,)))
                          for key, value in coefficients.items()}, source=source))

    def eq(name, coefficients):
        add(name, "=", coefficients, "actual face/link incidence")

    def ge(name, coefficients, source="actual face/link incidence"):
        add(name, ">=", coefficients, source)

    eq("edge_three_face_marginal", {"E3": 1, "f1": (-1, F(3, 2), F(-1, 2))})
    eq("two_three_face_marginal", {"K3": 1, "f2": (2, -1)})
    eq("three_four_face_marginal", {"K4": 1, "f3": (3, -1)})
    eq("bad_three_four_face_marginal", {"B4": 1, "B": (3, -1)})
    eq("four_face_count", {"f4": 1, "h_large": -1,
                           **{f"h{k}": -1 for k in range(5, 10)}})
    eq("four_face_facet_count", {"K4": 1, "K_large": -1,
                                 **{f"h{k}": -k for k in range(5, 10)}})
    eq("four_face_bad_count", {"B4": 1, "b_large": -1,
                               **{f"b{k}": -1 for k in range(5, 10)}})
    eq("three_face_edge_link_marginal", {"E3": 2, "P": -3, "Q": -4})
    eq("three_face_euler_link_marginal", {"K3": 2, "f3": -4, "P": -1, "Q": -2})
    eq("no_bad_in_five_facet_four_face", {"b5": 1})
    eq("no_bad_in_six_facet_four_face", {"b6": 1})
    for k, cap in ((7, 3), (8, 6), (9, 7)):
        ge(f"bad_cap_at_{k}", {f"h{k}": cap, f"b{k}": -1},
           "Prokhorov 1987 p405, Lemma 1 case proof")
    ge("large_four_face_facets", {"K_large": 1, "h_large": -10})
    ge("large_four_face_bad_cap", {"K_large": 1, "b_large": -1})
    ge("bad_subset_of_all", {"f3": 1, "B": -1})
    ge("finite_apex_capacity", {"A": 1, "B": -2})
    ge("ideal_equator_capacity", {"Q": 1, "B": -3})
    ge("right_trihedra_are_triangular_links", {"P": 1, "A": -1})
    ge("minimum_three_face_facets", {"K3": 1, "f3": -4, "B": -2})
    ge("minimum_three_face_edges", {"E3": 1, "f3": -6, "B": -3})
    ge("prokhorov_geometric_distance_charge", {"S": 1, "K3": 1, "B": 1, "f3": -7},
       "Prokhorov 1987 pp406,408-411, fixed h(d), not arbitrary weights")
    if parity == "odd":
        edge, facets, cofaces = ((0, 12), (3, -1)), ((-6, 6), (3, -1)), ((-16, 8), (5, -1))
    else:
        edge, facets, cofaces = ((-12, 12), (4, -1)), ((-12, 6), (4, -1)), ((-24, 8), (6, -1))
    for name, count, denominator, cap, total in (
            ("khovanskii_edges_of_three_faces", "E3", edge[1], edge[0], "f3"),
            ("khovanskii_facets_of_three_faces", "K3", facets[1], facets[0], "f3"),
            ("khovanskii_facets_of_four_faces", "K4", cofaces[1], cofaces[0], "f4")):
        ge(name, {count: denominator, total: cap},
           "weak closure of Prokhorov 1987 pp405-406 Khovanskii comparisons")
    return rows
