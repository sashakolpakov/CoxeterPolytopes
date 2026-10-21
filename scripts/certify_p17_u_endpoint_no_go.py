#!/usr/bin/env python3
"""Exact P17 rank-two U-endpoint profile and scoped carrier obstruction.

The 22 facet vectors of Vinberg's odd-unimodular P17 chamber admit eight
rank-two contraction endpoints whose surviving lattice is U.  This replay
recovers, without a chosen marking, the same six surviving curve profiles at
every endpoint.  The profiles force a local jet pattern which is impossible
both on an elliptic K3 surface with NS=U and on a geometrically ruled surface
over any smooth base in characteristic zero.

The geometric contradiction is proved in the accompanying report.  This
script freezes its exact lattice and vanishing-order inputs; it does not claim
nonexistence of the abstract P17 Coxeter polytope or of carriers on arbitrary
algebraic surfaces.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Sequence

try:
    from scripts.certify_rational_i1n_inverse_design import enumerate_simple_roots
    from scripts.certify_rational_p17_mmp import (
        contract_state,
        enumerate_contraction_levels,
        surface_pair,
        vector_add,
        vector_sum,
    )
except ModuleNotFoundError:  # direct ``python scripts/...`` replay
    from certify_rational_i1n_inverse_design import enumerate_simple_roots
    from certify_rational_p17_mmp import (
        contract_state,
        enumerate_contraction_levels,
        surface_pair,
        vector_add,
        vector_sum,
    )


Vector = tuple[int, ...]

LABEL_BY_BASE_CLASS = {
    (0, 1): "s",
    (1, 0): "f",
    (4, 2): "C",
    (5, 2): "Q",
    (6, 3): "A",
    (8, 4): "S",
}

EXPECTED_MULTIPLICITIES = {
    "s": (0,) * 16,
    "f": (1, 1) + (0,) * 14,
    "C": (1,) * 10 + (0,) * 6,
    "Q": (1,) * 14 + (0,) * 2,
    "A": (2,) + (1,) * 15,
    "S": (2,) * 6 + (1,) * 10,
}


def endpoint_image(curve: Vector, exceptional_basis: Sequence[Vector]) -> Vector:
    """Pull the endpoint image back to the original rank-18 lattice."""

    image = curve
    for exceptional in exceptional_basis:
        image = vector_add(
            image, exceptional, surface_pair(image, exceptional)
        )
    return image


def exceptional_chain_order(
    roots: Sequence[Vector], exceptional_basis: Sequence[Vector]
) -> tuple[Vector, ...]:
    """Recover E_1,...,E_16 from the contracted P17 facet chain.

    In an orthogonal total-transform basis, a contracted chain wall has
    coordinates E_i-E_(i+1), and the terminal wall has coordinate E_16.
    The calculation below discovers this path rather than assuming facet
    labels from one preferred endpoint.
    """

    vanished_coordinates: list[tuple[int, ...]] = []
    for root in roots:
        if not any(endpoint_image(root, exceptional_basis)):
            # If root=sum(c_i E_i) and E_i^2=-1, then c_i=-root.E_i.
            vanished_coordinates.append(
                tuple(-surface_pair(root, exceptional) for exceptional in exceptional_basis)
            )

    if len(vanished_coordinates) != 16:
        raise AssertionError("a U endpoint must contract exactly 16 facet curves")

    adjacency = {index: [] for index in range(16)}
    terminal_candidates: list[int] = []
    for coordinates in vanished_coordinates:
        support = [(index, value) for index, value in enumerate(coordinates) if value]
        if len(support) == 1:
            index, value = support[0]
            if value != 1:
                raise AssertionError("terminal exceptional wall is not primitive positive")
            terminal_candidates.append(index)
        elif len(support) == 2:
            (left, left_value), (right, right_value) = support
            if sorted((left_value, right_value)) != [-1, 1]:
                raise AssertionError("contracted wall is not an exceptional difference")
            adjacency[left].append(right)
            adjacency[right].append(left)
        else:
            raise AssertionError("contracted facets do not form an exceptional chain")

    if len(terminal_candidates) != 1:
        raise AssertionError("exceptional chain has no unique terminal wall")

    # Walk from E_16 to E_1, then reverse into blow-up order.
    path = [terminal_candidates[0]]
    previous: int | None = None
    current = path[0]
    while True:
        following = [index for index in adjacency[current] if index != previous]
        if not following:
            break
        if len(following) != 1:
            raise AssertionError("contracted exceptional graph branches")
        previous, current = current, following[0]
        path.append(current)
    if len(path) != 16:
        raise AssertionError("contracted exceptional graph is disconnected")
    return tuple(exceptional_basis[index] for index in reversed(path))


def endpoint_profile(
    roots: Sequence[Vector],
    state: Sequence[Vector],
    exceptional_basis: Sequence[Vector],
) -> dict[str, object]:
    negative = [curve for curve in state if surface_pair(curve, curve) < 0]
    isotropic = [curve for curve in state if surface_pair(curve, curve) == 0]
    if len(negative) != 1 or len(isotropic) != 1:
        raise AssertionError("unexpected rank-two endpoint cone")
    section = negative[0]
    fiber = isotropic[0]
    if surface_pair(section, section) != -2 or surface_pair(section, fiber) != 1:
        raise AssertionError("endpoint lattice is not U in section-fiber form")

    chain = exceptional_chain_order(roots, exceptional_basis)
    profiles: dict[str, dict[str, object]] = {}
    surviving_root_indices: list[int] = []
    for root_index, root in enumerate(roots):
        image = endpoint_image(root, exceptional_basis)
        if not any(image):
            continue
        surviving_root_indices.append(root_index)
        section_coefficient = surface_pair(image, fiber)
        fiber_coefficient = surface_pair(image, section) + 2 * section_coefficient
        base_class = (fiber_coefficient, section_coefficient)
        if base_class not in LABEL_BY_BASE_CLASS:
            raise AssertionError(f"unexpected surviving base class {base_class}")
        label = LABEL_BY_BASE_CLASS[base_class]
        multiplicities = tuple(surface_pair(root, exceptional) for exceptional in chain)
        if multiplicities != EXPECTED_MULTIPLICITIES[label]:
            raise AssertionError(f"unexpected {label} cluster {multiplicities}")
        profiles[label] = {
            "facet_index_zero_based": root_index,
            "base_class_f_s": list(base_class),
            "multiplicities_p1_through_p16": list(multiplicities),
            "strict_transform_square": surface_pair(root, root),
        }

    if set(profiles) != set(EXPECTED_MULTIPLICITIES):
        raise AssertionError("the six required U-endpoint curves were not recovered")

    canonical = vector_sum(exceptional_basis)
    if surface_pair(canonical, canonical) != -16:
        raise AssertionError("a 16-fold K3 blow-up must have K^2=-16")
    arithmetic_genera = []
    for root in roots:
        numerator = surface_pair(root, root) + surface_pair(canonical, root)
        if numerator % 2:
            raise AssertionError("adjunction parity failure")
        genus = 1 + numerator // 2
        if genus < 0:
            raise AssertionError("negative facet arithmetic genus")
        arithmetic_genera.append(genus)

    return {
        "surviving_facet_indices_zero_based": surviving_root_indices,
        "profiles": profiles,
        "blown_up_k3_canonical_square": -16,
        "facet_arithmetic_genera_for_k3_endpoint": arithmetic_genera,
    }


def u_pair(left: Sequence[int], right: Sequence[int]) -> int:
    """Pair classes recorded as (fiber coefficient, section coefficient)."""

    left_f, left_s = left
    right_f, right_s = right
    return left_f * right_s + left_s * right_f - 2 * left_s * right_s


def build_certificate() -> dict[str, object]:
    roots = tuple(tuple(root) for root in enumerate_simple_roots(17))
    levels, _ = enumerate_contraction_levels(roots)

    u_endpoints: list[dict[str, object]] = []
    rank_two_f1_basis_count = 0
    for state, bases in levels[2].items():
        negative = [curve for curve in state if surface_pair(curve, curve) < 0]
        isotropic = [curve for curve in state if surface_pair(curve, curve) == 0]
        if len(negative) != 1 or len(isotropic) != 1:
            raise AssertionError("unexpected rank-two state")
        section_square = surface_pair(negative[0], negative[0])
        if section_square == -1:
            rank_two_f1_basis_count += len(bases)
            continue
        if section_square != -2:
            raise AssertionError("rank-two endpoint is neither F1-type nor U-type")
        for basis in bases:
            u_endpoints.append(
                endpoint_profile(roots, state, tuple(sorted(basis)))
            )

    if len(u_endpoints) != 8 or rank_two_f1_basis_count != 8:
        raise AssertionError("expected eight U and eight F1-type endpoint markings")

    # All endpoints have identical unlabeled base classes and clusters.
    normalized_profiles = {
        json.dumps(endpoint["profiles"], sort_keys=True)
        for endpoint in u_endpoints
    }
    # Facet indices differ under diagram automorphisms, so compare without them.
    unlabeled_profiles = {
        tuple(
            (
                label,
                tuple(profile["base_class_f_s"]),
                tuple(profile["multiplicities_p1_through_p16"]),
                profile["strict_transform_square"],
            )
            for label, profile in sorted(endpoint["profiles"].items())
        )
        for endpoint in u_endpoints
    }
    if len(unlabeled_profiles) != 1:
        raise AssertionError("the eight U endpoints do not have one common profile")

    common = next(iter(unlabeled_profiles))
    common_dict = {
        label: {
            "base_class_f_s": list(base_class),
            "multiplicities_p1_through_p16": list(multiplicities),
            "strict_transform_square": square,
        }
        for label, base_class, multiplicities, square in common
    }
    c_class = common_dict["C"]["base_class_f_s"]
    a_class = common_dict["A"]["base_class_f_s"]
    f_class = common_dict["f"]["base_class_f_s"]
    if u_pair(c_class, f_class) != 2 or u_pair(c_class, a_class) != 12:
        raise AssertionError("unexpected base intersection numbers")
    local_c_a = sum(
        left * right
        for left, right in zip(
            common_dict["C"]["multiplicities_p1_through_p16"],
            common_dict["A"]["multiplicities_p1_through_p16"],
        )
    )
    if local_c_a != 11:
        raise AssertionError("the required C-A local intersection is not 11")

    certificate: dict[str, object] = {
        "schema": "p17_u_endpoint_no_go_v1",
        "characteristic": 0,
        "abstract_chamber": "classical Vinberg P17",
        "facet_count": len(roots),
        "rank_two_state_count": len(levels[2]),
        "u_endpoint_count": len(u_endpoints),
        "f1_type_endpoint_count": rank_two_f1_basis_count,
        "all_u_endpoints_have_one_common_profile": len(unlabeled_profiles) == 1,
        "common_u_endpoint_profile": common_dict,
        "endpoint_facet_labelings": u_endpoints,
        "forced_intersections": {
            "C_dot_f": 2,
            "C_dot_A": 12,
            "cluster_local_Ip_C_A": local_c_a,
        },
        "elliptic_k3_obstruction": {
            "D": "s+2f",
            "C": "2D",
            "A": "3D",
            "ramification_from_C_cap_f": "C|f=2p forces y(p)=0",
            "irreducible_A_term": "the y coefficient in H0(3D) is nonzero",
            "local_derivative": "in coordinates (t,y), partial(A)/partial(y) at p equals the nonzero y coefficient",
            "contradiction": "the cluster requires mult_p(A)=2 but every irreducible A is smooth at p",
        },
        "geometrically_ruled_obstruction": {
            "splitting": "the barycenter of the finite flat bisection C splits the affine-line torsor",
            "split_model": "P_B(O + L^-1) with deg(L)=2",
            "centered_bisection": "v^2=delta",
            "restricted_cubic": "A|C=P+vQ with P in H0(B,L^3), Q in H0(B,L^2)",
            "required_local_order": 11,
            "required_base_orders_for_order_at_least_11": {
                "P": 6,
                "Q": 5,
            },
            "degree_bounds": {"deg_L3": 6, "deg_L2": 4},
            "conclusion": "Q must vanish identically; P then has local pullback order 12 or vanishes identically, never exactly 11",
        },
        "geometric_conclusion": (
            "No one of the eight P17 U-endpoint markings is realized by a "
            "16-fold blow-up of an elliptic K3 surface with NS=U or of a "
            "geometrically ruled surface over a smooth base in characteristic zero."
        ),
        "not_claimed": (
            "nonexistence of the abstract P17 finite-volume Coxeter polytope, "
            "or nonexistence of P17 carriers on all algebraic surfaces"
        ),
    }
    # Silence a lint false-positive while retaining the stronger labeled replay
    # above in the certificate.
    if not normalized_profiles:
        raise AssertionError("no labeled endpoint profiles")
    payload = json.dumps(certificate, sort_keys=True, separators=(",", ":")).encode()
    certificate["certificate_payload_sha256"] = hashlib.sha256(payload).hexdigest()
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
