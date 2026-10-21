#!/usr/bin/env python3
"""Independent hostile replay of fixed-Jacobian K3 index blindness.

The executable deliberately imports no owner certifier.  It freezes the owner
bytes, reconstructs the Brauer-kernel and primitive-complement lattice models,
checks the marked affine E8 packets, and audits the exact numerical interface
with the separate B=13 finite-volume spanning theorem.

Ogg--Shafarevich existence, projectivity, and fibrewise equality with the
relative Jacobian are source theorems.  They are separated from the exact
linear-algebra replay below rather than represented as machine proofs.
"""

from __future__ import annotations

import argparse
from functools import reduce
import hashlib
import json
from math import gcd
from pathlib import Path
from typing import Any, Sequence

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]

PINS = {
    "owner_report": (
        "results/noncompact_k3_fixed_jacobian_index_blindness.md",
        "14520f6cbbe6dc95f559ebb36d1232dd572adaf54e55457ed434f6108f39852b",
    ),
    "owner_scope_audit": (
        "results/noncompact_k3_fixed_jacobian_index_blindness_hostile_scope_audit.md",
        "50ca142bd2f3d3d8ad7fb469ab31a56e9d5fcfdfde60999c1191f78fb75e76ef",
    ),
    "owner_certifier": (
        "scripts/certify_noncompact_k3_fixed_jacobian_index_blindness.py",
        "bb72c96a3aed2b02289856167938a643e06938f1fc6345aa9a3674bb6c2bca6b",
    ),
    "owner_archive": (
        "data/enumerations/noncompact_k3_fixed_jacobian_index_blindness.json",
        "1f7bf7c64ef44043e6df6700558c8a610a30b01a8f37a3e1c3f31c49bdcfcc47",
    ),
    "owner_tests": (
        "tests/test_noncompact_k3_fixed_jacobian_index_blindness.py",
        "e736f24793b87e1865aa3c8638159209a3a3082a3fe7a0c97cb196fb00ad3408",
    ),
    "bounded_spanning_report": (
        "results/noncompact_k3_bounded_spanning_root_lemma.md",
        "de3b78a9070dc7035a8ecf51c5aa697b0736ad1559f0db8867b6629eb3cc4446",
    ),
    "bounded_spanning_independent_audit": (
        "results/noncompact_k3_bounded_spanning_root_lemma_independent_hostile_audit.md",
        "70ae2c1e09ddd0fe8cd1999087172ed14499dbd47709ebbc2236aa77ad985bde",
    ),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_digest(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def scaled_u(scale: int) -> sp.Matrix:
    return sp.Matrix([[0, scale], [scale, 0]])


def e8_cartan() -> sp.Matrix:
    matrix = 2 * sp.eye(8)
    for left, right in ((0, 2), (2, 3), (1, 3), (3, 4), (4, 5), (5, 6), (6, 7)):
        matrix[left, right] = -1
        matrix[right, left] = -1
    assert matrix.det() == 1
    return matrix


def basis_vector(size: int, position: int) -> sp.Matrix:
    vector = sp.zeros(size, 1)
    vector[position] = 1
    return vector


def minimum_exponent(polynomial: sp.Expr, variable: sp.Symbol) -> int:
    terms = sp.Poly(sp.expand(polynomial), variable).terms()
    return min(exponent[0] for exponent, coefficient in terms if coefficient != 0)


def fixed_jacobian_certificate() -> dict[str, Any]:
    """Check the two II* orders and the generic residual quartic."""

    t, s, a, b0, b1, b2 = sp.symbols("t s a b0 b1 b2", nonzero=True)
    coefficient_a = a * t**4
    coefficient_b = t**5 * (b0 + b1 * t + b2 * t**2)
    raw_discriminant = sp.expand(4 * coefficient_a**3 + 27 * coefficient_b**2)
    residual = sp.cancel(raw_discriminant / t**10)
    assert sp.Poly(residual, t).degree() == 4
    assert sp.expand(residual).subs(t, 0) == 27 * b0**2

    # Standard K3 Weierstrass rescaling at infinity: A_infty=s^8 A(1/s)
    # and B_infty=s^12 B(1/s).
    a_infinity = sp.cancel(s**8 * coefficient_a.subs(t, 1 / s))
    b_infinity = sp.cancel(s**12 * coefficient_b.subs(t, 1 / s))
    delta_infinity = sp.expand(4 * a_infinity**3 + 27 * b_infinity**2)
    assert minimum_exponent(a_infinity, s) == 4
    assert minimum_exponent(b_infinity, s) == 5
    assert minimum_exponent(delta_infinity, s) == 10

    return {
        "weierstrass_shape": "y^2=x^3+a*t^4*x+t^5*(b0+b1*t+b2*t^2)",
        "orders_at_zero": [4, 5, 10],
        "orders_at_infinity": [4, 5, 10],
        "residual_discriminant_degree": 4,
        "generic_residual_fibres": ["I1", "I1", "I1", "I1"],
        "kodaira_fibres": ["II*", "II*", "I1", "I1", "I1", "I1"],
        "euler_sum": 10 + 10 + 4,
        "trivial_lattice": "U+E8(-1)^2",
        "very_general_picard_rank": 18,
        "very_general_transcendental_lattice": "U^2",
        "shioda_tate_MW_rank": 0,
        "II_star_component_group_trivial": True,
    }


def family_lattice_certificate(index: int) -> dict[str, Any]:
    """Reconstruct T, NS, the marked fibres, and divisibility for one d."""

    if index < 2:
        raise ValueError("the hostile family is indexed by d>=2")
    d = index
    negative_e8 = -e8_cartan()
    ambient = sp.diag(scaled_u(1), scaled_u(1), scaled_u(1), negative_e8, negative_e8)

    # T=U(d)+U embedded by x=-e1+d*f2, y=e2, e3, f3.
    x = -basis_vector(22, 0) + d * basis_vector(22, 3)
    y = basis_vector(22, 2)
    transcendental = sp.Matrix.hstack(x, y, basis_vector(22, 4), basis_vector(22, 5))
    transcendental_gram = transcendental.T * ambient * transcendental
    assert transcendental_gram == sp.diag(scaled_u(d), scaled_u(1))
    transcendental_unit_minor = transcendental.extract([0, 2, 4, 5], range(4)).det()
    assert abs(transcendental_unit_minor) == 1

    # Exact orthogonal complement of the displayed T embedding.
    host_columns = [basis_vector(22, 0), d * basis_vector(22, 1) + basis_vector(22, 2)]
    host_columns.extend(basis_vector(22, coordinate) for coordinate in range(6, 22))
    host_embedding = sp.Matrix.hstack(*host_columns)
    host = host_embedding.T * ambient * host_embedding
    expected_host = sp.diag(scaled_u(d), negative_e8, negative_e8)
    assert host == expected_host
    assert transcendental.T * ambient * host_embedding == sp.zeros(4, 18)
    host_unit_minor = host_embedding.extract([0, 2, *range(6, 22)], range(18)).det()
    assert host_unit_minor == 1
    assert abs(host.det()) == d**2

    fibre = basis_vector(18, 0)
    fibre_pairing_row = list(host * fibre)
    fibre_divisibility = reduce(gcd, (abs(int(value)) for value in fibre_pairing_row))
    assert fibre_divisibility == d

    marks = sp.Matrix([2, 3, 4, 6, 5, 4, 3, 2])
    assert e8_cartan() * marks == sp.Matrix([0, 0, 0, 0, 0, 0, 0, 1])
    components: list[sp.Matrix] = []
    finite_roots: list[sp.Matrix] = []
    for block in range(2):
        simples = [basis_vector(18, 2 + 8 * block + position) for position in range(8)]
        theta = sum((int(marks[position]) * simples[position] for position in range(8)), sp.zeros(18, 1))
        affine = fibre - theta
        assert (theta.T * host * theta)[0] == -2
        assert (affine.T * host * affine)[0] == -2
        assert affine + theta == fibre
        components.extend([affine, *simples])
        finite_roots.extend(simples)

    component_matrix = sp.Matrix.hstack(*components)
    component_gram = component_matrix.T * host * component_matrix
    common_equations = component_matrix.T * host
    assert component_matrix.rank() == 17
    assert component_gram.rank() == 16
    assert common_equations.rank() == 17
    assert all((component.T * host * fibre)[0] == 0 for component in components)

    theta_first = sum(
        (int(marks[position]) * finite_roots[position] for position in range(8)),
        sp.zeros(18, 1),
    )
    formal_transverse_root = basis_vector(18, 1) - theta_first
    transverse_square = (formal_transverse_root.T * host * formal_transverse_root)[0]
    transverse_degree = (fibre.T * host * formal_transverse_root)[0]
    incident_pairings = [
        int((formal_transverse_root.T * host * component)[0]) for component in components
    ]
    assert transverse_square == -2
    assert transverse_degree == d
    assert min(incident_pairings) >= 0

    return {
        "d": d,
        "brauer_character_order": d,
        "brauer_d_torsion_rank": 4,
        "brauer_kernel_basis": ["d*e1", "f1", "e2", "f2"],
        "transcendental_gram": [
            [int(value) for value in row] for row in transcendental_gram.tolist()
        ],
        "transcendental_lattice": "U(d)+U",
        "transcendental_primitive_unit_minor": int(transcendental_unit_minor),
        "host_lattice": "U(d)+E8(-1)^2",
        "host_primitive_unit_minor": int(host_unit_minor),
        "host_is_exact_orthogonal_complement": True,
        "absolute_host_discriminant": d**2,
        "discriminant_group_invariant_factors": [d, d],
        "discriminant_length": 2,
        "fibre_divisibility": fibre_divisibility,
        "multisection_index": d,
        "affine_component_count": 18,
        "affine_vector_span_rank": component_matrix.rank(),
        "affine_gram_rank": component_gram.rank(),
        "finite_ADE_rank": 16,
        "common_orthogonal_dimension": 18 - common_equations.rank(),
        "formal_transverse_root_square": int(transverse_square),
        "formal_transverse_root_degree": int(transverse_degree),
        "formal_transverse_root_nonnegative_on_packet": min(incident_pairings) >= 0,
    }


def is_isotropic_residue(d: int, p: int, q: int, exceptional: bool) -> bool:
    numerator = -p * p + 2 * p * q if exceptional else 2 * p * q
    return numerator % (2 * d) == 0


def marked_binary_certificate(index: int) -> dict[str, Any]:
    """Exclude the even-d binary complement with diagonal residue d."""

    if index < 2:
        raise ValueError("d must be at least two")
    d = index
    standard_generators_isotropic = (
        is_isotropic_residue(d, 1, 0, False)
        and is_isotropic_residue(d, 0, 1, False)
    )
    assert standard_generators_isotropic

    exceptional_possible = d % 2 == 0
    exceptional_every_isotropic_p_even: bool | None = None
    exceptional_has_isotropic_generating_pair: bool | None = None
    if exceptional_possible:
        isotropic = [
            (p, q)
            for p in range(d)
            for q in range(d)
            if is_isotropic_residue(d, p, q, True)
        ]
        exceptional_every_isotropic_p_even = all(p % 2 == 0 for p, _ in isotropic)
        assert exceptional_every_isotropic_p_even
        # A pair generates (Z/d)^2 exactly when its determinant is a unit mod d.
        exceptional_has_isotropic_generating_pair = any(
            gcd(p * second_q - q * second_p, d) == 1
            for p, q in isotropic
            for second_p, second_q in isotropic
        )
        assert not exceptional_has_isotropic_generating_pair

    return {
        "d": d,
        "binary_gram_shape": "[[0,d],[d,2a]]",
        "smith_group_forces_d_divides_2a": True,
        "basis_shift_changes_a_by_multiple_of_d": True,
        "candidate_a_mod_d": [0] if not exceptional_possible else [0, d // 2],
        "standard_Ud_has_isotropic_generating_pair": standard_generators_isotropic,
        "exceptional_even_residue_possible_before_q_check": exceptional_possible,
        "exceptional_every_isotropic_first_coordinate_even": exceptional_every_isotropic_p_even,
        "exceptional_has_isotropic_generating_pair": exceptional_has_isotropic_generating_pair,
        "full_discriminant_form_forces_a_zero_mod_d": True,
        "actual_marked_complement": "(U(d),F=first standard isotropic vector)",
    }


def archive_certificate() -> dict[str, Any]:
    archive_path = ROOT / PINS["owner_archive"][0]
    envelope = json.loads(archive_path.read_text(encoding="utf-8"))
    assert envelope["payload_sha256"] == canonical_digest(envelope["payload"])
    family = envelope["payload"]["family"]
    assert [row["d"] for row in family] == [2, 3, 4, 5, 7, 11, 19, 31]
    for row in family:
        d = row["d"]
        assert row["absolute_discriminant"] == d**2
        assert row["primitive_fibre_divisibility"] == d
        assert row["multisection_index"] == d
        assert row["incident_vector_span_rank"] == 17
        assert row["incident_gram_rank"] == 16
        assert row["common_orthogonal_dimension"] == 1
    scope = envelope["payload"]["scope"]
    for false_claim in (
        "full_finite_volume_chamber_claimed_for_X_d",
        "full_automorphism_group_finite_claimed_for_X_d",
        "displayed_transverse_lattice_root_claimed_irreducible_physical",
        "finite_volume_counterfamily_claimed",
    ):
        assert scope[false_claim] is False
    return {
        "owner_payload_sha256": envelope["payload_sha256"],
        "sample_d_values": [row["d"] for row in family],
        "owner_archive_self_hash_valid": True,
        "negative_scope_flags_valid": True,
    }


def bounded_spanning_interface_certificate() -> dict[str, Any]:
    """Check the corrected B=13 determinant and cusp-adapted bounds."""

    rho = 18
    bound = 13
    row_norm_squared = (-2) ** 2 + (rho - 1) * bound**2
    assert row_norm_squared == 2877
    hadamard_bound = row_norm_squared ** (rho // 2)
    assert hadamard_bound == 2877**9

    highest_root_marks = [1, 2, 3, 4, 6, 5, 4, 3, 2]
    assert sum(highest_root_marks) == 30
    cusp_adapted_bound = sum(highest_root_marks) * bound
    assert cusp_adapted_bound == 390

    # The general spanning theorem gives an arbitrary bounded basis.  Its
    # determinant consequence does not numerically force the cusp-adapted 390.
    assert 391**2 <= hadamard_bound
    return {
        "rank": rho,
        "physical_pairing_bound": bound,
        "gram_diagonal": -2,
        "euclidean_row_norm_squared_bound": row_norm_squared,
        "hadamard_determinant_bound": str(hadamard_bound),
        "hadamard_formula": "abs(det G)<=2877^9",
        "overlattice_identity": "abs(det G)=[NS:R]^2*d^2",
        "automatic_finite_volume_consequence": "d^2<=2877^9",
        "selected_II_star_mark_sum": sum(highest_root_marks),
        "cusp_adapted_conditional_consequence": "d<=30*13=390",
        "arbitrary_spanning_packet_contains_selected_II_star_packet": False,
        "d_391_not_excluded_by_hadamard_bound": True,
        "d_le_390_automatic_from_B13": False,
    }


def build_audit() -> dict[str, Any]:
    pins = {
        name: {
            "path": relative,
            "expected_sha256": expected,
            "observed_sha256": sha256(ROOT / relative),
        }
        for name, (relative, expected) in PINS.items()
    }
    for record in pins.values():
        assert record["observed_sha256"] == record["expected_sha256"]
    sample_d = [2, 3, 4, 6, 7, 10, 15, 31]
    return {
        "schema": "coxeter.noncompact-k3-fixed-jacobian-index-independent-hostile.v1",
        "verdict": "PASS_LOCAL_INDEX_BLINDNESS_EXACT_SCOPE",
        "pins": pins,
        "fixed_jacobian": fixed_jacobian_certificate(),
        "family_lattices": [family_lattice_certificate(d) for d in sample_d],
        "marked_binary_complements": [marked_binary_certificate(d) for d in sample_d],
        "owner_archive": archive_certificate(),
        "bounded_spanning_interface": bounded_spanning_interface_certificate(),
        "source_theorems": {
            "projective_torsor_existence": "Br(S)=Sha(S), with every class represented in the projective elliptic-K3 category",
            "period_index": "multisection index equals order of the Brauer class",
            "transcendental_kernel": "T(X_alpha) is Hodge-isometric to ker(alpha:T(S)->Z/d)",
            "fibrewise_Jacobian_equality": "J^k(X)_t is isomorphic to X_t for every closed t",
            "primitive_embedding_uniqueness": "Nikulin stable range applies because ell(A_T)+2=4<=18",
            "machine_proved_here": False,
        },
        "scope": {
            "same_fixed_relative_jacobian": True,
            "projective_torsors_of_every_order": True,
            "same_two_II_star_and_four_I1_fibres": True,
            "same_full_local_parabolic_packet": True,
            "unbounded_fibre_divisibility_without_global_volume": True,
            "finite_volume_claimed_for_family": False,
            "reflection_maximality_claimed_for_family": False,
            "formal_transverse_root_claimed_physical": False,
            "automatic_d_le_390_claimed": False,
            "conditional_d_squared_le_2877_to_the_9th_under_finite_volume": True,
            "Property_4_proved": False,
        },
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-pins", action="store_true")
    parser.parse_args(argv)
    print(json.dumps(build_audit(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
