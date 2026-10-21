from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_bugaenko_hodge_coherent_carrier import analyze


ROOT = Path(__file__).resolve().parents[1]


def test_exact_hodge_coherent_complement_identities() -> None:
    result = analyze()
    assert result["sqrt5_square_is_5I"] is True
    assert result["sqrt5_self_adjoint"] is True
    assert result["tau_polynomial_is_tau_plus_I"] is True
    assert result["tau_self_adjoint"] is True
    assert result["sqrt5_twisted_leading_minors"] == [6, 20, 20, 100]
    assert result["sqrt5_twisted_gram_positive_definite"] is True
    assert result["eigenspace_signatures"] == {
        "+sqrt5": [2, 0],
        "-sqrt5": [0, 2],
    }
    assert result["combined_N_plus_T_basis_index_in_k3_lattice"] == 4
    assert result["sqrt5_discriminant_action_on_N_is_identity"] is True
    assert result["sqrt5_discriminant_action_on_T_is_identity"] is True
    assert result["full_k3_sqrt5_action_integral"] is True
    assert result["full_k3_sqrt5_square_is_5I"] is True
    assert result["full_k3_sqrt5_self_adjoint"] is True
    assert result["sqrt5_extends_over_unimodular_glue"] is True
    assert result["full_k3_tau_action_integral"] is False
    assert result["full_k3_tau_polynomial_is_tau_plus_I"] is True
    assert result["full_k3_tau_self_adjoint"] is True
    assert result["combined_basis_gram_is_N_orthogonal_sum_T"] is True
    assert result["period_Bu_is_sqrt5_u"] is True
    assert result["period_Bv_is_sqrt5_v"] is True
    assert result["period_Ju_is_minus_sqrt20_v"] is True
    assert result["period_Jv_is_sqrt20_u"] is True
    assert result["period_omega_square_is_zero"] is True
    assert result["period_omega_dot_bar_omega"] == [40, 24]
    assert result["period_omega_dot_bar_omega_positive"] is True
    assert result["cm_i_action_square_is_minus_I"] is True
    assert result["buskin_cayley_z_is_T_isometry"] is True
    assert result["buskin_cayley_z_inverse_is_T_isometry"] is True
    assert result["buskin_sqrt5_span_identity"] is True


def test_frozen_certificate_matches_recomputation() -> None:
    frozen = json.loads(
        (ROOT / "data/candidates/bugaenko_hodge_coherent_carrier.json").read_text()
    )
    assert frozen == analyze()
    assert frozen["rank20_k3_h9_hodge_coherent_gate"] == "empty"
    assert frozen["k3_hilbert_scheme_full_h2_gate"] == "empty_in_every_dimension"
