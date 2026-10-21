from __future__ import annotations

from functools import lru_cache
import json

import sympy as sp

from scripts.certify_h9_14_t1_strong_k3_nef_obstruction import (
    OUTPUT,
    build_certificate,
)
from scripts.probe_h9_14_t1_strong_k3_nef_promotion import reconstruct


@lru_cache(maxsize=1)
def _payload() -> dict:
    return build_certificate()


def test_archive_replays_and_scope_is_exact() -> None:
    payload = _payload()
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["all_checks_pass"] is True
    assert all(payload["checks"].values())
    assert payload["scope"]["weak_K3_lattice_action_carrier"] is True
    assert payload["scope"]["strong_invariant_nef_realization"] is False
    assert payload["scope"]["all_K3_realizations_of_the_abstract_T1_polytope_excluded"] is False


def test_first_hidden_pair_and_crossing_rays_literal_equalities() -> None:
    payload = _payload()
    package = reconstruct()
    gram = package["n_gram"]
    action = package["involution"]
    fixed_gram = package["fixed_gram"]
    fixed_roots = package["fixed_roots"]
    fixed_embedding = package["fixed_embedding"]
    witness = payload["hidden_root_witnesses"][0]
    fixed = sp.Matrix(witness["fixed_normal"])
    plus = sp.Matrix(witness["lifted_root_plus"])
    minus = sp.Matrix(witness["lifted_root_minus"])

    assert (fixed.T * fixed_gram * fixed)[0] == -4
    assert (plus.T * gram * plus)[0] == -2
    assert (minus.T * gram * minus)[0] == -2
    assert (plus.T * gram * minus)[0] == 0
    assert action * plus == minus
    assert action * minus == plus
    assert plus + minus == fixed_embedding * fixed

    pairings = []
    for key in ("negative_extreme_ray", "positive_extreme_ray"):
        record = witness[key]
        ray = sp.Matrix(record["coordinates"])
        slacks = fixed_roots.T * fixed_gram * ray
        active = [index for index, value in enumerate(slacks) if value == 0]
        assert list(map(int, slacks)) == record["facet_slacks"]
        assert all(value >= 0 for value in slacks)
        assert fixed_roots[:, active].rank() == 9
        assert (ray.T * fixed_gram * ray)[0] == 10
        pairings.append(int((fixed.T * fixed_gram * ray)[0]))
        lifted = fixed_embedding * ray
        assert (lifted.T * gram * plus)[0] == pairings[-1] // 2
        assert (lifted.T * gram * minus)[0] == pairings[-1] // 2
    assert pairings == [-6, 6]


def test_all_five_hidden_pairs_cross_and_source_bridge_really_fails() -> None:
    payload = _payload()
    witnesses = payload["hidden_root_witnesses"]
    assert len(witnesses) == 5
    assert all(all(row["checks"].values()) for row in witnesses)
    assert all(
        row["negative_extreme_ray"]["hidden_normal_pairing"] < 0
        < row["positive_extreme_ray"]["hidden_normal_pairing"]
        for row in witnesses
    )

    boundary = payload["marked_bridge_boundary"]
    assert boundary["labelled_14_facet_normalized_Coxeter_gram_identity"] is True
    assert boundary["literal_equivariant_24_source_wall_bridge"] is False
    assert boundary["factory_source_package"]["ambient_rank"] == 14
    assert boundary["compiled_source_package"]["ambient_rank"] == 18
    assert boundary["factory_source_package"]["source_action_orbit_histogram"] == {
        "1": 8,
        "2": 8,
    }
    assert boundary["compiled_source_package"]["source_action_orbit_histogram"] == {
        "1": 4,
        "2": 10,
    }


def test_weak_period_gate_survives_but_strong_C_is_negative() -> None:
    payload = _payload()
    weak = payload["weak_K3_period_and_ample_gate"]
    assert weak["complement"] == "U+U(2)"
    assert weak["period_square"] == "0"
    assert weak["period_positive_norm"] == "8*sqrt(6)"
    assert weak["rational_orthogonal_is_zero"] is True
    assert weak["rootless_coinvariant_E8_minus2"] is True
    assert weak["genuine_symplectic_Nikulin_automorphism_exists"] is True
    assert weak["strict_compiled_interior"]["square"] == 680
    assert weak["strict_compiled_interior"]["strict_inside_all_14_compiled_facets"] is True

    obstruction = payload["effectivity_and_cone_obstruction"]
    assert obstruction["irreducibility_needed"] is False
    assert obstruction["candidate_cone_not_contained_in_any_invariant_nef_section"] is True
    assert obstruction["candidate_cone_cannot_equal_an_invariant_nef_section"] is True
    assert payload["URCV_boundary"]["C"].startswith("closed negatively")
