from __future__ import annotations

from functools import lru_cache
import json

from scripts.certify_bugaenko_h8_determinant_line_ag import OUTPUT, analyze


def _payload() -> dict[str, object]:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


@lru_cache(maxsize=1)
def _replay() -> dict[str, object]:
    return analyze()


def test_certificate_replays_exactly_from_archived_polytope_data() -> None:
    payload = _payload()
    assert _replay() == payload
    assert payload["all_checks_pass"] is True
    assert payload["source"] == {
        "compact_vertex_count": 41,
        "cusp_count": 0,
        "dimension": 8,
        "f_vector": [41, 164, 316, 374, 294, 156, 54, 11, 1],
        "facet_count": 11,
        "only_polytope_lattice_input": True,
        "path": "data/known_examples/bugaenko_compact_h8.json",
        "sha256": "c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c",
        "vinberg_completeness_status": "verified",
    }


def test_canonical_determinant_line_has_the_required_two_place_signatures() -> None:
    payload = _payload()
    assert payload["normalization"] == {
        "determinant": "2*phi",
        "even_integral": True,
        "lattice": "L=<2*phi> + E8(-1) over O=Z[phi]",
        "operation": "multiply the archived bilinear form by -1",
        "rank": 9,
        "signatures": {
            "identity": {"positive": 1, "negative": 8, "zero": 0},
            "conjugate": {"positive": 0, "negative": 9, "zero": 0},
        },
    }
    line = payload["determinant_line"]
    assert line["coordinate_free_module"] == "K_L=det_O(L)"
    assert line["generator_square"] == "4*phi"
    assert (line["sign_identity_embedding"], line["sign_conjugate_embedding"]) == (
        1,
        -1,
    )
    extension = payload["canonical_extension"]
    assert extension["signatures"] == {
        "identity": {"positive": 2, "negative": 8, "zero": 0},
        "conjugate": {"positive": 0, "negative": 10, "zero": 0},
    }
    assert extension["determinant_ideal"] == "(8*phi^2)=(8)"
    assert extension["absolute_discriminant_module_order"] == 64
    assert extension["anisotropic_over_Q_sqrt_5"] is True


def test_all_even_discriminant_gluings_are_exhausted() -> None:
    gluing = _payload()["discriminant_gluing"]
    assert gluing["element_count"] == 64
    assert gluing["underlying_abelian_group"] == "(Z/2)^2 + (Z/4)^2"
    assert gluing["quadratic_module"] == {
        "module": "O/(2) * u + O/(4) * v",
        "q_u": "1/(4*phi) mod O",
        "q_v": "1/(8*phi) mod O",
        "formula": "q(a*u+b*v)=(2*a^2+b^2)/(8*phi) mod O",
    }
    assert gluing["isotropic_elements"] == [
        {"a_mod_2": [0, 0], "b_mod_4": [0, 0]}
    ]
    assert gluing["isotropic_element_count"] == 1
    assert gluing["totally_isotropic_O_submodule_count"] == 1
    assert gluing["reflection_stable_even_overlattice_count"] == 1
    assert gluing["reflection_stable_even_overlattices"] == ["Lhat_0"]
    assert gluing["enumeration_sha256"] == (
        "bbb3975b9a818fe3b4a5bbecb88050aade28009a330aa8450282e7fe64356047"
    )


def test_all_eleven_roots_extend_and_have_one_heegner_type() -> None:
    payload = _payload()
    roots = payload["root_extensions"]
    assert roots["root_count"] == 11
    assert roots["wall_gram_matches_normalized_source"] is True
    assert roots["all_primitive_square_minus_two_divisibility_one"] is True
    assert roots["all_reflections_integral_involutive_isometries"] is True
    assert roots["discriminant_action"] == "identity"
    assert len(roots["reflection_matrix_sha256"]) == 11
    assert all(record["square_in_Lhat"] == "-2" for record in roots["roots"])
    assert all(record["divisibility_ideal"] == "O" for record in roots["roots"])

    heegner = payload["facet_heegner_types"]
    assert heegner == {
        "all_same_numerical_type": True,
        "bruinier_sign_type": {"m": 1, "mu": "0"},
        "facet_count": 11,
        "global_orbit_or_irreducible_component_equality_claimed": False,
        "signature_2_8_type": {"coset": "0", "divisibility": 1, "norm": "-2"},
    }


def test_stack_recovers_every_ridge_label_and_excludes_ultraparallel_pair() -> None:
    stack = _payload()["stacky_coxeter_labels"]
    assert stack["ridge_count"] == 54
    assert stack["finite_pair_equals_ridge_set"] is True
    assert stack["distinguished_dihedral_inertia_counts"] == {
        "I2(2)_order_4": 44,
        "I2(3)_order_6": 8,
        "I2(5)_order_10": 2,
    }
    assert stack["ultraparallel_pairs"] == [[9, 10]]
    assert sum(record["is_ridge"] for record in stack["pairs"]) == 54


def test_hidden_sheet_gate_has_exact_11_309_0_split() -> None:
    payload = _payload()
    hidden = payload["hidden_slice_gate"]
    assert hidden["classification"] == {
        "facet_mirrors": 11,
        "nonfacet_lower_boundary_contacts": 309,
        "interior_cutting_mirrors": 0,
    }
    assert hidden["distinct_mirror_count"] == 320
    assert hidden["nonfacet_contact_dimension_distribution"] == {
        "0": 88,
        "1": 43,
        "2": 26,
        "3": 17,
        "4": 96,
        "5": 25,
        "6": 14,
    }
    assert hidden["mirror_representative_sha256"] == (
        "815c0c4cad8e4ab9da32ef491c104ac71a781d4f2b4300f988ecb484506da69e"
    )
    assert hidden["closed_slice_is_facet_only"] is False
    assert hidden["open_chamber_interior_is_uncut_by_extra_shell_mirrors"] is True
    archive = payload["mirror_archive_cross_check"]
    assert archive["source_only_replay_matches_archive_exactly"] is True
    assert archive["sha256"] == (
        "e77bb81bcc738f4f542676748dbcf09418cab830da579a600fd1e309f03de782"
    )


def test_transfer_claims_are_kept_separate_from_missing_principality() -> None:
    payload = _payload()
    real_slice = payload["canonical_real_slice"]
    assert real_slice["facet_mirror_pullback"] == (
        "j_e^{-1}(D_(r_i,0))=r_i^perp"
    )
    assert real_slice["proper_projective_coarse_quotient"] is True
    assert real_slice["surface_is_canonical_without_auxiliary_choices"] is False

    conclusion = payload["AG_arrangement_conclusion"]
    assert conclusion["associated_projective_AG_object_exists"] is True
    assert conclusion["associated_surface_exists_after_neat_cover_and_complete_intersection"] is True
    assert conclusion["Coxeter_angles_recovered_stack_theoretically"] is True
    assert conclusion["full_numerical_Heegner_type_is_facet_only_on_closed_P"] is False
    assert conclusion["full_numerical_Heegner_type_cuts_open_P"] is False
    assert conclusion["automorphic_principality_or_Borcherds_product_claimed"] is False
    assert conclusion["new_abstract_Coxeter_polytope_claimed"] is False
