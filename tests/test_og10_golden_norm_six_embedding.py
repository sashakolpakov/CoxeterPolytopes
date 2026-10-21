import hashlib
import json

from scripts.certify_og10_golden_norm_six_embedding import (
    CERTIFICATE_PATH,
    build_certificate,
    verify_certificate,
)


EXPECTED_CERTIFICATE_FILE_SHA256 = (
    "e9b8c240e8a04def285e48486d7cb3daf7375b97e6159a7f1cf1f1b1d8a78a56"
)
EXPECTED_CERTIFICATE_PAYLOAD_SHA256 = (
    "2042bd2cbcd9d003ed226a20ebc7631829a7139f951a32a9a1e009e500c666be"
)


def test_norm_six_embedding_recomputes_exactly() -> None:
    archived = json.loads(CERTIFICATE_PATH.read_text())
    assert hashlib.sha256(CERTIFICATE_PATH.read_bytes()).hexdigest() == (
        EXPECTED_CERTIFICATE_FILE_SHA256
    )
    assert archived["certificate_payload_sha256"] == (
        EXPECTED_CERTIFICATE_PAYLOAD_SHA256
    )
    assert archived == build_certificate()
    assert all(verify_certificate(archived).values())
    assert all(archived["checks"].values())


def test_norm_six_discriminant_gluing_is_the_selected_order_twelve_graph() -> None:
    payload = json.loads(CERTIFICATE_PATH.read_text())
    graph = payload["discriminant_form_gluing"]

    assert payload["field_carrier"] == {
        "definition": "H6(k)=[[F,3k],[3k^t,6]], F=<-2phi>+E8+<2>",
        "field_determinant": "-6*phi",
        "field_rank": 11,
        "field_signatures": {"conjugate": [11, 0], "identity": [10, 1]},
        "normalized_covector": "k=(0^9,1)",
    }
    assert payload["trace_lattice"]["signature"] == [1, 21]
    assert payload["trace_lattice"]["determinant"] == -36
    assert payload["trace_lattice"]["smith_diagonal"] == [1] * 20 + [6, 6]

    assert graph["n_discriminant_group"] == "Z/6 + Z/6"
    assert graph["t_discriminant_group"] == "Z/2 + Z/6"
    assert graph["all_injective_anti_isometry_generator_pairs"] == 4
    assert graph["pairs_preserving_selected_D_over_3"] == 2
    assert graph["selected_pairs_form_one_complement_sign_orbit"] is True
    assert graph["chosen_generator_orders"] == [2, 6]
    assert graph["chosen_generator_squares"] == [2, 0]
    assert graph["chosen_generator_pairing"] == 1
    assert graph["chosen_graph_image_orders"] == {
        "N": 12,
        "N_plus_T": 12,
        "T": 12,
    }
    assert graph["chosen_graph_isotropic"] is True
    assert graph["selected_D_over_3_square_mod_2Z"] == "4/3"
    assert graph["selected_D_over_3_orthogonal_to_graph"] is True

    assert payload["gluing"]["complement"] == "<2>+<6>"
    assert payload["gluing"]["overlattice_index"] == 12
    assert payload["gluing"]["ambient_determinant"] == -3
    assert payload["gluing"]["ambient_smith_diagonal"] == [1] * 23 + [3]


def test_embedding_primitivity_selected_wall_and_literal_isometry_are_frozen() -> None:
    payload = json.loads(CERTIFICATE_PATH.read_text())

    embeddings = payload["primitive_embeddings"]
    assert embeddings["N_embedding_smith_diagonal"] == [1] * 22
    assert embeddings["T_embedding_smith_diagonal"] == [1, 1]
    assert embeddings["orthogonal_complement_gram"] == [[2, 0], [0, 6]]

    selected = payload["selected_trace_lifts"]
    assert selected["field_reflection"]["formula"] == (
        "s_e(x)=x-H6(x,e)*e/3"
    )
    d_class = selected["D_equals_e"]
    e_class = selected["E_equals_phi_minus_1_times_e"]
    assert (d_class["square"], d_class["ambient_divisibility"]) == (-6, 3)
    assert d_class["primitive"] is True
    assert d_class["discriminant_value"] == "-2/3 = 4/3 mod 2Z"
    assert (e_class["square"], e_class["ambient_divisibility"]) == (-6, 1)
    assert e_class["primitive"] is True

    matrices = payload["matrices"]
    assert len(matrices["ambient_gram_24"]) == 24
    assert len(matrices["ambient_to_standard_og10_isometry_24"]) == 24
    assert len(matrices["primitive_N_embedding_24x22"][0]) == 22
    assert len(matrices["primitive_T_embedding_24x2"][0]) == 2
    assert payload["matrix_sha256"] == {
        "ambient_gram_24": (
            "44f50ba5b4bce1cab90ec864b6d35c35888ae1f7476d2a969a4518a0dd1959c0"
        ),
        "ambient_to_standard_og10_isometry_24": (
            "a872feb3962af243d90e15f33b79a1ff97a48023752504370f4bb194f3b46ab9"
        ),
        "combined_standard_basis_in_N_plus_T_24": (
            "06a1851b4568d28eea1c5aefdbe566307d27c609e4170c8b5cc3949e13984c76"
        ),
        "full_gluing_basis_in_N_plus_T_24": (
            "8fa0efc379c4213597e6780ee2186635b867baa3862bf3980d3580f8e9a32631"
        ),
        "primitive_N_embedding_24x22": (
            "f137c784f38aabf08c5398ef60744fd85fd788904f4a4c74a1793f577c1021e1"
        ),
        "primitive_T_embedding_24x2": (
            "2fd141fed3ccf6a1bbe2c1d653838bf3dfa95b9542d80b0c90a70a305e13d1f7"
        ),
        "trace_gram_22": (
            "26eb861fae5fc3e5af609be8cb6a128cfb63e9e1f475ea01193399baaa2a6cc0"
        ),
    }

    geometry = payload["period_projectivity_and_MBM"]
    assert geometry["selected_D_is_an_OG10_MBM_wall_class"] is True
    assert geometry["positive_class"]["square"] == 2
    assert "full Picard lattice" in geometry["full_NS_statement"]
    assert payload["scope"].startswith(
        "This is a positive algebro-geometric carrier and wall theorem."
    )
