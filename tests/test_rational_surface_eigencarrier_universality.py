import json
from pathlib import Path

import pytest
import sympy as sp

from coxeter_polytopes.geometric.rational_surface_carrier import (
    embed_lorentzian_lattice,
    extend_self_adjoint_operator,
)
from scripts.certify_rational_surface_eigencarrier_universality import analyze


ROOT = Path(__file__).resolve().parents[1]


def test_generic_lorentzian_embedding_with_isotropic_basis() -> None:
    gram = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, -2]])
    result = embed_lorentzian_lattice(gram)
    embedding = sp.Matrix(result.embedding)
    ambient = sp.Matrix(result.ambient_gram)
    assert result.blowup_count == 9
    assert embedding.rank() == 3
    assert embedding.T * ambient * embedding == result.vector_scale**2 * gram


@pytest.mark.parametrize(
    "gram",
    [
        sp.Matrix([[7]]),
        sp.Matrix([[0, 1], [1, 0]]),
        sp.Matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, -2, 1], [0, 0, 1, -3]]),
        sp.Matrix([[9, 6, 3], [6, 1, 2], [3, 2, -2]]),
    ],
)
def test_adversarial_lorentzian_embeddings(gram: sp.Matrix) -> None:
    result = embed_lorentzian_lattice(gram)
    embedding = sp.Matrix(result.embedding)
    ambient = sp.Matrix(result.ambient_gram)
    assert embedding.rank() == gram.rows
    assert embedding.T * ambient * embedding == result.vector_scale**2 * gram
    source_gram = embedding.T * ambient * embedding
    projector = embedding * source_gram.inv() * embedding.T * ambient
    assert projector**2 == projector
    assert projector.T * ambient == ambient * projector
    complement = sp.Matrix.hstack(*(embedding.T * ambient).nullspace())
    complement_gram = complement.T * ambient * complement
    assert complement_gram.is_negative_definite is True


def test_embedding_rejects_forms_outside_exact_scope() -> None:
    rejected = (
        sp.Matrix([[1, 2], [0, -1]]),
        sp.Matrix([[1, 0], [0, 0]]),
        sp.diag(1, 1, -1),
        sp.Matrix([[sp.Rational(1, 2), 0], [0, -1]]),
    )
    for gram in rejected:
        with pytest.raises(ValueError):
            embed_lorentzian_lattice(gram)


def test_self_adjoint_extension_rejects_bad_ambient_and_operator() -> None:
    embedding = sp.eye(2)
    ambient = sp.diag(1, -1)
    with pytest.raises(ValueError):
        extend_self_adjoint_operator(
            embedding, sp.Matrix([[1, 1], [0, -1]]), sp.eye(2)
        )
    with pytest.raises(ValueError):
        extend_self_adjoint_operator(
            embedding, ambient, sp.Matrix([[0, 1], [0, 0]])
        )


def test_compact_bugaenko_rational_surface_carrier_certificate() -> None:
    result = analyze()
    frozen = json.loads(
        (ROOT / "data/enumerations/rational_surface_eigencarrier_universality.json").read_text()
    )
    assert result == frozen
    compact = result["compact_calibration"]
    assert compact["vinberg_completion_verified"] is True
    assert compact["finite_vertex_count"] == 41
    assert compact["cusp_count"] == 0
    assert compact["hyperideal_ray_count"] == 0
    assert compact["embedding_identity_verified"] is True
    assert compact["facet_divisors_integral"] is True
    assert compact["facet_divisor_square_set"] == [-8]
    assert compact["facet_pairing_transport_verified"] is True
    assert compact["field_root_trace_pairing_verified"] is True
    assert compact["identity_eigenspace_projected_wall_gram_verified"] is True
    assert compact["oriented_identity_component_verified"] is True
    assert compact["divisor_product_cycle_action_verified"] is True
    assert compact["divisor_product_cycle_transpose_invariant"] is True
    assert compact["operator_self_adjoint_verified"] is True
    assert compact["operator_characteristic_polynomial_verified"] is True
    assert compact[
        "trace_chow_projector_idempotent_under_cycle_composition"
    ] is True
    assert compact[
        "golden_chow_action_relations_under_cycle_composition"
    ] is True
    assert compact["rational_chow_action_cycles_transpose_invariant"] is True
    assert compact["integral_generator_is_scaled_rational_chow_action"] is True
    assert compact["trace_chow_projector_rank"] == 18
    assert compact["real_eigenspace_signatures_verified"] is True
    assert compact["points_distinct_verified"] is True
    assert compact["section_is_exactly_bugaenko_h8"] is True
