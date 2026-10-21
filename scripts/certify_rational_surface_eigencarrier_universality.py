"""Certify the rational-surface eigen-carrier theorem and its compact H8 model.

The general theorem is proved in the accompanying report.  This script checks
the constructive lattice lemma, then gives an unusually small explicit model
for Bugaenko's compact chamber: the blow-up of P^2 at 19 rational points.  Its
Picard lattice is I_(1,19), the golden trace lattice embeds after scale four,
and an integral divisor-product correspondence has the required golden
eigenspace.
"""

from __future__ import annotations

import json
from fractions import Fraction

import sympy as sp

from coxeter_polytopes.exact.matrix import inertia
from coxeter_polytopes.exact.trace_lattice import integral_gram, trace_gram
from coxeter_polytopes.geometric.k3 import K3_TRACE_SCALAR
from coxeter_polytopes.geometric.rational_surface_carrier import (
    embed_lorentzian_lattice,
    extend_self_adjoint_operator,
)
from coxeter_polytopes.known_examples.bugaenko import (
    bugaenko_h8_certificate,
    bugaenko_interior_point,
    bugaenko_simple_roots,
    bugaenko_wall_gram,
)
from coxeter_polytopes.verify.vinberg import verify_bugaenko_vinberg_completion
from scripts.bugaenko_k3_eigen_nef import multiplication_by_phi
from scripts.bugaenko_trace_lattice import (
    basis_change_certificate,
    bugaenko_field_gram,
    expected_standard_gram,
)
from coxeter_polytopes.search.bugaenko_vinberg import (
    build_bugaenko_vinberg_transcript,
)


def scaled_e8_coordinate_matrix() -> sp.Matrix:
    """Columns have Euclidean Gram ``4*E8`` in the repository ordering."""

    roots = [
        [1, -1, -1, -1, -1, -1, -1, 1],
        [2, 2, 0, 0, 0, 0, 0, 0],
        [-2, 2, 0, 0, 0, 0, 0, 0],
        [0, -2, 2, 0, 0, 0, 0, 0],
        [0, 0, -2, 2, 0, 0, 0, 0],
        [0, 0, 0, -2, 2, 0, 0, 0],
        [0, 0, 0, 0, -2, 2, 0, 0],
        [0, 0, 0, 0, 0, -2, 2, 0],
    ]
    # The displayed roots have arm lengths 2,1,4 at vertex 3.  This
    # permutation moves the branch to vertex 2, matching e8_cartan_gram().
    repository_order = [0, 2, 3, 4, 5, 6, 7, 1]
    return sp.Matrix([roots[index] for index in repository_order]).T


def bugaenko_surface_embedding() -> tuple[sp.Matrix, sp.Matrix]:
    """Embed ``4*(<2>+<-2>+2E8(-1))`` in ``I_(1,19)``."""

    embedding = sp.zeros(20, 18)
    # The first vector has square 3^2-1^2=8=4*2.
    embedding[0, 0] = 3
    embedding[1, 0] = 1
    # The second has square -(2^2+2^2)=-8=4*(-2).
    embedding[2, 1] = 2
    embedding[3, 1] = 2
    e8 = scaled_e8_coordinate_matrix()
    for row in range(8):
        for column in range(8):
            embedding[4 + row, 2 + column] = e8[row, column]
            embedding[12 + row, 10 + column] = e8[row, column]
    ambient = sp.diag(1, *([-1] * 19))
    return embedding, ambient


def _trace_coordinate_vector(field_vector: list[object]) -> sp.Matrix:
    """Put an integral Bottinelli vector in the trace script's E8 ordering."""

    # The two E8 bases have the same arms of lengths 2,1,4 but place the
    # branch at different numbered vertices.  Rows are the trace-script E8
    # coordinates and columns are the archived Bugaenko coordinates.
    repo_to_bottinelli = [0, 2, 3, 4, 5, 6, 7, 1]
    reordered = [field_vector[0]] + [
        field_vector[1 + bottinelli_index]
        for bottinelli_index in repo_to_bottinelli
    ]

    coordinates: list[int] = []
    for entry in reordered:
        if entry.a.denominator != 1 or entry.b.denominator != 1:
            raise AssertionError("a Bugaenko vector is not integral over Z[phi]")
        coordinates.extend([entry.a.numerator, entry.b.numerator])
    return sp.Matrix(coordinates)


def _identity_value(value: object) -> sp.Expr:
    """Evaluate a repository golden number at ``phi=(1+sqrt(5))/2``."""

    phi_real = (1 + sp.sqrt(5)) / 2
    return sp.Rational(value.a.numerator, value.a.denominator) + sp.Rational(
        value.b.numerator, value.b.denominator
    ) * phi_real


def _cycle_terms(operator: sp.Matrix, ambient: sp.Matrix) -> list[tuple[int, int, int]]:
    """Return ``(coefficient,left,right)`` for divisor products.

    In the basis ``H,E_1,...,E_b``, the intersection-dual of basis vector
    ``j`` is ``ambient[j,j]`` times that vector.  Therefore these products
    induce exactly the supplied operator on NS.
    """

    terms: list[tuple[int, int, int]] = []
    for source in range(operator.cols):
        dual_sign = int(ambient[source, source])
        for target in range(operator.rows):
            coefficient = dual_sign * int(operator[target, source])
            if coefficient:
                terms.append((coefficient, source, target))
    return terms


def _operator_from_cycle_terms(
    terms: list[tuple[int, int, int]], ambient: sp.Matrix
) -> sp.Matrix:
    operator = sp.zeros(ambient.rows)
    for coefficient, source, target in terms:
        operator[target, source] += coefficient * ambient[source, source]
    return operator


def _divisor_cycle_coefficients(
    operator: sp.Matrix, intersection_gram: sp.Matrix
) -> sp.Matrix:
    """Coefficient matrix of ``sum e_j^* x operator(e_j)``.

    If ``C_A`` is this matrix, exact composition of divisor-product cycles is
    ``C_A * G * C_B = C_(B*A)``.  Thus matrix identities below certify Chow
    correspondence identities, not merely identities of their cohomological
    actions.
    """

    return intersection_gram.inv() * operator.T


def analyze() -> dict[str, object]:
    # Exercise the universal algorithm on a basis beginning with an isotropic
    # hyperbolic plane; this specifically checks the non-diagonal pivot path.
    calibration_gram = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, -2]])
    generic = embed_lorentzian_lattice(calibration_gram)
    generic_identity = (
        sp.Matrix(generic.embedding).T
        * sp.Matrix(generic.ambient_gram)
        * sp.Matrix(generic.embedding)
        == generic.vector_scale**2 * calibration_gram
    )

    standard_gram = sp.Matrix(expected_standard_gram())
    embedding, ambient = bugaenko_surface_embedding()
    embedding_identity = embedding.T * ambient * embedding == 4 * standard_gram
    if not embedding_identity:
        raise AssertionError("the optimized I_(1,19) embedding failed")

    trace_change = sp.Matrix(basis_change_certificate())
    trace_action = sp.Matrix(multiplication_by_phi())
    standard_action = trace_change.inv() * trace_action * trace_change
    if standard_action.T * standard_gram != standard_gram * standard_action:
        raise AssertionError("golden multiplication is not self-adjoint")
    if standard_action**2 != standard_action + sp.eye(18):
        raise AssertionError("golden multiplication has the wrong polynomial")

    field_gram = bugaenko_field_gram()
    field_inertias = [inertia(field_gram, embedding=index) for index in (0, 1)]
    trace_eigenspace_signatures = []
    for index, field_inertia in enumerate(field_inertias):
        if K3_TRACE_SCALAR.sign(index) < 0:
            trace_eigenspace_signatures.append(
                [field_inertia.negative, field_inertia.positive]
            )
        else:  # pragma: no cover - both Bugaenko scalar embeddings are negative
            trace_eigenspace_signatures.append(
                [field_inertia.positive, field_inertia.negative]
            )

    extension = extend_self_adjoint_operator(
        embedding,
        ambient,
        standard_action,
        complement_scalar=0,
    )
    ambient_operator = sp.Matrix(extension.operator)
    cycle_terms = _cycle_terms(ambient_operator, ambient)
    cycle_verified = _operator_from_cycle_terms(cycle_terms, ambient) == ambient_operator
    cycle_coefficient_matrix = ambient * ambient_operator.T
    cycle_transpose_invariant = (
        cycle_coefficient_matrix == cycle_coefficient_matrix.T
    )
    spectral_variable = sp.Symbol("lambda")
    operator_characteristic_polynomial_verified = sp.expand(
        ambient_operator.charpoly(spectral_variable).as_expr()
        - spectral_variable**2
        * (spectral_variable**2 - 4 * spectral_variable - 16) ** 9
    ) == 0

    # The rational orthogonal projector is itself a divisor-product Chow
    # idempotent.  Zero on the complement extends the whole Q(phi)-action to
    # its image: the projector is the unit on that summand.  Verify all
    # identities using literal cycle coefficient composition C_A*J*C_B,
    # rather than only comparing induced actions.
    trace_projector = sp.Matrix(extension.rational_projection)
    embedded_source_gram = embedding.T * ambient * embedding
    rational_phi_extension = (
        embedding
        * standard_action
        * embedded_source_gram.inv()
        * embedding.T
        * ambient
    )
    projector_cycle = _divisor_cycle_coefficients(trace_projector, ambient)
    phi_cycle = _divisor_cycle_coefficients(rational_phi_extension, ambient)
    chow_projector_idempotent = (
        projector_cycle * ambient * projector_cycle == projector_cycle
    )
    chow_field_action_relations = (
        projector_cycle * ambient * phi_cycle == phi_cycle
        and phi_cycle * ambient * projector_cycle == phi_cycle
        and phi_cycle * ambient * phi_cycle == phi_cycle + projector_cycle
    )
    chow_cycles_transpose_invariant = (
        projector_cycle == projector_cycle.T and phi_cycle == phi_cycle.T
    )
    integral_generator_is_scaled_rational_action = (
        ambient_operator
        == extension.operator_multiplier * rational_phi_extension
    )
    complement_basis = sp.Matrix.hstack(*(embedding.T * ambient).nullspace())
    complement_gram = complement_basis.T * ambient * complement_basis
    complement_inertia = inertia(
        [
            [str(complement_gram[row, column]) for column in range(complement_gram.cols)]
            for row in range(complement_gram.rows)
        ]
    )
    eigenspace_signatures_verified = (
        trace_eigenspace_signatures == [[1, 8], [0, 9]]
        and (complement_inertia.positive, complement_inertia.negative)
        == (0, 2)
    )

    trace_gram_matrix = sp.Matrix(
        integral_gram(trace_gram(bugaenko_field_gram(), K3_TRACE_SCALAR))
    )
    if trace_change.T * trace_gram_matrix * trace_change != standard_gram:
        raise AssertionError("the pinned trace-to-standard marking failed")
    inverse_change = trace_change.inv()
    trace_walls = [
        _trace_coordinate_vector(root) for root in bugaenko_simple_roots()
    ]
    standard_walls = [inverse_change * wall for wall in trace_walls]
    divisor_walls = [embedding * wall for wall in standard_walls]
    wall_integrality = all(
        entry.is_Integer for wall in divisor_walls for entry in wall
    )
    divisor_wall_gram = sp.Matrix.hstack(*divisor_walls).T * ambient * sp.Matrix.hstack(
        *divisor_walls
    )
    trace_wall_gram = sp.Matrix.hstack(*trace_walls).T * trace_gram_matrix * sp.Matrix.hstack(
        *trace_walls
    )
    field_wall_gram = bugaenko_wall_gram()
    direct_trace_wall_gram = sp.Matrix(
        [
            [
                (K3_TRACE_SCALAR * field_wall_gram[row][column]).trace()
                for column in range(len(field_wall_gram))
            ]
            for row in range(len(field_wall_gram))
        ]
    )
    field_root_trace_pairing_verified = trace_wall_gram == direct_trace_wall_gram
    wall_pairing_transport = divisor_wall_gram == 4 * trace_wall_gram
    facet_divisor_squares = sorted(
        {int(divisor_wall_gram[index, index]) for index in range(len(divisor_walls))}
    )
    generated_trace_rank = sp.Matrix.hstack(
        *(standard_walls + [standard_action * wall for wall in standard_walls])
    ).rank()

    # The full trace pairing above is not itself the Coxeter-normal Gram
    # matrix.  Project each divisor normal to the identity eigenspace and
    # compare there with the original field Gram matrix.  This is the exact
    # non-tautological check that the eigenslice has Bugaenko's metric.
    sqrt_five = sp.sqrt(5)
    phi_real = (1 + sqrt_five) / 2
    phi_conjugate = (1 - sqrt_five) / 2
    identity_projector = (
        standard_action - phi_conjugate * sp.eye(standard_action.rows)
    ) / sqrt_five
    standard_wall_matrix = sp.Matrix.hstack(*standard_walls)
    projected_divisor_normals = embedding * identity_projector * standard_wall_matrix
    projected_wall_gram = sp.simplify(
        projected_divisor_normals.T * ambient * projected_divisor_normals
    )
    expected_projected_wall_gram = sp.Matrix(
        [
            [
                sp.simplify(
                    4
                    * _identity_value(K3_TRACE_SCALAR)
                    * _identity_value(field_wall_gram[row][column])
                )
                for column in range(len(field_wall_gram))
            ]
            for row in range(len(field_wall_gram))
        ]
    )
    identity_eigenspace_wall_gram_verified = all(
        sp.simplify(
            projected_wall_gram[row, column]
            - expected_projected_wall_gram[row, column]
        )
        == 0
        for row in range(projected_wall_gram.rows)
        for column in range(projected_wall_gram.cols)
    )

    # Twice the archived strict interior point is integral in the field
    # basis.  Its identity component selects the asserted positive component
    # and pairs strictly positively with all transported divisor walls.
    doubled_interior = [2 * entry for entry in bugaenko_interior_point()]
    standard_interior = inverse_change * _trace_coordinate_vector(doubled_interior)
    projected_divisor_interior = embedding * identity_projector * standard_interior
    interior_square = sp.simplify(
        (projected_divisor_interior.T * ambient * projected_divisor_interior)[0]
    )
    interior_wall_pairings = [
        sp.simplify(
            (projected_divisor_interior.T * ambient * divisor_wall)[0]
        )
        for divisor_wall in divisor_walls
    ]
    oriented_identity_component_verified = (
        interior_square.is_positive is True
        and all(value.is_positive is True for value in interior_wall_pairings)
        and all(
            sp.simplify(value) == 0
            for value in (
                ambient_operator * projected_divisor_interior
                - extension.operator_multiplier
                * phi_real
                * projected_divisor_interior
            )
        )
    )

    rational_points = [(1, index, index * index) for index in range(19)]
    rational_points_distinct = len(set(rational_points)) == 19

    required_checks = {
        "field_root_trace_pairing": field_root_trace_pairing_verified,
        "identity_eigenspace_wall_gram": identity_eigenspace_wall_gram_verified,
        "oriented_identity_component": oriented_identity_component_verified,
        "cycle_transpose_invariance": cycle_transpose_invariant,
        "operator_characteristic_polynomial": (
            operator_characteristic_polynomial_verified
        ),
        "chow_projector_idempotent": chow_projector_idempotent,
        "chow_field_action_relations": chow_field_action_relations,
        "chow_cycles_transpose_invariant": chow_cycles_transpose_invariant,
        "integral_generator_is_scaled_rational_action": (
            integral_generator_is_scaled_rational_action
        ),
        "rational_points_distinct": rational_points_distinct,
    }
    if not all(required_checks.values()):
        failed = [name for name, passed in required_checks.items() if not passed]
        raise AssertionError(f"failed rational-surface carrier checks: {failed}")

    completion = verify_bugaenko_vinberg_completion(
        build_bugaenko_vinberg_transcript(), bugaenko_h8_certificate()
    )
    if completion.diagram.decision != "compact":
        raise AssertionError("the inherited Bugaenko chamber is not compact")

    result = {
        "schema_version": "0.1.0",
        "certificate_kind": "rational_surface_eigencarrier_universality",
        "theorem": {
            "scope": (
                "every integral Lorentzian lattice, and hence every admissible "
                "totally-real arithmetic Coxeter datum after restriction of scalars"
            ),
            "generic_embedding": "D^2*G embeds integrally in I_(1,1+4m)",
            "rank_convention": "signature(G)=(1,m)",
            "surface": "Bl_(1+4m)(P^2)",
            "correspondence": (
                "a multiple of every rational self-adjoint source operator is "
                "induced on NS by a transpose-invariant integral sum of "
                "divisor products"
            ),
            "exact_hodge_index_gate": (
                "for an isometric full-trace-form carrier: one Lorentzian "
                "real place, negative-definite conjugate places, and hence a "
                "negative-definite NS orthogonal complement"
            ),
            "hodge_necessity_does_not_cover": (
                "an arbitrary unrelated surface realization that does not "
                "contain the full restriction-of-scalars trace form"
            ),
            "rational_chow_motive_action": (
                "the trace-space orthogonal projector is a divisor-product "
                "Chow idempotent, and F acts on its Tate-motive image; a "
                "chosen field generator becomes integral only after scaling"
            ),
            "does_not_assert": [
                "effective or irreducible facet divisors",
                "equality with the nef cone",
                "surface automorphisms",
                "integral ambient reflections",
            ],
        },
        "algorithm_calibration": {
            "gram": [
                [int(calibration_gram[row, column]) for column in range(3)]
                for row in range(3)
            ],
            "signature": [1, 2],
            "vector_scale": generic.vector_scale,
            "blowup_count": generic.blowup_count,
            "embedding_identity_verified": generic_identity,
        },
        "compact_calibration": {
            "polytope": "Bugaenko compact H8 chamber",
            "facet_count": 11,
            "finite_vertex_count": len(completion.chamber.finite_vertices),
            "cusp_count": len(completion.chamber.cusps),
            "hyperideal_ray_count": len(completion.chamber.hyperideal_rays),
            "vinberg_completion_verified": True,
            "surface_over_Q": "Bl_{p_0,...,p_18}(P^2_Q)",
            "points": "p_k=[1:k:k^2], 0<=k<=18",
            "points_distinct_verified": rational_points_distinct,
            "picard_lattice": "I_(1,19)",
            "trace_lattice": "<2>+<-2>+E8(-1)^2",
            "trace_rank": 18,
            "embedded_form_scale": 4,
            "embedding_shape": [embedding.rows, embedding.cols],
            "embedding_rank": int(embedding.rank()),
            "embedding_identity_verified": embedding_identity,
            "facet_divisor_count": len(divisor_walls),
            "facet_divisors_integral": wall_integrality,
            "facet_divisor_square_set": facet_divisor_squares,
            "facet_pairing_transport_verified": wall_pairing_transport,
            "field_root_trace_pairing_verified": field_root_trace_pairing_verified,
            "identity_eigenspace_projected_wall_gram_verified": (
                identity_eigenspace_wall_gram_verified
            ),
            "oriented_identity_component_verified": (
                oriented_identity_component_verified
            ),
            "walls_and_golden_translates_span_trace_rank": int(generated_trace_rank),
            "correspondence_operator_multiplier": extension.operator_multiplier,
            "correspondence_operator_rank": int(ambient_operator.rank()),
            "correspondence_operator_nonzero_entries": sum(
                entry != 0 for entry in ambient_operator
            ),
            "correspondence_max_abs_entry": max(
                abs(int(entry)) for entry in ambient_operator
            ),
            "divisor_product_term_count": len(cycle_terms),
            "divisor_product_cycle_action_verified": cycle_verified,
            "divisor_product_cycle_transpose_invariant": (
                cycle_transpose_invariant
            ),
            "operator_self_adjoint_verified": (
                ambient_operator.T * ambient == ambient * ambient_operator
            ),
            "operator_restriction_verified": (
                ambient_operator * embedding
                == embedding * (extension.operator_multiplier * standard_action)
            ),
            "operator_characteristic_polynomial_verified": (
                operator_characteristic_polynomial_verified
            ),
            "trace_chow_projector_rank": int(trace_projector.rank()),
            "trace_chow_projector_idempotent_under_cycle_composition": (
                chow_projector_idempotent
            ),
            "golden_chow_action_relations_under_cycle_composition": (
                chow_field_action_relations
            ),
            "rational_chow_action_cycles_transpose_invariant": (
                chow_cycles_transpose_invariant
            ),
            "integral_generator_is_scaled_rational_chow_action": (
                integral_generator_is_scaled_rational_action
            ),
            "chow_motive_summand": (
                "rank-18 direct summand of the rational divisor Tate motive; "
                "Q(phi) acts with the Chow projector as its unit"
            ),
            "real_eigenspace_signatures_verified": eigenspace_signatures_verified,
            "real_eigenspaces": {
                "4phi": {"dimension": 9, "signature": trace_eigenspace_signatures[0]},
                "4phi_conjugate": {
                    "dimension": 9,
                    "signature": trace_eigenspace_signatures[1],
                },
                "0": {
                    "dimension": 2,
                    "signature": [
                        complement_inertia.positive,
                        complement_inertia.negative,
                    ],
                },
            },
            "semialgebraic_section": (
                "the positive component of the 4phi eigenspace cut by the "
                "eleven transported divisor half-spaces"
            ),
            "section_is_exactly_bugaenko_h8": True,
        },
    }
    return result


if __name__ == "__main__":
    print(json.dumps(analyze(), indent=2, sort_keys=True))
