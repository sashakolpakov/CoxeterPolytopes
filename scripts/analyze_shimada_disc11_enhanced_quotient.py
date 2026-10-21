"""Exact enhanced-automorphism audit for the Shimada K3 traversal.

The global induced-chamber traversal repeatedly encounters physical copies of
the same Shimada representative.  A representative index alone does not say
whether the copies are related by an automorphism preserving the transported
golden action.  This module compares the cumulative row transports in the
marked Neron--Severi lattice.

It proves that the simultaneous returns at D68 and D27 give the *same* global
involution, determines its reflection action on the golden eigenspace, and
audits its quotient on a bounded exact traversal.  The bounded quotient is an
obstruction result, not a global-completeness assertion.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict, deque
import hashlib
import json
from pathlib import Path
import sys
from typing import Callable

from sympy import Matrix, Poly, cyclotomic_poly, symbols

from coxeter_polytopes.exact.coxeter import classify_root_pair
from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.exact.matrix import inertia, rank
from coxeter_polytopes.geometric.golden_halfspaces import (
    golden_dot,
    inverse_golden_matrix,
    positive_proportional_ratio,
)
from coxeter_polytopes.verify.chamber import diagnose_chamber_cone
from scripts.analyze_shimada_disc11_eigen_nef import (
    _induced_eigen_action,
    _restricted_coefficients,
)
from scripts.coupled_discriminant11_candidate import field_gram
from scripts.explore_shimada_disc11_global_eigen_nef import (
    _cross_internal_facet,
    _physical_state_key,
    _ray_key,
    _reduce_section,
    _root_and_internal_facets,
)
from scripts.fold_h11_seed_via_shimada import fold_h11_seed
from scripts.shimada_disc11_atlas import (
    ShimadaChamberState,
    ShimadaDisc11Atlas,
)
from scripts.shimada_disc11_conjugate import (
    AMPLE_FIELD_COEFFICIENTS,
    DISC11_SOURCE_SHA256,
    conjugated_phi_action_shimada,
    expected_shimada_gram,
    transported_identity_eigenbasis,
)


GoldenVector = tuple[Golden, ...]
PhysicalKey = tuple[tuple[int, ...], ...]
TraversalProgress = Callable[[str, int, ShimadaChamberState, int, int], None]

D68_RETURN_PATH = ((9, 0), (12, 0), (8, 1), (15, 0))
D27_RETURN_PATH = ((19, 0), (24, 0), (20, 1), (26, 0))
RETURN_TRANSPORT_SHA256 = (
    "1944304c8501f25118eaee6304aaf4038e690199a6e1a19847b2258645a65b4d"
)
RETURN_EIGEN_ACTION_SHA256 = (
    "5b6e9046cbbaf4cbffaf954c98a0332b169e3d246b851af9adeb43a8468cdefe"
)

# If an integral 20 by 20 matrix has finite order, every eigenvalue has order
# in this complete list (the orders n with Euler phi(n) <= 20).  Their lcm is
# a universal torsion exponent for GL(20,Z).  We first inspect the exact
# characteristic polynomial, so the large exponent is used only in the
# cyclotomic/unipotent edge case and cannot cause Salem-entry explosion.
_GL20_CYCLOTOMIC_ORDERS = (
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
    21, 22, 24, 25, 26, 27, 28, 30, 32, 33, 34, 36, 38, 40, 42, 44, 48, 50,
    54, 60, 66,
)
_GL20_TORSION_EXPONENT = 6_983_776_800
_GL20_TORSION_PRIMES = (2, 3, 5, 7, 11, 13, 17, 19)


def _canonical_sha(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def _integer_tuple(matrix: Matrix) -> tuple[tuple[int, ...], ...]:
    if any(not entry.is_Integer for entry in matrix):
        raise ValueError("an enhanced transport became nonintegral")
    return tuple(tuple(int(entry) for entry in row) for row in matrix.tolist())


def _decorated_state_record(
    key: PhysicalKey,
    state: ShimadaChamberState,
    interior: GoldenVector,
    depth: int,
    *,
    include_interior: bool = True,
) -> dict[str, object]:
    """Return a compact, canonical record for traversal-prefix hashing."""

    record: dict[str, object] = {
        "physical_state_sha256": _canonical_sha(key),
        "representative": state.index,
        "row_transport_sha256": _canonical_sha(state.row_isometry),
        "depth": depth,
    }
    if include_interior:
        record["interior"] = [entry.to_json() for entry in interior]
    return record


def _exact_matrix_order(matrix: Matrix) -> tuple[bool, int | None, str]:
    """Decide finite order for an integral matrix of rank at most 20."""

    if matrix.rows != matrix.cols or matrix.rows > 20:
        raise ValueError("the exact torsion test is scoped to square rank <= 20")
    variable = symbols("t")
    allowed = {
        Poly(cyclotomic_poly(order, variable), variable).as_expr()
        for order in _GL20_CYCLOTOMIC_ORDERS
    }
    factors = matrix.charpoly(variable).as_poly().factor_list()[1]
    noncyclotomic = [
        factor.monic().as_expr()
        for factor, _ in factors
        if factor.monic().as_expr() not in allowed
    ]
    if noncyclotomic:
        return False, None, "noncyclotomic_characteristic_factor"

    identity = Matrix.eye(matrix.rows)
    if matrix ** _GL20_TORSION_EXPONENT != identity:
        return False, None, "nonsemisimple_cyclotomic_part"

    order = _GL20_TORSION_EXPONENT
    for prime in _GL20_TORSION_PRIMES:
        while order % prime == 0 and matrix ** (order // prime) == identity:
            order //= prime
    return True, order, "universal_GL20_torsion_exponent"


def _walk(atlas, state, path):
    states = [state]
    for orbit, member in path:
        neighbor = atlas.neighboring_state(state, orbit, member)
        if neighbor is None:
            raise ValueError("an asserted enhanced return hit a K3 root wall")
        state = neighbor
        states.append(state)
    return state, tuple(states)


def _right_action(
    state: ShimadaChamberState,
    matrix: Matrix,
) -> ShimadaChamberState:
    product = Matrix(state.row_isometry) * matrix
    return ShimadaChamberState(state.index, _integer_tuple(product))


def _retained_traversal(
    atlas: ShimadaDisc11Atlas,
    *,
    max_states: int | None = None,
    max_depth: int | None = None,
    progress: TraversalProgress | None = None,
):
    """Run the exact global kernel while retaining marked state transports."""

    if max_states is None and max_depth is None:
        raise ValueError("a state or depth bound is required")
    if max_states is not None and max_states < 1:
        raise ValueError("max_states must be positive")
    if max_depth is not None and max_depth < 0:
        raise ValueError("max_depth must be nonnegative")

    start = atlas.representative_state(68)
    start_key = _physical_state_key(atlas, start)
    queue = deque([(start, tuple(AMPLE_FIELD_COEFFICIENTS), 0)])
    queued = {start_key}
    processed: dict[PhysicalKey, tuple[ShimadaChamberState, GoldenVector]] = {}
    depths: dict[PhysicalKey, int] = {start_key: 0}
    visible_root_rays: dict[GoldenVector, int] = {}
    edges: list[tuple[PhysicalKey, PhysicalKey]] = []

    while queue:
        if max_states is not None and len(processed) >= max_states:
            break
        if max_depth is not None and queue[0][2] > max_depth:
            break
        state, interior, depth = queue.popleft()
        key = _physical_state_key(atlas, state)
        if key in processed:
            continue
        if progress is not None:
            progress("start", len(processed) + 1, state, depth, len(queue))
        reduction = _reduce_section(atlas, state, interior)
        root_facets, internal_facets = _root_and_internal_facets(
            atlas, state, reduction
        )
        processed[key] = (state, interior)
        for wall_class in root_facets:
            ray = _ray_key(reduction.rows[wall_class[0]])
            visible_root_rays[ray] = visible_root_rays.get(ray, 0) + 1
        for wall_class in internal_facets:
            neighbor, opposite, _, _ = _cross_internal_facet(
                atlas,
                state,
                interior,
                reduction.witnesses[wall_class[0]],
                wall_class,
            )
            neighbor_key = _physical_state_key(atlas, neighbor)
            edges.append((key, neighbor_key))
            if neighbor_key not in queued and neighbor_key not in processed:
                queued.add(neighbor_key)
                depths[neighbor_key] = depth + 1
                queue.append((neighbor, opposite, depth + 1))
        if progress is not None:
            progress("finish", len(processed), state, depth, len(queue))
    return processed, tuple(queue), tuple(visible_root_rays), tuple(edges), depths


def _centralizer_ball_certificate(
    atlas: ShimadaDisc11Atlas,
    processed: dict[PhysicalKey, tuple[ShimadaChamberState, GoldenVector]],
    queue,
    depths: dict[PhysicalKey, int],
    tau: Matrix,
) -> dict[str, object]:
    """Recognize enhanced automorphisms in a one-witness traversal sample.

    A physical state ``(68,g)`` represents the chamber ``D_68 g``.  If an
    automorphism ``gamma`` maps the base chamber to that state, then
    ``g gamma^-1`` lies in the finite representative stabilizer ``H_68``;
    hence ``gamma=a^-1 g`` for some ``a in H_68``.  Enumerating that finite
    stabilizer is essential: the stored atlas transport alone is not a
    canonical automorphism label.

    The traversal supplying ``processed`` chooses one point on each restricted
    facet.  The D35 regression proves that this can omit ambient subdivision
    cells, so its breadth-first layers are *not* metric balls in the complete
    carrier.  The return labels found in the sample remain exact.
    """

    if queue:
        first_unfinished_depth = min(depth for _, _, depth in queue)
        complete_radius = first_unfinished_depth - 1
        traversal_closed = False
    else:
        complete_radius = max(depths.values(), default=0)
        traversal_closed = True

    layer_counts = Counter(
        depths[key] for key in processed if depths[key] <= complete_radius
    )
    base_stabilizer = tuple(Matrix(item) for item in atlas.stabilizer_group(68))
    centralizing: dict[tuple[tuple[int, ...], ...], dict[str, object]] = {}
    tested_returns: list[dict[str, object]] = []
    complete_state_records = []
    for key, (state, interior) in processed.items():
        depth = depths[key]
        if depth <= complete_radius:
            complete_state_records.append(
                _decorated_state_record(
                    key, state, interior, depth, include_interior=False
                )
            )
    complete_state_records.sort(
        key=lambda item: (item["depth"], item["physical_state_sha256"])
    )
    base_states = 0
    for key, (state, _) in processed.items():
        depth = depths[key]
        if depth > complete_radius or state.index != 68:
            continue
        base_states += 1
        transport = Matrix(state.row_isometry)
        for stabilizer_index, stabilizer in enumerate(base_stabilizer):
            candidate = stabilizer.inv() * transport
            candidate_tuple = _integer_tuple(candidate)
            finite_order, exact_order, order_reason = _exact_matrix_order(candidate)
            commutes = candidate.T * tau == tau * candidate.T
            return_record = {
                "one_witness_graph_depth": depth,
                "physical_state_sha256": _canonical_sha(key),
                "target_stabilizer_index": stabilizer_index,
                "commutes_with_tau": commutes,
                "determinant": str(candidate.det()),
                "characteristic_polynomial": str(
                    candidate.charpoly().as_expr().factor()
                ),
                "finite_order": finite_order,
                "exact_order": exact_order,
                "order_certificate": order_reason,
                "matrix_sha256": _canonical_sha(candidate_tuple),
            }
            tested_returns.append(return_record)
            if not commutes:
                continue
            centralizing.setdefault(
                candidate_tuple,
                return_record,
            )

    infinite_order_found = any(
        not item["finite_order"] for item in centralizing.values()
    )
    return {
        "facet_star_complete": False,
        "complete_restricted_facet_radius": None,
        "one_witness_prefix_radius": complete_radius,
        "one_witness_prefix_layer_counts": {
            str(depth): layer_counts[depth]
            for depth in sorted(layer_counts)
        },
        "states_in_one_witness_prefix_layers": sum(layer_counts.values()),
        "one_witness_prefix_state_records_sha256": _canonical_sha(
            complete_state_records
        ),
        "D68_type_states_in_one_witness_prefix": base_states,
        "D68_representative_stabilizer_order": len(base_stabilizer),
        "D68_return_candidates_tested": len(tested_returns),
        "D68_return_records": sorted(
            tested_returns,
            key=lambda item: (
                item["one_witness_graph_depth"],
                item["physical_state_sha256"],
                item["target_stabilizer_index"],
            ),
        ),
        "D68_return_records_sha256": _canonical_sha(
            sorted(
                tested_returns,
                key=lambda item: (
                    item["one_witness_graph_depth"],
                    item["physical_state_sha256"],
                    item["target_stabilizer_index"],
                ),
            )
        ),
        "centralizing_automorphisms_in_sample": len(centralizing),
        "centralizing_records": sorted(
            centralizing.values(),
            key=lambda item: (
                item["one_witness_graph_depth"],
                item["matrix_sha256"],
            ),
        ),
        "only_identity_and_A_in_one_witness_sample": (
            len(centralizing) == 2
            and sorted(item["exact_order"] for item in centralizing.values())
            == [1, 2]
            and sorted(
                item["characteristic_polynomial"]
                for item in centralizing.values()
            )
            == [
                "(lambda - 1)**18*(lambda + 1)**2",
                "(lambda - 1)**20",
            ]
        ),
        "centralizer_in_complete_ball_is_exactly_C2": False,
        "new_centralizer_requires_distance_at_least": None,
        "one_witness_queue_closed": traversal_closed,
        "traversal_closed": False,
        "infinite_order_return_found": infinite_order_found,
        "decision_status": (
            "infinite" if infinite_order_found
            else "pending_exact_ambient_carrier_enumeration"
        ),
        "withdrawal_reason": (
            "D35_has_two_distinct_opposite_induced_cells_across_one_"
            "restricted_facet"
        ),
        "decision_procedure": {
            "finite_stop": (
                "the_exact_ambient_carrier_queue_closes; then only_finitely_"
                "many_D68_images_and_finite_local_stabilizers_remain"
            ),
            "infinite_stop": (
                "a_D68_return_label_centralizes_tau_and_has_exact_infinite_order"
            ),
            "why_one_stop_occurs": (
                "finite_enhanced_centralizer_implies_a_compact_eigen_nef_"
                "chamber_and_hence_finitely_many_locally_finite_induced_cells;_"
                "an_infinite_finitely_generated_linear_centralizer_contains_an_"
                "infinite_order_element_whose_D68_image_occurs_at_finite_graph_"
                "distance"
            ),
            "implementation_totality": (
                "the_replacement_ambient_carrier_crosses_every_auxiliary_"
                "wall_member_and_uses_exact_closed_Q(phi)_phase_I_with_a_"
                "no_Slater_Lorentzian_KKT_fallback"
            ),
        },
        "scope": (
            "exact_return_labels_in_a_deterministic_one_witness_sample;_not_"
            "a_complete_metric_ball_or_displacement_bound"
        ),
    }


def _orbits(values: set[object], permutation: dict[object, object]):
    seen = set()
    result = []
    for value in values:
        if value in seen:
            continue
        image = permutation[value]
        orbit = frozenset((value, image))
        seen.update(orbit)
        result.append(orbit)
    return tuple(result)


def _active_set_orbits(values, facet_permutation):
    sets = {tuple(value) for value in values}
    permutation = {
        value: tuple(sorted(facet_permutation[index] for index in value))
        for value in sets
    }
    closed = set(permutation.values()) == sets
    orbits = _orbits(sets, permutation) if closed else ()
    return {
        "count": len(sets),
        "invariant": closed,
        "fixed": sum(value == permutation[value] for value in sets) if closed else 0,
        "quotient_orbits": len(orbits),
    }


def _root_quotient(rays, eigen_action):
    ray_index = {ray: index for index, ray in enumerate(rays)}
    permutation = []
    scales = []
    for ray in rays:
        image = tuple(
            sum(
                (eigen_action[i][j] * ray[j] for j in range(10)),
                Golden(0),
            )
            for i in range(10)
        )
        scale = golden_dot(image, AMPLE_FIELD_COEFFICIENTS)
        normalized = _ray_key(image)
        permutation.append(ray_index.get(normalized, -1))
        scales.append(scale)

    gram = tuple(tuple(entry for entry in row) for row in field_gram())
    inverse = inverse_golden_matrix(gram)
    roots = tuple(
        tuple(
            -sum((inverse[i][j] * functional[j] for j in range(10)), Golden(0))
            for i in range(10)
        )
        for functional in rays
    )

    def pair(left, right):
        return sum(
            (
                left[i] * gram[i][j] * right[j]
                for i in range(10)
                for j in range(10)
            ),
            Golden(0),
        )

    wall_gram = tuple(
        tuple(pair(left, right) for right in roots) for left in roots
    )
    relation_counts = Counter()
    all_coxeter = True
    for i in range(len(roots)):
        for j in range(i + 1, len(roots)):
            product = wall_gram[i][j]
            if product.sign() > 0:
                all_coxeter = False
                relation_counts[("nonacute", "none", str(product))] += 1
                continue
            relation = classify_root_pair(
                wall_gram[i][i], product, wall_gram[j][j]
            )
            all_coxeter &= relation.kind != "finite_unrecognized"
            relation_counts[
                (
                    relation.kind,
                    str(relation.label) if relation.label is not None else "none",
                    str(relation.q),
                )
            ] += 1

    valid_permutation = (
        -1 not in permutation
        and all(permutation[permutation[index]] == index for index in range(len(rays)))
    )
    scaled_gram_invariant = valid_permutation and all(
        wall_gram[i][j]
        == scales[i] * scales[j] * wall_gram[permutation[i]][permutation[j]]
        for i in range(len(rays))
        for j in range(len(rays))
    )
    geometry = diagnose_chamber_cone(
        roots, gram, tuple(AMPLE_FIELD_COEFFICIENTS), 9
    )
    finite = _active_set_orbits(geometry.finite_vertices, permutation)
    cusps = _active_set_orbits(geometry.cusps, permutation)
    hyperideal = _active_set_orbits(geometry.hyperideal_rays, permutation)
    opposite = _active_set_orbits(geometry.opposite_sheet_rays, permutation)
    return {
        "visible_root_rays": len(rays),
        "permutation": permutation,
        "positive_rescaling_factors": [value.to_json() for value in scales],
        "permutation_is_involution": valid_permutation,
        "fixed_root_rays": sum(index == image for index, image in enumerate(permutation)),
        "transposed_root_ray_pairs": [
            [index, image]
            for index, image in enumerate(permutation)
            if index < image
        ],
        "scaled_wall_gram_invariant": scaled_gram_invariant,
        "wall_gram_rank": rank(wall_gram),
        "wall_gram_inertia": {
            "positive": inertia(wall_gram).positive,
            "negative": inertia(wall_gram).negative,
            "zero": inertia(wall_gram).zero,
        },
        "all_pairs_coxeter": all_coxeter,
        "coxeter_relation_counts": [
            {
                "kind": key[0],
                "label": key[1],
                "q": key[2],
                "count": count,
            }
            for key, count in sorted(relation_counts.items())
        ],
        "outer_cone": {
            "finite_vertices": finite,
            "cusps": cusps,
            "hyperideal_rays": hyperideal,
            "opposite_sheet_rays": opposite,
            "quotient_still_has_hyperideal_rays": hyperideal["quotient_orbits"] > 0,
        },
    }


def analyze_centralizer_ball(
    path: Path,
    *,
    traversal_states: int = 50,
    complete_radius: int | None = None,
    progress: TraversalProgress | None = None,
) -> dict[str, object]:
    """Run the historical one-witness decorated traversal diagnostic.

    The local reductions and every reported return are exact, but D35 proves
    that the sampled restricted-facet graph is not the complete ambient
    carrier.  ``complete_radius`` is retained as a compatibility spelling for
    the requested one-witness depth and never produces a completeness claim.
    """

    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != DISC11_SOURCE_SHA256:
        raise ValueError("unexpected disc11.txt checksum")
    if traversal_states < 1:
        raise ValueError("traversal_states must be positive")
    if complete_radius is not None and complete_radius < 0:
        raise ValueError("complete_radius must be nonnegative")

    atlas = ShimadaDisc11Atlas(path)
    tau = Matrix(conjugated_phi_action_shimada())
    processed, queue, root_rays, edges, depths = _retained_traversal(
        atlas,
        max_states=None if complete_radius is not None else traversal_states,
        max_depth=complete_radius,
        progress=progress,
    )
    centralizer_ball = _centralizer_ball_certificate(
        atlas, processed, queue, depths, tau
    )

    processed_checkpoint_records = sorted(
        (
            _decorated_state_record(key, state, interior, depths[key])
            for key, (state, interior) in processed.items()
        ),
        key=lambda item: (item["depth"], item["physical_state_sha256"]),
    )
    processed_records = [
        {key: value for key, value in item.items() if key != "interior"}
        for item in processed_checkpoint_records
    ]
    queue_checkpoint_records = [
        _decorated_state_record(
            _physical_state_key(atlas, state), state, interior, depth
        )
        for state, interior, depth in queue
    ]
    queue_records = [
        {key: value for key, value in item.items() if key != "interior"}
        for item in queue_checkpoint_records
    ]
    edge_records = sorted(
        (
            {
                "source_physical_state_sha256": _canonical_sha(source),
                "target_physical_state_sha256": _canonical_sha(target),
            }
            for source, target in edges
        ),
        key=lambda item: (
            item["source_physical_state_sha256"],
            item["target_physical_state_sha256"],
        ),
    )
    root_ray_records = sorted(
        ([entry.to_json() for entry in ray] for ray in root_rays),
        key=_canonical_sha,
    )
    geometric_prefix_payload = {
        "processed": processed_records,
        "queue": queue_records,
        "edges": edge_records,
        "root_rays": root_ray_records,
    }
    checkpoint_prefix_payload = {
        "processed": processed_checkpoint_records,
        "queue": queue_checkpoint_records,
        "edges": edge_records,
        "root_rays": root_ray_records,
    }
    return {
        "analysis_id": "shimada_disc11_enhanced_one_witness_sample",
        "verified": False,
        "verified_exact_local_certificates": True,
        "source_sha256": digest,
        "requested_complete_radius": None,
        "requested_one_witness_radius": complete_radius,
        "traversal_state_cutoff": (
            None if complete_radius is not None else traversal_states
        ),
        "processed_physical_states": len(processed),
        "queued_unprocessed_states": len(queue),
        "discovered_physical_states": len(depths),
        "internal_crossings_from_processed_states": len(edges),
        "visible_root_rays": len(root_rays),
        "processed_records_sha256": _canonical_sha(processed_records),
        "queue_records_sha256": _canonical_sha(queue_records),
        "processed_checkpoint_records_sha256": _canonical_sha(
            processed_checkpoint_records
        ),
        "queue_checkpoint_records_sha256": _canonical_sha(
            queue_checkpoint_records
        ),
        "edge_records_sha256": _canonical_sha(edge_records),
        "root_ray_records_sha256": _canonical_sha(root_ray_records),
        "traversal_prefix_sha256": _canonical_sha(geometric_prefix_payload),
        "checkpoint_prefix_sha256": _canonical_sha(checkpoint_prefix_payload),
        "centralizer_decision_ball": centralizer_ball,
        "omitted_independent_diagnostics": [
            "Pell_seed_factorization",
            "root_wall_outer_cone",
            "C2_closed_bounded_state_quotient",
        ],
        "scope": (
            "exact_local_certificates_and_return_labels_in_a_one_witness_"
            "sample;_D35_proves_this_is_not_a_complete_metric_ball"
        ),
    }


def analyze_enhanced_quotient(
    path: Path,
    *,
    traversal_states: int = 50,
) -> dict[str, object]:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != DISC11_SOURCE_SHA256:
        raise ValueError("unexpected disc11.txt checksum")
    if traversal_states < 1:
        raise ValueError("traversal_states must be positive")

    atlas = ShimadaDisc11Atlas(path)
    d68 = atlas.representative_state(68)
    d68_return, d68_gallery = _walk(atlas, d68, D68_RETURN_PATH)
    if d68_return.index != 68:
        raise ValueError("the D68 simultaneous gallery did not return to D68")
    return_matrix = Matrix(d68_return.row_isometry)
    return_tuple = _integer_tuple(return_matrix)

    d86 = atlas.neighboring_state(d68, 11, 1)
    if d86 is None or d86.index != 86:
        raise ValueError("D68 wall (11,1) did not enter D86")
    d27 = atlas.neighboring_state(d86, 24, 0)
    if d27 is None or d27.index != 27:
        raise ValueError("D86 wall (24,0) did not enter D27")
    d27_return, d27_gallery = _walk(atlas, d27, D27_RETURN_PATH)
    if d27_return.index != 27:
        raise ValueError("the D27 simultaneous gallery did not return to D27")
    d27_base_matrix = Matrix(d27.row_isometry)
    d27_return_matrix = Matrix(d27_return.row_isometry)
    d27_global_relative = d27_base_matrix.inv() * d27_return_matrix
    d27_local_left = d27_return_matrix * d27_base_matrix.inv()

    gram = Matrix(expected_shimada_gram())
    tau = Matrix(conjugated_phi_action_shimada())
    discriminant_generator = gram.inv()[:, 2]
    basis = transported_identity_eigenbasis()
    eigen_action = _induced_eigen_action(return_tuple, basis)
    identity = [[Golden(int(i == j)) for j in range(10)] for i in range(10)]
    expected_eigen_action = [row[:] for row in identity]
    expected_eigen_action[1][0] = Golden(1, -1)
    expected_eigen_action[1][1] = Golden(-1)
    eigen_action_sha = _canonical_sha(
        [[entry.to_json() for entry in row] for row in eigen_action]
    )
    eigen_minus = [
        [eigen_action[i][j] - identity[i][j] for j in range(10)]
        for i in range(10)
    ]
    eigen_plus = [
        [eigen_action[i][j] + identity[i][j] for j in range(10)]
        for i in range(10)
    ]

    d68_rows = atlas.state_wall_orbits(d68)
    d27_rows = atlas.state_wall_orbits(d27)
    mirror_d68 = _restricted_coefficients(d68_rows[9][0], basis)
    mirror_d68_copy = _restricted_coefficients(d68_rows[16][0], basis)
    mirror_d27 = _restricted_coefficients(d27_rows[19][0], basis)
    mirror_d27_copy = _restricted_coefficients(d27_rows[31][0], basis)
    mirror_image = tuple(
        sum(
            (eigen_action[i][j] * mirror_d68[j] for j in range(10)),
            Golden(0),
        )
        for i in range(10)
    )
    mirrors_equal = (
        positive_proportional_ratio(mirror_d68, mirror_d68_copy) == Golden(0, 1)
        and positive_proportional_ratio(mirror_d68, mirror_d27) == Golden(1)
        and positive_proportional_ratio(mirror_d68, mirror_d27_copy)
        == Golden(0, 1)
    )

    processed, queue, root_rays, edges, depths = _retained_traversal(
        atlas, max_states=traversal_states
    )
    known_states = {key: state for key, (state, _) in processed.items()}
    for state, _, _ in queue:
        known_states[_physical_state_key(atlas, state)] = state
    action_on_known: dict[PhysicalKey, PhysicalKey] = {}
    image_states: dict[PhysicalKey, ShimadaChamberState] = {}
    for key, state in known_states.items():
        image = _right_action(state, return_matrix)
        image_key = _physical_state_key(atlas, image)
        action_on_known[key] = image_key
        image_states[image_key] = image
    known_closed = set(action_on_known.values()) == set(known_states)
    # A hard BFS cutoff need not contain both members of each symmetry orbit.
    # Close the retained set under A before counting quotient states.
    closed_states = dict(known_states)
    closed_states.update(image_states)
    action_on_closed: dict[PhysicalKey, PhysicalKey] = {}
    for key, state in closed_states.items():
        image = _right_action(state, return_matrix)
        action_on_closed[key] = _physical_state_key(atlas, image)
    closure_invariant = set(action_on_closed.values()) == set(closed_states)
    closed_orbits = (
        _orbits(set(closed_states), action_on_closed)
        if closure_invariant
        else ()
    )
    processed_keys = set(processed)
    queue_keys = {
        _physical_state_key(atlas, state) for state, _, _ in queue
    }
    processed_orbits = sum(bool(orbit & processed_keys) for orbit in closed_orbits)
    frontier_only_orbits = sum(
        not bool(orbit & processed_keys) for orbit in closed_orbits
    )

    edge_pairs = set(edges)
    equivariant_edges_tested = 0
    equivariant_edges_verified = True
    for source, target in edges:
        image_source = action_on_closed[source]
        image_target = action_on_closed[target]
        if image_source in processed_keys:
            equivariant_edges_tested += 1
            equivariant_edges_verified &= (image_source, image_target) in edge_pairs

    by_representative = defaultdict(list)
    for key, (state, _) in processed.items():
        by_representative[state.index].append((key, state))
    repeated_pair_comparisons = 0
    centralizing_pair_count = 0
    centralizing_relatives: list[Matrix] = []
    for records in by_representative.values():
        for i in range(len(records)):
            for j in range(i + 1, len(records)):
                repeated_pair_comparisons += 1
                left = Matrix(records[i][1].row_isometry)
                right = Matrix(records[j][1].row_isometry)
                relative = left.inv() * right
                if relative.T * tau != tau * relative.T:
                    continue
                centralizing_pair_count += 1
                if not any(relative == previous for previous in centralizing_relatives):
                    centralizing_relatives.append(relative)

    root_quotient = _root_quotient(root_rays, eigen_action)
    centralizer_ball = _centralizer_ball_certificate(
        atlas, processed, queue, depths, tau
    )
    pell_fold_records = []
    pell_fold_verified = True
    for exponent, expected in ((1, eigen_action), (2, identity)):
        fold = fold_h11_seed(
            path,
            seed_power=exponent,
            max_root_reflections=4,
            max_gallery=2000,
        )
        folded = tuple(
            tuple(Golden.from_json(entry) for entry in row)
            for row in fold["folded_row_action"]
        )
        expected_tuple = tuple(tuple(row) for row in expected)
        item_verified = (
            fold["terminal_same_chamber"]
            and fold["terminal_representative"] == 68
            and fold["root_reflections_used"] == exponent
            and folded == expected_tuple
            and not fold["verified_infinite_order_same_chamber_obstruction"]
        )
        pell_fold_verified &= item_verified
        pell_fold_records.append(
            {
                "seed_power": exponent,
                "root_reflections_used": fold["root_reflections_used"],
                "terminal_representative": fold["terminal_representative"],
                "residual": "E" if exponent % 2 else "identity",
                "verified": item_verified,
            }
        )
    discriminant_minus_identity = all(
        entry.is_Integer
        for entry in return_matrix.T * discriminant_generator
        + discriminant_generator
    )
    discriminant_identity = all(
        entry.is_Integer
        for entry in return_matrix.T * discriminant_generator
        - discriminant_generator
    )
    direct_verified = (
        _canonical_sha(return_tuple) == RETURN_TRANSPORT_SHA256
        and return_matrix.det() == 1
        and return_matrix * gram * return_matrix.T == gram
        and return_matrix * return_matrix == Matrix.eye(20)
        and return_matrix.T * tau == tau * return_matrix.T
        and not discriminant_identity
        and discriminant_minus_identity
        and d27_global_relative == return_matrix
        and d27_return_matrix == d27_base_matrix * return_matrix
        and d27_local_left.T * tau != tau * d27_local_left.T
        and eigen_action == expected_eigen_action
        and eigen_action_sha == RETURN_EIGEN_ACTION_SHA256
        and rank(eigen_minus) == 1
        and rank(eigen_plus) == 9
        and mirrors_equal
        and mirror_image == tuple(-entry for entry in mirror_d68)
    )
    quotient_verified = (
        closure_invariant
        and all(
            action_on_closed[action_on_closed[key]] == key
            for key in closed_states
        )
        and not any(action_on_closed[key] == key for key in closed_states)
        and len(centralizing_relatives) == 1
        and centralizing_relatives[0] == return_matrix
        and root_quotient["permutation_is_involution"]
        and root_quotient["scaled_wall_gram_invariant"]
        and root_quotient["outer_cone"]["quotient_still_has_hyperideal_rays"]
        and centralizer_ball["only_identity_and_A_in_one_witness_sample"]
    )
    verified = direct_verified and quotient_verified and pell_fold_verified
    return {
        "verified": verified,
        "source_sha256": digest,
        "enhanced_return": {
            "D68_gallery_representatives": [state.index for state in d68_gallery],
            "D27_gallery_representatives": [state.index for state in d27_gallery],
            "transport_sha256": RETURN_TRANSPORT_SHA256,
            "determinant": str(return_matrix.det()),
            "NS_characteristic_polynomial": str(
                return_matrix.charpoly().as_expr().factor()
            ),
            "NS_isometry": return_matrix * gram * return_matrix.T == gram,
            "involution": return_matrix * return_matrix == Matrix.eye(20),
            "centralizes_tau": return_matrix.T * tau == tau * return_matrix.T,
            "discriminant_action": "minus_identity",
            "discriminant_action_verified": (
                discriminant_minus_identity and not discriminant_identity
            ),
            "K3_extension": (
                "extend_by_minus_identity_on_the_oriented_rank2_"
                "transcendental_lattice_T11"
            ),
            "D27_global_relative_equals_D68_return": (
                d27_global_relative == return_matrix
            ),
            "D27_return_state_equals_base_times_D68_return": (
                d27_return_matrix == d27_base_matrix * return_matrix
            ),
            "D27_local_left_matrix_centralizes_tau": (
                d27_local_left.T * tau == tau * d27_local_left.T
            ),
            "golden_eigen_action_sha256": RETURN_EIGEN_ACTION_SHA256,
            "golden_eigen_action": [
                [entry.to_json() for entry in row] for row in eigen_action
            ],
            "golden_fixed_rank": 9,
            "golden_negative_rank": 1,
            "common_mirror_functional": [
                entry.to_json() for entry in mirror_d68
            ],
            "D68_and_D27_return_wall_classes_are_same_mirror": mirrors_equal,
            "mirror_is_negated": mirror_image
            == tuple(-entry for entry in mirror_d68),
            "generated_Coxeter_subgroup": "W(A1)=C2",
            "generated_subgroup_order": 2,
            "direct_certificate_verified": direct_verified,
        },
        "bounded_state_quotient": {
            "processed_states": len(processed),
            "queued_states": len(queue),
            "known_states": len(known_states),
            "known_state_set_A_invariant": known_closed,
            "A_images_missing_at_hard_cutoff": len(set(image_states) - set(known_states)),
            "A_closed_state_set": len(closed_states),
            "A_closure_invariant": closure_invariant,
            "A_closed_orbits": len(closed_orbits),
            "processed_orbits": processed_orbits,
            "frontier_only_orbits": frontier_only_orbits,
            "fixed_physical_states": sum(
                action_on_closed[key] == key for key in closed_states
            ),
            "internal_crossings_from_processed_states": len(edges),
            "equivariant_edges_tested": equivariant_edges_tested,
            "equivariant_edges_verified": equivariant_edges_verified,
            "representative_counts": dict(
                sorted(
                    Counter(
                        state.index for state, _ in processed.values()
                    ).items()
                )
            ),
            "same_representative_pair_comparisons": repeated_pair_comparisons,
            "centralizing_same_representative_pairs": centralizing_pair_count,
            "distinct_centralizing_relative_transports": len(
                centralizing_relatives
            ),
            "only_centralizing_relative_is_D68_involution": (
                len(centralizing_relatives) == 1
                and centralizing_relatives[0] == return_matrix
            ),
            "quotient_closed": not queue,
            "obstruction": (
                "the_return_generated_C2_halves_the_corridor_but_leaves_"
                "frontier_orbits_and_supplies_no_longitudinal_generator"
            ),
        },
        "centralizer_decision_ball": centralizer_ball,
        "Pell_seed_factorization": {
            "records": pell_fold_records,
            "verified": pell_fold_verified,
            "conclusion": (
                "the_infinite_Pell_stream_folds_modulo_root_reflections_"
                "to_E_power_k_and_supplies_no_second_stabilizer_generator"
            ),
        },
        "root_wall_quotient": root_quotient,
        "scope": (
            "exact_enhanced_return_and_bounded_quotient_obstruction_"
            "not_global_centralizer_or_wall_completeness"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("disc11", type=Path)
    parser.add_argument("--traversal-states", type=int, default=50)
    parser.add_argument(
        "--complete-radius",
        type=int,
        help=(
            "legacy spelling: process the one-witness sample through this "
            "depth; this is not a complete carrier radius"
        ),
    )
    parser.add_argument(
        "--centralizer-ball-only",
        action="store_true",
        help="omit the independent Pell, root-cone, and C2-quotient audits",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        help="write every Nth state start/finish record to stderr",
    )
    args = parser.parse_args()
    if args.complete_radius is not None and not args.centralizer_ball_only:
        parser.error("--complete-radius requires --centralizer-ball-only")
    if args.progress_every is not None and args.progress_every < 1:
        parser.error("--progress-every must be positive")

    def progress(event, count, state, depth, queued):
        if args.progress_every is None or count % args.progress_every:
            return
        print(
            f"traversal {event} state={count} depth={depth} "
            f"representative=D{state.index} queued={queued}",
            file=sys.stderr,
            flush=True,
        )

    analysis = (
        analyze_centralizer_ball(
            args.disc11,
            traversal_states=args.traversal_states,
            complete_radius=args.complete_radius,
            progress=progress if args.progress_every is not None else None,
        )
        if args.centralizer_ball_only
        else analyze_enhanced_quotient(
            args.disc11,
            traversal_states=args.traversal_states,
        )
    )
    print(
        json.dumps(
            analysis,
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
