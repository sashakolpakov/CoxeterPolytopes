"""Exact ambient-chamber carrier for the discriminant-11 eigensection.

The old traversal crossed one sample point on each *restricted* facet.  A
restricted facet can itself be subdivided by ambient Shimada walls, so that
graph is not complete.  The carrier used here keeps the ambient arrangement:

* a vertex is a decorated physical Shimada chamber whose closed restriction
  to the golden eigenspace contains a future-timelike point;
* every member of every auxiliary ambient wall orbit is crossed; and
* ``_minustwo`` walls are nef boundaries and are never crossed.

Closed feasibility is essential.  A chamber can belong to the carrier even
when its eigensection is lower-dimensional.  The generic exact gate uses
phase I and a no-Slater Lorentzian KKT certificate.
"""

from __future__ import annotations

import argparse
from collections import Counter, deque
from dataclasses import dataclass
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
from typing import Callable, Iterator

from sympy import Matrix

from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.geometric.golden_halfspaces import (
    ClosedLorentzianConeDecision,
    decide_lorentzian_closed_cone,
)
from scripts.certify_shimada_disc11_d35_facet_subdivision import (
    DEFAULT_FIXTURE,
)
from scripts.certify_shimada_disc11_bugaenko_v4 import (
    EXPECTED_FIELD_MATRIX_SHA256,
)
from scripts.analyze_shimada_disc11_eigen_nef import _induced_eigen_action
from scripts.analyze_shimada_disc11_enhanced_quotient import (
    _exact_matrix_order,
    _integer_tuple,
)
from scripts.explore_shimada_disc11_global_eigen_nef import (
    _Q_GRAM,
    _cross_internal_facet,
    _physical_state_key,
    _restricted_rows,
)
from scripts.shimada_disc11_atlas import ShimadaChamberState, ShimadaDisc11Atlas
from scripts.shimada_disc11_conjugate import (
    AMPLE_FIELD_COEFFICIENTS,
    DISC11_SOURCE_SHA256,
    conjugated_phi_action_shimada,
    transported_identity_eigenbasis,
)


PhysicalKey = tuple[tuple[int, ...], ...]
Progress = Callable[[str, int, ShimadaChamberState, int, int], None]
_CHECKPOINT_ALGORITHM = "ambient_carrier_closed_phaseI_KKT_v1"


def _canonical_sha(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


@dataclass(frozen=True, slots=True)
class AmbientCarrierEdge:
    source: PhysicalKey
    target: PhysicalKey
    source_representative: int
    target_representative: int
    orbit: int
    member: int
    target_retained: bool


def _state_to_json(state: ShimadaChamberState, depth: int) -> dict[str, object]:
    return {
        "representative": state.index,
        "row_isometry": [list(row) for row in state.row_isometry],
        "depth": depth,
    }


def _state_from_json(record: dict[str, object]) -> tuple[ShimadaChamberState, int]:
    raw_rows = record["row_isometry"]
    if not isinstance(raw_rows, list):
        raise ValueError("invalid carrier checkpoint row isometry")
    state = ShimadaChamberState(
        int(record["representative"]),
        tuple(tuple(int(entry) for entry in row) for row in raw_rows),
    )
    return state, int(record["depth"])


def _key_to_json(key: PhysicalKey) -> list[list[int]]:
    return [list(row) for row in key]


def _key_from_json(rows: object) -> PhysicalKey:
    if not isinstance(rows, list):
        raise ValueError("invalid carrier checkpoint physical key")
    return tuple(tuple(int(entry) for entry in row) for row in rows)


def _edge_to_json(edge: AmbientCarrierEdge) -> dict[str, object]:
    return {
        "source": _key_to_json(edge.source),
        "target": _key_to_json(edge.target),
        "source_representative": edge.source_representative,
        "target_representative": edge.target_representative,
        "orbit": edge.orbit,
        "member": edge.member,
        "target_retained": edge.target_retained,
    }


def _edge_from_json(record: dict[str, object]) -> AmbientCarrierEdge:
    return AmbientCarrierEdge(
        _key_from_json(record["source"]),
        _key_from_json(record["target"]),
        int(record["source_representative"]),
        int(record["target_representative"]),
        int(record["orbit"]),
        int(record["member"]),
        bool(record["target_retained"]),
    )


def _checkpoint_payload(
    atlas: ShimadaDisc11Atlas,
    start_key: PhysicalKey,
    processed: dict[
        PhysicalKey,
        tuple[ShimadaChamberState, ClosedLorentzianConeDecision],
    ],
    queue: deque[tuple[ShimadaChamberState, int]],
    edges: list[AmbientCarrierEdge],
    depths: dict[PhysicalKey, int],
    root_boundary_members: int,
) -> dict[str, object]:
    return {
        "schema_version": "0.1.0",
        "algorithm": _CHECKPOINT_ALGORITHM,
        "source_sha256": hashlib.sha256(atlas.path.read_bytes()).hexdigest(),
        "start_physical_key": _key_to_json(start_key),
        "processed": [
            _state_to_json(state, depths[key])
            for key, (state, _) in processed.items()
        ],
        "queue": [_state_to_json(state, depth) for state, depth in queue],
        "edges": [_edge_to_json(edge) for edge in edges],
        "root_boundary_wall_members": root_boundary_members,
    }


def _write_checkpoint(path: Path, payload: dict[str, object]) -> str:
    digest = _canonical_sha(payload)
    document = dict(payload)
    document["payload_sha256"] = digest
    temporary = path.with_name(path.name + ".tmp")
    temporary.parent.mkdir(parents=True, exist_ok=True)
    temporary.write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)
    return digest


def _read_checkpoint(
    atlas: ShimadaDisc11Atlas,
    path: Path,
    start_key: PhysicalKey,
) -> tuple[
    dict[PhysicalKey, tuple[ShimadaChamberState, ClosedLorentzianConeDecision]],
    deque[tuple[ShimadaChamberState, int]],
    set[PhysicalKey],
    dict[PhysicalKey, int],
    list[AmbientCarrierEdge],
    int,
    str,
]:
    document = json.loads(path.read_text())
    claimed_digest = document.pop("payload_sha256", None)
    actual_digest = _canonical_sha(document)
    if claimed_digest != actual_digest:
        raise ValueError("ambient-carrier checkpoint payload digest changed")
    if document.get("algorithm") != _CHECKPOINT_ALGORITHM:
        raise ValueError("ambient-carrier checkpoint algorithm changed")
    source_digest = hashlib.sha256(atlas.path.read_bytes()).hexdigest()
    if document.get("source_sha256") != source_digest:
        raise ValueError("ambient-carrier checkpoint belongs to another atlas")
    if _key_from_json(document["start_physical_key"]) != start_key:
        raise ValueError("ambient-carrier checkpoint has another base state")

    processed: dict[
        PhysicalKey,
        tuple[ShimadaChamberState, ClosedLorentzianConeDecision],
    ] = {}
    depths: dict[PhysicalKey, int] = {}
    for raw_record in document["processed"]:
        state, depth = _state_from_json(raw_record)
        key = _physical_state_key(atlas, state)
        decision = carrier_decision(atlas, state)
        if not decision.has_future_timelike_point:
            raise ValueError("checkpoint contains a rejected processed state")
        processed[key] = (state, decision)
        depths[key] = depth
    queue: deque[tuple[ShimadaChamberState, int]] = deque()
    for raw_record in document["queue"]:
        state, depth = _state_from_json(raw_record)
        key = _physical_state_key(atlas, state)
        if key in processed or key in depths:
            raise ValueError("checkpoint repeats a retained carrier state")
        queue.append((state, depth))
        depths[key] = depth
    discovered = set(depths)
    edges = [_edge_from_json(record) for record in document["edges"]]
    root_members = int(document["root_boundary_wall_members"])
    return (
        processed,
        queue,
        discovered,
        depths,
        edges,
        root_members,
        actual_digest,
    )


@lru_cache(maxsize=None)
def carrier_decision(
    atlas: ShimadaDisc11Atlas,
    state: ShimadaChamberState,
) -> ClosedLorentzianConeDecision:
    """Return the exact closed eigensection decision for one physical state."""

    return decide_lorentzian_closed_cone(
        _restricted_rows(atlas, state),
        quadratic_gram=_Q_GRAM,
        future_reference=AMPLE_FIELD_COEFFICIENTS,
    )


def iter_auxiliary_ambient_neighbors(
    atlas: ShimadaDisc11Atlas,
    state: ShimadaChamberState,
) -> Iterator[tuple[int, int, ShimadaChamberState]]:
    """Cross every physical member of every non-root ambient wall orbit."""

    chamber = atlas.chamber(state.index)
    for orbit, walls in enumerate(chamber.wall_orbits):
        if chamber.adjacencies[orbit].kind == "_minustwo":
            continue
        for member in range(len(walls)):
            neighbor = atlas.neighboring_state(state, orbit, member)
            if neighbor is None:
                raise ArithmeticError("an auxiliary ambient wall became a root wall")
            yield orbit, member, neighbor


def traverse_ambient_carrier(
    atlas: ShimadaDisc11Atlas,
    *,
    start: ShimadaChamberState | None = None,
    max_depth: int | None = None,
    max_states: int | None = None,
    progress: Progress | None = None,
    checkpoint_path: Path | None = None,
    resume: bool = True,
    checkpoint_every: int = 10,
) -> tuple[
    dict[PhysicalKey, tuple[ShimadaChamberState, ClosedLorentzianConeDecision]],
    tuple[tuple[ShimadaChamberState, int], ...],
    tuple[AmbientCarrierEdge, ...],
    dict[PhysicalKey, int],
    int,
    dict[str, object],
]:
    """Breadth-first traverse a bounded exact prefix of the ambient carrier."""

    if max_depth is None and max_states is None:
        raise ValueError("a depth or state bound is required")
    if max_depth is not None and max_depth < 0:
        raise ValueError("max_depth must be nonnegative")
    if max_states is not None and max_states < 1:
        raise ValueError("max_states must be positive")
    if checkpoint_every < 1:
        raise ValueError("checkpoint_every must be positive")

    base = start or atlas.representative_state(68)
    base_decision = carrier_decision(atlas, base)
    if not base_decision.has_future_timelike_point:
        raise ValueError("the requested carrier base state is not retained")
    base_key = _physical_state_key(atlas, base)
    resumed = bool(checkpoint_path and checkpoint_path.is_file() and resume)
    checkpoint_digest: str | None = None
    if resumed:
        assert checkpoint_path is not None
        (
            processed,
            queue,
            discovered,
            depths,
            edges,
            root_boundary_members,
            checkpoint_digest,
        ) = _read_checkpoint(atlas, checkpoint_path, base_key)
        if max_depth is not None and any(
            depths[key] > max_depth for key in processed
        ):
            raise ValueError("checkpoint already exceeds the requested depth")
        if max_states is not None and len(processed) > max_states:
            raise ValueError("checkpoint already exceeds the requested state bound")
    else:
        queue = deque([(base, 0)])
        discovered = {base_key}
        depths = {base_key: 0}
        processed = {}
        edges = []
        root_boundary_members = 0
    initial_processed = len(processed)

    while queue:
        if max_states is not None and len(processed) >= max_states:
            break
        if max_depth is not None and queue[0][1] > max_depth:
            break
        state, depth = queue.popleft()
        key = _physical_state_key(atlas, state)
        if key in processed:
            continue
        decision = carrier_decision(atlas, state)
        if not decision.has_future_timelike_point:
            raise ArithmeticError("a queued ambient carrier state became infeasible")
        if progress is not None:
            progress("start", len(processed) + 1, state, depth, len(queue))
        processed[key] = (state, decision)

        chamber = atlas.chamber(state.index)
        root_boundary_members += sum(
            len(chamber.wall_orbits[orbit])
            for orbit, adjacency in enumerate(chamber.adjacencies)
            if adjacency.kind == "_minustwo"
        )
        for orbit, member, neighbor in iter_auxiliary_ambient_neighbors(
            atlas, state
        ):
            neighbor_key = _physical_state_key(atlas, neighbor)
            neighbor_decision = carrier_decision(atlas, neighbor)
            retained = neighbor_decision.has_future_timelike_point
            edges.append(
                AmbientCarrierEdge(
                    key,
                    neighbor_key,
                    state.index,
                    neighbor.index,
                    orbit,
                    member,
                    retained,
                )
            )
            if retained and neighbor_key not in discovered:
                discovered.add(neighbor_key)
                depths[neighbor_key] = depth + 1
                queue.append((neighbor, depth + 1))
        if progress is not None:
            progress("finish", len(processed), state, depth, len(queue))
        if checkpoint_path is not None and len(processed) % checkpoint_every == 0:
            checkpoint_digest = _write_checkpoint(
                checkpoint_path,
                _checkpoint_payload(
                    atlas,
                    base_key,
                    processed,
                    queue,
                    edges,
                    depths,
                    root_boundary_members,
                ),
            )

    if checkpoint_path is not None:
        checkpoint_digest = _write_checkpoint(
            checkpoint_path,
            _checkpoint_payload(
                atlas,
                base_key,
                processed,
                queue,
                edges,
                depths,
                root_boundary_members,
            ),
        )
    return (
        processed,
        tuple(queue),
        tuple(edges),
        depths,
        root_boundary_members,
        {
            "enabled": checkpoint_path is not None,
            "path": str(checkpoint_path) if checkpoint_path is not None else None,
            "resumed": resumed,
            "initial_processed_states": initial_processed,
            "payload_sha256": checkpoint_digest,
            "algorithm": _CHECKPOINT_ALGORITHM,
            "checkpoint_every": checkpoint_every,
        },
    )


def _decision_record(
    key: PhysicalKey,
    state: ShimadaChamberState,
    decision: ClosedLorentzianConeDecision,
    depth: int,
) -> dict[str, object]:
    witness = decision.witness
    return {
        "physical_state_sha256": _canonical_sha(key),
        "representative": state.index,
        "row_transport_sha256": _canonical_sha(state.row_isometry),
        "ambient_carrier_depth": depth,
        "closed_section_reason": decision.reason,
        "phase_one_certificate_passes": decision.closed_linear_certificate.passes,
        "witness_sha256": (
            _canonical_sha([entry.to_json() for entry in witness])
            if witness is not None
            else None
        ),
        "maximum_quadratic_value": (
            decision.maximum_quadratic_value.to_json()
            if decision.maximum_quadratic_value is not None
            else None
        ),
    }


def carrier_centralizer_certificate(
    atlas: ShimadaDisc11Atlas,
    processed: dict[
        PhysicalKey,
        tuple[ShimadaChamberState, ClosedLorentzianConeDecision],
    ],
    queue: tuple[tuple[ShimadaChamberState, int], ...],
    depths: dict[PhysicalKey, int],
) -> dict[str, object]:
    """Audit all enhanced D68 returns in the complete carrier layers."""

    if queue:
        complete_radius = min(depth for _, depth in queue) - 1
        carrier_closed = False
    else:
        complete_radius = max(depths.values(), default=0)
        carrier_closed = True
    complete_keys = {
        key for key in processed if depths[key] <= complete_radius
    }
    tau = Matrix(conjugated_phi_action_shimada())
    eigenbasis = transported_identity_eigenbasis()
    base_stabilizer = tuple(Matrix(item) for item in atlas.stabilizer_group(68))
    tested = []
    centralizing: dict[tuple[tuple[int, ...], ...], dict[str, object]] = {}
    for key in complete_keys:
        state = processed[key][0]
        if state.index != 68:
            continue
        transport = Matrix(state.row_isometry)
        for stabilizer_index, stabilizer in enumerate(base_stabilizer):
            candidate = stabilizer.inv() * transport
            candidate_tuple = _integer_tuple(candidate)
            finite_order, exact_order, order_reason = _exact_matrix_order(candidate)
            commutes = candidate.T * tau == tau * candidate.T
            field_action_sha256 = None
            if commutes:
                field_action = _induced_eigen_action(candidate_tuple, eigenbasis)
                field_action_sha256 = _canonical_sha(
                    [
                        [entry.to_json() for entry in row]
                        for row in field_action
                    ]
                )
            record = {
                "ambient_carrier_distance": depths[key],
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
                "field_row_action_sha256": field_action_sha256,
            }
            tested.append(record)
            if commutes:
                centralizing.setdefault(candidate_tuple, record)

    tested.sort(
        key=lambda item: (
            item["ambient_carrier_distance"],
            item["physical_state_sha256"],
            item["target_stabilizer_index"],
        )
    )
    centralizing_records = sorted(
        centralizing.values(),
        key=lambda item: (
            item["ambient_carrier_distance"],
            item["matrix_sha256"],
        ),
    )
    infinite_order_found = any(
        not item["finite_order"] for item in centralizing_records
    )
    identity_field_sha = next(
        item["field_row_action_sha256"]
        for item in centralizing_records
        if item["exact_order"] == 1
    )
    known_v4_by_hash = {
        identity_field_sha: "identity",
        EXPECTED_FIELD_MATRIX_SHA256["A"]: "A",
        EXPECTED_FIELD_MATRIX_SHA256["B_plus"]: "B_plus",
        EXPECTED_FIELD_MATRIX_SHA256["B_minus"]: "B_minus",
    }
    centralizing_field_hashes = {
        item["field_row_action_sha256"] for item in centralizing_records
    }
    only_identity_and_a = centralizing_field_hashes == {
        identity_field_sha,
        EXPECTED_FIELD_MATRIX_SHA256["A"],
    }
    exactly_known_v4 = centralizing_field_hashes == set(known_v4_by_hash)
    missing_known_v4 = sorted(
        known_v4_by_hash[value]
        for value in set(known_v4_by_hash) - centralizing_field_hashes
    )
    if carrier_closed and missing_known_v4:
        raise ArithmeticError(
            "a closed carrier omitted a certified known V4 automorphism"
        )
    return {
        "complete_ambient_carrier_radius": complete_radius,
        "complete_layer_counts": {
            str(depth): count
            for depth, count in sorted(
                Counter(depths[key] for key in complete_keys).items()
            )
        },
        "states_in_complete_carrier_ball": len(complete_keys),
        "D68_representative_stabilizer_order": len(base_stabilizer),
        "D68_return_candidates_tested": len(tested),
        "D68_return_records": tested,
        "D68_return_records_sha256": _canonical_sha(tested),
        "centralizing_automorphisms_in_complete_carrier_ball": len(
            centralizing_records
        ),
        "centralizing_records": centralizing_records,
        "known_V4_elements_seen": sorted(
            known_v4_by_hash[value]
            for value in centralizing_field_hashes
            if value in known_v4_by_hash
        ),
        "known_V4_elements_not_seen": missing_known_v4,
        "missing_known_V4_elements_require_ambient_carrier_distance_at_least": (
            complete_radius + 1 if missing_known_v4 and not carrier_closed
            else None
        ),
        "only_identity_and_A_in_complete_carrier_ball": only_identity_and_a,
        "any_centralizer_element_other_than_identity_or_A_requires_ambient_carrier_distance_at_least": (
            complete_radius + 1 if only_identity_and_a and not carrier_closed
            else None
        ),
        "centralizing_set_in_complete_carrier_ball_is_exactly_known_V4": (
            exactly_known_v4
        ),
        "new_centralizer_beyond_known_V4_requires_ambient_carrier_distance_at_least": (
            complete_radius + 1 if exactly_known_v4 and not carrier_closed
            else None
        ),
        "carrier_queue_closed": carrier_closed,
        "infinite_order_return_found": infinite_order_found,
        "decision_status": (
            "infinite"
            if infinite_order_found
            else "finite"
            if carrier_closed
            else "pending_larger_exact_ambient_carrier_ball"
        ),
    }


def analyze_ambient_carrier(
    path: Path,
    *,
    max_depth: int | None = 1,
    max_states: int | None = None,
    progress: Progress | None = None,
    checkpoint_path: Path | None = None,
    resume: bool = True,
    checkpoint_every: int = 10,
) -> dict[str, object]:
    """Return a compact certificate for a bounded ambient-carrier prefix."""

    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != DISC11_SOURCE_SHA256:
        raise ValueError("unexpected disc11.txt checksum")
    atlas = ShimadaDisc11Atlas(path)
    (
        processed,
        queue,
        edges,
        depths,
        root_members,
        checkpoint,
    ) = traverse_ambient_carrier(
        atlas,
        max_depth=max_depth,
        max_states=max_states,
        progress=progress,
        checkpoint_path=checkpoint_path,
        resume=resume,
        checkpoint_every=checkpoint_every,
    )
    centralizer = carrier_centralizer_certificate(
        atlas, processed, queue, depths
    )
    state_records = sorted(
        (
            _decision_record(key, state, decision, depths[key])
            for key, (state, decision) in processed.items()
        ),
        key=lambda item: (
            item["ambient_carrier_depth"],
            item["physical_state_sha256"],
        ),
    )
    queue_records = sorted(
        (
            {
                "physical_state_sha256": _canonical_sha(
                    _physical_state_key(atlas, state)
                ),
                "representative": state.index,
                "ambient_carrier_depth": depth,
            }
            for state, depth in queue
        ),
        key=lambda item: (
            item["ambient_carrier_depth"],
            item["physical_state_sha256"],
        ),
    )
    edge_records = sorted(
        (
            {
                "source_physical_state_sha256": _canonical_sha(edge.source),
                "target_physical_state_sha256": _canonical_sha(edge.target),
                "source_representative": edge.source_representative,
                "target_representative": edge.target_representative,
                "ambient_wall": [edge.orbit, edge.member],
                "target_retained": edge.target_retained,
            }
            for edge in edges
        ),
        key=lambda item: (
            item["source_physical_state_sha256"],
            item["ambient_wall"],
        ),
    )
    traversal_closed = not queue
    return {
        "analysis_id": "shimada_disc11_exact_ambient_carrier",
        "verified_exact_closed_section_decisions": all(
            decision.closed_linear_certificate.passes
            for _, decision in processed.values()
        ),
        "source_sha256": digest,
        "atlas_chambers": atlas.chamber_count,
        "requested_max_depth": max_depth,
        "requested_max_states": max_states,
        "processed_carrier_states": len(processed),
        "queued_carrier_states": len(queue),
        "processed_depth_counts": {
            str(depth): count
            for depth, count in sorted(
                Counter(depths[key] for key in processed).items()
            )
        },
        "auxiliary_ambient_wall_crossings_tested": len(edges),
        "retained_crossings": sum(edge.target_retained for edge in edges),
        "rejected_crossings": sum(not edge.target_retained for edge in edges),
        "root_boundary_wall_members_not_crossed": root_members,
        "carrier_traversal_closed": traversal_closed,
        "checkpoint": checkpoint,
        "centralizer_decision_ball": centralizer,
        "state_records_sha256": _canonical_sha(state_records),
        "queue_records_sha256": _canonical_sha(queue_records),
        "edge_records_sha256": _canonical_sha(edge_records),
        "states": state_records,
        "queue": queue_records,
        "edges": edge_records,
        "scope": (
            "exact_bounded_ambient_carrier_prefix;_a_nonclosed_prefix_is_not_"
            "a_global_centralizer_or_finite_volume_decision"
        ),
    }


def certify_d35_carrier_regression(
    path: Path,
    fixture_path: Path = DEFAULT_FIXTURE,
) -> dict[str, object]:
    """Verify that the ambient carrier contains both D35 opposite galleries."""

    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    fixture = json.loads(fixture_path.read_text())
    if digest != DISC11_SOURCE_SHA256 or fixture["source_sha256"] != digest:
        raise ValueError("unexpected disc11.txt checksum")
    atlas = ShimadaDisc11Atlas(path)
    state = ShimadaChamberState(
        int(fixture["state_index"]),
        tuple(
            tuple(int(entry) for entry in row)
            for row in fixture["row_isometry"]
        ),
    )
    interior = tuple(Golden.from_json(entry) for entry in fixture["interior"])
    wall_class = tuple(
        (int(label[0]), int(label[1])) for label in fixture["wall_class"]
    )

    records = []
    target_hashes = set()
    for expected in fixture["witnesses"]:
        witness = tuple(Golden.from_json(entry) for entry in expected["witness"])
        neighbor, _, gallery, _ = _cross_internal_facet(
            atlas, state, interior, witness, wall_class
        )
        current = state
        states = [current]
        crossing_records = []
        for representative, orbit, member in gallery:
            if current.index != representative:
                raise AssertionError("the D35 ambient gallery source changed")
            available = {
                (candidate_orbit, candidate_member): candidate
                for candidate_orbit, candidate_member, candidate
                in iter_auxiliary_ambient_neighbors(atlas, current)
            }
            if (orbit, member) not in available:
                raise AssertionError("a D35 gallery wall left the auxiliary graph")
            current = available[(orbit, member)]
            crossing_records.append(
                {
                    "source_representative": representative,
                    "ambient_wall": [orbit, member],
                    "target_representative": current.index,
                }
            )
            states.append(current)
        if _physical_state_key(atlas, current) != _physical_state_key(
            atlas, neighbor
        ):
            raise AssertionError("the D35 ambient gallery target changed")

        decisions = [carrier_decision(atlas, item) for item in states]
        state_hashes = [
            _canonical_sha(_physical_state_key(atlas, item)) for item in states
        ]
        target_hash = state_hashes[-1]
        target_hashes.add(target_hash)
        records.append(
            {
                "representatives": [item.index for item in states],
                "physical_state_sha256": state_hashes,
                "ambient_crossings": crossing_records,
                "closed_section_reasons": [item.reason for item in decisions],
                "all_states_retained": all(
                    item.has_future_timelike_point for item in decisions
                ),
                "target_physical_state_sha256": target_hash,
            }
        )
        if target_hash != expected["target_physical_state_sha256"]:
            raise AssertionError("a D35 carrier target digest changed")

    verified = (
        len(records) == 2
        and len(target_hashes) == 2
        and all(record["all_states_retained"] for record in records)
    )
    return {
        "analysis_id": "shimada_disc11_D35_ambient_carrier_regression",
        "verified": verified,
        "source_sha256": digest,
        "source_physical_state_sha256": fixture["physical_state_sha256"],
        "wall_class": [list(label) for label in wall_class],
        "distinct_opposite_carrier_targets": len(target_hashes),
        "galleries": records,
        "conclusion": (
            "crossing_every_auxiliary_ambient_wall_retains_both_subdivision_"
            "branches_missed_by_a_one_witness_restricted_facet_graph"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("disc11", type=Path)
    parser.add_argument("--max-depth", type=int, default=1)
    parser.add_argument("--max-states", type=int)
    parser.add_argument("--progress-every", type=int)
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--checkpoint-every", type=int, default=10)
    parser.add_argument(
        "--fresh-checkpoint",
        action="store_true",
        help="ignore and replace an existing checkpoint",
    )
    parser.add_argument("--d35-regression-only", action="store_true")
    parser.add_argument("--d35-fixture", type=Path, default=DEFAULT_FIXTURE)
    args = parser.parse_args()
    if args.progress_every is not None and args.progress_every < 1:
        parser.error("--progress-every must be positive")
    if args.checkpoint_every < 1:
        parser.error("--checkpoint-every must be positive")
    if args.fresh_checkpoint and args.checkpoint is None:
        parser.error("--fresh-checkpoint requires --checkpoint")

    def progress(event, count, state, depth, queued):
        if args.progress_every is None or count % args.progress_every:
            return
        print(
            f"carrier {event} state={count} depth={depth} "
            f"representative=D{state.index} queued={queued}",
            file=sys.stderr,
            flush=True,
        )

    result = (
        certify_d35_carrier_regression(args.disc11, args.d35_fixture)
        if args.d35_regression_only
        else analyze_ambient_carrier(
            args.disc11,
            max_depth=args.max_depth,
            max_states=args.max_states,
            progress=progress if args.progress_every is not None else None,
            checkpoint_path=args.checkpoint,
            resume=not args.fresh_checkpoint,
            checkpoint_every=args.checkpoint_every,
        )
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
