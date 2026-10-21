"""Exact global traversal search for the discriminant-11 K3 eigen-nef cell.

Shimada's induced Borcherds chambers subdivide the nef cone.  Their non-root
walls are auxiliary: the actual eigen-nef chamber is obtained by gluing every
induced chamber met by the fixed golden eigenspace until a ``-2`` wall is
reached.  This script performs that gluing as a breadth-first traversal.

Face feasibility and cone membership are decided by a finite phase-I simplex
implemented directly over ``Q(phi)``.  If a strict linear face witness is not
timelike, an exact Lorentzian active-set calculation decides whether the face
reaches the positive cone.  A chamber enters the certified graph only after
all of the following have been checked over ``Q(phi)``:

* the complete positive-proportionality partition of its primitive walls;
* a timelike strict face point for every asserted facet;
* a positive exact dependence on facet rows for every redundant class; and
* a strict timelike interior point in every chamber reached across an
  auxiliary facet.

Consequently a completed run is an exact certificate relative to the
SHA-pinned Shimada atlas.  A state limit is reported as incomplete, never as
nonexistence or nonreflectivity.
"""

from __future__ import annotations

import argparse
from collections import deque
from dataclasses import dataclass
from functools import lru_cache
import hashlib
import json
from pathlib import Path
from typing import Mapping

from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.geometric.golden_halfspaces import (
    DependenceTerms,
    GoldenHalfspaceCertificateAudit,
    LorentzianFaceDecision,
    PositiveProportionalClass,
    audit_golden_halfspace_certificate,
    decide_lorentzian_strict_face,
    find_nonnegative_span_phase_one,
    golden_dot,
    golden_quadratic_value,
    inverse_golden_matrix,
    positive_proportional_classes,
)
from coxeter_polytopes.geometric.k3 import K3_TRACE_SCALAR
from scripts.analyze_shimada_disc11_eigen_nef import _restricted_coefficients
from scripts.coupled_discriminant11_candidate import field_gram
from scripts.shimada_disc11_atlas import ShimadaChamberState, ShimadaDisc11Atlas
from scripts.shimada_disc11_conjugate import (
    AMPLE_FIELD_COEFFICIENTS,
    DISC11_SOURCE_SHA256,
    transported_identity_eigenbasis,
)


WallLabel = tuple[int, int]
WallClass = tuple[WallLabel, ...]
GoldenVector = tuple[Golden, ...]


_Q_GRAM = tuple(
    tuple(K3_TRACE_SCALAR * entry for entry in row) for row in field_gram()
)
_Q_GRAM_INVERSE = inverse_golden_matrix(_Q_GRAM)


@lru_cache(maxsize=None)
def _restricted_rows(
    atlas: ShimadaDisc11Atlas,
    state: ShimadaChamberState,
) -> dict[WallLabel, GoldenVector]:
    basis = transported_identity_eigenbasis()
    return {
        (orbit_index, member_index): _restricted_coefficients(wall, basis)
        for orbit_index, orbit in enumerate(atlas.state_wall_orbits(state))
        for member_index, wall in enumerate(orbit)
    }


@lru_cache(maxsize=None)
def _physical_state_key(
    atlas: ShimadaDisc11Atlas,
    state: ShimadaChamberState,
) -> tuple[tuple[int, ...], ...]:
    """Canonicalize a physical chamber by its complete primitive wall set."""

    return tuple(
        sorted(
            wall
            for orbit in atlas.state_wall_orbits(state)
            for wall in orbit
        )
    )


def _state_digest(key: tuple[tuple[int, ...], ...]) -> str:
    payload = json.dumps(key, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def _ray_key(row: GoldenVector) -> GoldenVector:
    # Normalize by the fixed ample class, not by an arbitrary coordinate:
    # dividing by a negative pivot would reverse the oriented half-space and
    # corrupt the Coxeter Gram signs.  Every effective wall is strictly
    # positive on this ample point.
    scale = golden_dot(row, AMPLE_FIELD_COEFFICIENTS)
    if scale.sign() <= 0:
        raise ValueError("a visible effective wall is not positive on the ample class")
    return tuple(entry / scale for entry in row)


def _lorentz_pair(left: GoldenVector, right: GoldenVector) -> Golden:
    return sum(
        (
            left[i] * _Q_GRAM[i][j] * right[j]
            for i in range(len(left))
            for j in range(len(right))
        ),
        Golden(0),
    )


def _face_is_exact(
    rows: Mapping[WallLabel, GoldenVector],
    wall_class: WallClass,
    point: GoldenVector,
    interior: GoldenVector,
) -> bool:
    target = set(wall_class)
    time_pairing = _lorentz_pair(interior, point)
    return (
        all(golden_dot(rows[label], point) == 0 for label in target)
        and all(
            golden_dot(row, point).sign() > 0
            for label, row in rows.items()
            if label not in target
        )
        and golden_quadratic_value(_Q_GRAM, point).sign() > 0
        and time_pairing.sign() > 0
    )


def _project_to_wall(
    interior: GoldenVector, functional: GoldenVector
) -> GoldenVector | None:
    """Exact Lorentzian wall projection using the module-cached inverse."""

    normal = tuple(
        sum(
            (
                _Q_GRAM_INVERSE[row][column] * functional[column]
                for column in range(len(functional))
            ),
            Golden(0),
        )
        for row in range(len(functional))
    )
    normal_square = golden_dot(functional, normal)
    if not normal_square:
        return None
    scale = golden_dot(functional, interior) / normal_square
    projected = tuple(
        interior[index] - scale * normal[index]
        for index in range(len(interior))
    )
    if golden_dot(functional, projected):
        raise ArithmeticError("cached exact wall projection failed")
    return projected


def _discover_face_witness(
    rows: Mapping[WallLabel, GoldenVector],
    wall_class: WallClass,
    interior: GoldenVector,
) -> LorentzianFaceDecision:
    target = rows[wall_class[0]]
    target_members = set(wall_class)
    other_rows = {
        label: row for label, row in rows.items() if label not in target_members
    }
    return decide_lorentzian_strict_face(
        target,
        other_rows,
        quadratic_gram=_Q_GRAM,
        interior=interior,
    )


def _exact_redundancy(
    target: GoldenVector,
    generators: Mapping[WallLabel, GoldenVector],
) -> DependenceTerms | None:
    dependence, certificate = find_nonnegative_span_phase_one(target, generators)
    if not certificate.passes:
        raise ArithmeticError("exact redundancy phase-I certificate failed")
    return dependence


@dataclass(frozen=True, slots=True)
class ExactSectionReduction:
    rows: Mapping[WallLabel, GoldenVector]
    classes: tuple[PositiveProportionalClass, ...]
    linear_facet_classes: tuple[WallClass, ...]
    facet_classes: tuple[WallClass, ...]
    witnesses: Mapping[WallLabel, GoldenVector]
    redundancies: Mapping[WallLabel, DependenceTerms]
    causal_inactivities: Mapping[WallLabel, LorentzianFaceDecision]
    audit: GoldenHalfspaceCertificateAudit

    @property
    def root_facets(self) -> tuple[WallClass, ...]:
        raise AttributeError("root classification depends on the Shimada chamber")


def _reduce_section(
    atlas: ShimadaDisc11Atlas,
    state: ShimadaChamberState,
    interior: GoldenVector,
) -> ExactSectionReduction:
    rows = _restricted_rows(atlas, state)
    classes = positive_proportional_classes(rows)
    class_rows = {
        item.members[0]: rows[item.members[0]] for item in classes
    }
    linear_facets: list[WallClass] = []

    # First compute the extreme rays of the dual cone.  A class representative
    # is a linearly redundant inequality exactly when it lies in the
    # nonnegative span of all the other representatives.  Only the remaining
    # extreme rays can support a facet after intersecting with the future
    # Lorentz cone.  This exact phase-I pre-reduction is decisive for large
    # Shimada chambers: it prevents one expensive Lorentz/KKT decision for
    # every manifestly non-extreme row.
    for item in classes:
        wall_class = tuple(item.members)
        target = wall_class[0]
        other_rows = {
            label: row for label, row in class_rows.items() if label != target
        }
        if _exact_redundancy(rows[target], other_rows) is None:
            linear_facets.append(wall_class)

    facets: list[WallClass] = []
    witnesses: dict[WallLabel, GoldenVector] = {}
    unresolved: list[WallClass] = []

    # Exact projections are a particularly cheap certificate on the Shimada
    # atlas: when one is timelike and strict it already proves the facet.  No
    # completeness rests on this fast path; unresolved classes enter the
    # finite phase-I/KKT decision below.
    for wall_class in linear_facets:
        projected = _project_to_wall(interior, rows[wall_class[0]])
        if projected is not None and _face_is_exact(
            rows, wall_class, projected, interior
        ):
            facets.append(wall_class)
            witnesses[wall_class[0]] = projected
        else:
            unresolved.append(wall_class)

    linear_rows = {
        wall_class[0]: rows[wall_class[0]] for wall_class in linear_facets
    }
    for wall_class in unresolved:
        target = wall_class[0]
        decision = decide_lorentzian_strict_face(
            rows[target],
            {
                label: row
                for label, row in linear_rows.items()
                if label != target
            },
            quadratic_gram=_Q_GRAM,
            interior=interior,
        )
        if decision.is_facet:
            assert decision.witness is not None
            if not _face_is_exact(
                rows, wall_class, decision.witness, interior
            ):
                raise ArithmeticError(
                    "a reduced Lorentzian witness failed the full wall audit"
                )
            facets.append(wall_class)
            witnesses[target] = decision.witness

    # Reconstruct the final certificate only from actual timelike facets.
    # Usually every remaining row is now a direct positive span.  If a
    # linearly exposed facet is causally inactive, recompute that exceptional
    # decision against the full primitive system so the generic independent
    # audit below remains unchanged.
    facet_generators = {
        wall_class[0]: rows[wall_class[0]] for wall_class in facets
    }
    facet_sets = {frozenset(wall_class) for wall_class in facets}
    redundancies: dict[WallLabel, DependenceTerms] = {}
    causal_inactivities: dict[WallLabel, LorentzianFaceDecision] = {}
    for item in classes:
        if frozenset(item.members) in facet_sets:
            continue
        wall_class = tuple(item.members)
        target = wall_class[0]
        dependence = _exact_redundancy(rows[target], facet_generators)
        if dependence is not None:
            redundancies[target] = dependence
            continue
        decision = _discover_face_witness(rows, wall_class, interior)
        if decision.is_facet:
            raise ArithmeticError(
                "dual-cone pre-reduction omitted a timelike facet"
            )
        causal_inactivities[target] = decision

    audit = audit_golden_halfspace_certificate(
        rows,
        facet_classes=facets,
        facet_witnesses=witnesses,
        redundancies=redundancies,
        causal_inactivities=causal_inactivities,
        interior=interior,
        quadratic_gram=_Q_GRAM,
    )
    if not audit.passes:
        raise RuntimeError(f"state D{state.index}: exact half-space audit failed")
    return ExactSectionReduction(
        rows,
        classes,
        tuple(linear_facets),
        tuple(facets),
        witnesses,
        redundancies,
        causal_inactivities,
        audit,
    )


def _root_and_internal_facets(
    atlas: ShimadaDisc11Atlas,
    state: ShimadaChamberState,
    reduction: ExactSectionReduction,
) -> tuple[tuple[WallClass, ...], tuple[WallClass, ...]]:
    chamber = atlas.chamber(state.index)
    roots = []
    internals = []
    for wall_class in reduction.facet_classes:
        kinds = {chamber.adjacencies[orbit].kind for orbit, _ in wall_class}
        if "_minustwo" in kinds:
            roots.append(wall_class)
        else:
            internals.append(wall_class)
    return tuple(roots), tuple(internals)


def _strict_state_record(
    atlas: ShimadaDisc11Atlas,
    state: ShimadaChamberState,
    point: GoldenVector,
) -> tuple[bool, int, int]:
    rows = _restricted_rows(atlas, state)
    signs = tuple(golden_dot(row, point).sign() for row in rows.values())
    return (
        all(sign > 0 for sign in signs)
        and golden_quadratic_value(_Q_GRAM, point).sign() > 0
        and _lorentz_pair(tuple(AMPLE_FIELD_COEFFICIENTS), point).sign() > 0,
        sum(sign < 0 for sign in signs),
        sum(sign == 0 for sign in signs),
    )


def _locate_opposite_chamber(
    atlas: ShimadaDisc11Atlas,
    start: ShimadaChamberState,
    point: GoldenVector,
) -> tuple[ShimadaChamberState, tuple[tuple[int, int, int], ...]]:
    """Walk separating auxiliary walls until ``point`` is strictly interior.

    There is no step cutoff.  The segment from the starting interior to the
    fixed timelike target is compact, Shimada's induced arrangement is locally
    finite, and each crossed negative wall removes one separating wall.  Thus
    the exact walk is finite.  A repeated physical state is an arithmetic or
    atlas inconsistency and is rejected rather than reported as a timeout.
    """

    current = start
    path = []
    seen = set()
    while True:
        strict, negatives, zeros = _strict_state_record(atlas, current, point)
        if strict:
            return current, tuple(path)
        if zeros and not negatives:
            raise RuntimeError("the proposed opposite point lies on an auxiliary wall")
        key = _physical_state_key(atlas, current)
        if key in seen:
            raise RuntimeError("the violated-wall gallery cycled")
        seen.add(key)
        rows = _restricted_rows(atlas, current)
        negative_labels = [
            label for label, row in rows.items() if golden_dot(row, point).sign() < 0
        ]
        if not negative_labels:
            raise RuntimeError("point location stalled without a negative wall")
        orbit, member = negative_labels[0]
        neighbor = atlas.neighboring_state(current, orbit, member)
        if neighbor is None:
            raise RuntimeError("an auxiliary crossing unexpectedly hit a -2 boundary")
        path.append((current.index, orbit, member))
        current = neighbor


def _cross_internal_facet(
    atlas: ShimadaDisc11Atlas,
    state: ShimadaChamberState,
    interior: GoldenVector,
    witness: GoldenVector,
    wall_class: WallClass,
) -> tuple[ShimadaChamberState, GoldenVector, tuple[tuple[int, int, int], ...], int]:
    # The witness is in the relative interior of the selected facet.  Local
    # finiteness gives an open transverse interval meeting no unrelated wall,
    # so the unbounded dyadic sequence eventually lies in that interval.
    denominator = 2
    while True:
        epsilon = Golden(1) / denominator
        opposite = tuple(
            witness[index] - epsilon * (interior[index] - witness[index])
            for index in range(len(interior))
        )
        if (
            golden_quadratic_value(_Q_GRAM, opposite).sign() <= 0
            or _lorentz_pair(
                tuple(AMPLE_FIELD_COEFFICIENTS), opposite
            ).sign() <= 0
        ):
            denominator *= 2
            continue
        current_rows = _restricted_rows(atlas, state)
        negative_labels = {
            label
            for label, row in current_rows.items()
            if golden_dot(row, opposite).sign() < 0
        }
        # Stay in the open star of the selected facet.  Without this check a
        # coarse perturbation can cross several unrelated facets and skip the
        # immediate induced chambers needed for a completeness proof.
        if negative_labels != set(wall_class):
            denominator *= 2
            continue
        try:
            neighbor, gallery = _locate_opposite_chamber(atlas, state, opposite)
        except RuntimeError:
            denominator *= 2
            continue
        if _physical_state_key(atlas, neighbor) == _physical_state_key(atlas, state):
            # This can be a genuine enhanced-surface reflection, as at D68,
            # but the physical chambers on the two sides must still differ.
            if neighbor.row_isometry == state.row_isometry:
                denominator *= 2
                continue
        return neighbor, opposite, gallery, denominator


def explore_global_cell(path: Path, *, max_states: int = 200) -> dict[str, object]:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != DISC11_SOURCE_SHA256:
        raise ValueError("unexpected disc11.txt checksum")
    if max_states < 1:
        raise ValueError("max_states must be positive")

    atlas = ShimadaDisc11Atlas(path)
    start = atlas.representative_state(68)
    start_key = _physical_state_key(atlas, start)
    queue = deque([(start, tuple(AMPLE_FIELD_COEFFICIENTS))])
    queued = {start_key}
    processed: dict[tuple[tuple[int, ...], ...], dict[str, object]] = {}
    visible_root_rays: dict[GoldenVector, list[dict[str, object]]] = {}
    edges = []
    state_limit_hit = False

    while queue:
        if len(processed) >= max_states:
            state_limit_hit = True
            break
        state, interior = queue.popleft()
        key = _physical_state_key(atlas, state)
        if key in processed:
            continue
        reduction = _reduce_section(atlas, state, interior)
        root_facets, internal_facets = _root_and_internal_facets(
            atlas, state, reduction
        )
        digest_key = _state_digest(key)
        processed[key] = {
            "representative": state.index,
            "physical_state_sha256": digest_key,
            "primitive_walls": len(reduction.rows),
            "restricted_classes": len(reduction.classes),
            "linear_facet_classes": len(reduction.linear_facet_classes),
            "facets": len(reduction.facet_classes),
            "root_facets": len(root_facets),
            "internal_facets": len(internal_facets),
            "causal_inactive_classes": len(reduction.causal_inactivities),
            "exact_halfspace_audit": reduction.audit.passes,
            "feasibility_engine": (
                "exact_Q(phi)_phase_I_with_Lorentzian_KKT_fallback"
            ),
        }

        for wall_class in root_facets:
            ray = _ray_key(reduction.rows[wall_class[0]])
            visible_root_rays.setdefault(ray, []).append(
                {
                    "state": digest_key,
                    "representative": state.index,
                    "wall_class": [list(label) for label in wall_class],
                }
            )

        for wall_class in internal_facets:
            witness = reduction.witnesses[wall_class[0]]
            neighbor, opposite, gallery, denominator = _cross_internal_facet(
                atlas, state, interior, witness, wall_class
            )
            neighbor_key = _physical_state_key(atlas, neighbor)
            neighbor_digest = _state_digest(neighbor_key)
            edges.append(
                {
                    "from": digest_key,
                    "from_representative": state.index,
                    "wall_class": [list(label) for label in wall_class],
                    "gallery": [
                        {
                            "representative": representative,
                            "orbit": orbit,
                            "member": member,
                        }
                        for representative, orbit, member in gallery
                    ],
                    "crossing_denominator": denominator,
                    "to": neighbor_digest,
                    "to_representative": neighbor.index,
                }
            )
            if neighbor_key not in queued and neighbor_key not in processed:
                queued.add(neighbor_key)
                queue.append((neighbor, opposite))

    complete = not queue and not state_limit_hit
    return {
        "verified_exact_local_reductions": all(
            record["exact_halfspace_audit"] for record in processed.values()
        ),
        "bounded_discovery_cutoffs": False,
        "global_traversal_complete": complete,
        "source_sha256": digest,
        "atlas_chambers": atlas.chamber_count,
        "max_states": max_states,
        "state_limit_hit": state_limit_hit,
        "processed_physical_states": len(processed),
        "queued_unprocessed_states": len(queue),
        "unique_visible_root_rays": len(visible_root_rays),
        "visible_root_rows": [
            {
                "normalized_restricted_functional": [
                    entry.to_json() for entry in ray
                ],
                "occurrence_count": len(occurrences),
                "first_occurrence": occurrences[0],
            }
            for ray, occurrences in visible_root_rays.items()
        ],
        "internal_crossings": len(edges),
        "states": list(processed.values()),
        "edges": edges,
        "scope": (
            "exact_global_eigen_nef_wall_completeness_relative_to_Shimada_atlas"
            if complete
            else "partial_exact_traversal_no_global_volume_or_compactness_claim"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("disc11", type=Path)
    parser.add_argument("--max-states", type=int, default=200)
    args = parser.parse_args()
    print(
        json.dumps(
            explore_global_cell(args.disc11, max_states=args.max_states),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
