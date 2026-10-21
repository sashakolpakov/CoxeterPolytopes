"""Unrestricted compact-label obstruction for record 43's second nerve.

The whole-shadow incidence census leaves two exact 90-vertex nerves.  The
first is the previously frozen multiwedge selection.  This module treats the
other one.  Two of its order-five minimal nonfaces are disjoint Lanner
supports.  Each support has one distinguished vertex and four interchangeable
graph isolates.  We verify the resulting ``S4 x S4`` subgroup on the entire
minimal-nonface clutter, reduce all 420 labelled placements on each support
to 20 exact orbit representatives, and enumerate every one of the 25 cross
pairs in the full Burcroff alphabet ``{2,3,4,5}``.

Local ellipticity leaves 55 exceptional label-three orbit blocks.  Labels
four and five only decrease the Gram entries and are rejected exactly; no two
exceptional edges can coexist.  The 55 blocks lie in two weighted-graph
classes, whose exact principal minors give inertia ``(8,2,0)``.  Interlacing
therefore excludes an ambient compact H9 Gram of inertia ``(9,1,4)``.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations, permutations

from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import Pair, pair
from coxeter_polytopes.combinatorics.h9_14_disjoint_lanner_gate import (
    DisjointLannerCrossGateAudit,
    audit_disjoint_lanner_cross_gate,
)
from coxeter_polytopes.combinatorics.h9_14_record43 import (
    TARGET_MISSING_EDGES,
)
from coxeter_polytopes.combinatorics.h9_14_record43_all_completions import (
    SECOND_TARGET_MISSING_FIVES,
)
from coxeter_polytopes.combinatorics.h9_14_selected_nerve_coxeter_no_go import (
    ExactIndexTwoWitness,
    _a5_witnesses,
    _all_lanner_states,
    _canonical_lanner_states,
    _double_edge_gate,
    _isomorphic_blocks,
    _labelled_block_graph,
    _stronger_single_edge_survivors,
    _vertices_from_nonfaces,
)


def _distinguished_support_states(
    support: tuple[int, ...],
) -> tuple[tuple[dict[Pair, int], ...], tuple[tuple[int, int], ...]]:
    """Quotient all 420 A5 placements by the S4 fixing ``support[0]``."""

    distinguished = support[0]
    interchangeable = support[1:]
    edges = tuple(combinations(sorted(support), 2))
    orbit_counts: Counter[tuple[int, ...]] = Counter()
    for state in _all_lanner_states(support):
        orbit_codes = []
        for image_order in permutations(interchangeable):
            action = {
                distinguished: distinguished,
                **dict(zip(interchangeable, image_order, strict=True)),
            }
            inverse = {image: source for source, image in action.items()}
            orbit_codes.append(
                tuple(state[pair(inverse[left], inverse[right])] for left, right in edges)
            )
        orbit_counts[min(orbit_codes)] += 1
    representatives = tuple(
        dict(zip(edges, code, strict=True)) for code in sorted(orbit_counts)
    )
    return representatives, tuple(sorted(Counter(orbit_counts.values()).items()))


def _nerve_automorphism_subgroup_verified() -> bool:
    """Check all 24^2 permutations preserve the full ten-member clutter."""

    left_isolates = SECOND_TARGET_MISSING_FIVES[0][1:]
    right_isolates = SECOND_TARGET_MISSING_FIVES[1][1:]
    clutter = {
        frozenset(face)
        for face in (*TARGET_MISSING_EDGES, *SECOND_TARGET_MISSING_FIVES)
    }
    for left_order in permutations(left_isolates):
        for right_order in permutations(right_isolates):
            action = dict(zip(left_isolates, left_order, strict=True))
            action.update(dict(zip(right_isolates, right_order, strict=True)))
            image = {
                frozenset(action.get(vertex, vertex) for vertex in face)
                for face in clutter
            }
            if image != clutter:
                return False
    return True


@dataclass(frozen=True, slots=True)
class Record43SecondNerveLabelAudit:
    vertex_count: int
    left_support: tuple[int, ...]
    right_support: tuple[int, ...]
    verified_nerve_automorphism_count: int
    full_left_placement_count: int
    full_right_placement_count: int
    left_orbit_size_histogram: tuple[tuple[int, int], ...]
    right_orbit_size_histogram: tuple[tuple[int, int], ...]
    cross_gate: DisjointLannerCrossGateAudit
    stronger_single_edge_survivors: tuple[tuple[int, int], ...]
    raw_double_edge_trial_count: int
    monotonicity_residual_double_trial_count: int
    double_edge_passing_count: int
    exceptional_class_histogram: tuple[tuple[str, int], ...]
    exact_index_two_witnesses: tuple[ExactIndexTwoWitness, ...]

    @property
    def unrestricted_compact_coxeter_no_go(self) -> bool:
        return (
            self.vertex_count == 90
            and self.left_support == SECOND_TARGET_MISSING_FIVES[0]
            and self.right_support == SECOND_TARGET_MISSING_FIVES[1]
            and self.verified_nerve_automorphism_count == 576
            and self.full_left_placement_count == 420
            and self.full_right_placement_count == 420
            and self.left_orbit_size_histogram == ((12, 5), (24, 15))
            and self.right_orbit_size_histogram == ((12, 5), (24, 15))
            and self.cross_gate.left_state_count == 20
            and self.cross_gate.right_state_count == 20
            and self.cross_gate.sparse_state_count == 400
            and len(self.cross_gate.cross_pairs) == 25
            and self.cross_gate.unique_actual_restriction_size_histogram
            == ((6, 16), (7, 24), (8, 25))
            and self.cross_gate.label_three_test_count == 10_000
            and len(self.cross_gate.passing_trials) == 55
            and tuple(
                sorted(
                    Counter(
                        count
                        for _edge, count in (
                            self.cross_gate.passing_count_by_cross_pair
                        )
                    ).items()
                )
            )
            == ((0, 4), (1, 9), (2, 2), (3, 7), (6, 2), (9, 1))
            and self.stronger_single_edge_survivors == ((4, 0), (5, 0))
            and self.raw_double_edge_trial_count == 120_000
            and self.monotonicity_residual_double_trial_count == 0
            and self.double_edge_passing_count == 0
            and self.exceptional_class_histogram
            == (("A5_A5_type_1_1", 25), ("A5_A5_type_4_1", 30))
            and tuple(witness.inertia for witness in self.exact_index_two_witnesses)
            == ((8, 2, 0), (8, 2, 0))
        )

    def to_record(self) -> dict[str, object]:
        return {
            "rank_nine_vertex_count": self.vertex_count,
            "second_target_missing_fives": [
                list(face) for face in SECOND_TARGET_MISSING_FIVES
            ],
            "left_support": list(self.left_support),
            "right_support": list(self.right_support),
            "verified_nerve_automorphism_subgroup": (
                "S4 on {3,4,5,6} x S4 on {7,8,9,10}"
            ),
            "verified_nerve_automorphism_count": (
                self.verified_nerve_automorphism_count
            ),
            "full_left_placement_count": self.full_left_placement_count,
            "full_right_placement_count": self.full_right_placement_count,
            "left_orbit_size_histogram": {
                str(size): count for size, count in self.left_orbit_size_histogram
            },
            "right_orbit_size_histogram": {
                str(size): count for size, count in self.right_orbit_size_histogram
            },
            "orbit_coverage_argument": (
                "every full placement pair is carried to one of the 20 x 20 "
                "representative pairs by the verified independent S4 "
                "actions; all 25 cross edges and all C(25,2) edge pairs are "
                "then enumerated, so distinguished support vertices are not "
                "silently permuted with isolates"
            ),
            "support_gate": self.cross_gate.to_record(),
            "finite_label_domain": [2, 3, 4, 5],
            "stronger_single_edge_survivors": {
                str(label): count
                for label, count in self.stronger_single_edge_survivors
            },
            "raw_double_edge_trial_count": self.raw_double_edge_trial_count,
            "monotonicity_residual_double_trial_count": (
                self.monotonicity_residual_double_trial_count
            ),
            "double_edge_passing_count": self.double_edge_passing_count,
            "exceptional_class_histogram": {
                name: count for name, count in self.exceptional_class_histogram
            },
            "exact_index_two_witnesses": [
                witness.to_record() for witness in self.exact_index_two_witnesses
            ],
            "monotonicity_argument": (
                "testing label 3 is extremal for each non-right cross edge; "
                "labels 4 and 5 are also enumerated on all 55 local orbit "
                "survivors, and any multi-edge state would survive after "
                "lowering two of its edges to label 3"
            ),
            "interlacing_argument": (
                "zero cross edges give the orthogonal sum of two Lanner "
                "blocks; one exceptional edge gives an exact index-two "
                "ten-facet principal block; two edges are locally impossible"
            ),
            "unrestricted_compact_coxeter_no_go": (
                self.unrestricted_compact_coxeter_no_go
            ),
        }


@lru_cache(maxsize=1)
def audit_record43_second_nerve_labels() -> Record43SecondNerveLabelAudit:
    nonfaces = tuple(TARGET_MISSING_EDGES) + SECOND_TARGET_MISSING_FIVES
    vertices = _vertices_from_nonfaces(nonfaces)
    left = SECOND_TARGET_MISSING_FIVES[0]
    right = SECOND_TARGET_MISSING_FIVES[1]
    left_states, left_orbits = _distinguished_support_states(left)
    right_states, right_orbits = _distinguished_support_states(right)
    gate = audit_disjoint_lanner_cross_gate(
        vertices, left, right, left_states, right_states
    )
    stronger = _stronger_single_edge_survivors(
        gate, vertices, left_states, right_states
    )
    raw_double, residual_double, double_passing = _double_edge_gate(
        gate, vertices, left_states, right_states
    )

    witnesses = _a5_witnesses(left, right)
    # Exact witness reference graphs use the five canonical Lanner types;
    # this is classification only, not a quotient of the placement census.
    canonical_left = _canonical_lanner_states(left)
    canonical_right = _canonical_lanner_states(right)
    reference_graphs = (
        _labelled_block_graph(
            left, right, canonical_left[1], canonical_right[1], pair(left[0], right[0])
        ),
        _labelled_block_graph(
            left, right, canonical_left[4], canonical_right[1], pair(left[4], right[0])
        ),
    )
    class_histogram: Counter[str] = Counter()
    for left_index, right_index, edge in gate.passing_trials:
        graph = _labelled_block_graph(
            left,
            right,
            left_states[left_index],
            right_states[right_index],
            edge,
        )
        matches = tuple(
            index
            for index, reference in enumerate(reference_graphs)
            if _isomorphic_blocks(graph, reference)
        )
        if len(matches) != 1:
            raise AssertionError(
                "a second-nerve exceptional A5 block has no unique class"
            )
        class_histogram[witnesses[matches[0]].name] += 1

    audit = Record43SecondNerveLabelAudit(
        vertex_count=len(vertices),
        left_support=left,
        right_support=right,
        verified_nerve_automorphism_count=(
            576 if _nerve_automorphism_subgroup_verified() else 0
        ),
        full_left_placement_count=len(_all_lanner_states(left)),
        full_right_placement_count=len(_all_lanner_states(right)),
        left_orbit_size_histogram=left_orbits,
        right_orbit_size_histogram=right_orbits,
        cross_gate=gate,
        stronger_single_edge_survivors=stronger,
        raw_double_edge_trial_count=raw_double,
        monotonicity_residual_double_trial_count=residual_double,
        double_edge_passing_count=double_passing,
        exceptional_class_histogram=tuple(sorted(class_histogram.items())),
        exact_index_two_witnesses=witnesses,
    )
    if not audit.unrestricted_compact_coxeter_no_go:
        raise AssertionError("record-43 second-nerve compact label gate changed")
    return audit


__all__ = [
    "Record43SecondNerveLabelAudit",
    "audit_record43_second_nerve_labels",
]
