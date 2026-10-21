"""Unrestricted compact-Coxeter obstruction for the selected H9/14 record 17.

The exact ordinary-polytopal record-17 nerve has three order-five minimal
nonfaces with a common four-node core::

    {0,3,4,5,6}, {1,3,4,5,6}, {2,3,4,5,6}.

For a compact Coxeter realization each support must carry a classical Lanner
diagram.  Joining all 420 labelled order-five placements on the shared core
gives 4,620 compatible states.  Exact finite-Coxeter classification against
the complete list of 110 actual rank-nine vertices leaves 80 states after the
rank-five restrictions and none after four rank-six restrictions.

This is unrestricted in all remaining finite labels.  In normalized Coxeter
Gram convention ``G = I-A``, changing a right angle to any finite non-right
angle increases the symmetric nonnegative matrix ``A`` entrywise.  Perron--
Frobenius monotonicity therefore cannot repair a principal diagram which is
already non-elliptic in the all-right external shadow.

The result is scoped to this selected nerve and compact realizations.  It does
not address a different completion of the graph shadow or a finite-volume
noncompact chamber with ideal faces.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations

from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
    lanner_placement_codes,
    pair,
)
from coxeter_polytopes.combinatorics.h9_14_record17_polytopality import (
    TARGET_MISSING_EDGES,
    TARGET_MISSING_FIVES,
)
from coxeter_polytopes.combinatorics.h9_14_record35_coxeter_labels import (
    CommonFourCoreLannerShadowAudit,
    audit_common_four_core_lanner_shadow,
)


FACET_COUNT = 14
DIMENSION = 9
COMMON_CORE = (3, 4, 5, 6)
OUTER_FACETS = (0, 1, 2)
COMMON_CORE_SUPPORTS = tuple(
    tuple(sorted((*COMMON_CORE, outer))) for outer in OUTER_FACETS
)
MINIMAL_NONFACES = tuple(TARGET_MISSING_EDGES) + tuple(TARGET_MISSING_FIVES)


def record17_vertex_facets() -> tuple[tuple[int, ...], ...]:
    """Return all 110 rank-nine vertices of the selected record-17 nerve."""

    forbidden = tuple(frozenset(face) for face in MINIMAL_NONFACES)
    return tuple(
        facets
        for facets in combinations(range(FACET_COUNT), DIMENSION)
        if not any(face <= set(facets) for face in forbidden)
    )


@dataclass(frozen=True, slots=True)
class Record17CoxeterLabelAudit:
    facet_count: int
    dimension: int
    vertex_count: int
    meeting_pair_count: int
    dotted_pairs_in_cluster_union: tuple[tuple[int, int], ...]
    lanner_placements_per_support: int
    gate: CommonFourCoreLannerShadowAudit

    @property
    def unrestricted_compact_coxeter_no_go(self) -> bool:
        return (
            not self.dotted_pairs_in_cluster_union
            and self.gate.unrestricted_external_label_no_go
        )

    @property
    def verified(self) -> bool:
        stages = tuple(
            (
                stage.order,
                len(stage.restrictions),
                stage.passing_state_count,
            )
            for stage in self.gate.stages
        )
        return (
            self.facet_count == 14
            and self.dimension == 9
            and self.vertex_count == 110
            and self.meeting_pair_count == 86
            and not self.dotted_pairs_in_cluster_union
            and self.lanner_placements_per_support == 420
            and self.gate.outer_facets == OUTER_FACETS
            and self.gate.common_core == COMMON_CORE
            and self.gate.supports == COMMON_CORE_SUPPORTS
            and self.gate.join_counts == (420, 1_332, 4_620)
            and self.gate.core_extension_count_histogram
            == ((1, 36), (2, 30), (3, 40), (4, 51))
            and tuple(
                sorted(Counter(map(len, self.gate.actual_vertex_restrictions)).items())
            )
            == ((4, 9), (5, 16), (6, 4))
            and stages == ((4, 9, 4_620), (5, 16, 80), (6, 4, 0))
            and self.gate.last_nonempty_stage_order == 5
            and self.gate.sparse_shadow_passing_state_count == 0
            and self.unrestricted_compact_coxeter_no_go
        )

    def to_record(self) -> dict[str, object]:
        return {
            "facet_count": self.facet_count,
            "dimension": self.dimension,
            "rank_nine_vertex_count": self.vertex_count,
            "meeting_pair_count": self.meeting_pair_count,
            "minimal_nonfaces": {
                "dotted_pairs": [list(face) for face in TARGET_MISSING_EDGES],
                "order_five_lanner_supports": [
                    list(face) for face in TARGET_MISSING_FIVES
                ],
            },
            "decisive_common_four_core_gate": self.gate.to_record(),
            "dotted_pairs_in_decisive_cluster_union": [
                list(edge) for edge in self.dotted_pairs_in_cluster_union
            ],
            "lanner_placements_per_support": self.lanner_placements_per_support,
            "unrestricted_compact_coxeter_no_go": {
                "proved": self.unrestricted_compact_coxeter_no_go,
                "local_obstruction": (
                    "All 4,620 compatible placements on the three common-core "
                    "Lanner supports fail ellipticity on at least one actual "
                    "rank-nine vertex restriction; 80 survive through rank "
                    "five and zero survive the four rank-six restrictions."
                ),
                "monotonicity": (
                    "For G=I-A, changing any remaining m=2 meeting pair to "
                    "m>=3 increases the symmetric nonnegative matrix A "
                    "entrywise, so Perron--Frobenius monotonicity cannot turn "
                    "a non-elliptic principal diagram into an elliptic one."
                ),
                "gram_search_needed": False,
            },
            "scope": {
                "selected_nerve_only": True,
                "compact_only": True,
                "other_completions_of_graph_shadow_excluded": False,
                "finite_volume_noncompact_excluded": False,
            },
            "verified": self.verified,
        }


@lru_cache(maxsize=1)
def audit_record17_coxeter_labels() -> Record17CoxeterLabelAudit:
    """Exhaust the common-core Lanner joins and certify the compact no-go."""

    vertices = record17_vertex_facets()
    cluster_union = set((*OUTER_FACETS, *COMMON_CORE))
    dotted = tuple(
        sorted(
            pair(left, right)
            for left, right in TARGET_MISSING_EDGES
            if left in cluster_union and right in cluster_union
        )
    )
    gate = audit_common_four_core_lanner_shadow(
        OUTER_FACETS, COMMON_CORE, vertices
    )
    audit = Record17CoxeterLabelAudit(
        facet_count=FACET_COUNT,
        dimension=DIMENSION,
        vertex_count=len(vertices),
        meeting_pair_count=len(tuple(combinations(range(FACET_COUNT), 2)))
        - len(TARGET_MISSING_EDGES),
        dotted_pairs_in_cluster_union=dotted,
        lanner_placements_per_support=len(lanner_placement_codes(5)),
        gate=gate,
    )
    if not audit.verified:
        raise AssertionError("record-17 compact-Coxeter label gate changed")
    return audit


__all__ = [
    "COMMON_CORE",
    "COMMON_CORE_SUPPORTS",
    "OUTER_FACETS",
    "Record17CoxeterLabelAudit",
    "audit_record17_coxeter_labels",
    "record17_vertex_facets",
]
