"""Exact compact-Coxeter label obstruction for H9/14 record 35.

Record 35 has five minimal nonfaces of order five.  In a compact hyperbolic
Coxeter polytope each must support one of the classical order-five Lanner
diagrams, while every nine-facet vertex diagram must be elliptic.

Three of the missing fives have a common four-node core::

    {0,5,7,9,10}, {1,5,7,9,10}, {4,5,7,9,10}.

This module joins all 420 labelled Lanner placements on those supports.  The
12 relevant rank-five vertex restrictions leave 96 assignments.  Modulo the
natural ``S3 x S4`` action they form two orbits, and their representatives
contain respectively ``affine D4 + A1`` and ``affine D5`` rank-six principal
diagrams inside actual vertices.  Hence no all-right external completion is
locally elliptic.

The conclusion is unrestricted in the labels outside the forced Lanner
supports.  Write a normalized Coxeter Gram matrix as ``I - A``, where
``A_ij = cos(pi/m_ij)`` is nonnegative.  Replacing an external right angle by
any finite non-right Coxeter angle increases ``A`` entrywise.  Perron--
Frobenius monotonicity therefore cannot turn a sparse-shadow matrix with
spectral radius at least one into an elliptic matrix.  Thus record 35 is not
the nerve of a compact hyperbolic Coxeter polytope.  No global Gram or volume
claim is made for any other nerve.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations, permutations
from typing import Iterable, Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
    Pair,
    is_finite_coxeter_labeling,
    lanner_placement_codes,
    pair,
)
from coxeter_polytopes.combinatorics.h9_14_record35 import (
    FULL_MASK,
    TARGET_MISSING_EDGES,
    TARGET_MISSING_FIVES,
    record35_minimal_transversals,
    vertices,
)


FACET_COUNT = 14
ALL_PAIRS = tuple(combinations(range(FACET_COUNT), 2))

C_SUPPORTS = (TARGET_MISSING_FIVES[0], TARGET_MISSING_FIVES[4])
C_LEFT = (2, 4)
C_CORE = (3, 6, 8)
C_RIGHT = (11, 12)
C_UNION = tuple(sorted((*C_LEFT, *C_CORE, *C_RIGHT)))

D_SUPPORTS = TARGET_MISSING_FIVES[1:4]
D_OUTER = (0, 1, 4)
D_CORE = (5, 7, 9, 10)
D_UNION = tuple(sorted((*D_OUTER, *D_CORE)))
D_EXTERNAL_PAIRS = tuple(combinations(D_OUTER, 2))

LabelState = tuple[tuple[Pair, int], ...]


def record35_vertex_facets() -> tuple[tuple[int, ...], ...]:
    """Return the complete deterministic list of 102 rank-nine vertices."""

    return tuple(
        sorted(
            vertices(FULL_MASK ^ complement)
            for complement in record35_minimal_transversals()
        )
    )


@lru_cache(maxsize=None)
def _placement_options(support: tuple[int, ...]) -> tuple[LabelState, ...]:
    support = tuple(sorted(support))
    if len(support) != 5 or len(set(support)) != 5:
        raise ValueError("a Lanner support must have five distinct facets")
    edges = tuple(combinations(support, 2))
    return tuple(
        tuple(zip(edges, code, strict=True))
        for code in lanner_placement_codes(5)
    )


def _merge(left: LabelState, right: LabelState) -> LabelState | None:
    labels = dict(left)
    for edge, label in right:
        previous = labels.get(edge)
        if previous is not None and previous != label:
            return None
        labels[edge] = label
    return tuple(sorted(labels.items()))


def _join_tables(
    supports: Sequence[tuple[int, ...]],
) -> tuple[tuple[LabelState, ...], tuple[int, ...]]:
    states: tuple[LabelState, ...] = ((),)
    layer_counts: list[int] = []
    for support in supports:
        states = tuple(
            sorted(
                {
                    merged
                    for state in states
                    for option in _placement_options(tuple(support))
                    if (merged := _merge(state, option)) is not None
                }
            )
        )
        layer_counts.append(len(states))
    return states, tuple(layer_counts)


def _full_labels(union: Sequence[int], state: LabelState) -> dict[Pair, int]:
    labels = {
        pair(left, right): 2 for left, right in combinations(sorted(union), 2)
    }
    labels.update(state)
    return labels


def _finite(facets: Iterable[int], labels: Mapping[Pair, int]) -> bool:
    return is_finite_coxeter_labeling(tuple(sorted(facets)), labels)


def _cluster_restrictions(
    cluster_union: Sequence[int],
) -> tuple[tuple[int, ...], ...]:
    support = set(cluster_union)
    return tuple(
        sorted(
            {
                tuple(facet for facet in vertex if facet in support)
                for vertex in record35_vertex_facets()
            },
            key=lambda item: (len(item), item),
        )
    )


def _partition_canonical_code(
    labels: Mapping[Pair, int], partitions: Sequence[tuple[int, ...]]
) -> tuple[int, ...]:
    union = tuple(sorted(facet for block in partitions for facet in block))
    edges = tuple(combinations(union, 2))
    best: tuple[int, ...] | None = None

    def visit(index: int, image: dict[int, int]) -> None:
        nonlocal best
        if index == len(partitions):
            code = tuple(
                labels[pair(image[left], image[right])] for left, right in edges
            )
            if best is None or code < best:
                best = code
            return
        block = partitions[index]
        for permuted in permutations(block):
            visit(index + 1, image | dict(zip(block, permuted, strict=True)))

    visit(0, {})
    if best is None:
        raise AssertionError("partition canonicalization produced no image")
    return best


def _labels_from_code(
    union: Sequence[int], code: Sequence[int]
) -> dict[Pair, int]:
    return dict(zip(combinations(sorted(union), 2), code, strict=True))


def _nonright_labels(labels: Mapping[Pair, int]) -> tuple[tuple[int, int, int], ...]:
    return tuple(
        (left, right, label)
        for (left, right), label in sorted(labels.items())
        if label > 2
    )


def _contains_in_vertex(facets: Sequence[int]) -> tuple[int, ...]:
    target = set(facets)
    return next(
        vertex for vertex in record35_vertex_facets() if target <= set(vertex)
    )


def _simply_laced_cartan(
    facets: Sequence[int], labels: Mapping[Pair, int]
) -> tuple[tuple[int, ...], ...]:
    ordered = tuple(facets)
    if any(
        labels[pair(left, right)] not in (2, 3)
        for left, right in combinations(ordered, 2)
    ):
        raise ValueError("the affine witness is not simply laced")
    return tuple(
        tuple(
            2
            if left == right
            else -1
            if labels[pair(left, right)] == 3
            else 0
            for right in ordered
        )
        for left in ordered
    )


def _matrix_vector(
    matrix: Sequence[Sequence[int]], vector: Sequence[int]
) -> tuple[int, ...]:
    return tuple(
        sum(entry * coefficient for entry, coefficient in zip(row, vector, strict=True))
        for row in matrix
    )


@dataclass(frozen=True, slots=True)
class CommonFourCoreRestrictionStage:
    """One cumulative vertex-restriction stage in the reusable core gate."""

    order: int
    restrictions: tuple[tuple[int, ...], ...]
    passing_state_count: int

    def to_record(self) -> dict[str, object]:
        return {
            "order": self.order,
            "restrictions": [list(facets) for facets in self.restrictions],
            "restriction_count": len(self.restrictions),
            "cumulative_passing_state_count": self.passing_state_count,
        }


@dataclass(frozen=True, slots=True)
class CommonFourCoreLannerShadowAudit:
    """Reusable local gate for three order-five Lanners sharing four nodes.

    The conclusion is tied to ``actual_vertex_restrictions``.  It may be used
    for another nerve only after supplying that nerve's complete vertex list.
    If ``sparse_shadow_passing_state_count`` is zero, Perron--Frobenius
    monotonicity promotes the shadow failure to all external non-right labels.
    """

    outer_facets: tuple[int, ...]
    common_core: tuple[int, ...]
    supports: tuple[tuple[int, ...], ...]
    join_counts: tuple[int, ...]
    core_extension_count_histogram: tuple[tuple[int, int], ...]
    actual_vertex_restrictions: tuple[tuple[int, ...], ...]
    stages: tuple[CommonFourCoreRestrictionStage, ...]
    natural_partition_symmetry_preserves_restrictions: bool
    last_nonempty_stage_order: int
    last_nonempty_partition_orbit_size_histogram: tuple[tuple[int, int], ...]
    last_nonempty_partition_orbit_representative_nonright_labels: tuple[
        tuple[tuple[int, int, int], ...], ...
    ]
    sparse_shadow_passing_state_count: int

    @property
    def unrestricted_external_label_no_go(self) -> bool:
        return self.sparse_shadow_passing_state_count == 0

    def to_record(self) -> dict[str, object]:
        return {
            "outer_facets": list(self.outer_facets),
            "common_core": list(self.common_core),
            "supports": [list(face) for face in self.supports],
            "join_counts": list(self.join_counts),
            "core_extension_count_histogram": {
                str(extension_count): core_code_count
                for extension_count, core_code_count
                in self.core_extension_count_histogram
            },
            "actual_vertex_restrictions": [
                list(facets) for facets in self.actual_vertex_restrictions
            ],
            "cumulative_restriction_stages": [
                stage.to_record() for stage in self.stages
            ],
            "natural_S3_x_S4_symmetry_preserves_restrictions": (
                self.natural_partition_symmetry_preserves_restrictions
            ),
            "last_nonempty_stage_order": self.last_nonempty_stage_order,
            "last_nonempty_partition_orbit_size_histogram": {
                str(size): count
                for size, count
                in self.last_nonempty_partition_orbit_size_histogram
            },
            "last_nonempty_partition_orbit_representative_nonright_labels": [
                [list(edge) for edge in representative]
                for representative
                in self.last_nonempty_partition_orbit_representative_nonright_labels
            ],
            "sparse_shadow_passing_state_count": (
                self.sparse_shadow_passing_state_count
            ),
            "unrestricted_external_label_no_go": (
                self.unrestricted_external_label_no_go
            ),
            "reuse_scope": (
                "The no-go is valid only for a nerve whose supplied actual "
                "vertex list contains these restrictions; no abstract support "
                "pattern alone is asserted to be obstructed."
            ),
        }


def _partition_symmetry_preserves_restrictions(
    restrictions: Sequence[tuple[int, ...]],
    partitions: Sequence[tuple[int, ...]],
) -> bool:
    restriction_set = {frozenset(facets) for facets in restrictions}
    generators: list[dict[int, int]] = []
    for block in partitions:
        for index in range(len(block) - 1):
            image = {facet: facet for part in partitions for facet in part}
            image[block[index]], image[block[index + 1]] = (
                block[index + 1],
                block[index],
            )
            generators.append(image)
    return all(
        {
            frozenset(image[facet] for facet in restriction)
            for restriction in restrictions
        }
        == restriction_set
        for image in generators
    )


def _common_four_core_stage_survivors(
    outer_facets: tuple[int, ...],
    common_core: tuple[int, ...],
    vertex_facets: Sequence[Sequence[int]],
) -> tuple[
    tuple[LabelState, ...],
    tuple[int, ...],
    tuple[tuple[int, ...], ...],
    tuple[tuple[int, tuple[dict[Pair, int], ...]], ...],
]:
    supports = tuple(
        tuple(sorted((*common_core, outer))) for outer in outer_facets
    )
    states, join_counts = _join_tables(supports)
    union = tuple(sorted((*outer_facets, *common_core)))
    union_set = set(union)
    restrictions = tuple(
        sorted(
            {
                tuple(facet for facet in vertex if facet in union_set)
                for vertex in vertex_facets
            },
            key=lambda item: (len(item), item),
        )
    )
    labels = tuple(_full_labels(union, state) for state in states)
    staged: list[tuple[int, tuple[dict[Pair, int], ...]]] = []
    survivors = labels
    for order in sorted(set(map(len, restrictions))):
        tests = tuple(facets for facets in restrictions if len(facets) == order)
        survivors = tuple(
            labeling
            for labeling in survivors
            if all(_finite(facets, labeling) for facets in tests)
        )
        staged.append((order, survivors))
    return states, join_counts, restrictions, tuple(staged)


def audit_common_four_core_lanner_shadow(
    outer_facets: Iterable[int],
    common_core: Iterable[int],
    vertex_facets: Sequence[Sequence[int]],
) -> CommonFourCoreLannerShadowAudit:
    """Audit three Lanner supports ``K union {a,b,c}`` against real vertices.

    This is the reusable form of the record-35 obstruction.  The function
    does not infer vertex restrictions from the support pattern: callers must
    provide the complete actual vertex list of the nerve under study.
    """

    outer = tuple(sorted(set(outer_facets)))
    core = tuple(sorted(set(common_core)))
    if len(outer) != 3 or len(core) != 4 or set(outer) & set(core):
        raise ValueError("the gate needs three outer facets and a disjoint four-core")
    vertices_supplied = tuple(tuple(sorted(set(vertex))) for vertex in vertex_facets)
    if not vertices_supplied:
        raise ValueError("the gate needs the complete nonempty actual vertex list")
    union = set((*outer, *core))
    if any(not set(vertex) <= set(range(FACET_COUNT)) for vertex in vertices_supplied):
        raise ValueError("a supplied vertex lies outside the H9/14 facet set")
    if not union <= set().union(*(set(vertex) for vertex in vertices_supplied)):
        raise ValueError("some common-core cluster facet occurs in no actual vertex")

    states, join_counts, restrictions, staged = _common_four_core_stage_survivors(
        outer, core, vertices_supplied
    )
    supports = tuple(tuple(sorted((*core, facet))) for facet in outer)
    extensions: defaultdict[tuple[int, ...], int] = defaultdict(int)
    for option in _placement_options(supports[0]):
        labels = dict(option)
        core_code = tuple(
            labels[pair(left, right)] for left, right in combinations(core, 2)
        )
        extensions[core_code] += 1
    extension_histogram = tuple(sorted(Counter(extensions.values()).items()))
    if (
        sum(extensions.values()),
        sum(count**2 for count in extensions.values()),
        sum(count**3 for count in extensions.values()),
    ) != join_counts:
        raise AssertionError("common-core extension moments do not match joins")

    stages = tuple(
        CommonFourCoreRestrictionStage(
            order=order,
            restrictions=tuple(
                facets for facets in restrictions if len(facets) == order
            ),
            passing_state_count=len(survivors),
        )
        for order, survivors in staged
    )
    last_nonempty = next(
        ((order, survivors) for order, survivors in reversed(staged) if survivors),
        (0, ()),
    )
    symmetry = _partition_symmetry_preserves_restrictions(
        restrictions, (outer, core)
    )
    orbit_histogram: tuple[tuple[int, int], ...] = ()
    orbit_representatives: tuple[tuple[tuple[int, int, int], ...], ...] = ()
    if symmetry and last_nonempty[1]:
        orbit_codes = Counter(
            _partition_canonical_code(labels, (outer, core))
            for labels in last_nonempty[1]
        )
        orbit_histogram = tuple(sorted(Counter(orbit_codes.values()).items()))
        union_order = tuple(sorted(union))
        orbit_representatives = tuple(
            _nonright_labels(_labels_from_code(union_order, code))
            for code in sorted(orbit_codes)
        )

    return CommonFourCoreLannerShadowAudit(
        outer_facets=outer,
        common_core=core,
        supports=supports,
        join_counts=join_counts,
        core_extension_count_histogram=extension_histogram,
        actual_vertex_restrictions=restrictions,
        stages=stages,
        natural_partition_symmetry_preserves_restrictions=symmetry,
        last_nonempty_stage_order=last_nonempty[0],
        last_nonempty_partition_orbit_size_histogram=orbit_histogram,
        last_nonempty_partition_orbit_representative_nonright_labels=(
            orbit_representatives
        ),
        sparse_shadow_passing_state_count=(
            len(staged[-1][1]) if staged else len(states)
        ),
    )


@dataclass(frozen=True, slots=True)
class Record35DOrbit:
    orbit_size: int
    nonright_labels: tuple[tuple[int, int, int], ...]
    affine_type: str
    witness_facets: tuple[int, ...]
    containing_vertex: tuple[int, ...]
    cartan_null_vector: tuple[int, ...]
    failed_rank_six_restriction_count: int

    def to_record(self) -> dict[str, object]:
        return {
            "orbit_size": self.orbit_size,
            "representative_nonright_labels": [
                list(edge) for edge in self.nonright_labels
            ],
            "affine_witness": {
                "type": self.affine_type,
                "facets": list(self.witness_facets),
                "containing_rank_nine_vertex": list(self.containing_vertex),
                "simply_laced_cartan_null_vector": list(self.cartan_null_vector),
            },
            "failed_rank_six_restriction_count": (
                self.failed_rank_six_restriction_count
            ),
        }


@dataclass(frozen=True, slots=True)
class Record35CoxeterLabelAudit:
    facet_count: int
    vertex_count: int
    meeting_pair_count: int
    lanner_placements_per_support: int
    c_join_counts: tuple[int, ...]
    c_restriction_size_histogram: tuple[tuple[int, int], ...]
    c_passing_state_count: int
    c_orbit_size_histogram: tuple[tuple[int, int], ...]
    c_representative_nonright_labels: tuple[tuple[int, int, int], ...]
    d_common_four_core_gate: CommonFourCoreLannerShadowAudit
    d_core_extension_count_histogram: tuple[tuple[int, int], ...]
    d_join_counts: tuple[int, ...]
    d_rank_four_restriction_count: int
    d_rank_four_passing_state_count: int
    d_rank_five_restriction_count: int
    d_rank_five_passing_state_count: int
    d_rank_six_restriction_count: int
    d_rank_six_failure_multiplicity_histogram: tuple[tuple[int, int], ...]
    d_survivor_orbits: tuple[Record35DOrbit, ...]
    compatible_five_lanner_state_count: int
    sparse_passing_state_count: int
    one_external_nonright_state_count: int
    one_external_nonright_passing_state_count: int
    unrestricted_external_label_no_go: bool

    @property
    def verified(self) -> bool:
        return (
            self.facet_count == 14
            and self.vertex_count == 102
            and self.meeting_pair_count == 86
            and self.lanner_placements_per_support == 420
            and self.c_join_counts == (420, 11_634)
            and self.c_restriction_size_histogram == ((4, 10), (5, 16), (6, 3))
            and self.c_passing_state_count == 24
            and self.c_orbit_size_histogram == ((24, 1),)
            and self.d_common_four_core_gate.unrestricted_external_label_no_go
            and self.d_common_four_core_gate.last_nonempty_stage_order == 5
            and self.d_common_four_core_gate
            .last_nonempty_partition_orbit_size_histogram
            == ((24, 1), (72, 1))
            and self.d_core_extension_count_histogram
            == ((1, 36), (2, 30), (3, 40), (4, 51))
            and self.d_join_counts == (420, 1_332, 4_620)
            and self.d_rank_four_restriction_count == 13
            and self.d_rank_four_passing_state_count == 4_620
            and self.d_rank_five_restriction_count == 12
            and self.d_rank_five_passing_state_count == 96
            and self.d_rank_six_restriction_count == 4
            and self.d_rank_six_failure_multiplicity_histogram
            == ((1, 72), (2, 24))
            and tuple(orbit.orbit_size for orbit in self.d_survivor_orbits)
            == (24, 72)
            and tuple(orbit.affine_type for orbit in self.d_survivor_orbits)
            == ("affine_D4_plus_A1", "affine_D5")
            and self.compatible_five_lanner_state_count == 53_749_080
            and self.sparse_passing_state_count == 0
            and self.one_external_nonright_state_count == 41_580
            and self.one_external_nonright_passing_state_count == 0
            and self.unrestricted_external_label_no_go
        )

    def to_record(self) -> dict[str, object]:
        return {
            "facet_count": self.facet_count,
            "rank_nine_vertex_count": self.vertex_count,
            "meeting_pair_count": self.meeting_pair_count,
            "minimal_nonfaces": {
                "dotted_pairs": [list(edge) for edge in TARGET_MISSING_EDGES],
                "order_five_lanner_supports": [
                    list(face) for face in TARGET_MISSING_FIVES
                ],
            },
            "lanner_placements_per_support": self.lanner_placements_per_support,
            "two_support_cluster": {
                "supports": [list(face) for face in C_SUPPORTS],
                "join_counts": list(self.c_join_counts),
                "induced_vertex_restriction_size_histogram": {
                    str(order): count
                    for order, count in self.c_restriction_size_histogram
                },
                "passing_state_count": self.c_passing_state_count,
                "symmetry_group": "S2 x S3 x S2",
                "orbit_size_histogram": {
                    str(size): count
                    for size, count in self.c_orbit_size_histogram
                },
                "representative_nonright_labels": [
                    list(edge) for edge in self.c_representative_nonright_labels
                ],
            },
            "three_support_common_core_cluster": {
                "supports": [list(face) for face in D_SUPPORTS],
                "outer_facets": list(D_OUTER),
                "common_core": list(D_CORE),
                "common_core_extension_count_histogram": {
                    str(extension_count): core_code_count
                    for extension_count, core_code_count
                    in self.d_core_extension_count_histogram
                },
                "join_counts": list(self.d_join_counts),
                "rank_four_restrictions": {
                    "count": self.d_rank_four_restriction_count,
                    "passing_states": self.d_rank_four_passing_state_count,
                },
                "rank_five_restrictions": {
                    "count": self.d_rank_five_restriction_count,
                    "passing_states": self.d_rank_five_passing_state_count,
                },
                "rank_six_restrictions": {
                    "count": self.d_rank_six_restriction_count,
                    "failure_multiplicity_histogram": {
                        str(failures): states
                        for failures, states
                        in self.d_rank_six_failure_multiplicity_histogram
                    },
                },
                "rank_five_survivor_symmetry_group": "S3 x S4",
                "rank_five_survivor_orbits": [
                    orbit.to_record() for orbit in self.d_survivor_orbits
                ],
                "reusable_common_four_core_gate": (
                    self.d_common_four_core_gate.to_record()
                ),
            },
            "sparse_ansatz": {
                "definition": (
                    "choose an arbitrary classical order-five Lanner placement "
                    "on each support, require agreement on overlaps, and label "
                    "every other meeting pair by 2"
                ),
                "compatible_five_lanner_state_count": (
                    self.compatible_five_lanner_state_count
                ),
                "passing_state_count": self.sparse_passing_state_count,
            },
            "one_extra_nonright_pair_shell": {
                "eligible_pairs_inside_the_obstructing_cluster": [
                    list(edge) for edge in D_EXTERNAL_PAIRS
                ],
                "labels": [3, 4, 5],
                "state_count": self.one_external_nonright_state_count,
                "passing_state_count": (
                    self.one_external_nonright_passing_state_count
                ),
            },
            "unrestricted_compact_coxeter_no_go": {
                "proved": self.unrestricted_external_label_no_go,
                "monotonicity_argument": (
                    "For normalized Gram I-A, every added non-right label "
                    "increases the nonnegative cosine adjacency A entrywise. "
                    "Perron-Frobenius monotonicity cannot reduce spectral "
                    "radius from at least 1 to below 1, so no external labels "
                    "can repair a non-elliptic sparse-shadow principal diagram."
                ),
                "scope": (
                    "rules out record 35 as the nerve of a compact hyperbolic "
                    "Coxeter polytope; it does not assert a global Gram or "
                    "finite-volume result for any other nerve"
                ),
            },
            "verified": self.verified,
        }


def _build_d_orbits(
    survivors: Sequence[Mapping[Pair, int]],
    rank_six: Sequence[tuple[int, ...]],
) -> tuple[Record35DOrbit, ...]:
    orbit_codes = Counter(
        _partition_canonical_code(labels, (D_OUTER, D_CORE))
        for labels in survivors
    )
    if tuple(sorted(orbit_codes.values())) != (24, 72):
        raise AssertionError("the D-cluster survivor orbit census changed")

    result: list[Record35DOrbit] = []
    for code, orbit_size in sorted(orbit_codes.items(), key=lambda item: item[1]):
        labels = _labels_from_code(D_UNION, code)
        bad = tuple(facets for facets in rank_six if not _finite(facets, labels))
        if orbit_size == 24:
            affine_type = "affine_D4_plus_A1"
            witness = (0, 1, 4, 5, 7, 10)
            null_vector = (1, 1, 1, 1, 0, 2)
        elif orbit_size == 72:
            affine_type = "affine_D5"
            witness = (0, 1, 4, 7, 9, 10)
            null_vector = (1, 1, 1, 1, 2, 2)
        else:
            raise AssertionError("unexpected D-cluster orbit size")
        if witness not in bad:
            raise AssertionError("the canonical affine witness changed")
        cartan = _simply_laced_cartan(witness, labels)
        if _matrix_vector(cartan, null_vector) != (0,) * len(witness):
            raise AssertionError("the declared affine null vector is not null")
        result.append(
            Record35DOrbit(
                orbit_size=orbit_size,
                nonright_labels=_nonright_labels(labels),
                affine_type=affine_type,
                witness_facets=witness,
                containing_vertex=_contains_in_vertex(witness),
                cartan_null_vector=null_vector,
                failed_rank_six_restriction_count=len(bad),
            )
        )
    return tuple(result)


@lru_cache(maxsize=1)
def audit_record35_coxeter_labels() -> Record35CoxeterLabelAudit:
    """Exhaust the forced Lanner joins and certify the compact no-go."""

    vertex_facets = record35_vertex_facets()
    if len(vertex_facets) != 102 or any(len(vertex) != 9 for vertex in vertex_facets):
        raise AssertionError("record 35 no longer has 102 rank-nine vertices")
    dotted = {pair(*edge) for edge in TARGET_MISSING_EDGES}
    meeting_pairs = tuple(edge for edge in ALL_PAIRS if edge not in dotted)
    common_four_core_gate = audit_common_four_core_lanner_shadow(
        D_OUTER, D_CORE, vertex_facets
    )

    c_states, c_join_counts = _join_tables(C_SUPPORTS)
    c_restrictions = _cluster_restrictions(C_UNION)
    c_passing: list[dict[Pair, int]] = []
    for state in c_states:
        labels = _full_labels(C_UNION, state)
        if all(_finite(facets, labels) for facets in c_restrictions):
            c_passing.append(labels)
    c_orbits = Counter(
        _partition_canonical_code(labels, (C_LEFT, C_CORE, C_RIGHT))
        for labels in c_passing
    )
    if len(c_orbits) != 1:
        raise AssertionError("the C-cluster survivor orbit census changed")
    c_representative = _labels_from_code(C_UNION, next(iter(c_orbits)))

    d_states, d_join_counts = _join_tables(D_SUPPORTS)
    extensions: defaultdict[tuple[int, ...], int] = defaultdict(int)
    for option in _placement_options(D_SUPPORTS[0]):
        labels = dict(option)
        core_code = tuple(
            labels[pair(left, right)] for left, right in combinations(D_CORE, 2)
        )
        extensions[core_code] += 1
    extension_histogram = tuple(sorted(Counter(extensions.values()).items()))
    if (
        sum(extensions.values()),
        sum(count**2 for count in extensions.values()),
        sum(count**3 for count in extensions.values()),
    ) != d_join_counts:
        raise AssertionError("common-core extension moments do not match joins")

    d_restrictions = _cluster_restrictions(D_UNION)
    d_by_order = {
        order: tuple(facets for facets in d_restrictions if len(facets) == order)
        for order in (4, 5, 6)
    }
    d_labels = tuple(_full_labels(D_UNION, state) for state in d_states)
    rank_four_passing = tuple(
        labels
        for labels in d_labels
        if all(_finite(facets, labels) for facets in d_by_order[4])
    )
    rank_five_passing = tuple(
        labels
        for labels in rank_four_passing
        if all(_finite(facets, labels) for facets in d_by_order[5])
    )
    rank_six_failure_histogram = tuple(
        sorted(
            Counter(
                sum(
                    not _finite(facets, labels)
                    for facets in d_by_order[6]
                )
                for labels in rank_five_passing
            ).items()
        )
    )
    if any(
        all(_finite(facets, labels) for facets in d_by_order[6])
        for labels in rank_five_passing
    ):
        raise AssertionError("a D-cluster sparse state passed every restriction")
    d_orbits = _build_d_orbits(rank_five_passing, d_by_order[6])

    one_edge_states = 0
    one_edge_passing = 0
    for labels in d_labels:
        for edge in D_EXTERNAL_PAIRS:
            if labels[edge] != 2:
                raise AssertionError("an external D-cluster pair lies in a Lanner support")
            for label in (3, 4, 5):
                one_edge_states += 1
                extension = dict(labels)
                extension[edge] = label
                if all(
                    _finite(facets, extension)
                    for order in (4, 5, 6)
                    for facets in d_by_order[order]
                ):
                    one_edge_passing += 1

    audit = Record35CoxeterLabelAudit(
        facet_count=FACET_COUNT,
        vertex_count=len(vertex_facets),
        meeting_pair_count=len(meeting_pairs),
        lanner_placements_per_support=len(_placement_options(D_SUPPORTS[0])),
        c_join_counts=c_join_counts,
        c_restriction_size_histogram=tuple(
            sorted(Counter(map(len, c_restrictions)).items())
        ),
        c_passing_state_count=len(c_passing),
        c_orbit_size_histogram=tuple(
            sorted(Counter(c_orbits.values()).items())
        ),
        c_representative_nonright_labels=_nonright_labels(c_representative),
        d_common_four_core_gate=common_four_core_gate,
        d_core_extension_count_histogram=extension_histogram,
        d_join_counts=d_join_counts,
        d_rank_four_restriction_count=len(d_by_order[4]),
        d_rank_four_passing_state_count=len(rank_four_passing),
        d_rank_five_restriction_count=len(d_by_order[5]),
        d_rank_five_passing_state_count=len(rank_five_passing),
        d_rank_six_restriction_count=len(d_by_order[6]),
        d_rank_six_failure_multiplicity_histogram=rank_six_failure_histogram,
        d_survivor_orbits=d_orbits,
        compatible_five_lanner_state_count=len(c_states) * len(d_states),
        sparse_passing_state_count=0,
        one_external_nonright_state_count=one_edge_states,
        one_external_nonright_passing_state_count=one_edge_passing,
        unrestricted_external_label_no_go=True,
    )
    if not audit.verified:
        raise AssertionError("the record-35 compact Coxeter label gate changed")
    return audit


__all__ = [
    "C_SUPPORTS",
    "CommonFourCoreLannerShadowAudit",
    "CommonFourCoreRestrictionStage",
    "D_CORE",
    "D_EXTERNAL_PAIRS",
    "D_OUTER",
    "D_SUPPORTS",
    "Record35CoxeterLabelAudit",
    "Record35DOrbit",
    "audit_common_four_core_lanner_shadow",
    "audit_record35_coxeter_labels",
    "record35_vertex_facets",
]
