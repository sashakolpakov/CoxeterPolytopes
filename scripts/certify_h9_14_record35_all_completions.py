#!/usr/bin/env python3
"""Freeze the record-35 whole-shadow compact-Coxeter no-go theorem."""

from __future__ import annotations

import argparse
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record35 import (
    TARGET_GRAPH6,
    TARGET_MISSING_EDGES,
    TARGET_RECORD_INDEX,
)
from coxeter_polytopes.combinatorics.h9_14_record35_all_completions import (
    EXPECTED_INDEPENDENT_ACTIVE_MASKS,
    EXPECTED_SCHEDULES,
    audit_record35_all_completions,
)
from coxeter_polytopes.combinatorics.nerve_support import (
    verify_support_manifest,
)
from scripts.independent_h9_14_record35_all_completions_audit import (
    EXPECTED_SUMMARY_SHA256 as EXPECTED_HOSTILE_SUMMARY_SHA256,
    replay as independent_replay,
)


ROOT = Path(__file__).resolve().parents[1]
SUPPORT = ROOT / "data/enumerations/h9_14_nerve_support.json"
OUTPUT = ROOT / "data/enumerations/h9_14_record35_all_completions.json"

PROOF_SOURCE = (
    ROOT
    / "src/coxeter_polytopes/combinatorics/h9_14_record35_all_completions.py"
)
RECORD35_SOURCE = (
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record35.py"
)
FACET_AUDIT_SOURCE = (
    ROOT / "src/coxeter_polytopes/combinatorics/facet_completion.py"
)
HOMOLOGY_SOURCE = (
    ROOT / "src/coxeter_polytopes/combinatorics/simplicial_homology.py"
)
COXETER_LABEL_SOURCE = (
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_coxeter_labels.py"
)
EXACT_WITNESS_SOURCE = (
    ROOT
    / "src/coxeter_polytopes/combinatorics/"
    "h9_14_selected_nerve_coxeter_no_go.py"
)
RECORD46_SOURCE = (
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record46.py"
)
RECORD92_SOURCE = (
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record92.py"
)
LOCAL_LINK_SOURCE = (
    ROOT / "src/coxeter_polytopes/combinatorics/local_link.py"
)
MISSING_EDGES_SOURCE = (
    ROOT / "src/coxeter_polytopes/combinatorics/missing_edges.py"
)
NERVE_SUPPORT_SOURCE = (
    ROOT / "src/coxeter_polytopes/combinatorics/nerve_support.py"
)
INDEPENDENT_SOURCE = (
    ROOT / "scripts/independent_h9_14_record35_all_completions_audit.py"
)
CERTIFIER_SOURCE = Path(__file__).resolve()
TEST_SOURCE = ROOT / "tests/test_h9_14_record35_all_completions.py"

EXPECTED_SUPPORT_FILE_SHA256 = (
    "8121f869c182b1780f3d812d34823f1579a5c37b085a4a8ecc313ea94937c47d"
)
EXPECTED_CERTIFICATE_SHA256 = "4ebb69e2415ff97c37d62bd8347a2e44b19288cd2d3afa3893b9d0bd763227f8"
EXPECTED_ARCHIVE_FILE_SHA256 = "4b903a7010bfab088b5fd5f0e3fb80e809e0450602fa329d6ba76b58325e8c47"


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    # The in-memory audits contain integer histogram keys, whereas a parsed
    # JSON archive contains their string forms. Normalize through JSON before
    # sorting so write-time and read-time closure are literally identical.
    payload = json.loads(
        json.dumps(payload, sort_keys=True, separators=(",", ":"))
    )
    return sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _projection_digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _json_normalized(value: object) -> object:
    """Normalize integer JSON-object keys before replay equality checks."""

    return json.loads(json.dumps(value, sort_keys=True, separators=(",", ":")))


def _path_record(path: Path) -> dict[str, str]:
    return {
        "path": str(path.relative_to(ROOT)),
        "sha256": _file_sha256(path),
    }


def _normalized_certifier_source(source: str) -> str:
    """Normalize exactly the quoted RHS of both top-level outer freezes."""

    normalized = source
    for name in (
        "EXPECTED_CERTIFICATE_SHA256",
        "EXPECTED_ARCHIVE_FILE_SHA256",
    ):
        pattern = re.compile(
            rf'^(?P<prefix>{name}[ \t]*=[ \t]*)"[^"\n]*"(?P<suffix>[ \t]*)$',
            re.MULTILINE,
        )
        normalized, count = pattern.subn(
            lambda match: (
                f'{match.group("prefix")}"<OUTER_FREEZE>"'
                f'{match.group("suffix")}'
            ),
            normalized,
        )
        if count != 1:
            raise ValueError(
                f"certifier semantic source expected one top-level {name} assignment, got {count}"
            )
    return normalized


def _certifier_semantic_source_record() -> dict[str, str]:
    """Hash this certifier while excluding its two unavoidable outer refs."""

    normalized = _normalized_certifier_source(
        CERTIFIER_SOURCE.read_text(encoding="utf-8")
    )
    return {
        "path": str(CERTIFIER_SOURCE.relative_to(ROOT)),
        "sha256_excluding_outer_freeze_assignments": sha256(
            normalized.encode()
        ).hexdigest(),
    }


def _proof_sources() -> list[dict[str, str]]:
    return [
        _path_record(path)
        for path in (
            PROOF_SOURCE,
            RECORD35_SOURCE,
            FACET_AUDIT_SOURCE,
            HOMOLOGY_SOURCE,
            COXETER_LABEL_SOURCE,
            EXACT_WITNESS_SOURCE,
            RECORD46_SOURCE,
            RECORD92_SOURCE,
            LOCAL_LINK_SOURCE,
            MISSING_EDGES_SOURCE,
            NERVE_SUPPORT_SOURCE,
            INDEPENDENT_SOURCE,
            TEST_SOURCE,
        )
    ]


def _support_record() -> tuple[dict[str, object], str]:
    raw = SUPPORT.read_bytes()
    file_digest = sha256(raw).hexdigest()
    if file_digest != EXPECTED_SUPPORT_FILE_SHA256:
        raise ValueError("record-35 support archive file hash changed")
    source = json.loads(raw)
    verify_support_manifest(source)
    survivors = [
        row for row in source["records"] if row["survives_second_layer"]
    ]
    if TARGET_RECORD_INDEX >= len(survivors):
        raise ValueError("record-35 survivor index escaped the support archive")
    record = survivors[TARGET_RECORD_INDEX]
    if record.get("graph6") != TARGET_GRAPH6:
        raise ValueError("second-layer record-35 graph changed")
    if record.get("edges") != len(TARGET_MISSING_EDGES):
        raise ValueError("second-layer record-35 missing-edge count changed")
    support_core = record.get("ridge_support_core_masks")
    schedules = record.get("dehn_sommerville_schedule")
    if (
        not isinstance(support_core, list)
        or len(support_core) != 226
        or record.get("ridge_support_core_count") != len(support_core)
        or not isinstance(schedules, list)
        or len(schedules) != 5
    ):
        raise ValueError("record-35 support core or schedules changed")
    return record, file_digest


@lru_cache(maxsize=1)
def _hostile_replay() -> dict[str, object]:
    """Avoid repeating the expensive independent replay within one process."""

    return independent_replay()


def _validate_production(census: object) -> None:
    """Assert the exact counts that carry the compact whole-shadow theorem."""

    if not census.verified or census.exact_residual_incidence_orbit_count != 4:
        raise AssertionError("record-35 whole-shadow production audit failed")
    if (
        census.schedule.rows != EXPECTED_SCHEDULES
        or census.schedule.surviving_a != (0, 1, 2)
        or census.schedule.excluded_a != (3, 4)
    ):
        raise AssertionError("record-35 schedule partition changed")

    positive = census.positive_branches
    a1_censuses = positive.a1_completion_censuses
    if (
        not positive.verified
        or positive.individually_interval_feasible_a4_count != 490
        or positive.a1_viable_a4_count != 112
        or positive.a1_orbit_size_histogram != ((56, 2),)
        or len(a1_censuses) != 2
        or a1_censuses[0].scheduled_target_count != 1
        or a1_censuses[0].no_small_transversal_count != 1
        or a1_censuses[1].scheduled_target_count != 0
        or positive.a2_viable_pair_count != 28
        or positive.a2_orbit_size_histogram != ((28, 1),)
        or positive.a2_completion_census.no_small_transversal_count != 0
    ):
        raise AssertionError("record-35 positive-schedule census changed")

    a0 = census.a0
    if (
        not a0.verified
        or a0.independent_active_masks != EXPECTED_INDEPENDENT_ACTIVE_MASKS
        or a0.independent_active_mask_count != 18
        or a0.active_row_configuration_count != 1_221
        or a0.column_distribution_problem_count != 235
        or a0.total_distribution_path_count != 375_582
        or a0.repeated_row_path_count != 7_239
        or a0.small_transversal_path_count != 0
        or a0.blocker_valid_path_count != 368_343
        or a0.target_path_count != 20
        or a0.exact_target_family_count != 18
        or a0.topology_valid_family_indices != (0, 1, 2, 3, 4)
        or a0.all_target_full_graph_orbit_count != 11
        or a0.topology_valid_full_graph_orbit_count != 3
        or len(a0.label_audits) != 18
        or not all(row.compact_no_go for row in a0.label_audits)
    ):
        raise AssertionError("record-35 lossless a=0 census changed")

    a1_incidence = dict(census.a1_incidence_audit)
    if (
        a1_incidence.get("minimal_transversal_size_histogram") != {"5": 100}
        or a1_incidence.get("facet_count") != 100
        or a1_incidence.get("ridge_degree_histogram") != {"2": 450}
        or a1_incidence.get("face_vector")
        != [14, 86, 309, 725, 1165, 1295, 970, 450, 100]
        or a1_incidence.get("passes_connected_homology_audit") is not True
    ):
        raise AssertionError("record-35 a=1 incidence audit changed")

    a1_label = census.a1_label
    if (
        not a1_label.unrestricted_compact_no_go
        or a1_label.rank_nine_vertex_count != 100
        or (
            a1_label.order_four_placement_count,
            a1_label.overlapping_order_five_placement_count,
            a1_label.disjoint_order_five_placement_count,
        )
        != (108, 420, 420)
        or a1_label.left_cluster_passing_state_count != 84
        or a1_label.sparse_block_state_count != 35_280
        or a1_label.label_three_test_count != 1_058_400
        or a1_label.label_three_passing_count != 5_040
        or a1_label.stronger_single_edge_survivors != ((4, 0), (5, 0))
        or a1_label.exceptional_class_histogram
        != (("A5_A5_type_1_1", 4_320), ("A5_A5_type_4_1", 720))
        or tuple(row.inertia for row in a1_label.exact_index_two_witnesses)
        != ((8, 2, 0), (8, 2, 0))
        or tuple(
            row.multiplicity_per_support_gate
            for row in a1_label.exact_index_two_witnesses
        )
        != (4_320, 720)
    ):
        raise AssertionError("record-35 unrestricted a=1 label gate changed")


def _validate_hostile(hostile: Mapping[str, object]) -> None:
    if hostile.get("summary_sha256") != EXPECTED_HOSTILE_SUMMARY_SHA256:
        raise AssertionError("record-35 independent replay digest changed")
    positive = hostile.get("positive_branches")
    a0 = hostile.get("a0_census")
    a0_labels = hostile.get("a0_all_target_label_audits")
    a1_incidence = hostile.get("a1_incidence_audit")
    a1_label = hostile.get("a1_label_audit")
    theorem = hostile.get("theorem")
    if not all(
        isinstance(row, Mapping)
        for row in (positive, a0, a1_incidence, a1_label, theorem)
    ) or not isinstance(a0_labels, list):
        raise AssertionError("record-35 hostile replay is malformed")
    assert isinstance(positive, Mapping)
    assert isinstance(a0, Mapping)
    assert isinstance(a1_incidence, Mapping)
    assert isinstance(a1_label, Mapping)
    assert isinstance(theorem, Mapping)
    a0_topology = a0.get("full_topology_audits")
    a1_topology = a1_incidence.get("full_topology_audit")
    if (
        hostile.get("strict_schedule_survivors") != [0, 1, 2]
        or positive.get("a4_count") != 490
        or positive.get("a1_orbit_sizes") != [56, 56]
        or positive.get("a2_viable_count") != 28
        or positive.get("a2_orbit_sizes") != [28]
        or a0.get("independent_active_masks")
        != list(EXPECTED_INDEPENDENT_ACTIVE_MASKS)
        or a0.get("independent_active_mask_size_histogram")
        != {0: 1, 1: 5, 2: 8, 3: 4}
        or a0.get("independent_active_triples")
        != [[0, 1, 2], [0, 2, 12], [1, 2, 11], [2, 11, 12]]
        or a0.get("independent_active_mask_count") != 18
        or a0.get("active_row_configuration_count") != 1_221
        or a0.get("column_distribution_problem_count") != 235
        or a0.get("total_distribution_path_count") != 375_582
        or a0.get("repeated_row_path_count") != 7_239
        or a0.get("small_transversal_path_count") != 0
        or a0.get("blocker_valid_path_count") != 368_343
        or a0.get("target_path_count") != 20
        or a0.get("exact_target_family_count") != 18
        or a0.get("topology_valid_family_indices") != [0, 1, 2, 3, 4]
        or a0.get("all_target_orbit_count") != 11
        or a0.get("topology_valid_orbit_count") != 3
        or not isinstance(a0_topology, list)
        or len(a0_topology) != 5
        or any(
            not isinstance(row, Mapping)
            or row.get("topology_valid") is not True
            or row.get("all_links_have_spherical_euler_characteristic") is not True
            or row.get("all_links_are_mod2_homology_spheres") is not True
            for row in a0_topology
        )
        or len(a0_labels) != 18
        or not all(
            isinstance(row, Mapping) and row.get("compact_no_go") is True
            for row in a0_labels
        )
        or a1_incidence.get("minimal_transversal_size_histogram") != {5: 100}
        or a1_incidence.get("ridge_degree_histogram") != {2: 450}
        or not isinstance(a1_topology, Mapping)
        or a1_topology.get("topology_valid") is not True
        or a1_topology.get("all_links_have_spherical_euler_characteristic") is not True
        or a1_topology.get("all_links_are_mod2_homology_spheres") is not True
        or a1_label.get("placement_counts") != [108, 420, 420]
        or a1_label.get("left_cluster_passing_count") != 84
        or a1_label.get("sparse_state_count") != 35_280
        or a1_label.get("sparse_failure_count") != 0
        or a1_label.get("label_three_test_count") != 1_058_400
        or a1_label.get("label_three_passing_count") != 5_040
        or a1_label.get("stronger_survivors") != {4: 0, 5: 0}
        or a1_label.get("exceptional_class_histogram")
        != {"A5_A5_type_1_1": 4_320, "A5_A5_type_4_1": 720}
        or [row.get("inertia") for row in a1_label.get("exact_witnesses", [])]
        != [[8, 2, 0], [8, 2, 0]]
        or a1_label.get("compact_no_go") is not True
        or theorem.get("whole_graph_shadow_compact_coxeter_excluded") is not True
        or theorem.get("exact_topology_valid_residual_incidence_orbit_count") != 4
        or theorem.get("finite_volume_noncompact_chamber_excluded") is not False
    ):
        raise AssertionError("record-35 hostile replay did not close the theorem")


def build() -> dict[str, object]:
    support_record, support_digest = _support_record()
    schedules = support_record["dehn_sommerville_schedule"]
    support_core = support_record["ridge_support_core_masks"]
    schedules_json = json.dumps(schedules, sort_keys=True, separators=(",", ":"))
    support_core_json = json.dumps(
        support_core, sort_keys=True, separators=(",", ":")
    )
    census = audit_record35_all_completions(schedules_json, support_core_json)
    _validate_production(census)

    hostile = _hostile_replay()
    _validate_hostile(hostile)
    hostile_a0 = hostile["a0_census"]
    hostile_a1 = hostile["a1_label_audit"]
    assert isinstance(hostile_a0, Mapping)
    assert isinstance(hostile_a1, Mapping)

    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record35_whole_shadow_compact_coxeter_no_go",
        "claim": (
            "Under the compact H9/14 missing-face reductions, all five "
            "Dehn--Sommerville schedules for record 35 are exhaustive. "
            "Strict vertex-link positivity removes a=3,4; complete missing-"
            "four and forced-shadow censuses remove a=2 and leave one a=1 "
            "incidence orbit. A lossless all-cap-feasible-independent-active-"
            "row census "
            "leaves three a=0 incidence orbits passing the full necessary "
            "topology audit. Complete "
            "Lanner-placement joins reject all eighteen a=0 facet-count "
            "targets, while a monotonicity and exact-index-two gate rejects "
            "the a=1 orbit. Hence the whole graph shadow is not compact "
            "Coxeter."
        ),
        "scope": {
            "record_index": TARGET_RECORD_INDEX,
            "whole_graph_shadow_compact_coxeter_excluded": True,
            "all_five_dehn_sommerville_schedules_covered": True,
            "exact_topology_valid_residual_incidence_orbits": 4,
            "exact_topology_count_independently_replayed": True,
            "topology_audit_is_necessary_not_a_PL_sphere_certificate": True,
            "a0_topology_valid_residual_orbits": 3,
            "a1_topology_valid_residual_orbits": 1,
            "compact_only": True,
            "finite_volume_noncompact_chamber_excluded": False,
            "ideal_or_parabolic_faces_enumerated": False,
        },
        "nonclaims": [
            "No finite-volume noncompact chamber is excluded.",
            "No ideal-vertex or maximal-parabolic census is performed.",
            "No other H9/14 graph shadow is excluded by this certificate.",
            "The mod-two all-link homology audit is not a PL-sphere certificate and does not detect odd torsion or simple connectivity.",
            "Arbitrary ordinary completions with minimal nonfaces above the compact Coxeter order bound are not classified.",
        ],
        "proof_dependencies": [
            "Burcroff's no-missing-triples theorem for compact H9/14 nerves",
            "every compact Coxeter missing face is Lanner of order at most five",
            "Dehn--Sommerville equations for the simplicial polar boundary",
            "strict generalized-lower-bound and Macaulay inequalities in every vertex link",
            "facet complements are the blocker of the full minimal-nonface clutter",
            "complete enumeration of every independent active row subset allowed by the exact incidence caps and exact isolate-column multiplicities",
            "complete Lanner placement joins tested on all actual vertices",
            "Perron--Frobenius monotonicity, Cauchy interlacing, and exact negative-index-two principal blocks",
        ],
        "support_source": {
            "path": str(SUPPORT.relative_to(ROOT)),
            "file_sha256": support_digest,
            "second_layer_survivor_index": TARGET_RECORD_INDEX,
            "graph6": TARGET_GRAPH6,
            "missing_edges": [list(edge) for edge in TARGET_MISSING_EDGES],
            "ridge_support_core_count": len(support_core),
            "ridge_support_core_sha256": _projection_digest(support_core),
            "dehn_sommerville_schedule_sha256": _projection_digest(schedules),
        },
        "proof_sources": _proof_sources(),
        "certifier_source": _certifier_semantic_source_record(),
        "whole_shadow_production_audit": census.to_record(),
        "independent_hostile_replay": {
            "path": str(INDEPENDENT_SOURCE.relative_to(ROOT)),
            "source_sha256": _file_sha256(INDEPENDENT_SOURCE),
            "summary_sha256": hostile["summary_sha256"],
            "complete_active_row_scope": hostile_a0[
                "complete_active_row_scope"
            ],
            "independent_active_masks": hostile_a0[
                "independent_active_masks"
            ],
            "independent_active_mask_size_histogram": hostile_a0[
                "independent_active_mask_size_histogram"
            ],
            "independent_active_triples": hostile_a0[
                "independent_active_triples"
            ],
            "a0_total_distribution_path_count": hostile_a0[
                "total_distribution_path_count"
            ],
            "a0_exact_target_family_count": hostile_a0[
                "exact_target_family_count"
            ],
            "a0_topology_valid_orbit_count": hostile_a0[
                "topology_valid_orbit_count"
            ],
            "a0_full_topology_audit_sha256": _projection_digest(
                hostile_a0["full_topology_audits"]
            ),
            "a1_full_topology_audit_sha256": _projection_digest(
                hostile["a1_incidence_audit"]["full_topology_audit"]
            ),
            "exact_topology_valid_residual_incidence_orbit_count": 4,
            "a1_label_three_test_count": hostile_a1[
                "label_three_test_count"
            ],
            "a1_label_three_passing_count": hostile_a1[
                "label_three_passing_count"
            ],
            "a1_exact_witness_inertias": [
                row["inertia"] for row in hostile_a1["exact_witnesses"]
            ],
            "whole_graph_shadow_compact_coxeter_excluded": True,
            "finite_volume_noncompact_chamber_excluded": False,
        },
        "theorem": {
            "dimension": 9,
            "facet_count": 14,
            "record_index": TARGET_RECORD_INDEX,
            "whole_graph_shadow_compact_coxeter_excluded": True,
            "finite_volume_noncompact_chamber_excluded": False,
        },
    }
    manifest["certificate_sha256"] = _canonical_digest(manifest)
    return manifest


def _verify_semantics(document: Mapping[str, object]) -> None:
    if document.get("schema_version") != "1.0.0" or document.get("kind") != (
        "h9_14_record35_whole_shadow_compact_coxeter_no_go"
    ):
        raise ValueError("wrong record-35 whole-shadow certificate kind")
    scope = document.get("scope")
    theorem = document.get("theorem")
    production = document.get("whole_shadow_production_audit")
    hostile = document.get("independent_hostile_replay")
    if not all(
        isinstance(row, Mapping)
        for row in (scope, theorem, production, hostile)
    ):
        raise ValueError("record-35 whole-shadow certificate is malformed")
    assert isinstance(scope, Mapping)
    assert isinstance(theorem, Mapping)
    assert isinstance(production, Mapping)
    assert isinstance(hostile, Mapping)
    if (
        scope.get("record_index") != 35
        or scope.get("all_five_dehn_sommerville_schedules_covered") is not True
        or scope.get("exact_topology_valid_residual_incidence_orbits") != 4
        or scope.get("exact_topology_count_independently_replayed") is not True
        or scope.get("topology_audit_is_necessary_not_a_PL_sphere_certificate")
        is not True
        or scope.get("a0_topology_valid_residual_orbits") != 3
        or scope.get("a1_topology_valid_residual_orbits") != 1
        or scope.get("compact_only") is not True
        or scope.get("finite_volume_noncompact_chamber_excluded") is not False
        or scope.get("ideal_or_parabolic_faces_enumerated") is not False
        or theorem
        != {
            "dimension": 9,
            "facet_count": 14,
            "record_index": 35,
            "whole_graph_shadow_compact_coxeter_excluded": True,
            "finite_volume_noncompact_chamber_excluded": False,
        }
        or production.get("verified") is not True
        or production.get("exact_topology_valid_residual_incidence_orbit_count")
        != 4
        or production.get("whole_graph_shadow_compact_coxeter_excluded")
        is not True
        or production.get("finite_volume_noncompact_chamber_excluded") is not False
        or hostile.get("whole_graph_shadow_compact_coxeter_excluded") is not True
        or hostile.get("exact_topology_valid_residual_incidence_orbit_count") != 4
        or hostile.get("finite_volume_noncompact_chamber_excluded") is not False
    ):
        raise ValueError("record-35 certificate theorem semantics changed")
    if document.get("proof_sources") != _proof_sources():
        raise ValueError("record-35 proof-source closure changed")
    if document.get("certifier_source") != _certifier_semantic_source_record():
        raise ValueError("record-35 certifier semantic-source closure changed")
    support_record, support_digest = _support_record()
    expected_support = {
        "path": str(SUPPORT.relative_to(ROOT)),
        "file_sha256": support_digest,
        "second_layer_survivor_index": TARGET_RECORD_INDEX,
        "graph6": TARGET_GRAPH6,
        "missing_edges": [list(edge) for edge in TARGET_MISSING_EDGES],
        "ridge_support_core_count": len(
            support_record["ridge_support_core_masks"]
        ),
        "ridge_support_core_sha256": _projection_digest(
            support_record["ridge_support_core_masks"]
        ),
        "dehn_sommerville_schedule_sha256": _projection_digest(
            support_record["dehn_sommerville_schedule"]
        ),
    }
    if document.get("support_source") != expected_support:
        raise ValueError("record-35 support-source closure changed")


def verify(document: Mapping[str, object]) -> dict[str, object]:
    _verify_semantics(document)
    digest = _canonical_digest(document)
    if document.get("certificate_sha256") != digest:
        raise ValueError("record-35 whole-shadow certificate digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 not in ("TO_BE_FROZEN", digest):
        raise ValueError("record-35 whole-shadow certificate is not frozen")
    rebuilt = build()
    if _json_normalized(document) != _json_normalized(rebuilt):
        raise ValueError("record-35 whole-shadow certificate does not replay")
    return {
        "record_index": TARGET_RECORD_INDEX,
        "exact_residual_incidence_orbits": rebuilt["scope"][
            "exact_topology_valid_residual_incidence_orbits"
        ],
        "whole_graph_shadow_compact_coxeter_excluded": rebuilt["theorem"][
            "whole_graph_shadow_compact_coxeter_excluded"
        ],
        "finite_volume_noncompact_chamber_excluded": rebuilt["theorem"][
            "finite_volume_noncompact_chamber_excluded"
        ],
        "certificate_sha256": digest,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-certificate", action="store_true")
    args = parser.parse_args(argv)
    if args.write_certificate:
        document = build()
        OUTPUT.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    else:
        document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    if EXPECTED_ARCHIVE_FILE_SHA256 not in (
        "TO_BE_FROZEN",
        _file_sha256(OUTPUT),
    ):
        raise ValueError("record-35 certificate archive file hash changed")
    print(json.dumps({"output": str(OUTPUT), **verify(document)}, indent=2))


if __name__ == "__main__":
    main()
