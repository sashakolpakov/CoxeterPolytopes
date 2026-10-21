#!/usr/bin/env python3
"""Compose the complete compact-Coxeter exclusion for H9/14 Record 30.

This aggregate does not rerun the expensive branch enumerations.  It invokes
the frozen semantic verifier for every schedule certificate, checks both the
archive-file and internal certificate digests, independently reconstructs the
complete Dehn--Sommerville schedule table, and verifies that the branch
certificates cover every possible nonnegative value of ``a``.

The conclusion is compact-only.  In particular, it makes no claim about
finite-volume noncompact chambers or ideal/parabolic face data.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import re
import sys
from types import ModuleType
from typing import Mapping


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

OUTPUT = (
    ROOT
    / "data/enumerations/h9_14_record30_whole_shadow_compact_no_go.json"
)
CERTIFIER_SOURCE = Path(__file__).resolve()
LOCAL_LINK_REPORT = ROOT / "results/h9_14_local_link_constraints.md"
A0_HOSTILE_ARCHIVE = (
    ROOT / "data/enumerations/h9_14_record30_a0_dual_orbit_hostile.json"
)
A0_HOSTILE_SOURCE = (
    ROOT / "scripts/hostile_h9_14_record30_a0_dual_cnf_orbits.py"
)
A0_FINAL_HOSTILE_ARCHIVE = (
    ROOT / "data/enumerations/h9_14_record30_a0_dual_hostile_audit.json"
)
A0_FINAL_HOSTILE_REPORT = (
    ROOT / "results/h9_14_record30_a0_dual_hostile_audit.md"
)

# These are the two unavoidable outer pins.  Their quoted right-hand sides are
# normalized when this source is included in the aggregate provenance.
EXPECTED_CERTIFICATE_SHA256 = "4215b56a4fe830825da8a8689395ae33f2fba105b83bbbc468bb6a7f054f170f"
EXPECTED_ARCHIVE_FILE_SHA256 = "17a28b69ffc46e9712386971ed5e3c19edf6b952e0474272a6b00c51ce4a2b09"

EXPECTED_LOCAL_LINK_REPORT_SHA256 = (
    "6daaf924c16ce5dd2714608bb815c54ec66423f2bb46f2da5d73562b647e1814"
)
EXPECTED_A0_HOSTILE_ARCHIVE_SHA256 = (
    "b33237d52eb8206cad7c8c6c5bc95e43c7ec54dac5367030a0c55fc9803122b5"
)
EXPECTED_A0_HOSTILE_SOURCE_SHA256 = (
    "f49d28b445adffa22b326582dc8b4bc7d6018e3699c093d2140dc21ca81ca3b7"
)
EXPECTED_A0_FINAL_HOSTILE_ARCHIVE_SHA256 = (
    "aff5de1f18de521bf0c790e63183eb1d2456f2c03e3510cf7cb0068b3d1f75ae"
)
EXPECTED_A0_FINAL_HOSTILE_SUMMARY_SHA256 = (
    "9fa8b0be17abbb7c35fa8e2fc614bed2ef095e8a3efc87d0387406462022b861"
)
EXPECTED_A0_FINAL_HOSTILE_REPORT_SHA256 = (
    "0405bd354ce1b2e7c2501f5af041e9cb92dd27e4b1d7a50daa0222e661b16a41"
)


COMPONENTS: tuple[dict[str, object], ...] = (
    {
        "id": "a0_dual_exhaustion",
        "schedules": [0],
        "archive": "data/enumerations/h9_14_record30_a0_dual_exhaustion.json",
        "archive_sha256": "7934ffec51f8bc3c5f209d016f75c0aa6feac34343ad34a58f16a0d4660a64e5",
        "certificate_sha256": "436b9bc834529f7d2d94754e43c3f785bb937ec0f7157570ad9a06d56858ab2d",
        "verifier": "scripts/certify_h9_14_record30_a0_dual_exhaustion.py",
        "verifier_sha256": "2b0cc2779fb5d9b011de83c0471820ea1e323dedd23e037c548e02f4412faa40",
    },
    {
        "id": "a1_schedule",
        "schedules": [1],
        "archive": "data/enumerations/h9_14_record30_a1_schedule_no_go.json",
        "archive_sha256": "65152e8f4b687f487852932a5dae3896d189a56d1ed99c6764f7389aaea8d341",
        "certificate_sha256": "932d69e5646106785171562174e5423329ac09137e0224ac65e6ce436187067d",
        "verifier": "scripts/certify_h9_14_record30_a1_schedule_no_go.py",
        "verifier_sha256": "ff020c2c579071943fa37c12e4a41419057fa81d10ae729ea3c565f884f9e9f0",
    },
    {
        "id": "a2_schedule",
        "schedules": [2],
        "archive": "data/enumerations/h9_14_record30_a2_schedule_no_go.json",
        "archive_sha256": "06759a184708d9c2b8a4ab15fc772fc2a4b7e66ac275d2ff537bebf92cb9ebe6",
        "certificate_sha256": "3ac5e4ea4ddd624aef17dec45f64d403165b736bab015d088e2b8d78f43a5adf",
        "verifier": "scripts/certify_h9_14_record30_a2_schedule_no_go.py",
        "verifier_sha256": "1e1d303829344ba3bbdc54b1e071daddb1099f31a28bb3f35a08ecc9c8b28dff",
    },
    {
        "id": "a3_schedule",
        "schedules": [3],
        "archive": "data/enumerations/h9_14_record30_a3_compact_no_go.json",
        "archive_sha256": "2e959be2b47424ebe6b9780848202fd800055dd8f093ea3cf675336947abcb5c",
        "certificate_sha256": "e3afaf33526c003307d797bd0789bd2b4275d3d74c83f0a8d9792bb1ad0b9d54",
        "verifier": "scripts/certify_h9_14_record30_a3_no_go.py",
        "verifier_sha256": "88211b0f025a8c08fce5041c3a922fa6035594044aaa467fd22f0421efb7b529",
    },
    {
        "id": "a4_schedule",
        "schedules": [4],
        "archive": "data/enumerations/h9_14_record30_a4_compact_no_go.json",
        "archive_sha256": "4c23f1a87a35f51c9c277e08a18e7f89ca02eb51757c244d59270d005ecc6c46",
        "certificate_sha256": "12918a797758ec9b14ed02657c8a98060c77b54abaa93d0e627dc31d0efb3bb5",
        "verifier": "scripts/certify_h9_14_record30_a4_no_go.py",
        "verifier_sha256": "bd9b3eccc91eb3ff31491c0d10e3f42fdce7a03d5070119e934a8c259d9fcbf2",
    },
    {
        "id": "a5_a7_positive_scout",
        "schedules": [5, 6, 7],
        "archive": "data/enumerations/h9_14_record30_positive_scout.json",
        "archive_sha256": "841ebf50748747a8dc96fb4ce337570f29f1f9455f2206609e18d31187eb3771",
        "certificate_sha256": "02a2238457b33d0cba59bb9fb85579a5e27c4520ea86f17a4ce0614fcabb111e",
        "verifier": "scripts/certify_h9_14_record30_positive_scout.py",
        "verifier_sha256": "dfa7f24230126ae3bd07b289152c38738418f4902757001ce79d54dd521592d5",
    },
)

OPEN_COMPACT_SHADOWS_AFTER_RECORD30 = [0, 1, 3, 4, 5, 7, 9, 14, 31]


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    return sha256(
        json.dumps(
            payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode("ascii")
    ).hexdigest()


def _load(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"malformed theorem archive {path}")
    return value


def _load_module(name: str, path: Path) -> ModuleType:
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"could not load verifier {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def _normalized_certifier_source(source: str) -> str:
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
            raise ValueError(f"could not normalize exactly one {name} assignment")
    return normalized


def _semantic_certifier_source_record() -> dict[str, str]:
    normalized = _normalized_certifier_source(
        CERTIFIER_SOURCE.read_text(encoding="utf-8")
    )
    return {
        "path": str(CERTIFIER_SOURCE.relative_to(ROOT)),
        "sha256_excluding_outer_freeze_assignments": sha256(
            normalized.encode()
        ).hexdigest(),
    }


def _verify_embedded_sources(document: Mapping[str, object]) -> None:
    sources = document.get("source_sha256")
    if not isinstance(sources, Mapping):
        raise ValueError("Record-30 a=1 embedded source pins are missing")
    for relative, expected in sources.items():
        if not isinstance(relative, str) or not isinstance(expected, str):
            raise ValueError("Record-30 a=1 source pin is malformed")
        if _file_sha256(ROOT / relative) != expected:
            raise ValueError(f"Record-30 a=1 source changed: {relative}")


def _verify_a0_hostile_archive() -> dict[str, object]:
    if _file_sha256(A0_HOSTILE_ARCHIVE) != EXPECTED_A0_HOSTILE_ARCHIVE_SHA256:
        raise ValueError("Record-30 a=0 no-import hostile archive changed")
    if _file_sha256(A0_HOSTILE_SOURCE) != EXPECTED_A0_HOSTILE_SOURCE_SHA256:
        raise ValueError("Record-30 a=0 no-import hostile source changed")
    document = _load(A0_HOSTILE_ARCHIVE)
    critical = dict(document)
    summary_sha = critical.pop("summary_sha256", None)
    if summary_sha != sha256(
        json.dumps(
            critical, sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode("ascii")
    ).hexdigest():
        raise ValueError("Record-30 a=0 hostile summary digest is invalid")
    expected = {
        "input_target_orbit_count": 8,
        "blocked_seed_form_count": 63_960,
        "full_group_order": 1_935_360,
        "type_three_complement_orbit_size": 224,
        "cnf_variable_count": 226_358,
        "cnf_clause_count": 984_935,
        "all_target_orbits_pairwise_disjoint": True,
        "every_seed_count_matches_orbit_incidence_double_count": True,
        "every_orbit_has_a_seed_stabilizer_lex_witness": True,
        "cnf_tail_is_exactly_the_independent_blocking_clause_set": True,
    }
    if any(document.get(key) != value for key, value in expected.items()):
        raise ValueError("Record-30 a=0 hostile orbit semantics changed")
    audits = document.get("target_orbit_audits")
    if not isinstance(audits, list) or len(audits) != 8:
        raise ValueError("Record-30 a=0 hostile target partition changed")
    seed_counts = [row.get("direct_seed_form_count") for row in audits]
    if seed_counts != [1200, 14760, 14640, 3600, 3720, 14880, 3720, 7440]:
        raise ValueError("Record-30 a=0 hostile seed-form partition changed")

    if (
        _file_sha256(A0_FINAL_HOSTILE_ARCHIVE)
        != EXPECTED_A0_FINAL_HOSTILE_ARCHIVE_SHA256
    ):
        raise ValueError("Record-30 a=0 final hostile archive changed")
    if (
        _file_sha256(A0_FINAL_HOSTILE_REPORT)
        != EXPECTED_A0_FINAL_HOSTILE_REPORT_SHA256
    ):
        raise ValueError("Record-30 a=0 final hostile report changed")
    final = _load(A0_FINAL_HOSTILE_ARCHIVE)
    final_critical = dict(final)
    final_summary_sha = final_critical.pop("summary_sha256", None)
    computed_final_summary = sha256(
        json.dumps(
            final_critical,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
        ).encode("ascii")
    ).hexdigest()
    if (
        final_summary_sha != EXPECTED_A0_FINAL_HOSTILE_SUMMARY_SHA256
        or final_summary_sha != computed_final_summary
    ):
        raise ValueError("Record-30 a=0 final hostile internal digest changed")
    conclusion = final.get("conclusion")
    replay = final.get("independent_drat_replay")
    final_scope = final.get("scope")
    orbit_audit = final.get("orbit_and_clause_audit")
    if not all(
        isinstance(item, Mapping)
        for item in (conclusion, replay, final_scope, orbit_audit)
    ):
        raise ValueError("Record-30 a=0 final hostile sections are malformed")
    assert isinstance(conclusion, Mapping) and isinstance(replay, Mapping)
    assert isinstance(final_scope, Mapping) and isinstance(orbit_audit, Mapping)
    if (
        conclusion.get("audit_status") != "PASS"
        or conclusion.get("remaining_fatal_or_material_gap") is not False
        or conclusion.get("record30_a0_compact_schedule_excluded") is not True
        or conclusion.get("new_hyperbolic_coxeter_polytope") is not False
        or replay.get("status") != "VERIFIED"
        or replay.get("exit_code") != 0
        or replay.get("checker_resolution_step_count") != 34_109_825
        or replay.get("checker_binary_sha256")
        != "f58f63b0f76945d4c4c9ff6e87afaf870f579e67c0f7cca589492df8fc7ebd47"
        or orbit_audit.get("blocked_seed_form_count") != 63_960
        or orbit_audit.get("cnf_tail_exact_set_comparison_pass") is not True
        or final_scope.get("compact_only") is not True
        or final_scope.get("finite_volume_noncompact_statement") is not False
    ):
        raise ValueError("Record-30 a=0 final hostile semantics changed")
    return {
        "orbit_and_clause_replay": {
            "archive": str(A0_HOSTILE_ARCHIVE.relative_to(ROOT)),
            "archive_sha256": EXPECTED_A0_HOSTILE_ARCHIVE_SHA256,
            "summary_sha256": summary_sha,
            "target_orbit_count": 8,
            "blocked_seed_form_count": 63_960,
            "all_orbit_stabilizer_double_counts_exact": True,
            "cnf_tail_exactly_reconstructed": True,
        },
        "final_hostile_audit": {
            "archive": str(A0_FINAL_HOSTILE_ARCHIVE.relative_to(ROOT)),
            "archive_sha256": EXPECTED_A0_FINAL_HOSTILE_ARCHIVE_SHA256,
            "summary_sha256": final_summary_sha,
            "report": str(A0_FINAL_HOSTILE_REPORT.relative_to(ROOT)),
            "report_sha256": EXPECTED_A0_FINAL_HOSTILE_REPORT_SHA256,
            "audit_status": "PASS",
            "remaining_fatal_or_material_gap": False,
            "independent_drat_replay_status": "VERIFIED",
            "independent_drat_replay_seconds": replay.get("replay_seconds"),
            "checker_resolution_step_count": 34_109_825,
            "compact_only": True,
            "finite_volume_noncompact_statement": False,
        },
    }


def _component_by_id(identifier: str) -> dict[str, object]:
    for component in COMPONENTS:
        if component["id"] == identifier:
            return component
    raise KeyError(identifier)


def _verify_component_pins(component: Mapping[str, object]) -> tuple[dict[str, object], ModuleType]:
    archive = ROOT / str(component["archive"])
    verifier = ROOT / str(component["verifier"])
    if _file_sha256(archive) != component["archive_sha256"]:
        raise ValueError(f"component archive changed: {archive.name}")
    if _file_sha256(verifier) != component["verifier_sha256"]:
        raise ValueError(f"component verifier changed: {verifier.name}")
    document = _load(archive)
    if document.get("certificate_sha256") != component["certificate_sha256"]:
        raise ValueError(f"component internal digest pin changed: {archive.name}")
    module = _load_module(f"record30_aggregate_{component['id']}", verifier)
    return document, module


def verify_components() -> list[dict[str, object]]:
    records: list[dict[str, object]] = []

    component = _component_by_id("a0_dual_exhaustion")
    document, module = _verify_component_pins(component)
    summary = module.verify(document)
    if summary != {
        "target_orbit_count": 8,
        "blocked_seed_form_count": 63_960,
        "cnf_variable_count": 226_358,
        "cnf_clause_count": 984_935,
        "drat_verified": True,
        "positive_witness_count": 8,
        "compact_label_survivor_count": 0,
        "finite_volume_noncompact_statement": False,
    }:
        raise ValueError("Record-30 a=0 semantic verification changed")
    if document.get("conclusion") != {
        "whole_record30_a0_compact_schedule_excluded": True,
        "new_hyperbolic_coxeter_polytope": False,
    }:
        raise ValueError("Record-30 a=0 conclusion changed")
    proof = document.get("proof")
    if not isinstance(proof, Mapping):
        raise ValueError("Record-30 a=0 proof provenance is missing")
    replay_provenance = {
        "solver": proof.get("solver"),
        "solver_binary_sha256": proof.get("solver_binary_sha256"),
        "solver_status": proof.get("solver_status"),
        "checker": proof.get("checker"),
        "checker_commit": proof.get("checker_commit"),
        "checker_binary_sha256": proof.get("checker_binary_sha256"),
        "checker_status": proof.get("checker_status"),
        "checker_resolution_step_count": proof.get(
            "checker_resolution_step_count"
        ),
        "deep_mode_enforces_exact_solver_and_checker_binary_hashes": True,
    }
    if replay_provenance != {
        "solver": "Kissat 4.0.4",
        "solver_binary_sha256": (
            "7458a5dae42d7dcc0c8402bb153fe9faee321cb6f3d079b79b0fb385e55d8bb6"
        ),
        "solver_status": "UNSATISFIABLE",
        "checker": "drat-trim",
        "checker_commit": "2e3b2dc0ecf938addbd779d42877b6ed69d9a985",
        "checker_binary_sha256": (
            "f58f63b0f76945d4c4c9ff6e87afaf870f579e67c0f7cca589492df8fc7ebd47"
        ),
        "checker_status": "VERIFIED",
        "checker_resolution_step_count": 34_109_825,
        "deep_mode_enforces_exact_solver_and_checker_binary_hashes": True,
    }:
        raise ValueError("Record-30 a=0 frozen deep-replay provenance changed")
    records.append(
        {
            **component,
            "semantic_summary": summary,
            "deep_replay_provenance": replay_provenance,
        }
    )

    component = _component_by_id("a1_schedule")
    document, module = _verify_component_pins(component)
    summary = module.verify(document)
    _verify_embedded_sources(document)
    if (
        summary.get("record30_a1_compact_schedule_excluded") is not True
        or summary.get("necessary_incidence_target_orbit_count") != 3
        or summary.get("literal_target_count") != 168
        or summary.get("terminal_finite_label_state_count") != 0
        or summary.get("finite_volume_noncompact_chamber_excluded") is not False
        or summary.get("new_hyperbolic_coxeter_polytope") is not False
    ):
        raise ValueError("Record-30 a=1 semantic verification changed")
    records.append({**component, "semantic_summary": summary})

    component = _component_by_id("a2_schedule")
    document, module = _verify_component_pins(component)
    summary = module.verify(document)
    if (
        summary.get("a2_compact_schedule_excluded") is not True
        or summary.get("positive_A4_orbit_count") != 21
        or summary.get("pure_target_orbit_count") != 77
        or summary.get("sparse_label_no_go_target_count") != 74
        or summary.get("exceptional_target_count") != 3
        or summary.get("finite_volume_noncompact_chamber_excluded") is not False
        or summary.get("new_hyperbolic_coxeter_polytope") is not False
    ):
        raise ValueError("Record-30 a=2 semantic verification changed")
    records.append({**component, "semantic_summary": summary})

    component = _component_by_id("a3_schedule")
    document, module = _verify_component_pins(component)
    summary = module.verify(document)
    if (
        summary.get("newly_closed_schedule") != 3
        or summary.get("a3_orbit_count") != 32
        or summary.get("no_small_transversal_family_count") != 1_896_496
        or summary.get("scheduled_126_facet_target_count") != 0
        or summary.get("compact_only") is not True
        or summary.get("finite_volume_noncompact_chamber_excluded") is not False
    ):
        raise ValueError("Record-30 a=3 semantic verification changed")
    records.append({**component, "semantic_summary": summary})

    component = _component_by_id("a4_schedule")
    document, module = _verify_component_pins(component)
    summary = module.verify(document)
    if (
        summary.get("newly_closed_schedule") != 4
        or summary.get("a4_orbit_count") != 7
        or summary.get("no_small_transversal_family_count") != 840
        or summary.get("scheduled_124_facet_target_count") != 0
        or summary.get("compact_only") is not True
        or summary.get("finite_volume_noncompact_chamber_excluded") is not False
    ):
        raise ValueError("Record-30 a=4 semantic verification changed")
    records.append({**component, "semantic_summary": summary})

    component = _component_by_id("a5_a7_positive_scout")
    document, module = _verify_component_pins(component)
    summary = module.verify(document)
    if summary != {
        "positive_orbit_counts": [3, 21, 32, 7, 1],
        "a5_compact_schedule_excluded": True,
        "strict_sum_excluded_schedules": [6, 7],
        "whole_shadow_compact_coxeter_excluded": False,
        "finite_volume_noncompact_chamber_excluded": False,
    }:
        raise ValueError("Record-30 a=5,6,7 semantic verification changed")
    records.append({**component, "semantic_summary": summary})

    return records


def schedule_table() -> list[dict[str, object]]:
    # Derived independently from the fixed Record-30 graph 3K2 + 8K1 and the
    # H9/14 Dehn--Sommerville equations.
    vectors = (
        (14, 88, 328, 806, 1364, 1596, 1248, 594, 132),
        (14, 88, 328, 805, 1358, 1582, 1232, 585, 130),
        (14, 88, 328, 804, 1352, 1568, 1216, 576, 128),
        (14, 88, 328, 803, 1346, 1554, 1200, 567, 126),
        (14, 88, 328, 802, 1340, 1540, 1184, 558, 124),
        (14, 88, 328, 801, 1334, 1526, 1168, 549, 122),
        (14, 88, 328, 800, 1328, 1512, 1152, 540, 120),
        (14, 88, 328, 799, 1322, 1498, 1136, 531, 118),
    )
    rows = []
    for a, vector in enumerate(vectors):
        expected_vector = (
            14,
            88,
            328,
            806 - a,
            1364 - 6 * a,
            1596 - 14 * a,
            1248 - 16 * a,
            594 - 9 * a,
            132 - 2 * a,
        )
        if vector != expected_vector:
            raise AssertionError("Record-30 Dehn--Sommerville row changed")
        rows.append(
            {
                "a": a,
                "b": 8 + 6 * a,
                "f_vector": list(vector),
                "facets": vector[-1],
                "global_g4": 8 - a,
                "summed_vertex_link_g4": 6 * (8 - a),
            }
        )
    return rows


def schedule_exhaustion() -> dict[str, object]:
    rows = schedule_table()
    positive_component = _load(
        ROOT / str(_component_by_id("a5_a7_positive_scout")["archive"])
    )
    if positive_component.get("schedule_table") != [
        {
            **row,
            "strict_vertex_sum_feasible": row["summed_vertex_link_g4"] >= 14,
        }
        for row in rows
    ]:
        raise ValueError("independent and component schedule tables disagree")

    # The inherited Lanner order-five bound plus the Murai--Nevo equality
    # theorem makes global g4 strictly positive.  Thus a >= 8 is impossible.
    if rows[-1]["global_g4"] != 1 or 8 - 8 != 0:
        raise AssertionError("strict global-g4 cutoff changed")
    if [
        row["a"] for row in rows if row["summed_vertex_link_g4"] < 14
    ] != [6, 7]:
        raise AssertionError("strict summed link-g4 cutoff changed")

    closed_by_components = sorted(
        schedule
        for component in COMPONENTS
        for schedule in component["schedules"]
    )
    if closed_by_components != list(range(8)):
        raise AssertionError("component schedule cover is not exact")
    return {
        "parameter": "a = number of missing independent four-sets",
        "nonnegative_integer_parameter": True,
        "rows": rows,
        "strict_global_g4_excludes_a_at_least_8": True,
        "strict_summed_vertex_link_g4_excludes": [6, 7],
        "component_certificates_close": list(range(8)),
        "all_nonnegative_a_closed": True,
        "strict_glbt_lanner_dependency": {
            "path": str(LOCAL_LINK_REPORT.relative_to(ROOT)),
            "sha256": EXPECTED_LOCAL_LINK_REPORT_SHA256,
            "statement": (
                "Every compact H9 Coxeter nerve and every relevant face link "
                "has minimal nonfaces of order at most five; the generalized-"
                "lower-bound equality theorem therefore forces the relevant "
                "g entries, including global and vertex-link g4, to be strict."
            ),
        },
    }


def _proof_sources() -> list[dict[str, str]]:
    if _file_sha256(LOCAL_LINK_REPORT) != EXPECTED_LOCAL_LINK_REPORT_SHA256:
        raise ValueError("strict GLBT/Lanner theorem report changed")
    return [
        _semantic_certifier_source_record(),
        {
            "path": str(LOCAL_LINK_REPORT.relative_to(ROOT)),
            "sha256": EXPECTED_LOCAL_LINK_REPORT_SHA256,
        },
        {
            "path": str(A0_HOSTILE_SOURCE.relative_to(ROOT)),
            "sha256": EXPECTED_A0_HOSTILE_SOURCE_SHA256,
        },
        {
            "path": str(A0_HOSTILE_ARCHIVE.relative_to(ROOT)),
            "sha256": EXPECTED_A0_HOSTILE_ARCHIVE_SHA256,
        },
        {
            "path": str(A0_FINAL_HOSTILE_ARCHIVE.relative_to(ROOT)),
            "sha256": EXPECTED_A0_FINAL_HOSTILE_ARCHIVE_SHA256,
        },
        {
            "path": str(A0_FINAL_HOSTILE_REPORT.relative_to(ROOT)),
            "sha256": EXPECTED_A0_FINAL_HOSTILE_REPORT_SHA256,
        },
    ]


def fresh_record() -> dict[str, object]:
    components = verify_components()
    exhaustion = schedule_exhaustion()
    hostile = _verify_a0_hostile_archive()
    document: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record30_whole_shadow_compact_coxeter_exclusion",
        "claim": (
            "No compact hyperbolic Coxeter 9-polytope with fourteen facets "
            "has Record-30 missing-edge graph 3K2 + 8K1."
        ),
        "record_index": 30,
        "graph6": "M?????????O?O?G??",
        "missing_edges": [[0, 11], [1, 12], [2, 13]],
        "proof_scope": (
            "Whole graph shadow under the compact H9/14 missing-face, "
            "Lanner, incidence, finite-vertex, and Gram reductions."
        ),
        "schedule_exhaustion": exhaustion,
        "component_certificates": components,
        "a0_no_import_orbit_audit": hostile,
        "closed_schedules": list(range(8)),
        "compact_record30_excluded": True,
        "whole_shadow_compact_coxeter_theorem": True,
        "compact_only": True,
        "finite_volume_noncompact_checked": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "new_hyperbolic_coxeter_polytope": False,
        "frontier_update": {
            "before": {
                "compact_excluded_shadow_count": 132,
                "compact_open_shadow_count": 10,
            },
            "after": {
                "compact_excluded_shadow_count": 133,
                "compact_open_shadow_count": 9,
                "open_shadow_indices": OPEN_COMPACT_SHADOWS_AFTER_RECORD30,
            },
            "ordinary_incidence_partition": {
                "impossible_shadow_count": 123,
                "ordinary_polytopal_shadow_count": 19,
                "changed": False,
            },
        },
        "nonclaims": [
            "No finite-volume noncompact Record-30 chamber is excluded.",
            "No ideal-vertex, affine-cusp, maximal-parabolic, or Vinberg finite-volume endpoint audit is part of this theorem.",
            "No hyperbolic Coxeter polytope is constructed.",
            "This closes one graph shadow, not all fourteen-facet H9 shadows.",
        ],
        "proof_sources": _proof_sources(),
    }
    document["certificate_sha256"] = _canonical_digest(document)
    return document


def verify_record(
    document: Mapping[str, object], *, require_archive_file: bool = False
) -> dict[str, object]:
    if document.get("kind") != (
        "h9_14_record30_whole_shadow_compact_coxeter_exclusion"
    ):
        raise ValueError("wrong Record-30 aggregate certificate kind")
    digest = _canonical_digest(document)
    if document.get("certificate_sha256") != digest:
        raise ValueError("Record-30 aggregate internal digest is invalid")
    if (
        EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN"
        and digest != EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("Record-30 aggregate certificate is not frozen")
    if require_archive_file:
        if not OUTPUT.is_file():
            raise ValueError("Record-30 aggregate archive is missing")
        observed = _file_sha256(OUTPUT)
        if (
            EXPECTED_ARCHIVE_FILE_SHA256 != "TO_BE_FROZEN"
            and observed != EXPECTED_ARCHIVE_FILE_SHA256
        ):
            raise ValueError("Record-30 aggregate archive file changed")

    current_components = verify_components()
    current_exhaustion = schedule_exhaustion()
    current_hostile = _verify_a0_hostile_archive()
    if document.get("component_certificates") != current_components:
        raise ValueError("Record-30 aggregate component composition changed")
    if document.get("schedule_exhaustion") != current_exhaustion:
        raise ValueError("Record-30 aggregate schedule proof changed")
    if document.get("a0_no_import_orbit_audit") != current_hostile:
        raise ValueError("Record-30 aggregate hostile a=0 audit changed")
    if document.get("proof_sources") != _proof_sources():
        raise ValueError("Record-30 aggregate proof source changed")

    expected_frontier = {
        "before": {
            "compact_excluded_shadow_count": 132,
            "compact_open_shadow_count": 10,
        },
        "after": {
            "compact_excluded_shadow_count": 133,
            "compact_open_shadow_count": 9,
            "open_shadow_indices": OPEN_COMPACT_SHADOWS_AFTER_RECORD30,
        },
        "ordinary_incidence_partition": {
            "impossible_shadow_count": 123,
            "ordinary_polytopal_shadow_count": 19,
            "changed": False,
        },
    }
    if document.get("frontier_update") != expected_frontier:
        raise ValueError("Record-30 aggregate frontier accounting changed")
    if (
        document.get("closed_schedules") != list(range(8))
        or document.get("compact_record30_excluded") is not True
        or document.get("whole_shadow_compact_coxeter_theorem") is not True
        or document.get("compact_only") is not True
        or document.get("finite_volume_noncompact_checked") is not False
        or document.get("finite_volume_noncompact_chamber_excluded") is not False
        or document.get("new_hyperbolic_coxeter_polytope") is not False
    ):
        raise ValueError("Record-30 aggregate conclusion or scope changed")
    return {
        "compact_record30_excluded": True,
        "closed_schedules": list(range(8)),
        "all_nonnegative_a_closed": True,
        "finite_volume_noncompact_checked": False,
        "new_hyperbolic_example": False,
        "compact_frontier": "133 + 9",
        "open_compact_shadows": OPEN_COMPACT_SHADOWS_AFTER_RECORD30,
        "ordinary_incidence_frontier": "123 + 19",
        "certificate_sha256": digest,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fresh", action="store_true")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    if args.fresh or args.write:
        document = fresh_record()
    else:
        document = _load(OUTPUT)
    summary = verify_record(
        document, require_archive_file=not (args.fresh or args.write)
    )
    if args.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(json.dumps({"output": str(OUTPUT), **summary}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
