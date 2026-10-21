#!/usr/bin/env python3
"""Verify the no-project-import hostile replay for H9/14 record 32.

The default mode validates the frozen independent summary and all source
digests.  ``--positive`` recompiles and reruns the complete a=1,2,3 positive
A4 orbit census (about twenty seconds on the audit machine).  ``--full`` also
runs the 179,785,655-path independent a=0 stream and the exhaustive full-
action orbit comparator (a few minutes).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ARCHIVE = (
    ROOT
    / "data"
    / "enumerations"
    / "h9_14_record32_whole_shadow_hostile_audit.json"
)

EXPECTED_BRANCHES = [
    (1, 602, 602, 322, 4, {"28": 1, "70": 1, "112": 2}),
    (2, 180901, 172039, 2156, 6,
     {"56": 2, "112": 1, "420": 1, "672": 1, "840": 1}),
    (3, 36180200, 30954448, 56, 1, {"56": 1}),
]

EXPECTED_A0 = {
    "active_mask_count": 24,
    "active_after_caps": 22513,
    "active_after_cover": 14731,
    "column_paths": 179785655,
    "repeated_row_paths": 3513588,
    "literal_small_transversal_paths": 0,
    "blocker_valid_paths": 176272067,
    "target_paths": 1263,
    "target_family_count": 969,
}

EXPECTED_ORBIT_COMPARISON = {
    "production_literal_representative_count": 1025,
    "independent_literal_representative_count": 969,
    "production_orbit_count": 292,
    "independent_orbit_count": 292,
    "production_only_orbit_count": 0,
    "independent_only_orbit_count": 0,
    "target_orbit_sets_equal": True,
    "orbit_key_sha256": (
        "97f630bd464ec01d77087019c5dd02b549f0cbaebd466d5e621ae065797dc848"
    ),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_sha256(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def validate_archive(archive: dict) -> None:
    require(
        archive.get("kind") == "h9_14_record32_whole_shadow_hostile_audit",
        "unexpected archive kind",
    )
    status = archive["status"]
    require(status["independent_replay_complete"], "replay is not complete")
    require(
        status["whole_graph_shadow_compact_coxeter_no_go"],
        "compact whole-shadow verdict is absent",
    )
    scope = archive["scope"]
    require(scope["compact_only"], "scope lost compact-only qualifier")
    require(
        scope["abstract_targets_are_only_necessary_incidence_objects"],
        "target semantics were broadened",
    )
    require(
        not scope["ordinary_polytopality_of_abstract_targets_claimed"]
        and not scope["PL_sphere_status_of_abstract_targets_claimed"]
        and not scope["finite_volume_noncompact_chamber_excluded"],
        "archive makes an unsupported topology/noncompact claim",
    )

    for source in archive["independent_sources"].values():
        path = ROOT / source["path"]
        require(path.is_file(), f"missing hostile source {path}")
        require(sha256(path) == source["sha256"], f"source digest mismatch: {path}")

    positive = archive["positive_shadow_gate"]
    require(positive["independent_A4_count"] == 762, "wrong independent A4 count")
    require(positive["interval_A4_count"] == 602, "wrong interval A4 count")
    require(positive["all_local_lower_endpoints_zero"], "lower endpoint changed")
    require(len(positive["branches"]) == 3, "wrong positive-branch count")
    for branch, expected in zip(positive["branches"], EXPECTED_BRANCHES):
        summary = (
            branch["a"],
            branch["raw"],
            branch["interval_feasible"],
            branch["strict_aggregate_viable"],
            branch["orbit_count"],
            branch["orbit_size_histogram"],
        )
        require(summary == expected, f"positive branch a={branch['a']} changed")
        require(
            sum(row["count"] for row in branch["gate_histogram"])
            == branch["interval_feasible"],
            f"gate histogram does not exhaust a={branch['a']}",
        )
        require(
            sum(
                row["count"]
                for row in branch["gate_histogram"]
                if row["strict"] and row["aggregate"]
            )
            == branch["strict_aggregate_viable"],
            f"viable histogram mismatch for a={branch['a']}",
        )
        require(
            sum(int(size) * count
                for size, count in branch["orbit_size_histogram"].items())
            == branch["strict_aggregate_viable"],
            f"orbit sizes do not exhaust a={branch['a']}",
        )
        require(
            len(branch["orbit_keys_and_sizes"]) == branch["orbit_count"]
            and len({tuple(row["key"]) for row in branch["orbit_keys_and_sizes"]})
            == branch["orbit_count"]
            and sum(row["size"] for row in branch["orbit_keys_and_sizes"])
            == branch["strict_aggregate_viable"],
            f"orbit keys do not exhaust a={branch['a']}",
        )

    a0 = archive["a0"]
    require(a0["independent_replay"] == EXPECTED_A0, "a0 replay transcript changed")
    require(
        a0["full_action_orbit_comparison"] == EXPECTED_ORBIT_COMPARISON,
        "a0 orbit comparison changed",
    )
    target = a0["production_target_archive"]
    target_path = ROOT / target["path"]
    require(target_path.is_file(), "missing production target archive")
    require(sha256(target_path) == target["sha256"], "target archive digest mismatch")
    require(
        EXPECTED_A0["column_paths"] - EXPECTED_A0["repeated_row_paths"]
        == EXPECTED_A0["blocker_valid_paths"],
        "a0 path partition changed",
    )
    require(
        a0["production_literal_families"] == 1025
        and a0["pure_blocker_families"] == 941,
        "a0 literal/blocker census changed",
    )
    require(
        sum(a0["impure_blocker_histogram"].values()) == 84,
        "a0 impure-blocker histogram changed",
    )
    require(
        len(a0["ridge_valid_indices"]) == 12
        and a0["ridge_histogram"] == {"2": 513}
        and a0["ridge_valid_unrestricted_terminal_state_counts"] == [0] * 12,
        "a0 ridge/label terminal changed",
    )

    a1 = archive["a1"]
    require(len(a1["fixed_A4_transcripts"]) == 4, "wrong a1 orbit count")
    require(a1["ridge_valid_target_count"] == 2, "wrong a1 ridge count")
    require(a1["ridge_histogram"] == {"2": 504}, "wrong a1 ridge histogram")
    require(len(a1["ridge_valid_targets"]) == 2, "wrong derived a1 target identities")
    for gate in a1["cross_gates"]:
        require(
            gate["single_m3_tests"] == gate["sparse_states"] * gate["cross_pairs"]
            and gate["passing"] == 0,
            "a1 unrestricted cross gate changed",
        )

    a2 = archive["a2"]
    require(
        len(a2["representative_transcripts"]) == 6
        and sum(row["targets"] for row in a2["representative_transcripts"]) == 1,
        "a2 representative census changed",
    )
    require(
        a2["sole_target_is_older_selected_record32_nerve"]
        and a2["sole_target_uses_older_unrestricted_compact_label_no_go"],
        "a2 selected-nerve dependency changed",
    )
    selected_certificate = a2["older_selected_nerve_certificate"]
    selected_path = ROOT / selected_certificate["path"]
    require(
        selected_path.is_file()
        and sha256(selected_path) == selected_certificate["file_sha256"],
        "selected-a2 certificate file changed",
    )
    a3 = archive["a3"]
    require(
        a3["strict_shadow_family_count"] == 56
        and a3["strict_shadow_orbit_count"] == 1
        and a3["raw"] == 1540
        and a3["cap_feasible"] == 540
        and a3["no_small"] == 0
        and a3["cap_zero_one_two_generator_lossless"],
        "a3 lossless gate changed",
    )


def compile_source(cxx: str, source: Path, output: Path, pthread: bool = False) -> None:
    command = [cxx, "-std=c++20", "-O3"]
    if pthread:
        command.append("-pthread")
    command.extend([str(source), "-o", str(output)])
    subprocess.run(command, cwd=ROOT, check=True)


def parse_json_output(command: list[str]) -> dict:
    result = subprocess.run(
        command, cwd=ROOT, check=True, text=True, stdout=subprocess.PIPE
    )
    return json.loads(result.stdout)


def replay_positive(archive: dict, cxx: str, temporary: Path) -> dict:
    source = ROOT / archive["independent_sources"]["positive_orbit_replay"]["path"]
    executable = temporary / "positive_orbits"
    compile_source(cxx, source, executable)
    actual = parse_json_output([str(executable)])
    expected = archive["positive_shadow_gate"]
    require(
        actual["independent_A4_count"] == expected["independent_A4_count"]
        and actual["interval_A4_count"] == expected["interval_A4_count"]
        and actual["branches"] == expected["branches"],
        "fresh positive-orbit replay disagrees with archive",
    )
    return actual


def replay_additions(archive: dict, positive: dict) -> dict:
    source = archive["independent_sources"]["additions_replay"]
    actual = parse_json_output(
        [shutil.which("python3") or "python3", str(ROOT / source["path"])]
    )
    require(
        canonical_sha256(actual) == source["canonical_output_sha256"],
        "fresh additions replay disagrees with frozen full transcript",
    )
    branch_names = {1: "a1", 2: "a2", 3: "a3"}
    for branch in positive["branches"]:
        name = branch_names[branch["a"]]
        expected = {
            tuple(row["key"]): row["size"] for row in branch["orbit_keys_and_sizes"]
        }
        observed = dict(zip(
            map(tuple, actual["positive_representative_orbit_keys"][name]),
            actual["positive_representative_orbit_sizes"][name],
            strict=True,
        ))
        require(
            observed == expected,
            f"fresh {name} representative keys do not equal the complete positive orbits",
        )

    archived_a1 = archive["a1"]
    for actual_row, frozen in zip(
        actual["a1"]["representatives"], archived_a1["fixed_A4_transcripts"], strict=True
    ):
        require(actual_row["A4"] == frozen["A4"], "a1 representative changed")
        for key in (
            "active_configurations", "column_paths", "distinct_row_paths",
            "no_small_paths", "literal_targets", "target_orbits",
        ):
            require(actual_row[key] == frozen[key], f"a1 {key} transcript changed")
    observed_ridge = [
        {"A4": row["A4"], "A5s": row["additions"]}
        for row in actual["a1"]["ridge_targets"]
    ]
    require(
        observed_ridge == archived_a1["ridge_valid_targets"],
        "fresh a1 target derivation changed",
    )

    for actual_row, frozen in zip(
        actual["a2"]["representatives"], archive["a2"]["representative_transcripts"],
        strict=True,
    ):
        projection = {
            "shadow": actual_row["shadow"],
            "additional": actual_row["additional"],
            "candidates": actual_row["candidates"],
            "raw": actual_row["raw"],
            "cap_feasible": actual_row["cap_feasible"],
            "no_small": actual_row["no_small"],
            "targets": len(actual_row["targets"]),
        }
        require(projection == frozen, "fresh a2 representative transcript changed")
    require(
        actual["a2"]["selected_target"] == archive["a2"]["sole_target"],
        "fresh a2 target is not the selected record-32 nerve",
    )
    actual_a3 = actual["a3"]
    frozen_a3 = archive["a3"]
    require(
        (
            actual_a3["shadow"], actual_a3["additional"], actual_a3["candidates"],
            actual_a3["raw"], actual_a3["cap_feasible"], actual_a3["no_small"],
            len(actual_a3["targets"]),
        )
        == (
            frozen_a3["representative_shadow"], frozen_a3["representative_additional"],
            frozen_a3["candidates"], frozen_a3["raw"], frozen_a3["cap_feasible"],
            frozen_a3["no_small"], frozen_a3["targets"],
        ),
        "fresh a3 addition transcript changed",
    )
    return actual


def replay_downstream(archive: dict) -> dict:
    source = archive["independent_sources"]["downstream_replay"]
    actual = parse_json_output(
        [shutil.which("python3") or "python3", str(ROOT / source["path"])]
    )
    require(
        canonical_sha256(actual) == source["canonical_output_sha256"],
        "fresh blocker/Lanner/cross replay disagrees with frozen full transcript",
    )
    archived_a0 = archive["a0"]
    actual_a0 = actual["a0"]
    require(
        actual_a0["target_family_count"] == archived_a0["production_literal_families"]
        and actual_a0["pure_blocker_count"] == archived_a0["pure_blocker_families"]
        and actual_a0["impure_blocker_histogram"] == archived_a0["impure_blocker_histogram"]
        and actual_a0["ridge_valid_indices"] == archived_a0["ridge_valid_indices"]
        and actual_a0["ridge_histogram"] == archived_a0["ridge_histogram"]
        and [row["terminal_state_count"] for row in actual_a0["label_audits"]]
        == archived_a0["ridge_valid_unrestricted_terminal_state_counts"],
        "fresh a0 blocker/ridge/Lanner replay changed",
    )
    for row, frozen in zip(
        actual["a1"]["target_audits"], archive["a1"]["cross_gates"], strict=True
    ):
        projection = {
            "sparse_states": row["sparse_state_count"],
            "cross_pairs": row["cross_pair_count"],
            "single_m3_tests": row["single_m3_test_count"],
            "passing": row["single_m3_passing_count"],
            "first_bad_restriction_size_histogram": row[
                "first_bad_restriction_size_histogram"
            ],
        }
        require(projection == frozen, "fresh a1 cross gate changed")
    dependency = actual["selected_a2_dependency"]
    frozen_dependency = archive["a2"]["older_selected_nerve_certificate"]
    require(
        dependency["certificate_path"] == frozen_dependency["path"]
        and dependency["file_sha256"] == frozen_dependency["file_sha256"]
        and dependency["certificate_sha256"] == frozen_dependency["certificate_sha256"]
        and dependency["order_four_minimal_faces"]
        == archive["a2"]["sole_target"]["order_four_minimal_faces"]
        and dependency["order_five_minimal_face"]
        == archive["a2"]["sole_target"]["order_five_minimal_face"]
        and dependency["unrestricted_compact_coxeter_no_go"],
        "fresh selected-a2 dependency identity changed",
    )
    return actual


def replay_a0(archive: dict, cxx: str, temporary: Path) -> None:
    sources = archive["independent_sources"]
    source = ROOT / sources["a0_replay"]["path"]
    executable = temporary / "a0_replay"
    targets = temporary / "a0_targets.txt"
    compile_source(cxx, source, executable, pthread=True)
    actual = parse_json_output([str(executable), "--target-output", str(targets)])
    actual.pop("kind", None)
    require(actual == archive["a0"]["independent_replay"], "fresh a0 replay changed")

    comparator = ROOT / sources["a0_orbit_comparator"]["path"]
    production = ROOT / archive["a0"]["production_target_archive"]["path"]
    comparison = parse_json_output(
        [shutil.which("python3") or "python3", str(comparator), str(production), str(targets)]
    )
    require(
        comparison == archive["a0"]["full_action_orbit_comparison"],
        "fresh a0 full-action orbit comparison changed",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", nargs="?", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument(
        "--positive", action="store_true",
        help="compile and run the independent positive-orbit replay",
    )
    parser.add_argument(
        "--full", action="store_true",
        help="run every hostile replay, including the full independent a0 generator",
    )
    parser.add_argument(
        "--downstream", action="store_true",
        help="run positive, additions, and blocker/Lanner/cross replays (skip a0 generator)",
    )
    parser.add_argument("--cxx", help="C++ compiler (default: clang++ or c++)")
    args = parser.parse_args()

    archive_path = args.archive.resolve()
    archive = json.loads(archive_path.read_text(encoding="utf-8"))
    validate_archive(archive)
    cxx = args.cxx or shutil.which("clang++") or shutil.which("c++")
    if (args.positive or args.downstream or args.full) and not cxx:
        raise RuntimeError("no C++ compiler found")

    with tempfile.TemporaryDirectory(prefix="record32-hostile-") as name:
        temporary = Path(name)
        positive = None
        if args.positive or args.downstream or args.full:
            positive = replay_positive(archive, cxx, temporary)
        if args.downstream or args.full:
            replay_additions(archive, positive)
            replay_downstream(archive)
        if args.full:
            replay_a0(archive, cxx, temporary)

    output = {
        "kind": "h9_14_record32_whole_shadow_hostile_certificate",
        "archive_sha256": sha256(archive_path),
        "source_digests_verified": True,
        "semantic_scope_verified": True,
        "positive_replay_ran": bool(args.positive or args.downstream or args.full),
        "additions_replay_ran": bool(args.downstream or args.full),
        "blocker_lanner_cross_replay_ran": bool(args.downstream or args.full),
        "a0_full_replay_ran": bool(args.full),
        "whole_graph_shadow_compact_coxeter_no_go": True,
        "finite_volume_noncompact_inference": False,
        "status": "pass",
    }
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
