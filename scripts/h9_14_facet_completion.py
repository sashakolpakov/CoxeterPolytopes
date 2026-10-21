"""Solve or verify the exact H9/14 global facet-selection archive."""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.facet_completion import (
    audit_global_facet_selection,
    solve_facet_selection,
)
from coxeter_polytopes.combinatorics.local_link import analyze_local_link_graph
from coxeter_polytopes.combinatorics.nerve_support import verify_support_manifest


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT / "data" / "enumerations" / "h9_14_nerve_support.json"
DEFAULT_OUTPUT = ROOT / "data" / "enumerations" / "h9_14_facet_completion.json"
EXPECTED_SOURCE_SHA256 = (
    "8121f869c182b1780f3d812d34823f1579a5c37b085a4a8ecc313ea94937c47d"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "2789e320b0b23122d53fcf27b68d9198d575af64377a485dde9b973af77a8f0c"
)
EXPECTED_BASELINE_COUNTS = {
    "input_graph_shadows": 142,
    "sat": 0,
    "unsat": 115,
    "unknown": 27,
}


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _repository_path(path: Path) -> str:
    """Return stable repository-relative provenance when possible."""

    resolved = path.resolve()
    try:
        return resolved.relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return str(resolved)


def _problem_record(source_record: Mapping[str, object]) -> dict[str, object]:
    graph6 = source_record.get("graph6")
    core = source_record.get("ridge_support_core_masks")
    schedules = source_record.get("dehn_sommerville_schedule")
    if not isinstance(graph6, str):
        raise ValueError("source record lacks graph6")
    if not isinstance(core, list) or not all(isinstance(mask, int) for mask in core):
        raise ValueError("source record has invalid ridge-support core")
    if not isinstance(schedules, list) or not all(
        isinstance(schedule, Mapping) for schedule in schedules
    ):
        raise ValueError("source record has invalid DS schedules")
    return {
        "graph6": graph6,
        "edges": source_record.get("edges"),
        "ridge_support_core_count": len(core),
        "ridge_support_core_masks": core,
        "dehn_sommerville_schedule": schedules,
    }


def _certificate_digest(records: Sequence[Mapping[str, object]]) -> str:
    critical = [
        {
            "graph6": record["graph6"],
            "status": record["solver"]["status"],
            "selected_facet_complement_masks": record["solver"][
                "selected_facet_complement_masks"
            ],
        }
        for record in records
    ]
    payload = json.dumps(
        critical, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return sha256(payload).hexdigest()


def _summary(records: Sequence[Mapping[str, object]]) -> dict[str, object]:
    statuses = Counter(record["solver"]["status"] for record in records)
    by_edges: dict[str, dict[str, int]] = {}
    for edge_count in range(2, 10):
        subset = [record for record in records if record["edges"] == edge_count]
        edge_statuses = Counter(record["solver"]["status"] for record in subset)
        by_edges[str(edge_count)] = {
            "input": len(subset),
            "sat": edge_statuses["sat"],
            "unsat": edge_statuses["unsat"],
            "unknown": edge_statuses["unknown"],
        }
    return {
        "input_graph_shadows": len(records),
        "sat": statuses["sat"],
        "unsat": statuses["unsat"],
        "unknown": statuses["unknown"],
        "by_edges": by_edges,
    }


def generate(
    source_path: Path,
    output_path: Path,
    *,
    timeout_ms: int | None,
) -> dict[str, object]:
    source = json.loads(source_path.read_text(encoding="utf-8"))
    verify_support_manifest(source)
    records_by_index: dict[int, dict[str, object]] = {}
    survivors = [
        record for record in source["records"] if record["survives_second_layer"]
    ]
    # Exact local-link rejections are essentially free.  Among the remaining
    # 78 records, solving the smallest support cores first gives useful exact
    # progress even when a bounded run leaves larger instances undecided.
    work_order = sorted(
        enumerate(survivors),
        key=lambda item: (
            analyze_local_link_graph(item[1]["graph6"]).passes_edge_and_triple_feasibility,
            item[1]["ridge_support_core_count"],
            item[0],
        ),
    )
    for completed, (index, source_record) in enumerate(work_order, start=1):
        problem = _problem_record(source_record)
        result = solve_facet_selection(
            problem["graph6"],
            problem["ridge_support_core_masks"],
            problem["dehn_sommerville_schedule"],
            timeout_ms=timeout_ms,
        )
        record = {
            "index": index,
            "graph6": problem["graph6"],
            "edges": problem["edges"],
            "ridge_support_core_count": problem["ridge_support_core_count"],
            "solver": result.to_record(),
        }
        if result.status == "sat":
            audit = audit_global_facet_selection(
                problem["graph6"],
                problem["ridge_support_core_masks"],
                problem["dehn_sommerville_schedule"],
                result.selected_facet_complement_masks,
            )
            if not audit.passes_selection_model:
                raise AssertionError("search backend returned an invalid witness")
            record["witness_audit"] = audit.to_record()
        else:
            record["witness_audit"] = None
        records_by_index[index] = record
        print(
            f"[{completed}/{len(survivors)} index={index}] e={problem['edges']} "
            f"core={problem['ridge_support_core_count']} "
            f"status={result.status} seconds={result.elapsed_seconds:.3f}",
            flush=True,
        )

    records = [records_by_index[index] for index in range(len(survivors))]

    manifest = {
        "schema_version": "1.0.0",
        "object": "global facet-selection decisions, not polytopality certificates",
        "constraint_scope": [
            "selected five-set complements lie in the ridge-support core",
            "every six-set contains zero or two selected complements",
            "every independent triple is contained in a selected facet",
            "all face indicators are exactly the OR of containing selected facets",
            "one complete Dehn--Sommerville f-vector schedule holds",
            "every minimal nonface has order at most five",
            "all vertex, edge, and triangle figures are nonsimplex and obey the strict local g inequalities",
            "all enforced local g-vectors obey their exact Macaulay upper bounds",
            "missing four-face incidences obey every vertex, edge, and triangle interval",
        ],
        "nonclaims": [
            "A Z3 UNSAT disposition is exact but is not yet an independently checked proof certificate.",
            "Local-link precheck rejections are independently verified by the local-link archive.",
            "A SAT witness is not automatically a homology sphere or PL sphere.",
            "No oriented-matroid or polytopality certificate is supplied.",
            "No Coxeter labeling or Lorentzian Gram matrix is supplied.",
        ],
        "source": {
            "path": _repository_path(source_path),
            "sha256": _file_sha256(source_path),
            "schema_version": source["schema_version"],
            "second_layer_survivors": len(survivors),
        },
        "records": records,
    }
    manifest["counts"] = _summary(records)
    manifest["certificate_sha256"] = _certificate_digest(records)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def verify(
    source_path: Path,
    output_path: Path,
    *,
    rerun_unsat: bool,
    timeout_ms: int | None,
) -> dict[str, object]:
    source = json.loads(source_path.read_text(encoding="utf-8"))
    verify_support_manifest(source)
    manifest = json.loads(output_path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != "1.0.0":
        raise ValueError("unsupported facet-completion schema")
    if manifest.get("source", {}).get("path") != (
        "data/enumerations/h9_14_nerve_support.json"
    ):
        raise ValueError("facet-completion source path is not repository-relative")
    if manifest.get("source", {}).get("sha256") != _file_sha256(source_path):
        raise ValueError("facet-completion source hash mismatch")
    if _file_sha256(source_path) != EXPECTED_SOURCE_SHA256:
        raise ValueError("facet-completion source is not the frozen support archive")
    survivors = [
        record for record in source["records"] if record["survives_second_layer"]
    ]
    records = manifest.get("records")
    if not isinstance(records, list) or len(records) != len(survivors):
        raise ValueError("facet-completion record coverage mismatch")
    for index, (source_record, record) in enumerate(zip(survivors, records)):
        problem = _problem_record(source_record)
        if record.get("index") != index or record.get("graph6") != problem["graph6"]:
            raise ValueError(f"record {index} does not match the source ordering")
        solver_record = record.get("solver")
        if not isinstance(solver_record, Mapping):
            raise ValueError(f"record {index} lacks a solver disposition")
        status = solver_record.get("status")
        if status not in {"sat", "unsat", "unknown"}:
            raise ValueError(f"record {index} has invalid solver status")
        raw_selected = solver_record.get("selected_facet_complement_masks")
        if not isinstance(raw_selected, list) or not all(
            isinstance(mask, int) for mask in raw_selected
        ):
            raise ValueError(f"record {index} has invalid selected facets")
        backend = solver_record.get("backend")
        if backend == "exact_local_link_precheck":
            local = analyze_local_link_graph(problem["graph6"])
            if status != "unsat" or local.passes_edge_and_triple_feasibility:
                raise ValueError(
                    f"record {index} has an invalid local-link rejection"
                )
            if solver_record.get("assertion_count") != 0:
                raise ValueError(f"record {index} precheck carries Z3 assertions")
        elif backend != "z3":
            raise ValueError(f"record {index} has an unknown solver backend")
        if status == "sat":
            audit = audit_global_facet_selection(
                problem["graph6"],
                problem["ridge_support_core_masks"],
                problem["dehn_sommerville_schedule"],
                raw_selected,
            )
            if not audit.passes_selection_model:
                raise ValueError(f"record {index} SAT witness fails exact audit")
            if record.get("witness_audit") != audit.to_record():
                raise ValueError(f"record {index} witness audit mismatch")
        elif raw_selected:
            raise ValueError(f"record {index} non-SAT result carries a witness")
        if status == "unknown" and not solver_record.get("reason_unknown"):
            raise ValueError(f"record {index} unknown lacks a backend reason")
        if rerun_unsat and status == "unsat" and backend == "z3":
            rerun = solve_facet_selection(
                problem["graph6"],
                problem["ridge_support_core_masks"],
                problem["dehn_sommerville_schedule"],
                timeout_ms=timeout_ms,
            )
            if rerun.status != "unsat":
                raise ValueError(
                    f"record {index} UNSAT disposition did not reproduce: {rerun.status}"
                )
    if manifest.get("counts") != _summary(records):
        raise ValueError("facet-completion summary mismatch")
    if manifest.get("certificate_sha256") != _certificate_digest(records):
        raise ValueError("facet-completion certificate digest mismatch")
    if manifest.get("certificate_sha256") != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("facet-completion certificate is not the frozen baseline")
    counts = manifest["counts"]
    if any(counts.get(key) != value for key, value in EXPECTED_BASELINE_COUNTS.items()):
        raise ValueError("facet-completion baseline count vector mismatch")
    return manifest["counts"]


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--timeout-ms",
        type=int,
        default=None,
        help="preserve undecided if Z3 exceeds this per-record limit",
    )
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument(
        "--rerun-unsat",
        action="store_true",
        help="during verification, re-solve every archived UNSAT instance",
    )
    args = parser.parse_args(argv)
    if args.verify_only:
        counts = verify(
            args.source,
            args.output,
            rerun_unsat=args.rerun_unsat,
            timeout_ms=args.timeout_ms,
        )
        mode = "verified"
    else:
        manifest = generate(
            args.source,
            args.output,
            timeout_ms=args.timeout_ms,
        )
        counts = manifest["counts"]
        mode = "generated"
    print(json.dumps({"mode": mode, "counts": counts}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
