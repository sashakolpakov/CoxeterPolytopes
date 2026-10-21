"""Apply solver-independent refinements to the bounded H9/14 facet run."""

from __future__ import annotations

from collections import Counter
import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.facet_completion import (
    audit_global_facet_selection,
    explicit_product_facet_selection,
    forced_small_transversal_obstruction,
)
from coxeter_polytopes.combinatorics.nerve_support import verify_support_manifest
from h9_14_facet_completion import verify as verify_baseline


ROOT = Path(__file__).resolve().parents[1]
SUPPORT = ROOT / "data" / "enumerations" / "h9_14_nerve_support.json"
BASELINE = ROOT / "data" / "enumerations" / "h9_14_facet_completion.json"
OUTPUT = ROOT / "data" / "enumerations" / "h9_14_facet_refinement.json"
EXPECTED_BASELINE_SHA256 = (
    "2d329439e5baeffb574932ff7bbab2dd78daf706b2917a0674bdc9f01a2f5123"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "c976ac262bb78dc4020210e3eab8ed662d5719de0d058f1a94089c60c7f4ae25"
)
EXPECTED_COUNTS = {
    "input": 142,
    "certified_unsat": 113,
    "certified_explicit_witness": 2,
    "theorem_undecided": 27,
    "bounded_timeout_unknown": 25,
    "backend_unsat_pending_independent_certificate": 2,
}


def _sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _critical_digest(records: Sequence[Mapping[str, object]]) -> str:
    payload = [
        {
            "index": record["index"],
            "graph6": record["graph6"],
            "disposition": record["disposition"],
            "theorem_status": record["theorem_status"],
            "selected_facet_complement_masks": record.get(
                "selected_facet_complement_masks", []
            ),
            "small_transversal_certificate": record.get(
                "small_transversal_certificate"
            ),
        }
        for record in records
    ]
    return sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _counts(records: Sequence[Mapping[str, object]]) -> dict[str, object]:
    dispositions = Counter(record["disposition"] for record in records)
    theorem_statuses = Counter(record["theorem_status"] for record in records)
    return {
        "input": len(records),
        "certified_unsat": theorem_statuses["certified_unsat"],
        "certified_explicit_witness": theorem_statuses[
            "certified_explicit_witness"
        ],
        "theorem_undecided": theorem_statuses["undecided"],
        "bounded_timeout_unknown": dispositions["bounded_z3_unknown"],
        "backend_unsat_pending_independent_certificate": dispositions[
            "z3_unsat_not_independently_certified"
        ],
    }


def build() -> dict[str, object]:
    verify_baseline(
        SUPPORT, BASELINE, rerun_unsat=False, timeout_ms=None
    )
    support = json.loads(SUPPORT.read_text(encoding="utf-8"))
    verify_support_manifest(support)
    baseline = json.loads(BASELINE.read_text(encoding="utf-8"))
    source_records = [
        record for record in support["records"] if record["survives_second_layer"]
    ]
    records = []
    for index, (source, old) in enumerate(zip(source_records, baseline["records"])):
        graph6 = source["graph6"]
        witness = explicit_product_facet_selection(graph6)
        transversal = forced_small_transversal_obstruction(
            graph6, source["dehn_sommerville_schedule"]
        )
        old_status = old["solver"]["status"]
        old_backend = old["solver"]["backend"]
        record: dict[str, object] = {
            "index": index,
            "graph6": graph6,
            "edges": source["edges"],
            "selected_facet_complement_masks": [],
            "small_transversal_certificate": None,
            "witness_audit": None,
        }
        if old_backend == "exact_local_link_precheck":
            record.update(
                status="unsat",
                theorem_status="certified_unsat",
                disposition="independent_local_link_unsat",
            )
        elif transversal is not None:
            record.update(
                status="unsat",
                theorem_status="certified_unsat",
                disposition="independent_small_transversal_unsat",
                small_transversal_certificate=transversal,
            )
        elif witness is not None:
            audit = audit_global_facet_selection(
                graph6,
                source["ridge_support_core_masks"],
                source["dehn_sommerville_schedule"],
                witness,
            )
            if not audit.passes_connected_homology_audit:
                raise AssertionError("explicit product witness failed exact audit")
            record.update(
                status="sat",
                theorem_status="certified_explicit_witness",
                disposition="explicit_product_witness",
                selected_facet_complement_masks=list(witness),
                witness_audit=audit.to_record(),
            )
        elif old_status == "unsat":
            record.update(
                status="unsat",
                theorem_status="undecided",
                disposition="z3_unsat_not_independently_certified",
            )
        else:
            record.update(
                status="unknown",
                theorem_status="undecided",
                disposition="bounded_z3_unknown",
            )
        records.append(record)

    manifest = {
        "schema_version": "1.0.0",
        "object": "solver-independent refinement of bounded global facet selection",
        "source": {
            "support_path": "data/enumerations/h9_14_nerve_support.json",
            "support_sha256": _sha256(SUPPORT),
            "baseline_path": "data/enumerations/h9_14_facet_completion.json",
            "baseline_sha256": _sha256(BASELINE),
        },
        "scope": [
            "Local-link rejections are recomputable exact certificates.",
            "The small-transversal certificate is a human-checkable exact contradiction.",
            "Explicit product witnesses pass full face, link-Euler, and mod-two homology audits.",
            "The two remaining Z3 UNSAT results are not promoted to independent certificates.",
            "Unknown means undecided, never rejected.",
        ],
        "records": records,
        "counts": _counts(records),
    }
    manifest["certificate_sha256"] = _critical_digest(records)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    rebuilt = build()
    if manifest.get("source") != rebuilt["source"]:
        raise ValueError("facet refinement source provenance mismatch")
    if manifest.get("records") != rebuilt["records"]:
        raise ValueError("facet refinement records do not verify")
    if manifest.get("counts") != rebuilt["counts"]:
        raise ValueError("facet refinement counts do not verify")
    if manifest.get("certificate_sha256") != rebuilt["certificate_sha256"]:
        raise ValueError("facet refinement digest does not verify")
    if manifest.get("source", {}).get("baseline_sha256") != (
        EXPECTED_BASELINE_SHA256
    ):
        raise ValueError("facet refinement baseline is not frozen")
    if manifest.get("certificate_sha256") != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("facet refinement certificate is not frozen")
    if manifest.get("counts") != EXPECTED_COUNTS:
        raise ValueError("facet refinement frozen count vector mismatch")
    return rebuilt["counts"]


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args(argv)
    if args.verify_only:
        manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    else:
        manifest = build()
        OUTPUT.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    counts = verify(manifest)
    print(json.dumps({"output": str(OUTPUT), "counts": counts}, indent=2))


if __name__ == "__main__":
    main()
