#!/usr/bin/env python3
"""Audit the CoxIter header denominator used by the positive-example ledger.

This is intentionally a *negative-control* certificate, not a novelty
classifier.  It scans the two pinned CoxIter lanes and records the number of
files with each ``(hyperbolic dimension, facet count)`` header occurring in
the maintained positive ledger.  Several chambers whose prior lineage is
already closed by exact, independent identifications have no such header in
either lane.  Therefore zero header matches can only be reported as bounded
corpus evidence; they cannot imply worldwide priority.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import tarfile
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "results/higher_dimensional_polytope_lineage_ledger.md"
OUTPUT = ROOT / "data/enumerations/positive_example_lineage_header_controls.json"
COXITER_DIRECTORY = Path("/private/tmp/coxiter-current/CoxIter-master")
COXITER_ARCHIVE = Path("/private/tmp/coxiter-current.zip")
DOUBLING_RELATIVE = Path("graphs/testing_index2.tar.gz")
X2_LINEAGE = (
    ROOT / "data/enumerations/shimada_aut3singk3_x2_h10_37_global_lineage.json"
)

SCHEMA = "coxeter.positive_example_lineage.header_negative_controls.v1"
COXITER_COMMIT = "b800be48d8240fbabd2d1c9542e61054584ad7fe"
COXITER_ARCHIVE_SHA256 = (
    "ae04648b291c1e5f3d0e47d666c846dc805b244469a20486bca59674b4d171c3"
)
DOUBLING_ARCHIVE_SHA256 = (
    "3330ca41fb3e047123ace718f08542bc1dd0019b35922d4fa836a217647e7036"
)
X2_LINEAGE_SHA256 = (
    "52b4800b07608804b7834aad65a4671b3255ef918c07779497742f11091a984d"
)

# Each control is independently identified with a source predating its
# insertion in the present ledger.  This script checks only the orthogonal
# fact that its header is absent from both pinned CoxIter lanes.
KNOWN_ZERO_HEADER_CONTROLS = {
    "R11-D4A1x5-NEF-H10-90": {
        "prior_source": "Roulleau, An atlas of K3 surfaces with finite automorphism group",
        "source_url": "https://arxiv.org/abs/2003.08985",
        "exact_local_lineage": "results/roulleau_rank11_direct_graph_denominator.md",
    },
    "K3-D4-H12": {
        "prior_source": "Scharlau--Kirschmer reflective Lorentzian census, List12 entry 6",
        "source_url": "https://www.math.rwth-aachen.de/~Markus.Kirschmer/lorentz/index.html",
        "exact_local_lineage": "results/k3_rank15_order3_h12_scharlau_kirschmer_lineage.md",
    },
    "SS2-2C-H14-98": {
        "prior_source": "Alexeev, Reflective hyperbolic 2-elementary lattices",
        "source_url": "https://arxiv.org/abs/2209.09110",
        "exact_local_lineage": "results/borcherds_char2_h13_h14_lineage.md",
    },
    "D16-D5-H14": {
        "prior_source": "Alexeev, Reflective hyperbolic 2-elementary lattices",
        "source_url": "https://arxiv.org/abs/2209.09110",
        "exact_local_lineage": "results/shimada_disc16_d5_h14_literature_identification.md",
    },
    "X3-H19": {
        "prior_source": "Vinberg, The Two Most Algebraic K3 Surfaces",
        "source_url": "https://doi.org/10.1007/BF01456933",
        "exact_local_lineage": "results/vinberg_x3_h19.md",
    },
    "BOR-D4-H21-210": {
        "prior_source": "Borcherds, Automorphism groups of Lorentzian lattices",
        "source_url": "https://doi.org/10.1016/0021-8693(87)90245-6",
        "exact_local_lineage": "results/higher_dimensional_algebraic_reflection_examples.md",
    },
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _parse_header(text: str) -> tuple[int, int]:
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        facets, dimension = map(int, line.split()[:2])
        return dimension, facets
    raise ArithmeticError("CoxIter file has no header")


def _parse_ledger() -> list[dict[str, Any]]:
    text = LEDGER.read_text(encoding="utf-8")
    table = text.split("## Current catch ledger", 1)[1].split(
        "The five complete K3 nef chambers", 1
    )[0]
    rows: list[dict[str, Any]] = []
    for line in table.splitlines():
        if not line.startswith("| `"):
            continue
        fields = [field.strip() for field in line.strip().strip("|").split("|")]
        if len(fields) != 9:
            raise ArithmeticError(f"unexpected ledger row: {line}")
        rows.append(
            {
                "id": fields[0].strip("`"),
                "dimension": int(fields[1].strip("`")[1:]),
                "facets": int(fields[2]),
                "volume_type": fields[3],
            }
        )
    return rows


def _scan_corpus(corpus: Path) -> tuple[Counter[tuple[int, int]], Counter[tuple[int, int]], int, int]:
    ordinary: Counter[tuple[int, int]] = Counter()
    ordinary_count = 0
    for path in sorted(corpus.rglob("*.coxiter")):
        ordinary[_parse_header(path.read_text(encoding="utf-8"))] += 1
        ordinary_count += 1

    doubling: Counter[tuple[int, int]] = Counter()
    doubling_count = 0
    with tarfile.open(corpus / DOUBLING_RELATIVE, "r:gz") as archive:
        for member in archive:
            if not member.isfile() or not member.name.endswith(".coxiter"):
                continue
            handle = archive.extractfile(member)
            if handle is None:
                raise ArithmeticError(f"failed to read {member.name}")
            doubling[_parse_header(handle.read().decode("utf-8"))] += 1
            doubling_count += 1
    return ordinary, doubling, ordinary_count, doubling_count


def build_certificate(corpus: Path = COXITER_DIRECTORY) -> dict[str, Any]:
    if _sha256(COXITER_ARCHIVE) != COXITER_ARCHIVE_SHA256:
        raise ArithmeticError("pinned ordinary CoxIter archive changed")
    doubling_path = corpus / DOUBLING_RELATIVE
    if _sha256(doubling_path) != DOUBLING_ARCHIVE_SHA256:
        raise ArithmeticError("pinned CoxIter doubling archive changed")

    rows = _parse_ledger()
    ordinary, doubling, ordinary_count, doubling_count = _scan_corpus(corpus)
    row_by_id = {row["id"]: row for row in rows}
    if len(row_by_id) != len(rows):
        raise ArithmeticError("duplicate positive-ledger identifier")

    audited_rows = []
    for row in rows:
        header = (row["dimension"], row["facets"])
        audited_rows.append(
            {
                **row,
                "ordinary_header_count": ordinary[header],
                "doubling_header_count": doubling[header],
                "zero_header_in_both_pinned_lanes": ordinary[header] + doubling[header] == 0,
            }
        )

    controls = []
    for identifier, metadata in KNOWN_ZERO_HEADER_CONTROLS.items():
        row = row_by_id[identifier]
        header = (row["dimension"], row["facets"])
        controls.append(
            {
                "id": identifier,
                "dimension": row["dimension"],
                "facets": row["facets"],
                "ordinary_header_count": ordinary[header],
                "doubling_header_count": doubling[header],
                **metadata,
            }
        )

    unique_headers = {(row["dimension"], row["facets"]) for row in rows}
    zero_unique_headers = sum(
        ordinary[header] + doubling[header] == 0 for header in unique_headers
    )
    dimension_histogram = Counter(row["dimension"] for row in rows)

    x2 = json.loads(X2_LINEAGE.read_text(encoding="utf-8"))
    checks = {
        "ledger_has_56_positive_rows": len(rows) == 56,
        "ledger_has_43_distinct_headers": len(unique_headers) == 43,
        "dimension_histogram_matches": dimension_histogram
        == Counter({10: 29, 11: 14, 12: 6, 13: 2, 14: 3, 19: 1, 21: 1}),
        "ordinary_record_count_1214": ordinary_count == 1214,
        "doubling_record_count_2959": doubling_count == 2959,
        "thirty_of_43_headers_are_absent": zero_unique_headers == 30,
        "thirty_five_of_56_rows_have_absent_headers": sum(
            row["zero_header_in_both_pinned_lanes"] for row in audited_rows
        )
        == 35,
        "all_six_prior_controls_have_zero_headers": all(
            control["ordinary_header_count"] == 0
            and control["doubling_header_count"] == 0
            for control in controls
        ),
        "all_six_control_lineage_reports_exist": all(
            (ROOT / control["exact_local_lineage"]).is_file()
            for control in controls
        ),
        "x2_strengthened_certificate_hash_verified": _sha256(X2_LINEAGE)
        == X2_LINEAGE_SHA256,
        "x2_all_45_same_node_records_excluded": bool(
            x2["checks"][
                "all_45_37_node_records_have_different_relation_histogram"
            ]
        ),
        "x2_worldwide_priority_still_open": x2["claim"]["worldwide_priority"]
        == "open",
    }
    if not all(checks.values()):
        raise ArithmeticError(f"positive-lineage header audit failed: {checks}")

    return {
        "schema": SCHEMA,
        "status_date": "2026-08-07",
        "claim": {
            "theorem": (
                "The two pinned CoxIter lanes are not a worldwide-publication "
                "denominator: six chambers with independently closed prior "
                "lineage have zero matching (dimension, facet-count) headers."
            ),
            "valid_use": (
                "Header absence, or an exact diagram nonmatch within these lanes, "
                "is bounded-corpus evidence only."
            ),
            "invalid_inference": (
                "No worldwide novelty or first-publication claim follows from "
                "absence in the pinned CoxIter lanes."
            ),
            "x2_exception_scope": (
                "X2-RAT-C2-H10-37 has the stronger exact exclusion of all 45 "
                "37-node records in the pinned lanes, regardless of declared "
                "dimension; that stronger result is still bounded and leaves "
                "worldwide priority open."
            ),
        },
        "ledger": {
            "path": LEDGER.relative_to(ROOT).as_posix(),
            "sha256": _sha256(LEDGER),
            "positive_row_count": len(rows),
            "dimension_histogram": {
                f"H{dimension}": dimension_histogram[dimension]
                for dimension in sorted(dimension_histogram)
            },
            "distinct_header_count": len(unique_headers),
            "zero_header_distinct_count": zero_unique_headers,
            "zero_header_row_count": sum(
                row["zero_header_in_both_pinned_lanes"] for row in audited_rows
            ),
            "rows": audited_rows,
        },
        "corpus": {
            "repository": "https://github.com/rgugl/CoxIter",
            "commit": COXITER_COMMIT,
            "ordinary_archive_sha256": COXITER_ARCHIVE_SHA256,
            "doubling_archive_relative_path": DOUBLING_RELATIVE.as_posix(),
            "doubling_archive_sha256": DOUBLING_ARCHIVE_SHA256,
            "ordinary_record_count": ordinary_count,
            "doubling_record_count": doubling_count,
        },
        "known_zero_header_controls": controls,
        "x2_strengthened_bounded_audit": {
            "path": X2_LINEAGE.relative_to(ROOT).as_posix(),
            "sha256": X2_LINEAGE_SHA256,
            "same_node_records_checked": x2["corpus"]["all_37_node_record_count"],
            "histogram_match_count": x2["corpus"]["candidate_histogram_match_count"],
            "worldwide_priority": x2["claim"]["worldwide_priority"],
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate()
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.exists():
            raise SystemExit(f"missing archive: {args.output}")
        archived = json.loads(args.output.read_text(encoding="utf-8"))
        # The ledger SHA is expected to change when wording is clarified; it
        # is provenance, not mathematical content of this header census.
        archived["ledger"]["sha256"] = certificate["ledger"]["sha256"]
        if archived != certificate:
            raise SystemExit("archived positive-lineage header audit differs")
        print("archived positive-lineage header audit matches")
        return
    args.output.write_text(rendered, encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
