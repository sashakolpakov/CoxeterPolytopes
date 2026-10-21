"""Bounded exact-match audit against the official CoxIter graph corpus.

This scan uses only format invariants that every literal chamber match must
preserve: hyperbolic dimension and number of facets (Coxeter graph vertices).
It therefore certifies absence of this exact 60-facet H14 chamber from the
specified corpus snapshot.  It does not test whether the underlying lattice,
some commensurable group, or another chamber for it occurs in the literature.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CORPUS_ROOT = Path("/private/tmp/coxiter-current/CoxIter-master")
DEFAULT_ARCHIVE = Path("/private/tmp/coxiter-current.zip")
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_d5_h14_coxiter_corpus_audit.json"
)
COXITER_COMMIT = "b800be48d8240fbabd2d1c9542e61054584ad7fe"
COXITER_COMMIT_DATE = "2025-07-19T07:15:31Z"
MASTER_VERIFIED_DATE = "2026-08-06"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _first_nonempty_line(path: Path) -> str:
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.strip():
                return line.strip()
    raise ValueError("empty CoxIter file")


def build_audit(corpus_root: Path, archive: Path | None = None) -> dict[str, Any]:
    records: list[dict[str, Any]] = []
    errors: list[dict[str, str]] = []
    for path in sorted(corpus_root.rglob("*.coxiter")):
        relative = path.relative_to(corpus_root).as_posix()
        try:
            tokens = _first_nonempty_line(path).split()
            if len(tokens) < 2:
                raise ValueError("header has fewer than two fields")
            facet_count, dimension = map(int, tokens[:2])
            if facet_count <= 0 or dimension <= 0:
                raise ValueError("nonpositive facet count or dimension")
            records.append(
                {
                    "path": relative,
                    "facet_count": facet_count,
                    "hyperbolic_dimension": dimension,
                }
            )
        except (OSError, ValueError) as error:
            errors.append({"path": relative, "error": str(error)})

    if not records or errors:
        raise ArithmeticError(
            f"expected a nonempty fully parseable corpus; records={len(records)}, "
            f"errors={len(errors)}"
        )

    pair_histogram = Counter(
        (record["hyperbolic_dimension"], record["facet_count"])
        for record in records
    )
    dimension_14 = [
        record for record in records if record["hyperbolic_dimension"] == 14
    ]
    maximum_facets = max(record["facet_count"] for record in records)
    maximum_facet_records = [
        record for record in records if record["facet_count"] == maximum_facets
    ]
    candidate_dimension = 14
    candidate_facets = 60
    literal_count_matches = [
        record
        for record in records
        if record["hyperbolic_dimension"] == candidate_dimension
        and record["facet_count"] == candidate_facets
    ]

    result: dict[str, Any] = {
        "schema": "shimada_disc16_d5_h14_coxiter_corpus_audit_v1",
        "candidate": {
            "hyperbolic_dimension": candidate_dimension,
            "facet_count": candidate_facets,
            "intrinsic_lattice": "U(2) + E8(-1) + A1(-1)^5",
        },
        "corpus_snapshot": {
            "repository": "https://github.com/rgugliel/CoxIter",
            "commit": COXITER_COMMIT,
            "commit_url": f"https://github.com/rgugliel/CoxIter/commit/{COXITER_COMMIT}",
            "commit_date": COXITER_COMMIT_DATE,
            "master_commit_verified_on": MASTER_VERIFIED_DATE,
            "master_verification_method": "GitHub REST API GET /repos/rgugliel/CoxIter/commits/master",
            "archive_sha256": _sha256(archive) if archive is not None else None,
        },
        "scan": {
            "coxiter_file_count": len(records),
            "parse_error_count": len(errors),
            "dimension_facet_histogram": {
                f"H{dimension}_facets_{facets}": count
                for (dimension, facets), count in sorted(pair_histogram.items())
            },
            "global_maximum_facet_count": maximum_facets,
            "global_maximum_facet_records": maximum_facet_records,
            "h14_record_count": len(dimension_14),
            "h14_records": dimension_14,
            "h14_facet_histogram": {
                str(facets): count
                for facets, count in sorted(
                    Counter(record["facet_count"] for record in dimension_14).items()
                )
            },
            "candidate_dimension_and_facet_count_matches": literal_count_matches,
        },
        "decision": {
            "exact_60_facet_h14_chamber_absent_from_snapshot": not literal_count_matches,
            "reason": (
                "a literal Coxeter-chamber match must preserve hyperbolic dimension "
                "and the number of graph vertices/facets"
            ),
            "does_not_claim_underlying_lattice_absent": True,
            "does_not_claim_global_literature_novelty": True,
        },
    }
    if not (
        len(records) == 1214
        and maximum_facets == 50
        and len(dimension_14) == 3
        and Counter(record["facet_count"] for record in dimension_14)
        == Counter({17: 2, 20: 1})
        and not literal_count_matches
    ):
        raise ArithmeticError("the pinned CoxIter corpus census changed")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus-root", type=Path, default=DEFAULT_CORPUS_ROOT)
    parser.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    result = build_audit(arguments.corpus_root, arguments.archive)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file() or arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archived artifact differs: {arguments.output}")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
    print(
        json.dumps(
            {
                "coxiter_file_count": result["scan"]["coxiter_file_count"],
                "global_maximum_facet_count": result["scan"]["global_maximum_facet_count"],
                "h14_records": result["scan"]["h14_records"],
                "exact_match_absent": result["decision"]["exact_60_facet_h14_chamber_absent_from_snapshot"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
