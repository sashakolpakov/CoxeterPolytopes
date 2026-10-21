#!/usr/bin/env python3
"""Bounded literal-lineage audit for the X_2 rational-involution H10/37.

The composite geometric certificate already records that no file with header
``37 10`` occurs in the pinned CoxIter ordinary and index-two-doubling lanes.
This script makes the comparison stronger without pretending to settle
worldwide priority: it scans *every* CoxIter record having 37 vertices,
regardless of the dimension declared in its header.  The complete histogram
of CoxIter relation labels is an isomorphism invariant, so a histogram
mismatch rules out an exact relation-diagram match before graph isomorphism.

CoxIter uses label ``0`` for parallel walls and label ``1`` for a generic
ultraparallel pair.  Thus this audit compares the complete CoxIter relation
matrix, not the finer metric carried by a dotted edge with an explicit
distance parameter.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import tarfile
from typing import Any, Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = (
    ROOT
    / "data/enumerations/shimada_aut3singk3_x2_rational_involution_h10_37.json"
)
OUTPUT = (
    ROOT
    / "data/enumerations/shimada_aut3singk3_x2_h10_37_global_lineage.json"
)
COXITER_DIRECTORY = Path("/private/tmp/coxiter-current/CoxIter-master")
COXITER_ARCHIVE = Path("/private/tmp/coxiter-current.zip")
DOUBLING_RELATIVE = Path("graphs/testing_index2.tar.gz")

SCHEMA = "coxeter.shimada.aut3singk3.x2.h10_37.global_lineage.v1"
COXITER_COMMIT = "b800be48d8240fbabd2d1c9542e61054584ad7fe"
COXITER_ARCHIVE_SHA256 = (
    "ae04648b291c1e5f3d0e47d666c846dc805b244469a20486bca59674b4d171c3"
)
DOUBLING_ARCHIVE_SHA256 = (
    "3330ca41fb3e047123ace718f08542bc1dd0019b35922d4fa836a217647e7036"
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_sha256(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _parse_coxiter_text(text: str) -> tuple[int, int, list[list[str]]]:
    lines = [
        content
        for line in text.splitlines()
        if (content := line.split("#", 1)[0].strip())
    ]
    facets, dimension = map(int, lines[0].split()[:2])
    matrix = [
        ["1" if left == right else "2" for right in range(facets)]
        for left in range(facets)
    ]
    vertex_line = next(
        (line for line in lines[1:] if line.startswith("vertices labels:")),
        None,
    )
    labels = (
        vertex_line.split(":", 1)[1].split()
        if vertex_line is not None
        else [str(index + 1) for index in range(facets)]
    )
    if len(labels) != facets or len(set(labels)) != facets:
        raise ArithmeticError("invalid CoxIter vertex labels")
    index_of = {label: index for index, label in enumerate(labels)}
    for line in lines[1:]:
        if line.startswith("vertices labels:"):
            continue
        tokens = line.split()
        if len(tokens) != 3:
            raise ArithmeticError(f"unsupported CoxIter line: {line}")
        left, right = index_of[tokens[0]], index_of[tokens[1]]
        matrix[left][right] = matrix[right][left] = tokens[2]
    return facets, dimension, matrix


def _relation_histogram(matrix: Sequence[Sequence[str]]) -> dict[str, int]:
    counts: Counter[str] = Counter()
    for left in range(len(matrix)):
        for right in range(left + 1, len(matrix)):
            counts[matrix[left][right]] += 1
    return {label: counts[label] for label in sorted(counts)}


def _candidate_relations(payload: Sequence[Sequence[str]]) -> list[list[str]]:
    # The archived facet Gram has negative diagonal and nonnegative off-diagonal.
    gram = [[-Fraction(entry) for entry in row] for row in payload]
    size = len(gram)
    labels = {
        Fraction(1, 4): "3",
        Fraction(1, 2): "4",
        Fraction(3, 4): "6",
        Fraction(1): "0",
    }
    result = [
        ["1" if left == right else "2" for right in range(size)]
        for left in range(size)
    ]
    for left in range(size):
        if len(gram[left]) != size or gram[left][left] <= 0:
            raise ArithmeticError("invalid candidate Gram matrix")
        for right in range(left + 1, size):
            product = gram[left][right]
            if product != gram[right][left] or product > 0:
                raise ArithmeticError("candidate Gram is not symmetric acute")
            if product == 0:
                relation = "2"
            else:
                ratio = product * product / (
                    gram[left][left] * gram[right][right]
                )
                relation = labels.get(ratio, "1" if ratio > 1 else "")
                if not relation:
                    raise ArithmeticError(f"unsupported Coxeter ratio {ratio}")
            result[left][right] = result[right][left] = relation
    return result


def _histogram_classes(records: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    classes: dict[tuple[int, tuple[tuple[str, int], ...]], list[str]] = {}
    for record in records:
        key = (
            record["declared_dimension"],
            tuple(sorted(record["relation_histogram"].items())),
        )
        classes.setdefault(key, []).append(record["path"])
    return [
        {
            "declared_dimension": dimension,
            "relation_histogram": dict(histogram),
            "record_count": len(paths),
            "paths": sorted(paths),
        }
        for (dimension, histogram), paths in sorted(classes.items())
    ]


def build_certificate(
    candidate_path: Path = CANDIDATE,
    corpus: Path = COXITER_DIRECTORY,
    corpus_archive: Path = COXITER_ARCHIVE,
) -> dict[str, Any]:
    if _sha256(corpus_archive) != COXITER_ARCHIVE_SHA256:
        raise ArithmeticError("the pinned CoxIter ordinary archive changed")
    doubling_path = corpus / DOUBLING_RELATIVE
    if _sha256(doubling_path) != DOUBLING_ARCHIVE_SHA256:
        raise ArithmeticError("the pinned CoxIter doubling archive changed")

    candidate = json.loads(candidate_path.read_text(encoding="utf-8"))
    relations = _candidate_relations(
        candidate["rational_involution_fold"]["facet_gram"]
    )
    candidate_histogram = _relation_histogram(relations)

    ordinary_count = 0
    ordinary_37: list[dict[str, Any]] = []
    for path in sorted(corpus.rglob("*.coxiter")):
        facets, dimension, matrix = _parse_coxiter_text(
            path.read_text(encoding="utf-8")
        )
        ordinary_count += 1
        if facets == 37:
            ordinary_37.append(
                {
                    "path": path.relative_to(corpus).as_posix(),
                    "declared_dimension": dimension,
                    "relation_histogram": _relation_histogram(matrix),
                }
            )

    doubling_count = 0
    doubling_37: list[dict[str, Any]] = []
    with tarfile.open(doubling_path, "r:gz") as archive:
        for member in archive:
            if not member.isfile() or not member.name.endswith(".coxiter"):
                continue
            handle = archive.extractfile(member)
            if handle is None:
                raise ArithmeticError(f"failed to read {member.name}")
            facets, dimension, matrix = _parse_coxiter_text(
                handle.read().decode("utf-8")
            )
            doubling_count += 1
            if facets == 37:
                doubling_37.append(
                    {
                        "path": member.name,
                        "declared_dimension": dimension,
                        "relation_histogram": _relation_histogram(matrix),
                    }
                )

    all_37 = [
        *({"lane": "ordinary", **record} for record in ordinary_37),
        *({"lane": "doubling", **record} for record in doubling_37),
    ]
    histogram_matches = [
        record
        for record in all_37
        if record["relation_histogram"] == candidate_histogram
    ]
    checks = {
        "pinned_ordinary_archive_verified": True,
        "pinned_doubling_archive_verified": True,
        "ordinary_record_count_1214": ordinary_count == 1214,
        "doubling_record_count_2959": doubling_count == 2959,
        "candidate_has_37_facets": len(relations) == 37,
        "ordinary_37_node_record_count_2": len(ordinary_37) == 2,
        "doubling_37_node_record_count_43": len(doubling_37) == 43,
        "all_45_37_node_records_have_different_relation_histogram": not histogram_matches,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"bounded lineage audit failed: {checks}")

    return {
        "schema": SCHEMA,
        "id": "SHIMADA-X2-RATIONAL-C2-H10-37-BOUNDED-LITERAL-LINEAGE",
        "status_date": "2026-08-07",
        "claim": {
            "bounded_result": (
                "No 37-node record in the two pinned CoxIter lanes can be "
                "isomorphic to the candidate CoxIter relation diagram: all 45 "
                "are already excluded by their complete relation-label histogram."
            ),
            "worldwide_priority": "open",
            "scope_warning": (
                "This is an exact bounded-corpus exclusion, not a theorem that "
                "the abstract H10/37 polytope was never published elsewhere."
            ),
        },
        "candidate": {
            "composite_certificate_path": candidate_path.relative_to(ROOT).as_posix(),
            "composite_certificate_sha256": _sha256(candidate_path),
            "facet_count": 37,
            "declared_dimension": 10,
            "relation_histogram": candidate_histogram,
            "relation_histogram_semantics": {
                "0": "parallel",
                "1": "generic ultraparallel",
                "2": "m=2",
                "3": "m=3",
                "4": "m=4",
            },
        },
        "corpus": {
            "repository": "https://github.com/rgugl/CoxIter",
            "commit": COXITER_COMMIT,
            "ordinary_archive_sha256": COXITER_ARCHIVE_SHA256,
            "doubling_archive_relative_path": DOUBLING_RELATIVE.as_posix(),
            "doubling_archive_sha256": DOUBLING_ARCHIVE_SHA256,
            "ordinary_record_count": ordinary_count,
            "doubling_record_count": doubling_count,
            "ordinary_37_node_records": ordinary_37,
            "doubling_37_node_histogram_classes": _histogram_classes(doubling_37),
            "all_37_node_record_count": len(all_37),
            "all_37_node_records_canonical_sha256": _canonical_sha256(all_37),
            "candidate_histogram_match_count": len(histogram_matches),
            "candidate_histogram_matches": histogram_matches,
        },
        "method": (
            "Scan every .coxiter record, ignore declared dimension when selecting "
            "37-node diagrams, and compare the complete multiset of all 666 "
            "off-diagonal CoxIter relation labels.  Equality of diagrams implies "
            "equality of this histogram, so mismatch is a rigorous obstruction."
        ),
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=Path, default=CANDIDATE)
    parser.add_argument("--corpus", type=Path, default=COXITER_DIRECTORY)
    parser.add_argument("--corpus-archive", type=Path, default=COXITER_ARCHIVE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()

    certificate = build_certificate(
        candidate_path=args.candidate,
        corpus=args.corpus,
        corpus_archive=args.corpus_archive,
    )
    if args.check_archive:
        archived = json.loads(args.output.read_text(encoding="utf-8"))
        if archived != certificate:
            raise SystemExit("archived bounded-lineage certificate differs")
        print("archived bounded-lineage certificate matches")
        return
    print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
