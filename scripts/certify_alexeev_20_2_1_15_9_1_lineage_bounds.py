#!/usr/bin/env python3
"""Bounded lineage and rational-class audit for three Alexeev sections.

This deliberately proves less than worldwide novelty.  It identifies the
one 2-elementary fixed lattice by Nikulin's invariants, compares all three
rational quadratic spaces by exact local Hasse invariants, and checks whether
the literal facet counts occur anywhere in the pinned ordinary or index-two
doubling CoxIter snapshots.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import tarfile
from typing import Any, Iterable

import sympy as sp

from scripts.certify_higher_dimensional_rational_lineage import (
    _compare,
    _e8_cartan,
)
from scripts.screen_alexeev_20_2_1_15_9_1_fixed_sections import (
    OUTPUT as SECTION_CENSUS,
    ROOT,
    _matrix,
)


COXITER_ROOT = Path("/private/tmp/coxiter-current/CoxIter-master")
COXITER_ARCHIVE = Path("/private/tmp/coxiter-current.zip")
DOUBLING_ARCHIVE = COXITER_ROOT / "graphs" / "testing_index2.tar.gz"
COXITER_COMMIT = "b800be48d8240fbabd2d1c9542e61054584ad7fe"
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_20_2_1_15_9_1_lineage_bounds.json"
)
SCHEMA = "coxeter.alexeev_20_2_1_15_9_1_lineage_bounds.v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _header(lines: Iterable[str]) -> tuple[int, int]:
    content = next(
        line.split("#", 1)[0].strip()
        for line in lines
        if line.split("#", 1)[0].strip()
    )
    facets, dimension = map(int, content.split()[:2])
    return dimension, facets


def _scan_coxiter(
    targets: set[tuple[int, int]],
) -> dict[str, Any]:
    ordinary_histogram: Counter[tuple[int, int]] = Counter()
    ordinary_matches = {target: [] for target in targets}
    for path in sorted(COXITER_ROOT.rglob("*.coxiter")):
        dimension, facets = _header(
            path.read_text(encoding="utf-8", errors="replace").splitlines()
        )
        key = (dimension, facets)
        ordinary_histogram[key] += 1
        if key in ordinary_matches:
            ordinary_matches[key].append(path.relative_to(COXITER_ROOT).as_posix())

    doubling_histogram: Counter[tuple[int, int]] = Counter()
    doubling_matches = {target: [] for target in targets}
    with tarfile.open(DOUBLING_ARCHIVE, "r:gz") as archive:
        members = [
            member
            for member in archive
            if member.isfile() and member.name.endswith(".coxiter")
        ]
        for member in members:
            handle = archive.extractfile(member)
            if handle is None:
                raise ArithmeticError(f"cannot read {member.name}")
            dimension, facets = _header(
                raw.decode("utf-8", errors="replace") for raw in handle
            )
            key = (dimension, facets)
            doubling_histogram[key] += 1
            if key in doubling_matches:
                doubling_matches[key].append(member.name)

    target_rows = {}
    for dimension, facets in sorted(targets):
        key = (dimension, facets)
        target_rows[f"H{dimension}_facets_{facets}"] = {
            "ordinary_header_match_count": ordinary_histogram[key],
            "ordinary_header_matches": ordinary_matches[key],
            "doubling_header_match_count": doubling_histogram[key],
            "doubling_header_matches": doubling_matches[key],
            "absent_from_both_snapshots_by_dimension_and_facet_count": not (
                ordinary_histogram[key] or doubling_histogram[key]
            ),
        }
    return {
        "repository": "https://github.com/rgugl/CoxIter",
        "commit": COXITER_COMMIT,
        "archive_path": str(COXITER_ARCHIVE),
        "archive_sha256": _sha256(COXITER_ARCHIVE),
        "ordinary_coxiter_file_count": sum(ordinary_histogram.values()),
        "doubling_archive_path": str(DOUBLING_ARCHIVE),
        "doubling_archive_sha256": _sha256(DOUBLING_ARCHIVE),
        "doubling_coxiter_file_count": sum(doubling_histogram.values()),
        "targets": target_rows,
    }


def _delta(lattice: sp.Matrix) -> int | None:
    """Return Nikulin delta for a 2-elementary even lattice, else None."""

    inverse = lattice.inv()
    if not all(2 * value in sp.ZZ for value in inverse):
        return None
    return int(any(value.q == 2 for value in inverse.diagonal()))


def _rational_comparisons(lattice: sp.Matrix) -> dict[str, Any]:
    dimension = lattice.rows - 1
    references = {
        "f1_odd_unimodular": sp.diag(1, *([-1] * dimension)),
        "f2": sp.diag(2, *([-1] * dimension)),
        "f3_McLeod": sp.diag(3, *([-1] * dimension)),
    }
    rows = {}
    for name, reference in references.items():
        scalar = Fraction(int(lattice.det()), int(reference.det()))
        rows[name] = _compare(lattice, reference, scalar)
    return rows


def build() -> dict[str, Any]:
    sections = json.loads(SECTION_CENSUS.read_text(encoding="utf-8"))
    cases = []
    targets = set()
    for row in sections["cases"]:
        lattice = _matrix(row["saturated_fixed_lattice"]["gram"])
        dimension = int(row["case"]["hyperbolic_dimension"])
        facets = int(row["exact_reduction"]["active_hyperbolic_facet_count"])
        targets.add((dimension, facets))
        comparisons = _rational_comparisons(lattice)
        delta = _delta(lattice)
        two_elementary = delta is not None
        classification = None
        if row["id"] == "A20_2_1_C02":
            if not (
                two_elementary
                and delta == 1
                and lattice.rows == 12
                and abs(int(lattice.det())) == 2**10
            ):
                raise ArithmeticError("the H11 2-elementary invariants changed")
            standard = sp.diag(
                sp.Matrix([[0, 1], [1, 0]]),
                -2 * _e8_cartan(),
                sp.Matrix([[-2]]),
                sp.Matrix([[-2]]),
            )
            if not (
                standard.det() == lattice.det()
                and _delta(standard) == delta
            ):
                raise ArithmeticError("standard (12,10,1) lattice check failed")
            classification = {
                "nikulin_invariants": [12, 10, 1],
                "isometry_type": "U + E8(-2) + <-2>^2",
                "identification_method": (
                    "uniqueness of indefinite even 2-elementary lattices from "
                    "signature and (r,a,delta)"
                ),
                "known_complete_reflection_chamber": {
                    "source": "Alexeev arXiv:2209.09110v4, case (12,10,1)",
                    "facet_count": 14,
                },
                "section_group_has_finite_index_in_known_full_reflection_group": True,
                "index_value": None,
            }
        cases.append(
            {
                "id": row["id"],
                "hyperbolic_dimension": dimension,
                "facet_count": facets,
                "fixed_lattice_determinant": int(lattice.det()),
                "fixed_lattice_smith_diagonal_absolute": row[
                    "saturated_fixed_lattice"
                ]["smith_diagonal_absolute"],
                "all_facet_reflections_integral_on_fixed_lattice": row[
                    "restricted_mirrors"
                ]["all_active_reflections_integral_on_saturated_fixed_lattice"],
                "two_elementary_delta": delta,
                "integral_classification": classification,
                "rational_similarity_comparisons": comparisons,
                "arithmetic_commensurability_conclusion": (
                    "odd-unimodular/f1 rational class"
                    if comparisons["f1_odd_unimodular"]["rationally_similar"]
                    else "not the odd-unimodular/f1 rational class"
                ),
            }
        )

    coxiter = _scan_coxiter(targets)
    checks = {
        "three_cases": len(cases) == 3,
        "all_reflections_integral": all(
            row["all_facet_reflections_integral_on_fixed_lattice"] for row in cases
        ),
        "all_three_in_f1_rational_class": all(
            row["rational_similarity_comparisons"]["f1_odd_unimodular"][
                "rationally_similar"
            ]
            for row in cases
        ),
        "h11_is_known_12_10_1_lattice": next(
            row for row in cases if row["id"] == "A20_2_1_C02"
        )["integral_classification"]
        is not None,
        "all_literal_facet_counts_absent_from_pinned_coxiter": all(
            row["absent_from_both_snapshots_by_dimension_and_facet_count"]
            for row in coxiter["targets"].values()
        ),
        "ordinary_snapshot_has_1214_files": coxiter[
            "ordinary_coxiter_file_count"
        ]
        == 1214,
        "doubling_snapshot_has_2959_files": coxiter[
            "doubling_coxiter_file_count"
        ]
        == 2959,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"lineage-bound checks failed: {checks}")
    return {
        "schema": SCHEMA,
        "source": {
            "fixed_section_census": str(SECTION_CENSUS.relative_to(ROOT)),
            "fixed_section_census_sha256": _sha256(SECTION_CENSUS),
        },
        "scope": {
            "positive_statement": (
                "exact integral/rational lattice lineage and literal absence "
                "from two pinned CoxIter snapshots"
            ),
            "not_claimed": (
                "worldwide novelty, absence from all literature, or a new "
                "arithmetic commensurability class"
            ),
        },
        "cases": cases,
        "coxiter_bounded_comparison": coxiter,
        "checks": checks,
        "all_checks_pass": True,
    }


def main() -> None:
    payload = build()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(OUTPUT)


if __name__ == "__main__":
    main()
