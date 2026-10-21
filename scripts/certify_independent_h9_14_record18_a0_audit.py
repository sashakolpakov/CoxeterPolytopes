#!/usr/bin/env python3
"""Freeze and cross-check the corrected record-18 ``a=0`` survivor census.

The independent census and the production census use different recursive
searches.  This certifier compares every literal survivor family and also
compares the multisets obtained from an explicit incidence-matrix
canonicalizer.  It does not use graph hashes or graph-isomorphism packages.

The expensive independent enumeration may be supplied as JSON with
``--independent-json``.  With ``--fresh-independent`` it is rebuilt and run
from source.  The production census is always rebuilt and rerun.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import re
import subprocess
import tempfile
from typing import Sequence


ROOT = Path(__file__).resolve().parents[1]
INDEPENDENT_SOURCE = (
    ROOT / "scripts/independent_h9_14_record18_a0_balanced_cut_audit.cpp"
)
PRODUCTION_SOURCE = ROOT / "scripts/probe_h9_14_record18_a0_balanced_cuts.cpp"
DEFAULT_OUTPUT = (
    ROOT / "data/enumerations/h9_14_record18_a0_independent_audit.json"
)
SURVIVOR_RE = re.compile(r"^SURVIVOR key=([^ ]+) rows=([0-9,]+)$")
EXPECTED_KEYS = {
    (1, 2, 4, 5, 5, 6, 8, 10, 29, 46, 46, 81, 81, 166, 193, 194): 36,
    (1, 2, 4, 5, 5, 8, 10, 10, 14, 14, 57, 70, 70, 177, 177, 225): 144,
    (1, 2, 4, 5, 6, 8, 9, 10, 13, 13, 22, 22, 105, 178, 225, 226): 16,
}


def compile_cpp(source: Path, output: Path) -> None:
    subprocess.run(
        ["c++", "-std=c++20", "-O3", str(source), "-o", str(output)],
        cwd=ROOT,
        check=True,
    )


def independent_document(binary: Path, supplied: Path | None) -> dict[str, object]:
    if supplied is not None:
        return json.loads(supplied.read_text(encoding="utf-8"))
    completed = subprocess.run(
        [str(binary)], cwd=ROOT, check=True, capture_output=True, text=True
    )
    return json.loads(completed.stdout)


def production_families(binary: Path) -> tuple[tuple[int, ...], ...]:
    completed = subprocess.run(
        [str(binary)], cwd=ROOT, check=True, capture_output=True, text=True
    )
    families = []
    for line in completed.stdout.splitlines():
        match = SURVIVOR_RE.match(line)
        if match is None:
            continue
        rows = tuple(
            sorted(int(value) for value in match.group(2).split(",") if value)
        )
        if len(rows) != 8 or len(set(rows)) != 8:
            raise AssertionError("malformed production survivor family")
        families.append(rows)
    if len(families) != 196 or len(set(families)) != 196:
        raise AssertionError("production survivor family count changed")
    return tuple(families)


def canonical_keys(
    binary: Path, families: Sequence[Sequence[int]]
) -> tuple[tuple[int, ...], ...]:
    arguments = [
        str(binary),
        "--canonicalize",
        *(",".join(map(str, family)) for family in families),
    ]
    completed = subprocess.run(
        arguments, cwd=ROOT, check=True, capture_output=True, text=True
    )
    return tuple(tuple(map(int, key)) for key in json.loads(completed.stdout))


def key_record(counter: Counter[tuple[int, ...]]) -> list[dict[str, object]]:
    return [
        {"canonical_family_key": list(key), "multiplicity": multiplicity}
        for key, multiplicity in sorted(counter.items())
    ]


def build_record(
    independent: dict[str, object],
    production: Sequence[Sequence[int]],
    canonicalizer: Path,
) -> dict[str, object]:
    raw_survivors = independent.get("survivors")
    if not isinstance(raw_survivors, list):
        raise AssertionError("independent survivor archive is missing")
    independent_families = tuple(
        tuple(sorted(map(int, survivor["missing_five_masks"])))
        for survivor in raw_survivors
    )
    if len(independent_families) != 196 or len(set(independent_families)) != 196:
        raise AssertionError("independent survivor family count changed")

    independent_keys = canonical_keys(canonicalizer, independent_families)
    production_keys = canonical_keys(canonicalizer, production)
    independent_counter = Counter(independent_keys)
    production_counter = Counter(production_keys)
    if independent_counter != Counter(EXPECTED_KEYS):
        raise AssertionError("independent canonical orbit census changed")
    if production_counter != independent_counter:
        raise AssertionError("production/independent normalized-family mismatch")

    literal_independent = set(independent_families)
    literal_production = set(map(tuple, production))
    if literal_independent != literal_production:
        raise AssertionError("production/independent literal-family mismatch")

    independent_rows = [
        {
            "missing_five_masks": list(family),
            "canonical_family_key": list(key),
        }
        for family, key in sorted(
            zip(independent_families, independent_keys, strict=True)
        )
    ]
    production_rows = [
        {
            "missing_five_masks": list(family),
            "canonical_family_key": list(key),
        }
        for family, key in sorted(zip(production, production_keys, strict=True))
    ]
    representatives = []
    for key in sorted(independent_counter):
        representative = min(
            family
            for family, observed in zip(
                independent_families, independent_keys, strict=True
            )
            if observed == key
        )
        representatives.append(
            {
                "canonical_family_key": list(key),
                "multiplicity": independent_counter[key],
                "representative_missing_five_masks": list(representative),
            }
        )

    totals = independent.get("totals")
    if not isinstance(totals, dict) or totals.get("ridge_valid") != 196:
        raise AssertionError("independent corrected census totals changed")
    return {
        "kind": "independent_h9_14_record18_a0_exact_family_audit",
        "compact_only": True,
        "finite_volume_noncompact_checked": False,
        "theorem_certificate": False,
        "independent_search_totals": totals,
        "independent_survivor_count": len(independent_families),
        "production_survivor_count": len(production),
        "literal_family_sets_equal": True,
        "canonical_normalized_multisets_equal": True,
        "canonical_orbit_count": len(independent_counter),
        "canonical_orbits": key_record(independent_counter),
        "canonical_representatives": representatives,
        "independent_survivors": independent_rows,
        "production_survivors": production_rows,
        "scope": (
            "Exact corrected a=0 incidence census comparison. This freezes all "
            "196 literal families and their explicit incidence-matrix normal "
            "forms. Topology and Coxeter labels remain separate gates."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--independent-json", type=Path)
    parser.add_argument("--fresh-independent", action="store_true")
    parser.add_argument("--write", type=Path)
    args = parser.parse_args()
    if args.independent_json is None and not args.fresh_independent:
        parser.error("choose --independent-json or --fresh-independent")
    if args.independent_json is not None and args.fresh_independent:
        parser.error("the independent input modes are mutually exclusive")

    with tempfile.TemporaryDirectory(prefix="record18-a0-audit-") as directory:
        temporary = Path(directory)
        independent_binary = temporary / "independent"
        production_binary = temporary / "production"
        compile_cpp(INDEPENDENT_SOURCE, independent_binary)
        compile_cpp(PRODUCTION_SOURCE, production_binary)
        independent = independent_document(
            independent_binary,
            None if args.fresh_independent else args.independent_json,
        )
        production = production_families(production_binary)
        record = build_record(independent, production, independent_binary)

    rendered = json.dumps(record, indent=2, sort_keys=True) + "\n"
    if args.write is not None:
        args.write.parent.mkdir(parents=True, exist_ok=True)
        args.write.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
