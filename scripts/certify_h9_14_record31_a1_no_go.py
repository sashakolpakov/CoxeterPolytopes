#!/usr/bin/env python3
"""Freeze and verify the compact-only Record 31 a=1 incidence no-go."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import tempfile
from typing import Iterable, Mapping, Sequence

from scripts.generate_h9_14_record31_a1_addition_engine import generate


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/h9_14_record31_a1_compact_no_go.json"
POSITIVE_ARCHIVE = ROOT / "data/enumerations/h9_14_record31_positive_a5.json"
TEMPLATE = ROOT / "scripts/explore_h9_14_record31_a2_additions.cpp"
GENERATOR = ROOT / "scripts/generate_h9_14_record31_a1_addition_engine.py"
INDEPENDENT_LINK = ROOT / "scripts/independent_h9_14_record31_a1_link_audit.py"
CERTIFIER = Path(__file__).resolve()

GRAPH6 = "M?????????O?O?o??"
EXPECTED_F_VECTOR = [14, 87, 319, 770, 1281, 1477, 1142, 540, 120]
EXPECTED_MASKS = (60, 2076, 6156, 8220, 10252, 14340)
EXPECTED_ADDITIONS = (5, 6, 7, 7, 7, 7)
EXPECTED_FAILURES = (1,)


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"Record 31 a=1 {name} is malformed")
    return value


def _source(path: Path) -> dict[str, str]:
    return {
        "path": str(path.relative_to(ROOT)),
        "sha256": sha256(path.read_bytes()).hexdigest(),
    }


def _digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    return sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _mask(vertices: Iterable[int]) -> int:
    return sum(1 << int(vertex) for vertex in vertices)


def _numerical_witness_signature(raw: object) -> str:
    witness = _mapping(raw, "strict-link witness")
    payload = {key: value for key, value in witness.items() if key != "subset"}
    return json.dumps(payload, sort_keys=True, separators=(",", ":"))


def _rows(census: Mapping[str, object]) -> list[Mapping[str, object]]:
    raw_rows = census.get("orbits")
    if not isinstance(raw_rows, list) or len(raw_rows) != 6:
        raise ValueError("Record 31 a=1 incidence orbit coverage changed")
    return [_mapping(row, "incidence orbit") for row in raw_rows]


def combine_single_orbit_documents(
    documents: Sequence[Mapping[str, object]], mode: str
) -> dict[str, object]:
    """Combine deterministic one-orbit runs into the ordinary six-row census."""

    if len(documents) != 6:
        raise ValueError("Record 31 a=1 parallel branch coverage changed")
    rows: list[Mapping[str, object]] = []
    for expected_index, document in enumerate(documents):
        if (
            document.get("kind")
            != "h9_14_record31_a1_pivoted_addition_scout"
            or document.get("compact_only") is not True
            or document.get("enumeration_mode") != mode
            or document.get("orbit_count") != 1
        ):
            raise ValueError("Record 31 a=1 single-orbit transcript changed")
        raw_rows = document.get("orbits")
        if not isinstance(raw_rows, list) or len(raw_rows) != 1:
            raise ValueError("Record 31 a=1 single-orbit row changed")
        row = _mapping(raw_rows[0], "single incidence orbit")
        if row.get("orbit_index") != expected_index:
            raise ValueError("Record 31 a=1 parallel row order changed")
        rows.append(row)
    pure_total = sum(int(row["pure_120_facet_target_count"]) for row in rows)
    return {
        "kind": "h9_14_record31_a1_pivoted_addition_scout",
        "compact_only": True,
        "enumeration_mode": mode,
        "orbit_count": 6,
        "orbits": rows,
        "pure_120_facet_target_count_across_representatives": pure_total,
        "scope": (
            "Exact necessary compact incidence census for Record 31 a=1; "
            "no finite-volume inference."
        ),
    }


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record31_a1_compact_schedule_no_go":
        raise ValueError("wrong Record 31 a=1 certificate kind")
    if document.get("graph6") != GRAPH6:
        raise ValueError("Record 31 graph source changed")
    if document.get("schedule") != {
        "a": 1,
        "missing_five_count": 15,
        "facet_count": 120,
        "f_vector": EXPECTED_F_VECTOR,
    }:
        raise ValueError("Record 31 a=1 schedule changed")
    expected_scope = {
        "compact_only": True,
        "entire_a1_schedule_excluded": True,
        "whole_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "new_hyperbolic_coxeter_polytope": False,
    }
    if document.get("scope") != expected_scope:
        raise ValueError("Record 31 a=1 scope changed or was over-promoted")

    production = _mapping(document.get("production"), "production census")
    independent = _mapping(
        document.get("independent_incidence"), "independent census"
    )
    if production.get("enumeration_mode") != "nauty_partial_quotient":
        raise ValueError("Record 31 a=1 production quotient changed")
    if (
        independent.get("enumeration_mode")
        != "direct_binary_matrix_partial_quotient"
    ):
        raise ValueError("Record 31 a=1 independent quotient changed")

    census_rows = (_rows(production), _rows(independent))
    for census, rows in zip((production, independent), census_rows, strict=True):
        if (
            census.get("kind") != "h9_14_record31_a1_pivoted_addition_scout"
            or census.get("compact_only") is not True
            or int(census.get("orbit_count", -1)) != 6
        ):
            raise ValueError("Record 31 a=1 incidence census header changed")
        for index, row in enumerate(rows):
            if row.get("orbit_index") != index:
                raise ValueError("Record 31 a=1 incidence row order changed")
            masks = row.get("positive_A4_masks")
            if (
                not isinstance(masks, list)
                or len(masks) != 1
                or _mask(masks[0]) != EXPECTED_MASKS[index]
            ):
                raise ValueError("Record 31 a=1 positive-A4 source changed")
            if int(row["additional_A5_count"]) != EXPECTED_ADDITIONS[index]:
                raise ValueError("Record 31 a=1 addition schedule changed")
            witnesses = row.get("fixed_local_failure_witnesses")
            if not isinstance(witnesses, list):
                raise ValueError("Record 31 a=1 link witnesses changed")
            if bool(row.get("fixed_local_failure")) != (index in EXPECTED_FAILURES):
                raise ValueError("Record 31 a=1 strict-link partition changed")
            if bool(witnesses) != bool(row.get("fixed_local_failure")):
                raise ValueError("Record 31 a=1 link witness flag changed")
            if index in EXPECTED_FAILURES:
                if len(witnesses) != 4:
                    raise ValueError("Record 31 a=1 witness count changed")
                for raw_witness in witnesses:
                    witness = _mapping(raw_witness, "incidence link witness")
                    if (
                        witness.get("additional_A5_incidence_interval") != [0, -1]
                        or witness.get("g_prefix", [])[-1] != 0
                    ):
                        raise ValueError("Record 31 a=1 link contradiction changed")
            if index not in EXPECTED_FAILURES and int(row["search_node_count"]) <= 0:
                raise ValueError("Record 31 a=1 residual branch was not searched")
            if (
                int(row["ridge_valid_pure_target_count"]) != 0
                or int(row["ridge_valid_target_orbit_count"]) != 0
                or row.get("ridge_valid_target_minimal_nonfaces") != []
            ):
                raise ValueError("a Record 31 a=1 ridge target survived")
        pure_total = sum(int(row["pure_120_facet_target_count"]) for row in rows)
        if int(census["pure_120_facet_target_count_across_representatives"]) != pure_total:
            raise ValueError("Record 31 a=1 pure-target total changed")

    production_nodes = [int(row["search_node_count"]) for row in census_rows[0]]
    independent_nodes = [int(row["search_node_count"]) for row in census_rows[1]]
    if production_nodes == independent_nodes:
        raise ValueError("Record 31 a=1 independent traversal collapsed")

    link = _mapping(document.get("independent_link_audit"), "link audit")
    if (
        link.get("kind") != "independent_h9_14_record31_a1_strict_face_link_audit"
        or link.get("schedule") != document.get("schedule")
        or link.get("positive_A4_orbit_count") != 6
        or link.get("constraints_per_orbit") != 420
        or link.get("total_constraints") != 2520
        or link.get("strict_face_link_excluded_orbit_indices") != [1]
        or link.get("residual_orbit_indices") != [0, 2, 3, 4, 5]
    ):
        raise ValueError("independent Record 31 a=1 link replay changed")
    link_rows = link.get("orbits")
    if not isinstance(link_rows, list) or len(link_rows) != 6:
        raise ValueError("independent Record 31 a=1 link rows changed")
    for index in range(6):
        link_row = _mapping(link_rows[index], "independent link row")
        independent_witnesses = link_row.get("witnesses")
        if not isinstance(independent_witnesses, list):
            raise ValueError("independent Record 31 a=1 witnesses changed")
        production_witnesses = census_rows[0][index][
            "fixed_local_failure_witnesses"
        ]
        if sorted(map(_numerical_witness_signature, independent_witnesses)) != sorted(
            map(_numerical_witness_signature, production_witnesses)
        ):
            raise ValueError("Record 31 a=1 numerical link witnesses differ")

    generated_digest = sha256(generate().encode()).hexdigest()
    if document.get("generated_cpp_sha256") != generated_digest:
        raise ValueError("Record 31 a=1 generated C++ source changed")
    expected_provenance = [
        _source(path)
        for path in (TEMPLATE, GENERATOR, INDEPENDENT_LINK, POSITIVE_ARCHIVE, CERTIFIER)
    ]
    if document.get("source_provenance") != expected_provenance:
        raise ValueError("Record 31 a=1 source provenance changed")
    if document.get("certificate_sha256") != _digest(document):
        raise ValueError("Record 31 a=1 certificate digest changed")
    return {
        "strict_face_link_excluded_orbit_count": 1,
        "enumerated_residual_orbit_count": 5,
        "ridge_target_orbit_count": 0,
        "certificate_sha256": document["certificate_sha256"],
    }


def assemble(
    production: Mapping[str, object],
    independent: Mapping[str, object],
    link: Mapping[str, object],
) -> dict[str, object]:
    document: dict[str, object] = {
        "kind": "h9_14_record31_a1_compact_schedule_no_go",
        "claim": (
            "Every completion of H9/14 graph-shadow Record 31 in the compact "
            "Dehn--Sommerville schedule a=1 fails an exact necessary compact "
            "incidence condition."
        ),
        "graph6": GRAPH6,
        "schedule": {
            "a": 1,
            "missing_five_count": 15,
            "facet_count": 120,
            "f_vector": EXPECTED_F_VECTOR,
        },
        "scope": {
            "compact_only": True,
            "entire_a1_schedule_excluded": True,
            "whole_shadow_compact_coxeter_theorem": False,
            "finite_volume_noncompact_chamber_excluded": False,
            "new_hyperbolic_coxeter_polytope": False,
        },
        "production": production,
        "independent_incidence": independent,
        "independent_link_audit": link,
        "generated_cpp_sha256": sha256(generate().encode()).hexdigest(),
        "source_provenance": [
            _source(path)
            for path in (TEMPLATE, GENERATOR, INDEPENDENT_LINK, POSITIVE_ARCHIVE, CERTIFIER)
        ],
        "nonclaims": [
            "No noncompact finite-volume chamber is excluded.",
            "No parabolic or Vinberg endpoint census is performed.",
            "No topology or Coxeter-label gate is needed because no ridge-valid incidence target survives.",
            "No new hyperbolic Coxeter polytope is constructed.",
            "The compact schedule a=0 of Record 31 remains open.",
        ],
    }
    document["certificate_sha256"] = _digest(document)
    verify(document)
    return document


def _run_mode(binary: Path, mode: str) -> dict[str, object]:
    def run_one(pair: tuple[int, int]) -> Mapping[str, object]:
        index, mask = pair
        result = subprocess.run(
            [str(binary), mode],
            input=f"1\n{index} {mask}\n",
            text=True,
            capture_output=True,
            check=True,
            cwd=ROOT,
        )
        return _mapping(json.loads(result.stdout), "fresh single-orbit output")

    with ThreadPoolExecutor(max_workers=6) as executor:
        documents = list(executor.map(run_one, enumerate(EXPECTED_MASKS)))
    mode_name = {
        "--quotient": "nauty_partial_quotient",
        "--matrix-quotient": "direct_binary_matrix_partial_quotient",
    }[mode]
    return combine_single_orbit_documents(documents, mode_name)


def fresh() -> dict[str, object]:
    with tempfile.TemporaryDirectory(prefix="record31-a1-") as directory:
        temporary = Path(directory)
        source = temporary / "enumerate.cpp"
        binary = temporary / "enumerate"
        source.write_text(generate(), encoding="utf-8")
        subprocess.run(
            [
                "c++", "-std=c++20", "-O3", "-DNDEBUG",
                "-I/opt/homebrew/include/nauty", str(source),
                "-L/opt/homebrew/lib", "-lnauty", "-o", str(binary),
            ],
            check=True,
            cwd=ROOT,
        )
        production = _run_mode(binary, "--quotient")
        independent = _run_mode(binary, "--matrix-quotient")
        link_result = subprocess.run(
            ["python3", str(INDEPENDENT_LINK), str(POSITIVE_ARCHIVE)],
            text=True,
            capture_output=True,
            check=True,
            cwd=ROOT,
        )
        link = _mapping(json.loads(link_result.stdout), "fresh link output")
        return assemble(production, independent, link)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fresh", action="store_true")
    parser.add_argument("--assemble", action="store_true")
    parser.add_argument("--production", type=Path)
    parser.add_argument("--independent", type=Path)
    parser.add_argument("--link", type=Path)
    args = parser.parse_args()
    if args.fresh:
        document = fresh()
        ARCHIVE.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    elif args.assemble:
        if args.production is None or args.independent is None or args.link is None:
            parser.error("--assemble requires --production, --independent, and --link")
        document = assemble(
            _mapping(json.loads(args.production.read_text(encoding="utf-8")), "production input"),
            _mapping(json.loads(args.independent.read_text(encoding="utf-8")), "independent input"),
            _mapping(json.loads(args.link.read_text(encoding="utf-8")), "link input"),
        )
        ARCHIVE.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    else:
        document = _mapping(
            json.loads(ARCHIVE.read_text(encoding="utf-8")), "archive"
        )
    print(json.dumps(verify(document), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
