#!/usr/bin/env python3
"""Freeze and verify the exact mixed-E8 projected-wall shell-3 no-go."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Iterable, Iterator

from scripts.independent_og10_golden_projected_wall_mixed_shell3 import (
    add as independent_add,
    coefficient_records as independent_coefficient_records,
    is_signature_admissible as independent_is_signature_admissible,
    iter_independent_mixed_records,
    run_independent_replay,
)
from scripts.search_og10_golden_projected_wall import (
    _add as production_add,
    _archived_vectors_match_source,
    _coefficient_records as production_coefficient_records,
    _is_signature_admissible as production_is_signature_admissible,
    _mixed_e8_records,
    _source_hashes as production_source_hashes,
    audit_mixed_e8_family,
)


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE_PATH = (
    ROOT
    / "data/candidates/og10_golden_projected_wall_mixed_shell3_no_go.json"
)
INDEPENDENT_SCRIPT_PATH = (
    ROOT
    / "scripts/independent_og10_golden_projected_wall_mixed_shell3.py"
)
COEFFICIENT_BOUND = 30
B_SHELL_BOUND = 3

EXPECTED_ELIGIBLE_COUNT = 569_760
EXPECTED_ELIGIBLE_HASH = (
    "f3437fc0d251f3c0d0fb0729bb3e317e72676ba52988002b1f4c1aa759608ecf"
)
EXPECTED_WALL_TEST_COUNT = 1_472_640
EXPECTED_WALL_HASH = (
    "851270ff6c0ab5ea5afefb7e7890e957c20ca31bac6fe2dcfa3d3114fe08ad73"
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_json_sha256(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def local_source_hashes() -> dict[str, object]:
    return {
        "production_dependencies": production_source_hashes(),
        "independent_replay": sha256_file(INDEPENDENT_SCRIPT_PATH),
        "certifier": sha256_file(Path(__file__).resolve()),
    }


def selected_agreement(
    production: dict[str, object], independent: dict[str, object]
) -> dict[str, bool]:
    fields = (
        "b_vector_count",
        "eligible_w_count",
        "eligible_w_square_class_count",
        "eligible_w_records_sha256",
        "signature_admissible_parameter_pair_count_summed_over_w_square_classes",
        "exact_wall_test_count",
        "canonical_wall_records_sha256",
        "cutting_covector_count",
        "cutting_records",
    )
    return {field: production[field] == independent[field] for field in fields}


def iter_wall_records(
    records: Iterable[object],
    coefficients: tuple[object, ...],
    squares: tuple[object, ...],
    add: object,
    is_signature_admissible: object,
) -> Iterator[object]:
    admissible_by_square = {
        square: tuple(
            data
            for data in coefficients
            if is_signature_admissible(add(data[2], square))
        )
        for square in squares
    }
    for record in records:
        square, _, e8_interior, e8_ray, a_vector, b_vector = record
        for A, C, non_e8_square, non_e8_interior, non_e8_ray in (
            admissible_by_square[square]
        ):
            yield (
                a_vector,
                b_vector,
                A,
                C,
                add(non_e8_square, square),
                add(non_e8_interior, e8_interior),
                add(non_e8_ray, e8_ray),
            )


def recompute() -> tuple[dict[str, object], dict[str, object], dict[str, object]]:
    production = audit_mixed_e8_family(COEFFICIENT_BOUND, B_SHELL_BOUND)
    independent = run_independent_replay(COEFFICIENT_BOUND)
    production_records = tuple(
        sorted(_mixed_e8_records(B_SHELL_BOUND), key=lambda record: (record[5], record[4]))
    )

    literal_eligible_count = 0
    for index, (left, right) in enumerate(
        zip(production_records, iter_independent_mixed_records(), strict=True)
    ):
        if left != right:
            raise AssertionError(f"first production/independent mismatch at {index}")
        literal_eligible_count += 1
    if literal_eligible_count != EXPECTED_ELIGIBLE_COUNT:
        raise AssertionError(
            f"literal eligible comparison covered {literal_eligible_count} records"
        )

    production_coefficients = production_coefficient_records(COEFFICIENT_BOUND)
    independent_coefficients = independent_coefficient_records(COEFFICIENT_BOUND)
    coefficient_literal_equality = production_coefficients == independent_coefficients
    if not coefficient_literal_equality:
        raise AssertionError("production/independent coefficient records differ")

    squares = tuple(sorted({record[0] for record in production_records}))
    if len(squares) != 8:
        raise AssertionError(f"expected eight w-square classes, found {len(squares)}")
    production_walls = iter_wall_records(
        production_records,
        production_coefficients,
        squares,
        production_add,
        production_is_signature_admissible,
    )
    independent_walls = iter_wall_records(
        iter_independent_mixed_records(),
        independent_coefficients,
        squares,
        independent_add,
        independent_is_signature_admissible,
    )
    literal_wall_count = 0
    for index, (left, right) in enumerate(
        zip(production_walls, independent_walls, strict=True)
    ):
        if left != right:
            raise AssertionError(f"first production/independent wall mismatch at {index}")
        literal_wall_count += 1
    if literal_wall_count != EXPECTED_WALL_TEST_COUNT:
        raise AssertionError(
            f"literal wall comparison covered {literal_wall_count} records"
        )

    agreement = selected_agreement(production, independent)
    if not all(agreement.values()):
        raise AssertionError(f"production/independent summaries differ: {agreement}")
    comparison = {
        "literal_eligible_record_equality_before_hash": True,
        "literal_eligible_record_count": literal_eligible_count,
        "literal_coefficient_record_equality_before_hash": coefficient_literal_equality,
        "literal_coefficient_record_count": len(production_coefficients),
        "literal_wall_record_equality_before_hash": True,
        "literal_wall_record_count": literal_wall_count,
        "summary_field_agreement": agreement,
        "hashes_are_checksums_after_literal_equality_not_substitutes_for_it": True,
    }
    return production, independent, comparison


def build_certificate() -> dict[str, object]:
    production, independent, comparison = recompute()
    if not _archived_vectors_match_source():
        raise AssertionError("archived interior/ray coordinates changed")
    pinned = {
        "b_vector_count": production["b_vector_count"] == 9120,
        "eligible_w_count": production["eligible_w_count"]
        == EXPECTED_ELIGIBLE_COUNT,
        "eligible_w_square_class_count": production["eligible_w_square_class_count"]
        == 8,
        "eligible_w_records_hash": production["eligible_w_records_sha256"]
        == EXPECTED_ELIGIBLE_HASH,
        "signature_parameter_pair_count": production[
            "signature_admissible_parameter_pair_count_summed_over_w_square_classes"
        ]
        == 144,
        "wall_test_count": production["exact_wall_test_count"]
        == EXPECTED_WALL_TEST_COUNT,
        "wall_records_hash": production["canonical_wall_records_sha256"]
        == EXPECTED_WALL_HASH,
        "no_cutting_covectors": production["cutting_covector_count"] == 0,
        "independent_sign_partition_has_no_opposite_or_zero_pair": independent[
            "separation_sign_histogram"
        ]
        == {"-1,-1": 736320, "1,1": 736320},
    }
    if not all(pinned.values()):
        raise AssertionError(f"shell3 pinned theorem data failed: {pinned}")

    payload: dict[str, object] = {
        "schema_version": "1.0.0",
        "candidate_id": "og10_golden_projected_wall_mixed_e8_shell3_no_go",
        "certificate_kind": "exact_replayable_finite_box_no_go",
        "field": "Q(sqrt(5))",
        "ring": "Z[phi]",
        "fixed_lattice": "F=<-2*phi>+E8+<2>",
        "family": "w=a+phi*b with a,b in E8",
        "box": {
            "A_and_C_integral_coefficient_interval": [-30, 30],
            "b_square_values": [2, 4, 6],
            "b_shell_bound_cli_value": 3,
            "conjugate_E8_gate": "0 < (w')^2 < 2",
            "signature_gate": "0 < q=2*p^2 < 4 and q' < 4",
            "cut_gate": "ell(interior)*ell(v61)<0",
        },
        "production_audit": production,
        "independent_replay": independent,
        "comparison": comparison,
        "completeness": {
            "production": (
                "For every theta-pinned E8 b-vector with b^2<=6, exact "
                "prefix Schur complements enumerate all integral a with "
                "(a+phi'*b)^2<2."
            ),
            "independent": (
                "The inequality a^2<4+12*(phi')^2<10 and E8 evenness "
                "reduce a to complete theta-pinned shells 0,2,4,6,8. The "
                "single Weyl orbits of b-square 2,4,6 transport the complete "
                "eligible seed lists bijectively to all 9,120 b-vectors."
            ),
        },
        "checks": {
            "archived_vectors_match_fixed_slice": True,
            "all_signature_and_separation_tests_are_exact": True,
            "production_and_independent_records_are_literally_equal": True,
            "pinned_counts_and_hashes_match": pinned,
        },
        "source_sha256": local_source_hashes(),
        "commands": [
            "PYTHONPATH=src:. python3 scripts/certify_og10_golden_projected_wall_mixed_shell3.py --recompute",
            "python3 scripts/independent_og10_golden_projected_wall_mixed_shell3.py --compare-certificate data/candidates/og10_golden_projected_wall_mixed_shell3_no_go.json",
            "PYTHONPATH=src:. python3 scripts/certify_og10_golden_projected_wall_mixed_shell3.py --verify-certificate",
        ],
        "conclusion": (
            "No covector in the exact mixed-E8 b^2<=6, coefficient-bound-30 "
            "box separates the archived invariant interior from v61."
        ),
        "scope": (
            "Finite-box theorem only. It does not cover b^2>6, coefficients "
            "outside [-30,30], other fixed carriers, a relaxed reflection "
            "involution, or multi-root/rank-two coupling; it is neither a "
            "global OG10 obstruction nor a new hyperbolic Coxeter example."
        ),
    }
    payload["certificate_payload_sha256"] = canonical_json_sha256(payload)
    return payload


def verify_certificate(
    payload: dict[str, object], *, full_recompute: bool = False
) -> dict[str, bool]:
    body = dict(payload)
    digest = body.pop("certificate_payload_sha256", None)
    comparison = payload["comparison"]
    production = payload["production_audit"]
    independent = payload["independent_replay"]
    checks = {
        "payload_digest_matches": digest == canonical_json_sha256(body),
        "source_hashes_match": payload.get("source_sha256") == local_source_hashes(),
        "archived_vectors_match_fixed_slice": _archived_vectors_match_source(),
        "literal_equality_was_recorded": comparison[
            "literal_eligible_record_equality_before_hash"
        ]
        is True
        and comparison["literal_coefficient_record_equality_before_hash"] is True
        and comparison["literal_wall_record_equality_before_hash"] is True,
        "summary_fields_agree": all(selected_agreement(production, independent).values()),
        "eligible_count_is_pinned": production["eligible_w_count"]
        == EXPECTED_ELIGIBLE_COUNT,
        "eligible_hash_is_pinned": production["eligible_w_records_sha256"]
        == EXPECTED_ELIGIBLE_HASH,
        "wall_count_is_pinned": production["exact_wall_test_count"]
        == EXPECTED_WALL_TEST_COUNT,
        "wall_hash_is_pinned": production["canonical_wall_records_sha256"]
        == EXPECTED_WALL_HASH,
        "zero_cutters_is_pinned": production["cutting_covector_count"] == 0,
    }
    if full_recompute:
        replayed_production, replayed_independent, replayed_comparison = recompute()
        checks.update(
            {
                "production_full_replay_matches": replayed_production == production,
                "independent_full_replay_matches": replayed_independent == independent,
                "literal_record_equality_replays": replayed_comparison == comparison,
            }
        )
    if not all(checks.values()):
        raise AssertionError(f"shell3 certificate verification failed: {checks}")
    return checks


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-certificate", action="store_true")
    parser.add_argument("--verify-certificate", action="store_true")
    parser.add_argument("--recompute", action="store_true")
    args = parser.parse_args()
    if args.write_certificate:
        payload = build_certificate()
        CERTIFICATE_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(f"wrote {CERTIFICATE_PATH}")
        print(payload["certificate_payload_sha256"])
        return
    payload = json.loads(CERTIFICATE_PATH.read_text())
    checks = verify_certificate(payload, full_recompute=args.recompute)
    print(json.dumps(checks, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
