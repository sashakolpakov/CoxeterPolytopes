#!/usr/bin/env python3
"""Exact CC tower certificate for five-black WW partition (3^2,1^23).

For the white partition (3^2,1^23), the CC marking fixes the two critical
vertices at the normalized endpoints.  Hence K=x(x-1) and

    D=29*K^2,
    G(0)=G(1).

The four coefficient equations have degrees 4,3,2,1 and define a finite
base algebra of dimension 24.  The endpoint-value equation is then adjoined
inside that exact algebra.  After the predetermined 41-record tranche with
17<=a_lambda<=20 closed, the same engine certifies all 1,292 records.  The
12 records with a_lambda>=21 were already eliminated for every endpoint
status, so this full CC theorem removes 1,280 new marked passports.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
from typing import Any

from certify_gram_belyi_m5_degree29_reduction import lambda_mod
from certify_gram_belyi_m5_degree29_five_black_quotient import canonical_black_records
from certify_gram_belyi_m5_degree29_five_black_ww_partition5 import residual_derivative


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data/enumerations/gram_belyi_m5_degree29_five_black_ww_partition33_cc_tower.json"
)
SCHEMA = "coxeter.gram_belyi.m5_degree29_five_black_ww_partition33_cc_tower.v1"
SPLIT_PRIMES = (101, 109, 131, 139, 149, 151, 179, 181, 191, 199)
EXPECTED_RECORDS = 1292
PREDETERMINED_TRANCHE_RECORDS = 41
BASE_QUOTIENT_DIMENSION = 24
STRUCTURAL_HOMOGENIZATION_DEGREES = (4, 3, 2, 1, 29)


def equation_strings(weights, derivative) -> list[str]:
    a, b, c, d, e = weights
    base = [
        str(derivative[0]),
        str(derivative[1]),
        f"({derivative[2]})-29",
        f"({derivative[3]})+58",
    ]
    g0 = f"(-L)^{a}*(-u)^{b}*(-v)^{c}*(-w)^{d}*(-z)^{e}"
    g1 = f"(1-L)^{a}*(1-u)^{b}*(1-v)^{c}*(1-w)^{d}*(1-z)^{e}"
    return base + [f"({g0})-({g1})"]


def successive_basis_certificate(prime: int, lam: int, equations):
    program = "\n".join(
        [
            f"ring R={prime},(u,v,w,z),dp;",
            "option(redSB);",
            f"number L={lam};",
            *(f"poly f{i}={value};" for i, value in enumerate(equations)),
            "ideal I=f0,f1,f2,f3; ideal B=std(I);",
            'print("BASE_DIM"); print(dim(B)); print(vdim(B));',
            "ideal J=std(B+ideal(f4));",
            'if (size(J)==1 && string(J[1])=="1") { print("UNIT"); } else { print("NONUNIT"); print(dim(J)); print(vdim(J)); }',
            "quit;",
        ]
    )
    completed = subprocess.run(
        ["Singular", "-q"], input=program, text=True, capture_output=True, check=True
    )
    stdout = completed.stdout.strip()
    return stdout == "BASE_DIM\n0\n24\nUNIT", stdout, completed.stderr.strip()


def certify_one(prime: int, weights) -> dict[str, Any] | None:
    sqrt5, lam = lambda_mod(prime)
    _context, derivative = residual_derivative(prime, lam, weights)
    equations = equation_strings(weights, derivative)
    success, stdout, stderr = successive_basis_certificate(prime, lam, equations)
    if stderr:
        raise ArithmeticError((weights, prime, "Singular stderr", stderr))
    if not success:
        return None
    return {
        "a_lambda": weights[0],
        "b": weights[1],
        "c": weights[2],
        "d": weights[3],
        "e": weights[4],
        "prime": prime,
        "sqrt5_mod_prime": sqrt5,
        "lambda_mod_prime": lam,
        "base_dimension": 0,
        "base_quotient_dimension": BASE_QUOTIENT_DIMENSION,
        "successive_endpoint_standard_basis": [1],
        "singular_stdout": stdout,
        "structural_homogenization_degrees": list(
            STRUCTURAL_HOMOGENIZATION_DEGREES
        ),
        "uniform_partition_infinity_lemma_applied": True,
        "equation_sha256": hashlib.sha256("\n".join(equations).encode()).hexdigest(),
    }


def histogram(rows) -> dict[str, int]:
    return {
        str(prime): count
        for prime, count in sorted(Counter(row["prime"] for row in rows).items())
    }


def build_archive(progress_every: int = 10) -> dict[str, Any]:
    records = [tuple(value) for value in canonical_black_records(0)]
    if len(records) != EXPECTED_RECORDS:
        raise ArithmeticError(("full CC record count", len(records)))
    rows = []
    started = time.monotonic()
    for index, weights in enumerate(records, 1):
        found = None
        for prime in SPLIT_PRIMES:
            found = certify_one(prime, weights)
            if found is not None:
                break
        if found is None:
            raise ArithmeticError(("unresolved partition-33 CC tower", weights))
        rows.append(found)
        if progress_every and index % progress_every == 0:
            print(
                json.dumps(
                    {
                        "progress": index,
                        "total": len(records),
                        "elapsed_seconds": round(time.monotonic() - started, 3),
                    },
                    sort_keys=True,
                ),
                file=sys.stderr,
                flush=True,
            )

    return {
        "schema": SCHEMA,
        "status": "EXACT_FIVE_BLACK_WW_PARTITION_33_CC_TOWER_NO_GO",
        "global_degree_29_status": "OPEN",
        "five_black_ww_status": "PARTIAL",
        "white_partition": "(3^2,1^23)",
        "endpoint_status": "CC",
        "marked_black_degree_condition": "all canonical records",
        "predetermined_tranche": {
            "condition": "17<=a_lambda<=20",
            "records": PREDETERMINED_TRANCHE_RECORDS,
            "passed_before_full_CC_run": True,
        },
        "CC_black_degree_records_certified": EXPECTED_RECORDS,
        "CC_marked_passports_certified": EXPECTED_RECORDS,
        "previously_eliminated_CC_overlap_from_high_a_theorem": 12,
        "new_marked_passports_eliminated": 1280,
        "base_quotient_dimension_every_row": BASE_QUOTIENT_DIMENSION,
        "partition_black_degree_systems_remaining": 1280,
        "partition_marked_passports_before": 5120,
        "partition_marked_passports_remaining": 3840,
        "five_black_status_free_necessary_systems_remaining": 9742,
        "five_black_binary_status_master_systems_remaining": 22068,
        "five_black_ww_marked_passports_remaining": 28388,
        "five_black_marked_passports_remaining": 129788,
        "global_marked_passports_eliminated": 108602,
        "global_marked_passports_remaining": 42502076,
        "prime_histogram": histogram(rows),
        "successive_exact_algebra": {
            "critical_quadratic_factor": "K=x*(x-1)",
            "base_equations": "D=29*K^2",
            "base_equation_degrees": [4, 3, 2, 1],
            "base_quotient_dimension": BASE_QUOTIENT_DIMENSION,
            "adjoined_value_equation": "G(0)-G(1)=0",
            "final_ideal": [1],
        },
        "projective_infinity_certificate": {
            "structural_degrees": list(STRUCTURAL_HOMOGENIZATION_DEGREES),
            "lemma": (
                "Homogenizing the fixed CC conditions t=0,s=1 as t=0,s=h "
                "gives t=s=0 at h=0. Then D_infinity=29*x^4; the collided-root "
                "factorization and p>29 force all five black labels to collide "
                "at the fixed root 0, the forbidden zero projective point."
            ),
            "valid_for_all_selected_primes": True,
        },
        "rows": rows,
        "complete_CC_no_go_claimed": True,
        "complete_partition_no_go_claimed": False,
        "complete_five_black_ww_no_go_claimed": False,
        "complete_degree_29_no_go_claimed": False,
        "degree_29_witness_found": False,
        "broad_plane_tree_search_used": False,
    }


def verify_archive(archive: dict[str, Any]) -> None:
    checks = {
        "CC_black_degree_records_certified": 1292,
        "CC_marked_passports_certified": 1292,
        "new_marked_passports_eliminated": 1280,
        "partition_black_degree_systems_remaining": 1280,
        "partition_marked_passports_remaining": 3840,
        "five_black_ww_marked_passports_remaining": 28388,
        "five_black_marked_passports_remaining": 129788,
        "global_marked_passports_eliminated": 108602,
        "global_marked_passports_remaining": 42502076,
    }
    for key, expected in checks.items():
        if archive[key] != expected:
            raise ArithmeticError((key, archive[key], expected))
    if len(archive["rows"]) != EXPECTED_RECORDS:
        raise ArithmeticError("CC tranche row count changed")
    if sum(archive["prime_histogram"].values()) != EXPECTED_RECORDS:
        raise ArithmeticError("prime histogram changed")
    for row in archive["rows"]:
        if row["base_quotient_dimension"] != BASE_QUOTIENT_DIMENSION:
            raise ArithmeticError(("base quotient dimension", row))
    for key in (
        "complete_partition_no_go_claimed",
        "complete_five_black_ww_no_go_claimed",
        "complete_degree_29_no_go_claimed",
        "degree_29_witness_found",
        "broad_plane_tree_search_used",
    ):
        if archive[key]:
            raise ArithmeticError((key, "must remain false"))
    if not archive["complete_CC_no_go_claimed"]:
        raise ArithmeticError("complete CC theorem flag must be true")


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--progress-every", type=int, default=10)
    arguments = parser.parse_args()
    archive = build_archive(arguments.progress_every)
    verify_archive(archive)
    encoded = canonical_bytes(archive)
    digest = hashlib.sha256(encoded).hexdigest()
    if arguments.check:
        if not arguments.output.exists() or arguments.output.read_bytes() != encoded:
            raise SystemExit("partition-33 CC tower archive is stale")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_bytes(encoded)
    print(
        json.dumps(
            {
                "status": archive["status"],
                "systems_certified": archive["CC_black_degree_records_certified"],
                "new_marked_passports_eliminated": archive[
                    "new_marked_passports_eliminated"
                ],
                "marked_passports_remaining": archive["partition_marked_passports_remaining"],
                "prime_histogram": archive["prime_histogram"],
                "sha256": digest,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
