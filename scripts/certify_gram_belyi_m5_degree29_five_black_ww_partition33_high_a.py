#!/usr/bin/env python3
"""Exact high-marked-degree no-go inside WW partition (3^2,1^23).

For the 12 canonical black-degree records (a,b,c,d,e) with a>=21, this
script checks the exact factorized system

    D=29(x-t)^2(x-s)^2,
    G(s)=G(t), G(0)=G(t), G(1)=G(t).

One system covers SS,SC,CS,CC.  A uniform collision argument excludes every
projective point at infinity in characteristic p>29.  This is only a strict
sublocus theorem; the remaining 1,280 black-degree systems stay open.
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
from certify_gram_belyi_m5_degree29_five_black_quotient import (
    canonical_black_records,
)
from certify_gram_belyi_m5_degree29_five_black_ww_partition5 import (
    residual_derivative,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data/enumerations/gram_belyi_m5_degree29_five_black_ww_partition33_high_a.json"
)
SCHEMA = "coxeter.gram_belyi.m5_degree29_five_black_ww_partition33_high_a.v1"
SPLIT_PRIMES = (101, 109, 131, 139, 149, 151, 179, 181, 191, 199)
MINIMUM_MARKED_DEGREE = 21
EXPECTED_RECORDS = 12
ENDPOINT_MARKINGS_PER_RECORD = 4
HOMOGENIZATION_DEGREES = (4, 3, 2, 1, 29, 29, 29)


def equation_strings(weights, derivative) -> list[str]:
    """Return D-factor, common-value, and endpoint equations."""

    a, b, c_weight, d, e = weights
    desired = (
        "29*t^2*s^2",
        "-58*t*s*(t+s)",
        "29*(t^2+4*t*s+s^2)",
        "-58*(t+s)",
    )
    critical = [f"({derivative[index]})-({desired[index]})" for index in range(4)]
    g0 = f"(-L)^{a}*(-u)^{b}*(-v)^{c_weight}*(-w)^{d}*(-z)^{e}"
    g1 = (
        f"(1-L)^{a}*(1-u)^{b}*(1-v)^{c_weight}*"
        f"(1-w)^{d}*(1-z)^{e}"
    )
    gt = (
        f"(t-L)^{a}*(t-u)^{b}*(t-v)^{c_weight}*"
        f"(t-w)^{d}*(t-z)^{e}"
    )
    gs = (
        f"(s-L)^{a}*(s-u)^{b}*(s-v)^{c_weight}*"
        f"(s-w)^{d}*(s-z)^{e}"
    )
    return critical + [f"({gs})-({gt})", f"({g0})-({gt})", f"({g1})-({gt})"]


def affine_standard_basis_is_unit(prime: int, lam: int, equations):
    program = "\n".join(
        [
            f"ring R={prime},(u,v,w,z,t,s),dp;",
            "option(redSB);",
            f"number L={lam};",
            *(f"poly f{i}={value};" for i, value in enumerate(equations)),
            "ideal I=f0,f1,f2,f3,f4,f5,f6; ideal J=std(I);",
            'if (size(J)==1 && string(J[1])=="1") { print("UNIT"); } else { print("NONUNIT"); print(dim(J)); print(vdim(J)); print(J); }',
            "quit;",
        ]
    )
    completed = subprocess.run(
        ["Singular", "-q"], input=program, text=True, capture_output=True, check=True
    )
    stdout = completed.stdout.strip()
    return stdout == "UNIT", stdout, completed.stderr.strip()


def certify_one(prime: int, weights) -> dict[str, Any] | None:
    sqrt5, lam = lambda_mod(prime)
    _context, derivative = residual_derivative(prime, lam, weights)
    equations = equation_strings(weights, derivative)
    affine_unit, stdout, stderr = affine_standard_basis_is_unit(prime, lam, equations)
    if stderr:
        raise ArithmeticError((weights, prime, "Singular stderr", stderr))
    if not affine_unit:
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
        "homogenization_degrees": list(HOMOGENIZATION_DEGREES),
        "affine_standard_basis": [1],
        "singular_stdout": stdout,
        "uniform_infinity_lemma_applied": True,
        "equation_sha256": hashlib.sha256("\n".join(equations).encode()).hexdigest(),
    }


def histogram(rows) -> dict[str, int]:
    return {
        str(prime): count
        for prime, count in sorted(Counter(row["prime"] for row in rows).items())
    }


def build_archive(progress_every: int = 5) -> dict[str, Any]:
    records = [
        tuple(value)
        for value in canonical_black_records(0)
        if value[0] >= MINIMUM_MARKED_DEGREE
    ]
    if len(records) != EXPECTED_RECORDS:
        raise ArithmeticError(("high-a record count", len(records)))
    rows = []
    started = time.monotonic()
    for index, weights in enumerate(records, 1):
        found = None
        for prime in SPLIT_PRIMES:
            found = certify_one(prime, weights)
            if found is not None:
                break
        if found is None:
            raise ArithmeticError(("unresolved partition-33 high-a system", weights))
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

    eliminated = EXPECTED_RECORDS * ENDPOINT_MARKINGS_PER_RECORD
    return {
        "schema": SCHEMA,
        "status": "EXACT_FIVE_BLACK_WW_PARTITION_33_HIGH_A_NO_GO",
        "global_degree_29_status": "OPEN",
        "five_black_ww_status": "PARTIAL",
        "white_partition": "(3^2,1^23)",
        "marked_black_degree_condition": "a_lambda>=21",
        "ordered_endpoint_markings_covered": ["SS", "SC", "CS", "CC"],
        "black_degree_records_eliminated": EXPECTED_RECORDS,
        "marked_passports_per_black_record": ENDPOINT_MARKINGS_PER_RECORD,
        "marked_passports_eliminated": eliminated,
        "partition_black_degree_systems_before": 1292,
        "partition_black_degree_systems_remaining": 1292 - EXPECTED_RECORDS,
        "partition_marked_passports_before": 5168,
        "partition_marked_passports_remaining": 5168 - eliminated,
        "five_black_status_free_necessary_systems_remaining": 9742,
        "five_black_binary_status_master_systems_remaining": 22068,
        "five_black_ww_marked_passports_remaining": 29716 - eliminated,
        "five_black_marked_passports_remaining": 131116 - eliminated,
        "global_marked_passports_eliminated": 107274 + eliminated,
        "global_marked_passports_remaining": 42503404 - eliminated,
        "prime_histogram": histogram(rows),
        "projective_infinity_certificate": {
            "homogenized_variables": ["u", "v", "w", "z", "t", "s", "h"],
            "equation_degrees": list(HOMOGENIZATION_DEGREES),
            "lemma": (
                "At h=0 the value equations make t,s black-root clusters. "
                "For t!=s, D=29(x-t)^2(x-s)^2 and E(alpha)!=0 force "
                "n_t=n_s=3, impossible with five labels. For t=s, order four "
                "forces all five labels to collide with the fixed root 0, giving "
                "the forbidden zero projective point."
            ),
            "valid_for_all_selected_primes": True,
        },
        "rows": rows,
        "complete_partition_no_go_claimed": False,
        "complete_five_black_ww_no_go_claimed": False,
        "complete_five_black_no_go_claimed": False,
        "complete_degree_29_no_go_claimed": False,
        "degree_29_witness_found": False,
        "broad_plane_tree_search_used": False,
    }


def verify_archive(archive: dict[str, Any]) -> None:
    checks = {
        "black_degree_records_eliminated": 12,
        "marked_passports_eliminated": 48,
        "partition_black_degree_systems_remaining": 1280,
        "partition_marked_passports_remaining": 5120,
        "five_black_ww_marked_passports_remaining": 29668,
        "five_black_marked_passports_remaining": 131068,
        "global_marked_passports_remaining": 42503356,
    }
    for key, expected in checks.items():
        if archive[key] != expected:
            raise ArithmeticError((key, archive[key], expected))
    if len(archive["rows"]) != EXPECTED_RECORDS:
        raise ArithmeticError("row count changed")
    if sum(archive["prime_histogram"].values()) != EXPECTED_RECORDS:
        raise ArithmeticError("prime histogram changed")
    for key in (
        "complete_partition_no_go_claimed",
        "complete_five_black_ww_no_go_claimed",
        "complete_five_black_no_go_claimed",
        "complete_degree_29_no_go_claimed",
        "degree_29_witness_found",
        "broad_plane_tree_search_used",
    ):
        if archive[key]:
            raise ArithmeticError((key, "must remain false"))


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--progress-every", type=int, default=5)
    arguments = parser.parse_args()
    archive = build_archive(arguments.progress_every)
    verify_archive(archive)
    encoded = canonical_bytes(archive)
    digest = hashlib.sha256(encoded).hexdigest()
    if arguments.check:
        if not arguments.output.exists() or arguments.output.read_bytes() != encoded:
            raise SystemExit("partition-33 high-a archive is stale")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_bytes(encoded)
    print(
        json.dumps(
            {
                "status": archive["status"],
                "systems_eliminated": archive["black_degree_records_eliminated"],
                "marked_passports_eliminated": archive["marked_passports_eliminated"],
                "systems_remaining": archive["partition_black_degree_systems_remaining"],
                "marked_passports_remaining": archive["partition_marked_passports_remaining"],
                "prime_histogram": archive["prime_histogram"],
                "sha256": digest,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
