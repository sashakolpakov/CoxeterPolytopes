#!/usr/bin/env python3
"""Exact no-go certificate for the five-black WW partition (5,1^24).

Write

    G(x)=(x-lambda)^a(x-u)^b(x-v)^c(x-w)^d(x-z)^e.

For the white partition (5,1^24), the residual derivative quartic is
29(x-t)^4.  The two normalized source endpoints are white precisely only if
G(0)=G(t)=G(1), with t=0 or t=1 allowed for the CS and SC markings.  Thus one
closed six-equation system covers the three ordered endpoint markings SS,
CS,SC for each of the 1,292 canonical black-degree records.

At a good split prime this script proves the affine ideal is the unit ideal.
A uniform leading-form argument excludes projective points at infinity: the
endpoint equations force t to be a black root, while the factorization of the
residual derivative at collided black roots makes a fourth-order zero at t
possible only when all five labelled roots coincide.  Since the marked root
is 0 at infinity, that would be the forbidden zero projective point.

This is a bounded coefficient census.  It proves no complete five-black WW
or degree-29 no-go.
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

from flint import nmod_mpoly_ctx

from certify_gram_belyi_m5_degree29_reduction import lambda_mod
from certify_gram_belyi_m5_degree29_five_black_quotient import (
    canonical_black_records,
)
from certify_gram_belyi_m5_degree29_four_black_ww_simple import (
    outer_add,
    outer_mul,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data/enumerations/gram_belyi_m5_degree29_five_black_ww_partition5.json"
)
SCHEMA = "coxeter.gram_belyi.m5_degree29_five_black_ww_partition5.v1"
SPLIT_PRIMES = (101, 109, 131, 139, 149, 151, 179, 181, 191, 199)
DERIVATIVE_HOMOGENIZATION_DEGREES = (4, 3, 2, 1)
ENDPOINT_HOMOGENIZATION_DEGREES = (29, 29)
EXPECTED_RECORDS = 1292
ENDPOINT_MARKINGS_PER_RECORD = 3  # SS, SC, CS for (5,1^24)


def residual_derivative(
    prime: int, lam: int, weights: tuple[int, int, int, int, int]
):
    """Return the residual derivative quartic coefficients over F_p."""

    context = nmod_mpoly_ctx.get(("u", "v", "w", "z"), prime, "degrevlex")
    u, v, w, z = context.gens()
    roots = (context.constant(lam), u, v, w, z)
    factors = [[-root, context.constant(1)] for root in roots]
    derivative = [context.constant(0)]
    for omitted, weight in enumerate(weights):
        term = [context.constant(1)]
        for index, factor in enumerate(factors):
            if index != omitted:
                term = outer_mul(term, factor)
        derivative = outer_add(
            derivative, [weight * coefficient for coefficient in term]
        )
    if len(derivative) != 5 or derivative[-1] != context.constant(29):
        raise ArithmeticError((weights, "residual derivative is not monic times 29"))
    return context, derivative


def equation_strings(
    lam: int,
    weights: tuple[int, int, int, int, int],
    derivative,
) -> tuple[list[str], list[str]]:
    """Return derivative-factor and endpoint equations for Singular."""

    a, b, c_weight, d, e = weights
    desired = ("29*t^4", "-116*t^3", "174*t^2", "-116*t")
    derivative_equations = [
        f"({derivative[index]})-({desired[index]})" for index in range(4)
    ]
    g0 = f"(-L)^{a}*(-u)^{b}*(-v)^{c_weight}*(-w)^{d}*(-z)^{e}"
    g1 = (
        f"(1-L)^{a}*(1-u)^{b}*(1-v)^{c_weight}*"
        f"(1-w)^{d}*(1-z)^{e}"
    )
    gt = (
        f"(t-L)^{a}*(t-u)^{b}*(t-v)^{c_weight}*"
        f"(t-w)^{d}*(t-z)^{e}"
    )
    return derivative_equations, [f"({g0})-({gt})", f"({g1})-({gt})"]


def affine_standard_basis_is_unit(
    prime: int,
    lam: int,
    derivative_equations: list[str],
    endpoint_equations: list[str],
) -> tuple[bool, str, str]:
    """Check the exact six-equation ideal over F_p by Singular."""

    equations = derivative_equations + endpoint_equations
    program = "\n".join(
        [
            f"ring r={prime},(u,v,w,z,t),dp;",
            "option(redSB);",
            f"number L={lam};",
            *(f"poly f{index}={value};" for index, value in enumerate(equations)),
            "ideal I=f0,f1,f2,f3,f4,f5;",
            "ideal J=std(I);",
            'if (size(J)==1 && string(J[1])=="1") { print("UNIT"); } else { print("NONUNIT"); print(dim(J)); print(vdim(J)); print(J); }',
            "quit;",
        ]
    )
    completed = subprocess.run(
        ["Singular", "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    stdout = completed.stdout.strip()
    return stdout == "UNIT", stdout, completed.stderr.strip()


def certify_one(
    prime: int, weights: tuple[int, int, int, int, int]
) -> dict[str, Any] | None:
    sqrt5, lam = lambda_mod(prime)
    _context, derivative = residual_derivative(prime, lam, weights)
    derivative_equations, endpoint_equations = equation_strings(
        lam, weights, derivative
    )

    actual_derivative_degrees = tuple(
        max(int(derivative[index].total_degree()), 4 - index)
        for index in range(4)
    )
    if any(
        actual > bound
        for actual, bound in zip(
            actual_derivative_degrees, DERIVATIVE_HOMOGENIZATION_DEGREES
        )
    ):
        raise ArithmeticError(
            (weights, "derivative homogenization bound", actual_derivative_degrees)
        )

    affine_unit, stdout, stderr = affine_standard_basis_is_unit(
        prime, lam, derivative_equations, endpoint_equations
    )
    if stderr:
        raise ArithmeticError((weights, prime, "Singular stderr", stderr))
    if not affine_unit:
        return None

    encoded_equations = "\n".join(derivative_equations + endpoint_equations).encode()
    return {
        "a_lambda": weights[0],
        "b": weights[1],
        "c": weights[2],
        "d": weights[3],
        "e": weights[4],
        "prime": prime,
        "sqrt5_mod_prime": sqrt5,
        "lambda_mod_prime": lam,
        "derivative_homogenization_degrees": list(
            DERIVATIVE_HOMOGENIZATION_DEGREES
        ),
        "endpoint_homogenization_degrees": list(
            ENDPOINT_HOMOGENIZATION_DEGREES
        ),
        "affine_standard_basis": [1],
        "singular_stdout": stdout,
        "equation_sha256": hashlib.sha256(encoded_equations).hexdigest(),
        "uniform_infinity_lemma_applied": True,
    }


def histogram(rows: list[dict[str, Any]]) -> dict[str, int]:
    return {
        str(prime): count
        for prime, count in sorted(Counter(row["prime"] for row in rows).items())
    }


def build_archive(progress_every: int = 100) -> dict[str, Any]:
    records = [tuple(value) for value in canonical_black_records(0)]
    if len(records) != EXPECTED_RECORDS:
        raise ArithmeticError(("five-black M0 record count", len(records)))

    rows = []
    started = time.monotonic()
    for index, weights in enumerate(records, 1):
        found = None
        for prime in SPLIT_PRIMES:
            found = certify_one(prime, weights)
            if found is not None:
                break
        if found is None:
            raise ArithmeticError(("unresolved five-black WW partition-5 system", weights))
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

    eliminated_passports = EXPECTED_RECORDS * ENDPOINT_MARKINGS_PER_RECORD
    return {
        "schema": SCHEMA,
        "status": "EXACT_FIVE_BLACK_WW_PARTITION_5_NO_GO",
        "global_degree_29_status": "OPEN",
        "five_black_ww_status": "PARTIAL",
        "white_partition": "(5,1^24)",
        "ordered_endpoint_markings_covered": ["SS", "SC", "CS"],
        "black_degree_records": EXPECTED_RECORDS,
        "exact_partition_systems_eliminated": EXPECTED_RECORDS,
        "marked_passports_per_black_record": ENDPOINT_MARKINGS_PER_RECORD,
        "eliminated_five_black_binary_status_master_systems": 0,
        "ww_binary_status_master_systems_before_collapse": 5168,
        "ww_status_free_master_systems_after_collapse": 1292,
        "status_free_master_systems_after_collapse_by_endpoint_block": {
            "BW": 4225,
            "WB": 4225,
            "WW": 1292,
            "total": 9742,
        },
        "remaining_five_black_status_free_necessary_systems": 9742,
        "eliminated_five_black_marked_passports": eliminated_passports,
        "remaining_five_black_binary_status_master_systems": 22068,
        "remaining_five_black_marked_passports": 134992 - eliminated_passports,
        "remaining_five_black_ww_marked_passports": 33592 - eliminated_passports,
        "global_marked_passports_eliminated": 103398 + eliminated_passports,
        "global_marked_passports_remaining": 42507280 - eliminated_passports,
        "prime_histogram": histogram(rows),
        "projective_infinity_certificate": {
            "homogenized_variables": ["u", "v", "w", "z", "t", "s"],
            "marked_root_at_infinity": 0,
            "derivative_equation_degrees": list(
                DERIVATIVE_HOMOGENIZATION_DEGREES
            ),
            "endpoint_equation_degrees": list(
                ENDPOINT_HOMOGENIZATION_DEGREES
            ),
            "lemma": (
                "G_infinity(t)=0 makes t a black root; grouping collided labelled "
                "roots alpha with total weights A_alpha gives "
                "D_infinity=product(x-alpha)^(n_alpha-1) E and "
                "E(alpha)=A_alpha*product_{beta!=alpha}(alpha-beta)!=0 for p>29. "
                "Therefore ord_t(D_infinity)=4 forces n_t=5; the fixed root 0 "
                "then makes every projective coordinate zero."
            ),
            "valid_for_all_selected_primes": True,
        },
        "remaining_exact_ww_partition_rows": {
            "(4,2,1^23)": 9044,
            "(3^2,1^23)": 5168,
            "(3,2^2,1^22)": 10336,
            "(2^4,1^21)": 5168,
            "total": 29716,
        },
        "rows": rows,
        "degree_29_witness_found": False,
        "complete_five_black_ww_no_go_claimed": False,
        "complete_five_black_no_go_claimed": False,
        "complete_degree_29_no_go_claimed": False,
        "broad_plane_tree_search_used": False,
    }


def verify_archive(archive: dict[str, Any]) -> None:
    if archive["exact_partition_systems_eliminated"] != 1292:
        raise ArithmeticError("partition-5 system count changed")
    if archive["eliminated_five_black_marked_passports"] != 3876:
        raise ArithmeticError("partition-5 marked-passport count changed")
    if archive["remaining_five_black_binary_status_master_systems"] != 22068:
        raise ArithmeticError("binary master count must remain conservative")
    if archive["ww_status_free_master_systems_after_collapse"] != 1292:
        raise ArithmeticError("WW status-free system quotient changed")
    if archive["remaining_five_black_status_free_necessary_systems"] != 9742:
        raise ArithmeticError("five-black status-free system remainder changed")
    if archive["status_free_master_systems_after_collapse_by_endpoint_block"] != {
        "BW": 4225,
        "WB": 4225,
        "WW": 1292,
        "total": 9742,
    }:
        raise ArithmeticError("five-black status-free block quotient changed")
    if archive["remaining_five_black_marked_passports"] != 131116:
        raise ArithmeticError("five-black raw remainder changed")
    if archive["remaining_five_black_ww_marked_passports"] != 29716:
        raise ArithmeticError("five-black WW raw remainder changed")
    if archive["global_marked_passports_remaining"] != 42503404:
        raise ArithmeticError("global marked-passport remainder changed")
    if sum(archive["prime_histogram"].values()) != 1292:
        raise ArithmeticError("prime histogram does not cover every row")
    if len(archive["rows"]) != 1292:
        raise ArithmeticError("row count changed")
    if archive["complete_five_black_ww_no_go_claimed"]:
        raise ArithmeticError("five-black WW remains open")
    if archive["complete_five_black_no_go_claimed"]:
        raise ArithmeticError("five-black stratum remains open")
    if archive["complete_degree_29_no_go_claimed"]:
        raise ArithmeticError("degree 29 remains open")


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--progress-every", type=int, default=100)
    arguments = parser.parse_args()
    archive = build_archive(arguments.progress_every)
    verify_archive(archive)
    encoded = canonical_bytes(archive)
    digest = hashlib.sha256(encoded).hexdigest()
    if arguments.check:
        if not arguments.output.exists() or arguments.output.read_bytes() != encoded:
            raise SystemExit("five-black WW partition-5 archive is stale")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_bytes(encoded)
    print(
        json.dumps(
            {
                "status": archive["status"],
                "eliminated_exact_partition_systems": archive[
                    "exact_partition_systems_eliminated"
                ],
                "eliminated_marked_passports": archive[
                    "eliminated_five_black_marked_passports"
                ],
                "remaining_five_black_ww_marked_passports": archive[
                    "remaining_five_black_ww_marked_passports"
                ],
                "prime_histogram": archive["prime_histogram"],
                "sha256": digest,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
