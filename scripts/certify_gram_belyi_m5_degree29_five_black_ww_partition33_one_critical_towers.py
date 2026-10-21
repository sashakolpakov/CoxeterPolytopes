#!/usr/bin/env python3
"""Exact CS theorem and bounded SC tranche in WW partition (3^2,1^23).

Endpoint statuses are ordered as (x=0,x=1).  The two separately formulated
one-critical-endpoint orientations are:

    CS: K=x(x-r),       G(r)=G(0), G(1)=G(0);
    SC: K=(x-1)(x-r),   G(r)=G(1), G(0)=G(1).

The full 1,280-record CS orientation and the predetermined 41-record SC
tranche 17<=a_lambda<=20 are certified.  A strategic stop forbids completing
the remaining SC census here.  The derivative ideal is one-dimensional.
Adjoining the critical-value equation gives an exact zero-dimensional algebra
of dimension 696; adjoining the simple-endpoint value equation then gives the
unit ideal.  Reflection is not used: lambda and 1-lambda do not satisfy the
same marked quadratic.
"""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
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
    / "data/enumerations/gram_belyi_m5_degree29_five_black_ww_partition33_one_critical_towers.json"
)
SCHEMA = "coxeter.gram_belyi.m5_degree29_five_black_ww_partition33_one_critical_towers.v1"
SPLIT_PRIMES = (101, 109, 131, 139, 149, 151, 179, 181, 191, 199)
ORIENTATIONS = ("CS", "SC")
MAXIMUM_UNRESOLVED_MARKED_DEGREE = 20
EXPECTED_RECORDS_PER_ORIENTATION = 1280
PREDETERMINED_TRANCHE_MINIMUM_A = 17
PREDETERMINED_TRANCHE_RECORDS_PER_ORIENTATION = 41
CRITICAL_VALUE_QUOTIENT_DIMENSION = 696
STRUCTURAL_HOMOGENIZATION_DEGREES = (4, 3, 2, 1, 29, 29)


def equation_strings(weights, derivative, orientation: str) -> list[str]:
    a, b, c, d, e = weights
    g0 = f"(-L)^{a}*(-u)^{b}*(-v)^{c}*(-w)^{d}*(-z)^{e}"
    g1 = f"(1-L)^{a}*(1-u)^{b}*(1-v)^{c}*(1-w)^{d}*(1-z)^{e}"
    gr = f"(r-L)^{a}*(r-u)^{b}*(r-v)^{c}*(r-w)^{d}*(r-z)^{e}"

    if orientation == "CS":
        derivative_equations = [
            str(derivative[0]),
            str(derivative[1]),
            f"({derivative[2]})-29*r^2",
            f"({derivative[3]})+58*r",
        ]
        value_equations = [f"({gr})-({g0})", f"({g1})-({g0})"]
    elif orientation == "SC":
        derivative_equations = [
            f"({derivative[0]})-29*r^2",
            f"({derivative[1]})+58*r*(1+r)",
            f"({derivative[2]})-29*(1+4*r+r^2)",
            f"({derivative[3]})+58*(1+r)",
        ]
        value_equations = [f"({gr})-({g1})", f"({g0})-({g1})"]
    else:
        raise ValueError(orientation)
    return derivative_equations + value_equations


def successive_basis_certificate(prime: int, lam: int, equations):
    program = "\n".join(
        [
            f"ring R={prime},(r,u,v,w,z),dp;",
            "option(redSB);",
            f"number L={lam};",
            *(f"poly f{i}={value};" for i, value in enumerate(equations)),
            "ideal I=f0,f1,f2,f3; ideal B=std(I);",
            'print("DERIVATIVE_BASE"); print(dim(B));',
            "ideal C=std(B+ideal(f4));",
            'print("CRITICAL_VALUE_BASE"); print(dim(C)); print(vdim(C));',
            "ideal J=std(C+ideal(f5));",
            'if (size(J)==1 && string(J[1])=="1") { print("UNIT"); } else { print("NONUNIT"); print(dim(J)); print(vdim(J)); }',
            "quit;",
        ]
    )
    completed = subprocess.run(
        ["Singular", "-q"], input=program, text=True, capture_output=True, check=True
    )
    stdout = completed.stdout.strip()
    expected = "DERIVATIVE_BASE\n1\nCRITICAL_VALUE_BASE\n0\n696\nUNIT"
    return stdout == expected, stdout, completed.stderr.strip()


def certify_one(prime: int, weights, orientation: str) -> dict[str, Any] | None:
    sqrt5, lam = lambda_mod(prime)
    _context, derivative = residual_derivative(prime, lam, weights)
    equations = equation_strings(weights, derivative, orientation)
    success, stdout, stderr = successive_basis_certificate(prime, lam, equations)
    if stderr:
        raise ArithmeticError((orientation, weights, prime, "Singular stderr", stderr))
    if not success:
        return None
    return {
        "orientation": orientation,
        "a_lambda": weights[0],
        "b": weights[1],
        "c": weights[2],
        "d": weights[3],
        "e": weights[4],
        "prime": prime,
        "sqrt5_mod_prime": sqrt5,
        "lambda_mod_prime": lam,
        "derivative_base_dimension": 1,
        "critical_value_base_dimension": 0,
        "critical_value_quotient_dimension": CRITICAL_VALUE_QUOTIENT_DIMENSION,
        "successive_simple_endpoint_standard_basis": [1],
        "singular_stdout": stdout,
        "structural_homogenization_degrees": list(
            STRUCTURAL_HOMOGENIZATION_DEGREES
        ),
        "uniform_partition_infinity_lemma_applied": True,
        "equation_sha256": hashlib.sha256("\n".join(equations).encode()).hexdigest(),
    }


def certify_task(task) -> dict[str, Any]:
    orientation, weights = task
    for prime in SPLIT_PRIMES:
        found = certify_one(prime, weights, orientation)
        if found is not None:
            return found
    raise ArithmeticError(
        ("unresolved partition-33 one-critical tower", orientation, weights)
    )


def records(tranche_only: bool):
    answer = [
        tuple(value)
        for value in canonical_black_records(0)
        if value[0] <= MAXIMUM_UNRESOLVED_MARKED_DEGREE
    ]
    if len(answer) != EXPECTED_RECORDS_PER_ORIENTATION:
        raise ArithmeticError(("one-critical record count", len(answer)))
    if tranche_only:
        answer = [
            value for value in answer if value[0] >= PREDETERMINED_TRANCHE_MINIMUM_A
        ]
        if len(answer) != PREDETERMINED_TRANCHE_RECORDS_PER_ORIENTATION:
            raise ArithmeticError(("predetermined tranche count", len(answer)))
    return answer


def histogram(rows) -> dict[str, int]:
    return {
        str(prime): count
        for prime, count in sorted(Counter(row["prime"] for row in rows).items())
    }


def orientation_histograms(rows) -> dict[str, dict[str, int]]:
    return {
        orientation: histogram(
            [row for row in rows if row["orientation"] == orientation]
        )
        for orientation in ORIENTATIONS
    }


def build_archive(
    progress_every: int = 10, tranche_only: bool = False, workers: int = 1
) -> dict[str, Any]:
    black_records = records(tranche_only)
    rows = []
    started = time.monotonic()
    if tranche_only:
        tasks = [
            (orientation, weights)
            for orientation in ORIENTATIONS
            for weights in black_records
        ]
    else:
        sc_tranche = [
            value
            for value in black_records
            if value[0] >= PREDETERMINED_TRANCHE_MINIMUM_A
        ]
        if len(sc_tranche) != PREDETERMINED_TRANCHE_RECORDS_PER_ORIENTATION:
            raise ArithmeticError(("SC frozen tranche count", len(sc_tranche)))
        tasks = [("CS", weights) for weights in black_records] + [
            ("SC", weights) for weights in sc_tranche
        ]
    total = len(tasks)
    if workers == 1:
        iterator = map(certify_task, tasks)
        executor = None
    else:
        # The heavy algebra runs in isolated Singular subprocesses.  Threads
        # only orchestrate those children and avoid platform semaphore probes
        # that are unavailable in the restricted replay environment.
        executor = ThreadPoolExecutor(max_workers=workers)
        iterator = executor.map(certify_task, tasks, chunksize=1)
    try:
        for found in iterator:
            rows.append(found)
            if progress_every and len(rows) % progress_every == 0:
                print(
                    json.dumps(
                        {
                            "progress": len(rows),
                            "total": total,
                            "elapsed_seconds": round(time.monotonic() - started, 3),
                        },
                        sort_keys=True,
                    ),
                    file=sys.stderr,
                    flush=True,
                )
    finally:
        if executor is not None:
            executor.shutdown()

    if tranche_only:
        return {
            "schema": SCHEMA + ".preflight_tranche",
            "status": "EXACT_PREDETERMINED_ONE_CRITICAL_TRANCHE_NO_GO",
            "condition": "17<=a_lambda<=20",
            "records_per_orientation": len(black_records),
            "orientations": list(ORIENTATIONS),
            "rows": rows,
            "prime_histogram": histogram(rows),
            "orientation_prime_histograms": orientation_histograms(rows),
            "broader_claim_made": False,
        }

    return {
        "schema": SCHEMA,
        "status": "EXACT_PARTITION_33_CS_NO_GO_AND_SC_TRANCHE_NO_GO",
        "global_degree_29_status": "OPEN",
        "five_black_ww_status": "PARTIAL",
        "white_partition": "(3^2,1^23)",
        "endpoint_statuses": {
            "CS": "complete after adjoining the preceding a_lambda>=21 theorem",
            "SC": "only the frozen 17<=a_lambda<=20 tranche",
        },
        "reflection_symmetry_used": False,
        "marked_black_degree_condition": "a_lambda<=20",
        "preceding_all_status_theorem_condition": "a_lambda>=21",
        "predetermined_tranche": {
            "condition": "17<=a_lambda<=20",
            "records_per_orientation": PREDETERMINED_TRANCHE_RECORDS_PER_ORIENTATION,
            "passed_before_full_run": True,
        },
        "CS_black_degree_records_certified_here": len(black_records),
        "CS_preceding_high_a_records": 12,
        "SC_black_degree_records_certified_here": len(sc_tranche),
        "SC_black_degree_records_remaining": 1239,
        "one_critical_marked_passports_certified": len(rows),
        "new_marked_passports_eliminated": len(rows),
        "critical_value_quotient_dimension_every_row": (
            CRITICAL_VALUE_QUOTIENT_DIMENSION
        ),
        "partition_black_degree_systems_remaining": 1280,
        "partition_marked_passports_before": 3840,
        "partition_marked_passports_remaining": 2519,
        "partition_remaining_endpoint_status_counts": {"SS": 1280, "SC": 1239},
        "five_black_status_free_necessary_systems_remaining": 9742,
        "five_black_binary_status_master_systems_remaining": 22068,
        "five_black_ww_marked_passports_remaining": 27067,
        "five_black_marked_passports_remaining": 128467,
        "global_marked_passports_eliminated": 109923,
        "global_marked_passports_remaining": 42500755,
        "prime_histogram": histogram(rows),
        "orientation_prime_histograms": orientation_histograms(rows),
        "successive_exact_algebra": {
            "CS_critical_quadratic_factor": "K=x*(x-r)",
            "SC_critical_quadratic_factor": "K=(x-1)*(x-r)",
            "derivative_equations": "D=29*K^2",
            "derivative_base_dimension": 1,
            "critical_value_quotient_dimension": (
                CRITICAL_VALUE_QUOTIENT_DIMENSION
            ),
            "adjoined_simple_endpoint_equation_final_ideal": [1],
        },
        "projective_infinity_certificate": {
            "structural_degrees": list(STRUCTURAL_HOMOGENIZATION_DEGREES),
            "lemma": (
                "At h=0 the fixed marked black root and the fixed endpoints "
                "coalesce at 0. The critical-value equations force both roots "
                "of K to be black-root clusters. If r!=0, D_infinity has two "
                "order-two roots and forces disjoint clusters of three labels "
                "each, impossible with five labels. If r=0, its order-four "
                "root forces all five labels to collide at 0, the forbidden "
                "zero projective point."
            ),
            "valid_for_all_selected_primes": True,
        },
        "rows": rows,
        "complete_CS_no_go_claimed": True,
        "complete_SC_no_go_claimed": False,
        "complete_SC_tranche_no_go_claimed": True,
        "complete_partition_no_go_claimed": False,
        "complete_five_black_ww_no_go_claimed": False,
        "complete_degree_29_no_go_claimed": False,
        "degree_29_witness_found": False,
        "broad_plane_tree_search_used": False,
    }


def verify_archive(archive: dict[str, Any], tranche_only: bool = False) -> None:
    expected_rows = (
        PREDETERMINED_TRANCHE_RECORDS_PER_ORIENTATION * len(ORIENTATIONS)
        if tranche_only
        else EXPECTED_RECORDS_PER_ORIENTATION
        + PREDETERMINED_TRANCHE_RECORDS_PER_ORIENTATION
    )
    if len(archive["rows"]) != expected_rows:
        raise ArithmeticError(("one-critical row count", len(archive["rows"])))
    if sum(archive["prime_histogram"].values()) != expected_rows:
        raise ArithmeticError("prime histogram changed")
    expected_orientation_counts = (
        {"CS": 41, "SC": 41} if tranche_only else {"CS": 1280, "SC": 41}
    )
    for orientation, expected_count in expected_orientation_counts.items():
        count = sum(row["orientation"] == orientation for row in archive["rows"])
        if count != expected_count:
            raise ArithmeticError((orientation, count, expected_count))
    for row in archive["rows"]:
        if row["critical_value_quotient_dimension"] != 696:
            raise ArithmeticError(("critical-value quotient dimension", row))
    if tranche_only:
        if archive["broader_claim_made"]:
            raise ArithmeticError("tranche archive cannot make a broader claim")
        return

    checks = {
        "CS_black_degree_records_certified_here": 1280,
        "SC_black_degree_records_certified_here": 41,
        "one_critical_marked_passports_certified": 1321,
        "new_marked_passports_eliminated": 1321,
        "partition_marked_passports_remaining": 2519,
        "five_black_ww_marked_passports_remaining": 27067,
        "five_black_marked_passports_remaining": 128467,
        "global_marked_passports_eliminated": 109923,
        "global_marked_passports_remaining": 42500755,
    }
    for key, expected in checks.items():
        if archive[key] != expected:
            raise ArithmeticError((key, archive[key], expected))
    for key in (
        "complete_partition_no_go_claimed",
        "complete_five_black_ww_no_go_claimed",
        "complete_degree_29_no_go_claimed",
        "degree_29_witness_found",
        "broad_plane_tree_search_used",
        "reflection_symmetry_used",
        "complete_SC_no_go_claimed",
    ):
        if archive[key]:
            raise ArithmeticError((key, "must remain false"))
    if not archive["complete_CS_no_go_claimed"]:
        raise ArithmeticError("CS theorem must be complete")
    if not archive["complete_SC_tranche_no_go_claimed"]:
        raise ArithmeticError("the frozen SC tranche must be complete")


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--progress-every", type=int, default=10)
    parser.add_argument("--tranche-only", action="store_true")
    parser.add_argument("--workers", type=int, default=1)
    arguments = parser.parse_args()
    if not 1 <= arguments.workers <= 4:
        raise SystemExit("--workers must lie in [1,4] under the 0.6 GiB guard")
    archive = build_archive(
        arguments.progress_every, arguments.tranche_only, arguments.workers
    )
    verify_archive(archive, arguments.tranche_only)
    encoded = canonical_bytes(archive)
    digest = hashlib.sha256(encoded).hexdigest()

    if arguments.tranche_only:
        if arguments.check:
            raise SystemExit("--check is not defined for the transient tranche preflight")
    elif arguments.check:
        if not arguments.output.exists() or arguments.output.read_bytes() != encoded:
            raise SystemExit("partition-33 one-critical tower archive is stale")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_bytes(encoded)

    print(
        json.dumps(
            {
                "status": archive["status"],
                "rows_certified": len(archive["rows"]),
                "prime_histogram": archive["prime_histogram"],
                "sha256": digest,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
