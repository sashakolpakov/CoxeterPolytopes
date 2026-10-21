#!/usr/bin/env python3
"""Freeze the exact H14 first-packet prefix for silver Lambda22."""

from __future__ import annotations

import json
import shutil
import subprocess
from fractions import Fraction
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_lambda22 import (
    FIRST_PACKET_CONTROL,
    FIRST_PACKET_DYADIC_SHELL_NORM,
    FIRST_PACKET_HEIGHT,
    FIRST_PACKET_NORM,
    FIRST_PACKET_SHELL_NORM,
    FIRST_PACKET_TARGET,
    GRAM_FILE_SHA256,
    GRAM_STREAM_SHA256,
    LAMBDA22_COEFFICIENT,
    LAMBDA22_MIRROR_NORMS,
    LAMBDA22_NORM10_SHELL_COUNT,
    LAMBDA22_NORM_GATE_PROOF,
    S2_MINIMUM,
    S2_NORM20_SHELL_COUNT,
    S3_MINIMUM,
    S6_MINIMUM,
    SOURCE_HTML_SHA256,
    SOURCE_URL,
    cross_packet_pairing_from_dot,
    dyadic_packet_pairing_from_dot,
    dyadic_witness_pairing_from_dot,
    enumerate_lambda22_controls,
    first_packet_pairing_from_dot,
    lambda22_basis_audit,
    lambda22_gram,
    lambda22_gram_file_sha256,
    lambda22_gram_stream_sha256,
    lambda22_pairing_radical_gram,
    lambda22_smith_invariants,
    necessary_components_through_first_packet,
    witness_norm,
    witness_pairing_from_dot,
)


ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = ROOT / "data/enumerations/silver_lambda22_h14.json"
REPORT_PATH = ROOT / "results/silver_lambda22_h14.md"


def value_json(value: Silver) -> dict[str, object]:
    return {"exact": value.to_json(), "display": str(value), "signs": [value.sign(0), value.sign(1)]}


def pari_sublattice_audit() -> dict[str, object]:
    gp = shutil.which("gp")
    if gp is None:
        raise RuntimeError("PARI/GP is required")
    rows = []
    for modulus, bound in ((1, 4), (2, 20), (3, 24), (6, 96)):
        if modulus == 2:
            # The full flag-2 matrix has 713,460 columns.  Its cumulative
            # exact result and shell differences are independently frozen;
            # materialising it in every routine audit costs about 2 GB.
            rows.append(
                {
                    "pairing_modulus": 2,
                    "bound": 20,
                    "signed_vector_count_at_or_below_bound": 1_426_920,
                    "maximum_norm_returned": 20,
                    "stored_representatives_mod_sign_shape": [22, 713_460],
                    "exact_shell_histogram": {"12": 8_448, "16": 49_896, "20": 1_368_576},
                    "replay_note": "full PARI materialization is pinned but skipped in the routine generator",
                }
            )
            continue
        gram = lambda22_pairing_radical_gram(modulus)
        matrix = "[" + ";".join(",".join(map(str, row)) for row in gram) + "]"
        program = (
            f"G={matrix};R=qfminim(G,{bound},,2);"
            "print(R[1]);print(R[2]);print(matsize(R[3]));quit()\n"
        )
        process = subprocess.run(
            [gp, "-q", "-s", "256000000"],
            input=program,
            text=True,
            capture_output=True,
            timeout=180,
            check=True,
        )
        lines = [line.strip() for line in process.stdout.splitlines() if line.strip()]
        if len(lines) != 3:
            raise AssertionError(process.stdout)
        shape = [int(item.strip()) for item in lines[2].strip("[]").split(",")]
        rows.append(
            {
                "pairing_modulus": modulus,
                "bound": bound,
                "signed_vector_count_at_or_below_bound": int(lines[0]),
                "maximum_norm_returned": int(lines[1].split(".", 1)[0]),
                "stored_representatives_mod_sign_shape": shape,
            }
        )
    expected = [
        (1, 49_896, 4),
        (2, 1_426_920, 20),
        (3, 1_782, 24),
        (6, 1_782, 96),
    ]
    if [
        (row["pairing_modulus"], row["signed_vector_count_at_or_below_bound"], row["maximum_norm_returned"])
        for row in rows
    ] != expected:
        raise AssertionError(f"PARI sublattice census drifted: {rows!r}")
    return {"algorithm": "PARI qfminim(HNF-gated Gram,bound,,flag=2)", "rows": rows}


def build_payload() -> dict[str, object]:
    if not lambda22_basis_audit():
        raise AssertionError("the vendored Lambda22 basis audit failed")
    controls = enumerate_lambda22_controls()
    components = necessary_components_through_first_packet()
    if len(controls) != 25:
        raise AssertionError("the H14 control count drifted")
    if lambda22_smith_invariants() != (1,) * 20 + (2, 6):
        raise AssertionError("the Lambda22 Smith form drifted")
    pari = pari_sublattice_audit()
    nonempty = [component for component in components if component.triples]
    if [component.level.height for component in nonempty] != [
        Silver(2, Fraction(3, 2)),
        Silver(2, 2),
        Silver(6, Fraction(9, 2)),
        Silver(7, 5),
        Silver(12, Fraction(17, 2)),
        Silver(12, 9),
        FIRST_PACKET_HEIGHT,
        FIRST_PACKET_HEIGHT,
    ]:
        raise AssertionError("the H14 necessary-component heights drifted")
    d = (1,) + (0,) * 21
    if witness_norm(d) != Silver(13, 0) / 9 - Silver(0, 1):
        raise AssertionError("the hyperideal witness norm drifted")
    if witness_pairing_from_dot(5).sign(0) >= 0:
        raise AssertionError("the ordinary packet no longer leaves a strict witness")
    if dyadic_witness_pairing_from_dot(8).sign(0) >= 0:
        raise AssertionError("the dyadic packet no longer leaves a strict witness")

    return {
        "schema_version": "0.1.0",
        "seed_id": "silver-H22-Lambda22-a-sqrt2-1",
        "scope": (
            "Complete exact Vinberg prefix through H=14+10sqrt(2). The strict "
            "hyperideal witness proves only that this prefix has infinite volume; "
            "later mirrors may cut it."
        ),
        "field": "Q(sqrt(2))",
        "ring": "Z[sqrt(2)]",
        "hyperbolic_dimension": 22,
        "lattice": "<-(sqrt(2)-1)> + Lambda22",
        "timelike_coefficient": value_json(LAMBDA22_COEFFICIENT),
        "lambda22_source": {
            "url": SOURCE_URL,
            "source_html_sha256": SOURCE_HTML_SHA256,
            "gram_file_sha256": lambda22_gram_file_sha256(),
            "expected_gram_file_sha256": GRAM_FILE_SHA256,
            "gram_stream_sha256": lambda22_gram_stream_sha256(),
            "expected_gram_stream_sha256": GRAM_STREAM_SHA256,
            "dimension": 22,
            "determinant": 12,
            "smith_invariants": list(lambda22_smith_invariants()),
            "minimum": 4,
            "kissing_number": 49_896,
        },
        "local_mirror_norms": {
            "canonical_norms": [value_json(norm) for norm in LAMBDA22_MIRROR_NORMS],
            "dyadic_valuations": [0, 2, 3, 4],
            "ternary_valuations": [0, 1],
            "proof": LAMBDA22_NORM_GATE_PROOF,
            "gates": {
                "v_pi_0_or_2": "no spatial pairing gate",
                "v_pi_3": "G*u even for y=u+sqrt(2)v",
                "v_pi_4": "G*u and G*v even",
                "v_3_1": "G*u and G*v divisible by 3",
            },
        },
        "gated_sublattices": {
            "S2_minimum": S2_MINIMUM,
            "S3_minimum": S3_MINIMUM,
            "S6_minimum": S6_MINIMUM,
            "pari": pari,
        },
        "height_bound": value_json(FIRST_PACKET_HEIGHT),
        "control_level_count": len(controls),
        "controls": [
            {
                "height": value_json(row.height),
                "norm": value_json(row.norm),
                "control": value_json(row.control),
                "target": value_json(row.target),
                "dyadic_valuation": row.dyadic_valuation,
                "ternary_gate": row.ternary_gate,
            }
            for row in controls
        ],
        "necessary_components": [
            {
                "height": value_json(row.level.height),
                "norm": value_json(row.level.norm),
                "control": value_json(row.level.control),
                "target": value_json(row.level.target),
                "dyadic_valuation": row.level.dyadic_valuation,
                "ternary_gate": row.level.ternary_gate,
                "gram_triples": [list(triple) for triple in row.triples],
                "resolution": row.exclusion,
            }
            for row in components
            if row.triples
        ],
        "first_packet": {
            "height": value_json(FIRST_PACKET_HEIGHT),
            "ordinary_component": {
                "norm": value_json(FIRST_PACKET_NORM),
                "control": value_json(FIRST_PACKET_CONTROL),
                "target": value_json(FIRST_PACKET_TARGET),
                "roots": "(6+4sqrt(2),(1+sqrt(2))*w), w in full Lambda22 norm-10 shell",
                "shell_norm": FIRST_PACKET_SHELL_NORM,
                "shell_count": LAMBDA22_NORM10_SHELL_COUNT,
                "count_proof": "exact PARI qflllgram + count-only qfminim cumulative shell difference",
            },
            "dyadic_component": {
                "norm": value_json(Silver(4)),
                "control": value_json(Silver(8, 6)),
                "target": value_json(Silver(60, 40)),
                "roots": "(8+6sqrt(2),(1+sqrt(2))*u), u in full S2 norm-20 shell",
                "shell_norm": FIRST_PACKET_DYADIC_SHELL_NORM,
                "shell_count": S2_NORM20_SHELL_COUNT,
                "primitivity": "automatic: norm 20 is 4 mod 8, hence u is in a nonzero mod-2 radical class",
            },
            "total_root_count": LAMBDA22_NORM10_SHELL_COUNT + S2_NORM20_SHELL_COUNT,
            "same_height_acuteness": {
                "ordinary_distinct_dot_upper": 8,
                "ordinary_pairing_at_upper": value_json(first_packet_pairing_from_dot(8)),
                "dyadic_distinct_dot_upper": 14,
                "dyadic_pairing_at_upper": value_json(dyadic_packet_pairing_from_dot(14)),
                "cross_dot_upper": 13,
                "cross_pairing_at_upper": value_json(cross_packet_pairing_from_dot(13)),
                "cross_dot_14_excluded": "(u-w)^2=2 contradicts min(Lambda22)=4",
            },
        },
        "strict_hyperideal_witness": {
            "vector": "z=(1,d/3), d any Lambda22 minimal vector",
            "norm": value_json(witness_norm(d)),
            "ordinary_dot_bound": 5,
            "ordinary_max_pairing": value_json(witness_pairing_from_dot(5)),
            "dyadic_dot_bound": 8,
            "dyadic_max_pairing": value_json(dyadic_witness_pairing_from_dot(8)),
            "proof": (
                "For w^2=10,d^2=4, |d.w|<=5 because |d.w|>=6 makes one "
                "of (w+-d)^2 at most 2. For u^2=20,d^2=4, Cauchy gives "
                "|d.u|<=8. Both exact wall pairings are strictly negative."
            ),
        },
        "conclusion": "H14 prefix has infinite hyperbolic volume",
        "discovery": None,
        "next_target": (
            "Continue the same coefficient to the first unresolved possible "
            "tangency at H=36+51sqrt(2)/2 in the mixed norm 12+6sqrt(2) class."
        ),
    }


def report(payload: dict[str, object]) -> str:
    return f"""# Silver Lambda22 exact H14 prefix

The exact host

`<-(sqrt(2)-1)> + Lambda22`

in hyperbolic dimension 22 has a complete Vinberg prefix through
`H=14+10sqrt(2)`.  It is **not** finite-volume at this cutoff: for any
minimal `d in Lambda22`, the exact ray `z=(1,d/3)` has

`z^2 = 13/9-sqrt(2) > 0`

and satisfies every frozen wall strictly.  This is a prefix result, not a
global nonreflectivity theorem.

## Exact lattice and local gates

The vendored Nebe--Sloane Gram matrix has determinant 12, minimum 4,
kissing number 49,896, and Smith form `1^20,2,6`.  The complete primitive
mirror-norm alphabet is

`1, 2, 4+2sqrt(2), 4, 3, 6, 12+6sqrt(2), 12`.

The dyadic valuations are `0,2,3,4`; exponent 1 is excluded by evenness.
For `y=u+sqrt(2)v`, valuation 3 requires `G*u` even, valuation 4 requires
both `G*u,G*v` even, and the 3-primary classes require both divisible by 3.
Independent HNF+PARI checks give minima `min(S2)=12`, `min(S3)=24`, and
`min(S6)=96`.

## First simultaneous packet

At `H=14+10sqrt(2)` there are exactly two shell descriptions:

1. `q=2`, `x=6+4sqrt(2)`, and `y=(1+sqrt(2))w` for the complete ordinary
   `Lambda22` norm-10 shell.  Its exact cardinality is 467,596,800.
2. `q=4`, `x=8+6sqrt(2)`, and `y=(1+sqrt(2))u` for the complete norm-20
   shell of `S2={{u:G*u even}}`.  Its exact signed cardinality is 1,368,576.

The simultaneous packet therefore contains 468,965,376 roots.

All other necessary Gram types through the same height are killed by the
stated gate minima or mod-8 radical residues.  The two full components are
mutually Vinberg-acute.  In particular, cross dot 14 would give a forbidden
norm-2 vector, so the cross dot is at most 13.

For the witness, ordinary-shell pairings use `|d.w|<=5`; dyadic-shell
pairings use `|d.u|<=8`.  The respective worst exact pairings are
`-(1+sqrt(2))/3` and `-4/3+2sqrt(2)/3`, both negative.

## Frontier

The next sharp target is the first unresolved possible tangency at
`H=36+51sqrt(2)/2`, in the mixed norm `12+6sqrt(2)` class.  That continuation
is not part of this frozen H14 certificate.

## Reproduction

```sh
PYTHONPATH=src python3 scripts/audit_silver_lambda22_h14.py
PYTHONPATH=src pytest -q tests/test_silver_lambda22.py
```

Source: {SOURCE_URL}
"""


def main() -> None:
    payload = build_payload()
    JSON_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    JSON_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    REPORT_PATH.write_text(report(payload))
    print(JSON_PATH)
    print(REPORT_PATH)


if __name__ == "__main__":
    main()
