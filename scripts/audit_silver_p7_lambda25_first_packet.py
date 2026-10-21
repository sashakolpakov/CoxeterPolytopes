#!/usr/bin/env python3
"""Freeze the exact first packet of the silver p7 Lambda25 host."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from collections import Counter
from pathlib import Path

from sympy import Matrix

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_lambda25 import (
    COMPLETE_P7_LAMBDA25_NORM_PROOF,
    FIRST_PACKET_CONTROL,
    FIRST_PACKET_HEIGHT,
    FIRST_PACKET_NORM,
    FIRST_PACKET_SHELL_NORM,
    FIRST_PACKET_TARGET,
    FIRST_PACKET_WALL_COUNT,
    FULL_MINIMAL_SHELL_COUNT,
    HEIGHT_ZERO_VECTOR,
    LAMBDA25_SOURCE_URL,
    P7_GENERATOR,
    P7_LAMBDA25_COEFFICIENT,
    P7_LAMBDA25_NORM_CLASSES,
    PAIRING_RADICAL_MINIMAL_COUNT,
    PAIRING_RADICAL_MINIMUM,
    enumerate_p7_lambda25_controls,
    first_packet_pairing_from_dot,
    first_packet_reflection_integral,
    gated_gram_triples,
    height_zero_gram_triples,
    height_zero_vector_audit,
    lambda25_catalogue_audit,
    lambda25_gram,
    lambda25_pairing_radical_basis,
    lambda25_pairing_radical_gram,
    lambda25_smith_invariants,
    reduced_spatial_target,
    strict_witness_first_packet_pairing,
    strict_witness_height_zero_pairing,
    strict_witness_norm,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "data/enumerations/silver_p7_lambda25_first_packet.json"
DEFAULT_REPORT = ROOT / "results/silver_p7_lambda25_first_packet.md"


def value_json(value: Silver) -> dict[str, object]:
    return {
        "exact": value.to_json(),
        "display": str(value),
        "signs": [value.sign(0), value.sign(1)],
    }


def _gp_matrix(matrix: tuple[tuple[int, ...], ...]) -> str:
    return "[" + ";".join(",".join(map(str, row)) for row in matrix) + "]"


def pari_shell_audit() -> dict[str, object]:
    """Replay all shell counts and the decisive orthogonality in PARI/GP."""

    gp = shutil.which("gp")
    if gp is None:
        raise RuntimeError("PARI/GP is required to regenerate this certificate")
    gram = lambda25_gram()
    pieces = [f"G={_gp_matrix(gram)}"]
    for modulus in (1, 2, 4, 8, 16):
        basis = lambda25_pairing_radical_basis(modulus)
        pieces.extend(
            (
                f"B{modulus}={_gp_matrix(basis)}",
                f"Q{modulus}=B{modulus}~*G*B{modulus}",
                f"U{modulus}=qflllgram(Q{modulus})",
                f"S{modulus}=U{modulus}~*Q{modulus}*U{modulus}",
                f"R{modulus}=qfminim(S{modulus},,1)",
                f'print("radical_{modulus}_signed_count=",R{modulus}[1])',
                f'print("radical_{modulus}_minimum=",floor(R{modulus}[2]+1/2))',
            )
        )
    pieces.extend(
        (
            # R1 already gives the complete signed minimal-shell count.  Do
            # not materialise the irrelevant 98,328-column full-shell matrix.
            'print("full_signed_count=",R1[1])',
            'print("full_minimum=",floor(R1[2]+1/2))',
            "T=qfminim(S2,8,,2)",
            "W=B2*U2*T[3]",
            "nt=matsize(W)[2]",
            "MT=vector(25,i,vecmax(vector(nt,j,abs(W[i,j]))))",
            'print("s2_signed_count=",T[1])',
            'print("s2_representatives_mod_sign=",nt)',
            'print("s2_coordinate_maxima=",MT)',
            "Z=qfminim(S8,16,,2)",
            "W0=B8*U8*Z[3][,1]",
            "if(W0[25]>0,W0=-W0)",
            'print("height_zero_vector=",Vec(W0))',
            'print("height_zero_square=",W0~*G*W0)',
            'print("height_zero_gram_image=",Vec(G*W0))',
            "D=W0~*G*W",
            'print("s2_nonorthogonal_count_mod_sign=",sum(j=1,nt,D[j]!=0))',
            'print("s2_last_coordinate_nonzero_count_mod_sign=",sum(j=1,nt,W[25,j]!=0))',
            "quit()",
        )
    )
    process = subprocess.run(
        [gp, "-q", "-s", "4000000000"],
        input=";".join(pieces) + "\n",
        text=True,
        capture_output=True,
        timeout=300,
        check=True,
    )
    parsed: dict[str, object] = {}
    for raw in process.stdout.splitlines():
        key, text = raw.strip().split("=", 1)
        parsed[key] = json.loads(text)

    expected = {
        "radical_1_signed_count": 196_656,
        "radical_1_minimum": 8,
        "radical_2_signed_count": 98_256,
        "radical_2_minimum": 8,
        "radical_4_signed_count": 2,
        "radical_4_minimum": 16,
        "radical_8_signed_count": 2,
        "radical_8_minimum": 16,
        "radical_16_signed_count": 2,
        "radical_16_minimum": 64,
        "full_signed_count": 196_656,
        "full_minimum": 8,
        "s2_signed_count": 98_256,
        "s2_representatives_mod_sign": 49_128,
        "s2_coordinate_maxima": [
            14, 4, 5, 4, 5, 2, 2, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0,
        ],
        "height_zero_vector": list(HEIGHT_ZERO_VECTOR),
        "height_zero_square": 16,
        "height_zero_gram_image": [0] * 24 + [-8],
        "s2_nonorthogonal_count_mod_sign": 0,
        "s2_last_coordinate_nonzero_count_mod_sign": 0,
    }
    if parsed != expected:
        raise AssertionError(f"the Lambda25 shell transcript drifted: {parsed!r}")
    return {
        "engine": "PARI/GP qflllgram and qfminim with exact integer replay",
        **parsed,
    }


def build_payload() -> dict[str, object]:
    if not lambda25_catalogue_audit():
        raise AssertionError("the catalogue Gram audit failed")
    if lambda25_smith_invariants() != (1, 1) + (2,) * 22 + (16,):
        raise AssertionError("the Lambda25 Smith form drifted")
    if not height_zero_vector_audit():
        raise AssertionError("the pinned height-zero vector drifted")

    height_zero = tuple(
        (row, height_zero_gram_triples(row))
        for row in P7_LAMBDA25_NORM_CLASSES
        if height_zero_gram_triples(row)
    )
    if len(height_zero) != 1:
        raise AssertionError("the unique height-zero component changed")
    if height_zero[0][0].norm != 16 or [
        (row.u_square, row.v_square, row.u_dot_v) for row in height_zero[0][1]
    ] != [(16, 0, 0)]:
        raise AssertionError("the height-zero Gram type changed")

    controls = enumerate_p7_lambda25_controls()
    surviving = tuple(row for row in controls if gated_gram_triples(row))
    if len(controls) != 17 or len({row.height for row in controls}) != 9:
        raise AssertionError("the prefix control census changed")
    if len(surviving) != 1:
        raise AssertionError("the first realized packet ceased to be unique")
    first = surviving[0]
    if (
        first.height != FIRST_PACKET_HEIGHT
        or first.norm_class.norm != FIRST_PACKET_NORM
        or first.control != FIRST_PACKET_CONTROL
        or first.target != FIRST_PACKET_TARGET
        or [
            (row.u_square, row.v_square, row.u_dot_v)
            for row in gated_gram_triples(first)
        ]
        != [(8, 8, 8)]
    ):
        raise AssertionError("the first packet changed")

    explicit_minimal = (0, 1) + (0,) * 23
    if not first_packet_reflection_integral(explicit_minimal):
        raise AssertionError("the explicit S2 minimal root ceased to reflect")
    if first_packet_pairing_from_dot(4).sign(0) >= 0:
        raise AssertionError("same-height acuteness bound failed")
    if strict_witness_norm().sign(0) <= 0:
        raise AssertionError("the hyperideal witness ceased to be spacelike")
    if strict_witness_first_packet_pairing().sign(0) >= 0:
        raise AssertionError("the hyperideal witness lost strict packet margin")

    shell = pari_shell_audit()
    height_counts = Counter(str(row.height) for row in controls)
    radical_indices = {
        str(modulus): abs(int(Matrix(lambda25_pairing_radical_basis(modulus)).det()))
        for modulus in (1, 2, 4, 8, 16)
    }
    return {
        "schema_version": "0.1.0",
        "seed_id": "silver-H25-Lambda25-a-2sqrt2-minus1",
        "field": "Q(sqrt(2))",
        "ring": "Z[sqrt(2)]",
        "hyperbolic_dimension": 25,
        "lattice": "<-(2sqrt(2)-1)> + Lambda25",
        "scope": (
            "Complete exact Vinberg prefix through the first positive packet "
            "H=5+4sqrt(2), including the complete height-zero packet. The "
            "strict spacelike witness proves only that this prefix has infinite "
            "volume; later mirrors may cut it."
        ),
        "timelike_coefficient": value_json(P7_LAMBDA25_COEFFICIENT),
        "prime_seven_generator": value_json(P7_GENERATOR),
        "lambda25_source": {
            "url": LAMBDA25_SOURCE_URL,
            "dimension": 25,
            "determinant": 2**26,
            "smith_invariants": list(lambda25_smith_invariants()),
            "minimum": 8,
            "kissing_number": FULL_MINIMAL_SHELL_COUNT,
        },
        "local_mirror_norms": {
            "class_count": len(P7_LAMBDA25_NORM_CLASSES),
            "proof": COMPLETE_P7_LAMBDA25_NORM_PROOF,
            "classes": [
                {
                    "norm": value_json(row.norm),
                    "dyadic_valuation": row.dyadic_valuation,
                    "p7_gate": row.p7_gate,
                    "coefficient_pairing_moduli": list(row.coefficient_pairing_moduli),
                }
                for row in P7_LAMBDA25_NORM_CLASSES
            ],
            "pairing_radicals": {
                str(modulus): {
                    "index": radical_indices[str(modulus)],
                    "minimum": PAIRING_RADICAL_MINIMUM[modulus],
                    "minimal_signed_count": PAIRING_RADICAL_MINIMAL_COUNT[modulus],
                }
                for modulus in (1, 2, 4, 8, 16)
            },
        },
        "height_zero": {
            "unique_norm": value_json(height_zero[0][0].norm),
            "unique_gram_type": [16, 0, 0],
            "root_vectors_mod_sign": 1,
            "wall_count": 1,
            "oriented_root": list(HEIGHT_ZERO_VECTOR),
            "gram_image": [0] * 24 + [-8],
            "square": 16,
        },
        "positive_prefix": {
            "height_bound": value_json(FIRST_PACKET_HEIGHT),
            "control_count": len(controls),
            "distinct_height_count": len(height_counts),
            "earlier_control_count": sum(row.height < FIRST_PACKET_HEIGHT for row in controls),
            "earlier_fully_gated_gram_type_count": 0,
            "controls": [
                {
                    "height": value_json(row.height),
                    "norm": value_json(row.norm_class.norm),
                    "control": value_json(row.control),
                    "target": value_json(row.target),
                    "reduced_target_after_p7_gate": value_json(reduced_spatial_target(row)),
                    "dyadic_valuation": row.norm_class.dyadic_valuation,
                    "p7_gate": row.norm_class.p7_gate,
                    "coefficient_pairing_moduli": list(
                        row.norm_class.coefficient_pairing_moduli
                    ),
                    "fully_gated_gram_triples": [
                        [triple.u_square, triple.v_square, triple.u_dot_v]
                        for triple in gated_gram_triples(row)
                    ],
                }
                for row in controls
            ],
        },
        "first_packet": {
            "height": value_json(FIRST_PACKET_HEIGHT),
            "norm": value_json(FIRST_PACKET_NORM),
            "control": value_json(FIRST_PACKET_CONTROL),
            "target": value_json(FIRST_PACKET_TARGET),
            "gram_triple": [8, 8, 8],
            "roots": (
                "r_b=(2+2sqrt(2),(1+sqrt(2))*b), where b^2=8 and "
                "G*b=0 mod 2"
            ),
            "normalized_walls": "(2,b)",
            "full_minimal_shell_count": FULL_MINIMAL_SHELL_COUNT,
            "reflecting_s2_shell_count": FIRST_PACKET_WALL_COUNT,
            "excluded_by_reflection_gate": FULL_MINIMAL_SHELL_COUNT
            - FIRST_PACKET_WALL_COUNT,
            "same_height_distinct_dot_upper": 4,
            "same_height_pairing_at_upper": value_json(
                first_packet_pairing_from_dot(4)
            ),
            "shell_audit": shell,
        },
        "strict_hyperideal_witness": {
            "vector": "z=(1,-w0)",
            "spatial_coordinates": [-coordinate for coordinate in HEIGHT_ZERO_VECTOR],
            "norm": value_json(strict_witness_norm()),
            "height_zero_pairing": value_json(strict_witness_height_zero_pairing()),
            "first_packet_pairing": value_json(strict_witness_first_packet_pairing()),
            "orthogonality_certificate": (
                "Every one of the 49,128 S2 minimal representatives has "
                "b_24=0, while G*w0=-8*e_24; hence (w0,b)=0."
            ),
            "strict": True,
        },
        "discarded_full_shell_polar": {
            "reason": (
                "The 196,656-vector full-shell contact polar ignores the exact "
                "q=4 reflection gate G*b=0 mod 2. Only 98,256 shell vectors "
                "are walls, so conclusions from the full-shell polar do not "
                "describe the Vinberg cone."
            )
        },
        "conclusion": "complete first positive packet prefix has infinite hyperbolic volume",
        "discovery": None,
    }


def report() -> str:
    return """# Silver p7 Lambda25: the exact first packet is infinite-volume

For the arithmetic host

`<-(2sqrt(2)-1)> + Lambda25`

in hyperbolic dimension 25, the complete prefix through the first positive
packet does **not** close a finite-volume chamber.

## Exact local and chronological screen

The catalogue Gram matrix has determinant `2^26`, Smith form
`1,1,2^22,16`, minimum 8, and 196,656 minimal vectors.  Its complete
primitive mirror alphabet has the ten dyadic valuations
`0,2,3,4,5,6,7,8,9,10`, with an independent optional factor `3+sqrt(2)`
above seven: twenty classes total.  Exact pairing-radical minima for moduli
`1,2,4,8,16` are respectively `8,8,16,16,64`.

There is exactly one height-zero wall.  With the chosen orientation it is
represented by

```
w0=(-2,1,0,1,2,1,1,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2),
```

and exact arithmetic gives `w0^2=16` and `G*w0=-8*e_24`.

There are 17 arithmetic controls at nine heights through
`H=5+4sqrt(2)`.  All 14 earlier controls are empty.  At the final height the
norm-one and norm-two controls are empty; the unique realized component is

`r_b=(2+2sqrt(2),(1+sqrt(2))*b)`

with `b^2=8` and `G*b=0 mod 2`.

## Why the full-shell compactness signal was false

The norm-four reflection gate is essential.  The full Lambda25 minimal shell
has 196,656 vectors, but PARI `qfminim` finds only **98,256** minimal vectors
in `S2={b:G*b=0 mod 2}`.  Therefore the first packet has 98,256 walls, not
196,656.  A contact polar computed from the full shell is not a polar of the
Vinberg cone.

Moreover, every reflecting minimal vector has `b_24=0`.  Since
`G*w0=-8*e_24`, all 98,256 first-packet walls are exactly orthogonal to the
height-zero wall.

## Strict hyperideal certificate

Choose the height-zero side `(z,w0)<=0` and put `z=(1,-w0)`.  Then

`z^2=17-2sqrt(2)>0`,

its pairing with the height-zero wall is `-16`, and its pairing with every
first-packet root is `-6-2sqrt(2)<0`.  Thus a strict open hyperideal patch
survives the complete prefix, proving that prefix has infinite hyperbolic
volume.  Later mirrors may still cut it; this is not a global
nonreflectivity result.

## Reproduction

```sh
PYTHONPATH=src python3 scripts/audit_silver_p7_lambda25_first_packet.py
PYTHONPATH=src pytest -q tests/test_silver_p7_lambda25_first_packet.py
```

Source: https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/LAMBDA25.html
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()
    payload = build_payload()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(report())
    print(args.output)
    print(args.report)


if __name__ == "__main__":
    main()
