#!/usr/bin/env python3
"""Freeze the exact first packet of the silver p7 Lambda22 host."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from collections import Counter
from pathlib import Path

from coxeter_polytopes.exact.silver import SQRT2, Silver
from coxeter_polytopes.search.cubic7_lambda22 import (
    PARI_QFMINIM_CERTIFICATES,
    lambda22_gram,
)
from coxeter_polytopes.search.silver_p7_lambda22 import (
    COMPLETE_P7_LAMBDA22_NORM_PROOF,
    FIRST_PACKET_CONTROL,
    FIRST_PACKET_HEIGHT,
    FIRST_PACKET_NORM,
    FIRST_PACKET_SHELL_NORM,
    FIRST_PACKET_TARGET,
    FIRST_PACKET_WALL_COUNT,
    GRAM_FILE_SHA256,
    GRAM_STREAM_SHA256,
    P7_DUAL_DIRECTION,
    P7_GENERATOR,
    P7_LAMBDA22_COEFFICIENT,
    P7_LAMBDA22_NORM_CLASSES,
    PAIRING_RADICAL_MINIMUM,
    SOURCE_HTML_SHA256,
    SOURCE_URL,
    ambient_gram_triples,
    dual_direction_audit,
    enumerate_p7_lambda22_controls,
    first_packet_pairing_from_dot,
    first_packet_reflection_integral,
    gated_gram_triples,
    height_zero_gram_triples,
    lambda22_basis_audit,
    lambda22_gram_file_sha256,
    lambda22_gram_stream_sha256,
    lambda22_smith_invariants,
    p7_lambda22_factorization_audit,
    reduced_spatial_target,
    strict_witness_norm,
    strict_witness_pairing_from_last_coordinate,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "data/enumerations/silver_p7_lambda22_first_packet.json"
DEFAULT_REPORT = ROOT / "results/silver_p7_lambda22_first_packet.md"


def value_json(value: Silver) -> dict[str, object]:
    return {
        "exact": value.to_json(),
        "display": str(value),
        "signs": [value.sign(0), value.sign(1)],
    }


def pari_minimal_shell_audit() -> dict[str, object]:
    """Replay the complete minimal shell and its final-coordinate spectrum."""

    gp = shutil.which("gp")
    if gp is None:
        raise RuntimeError("PARI/GP is required to regenerate this certificate")
    gram = lambda22_gram()
    matrix = "[" + ";".join(",".join(map(str, row)) for row in gram) + "]"
    program = ";".join(
        (
            f"G={matrix}",
            "R=qfminim(G,4,,2)",
            "V=R[3]",
            "m=matsize(V)[2]",
            'print("signed_count=",R[1])',
            'print("maximum=",floor(R[2]+1/2))',
            'print("representatives_mod_sign=",m)',
            'print("coordinate_abs_values=",Set(vector(m,j,abs(V[22,j]))))',
            'print("coordinate_abs_0=",sum(j=1,m,abs(V[22,j])==0))',
            'print("coordinate_abs_1=",sum(j=1,m,abs(V[22,j])==1))',
            'print("coordinate_abs_2=",sum(j=1,m,abs(V[22,j])==2))',
            "quit()",
        )
    )
    process = subprocess.run(
        [gp, "-q", "-s", "512000000"],
        input=program + "\n",
        text=True,
        capture_output=True,
        timeout=180,
        check=True,
    )
    parsed: dict[str, object] = {}
    for raw in process.stdout.splitlines():
        key, text = raw.strip().split("=", 1)
        parsed[key] = json.loads(text)
    expected = {
        "signed_count": 49_896,
        "maximum": 4,
        "representatives_mod_sign": 24_948,
        "coordinate_abs_values": [0, 1, 2],
        "coordinate_abs_0": 13_860,
        "coordinate_abs_1": 10_752,
        "coordinate_abs_2": 336,
    }
    if parsed != expected:
        raise AssertionError(f"the Lambda22 minimal-shell census drifted: {parsed!r}")
    ordinary = PARI_QFMINIM_CERTIFICATES[1]
    return {
        "algorithm": "PARI qfminim(G,4,,flag=2)",
        **parsed,
        "representative_stream_sha256": ordinary[
            "ordinary_vector_stream_sha256"
        ],
        "coordinate_note": (
            "Counts by |d_22| are for the one-representative-per-sign matrix; "
            "the absolute spectrum is sign invariant."
        ),
    }


def build_payload() -> dict[str, object]:
    if not lambda22_basis_audit():
        raise AssertionError("the authoritative Lambda22 basis audit failed")
    if not p7_lambda22_factorization_audit():
        raise AssertionError("the silver prime-seven factorization audit failed")
    if not dual_direction_audit():
        raise AssertionError("the pinned dual direction lost its exact identities")
    if lambda22_smith_invariants() != (1,) * 20 + (2, 6):
        raise AssertionError("the Lambda22 Smith form drifted")

    controls = enumerate_p7_lambda22_controls()
    if len(controls) != 16 or len({row.height for row in controls}) != 10:
        raise AssertionError("the first-packet control census drifted")
    first_rows = tuple(row for row in controls if row.height == FIRST_PACKET_HEIGHT)
    if tuple(row.norm_class.norm for row in first_rows) != (
        Silver(1),
        Silver(2),
        Silver(4),
    ):
        raise AssertionError("the first simultaneous packet changed")
    surviving = tuple(row for row in controls if gated_gram_triples(row))
    if len(surviving) != 1 or surviving[0].norm_class.norm != FIRST_PACKET_NORM:
        raise AssertionError("the first realized packet changed")
    if gated_gram_triples(surviving[0]) != ((4, 4, 4),):
        raise AssertionError("the first packet ceased to be the minimal shell")
    if any(height_zero_gram_triples(row) for row in P7_LAMBDA22_NORM_CLASSES):
        raise AssertionError("a height-zero reflecting wall appeared")

    minimal = (1,) + (0,) * 21
    if not first_packet_reflection_integral(minimal):
        raise AssertionError("the explicit first-packet representative ceased to reflect")
    shell = pari_minimal_shell_audit()
    if shell["signed_count"] != FIRST_PACKET_WALL_COUNT:
        raise AssertionError("the first-packet wall count drifted")
    witness_norm = strict_witness_norm()
    witness_maximum = strict_witness_pairing_from_last_coordinate(2)
    if witness_norm.sign(0) <= 0 or witness_maximum.sign(0) >= 0:
        raise AssertionError("the strict hyperideal witness ceased to be strict")

    height_counts = Counter(str(row.height) for row in controls)
    return {
        "schema_version": "0.1.0",
        "seed_id": "silver-H22-Lambda22-a-2sqrt2-minus1",
        "field": "Q(sqrt(2))",
        "ring": "Z[sqrt(2)]",
        "hyperbolic_dimension": 22,
        "lattice": "<-(2sqrt(2)-1)> + Lambda22",
        "scope": (
            "Complete exact Vinberg prefix through the first nonempty packet "
            "H=5+4sqrt(2). The strict spacelike witness proves only that this "
            "bounded prefix has infinite volume; later mirrors may cut it. No "
            "eventual reflectivity, finite-volume, compactness, or global "
            "nonreflectivity conclusion is made."
        ),
        "timelike_coefficient": value_json(P7_LAMBDA22_COEFFICIENT),
        "prime_seven_generator": value_json(P7_GENERATOR),
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
            "class_count": len(P7_LAMBDA22_NORM_CLASSES),
            "proof": COMPLETE_P7_LAMBDA22_NORM_PROOF,
            "classes": [
                {
                    "norm": value_json(row.norm),
                    "dyadic_valuation": row.dyadic_valuation,
                    "ternary_gate": row.ternary_gate,
                    "p7_gate": row.p7_gate,
                    "coefficient_pairing_moduli": list(
                        row.coefficient_pairing_moduli
                    ),
                    "coefficient_minima": [
                        PAIRING_RADICAL_MINIMUM[modulus]
                        for modulus in row.coefficient_pairing_moduli
                    ],
                }
                for row in P7_LAMBDA22_NORM_CLASSES
            ],
        },
        "height_zero": {
            "class_count": len(P7_LAMBDA22_NORM_CLASSES),
            "fully_gated_gram_type_count": 0,
            "reflecting_wall_count": 0,
            "statement": (
                "All sixteen local classes are empty after the exact p7 target "
                "division and the Lambda22 pairing-radical minima."
            ),
        },
        "positive_prefix": {
            "height_bound": value_json(FIRST_PACKET_HEIGHT),
            "control_count": len(controls),
            "distinct_height_count": len(height_counts),
            "earlier_control_count": sum(
                row.height < FIRST_PACKET_HEIGHT for row in controls
            ),
            "earlier_fully_gated_gram_type_count": 0,
            "controls": [
                {
                    "height": value_json(row.height),
                    "norm": value_json(row.norm_class.norm),
                    "control": value_json(row.control),
                    "target": value_json(row.target),
                    "reduced_target_after_p7_gate": value_json(
                        reduced_spatial_target(row)
                    ),
                    "dyadic_valuation": row.norm_class.dyadic_valuation,
                    "ternary_gate": row.norm_class.ternary_gate,
                    "p7_gate": row.norm_class.p7_gate,
                    "ambient_gram_triples": [
                        list(triple) for triple in ambient_gram_triples(row)
                    ],
                    "fully_gated_gram_triples": [
                        list(triple) for triple in gated_gram_triples(row)
                    ],
                }
                for row in controls
            ],
        },
        "first_packet": {
            "height": value_json(FIRST_PACKET_HEIGHT),
            "simultaneous_control_norms": [
                value_json(row.norm_class.norm) for row in first_rows
            ],
            "unique_primitive_component": {
                "norm": value_json(FIRST_PACKET_NORM),
                "control": value_json(FIRST_PACKET_CONTROL),
                "target": value_json(FIRST_PACKET_TARGET),
                "gram_triple": [4, 4, 4],
                "roots": (
                    "r_d=(2+sqrt(2),(1+sqrt(2))*d), over every signed "
                    "minimal d in Lambda22"
                ),
                "shell_norm": FIRST_PACKET_SHELL_NORM,
                "wall_count": FIRST_PACKET_WALL_COUNT,
            },
            "other_simultaneous_controls": {
                "norm_1": "no ambient minimum-four Gram type",
                "norm_4": (
                    "ambient types (16,4,8),(8,8,8) are killed by the "
                    "dyadic gate minimum min(S2)=12"
                ),
            },
            "minimal_shell_audit": shell,
            "same_height_acuteness": {
                "distinct_minimal_dot_upper": 2,
                "pairing_at_upper": value_json(first_packet_pairing_from_dot(2)),
                "proof": (
                    "For distinct d,e of norm four, (d-e)^2>=4 gives d.e<=2."
                ),
            },
        },
        "strict_hyperideal_witness": {
            "vector": "z=(1,(4-sqrt(2))*k/8)",
            "dual_integer_vector_k": list(P7_DUAL_DIRECTION),
            "dual_identities": "G*k=3*e_22 and k^2=24",
            "shell_coordinate_bound": "|d_22|<=2 for all 49,896 minimal d",
            "norm": value_json(witness_norm),
            "maximum_wall_pairing": value_json(witness_maximum),
            "strict": True,
            "proof": (
                "Since (k,d)=3*d_22 and the exact shell census gives "
                "|d_22|<=2, the largest pairing is attained at d_22=2 and "
                "equals -(2+3sqrt(2))/4. The witness norm is "
                "31/4-5sqrt(2)>0."
            ),
        },
        "conclusion": "complete first-packet prefix has infinite hyperbolic volume",
        "discovery": None,
    }


def report(payload: dict[str, object]) -> str:
    return """# Silver p7 Lambda22: exact first packet is infinite-volume

For the arithmetic host

`<-(2sqrt(2)-1)> + Lambda22`

in hyperbolic dimension 22, the complete first nonempty Vinberg packet does
**not** close a finite-volume chamber.  An exact strict hyperideal witness
survives every wall in the packet.

## Complete local and chronological screen

The Smith form of the authoritative Nebe--Sloane `Lambda22` Gram matrix is
`1^20,2,6`.  Combining its dyadic valuations `0,2,3,4`, the optional inert
prime-three factor, and the optional prime-seven factor `3+sqrt(2)` gives
exactly 16 primitive mirror-norm classes.  For a prime-seven class the
spatial vector is first divided by `3+sqrt(2)`; the dyadic and ternary
pairing gates are then applied to the quotient.

There are 16 reflective controls at 10 distinct heights through
`H=5+4sqrt(2)`.  All 13 earlier controls are empty after exact target
division, minimum-four Gram enumeration, and the pairing-radical minima
`min(S2)=12`, `min(S3)=24`, `min(S6)=96`.  There are no height-zero walls.

At `H=5+4sqrt(2)`, the simultaneous norm-1 and norm-4 controls are empty.
The unique primitive component is

`r_d=(2+sqrt(2),(1+sqrt(2))*d),  d^2=4 in Lambda22`,

of norm 2.  PARI `qfminim` enumerates the full signed minimal shell, so this
packet contains exactly **49,896 walls**.  Distinct minimal vectors satisfy
`d.e<=2`, which makes all distinct same-height wall pairings strictly
negative.

## Strict hyperideal certificate

Let

```
k=(-4490,27287,17909,13568,8412,18464,12340,6324,3634,1942,
   2575,1337,544,272,408,214,32,64,38,22,12,8).
```

For the vendored Gram matrix, exact integer arithmetic gives
`G*k=3*e_22` and `k^2=24`.  The complete minimal-shell enumeration gives
`|d_22|<=2`.  Therefore

`z=(1,(4-sqrt(2))*k/8)`

has

`z^2=31/4-5sqrt(2)>0`,

while its largest pairing with any of the 49,896 walls is

`-(2+3sqrt(2))/4<0`.

Thus the complete first-packet prefix has infinite hyperbolic volume.  This
is not an eventual nonreflectivity result: a later mirror may cut the
displayed witness, and no global compactness or reflectivity claim is made.

## Reproduction

```sh
PYTHONPATH=src python3 scripts/audit_silver_p7_lambda22_first_packet.py
PYTHONPATH=src pytest -q tests/test_silver_p7_lambda22_first_packet.py
```

Source: https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/LAMBDA22.html
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
    args.report.write_text(report(payload))
    print(args.output)
    print(args.report)


if __name__ == "__main__":
    main()
