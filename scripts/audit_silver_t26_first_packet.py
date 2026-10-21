#!/usr/bin/env python3
"""Freeze the exact full-shell polar obstruction for silver p7 T26."""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

from sympy import Matrix
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_t26 import (
    COMPLETE_NORM_CLASS_PROOF,
    FIRST_PACKET_CONTROL,
    FIRST_PACKET_HEIGHT,
    FIRST_PACKET_NORM,
    FIRST_PACKET_TARGET,
    FIRST_PACKET_WALL_COUNT,
    NORM_CLASSES,
    P7_GENERATOR,
    POLAR_BOUNDARY_NORM,
    POLAR_STRICT_FIRST_PACKET_MARGIN,
    POLAR_STRICT_NORM,
    POLAR_STRICT_SCALE,
    POLAR_VERTEX,
    POLAR_VERTEX_CONTACT_COUNT,
    POLAR_VERTEX_SHELL_MAXIMUM,
    POLAR_VERTEX_SQUARE,
    POLAR_WALL_BOUND,
    RHO,
    T26_KISSING_NUMBER,
    T26_MINIMUM,
    T26_RADICAL_THREE_KISSING_NUMBER,
    T26_RADICAL_THREE_MINIMUM,
    T26_SMITH_INVARIANTS,
    T26_SOURCE_URL,
    enumerate_controls,
    gated_gram_triples,
    height_zero_gated_gram_triples,
    lorentz_pair,
    polar_strict_witness,
    polar_strict_witness_audit,
    polar_vertex_audit,
    reduced_target,
    t26_catalogue_audit,
    t26_gram,
    t26_radical_three_basis,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "data/enumerations/silver_t26_first_packet.json"
DEFAULT_REPORT = ROOT / "results/silver_t26_first_packet.md"
MINIMAL_REPRESENTATIVE_STREAM_SHA256 = (
    "69cdb97ab9d5b89e24db1b1fe94b83f5488e68f14e71e65cf04a8d991b8c4870"
)


def value_json(value: Silver) -> dict[str, object]:
    return {
        "exact": value.to_json(),
        "display": str(value),
        "signs": [value.sign(0), value.sign(1)],
    }


def _gp_matrix(matrix) -> str:
    return "[" + ";".join(",".join(map(str, row)) for row in matrix) + "]"


def pari_shell_audit() -> dict[str, object]:
    """Enumerate, hash, and maximize over every minimal sign pair."""

    gp = shutil.which("gp")
    if gp is None:
        raise RuntimeError("PARI/GP is required to regenerate this certificate")
    gram = _gp_matrix(t26_gram())
    radical_basis = _gp_matrix(t26_radical_three_basis())
    polar = "[" + ",".join(
        f"({value.numerator}/{value.denominator})" for value in POLAR_VERTEX
    ) + "]~"
    program = ";".join(
        (
            f"G={gram}",
            f"B3={radical_basis}",
            "Q3=B3~*G*B3",
            "U3=qflllgram(Q3)",
            "S3=U3~*Q3*U3",
            "R3=qfminim(S3,,1)",
            f"Y={polar}",
            "R=qfminim(G,4,,2)",
            "W=R[3]",
            "n=matsize(W)[2]",
            "M=vector(26,i,vecmax(vector(n,j,abs(W[i,j]))))",
            "mx=0;contacts=0",
            "for(j=1,n,d=abs(Y~*G*W[,j]);if(d>mx,mx=d;contacts=1,if(d==mx,contacts++)))",
            'print("SUMMARY|",R[1],"|",floor(R[2]+1/2),"|",n,"|",mx,"|",contacts)',
            'print("MAXIMA|",M)',
            'print("RADICAL3|",R3[1],"|",floor(R3[2]+1/2))',
            'for(j=1,n,print("V|",Vec(W[,j])))',
            "quit()",
        )
    )
    completed = subprocess.run(
        [gp, "-q", "-s", "4000000000"],
        input=program + "\n",
        text=True,
        capture_output=True,
        timeout=300,
        check=True,
    )
    lines = completed.stdout.splitlines()
    if len(lines) != 3 + T26_KISSING_NUMBER // 2:
        raise AssertionError(f"the PARI shell transcript was truncated: {len(lines)}")
    summary = lines[0].split("|")
    maxima = ast.literal_eval(lines[1].split("|", 1)[1])
    radical = lines[2].split("|")
    expected_summary = ["SUMMARY", "117936", "4", "58968", "1", "66"]
    expected_maxima = [
        6, 14, 5, 4, 7, 4, 5, 4, 4, 4, 3, 4, 4,
        4, 4, 3, 4, 4, 2, 2, 2, 2, 2, 4, 2, 2,
    ]
    if summary != expected_summary or maxima != expected_maxima:
        raise AssertionError(
            f"the T26 minimal shell changed: {summary!r}, {maxima!r}"
        )
    if radical != ["RADICAL3", "1638", "24"]:
        raise AssertionError(f"the T26 S3 shell changed: {radical!r}")

    digest = hashlib.sha256()
    for line in lines[3:]:
        if not line.startswith("V|"):
            raise AssertionError("malformed PARI minimal-vector stream")
        vector = ast.literal_eval(line[2:])
        if len(vector) != 26:
            raise AssertionError("a PARI minimal vector has the wrong dimension")
        digest.update((" ".join(map(str, vector)) + "\n").encode("ascii"))
    stream_sha = digest.hexdigest()
    if stream_sha != MINIMAL_REPRESENTATIVE_STREAM_SHA256:
        raise AssertionError(f"the T26 vector-stream digest changed: {stream_sha}")
    return {
        "engine": "PARI/GP qfminim(G,4,,flag=2), exact rational maximization",
        "signed_count": int(summary[1]),
        "minimum": int(summary[2]),
        "representatives_mod_sign": int(summary[3]),
        "polar_absolute_maximum": summary[4],
        "polar_contact_sign_pairs": int(summary[5]),
        "coordinate_absolute_maxima": maxima,
        "representative_stream_serialization": (
            "one PARI column per line, decimal coordinates separated by one space"
        ),
        "representative_stream_sha256": stream_sha,
        "radical_three_minimum": int(radical[2]),
        "radical_three_kissing_number": int(radical[1]),
    }


def build_payload() -> dict[str, object]:
    if not t26_catalogue_audit():
        raise AssertionError("the catalogue Gram audit failed")
    smith = tuple(
        abs(int(value))
        for value in smith_normal_form(Matrix(t26_gram()), domain=ZZ).diagonal()
    )
    if smith != T26_SMITH_INVARIANTS:
        raise AssertionError("the T26 Smith form drifted")
    if not polar_vertex_audit() or not polar_strict_witness_audit():
        raise AssertionError("the exact polar witness failed")

    controls = enumerate_controls()
    nonempty = tuple(row for row in controls if gated_gram_triples(row))
    if len(controls) != 10 or len({row.height for row in controls}) != 7:
        raise AssertionError("the first-packet control census changed")
    if len(nonempty) != 1:
        raise AssertionError("the first realized component changed")
    first = nonempty[0]
    if (
        first.height != FIRST_PACKET_HEIGHT
        or first.norm_class.norm != FIRST_PACKET_NORM
        or first.control != FIRST_PACKET_CONTROL
        or first.target != FIRST_PACKET_TARGET
        or [
            (item.u_square, item.v_square, item.u_dot_v)
            for item in gated_gram_triples(first)
        ]
        != [(4, 4, 4)]
    ):
        raise AssertionError("the first packet changed")
    if any(height_zero_gated_gram_triples(row) for row in NORM_CLASSES):
        raise AssertionError("a height-zero wall appeared")

    shell = pari_shell_audit()
    if shell["signed_count"] != FIRST_PACKET_WALL_COUNT:
        raise AssertionError("the first-packet wall count changed")
    witness = polar_strict_witness()
    if lorentz_pair(witness, witness) != POLAR_STRICT_NORM:
        raise AssertionError("the strict witness norm changed")

    return {
        "schema_version": "0.1.0",
        "seed_id": "silver-H26-T26-a-2sqrt2-minus1",
        "field": "Q(sqrt(2))",
        "ring": "Z[sqrt(2)]",
        "hyperbolic_dimension": 26,
        "lattice": "<-(2sqrt(2)-1)> + T26",
        "scope": (
            "Complete exact Vinberg prefix through the first nonempty packet "
            "H=5+4sqrt(2). The strict full-shell-polar witness proves only "
            "that this prefix has infinite volume. Later mirrors may cut it; "
            "no eventual reflectivity or global nonreflectivity claim is made."
        ),
        "timelike_coefficient": value_json(RHO),
        "prime_seven_generator": value_json(P7_GENERATOR),
        "t26_source": {
            "url": T26_SOURCE_URL,
            "dimension": 26,
            "determinant": 3,
            "smith_invariants": list(smith),
            "minimum": T26_MINIMUM,
            "kissing_number": T26_KISSING_NUMBER,
            "radical_three_minimum": T26_RADICAL_THREE_MINIMUM,
            "radical_three_kissing_number": T26_RADICAL_THREE_KISSING_NUMBER,
        },
        "local_mirror_norms": {
            "class_count": len(NORM_CLASSES),
            "proof": COMPLETE_NORM_CLASS_PROOF,
            "classes": [
                {
                    "norm": value_json(row.norm),
                    "dyadic": row.dyadic,
                    "p7_gate": row.p7,
                    "discriminant_three_gate": row.discriminant_three,
                }
                for row in NORM_CLASSES
            ],
        },
        "height_zero": {
            "class_count": len(NORM_CLASSES),
            "fully_gated_gram_type_count": 0,
            "reflecting_wall_count": 0,
        },
        "positive_prefix": {
            "height_bound": value_json(FIRST_PACKET_HEIGHT),
            "control_count": len(controls),
            "distinct_height_count": len({row.height for row in controls}),
            "earlier_control_count": sum(row.height < FIRST_PACKET_HEIGHT for row in controls),
            "earlier_fully_gated_gram_type_count": 0,
            "controls": [
                {
                    "height": value_json(row.height),
                    "norm": value_json(row.norm_class.norm),
                    "control": value_json(row.control),
                    "target": value_json(row.target),
                    "reduced_target_after_p7_gate": value_json(reduced_target(row)),
                    "dyadic": row.norm_class.dyadic,
                    "p7_gate": row.norm_class.p7,
                    "discriminant_three_gate": row.norm_class.discriminant_three,
                    "fully_gated_gram_triples": [
                        [item.u_square, item.v_square, item.u_dot_v]
                        for item in gated_gram_triples(row)
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
            "gram_triple": [4, 4, 4],
            "roots": "r_b=(2+sqrt(2),(1+sqrt(2))*b), b^2=4 in T26",
            "normalized_walls": "(sqrt(2),b)",
            "wall_count": FIRST_PACKET_WALL_COUNT,
            "same_height_distinct_minimal_dot_upper": 2,
            "same_height_maximum_root_pairing": value_json(Silver(-4, -4)),
            "minimal_shell_audit": shell,
        },
        "full_shell_polar": {
            "rational_vertex_coordinates": [str(value) for value in POLAR_VERTEX],
            "square": str(POLAR_VERTEX_SQUARE),
            "wall_bound": value_json(POLAR_WALL_BOUND),
            "absolute_shell_support": str(POLAR_VERTEX_SHELL_MAXIMUM),
            "contact_sign_pairs": POLAR_VERTEX_CONTACT_COUNT,
            "boundary_witness_norm": value_json(POLAR_BOUNDARY_NORM),
        },
        "strict_hyperideal_witness": {
            "vector": [value.to_json() for value in witness],
            "construction": "z=(1,(999/1000)*(4-sqrt(2))*y_polar)",
            "shrink_scale": str(POLAR_STRICT_SCALE),
            "norm": value_json(POLAR_STRICT_NORM),
            "maximum_first_packet_pairing": value_json(
                POLAR_STRICT_FIRST_PACKET_MARGIN
            ),
            "strict": True,
        },
        "conclusion": "complete first-packet prefix has infinite hyperbolic volume",
        "discovery": None,
    }


def report() -> str:
    return """# Silver p7 T26: exact full-shell polar leaves a hyperideal patch

For the arithmetic host

`<-(2sqrt(2)-1)> + T26`

in hyperbolic dimension 26, the complete first Vinberg packet has infinite
hyperbolic volume.  This replaces the earlier, weak coordinate-ray estimate
by an exact vertex of the complete 117,936-wall shell polar.

## Complete first packet

The primary catalogue Gram matrix is even of rank 26, determinant 3, Smith
form `1^25,3`, minimum 4, and kissing number 117,936.  Its pairing radical
modulo three has minimum 24 and kissing number 1,638.

The complete primitive mirror-norm alphabet has eight classes: independent
optional dyadic, coefficient-prime-seven, and discriminant-three factors.
There are ten arithmetic controls at seven heights through
`H=5+4sqrt(2)`; all eight earlier controls are empty, and there are no
height-zero walls.  The unique realized component is

`r_b=(2+sqrt(2),(1+sqrt(2))*b),  b^2=4 in T26`,

of norm 2.  Hence the first packet is the complete 117,936-vector minimal
shell, with normalized walls `(sqrt(2),b)`.

## Exact 66-contact polar vertex

Let `y` be `1/30` times

```
(43,119,50,-11,1,9,60,5,-27,-3,-20,30,0,-18,-43,-35,
 6,30,4,-11,2,-1,-16,30,-19,14).
```

Exact arithmetic gives `y^2=97/150`.  A full PARI enumeration of all 58,968
minimal representatives modulo sign proves

`max_{b^2=4}|(y,b)|=1`,

with exactly 66 contacting sign pairs.  The canonical representative stream
has SHA-256
`69cdb97ab9d5b89e24db1b1fe94b83f5488e68f14e71e65cf04a8d991b8c4870`.

The boundary polar vector has spatial part `(4-sqrt(2))*y` and Lorentz norm

`316/25-(538/75)sqrt(2)>0`.

After shrinking by `999/1000`, the strict witness

`z=(1,(999/1000)*(4-sqrt(2))*y)`

still has norm

`315418291/25000000-(44768699/6250000)sqrt(2)>0`,

while its maximum pairing with every first-packet wall is

`-1/500-(3/1000)sqrt(2)<0`.

Thus a strict open hyperideal patch survives the complete first packet.  This
is a bounded-prefix statement only: later mirrors may cut it.

## Reproduction

```sh
PYTHONPATH=src python3 scripts/audit_silver_t26_first_packet.py
PYTHONPATH=src pytest -q tests/test_silver_t26_first_packet.py
```

Source: https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/neb4.html
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
