#!/usr/bin/env python3
"""Freeze the first silver Vinberg packet on the rank-28 J2dim14 lattice."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_j2dim14 import (
    FIVE_P7_SEED,
    FIVE_SQRT2_SEED,
    FIVE_TWO_SQRT2_SEED,
    FIVE_UNIT_SEED,
    J2DIM14_DETERMINANT,
    J2DIM14_GROUP_NAME,
    J2DIM14_KISSING_NUMBER,
    J2DIM14_MINIMUM,
    J2DIM14_NORM_TEN_SHELL_COUNT,
    J2DIM14_SMITH_INVARIANTS,
    J2DIM14_SOURCE_URL,
    P7_SEED,
    SQRT2_SEED,
    TWO_SQRT2_SEED,
    UNIT_FIRST_HEIGHT,
    UNIT_SEED,
    enumerate_controls,
    gated_gram_triples,
    height_zero_gated_gram_triples,
    j2dim14_catalogue_audit,
    j2dim14_gram,
    j2dim14_radical_five_gram,
    norm_classes,
    reduced_target,
    unit_first_packet_audit,
)


ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = ROOT / "data/enumerations/silver_j2dim14_h28_first_packet.json"
REPORT_PATH = ROOT / "results/silver_j2dim14_h28_first_packet.md"


def silver(value: Silver) -> dict[str, object]:
    return {
        "exact": value.to_json(),
        "display": str(value),
        "signs": [value.sign(0), value.sign(1)],
    }


def matrix_text(matrix) -> str:
    return "\n".join(" ".join(map(str, row)) for row in matrix) + "\n"


def gp_matrix(matrix) -> str:
    return "[" + ";".join(",".join(map(str, row)) for row in matrix) + "]"


def pari_shell_audit() -> dict[str, object]:
    executable = shutil.which("gp")
    if executable is None:
        raise RuntimeError("PARI/GP is required to regenerate the certificate")
    gram = gp_matrix(j2dim14_gram())
    radical = gp_matrix(
        tuple(tuple(value // 5 for value in row) for row in j2dim14_radical_five_gram())
    )
    program = (
        f"G={gram};R={radical};"
        "print(\"det|\",matdet(G));print(\"snf|\",matsnf(G));"
        "for(i=1,2,B=[8,10][i];Q=qfminim(G,B,,2);"
        "print(\"shell|\",B,\"|\",Q[1],\"|\",floor(Q[2]+1/2),\"|\",matsize(Q[3])[2]);"
        "if(B==10,V=Q[3];M=vector(28,k,0);for(j=1,matsize(V)[2],v=V[,j];"
        "if(v~*G*v==10,for(k=1,28,M[k]=max(M[k],abs(v[k])))));print(\"max|\",M)));"
        "Q=qfminim(R,8,,2);print(\"rad|\",matdet(R),\"|\",Q[1],\"|\",floor(Q[2]+1/2));quit()"
    )
    process = subprocess.run(
        [executable, "-q", "-s", "2000000000"],
        input=program + "\n",
        text=True,
        capture_output=True,
        timeout=600,
        check=True,
    )
    lines = [line.strip() for line in process.stdout.splitlines() if line.strip()]
    det = int(next(line.split("|")[1] for line in lines if line.startswith("det|")))
    shell_lines = [line for line in lines if line.startswith("shell|")]
    shell = {}
    for line in shell_lines:
        _, bound, count, minimum, representatives = line.split("|")
        shell[int(bound)] = (int(count), int(minimum), int(representatives))
    maxima_line = next(line for line in lines if line.startswith("max|"))
    maxima = tuple(map(int, re.findall(r"-?\d+", maxima_line.split("|", 1)[1])))
    radical_line = next(line for line in lines if line.startswith("rad|"))
    _, radical_det, radical_count, radical_minimum = radical_line.split("|")
    expected_maxima = (5,5,3,14,13,7,6,4,9,10,10,7,4,2,7,5,4,2,4,6,4,4,2,2,2,2,2,2)
    if det != J2DIM14_DETERMINANT:
        raise AssertionError("PARI determinant drifted")
    if shell != {8: (23_310, 8, 11_655), 10: (398_286, 10, 199_143)}:
        raise AssertionError(f"PARI shell transcript drifted: {shell!r}")
    if shell[10][0] - shell[8][0] != J2DIM14_NORM_TEN_SHELL_COUNT:
        raise AssertionError("the exact norm-ten shell count drifted")
    if maxima != expected_maxima:
        raise AssertionError(f"the norm-ten coordinate maxima drifted: {maxima!r}")
    if (int(radical_det), int(radical_count), int(radical_minimum)) != (
        J2DIM14_DETERMINANT,
        23_310,
        8,
    ):
        raise AssertionError("the scaled radical-five shell transcript drifted")
    return {
        "engine": "PARI/GP qfminim flag=2",
        "cumulative_signed_counts": {"8": 23_310, "10": 398_286},
        "exact_signed_shell_counts": {"8": 23_310, "10": 374_976},
        "representatives_mod_sign_through_10": 199_143,
        "norm_ten_coordinate_absolute_maxima": list(maxima),
        "witness_coordinate_zero_based": 13,
        "witness_coordinate_maximum": maxima[13],
        "radical_five_scaled_form": {
            "determinant": int(radical_det),
            "cumulative_signed_count_through_8": int(radical_count),
            "minimum": int(radical_minimum),
        },
    }


def lane_summary(seed, bound: int) -> dict[str, object]:
    controls = enumerate_controls(seed, Silver(bound))
    live = tuple(row for row in controls if gated_gram_triples(row))
    first_height = live[0].height if live else None
    first = tuple(row for row in live if row.height == first_height)
    return {
        "seed_id": seed.seed_id,
        "timelike_coefficient": silver(seed.coefficient),
        "complete_norm_classes": [
            {
                "norm": silver(row.norm),
                "inert_five_gate": row.five,
                "prime_seven_gate": row.p7,
                "dyadic_spatial_gate": row.dyadic_spatial_gate,
            }
            for row in norm_classes(seed)
        ],
        "height_zero_empty": all(
            not height_zero_gated_gram_triples(row) for row in norm_classes(seed)
        ),
        "screen_bound": bound,
        "control_count": len(controls),
        "distinct_height_count": len({row.height for row in controls}),
        "necessary_nonempty_control_count": len(live),
        "first_necessary_height": silver(first_height) if first_height else None,
        "first_necessary_rows": [
            {
                "norm": silver(row.norm_class.norm),
                "control": silver(row.control),
                "target": silver(row.target),
                "reduced_target": silver(reduced_target(row)),
                "gram_types": [
                    [triple.u_square, triple.v_square, triple.u_dot_v]
                    for triple in gated_gram_triples(row)
                ],
            }
            for row in first
        ],
    }


def build_payload() -> dict[str, object]:
    if not j2dim14_catalogue_audit():
        raise AssertionError("the exact J2dim14 catalogue audit failed")
    shell = pari_shell_audit()
    audit = unit_first_packet_audit(
        maximum_shell_pairing=shell["witness_coordinate_maximum"]
    )
    seeds = (
        UNIT_SEED,
        SQRT2_SEED,
        TWO_SQRT2_SEED,
        P7_SEED,
        FIVE_UNIT_SEED,
        FIVE_SQRT2_SEED,
        FIVE_TWO_SQRT2_SEED,
        FIVE_P7_SEED,
    )
    lanes = [lane_summary(seed, 300) for seed in seeds]
    if not all(row["height_zero_empty"] for row in lanes):
        raise AssertionError("a height-zero J2dim14 packet survived")
    return {
        "schema_version": "0.1",
        "seed_id": UNIT_SEED.seed_id,
        "scope": (
            "Exact catalogue lattice and complete local mirror norms for eight "
            "coefficient lanes; necessary chronology through H=300 for all lanes; "
            "exact complete first packet and strict hyperideal witness for a=1+sqrt(2). "
            "This is bounded and does not assert global nonreflectivity."
        ),
        "field": "Q(sqrt(2))",
        "ring": "Z[sqrt(2)]",
        "hyperbolic_dimension": 28,
        "host": "<-(1+sqrt(2))> + J2dim14",
        "spatial_lattice": {
            "name": "extremal 5-modular J2dim14",
            "source": J2DIM14_SOURCE_URL,
            "rank": 28,
            "determinant": J2DIM14_DETERMINANT,
            "smith_invariants": list(J2DIM14_SMITH_INVARIANTS),
            "minimum": J2DIM14_MINIMUM,
            "kissing_number": J2DIM14_KISSING_NUMBER,
            "group_name": J2DIM14_GROUP_NAME,
            "gram_sha256": hashlib.sha256(matrix_text(j2dim14_gram()).encode()).hexdigest(),
            "pari_shell_audit": shell,
            "five_modular_audit": (
                "The preimage of the mod-five radical has index 5^14; dividing "
                "its Gram by five gives an integral lattice isometric to J2dim14, "
                "with determinant 5^14, minimum 8, and the same kissing number."
            ),
        },
        "local_mirror_norm_proof": (
            "For a primitive r=(x,y), q(r) divides 2(a*x,(y,J)). Away from "
            "the dyadic prime, inert 5, and coefficient primes, primitivity and "
            "unimodularity exclude q-primes. At inert 5, Smith form 1^14,5^14 "
            "allows exponent at most one: the 5-class forces y into the pairing "
            "preimage and, unless 5 divides a, x into 5O. At the dyadic prime, "
            "evenness gives the standard unit/pi/pi^3 alphabets. A p7 factor is "
            "optional exactly when p7 divides a and forces y in p7*J. Class number "
            "one and U+=U^2 yield the displayed canonical generators."
        ),
        "coefficient_lane_screen": lanes,
        "unit_first_packet": {
            "height": silver(UNIT_FIRST_HEIGHT),
            "control_count_through_packet": audit.control_count,
            "distinct_height_count_through_packet": audit.distinct_height_count,
            "earlier_control_count": audit.earlier_control_count,
            "earlier_necessary_survivor_count": 0,
            "norm": silver(audit.first_control.norm_class.norm),
            "control": silver(audit.first_control.control),
            "target": silver(audit.first_control.target),
            "gram_type": [
                audit.first_triple.u_square,
                audit.first_triple.v_square,
                audit.first_triple.u_dot_v,
            ],
            "description": "r_d=(2+2sqrt(2),(1+sqrt(2))*d), d^2=10",
            "wall_count": audit.wall_count,
        },
        "geometry": {
            "dual_direction": "h=G^(-1)e_13 (zero based), so h.d=d_13",
            "dual_square": str(audit.dual_square),
            "norm_ten_pairing_maximum": audit.maximum_shell_pairing,
            "witness": "z=(1,(3/2)h)",
            "witness_scale": str(audit.witness_scale),
            "witness_spatial_square": str(audit.witness_spatial_square),
            "witness_lorentz_norm": silver(audit.witness_lorentz_norm),
            "normalized_packet_slab_bound": silver(2 + 2 * Silver(0, 1)),
            "strict_wall_margin": silver(audit.wall_margin),
            "strict_hyperideal_open_patch": audit.strict_hyperideal,
            "finite_volume_through_prefix": False,
            "compact_through_prefix": False,
        },
        "conclusion": (
            "No new compact or finite-volume example: the exact 374976-wall first "
            "packet has a strict hyperideal point. Among the screened coefficients, "
            "the unit lane is earliest; sqrt(2) and 2sqrt(2) first survive only at "
            "H=96+68sqrt(2), and p7 only at H=124+88sqrt(2)."
        ),
    }


def render_report(payload: dict[str, object]) -> str:
    packet = payload["unit_first_packet"]
    geometry = payload["geometry"]
    spatial = payload["spatial_lattice"]
    rows = []
    for lane in payload["coefficient_lane_screen"]:
        first = lane["first_necessary_height"]
        rows.append(
            f"| `{lane['timelike_coefficient']['display']}` | {len(lane['complete_norm_classes'])} | "
            f"{lane['control_count']} | {lane['necessary_nonempty_control_count']} | "
            f"`{first['display'] if first else '--'}` |"
        )
    return rf"""# Silver J2dim14 host in H28: exact first packet

## Outcome

The new rank-28 host

\[
L=\langle-(1+\sqrt2)\rangle\oplus J2dim14
\]

has no height-zero walls. Its exact first nonempty Vinberg packet is at
`H={packet['height']['display']}` and contains **{packet['wall_count']:,}**
primitive walls. A strict hyperideal point survives the entire packet, so
this prefix has infinite volume and is not a compact example. This is a
bounded-prefix statement; later walls may cut the point.

## Exact lattice

The literal Nebe--Sloane catalogue Gram has determinant `{spatial['determinant']}`
and Smith invariants `1^14 5^14`. PARI recovers minimum eight and exactly
{spatial['kissing_number']:,} minimal vectors. The full preimage of the
mod-five radical has index `5^14`; its Gram divided by five is integral and
isometric to the original lattice, independently confirming the 5-modular
local gate.

PARI gives cumulative signed shell counts `23,310` through norm eight and
`398,286` through norm ten. Thus the norm-ten shell has exactly
**{J2DIM14_NORM_TEN_SHELL_COUNT:,}** vectors.

## Eight coefficient lanes through H=300

All local alphabets and all height-zero exclusions are complete. Nonempty
counts here are necessary Gram packets; only the unit first packet is also
proved exactly realized.

| coefficient | norm classes | controls | necessary nonempty | first height |
|:---|---:|---:|---:|:---|
{chr(10).join(rows)}

The four coefficients divisible by five reproduce their unscaled normalized
chronologies in the mirror-norm class carrying the inert-five factor.

## Exact first packet

There are {packet['control_count_through_packet']} controls at
{packet['distinct_height_count_through_packet']} heights through the endpoint.
All {packet['earlier_control_count']} earlier controls have empty necessary
minimum-eight Gram censuses. The sole surviving row is

\[
q=2,\quad x=2+2\sqrt2,\quad
y=(1+\sqrt2)d,\qquad d^2=10.
\]

The binary Gram `(10,10,10)` forces the scalar form, so the packet is the
complete norm-ten shell and has {packet['wall_count']:,} walls.

## Strict infinite-volume certificate

Let `h=G^(-1)e_13`, using zero-based index 13. Then `h^2=8/5` and
`h.d=d_13`. Direct enumeration of the complete norm-ten shell gives the
coordinate-absolute maxima

`5,5,3,14,13,7,6,4,9,10,10,7,4,2,7,5,4,2,4,6,4,4,2,2,2,2,2,2`,

hence `|h.d|<=2`. Put `z=(1,(3/2)h)`. Its Lorentz excess is

\[
\frac94\frac85-(1+\sqrt2)=\frac{{13}}5-\sqrt2>0.
\]

After dividing a packet root by `1+sqrt(2)`, the slab bound is
`2+2sqrt(2)`. The witness pairing is at most three, leaving exact margin

\[
-1+2\sqrt2>0.
\]

Thus an open hyperideal patch survives every first-packet wall.

## Reproduction

```bash
PYTHONPATH=src python3 scripts/audit_silver_j2dim14_first_packet.py --verify
PYTHONPATH=src pytest -q tests/test_silver_j2dim14_first_packet.py
```
"""


def write_artifacts() -> None:
    payload = build_payload()
    JSON_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    JSON_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    REPORT_PATH.write_text(render_report(payload))


def verify_artifacts() -> None:
    expected_json = json.dumps(build_payload(), indent=2, sort_keys=True) + "\n"
    expected_report = render_report(json.loads(expected_json))
    if JSON_PATH.read_text() != expected_json:
        raise AssertionError(f"stale certificate: {JSON_PATH}")
    if REPORT_PATH.read_text() != expected_report:
        raise AssertionError(f"stale report: {REPORT_PATH}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.verify:
        verify_artifacts()
    else:
        write_artifacts()
    print(JSON_PATH)
    print(REPORT_PATH)


if __name__ == "__main__":
    main()
