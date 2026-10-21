#!/usr/bin/env python3
"""Freeze the first cutter and replacement ray for the silver J2dim14 host."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path

from coxeter_polytopes.exact.silver import SQRT2, Silver
from coxeter_polytopes.search.silver_j2dim14 import (
    EPSILON,
    J2DIM14_NORM_TEN_SHELL_COUNT,
    UNIT_CONTINUATION_CONTROL,
    UNIT_CONTINUATION_HEIGHT,
    UNIT_CUTTER_DIFFERENCE,
    UNIT_CUTTER_U,
    UNIT_CUTTER_V,
    UNIT_FIRST_PACKET_MAXIMIZER,
    UNIT_FIRST_HEIGHT,
    UNIT_REPLACEMENT_SCALE,
    UNIT_SEED,
    UNIT_WITNESS_SCALE,
    coordinate_slice_vectors,
    gated_gram_triples,
    j2_pair,
    j2dim14_gram,
    unit_cutter_root_audit,
    unit_h137_continuation_audit,
)


ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = ROOT / "data/enumerations/silver_j2dim14_h137_replacement.json"
REPORT_PATH = ROOT / "results/silver_j2dim14_h137_replacement.md"
FIRST_JSON_PATH = ROOT / "data/enumerations/silver_j2dim14_h28_first_packet.json"


def silver(value: Silver) -> dict[str, object]:
    return {
        "exact": value.to_json(),
        "display": str(value),
        "signs": [value.sign(0), value.sign(1)],
    }


def gp_matrix(matrix) -> str:
    return "[" + ";".join(",".join(map(str, row)) for row in matrix) + "]"


def gp_column(vector) -> str:
    return "[" + ",".join(map(str, vector)) + "]~"


def pari_first_packet_pairing_audit() -> dict[str, object]:
    executable = shutil.which("gp")
    if executable is None:
        raise RuntimeError("PARI/GP is required to regenerate the continuation")
    gram = gp_matrix(j2dim14_gram())
    u = gp_column(UNIT_CUTTER_U)
    v = gp_column(UNIT_CUTTER_V)
    program = (
        f"G={gram};U={u};V={v};Q=qfminim(G,10,,2);D=Q[3];"
        "gu=G*U;gv=G*V;c8=0;c10=0;m8=0;m10=0;"
        "ba=-10^9;bb=0;bj=0;bs=1;"
        "for(j=1,matsize(D)[2],d=D[,j];n=d~*G*d;"
        "if(n==8,c8++;m8=max(m8,abs(d[14])));"
        "if(n==10,c10++;m10=max(m10,abs(d[14]));aa=d~*gu;bb0=d~*gv;"
        "forstep(s=-1,1,2,a=s*aa;b=s*bb0;"
        "if(a+b*sqrt(2)>ba+bb*sqrt(2),ba=a;bb=b;bj=j;bs=s))));"
        "print(\"count|\",Q[1],\"|\",c8,\"|\",c10);"
        "print(\"coord|\",m8,\"|\",m10);"
        "print(\"pair|\",ba,\"|\",bb,\"|\",bj,\"|\",bs,\"|\",Vec(bs*D[,bj]));quit()"
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
    if not any(line.startswith("count|") for line in lines):
        raise RuntimeError(
            "PARI continuation transcript is missing:\nSTDOUT:\n"
            + process.stdout
            + "\nSTDERR:\n"
            + process.stderr
        )
    count = tuple(map(int, next(line for line in lines if line.startswith("count|")).split("|")[1:]))
    coordinate = tuple(map(int, next(line for line in lines if line.startswith("coord|")).split("|")[1:]))
    pairing_parts = next(line for line in lines if line.startswith("pair|")).split("|", 5)
    pairing = (int(pairing_parts[1]), int(pairing_parts[2]))
    representative_index = int(pairing_parts[3])
    sign = int(pairing_parts[4])
    maximizer = tuple(map(int, re.findall(r"-?\d+", pairing_parts[5])))
    if count != (398_286, 11_655, 187_488):
        raise AssertionError(f"the PARI first-shell counts drifted: {count!r}")
    if coordinate != (2, 2):
        raise AssertionError(f"the short-shell coordinate bounds drifted: {coordinate!r}")
    if pairing != (27, 20) or maximizer != UNIT_FIRST_PACKET_MAXIMIZER:
        raise AssertionError(
            f"the cutter-to-first-packet maximum drifted: {pairing!r}, {maximizer!r}"
        )
    return {
        "engine": (
            "PARI/GP qfminim flag=2; the selected maximum is replayed "
            "exactly in Q(sqrt(2))"
        ),
        "cumulative_signed_count_through_10": count[0],
        "representatives_mod_sign": {"norm_8": count[1], "norm_10": count[2]},
        "exact_signed_norm_ten_count": 2 * count[2],
        "coordinate_13_absolute_maxima": {"norm_8": coordinate[0], "norm_10": coordinate[1]},
        "cutter_pairing_maximum": silver(Silver(*pairing)),
        "maximizing_norm_ten_vector": list(maximizer),
        "representative_index_one_based": representative_index,
        "representative_sign": sign,
    }


def build_payload() -> dict[str, object]:
    if not unit_cutter_root_audit():
        raise AssertionError("the explicit H137 cutter failed its arithmetic audit")
    pari = pari_first_packet_pairing_audit()
    continuation = unit_h137_continuation_audit(
        cutter_first_packet_maximum=Silver(27, 20)
    )
    slice12 = coordinate_slice_vectors(12, 4)
    slice14 = coordinate_slice_vectors(14, 4)
    if not continuation.replacement_all_margins_strict:
        raise AssertionError("the J2dim14 replacement ray ceased to be strict")
    first_hash = hashlib.sha256(FIRST_JSON_PATH.read_bytes()).hexdigest()
    rows = []
    for row in continuation.rows:
        rows.append(
            {
                "height": silver(row.control.height),
                "norm": silver(row.control.norm_class.norm),
                "control": silver(row.control.control),
                "target": silver(row.control.target),
                "gram_type": [
                    row.triple.u_square,
                    row.triple.v_square,
                    row.triple.u_dot_v,
                ],
                "retained_augmented_coordinate_pairs": row.pairing_bound.retained_count,
                "maximum_h_pairing": silver(row.pairing_bound.maximum),
                "maximizing_coordinates": [
                    row.pairing_bound.maximizing_u_coordinate,
                    row.pairing_bound.maximizing_v_coordinate,
                ],
                "old_witness_margin": silver(row.old_witness_margin),
                "replacement_margin": silver(row.replacement_margin),
            }
        )
    return {
        "schema_version": "0.1",
        "seed_id": UNIT_SEED.seed_id,
        "scope": (
            "Complete unit-coefficient arithmetic controls through H=137+97sqrt(2), "
            "exact short-shell/CVP coordinate gates, an explicit chronologically "
            "accepted cutter, and a strict replacement hyperideal ray. This is a "
            "bounded prefix, not global nonreflectivity or Vinberg termination."
        ),
        "field": "Q(sqrt(2))",
        "hyperbolic_dimension": 28,
        "host": "<-(1+sqrt(2))> + J2dim14",
        "first_packet_dependency": {
            "path": str(FIRST_JSON_PATH.relative_to(ROOT)),
            "sha256": first_hash,
            "height": silver(UNIT_FIRST_HEIGHT),
            "wall_count": J2DIM14_NORM_TEN_SHELL_COUNT,
        },
        "short_shell_coordinate_certificate": {
            "pari": pari,
            "logic": (
                "For h=G^-1 e_13, h.d=d_13 and h^2=8/5. Cauchy gives "
                "|d_13|<=4 at norms 12 and 14. Exact rational LDL affine-slice "
                "enumeration finds two vectors at norm 12, coordinate +4, and no "
                "vector at norm 14, coordinate +4. Therefore the safe bounds at "
                "norms 8,10,12,14 are respectively 2,2,4,3."
            ),
            "norm_12_coordinate_4": {
                "node_count": slice12.node_count,
                "vector_count": len(slice12.vectors),
                "vectors": [list(vector) for vector in slice12.vectors],
            },
            "norm_14_coordinate_4": {
                "node_count": slice14.node_count,
                "vector_count": len(slice14.vectors),
            },
        },
        "chronology": {
            "height_bound": silver(UNIT_CONTINUATION_HEIGHT),
            "control_count": continuation.control_count,
            "distinct_height_count": continuation.distinct_height_count,
            "nonempty_control_count": continuation.nonempty_control_count,
            "gated_type_count": continuation.gated_type_count,
            "rows": rows,
            "first_old_witness_threat_height": silver(
                continuation.first_old_witness_threat_height
            ),
            "no_earlier_cutter_logic": (
                "The only later necessary rows before H137 are H82 and H103. "
                "Their difference vectors have norm 8 or 10, so the exact "
                "coordinate bound two reduces their maximum pairings to 8+6sqrt(2); "
                "the old witness margins remain respectively 2+sqrt(2) and "
                "4+2sqrt(2)."
            ),
        },
        "explicit_cutter": {
            "height": silver(UNIT_CONTINUATION_HEIGHT),
            "norm": silver(Silver(2)),
            "control": silver(UNIT_CONTINUATION_CONTROL),
            "u": list(UNIT_CUTTER_U),
            "v": list(UNIT_CUTTER_V),
            "difference_u_minus_v": list(UNIT_CUTTER_DIFFERENCE),
            "gram_u_v": [
                int(j2_pair(UNIT_CUTTER_U, UNIT_CUTTER_U)),
                int(j2_pair(UNIT_CUTTER_V, UNIT_CUTTER_V)),
                int(j2_pair(UNIT_CUTTER_U, UNIT_CUTTER_V)),
            ],
            "description": "r=(8+5sqrt(2), u+sqrt(2)v)",
            "primitive": True,
            "reflection_integral": True,
            "old_witness_pairing": silver(
                continuation.cutter_old_witness_pairing
            ),
            "cuts_old_witness": continuation.cutter_old_witness_pairing.sign(0) > 0,
            "first_packet_pairing_maximum": pari["cutter_pairing_maximum"],
            "first_packet_acute_threshold": silver(
                2 * EPSILON * UNIT_CONTINUATION_CONTROL
            ),
            "first_packet_acute_margin": silver(
                continuation.cutter_first_packet_acute_margin
            ),
            "chronologically_accepted": (
                continuation.cutter_first_packet_acute_margin.sign(0) > 0
            ),
        },
        "replacement_hyperideal_ray": {
            "vector": "z'=(1,(5/4)h), h=G^-1e_13",
            "scale": str(UNIT_REPLACEMENT_SCALE),
            "lorentz_norm": silver(continuation.replacement_lorentz_norm),
            "first_packet_margin_after_normalizing_by_1+sqrt(2)": silver(
                continuation.replacement_first_packet_margin
            ),
            "worst_h137_margin": silver(
                min(row.replacement_margin for row in continuation.rows)
            ),
            "all_possible_walls_through_h137_strict": (
                continuation.replacement_all_margins_strict
            ),
            "finite_volume_through_prefix": False,
            "compact_through_prefix": False,
        },
        "conclusion": (
            "The first accepted wall that cuts the 3h/2 ray occurs at "
            "H=137+97sqrt(2). The smaller 5h/4 ray is strictly hyperideal and "
            "survives every possible wall through that complete height, so no "
            "compact or finite-volume example is obtained at this prefix."
        ),
    }


def render_report(payload: dict[str, object]) -> str:
    chronology = payload["chronology"]
    cutter = payload["explicit_cutter"]
    replacement = payload["replacement_hyperideal_ray"]
    table = []
    for row in chronology["rows"]:
        table.append(
            f"| `{row['height']['display']}` | `{row['gram_type']}` | "
            f"`{row['maximum_h_pairing']['display']}` | "
            f"`{row['old_witness_margin']['display']}` | "
            f"`{row['replacement_margin']['display']}` |"
        )
    return rf"""# Silver J2dim14 H28 continuation through H=137+97sqrt(2)

## Outcome

The old strict hyperideal ray `z=(1,(3/2)h)`, with `h=G^(-1)e_13`,
survives the complete unit-coefficient chronology through `H=103+73sqrt(2)`.
Its first genuine accepted cutter occurs at

\[
H=137+97\sqrt2.
\]

The cutter does not close the chamber: `z'=(1,(5/4)h)` is strictly
hyperideal and survives every possible wall through that height. Thus this
bounded prefix remains infinite-volume and supplies no compact example.

## Exact short-shell and CVP gate

PARI enumerates the complete norm-eight and norm-ten shells and gives
`|d_13|<=2` on each. For norms twelve and fourteen, Cauchy permits at most
four. An exact rational-LDL affine-slice enumeration visits respectively
{payload['short_shell_coordinate_certificate']['norm_12_coordinate_4']['node_count']}
and {payload['short_shell_coordinate_certificate']['norm_14_coordinate_4']['node_count']}
nodes. It finds exactly two norm-twelve vectors with coordinate `+4`, but no
norm-fourteen vector with coordinate `+4`. Hence the certified bounds are

| norm | 8 | 10 | 12 | 14 |
|:---|---:|---:|---:|---:|
| maximum `|d_13|` | 2 | 2 | 4 | 3 |

These bounds are imposed on every short combination of the binary Gram
vectors, together with exact positive semidefiniteness of the augmented Gram.
The resulting census is a necessary superset of all actual walls.

## Complete chronology

There are {chronology['control_count']} arithmetic controls at
{chronology['distinct_height_count']} heights through the endpoint. Exactly
{chronology['nonempty_control_count']} controls retain {chronology['gated_type_count']}
necessary binary types.

| height | binary Gram `(u^2,v^2,u.v)` | max `h.(u+sqrt(2)v)` | margin at `3h/2` | margin at `5h/4` |
|:---|:---|:---|:---|:---|
{chr(10).join(table)}

At `H=82+58sqrt(2)` and `H=103+73sqrt(2)`, the difference `u-v` has
norm eight or ten. The shell bound eliminates the raw Cauchy threats and
leaves strictly positive old-ray margins. The first surviving threat is the
Gram type `(136,70,97)` at `H=137+97sqrt(2)`.

## Explicit accepted cutter

Take the integral vectors `u,v` frozen in the certificate. They satisfy

\[
u^2=136,\quad v^2=70,\quad u\mathbin{{\cdot}}v=97,
\quad (u-v)^2=12,\quad (u_{{13}},v_{{13}})=(14,10).
\]

Then

\[
r=(8+5\sqrt2,u+\sqrt2v),\quad r^2=2.
\]

The root is primitive and reflection-integral. It cuts the old ray by
`{cutter['old_witness_pairing']['display']}>0`.

It is nevertheless chronologically accepted. Exhaustion of all
{J2DIM14_NORM_TEN_SHELL_COUNT:,} first-packet walls gives maximum spatial
pairing `{cutter['first_packet_pairing_maximum']['display']}`, below the acute
threshold `{cutter['first_packet_acute_threshold']['display']}` by the strict
margin `{cutter['first_packet_acute_margin']['display']}`.

## Replacement ray

For `z'=(1,(5/4)h)`,

\[
(z')^2=5/2-(1+\sqrt2)=3/2-\sqrt2>0.
\]

The normalized first-packet margin is
`{replacement['first_packet_margin_after_normalizing_by_1+sqrt(2)']['display']}`.
The worst margin over every necessary wall type through `H=137+97sqrt(2)` is
`{replacement['worst_h137_margin']['display']}`. Both are strictly positive,
so an open hyperideal patch remains.

## Reproduction

```bash
PYTHONPATH=src python3 scripts/audit_silver_j2dim14_h137_replacement.py --verify
PYTHONPATH=src pytest -q tests/test_silver_j2dim14_h137_replacement.py
```
"""


def write_artifacts() -> None:
    payload = build_payload()
    JSON_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    JSON_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    REPORT_PATH.write_text(render_report(payload))


def verify_artifacts() -> None:
    payload = build_payload()
    expected_json = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if JSON_PATH.read_text() != expected_json:
        raise AssertionError(f"stale certificate: {JSON_PATH}")
    if REPORT_PATH.read_text() != render_report(payload):
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
