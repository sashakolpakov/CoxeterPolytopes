#!/usr/bin/env python3
"""Certify the first cutter and replacement for the silver K12^2 host."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from fractions import Fraction
from pathlib import Path

from coxeter_polytopes.exact.silver import SILVER_ONE, SQRT2, Silver
from coxeter_polytopes.search.silver_coxeter_todd_squared import (
    COEFFICIENT,
    FIRST_CUTTING_CONTROL,
    FIRST_CUTTING_D,
    FIRST_CUTTING_E,
    FIRST_CUTTING_HEIGHT,
    FIRST_CUTTING_NORM,
    FIRST_CUTTING_U,
    FIRST_CUTTING_V,
    FIRST_PACKET_HEIGHT,
    K12_GRAM,
    REPLACEMENT_WITNESS_SCALE,
    WITNESS_SCALE,
    cross_gram_maximum_against_first_cutter,
    dual_coordinate_coset_minimum,
    dual_direction_level_bound,
    enumerate_controls,
    first_cutting_root,
    first_cutting_root_audit,
    gated_gram_triples,
    lorentz_pair,
    replacement_hyperideal_witness,
    replacement_witness_norm,
    strict_hyperideal_witness,
    transformed_difference_gram,
)


ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = ROOT / "data/enumerations/silver_coxeter_todd_squared_h393_replacement.json"
REPORT_PATH = ROOT / "results/silver_coxeter_todd_squared_h393_replacement.md"
FIRST_PACKET_JSON = (
    ROOT / "data/enumerations/silver_coxeter_todd_squared_first_packet.json"
)


def value_json(value: Silver) -> dict[str, object]:
    return {
        "exact": value.to_json(),
        "display": str(value),
        "signs": [value.sign(0), value.sign(1)],
    }


def _gp_matrix(matrix) -> str:
    return "[" + ";".join(",".join(map(str, row)) for row in matrix) + "]"


def _run_gp(program: str) -> list[str]:
    executable = shutil.which("gp")
    if executable is None:
        raise RuntimeError("PARI/GP is required to regenerate this certificate")
    completed = subprocess.run(
        [executable, "-q", "-s", "1500000000"],
        input=program + "\n",
        text=True,
        capture_output=True,
        timeout=180,
        check=True,
    )
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    if not lines and completed.stderr.strip():
        raise AssertionError(f"PARI emitted no output: {completed.stderr.strip()}")
    return lines


def pari_direction_and_shell_audit() -> dict[str, object]:
    """Replay the dual cosets and explicit cutter's first-packet maxima."""

    gram = _gp_matrix(K12_GRAM)
    u = "[" + ",".join(map(str, FIRST_CUTTING_U)) + "]~"
    v = "[" + ",".join(map(str, FIRST_CUTTING_V)) + "]~"
    lines = _run_gp(
        f"G={gram};U={u};V0={v};"
        "cmp2(a,b,c,d)={my(A=a-c,B=b-d,N);"
        "if(A==0,return(sign(B)));if(B==0,return(sign(A)));"
        "if(sign(A)==sign(B),return(sign(A)));N=A*A-2*B*B;"
        "if(N>0,return(sign(A)),return(sign(B)))};\n"
        "R=qfminim(G,10,,2);W=R[3];m1=10^9;m2=10^9;m3=10^9;"
        "c1=0;c2=0;c3=0;"
        "for(j=1,matsize(W)[2],forstep(s=-1,1,2,w=s*W[,j];q=w~*G*w;k=w[1];"
        "if(k==1,if(q<m1,m1=q;c1=1,if(q==m1,c1++)));"
        "if(k==2,if(q<m2,m2=q;c2=1,if(q==m2,c2++)));"
        "if(k==3,if(q<m3,m3=q;c3=1,if(q==m3,c3++)))));"
        'print("COSET|",R[1],"|",matsize(W)[2],"|",m1,"|",c1,"|",m2,"|",c2,"|",m3,"|",c3);'
        "S=qfminim(G,8,,2);W=S[3];a4=-10^9;b4=0;a8=-10^9;b8=0;n8=0;"
        "for(j=1,matsize(W)[2],forstep(s=-1,1,2,w=s*W[,j];q=w~*G*w;"
        "A=U~*G*w;B=V0~*G*w;"
        "if(q==4&&cmp2(A,B,a4,b4)>0,a4=A;b4=B);"
        "if(q<=8&&cmp2(A,B,a8,b8)>0,a8=A;b8=B;n8=q)));"
        'print("MAX4|",a4,"|",b4);print("MAX8|",a8,"|",b8,"|",n8);quit()'
    )
    if len(lines) != 3:
        raise AssertionError(f"unexpected PARI transcript: {lines!r}")
    coset = lines[0].split("|")
    max4 = lines[1].split("|")
    max8 = lines[2].split("|")
    if coset != ["COSET", "85680", "42840", "4", "160", "4", "2", "10", "160"]:
        raise AssertionError(f"dual-coset replay changed: {coset!r}")
    if max4 != ["MAX4", "36", "26"] or max8 != ["MAX8", "40", "30", "8"]:
        raise AssertionError(f"first-packet maximum changed: {max4!r}, {max8!r}")
    return {
        "engine": "PARI/GP qfminim with exact integer Q(sqrt(2)) comparison",
        "dual_coordinate_base_residues": {
            "qfminim_signed_count_through_norm_10": int(coset[1]),
            "representatives_mod_sign": int(coset[2]),
            "rows": [
                {"coordinate": 0, "minimum": 0, "minimizer_count": 1},
                {"coordinate": 1, "minimum": int(coset[3]), "minimizer_count": int(coset[4])},
                {"coordinate": 2, "minimum": int(coset[5]), "minimizer_count": int(coset[6])},
                {"coordinate": 3, "minimum": int(coset[7]), "minimizer_count": int(coset[8])},
            ],
            "translation": "3*h is integral, has coordinate 4 and norm 12",
            "formula": "mu(k)=3*k^2/4 + rho(k mod 4), rho=(0,13/4,1,13/4)",
        },
        "explicit_cutter_first_packet_maxima": {
            "norm_four_K12_shell": value_json(Silver(int(max4[1]), int(max4[2]))),
            "norm_at_most_eight_K12_shell": value_json(Silver(int(max8[1]), int(max8[2]))),
            "maximizing_norm_for_second_row": int(max8[3]),
        },
    }


def _control(height: Silver, norm: Silver):
    rows = tuple(
        row
        for row in enumerate_controls(FIRST_CUTTING_HEIGHT)
        if row.height == height and row.norm_class.norm == norm
    )
    if len(rows) != 1:
        raise AssertionError(f"control lookup changed for {height}, {norm}: {rows!r}")
    return rows[0]


def _bound_json(bound, old_margin: Silver, replacement_margin: Silver) -> dict[str, object]:
    return {
        "height": value_json(bound.height),
        "norm": value_json(bound.norm),
        "control": value_json(bound.control),
        "gated_gram_type_count": bound.gated_gram_type_count,
        "finite_directional_pair_count": bound.directional_pair_count,
        "maximum_spatial_pairing": value_json(bound.maximum_spatial_pairing),
        "maximizing_gram_type": list(bound.maximizing_gram_type),
        "maximizing_pair_hu_hv": list(bound.maximizing_pair),
        "old_9_over_5_witness_margin": value_json(old_margin),
        "replacement_7_over_4_witness_margin": value_json(replacement_margin),
    }


def _cross_json(control, triple: tuple[int, int, int]) -> dict[str, object]:
    transformed = transformed_difference_gram(triple)
    maximum = cross_gram_maximum_against_first_cutter(transformed)
    root_margin = (
        control.norm_class.odd_spatial_divisor * maximum.maximum
        - COEFFICIENT * FIRST_CUTTING_CONTROL * control.control
    )
    return {
        "height": value_json(control.height),
        "norm": value_json(control.norm_class.norm),
        "original_gram_type": list(triple),
        "difference_gram": [list(row) for row in transformed],
        "first_column_count": maximum.first_column_count,
        "second_column_count": maximum.second_column_count,
        "psd_cross_gram_count": maximum.feasible_cross_gram_count,
        "maximum_unscaled_spatial_pairing": value_json(maximum.maximum),
        "maximizing_cross_entries_A_B_C_D": list(maximum.maximizing_cross_entries),
        "maximum_root_pairing": value_json(root_margin),
    }


def build_payload() -> dict[str, object]:
    if not FIRST_PACKET_JSON.exists():
        raise AssertionError("the frozen first-packet dependency is missing")
    if not first_cutting_root_audit():
        raise AssertionError("the explicit first cutting root failed replay")
    if replacement_witness_norm() != Silver(Fraction(205, 12), -12):
        raise AssertionError("the replacement norm changed")
    if replacement_witness_norm().sign(0) <= 0:
        raise AssertionError("the replacement is not hyperideal")

    pari = pari_direction_and_shell_audit()
    controls = enumerate_controls(FIRST_CUTTING_HEIGHT)
    nonempty = tuple(row for row in controls if gated_gram_triples(row))
    expected = (
        (FIRST_PACKET_HEIGHT, Silver(1), 1),
        (FIRST_PACKET_HEIGHT, Silver(2), 2),
        (Silver(356, 252), Silver(6, 2), 2),
        (Silver(362, 256), Silver(10, -4), 1),
        (FIRST_CUTTING_HEIGHT, Silver(1), 4),
        (FIRST_CUTTING_HEIGHT, Silver(2), 11),
    )
    observed = tuple(
        (row.height, row.norm_class.norm, len(gated_gram_triples(row)))
        for row in nonempty
    )
    if observed != expected:
        raise AssertionError(f"the H393 necessary frontier changed: {observed!r}")
    if len(controls) != 128 or len({row.height for row in controls}) != 87:
        raise AssertionError("the complete H393 control chronology changed")

    bounds = []
    for row in nonempty:
        bound = dual_direction_level_bound(row)
        old_margin = WITNESS_SCALE * bound.maximum_spatial_pairing - COEFFICIENT * row.control
        replacement_margin = (
            REPLACEMENT_WITNESS_SCALE * bound.maximum_spatial_pairing
            - COEFFICIENT * row.control
        )
        bounds.append(_bound_json(bound, old_margin, replacement_margin))
    post_first = [row for row in bounds if Silver.from_json(row["height"]["exact"]) > FIRST_PACKET_HEIGHT]
    if [row["old_9_over_5_witness_margin"]["signs"][0] for row in post_first] != [
        -1,
        -1,
        -1,
        1,
    ]:
        raise AssertionError("the first-cutting direction chronology changed")
    if any(row["replacement_7_over_4_witness_margin"]["signs"][0] >= 0 for row in bounds):
        raise AssertionError("the replacement ceased to be strict")

    cutter = first_cutting_root()
    old_pairing = lorentz_pair(strict_hyperideal_witness(), cutter)
    replacement_pairing = lorentz_pair(replacement_hyperideal_witness(), cutter)
    if old_pairing != Silver(Fraction(8, 5), Fraction(-1, 5)) or old_pairing.sign(0) <= 0:
        raise AssertionError("the explicit wall no longer cuts the old witness")
    if replacement_pairing != Silver(Fraction(1, 2), -1) or replacement_pairing.sign(0) >= 0:
        raise AssertionError("the replacement is not strict at the explicit cutter")

    # Exact shell maxima against the 613,872-wall packet.  The candidate has
    # no second-factor component, so K12^2(4) reduces to K12(4), while
    # K12^2(8) reduces to K12(<=8).
    max4 = Silver(36, 26)
    max8 = Silver(40, 30)
    first_q1_margin = (Silver(3, 2)) * (
        max4 - COEFFICIENT * FIRST_CUTTING_CONTROL
    )
    first_q2_margin = Silver(3, 2) * (
        max8 - SQRT2 * COEFFICIENT * FIRST_CUTTING_CONTROL
    )
    if first_q1_margin.sign(0) >= 0 or first_q2_margin.sign(0) >= 0:
        raise AssertionError("the explicit cutter is not acute to the first packet")

    earlier_cross = []
    for height, norm in (
        (Silver(356, 252), Silver(6, 2)),
        (Silver(362, 256), Silver(10, -4)),
    ):
        row = _control(height, norm)
        earlier_cross.extend(_cross_json(row, triple) for triple in gated_gram_triples(row))
    if [row["maximum_root_pairing"]["display"] for row in earlier_cross] != [
        "0",
        "0",
        "-7 - 4*sqrt(2)",
    ]:
        raise AssertionError("the earlier-packet cross-Gram maxima changed")

    q1_control = _control(FIRST_CUTTING_HEIGHT, SILVER_ONE)
    same_height_q1 = tuple(
        _cross_json(q1_control, triple) for triple in gated_gram_triples(q1_control)
    )
    if any(row["maximum_root_pairing"]["signs"][0] >= 0 for row in same_height_q1):
        raise AssertionError("the cutter ceased to be acute to a simultaneous q=1 root")

    return {
        "schema_version": "0.1",
        "seed_id": "silver-H24-K12-squared-a-12sqrt2-minus13-H393-replacement",
        "scope": (
            "Complete primitive control and necessary binary-Gram chronology through "
            "H=393+278*sqrt(2). Determines the first later height capable of cutting "
            "the 9/5 dual witness, realizes an accepted cutter there, and certifies a "
            "strict 7/4 replacement through that bound. This is not Vinberg termination."
        ),
        "host": "<-(12*sqrt(2)-13)> + K12 + K12",
        "hyperbolic_dimension": 24,
        "first_packet_dependency": {
            "path": str(FIRST_PACKET_JSON.relative_to(ROOT)),
            "sha256": hashlib.sha256(FIRST_PACKET_JSON.read_bytes()).hexdigest(),
        },
        "chronology": {
            "height_bound": value_json(FIRST_CUTTING_HEIGHT),
            "control_count": len(controls),
            "distinct_height_count": len({row.height for row in controls}),
            "controls_through_first_packet": sum(row.height <= FIRST_PACKET_HEIGHT for row in controls),
            "controls_strictly_between_packets": sum(
                FIRST_PACKET_HEIGHT < row.height < FIRST_CUTTING_HEIGHT for row in controls
            ),
            "controls_at_cutting_height": sum(row.height == FIRST_CUTTING_HEIGHT for row in controls),
            "nonempty_necessary_rows": bounds,
        },
        "dual_coordinate_gate": pari["dual_coordinate_base_residues"],
        "first_cutting_wall": {
            "height": value_json(FIRST_CUTTING_HEIGHT),
            "norm": value_json(FIRST_CUTTING_NORM),
            "control": value_json(FIRST_CUTTING_CONTROL),
            "difference_vectors_D_E": [list(FIRST_CUTTING_D), list(FIRST_CUTTING_E)],
            "coefficient_vectors_U_V": [list(FIRST_CUTTING_U), list(FIRST_CUTTING_V)],
            "difference_gram": [[28, 46], [46, 84]],
            "coefficient_gram": [380, 204, 278],
            "h_pairings_U_V": [22, 16],
            "root_vector": [value.to_json() for value in cutter],
            "reflection_integral": True,
            "primitive": True,
            "primitive_proof": (
                "A nonunit divisor would have square dividing norm 2, hence is "
                "associated to sqrt(2); but U has an odd coefficient, so the "
                "spatial vector is not divisible by sqrt(2)."
            ),
            "old_witness_pairing": value_json(old_pairing),
            "first_packet_shell_replay": pari["explicit_cutter_first_packet_maxima"],
            "first_packet_q1_maximum_root_pairing": value_json(first_q1_margin),
            "first_packet_q2_maximum_root_pairing": value_json(first_q2_margin),
            "earlier_packet_cross_gram_bounds": earlier_cross,
            "same_height_q1_cross_gram_bounds": list(same_height_q1),
            "same_height_q2_argument": (
                "For two q=2 roots at the same control, g=2-(y-y')^2/2. "
                "At the identity embedding g<=2, while positive-definite Cauchy "
                "at the conjugate embedding gives |g'|<=2. Thus positive g in "
                "Z[sqrt(2)] is 1, sqrt(2), or 2. Values 1 and sqrt(2) contradict "
                "the minimum-four even K12^2 trace equation; g=2 forces equal roots."
            ),
            "accepted": True,
        },
        "replacement_hyperideal_ray": {
            "vector": [value.to_json() for value in replacement_hyperideal_witness()],
            "construction": "z'=(1,(7/4)h,0), h=G_K12^-1*e1",
            "scale": str(REPLACEMENT_WITNESS_SCALE),
            "norm": value_json(replacement_witness_norm()),
            "strict_against_every_necessary_row_through_bound": True,
            "explicit_cutter_pairing": value_json(replacement_pairing),
            "finite_volume_prefix": False,
            "compact_prefix": False,
        },
        "conclusion": {
            "first_later_cutting_height": value_json(FIRST_CUTTING_HEIGHT),
            "old_witness_is_cut": True,
            "replacement_is_strict_hyperideal": True,
            "prefix_through_bound_has_infinite_volume": True,
            "vinberg_terminated": False,
            "compact_example": False,
            "finite_volume_example": False,
        },
    }


def render_report(payload: dict[str, object]) -> str:
    chronology = payload["chronology"]
    cutter = payload["first_cutting_wall"]
    replacement = payload["replacement_hyperideal_ray"]
    return rf"""# Silver Coxeter--Todd-squared host through the first cutter

## Outcome

For

\[
L=\langle-(12\sqrt2-13)\rangle\oplus K_{{12}}\oplus K_{{12}}
\subset \mathbf H^{{24}},
\]

the old `z=(1,(9/5)h,0)` certificate first becomes cut at

\[
H={cutter['height']['display']}.
\]

An explicit primitive norm-two wall at that height is accepted by every
earlier packet and pairs with the old witness by
`{cutter['old_witness_pairing']['display']}>0`.  Replacing it by

\[
z'=(1,(7/4)h,0)
\]

gives the strict hyperideal norm
`{replacement['norm']['display']}>0` and survives every possible wall through
the same height.  Thus the complete bounded prefix still has infinite volume;
this is **not** a compact or finite-volume example and is not a termination
claim.

## Complete chronology through the bound

The 16 certified mirror-norm classes give exactly
**{chronology['control_count']} controls at {chronology['distinct_height_count']} heights**.
There are {chronology['controls_strictly_between_packets']} controls strictly
between the 613,872-wall first packet and the cutting height.  Only the
following necessary Gram rows survive:

| height | norm | Gram types | directional pairs | max `h.y` | old margin | replacement margin |
|:---|:---|---:|---:|:---|:---|:---|
""" + "\n".join(
        f"| `{row['height']['display']}` | `{row['norm']['display']}` | "
        f"{row['gated_gram_type_count']} | {row['finite_directional_pair_count']} | "
        f"`{row['maximum_spatial_pairing']['display']}` | "
        f"`{row['old_9_over_5_witness_margin']['display']}` | "
        f"`{row['replacement_7_over_4_witness_margin']['display']}` |"
        for row in chronology["nonempty_necessary_rows"]
    ) + rf"""

The exact fixed-coordinate formula

\[
\min\{{x^2:x_1=k\}}=\frac34k^2+\rho(k\bmod4),\qquad
\rho=(0,13/4,1,13/4),
\]

is replayed by PARI on the four base cosets.  It turns the apparent earlier
`H=356+252sqrt(2)` cutter into the forbidden combination of norm eight and
dual coordinate three.  The old witness therefore remains strict at both
intermediate necessary heights.  At `H=393+278sqrt(2)`, the norm-one row is
still strict, while the norm-two row has positive old margin.

## Explicit accepted cutter

Set `D=U-V` and `E=2V-U`.  The displayed vectors in the JSON certificate obey

\[
\operatorname{{Gram}}(D,E)=
\begin{{pmatrix}}28&46\\46&84\end{{pmatrix}},\qquad
\operatorname{{Gram}}(U,V)=(380,204,278),
\]

with `(h.U,h.V)=(22,16)`.  Hence

\[
r=(10+7\sqrt2;\ U+\sqrt2V,0),\qquad r^2=2.
\]

Exact `qfminim` shell maximization gives `36+26sqrt(2)` against the norm-four
shell and `40+30sqrt(2)` against the norm-at-most-eight shell.  Both are
strictly below the relevant first-packet thresholds.

All three possible intermediate binary Gram types are handled without
enumerating their walls: an exhaustive integral Schur-complement census of
the four-vector Gram matrix gives maximum cutter pairing `0, 0`, and
`-7-4sqrt(2)`.  The four simultaneous norm-one types are strictly negative.
For simultaneous norm-two roots, equal height and the two real Cauchy bounds
leave positive inner products only `1,sqrt(2),2`; the first two contradict the
even minimum-four trace equation and the last identifies the same root.

## Replacement and bounded conclusion

The largest necessary pairing at the cutting norm-two row is `22+16sqrt(2)`.
Thus the replacement margin is

\[
\frac74(22+16\sqrt2)-(38+29\sqrt2)=\frac12-\sqrt2<0.
\]

Every other row in the table is also strict.  Since the basepoint is strict,
the segment to the spacelike `z'` crosses the light cone away from every facet,
so the chamber prefix through `H=393+278sqrt(2)` has an open ideal patch and
infinite volume.  Later walls may still cut this replacement.

## Reproduction

```bash
PYTHONPATH=src python3 scripts/audit_silver_coxeter_todd_squared_h393_replacement.py --verify
PYTHONPATH=src pytest -q tests/test_silver_coxeter_todd_squared_h393_replacement.py
```
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    arguments = parser.parse_args()
    payload = build_payload()
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    report = render_report(payload)
    if arguments.verify:
        if JSON_PATH.read_text() != serialized:
            raise SystemExit(f"artifact drift: {JSON_PATH}")
        if REPORT_PATH.read_text() != report:
            raise SystemExit(f"report drift: {REPORT_PATH}")
        print(f"PASS {JSON_PATH.relative_to(ROOT)}")
        print(f"PASS {REPORT_PATH.relative_to(ROOT)}")
        return
    JSON_PATH.write_text(serialized)
    REPORT_PATH.write_text(report)
    print(JSON_PATH)
    print(REPORT_PATH)


if __name__ == "__main__":
    main()
