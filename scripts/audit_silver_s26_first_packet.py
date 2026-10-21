#!/usr/bin/env python3
"""Freeze the first packet of ``<-(2sqrt(2)-1)> + S26``."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_s26 import (
    COMPLETE_NORM_CLASS_PROOF,
    FIRST_PACKET_HEIGHT,
    NORM_CLASSES,
    RHO,
    S26_GROUP_ORDER,
    S26_KISSING_NUMBER,
    S26_MINIMUM,
    S26_SHELL_COUNTS,
    S26_SOURCE_URL,
    enumerate_controls,
    first_packet_audit,
    first_packet_root_is_reflective,
    gated_gram_triples,
    h1000_minimal_direction_audit,
    height_zero_gated_gram_triples,
    s26_catalogue_audit,
    s26_gram,
)


ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = ROOT / "data/enumerations/silver_s26_h26_first_packet.json"
REPORT_PATH = ROOT / "results/silver_s26_h26_first_packet.md"


def value_json(value: Silver) -> dict[str, object]:
    return {
        "exact": value.to_json(),
        "display": str(value),
        "signs": [value.sign(0), value.sign(1)],
    }


def gram_sha256() -> str:
    stream = "\n".join(" ".join(map(str, row)) for row in s26_gram()) + "\n"
    return hashlib.sha256(stream.encode()).hexdigest()


def _gp_matrix(matrix) -> str:
    return "[" + ";".join(",".join(map(str, row)) for row in matrix) + "]"


def pari_shell_audit() -> dict[str, object]:
    executable = shutil.which("gp")
    if executable is None:
        raise RuntimeError("PARI/GP is required to regenerate the S26 certificate")
    gram = _gp_matrix(s26_gram())
    program = (
        f"G={gram};"
        "for(i=1,2,B=[3,4][i];R=qfminim(G,B,,2);"
        'print(B,"|",R[1],"|",floor(R[2]+1/2),"|",matsize(R[3])[2]));quit()'
    )
    process = subprocess.run(
        [executable, "-q", "-s", "512000000"],
        input=program + "\n",
        text=True,
        capture_output=True,
        timeout=240,
        check=True,
    )
    lines = [line.strip() for line in process.stdout.splitlines() if line.strip()]
    expected = {
        3: (3_120, 3, 1_560),
        4: (105_300, 4, 52_650),
    }
    parsed = {}
    for line in lines:
        bound, count, minimum, representatives = map(int, line.split("|"))
        parsed[bound] = (count, minimum, representatives)
    if parsed != expected:
        raise AssertionError(f"the S26 PARI shell transcript drifted: {parsed!r}")
    exact = {
        3: parsed[3][0],
        4: parsed[4][0] - parsed[3][0],
    }
    if exact != {3: S26_SHELL_COUNTS[3], 4: S26_SHELL_COUNTS[4]}:
        raise AssertionError(f"the S26 exact shell counts drifted: {exact!r}")
    return {
        "engine": "PARI/GP qfminim with flag=2",
        "cumulative_signed_counts": {str(k): parsed[k][0] for k in parsed},
        "exact_signed_shell_counts": {str(k): value for k, value in exact.items()},
        "representatives_mod_sign": {str(k): parsed[k][2] for k in parsed},
    }


def build_payload() -> dict[str, object]:
    if not s26_catalogue_audit():
        raise AssertionError("the exact S26 catalogue audit failed")
    shell_audit = pari_shell_audit()
    audit = first_packet_audit(shell_four_count=S26_SHELL_COUNTS[4])
    h1000 = h1000_minimal_direction_audit()
    explicit_d = (-1, 0, 0, 0, 1) + (0,) * 21
    if not first_packet_root_is_reflective(explicit_d):
        raise AssertionError("the explicit first-packet root is not reflective")
    controls = enumerate_controls()
    return {
        "schema_version": "0.1",
        "seed_id": "silver-H26-S26-a-2sqrt2-minus1",
        "scope": (
            "Complete primitive mirror-norm alphabet, all height-zero exclusions, "
            "and every positive Vinberg control through H=1000. The first nonempty "
            "packet is exact; later necessary Gram supersets bound every possible "
            "wall. A strict hyperideal point proves this bounded prefix has infinite "
            "volume; walls above H=1000 may cut it."
        ),
        "field": "Q(sqrt(2))",
        "ring": "Z[sqrt(2)]",
        "hyperbolic_dimension": 26,
        "host": "<-(2*sqrt(2)-1)> + S26",
        "timelike_coefficient": value_json(RHO),
        "coefficient_norm": str(RHO.norm()),
        "spatial_lattice": {
            "name": "Conway--Borcherds S26",
            "source": S26_SOURCE_URL,
            "rank": 26,
            "determinant": 1,
            "odd_unimodular": True,
            "minimum": S26_MINIMUM,
            "rootless": True,
            "kissing_number": S26_KISSING_NUMBER,
            "automorphism_group_order": S26_GROUP_ORDER,
            "gram_stream_sha256": gram_sha256(),
            "pari_shell_audit": shell_audit,
        },
        "local_mirror_norms": {
            "class_count": len(NORM_CLASSES),
            "proof": COMPLETE_NORM_CLASS_PROOF,
            "classes": [
                {
                    "norm": value_json(row.norm),
                    "dyadic_exponent": row.dyadic_exponent,
                    "prime_7_gate": row.p7,
                }
                for row in NORM_CLASSES
            ],
        },
        "height_zero": {
            "class_count": len(NORM_CLASSES),
            "all_empty": all(not height_zero_gated_gram_triples(row) for row in NORM_CLASSES),
            "reflecting_wall_count": 0,
            "reason": (
                "After the p7 divisibility gate, every target is excluded by "
                "nonintegrality or the exact minimum-three Gram gate."
            ),
        },
        "chronology": {
            "control_count": audit.control_count,
            "distinct_height_count": audit.distinct_height_count,
            "earlier_control_count": audit.earlier_control_count,
            "earlier_surviving_control_count": audit.earlier_surviving_control_count,
            "controls": [
                {
                    "height": value_json(row.height),
                    "norm": value_json(row.norm_class.norm),
                    "control": value_json(row.control),
                    "target": value_json(row.target),
                    "prime_7_gate": row.norm_class.p7,
                    "gated_gram_triples": [
                        [triple.u_square, triple.v_square, triple.u_dot_v]
                        for triple in gated_gram_triples(row)
                    ],
                }
                for row in controls
            ],
        },
        "first_packet": {
            "height": value_json(FIRST_PACKET_HEIGHT),
            "norm": value_json(audit.first_norm),
            "control": value_json(audit.first_control),
            "target": value_json(audit.first_target),
            "surviving_gram_triple": [
                audit.surviving_triple.u_square,
                audit.surviving_triple.v_square,
                audit.surviving_triple.u_dot_v,
            ],
            "description": "r_d=(2+sqrt(2), (1+sqrt(2))*d), d in S26, d^2=4",
            "spatial_shell_norm": audit.scalar_shell_norm,
            "wall_count": audit.wall_count,
            "explicit_d_coordinates": list(explicit_d),
            "explicit_root_reflective": True,
        },
        "geometry": {
            "minimal_direction": "choose any e in S26 with e^2=3 and put z=(1,e)",
            "slab_bound_after_dividing_by_1+sqrt(2)": value_json(audit.slab_bound),
            "integral_pairing_bound": audit.pairing_bound,
            "witness_scale": value_json(audit.witness_scale),
            "witness_spatial_square": value_json(audit.witness_square),
            "witness_lorentz_norm": value_json(audit.witness_lorentz_norm),
            "wall_margin": value_json(audit.wall_margin),
            "strict_hyperideal_open_patch": audit.strict_hyperideal,
            "finite_volume_through_prefix": audit.finite_volume,
            "compact_through_prefix": audit.compact,
            "logic": (
                "For d^2=4 and fixed e^2=3, integrality and minimum three "
                "applied to e+-d give |e.d|<=2. Hence |e.d| < 4-sqrt(2). "
                "The point is spacelike because 3-(2sqrt(2)-1)>0; the segment from "
                "the basepoint crosses the light cone away from every facet."
            ),
        },
        "h1000_persistence": {
            "complete_height_bound": value_json(h1000.complete_height_bound),
            "control_count": h1000.control_count,
            "distinct_height_count": h1000.distinct_height_count,
            "nonempty_control_count": h1000.nonempty_control_count,
            "gated_binary_type_count": h1000.gated_binary_type_count,
            "nonempty_height_count": h1000.nonempty_height_count,
            "binary_gate_scope": (
                "For y=u+sqrt(2)v, all integral binary Gram triples satisfying "
                "positive semidefiniteness and the minimum-three short-combination "
                "gate are retained; this is a necessary superset of actual S26 pairs."
            ),
            "augmented_gate_scope": (
                "For fixed d^2=3, every integral (d.u,d.v) satisfying exact 3x3 "
                "positive semidefiniteness and excluding norms 1,2 on [-3,3]^3 "
                "is retained. This bounded gate is necessary, not a proof of true "
                "ternary minimum; maximising over its superset is safe."
            ),
            "witness": "z=(1,(43/55)d), d^2=3",
            "witness_scale": value_json(h1000.witness_scale),
            "witness_spatial_square": value_json(h1000.witness_square),
            "witness_lorentz_norm": value_json(h1000.witness_lorentz_norm),
            "all_wall_margins_strict": h1000.all_margins_strict,
            "strict_hyperideal_open_patch": h1000.strict_hyperideal,
            "worst_row": {
                "height": value_json(h1000.worst_row.control.height),
                "norm": value_json(h1000.worst_row.control.norm_class.norm),
                "control": value_json(h1000.worst_row.control.control),
                "prime_7_gate": h1000.worst_row.control.norm_class.p7,
                "gram_triple": [
                    h1000.worst_row.triple.u_square,
                    h1000.worst_row.triple.v_square,
                    h1000.worst_row.triple.u_dot_v,
                ],
                "maximizing_pair": [
                    h1000.worst_row.pairing.maximizing_u_pairing,
                    h1000.worst_row.pairing.maximizing_v_pairing,
                ],
                "maximum_pairing": value_json(
                    h1000.worst_row.pairing.maximum_pairing
                ),
                "retained_extension_count": (
                    h1000.worst_row.pairing.retained_extension_count
                ),
                "normalized_wall_bound": value_json(
                    h1000.worst_row.normalized_wall_bound
                ),
                "witness_margin": value_json(h1000.worst_row.witness_margin),
            },
            "finite_volume_through_prefix": h1000.finite_volume,
            "compact_through_prefix": h1000.compact,
        },
        "conclusion": (
            "This genuinely new H26 rootless host has an exact 102180-wall "
            "first packet, and a strict hyperideal point survives every possible "
            "wall through H=1000, but this bounded prefix has infinite volume. "
            "No compact or finite-volume example is obtained; walls above H=1000 "
            "remain open."
        ),
    }


def render_report(payload: dict[str, object]) -> str:
    g = payload["geometry"]
    f = payload["first_packet"]
    c = payload["chronology"]
    p = payload["h1000_persistence"]
    return f"""# Silver S26 host in H26: exact first packet

## Outcome

The genuinely new rootless host

\\[
 L=\\langle-(2\\sqrt2-1)\\rangle\\oplus S_{{26}}
 \\quad\\text{{in }}\\mathbf H^{{26}}
\\]

has no height-zero walls. Its complete first nonempty Vinberg packet is at
`H={FIRST_PACKET_HEIGHT}` and contains **{f['wall_count']:,}** primitive
walls. An exact strict hyperideal point survives the whole packet and, after
replacement, every possible wall through `H<=1000`. This bounded prefix has
infinite volume and is not a compact example. Walls above the frozen cutoff
may cut the point; this is not a nonreflectivity theorem.

## Exact host and local arithmetic

The spatial lattice is Conway--Borcherds' odd unimodular rank-26 lattice
`S26`, of minimum three and kissing number 3,120. The catalogue lower-triangle
Gram matrix is stored literally; exact Bareiss elimination gives determinant
one. PARI independently recovers cumulative signed counts 3,120 and 105,300
through norms three and four, hence the exact norm-four shell has 102,180
vectors.

For `rho=2sqrt(2)-1`, the primitive pairing-ideal argument leaves exactly six
norm classes

`(1, 2+sqrt(2), or 2) * (1 or 3+sqrt(2))`.

The second factor is accompanied by the indispensable spatial divisibility
`y in rho*S26`. All six height-zero classes are empty.

## Complete chronology through the first packet

There are exactly **{c['control_count']}** arithmetic controls at
**{c['distinct_height_count']}** distinct heights through the first packet.
All **{c['earlier_control_count']}** earlier controls have empty exact
minimum-three Gram censuses.

At `H={FIRST_PACKET_HEIGHT}` the sole surviving component is

\\[
 q=2,\\quad x=2+\\sqrt2,\\quad
 y=(1+\\sqrt2)d,\\qquad d^2=4.
\\]

The Gram triple `(u^2,v^2,u.v)=(4,4,4)` forces this scalar form, so every
norm-four vector occurs and the packet has exactly 102,180 walls.

## Exact infinite-volume certificate

Fix a minimal vector `e` of `S26`, so `e^2=3`, and put

\\[
 z=(1,e).
\\]

For every `d^2=4`, integrality and minimum three applied to `e+-d` give
`|e.d|<=2`. After dividing a
packet root by `1+sqrt(2)`, its section inequality is
`z.d <= 4-sqrt(2)`, and

\\[
 |e.d|\\le2<4-\\sqrt2.
\\]

The wall margin is `{g['wall_margin']['display']}`. Meanwhile

\\[
 z^2=3,\\qquad z^2-(2\\sqrt2-1)
 ={g['witness_lorentz_norm']['display']}>0.
\\]

Thus `z` is strictly spacelike and strictly inside every first-packet slab.
The segment from the Vinberg basepoint to `z` crosses the light cone away
from every facet, producing an open ideal patch. The complete first-packet
chamber has infinite volume.

## Persistence through every possible wall at H<=1000

The fixed minimal direction can be continued much farther. Put

\\[
 z=(1,(43/55)e),\\qquad e^2=3.
\\]

The exhaustive arithmetic census contains **{p['control_count']}** controls
at **{p['distinct_height_count']}** distinct heights. Exactly
**{p['nonempty_control_count']}** controls at **{p['nonempty_height_count']}**
heights retain a necessary binary Gram type, for
**{p['gated_binary_type_count']}** types in total.

For each type `(u^2,v^2,u.v)`, the audit maximises `e.(u+sqrt(2)v)` over
every integral augmented Gram that is positive semidefinite and has no norm
one or two among combinations in `[-3,3]^3`. This local test is explicitly
only a necessary gate, not a proof of true ternary minimum. Retaining the
larger superset makes the resulting upper bound valid for every actual S26
root.

The tight row occurs at

\\[
 H={p['worst_row']['height']['display']},\\quad q=2,\\quad
 x={p['worst_row']['control']['display']},\\quad
 (u^2,v^2,u.v)=(460,226,322).
\\]

Its maximum is `37+26sqrt(2)`, attained at `(e.u,e.v)=(37,26)`, while the
normalised wall bound is `28+21sqrt(2)`. The exact remaining margin is

\\[
 {p['worst_row']['witness_margin']['display']}>0.
\\]

Every other row has larger margin. The witness remains spacelike by

\\[
 3(43/55)^2-(2\\sqrt2-1)
 ={p['witness_lorentz_norm']['display']}>0.
\\]

Consequently an open hyperideal patch survives every possible wall, hence
every Vinberg-accepted wall, through the complete `H<=1000` prefix.

## Reproduction

```bash
PYTHONPATH=src python3 scripts/audit_silver_s26_first_packet.py --verify
PYTHONPATH=src pytest -q tests/test_silver_s26_first_packet.py
```
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    payload = build_payload()
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    report = render_report(payload)
    if args.verify:
        if JSON_PATH.read_text() != text:
            raise SystemExit(f"certificate drift: {JSON_PATH}")
        if REPORT_PATH.read_text() != report:
            raise SystemExit(f"report drift: {REPORT_PATH}")
        print("silver S26 H26 certificate verified")
        return
    JSON_PATH.write_text(text)
    REPORT_PATH.write_text(report)
    print(JSON_PATH)
    print(REPORT_PATH)


if __name__ == "__main__":
    main()
