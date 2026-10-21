#!/usr/bin/env python3
"""Freeze the first packet of ``<-(12sqrt(2)-13)> + K12^2``."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_coxeter_todd_squared import (
    COEFFICIENT,
    COMPLETE_NORM_CLASS_PROOF,
    FIRST_PACKET_HEIGHT,
    FIRST_PACKET_WALL_COUNT,
    K12_GRAM,
    K12_GRAM_STREAM_SHA256,
    K12_RADICAL_THREE_MINIMUM,
    K12_SHELL_COUNTS,
    K12_SMITH_INVARIANTS,
    K12_SOURCE_URL,
    NORM_CLASSES,
    P17_GENERATOR,
    P7_GENERATOR,
    WITNESS_SCALE,
    coefficient_factorization_audit,
    dual_direction_audit,
    enumerate_controls,
    first_packet_chronology_audit,
    first_packet_components,
    first_packet_root,
    first_packet_root_is_reflective,
    gated_gram_triples,
    height_zero_gated_gram_triples,
    k12_basis_audit,
    k12_radical_three_gram,
    lorentz_pair,
    same_height_acuteness_bounds,
    strict_hyperideal_audit,
    strict_hyperideal_witness,
    strict_witness_norm,
    strict_witness_pairing_bounds,
)


ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = ROOT / "data/enumerations/silver_coxeter_todd_squared_first_packet.json"
REPORT_PATH = ROOT / "results/silver_coxeter_todd_squared_first_packet.md"


def value_json(value: Silver) -> dict[str, object]:
    return {
        "exact": value.to_json(),
        "display": str(value),
        "signs": [value.sign(0), value.sign(1)],
    }


def vector_json(vector) -> list[dict[str, str]]:
    return [value.to_json() for value in vector]


def _gp_matrix(matrix) -> str:
    return "[" + ";".join(",".join(map(str, row)) for row in matrix) + "]"


def _run_gp(program: str) -> list[str]:
    executable = shutil.which("gp")
    if executable is None:
        raise RuntimeError("PARI/GP is required to regenerate this certificate")
    process = subprocess.run(
        [executable, "-q", "-s", "512000000"],
        input=program + "\n",
        text=True,
        capture_output=True,
        timeout=180,
        check=True,
    )
    return [line.strip() for line in process.stdout.splitlines() if line.strip()]


def pari_shell_audit() -> dict[str, object]:
    gram = _gp_matrix(K12_GRAM)
    lines = _run_gp(
        f"G={gram};"
        "for(i=1,2,B=[4,8][i];R=qfminim(G,B,,2);"
        'print(B,"|",R[1],"|",matsize(R[3])[2]));'
        "quit()"
    )
    # The loop intentionally visits 4 and 8; query six separately to keep
    # the GP source legible and deterministic.
    six = _run_gp(f"R=qfminim({gram},6,,2);print(R[1]);print(matsize(R[3])[2]);quit()")
    if len(lines) != 2 or len(six) != 2:
        raise AssertionError(f"unexpected qfminim transcript: {lines!r}, {six!r}")
    cumulative = {
        4: int(lines[0].split("|")[1]),
        6: int(six[0]),
        8: int(lines[1].split("|")[1]),
    }
    representatives = {
        4: int(lines[0].split("|")[2]),
        6: int(six[1]),
        8: int(lines[1].split("|")[2]),
    }
    exact = {
        4: cumulative[4],
        6: cumulative[6] - cumulative[4],
        8: cumulative[8] - cumulative[6],
    }
    if exact != K12_SHELL_COUNTS or representatives != {
        4: 378,
        6: 2_394,
        8: 12_600,
    }:
        raise AssertionError(
            f"the Coxeter--Todd shell census drifted: {exact!r}, {representatives!r}"
        )

    radical = k12_radical_three_gram()
    radical_lines = _run_gp(
        f"R=qfminim({_gp_matrix(radical)},12,,2);"
        "print(R[1]);print(floor(R[2]+1/2));print(matsize(R[3])[2]);quit()"
    )
    if radical_lines != ["756", "12", "378"]:
        raise AssertionError(f"the ternary radical audit drifted: {radical_lines!r}")
    return {
        "engine": "PARI/GP qfminim with flag=2",
        "ordinary_K12": {
            "cumulative_signed_counts": {str(key): value for key, value in cumulative.items()},
            "exact_signed_shell_counts": {str(key): value for key, value in exact.items()},
            "representatives_mod_sign": {
                str(key): value for key, value in representatives.items()
            },
        },
        "pairing_radical_mod_3": {
            "basis_index": 729,
            "minimum": 12,
            "minimal_shell_signed_count": 756,
            "representatives_mod_sign": 378,
        },
    }


def gram_stream_sha256() -> str:
    stream = "\n".join(" ".join(map(str, row)) for row in K12_GRAM) + "\n"
    return hashlib.sha256(stream.encode()).hexdigest()


def build_payload() -> dict[str, object]:
    if not k12_basis_audit():
        raise AssertionError("the exact K12 basis audit failed")
    if gram_stream_sha256() != K12_GRAM_STREAM_SHA256:
        raise AssertionError("the normalized K12 Gram stream changed")
    if not coefficient_factorization_audit():
        raise AssertionError("the norm-119 coefficient factorization changed")
    if not first_packet_chronology_audit():
        raise AssertionError("the first-packet chronology audit failed")
    if not strict_hyperideal_audit():
        raise AssertionError("the strict hyperideal certificate failed")

    controls = enumerate_controls()
    components = first_packet_components()
    shell_audit = pari_shell_audit()
    q1_example = first_packet_root(4, (1,) + (0,) * 23)
    q2_example = first_packet_root(8, (1, 1) + (0,) * 22)
    if not all(first_packet_root_is_reflective(root) for root in (q1_example, q2_example)):
        raise AssertionError("an explicit first-packet root lost reflectivity")

    return {
        "schema_version": "0.1",
        "seed_id": "silver-H24-K12-squared-a-12sqrt2-minus13",
        "scope": (
            "Complete primitive mirror-norm alphabet, all height-zero exclusions, "
            "and every positive Vinberg control through the first nonempty packet. "
            "The strict hyperideal point proves that this prefix has infinite volume; "
            "later walls may cut it."
        ),
        "field": "Q(sqrt(2))",
        "ring": "Z[sqrt(2)]",
        "hyperbolic_dimension": 24,
        "host": "<-(12*sqrt(2)-13)> + K12 + K12",
        "timelike_coefficient": value_json(COEFFICIENT),
        "coefficient_norm": str(COEFFICIENT.norm()),
        "coefficient_squareclass_novelty": (
            "the coefficient is unit-square-equivalent to the archived norm-119 "
            "Leech direction 9+10sqrt(2), but its rootless 3-modular spatial host "
            "K12^2 is new and changes both the ternary gates and the first packet"
        ),
        "factorization": {
            "identity": "a=(sqrt(2)-1)*(3+sqrt(2))*(5-2sqrt(2))",
            "prime_7_generator": value_json(P7_GENERATOR),
            "prime_17_generator": value_json(P17_GENERATOR),
        },
        "spatial_lattice": {
            "name": "K12^2 (two Coxeter--Todd factors)",
            "source": K12_SOURCE_URL,
            "rank": 24,
            "determinant": 3**12,
            "minimum": 4,
            "rootless": True,
            "K12_determinant": 729,
            "K12_smith_invariants": list(K12_SMITH_INVARIANTS),
            "K12_gram_stream_sha256": K12_GRAM_STREAM_SHA256,
            "ternary_pairing_radical_minimum": K12_RADICAL_THREE_MINIMUM,
            "pari_shell_audit": shell_audit,
        },
        "local_mirror_norms": {
            "class_count": len(NORM_CLASSES),
            "proof": COMPLETE_NORM_CLASS_PROOF,
            "classes": [
                {
                    "norm": value_json(row.norm),
                    "dyadic": row.dyadic,
                    "ternary_gate": row.ternary,
                    "prime_7_gate": row.p7,
                    "prime_17_gate": row.p17,
                    "odd_spatial_divisor": value_json(row.odd_spatial_divisor),
                }
                for row in NORM_CLASSES
            ],
        },
        "height_zero": {
            "class_count": len(NORM_CLASSES),
            "all_empty": all(not height_zero_gated_gram_triples(row) for row in NORM_CLASSES),
            "reflecting_wall_count": 0,
            "reason": (
                "After the p7/p17 spatial divisibility gates, every target is "
                "excluded by integrality, evenness, K12 minimum four, or the "
                "minimum-twelve mod-three pairing radical."
            ),
        },
        "chronology": {
            "control_count": len(controls),
            "distinct_height_count": len({row.height for row in controls}),
            "earlier_control_count": sum(row.height < FIRST_PACKET_HEIGHT for row in controls),
            "earlier_surviving_control_count": sum(
                bool(gated_gram_triples(row))
                for row in controls
                if row.height < FIRST_PACKET_HEIGHT
            ),
            "controls": [
                {
                    "height": value_json(row.height),
                    "norm": value_json(row.norm_class.norm),
                    "control": value_json(row.control),
                    "target": value_json(row.target),
                    "reduced_target": value_json(
                        row.target
                        / (row.norm_class.odd_spatial_divisor**2)
                    ),
                    "gated_gram_triples": [list(triple) for triple in gated_gram_triples(row)],
                }
                for row in controls
            ],
        },
        "first_packet": {
            "height": value_json(FIRST_PACKET_HEIGHT),
            "component_count": len(components),
            "wall_count": FIRST_PACKET_WALL_COUNT,
            "components": [
                {
                    "norm": value_json(component.norm),
                    "control": value_json(component.control),
                    "target": value_json(component.target),
                    "raw_gram_triples": [list(row) for row in component.raw_gram_triples],
                    "imprimitive_gram_triples": [
                        list(row) for row in component.imprimitive_gram_triples
                    ],
                    "impossible_gram_triples": [
                        list(row) for row in component.impossible_gram_triples
                    ],
                    "primitive_gram_triples": [
                        list(row) for row in component.primitive_gram_triples
                    ],
                    "spatial_shell_norm": component.spatial_shell_norm,
                    "wall_count": component.wall_count,
                }
                for component in components
            ],
            "q1_decomposition": (
                "(36,16,24) has 2u=3v; d=u-v gives u=3d,v=2d,d^2=4, "
                "so r=(3+2sqrt(2))*(1,d)"
            ),
            "q2_decomposition": (
                "(64,36,48) has 3u=4v and is sqrt(2) times a q=1 root; "
                "(68,34,48) has (2u-3v)^2=2 and is impossible; "
                "(72,32,48) has 2u=3v and yields the primitive roots "
                "r=(sqrt(2)(3+2sqrt(2)),(3+2sqrt(2))f), f^2=8"
            ),
            "shell_count_identity": "1512 + (2*20412 + 756^2) = 613872",
            "explicit_root_norms": [
                value_json(lorentz_pair(q1_example, q1_example)),
                value_json(lorentz_pair(q2_example, q2_example)),
            ],
            "same_height_acuteness_bounds": [
                {"pair": name, "maximum": value_json(value)}
                for name, value in same_height_acuteness_bounds()
            ],
        },
        "strict_hyperideal_prefix_certificate": {
            "vector": vector_json(strict_hyperideal_witness()),
            "construction": (
                "z=(1,(9/5)h,0), where h=G_K12^-1 e1 has h^2=4/3 "
                "and integral pairing with K12"
            ),
            "scale": str(WITNESS_SCALE),
            "norm": value_json(strict_witness_norm()),
            "pairing_bounds": [
                {
                    "component": name,
                    "maximum_absolute_dual_pairing": maximum_dot,
                    "maximum_wall_pairing": value_json(value),
                }
                for name, maximum_dot, value in strict_witness_pairing_bounds()
            ],
            "dual_direction_audit": dual_direction_audit(),
            "strict": True,
            "finite_volume": False,
            "compact": False,
            "open_patch_reason": (
                "Both the basepoint and z satisfy every wall strictly, while z is "
                "spacelike. Their segment crosses the light cone away from all "
                "facets, yielding an open ideal patch and hence infinite volume."
            ),
        },
        "discovery": None,
    }


def render_report(payload: dict[str, object]) -> str:
    components = payload["first_packet"]["components"]
    pairings = payload["strict_hyperideal_prefix_certificate"]["pairing_bounds"]
    return f"""# Silver Coxeter--Todd-squared host: exact first packet

## Outcome

The genuinely new rootless host

\\[
  L=\\langle-(12\\sqrt2-13)\\rangle\\oplus K_{{12}}\\oplus K_{{12}}
  \\quad\\text{{in }}\\mathbf H^{{24}}
\\]

has a complete first nonempty Vinberg packet at
`H={FIRST_PACKET_HEIGHT}`.  It contains **{payload['first_packet']['wall_count']:,}**
primitive walls.  An exact strict hyperideal point survives the whole packet,
so this prefix has infinite volume and is not a compact example.  Later walls
may cut the point; this is not a nonreflectivity theorem.

## Why this lane is new

The spatial lattice is the rootless even 3-modular lattice `K12^2`, not a
Niemeier, Leech, odd-Leech-section, or `Lambda22` host.  The coefficient

\\[
 a=12\\sqrt2-13=(\\sqrt2-1)(3+\\sqrt2)(5-2\\sqrt2)
\\]

has norm `-119=-7*17`.  It is unit-square-equivalent to the archived Leech
direction `9+10sqrt(2)`, but the host is genuinely new: `K12^2` introduces
a rank-12 ternary pairing radical and a different complete first packet.

The authoritative Coxeter--Todd Gram matrix is from the
[Nebe--Sloane lattice catalogue]({K12_SOURCE_URL}).  Exact integral checks give
`det(K12)=729`, Smith form `1^6 3^6`, positive leading minors, minimum four,
and 756 minimal vectors.  PARI independently recovers the shell counts
`756, 4032, 20412` at norms `4,6,8`.  The mod-three pairing radical has
minimum 12 and 756 minimal vectors.

## Complete local and chronological search

The primitive pairing-ideal argument leaves exactly **16** norm classes:

`(1 or 2) * (1 or 3) * (1 or (3+sqrt(2))) * (1 or (5-2sqrt(2)))`.

All 16 height-zero classes are empty.  Through the first packet there are
exactly **{payload['chronology']['control_count']}** controls at
**{payload['chronology']['distinct_height_count']}** heights.  All
**{payload['chronology']['earlier_control_count']}** earlier controls are
empty after the exact odd-prime divisibility, evenness, minimum, Cauchy, and
short-combination gates.

## First packet

The simultaneous components are:

| root norm | control | spatial shell | primitive walls |
|---:|:---|---:|---:|
| `1` | `{components[0]['control']['display']}` | `K12^2(4)` | {components[0]['wall_count']:,} |
| `2` | `{components[1]['control']['display']}` | `K12^2(8)` | {components[1]['wall_count']:,} |

For norm one, the only Gram triple `(36,16,24)` forces
`u=3d,v=2d,d^2=4`.  For norm two, `(64,36,48)` is exactly a
`sqrt(2)`-multiple of a norm-one root, while `(68,34,48)` would give the
forbidden vector `(2u-3v)^2=2`.  The remaining triple `(72,32,48)` forces
`u=3f,v=2f,f^2=8` and is primitive.  Thus

\\[
  2\\cdot756+\\bigl(2\\cdot20412+756^2\\bigr)=613872.
\\]

Minimum four and integrality give strict same-height acuteness bounds
`q1/q1={payload['first_packet']['same_height_acuteness_bounds'][0]['maximum']['display']}`,
`q1/q2={payload['first_packet']['same_height_acuteness_bounds'][1]['maximum']['display']}`,
and `q2/q2={payload['first_packet']['same_height_acuteness_bounds'][2]['maximum']['display']}`.

## Exact infinite-volume certificate

Let `h=G_K12^-1 e1`.  Then `h` lies in `K12^*`, `h^2=4/3`, and `h.d` is
integral for every `d` in `K12`.  Put

\\[
 z=(1,(9/5)h,0).
\\]

Its Lorentz norm is

\\[
 z^2={payload['strict_hyperideal_prefix_certificate']['norm']['display']}>0.
\\]

Cauchy plus integrality gives `|h.d|<=2` on the norm-four shell and
`|h.f|<=3` on every first-factor projection of a norm-eight vector.  Hence
the exact maximum pairings with the two packet components are respectively
`{pairings[0]['maximum_wall_pairing']['display']}` and
`{pairings[1]['maximum_wall_pairing']['display']}`, both strictly negative.

The basepoint is also strict.  The segment from it to the spacelike `z`
therefore crosses the light cone away from every facet, producing an open
ideal patch.  The complete first-packet chamber has infinite volume.

## Reproduction

```bash
PYTHONPATH=src python3 scripts/audit_silver_coxeter_todd_squared_first_packet.py --verify
PYTHONPATH=src pytest -q tests/test_silver_coxeter_todd_squared_first_packet.py
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
