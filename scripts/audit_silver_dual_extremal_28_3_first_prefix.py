#!/usr/bin/env python3
"""Freeze the bounded silver ``DualExtremal(28,3)`` certificate."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from coxeter_polytopes.search.silver_dual_extremal_28_3 import (  # noqa: E402
    COMPLETE_NORM_CLASS_PROOF,
    DETERMINANT,
    DIMENSION,
    FIRST_NECESSARY_HEIGHT,
    KISSING_NUMBER,
    MINIMUM,
    NORM_CLASSES,
    RADICAL_THREE_INDEX,
    RADICAL_THREE_KISSING_NUMBER,
    RADICAL_THREE_MINIMUM,
    SMITH_INVARIANTS,
    SOURCE_URL,
    catalogue_audit,
    first_prefix_audit,
    gram,
)


ARTIFACT = ROOT / "data/enumerations/silver_dual_extremal_28_3_first_prefix.json"
REPORT = ROOT / "results/silver_dual_extremal_28_3_first_prefix.md"


def silver_text(value: object) -> str:
    return str(value).replace("sqrt(2)", "sqrt(2)")


def matrix_text(matrix: tuple[tuple[int, ...], ...]) -> str:
    return "[" + ";".join(",".join(map(str, row)) for row in matrix) + "]"


def pari_audit() -> dict[str, object]:
    executable = shutil.which("gp")
    if executable is None:
        raise RuntimeError("PARI/GP is required to regenerate this certificate")
    program = (
        f"G={matrix_text(gram())};"
        'print("DET|",matdet(G));'
        'print("SNF|",matsnf(G));'
        "M=qfminim(G,8,1,0);"
        'print("HOST|",M[1],"|",floor(M[2]+1/2));'
        "K=lift(matker(Mod(G,3)));"
        "B=mathnf(matconcat([3*matid(28),K]));"
        "Q=B~*G*B;U=qflllgram(Q);R=U~*Q*U;"
        "S=qfminim(R,12,1,0);"
        'print("S3|",matdet(B),"|",S[1],"|",floor(S[2]+1/2));quit()'
    )
    completed = subprocess.run(
        [executable, "-q", "-s", "3000000000"],
        input=program + "\n",
        text=True,
        capture_output=True,
        timeout=300,
        check=True,
    )
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    if len(lines) != 4:
        raise AssertionError(f"unexpected PARI output: {lines!r}; {completed.stderr!r}")
    determinant = int(lines[0].split("|")[1])
    smith = tuple(map(int, lines[1].split("|")[1].strip("[]").split(",")))
    host = tuple(map(int, lines[2].split("|")[1:]))
    radical = tuple(map(int, lines[3].split("|")[1:]))
    expected_smith = tuple(reversed(SMITH_INVARIANTS))
    if (
        determinant != DETERMINANT
        or smith != expected_smith
        or host != (KISSING_NUMBER, MINIMUM)
        or radical != (RADICAL_THREE_INDEX, RADICAL_THREE_KISSING_NUMBER, RADICAL_THREE_MINIMUM)
    ):
        raise AssertionError((determinant, smith, host, radical))
    return {
        "engine": "PARI/GP matdet, matsnf, qfminim, matker(Mod(G,3)), mathnf, qflllgram",
        "transcript": lines,
        "determinant": determinant,
        "smith_invariants_ascending": list(SMITH_INVARIANTS),
        "minimum": host[1],
        "minimum_shell_signed_count": host[0],
        "pairing_radical_preimage_index": radical[0],
        "pairing_radical_preimage_minimum": radical[2],
        "pairing_radical_preimage_minimum_shell_signed_count": radical[1],
    }


def build_payload() -> dict[str, object]:
    if not catalogue_audit():
        raise AssertionError("literal catalogue audit failed")
    arithmetic = pari_audit()
    audit = first_prefix_audit()
    components = []
    maximum_index = 0
    for norm, control, triples in audit.endpoint_components:
        rows = []
        for triple in triples:
            maximum = audit.maxima[maximum_index]
            margin = audit.margins[maximum_index]
            maximum_index += 1
            rows.append(
                {
                    "gram": [triple.u_square, triple.v_square, triple.u_dot_v],
                    "maximum_d_dot_y": silver_text(maximum.maximum),
                    "maximizing_pair": list(maximum.maximizing_pair),
                    "augmented_determinant": maximum.augmented_determinant,
                    "retained_local_extensions": maximum.retained_extensions,
                    "witness_margin": silver_text(margin),
                }
            )
        components.append(
            {
                "mirror_norm": silver_text(norm),
                "time_control": silver_text(control),
                "necessary_gram_types": rows,
            }
        )
    return {
        "schema_version": "0.1",
        "field": "Q(sqrt(2))",
        "hyperbolic_dimension": 28,
        "host": "<-(2*sqrt(2)-1)> + DualExtremal(28,3)",
        "source": SOURCE_URL,
        "scope": (
            "Complete arithmetic controls and a rigorous necessary spatial-Gram "
            "superset through H=93+66*sqrt(2). This is a bounded prefix theorem, "
            "not Vinberg termination or global nonreflectivity."
        ),
        "spatial_lattice": {
            "rank": DIMENSION,
            "determinant": str(DETERMINANT),
            "minimum": MINIMUM,
            "minimum_shell_signed_count": KISSING_NUMBER,
            "smith_invariants": list(SMITH_INVARIANTS),
            "pari": arithmetic,
        },
        "mirror_norms": {
            "complete_classes": [silver_text(row.norm) for row in NORM_CLASSES],
            "proof": COMPLETE_NORM_CLASS_PROOF,
            "height_zero_walls": 0,
            "height_zero_proof": (
                "Ungated norms 1 and 2 are below the host minimum 8; norms 3 and 6 "
                "use the 3-pairing radical of minimum 12; rho-gated targets are "
                "not integral after division by rho^2."
            ),
        },
        "controls": {
            "inclusive_height_bound": silver_text(FIRST_NECESSARY_HEIGHT),
            "count": audit.control_count,
            "distinct_height_count": audit.height_count,
            "earlier_nonempty_necessary_controls": audit.earlier_nonempty_controls,
            "endpoint_components": components,
        },
        "witness": {
            "description": "z=(1,d/2) for any spatial minimum vector d with d^2=8",
            "scale": silver_text(audit.witness_scale),
            "lorentz_norm": silver_text(audit.witness_excess),
            "strict": True,
            "proof_scope": (
                "Every actual root lies in the enumerated integral PSD/minimum-eight "
                "augmented-Gram superset; all six exact worst-case margins are positive."
            ),
        },
        "verdict": {
            "bounded_prefix_has_open_hyperideal_patch": True,
            "bounded_prefix_has_finite_volume": False,
            "compact_example_found": False,
            "eventual_chamber": "OPEN",
        },
    }


def render(payload: dict[str, object]) -> str:
    controls = payload["controls"]
    witness = payload["witness"]
    rows = []
    for component in controls["endpoint_components"]:
        for row in component["necessary_gram_types"]:
            rows.append(
                "| {norm} | {control} | {gram} | {maximum} | {margin} |".format(
                    norm=component["mirror_norm"],
                    control=component["time_control"],
                    gram=tuple(row["gram"]),
                    maximum=row["maximum_d_dot_y"],
                    margin=row["witness_margin"],
                )
            )
    table = "\n".join(rows)
    return f"""# Silver `DualExtremal(28,3)` first necessary prefix

## Result

For

```text
L = <-(2*sqrt(2)-1)> + DualExtremal(28,3)
```

in `H28`, the complete mirror-norm alphabet has eight classes. Exact controls
through `H={controls['inclusive_height_bound']}` contain {controls['count']}
controls on {controls['distinct_height_count']} heights. Every earlier control
is empty in the rigorous spatial Gram superset. The endpoint has the following
two norm components and six necessary Gram types.

| mirror norm | time control | `(u^2,v^2,u.v)` | max `d.y` | margin for `d/2` |
|---:|---:|---:|---:|---:|
{table}

The lattice has determinant `3^26`, minimum 8 and 1,260 minimum vectors. Its
mod-three pairing-radical preimage has index 9, exact minimum 12 and 65,520
minimum vectors; these data are independently replayed with PARI.

For any minimum vector `d`, the section point {witness['description']} has

```text
z^2 = {witness['lorentz_norm']} > 0.
```

Every displayed margin is strict even against the full necessary augmented
Gram superset. Therefore the complete prefix through the stated height has an
open hyperideal patch and infinite volume.

## Scope

This is a bounded negative screen. It neither realizes every endpoint Gram
type nor decides later Vinberg walls, termination, reflectivity, or the final
chamber. It is not a compact example and not a global nonreflectivity theorem.

Source: [Nebe--Sloane catalogue]({payload['source']}).
"""


def encoded(payload: dict[str, object]) -> bytes:
    return (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    payload = build_payload()
    artifact = encoded(payload)
    report = render(payload).encode()
    if args.verify:
        if ARTIFACT.read_bytes() != artifact or REPORT.read_bytes() != report:
            raise SystemExit("frozen silver DualExtremal certificate differs")
        print(f"PASS {hashlib.sha256(artifact).hexdigest()}")
        return
    ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    ARTIFACT.write_bytes(artifact)
    REPORT.write_bytes(report)
    print(hashlib.sha256(artifact).hexdigest())


if __name__ == "__main__":
    main()
