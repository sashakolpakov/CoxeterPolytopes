#!/usr/bin/env python3
"""Freeze the cubic-seven ``K12^2`` support screen and best first packets."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from collections import Counter
from pathlib import Path

from coxeter_polytopes.exact.cubic7 import Cubic7
from coxeter_polytopes.search.cubic7_coxeter_todd_squared import (
    LOCAL_NORM_PROOF,
    SUPPORTED_K12_SEEDS,
    TWO_THETA_K12_SEED,
    screen_cubic7_k12_seed,
    two_theta_k12_first_two_packets,
)
from coxeter_polytopes.search.silver_coxeter_todd_squared import K12_GRAM


ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = ROOT / "data/enumerations/cubic7_coxeter_todd_squared_h100.json"
REPORT_PATH = ROOT / "results/cubic7_coxeter_todd_squared_h100.md"


def _value(value: Cubic7) -> dict[str, object]:
    return {
        "exact": value.to_json(),
        "display": str(value),
        "signs": [value.sign(place) for place in range(3)],
    }


def _gram(gram) -> dict[str, object]:
    return {
        "matrix": gram.matrix,
        "determinant": gram.determinant,
        "rank": gram.rank,
    }


def _gp_matrix(matrix) -> str:
    return "[" + ";".join(",".join(map(str, row)) for row in matrix) + "]"


def minimal_shell_inner_product_audit() -> dict[str, object]:
    """Independently stream the K12 minimal shell with PARI/GP."""

    executable = shutil.which("gp")
    if executable is None:
        raise RuntimeError("PARI/GP is required to regenerate this artifact")
    program = (
        f"V=qfminim({_gp_matrix(K12_GRAM)},4,,2)[3];"
        "for(j=1,matsize(V)[2],print(vector(12,i,V[i,j])));quit()\n"
    )
    process = subprocess.run(
        [executable, "-q", "-s", "512000000"],
        input=program,
        text=True,
        capture_output=True,
        timeout=180,
        check=True,
    )
    import ast

    representatives = [
        tuple(ast.literal_eval(line))
        for line in process.stdout.splitlines()
        if line.strip()
    ]
    vectors = tuple(
        vector
        for representative in representatives
        for vector in (representative, tuple(-value for value in representative))
    )

    def pair(left, right) -> int:
        return sum(
            left[i] * K12_GRAM[i][j] * right[j]
            for i in range(12)
            for j in range(12)
            if K12_GRAM[i][j]
        )

    if len(vectors) != 756 or len(set(vectors)) != 756:
        raise AssertionError("the signed K12 minimal shell changed")
    distribution = Counter(pair(vectors[0], vector) for vector in vectors)
    expected = {-4: 1, -2: 82, -1: 160, 0: 270, 1: 160, 2: 82, 4: 1}
    if dict(sorted(distribution.items())) != expected:
        raise AssertionError("the K12 minimal-shell distribution changed")
    stream = "".join(
        ",".join(map(str, vector)) + "\n" for vector in sorted(vectors)
    )
    return {
        "engine": "PARI/GP qfminim, flag=2",
        "signed_vector_count": len(vectors),
        "representatives_mod_sign": len(representatives),
        "fixed_vector_inner_product_distribution": {
            str(key): value for key, value in expected.items()
        },
        "neighbors_at_inner_product_minus_one": expected[-1],
        "ordered_pairs_at_inner_product_minus_one": 756 * expected[-1],
        "sorted_signed_coordinate_stream_sha256": hashlib.sha256(
            stream.encode()
        ).hexdigest(),
    }


def build_payload() -> dict[str, object]:
    shell = minimal_shell_inner_product_audit()
    screens = tuple(screen_cubic7_k12_seed(seed, 100) for seed in SUPPORTED_K12_SEEDS)
    packets = two_theta_k12_first_two_packets()
    records = []
    for screen in screens:
        records.append(
            {
                "seed_id": screen.seed.seed_id,
                "radical_norm": screen.seed.radical_norm,
                "ideal_support": list(screen.seed.support),
                "support_element": _value(screen.seed.support_element),
                "timelike_coefficient": _value(screen.seed.coefficient),
                "canonical_mirror_norms": [
                    _value(norm) for norm in screen.seed.canonical_mirror_norms
                ],
                "control_level_count": len(screen.controls),
                "nonempty_necessary_level_count": len(screen.nonempty),
                "necessary_gram_type_count": sum(
                    len(types) for _level, types in screen.nonempty
                ),
                "nonempty_levels": [
                    {
                        "height": _value(level.height),
                        "norm": _value(level.norm),
                        "control": _value(level.control),
                        "target": _value(level.target),
                        "necessary_gram_types": [_gram(gram) for gram in types],
                    }
                    for level, types in screen.nonempty
                ],
            }
        )
    nonempty = [record for record in records if record["nonempty_necessary_level_count"]]
    return {
        "schema": "coxeter.cubic7_k12_squared_h100.v1",
        "field": "Q(theta), theta^3+theta^2-2theta-1=0",
        "hyperbolic_dimension": 24,
        "spatial_lattice": "K12 + K12",
        "spatial_determinant": 3**12,
        "spatial_minimum": 4,
        "local_arithmetic": {
            "proof": LOCAL_NORM_PROOF,
            "manifest_scope": (
                "all squarefree support products generated by inert (2), "
                "discriminant (3), ramified (2-theta), and at most one of "
                "the three primes above 13, radical norm at most 351"
            ),
            "seed_count": len(records),
        },
        "height_bound": 100,
        "seeds": records,
        "screen_summary": {
            "seed_count": len(records),
            "empty_necessary_seed_count": len(records) - len(nonempty),
            "nonempty_necessary_seed_count": len(nonempty),
            "nonempty_seed_ids": [record["seed_id"] for record in nonempty],
            "warning": (
                "emptiness is conclusive in this height window; a nonempty "
                "necessary Gram type need not be a realized K12^2 vector"
            ),
        },
        "selected_lane": {
            "seed_id": TWO_THETA_K12_SEED.seed_id,
            "selection_reason": (
                "earliest genuinely multidirectional realized packet among "
                "the screened squareclasses"
            ),
            "minimal_shell_audit": shell,
            "first_packet": {
                "height": _value(packets.first_level.height),
                "norm": _value(packets.first_level.norm),
                "control": _value(packets.first_level.control),
                "target": _value(packets.first_level.target),
                "gram_type": _gram(packets.first_gram_type),
                "realization": "y=(2theta+theta^2)M, M minimal in either K12 factor",
                "wall_count": packets.first_wall_count,
                "same_height_pairing_maximum": _value(
                    packets.first_same_height_maximum
                ),
            },
            "second_packet": {
                "height": _value(packets.second_level.height),
                "norm": _value(packets.second_level.norm),
                "control": _value(packets.second_level.control),
                "target": _value(packets.second_level.target),
                "gram_type": _gram(packets.second_gram_type),
                "realization": (
                    "y=(2theta+theta^2)(E-D)-D for ordered minimal "
                    "D,E in one K12 factor with D.E=-1"
                ),
                "minimal_vectors_per_factor": packets.minimal_vectors_per_factor,
                "neighbors_per_minimal_vector": packets.minimal_neighbors_inner_product_minus_one,
                "wall_count": packets.second_wall_count,
                "first_second_pairing_maximum_bound": _value(
                    packets.first_second_maximum_bound
                ),
                "same_height_pairing_maximum_bound": _value(
                    packets.second_same_height_maximum_bound
                ),
            },
        },
        "conclusion": {
            "new_compact_example": False,
            "new_finite_volume_example": False,
            "new_arithmetic_host_and_packets": True,
            "scope": (
                "exact H=100 necessary-packet screen for 19 coefficients and "
                "complete realization/acuteness of the first two packets of "
                "the selected 2theta lane; later walls remain open"
            ),
        },
    }


def render_report(payload: dict[str, object]) -> str:
    summary = payload["screen_summary"]
    selected = payload["selected_lane"]
    return f"""# Cubic-seven `K12^2` screen through H=100

## Outcome

No compact or finite-volume H24 chamber was found.  The search nevertheless
produced a genuinely new spatial host and two complete, very large Vinberg
packets for

```text
<-2*theta> + K12 + K12.
```

The first packet has **{selected['first_packet']['wall_count']:,}** walls and
the second has **{selected['second_packet']['wall_count']:,}** walls.  Both
are simultaneous acute packets.  This is a bounded prefix, not a reflective
or finite-volume conclusion.

## Finite coefficient screen

The manifest contains all squarefree support products generated by inert
`(2)`, discriminant `(3)`, ramified `(2-theta)`, and at most one of the three
split primes above 13, subject to radical norm at most 351.  It contains
{summary['seed_count']} coefficient squareclasses.  Through exact Vinberg
height 100, {summary['empty_necessary_seed_count']} have no admissible
minimum/congruence Gram type at all.  Only these three lanes survive the
necessary screen:

```text
{chr(10).join(summary['nonempty_seed_ids'])}
```

Emptiness here is a proof for the stated window.  Nonemptiness is only a
necessary condition unless an actual K12 realization is supplied.

## Selected `a=2*theta` packets

Put `s=2theta+theta^2`.  The first control is
`1+2theta+theta^2`, at height `5+14theta+6theta^2`.  Every wall is

```text
(1+2theta+theta^2, s*M),
```

where `M` is one of the 756 minimal vectors in either K12 factor.  Hence the
packet has `2*756=1,512` walls.

The second control is `1+4theta+2theta^2`, at height
`16+45theta+20theta^2`.  Its walls are indexed by ordered minimal pairs
`D.E=-1` in either factor:

```text
y = s*(E-D)-D.
```

PARI independently gives 160 such neighbors per fixed minimal vector, so the
packet has `2*756*160=241,920` walls.  Complete necessary Gram comparisons
give the exact upper pairings `-(3+7theta+3theta^2)` against packet one and
`-1` within packet two; both are strictly acute.

## Exact local scope

For primitive `r=(x,y)`, reflection integrality gives `q(r)|2I(r)`, where
`I=(a*x,(y,K12^2))`.  The Smith form `1^12,3^12`, evenness, and the
squarefree coefficient support leave precisely the archived canonical norm
products.  A factor 3 forces the spatial pairing-radical gate, whose minimum
is 12; an odd supported coefficient prime forces spatial divisibility by
that prime.  The archive records every canonical norm, control, target, and
surviving necessary Gram type.

## Reproducibility

- Backend: `src/coxeter_polytopes/search/cubic7_coxeter_todd_squared.py`
- Audit: `scripts/audit_cubic7_coxeter_todd_squared.py`
- Tests: `tests/test_cubic7_coxeter_todd_squared.py`
- Archive: `data/enumerations/cubic7_coxeter_todd_squared_h100.json`
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=JSON_PATH)
    parser.add_argument("--report", type=Path, default=REPORT_PATH)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    payload = build_payload()
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    report = render_report(payload)
    if args.verify:
        if args.output.read_text() != encoded or args.report.read_text() != report:
            raise SystemExit("frozen cubic-seven K12^2 artifact drifted")
        print("PASS")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(encoded)
    args.report.write_text(report)
    print(json.dumps(payload["screen_summary"], indent=2))


if __name__ == "__main__":
    main()
