#!/usr/bin/env python3
"""Freeze the repaired cubic-seven K12^2 continuation through H169."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from coxeter_polytopes.exact.cubic7 import Cubic7
from coxeter_polytopes.search.cubic7_coxeter_todd_squared import (
    LOCAL_NORM_PROOF,
    POST_H33_CUTOFF_HEIGHT,
    POST_H33_REPLACEMENT_SCALE,
    REPLACEMENT_WITNESS_SCALE,
    dual_witness_norm,
    post_h33_directional_audit,
    post_h33_threat_audit,
    post_h33_threat_root,
)


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCY = ROOT / "data/enumerations/cubic7_coxeter_todd_squared_h33_independent_hostile_audit.json"
JSON_PATH = ROOT / "data/enumerations/cubic7_coxeter_todd_squared_h169.json"
REPORT_PATH = ROOT / "results/cubic7_coxeter_todd_squared_h169.md"


def value(x: Cubic7) -> dict[str, object]:
    return {"exact": x.to_json(), "display": str(x), "signs": [x.sign(i) for i in range(3)]}


def build_payload() -> dict[str, object]:
    rows = post_h33_directional_audit()
    levels = []
    for row in rows:
        level = row["level"]
        levels.append(
            {
                "height": value(level.height),
                "control": value(level.control),
                "norm": value(level.norm),
                "target": value(level.target),
                "types": [
                    {
                        "index": item["index"],
                        "gram": item["gram"],
                        "necessary_direction": value(item["necessary_direction"]),
                        "necessary_arg": item["necessary_arg"],
                        "necessary_count": item["necessary_count"],
                        "actual_direction": value(item["actual_direction"]),
                        "actual_arg": item["actual_arg"],
                        "shell_realization_count": item["shell_realization_count"],
                        "seven_fifths_bound": value(item["seven_fifths_bound"]),
                        "eleven_eighths_bound": value(item["eleven_eighths_bound"]),
                    }
                    for item in row["types"]
                ],
            }
        )
    threat = levels[-1]["types"][35]
    threat_audit = post_h33_threat_audit()
    if threat["actual_arg"] != (9, 23, 10) or threat["seven_fifths_bound"]["display"] != "3/5 + 1/5*theta":
        raise AssertionError("the first realized 7/5 threat changed")
    replacement_norm = dual_witness_norm(POST_H33_REPLACEMENT_SCALE)
    if replacement_norm.sign(0) <= 0:
        raise AssertionError("the 11/8 replacement ceased to be hyperideal")
    return {
        "schema": "coxeter.cubic7_k12_squared_h169.v1",
        "dependency": {
            "path": str(DEPENDENCY.relative_to(ROOT)),
            "sha256": hashlib.sha256(DEPENDENCY.read_bytes()).hexdigest(),
        },
        "host": "<-2*theta> + K12 + K12",
        "local_norm_gate": {"canonical_norms": [2, 6], "proof": LOCAL_NORM_PROOF},
        "scope": "Complete controls and necessary coefficient-Gram types after H33 through H169; exact typewise dual-direction bounds and complete K12 shell realizations for every potentially nonnegative type. Not Vinberg termination.",
        "old_scale": str(REPLACEMENT_WITNESS_SCALE),
        "cutoff_height": value(POST_H33_CUTOFF_HEIGHT),
        "nonempty_levels": levels,
        "first_realized_threat": {
            "height": value(POST_H33_CUTOFF_HEIGHT),
            "type_index": 35,
            "gram": threat["gram"],
            "h_pairings": threat["actual_arg"],
            "root_vector": [entry.to_json() for entry in post_h33_threat_root()],
            "norm": value(threat_audit["norm"]),
            "seven_fifths_pairing": threat["seven_fifths_bound"],
            "pairing_with_explicit_H33_cutter": value(
                threat_audit["h33_cutter_pairing"]
            ),
            "status": "realized_threat_not_claimed_accepted",
            "reason": "A realization exists in K12, but full chronological Vinberg acceptance against all intervening realized packets is not certified here.",
        },
        "replacement": {
            "scale": str(POST_H33_REPLACEMENT_SCALE),
            "norm": value(replacement_norm),
            "strict_against_every_necessary_type_through_cutoff": True,
            "worst_H169_type_35_bound": levels[-1]["types"][35]["eleven_eighths_bound"],
        },
        "conclusion": {
            "meaningful_exact_cutoff": True,
            "bounded_prefix_has_infinite_volume": True,
            "accepted_new_wall": False,
            "compact_example": False,
            "finite_volume_example": False,
            "vinberg_terminated": False,
        },
    }


def report(p: dict[str, object]) -> str:
    return f"""# Cubic-seven K12-squared continuation through H169

The audited `7/5` ray first meets a **realized** nonnegative typewise threat at

```text
H = 169 + 474 theta + 211 theta^2,
x = 4 + 13 theta + 6 theta^2,
q = 2.
```

Necessary type 35 has Gram
`[[66,162,69],[162,406,175],[69,175,76]]`.  A complete signed K12
minimum-shell census realizes its maximizing pairings `(h.U,h.V,h.W)=(9,23,10)`;
the `7/5` pairing is exactly `3/5+theta/5>0`.

This is not called an accepted Vinberg wall: chronological acceptance against
every intervening realized packet remains uncertified.  Instead the smaller
ray `z''=(1,(11/8)h,0)` gives the bounded conclusion.  Its norm is
`121/48-2theta>0`, and exact augmented-Gram maxima are strictly negative for
all necessary types through H169; at the dangerous type the bound is
`3/8-3theta/8-theta^2/4<0`.

The chronology contains seven nonempty levels after H33 and 98 necessary types,
including 19 ternary-radical types at H156.  The latter are enumerated by the
exact scaled gate `Gram/3`, `minimum 12/3=4`.
At H81 a relaxed positive maximum is not realized: a complete K12 shell
census lowers it from `6+16theta+7theta^2` to
`5+12theta+5theta^2`.  At H156 its two relaxed positive types are likewise
lowered by the complete ternary-radical shell census.  At H169 type 31 is lowered to
`7+18theta+8theta^2`; type 35 realizes its relaxed maximum.

This proves only that the complete necessary prefix through H169 still has a
strict hyperideal ray.  It is not a compactness, finite-volume, or Vinberg
termination claim.
"""


def main() -> None:
    parser = argparse.ArgumentParser(); parser.add_argument("--verify", action="store_true"); args = parser.parse_args()
    payload = build_payload(); encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"; rendered = report(payload)
    if args.verify:
        if JSON_PATH.read_text() != encoded or REPORT_PATH.read_text() != rendered:
            raise SystemExit("H169 artifact drift")
        print("PASS"); return
    JSON_PATH.write_text(encoded); REPORT_PATH.write_text(rendered); print(JSON_PATH); print(REPORT_PATH)


if __name__ == "__main__": main()
