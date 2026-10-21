#!/usr/bin/env python3
"""Freeze the first cutter and replacement for silver DualExtremal(28,3)."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_dual_extremal_28_3 import (
    FIRST_CUTTING_CONTROL,
    FIRST_CUTTING_D,
    FIRST_CUTTING_HEIGHT,
    FIRST_CUTTING_NORM,
    FIRST_CUTTING_TRIPLE,
    FIRST_CUTTING_U,
    FIRST_CUTTING_V,
    FIRST_NECESSARY_HEIGHT,
    EARLIER_ACCEPTANCE_SPECS,
    REPLACEMENT_WITNESS_SCALE,
    RHO,
    controls,
    earlier_acceptance_bound,
    first_cutting_old_witness_pairing,
    first_cutting_root_audit,
    necessary_triples,
    norm_eight_pairing_maximum,
    replacement_margin,
    replacement_witness_norm,
)


ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = ROOT / "data/enumerations/silver_dual_extremal_28_3_h362_replacement.json"
REPORT_PATH = ROOT / "results/silver_dual_extremal_28_3_h362_replacement.md"
FIRST_PREFIX = ROOT / "data/enumerations/silver_dual_extremal_28_3_first_prefix.json"


def value(value: Silver) -> dict[str, object]:
    return {"exact": value.to_json(), "display": str(value), "signs": [value.sign(0), value.sign(1)]}


def chronology_digest(rows) -> str:
    stream = "".join(
        "|".join(
            (
                str(row.height.a), str(row.height.b), str(row.norm_class.norm.a),
                str(row.norm_class.norm.b), str(row.control.a), str(row.control.b),
                str(row.target.a), str(row.target.b), str(int(row.norm_class.dyadic)),
                str(int(row.norm_class.p7)), str(int(row.norm_class.discriminant_three)),
            )
        ) + "\n"
        for row in rows
    )
    return hashlib.sha256(stream.encode("ascii")).hexdigest()


def build_payload() -> dict[str, object]:
    if not FIRST_PREFIX.exists() or not first_cutting_root_audit():
        raise AssertionError("the first-prefix dependency or explicit cutter failed")
    rows = controls(FIRST_CUTTING_HEIGHT)
    nonempty = tuple((row, necessary_triples(row)) for row in rows if necessary_triples(row))
    expected = (
        (Silver(93, 66), Silver(3), 1),
        (Silver(93, 66), Silver(6), 5),
        (Silver(124, 88), Silver(2), 1),
        (Silver(181, 128), Silver(2), 6),
        (Silver(186, 132), Silver(6), 4),
        (Silver(362, 256), Silver(2), 8),
    )
    observed = tuple((row.height, row.norm_class.norm, len(packet)) for row, packet in nonempty)
    if observed != expected or len(rows) != 103 or len({row.height for row in rows}) != 70:
        raise AssertionError((len(rows), observed))

    levels = []
    first_nonpositive = None
    for row, packet in nonempty:
        types = []
        for triple in packet:
            maximum = norm_eight_pairing_maximum(triple)
            old_margin = RHO * row.control - maximum.maximum / 2
            new_margin = replacement_margin(row, triple)
            record = {
                "gram": [triple.u_square, triple.v_square, triple.u_dot_v],
                "maximum_d_dot_y": value(maximum.maximum),
                "maximizing_pair": list(maximum.maximizing_pair),
                "augmented_determinant": maximum.augmented_determinant,
                "retained_local_extensions": maximum.retained_extensions,
                "old_witness_margin": value(old_margin),
                "replacement_margin": value(new_margin),
            }
            types.append(record)
            if old_margin.sign(0) <= 0 and first_nonpositive is None:
                first_nonpositive = (row, triple, maximum, old_margin)
            if new_margin.sign(0) <= 0:
                raise AssertionError("replacement lost strictness")
        levels.append(
            {
                "height": value(row.height), "mirror_norm": value(row.norm_class.norm),
                "time_control": value(row.control), "target": value(row.target),
                "necessary_gram_types": types,
            }
        )
    if first_nonpositive is None:
        raise AssertionError("no cutter appeared")
    row, triple, maximum, old_margin = first_nonpositive
    if (
        row.height != FIRST_CUTTING_HEIGHT or row.norm_class.norm != FIRST_CUTTING_NORM
        or row.control != FIRST_CUTTING_CONTROL or triple != FIRST_CUTTING_TRIPLE
        or maximum.maximum != Silver(53, 37) or maximum.maximizing_pair != (53, 37)
        or maximum.augmented_determinant != 0
        or old_margin != Silver(-1, -1) / 2
    ):
        raise AssertionError("first threat changed")
    if first_cutting_old_witness_pairing() != Silver(1, 1) / 2:
        raise AssertionError("old-witness cutter pairing changed")
    if replacement_witness_norm() != Silver(1777) / 625 - Silver(0, 2):
        raise AssertionError("replacement norm changed")

    acceptance = []
    for spec in EARLIER_ACCEPTANCE_SPECS:
        bound = earlier_acceptance_bound(spec)
        if bound.admitted_maximum.sign(0) > 0:
            raise AssertionError("the cutter is not acute to an earlier necessary type")
        acceptance.append(
            {
                "height": value(spec.height),
                "norm": value(spec.norm),
                "original_gram": [spec.triple.u_square, spec.triple.v_square, spec.triple.u_dot_v],
                "reduced_gram": [list(row) for row in spec.reduced_gram],
                "reduced_field_coefficients": [str(spec.first_coefficient), str(spec.second_coefficient)],
                "cross_pairing_modulus": spec.pairing_modulus,
                "feasible_psd_cross_grams": bound.feasible_cross_gram_count,
                "relaxed_maximum_root_pairing": value(bound.relaxed_maximum),
                "relaxed_maximizing_entries": list(bound.relaxed_maximizing_entries),
                "positive_relaxed_extensions": bound.positive_relaxed_count,
                "forbidden_positive_extensions": [
                    {
                        "cross_entries": list(entries),
                        "root_pairing": value(pairing),
                        "short_combination": list(coefficients),
                        "short_square": square,
                    }
                    for entries, pairing, coefficients, square in bound.forbidden_positive_extensions
                ],
                "admitted_maximum_root_pairing": value(bound.admitted_maximum),
                "admitted_maximizing_entries": list(bound.admitted_maximizing_entries),
            }
        )
    if sum(row["positive_relaxed_extensions"] for row in acceptance) != 2:
        raise AssertionError("the positive cross-Gram exception census changed")
    return {
        "schema_version": "0.1",
        "host": "<-(2*sqrt(2)-1)> + DualExtremal(28,3)",
        "hyperbolic_dimension": 28,
        "scope": "Complete controls and exact necessary Gram chronology through the first accepted cutter H=362+256*sqrt(2), followed by a strict same-direction replacement. Not Vinberg termination.",
        "first_prefix_dependency": {"path": str(FIRST_PREFIX.relative_to(ROOT)), "sha256": hashlib.sha256(FIRST_PREFIX.read_bytes()).hexdigest()},
        "chronology": {
            "height_bound": value(FIRST_CUTTING_HEIGHT),
            "control_count": len(rows),
            "distinct_height_count": len({row.height for row in rows}),
            "sha256": chronology_digest(rows),
            "binary_minimum_gate": {
                "method": "Exact determinant bounds exhaust every coefficient pair capable of norm below the host/S3 minimum; rank-one forms use gcd(u^2,v^2,u.v).",
                "repaired_former_fixed_box_gate": True,
                "removed_spurious_types": [
                    {"height": "124+88*sqrt(2)", "gram": [118, 66, 88], "short_square": 6},
                    {"height": "181+128*sqrt(2)", "gram": [164, 100, 128], "short_square": 4},
                    {"height": "186+132*sqrt(2)", "gram": [534, 294, 396], "short_square": 6},
                    {"height": "248+176*sqrt(2)", "gram": [242, 128, 176], "short_square": 2},
                    {"height": "362+256*sqrt(2)", "gram": [342, 192, 256], "short_square": 6},
                    {"height": "362+256*sqrt(2)", "gram": [338, 194, 256], "short_square": 2},
                ],
            },
            "nonempty_levels": levels,
        },
        "first_cutting_wall": {
            "norm": value(FIRST_CUTTING_NORM), "time_control": value(FIRST_CUTTING_CONTROL),
            "gram": [FIRST_CUTTING_TRIPLE.u_square, FIRST_CUTTING_TRIPLE.v_square, FIRST_CUTTING_TRIPLE.u_dot_v],
            "d": list(FIRST_CUTTING_D), "u": list(FIRST_CUTTING_U), "v": list(FIRST_CUTTING_V),
            "relation": "u=2*d+v", "d_pairings": [53, 37],
            "old_witness_pairing": value(first_cutting_old_witness_pairing()),
            "reflective_and_primitive": True, "accepted": True,
            "earlier_packet_cross_gram_bounds": acceptance,
            "same_height_acuteness": "For distinct q=2 roots, g in O with 0<g<=2 and |g'|<=2 gives g in {1,sqrt(2),2}; g=1 or sqrt(2) would make the coefficient-vector difference have rational trace part 2 or 4 below min(L)=8, and g=2 is equality.",
        },
        "replacement": {
            "ray": "z'=(1,(12/25)d)", "scale": str(REPLACEMENT_WITNESS_SCALE),
            "norm": value(replacement_witness_norm()),
            "tight_margin": value(replacement_margin(row, triple)),
            "strict_against_every_necessary_type_through_bound": True,
        },
        "conclusion": {"old_witness_is_cut": True, "accepted_cutter_found": True, "prefix_through_bound_has_infinite_volume": True, "compact_example": False, "finite_volume_example": False},
    }


def render(payload: dict[str, object]) -> str:
    c = payload["chronology"]
    r = payload["replacement"]
    return f"""# Silver DualExtremal(28,3): first cutter at H362 and replacement

The complete eight-class chronology through `H={c['height_bound']['display']}` has
{c['control_count']} controls on {c['distinct_height_count']} heights. After the
audited first packet, the only nonempty necessary levels before the endpoint are
`124+88*sqrt(2)`, `181+128*sqrt(2)`, and `186+132*sqrt(2)`.
The continuation also repairs the former fixed `[-3,3]` binary minimum gate:
exact determinant coefficient bounds remove six spurious types, including the
otherwise nonempty-looking level `248+176*sqrt(2)`.

At the endpoint, the q=2 type `(362,182,256)` has exact augmented maximum
`d.(u+sqrt(2)v)=53+37*sqrt(2)` with determinant zero. The explicit vectors in
the JSON satisfy `u=2d+v`, so this is a realized primitive reflective root. Its
pairing with `(1,d/2)` is `{payload['first_cutting_wall']['old_witness_pairing']['display']}>0`.
Exact four-vector PSD cross-Gram bounds, with the two sole positive relaxed
extensions rejected by explicit norm-6 and norm-2 combinations, make this root
acute to every earlier necessary type. The elementary equal-height norm-two
difference lemma handles its simultaneous packet.

The cutter does not close the cone. The replacement `{r['ray']}` has exact norm
`{r['norm']['display']}>0`; its tight necessary-type margin is
`{r['tight_margin']['display']}>0`. Thus the accepted prefix through H362 still
has an open hyperideal patch and infinite volume.

This is a bounded result, not Vinberg termination or global nonreflectivity.
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    payload = build_payload()
    artifact = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()
    report = render(payload).encode()
    if args.verify:
        if JSON_PATH.read_bytes() != artifact or REPORT_PATH.read_bytes() != report:
            raise SystemExit("frozen H362 certificate differs")
        print("PASS", hashlib.sha256(artifact).hexdigest())
        return
    JSON_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    JSON_PATH.write_bytes(artifact)
    REPORT_PATH.write_bytes(report)
    print(hashlib.sha256(artifact).hexdigest())


if __name__ == "__main__":
    main()
