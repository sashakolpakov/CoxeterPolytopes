#!/usr/bin/env python3
"""Clean-room replay of the silver Lambda22 H10000--H25000 frontier.

This imports only the standalone arithmetic/enumeration kernel from the
independent H10000 replay; it imports no ``coxeter_polytopes`` project code.
The control box and both new obstruction reductions are derived here for the
larger cutoff.  The result is a bounded necessary-superset audit, not a
termination or global nonreflectivity theorem.
"""

from __future__ import annotations

import argparse
from functools import cmp_to_key
from hashlib import sha256
import json
from math import lcm
from pathlib import Path

from independent_hostile_audit_silver_lambda22_h10000_replacement import (
    A,
    K2,
    NORM_GATES,
    NORMS,
    candidate_grams,
    compare_controls,
    compare_k2,
    extensions,
    norm_gate,
    sign_integer_pair,
    Control,
)


ROOT = Path(__file__).resolve().parents[1]
PRODUCTION = ROOT / "data/enumerations/silver_lambda22_h25000_replacement.json"
OUTPUT = (
    ROOT
    / "data/enumerations/silver_lambda22_h25000_replacement_independent_hostile_audit.json"
)
START = K2.make(10_000)
BOUND = K2.make(25_000)
SCALE = K2.make("28/87")
COEFFICIENT_BOX = 561
EXPECTED_TYPE_STREAM = "21240a46b23221f1c167ffd5c6267c71eb328a57b0a067fa70e488b3955f2e5f"
EXPECTED_THREATS = (
    (
        K2.make(6392, 4520),
        K2.make(4),
        K2.make(176, 124),
        (12772, 6400, 9040),
        (226, 160, 0),
        K2.make("64/87", "-44/87"),
    ),
    (
        K2.make(9588, 6780),
        K2.make(6),
        K2.make(264, 186),
        (28734, 14400, 20340),
        (339, 240, 0),
        K2.make("32/29", "-22/29"),
    ),
)


def coefficient_box_proof() -> dict[str, object]:
    """Prove the exhaustive box for ``x=m+n*sqrt(2)`` at H<=25000."""

    # H<=25000 implies A*x_0^2<=25000*q_0.  The uniform bound |x_0|<1118
    # has strict slack for every canonical mirror norm.  Target positivity
    # at the conjugate place and q+9A<0 give |x_1|<3.
    for norm in NORMS:
        if norm.sign(0) <= 0 or norm.sign(1) <= 0:
            raise AssertionError("a mirror norm ceased to be totally positive")
        if (K2.make(1118 * 1118) * A - BOUND * norm).sign(0) <= 0:
            raise AssertionError("the identity control bound is insufficient")
        if (norm + K2.make(9) * A).sign(1) >= 0:
            raise AssertionError("the conjugate control bound is insufficient")
    # |m|=|x0+x1|/2<1121/2=560.5, so |m|<=560.  Also
    # |n|sqrt(2)<560.5 and 397sqrt(2)>560.5, so |n|<=396.
    if (K2.make(0, 397) - K2.make("1121/2")).sign(0) <= 0:
        raise AssertionError("the radical coefficient bound was not proved")
    return {
        "identity_embedding_absolute_bound": "<1118",
        "conjugate_embedding_absolute_bound": "<3",
        "derived_rational_coefficient_bound": "<=560",
        "derived_radical_coefficient_bound": "<=396",
        "enumerated_symmetric_coefficient_box": COEFFICIENT_BOX,
    }


def enumerate_controls() -> tuple[Control, ...]:
    """Independently reconstruct every control in 10000<H<=25000."""

    coefficient_box_proof()
    rows: set[Control] = set()
    for norm in NORMS:
        qa, qb = int(norm.a), int(norm.b)
        qnorm = qa * qa - 2 * qb * qb
        for rational in range(-COEFFICIENT_BOX, COEFFICIENT_BOX + 1):
            for radical in range(-COEFFICIENT_BOX, COEFFICIENT_BOX + 1):
                if sign_integer_pair(rational, radical) <= 0:
                    continue
                # 2*A*x=c+d*sqrt(2).  Exact division by q lies in Z[sqrt2]
                # iff both numerator coefficients are divisible by N(q).
                c = -2 * rational + 4 * radical
                d = 2 * rational - 2 * radical
                if (c * qa - 2 * d * qb) % qnorm:
                    continue
                if (d * qa - c * qb) % qnorm:
                    continue
                target_a = (
                    qa
                    - rational * rational
                    - 2 * radical * radical
                    + 4 * rational * radical
                )
                target_b = (
                    qb
                    + rational * rational
                    + 2 * radical * radical
                    - 2 * rational * radical
                )
                if sign_integer_pair(target_a, target_b, 1) < 0:
                    continue
                height_a = target_a - qa
                height_b = target_b - qb
                if sign_integer_pair(height_a - 10_000 * qa, height_b - 10_000 * qb) <= 0:
                    continue
                if sign_integer_pair(height_a - 25_000 * qa, height_b - 25_000 * qb) > 0:
                    continue
                rows.add(
                    Control(
                        norm,
                        K2.make(rational, radical),
                        K2.make(target_a, target_b),
                        K2.make(height_a, height_b) / norm,
                    )
                )
    return tuple(sorted(rows, key=cmp_to_key(compare_controls)))


def ray_margin(control: K2, p: int, q: int) -> K2:
    return SCALE * K2.make(p, q) - A * control


def scan_frontier() -> dict[str, object]:
    controls = enumerate_controls()
    controls_by_norm: dict[str, int] = {}
    retained_by_norm: dict[str, int] = {}
    extensions_by_norm: dict[str, int] = {}
    retained = extension_count = 0
    threats = []
    digest = sha256()
    for index, control in enumerate(controls, start=1):
        label = control.norm.display()
        controls_by_norm[label] = controls_by_norm.get(label, 0) + 1
        u_modulus, v_modulus = NORM_GATES[control.norm]
        for gram in candidate_grams(control.target):
            u_square, v_square, dot = gram
            if (
                not norm_gate(u_square, u_modulus)
                or not norm_gate(v_square, v_modulus)
                or dot % lcm(u_modulus, v_modulus)
            ):
                continue
            rows = extensions(gram, u_modulus, v_modulus)
            if not rows:
                continue
            retained += 1
            extension_count += len(rows)
            retained_by_norm[label] = retained_by_norm.get(label, 0) + 1
            extensions_by_norm[label] = extensions_by_norm.get(label, 0) + len(rows)
            digest.update(
                (
                    f"{control.norm.a},{control.norm.b};"
                    f"{control.height.a},{control.height.b};"
                    f"{control.x.a},{control.x.b};"
                    f"{u_square},{v_square},{dot};{len(rows)}\n"
                ).encode()
            )
            nonnegative = tuple(
                (p, q, determinant, ray_margin(control.x, p, q))
                for p, q, determinant in rows
                if ray_margin(control.x, p, q).sign(0) >= 0
            )
            if nonnegative:
                threats.append(
                    (
                        control.height,
                        control.norm,
                        control.x,
                        gram,
                        nonnegative,
                    )
                )
        if index % 50 == 0:
            print(
                "PROGRESS",
                index,
                len(controls),
                retained,
                extension_count,
                len(threats),
                flush=True,
            )
    if (
        len(controls),
        retained,
        extension_count,
        len(threats),
        digest.hexdigest(),
    ) != (301, 5_575, 2_663_705, 2, EXPECTED_TYPE_STREAM):
        raise AssertionError("the independent H25000 frontier changed")
    for replay, expected in zip(threats, EXPECTED_THREATS, strict=True):
        height, norm, control, gram, nonnegative = replay
        eh, en, ec, eg, extension, margin = expected
        if (height, norm, control, gram, nonnegative) != (
            eh,
            en,
            ec,
            eg,
            ((*extension, margin),),
        ):
            raise AssertionError("an independent threat row changed")
    return {
        "control_count": len(controls),
        "controls_by_norm": controls_by_norm,
        "first_control_height": controls[0].height.exact(),
        "last_control_height": controls[-1].height.exact(),
        "retained_type_count": retained,
        "retained_types_by_norm": retained_by_norm,
        "locally_admissible_extension_count_with_multiplicity": extension_count,
        "locally_admissible_extensions_by_norm": extensions_by_norm,
        "retained_type_stream_sha256": digest.hexdigest(),
        "nonnegative_type_count": len(threats),
        "threats": [
            {
                "height": height.exact(),
                "norm": norm.exact(),
                "control": control.exact(),
                "gram": list(gram),
                "nonnegative_extensions": [
                    {
                        "d_pairings": [p, q],
                        "extension_determinant": determinant,
                        "margin": margin.exact(),
                    }
                    for p, q, determinant, margin in nonnegative
                ],
            }
            for height, norm, control, gram, nonnegative in threats
        ],
    }


def critical_reductions() -> list[dict[str, object]]:
    output = []
    parameters = ((40, 56, 4, 2, 12), (60, 84, 6, 3, 24))
    for expected, (multiple, shift, short_norm, modulus, minimum) in zip(
        EXPECTED_THREATS, parameters, strict=True
    ):
        _, _, _, gram, extension, margin = expected
        p, q, determinant = extension
        u_square, v_square, dot = gram
        relation_norm = v_square + multiple * multiple * 4 - 2 * multiple * q
        relation_cross = dot - multiple * p
        short_square = u_square - 2 * shift * p + shift * shift * 4
        if (
            determinant,
            relation_norm,
            relation_cross,
            short_square,
            margin.sign(0),
        ) != (0, 0, 0, short_norm, 1):
            raise AssertionError("an independent short-vector reduction changed")
        output.append(
            {
                "forced_relation": f"v={multiple}d",
                "relation_norm": relation_norm,
                "relation_cross": relation_cross,
                "short_vector": f"u-{shift}d",
                "short_vector_norm": short_square,
                "gate_sublattice": f"S{modulus}",
                "gate_sublattice_minimum": minimum,
                "empty": True,
            }
        )
    return output


def compare_production(frontier: dict[str, object]) -> dict[str, object]:
    production = json.loads(PRODUCTION.read_text())
    frozen = production["frontier"]
    comparisons = {
        "control_count": "control_count",
        "controls_by_norm": "controls_by_norm",
        "retained_type_count": "retained_gate_compatible_type_count",
        "retained_types_by_norm": "retained_types_by_norm",
        "locally_admissible_extension_count_with_multiplicity": (
            "locally_admissible_extension_count_with_multiplicity"
        ),
        "locally_admissible_extensions_by_norm": "locally_admissible_extensions_by_norm",
        "retained_type_stream_sha256": "retained_type_stream_sha256",
        "nonnegative_type_count": "nonnegative_type_count",
    }
    for replay_key, frozen_key in comparisons.items():
        if frontier[replay_key] != frozen[frozen_key]:
            raise AssertionError(f"production disagrees at {replay_key}")
    return {
        "production_artifact": str(PRODUCTION.relative_to(ROOT)),
        "production_sha256": sha256(PRODUCTION.read_bytes()).hexdigest(),
        "all_frontier_counts_and_digest_agree": True,
    }


def build_payload() -> dict[str, object]:
    frontier = scan_frontier()
    source = compare_production(frontier)
    reductions = critical_reductions()
    return {
        "schema_version": "0.1.0",
        "verdict": "PASS",
        "scope": (
            "Independent exact necessary-superset replay strictly after rational "
            "identity height 10000 through inclusive height 25000."
        ),
        "independence": (
            "No project imports. Reuses only the standalone Q(sqrt2), Gram, and "
            "ellipse kernel of the independently audited H10000 replay; derives "
            "the larger control box and new obstruction reductions locally."
        ),
        "source": source,
        "cutoff_semantics": {
            "lower": "strictly greater than the rational field element 10000",
            "upper": "less than or equal to the rational field element 25000",
            "coefficient_box_proof": coefficient_box_proof(),
        },
        "frontier": frontier,
        "critical_reductions": reductions,
        "logical_scope": {
            "all_controls_reconstructed": True,
            "all_binary_Gram_types_reconstructed": True,
            "all_integral_PSD_extensions_enumerated": True,
            "both_nonnegative_types_empty": True,
            "bounded_prefix_only": True,
            "vinberg_termination": False,
            "global_nonreflectivity": False,
        },
        "conclusion": (
            "The production H25000 persistence certificate passes the clean-room "
            "replay. The two nonnegative necessary types are empty by the S2 and "
            "S3 minimum contradictions, so the H10000 ray remains strict through "
            "the stated cutoff."
        ),
    }


def encoded_payload() -> bytes:
    return (json.dumps(build_payload(), indent=2, sort_keys=True) + "\n").encode()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    arguments = parser.parse_args()
    encoded = encoded_payload()
    if arguments.verify:
        if not OUTPUT.exists() or OUTPUT.read_bytes() != encoded:
            raise SystemExit("FAIL: silver H25000 hostile replay differs from frozen JSON")
        status = "PASS"
    else:
        OUTPUT.write_bytes(encoded)
        status = "WROTE"
    payload = json.loads(encoded)
    print(
        json.dumps(
            {
                "status": status,
                "output": str(OUTPUT),
                "sha256": sha256(encoded).hexdigest(),
                "controls": payload["frontier"]["control_count"],
                "types": payload["frontier"]["retained_type_count"],
                "extensions": payload["frontier"][
                    "locally_admissible_extension_count_with_multiplicity"
                ],
                "threats": payload["frontier"]["nonnegative_type_count"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
