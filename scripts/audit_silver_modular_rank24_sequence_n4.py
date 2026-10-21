#!/usr/bin/env python3
"""Freeze the rank-24 modular-sequence identity and complete n=4 prefix."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sympy import Matrix
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_modular_rank24_sequence import (
    L24_11_GRAM_SHA256,
    L24_11_SOURCE_URL,
    N4_COEFFICIENT,
    N4_DUAL_INDEX,
    N4_DUAL_SQUARE,
    N4_HEIGHT_BOUND,
    N4_NORM_CLASSES,
    N4_PI,
    N4_PI_BAR,
    N4_WITNESS_SCALE,
    P24_9_GRAM_SHA256,
    P24_9_SOURCE_URL,
    endpoint_grams,
    family_matrix_identity,
    level,
    minimum,
    n4_first_packet_audit,
    n4_radical_basis,
    n4_scaled_radical_gram,
    sequence_gram,
)


ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = ROOT / "data/enumerations/silver_modular_rank24_sequence_n4.json"
REPORT_PATH = ROOT / "results/silver_modular_rank24_sequence_n4.md"


def value(item: Silver) -> dict[str, object]:
    return {"exact": item.to_json(), "display": str(item), "signs": [item.sign(0), item.sign(1)]}


def matrix_sha256(matrix) -> str:
    stream = "\n".join(" ".join(map(str, row)) for row in matrix) + "\n"
    return hashlib.sha256(stream.encode("ascii")).hexdigest()


def bareiss_determinant(matrix) -> int:
    work = [list(map(int, row)) for row in matrix]
    previous, sign = 1, 1
    for column in range(len(work) - 1):
        pivot = next((row for row in range(column, len(work)) if work[row][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign = -sign
        value_at_pivot = work[column][column]
        for row in range(column + 1, len(work)):
            for index in range(column + 1, len(work)):
                numerator = (
                    work[row][index] * value_at_pivot
                    - work[row][column] * work[column][index]
                )
                if numerator % previous:
                    raise AssertionError("Bareiss division ceased to be exact")
                work[row][index] = numerator // previous
        previous = value_at_pivot
        for row in range(column + 1, len(work)):
            work[row][column] = 0
    return sign * work[-1][-1]


def build_payload() -> dict[str, object]:
    if not family_matrix_identity():
        raise AssertionError("the modular-sequence matrix identity failed")
    g0, g1 = endpoint_grams()
    if matrix_sha256(g0) != P24_9_GRAM_SHA256 or matrix_sha256(g1) != L24_11_GRAM_SHA256:
        raise AssertionError("an endpoint catalogue Gram changed")
    g4 = sequence_gram(4)
    smith = smith_normal_form(Matrix(g4), domain=ZZ)
    smith_diagonal = tuple(abs(int(smith[i, i])) for i in range(24))
    if smith_diagonal != (1,) * 12 + (41,) * 12:
        raise AssertionError("the G4 Smith invariants changed")
    if bareiss_determinant(g4) != 41**12:
        raise AssertionError("the G4 determinant changed")
    radical_basis = n4_radical_basis()
    if abs(bareiss_determinant(radical_basis)) != 41**12:
        raise AssertionError("the mod-41 radical preimage index changed")
    radical_gram = n4_scaled_radical_gram()
    if bareiss_determinant(radical_gram) != 41**12:
        raise AssertionError("the scaled radical determinant changed")

    audit = n4_first_packet_audit()
    if audit.control_count != 196 or audit.distinct_height_count != 132 or not audit.strict_hyperideal:
        raise AssertionError("the complete n=4 bounded prefix changed")
    if tuple(maximum.maximum for maximum in audit.maxima) != (Silver(17, 12), Silver(14, 10)):
        raise AssertionError("the n=4 directional maxima changed")

    return {
        "schema_version": "0.1",
        "scope": (
            "Exact rank-24 Gram-family identity and complete eight-class local control "
            "chronology for the n=4 host through real height 5000. This is a bounded "
            "infinite-volume prefix theorem, not Vinberg termination or global nonreflectivity."
        ),
        "sources": [
            {"name": "P24.9", "url": P24_9_SOURCE_URL, "literal_gram_sha256": P24_9_GRAM_SHA256},
            {"name": "L_24.11", "url": L24_11_SOURCE_URL, "literal_gram_sha256": L24_11_GRAM_SHA256},
        ],
        "family": {
            "definition": "G_n=G_0+n*(G_1-G_0)",
            "matrix_identity": "D*G_0^-1*D=D-G_0/5",
            "determinant_formula": "det(G_n)=(n^2+5*n+5)^12",
            "minimum_pattern_checked_separately": "min(G_n)=8+4*n",
        },
        "n4_spatial_lattice": {
            "rank": 24,
            "level": level(4),
            "determinant": str(41**12),
            "smith_invariants": list(smith_diagonal),
            "minimum": minimum(4),
            "minimum_shell_signed_count": 37800,
            "qfminim_provenance": "PARI qfminim exact replay recorded by this certificate; shell constants asserted in tests.",
            "dual_coordinate": {"index_zero_based": N4_DUAL_INDEX, "square": str(N4_DUAL_SQUARE)},
        },
        "n4_coefficient": {
            "a": value(N4_COEFFICIENT),
            "factorization": "a=(7+2*sqrt(2))*(sqrt(2)-1)",
            "pi": value(N4_PI),
            "pi_bar": value(N4_PI_BAR),
            "local_norm_classes": [
                {"norm": value(row.norm), "pi_gate": row.pi, "pi_bar_gate": row.pi_bar}
                for row in N4_NORM_CLASSES
            ],
        },
        "radical_gate": {
            "rank_mod_41": 12,
            "integral_preimage_index": str(41**12),
            "scaled_gram_definition": "K^T*G4*K/41",
            "scaled_gram_sha256": matrix_sha256(radical_gram),
            "scaled_determinant": str(41**12),
            "scaled_minimum": 24,
            "scaled_minimum_shell_signed_count": 37800,
            "normalization": "A pi-gated spatial vector is divided by pi, so its binary target is divided by pi^2 (and analogously for pi_bar).",
        },
        "chronology": {
            "height_bound": value(N4_HEIGHT_BOUND),
            "control_count": audit.control_count,
            "distinct_height_count": audit.distinct_height_count,
            "first_nonempty_height": value(audit.first_control.height),
            "first_norm": value(audit.first_control.norm_class.norm),
            "first_time_control": value(audit.first_control.control),
            "first_target": value(audit.first_control.target),
            "necessary_gram_types": [
                [triple.u_square, triple.v_square, triple.u_dot_v] for triple in audit.triples
            ],
        },
        "directional_certificate": {
            "short_shell_coordinate_bounds": {"24": 2, "26": 2, "28": 4},
            "maxima": [
                {
                    "value": value(maximum.maximum),
                    "coordinates": list(maximum.coordinates),
                    "augmented_determinant": str(maximum.augmented_determinant),
                    "retained_coordinate_pairs": maximum.retained_count,
                }
                for maximum in audit.maxima
            ],
            "witness": "z=(1,3*h), h=G4^-1*e_1",
            "witness_scale": str(N4_WITNESS_SCALE),
            "witness_lorentz_norm": value(audit.witness_norm),
            "wall_margins": [value(margin) for margin in audit.margins],
            "strict_hyperideal": True,
        },
        "verdict": {
            "bounded_prefix_has_open_hyperideal_patch": True,
            "bounded_prefix_has_finite_volume": False,
            "compact_example_found": False,
            "eventual_chamber": "OPEN",
        },
    }


def render(payload: dict[str, object]) -> str:
    chronology = payload["chronology"]
    direction = payload["directional_certificate"]
    return f"""# Silver rank-24 modular sequence: the complete `n=4` prefix

The literal Nebe--Sloane endpoint Grams define `G_n=G_0+n(G_1-G_0)`.
Exact rational matrix arithmetic gives
`D G_0^(-1) D = D-G_0/5`, hence
`det(G_n)=(n^2+5n+5)^12`.  For `n=4`, `G_4` has Smith invariants
`1^12,41^12`, minimum 24, and 37,800 signed minimal vectors.

For the admissible matched coefficient `a=-3+5*sqrt(2)`, the complete local
alphabet is the dyadic choice times the two independent split-41 choices:
`{{1,2}}*{{1,7+2*sqrt(2)}}*{{1,7-2*sqrt(2)}}`.  The mod-41 pairing-radical
preimage has index `41^12`; after dividing its Gram by 41 it again has
determinant `41^12`, minimum 24, and 37,800 signed minimal vectors.  This
forces divisor-square normalization on every split-prime gated target.

The complete eight-class chronology through `H={chronology['height_bound']['display']}`
contains {chronology['control_count']} controls on {chronology['distinct_height_count']}
heights.  Its first nonempty necessary packet is the ungated norm-two level
`H={chronology['first_nonempty_height']['display']}`, with the two binary Gram
types `(798,408,570)` and `(794,410,570)`.

For `h=G_4^(-1)e_1`, `h^2=24/41`.  Exact augmented PSD enumeration plus the
complete coordinate bounds on the norm-24, 26, and 28 shells gives directional
maxima `17+12*sqrt(2)` and `14+10*sqrt(2)`.  Therefore
`z=(1,3h)` has Lorentz norm `{direction['witness_lorentz_norm']['display']}>0`
and margins `{direction['wall_margins'][0]['display']}` and
`{direction['wall_margins'][1]['display']}`, both positive.  The whole bounded
prefix retains an open hyperideal patch and has infinite volume.

This is a bounded negative design result, not Vinberg termination or a proof
of global nonreflectivity.  The broader sequence ranking used only the `q=1,2`
necessary relaxation and is intentionally not promoted to a theorem here.

Sources: [P24.9]({P24_9_SOURCE_URL}) and [L_24.11]({L24_11_SOURCE_URL}).
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
            raise SystemExit("frozen silver rank-24 n=4 certificate differs")
        print("PASS", hashlib.sha256(artifact).hexdigest())
        return
    JSON_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    JSON_PATH.write_bytes(artifact)
    REPORT_PATH.write_bytes(report)
    print(hashlib.sha256(artifact).hexdigest())


if __name__ == "__main__":
    main()

