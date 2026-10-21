#!/usr/bin/env python3
"""Freeze the exact H45-prefix audit for the p23 silver odd-Leech host."""

from __future__ import annotations

import json
import shutil
import subprocess
from fractions import Fraction
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_odd_leech import (
    HEIGHT_ZERO_CANONICAL_NORMS,
    HEIGHT_ZERO_TRACE_MINIMUM,
    ODD_LEECH_COEFFICIENT,
    ODD_LEECH_GRAM,
    ODD_LEECH_NORM4_WITNESS_SCALED,
    ODD_LEECH_NORM_GATE_PROOF,
    odd_leech_first_packet_audit,
    odd_leech_h45_prefix_audit,
)


ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = ROOT / "data/enumerations/silver_odd_leech_p23_first_packet.json"
REPORT_PATH = ROOT / "results/silver_odd_leech_p23_first_packet.md"


def rational_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def silver_json(value: Silver) -> dict[str, str]:
    return {"a": rational_text(value.a), "b": rational_text(value.b)}


def pari_qfminim_cross_check() -> dict[str, object]:
    """Independently enumerate the Gram form's norm-three vectors in PARI."""

    executable = shutil.which("gp")
    if executable is None:
        raise RuntimeError("PARI/GP is required for the independent qfminim cross-check")
    matrix = "[" + ";".join(
        ",".join(str(value) for value in row) for row in ODD_LEECH_GRAM
    ) + "]"
    program = (
        f"G={matrix};"
        "R=qfminim(G,3,,2);V=R[3];H=vector(5);ok=1;"
        "for(j=1,matsize(V)[2],v=V[,j];"
        "if(v~*G*v!=3,ok=0);p=(G*v)[2];"
        "if(abs(p)>2,ok=0,H[p+3]++;H[-p+3]++));"
        "print(R[1]);print(matsize(V)[1]);print(matsize(V)[2]);"
        "print(ok);print(H);quit()\n"
    )
    process = subprocess.run(
        [executable, "-q"],
        input=program,
        text=True,
        capture_output=True,
        timeout=180,
        check=True,
    )
    lines = [line.strip() for line in process.stdout.splitlines() if line.strip()]
    if len(lines) != 5:
        raise AssertionError(f"unexpected PARI output: {process.stdout!r}")
    histogram = [int(value.strip()) for value in lines[4].strip("[]").split(",")]
    result = {
        "algorithm": "PARI qfminim(G,3,,flag=2)",
        "full_norm_at_most_3_count": int(lines[0]),
        "stored_half_shell_matrix_shape": [int(lines[1]), int(lines[2])],
        "all_stored_vectors_have_exact_norm_3_and_pairing_bound": lines[3] == "1",
        "norm4_witness_full_shell_pairing_histogram": [
            [pairing, histogram[pairing + 2]] for pairing in range(-2, 3)
        ],
    }
    if result != {
        "algorithm": "PARI qfminim(G,3,,flag=2)",
        "full_norm_at_most_3_count": 4096,
        "stored_half_shell_matrix_shape": [24, 2048],
        "all_stored_vectors_have_exact_norm_3_and_pairing_bound": True,
        "norm4_witness_full_shell_pairing_histogram": [
            [-2, 0],
            [-1, 1024],
            [0, 2048],
            [1, 1024],
            [2, 0],
        ],
    }:
        raise AssertionError(f"the independent PARI audit drifted: {result!r}")
    return result


def build_payload() -> dict[str, object]:
    audit = odd_leech_first_packet_audit()
    h45 = odd_leech_h45_prefix_audit()
    model = audit.model
    pari = pari_qfminim_cross_check()
    return {
        "schema_version": "0.2",
        "seed_id": "silver-H24-odd-Leech-a--3+4sqrt2",
        "scope": (
            "Complete primitive local mirror-norm candidate classes and exact "
            "Vinberg necessary-wall prefix through H=45/2+16sqrt(2). "
            "The strict hyperideal witness proves only that this prefix has "
            "infinite volume; later mirrors may cut it."
        ),
        "field": "Q(sqrt(2))",
        "ring": "Z[sqrt(2)]",
        "hyperbolic_dimension": 24,
        "timelike_coefficient": silver_json(ODD_LEECH_COEFFICIENT),
        "coefficient_absolute_norm": 23,
        "odd_leech_model": {
            "source": "https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/oddLeech.html",
            "coordinate_scale": "sqrt(8)",
            "construction": (
                "(2*C+4*Z_+^24) union (1+2*C+4*Z_+^24), divided by sqrt(8)"
            ),
            "gram_determinant": model.determinant,
            "leading_principal_minors": list(model.leading_principal_minors),
            "golay_weight_distribution": [list(row) for row in model.golay_weight_distribution],
            "minimum": 3,
            "minimal_shell_count": model.minimal_shell_count,
            "minimal_shell_sha256": model.minimal_shell_digest,
            "minimal_shell_proof": (
                "The unshifted raw construction has covolume 2^37: "
                "4*Z_+^24 has covolume 2^49 and the 12-dimensional code "
                "reduces it by 2^12. Adding the shifted coset halves this to "
                "2^36, and division by sqrt(8) in rank 24 gives covolume one. "
                "Every catalogue basis row belongs to the construction and "
                "its Gram determinant is one, so equal covolumes prove that "
                "the basis and two-coset models coincide. "
                "In the shifted coset every raw coordinate is odd, hence raw "
                "square at least 24 and lattice norm at least 3. Equality "
                "forces every coordinate to be +/-1, and the mod-4 residues "
                "give exactly one sign vector for each of the 4096 Golay "
                "codewords. In the unshifted coset, a nonzero codeword has "
                "weight at least 8; for codeword zero, the even-sum condition "
                "forces z^2 at least 2. Both cases have norm at least 4."
            ),
            "norm4_witness_raw_scaled": list(ODD_LEECH_NORM4_WITNESS_SCALED),
            "norm4_witness_square": model.witness_square,
            "norm4_witness_pairing_histogram": [
                list(row) for row in model.witness_pairing_histogram
            ],
            "independent_pari_cross_check": pari,
        },
        "mirror_norm_gate_proof": ODD_LEECH_NORM_GATE_PROOF,
        "canonical_mirror_norms": [silver_json(value) for value in audit.mirror_norms],
        "p23_spatial_gate": (
            "For q in {5+sqrt(2),12+7sqrt(2),10+2sqrt(2)}, y lies in "
            "a*O24 and the shell target must be divided by a^2."
        ),
        "height_zero": {
            "wall_count": audit.height_zero_wall_count,
            "remaining_canonical_norms": [
                silver_json(value) for value in HEIGHT_ZERO_CANONICAL_NORMS
            ],
            "their_traces": [int(value.trace()) for value in HEIGHT_ZERO_CANONICAL_NORMS],
            "nonzero_spatial_trace_minimum": HEIGHT_ZERO_TRACE_MINIMUM,
        },
        "control_census": [
            {
                "height": silver_json(level.height),
                "norm": silver_json(level.norm),
                "control": silver_json(level.control),
                "spatial_square_target": silver_json(level.spatial_square_target),
                "p23_spatial_gate": level.p23_spatial_gate,
                "reduced_target": silver_json(level.reduced_target),
                "reduced_target_integral": level.reduced_target.is_integral,
                "necessary_gram_triples": [
                    {
                        "u_square": triple.u_square,
                        "v_square": triple.v_square,
                        "u_dot_v": triple.u_dot_v,
                    }
                    for triple in level.necessary_gram_triples
                ],
            }
            for level in audit.controls
        ],
        "first_packet": {
            "height": silver_json(audit.first_height),
            "norm": silver_json(audit.first_norm),
            "control": silver_json(audit.first_control),
            "target": silver_json(audit.first_target),
            "unique_gram_triple": {
                "u_square": audit.first_gram_triple.u_square,
                "v_square": audit.first_gram_triple.v_square,
                "u_dot_v": audit.first_gram_triple.u_dot_v,
            },
            "root_description": "r_b=(1+sqrt(2))*(1,b), b in O24, b^2=3",
            "root_count": audit.first_packet_root_count,
            "distinct_spatial_pairings_between_distinct_b": list(
                audit.distinct_spatial_pairings
            ),
            "distinct_root_pairings": [
                silver_json(value) for value in audit.distinct_root_pairings
            ],
            "all_distinct_root_pairings_negative": True,
        },
        "strict_hyperideal_witness": {
            "affine_point": "(1,d), d=(4,4,0^22)/sqrt(8)",
            "d_square": audit.witness_square,
            "maximum_b_dot_d": audit.witness_maximum_pairing,
            "minimum_wall_margin_a_minus_b_dot_d": silver_json(
                audit.witness_minimum_wall_margin
            ),
            "hyperideal_excess_d_square_minus_a": silver_json(
                audit.witness_hyperideal_excess
            ),
            "conclusion": (
                "The first-packet cone contains a strict spacelike point and "
                "therefore an open hyperideal patch; this prefix has infinite volume."
            ),
        },
        "h45_necessary_wall_prefix": {
            "endpoint": silver_json(h45.terminal_height),
            "complete_control_count": len(h45.controls),
            "nonempty_necessary_level_count": len(h45.nonempty_necessary_levels),
            "terminal_norm": silver_json(h45.terminal_norm),
            "terminal_control": silver_json(h45.terminal_control),
            "terminal_target": silver_json(h45.terminal_target),
            "terminal_necessary_gram_triples": [
                {
                    "u_square": triple.u_square,
                    "v_square": triple.v_square,
                    "u_dot_v": triple.u_dot_v,
                }
                for triple in h45.terminal_triples
            ],
            "fixed_direction_bounds": [
                {
                    "triple": {
                        "u_square": bound.original_triple.u_square,
                        "v_square": bound.original_triple.v_square,
                        "u_dot_v": bound.original_triple.u_dot_v,
                    },
                    "change_of_variables_and_psd_bound": bound.transformed_description,
                    "maximum_pairing": silver_json(
                        bound.fixed_direction_maximum_pairing
                    ),
                    "wall_bound": silver_json(bound.wall_bound),
                    "strict_margin": silver_json(bound.strict_margin),
                }
                for bound in h45.fixed_direction_bounds
            ],
            "unscaled_direction_withdrawal": (
                "The unscaled d is not claimed to survive H45: the joint PSD "
                "superset maxima for d are 10+6sqrt(2), 9+6sqrt(2), and "
                "9+6sqrt(2), all exceeding the wall bound. The "
                "certified replacement is z=(5/6)d."
            ),
            "replacement_witness_scale": rational_text(h45.witness_scale),
            "global_minimum_wall_margin": silver_json(h45.minimum_wall_margin),
            "hyperideal_excess": silver_json(h45.hyperideal_excess),
            "proof_scope": (
                "The 3x3 PSD pairing optimization is over a necessary superset "
                "of all realized odd-Leech pairs, so it certifies every possible "
                "wall through the endpoint without expanding the enormous shells."
            ),
        },
        "status": "H45_NECESSARY_PREFIX_INFINITE_VOLUME; EVENTUAL_CHAMBER_UNRESOLVED",
    }


def silver_text(value: dict[str, str]) -> str:
    a = value["a"]
    b = value["b"]
    if b == "0":
        return a
    return f"{a} + ({b}) sqrt(2)"


def build_report(payload: dict[str, object]) -> str:
    controls = payload["control_census"]
    assert isinstance(controls, list)
    table_rows = []
    for index, row in enumerate(controls, start=1):
        assert isinstance(row, dict)
        triples = row["necessary_gram_triples"]
        table_rows.append(
            "| {index} | `{height}` | `{norm}` | `{control}` | {gate} | "
            "`{reduced}` | {integral} | {triples} |".format(
                index=index,
                height=silver_text(row["height"]),
                norm=silver_text(row["norm"]),
                control=silver_text(row["control"]),
                gate="yes" if row["p23_spatial_gate"] else "no",
                reduced=silver_text(row["reduced_target"]),
                integral="yes" if row["reduced_target_integral"] else "no",
                triples=len(triples),
            )
        )
    table = "\n".join(table_rows)
    model = payload["odd_leech_model"]
    first = payload["first_packet"]
    witness = payload["strict_hyperideal_witness"]
    h45 = payload["h45_necessary_wall_prefix"]
    assert isinstance(model, dict) and isinstance(first, dict) and isinstance(witness, dict)
    assert isinstance(h45, dict)
    return f"""# Silver odd-Leech H24: exact H45 prefix for norm-23 coefficient

## Outcome

For

`L = <-(-3+4 sqrt(2))> + O_24`

over `Z[sqrt(2)]`, the first nonempty positive Vinberg packet occurs at

`H = 7/2 + 3 sqrt(2)`.

It consists of exactly **4096 primitive norm-two roots**
`r_b=(1+sqrt(2))(1,b)`, one for each norm-three vector `b` of the odd Leech
lattice.  All distinct roots in the packet have negative pairing, so the
whole packet is accepted.  The exact affine point

`(1,d)`, with `d=(4,4,0^22)/sqrt(8)` and `d^2=4`,

does **not** itself carry the H45 claim.  The exact replacement
`(1,z)` with `z=(5/6)d` satisfies every possible wall through
`H=45/2+16sqrt(2)` strictly and is
hyperideal.  Therefore this complete necessary-wall prefix has **infinite
hyperbolic volume**.  This is not a
global nonreflectivity result: later mirrors may cut the patch.

## Exact odd-Leech model

Use the [Nebe--Sloane catalogue basis]({model['source']}) on the `sqrt(8)`
raw-coordinate scale.  With the catalogue extended Golay code `C` and
`Z_+^24={{z: sum(z) even}}`, the lattice is

`((2C+4Z_+^24) union (1+2C+4Z_+^24))/sqrt(8)`.

The encoded 24-row basis has integral Gram determinant **{model['gram_determinant']}**;
all 24 leading principal minors are positive.  Its Golay generators have
dimension 12 and weight enumerator
`1 + 759 x^8 + 2576 x^12 + 759 x^16 + x^24`.

The construction and basis are exactly equal, rather than merely similar
models.  Raw `4Z_+^24` has covolume `2^49`; adjoining the 12-dimensional code
reduces this to `2^37`, and the shifted coset halves it to `2^36`.  Dividing
24 coordinates by `sqrt(8)` gives covolume one.  Every catalogue basis row
lies in that construction, and its Gram determinant is also one, so the
containment and equal covolumes force equality.

The count 4096 is proved, not guessed from a sign-pattern shortcut.  In the
shifted coset all raw coordinates are odd, so norm three forces all 24 to be
`+/-1`.  Their residues modulo four are exactly one word of `C`, giving
precisely `|C|=4096` vectors.  The unshifted coset has minimum at least four:
a nonzero codeword has weight at least eight, while codeword zero uses a
nonzero even-sum integer vector.  Independently, PARI
`qfminim(G,3,,flag=2)` returned **4096** vectors from the determinant-one Gram
matrix (a `24 x 2048` stored half-shell).

An arbitrary collection of 4096 signs, a wrongly ordered Golay code, or an LP
that has not first proved compatibility with this basis is not a valid odd-
Leech shell model.  The certificate here avoids that failure mode entirely.

## Complete local mirror-norm candidates

{payload['mirror_norm_gate_proof']}

At height zero the p23-supported cases contradict primitivity.  The remaining
canonical norms have traces `2,4,4`, whereas any nonzero
`y=u+sqrt(2)v` has `Tr(y^2)=2(u^2+2v^2)>=6`.  Hence there are no height-zero
walls.

## Exact control census through the first packet

The search uses the complete Minkowski box for `H<=8` and then applies the
exact endpoint `H<=7/2+3sqrt(2)`, crystallographic integrality, conjugate
positivity, and the p23 spatial gate.

| # | height | q | x | p23 gate | target after gate | integral? | necessary triples |
|---:|---|---|---|:---:|---|:---:|---:|
{table}

The eight earlier controls are all p23-gated, and division by `a^2` makes
their reduced spatial targets nonintegral.  The endpoint has

`q=2`, `x=1+sqrt(2)`, `y^2=9+6sqrt(2)`.

Writing `y=u+sqrt(2)v` gives the unique necessary Gram triple
`(u^2,v^2,u.v)=(3,3,3)`.  Positive definiteness then gives `u=v=b`, with
`b^2=3`; conversely every one of the 4096 minimal vectors realizes a root.
Distinct Golay sign vectors have pairings `1,0,-1,-3`, producing root
pairings `-4-4sqrt(2), -7-6sqrt(2), -10-8sqrt(2), -16-12sqrt(2)`, all
strictly negative.

## Strict hyperideal witness

For every minimal `b`, the explicit norm-four vector `d` has
`b.d in {{-1,0,1}}` (counts `1024,2048,1024`).  The first walls in the affine
section are `b.z<=a`.  At `z=d`, the minimum wall margin is

`a-1 = -4+4sqrt(2) > 0`,

while the Lorentzian excess is

`d^2-a = 7-4sqrt(2) > 0`.

Thus a neighborhood of the segment's sphere-crossing point remains inside
the first-prefix chamber, proving an open hyperideal patch and infinite
volume for this certified prefix.

## Continuation through `H=45/2+16sqrt(2)`

The complete normalized arithmetic census through this endpoint contains
**{h45['complete_control_count']} controls**.  Apart from the first packet,
every possible level below the endpoint is empty after the p23 gate,
integrality, and all universal consequences of odd-Leech minimum three.  In
particular, the apparent intermediate target `(9,4,6)` is rank one and would
force a forbidden norm-one primitive direction.

At the endpoint,

`q=2`, `x=3+2sqrt(2)`, `y^2=47+32sqrt(2)`,

and exactly three necessary binary Gram types remain:

- `(27,10,16)`, equivalently `(b^2,w^2,b.w)=(3,10,-4)` for `b=u-2v`
  with `y=b+(2+sqrt(2))w`;
- `(25,11,16)`, equivalently `(d'^2,v^2,d'.v)=(4,11,5)` with
  `y=d'+(1+sqrt(2))v`;
- `(23,12,16)`, equivalently `(b^2,w^2,b.w)=(3,12,4)` with
  `y=b+(1+sqrt(2))w`.

There is no need to expand their enormous ordinary shells.  Extend each
binary Gram matrix by the fixed `d^2=4` row and column, require its exact
`3 x 3` Gram matrix to be positive semidefinite, and enumerate the finitely
many **joint** integral pairings.  Optimizing the full pair
`(first-vector.d, second-vector.d)`—not combining incompatible coordinatewise
bounds—gives for the unscaled `d`, respectively,

`y.d <= 10+6sqrt(2), 9+6sqrt(2), 9+6sqrt(2)`.

Consequently the unscaled `d` is explicitly **withdrawn as an H45 witness**:
all three necessary-superset maxima exceed the wall bound.  Scale along
the same ray and use the replacement `z=(5/6)d`.  Its maxima are

`y.z <= 25/3+5sqrt(2), 15/2+5sqrt(2), 15/2+5sqrt(2)`.

The wall bound is `a*x=7+6sqrt(2)`, leaving exact margins

`-4/3+sqrt(2), -1/2+sqrt(2), -1/2+sqrt(2)`.

All are strictly positive.  The tightest clearance in the entire prefix is
therefore `-4/3+sqrt(2)>0`.  The replacement remains hyperideal because

`z^2-a = 52/9-4sqrt(2)>0`.

Because the Gram scan deliberately optimizes over a
superset of realized odd-Leech pairs, the conclusion applies a fortiori to
every actual Vinberg wall through H45.

## Scope

Status: **{payload['status']}**.  The mirror-norm gates, height-zero exclusion,
control census, first shell, and the necessary-superset hyperideal witness are
exact and complete only through the stated H45 endpoint.  No claim is made
that later Vinberg walls preserve the patch or that the eventual reflection
chamber has infinite volume.
"""


def main() -> None:
    payload = build_payload()
    JSON_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    REPORT_PATH.write_text(build_report(payload))
    print(JSON_PATH.relative_to(ROOT))
    print(REPORT_PATH.relative_to(ROOT))


if __name__ == "__main__":
    main()
