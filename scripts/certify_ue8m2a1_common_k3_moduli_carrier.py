#!/usr/bin/env python3
"""Certify the common-moduli K3 carrier criterion for H10 chambers.

The lattice

    F = U + E8(-2) + <-2>

is the very-general Neron--Severi lattice of an explicit nine-parameter K3
family.  Consequently every oriented integral Coxeter realization by
primitive negative vectors in F is a weak equation-level full-Picard divisor
carrier on that family.  This script also derives the exhaustive primitive
negative integral-reflection alphabet of F from the exponent-two
discriminant group and replays the 13-, 21-, and 37-wall calibrations already
transported to F.

This is a carrier theorem.  Finite volume is inherited only from each
calibration's separate complete Vinberg certificate, and full-nef equality is
subject to the K3 root-completeness gate.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form

from scripts.certify_a15_11_ue8m2a1_basechange_family import F_GRAM
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data/enumerations/ue8m2a1_common_k3_moduli_carrier.json"
)
SOURCES = {
    "family": ROOT
    / "data/enumerations/a15_11_ue8m2a1_basechange_family.json",
    "h10_13": ROOT
    / "data/enumerations/"
    "k3_s3_c2_h10_13_ue8m2a1_full_picard_carrier.json",
    "h10_21": ROOT
    / "data/enumerations/a15_11_h10_21_x2_k3_carrier.json",
    "h10_37": ROOT
    / "data/enumerations/x2_h10_37_ue8m2a1_full_picard_carrier.json",
}
H37_ORIGIN = (
    ROOT
    / "data/enumerations/"
    "shimada_aut3singk3_x2_rational_involution_h10_37.json"
)
H37_ORIGIN_SHA256 = (
    "12914344403b7d3d6b69a2866c4150f0a0420a64914e4bc8584fa372799599a7"
)
EXPECTED_FILE_SHA256 = {
    "family": "1a5902092e3150bd380af19916ec5a6733392404a6cf2ecdba523f82b488532c",
    "h10_13": "93cec35a3c395482f4f99be827b3f5a6ba44764fc7defa22bf23eb6ee951acf0",
    "h10_21": "5f241c3659a79cf674d09558f72a5bc075d50b9fc6d87444ce8d253f1a4ab90b",
    "h10_37": "da8179b121df58e3c735c60ce31258a772c84400f7423e16c66f622d3c3d278f",
}
EXPECTED_PAYLOAD_SHA256 = {
    "family": "83d4d64c5a4003c296efc99c55557bd7ec4170fe26d578b9e950c02edff45f6d",
    "h10_13": "c033a684aa04b2c6d0bfe80d1b5393f0c52f5a9c9c1039140ca674ea9cc67a34",
    "h10_21": "a2905e05f15d5538d7f5896314bac4f52228f79a0ba2609e2ea13e99ff2e963c",
    "h10_37": "50c48a6ef7dcaaa4719b21d91f303a6d1983d5ab15cd998e4832d498a7ae4e08",
}
SCHEMA = "coxeter.k3.ue8m2a1_common_moduli_carrier.v1"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _payload_sha256(payload: dict[str, Any]) -> str:
    unsealed = dict(payload)
    unsealed.pop("payload_sha256", None)
    return hashlib.sha256(
        json.dumps(unsealed, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _smith(matrix: sp.Matrix) -> list[int]:
    diagonal = smith_normal_form(matrix, domain=ZZ)
    return [
        abs(int(diagonal[index, index]))
        for index in range(min(diagonal.shape))
        if diagonal[index, index]
    ]


def _matrix(payload: list[list[int | str]]) -> sp.Matrix:
    return sp.Matrix(
        [[sp.Rational(value) for value in row] for row in payload]
    )


def _reflection_alphabet(discriminant_exponent: int) -> list[dict[str, int]]:
    """Arithmetic possibilities for primitive negative mirror vectors.

    If v is primitive and d=div(v), the class v/d has exact order d in the
    discriminant group.  Hence d divides the discriminant exponent.  The
    reflection in v is integral exactly when |v^2| divides 2d.
    """

    answer: list[dict[str, int]] = []
    for divisibility in range(1, discriminant_exponent + 1):
        if discriminant_exponent % divisibility:
            continue
        for absolute_square in range(2, 2 * divisibility + 1, 2):
            if 2 * divisibility % absolute_square == 0:
                answer.append(
                    {
                        "square": -absolute_square,
                        "divisibility": divisibility,
                    }
                )
    return answer


def _primitive_reflection_record(vector: sp.Matrix) -> dict[str, Any]:
    """Audit one primitive vector and its lattice reflection on F."""

    content = math.gcd(*(abs(int(value)) for value in vector))
    square = int((vector.T * F_GRAM * vector)[0])
    divisibility = math.gcd(
        *(abs(int(value)) for value in F_GRAM * vector)
    )
    reflection = sp.simplify(
        sp.eye(F_GRAM.rows)
        - sp.Rational(2, square) * vector * (vector.T * F_GRAM)
    )
    dual_class = vector / divisibility
    exact_order = next(
        order
        for order in range(1, divisibility + 1)
        if all(value.is_Integer for value in order * dual_class)
    )
    return {
        "coordinates": [int(value) for value in vector],
        "square": square,
        "divisibility": divisibility,
        "primitive": content == 1,
        "dual_class_exact_order": exact_order,
        "reflection_integral": all(value.is_Integer for value in reflection),
        "reflection_preserves_F": reflection.T * F_GRAM * reflection
        == F_GRAM,
        "hyperbolic_mirror_normal": square < 0,
    }


def _wall_span(records: list[dict[str, Any]]) -> dict[str, Any]:
    columns = [
        sp.Matrix([sp.Rational(value) for value in record["primitive_F_coordinates"]])
        for record in records
    ]
    matrix = sp.Matrix.hstack(*columns)
    invariants = _smith(matrix)
    rank = len(invariants)
    return {
        "rank": rank,
        "smith": invariants,
        "index_in_F": math.prod(invariants) if rank == F_GRAM.rows else None,
    }


def _load_sources() -> dict[str, dict[str, Any]]:
    sources = {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCES.items()
    }
    for key, path in SOURCES.items():
        if _sha256(path) != EXPECTED_FILE_SHA256[key]:
            raise ArithmeticError(f"pinned {key} file changed")
        if sources[key].get("all_checks_pass") is not True:
            raise ArithmeticError(f"pinned {key} certificate is not green")
    observed_payloads = {
        "family": sources["family"]["payload_sha256"],
        "h10_13": sources["h10_13"]["payload_sha256"],
        "h10_21": sources["h10_21"]["certificate_payload_sha256"],
        "h10_37": sources["h10_37"]["payload_sha256"],
    }
    if observed_payloads != EXPECTED_PAYLOAD_SHA256:
        raise ArithmeticError(
            f"pinned payload map changed: {observed_payloads}"
        )
    return sources


def build_certificate() -> dict[str, Any]:
    sources = _load_sources()
    if _sha256(H37_ORIGIN) != H37_ORIGIN_SHA256:
        raise ArithmeticError("pinned H10/37 lineage source changed")
    h37_origin = json.loads(H37_ORIGIN.read_text(encoding="utf-8"))
    if h37_origin.get("all_checks_pass") is not True:
        raise ArithmeticError("pinned H10/37 lineage source is not green")
    family = sources["family"]
    h13 = sources["h10_13"]
    h21 = sources["h10_21"]
    h37 = sources["h10_37"]

    inverse = F_GRAM.inv()
    smith = _smith(F_GRAM)
    exponent_two = all(value.q == 1 for value in 2 * inverse)
    exponent_not_one = any(value.q != 1 for value in inverse)
    alphabet = _reflection_alphabet(2)
    expected_alphabet = [
        {"square": -2, "divisibility": 1},
        {"square": -2, "divisibility": 2},
        {"square": -4, "divisibility": 2},
    ]
    explicit_vectors = {
        "minus2_div1": sp.Matrix([1, -1, *([0] * 9)]),
        "minus2_div2": sp.Matrix([*([0] * 10), 1]),
        "minus4_div2": sp.Matrix([0, 0, 1, *([0] * 8)]),
    }
    alphabet_witnesses = {
        key: _primitive_reflection_record(vector)
        for key, vector in explicit_vectors.items()
    }
    # e+f in the U summand is a primitive positive-norm integral lattice
    # reflection.  It is deliberately excluded from the hyperbolic mirror
    # alphabet, whose normals must have negative square.
    positive_norm_boundary = _primitive_reflection_record(
        sp.Matrix([1, 1, *([0] * 9)])
    )

    h13_chamber = h13["transported_H10_13_chamber"]
    h37_chamber = h37["transported_H10_37_chamber"]
    h21_walls = h21["walls"]
    h21_span = _wall_span(h21_walls["records"])
    full_lineage = h37_origin["arithmetic_and_literal_lineage"]
    full_chamber = full_lineage["published_full_reflection_chamber"]
    solomon = full_lineage["solomon"]
    chi_full = Fraction(solomon["full_Wr_chamber"]["rational_euler_characteristic"])
    chi_h21 = Fraction(h21["arithmetic_lineage"]["coxeter_euler_characteristic"])
    chi_h37 = Fraction(solomon["candidate"]["rational_euler_characteristic"])

    calibrations = {
        "K3-S3-C2-H10-13": {
            "facets": h13_chamber["facet_count"],
            "wall_types": h13_chamber["primitive_mirror_type_histogram"],
            "wall_span_index": h13_chamber["wall_span_index"],
            "finite_vertices": h13["finite_volume"]["finite_vertices"],
            "cusps": h13["finite_volume"]["cusps"],
            "elliptic_rank9_edges": h13["finite_volume"]["elliptic_rank9_edges"],
            "classification": h13["finite_volume"]["classification"],
            "index_in_full_reflection_group": 1,
            "rational_euler_characteristic": str(chi_full),
            "all_checks_pass": all(h13_chamber["checks"].values())
            and h13["finite_volume"]["literal_bad_edges"] == 0,
        },
        "A15_11-C2-H10-21": {
            "facets": h21_walls["count"],
            "wall_types": h21_walls["primitive_F_type_histogram"],
            "wall_span_index": h21_span["index_in_F"],
            "finite_vertices": h21["finite_volume"]["literal_vinberg"]["finite_vertices"],
            "cusps": h21["finite_volume"]["literal_vinberg"]["ideal_vertices"],
            "elliptic_rank9_edges": h21["finite_volume"]["literal_vinberg"]["elliptic_edges"],
            "classification": h21["finite_volume"]["literal_vinberg"]["classification"],
            "index_in_full_reflection_group": h21["arithmetic_lineage"]["index_in_full_reflection_group"],
            "rational_euler_characteristic": str(chi_h21),
            "all_checks_pass": h21_walls["all_21_integral_on_F"]
            and h21_span["rank"] == 11
            and all(
                all(record["checks"].values())
                for record in h21_walls["records"]
            ),
        },
        "X2-RAT-C2-H10-37": {
            "facets": h37_chamber["facet_count"],
            "wall_types": h37_chamber["primitive_mirror_type_histogram"],
            "wall_span_index": h37_chamber["wall_span_index"],
            "finite_vertices": h37["finite_volume"]["finite_vertices"],
            "cusps": h37["finite_volume"]["cusps"],
            "elliptic_rank9_edges": h37["finite_volume"]["elliptic_rank9_edges"],
            "classification": h37["finite_volume"]["classification"],
            "index_in_full_reflection_group": solomon["index_candidate_in_full_Wr"],
            "rational_euler_characteristic": str(chi_h37),
            "all_checks_pass": all(h37_chamber["checks"].values())
            and h37["finite_volume"]["bad_edges"] == 0,
        },
    }

    family_gram = _matrix(
        family["fixed_locus_and_invariant_lattice"]["gram"]
    )
    source_grams = {
        "h10_13": _matrix(
            h13["full_Picard_lattice"]["standard_geometric_gram"]
        ),
        "h10_21": _matrix(h21["lattice_bridge"]["standard_model_gram"]),
        "h10_37": _matrix(
            h37["full_Picard_lattice"]["standard_geometric_gram"]
        ),
    }

    checks = {
        "F_even_rank11_signature_1_10": F_GRAM.rows == 11
        and all(F_GRAM[index, index] % 2 == 0 for index in range(11))
        and exact_inertia(F_GRAM) == (1, 10, 0),
        "F_determinant_512": F_GRAM.det() == 512,
        "F_discriminant_is_2_elementary_length9": smith
        == [1, 1, *([2] * 9)],
        "F_discriminant_exponent_exactly2": exponent_two
        and exponent_not_one,
        "reflection_alphabet_exact": alphabet == expected_alphabet,
        "all_three_negative_types_have_explicit_F_witnesses": (
            {
                (row["square"], row["divisibility"])
                for row in alphabet_witnesses.values()
            }
            == {(row["square"], row["divisibility"]) for row in alphabet}
            and all(
                row["primitive"]
                and row["reflection_integral"]
                and row["reflection_preserves_F"]
                and row["hyperbolic_mirror_normal"]
                for row in alphabet_witnesses.values()
            )
        ),
        "primitive_divisibility_equals_dual_class_order_on_witnesses": all(
            row["dual_class_exact_order"] == row["divisibility"]
            for row in alphabet_witnesses.values()
        ),
        "positive_norm_integral_reflection_is_explicitly_out_of_scope": (
            positive_norm_boundary["square"] == 2
            and positive_norm_boundary["divisibility"] == 1
            and positive_norm_boundary["primitive"]
            and positive_norm_boundary["reflection_integral"]
            and positive_norm_boundary["reflection_preserves_F"]
            and not positive_norm_boundary["hyperbolic_mirror_normal"]
        ),
        "explicit_family_same_F": family_gram == F_GRAM,
        "all_three_sources_same_F": all(
            gram == F_GRAM for gram in source_grams.values()
        ),
        "family_is_nine_dimensional_and_very_general_full_Picard": (
            family["family"]["parameter_count"] == 9
            and family["generic_full_picard"]["parameter_moduli_dimension"]
            == 9
            and family["generic_full_picard"]["claim"]
            == "NS(X)=U+E8(-2)+<-2> for a very general parameter"
        ),
        "family_open_nonempty": bool(
            family["family"]["nonempty_open_witness"]
        ),
        "three_calibrations_green": all(
            row["all_checks_pass"] for row in calibrations.values()
        ),
        "all_observed_types_lie_in_alphabet": all(
            set(row["wall_types"])
            <= {f"{item['square']},{item['divisibility']}" for item in alphabet}
            for row in calibrations.values()
        ),
        "facet_counts_13_21_37": [
            row["facets"] for row in calibrations.values()
        ]
        == [13, 21, 37],
        "full_reflection_chamber_is_13_wall_nonuniform": (
            h13_chamber["known_full_reflection_chamber"] is True
            and full_chamber["id"] == "K3-S3-C2-H10-13"
            and full_chamber["facet_count"] == 13
            and full_chamber["cusp_count"] == 4
            and all(full_lineage["checks"].values())
            and full_lineage["checks"][
                "candidate_reflections_lie_in_full_root_reflection_group"
            ]
        ),
        "subgroup_indices_1_6720_241920": [
            calibrations[key]["index_in_full_reflection_group"]
            for key in [
                "K3-S3-C2-H10-13",
                "A15_11-C2-H10-21",
                "X2-RAT-C2-H10-37",
            ]
        ]
        == [1, 6720, 241920],
        "subgroup_indices_equal_exact_Euler_ratios": (
            chi_h21 / chi_full == 6720
            and chi_h37 / chi_full == 241920
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(
            [key for key, value in checks.items() if not value]
        )

    certificate: dict[str, Any] = {
        "schema": SCHEMA,
        "theorem": {
            "common_equation_family_criterion": (
                "Every genuine finite-sided oriented H10 Coxeter polyhedron "
                "admitting a realization by primitive negative vectors whose "
                "reflections preserve F=U+E8(-2)+<-2> is a weak equation-level "
                "full-Picard divisor carrier, after choosing an F-marking of "
                "a very general member of the displayed nine-parameter K3 family."
            ),
            "converse_for_this_marked_family": (
                "Every primitive negative Picard vector defining an integral "
                "hyperbolic mirror reflection on a very general marked member "
                "has type (-2,div1), (-2,div2), or (-4,div2). Positive-norm "
                "integral lattice reflections are outside this mirror alphabet."
            ),
            "finite_volume_subgroup_criterion": (
                "Let W(F) be the hyperbolic group generated by all primitive "
                "negative integral mirror reflections of F. Its complete "
                "13-wall chamber has finite volume and four cusps. For a "
                "genuine Coxeter polyhedron P which is the fundamental domain "
                "of its facet-reflection group Gamma<=W(F), vol(P) is finite "
                "exactly when [W(F):Gamma] is finite; every such P is noncompact."
            ),
            "compact_no_go_inside_F": (
                "No compact H10 Coxeter polytope admits an oriented primitive "
                "negative integral-mirror realization in F."
            ),
            "proof_inputs": [
                "the explicit family has very-general NS exactly F",
                "the primitive K3 embedding of F supplies the marked period space",
                "a primitive vector has divisibility dividing exp(A_F)=2",
                "reflection integrality is equivalent to |v^2| dividing 2*div(v)",
                "a genuine Coxeter polyhedron is a fundamental domain for its facet-reflection group",
            ],
        },
        "lattice": {
            "name": "F=U+E8(-2)+<-2>",
            "rank": 11,
            "signature": [1, 10, 0],
            "determinant": 512,
            "smith": smith,
            "discriminant_exponent": 2,
            "primitive_integral_reflection_alphabet": alphabet,
            "explicit_negative_type_witnesses": alphabet_witnesses,
            "positive_norm_integral_reflection_boundary_witness": positive_norm_boundary,
            "alphabet_scope": (
                "These are the exhaustive possibilities for primitive negative "
                "vectors defining hyperbolic mirrors. They do not classify "
                "positive-norm integral lattice reflections."
            ),
        },
        "equation_family": {
            "equation": family["family"]["k3_weierstrass_model"],
            "A": family["family"]["A"],
            "B": family["family"]["B"],
            "parameter_count": family["family"]["parameter_count"],
            "parameter_open": family["family"]["parameter_open"],
            "nonempty_open_witness": family["family"]["nonempty_open_witness"],
            "very_general_full_Picard": family["generic_full_picard"]["claim"],
            "moduli_dimension": family["generic_full_picard"]["parameter_moduli_dimension"],
            "weak_carrier_scope": (
                "Equation-level refers to the ambient K3 surface. An F-marking "
                "makes every transported wall an algebraic divisor class in "
                "full NS; it does not provide equations, effective representatives, "
                "or irreducible curves for arbitrary transported classes."
            ),
        },
        "calibrations": calibrations,
        "full_reflection_group": {
            "name": "W(F)=W_r(11,9,1)",
            "fundamental_chamber": "K3-S3-C2-H10-13",
            "fundamental_facets": full_chamber["facet_count"],
            "finite_vertices": full_chamber["finite_vertex_count"],
            "cusps": full_chamber["cusp_count"],
            "nonuniform": full_chamber["cusp_count"] > 0,
            "rational_euler_characteristic": solomon["full_Wr_chamber"]["rational_euler_characteristic"],
            "proof_of_subgroup_criterion": (
                "For a genuine Coxeter polyhedron P, the Coxeter-Poincare "
                "theorem makes P a fundamental domain of its facet group Gamma. "
                "Primitive reduction puts every facet reflection in W(F), so "
                "Gamma<=W(F). A Gamma domain is a disjoint union modulo "
                "boundaries of [W(F):Gamma] W(F)-chambers; hence its volume is "
                "finite iff the index is finite. A finite-index subgroup of "
                "the nonuniform W(F) contains powers of parabolic elements and "
                "is nonuniform, so P cannot be compact."
            ),
        },
        "promotion_boundary": {
            "minus_two_effectivity": (
                "For each chamber separately, a Weyl re-marking moves a "
                "generic strict interior point into the ample cone; K3 "
                "Riemann--Roch then makes every oriented (-2) facet effective."
            ),
            "reducibility": (
                "An effective (-2) class may be reducible and therefore need "
                "not define an actual nef facet."
            ),
            "full_nef_criterion": (
                "After the chamber-specific Weyl re-marking, full-nef equality "
                "holds exactly when every irredundant facet is a primitive "
                "(-2)-root hyperplane and no other (-2)-root hyperplane meets "
                "the chamber interior."
            ),
            "minus_four_obstruction": (
                "An irredundant primitive (-4,div2) facet cannot be a full "
                "K3 nef wall, though it remains an integral divisor reflection."
            ),
            "volume_gate": (
                "The family supplies no volume conclusion. Use either a "
                "separate complete Vinberg endpoint/parabolic audit or an "
                "independent proof of finite index in W(F)."
            ),
        },
        "scope": {
            "new_abstract_polytopes": "none",
            "new_result": (
                "a reusable equation-family carrier criterion and exhaustive "
                "wall alphabet, calibrated by three known H10 chambers"
            ),
            "surface_automorphisms": "not asserted",
            "simultaneous_ample_marking": (
                "not asserted; separate chambers may require separate Weyl markings"
            ),
            "equation_level_wall_representatives": (
                "not asserted; arbitrary transported walls are divisor classes under a marking"
            ),
            "full_nef_for_calibrations": "not asserted",
        },
        "sources": {
            key: {
                "path": str(path.relative_to(ROOT)),
                "file_sha256": EXPECTED_FILE_SHA256[key],
                "payload_sha256": EXPECTED_PAYLOAD_SHA256[key],
            }
            for key, path in SOURCES.items()
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    certificate["sources"]["h10_37_lineage"] = {
        "path": str(H37_ORIGIN.relative_to(ROOT)),
        "file_sha256": H37_ORIGIN_SHA256,
    }
    certificate["payload_sha256"] = _payload_sha256(certificate)
    return certificate


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()

    certificate = build_certificate()
    if args.write:
        OUTPUT.write_text(
            json.dumps(certificate, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    if args.check_archive:
        archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
        if archived != certificate:
            raise ArithmeticError("archived common-carrier certificate changed")
    print(json.dumps(certificate, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
