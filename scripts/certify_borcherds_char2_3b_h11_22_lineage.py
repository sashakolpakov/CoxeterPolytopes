#!/usr/bin/env python3
"""Certify the full-reflection-group lineage of ``SS2-3B-H11-22``.

The production certificate proves that the restricted 22-wall chamber is a
finite-volume noncompact Coxeter polytope.  This independent supplement pins
the matching Scharlau--Kirschmer class-number-one census entry and proves that
the 22 displayed mirrors are *all* the simple mirrors of the lattice.

The all-root step uses an integral control vector and PARI/GP ``qfminim``.
It exhausts the height-zero stabilizer and the first two admissible Vinberg
shells.  Those shells contain exactly the displayed 12 short and 10 long
roots.  Since their chamber is already finite-volume, Vinberg's stopping
criterion proves that the section group is the full reflection group.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter, deque
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import subprocess
from typing import Any
import zipfile

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form

from scripts.certify_borcherds_char2_3b_h11_22 import (
    _digest,
    _matrix_payload,
    _smith_decomposition,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
PRODUCTION_CERTIFICATE = (
    ROOT / "data/enumerations/borcherds_char2_3b_h11_22_certificate.json"
)
OUTPUT = (
    ROOT / "data/enumerations/borcherds_char2_3b_h11_22_lineage.json"
)
DEFAULT_SOURCE_ARCHIVE = Path("/private/tmp/Scharlau-Kirschmer-Lorentz.zip")
SCHEMA = "coxeter.borcherds_char2_3b_h11_22_lineage.v1"

SOURCE_PAGE_URL = (
    "https://www.math.rwth-aachen.de/~Markus.Kirschmer/lorentz/index.html"
)
SOURCE_ARCHIVE_URL = (
    "https://www.math.rwth-aachen.de/~Markus.Kirschmer/lorentz/Lorentz.zip"
)
SOURCE_ARCHIVE_SHA256 = (
    "74a4dd3e534f1e793cba45765ff33575dbd331d5cdce6f34cdd64340b3f0d4a5"
)
LIST11_SHA256 = "2ca0725232cb101e60380f14ad0ab7777478c60eff06991ba591b6d3c71e3184"
GRAM11_SHA256 = "af8d52984fc9caabc00ffe670ae646d2151e68d141082544e18f4827ba142540"
PRODUCTION_FACET_GRAM_SHA256 = (
    "2944c4ce6ebab9a5fb1b52c9a0f7c2a105656863036c3e2badbe52b15c9763cb"
)

# One-based entry 27.  The official table consists of class-number-one
# strongly squarefree Lorentzian lattices and records, in order, determinant,
# genus symbol, elliptic-subdiagram counts, and number of cusps.
FROZEN_LIST_ENTRY = (
    243,
    "II_{11,1}(3^{-5})",
    [22, 231, 1480, 6320, 18704, 38986,
     57182, 57785, 38140, 14580, 2304],
    72,
)

# Gram11 entry 27, reshaped from the source's flattened 12-by-12 record.
FROZEN_SOURCE_GRAM = [
    [0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 2, -1, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, -1, 2, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 2, -1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, -1, 2, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 2, 1, 1, 1, 1, -1],
    [0, 0, 0, 0, 0, 0, 1, 2, 1, 0, 0, -1],
    [0, 0, 0, 0, 0, 0, 1, 1, 2, 1, 1, -1],
    [0, 0, 0, 0, 0, 0, 1, 0, 1, 2, 1, -1],
    [0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 2, -1],
    [0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, 2],
]


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256_file(path: Path) -> str:
    return _sha256_bytes(path.read_bytes())


def _mod_one(value: sp.Rational) -> sp.Rational:
    value = sp.Rational(value)
    return value - sp.floor(value)


def _mod_two(value: sp.Rational) -> sp.Rational:
    value = sp.Rational(value)
    return value - 2 * sp.floor(value / 2)


def _discriminant_q_histogram(gram: sp.Matrix) -> dict[str, int]:
    """Enumerate the finite quadratic form ``q(x)=x^2 mod 2Z`` exactly."""

    inverse = gram.inv()
    generators = [
        tuple(_mod_one(inverse[row, column]) for row in range(gram.rows))
        for column in range(gram.cols)
    ]
    zero = tuple(sp.Rational(0) for _ in range(gram.rows))
    classes = {zero}
    queue: deque[tuple[sp.Rational, ...]] = deque([zero])
    while queue:
        current = queue.popleft()
        for generator in generators:
            candidate = tuple(
                _mod_one(left + right)
                for left, right in zip(current, generator, strict=True)
            )
            if candidate not in classes:
                classes.add(candidate)
                queue.append(candidate)
    histogram = Counter(
        str(_mod_two((sp.Matrix(vector).T * gram * sp.Matrix(vector))[0]))
        for vector in classes
    )
    return {key: histogram[key] for key in sorted(histogram, key=sp.Rational)}


def _standard_model() -> sp.Matrix:
    hyperbolic_plane = sp.Matrix([[0, 1], [1, 0]])
    a2_negative = sp.Matrix([[-2, 1], [1, -2]])
    return sp.diag(hyperbolic_plane, *([a2_negative] * 5))


def verify_source_archive(path: Path) -> dict[str, Any]:
    raw = path.read_bytes()
    if _sha256_bytes(raw) != SOURCE_ARCHIVE_SHA256:
        raise AssertionError("Lorentz.zip SHA-256 mismatch")
    with zipfile.ZipFile(path) as archive:
        list_raw = archive.read("List11")
        gram_raw = archive.read("Gram11")
        if _sha256_bytes(list_raw) != LIST11_SHA256:
            raise AssertionError("List11 SHA-256 mismatch")
        if _sha256_bytes(gram_raw) != GRAM11_SHA256:
            raise AssertionError("Gram11 SHA-256 mismatch")
        entries = ast.literal_eval(
            list_raw.decode("ascii").replace("<", "(").replace(">", ")")
        )
        grams = ast.literal_eval(gram_raw.decode("ascii"))
    if tuple(entries[26]) != FROZEN_LIST_ENTRY:
        raise AssertionError("List11 entry 27 changed")
    source_gram = sp.Matrix(12, 12, grams[26])
    if source_gram != sp.Matrix(FROZEN_SOURCE_GRAM):
        raise AssertionError("Gram11 entry 27 changed")
    return {
        "archive_path": str(path),
        "archive_sha256": SOURCE_ARCHIVE_SHA256,
        "list11_sha256": LIST11_SHA256,
        "gram11_sha256": GRAM11_SHA256,
        "entry_one_based": 27,
        "entry_reparsed_exactly": True,
    }


def _gp_matrix(matrix: sp.Matrix) -> str:
    return "[" + ";".join(
        ",".join(str(int(matrix[row, column])) for column in range(matrix.cols))
        for row in range(matrix.rows)
    ) + "]"


def _run_gp(program: str, timeout: int = 60) -> list[str]:
    completed = subprocess.run(
        ["gp", "-fq", "-s", "500000000"],
        input=program,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=True,
    )
    if completed.stderr.strip():
        raise AssertionError(f"PARI/GP error: {completed.stderr.strip()}")
    return completed.stdout.splitlines()


def _parse_tagged(lines: list[str]) -> dict[str, int]:
    values: dict[str, int] = {}
    for line in lines:
        if line.startswith("ROOT|"):
            continue
        key, value = line.split("|", 1)
        values[key] = int(value)
    return values


def _parse_root(line: str) -> tuple[int, ...]:
    payload = line.split("|", 1)[1].strip()
    if not (payload.startswith("[") and payload.endswith("]")):
        raise AssertionError(f"malformed PARI root row: {line}")
    return tuple(int(value) for value in payload[1:-1].split(","))


def _height_zero_audit(
    gram: sp.Matrix, kernel: sp.Matrix
) -> dict[str, Any]:
    positive_gram = -(kernel.T * gram * kernel)
    program = (
        f"Q={_gp_matrix(positive_gram)};"
        f"K={_gp_matrix(kernel)};G={_gp_matrix(gram)};"
        "R=qfminim(Q,6,,2);V=R[3];"
        "c2=0;c4=0;c6=0;c6d3=0;"
        "for(j=1,matsize(V)[2],z=V[,j];n=z~*Q*z;"
        "if(n==2,c2++);if(n==4,c4++);"
        "if(n==6,c6++;x=K*z;y=G*x;ok=1;"
        "for(k=1,12,if(y[k]%3,ok=0));if(ok,c6d3++)));"
        "print(\"TOTAL|\",R[1]);print(\"N2|\",c2);"
        "print(\"N4|\",c4);print(\"N6|\",c6);"
        "print(\"N6DIV3|\",c6d3);quit()"
    )
    counts = _parse_tagged(_run_gp(program))
    expected = {"TOTAL": 1226, "N2": 0, "N4": 111, "N6": 502, "N6DIV3": 0}
    if counts != expected:
        raise AssertionError(f"height-zero qfminim census changed: {counts}")
    return {
        "positive_kernel_gram": _matrix_payload(positive_gram),
        "positive_kernel_determinant": int(positive_gram.det()),
        "qfminim_total_signed_vectors_norm_at_most_6": counts["TOTAL"],
        "representative_pairs_by_norm": {
            "2": counts["N2"], "4": counts["N4"], "6": counts["N6"]
        },
        "norm_6_divisibility_3_pairs": counts["N6DIV3"],
        "reflective_root_pair_count": counts["N2"] + counts["N6DIV3"],
    }


def _affine_shell_audit(
    gram: sp.Matrix,
    control: sp.Matrix,
    kernel: sp.Matrix,
    bezout: sp.Matrix,
    negative_norm: int,
    pairing: int,
) -> dict[str, Any]:
    # For y=(z,t), x=Kz+t*pairing*b and
    # E(y)=-x^2+t^2.  The Schur complement is 1-pairing^2/p^2>0
    # for pairing 1 and 3, so E is positive definite.  Every target root has
    # t=1 and E-value negative_norm+1.
    affine_basis = kernel.row_join(pairing * bezout)
    positive_form = -(affine_basis.T * gram * affine_basis)
    positive_form[-1, -1] += 1
    if not positive_form.is_positive_definite:
        raise AssertionError("affine-shell bounding form is not positive definite")
    expected_divisibility = 1 if negative_norm == 2 else 3
    program = (
        f"Q={_gp_matrix(positive_form)};B={_gp_matrix(affine_basis)};"
        f"G={_gp_matrix(gram)};p={_gp_matrix(control)};"
        f"R=qfminim(Q,{negative_norm + 1},,2);V=R[3];"
        "c=0;reflective=0;"
        "for(j=1,matsize(V)[2],z=V[,j];"
        "if(abs(z[12])==1,x=B*z;if(z[12]==-1,x=-x);"
        f"if(x~*G*x=={-negative_norm} && (x~*G*p)[1]=={pairing},"
        "c++;y=G*x;d=0;for(k=1,12,d=gcd(d,abs(y[k])));"
        f"if(d=={expected_divisibility},reflective++;"
        "print(\"ROOT|\",Vec(x))))));"
        "print(\"TOTAL|\",R[1]);print(\"CANDIDATES|\",c);"
        "print(\"REFLECTIVE|\",reflective);quit()"
    )
    lines = _run_gp(program)
    roots = {_parse_root(line) for line in lines if line.startswith("ROOT|")}
    counts = _parse_tagged([line for line in lines if not line.startswith("ROOT|")])
    expected_counts = (
        {"TOTAL": 24, "CANDIDATES": 12, "REFLECTIVE": 12}
        if (negative_norm, pairing) == (2, 1)
        else {"TOTAL": 6534, "CANDIDATES": 1720, "REFLECTIVE": 10}
    )
    if counts != expected_counts or len(roots) != expected_counts["REFLECTIVE"]:
        raise AssertionError(
            f"Vinberg shell {(negative_norm, pairing)} changed: "
            f"{counts}, {len(roots)} roots"
        )
    return {
        "negative_norm": negative_norm,
        "control_pairing": pairing,
        "vinberg_height_pairing_squared_over_negative_norm": str(
            Fraction(pairing**2, negative_norm)
        ),
        "positive_bounding_gram": _matrix_payload(positive_form),
        "qfminim_total_signed_vectors_through_bound": counts["TOTAL"],
        "oriented_target_shell_candidate_count": counts["CANDIDATES"],
        "reflective_target_root_count": counts["REFLECTIVE"],
        "oriented_reflective_root_coordinates_sha256": _digest(sorted(roots)),
        "roots": roots,
    }


def _vinberg_all_root_audit(
    gram: sp.Matrix, roots: list[sp.Matrix], finite_volume: bool
) -> dict[str, Any]:
    root_sum = sum(roots, sp.zeros(gram.rows, 1))
    if any(int(value) % 3 for value in root_sum):
        raise AssertionError("one third of the root sum is no longer integral")
    control = root_sum / 3
    pairing = control.T * gram
    control_divisibility = math.gcd(*(abs(int(value)) for value in gram * control))
    if not (
        math.gcd(*(abs(int(value)) for value in control)) == 1
        and int((control.T * gram * control)[0]) == 14
        and control_divisibility == 1
    ):
        raise AssertionError("the primitive Vinberg control vector changed")

    diagonal, _, right = _smith_decomposition(pairing)
    if diagonal[0, 0] != 1 or any(diagonal[0, index] for index in range(1, 12)):
        raise AssertionError("control pairing row is no longer primitive")
    bezout = right[:, 0]
    kernel = right[:, 1:]
    if not (pairing * bezout == sp.ones(1, 1) and pairing * kernel == sp.zeros(1, 11)):
        raise AssertionError("Smith control decomposition failed")

    height_zero = _height_zero_audit(gram, kernel)
    short_shell = _affine_shell_audit(gram, control, kernel, bezout, 2, 1)
    long_shell = _affine_shell_audit(gram, control, kernel, bezout, 6, 3)

    expected_short = {
        tuple(int(value) for value in root)
        for root in roots
        if int((root.T * gram * root)[0]) == -2
    }
    expected_long = {
        tuple(int(value) for value in root)
        for root in roots
        if int((root.T * gram * root)[0]) == -6
    }
    shell_match = (
        short_shell.pop("roots") == expected_short
        and long_shell.pop("roots") == expected_long
    )
    if not shell_match:
        raise AssertionError("qfminim shells do not equal the displayed roots")
    if not finite_volume:
        raise AssertionError("Vinberg stopping requires the production volume theorem")
    pairings = Counter(int((control.T * gram * root)[0]) for root in roots)
    if pairings != Counter({1: 12, 3: 10}):
        raise AssertionError("simple-root control pairings changed")
    version = subprocess.run(
        ["gp", "--version-short"], capture_output=True, text=True,
        timeout=5, check=True,
    ).stdout.strip()
    return {
        "method": "exact PARI/GP shell enumeration plus Vinberg finite-volume stopping",
        "pari_gp_version": version,
        "control_vector": [int(value) for value in control],
        "control_square": 14,
        "control_divisibility": control_divisibility,
        "control_pairing_row": [int(value) for value in pairing],
        "simple_root_pairing_histogram": {str(k): v for k, v in sorted(pairings.items())},
        "saturated_control_orthogonal_kernel_basis": _matrix_payload(kernel),
        "height_zero": height_zero,
        "positive_shells": [short_shell, long_shell],
        "shell_roots_equal_displayed_simple_roots": shell_match,
        "finite_volume_input": finite_volume,
        "vinberg_stopping_applies": True,
        "full_reflection_chamber": True,
        "reflection_subgroup_index": 1,
    }


def build_certificate(
    production_path: Path = PRODUCTION_CERTIFICATE,
    source_archive: Path | None = None,
) -> dict[str, Any]:
    production = json.loads(production_path.read_text(encoding="utf-8"))
    chamber = production["restricted_chamber"]
    if not (
        production["id"] == "SS2-3B-H11-22"
        and production["all_checks_pass"] is True
        and production["vinberg"]["finite_volume"] is True
        and chamber["facet_gram_sha256"] == PRODUCTION_FACET_GRAM_SHA256
    ):
        raise AssertionError("production H11/22 certificate changed")

    gram = sp.Matrix(production["fixed_lattice"]["gram"])
    roots = [
        sp.Matrix(record["primitive_fixed_coordinates"])
        for record in chamber["primitive_facet_lattice_records"]
    ]
    model = _standard_model()
    source_gram = sp.Matrix(FROZEN_SOURCE_GRAM)
    smith = smith_normal_form(gram, domain=ZZ)
    smith_absolute = [abs(int(smith[index, index])) for index in range(12)]
    q_histogram = _discriminant_q_histogram(gram)
    model_q_histogram = _discriminant_q_histogram(model)
    source_q_histogram = _discriminant_q_histogram(-source_gram)
    expected_q_histogram = {"0": 81, "2/3": 72, "4/3": 90}
    lattice_checks = {
        "candidate_rank_12_signature_1_11": exact_inertia(gram) == (1, 11, 0),
        "candidate_even": all(int(gram[i, i]) % 2 == 0 for i in range(12)),
        "candidate_determinant_minus_3_power_5": int(gram.det()) == -243,
        "candidate_smith_1_power_7_3_power_5": smith_absolute == [1] * 7 + [3] * 5,
        "model_U_plus_A2_minus_power_5": (
            exact_inertia(model) == (1, 11, 0) and int(model.det()) == -243
        ),
        "finite_discriminant_forms_match": (
            q_histogram == model_q_histogram == source_q_histogram == expected_q_histogram
        ),
        "source_gram_has_opposite_signature": exact_inertia(source_gram) == (11, 1, 0),
    }
    if not all(lattice_checks.values()):
        raise AssertionError(f"entry-27 lattice match failed: {lattice_checks}")

    # For a primitive reflective root r with r^2=-2k and divisibility d,
    # integrality gives k|d.  The exponent-three discriminant group gives
    # d in {1,3}.  The only extra formal possibility, (-2,div 3), would put
    # r/3 in A_L with q=-2/9=16/9 mod 2, absent from the enumerated form.
    excluded_q = str(_mod_two(sp.Rational(-2, 9)))
    if excluded_q in q_histogram:
        raise AssertionError("the forbidden (-2,div 3) discriminant value appeared")
    root_classification = {
        "reflection_integrality_rule": "for r^2=-2k and div(r)=d, k divides d",
        "discriminant_group_exponent": 3,
        "formal_primitive_types_before_discriminant_filter": [
            [-2, 1], [-2, 3], [-6, 3]
        ],
        "excluded_type": [-2, 3],
        "excluded_discriminant_value_mod_2Z": excluded_q,
        "discriminant_q_values_present_mod_2Z": list(q_histogram),
        "allowed_primitive_reflective_types_square_divisibility": [[-2, 1], [-6, 3]],
    }

    all_root = _vinberg_all_root_audit(
        gram, roots, production["vinberg"]["finite_volume"]
    )
    source_verification = (
        verify_source_archive(source_archive) if source_archive is not None else None
    )
    census_match = {
        "entry_one_based": 27,
        "list_entry": list(FROZEN_LIST_ENTRY),
        "source_gram": FROZEN_SOURCE_GRAM,
        "source_gram_determinant": int(source_gram.det()),
        "source_gram_signature_positive_negative": [11, 1],
        "candidate_elliptic_face_vector_matches_entry": (
            production["vinberg"]["elliptic_face_vector_rank_1_through_11"]
            == FROZEN_LIST_ENTRY[2]
        ),
        "candidate_cusp_count_matches_entry": (
            production["vinberg"]["ideal_vertex_count"] == FROZEN_LIST_ENTRY[3]
        ),
        "class_number_one_source_semantics": True,
        "identified_lattice": "II_(11,1)(3^-5), equivalently U+A2(-1)^5",
    }
    if not all(
        census_match[key]
        for key in (
            "candidate_elliptic_face_vector_matches_entry",
            "candidate_cusp_count_matches_entry",
            "class_number_one_source_semantics",
        )
    ):
        raise AssertionError("Scharlau--Kirschmer chamber census mismatch")

    checks = {
        "production_finite_volume_certificate_pinned": True,
        "entry_27_lattice_genus_match": all(lattice_checks.values()),
        "entry_27_face_and_cusp_census_match": True,
        "all_primitive_reflective_root_types_classified": True,
        "height_zero_reflective_root_set_empty": (
            all_root["height_zero"]["reflective_root_pair_count"] == 0
        ),
        "first_two_vinberg_shells_exact": (
            all_root["shell_roots_equal_displayed_simple_roots"]
        ),
        "vinberg_finite_volume_stopping": all_root["vinberg_stopping_applies"],
        "full_reflection_group_index_one": all_root["reflection_subgroup_index"] == 1,
    }
    return {
        "schema": SCHEMA,
        "id": "SS2-3B-H11-22-LINEAGE",
        "production_certificate": {
            "path": str(production_path),
            "sha256": _sha256_file(production_path),
            "facet_gram_sha256": PRODUCTION_FACET_GRAM_SHA256,
        },
        "scharlau_kirschmer_source": {
            "page_url": SOURCE_PAGE_URL,
            "archive_url": SOURCE_ARCHIVE_URL,
            "archive_sha256": SOURCE_ARCHIVE_SHA256,
            "list11_sha256": LIST11_SHA256,
            "gram11_sha256": GRAM11_SHA256,
            "optional_archive_verification": source_verification,
        },
        "lattice": {
            "candidate_gram": _matrix_payload(gram),
            "standard_model_gram": _matrix_payload(model),
            "candidate_smith_diagonal_absolute": smith_absolute,
            "discriminant_q_histogram": q_histogram,
            "checks": lattice_checks,
        },
        "census_match": census_match,
        "reflective_root_classification": root_classification,
        "all_root_vinberg_audit": all_root,
        "decision": {
            "full_reflection_chamber": True,
            "reflection_subgroup_index": 1,
            "abstract_chamber_lineage": "known Scharlau--Kirschmer List11/Gram11 entry 27",
            "literal_k3_fixed_section_priority": "not located in bounded source audit",
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--production", type=Path, default=PRODUCTION_CERTIFICATE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--verify-source-archive", type=Path)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate(args.production, args.verify_source_archive)
    if args.check_archive:
        archived = json.loads(args.output.read_text(encoding="utf-8"))
        if payload != archived:
            raise AssertionError("rebuilt lineage certificate differs from archive")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "id": payload["id"],
        "entry": payload["census_match"]["entry_one_based"],
        "index": payload["decision"]["reflection_subgroup_index"],
        "height_zero_reflective_pairs": payload["all_root_vinberg_audit"]["height_zero"]["reflective_root_pair_count"],
        "shell_root_counts": [
            shell["reflective_target_root_count"]
            for shell in payload["all_root_vinberg_audit"]["positive_shells"]
        ],
        "all_checks_pass": payload["all_checks_pass"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
