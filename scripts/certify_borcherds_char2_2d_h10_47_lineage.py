#!/usr/bin/env python3
"""Certify the arithmetic lineage and exact index of ``SS2-2D-H10-47``.

This supplement is deliberately separate from the production Vinberg-volume
certificate.  It gives an explicit integral isometry from one half of the
fixed divisor form to Vinberg's determinant-two form

    f_2^10 = <2> + <-1>^10,

pins Scharlau--Kirschmer ``List10/Gram10`` entry 2, and independently replays
Solomon's sum over every finite special subgroup of both the 47-wall chamber
and the classical 13-wall full reflection chamber.  Their Euler-characteristic
ratio is the exact reflection-subgroup index.
"""

from __future__ import annotations

import argparse
import ast
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any
import zipfile

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form

from scripts.certify_borcherds_char2_3b_h11_22 import _matrix_payload
from scripts.certify_jacobian_kummer_h10_reflection_index import (
    _connected_masks_direct,
    _euler_from_connected_components,
    _relation_matrix,
)
from scripts.certify_shimada_disc16_d5_h10_30_vinberg_faces import (
    DiagramEnumerator,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
PRODUCTION_CERTIFICATE = (
    ROOT / "data/enumerations/borcherds_char2_2d_h10_47_certificate.json"
)
MAXIMAL_CHAMBER_CERTIFICATE = (
    ROOT / "data/enumerations/jacobian_kummer_h10_maximal_reflection_lineage.json"
)
OUTPUT = ROOT / "data/enumerations/borcherds_char2_2d_h10_47_lineage.json"
DEFAULT_SOURCE_ARCHIVE = Path("/private/tmp/Scharlau-Kirschmer-Lorentz.zip")
SCHEMA = "coxeter.borcherds_char2_2d_h10_47_lineage.v1"

SOURCE_PAGE_URL = (
    "https://www.math.rwth-aachen.de/~Markus.Kirschmer/lorentz/index.html"
)
SOURCE_ARCHIVE_URL = (
    "https://www.math.rwth-aachen.de/~Markus.Kirschmer/lorentz/Lorentz.zip"
)
SOURCE_ARCHIVE_SHA256 = (
    "74a4dd3e534f1e793cba45765ff33575dbd331d5cdce6f34cdd64340b3f0d4a5"
)
LIST10_SHA256 = "9e627c4e1303015fd96acdcf6ccbeb264ebadba49abfcd761705f54d211b2e21"
GRAM10_SHA256 = "b8186e1fa2771c517a6040998f3b234568c3edde3ec8b7bb8412865d7710af18"
PRODUCTION_FACET_GRAM_SHA256 = (
    "454e34ef2358603578a4360145ac6bdd7abfe6e0b539344323ffa5e354038881"
)
MAXIMAL_CERTIFICATE_SHA256 = (
    "3cb8b0e9120d13479f9b9e3a32cfb63922e03f41f2f0fc1dc501b0a9f6654c29"
)

FROZEN_LIST_ENTRY = (
    2,
    "I_{10,1}([1^{10}2^{1}]_{1})",
    [13, 76, 264, 606, 966, 1092, 876, 488, 175, 29],
    3,
)

# Columns give a determinant-minus-one basis change.  If H is half the
# production fixed Gram, then C^T H C = diag(2,-1,...,-1).  The last ten
# columns are pairwise orthogonal primitive norm-minus-one vectors; the first
# is the primitive norm-two generator of their orthogonal complement.
EXPLICIT_F2_ISOMETRY = [
    [404, 8, 6, 16, 9, 14, 264, 18, 14, 9, 102],
    [-150, -3, -2, -6, -3, -5, -98, -7, -5, -3, -38],
    [-277, -6, -4, -11, -6, -9, -181, -13, -9, -6, -70],
    [-372, -8, -6, -15, -8, -12, -243, -18, -12, -8, -94],
    [-188, -4, -3, -8, -4, -6, -123, -9, -6, -4, -47],
    [620, 12, 9, 25, 14, 21, 405, 28, 21, 14, 157],
    [836, 16, 12, 33, 19, 28, 546, 38, 28, 19, 212],
    [1290, 25, 18, 51, 29, 42, 843, 59, 43, 29, 326],
    [1524, 30, 21, 60, 34, 50, 996, 70, 50, 34, 385],
    [-450, -8, -6, -18, -11, -15, -294, -20, -15, -10, -114],
    [756, 15, 10, 30, 17, 25, 494, 35, 25, 17, 191],
]


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256_file(path: Path) -> str:
    return _sha256_bytes(path.read_bytes())


def verify_source_archive(path: Path) -> dict[str, Any]:
    raw = path.read_bytes()
    if _sha256_bytes(raw) != SOURCE_ARCHIVE_SHA256:
        raise AssertionError("Lorentz.zip SHA-256 mismatch")
    with zipfile.ZipFile(path) as archive:
        list_raw = archive.read("List10")
        gram_raw = archive.read("Gram10")
        if _sha256_bytes(list_raw) != LIST10_SHA256:
            raise AssertionError("List10 SHA-256 mismatch")
        if _sha256_bytes(gram_raw) != GRAM10_SHA256:
            raise AssertionError("Gram10 SHA-256 mismatch")
        entries = ast.literal_eval(
            list_raw.decode("ascii").replace("<", "(").replace(">", ")")
        )
        grams = ast.literal_eval(gram_raw.decode("ascii"))
    if tuple(entries[1]) != FROZEN_LIST_ENTRY:
        raise AssertionError("List10 entry 2 changed")
    source_gram = sp.Matrix(11, 11, grams[1])
    expected = sp.diag(*([1] * 10), -2)
    if source_gram != expected:
        raise AssertionError("Gram10 entry 2 changed")
    return {
        "archive_path": str(path),
        "archive_sha256": SOURCE_ARCHIVE_SHA256,
        "list10_sha256": LIST10_SHA256,
        "gram10_sha256": GRAM10_SHA256,
        "entry_one_based": 2,
        "entry_reparsed_exactly": True,
    }


def _subgroup_solomon(positive_gram: sp.Matrix) -> dict[str, Any]:
    enumerator = DiagramEnumerator(positive_gram, dimension=10)
    connected_masks: set[int] = set()
    direct_counts = [1]
    connected_counts: list[int] = []
    for rank in range(1, 11):
        finite = enumerator.finite_subdiagrams(rank)
        direct_counts.append(len(finite))
        connected = {mask for mask, signature in finite.items() if len(signature) == 1}
        connected_masks.update(connected)
        connected_counts.append(len(connected))
    relations = _relation_matrix(
        [
            [str(positive_gram[row, column]) for column in range(positive_gram.cols)]
            for row in range(positive_gram.rows)
        ]
    )
    solomon = _euler_from_connected_components(relations, connected_masks)
    if solomon["finite_special_subset_count_by_rank"] != direct_counts:
        raise AssertionError("Solomon recursion disagrees with direct finite diagrams")
    solomon["independent_direct_finite_counts_by_rank"] = direct_counts
    solomon["connected_finite_counts_by_rank_1_through_10"] = connected_counts
    return solomon


def _maximal_solomon(relations: list[list[str]]) -> dict[str, Any]:
    connected = _connected_masks_direct(relations)
    return _euler_from_connected_components(relations, connected)


def build_certificate(
    production_path: Path = PRODUCTION_CERTIFICATE,
    maximal_path: Path = MAXIMAL_CHAMBER_CERTIFICATE,
    source_archive: Path | None = None,
) -> dict[str, Any]:
    production = json.loads(production_path.read_text(encoding="utf-8"))
    maximal = json.loads(maximal_path.read_text(encoding="utf-8"))
    if not (
        production["id"] == "SS2-2D-H10-47"
        and production["all_checks_pass"] is True
        and production["vinberg"]["finite_volume"] is True
        and production["restricted_chamber"]["facet_gram_sha256"]
        == PRODUCTION_FACET_GRAM_SHA256
    ):
        raise AssertionError("production H10/47 certificate changed")
    if _sha256_file(maximal_path) != MAXIMAL_CERTIFICATE_SHA256:
        raise AssertionError("pinned 13-wall maximal-chamber certificate changed")
    if not maximal["all_checks_pass"]:
        raise AssertionError("the pinned maximal-chamber Vinberg audit failed")

    fixed_gram = sp.Matrix(production["fixed_lattice"]["gram"])
    if any(value.q != 1 for value in fixed_gram / 2):
        raise AssertionError("the fixed Gram is no longer divisible by two")
    half_gram = fixed_gram / 2
    model = sp.diag(2, *([-1] * 10))
    isometry = sp.Matrix(EXPLICIT_F2_ISOMETRY)
    half_smith = smith_normal_form(half_gram, domain=ZZ)
    half_smith_absolute = [
        abs(int(half_smith[index, index])) for index in range(11)
    ]
    lattice_checks = {
        "fixed_gram_entrywise_divisible_by_2": True,
        "half_form_integral_odd": (
            all(value.q == 1 for value in half_gram)
            and any(int(half_gram[index, index]) % 2 for index in range(11))
        ),
        "half_form_signature_1_10": exact_inertia(half_gram) == (1, 10, 0),
        "half_form_determinant_2": int(half_gram.det()) == 2,
        "half_form_smith_1_power_10_2": half_smith_absolute == [1] * 10 + [2],
        "explicit_basis_change_unimodular": abs(int(isometry.det())) == 1,
        "explicit_congruence_to_f2_10": isometry.T * half_gram * isometry == model,
        "source_Gram10_entry_2_is_negative_model": (
            sp.diag(*([1] * 10), -2) == -model.extract(
                list(range(1, 11)) + [0], list(range(1, 11)) + [0]
            )
        ),
    }
    if not all(lattice_checks.values()):
        raise AssertionError(f"f_2^10 lattice identification failed: {lattice_checks}")

    positive_facet_gram = -sp.Matrix(
        production["restricted_chamber"]["facet_gram"]
    )
    subgroup = _subgroup_solomon(positive_facet_gram)
    maximal_relations = maximal["vinberg_maximal_integral_reflection_chamber"][
        "normalized_coxeter_relation_matrix"
    ]
    full_group = _maximal_solomon(maximal_relations)
    subgroup_chi = Fraction(subgroup["rational_euler_characteristic"])
    full_chi = Fraction(full_group["rational_euler_characteristic"])
    index = subgroup_chi / full_chi
    volume_coefficient = -subgroup_chi * Fraction(32, 945)
    expected_subgroup_counts = [
        1, 47, 540, 3230, 12215, 31212,
        54582, 64050, 47760, 19910, 3136,
    ]
    expected_full_counts = [1, 13, 76, 264, 606, 966, 1092, 876, 488, 175, 29]
    full_face_vector = [
        expected_subgroup_counts[10] + production["vinberg"]["ideal_vertex_count"],
        *reversed(expected_subgroup_counts[1:10]),
        1,
    ]
    solomon_checks = {
        "subgroup_finite_counts_expected": (
            subgroup["finite_special_subset_count_by_rank"] == expected_subgroup_counts
        ),
        "subgroup_connected_counts_expected": (
            subgroup["connected_finite_counts_by_rank_1_through_10"]
            == [47, 50, 120, 240, 300, 360, 420, 420, 240, 120]
        ),
        "subgroup_chi_expected": subgroup_chi == Fraction(-527, 120_960),
        "full_group_finite_counts_match_SK_entry": (
            full_group["finite_special_subset_count_by_rank"] == expected_full_counts
            and expected_full_counts[1:] == FROZEN_LIST_ENTRY[2]
        ),
        "full_group_chi_expected": full_chi == Fraction(-31, 11_147_673_600),
        "euler_ratio_positive_integer": index.denominator == 1 and index > 0,
        "index_expected": index == 1_566_720,
        "gauss_bonnet_volume_expected": volume_coefficient == Fraction(527, 3_572_100),
        "full_face_vector_expected": full_face_vector == [
            3352, 19910, 47760, 64050, 54582, 31212,
            12215, 3230, 540, 47, 1,
        ],
    }
    if not all(solomon_checks.values()):
        raise AssertionError(f"H10/47 Solomon/index audit failed: {solomon_checks}")

    source_verification = (
        verify_source_archive(source_archive) if source_archive is not None else None
    )
    checks = {
        "production_finite_volume_certificate_pinned": True,
        "explicit_integral_f2_10_isometry": all(lattice_checks.values()),
        "SK_entry_2_pinned": True,
        "maximal_13_wall_vinberg_certificate_pinned": True,
        "complete_subgroup_solomon_replay": all(solomon_checks.values()),
        "reflection_subgroup_index_certified": index == 1_566_720,
    }
    return {
        "schema": SCHEMA,
        "id": "SS2-2D-H10-47-LINEAGE",
        "production_certificate": {
            "path": str(production_path),
            "sha256": _sha256_file(production_path),
            "facet_gram_sha256": PRODUCTION_FACET_GRAM_SHA256,
        },
        "scharlau_kirschmer_source": {
            "page_url": SOURCE_PAGE_URL,
            "archive_url": SOURCE_ARCHIVE_URL,
            "archive_sha256": SOURCE_ARCHIVE_SHA256,
            "list10_sha256": LIST10_SHA256,
            "gram10_sha256": GRAM10_SHA256,
            "entry_one_based": 2,
            "list_entry": list(FROZEN_LIST_ENTRY),
            "optional_archive_verification": source_verification,
        },
        "lattice_isometry": {
            "fixed_gram": _matrix_payload(fixed_gram),
            "rescaling": "divide the bilinear form by 2",
            "half_gram": _matrix_payload(half_gram),
            "target_model": "f_2^10 = <2> + <-1>^10",
            "target_gram": _matrix_payload(model),
            "basis_change_columns": EXPLICIT_F2_ISOMETRY,
            "basis_change_determinant": int(isometry.det()),
            "half_gram_smith_diagonal_absolute": half_smith_absolute,
            "checks": lattice_checks,
        },
        "maximal_chamber_certificate": {
            "path": str(maximal_path),
            "sha256": MAXIMAL_CERTIFICATE_SHA256,
            "facet_count": 13,
            "finite_vertex_count": 29,
            "ideal_vertex_count": 3,
            "vinberg_stopping_fired": maximal[
                "vinberg_maximal_integral_reflection_chamber"
            ]["vinberg_stopping_fired"],
        },
        "solomon_index": {
            "method": "complete Solomon sum over every finite special subgroup",
            "subgroup_47_walls": subgroup,
            "full_f2_10_group_13_walls": full_group,
            "index": int(index),
            "index_prime_factorization": {"2": 11, "3": 2, "5": 1, "17": 1},
            "hyperbolic_covolume": "pi^5 * 527/3572100",
            "full_face_vector": full_face_vector,
            "checks": solomon_checks,
        },
        "decision": {
            "ambient_lattice_and_full_reflection_group": (
                "known: Vinberg f_2^10; Scharlau--Kirschmer entry 2"
            ),
            "literal_47_wall_subgroup_chamber": (
                "not located in the bounded ordinary/doubling CoxIter audit; worldwide priority open"
            ),
            "reflection_subgroup_index": int(index),
            "finite_volume_noncompact": True,
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--production", type=Path, default=PRODUCTION_CERTIFICATE)
    parser.add_argument("--maximal", type=Path, default=MAXIMAL_CHAMBER_CERTIFICATE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--verify-source-archive", type=Path)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate(
        args.production, args.maximal, args.verify_source_archive
    )
    if args.check_archive:
        archived = json.loads(args.output.read_text(encoding="utf-8"))
        if payload != archived:
            raise AssertionError("rebuilt H10/47 lineage differs from archive")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "id": payload["id"],
        "lattice": payload["lattice_isometry"]["target_model"],
        "index": payload["solomon_index"]["index"],
        "chi": payload["solomon_index"]["subgroup_47_walls"]["rational_euler_characteristic"],
        "all_checks_pass": payload["all_checks_pass"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
