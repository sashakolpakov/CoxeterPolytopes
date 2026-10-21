#!/usr/bin/env python3
"""Certify the exact lineage of the 328-wall char-3 V4 fixed section.

The divisor-fixed lattice is not itself preserved by every induced mirror.
The upstream mirror certificate constructs its canonical minimal
reflection-stable overlattice.  Here an explicit unimodular basis map
identifies that stable lattice with entry 20 of Scharlau--Kirschmer's
List11/Gram11 reflective-lattice census.  The script also gives an exact
discriminant-form obstruction to a tempting but incorrect orthogonal
decomposition and compares the rational quadratic space with McLeod's
f_3^11 space.

The 328-wall chamber is not enumerated again here.  Its finite-volume gate is
read from the independent exact GMP cone certificate.
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
import tarfile
from typing import Any, Sequence
import zipfile

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

from scripts.certify_higher_dimensional_rational_lineage import _compare
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
MIRROR_CERTIFICATE = (
    ROOT / "data/enumerations/ss_fermat_quartic_char3_h11_mirror_lattices.json"
)
V4_SCREEN = (
    ROOT / "data/enumerations/ss_fermat_quartic_char3_v4_fold_screen.json"
)
EXACT_CONE = (
    ROOT
    / "data/enumerations/ss_fermat_quartic_char3_v4_h11_exact_cone_certificate.json"
)
OUTPUT = (
    ROOT / "data/enumerations/ss_fermat_quartic_char3_v4_h11_lineage.json"
)

DEFAULT_SK_ARCHIVE = Path("/private/tmp/Scharlau-Kirschmer-Lorentz.zip")
DEFAULT_COXITER_CORPUS = Path("/private/tmp/coxiter-current/CoxIter-master")
COXITER_ARCHIVE = Path("/private/tmp/coxiter-current.zip")

SK_SOURCE_PAGE = "https://www.math.rwth-aachen.de/~Markus.Kirschmer/lorentz/index.html"
SK_ARCHIVE_URL = (
    "https://www.math.rwth-aachen.de/~Markus.Kirschmer/lorentz/Lorentz.zip"
)
SK_ARCHIVE_SHA256 = (
    "74a4dd3e534f1e793cba45765ff33575dbd331d5cdce6f34cdd64340b3f0d4a5"
)
LIST11_SHA256 = (
    "2ca0725232cb101e60380f14ad0ab7777478c60eff06991ba591b6d3c71e3184"
)
GRAM11_SHA256 = (
    "af8d52984fc9caabc00ffe670ae646d2151e68d141082544e18f4827ba142540"
)
SK_GRAM20_SHA256 = (
    "c854fc82d4602c0e6cb3752d76e0a195e51194532754feaab85d2b92b99e455f"
)
STABLE_TO_SK_SHA256 = (
    "72068a77b06d861f257ce2ef486f5ba1af2e4fbcd8c689e0fa8d23068b23095a"
)

MCLEOD_SOURCE = "https://arxiv.org/abs/1007.2299"
MCLEOD_COXITER_RELATIVE = Path(
    "graphs/group units qf/11-3,1,1,1,1,1,1,1,1,1,1,1.coxiter"
)
MCLEOD_COXITER_SHA256 = (
    "6e53b322d06c6ed678e2af4317290c4e3a6fc879d3e00efd17f3b87a9a3761dd"
)

COXITER_COMMIT = "b800be48d8240fbabd2d1c9542e61054584ad7fe"
COXITER_ARCHIVE_SHA256 = (
    "ae04648b291c1e5f3d0e47d666c846dc805b244469a20486bca59674b4d171c3"
)
DOUBLING_ARCHIVE_SHA256 = (
    "3330ca41fb3e047123ace718f08542bc1dd0019b35922d4fa836a217647e7036"
)


FROZEN_LIST11_ENTRY20 = (
    48,
    "II_{11,1}(2^{4} 3^{-1})",
    [22, 204, 1089, 3789, 9126, 15660, 19295, 16851, 9973, 3564, 543],
    22,
)

# Scharlau--Kirschmer Gram11 entry 20, in their signature-(11,1)
# convention.  The candidate convention is its negative, of signature (1,11).
SK_GRAM20 = sp.Matrix(
    [
        [0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 2, -1, -1, -1, -1, -1, -1, 0, 1, -1],
        [0, 0, -1, 2, 1, 0, 1, 0, 1, 0, 0, 0],
        [0, 0, -1, 1, 2, 1, 0, 0, 1, 0, 0, 0],
        [0, 0, -1, 0, 1, 2, 0, 0, 0, 0, -1, 1],
        [0, 0, -1, 1, 0, 0, 2, 0, 0, 0, 0, 0],
        [0, 0, -1, 0, 0, 0, 0, 2, 1, 0, -1, 1],
        [0, 0, -1, 1, 1, 0, 0, 1, 2, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 1],
        [0, 0, 1, 0, 0, -1, 0, -1, 0, -1, 4, -3],
        [0, 0, -1, 0, 0, 1, 0, 1, 0, 1, -3, 4],
    ]
)

# Columns are the Scharlau--Kirschmer entry-20 basis in stable-lattice
# coordinates.  Thus R^T G_stable R = -Gram11[20].
STABLE_TO_SK = sp.Matrix(
    [
        [0, -1, -1, 2, 1, 0, 0, 1, 0, 0, -1, 1],
        [0, 1, 0, -1, -1, 0, 0, -1, 1, 0, 1, 0],
        [5, -6, 0, 0, 1, 1, -3, 3, -1, 2, -10, 6],
        [2, 0, 0, -1, 0, 1, -1, -1, 0, 1, -2, 2],
        [0, 1, 0, 1, 0, 0, 0, -1, 0, 0, 1, 1],
        [-2, 2, 0, 2, 0, -1, 1, -1, 0, -1, 4, -2],
        [2, -2, 0, -1, 1, 1, -1, 1, 0, 1, -4, 1],
        [0, 0, 0, -2, -3, 0, 1, -1, 0, -1, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0],
        [1, -2, 0, -2, -2, 0, 0, 1, 0, 1, -3, 2],
        [1, 0, 0, -1, -1, 0, 0, 0, 0, 0, -1, 1],
        [0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0],
    ]
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _matrix_payload(matrix: sp.Matrix) -> list[list[int]]:
    if any(value.is_Integer is not True for value in matrix):
        raise ArithmeticError("expected an integral matrix")
    return [[int(value) for value in row] for row in matrix.tolist()]


def _matrix_sha256(matrix: sp.Matrix) -> str:
    payload = ";".join(
        ",".join(str(int(matrix[row, column])) for column in range(matrix.cols))
        for row in range(matrix.rows)
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def _smith_diagonal(matrix: sp.Matrix) -> list[int]:
    smith = smith_normal_form(matrix, domain=sp.ZZ)
    return [abs(int(smith[index, index])) for index in range(matrix.rows)]


def _d8_cartan() -> sp.Matrix:
    matrix = 2 * sp.eye(8)
    for left, right in (
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 6),
        (5, 7),
    ):
        matrix[left, right] = matrix[right, left] = -1
    return matrix


def _tempting_decomposition_gram() -> sp.Matrix:
    return sp.diag(
        sp.Matrix([[0, 1], [1, 0]]),
        -_d8_cartan(),
        -2 * sp.Matrix([[2, -1], [-1, 2]]),
    )


def _mod_one(value: sp.Rational) -> sp.Rational:
    value = sp.Rational(value)
    return value - sp.floor(value)


def _mod_two(value: sp.Rational) -> sp.Rational:
    value = sp.Rational(value)
    return value - 2 * sp.floor(value / 2)


def _discriminant_q_histogram(gram: sp.Matrix) -> dict[str, int]:
    """Enumerate q on A_L = L^*/L exactly; only 48 classes occur here."""

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


def _root_type(gram: sp.Matrix, root: sp.Matrix) -> tuple[int, int]:
    square = int((root.T * gram * root)[0])
    divisibility = math.gcd(*(abs(int(value)) for value in gram * root))
    return square, divisibility


def _reflection_is_integral(gram: sp.Matrix, root: sp.Matrix) -> bool:
    square = int((root.T * gram * root)[0])
    reflection = sp.eye(gram.rows) - sp.Rational(2, square) * root * (
        root.T * gram
    )
    return (
        all(value.is_Integer is True for value in reflection)
        and reflection.T * gram * reflection == gram
        and reflection * reflection == sp.eye(gram.rows)
        and reflection * root == -root
    )


def _replay_328_mirrors(
    gram: sp.Matrix, records: Sequence[dict[str, Any]]
) -> tuple[list[dict[str, Any]], Counter[tuple[int, int]]]:
    inverse_isometry = STABLE_TO_SK.inv()
    source_gram_in_candidate_sign = -SK_GRAM20
    output: list[dict[str, Any]] = []
    types: Counter[tuple[int, int]] = Counter()
    for source in records:
        root = sp.Matrix(source["stable_lattice_coordinates"])
        source_coordinates = inverse_isometry * root
        if any(value.is_Integer is not True for value in source_coordinates):
            raise ArithmeticError("a mirror left the Scharlau--Kirschmer lattice")
        square, divisibility = _root_type(gram, root)
        if (
            math.gcd(*(abs(int(value)) for value in root)) != 1
            or [square, divisibility]
            != [int(source["square"]), int(source["divisibility"])]
            or not _reflection_is_integral(gram, root)
            or not _reflection_is_integral(
                source_gram_in_candidate_sign, source_coordinates
            )
            or (source_coordinates.T * source_gram_in_candidate_sign * source_coordinates)[
                0
            ]
            != square
        ):
            raise ArithmeticError("a stable-lattice mirror replay failed")
        types[(square, divisibility)] += 1
        output.append(
            {
                "stable_lattice_coordinates": [int(value) for value in root],
                "sk_entry20_coordinates": [
                    int(value) for value in source_coordinates
                ],
                "candidate_sign_square": square,
                "divisibility": divisibility,
                "integral_involutive_isometry_in_both_bases": True,
            }
        )
    return output, types


def _first_header(lines: Sequence[str]) -> tuple[int, int]:
    content = next(
        stripped
        for line in lines
        if (stripped := line.split("#", 1)[0].strip())
    )
    facets, dimension = map(int, content.split()[:2])
    return facets, dimension


def _verify_external_sources(sk_archive: Path, coxiter_corpus: Path) -> None:
    if _sha256(sk_archive) != SK_ARCHIVE_SHA256:
        raise ArithmeticError("Scharlau--Kirschmer archive hash mismatch")
    with zipfile.ZipFile(sk_archive) as archive:
        list_raw = archive.read("List11")
        gram_raw = archive.read("Gram11")
    if _sha256_bytes(list_raw) != LIST11_SHA256:
        raise ArithmeticError("List11 hash mismatch")
    if _sha256_bytes(gram_raw) != GRAM11_SHA256:
        raise ArithmeticError("Gram11 hash mismatch")
    entries = ast.literal_eval(
        list_raw.decode("ascii").replace("<", "(").replace(">", ")")
    )
    grams = ast.literal_eval(gram_raw.decode("ascii"))
    if tuple(entries[19]) != FROZEN_LIST11_ENTRY20:
        raise ArithmeticError("List11 entry 20 differs from the freeze")
    if sp.Matrix(12, 12, grams[19]) != SK_GRAM20:
        raise ArithmeticError("Gram11 entry 20 differs from the freeze")

    if _sha256(COXITER_ARCHIVE) != COXITER_ARCHIVE_SHA256:
        raise ArithmeticError("CoxIter archive hash mismatch")
    mcleod = coxiter_corpus / MCLEOD_COXITER_RELATIVE
    if _sha256(mcleod) != MCLEOD_COXITER_SHA256:
        raise ArithmeticError("McLeod H11 CoxIter file hash mismatch")
    if _first_header(mcleod.read_text(encoding="utf-8").splitlines()) != (15, 11):
        raise ArithmeticError("McLeod H11 CoxIter header changed")

    ordinary: Counter[tuple[int, int]] = Counter()
    ordinary_total = 0
    for path in coxiter_corpus.rglob("*.coxiter"):
        facets, dimension = _first_header(
            path.read_text(encoding="utf-8").splitlines()
        )
        ordinary[(dimension, facets)] += 1
        ordinary_total += 1
    if ordinary_total != 1_214 or Counter(
        {key: value for key, value in ordinary.items() if key[0] == 11}
    ) != Counter({(11, 12): 14, (11, 13): 6, (11, 15): 2}):
        raise ArithmeticError("ordinary CoxIter H11 census changed")

    doubling_path = coxiter_corpus / "graphs/testing_index2.tar.gz"
    if _sha256(doubling_path) != DOUBLING_ARCHIVE_SHA256:
        raise ArithmeticError("CoxIter doubling archive hash mismatch")
    doubling: Counter[tuple[int, int]] = Counter()
    doubling_names_22: list[str] = []
    doubling_total = 0
    with tarfile.open(doubling_path, "r:gz") as archive:
        for member in archive:
            if not member.isfile() or not member.name.endswith(".coxiter"):
                continue
            handle = archive.extractfile(member)
            if handle is None:
                raise ArithmeticError("failed to read a doubling record")
            facets, dimension = _first_header(
                [raw.decode("utf-8") for raw in handle]
            )
            doubling[(dimension, facets)] += 1
            doubling_total += 1
            if (dimension, facets) == (11, 22):
                doubling_names_22.append(member.name)
    expected_doubling = Counter(
        {
            (11, 15): 2,
            (11, 16): 2,
            (11, 18): 2,
            (11, 22): 2,
            (11, 30): 2,
            (11, 46): 2,
        }
    )
    expected_names = [
        "testing_index2/11-2,1,1,1,1,1,1,1,1,1,1,1-15-11_3.coxiter",
        "testing_index2/11-2,1,1,1,1,1,1,1,1,1,1,1-11-15_3.coxiter",
    ]
    if (
        doubling_total != 2_959
        or Counter(
            {key: value for key, value in doubling.items() if key[0] == 11}
        )
        != expected_doubling
        or doubling_names_22 != expected_names
    ):
        raise ArithmeticError("doubling CoxIter H11 census changed")


def build_certificate() -> dict[str, Any]:
    mirror_payload = json.loads(MIRROR_CERTIFICATE.read_text(encoding="utf-8"))
    v4 = mirror_payload["v4_h11"]
    stable = v4["minimal_reflection_stable_overlattice"]
    gram = sp.Matrix(stable["gram"])
    exact_cone = json.loads(EXACT_CONE.read_text(encoding="utf-8"))
    screen = json.loads(V4_SCREEN.read_text(encoding="utf-8"))
    h11_screen = next(
        row
        for row in screen["casts"]
        if int(row["subgroup"]["hyperbolic_dimension"]) == 11
    )

    mirrors, mirror_types = _replay_328_mirrors(
        gram, stable["mirrors_in_stable_lattice"]
    )
    target_guess = _tempting_decomposition_gram()
    stable_q_histogram = _discriminant_q_histogram(gram)
    target_q_histogram = _discriminant_q_histogram(target_guess)
    rational_mcleod = _compare(
        gram, sp.diag(3, *([-1] * 11)), Fraction(1)
    )

    list_faces = FROZEN_LIST11_ENTRY20[2]
    finite_face_vector_high_to_low = list(list_faces)
    full_face_vector_low_to_high = [
        list_faces[-1] + FROZEN_LIST11_ENTRY20[3],
        *reversed(list_faces[:-1]),
        1,
    ]
    expected_smith = [*([1] * 8), 2, 2, 2, 6]
    checks = {
        "upstream_mirror_certificate_passes": mirror_payload["all_checks_pass"],
        "fixed_divisor_lattice_has_36_nonintegral_reflections": (
            v4["nonintegral_reflection_count_on_fixed_lattice"] == 36
        ),
        "minimal_stable_lattice_is_index4_even_det48_H11": (
            stable["overlattice_index"] == 4
            and stable["even"]
            and stable["determinant_absolute"] == 48
            and stable["signature"] == [1, 11]
            and stable["smith_diagonal_absolute"] == expected_smith
            and exact_inertia(gram) == (1, 11, 0)
            and _smith_diagonal(gram) == expected_smith
        ),
        "frozen_sk_gram_hash_matches": (
            _matrix_sha256(SK_GRAM20) == SK_GRAM20_SHA256
        ),
        "explicit_unimodular_sk_entry20_isometry": (
            int(STABLE_TO_SK.det()) == -1
            and STABLE_TO_SK.T * gram * STABLE_TO_SK == -SK_GRAM20
            and _matrix_sha256(STABLE_TO_SK) == STABLE_TO_SK_SHA256
        ),
        "exactly_328_distinct_primitive_integral_reflections": (
            len(mirrors) == 328
            and len(
                {tuple(record["stable_lattice_coordinates"]) for record in mirrors}
            )
            == 328
            and mirror_types
            == Counter({(-6, 3): 276, (-2, 1): 28, (-4, 2): 24})
        ),
        "exact_cone_proves_finite_volume_noncompact": (
            exact_cone["decision"]
            == {
                "classification": "finite_volume_noncompact",
                "compact": False,
                "criterion": (
                    "complete pointed chamber has no spacelike/hyperideal extreme ray"
                ),
                "finite_volume": True,
            }
            and exact_cone["exact_double_description"]["finite_vertex_count"]
            == 225_408
            and exact_cone["exact_double_description"]["ideal_vertex_count"]
            == 11_554
            and exact_cone["exact_double_description"]["hyperideal_vertex_count"]
            == 0
        ),
        "strict_328_facet_coxeter_gate_passes": (
            h11_screen["strict_facet_gate"]["certified_strict_facet_count"] == 328
            and h11_screen["strict_facet_gate"]["bad_pair_count"] == 0
        ),
        "tempting_decomposition_rejected_by_discriminant_form": (
            abs(int(target_guess.det())) == 48
            and _smith_diagonal(target_guess) == expected_smith
            and stable_q_histogram
            == {"0": 10, "1/3": 12, "1": 6, "4/3": 20}
            and target_q_histogram
            == {"0": 6, "2/3": 12, "1": 10, "5/3": 20}
            and stable_q_histogram != target_q_histogram
        ),
        "rational_space_matches_mcleod_f3_H11": rational_mcleod[
            "rationally_similar"
        ],
        "alexeev_2_elementary_census_is_inapplicable": (
            stable["smith_diagonal_absolute"][-1] == 6
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"V4-H11 lineage checks failed: {checks}")

    return {
        "schema": "ss_fermat_quartic_char3_v4_h11_lineage_v1",
        "candidate": {
            "id": "F3-V4-H11-328",
            "dimension": 11,
            "facet_count": 328,
            "finite_vertex_count": 225_408,
            "ideal_vertex_count": 11_554,
            "hyperideal_vertex_count": 0,
            "volume_decision": "finite_volume_noncompact",
        },
        "divisor_fixed_lattice_warning": {
            "fixed_lattice_id": "char3-V4-H11",
            "determinant_absolute": 768,
            "nonintegral_facet_reflection_count": 36,
            "conclusion": (
                "the saturated divisor-fixed lattice itself is not the integral "
                "reflection lattice of all 328 mirrors"
            ),
        },
        "minimal_reflection_stable_lattice": {
            "construction": stable["construction"],
            "minimality_reason": stable["minimality_reason"],
            "index_over_divisor_fixed_lattice": 4,
            "gram": _matrix_payload(gram),
            "signature": [1, 11],
            "even": True,
            "determinant_absolute": 48,
            "smith_diagonal_absolute": expected_smith,
            "discriminant_group": "(Z/2)^3 + Z/6",
            "primitive_mirror_type_histogram": {
                f"{square},{divisibility}": count
                for (square, divisibility), count in sorted(mirror_types.items())
            },
        },
        "known_maximal_integral_reflection_chamber": {
            "source": {
                "authors": "Rudolf Scharlau and Markus Kirschmer",
                "maintained_census_page": SK_SOURCE_PAGE,
                "archive_url": SK_ARCHIVE_URL,
                "archive_sha256": SK_ARCHIVE_SHA256,
                "list11_sha256": LIST11_SHA256,
                "gram11_sha256": GRAM11_SHA256,
                "entry_one_based": 20,
                "genus_symbol": FROZEN_LIST11_ENTRY20[1],
                "source_scope": (
                    "strongly squarefree reflective Lorentzian lattices of "
                    "signature (n,1); the maintained page says each listed genus "
                    "has class number one"
                ),
            },
            "source_gram_signature_11_1": _matrix_payload(SK_GRAM20),
            "source_gram_sha256": _matrix_sha256(SK_GRAM20),
            "identity": "R^T * G_stable * R = -Gram11[20]",
            "explicit_unimodular_isometry_R": _matrix_payload(STABLE_TO_SK),
            "explicit_isometry_determinant": int(STABLE_TO_SK.det()),
            "explicit_isometry_sha256": _matrix_sha256(STABLE_TO_SK),
            "facet_count": list_faces[0],
            "finite_face_counts_from_facets_to_finite_vertices": (
                finite_face_vector_high_to_low
            ),
            "finite_vertex_count": list_faces[-1],
            "cusp_count": FROZEN_LIST11_ENTRY20[3],
            "full_face_vector_f0_through_polytope": full_face_vector_low_to_high,
            "maximality_scope": (
                "W(L_stable), the reflection group generated by every integral "
                "root reflection of this lattice"
            ),
        },
        "induced_328_wall_reflection_subgroup": {
            "facet_count": len(mirrors),
            "primitive_square_divisibility_histogram": {
                f"{square},{divisibility}": count
                for (square, divisibility), count in sorted(mirror_types.items())
            },
            "all_reflections_integral_involutive_isometries": True,
            "mirrors_in_stable_and_sk_bases": mirrors,
            "subgroup_conclusion": (
                "W_328 is a literal reflection subgroup of the known full "
                "22-wall Scharlau--Kirschmer entry-20 reflection group"
            ),
            "finite_index_reason": (
                "both the 328-wall chamber and the full 22-wall chamber have "
                "finite hyperbolic covolume"
            ),
            "numeric_index_status": (
                "not computed: in H11 both rational Euler characteristics vanish, "
                "so a Solomon-Euler ratio is 0/0; use explicit cosets/tiling or an "
                "odd-dimensional hyperbolic-volume comparison"
            ),
        },
        "near_miss_decomposition": {
            "candidate": "U + D8(-1) + A2(-2)",
            "same_rank_determinant_and_smith_diagonal": True,
            "candidate_gram": _matrix_payload(target_guess),
            "stable_discriminant_q_mod_2Z_histogram": stable_q_histogram,
            "candidate_discriminant_q_mod_2Z_histogram": target_q_histogram,
            "integrally_isometric": False,
            "obstruction": "the exact discriminant quadratic forms differ",
        },
        "rational_and_census_comparisons": {
            "mcleod_f3_H11": {
                "source": MCLEOD_SOURCE,
                "reference_form_candidate_sign": "diag(3,-1^11)",
                "maximal_reflection_chamber_facets": 15,
                "relationship": (
                    "same rational quadratic space and hence the same arithmetic "
                    "commensurability class, but not the same integral lattice "
                    "(determinants 48 and 3)"
                ),
                "exact_rational_comparison": rational_mcleod,
            },
            "alexeev_2_elementary_census": {
                "applicable": False,
                "reason": (
                    "the discriminant group has a nontrivial 3-primary part, so "
                    "the lattice is not 2-elementary"
                ),
            },
        },
        "bounded_weighted_diagram_audit": {
            "coxiter_snapshot": {
                "repository": "https://github.com/rgugliel/CoxIter",
                "commit": COXITER_COMMIT,
                "archive_sha256": COXITER_ARCHIVE_SHA256,
                "ordinary_file_count": 1_214,
                "ordinary_H11_facet_histogram": {"12": 14, "13": 6, "15": 2},
                "ordinary_global_maximum_facets": 50,
                "mcleod_H11_record": MCLEOD_COXITER_RELATIVE.as_posix(),
                "mcleod_H11_record_sha256": MCLEOD_COXITER_SHA256,
                "doubling_archive_sha256": DOUBLING_ARCHIVE_SHA256,
                "doubling_file_count": 2_959,
                "doubling_H11_facet_histogram": {
                    "15": 2,
                    "16": 2,
                    "18": 2,
                    "22": 2,
                    "30": 2,
                    "46": 2,
                },
                "doubling_H11_22_facet_records": [
                    (
                        "testing_index2/11-2,1,1,1,1,1,1,1,1,1,1,1-"
                        "15-11_3.coxiter"
                    ),
                    (
                        "testing_index2/11-2,1,1,1,1,1,1,1,1,1,1,1-"
                        "11-15_3.coxiter"
                    ),
                ],
                "doubling_22_facet_provenance_note": (
                    "these are index-two doubling descendants of the qf2 H11 "
                    "record, not source records for the determinant-48 SK lattice"
                ),
                "doubling_global_maximum_facets": 264,
            },
            "literal_328_facet_record_present": False,
            "scope_warning": (
                "facet-count exclusion proves absence only from these two pinned "
                "corpora; it is not a worldwide subgroup-chamber priority proof"
            ),
        },
        "lineage_decision": {
            "new_stable_integral_reflection_lattice": False,
            "new_maximal_reflection_group": False,
            "new_arithmetic_commensurability_class": False,
            "positive_algebraic_construction": True,
            "literal_328_wall_subgroup_chamber_priority": (
                "open after bounded corpus exclusion"
            ),
        },
        "upstream_artifacts": {
            str(path.relative_to(ROOT)): _sha256(path)
            for path in (MIRROR_CERTIFICATE, V4_SCREEN, EXACT_CONE)
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--verify-external", action="store_true")
    parser.add_argument("--sk-archive", type=Path, default=DEFAULT_SK_ARCHIVE)
    parser.add_argument(
        "--coxiter-corpus", type=Path, default=DEFAULT_COXITER_CORPUS
    )
    arguments = parser.parse_args()
    if arguments.verify_external:
        _verify_external_sources(arguments.sk_archive, arguments.coxiter_corpus)
    result = build_certificate()
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if (
            not arguments.output.is_file()
            or arguments.output.read_text(encoding="utf-8") != rendered
        ):
            raise SystemExit(f"archived artifact differs: {arguments.output}")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
    print(arguments.output)


if __name__ == "__main__":
    main()
