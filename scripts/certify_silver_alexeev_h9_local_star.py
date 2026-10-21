#!/usr/bin/env python3
"""Certify the complete local Weyl star at the silver bad ridge.

The frozen silver-Alexeev obstruction supplies a timelike point which is
orthogonal exactly to ambient walls 2, 4, and 5.  Those walls have spherical
type ``A1 + A2``.  Hence the twelve chambers in their finite standard
parabolic Weyl orbit are the complete ambient chamber star at that face.

This certificate restricts all twelve local chambers to the fixed silver
``sqrt(2)`` eigenspace.  Exact Gordan-alternative signs decide whether the
three transverse inequalities have a strict two-dimensional sector.  Every
one of the six full-dimensional sectors has a non-Coxeter ridge angle.
"""

from __future__ import annotations

import argparse
from collections import deque
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from coxeter_polytopes.exact.silver import SILVER_ZERO, Silver
from scripts.certify_silver_alexeev_h9_folding_obstruction import (
    _decode_integer_matrix,
    _decode_silver,
    _pair,
    _projected_pair,
    _reflection,
    verify as verify_base,
)


ROOT = Path(__file__).resolve().parents[1]
BASE_CERTIFICATE = ROOT / "data/candidates/silver_alexeev_h9_folding_obstruction.json"
OUTPUT = ROOT / "data/candidates/silver_alexeev_h9_local_star.json"
BASE_CERTIFICATE_SHA256 = (
    "d4a8502d51ea7b33225af999825d31a18506f0ac58ee86d96f4db1069daacb23"
)
EXPECTED_CERTIFICATE_SHA256 = "fe055505befad84e1cb0e9956b54fbdeb88184e645bbcc423a14e42ca655225e"
PARABOLIC_WALLS = (2, 4, 5)


def _canonical_sha256(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _certificate_digest(manifest: Mapping[str, object]) -> str:
    critical = dict(manifest)
    critical.pop("certificate_sha256", None)
    return _canonical_sha256(critical)


def _silver(value: Silver | int) -> dict[str, str]:
    return Silver.coerce(value).to_json()


def _silver_matrix(matrix: Sequence[Sequence[Silver]]) -> list[list[dict[str, str]]]:
    return [[_silver(value) for value in row] for row in matrix]


def _matrix_key(matrix: np.ndarray) -> tuple[int, ...]:
    return tuple(int(value) for value in matrix.flat)


def _word_matrix(
    generators: Mapping[int, np.ndarray],
    word: Sequence[int],
) -> np.ndarray:
    result = np.eye(20, dtype=object)
    for wall in word:
        result = result @ generators[wall]
    return result


def _local_group(
    generators: Mapping[int, np.ndarray],
) -> tuple[tuple[tuple[int, ...], np.ndarray], ...]:
    """Enumerate shortest lexicographic words in the finite local group."""

    identity = np.eye(20, dtype=object)
    seen: dict[tuple[int, ...], tuple[tuple[int, ...], np.ndarray]] = {
        _matrix_key(identity): ((), identity)
    }
    queue: deque[np.ndarray] = deque((identity,))
    while queue:
        element = queue.popleft()
        word = seen[_matrix_key(element)][0]
        for wall in PARABOLIC_WALLS:
            candidate = element @ generators[wall]
            key = _matrix_key(candidate)
            if key in seen:
                continue
            seen[key] = (word + (wall,), candidate)
            queue.append(candidate)
            if len(seen) > 12:
                raise AssertionError("the declared A1+A2 subgroup has order above 12")
    return tuple(sorted(seen.values(), key=lambda item: (len(item[0]), item[0])))


def _kernel_vector(matrix: Sequence[Sequence[Silver]]) -> tuple[Silver, Silver, Silver]:
    """Return a nonzero cross product of two rows of a rank-two 3x3 matrix."""

    for left, right in ((0, 1), (0, 2), (1, 2)):
        vector = (
            matrix[left][1] * matrix[right][2]
            - matrix[left][2] * matrix[right][1],
            matrix[left][2] * matrix[right][0]
            - matrix[left][0] * matrix[right][2],
            matrix[left][0] * matrix[right][1]
            - matrix[left][1] * matrix[right][0],
        )
        if any(vector):
            return vector
    raise AssertionError("the transverse projected Gram has rank below two")


def _kernel_checks(
    matrix: Sequence[Sequence[Silver]],
    kernel: Sequence[Silver],
) -> bool:
    annihilates = all(
        sum((matrix[row][column] * kernel[column] for column in range(3)), SILVER_ZERO)
        == 0
        for row in range(3)
    )
    negative_normals = all(matrix[index][index].sign() < 0 for index in range(3))
    rank_at_least_two = any(
        (
            matrix[left][left] * matrix[right][right]
            - matrix[left][right] ** 2
        ).sign()
        > 0
        for left in range(3)
        for right in range(left + 1, 3)
    )
    return annihilates and negative_normals and rank_at_least_two


def _sector_facets(kernel: Sequence[Silver]) -> tuple[bool, tuple[tuple[int, ...], ...]]:
    """Apply the strict Gordan alternative and identify the two extreme rays.

    The three projected normals span a negative-definite plane and have the
    displayed one-dimensional relation.  A strict transverse sector exists
    exactly when the relation has both signs.  In that case either a zero
    coefficient records two positively proportional duplicate walls, or the
    unique coefficient of one sign records the redundant middle normal.
    """

    signs = tuple(value.sign() for value in kernel)
    full = min(signs) < 0 < max(signs)
    if not full:
        return False, ()

    zero_indices = tuple(index for index, sign in enumerate(signs) if sign == 0)
    if zero_indices:
        if len(zero_indices) != 1:
            raise AssertionError("an unexpected transverse kernel zero pattern occurred")
        duplicate_group = tuple(
            PARABOLIC_WALLS[index]
            for index in range(3)
            if index not in zero_indices
        )
        other_group = (PARABOLIC_WALLS[zero_indices[0]],)
        return True, (duplicate_group, other_group)

    exceptional_sign = next(sign for sign in (-1, 1) if signs.count(sign) == 1)
    redundant = signs.index(exceptional_sign)
    return True, tuple(
        (PARABOLIC_WALLS[index],)
        for index in range(3)
        if index != redundant
    )


def _matrix_times_silver(
    matrix: np.ndarray,
    vector: Sequence[Silver],
) -> tuple[Silver, ...]:
    return tuple(
        sum(
            (Silver.coerce(int(matrix[row, column])) * vector[column] for column in range(20)),
            SILVER_ZERO,
        )
        for row in range(20)
    )


def _derive(base: Mapping[str, object]) -> dict[str, object]:
    base_result = verify_base(base)
    if base_result["certificate_sha256"] != BASE_CERTIFICATE_SHA256:
        raise ValueError("the local-star base certificate is not the frozen one")

    gram = _decode_integer_matrix(base.get("basis_gram"))
    roots = _decode_integer_matrix(base.get("wall_coordinates"))
    action = _decode_integer_matrix(base.get("framed_silver_action"))
    witness_value = base.get("ridge_witness")
    if not isinstance(witness_value, list) or len(witness_value) != 20:
        raise ValueError("the base ridge witness is malformed")
    witness = tuple(_decode_silver(value) for value in witness_value)

    generators = {
        wall: _reflection(
            roots[:, wall],
            gram,
            int(roots[:, wall] @ gram @ roots[:, wall]),
        )
        for wall in PARABOLIC_WALLS
    }
    group = _local_group(generators)
    sectors: list[dict[str, object]] = []
    all_kernel_checks = True
    all_words_replay = True
    all_witness_fixed = True
    all_face_pairings = True
    full_ratios: list[Silver] = []

    for word, element in group:
        all_words_replay &= np.array_equal(element, _word_matrix(generators, word))
        all_witness_fixed &= _matrix_times_silver(element, witness) == witness
        transformed_roots = element @ roots
        face_pairings = tuple(
            _pair(
                witness,
                gram,
                tuple(int(value) for value in transformed_roots[:, wall]),
            )
            for wall in range(transformed_roots.shape[1])
        )
        all_face_pairings &= (
            tuple(index for index, value in enumerate(face_pairings) if not value)
            == PARABOLIC_WALLS
            and all(
                value.sign() > 0
                for index, value in enumerate(face_pairings)
                if index not in PARABOLIC_WALLS
            )
        )
        projected = tuple(
            tuple(
                _projected_pair(
                    transformed_roots[:, left],
                    transformed_roots[:, right],
                    gram,
                    action,
                )
                for right in PARABOLIC_WALLS
            )
            for left in PARABOLIC_WALLS
        )
        kernel = _kernel_vector(projected)
        all_kernel_checks &= _kernel_checks(projected, kernel)
        full, facet_groups = _sector_facets(kernel)
        sector: dict[str, object] = {
            "word": list(word),
            "projected_gram_scaled_by_four": _silver_matrix(projected),
            "kernel_relation": [_silver(value) for value in kernel],
            "kernel_signs": [value.sign() for value in kernel],
            "full_dimensional_silver_sector": full,
            "facet_wall_groups": [list(group) for group in facet_groups],
        }
        if full:
            left = PARABOLIC_WALLS.index(facet_groups[0][0])
            right = PARABOLIC_WALLS.index(facet_groups[1][0])
            numerator = projected[left][right]
            ratio = numerator**2 / (projected[left][left] * projected[right][right])
            if numerator.sign() <= 0:
                raise AssertionError("a selected local facet pairing is not positive")
            if 0 < ratio < Silver(1, 0) / 4:
                interval = "between_m2_and_m3"
                lower_gap = ratio
                upper_gap = Silver(1, 0) / 4 - ratio
            elif Silver(1, 0) / 4 < ratio < Silver(1, 0) / 2:
                interval = "between_m3_and_m4"
                lower_gap = ratio - Silver(1, 0) / 4
                upper_gap = Silver(1, 0) / 2 - ratio
            else:
                raise AssertionError("a full local sector did not reproduce a bad angle")
            sector.update(
                {
                    "cosine_squared": _silver(ratio),
                    "coxeter_gap_interval": interval,
                    "lower_gap": _silver(lower_gap),
                    "upper_gap": _silver(upper_gap),
                }
            )
            full_ratios.append(ratio)
        sectors.append(sector)

    ratio_a = Silver(1, 0) / 4 + Silver(0, 1) / 8
    ratio_b = Silver(3, 0) / 8 - Silver(0, 1) / 4
    parabolic_gram = [
        [int(roots[:, left] @ gram @ roots[:, right]) for right in PARABOLIC_WALLS]
        for left in PARABOLIC_WALLS
    ]
    identity = np.eye(20, dtype=object)
    checks = {
        "parabolic_gram_is_A1_plus_A2": parabolic_gram
        == [[-2, 0, 0], [0, -2, 1], [0, 1, -2]],
        "parabolic_generator_involutions": all(
            np.array_equal(generator @ generator, identity)
            for generator in generators.values()
        ),
        "parabolic_group_order_12": len(group) == 12,
        "canonical_words_replay": all_words_replay,
        "common_ridge_witness_fixed_by_local_group": all_witness_fixed,
        "common_face_zero_set_and_other_wall_strictness": all_face_pairings,
        "all_projected_grams_have_exact_rank_two": all_kernel_checks,
        "six_full_and_six_nonfull_sectors": (
            sum(bool(sector["full_dimensional_silver_sector"]) for sector in sectors) == 6
        ),
        "all_six_full_sectors_have_non_coxeter_gap": len(full_ratios) == 6,
        "angle_histogram_is_four_plus_two": (
            full_ratios.count(ratio_a) == 4 and full_ratios.count(ratio_b) == 2
        ),
    }
    return {
        "parabolic_gram": parabolic_gram,
        "sectors": sectors,
        "checks": checks,
        "passes": all(checks.values()),
        "full_sector_count": 6,
        "nonfull_sector_count": 6,
        "angle_histogram": [
            {"cosine_squared": _silver(ratio_a), "count": 4},
            {"cosine_squared": _silver(ratio_b), "count": 2},
        ],
    }


def build() -> dict[str, object]:
    base = json.loads(BASE_CERTIFICATE.read_text(encoding="utf-8"))
    derived = _derive(base)
    if not derived["passes"]:
        raise AssertionError("the exact silver local-star checks failed")
    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "silver_alexeev_h9_local_star_obstruction",
        "claim": (
            "The complete finite A1+A2 standard-parabolic star at the certified silver ridge "
            "has twelve ambient chambers. Exactly six have full-dimensional "
            "H9 silver sectors, and all six sectors have a non-Coxeter ridge."
        ),
        "base_certificate": {
            "path": "data/candidates/silver_alexeev_h9_folding_obstruction.json",
            "certificate_sha256": BASE_CERTIFICATE_SHA256,
        },
        "parabolic_walls": list(PARABOLIC_WALLS),
        "parabolic_type": "A1 + A2",
        "parabolic_gram": derived["parabolic_gram"],
        "local_group_order": 12,
        "full_sector_count": derived["full_sector_count"],
        "nonfull_sector_count": derived["nonfull_sector_count"],
        "angle_histogram": derived["angle_histogram"],
        "sectors": derived["sectors"],
        "checks": derived["checks"],
        "scope": {
            "proved": (
                "All five alternative full-dimensional Weyl framings in the "
                "complete local standard-parabolic star, as well as the original "
                "framing, fail the Coxeter angle gate."
            ),
            "not_proved": [
                "Weyl framings outside this twelve-chamber local star are excluded.",
                "Every silver action on the (20,2,1) K3 lattice is excluded.",
                "Finite or infinite volume for any different silver root chamber.",
            ],
        },
    }
    manifest["certificate_sha256"] = _certificate_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    base = json.loads(BASE_CERTIFICATE.read_text(encoding="utf-8"))
    derived = _derive(base)
    expected_fields = {
        "parabolic_walls": list(PARABOLIC_WALLS),
        "parabolic_type": "A1 + A2",
        "parabolic_gram": derived["parabolic_gram"],
        "local_group_order": 12,
        "full_sector_count": derived["full_sector_count"],
        "nonfull_sector_count": derived["nonfull_sector_count"],
        "angle_histogram": derived["angle_histogram"],
        "sectors": derived["sectors"],
        "checks": derived["checks"],
    }
    for key, expected in expected_fields.items():
        if manifest.get(key) != expected:
            raise ValueError(f"the frozen local-star field {key} does not reproduce")
    if not derived["passes"]:
        raise ValueError("the exact silver local-star checks failed")
    digest = _certificate_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("silver local-star certificate digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("silver local-star certificate is not frozen")
    return {
        "local_group_order": 12,
        "full_sector_count": derived["full_sector_count"],
        "nonfull_sector_count": derived["nonfull_sector_count"],
        "angle_histogram": derived["angle_histogram"],
        "certificate_sha256": digest,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)
    if args.write:
        manifest = build()
        OUTPUT.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        print(json.dumps({"certificate_sha256": manifest["certificate_sha256"]}, indent=2))
        return
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    print(json.dumps(verify(manifest), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
