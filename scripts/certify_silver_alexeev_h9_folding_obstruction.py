#!/usr/bin/env python3
"""Freeze the exact non-Coxeter obstruction for the silver Alexeev slice.

The construction starts with Alexeev's 50-wall chamber for the unique
2-elementary lattice ``(20,2,1)``.  The silver action is the direct trace of
the saturated one-4 C5 Lanner block, F4, and A1.  After a 108-reflection Weyl
framing, two genuine restricted facets meet at an angle which is not pi/m.

Building the archive needs Alexeev's v4 ancillary file.  Verifying an already
frozen archive is self-contained: the basis Gram, all wall coordinates, both
silver actions, the framing word, and an exact ridge witness are stored.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
from sympy import Matrix

from coxeter_polytopes.exact.silver import SILVER_ZERO, SQRT2, Silver
from scripts.probe_silver_alexeev_h9_slice import (
    _alexeev_matrix,
    _reflection,
    _wall_coordinates,
    silver_action_in_alexeev_basis,
    weyl_frame,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/candidates/silver_alexeev_h9_folding_obstruction.json"
ALEXEEV_SOURCE_SHA256 = (
    "4cdf39f61f01238a883129551c0181427e5821b7fed2fabfea436cbf3627f457"
)
EXPECTED_CERTIFICATE_SHA256 = "d4a8502d51ea7b33225af999825d31a18506f0ac58ee86d96f4db1069daacb23"
PLUS_BASIS_COLUMNS = (0, 1, 2, 3, 4, 9, 10, 11, 12, 17)
RIDGE_WALLS = (2, 4)
DUPLICATE_WALL = 5
RIDGE_COEFFICIENTS = (
    Silver(-69, -49),
    Silver(-26, -19),
    Silver(-96, -68),
    Silver(-74, -53),
    Silver(-48, -35),
    Silver(121, 84),
    Silver(142, 100),
    Silver(90, 63),
    Silver(138, 97),
    Silver(-1, Fraction(-1, 2)),
)


def _canonical_sha256(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _certificate_digest(manifest: Mapping[str, object]) -> str:
    critical = dict(manifest)
    critical.pop("certificate_sha256", None)
    return _canonical_sha256(critical)


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _silver(value: Silver | int | Fraction) -> dict[str, str]:
    return Silver.coerce(value).to_json()


def _silver_matrix(matrix: Sequence[Sequence[Silver]]) -> list[list[dict[str, str]]]:
    return [[_silver(value) for value in row] for row in matrix]


def _integer_matrix(matrix: np.ndarray) -> list[list[int]]:
    return [[int(value) for value in row] for row in matrix]


def _decode_integer_matrix(value: object) -> np.ndarray:
    if not isinstance(value, list) or not value:
        raise ValueError("an integer matrix must be a nonempty row list")
    if any(not isinstance(row, list) for row in value):
        raise ValueError("an integer matrix must be a row list")
    width = len(value[0])
    if any(len(row) != width for row in value):
        raise ValueError("an integer matrix is ragged")
    if any(not isinstance(entry, int) for row in value for entry in row):
        raise ValueError("an integer matrix has a noninteger entry")
    return np.array(value, dtype=object)


def _decode_silver(value: object) -> Silver:
    if not isinstance(value, Mapping):
        raise ValueError("a silver entry must be an object")
    return Silver.from_json(value)


def _plus_basis(action: np.ndarray) -> tuple[tuple[Silver, ...], ...]:
    columns: list[tuple[Silver, ...]] = []
    for column in PLUS_BASIS_COLUMNS:
        columns.append(
            tuple(
                Silver(1 if row == column else 0, Fraction(int(action[row, column]), 2))
                for row in range(20)
            )
        )
    return tuple(columns)


def _ambient_witness(action: np.ndarray) -> tuple[Silver, ...]:
    basis = _plus_basis(action)
    return tuple(
        sum(
            (basis[column][row] * RIDGE_COEFFICIENTS[column] for column in range(10)),
            SILVER_ZERO,
        )
        for row in range(20)
    )


def _pair(
    left: Sequence[Silver],
    gram: np.ndarray,
    right: Sequence[int] | Sequence[Silver],
) -> Silver:
    return sum(
        (
            left[row]
            * int(gram[row, column])
            * Silver.coerce(right[column])
            for row in range(20)
            for column in range(20)
        ),
        SILVER_ZERO,
    )


def _framed_action(
    gram: np.ndarray,
    roots: np.ndarray,
    unframed: np.ndarray,
    word: Sequence[int],
) -> tuple[np.ndarray, np.ndarray]:
    weyl = np.eye(20, dtype=object)
    for wall in word:
        root = roots[:, wall]
        square = int(root @ gram @ root)
        weyl = _reflection(root, gram, square) @ weyl
    inverse = Matrix(weyl.tolist()).inv()
    if any(Fraction(entry).denominator != 1 for entry in inverse):
        raise AssertionError("the framing word is not unimodular")
    framed = Matrix(weyl.tolist()) * Matrix(unframed.tolist()) * inverse
    return (
        np.array([[int(entry) for entry in row] for row in framed.tolist()], dtype=object),
        weyl,
    )


def _projected_pair(
    left: np.ndarray,
    right: np.ndarray,
    gram: np.ndarray,
    action: np.ndarray,
) -> Silver:
    # Four times the Gram pairing of the positive-eigenvalue projections.
    q_value = int(left @ gram @ right)
    b_value = int(left @ gram @ action @ right)
    return Silver(2 * q_value, b_value)


def _derived_checks(
    gram: np.ndarray,
    roots: np.ndarray,
    unframed: np.ndarray,
    framed: np.ndarray,
    word: Sequence[int],
) -> dict[str, object]:
    identity = np.eye(20, dtype=object)
    witness = _ambient_witness(framed)
    pairings = tuple(
        _pair(witness, gram, tuple(int(value) for value in roots[:, wall]))
        for wall in range(roots.shape[1])
    )
    witness_square = _pair(witness, gram, witness)
    projected = tuple(
        tuple(
            _projected_pair(roots[:, left], roots[:, right], gram, framed)
            for right in RIDGE_WALLS
        )
        for left in RIDGE_WALLS
    )
    angle_ratio = projected[0][1] ** 2 / (projected[0][0] * projected[1][1])
    lower_gap = angle_ratio - Fraction(1, 4)
    upper_gap = Fraction(1, 2) - angle_ratio
    _, weyl = _framed_action(gram, roots, unframed, word)
    checks = {
        "basis_gram_symmetric": np.array_equal(gram, gram.T),
        "wall_gram_reconstructed": all(
            int(roots[:, index] @ gram @ roots[:, index]) in (-2, -4)
            for index in range(roots.shape[1])
        ),
        "unframed_action_square_two": np.array_equal(unframed @ unframed, 2 * identity),
        "unframed_action_self_adjoint": np.array_equal(unframed.T @ gram, gram @ unframed),
        "framing_word_length_108": len(word) == 108,
        "framing_isometry": np.array_equal(weyl.T @ gram @ weyl, gram),
        "framed_action_square_two": np.array_equal(framed @ framed, 2 * identity),
        "framed_action_self_adjoint": np.array_equal(framed.T @ gram, gram @ framed),
        "witness_is_positive_eigenvector": all(
            sum(
                (Silver.coerce(int(framed[row, column])) * witness[column] for column in range(20)),
                SILVER_ZERO,
            )
            == SQRT2 * witness[row]
            for row in range(20)
        ),
        "witness_is_timelike": witness_square.sign() > 0,
        "ridge_zero_set_is_2_4_5": tuple(
            index for index, value in enumerate(pairings) if not value
        )
        == (2, 4, 5),
        "all_other_wall_pairings_strict": all(
            value.sign() > 0
            for index, value in enumerate(pairings)
            if index not in (2, 4, 5)
        ),
        "walls_2_and_5_have_same_projected_mirror": (
            _projected_pair(roots[:, 2], roots[:, 2], gram, framed)
            * _projected_pair(roots[:, 5], roots[:, 5], gram, framed)
            == _projected_pair(roots[:, 2], roots[:, 5], gram, framed) ** 2
        ),
        "ridge_normals_are_independent": (
            projected[0][0] * projected[1][1] - projected[0][1] ** 2
        ).sign()
        > 0,
        "angle_ratio_strictly_between_m3_and_m4": (
            lower_gap.sign() > 0 and upper_gap.sign() > 0
        ),
    }
    return {
        "checks": checks,
        "passes": all(checks.values()),
        "witness": witness,
        "pairings": pairings,
        "witness_square": witness_square,
        "projected_gram_scaled_by_four": projected,
        "angle_ratio": angle_ratio,
        "lower_gap_from_cos2_pi_over_3": lower_gap,
        "upper_gap_to_cos2_pi_over_4": upper_gap,
    }


def build(source: Path) -> dict[str, object]:
    if _file_sha256(source) != ALEXEEV_SOURCE_SHA256:
        raise ValueError("Alexeev ancillary source digest mismatch")
    chamber = _alexeev_matrix(source)
    gram, unframed = silver_action_in_alexeev_basis(chamber)
    roots = _wall_coordinates(chamber)
    squares = [int(chamber[index, index]) for index in range(chamber.rows)]
    framed, word_array, _ = weyl_frame(
        gram,
        unframed,
        roots,
        squares,
        max_steps=1_000,
    )
    word = tuple(int(value) for value in word_array)
    replayed, _ = _framed_action(gram, roots, unframed, word)
    if not np.array_equal(framed, replayed):
        raise AssertionError("the stored framing word does not reproduce")
    derived = _derived_checks(gram, roots, unframed, framed, word)
    if not derived["passes"]:
        raise AssertionError("the silver folding obstruction failed")

    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "silver_alexeev_h9_folding_obstruction",
        "claim": (
            "The naturally framed silver eigenslice of Alexeev's (20,2,1) "
            "50-wall chamber is not Coxeter: genuine ridge walls 2 and 4 "
            "have cos(theta)^2=1/4+sqrt(2)/8, strictly between the m=3 and "
            "m=4 Coxeter values."
        ),
        "source": {
            "description": "Alexeev arXiv:2209.09110v4 ancillary all_cases_in_2el_paper.txt",
            "sha256": ALEXEEV_SOURCE_SHA256,
        },
        "basis_gram": _integer_matrix(gram),
        "wall_coordinates": _integer_matrix(roots),
        "unframed_silver_action": _integer_matrix(unframed),
        "framing_word_zero_based_wall_indices": list(word),
        "framed_silver_action": _integer_matrix(framed),
        "plus_basis_columns": list(PLUS_BASIS_COLUMNS),
        "ridge_walls": list(RIDGE_WALLS),
        "duplicate_projected_wall": DUPLICATE_WALL,
        "ridge_coefficients": [_silver(value) for value in RIDGE_COEFFICIENTS],
        "ridge_witness": [_silver(value) for value in derived["witness"]],
        "ridge_wall_pairings": [_silver(value) for value in derived["pairings"]],
        "ridge_witness_square": _silver(derived["witness_square"]),
        "projected_ridge_gram_scaled_by_four": _silver_matrix(
            derived["projected_gram_scaled_by_four"]
        ),
        "cosine_squared": _silver(derived["angle_ratio"]),
        "gap_above_cosine_squared_pi_over_3": _silver(
            derived["lower_gap_from_cos2_pi_over_3"]
        ),
        "gap_below_cosine_squared_pi_over_4": _silver(
            derived["upper_gap_to_cos2_pi_over_4"]
        ),
        "checks": derived["checks"],
        "scope": {
            "proved": (
                "This explicit Weyl framing of the Alexeev chamber has a "
                "non-Coxeter H9 eigenslice."
            ),
            "not_proved": [
                "Every Weyl framing or every silver chamber is non-Coxeter.",
                "The norm-two silver reflection group is nonreflective.",
                "Finite or infinite volume for a different silver root chamber.",
            ],
        },
    }
    manifest["certificate_sha256"] = _certificate_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    gram = _decode_integer_matrix(manifest.get("basis_gram"))
    roots = _decode_integer_matrix(manifest.get("wall_coordinates"))
    unframed = _decode_integer_matrix(manifest.get("unframed_silver_action"))
    framed = _decode_integer_matrix(manifest.get("framed_silver_action"))
    word_value = manifest.get("framing_word_zero_based_wall_indices")
    if not isinstance(word_value, list) or any(not isinstance(value, int) for value in word_value):
        raise ValueError("the framing word is malformed")
    word = tuple(word_value)
    if gram.shape != (20, 20) or roots.shape != (20, 50):
        raise ValueError("the frozen lattice or wall matrix has the wrong shape")
    replayed, _ = _framed_action(gram, roots, unframed, word)
    if not np.array_equal(replayed, framed):
        raise ValueError("the frozen framing word does not reproduce")
    derived = _derived_checks(gram, roots, unframed, framed, word)
    if not derived["passes"]:
        raise ValueError("the exact silver folding checks failed")
    expected_fields = {
        "ridge_witness": [_silver(value) for value in derived["witness"]],
        "ridge_wall_pairings": [_silver(value) for value in derived["pairings"]],
        "ridge_witness_square": _silver(derived["witness_square"]),
        "projected_ridge_gram_scaled_by_four": _silver_matrix(
            derived["projected_gram_scaled_by_four"]
        ),
        "cosine_squared": _silver(derived["angle_ratio"]),
        "gap_above_cosine_squared_pi_over_3": _silver(
            derived["lower_gap_from_cos2_pi_over_3"]
        ),
        "gap_below_cosine_squared_pi_over_4": _silver(
            derived["upper_gap_to_cos2_pi_over_4"]
        ),
        "checks": derived["checks"],
    }
    for key, expected in expected_fields.items():
        if manifest.get(key) != expected:
            raise ValueError(f"the frozen field {key} does not reproduce")
    digest = _certificate_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("silver folding certificate digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("silver folding certificate is not frozen")
    return {
        "framing_steps": len(word),
        "ridge_zero_walls": [2, 4, 5],
        "cosine_squared": expected_fields["cosine_squared"],
        "certificate_sha256": digest,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)
    if args.source is not None:
        manifest = build(args.source)
        if args.write:
            OUTPUT.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        print(json.dumps({"certificate_sha256": manifest["certificate_sha256"]}, indent=2))
        return
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    print(json.dumps(verify(manifest), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
