"""Search short exact words in Shimada's Aut(X) generators centralizing tau.

For row isometries ``M_i``, the induced column action is ``C_i=M_i^T``.
A two-letter row word ``M_i M_j`` centralizes the golden action exactly when

    C_j C_i tau = tau C_j C_i.

Instead of testing all 767^2 products, hash the two conjugacy lists

    C_i tau C_i^-1  and  C_j^-1 tau C_j.

Their intersections are precisely the centralizing two-letter words.  All
arithmetic is integral/rational and every reported product is rechecked
entrywise.  Absence is only a word-length-two statement, not a centralizer
classification.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from sympy import Matrix

from scripts.analyze_shimada_disc11 import read_shimada_data
from scripts.analyze_shimada_disc11_enhanced_quotient import (
    D68_RETURN_PATH,
    _walk,
)
from scripts.shimada_disc11_atlas import ShimadaDisc11Atlas
from scripts.shimada_disc11_conjugate import (
    DISC11_SOURCE_SHA256,
    conjugated_phi_action_shimada,
)


def _key(matrix: Matrix) -> tuple[int, ...]:
    if any(not entry.is_Integer for entry in matrix):
        raise ValueError("a conjugated golden action became nonintegral")
    return tuple(int(entry) for entry in matrix)


def _encoded_matrix(matrix: Matrix) -> list[list[int]]:
    if any(not entry.is_Integer for entry in matrix):
        raise ValueError("a centralizing word became nonintegral")
    return [[int(entry) for entry in row] for row in matrix.tolist()]


def _canonical_sha(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def _infinite_order_trace_power(
    matrix: Matrix, *, max_power: int = 40
) -> tuple[int | None, int | None]:
    """Prove infinite order if one power has absolute trace above rank."""

    power = Matrix.eye(matrix.rows)
    for exponent in range(1, max_power + 1):
        power = power * matrix
        trace = int(power.trace())
        if abs(trace) > matrix.rows:
            return exponent, trace
    return None, None


def _action_record(matrix: Matrix) -> dict[str, object]:
    trace_power, trace_value = _infinite_order_trace_power(matrix)
    return {
        "column_action_sha256": _canonical_sha(_encoded_matrix(matrix)),
        "is_identity": matrix == Matrix.eye(matrix.rows),
        "square_is_identity": matrix * matrix == Matrix.eye(matrix.rows),
        "trace": int(matrix.trace()),
        "characteristic_polynomial": str(matrix.charpoly().as_expr().factor()),
        "infinite_order_trace_power": trace_power,
        "infinite_order_trace_value": trace_value,
    }


def search_length_two(
    path: Path,
    *,
    include_match_pairs: bool = False,
    max_word_length: int = 3,
) -> dict[str, object]:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != DISC11_SOURCE_SHA256:
        raise ValueError("unexpected disc11.txt checksum")
    _, _, generators = read_shimada_data(path)
    columns = tuple(Matrix(generator).T for generator in generators)
    inverses = tuple(column.inv() for column in columns)
    tau = Matrix(conjugated_phi_action_shimada())
    identity = Matrix.eye(20)

    left: dict[tuple[int, ...], list[int]] = {}
    right: dict[tuple[int, ...], list[int]] = {}
    left_matrices: list[Matrix] = []
    right_matrices: list[Matrix] = []
    for index, (column, inverse) in enumerate(zip(columns, inverses, strict=True)):
        left_matrix = column * tau * inverse
        right_matrix = inverse * tau * column
        left_matrices.append(left_matrix)
        right_matrices.append(right_matrix)
        left.setdefault(_key(left_matrix), []).append(index)
        right.setdefault(_key(right_matrix), []).append(index)

    match_pairs: list[tuple[int, int]] = []
    distinct: dict[tuple[int, ...], Matrix] = {}
    for conjugate in sorted(set(left) & set(right)):
        for first in left[conjugate]:
            for second in right[conjugate]:
                # The row word is M_first M_second, hence its column action
                # is C_second C_first.
                column_word = columns[second] * columns[first]
                if column_word * tau != tau * column_word:
                    raise AssertionError("conjugacy hash produced a false match")
                word_key = _key(column_word)
                distinct[word_key] = column_word
                match_pairs.append((first, second))

    distinct_records = [_action_record(distinct[key]) for key in sorted(distinct)]

    length_three_words: list[tuple[int, int, int]] = []
    length_three_actions: dict[tuple[int, ...], Matrix] = {}
    if max_word_length >= 3:
        # A modular probe is a necessary filter only.  Every surviving word
        # is recomputed and accepted over the exact integers below, so modular
        # collisions cannot create a false result or hide a genuine word.
        prime = 1_000_003
        probes = np.array(
            [[i + 1, (i + 1) * (i + 3) + 7] for i in range(20)],
            dtype=np.int64,
        )

        def modular_array(matrix: Matrix) -> np.ndarray:
            return np.array(
                [[int(entry) % prime for entry in row] for row in matrix.tolist()],
                dtype=np.int64,
            )

        left_mod = np.stack([modular_array(matrix) for matrix in left_matrices])
        right_mod = np.stack([modular_array(matrix) for matrix in right_matrices])
        columns_mod = np.stack([modular_array(matrix) for matrix in columns])
        inverses_mod = np.stack([modular_array(matrix) for matrix in inverses])
        target_images = np.matmul(right_mod, probes) % prime
        targets: dict[bytes, list[int]] = {}
        for third in range(len(columns)):
            targets.setdefault(target_images[third].tobytes(), []).append(third)

        for second in range(len(columns)):
            pulled_probes = inverses_mod[second] @ probes % prime
            first_images = np.matmul(left_mod, pulled_probes) % prime
            pair_images = (
                np.einsum(
                    "ab,ibk->iak",
                    columns_mod[second],
                    first_images,
                    optimize=True,
                )
                % prime
            )
            for first in range(len(columns)):
                for third in targets.get(pair_images[first].tobytes(), ()):
                    column_word = columns[third] * columns[second] * columns[first]
                    if column_word * tau != tau * column_word:
                        continue
                    length_three_words.append((first, second, third))
                    length_three_actions[_key(column_word)] = column_word

    # Express the independently discovered D68 return directly in Shimada's
    # published generating set.  Along this gallery the only nonidentity
    # local factors are Gamma[39], Gamma[1], Gamma[507], accumulated on the
    # left; hence the returned row action is Gamma[507] Gamma[1] Gamma[39].
    atlas = ShimadaDisc11Atlas(path)
    returned, _ = _walk(atlas, atlas.representative_state(68), D68_RETURN_PATH)
    return_row = Matrix(returned.row_isometry)
    generator_rows = tuple(Matrix(generator) for generator in generators)
    length_three_indices = (507, 1, 39)
    length_three_row = (
        generator_rows[length_three_indices[0]]
        * generator_rows[length_three_indices[1]]
        * generator_rows[length_three_indices[2]]
    )
    if length_three_row != return_row:
        raise AssertionError("the archived D68 return factorization changed")
    length_three_column = length_three_row.T
    if length_three_column * tau != tau * length_three_column:
        raise AssertionError("the archived D68 return no longer centralizes tau")
    if max_word_length >= 3 and _key(length_three_column) not in length_three_actions:
        raise AssertionError("the exhaustive length-three search missed the D68 return")

    result = {
        "verified": True,
        "source_sha256": digest,
        "published_generator_count": len(columns),
        "individual_centralizing_generators": sum(
            column * tau == tau * column for column in columns
        ),
        "left_conjugacy_classes": len(left),
        "right_conjugacy_classes": len(right),
        "matching_conjugacy_classes": len(set(left) & set(right)),
        "centralizing_ordered_length_two_words": len(match_pairs),
        "distinct_centralizing_column_actions": len(distinct),
        "nonidentity_centralizing_length_two_actions": sum(
            matrix != identity for matrix in distinct.values()
        ),
        "length_two_match_pairs_sha256": _canonical_sha(match_pairs),
        "distinct_actions": distinct_records,
        "exhaustive_length_three": {
            "performed": max_word_length >= 3,
            "modular_filter_prime": 1_000_003 if max_word_length >= 3 else None,
            "probe_columns": 2 if max_word_length >= 3 else None,
            "centralizing_ordered_words": len(length_three_words),
            "word_triples_sha256": (
                _canonical_sha(length_three_words) if max_word_length >= 3 else None
            ),
            "distinct_centralizing_actions": len(length_three_actions),
            "distinct_actions": [
                _action_record(length_three_actions[key])
                for key in sorted(length_three_actions)
            ],
        },
        "known_nontrivial_length_three_return": {
            "row_generator_indices_zero_based": list(length_three_indices),
            "column_action_sha256": _canonical_sha(
                _encoded_matrix(length_three_column)
            ),
            "square_is_identity": length_three_column * length_three_column
            == identity,
            "trace": int(length_three_column.trace()),
            "characteristic_polynomial": str(
                length_three_column.charpoly().as_expr().factor()
            ),
        },
        "scope": (
            "exact_words_of_published_generator_length_at_most_three_only"
            if max_word_length >= 3
            else "exact_words_of_published_generator_length_at_most_two_only"
        ),
    }
    if include_match_pairs:
        result["length_two_match_pairs"] = [list(pair) for pair in match_pairs]
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("disc11", type=Path)
    parser.add_argument("--include-match-pairs", action="store_true")
    parser.add_argument("--max-word-length", type=int, choices=(2, 3), default=3)
    args = parser.parse_args()
    print(
        json.dumps(
            search_length_two(
                args.disc11,
                include_match_pairs=args.include_match_pairs,
                max_word_length=args.max_word_length,
            ),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
