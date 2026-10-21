"""Exhaustively search length-four Shimada words centralizing ``tau``.

The published row isometries are ``M_i`` and their column actions are
``C_i = M_i.T``.  Split a length-four column word as

    (C_l C_k) (C_j C_i) = B A.

It centralizes ``tau`` if and only if

    A tau A^-1 = B^-1 tau B.                         (1)

There are only ``767^2`` possibilities on each side of (1).  We index the
right-hand list by a modular probe fingerprint.  This fingerprint is only a
necessary filter: every candidate quadruple is reconstructed over ``ZZ``
using python-flint and its full ``20 x 20`` commutator with ``tau`` is checked
entrywise.  Hence modular or digest collisions can add exact checks but can
neither hide a solution nor create one.

The optional python-flint dependency is imported lazily so that merely
importing this module does not impose it on the rest of the project.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import struct
import sys
import time
from typing import Any

from sympy import Matrix

from scripts.analyze_shimada_disc11 import read_shimada_data
from scripts.search_shimada_disc11_centralizer_words import (
    _action_record,
)
from scripts.shimada_disc11_conjugate import (
    DISC11_SOURCE_SHA256,
    conjugated_phi_action_shimada,
)


DEFAULT_PRIOR_CERTIFICATE = Path(
    "data/candidates/shimada_disc11_centralizer_words3.json"
)
FINGERPRINT_PRIME = 1_000_003
FINGERPRINT_PROBE_COLUMNS = 4


def _require_flint() -> tuple[Any, Any]:
    try:
        from flint import fmpz_mat
    except ImportError as exc:  # pragma: no cover - depends on local tooling
        raise RuntimeError(
            "the exhaustive length-four search requires optional python-flint"
        ) from exc
    return fmpz_mat, struct.Struct(f"<{20 * FINGERPRINT_PROBE_COLUMNS}I")


def _fmpz_matrix(fmpz_mat: Any, rows: Any) -> Any:
    return fmpz_mat([[int(entry) for entry in row] for row in rows])


def _integral_inverse(matrix: Any) -> Any:
    numerator, denominator = matrix.inv().numer_denom()
    if denominator != 1:
        raise ValueError("a published automorphism matrix has nonintegral inverse")
    return numerator


def _probe_matrix(fmpz_mat: Any) -> Any:
    """The first two columns agree with the exhaustive length-three probe."""

    return fmpz_mat(
        [
            [
                row + 1,
                (row + 1) * (row + 3) + 7,
                (row + 2) ** 3 + 11,
                1 if row in (0, 7, 13, 19) else 0,
            ]
            for row in range(20)
        ]
    )


def _fingerprint(
    matrix: Any,
    probe: Any,
    packer: struct.Struct,
    *,
    prime: int = FINGERPRINT_PRIME,
) -> bytes:
    """Return a necessary-equality fingerprint, never a proof of equality."""

    image = matrix * probe
    payload = packer.pack(*(int(entry) % prime for entry in image.entries()))
    return hashlib.blake2b(
        payload,
        digest_size=16,
        person=b"disc11-word4",
    ).digest()


def _add_bucket(buckets: dict[bytes, int | list[int]], key: bytes, value: int) -> None:
    previous = buckets.get(key)
    if previous is None:
        buckets[key] = value
    elif isinstance(previous, int):
        buckets[key] = [previous, value]
    else:
        previous.append(value)


def _bucket_values(bucket: int | list[int]) -> tuple[int, ...] | list[int]:
    return (bucket,) if isinstance(bucket, int) else bucket


def _stream_json_word_digest_update(
    digest: Any,
    word: tuple[int, int, int, int],
    *,
    first_word: bool,
) -> None:
    if not first_word:
        digest.update(b",")
    digest.update(json.dumps(word, separators=(",", ":")).encode())


def _sympy_matrix(matrix: Any) -> Matrix:
    return Matrix(matrix.nrows(), matrix.ncols(), [int(x) for x in matrix.entries()])


def _load_prior_certificate(path: Path) -> tuple[dict[str, object], str]:
    encoded = path.read_bytes()
    certificate = json.loads(encoded)
    if certificate.get("verified") is not True:
        raise ValueError("the length-three certificate is not marked verified")
    if certificate.get("source_sha256") != DISC11_SOURCE_SHA256:
        raise ValueError("the length-three certificate uses a different source")
    if certificate.get("published_generator_count") != 767:
        raise ValueError("the length-three certificate has an unexpected alphabet")
    return certificate, hashlib.sha256(encoded).hexdigest()


def search_length_four(
    disc11_path: Path,
    *,
    prior_certificate_path: Path = DEFAULT_PRIOR_CERTIFICATE,
    progress: bool = False,
) -> dict[str, object]:
    """Exhaust all ordered words of exact published-generator length four."""

    started = time.perf_counter()
    source_digest = hashlib.sha256(disc11_path.read_bytes()).hexdigest()
    if source_digest != DISC11_SOURCE_SHA256:
        raise ValueError("unexpected disc11.txt checksum")
    prior, prior_digest = _load_prior_certificate(prior_certificate_path)

    fmpz_mat, packer = _require_flint()
    _, _, generator_rows = read_shimada_data(disc11_path)
    rank = 20
    count = len(generator_rows)
    if count != 767:
        raise ValueError("unexpected published generator count")

    columns = tuple(
        fmpz_mat(
            [
                [int(generator_rows[index][column][row]) for column in range(rank)]
                for row in range(rank)
            ]
        )
        for index in range(count)
    )
    inverses = tuple(_integral_inverse(column) for column in columns)
    identity = fmpz_mat(rank, rank)
    for diagonal in range(rank):
        identity[diagonal, diagonal] = 1
    for column, inverse in zip(columns, inverses, strict=True):
        if column * inverse != identity or inverse * column != identity:
            raise AssertionError("an exact generator inverse check failed")

    tau = _fmpz_matrix(fmpz_mat, conjugated_phi_action_shimada())
    # Exact column action of the independently known row factorization
    # g_507 g_1 g_39.  Comparing matrices (not merely their digests) makes
    # the final identification of the sole nonidentity action unconditional.
    known_return = columns[39] * columns[1] * columns[507]
    if known_return * tau != tau * known_return:
        raise AssertionError("the known D68/D27 return no longer centralizes tau")
    probe = _probe_matrix(fmpz_mat)
    left_single = tuple(
        column * tau * inverse
        for column, inverse in zip(columns, inverses, strict=True)
    )
    right_single = tuple(
        inverse * tau * column
        for column, inverse in zip(columns, inverses, strict=True)
    )

    # For B=C_l C_k, B^-1 tau B = C_k^-1 (C_l^-1 tau C_l) C_k.
    right_buckets: dict[bytes, int | list[int]] = {}
    for third in range(count):
        inverse_third = inverses[third]
        column_third = columns[third]
        for fourth in range(count):
            conjugate = inverse_third * right_single[fourth] * column_third
            key = _fingerprint(conjugate, probe, packer)
            _add_bucket(right_buckets, key, third * count + fourth)
        if progress and (third + 1) % 100 == 0:
            print(
                f"right fingerprints: {third + 1}/{count}",
                file=sys.stderr,
                flush=True,
            )

    right_collision_buckets = sum(
        isinstance(bucket, list) for bucket in right_buckets.values()
    )
    largest_right_bucket = max(
        len(bucket) if isinstance(bucket, list) else 1
        for bucket in right_buckets.values()
    )

    # Accepted words are visited lexicographically.  Hash the canonical compact
    # JSON list incrementally, avoiding a potentially large in-memory word list.
    word_digest = hashlib.sha256()
    word_digest.update(b"[")
    accepted_words = 0
    modular_candidates = 0
    distinct_actions: dict[tuple[int, ...], tuple[Any, int]] = {}
    identity_words = 0

    for first in range(count):
        for second in range(count):
            left_conjugate = columns[second] * left_single[first] * inverses[second]
            bucket = right_buckets.get(_fingerprint(left_conjugate, probe, packer))
            if bucket is None:
                continue
            candidates = _bucket_values(bucket)
            modular_candidates += len(candidates)
            first_half = columns[second] * columns[first]
            for pair_index in candidates:
                third, fourth = divmod(pair_index, count)
                second_half = columns[fourth] * columns[third]
                word_action = second_half * first_half
                if word_action * tau != tau * word_action:
                    continue

                word = (first, second, third, fourth)
                _stream_json_word_digest_update(
                    word_digest,
                    word,
                    first_word=accepted_words == 0,
                )
                accepted_words += 1
                if word_action == identity:
                    identity_words += 1
                    continue

                action_key = tuple(int(entry) for entry in word_action.entries())
                previous = distinct_actions.get(action_key)
                if previous is None:
                    distinct_actions[action_key] = (word_action, 1)
                else:
                    distinct_actions[action_key] = (previous[0], previous[1] + 1)
        if progress and (first + 1) % 100 == 0:
            print(
                f"left fingerprints: {first + 1}/{count}; "
                f"modular candidates={modular_candidates}; accepted={accepted_words}",
                file=sys.stderr,
                flush=True,
            )

    word_digest.update(b"]")
    nonidentity_records: list[dict[str, object]] = []
    for action_key in sorted(distinct_actions):
        action, multiplicity = distinct_actions[action_key]
        record = _action_record(_sympy_matrix(action))
        record["ordered_word_multiplicity"] = multiplicity
        record["is_known_D68_D27_involution"] = action == known_return
        nonidentity_records.append(record)

    prior_distinct = prior["distinct_length_three_actions"]
    if not isinstance(prior_distinct, list):
        raise ValueError("malformed length-three distinct-action record")

    result: dict[str, object] = {
        "verified": True,
        "source_sha256": source_digest,
        "published_generator_count": count,
        "prior_length_three_certificate": {
            "path": str(prior_certificate_path),
            "file_sha256": prior_digest,
            "centralizing_ordered_words": prior[
                "centralizing_ordered_length_three_words"
            ],
            "distinct_centralizing_actions": len(prior_distinct),
        },
        "exact_length_four": {
            "ordered_words_exhausted": count**4,
            "half_words_per_side": count**2,
            "split_equation": "A*tau*A^-1 = B^-1*tau*B",
            "fingerprint": {
                "role": "necessary_filter_only",
                "prime": FINGERPRINT_PRIME,
                "probe_columns": FINGERPRINT_PROBE_COLUMNS,
                "digest": "BLAKE2b-128",
                "right_buckets": len(right_buckets),
                "right_collision_buckets": right_collision_buckets,
                "largest_right_bucket": largest_right_bucket,
            },
            "modular_candidate_quadruples": modular_candidates,
            "exact_integral_commutator_checks": modular_candidates,
            "false_modular_candidates": modular_candidates - accepted_words,
            "centralizing_ordered_words": accepted_words,
            "word_quadruples_sha256": word_digest.hexdigest(),
            "word_digest_encoding": (
                "canonical compact JSON array in lexicographic word order"
            ),
            "identity_word_multiplicity": identity_words,
            "nonidentity_word_multiplicity": accepted_words - identity_words,
            "distinct_centralizing_actions": (
                (1 if identity_words else 0) + len(distinct_actions)
            ),
            "distinct_actions": [
                {
                    **_action_record(Matrix.eye(rank)),
                    "ordered_word_multiplicity": identity_words,
                },
                *nonidentity_records,
            ],
        },
        "elapsed_seconds": time.perf_counter() - started,
        "scope": (
            "exact_words_of_published_generator_length_four_only; "
            "combined_with_the_SHA_pinned_prior_certificate_through_length_four"
        ),
    }
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("disc11", type=Path)
    parser.add_argument(
        "--prior-certificate",
        type=Path,
        default=DEFAULT_PRIOR_CERTIFICATE,
    )
    parser.add_argument("--progress", action="store_true")
    args = parser.parse_args()
    print(
        json.dumps(
            search_length_four(
                args.disc11,
                prior_certificate_path=args.prior_certificate,
                progress=args.progress,
            ),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
