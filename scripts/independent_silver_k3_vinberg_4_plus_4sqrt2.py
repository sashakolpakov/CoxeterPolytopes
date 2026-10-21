#!/usr/bin/env python3
"""Independent exhaustive block census for ``x=4+4sqrt(2)``."""

from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
import json
from pathlib import Path
import sys
from typing import Sequence

from coxeter_polytopes.exact.silver import SQRT2, Silver
from coxeter_polytopes.search.silver_k3_vinberg import (
    Vector,
    _enumeration_data,
    _silver_dot,
    _silver_mat_vec,
    _silver_quadratic,
    _weight_alphabet,
    admissible_silver_control_coordinates,
)
from scripts.diagnose_silver_k3_726_global_vinberg import _decode_roots
from scripts.independent_silver_k3_vinberg_4_plus_3sqrt2 import (
    _direct_bounded_block,
    _vector_key,
)
from scripts.independent_silver_k3_vinberg_5_plus_3sqrt2 import (
    _canonical_sha256,
    _jsonable,
)


ROOT = Path(__file__).resolve().parents[1]
PRIOR_ARCHIVE = (
    ROOT / "data/candidates/silver_k3_norm2_vinberg_extension_5_plus_3sqrt2.json"
)
EXPECTED_PRIOR_CERTIFICATE_SHA256 = (
    "1de8fe6d72f76951410950c34ce6b7aee223f12dcdb3e78ca8f481c17c12e189"
)
EXPECTED_REPLAY_SHA256 = (
    "34cd4e858e2bc4181b06bbf20c1d445ec2b08fb83c0a4d0d39585712147ccc42"
)
DEFAULT_OUTPUT = Path("/private/tmp/silver_k3_x4_plus_4sqrt2_block_census.json")


def _phase(message: str) -> None:
    print(f"[independent x=4+4sqrt(2)] {message}", file=sys.stderr, flush=True)


@lru_cache(maxsize=1)
def build_census() -> dict[str, object]:
    prior_document = json.loads(PRIOR_ARCHIVE.read_text(encoding="utf-8"))
    if prior_document.get("certificate_sha256") != EXPECTED_PRIOR_CERTIFICATE_SHA256:
        raise ValueError("the literal 1,076-wall input certificate changed")
    prior_roots = _decode_roots(prior_document)
    if len(prior_roots) != 1076:
        raise ValueError("the literal prior prefix does not have 1,076 roots")
    _phase("loaded and checked the literal 1,076-wall prior prefix")

    gram, finite_inverse, cross, inverse_simple_transpose, weight_form, _ = (
        _enumeration_data()
    )
    coordinate = 4 + 4 * SQRT2
    target = 2 + (1 + SQRT2) * coordinate * coordinate
    alphabet_record = _weight_alphabet(target)
    alphabet = alphabet_record.values
    blocks = ((0, 1, 2, 3), (4, 5, 6, 7), (8,))
    block_of = {
        coordinate_index: block_index
        for block_index, block in enumerate(blocks)
        for coordinate_index in block
    }
    if any(
        weight_form[i][j]
        for i in range(9)
        for j in range(9)
        if block_of[i] != block_of[j]
    ):
        raise AssertionError("the finite weight form lost its orthogonal blocks")

    bounded: list[tuple[tuple[Vector, Silver], ...]] = []
    for block in blocks:
        block_form = tuple(
            tuple(weight_form[i][j] for j in block) for i in block
        )
        bounded.append(_direct_bounded_block(alphabet, block_form, target))
    _phase(
        "component boxes exhausted; bounded counts="
        + repr(tuple(len(vectors) for vectors in bounded))
    )

    middle_by_norm: dict[Silver, list[Vector]] = {}
    for vector, norm in bounded[1]:
        middle_by_norm.setdefault(norm, []).append(vector)
    weight_solutions: list[Vector] = []
    for left, left_norm in bounded[0]:
        for tail, tail_norm in bounded[2]:
            required = target - left_norm - tail_norm
            for middle in middle_by_norm.get(required, ()):
                weight = (*left, *middle, *tail)
                if _silver_quadratic(weight, weight_form) != target:
                    raise AssertionError("a block join has the wrong full norm")
                weight_solutions.append(weight)
    if len(set(weight_solutions)) != len(weight_solutions):
        raise AssertionError("the block join produced duplicate weights")
    _phase(f"norm join complete; exact weights={len(weight_solutions)}")

    reconstructed: list[Vector] = []
    nonintegral = 0
    accepted: list[Vector] = []
    rejected: list[tuple[Vector, tuple[int, ...]]] = []
    for weight in weight_solutions:
        finite_coordinates = _silver_mat_vec(
            finite_inverse,
            tuple(
                entry - cross[index] * coordinate
                for index, entry in enumerate(
                    _silver_mat_vec(inverse_simple_transpose, weight)
                )
            ),
        )
        if any(not entry.is_integral for entry in finite_coordinates):
            nonintegral += 1
            continue
        root = (coordinate, *finite_coordinates)
        if _silver_quadratic(root, gram) != 2:
            raise AssertionError("the block census reconstructed a non-root")
        reconstructed.append(root)
        gram_root = _silver_mat_vec(gram, root)
        blockers = tuple(
            index
            for index, wall in enumerate(prior_roots)
            if _silver_dot(wall, gram_root).sign(0) > 0
        )
        if blockers:
            rejected.append((root, blockers))
        else:
            accepted.append(root)

    accepted_sorted = tuple(sorted(accepted, key=_vector_key))
    rejected_sorted = tuple(
        sorted(rejected, key=lambda item: (_vector_key(item[0]), item[1]))
    )
    reconstructed_sorted = tuple(sorted(reconstructed, key=_vector_key))
    _phase(
        "reconstruction/blockers complete; "
        f"nonintegral={nonintegral}, rejected={len(rejected_sorted)}, "
        f"accepted={len(accepted_sorted)}"
    )

    accepted_images = tuple(_silver_mat_vec(gram, root) for root in accepted_sorted)
    same_height_pairings = tuple(
        _silver_dot(accepted_sorted[i], accepted_images[j])
        for i in range(len(accepted_sorted))
        for j in range(i)
    )
    if any(pairing.sign(0) > 0 for pairing in same_height_pairings):
        raise AssertionError("the accepted same-height cohort is not acute")
    _phase(f"same-height census complete; pairings={len(same_height_pairings)}")

    after = tuple(
        value
        for value in admissible_silver_control_coordinates(14)
        if (value - coordinate).sign(0) > 0
    )
    if not after:
        raise AssertionError("the control-coordinate search box is too small")
    blocker_histogram = Counter(len(blockers) for _root, blockers in rejected_sorted)
    payload: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "silver_k3_x4_plus_4sqrt2_independent_exhaustive_block_replay",
        "prior_certificate_sha256": EXPECTED_PRIOR_CERTIFICATE_SHA256,
        "prior_wall_count": len(prior_roots),
        "control_coordinate": coordinate,
        "height": (1 + SQRT2) * coordinate * coordinate / 2,
        "target": target,
        "alphabet": alphabet,
        "identity_upper": alphabet_record.identity_upper,
        "conjugate_upper": alphabet_record.conjugate_upper,
        "coefficient_bound": alphabet_record.coefficient_bound,
        "block_partition": blocks,
        "component_box_counts": tuple(len(alphabet) ** len(block) for block in blocks),
        "component_bounded_counts": tuple(len(vectors) for vectors in bounded),
        "norm_join_probe_count": len(bounded[0]) * len(bounded[2]),
        "exact_weight_solutions": len(weight_solutions),
        "nonintegral_reconstructions": nonintegral,
        "integral_reconstructions": len(reconstructed_sorted),
        "rejected_root_count": len(rejected_sorted),
        "blocker_set_size_histogram": tuple(sorted(blocker_histogram.items())),
        "accepted_root_count": len(accepted_sorted),
        "accepted_roots": accepted_sorted,
        "same_height_negative_pair_count": sum(
            pairing.sign(0) < 0 for pairing in same_height_pairings
        ),
        "same_height_zero_pair_count": sum(
            pairing.sign(0) == 0 for pairing in same_height_pairings
        ),
        "prefix_wall_count": len(prior_roots) + len(accepted_sorted),
        "next_admissible_control_coordinate": after[0],
        "all_reconstructed_roots_sha256": _canonical_sha256(
            _jsonable(reconstructed_sorted)
        ),
        "rejection_witnesses_sha256": _canonical_sha256(
            _jsonable(rejected_sorted)
        ),
        "accepted_roots_sha256": _canonical_sha256(_jsonable(accepted_sorted)),
        "stopping_status": "undecidable_incomplete_simple_roots",
    }
    payload["census_sha256"] = _canonical_sha256(_jsonable(payload))
    if (
        EXPECTED_REPLAY_SHA256 != "TO_BE_FROZEN"
        and payload["census_sha256"] != EXPECTED_REPLAY_SHA256
    ):
        raise AssertionError("the frozen x=4+4sqrt(2) block replay changed")
    return _jsonable(payload)


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    document = build_census()
    args.output.write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    omitted = {"alphabet", "accepted_roots"}
    print(
        json.dumps(
            {key: value for key, value in document.items() if key not in omitted},
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
