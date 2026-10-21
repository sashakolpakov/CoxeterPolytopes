#!/usr/bin/env python3
"""Certify a four-wall orthogonal doubling cube in the JK-H12 chamber."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
import math
from pathlib import Path
from typing import Any, Sequence

import sympy as sp

import scripts.certify_jacobian_kummer_h12_single_doublings as base


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/jacobian_kummer_h12_doubling_cube.json"
SINGLE_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h12_single_doublings.json"
SINGLE_BUILDER_SHA256 = (
    "0bd5672b54e091da680192be8873fa7397a4e9c162c592dff6cf1ebd6e15386b"
)
SINGLE_ARCHIVE_SHA256 = (
    "01eb89baef2cedcf78d3a20dd1c55931db9a9a2ca83bc5f4b7c1a1329cbdca98"
)
SCHEMA = "coxeter.jacobian_kummer.h12_doubling_cube.v1"

DOUBLING_WALLS = (30, 57, 58, 67)  # zero-based: 31,58,59,68
EXPECTED_COUNTS = {
    (30,): 296,
    (57,): 296,
    (58,): 296,
    (67,): 317,
    (30, 57): 523,
    (30, 58): 523,
    (30, 67): 559,
    (57, 58): 523,
    (57, 67): 559,
    (58, 67): 559,
    (30, 57, 58): 941,
    (30, 57, 67): 1004,
    (30, 58, 67): 1004,
    (57, 58, 67): 1004,
    (30, 57, 58, 67): 1827,
}
TYPE_REPRESENTATIVES = {
    296: (30,),
    317: (67,),
    523: (30, 57),
    559: (30, 67),
    941: (30, 57, 58),
    1004: (30, 57, 67),
    1827: (30, 57, 58, 67),
}


def _file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _root_digest(roots: Sequence[Sequence[int]]) -> str:
    return sha256(
        json.dumps(roots, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _apply_mask(
    gram: base.Matrix,
    seed_roots: Sequence[base.IntegralVector],
    selected: Sequence[int],
    mask: int,
    vector: Sequence[int | Fraction],
) -> base.Vector:
    image = tuple(base._fraction(value) for value in vector)
    for position, wall in enumerate(selected):
        if mask & (1 << position):
            image = base._reflection_apply(gram, seed_roots[wall], image)
    return image


def _construct(
    gram: base.Matrix,
    seed_roots: Sequence[base.IntegralVector],
    seed_interior: base.IntegralVector,
    seed_witnesses: Sequence[base.Vector],
    selected: tuple[int, ...],
    validate_witnesses: bool,
) -> dict[str, Any]:
    selected_set = set(selected)
    copy_count = 1 << len(selected)
    transformed_interiors = tuple(
        _apply_mask(gram, seed_roots, selected, mask, seed_interior)
        for mask in range(copy_count)
    )
    interior = tuple(
        sum((image[index] for image in transformed_interiors), Fraction(0))
        for index in range(len(seed_interior))
    )

    boundary: dict[base.IntegralVector, dict[str, Any]] = {}
    for mask in range(copy_count):
        for seed_index, seed_root in enumerate(seed_roots):
            if seed_index in selected_set:
                continue
            root = base._primitive(
                _apply_mask(gram, seed_roots, selected, mask, seed_root)
            )
            witness = _apply_mask(
                gram, seed_roots, selected, mask, seed_witnesses[seed_index]
            )
            key = base._unoriented_key(root)
            pairing = base._pair(gram, interior, root)
            if pairing == 0:
                raise ArithmeticError("cube interior lies on a boundary wall")
            if pairing < 0:
                root = tuple(-value for value in root)
            boundary.setdefault(
                key,
                {
                    "root": root,
                    "witness": witness,
                    "mask": mask,
                    "seed_facet": seed_index + 1,
                },
            )

    records = tuple(boundary[key] for key in sorted(boundary))
    roots = tuple(record["root"] for record in records)
    covectors = tuple(base._covector(gram, root) for root in roots)
    if not (
        sp.Matrix(roots).rank() == 13
        and base._pair(gram, interior, interior) > 0
        and all(base._pair(gram, interior, root) > 0 for root in roots)
        and all(base._content(root) == 1 for root in roots)
        and all(base._reflection_integral(gram, root) for root in roots)
    ):
        raise ArithmeticError("a cube boundary failed rank/orientation/integrality")

    witness_payload: dict[str, Any]
    if validate_witnesses:
        minimum_square: Fraction | None = None
        minimum_slack: Fraction | None = None
        for target, record in enumerate(records):
            witness = record["witness"]
            square = base._pair(gram, witness, witness)
            if not (
                square > 0
                and base._pair_from_covector(witness, covectors[target]) == 0
            ):
                raise ArithmeticError("a cube witness failed its target")
            minimum_square = square if minimum_square is None else min(minimum_square, square)
            for index, covector in enumerate(covectors):
                if index == target:
                    continue
                slack = base._pair_from_covector(witness, covector)
                if slack <= 0:
                    raise ArithmeticError("a cube boundary inequality is redundant")
                minimum_slack = slack if minimum_slack is None else min(minimum_slack, slack)
        witness_payload = {
            "validated_for_every_facet": True,
            "count": len(records),
            "minimum_square": base._payload_number(minimum_square),
            "minimum_off_target_slack": base._payload_number(minimum_slack),
        }
    else:
        witness_payload = {
            "validated_for_every_facet": False,
            "reason": "a separate representative with this facet count receives the exhaustive witness replay",
        }

    squares = tuple(
        base._pair_from_covector(root, covector)
        for root, covector in zip(roots, covectors, strict=True)
    )
    divisibilities = tuple(
        math.gcd(*(abs(int(value)) for value in covector if value))
        for covector in covectors
    )
    root_payload = [list(root) for root in roots]
    return {
        "selected_walls_one_based": [wall + 1 for wall in selected],
        "selected_wall_count": len(selected),
        "copy_count_and_volume_multiple": copy_count,
        "facet_count": len(roots),
        "boundary_root_digest": _root_digest(root_payload),
        "full_linear_rank": 13,
        "primitive_root_square_histogram": {
            str(key): value for key, value in sorted(Counter(squares).items())
        },
        "primitive_root_divisibility_histogram": {
            str(key): value for key, value in sorted(Counter(divisibilities).items())
        },
        "strict_interior_square": base._payload_number(
            base._pair(gram, interior, interior)
        ),
        "strict_facet_witnesses": witness_payload,
        "finite_volume": True,
        "compact": False,
        "all_exact_checks_pass": True,
    }


def build() -> dict[str, Any]:
    if _file_sha256(Path(base.__file__)) != SINGLE_BUILDER_SHA256:
        raise ArithmeticError("the pinned single-doubling builder changed")
    if _file_sha256(SINGLE_ARCHIVE) != SINGLE_ARCHIVE_SHA256:
        raise ArithmeticError("the pinned single-doubling archive changed")
    single = json.loads(SINGLE_ARCHIVE.read_text(encoding="utf-8"))
    if not single["all_checks_pass"]:
        raise ArithmeticError("the single-doubling theorem is not green")

    source, gram, roots, interior = base._load_source()
    relations = base._seed_relations(gram, roots)
    doubling = set(base._doubling_walls(relations))
    if not (
        set(DOUBLING_WALLS) <= doubling
        and all(
            base._pair(gram, roots[left], roots[right]) == 0
            for position, left in enumerate(DOUBLING_WALLS)
            for right in DOUBLING_WALLS[position + 1 :]
        )
    ):
        raise ArithmeticError("the claimed C2^4 doubling cube changed")

    seed_witnesses = tuple(
        base._seed_facet_witness(gram, interior, root) for root in roots
    )
    constructions = []
    for size in range(1, len(DOUBLING_WALLS) + 1):
        for selected in combinations(DOUBLING_WALLS, size):
            expected = EXPECTED_COUNTS[selected]
            validate = TYPE_REPRESENTATIVES.get(expected) == selected
            row = _construct(
                gram,
                roots,
                interior,
                seed_witnesses,
                selected,
                validate,
            )
            if row["facet_count"] != expected:
                raise ArithmeticError("a cube facet count changed")
            constructions.append(row)

    facet_histogram = Counter(row["facet_count"] for row in constructions)
    expected_histogram = Counter({296: 3, 317: 1, 523: 3, 559: 3, 941: 1, 1004: 3, 1827: 1})
    if facet_histogram != expected_histogram:
        raise ArithmeticError("the cube facet-count census changed")
    representative_rows = [
        row
        for row in constructions
        if TYPE_REPRESENTATIVES[row["facet_count"]]
        == tuple(wall - 1 for wall in row["selected_walls_one_based"])
    ]
    if not all(
        row["strict_facet_witnesses"]["validated_for_every_facet"]
        for row in representative_rows
    ):
        raise ArithmeticError("a diagram-type representative lacks witnesses")

    return {
        "schema": SCHEMA,
        "status_date": "2026-08-08",
        "all_checks_pass": True,
        "pinned_inputs": {
            "single_doubling_builder_sha256": SINGLE_BUILDER_SHA256,
            "single_doubling_archive_sha256": SINGLE_ARCHIVE_SHA256,
            "seed_id": source["chamber"]["id"],
        },
        "orthogonal_doubling_cube": {
            "walls_one_based": [wall + 1 for wall in DOUBLING_WALLS],
            "group": "C2^4",
            "faithful": True,
            "complete_nonempty_subset_count": len(constructions),
            "facet_count_histogram": {
                str(key): value for key, value in sorted(facet_histogram.items())
            },
            "pairwise_nonisomorphic_type_lower_bound": len(facet_histogram),
            "reason_for_lower_bound": "seven distinct facet counts",
        },
        "constructions": constructions,
        "decision": {
            "finite_volume_noncompact_H12_polytopes": True,
            "new_to_project_types_beyond_single_doublings": 5,
            "total_cube_type_lower_bound": 7,
            "largest_certified_descendant": "H12/1827",
            "volume_multipliers": [2, 4, 8, 16],
            "same_weak_full_Picard_K3_carrier": True,
            "worldwide_priority": "open",
        },
        "scope": {
            "all_15_marked_subsets_constructed": True,
            "one_representative_per_facet_count_has_all_facet_witnesses": True,
            "same_facet_count_labelled_isomorphism_not_claimed": True,
            "Coxeter_and_volume_proof": "iterated exact even-angle doubling along four commuting orthogonal walls",
            "literal_descendant_Vinberg_replay": False,
            "full_nef_or_surface_automorphism_claimed": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if OUTPUT.read_text(encoding="utf-8") != rendered:
            raise ArithmeticError("JK-H12 doubling-cube archive changed")
    else:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(rendered, encoding="utf-8")
    print(json.dumps(payload["decision"], sort_keys=True))
    print(OUTPUT)


if __name__ == "__main__":
    main()
