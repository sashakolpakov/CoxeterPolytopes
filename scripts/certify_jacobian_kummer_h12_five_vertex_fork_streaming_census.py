#!/usr/bin/env python3
"""Exact streaming census for invariant shards of the JK--H12 five-fork.

The marked tree has a central edge ``a``, its long-arm continuation ``d``,
and an unordered pair of short-arm colours ``{b,c}``.  Thus ``a != d``,
``a,b,c`` are pairwise distinct, while ``d`` may equal ``b`` or ``c``.

This certifier is deliberately shardable by the number of B-coloured edges.
It performs the exact seed-action quotient, rebuilds every rank-13 boundary,
and classifies the resulting ordinary Coxeter matrices.  SHA-256 is only a
one-way inequality filter: every equal-digest bucket is replayed literally.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
import gzip
import hashlib
import io
from itertools import combinations, product
import json
import math
from pathlib import Path
import resource
import sys
import tempfile
import time
from typing import Any, Iterable, Iterator, Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import scripts.certify_jacobian_kummer_h12_five_vertex_p5_b_le_2_census as p5
import scripts.certify_jacobian_kummer_h12_five_vertex_p5_streaming_shard as p5stream
import scripts.certify_jacobian_kummer_h12_five_vertex_star_census as star
import scripts.certify_jacobian_kummer_h12_four_vertex_subtree_census as four


SCHEMA = "coxeter.jacobian_kummer.h12_five_vertex_fork_streaming_census.v1"
STATUS_DATE = "2026-08-10"
RSS_CEILING_BYTES = 1 << 30

P5_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_p5_b_le_2_census.json"
STAR_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_star_census.json"
DENOMINATOR_ARCHIVE = ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_tree_denominators.json"
FORK_DENOMINATOR_OUTPUT = (
    ROOT
    / "data/enumerations/jacobian_kummer_h12_five_vertex_fork_pattern_denominator.json"
)
P5_ABBB_ARCHIVE = (
    ROOT / "data/enumerations/jacobian_kummer_h12_five_vertex_p5_abbb_census.json"
)
P5_BABB_ROWS = (
    ROOT
    / "data/enumerations/"
    "jacobian_kummer_h12_five_vertex_p5_babb_streaming_census_rows.jsonl.gz"
)
P5_BBBB_ROWS = (
    ROOT
    / "data/enumerations/"
    "jacobian_kummer_h12_five_vertex_p5_bbbb_streaming_census_rows.jsonl.gz"
)
HIGH_B_PRIOR_SHA256 = {
    "P5_ABBB_archive": "b535bbe5d1639358e5a336ce606a4d7eae8dedd0b55835f9ee057321ac394299",
    "P5_BABB_rows": "04c9ff621b6ff47e690902e872204790982040cc725548ed3abd808670219e22",
    "P5_BBBB_rows": "d367b7d208ad8c7fff294978a1857f7557d3755f33df8e4afff6689ad2aec1a7",
}

# (class(a), class(d), class(b), class(c)); the last two classes are sorted.
PATTERNS: tuple[tuple[str, tuple[str, str, str, str]], ...] = (
    ("AA|AA", ("A", "A", "A", "A")),
    ("BA|AA", ("B", "A", "A", "A")),
    ("AB|AA", ("A", "B", "A", "A")),
    ("AA|AB", ("A", "A", "A", "B")),
    ("BB|AA", ("B", "B", "A", "A")),
    ("BA|AB", ("B", "A", "A", "B")),
    ("AB|AB", ("A", "B", "A", "B")),
    ("AA|BB", ("A", "A", "B", "B")),
    ("BB|AB", ("B", "B", "A", "B")),
    ("BA|BB", ("B", "A", "B", "B")),
    ("AB|BB", ("A", "B", "B", "B")),
    ("BB|BB", ("B", "B", "B", "B")),
)
PATTERN_CLASSES = dict(PATTERNS)

ColorWord = tuple[int, int, int, int]
Permutation = tuple[int, ...]


@dataclass
class CompactRow:
    pattern: str
    colors: ColorWord
    orbit_size: int
    repetition_count: int
    facet_count: int
    degree_sha256: str
    gate_sha256: str
    boundary_root_sha256: str
    refinement_sha256: str | None = None
    refinement_discrete: bool | None = None
    refinement_final_cells: int | None = None


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def _stream_update(digest: Any, value: bytes) -> None:
    digest.update(len(value).to_bytes(8, "big"))
    digest.update(value)


def _rss_bytes() -> int:
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(peak if sys.platform == "darwin" else peak * 1024)


def _rss_guard() -> None:
    if _rss_bytes() > RSS_CEILING_BYTES:
        raise MemoryError("the fork census exceeded its explicit 1 GiB RSS ceiling")


def _class_range(letter: str) -> range:
    if letter == "A":
        return range(p5.A_SIZE)
    if letter == "B":
        return range(p5.A_SIZE, p5.STOCK_SIZE)
    raise ArithmeticError(f"unknown colour class {letter}")


def _class_size(letter: str) -> int:
    return p5.A_SIZE if letter == "A" else p5.B_SIZE


def _pattern_b_count(pattern: str) -> int:
    return PATTERN_CLASSES[pattern].count("B")


def _selected_patterns(min_b: int, max_b: int) -> tuple[str, ...]:
    if not (0 <= min_b <= max_b <= 4):
        raise ArithmeticError("fork B-range must satisfy 0 <= min <= max <= 4")
    return tuple(
        name for name, _classes in PATTERNS if min_b <= _pattern_b_count(name) <= max_b
    )


def _raw_pattern_count(pattern: str) -> int:
    ca, cd, cb, cc = PATTERN_CLASSES[pattern]
    a_count = _class_size(ca)
    d_count = _class_size(cd) - int(ca == cd)
    if cb == cc:
        short_count = math.comb(_class_size(cb) - int(ca == cb), 2)
    else:
        short_count = (
            (_class_size(cb) - int(ca == cb))
            * (_class_size(cc) - int(ca == cc))
        )
    return a_count * d_count * short_count


def _cycle_counts_by_class(action: Permutation, letter: str) -> Counter[int]:
    values = set(_class_range(letter))
    counts: Counter[int] = Counter()
    while values:
        start = min(values)
        current = start
        length = 0
        while current in values:
            values.remove(current)
            current = action[current]
            length += 1
        if current != start:
            raise ArithmeticError("the seed action does not preserve an A/B class")
        counts[length] += 1
    return counts


def _burnside_pattern_orbits(
    pattern: str, actions: Sequence[Permutation]
) -> dict[str, int]:
    ca, cd, cb, cc = PATTERN_CLASSES[pattern]
    fixed_sum = 0
    for action in actions:
        cycles = {
            letter: _cycle_counts_by_class(action, letter) for letter in ("A", "B")
        }
        fixed = {letter: cycles[letter][1] for letter in ("A", "B")}
        a_count = fixed[ca]
        if not a_count:
            continue
        d_count = fixed[cd] - int(ca == cd)
        if cb == cc:
            available_fixed = fixed[cb] - int(ca == cb)
            short_count = math.comb(available_fixed, 2) + cycles[cb][2]
        else:
            short_count = (
                (fixed[cb] - int(ca == cb))
                * (fixed[cc] - int(ca == cc))
            )
        fixed_sum += a_count * d_count * short_count
    if fixed_sum % len(actions):
        raise ArithmeticError("a fork Burnside numerator is not divisible")
    return {
        "Burnside_fixed_configuration_sum": fixed_sum,
        "Burnside_denominator": len(actions),
        "orbit_count": fixed_sum // len(actions),
    }


def _local_code(word: ColorWord, pattern: str) -> int:
    code = 0
    for color, letter in zip(word, PATTERN_CLASSES[pattern], strict=True):
        local = color if letter == "A" else color - p5.A_SIZE
        size = _class_size(letter)
        if not 0 <= local < size:
            raise ArithmeticError("a fork colour left its invariant class")
        code = code * size + local
    return code


def _pattern_color_orbits(
    pattern: str, actions: Sequence[Permutation]
) -> tuple[dict[str, Any], ...]:
    ca, cd, cb, cc = PATTERN_CLASSES[pattern]
    dense_size = math.prod(_class_size(letter) for letter in (ca, cd, cb, cc))
    seen = bytearray(dense_size)
    rows: list[dict[str, Any]] = []
    short_iterator: Iterable[tuple[int, int]]
    if cb == cc:
        short_iterator_factory = lambda: combinations(_class_range(cb), 2)
    else:
        short_iterator_factory = lambda: product(_class_range(cb), _class_range(cc))
    for a in _class_range(ca):
        for d in _class_range(cd):
            if a == d:
                continue
            for b, c in short_iterator_factory():
                if a == b or a == c:
                    continue
                word = (a, d, b, c)
                code = _local_code(word, pattern)
                if seen[code]:
                    continue
                orbit: set[ColorWord] = set()
                for action in actions:
                    image_b, image_c = sorted((action[b], action[c]))
                    image = (action[a], action[d], image_b, image_c)
                    orbit.add(image)
                for image in orbit:
                    seen[_local_code(image, pattern)] = 1
                rows.append(
                    {
                        "fork_class_pattern": pattern,
                        "colors_zero_based_stock": word,
                        "unrooted_coloring_orbit_size_under_seed_action": len(orbit),
                    }
                )
    raw_count = _raw_pattern_count(pattern)
    burnside = _burnside_pattern_orbits(pattern, actions)
    if not (
        sum(seen) == raw_count
        and sum(
            row["unrooted_coloring_orbit_size_under_seed_action"] for row in rows
        )
        == raw_count
        and len(rows) == burnside["orbit_count"]
    ):
        raise ArithmeticError(f"the exact {pattern} fork quotient failed")
    return tuple(rows)


class ExactForkBoundaryEngine(p5.ExactP5BoundaryEngine):
    """Boundary engine for transforms ``1,Ra,Rb,Rc,Ra Rd``."""

    def _transformations(self, stock_colors: ColorWord) -> tuple[np.ndarray, ...]:
        a, d, b, c = (color + p5.STOCK_START for color in stock_colors)
        transformations = (
            self.identity,
            self.reflections[a],
            self.reflections[b],
            self.reflections[c],
            self.reflections[a] @ self.reflections[d],
        )
        if len({p5._matrix_key(matrix) for matrix in transformations}) != 5:
            raise ArithmeticError("a proper fork word did not give five chambers")
        return transformations

    def boundary(self, stock_colors: ColorWord) -> tuple[np.ndarray, dict[str, Any]]:
        if not (
            len(stock_colors) == 4
            and all(0 <= color < p5.STOCK_SIZE for color in stock_colors)
        ):
            raise ArithmeticError("invalid fork colour word")
        a, d, b, c = stock_colors
        if a == d or len({a, b, c}) != 3 or b >= c:
            raise ArithmeticError("fork proper-colouring convention failed")
        absolute = tuple(color + p5.STOCK_START for color in stock_colors)
        absolute_a, absolute_d, absolute_b, absolute_c = absolute
        transformations = self._transformations(stock_colors)
        transform_index = {
            p5._matrix_key(matrix): index for index, matrix in enumerate(transformations)
        }
        adjacency: set[tuple[int, int, int]] = set()
        internal: set[tuple[int, ...]] = set()
        for chamber, transformation in enumerate(transformations):
            for color in set(absolute):
                neighbor = transformation @ self.reflections[color]
                other = transform_index.get(p5._matrix_key(neighbor))
                if other is None:
                    continue
                left_key = star._primitive_key(transformation @ self.roots[color])
                right_key = star._primitive_key(neighbor @ self.roots[color])
                if left_key != right_key:
                    raise ArithmeticError("the two sides of an internal fork mirror differ")
                internal.add(left_key)
                if chamber < other:
                    adjacency.add((chamber, other, color))
        expected_adjacency = {
            (0, 1, absolute_a),
            (1, 4, absolute_d),
            (0, 2, absolute_b),
            (0, 3, absolute_c),
        }
        if adjacency != expected_adjacency or len(internal) != 4:
            raise ArithmeticError("the transformations do not induce the marked fork")

        union_interior = sum(
            (transformation @ self.interior for transformation in transformations),
            start=np.zeros(p5.RANK, dtype=np.int64),
        )
        boundary: dict[tuple[int, ...], tuple[int, ...]] = {}
        maximum_orientation_pairing = 0
        maximum_transformed_root = 0
        for transformation in transformations:
            images = self.roots @ transformation.T
            image_maximum = int(np.max(np.abs(images)))
            maximum_transformed_root = max(maximum_transformed_root, image_maximum)
            orientation_bound = (
                p5.RANK**2
                * self.gram_maximum
                * image_maximum
                * int(np.max(np.abs(union_interior)))
            )
            if orientation_bound >= p5.INT64_MAX:
                raise ArithmeticError("a fork orientation pairing could overflow int64")
            pairings = images @ self.gram @ union_interior
            maximum_orientation_pairing = max(
                maximum_orientation_pairing, int(np.max(np.abs(pairings)))
            )
            for image, pairing in zip(images, pairings, strict=True):
                key = star._primitive_key(image)
                if key in internal:
                    continue
                if pairing == 0:
                    raise ArithmeticError("the summed fork interior lies on a boundary")
                oriented = image if pairing > 0 else -image
                value = tuple(int(entry) for entry in oriented)
                existing = boundary.setdefault(key, value)
                if existing != value:
                    raise ArithmeticError("a fork wall received opposite orientations")

        ordered_keys = tuple(sorted(boundary))
        if not self.core_keys <= set(ordered_keys):
            raise ArithmeticError("the fixed Picard core left a fork boundary")
        roots = np.array([boundary[key] for key in ordered_keys], dtype=np.int64)
        coordinate_maximum = int(np.max(np.abs(roots)))
        pair_bound = p5.RANK**2 * self.gram_maximum * coordinate_maximum**2
        interior_maximum = int(np.max(np.abs(union_interior)))
        interior_bound = p5.RANK**2 * self.gram_maximum * interior_maximum**2
        if max(pair_bound, interior_bound) >= p5.INT64_MAX:
            raise ArithmeticError("a fork boundary Gram product could overflow int64")
        covectors = roots @ self.gram
        squares = np.einsum("ij,ij->i", roots, covectors)
        if not np.all(squares < 0):
            raise ArithmeticError("a fork boundary root is not spacelike")
        primitive_integral = 0
        for root, square, covector in zip(roots, squares, covectors, strict=True):
            if not (
                star._content(root) == 1
                and all((2 * int(value)) % int(square) == 0 for value in covector)
            ):
                raise ArithmeticError("a fork boundary reflection is nonintegral")
            primitive_integral += 1
        interior_pairings = roots @ self.gram @ union_interior
        interior_square = int(union_interior @ self.gram @ union_interior)
        if not (interior_square > 0 and np.all(interior_pairings > 0)):
            raise ArithmeticError("a fork strict interior failed")
        return roots, {
            "facet_count": len(roots),
            "chamber_count": 5,
            "internal_mirror_count": len(internal),
            "marked_chamber_adjacency": [
                [left + 1, right + 1, color + 1]
                for left, right, color in sorted(adjacency)
            ],
            "strict_interior_square": interior_square,
            "minimum_strict_interior_facet_pairing": int(np.min(interior_pairings)),
            "primitive_integral_reflection_root_count": primitive_integral,
            "fixed_boundary_core_Smith_diagonal": self.core_smith,
            "full_linear_rank": p5.RANK,
            "marked_boundary_root_sha256": star._json_sha256(roots.tolist()),
            "maximum_boundary_coordinate_absolute": coordinate_maximum,
            "maximum_transformed_seed_root_coordinate_absolute": maximum_transformed_root,
            "maximum_orientation_pairing_absolute": maximum_orientation_pairing,
            "maximum_chamber_transform_coordinate_absolute": max(
                int(np.max(np.abs(matrix))) for matrix in transformations
            ),
        }


def _degree_replay(
    engine: Any, colors: ColorWord, expected_digest: str | None = None
) -> tuple[bytes, Any]:
    boundary, _boundary_data = engine.boundary(colors)
    labels, degree_key, _ordinary = engine.ordinary_data(boundary)
    if expected_digest is not None and hashlib.sha256(degree_key).hexdigest() != expected_digest:
        raise ArithmeticError("a degree-profile digest changed on exact replay")
    return degree_key, labels


def _refinement(
    engine: Any, colors: ColorWord, expected_degree: bytes
) -> tuple[bytes, dict[str, Any]]:
    degree, labels = _degree_replay(engine, colors)
    if degree != expected_degree:
        raise ArithmeticError("a literal degree profile changed before refinement")
    return four._canonical_color_refinement(labels)


def _gzip_archived_degree_hashes(path: Path, expected: int) -> tuple[str, ...]:
    values: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            row = json.loads(line)
            digest = row.get("ordinary_degree_profile_multiset_sha256")
            if not isinstance(digest, str) or len(digest) != 64:
                raise ArithmeticError(f"a compact P5 degree digest changed in {path}")
            values.append(digest)
    if len(values) != expected:
        raise ArithmeticError(f"the compact P5 row count changed for {path}")
    return tuple(values)


def _p5_prior_rows_for_b(
    target_b: int, actions: Sequence[Permutation]
) -> Iterator[tuple[str, int, ColorWord, str, str]]:
    """Yield every completed P5 row with exactly ``target_b`` B edges."""

    if target_b <= 2:
        hashes = p5stream._archived_degree_hashes(P5_ARCHIVE, 20_693)
        offset = 0
        for pattern in p5.PRODUCTION_PATTERNS:
            rows = p5._pattern_color_orbits(pattern, actions)
            pattern_hashes = hashes[offset : offset + len(rows)]
            if len(pattern_hashes) != len(rows):
                raise ArithmeticError("the P5 B<=2 archive/order correspondence changed")
            if pattern.count("B") == target_b:
                for index, (row, digest) in enumerate(
                    zip(rows, pattern_hashes, strict=True), start=1
                ):
                    yield (
                        f"P5_{pattern}",
                        index,
                        tuple(int(value) for value in row["colors_zero_based_stock"]),
                        digest,
                        "P5",
                    )
            offset += len(rows)
        if offset != len(hashes):
            raise ArithmeticError("the complete B<=2 P5 row count changed")
        return

    if target_b == 3:
        high_patterns = (
            ("ABBB", P5_ABBB_ARCHIVE, 36_247, False),
            ("BABB", P5_BABB_ROWS, 36_648, True),
        )
    elif target_b == 4:
        high_patterns = (("BBBB", P5_BBBB_ROWS, 108_159, True),)
    else:
        raise ArithmeticError("a P5 same-facet B count left 0..4")
    for pattern, archive, expected, compressed in high_patterns:
        rows = p5._pattern_color_orbits(pattern, actions)
        hashes = (
            _gzip_archived_degree_hashes(archive, expected)
            if compressed
            else p5stream._archived_degree_hashes(archive, expected)
        )
        if len(rows) != expected:
            raise ArithmeticError(f"the complete P5 {pattern} denominator changed")
        for index, (row, digest) in enumerate(
            zip(rows, hashes, strict=True), start=1
        ):
            yield (
                f"P5_{pattern}",
                index,
                tuple(int(value) for value in row["colors_zero_based_stock"]),
                digest,
                "P5",
            )


def _prior_rows(
    min_b: int, max_b: int, actions: Sequence[Permutation]
) -> Iterator[tuple[str, int, ColorWord, str, str]]:
    """Yield every completed P5/star row with exactly the selected headers."""

    for target_b in range(min_b, max_b + 1):
        yield from _p5_prior_rows_for_b(target_b, actions)

    star_rows = star._star_color_orbits(actions)
    star_hashes = p5stream._archived_degree_hashes(STAR_ARCHIVE, 16_872)
    for index, (row, digest) in enumerate(
        zip(star_rows, star_hashes, strict=True), start=1
    ):
        colors = tuple(int(value) for value in row["colors_zero_based_stock"])
        if min_b <= sum(color >= p5.A_SIZE for color in colors) <= max_b:
            yield ("K1,4_star", index, colors, digest, "star")


def _compact_payload(index: int, row: CompactRow) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "orbit_number": index,
        "fork_class_pattern": row.pattern,
        "edge_roles": ["central_long", "long_extension", "short_1", "short_2"],
        "edge_colors_one_based": [
            color + p5.STOCK_START + 1 for color in row.colors
        ],
        "unrooted_coloring_orbit_size_under_seed_action": row.orbit_size,
        "nonadjacent_color_repetition_count": row.repetition_count,
        "facet_count": row.facet_count,
        "marked_boundary_root_sha256": row.boundary_root_sha256,
        "ordinary_degree_profile_multiset_sha256": row.degree_sha256,
        "exact_gate_record_sha256": row.gate_sha256,
        "color_refinement_applied": row.refinement_sha256 is not None,
    }
    if row.refinement_sha256 is not None:
        payload.update(
            {
                "color_refinement_invariant_sha256": row.refinement_sha256,
                "color_refinement_discrete_canonical_form": row.refinement_discrete,
                "color_refinement_final_cell_count": row.refinement_final_cells,
            }
        )
    return payload


def _write_rows(path: Path, rows: Sequence[CompactRow]) -> tuple[str, str]:
    transcript = hashlib.sha256()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as compressed:
            for index, row in enumerate(rows, start=1):
                rendered = _json_bytes(_compact_payload(index, row)) + b"\n"
                transcript.update(rendered)
                compressed.write(rendered)
    return transcript.hexdigest(), _sha256(path)


def build(
    *, min_b: int, max_b: int, rows_output: Path, verbose: bool = False
) -> dict[str, Any]:
    patterns = _selected_patterns(min_b, max_b)
    small, pilot, denominator = p5._pinned_predecessors()
    if not (
        denominator["tree_shapes"]["degree_3_fork"][
            "seed_and_tree_automorphism_orbits"
        ]
        == 199_362
        and _sha256(DENOMINATOR_ARCHIVE)
        == p5.PINNED_SHA256["five_vertex_denominator_archive"]
        and _sha256(P5_ARCHIVE)
        == p5stream.PINNED_SHA256["P5_B_le_2_archive"]
        and _sha256(STAR_ARCHIVE)
        == p5stream.PINNED_SHA256["five_vertex_star_archive"]
    ):
        raise ArithmeticError("the exact fork denominator prerequisite changed")
    selected_high_b_inputs: dict[str, str] = {}
    if max_b >= 3:
        selected_high_b_inputs.update(
            {
                "P5_ABBB_archive": _sha256(P5_ABBB_ARCHIVE),
                "P5_BABB_rows": _sha256(P5_BABB_ROWS),
            }
        )
    if max_b >= 4:
        selected_high_b_inputs["P5_BBBB_rows"] = _sha256(P5_BBBB_ROWS)
    if selected_high_b_inputs != {
        key: HIGH_B_PRIOR_SHA256[key] for key in selected_high_b_inputs
    }:
        raise ArithmeticError(
            f"a frozen high-B same-facet P5 archive changed: {selected_high_b_inputs}"
        )
    pinned_input_hashes: dict[str, Any] = {
        "five_vertex_tree_denominator": _sha256(DENOMINATOR_ARCHIVE),
        "complete_P5_B_le_2_census": _sha256(P5_ARCHIVE),
        "complete_star_census": _sha256(STAR_ARCHIVE),
    }
    if selected_high_b_inputs:
        pinned_input_hashes["selected_high_B_same_facet_P5_archives"] = (
            selected_high_b_inputs
        )
    carrier, gram, roots, interior, bounded = p5.base._load_inputs()
    actions = star._restricted_actions()
    if not all(
        all(
            (color < p5.A_SIZE) == (action[color] < p5.A_SIZE)
            for color in range(p5.STOCK_SIZE)
        )
        for action in actions
    ):
        raise ArithmeticError("the seed action no longer preserves the A/B split")
    stock_gram_sha256 = star._root_gram_action_gate(gram, roots, actions)

    color_rows: list[dict[str, Any]] = []
    pattern_denominator: dict[str, dict[str, Any]] = {}
    for pattern in patterns:
        tick = time.perf_counter()
        shard = _pattern_color_orbits(pattern, actions)
        burnside = _burnside_pattern_orbits(pattern, actions)
        pattern_denominator[pattern] = {
            "class_roles": list(PATTERN_CLASSES[pattern]),
            "B_edge_occurrence_count": _pattern_b_count(pattern),
            "raw_unrooted_proper_colorings": _raw_pattern_count(pattern),
            **burnside,
            "explicit_orbit_cover_count": len(shard),
            "explicit_orbit_size_sum": sum(
                row["unrooted_coloring_orbit_size_under_seed_action"] for row in shard
            ),
        }
        color_rows.extend(shard)
        if verbose:
            print(
                f"orbit_shard={pattern} rows={len(shard)} "
                f"seconds={time.perf_counter()-tick:.2f}",
                flush=True,
            )
        _rss_guard()

    engine = ExactForkBoundaryEngine(gram, roots, interior)
    digest_groups: defaultdict[str, list[int]] = defaultdict(list)
    rows: list[CompactRow] = []
    relation_histograms: Counter[tuple[tuple[str, int], ...]] = Counter()
    facet_histogram: Counter[int] = Counter()
    raw_facet_histogram: Counter[int] = Counter()
    orbit_size_histogram: Counter[int] = Counter()
    repetition_histogram: Counter[int] = Counter()
    pattern_type_histogram: Counter[str] = Counter()
    gate_stream = hashlib.sha256()
    degree_stream = hashlib.sha256()
    coordinate_maximum = transformed_root_maximum = transform_maximum = pairing_maximum = 0
    minimum_interior_pairing: int | None = None
    boundary_started = time.perf_counter()
    for index, color_row in enumerate(color_rows):
        pattern = str(color_row["fork_class_pattern"])
        colors = tuple(int(value) for value in color_row["colors_zero_based_stock"])
        boundary, boundary_data = engine.boundary(colors)
        labels, degree_key, ordinary_data = engine.ordinary_data(boundary)
        del labels
        facets = int(boundary_data["facet_count"])
        expected_facets = 752 + 9 * _pattern_b_count(pattern)
        if facets != expected_facets:
            raise ArithmeticError("the five-fork facet formula changed")
        gate_data = {**boundary_data, **ordinary_data}
        if not (
            boundary_data["chamber_count"] == 5
            and boundary_data["internal_mirror_count"] == 4
            and boundary_data["strict_interior_square"] > 0
            and boundary_data["minimum_strict_interior_facet_pairing"] > 0
            and boundary_data["primitive_integral_reflection_root_count"] == facets
            and boundary_data["fixed_boundary_core_Smith_diagonal"] == [1] * p5.RANK
            and boundary_data["full_linear_rank"] == p5.RANK
            and ordinary_data["primitive_integral_reflection_root_count"] == facets
            and sum(ordinary_data["ordinary_Coxeter_relation_histogram"].values())
            == facets * (facets - 1) // 2
        ):
            raise ArithmeticError("a fork arithmetic/Coxeter gate failed")
        degree_digest = hashlib.sha256(degree_key).hexdigest()
        digest_groups[degree_digest].append(index)
        _stream_update(degree_stream, degree_key)
        gate_rendered = _json_bytes(gate_data)
        _stream_update(gate_stream, gate_rendered)
        orbit_size = int(color_row["unrooted_coloring_orbit_size_under_seed_action"])
        repetition = 4 - len(set(colors))
        rows.append(
            CompactRow(
                pattern=pattern,
                colors=colors,
                orbit_size=orbit_size,
                repetition_count=repetition,
                facet_count=facets,
                degree_sha256=degree_digest,
                gate_sha256=hashlib.sha256(gate_rendered).hexdigest(),
                boundary_root_sha256=boundary_data["marked_boundary_root_sha256"],
            )
        )
        relation_histograms[
            tuple(sorted(ordinary_data["ordinary_Coxeter_relation_histogram"].items()))
        ] += 1
        facet_histogram[facets] += 1
        raw_facet_histogram[facets] += orbit_size
        orbit_size_histogram[orbit_size] += 1
        repetition_histogram[repetition] += 1
        pattern_type_histogram[pattern] += 1
        coordinate_maximum = max(
            coordinate_maximum, boundary_data["maximum_boundary_coordinate_absolute"]
        )
        transformed_root_maximum = max(
            transformed_root_maximum,
            boundary_data["maximum_transformed_seed_root_coordinate_absolute"],
        )
        transform_maximum = max(
            transform_maximum,
            boundary_data["maximum_chamber_transform_coordinate_absolute"],
        )
        pairing_maximum = max(
            pairing_maximum, ordinary_data["maximum_boundary_pairing_absolute"]
        )
        value = int(boundary_data["minimum_strict_interior_facet_pairing"])
        minimum_interior_pairing = (
            value if minimum_interior_pairing is None else min(minimum_interior_pairing, value)
        )
        if (index + 1) % 100 == 0:
            _rss_guard()
        if verbose and ((index + 1) % 500 == 0 or index + 1 == len(color_rows)):
            print(
                f"boundary_progress={index+1}/{len(color_rows)} "
                f"elapsed_seconds={time.perf_counter()-boundary_started:.1f} "
                f"peak_rss_gib={_rss_bytes()/(1<<30):.3f}",
                flush=True,
            )

    digest_collision_groups = [
        indices for indices in digest_groups.values() if len(indices) > 1
    ]
    literal_replayed = 0
    literal_degree_group_count = 0
    literal_degree_member_count = 0
    refined_current_indices: set[int] = set()
    refinement_discrete = 0
    for digest_group in digest_collision_groups:
        exact: defaultdict[bytes, list[int]] = defaultdict(list)
        for index in digest_group:
            degree, _labels = _degree_replay(
                engine, rows[index].colors, rows[index].degree_sha256
            )
            exact[degree].append(index)
            literal_replayed += 1
        for degree, indices in exact.items():
            if len(indices) <= 1:
                continue
            literal_degree_group_count += 1
            literal_degree_member_count += len(indices)
            literal_refinements: list[bytes] = []
            for index in indices:
                invariant, data = _refinement(engine, rows[index].colors, degree)
                if invariant in literal_refinements:
                    raise ArithmeticError(
                        "UNKNOWN: exact refinement did not separate a fork degree group"
                    )
                literal_refinements.append(invariant)
                rows[index].refinement_sha256 = hashlib.sha256(invariant).hexdigest()
                rows[index].refinement_discrete = bool(
                    data["color_refinement_discrete_canonical_form"]
                )
                rows[index].refinement_final_cells = int(
                    data["color_refinement_final_cell_count"]
                )
                refinement_discrete += int(rows[index].refinement_discrete)
                refined_current_indices.add(index)
        _rss_guard()

    current_digests = set(digest_groups)
    p5_engine = p5.ExactP5BoundaryEngine(gram, roots, interior)
    star_engine = star.ExactStarBoundaryEngine(gram, roots, interior)
    prior_same_facet_count = 0
    prior_hash_candidates = 0
    prior_literal_replayed = 0
    cross_literal_matches: list[tuple[str, int, int, ColorWord, str]] = []
    for source, prior_index, colors, digest, engine_name in _prior_rows(
        min_b, max_b, actions
    ):
        prior_same_facet_count += 1
        if digest not in current_digests:
            continue
        prior_hash_candidates += 1
        prior_engine = p5_engine if engine_name == "P5" else star_engine
        prior_degree, _labels = _degree_replay(prior_engine, colors, digest)
        prior_literal_replayed += 1
        for current_index in digest_groups[digest]:
            current_degree, _labels = _degree_replay(
                engine, rows[current_index].colors, digest
            )
            if current_degree == prior_degree:
                cross_literal_matches.append(
                    (source, prior_index, current_index, colors, engine_name)
                )
        _rss_guard()

    for source, prior_index, current_index, prior_colors, engine_name in cross_literal_matches:
        prior_engine = p5_engine if engine_name == "P5" else star_engine
        current_degree, _labels = _degree_replay(
            engine, rows[current_index].colors, rows[current_index].degree_sha256
        )
        prior_degree, _labels = _degree_replay(
            prior_engine, prior_colors, rows[current_index].degree_sha256
        )
        if current_degree != prior_degree:
            raise ArithmeticError("a recorded fork/prior degree match changed")
        current_refinement, current_data = _refinement(
            engine, rows[current_index].colors, current_degree
        )
        prior_refinement, _prior_data = _refinement(
            prior_engine, prior_colors, prior_degree
        )
        if current_refinement == prior_refinement:
            raise ArithmeticError(
                f"UNKNOWN: fork row {current_index+1} and {source} row {prior_index} "
                "retain the same exact refinement invariant"
            )
        rows[current_index].refinement_sha256 = hashlib.sha256(
            current_refinement
        ).hexdigest()
        rows[current_index].refinement_discrete = bool(
            current_data["color_refinement_discrete_canonical_form"]
        )
        rows[current_index].refinement_final_cells = int(
            current_data["color_refinement_final_cell_count"]
        )
        refined_current_indices.add(current_index)
        _rss_guard()

    rows_transcript_sha256, rows_gzip_sha256 = _write_rows(rows_output, rows)
    raw_count = sum(_raw_pattern_count(pattern) for pattern in patterns)
    if sum(row.orbit_size for row in rows) != raw_count:
        raise ArithmeticError("fork transcript orbit sizes do not cover the raw shard")
    checks = {
        "frozen_seed_and_fork_denominator_replayed": True,
        "exact_action_preserves_stock_Gram_and_AB_split": True,
        "explicit_orbit_cover_and_Burnside_agree_pattern_by_pattern": True,
        "every_boundary_rebuilt_with_ordinary_Coxeter_matrix": True,
        "every_root_primitive_integral_and_reflective": True,
        "every_boundary_has_strict_timelike_interior": True,
        "every_boundary_contains_fixed_Smith_1_power_13_core": True,
        "different_hashes_used_only_as_exact_inequality_filters": True,
        "every_same_digest_bucket_replayed_with_literal_degree_bytes": True,
        "literal_degree_collisions_separated_by_exact_refinement": True,
        "all_rows_pairwise_nonisomorphic_within_selected_fork_shards": True,
        "all_same_facet_P5_and_star_hash_candidates_replayed_literally": True,
        "all_literal_cross_shape_matches_separated_by_exact_refinement": True,
        "Allcock_finite_volume_noncompact_multiplier_five_gate": True,
        "same_full_Picard_weak_K3_carrier_retained": True,
        "RSS_ceiling_at_most_1_GiB": RSS_CEILING_BYTES <= 1 << 30,
        "worldwide_priority_explicitly_open": True,
    }
    return {
        "schema": SCHEMA,
        "status_date": STATUS_DATE,
        "all_checks_pass": all(checks.values()),
        "selected_B_range": [min_b, max_b],
        "selected_patterns": list(patterns),
        "pinned_inputs": pinned_input_hashes,
        "seed": {
            "id": "JK-H12/172",
            "hyperbolic_dimension": 12,
            "facet_count": 172,
            "selected_pairwise_disjoint_doubling_stock": "W61--W172",
            "stock_size": 112,
            "folded_seed_action_order": len(actions),
            "stock_Gram_sha256": stock_gram_sha256,
            "full_Picard_lattice": "U(2)+D5(-1)^2+<-4>",
            "K3_moduli_dimension": 7,
            "source_carrier_theorem": carrier["theorem"],
        },
        "coloring_denominator": {
            "pattern_rows": pattern_denominator,
            "raw_unrooted_proper_colorings": raw_count,
            "explicit_orbit_cover_count": len(rows),
            "explicit_orbit_size_sum": sum(row.orbit_size for row in rows),
            "complete_fork_orbit_denominator": 199_362,
        },
        "theorem": {
            "scope": f"complete degree-three fork B={min_b}..{max_b} layer",
            "ordinary_Coxeter_type_count": len(rows),
            "pairwise_nonisomorphic_within_selected_fork_shards": True,
            "separated_from_all_completed_same_facet_P5_and_star_rows": True,
            "new_pairwise_nonisomorphic_increment_over_completed_five_vertex_shapes": len(rows),
            "type_facet_histogram": {
                str(key): value for key, value in sorted(facet_histogram.items())
            },
            "pattern_type_histogram": dict(sorted(pattern_type_histogram.items())),
            "raw_unrooted_coloring_facet_histogram": {
                str(key): value for key, value in sorted(raw_facet_histogram.items())
            },
            "reflection_index_and_volume_multiple": 5,
            "finite_volume": True,
            "compact": False,
        },
        "compact_exact_transcript": {
            "rows_path": str(rows_output),
            "row_count": len(rows),
            "uncompressed_canonical_JSONL_sha256": rows_transcript_sha256,
            "deterministic_gzip_sha256": rows_gzip_sha256,
            "exact_gate_record_stream_sha256": gate_stream.hexdigest(),
            "exact_degree_profile_byte_stream_sha256": degree_stream.hexdigest(),
            "hash_role": "different digests prove unequal bytes; equal digests are always replayed literally",
        },
        "ordinary_Coxeter_nonisomorphism_certificate": {
            "degree_digest_class_count": len(digest_groups),
            "degree_digest_collision_group_count": len(digest_collision_groups),
            "degree_digest_collision_member_count": sum(map(len, digest_collision_groups)),
            "literal_degree_replayed_member_count": literal_replayed,
            "literal_degree_collision_group_count": literal_degree_group_count,
            "literal_degree_collision_member_count": literal_degree_member_count,
            "color_refined_current_member_count": len(refined_current_indices),
            "color_refined_discrete_current_count": refinement_discrete,
            "all_literal_degree_collision_groups_refinement_distinct": True,
        },
        "cross_shape_exact_separation": {
            "prior_same_facet_P5_and_star_row_count": prior_same_facet_count,
            "archived_degree_digest_candidate_count": prior_hash_candidates,
            "prior_literal_degree_replay_count": prior_literal_replayed,
            "literal_cross_shape_degree_match_count": len(cross_literal_matches),
            "exact_refinement_separated_pair_count": len(cross_literal_matches),
            "unresolved_cross_shape_pair_count": 0,
        },
        "aggregate_arithmetic_gates": {
            "ordinary_relation_histogram_inventory": [
                {"histogram": dict(histogram), "type_count": count}
                for histogram, count in sorted(relation_histograms.items())
            ],
            "minimum_strict_interior_facet_pairing": minimum_interior_pairing,
            "maximum_boundary_coordinate_absolute": coordinate_maximum,
            "maximum_transformed_seed_root_coordinate_absolute": transformed_root_maximum,
            "maximum_chamber_transform_coordinate_absolute": transform_maximum,
            "maximum_boundary_pairing_absolute": pairing_maximum,
            "seed_orbit_size_histogram": {
                str(key): value for key, value in sorted(orbit_size_histogram.items())
            },
            "nonadjacent_color_repetition_type_histogram": {
                str(key): value for key, value in sorted(repetition_histogram.items())
            },
        },
        "resource_certificate": {
            "algorithm": "one boundary at a time; compact row digests; literal replay of every equality candidate",
            "RSS_ceiling_bytes": RSS_CEILING_BYTES,
            "in_process_peak_RSS_guard": True,
            "recommended_outer_watchdog": "python3 scripts/run_with_rss_limit.py --memory-gib 1.0 -- COMMAND",
            "all_Gram_matrices_materialized_together": False,
        },
        "Allcock_and_carrier": {
            "finite_subtree_theorem": "each row is a five-copy degree-three-fork fundamental polytope",
            "finite_volume": bool(pilot["result"]["finite_volume"]),
            "compact": bool(pilot["result"]["compact"]),
            "reflection_index_and_volume_multiplier": pilot["result"][
                "reflection_index_and_volume_multiplier"
            ],
            "same_full_Picard_K3_family": True,
            "fixed_index_one_divisor_lattice_core": True,
            "weak_divisor_reflection_carrier_only": True,
            "full_nef_effectivity_or_surface_action_claim": False,
        },
        "lineage_and_claim_boundary": {
            "direction": "SURFACE_TO_POLYTOPE",
            "selected_fork_layer_complete": True,
            "degree_three_fork_complete_census_claimed": min_b == 0 and max_b == 4,
            "all_five_vertex_shapes_complete_census_claimed": min_b == 0 and max_b == 4,
            "worldwide_priority": "open",
            "worldwide_novelty_claim": False,
            "bounded_CoxIter_commit": bounded["commit"],
            "weak_carrier_not_full_nef": True,
        },
        "checks": checks,
    }


def denominator_only(*, verbose: bool = False) -> dict[str, Any]:
    _small, _pilot, denominator = p5._pinned_predecessors()
    if (
        denominator["tree_shapes"]["degree_3_fork"][
            "seed_and_tree_automorphism_orbits"
        ]
        != 199_362
    ):
        raise ArithmeticError("the frozen aggregate fork denominator changed")
    actions = star._restricted_actions()
    rows: dict[str, Any] = {}
    total_raw = total_orbits = 0
    for pattern, _classes in PATTERNS:
        tick = time.perf_counter()
        orbits = _pattern_color_orbits(pattern, actions)
        burnside = _burnside_pattern_orbits(pattern, actions)
        raw = _raw_pattern_count(pattern)
        rows[pattern] = {
            "class_roles": list(PATTERN_CLASSES[pattern]),
            "B_edge_occurrence_count": _pattern_b_count(pattern),
            "raw_unrooted_proper_colorings": raw,
            "orbit_count": len(orbits),
            **burnside,
        }
        total_raw += raw
        total_orbits += len(orbits)
        if verbose:
            print(
                f"pattern={pattern} raw={raw} orbits={len(orbits)} "
                f"seconds={time.perf_counter()-tick:.2f}",
                flush=True,
            )
        _rss_guard()
    if total_raw != 75_897_360 or total_orbits != 199_362:
        raise ArithmeticError(
            f"the fork denominator partition failed: raw={total_raw}, orbits={total_orbits}"
        )
    by_b = Counter()
    raw_by_b = Counter()
    for row in rows.values():
        by_b[row["B_edge_occurrence_count"]] += row["orbit_count"]
        raw_by_b[row["B_edge_occurrence_count"]] += row[
            "raw_unrooted_proper_colorings"
        ]
    checks = {
        "all_twelve_AB_role_patterns_partition_the_fork": len(rows) == 12,
        "every_pattern_has_explicit_cover_equal_to_Burnside": True,
        "raw_pattern_sum_is_75897360": total_raw == 75_897_360,
        "orbit_pattern_sum_is_199362": total_orbits == 199_362,
        "frozen_aggregate_denominator_recovered": True,
    }
    return {
        "schema": "coxeter.jacobian_kummer.h12_five_vertex_fork_pattern_denominator.v1",
        "status_date": STATUS_DATE,
        "all_checks_pass": all(checks.values()),
        "pinned_aggregate_denominator_sha256": _sha256(DENOMINATOR_ARCHIVE),
        "patterns": rows,
        "raw_coloring_count_by_B": {
            str(key): value for key, value in sorted(raw_by_b.items())
        },
        "orbit_count_by_B": {
            str(key): value for key, value in sorted(by_b.items())
        },
        "total_raw": total_raw,
        "total_orbits": total_orbits,
        "checks": checks,
    }


def _default_paths(min_b: int, max_b: int) -> tuple[Path, Path]:
    stem = f"jacobian_kummer_h12_five_vertex_fork_b{min_b}_to_b{max_b}_streaming_census"
    return (
        ROOT / "data/enumerations" / f"{stem}.json",
        ROOT / "data/enumerations" / f"{stem}_rows.jsonl.gz",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-b", type=int, default=0)
    parser.add_argument("--max-b", type=int, default=1)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--rows", type=Path)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--denominator-only", action="store_true")
    parser.add_argument("--denominator-output", type=Path, default=FORK_DENOMINATOR_OUTPUT)
    parser.add_argument("--check-denominator-archive", action="store_true")
    arguments = parser.parse_args()
    if arguments.denominator_only:
        payload = denominator_only(verbose=True)
        rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
        if arguments.check_denominator_archive:
            if (
                not arguments.denominator_output.is_file()
                or arguments.denominator_output.read_text(encoding="utf-8") != rendered
            ):
                raise SystemExit("the fork pattern denominator archive differs")
            print("denominator_archive_check=PASS")
        else:
            arguments.denominator_output.parent.mkdir(parents=True, exist_ok=True)
            arguments.denominator_output.write_text(rendered, encoding="utf-8")
            print(f"wrote_denominator={arguments.denominator_output}")
        print(f"denominator_sha256={hashlib.sha256(rendered.encode()).hexdigest()}")
        return
    default_summary, default_rows = _default_paths(arguments.min_b, arguments.max_b)
    summary_path = arguments.summary or default_summary
    rows_path = arguments.rows or default_rows
    with tempfile.TemporaryDirectory(prefix="jk_h12_fork_", dir="/private/tmp") as temporary:
        temporary_rows = Path(temporary) / "rows.jsonl.gz"
        payload = build(
            min_b=arguments.min_b,
            max_b=arguments.max_b,
            rows_output=temporary_rows,
            verbose=True,
        )
        payload["compact_exact_transcript"]["rows_path"] = str(rows_path.relative_to(ROOT))
        rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
        if arguments.check_archive:
            if not summary_path.is_file() or not rows_path.is_file():
                raise SystemExit("the fork summary or row archive is missing")
            if summary_path.read_text(encoding="utf-8") != rendered:
                raise SystemExit("the fork summary differs from exact replay")
            if rows_path.read_bytes() != temporary_rows.read_bytes():
                raise SystemExit("the deterministic fork rows differ from exact replay")
            print("archive_check=PASS")
        else:
            summary_path.parent.mkdir(parents=True, exist_ok=True)
            rows_path.parent.mkdir(parents=True, exist_ok=True)
            summary_path.write_text(rendered, encoding="utf-8")
            temporary_rows.replace(rows_path)
            print(f"wrote_summary={summary_path}")
            print(f"wrote_rows={rows_path}")
    print("all_checks_pass=true")
    print(f"B_range={arguments.min_b}..{arguments.max_b}")
    print(f"type_count={payload['theorem']['ordinary_Coxeter_type_count']}")
    print(f"peak_rss_gib={_rss_bytes()/(1<<30):.3f}")
    print(f"summary_sha256={_sha256(summary_path)}")
    print(f"rows_sha256={_sha256(rows_path)}")


if __name__ == "__main__":
    main()
