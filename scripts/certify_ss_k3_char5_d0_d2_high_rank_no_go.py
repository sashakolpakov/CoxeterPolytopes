"""Certify the high-dimensional fixed-slice no-go for the char-5 D0--D2 atlas.

The input is the exact computation data published with Katsura--Kondo--
Shimada's study of the supersingular K3 surface of Artin invariant 1 in
characteristic 5.  The paper gives complete wall lists and factorizations of
the period-preserving groups ``Aut_X(D_i)`` for three highlighted induced
chambers ``D0``, ``D1`` and ``D2``.  It does *not* publish the full abstract
lattice stabilizers ``Aut(D_i)``, nor these data for every induced chamber of
the surface.  The scope of this program is exactly the three published groups.

For every element of each published stabilizer, reduction modulo 1009 gives
an exact upper bound on its rational fixed rank.  Only the identity and
rank-14 involutions survive above rank 10.  Pairwise modular row reduction
then proves that distinct surviving involutions have common fixed rank at
most 10.  Centralizers and conjugacy invariants reduce the surviving
involutions to seven H^13 casts.

For each cast the archived wall permutation is derived from the lattice
matrix, orbit-sum normals are formed in the fixed space, and two walls are
proved to be genuine facets by exact strict relative-interior witnesses.
Their normalized squared product is not a Coxeter value.  The local angle
failure is therefore decisive; no finite-volume test is run on a failed
Coxeter candidate.

NumPy/Numba are used to make the exhaustive finite-group pass practical.
The decisive implications are exact: integer matrix products, rank over a
finite field, and a final object-integer replay of every facet witness.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import sympy as sp

try:
    from numba import njit, prange
except ImportError:  # pragma: no cover - the exact fallback is simply slow
    def njit(*args: Any, **kwargs: Any):
        def decorate(function: Any) -> Any:
            return function
        return decorate

    prange = range


ROOT = Path(__file__).resolve().parents[1]
SOURCE_DIRECTORY = (
    ROOT / "data" / "external" / "katsura_kondo_shimada_char5"
)
BASIC_SOURCE = SOURCE_DIRECTORY / "compdataB.txt"
CHAMBER_SOURCE = SOURCE_DIRECTORY / "compdataChams.txt"
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "ss_k3_char5_d0_d2_high_rank_no_go.json"
)

RANK = 22
MODULUS = 1009
EXPECTED_HASHES = {
    "compdataB.txt": "956391d546707c958947d55c9bfc43ae1ae7952f73a01d4a4d1a0505169e5a1e",
    "compdataChams.txt": "5173c1c2974c9db1a78c3bf19f6cdc2ad4f08c10be3df9e99cfdd85c6dce7805",
}

# Zero-based selectors in the published factor lists.  The final artifact
# also prints one-based selectors and flattened positions for easy checking.
CLASS_REPRESENTATIVES = {
    0: ((0, 206),),
    1: ((0, 42), (1, 9)),
    2: ((0, 0), (0, 1), (0, 26), (12, 38)),
}


def _extract_assignment(text: str, name: str) -> str:
    """Return one top-level bracketed GAP assignment as Python text."""

    marker = f"{name}:="
    start = text.index(marker) + len(marker)
    while text[start].isspace():
        start += 1
    if text[start] != "[":
        raise ValueError(f"assignment {name} is not a bracketed list")
    depth = 0
    quote: str | None = None
    for end in range(start, len(text)):
        character = text[end]
        if quote is not None:
            if character == quote and text[end - 1] != "\\":
                quote = None
            continue
        if character in {'"', "'"}:
            quote = character
        elif character == "[":
            depth += 1
        elif character == "]":
            depth -= 1
            if depth == 0:
                return text[start : end + 1]
    raise ValueError(f"unterminated assignment {name}")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _integer_numpy(payload: Any) -> np.ndarray:
    return np.asarray(payload, dtype=np.int64)


@lru_cache(maxsize=4)
def parse_source(
    basic_path: Path = BASIC_SOURCE,
    chamber_path: Path = CHAMBER_SOURCE,
) -> dict[str, Any]:
    """Parse exactly the lattice, three wall stocks and group factors."""

    hashes = {
        basic_path.name: _sha256(basic_path),
        chamber_path.name: _sha256(chamber_path),
    }
    if basic_path == BASIC_SOURCE and chamber_path == CHAMBER_SOURCE:
        if hashes != EXPECTED_HASHES:
            raise ArithmeticError("the archived source hashes have changed")

    basic_text = basic_path.read_text(encoding="utf-8")
    chamber_text = chamber_path.read_text(encoding="utf-8")
    gram = sp.Matrix(ast.literal_eval(_extract_assignment(basic_text, "GramSX")))
    scaled_dual = 5 * gram.inv()
    if any(sp.Rational(value).q != 1 for value in scaled_dual):
        raise ArithmeticError("5*GramSX^-1 is not integral")

    chambers: list[dict[str, Any]] = []
    for index in range(3):
        wall_blocks_payload = ast.literal_eval(
            _extract_assignment(chamber_text, f"walls[{index}]")
        )
        wall_blocks = tuple(_integer_numpy(block) for block in wall_blocks_payload)
        walls = np.concatenate(wall_blocks, axis=0)
        s_factors = _integer_numpy(
            ast.literal_eval(_extract_assignment(chamber_text, f"SAutXD[{index}]"))
        )
        t_factors = _integer_numpy(
            ast.literal_eval(_extract_assignment(chamber_text, f"TAutXD[{index}]"))
        )
        weyl = _integer_numpy(
            ast.literal_eval(_extract_assignment(chamber_text, f"w[{index}]"))[:RANK]
        )
        polarization = _integer_numpy(
            ast.literal_eval(_extract_assignment(chamber_text, f"h[{index}]"))
        )
        if walls.shape[1] != RANK:
            raise ArithmeticError(f"D{index} wall width is not 22")
        if s_factors.shape[1:] != (RANK, RANK):
            raise ArithmeticError(f"D{index} S-factor shape is wrong")
        if t_factors.shape[1:] != (RANK, RANK):
            raise ArithmeticError(f"D{index} T-factor shape is wrong")
        chambers.append(
            {
                "index": index,
                "weyl_covector": weyl,
                "polarization": polarization,
                "wall_blocks": wall_blocks,
                "walls": walls,
                "s_factors": s_factors,
                "t_factors": t_factors,
            }
        )
    return {
        "gram_sympy": gram,
        "gram": _integer_numpy(gram.tolist()),
        "scaled_dual": _integer_numpy(scaled_dual.tolist()),
        "chambers": tuple(chambers),
        "hashes": hashes,
    }


@njit(cache=True)
def _inverse_mod(value: int, prime: int) -> int:
    result = 1
    base = value % prime
    exponent = prime - 2
    while exponent:
        if exponent & 1:
            result = (result * base) % prime
        base = (base * base) % prime
        exponent >>= 1
    return result


@njit(cache=True)
def _rank_mod(matrix: np.ndarray, prime: int) -> int:
    """Exact row rank over F_prime."""

    reduced = matrix.copy()
    rows, columns = reduced.shape
    rank = 0
    for column in range(columns):
        pivot = -1
        for row in range(rank, rows):
            if reduced[row, column] % prime:
                pivot = row
                break
        if pivot < 0:
            continue
        if pivot != rank:
            for position in range(column, columns):
                temporary = reduced[rank, position]
                reduced[rank, position] = reduced[pivot, position]
                reduced[pivot, position] = temporary
        inverse = _inverse_mod(reduced[rank, column], prime)
        for position in range(column, columns):
            reduced[rank, position] = (
                reduced[rank, position] * inverse
            ) % prime
        for row in range(rows):
            if row == rank:
                continue
            factor = reduced[row, column] % prime
            if factor:
                for position in range(column, columns):
                    reduced[row, position] = (
                        reduced[row, position]
                        - factor * reduced[rank, position]
                    ) % prime
        rank += 1
        if rank == columns:
            break
    return rank


@njit(cache=True, parallel=True)
def _fixed_rank_upper_bounds(products: np.ndarray, prime: int) -> np.ndarray:
    """Return 22-rank_Fp(g-I), an upper bound for dim_Q Fix(g)."""

    count = len(products)
    bounds = np.empty(count, dtype=np.int16)
    for index in prange(count):
        constraint = products[index].copy()
        for diagonal in range(RANK):
            constraint[diagonal, diagonal] -= 1
        bounds[index] = RANK - _rank_mod(constraint, prime)
    return bounds


@njit(cache=True, parallel=True)
def _pair_fixed_rank_upper_bounds(
    involutions: np.ndarray, prime: int
) -> np.ndarray:
    """Upper bounds for common fixed ranks of all unordered pairs."""

    count = len(involutions)
    bounds = np.empty(count * (count - 1) // 2, dtype=np.int16)
    pair_index = np.empty((len(bounds), 2), dtype=np.int32)
    cursor = 0
    for left in range(count):
        for right in range(left + 1, count):
            pair_index[cursor, 0] = left
            pair_index[cursor, 1] = right
            cursor += 1
    identity = np.eye(RANK, dtype=np.int64)
    for index in prange(len(bounds)):
        stacked = np.empty((2 * RANK, RANK), dtype=np.int64)
        left = pair_index[index, 0]
        right = pair_index[index, 1]
        stacked[:RANK] = involutions[left] - identity
        stacked[RANK:] = involutions[right] - identity
        bounds[index] = RANK - _rank_mod(stacked, prime)
    return bounds


def _group_chunks(chamber: dict[str, Any], chunk_size: int = 50):
    s_factors = chamber["s_factors"]
    t_factors = chamber["t_factors"]
    for start in range(0, len(s_factors), chunk_size):
        stop = min(start + chunk_size, len(s_factors))
        products = (
            s_factors[start:stop, None] @ t_factors[None]
        ).reshape(-1, RANK, RANK)
        yield start, products


def _matrix_at(chamber: dict[str, Any], selector: tuple[int, int]) -> np.ndarray:
    return chamber["s_factors"][selector[0]] @ chamber["t_factors"][selector[1]]


def _full_high_rank_census(chamber: dict[str, Any]) -> dict[str, Any]:
    """Exact modular screen of every published stabilizer element."""

    t_count = len(chamber["t_factors"])
    histogram: Counter[int] = Counter()
    survivors: list[tuple[int, int, int]] = []
    for s_start, products in _group_chunks(chamber):
        bounds = _fixed_rank_upper_bounds(products, MODULUS)
        histogram.update(map(int, bounds))
        for offset in np.flatnonzero(bounds > 10):
            s_index = s_start + int(offset) // t_count
            t_index = int(offset) % t_count
            survivors.append((s_index, t_index, int(bounds[offset])))

    identity = np.eye(RANK, dtype=np.int64)
    involution_selectors: list[tuple[int, int]] = []
    identity_selectors: list[tuple[int, int]] = []
    involutions: list[np.ndarray] = []
    for s_index, t_index, bound in survivors:
        matrix = _matrix_at(chamber, (s_index, t_index))
        if np.array_equal(matrix, identity):
            if bound != RANK:
                raise ArithmeticError("identity has wrong modular fixed rank")
            identity_selectors.append((s_index, t_index))
            continue
        if not np.array_equal(matrix @ matrix, identity) or int(np.trace(matrix)) != 6:
            raise ArithmeticError("an unclassified element survives above fixed rank 10")
        # For an involution in characteristic zero, trace 6 in rank 22 gives
        # dim Fix = (22+6)/2 = 14 exactly.
        if bound != 14:
            raise ArithmeticError("a trace-6 involution has wrong modular bound")
        involution_selectors.append((s_index, t_index))
        involutions.append(matrix)

    if len(identity_selectors) != 1:
        raise ArithmeticError("the factorization does not contain one identity")
    involution_array = np.asarray(involutions, dtype=np.int64)
    if len({matrix.tobytes() for matrix in involution_array}) != len(involution_array):
        raise ArithmeticError("the high-rank involution list contains duplicates")
    pair_bounds = _pair_fixed_rank_upper_bounds(involution_array, MODULUS)
    return {
        "modular_upper_bound_histogram": {
            str(key): value for key, value in sorted(histogram.items())
        },
        "survivor_count": len(survivors),
        "identity_selector": identity_selectors[0],
        "involution_selectors": tuple(involution_selectors),
        "involutions": involution_array,
        "pair_upper_bound_histogram": {
            str(key): value
            for key, value in sorted(Counter(map(int, pair_bounds)).items())
        },
        "maximum_distinct_pair_fixed_rank": int(pair_bounds.max()),
    }


def _centralizer_size(chamber: dict[str, Any], matrix: np.ndarray) -> int:
    count = 0
    for _, products in _group_chunks(chamber):
        left = products @ matrix
        right = matrix @ products
        count += int(np.all(left == right, axis=(1, 2)).sum())
    return count


def _wall_cycles(
    walls: np.ndarray, involution: np.ndarray
) -> tuple[tuple[tuple[int, ...], ...], np.ndarray]:
    """Derive the dual action and return wall cycles and orbit sums."""

    # The source matrices act on primal rows.  Wall rows are covectors, so
    # their right action is inverse-transpose.  Here g^-1=g.
    dual_action = involution.T
    lookup = {tuple(map(int, wall)): index for index, wall in enumerate(walls)}
    permutation: list[int] = []
    for wall in walls:
        image = tuple(map(int, wall @ dual_action))
        if image not in lookup:
            raise ArithmeticError("a stabilizer element does not permute the walls")
        permutation.append(lookup[image])
    unseen = set(range(len(walls)))
    cycles: list[tuple[int, ...]] = []
    normals: list[np.ndarray] = []
    while unseen:
        left = min(unseen)
        right = permutation[left]
        if permutation[right] != left:
            raise ArithmeticError("the derived wall action is not involutory")
        cycle = tuple(sorted({left, right}))
        unseen.difference_update(cycle)
        cycles.append(cycle)
        normals.append(walls[list(cycle)].sum(axis=0))
    return tuple(cycles), np.asarray(normals, dtype=np.int64)


def _fixed_counts_by_block(
    cycles: Iterable[tuple[int, ...]], block_sizes: tuple[int, ...]
) -> tuple[int, ...]:
    offsets = np.cumsum((0, *block_sizes))
    counts = []
    for block, size in enumerate(block_sizes):
        start, stop = int(offsets[block]), int(offsets[block + 1])
        counts.append(
            sum(
                len(cycle) == 1 and start <= cycle[0] < stop
                for cycle in cycles
            )
        )
    return tuple(counts)


def _is_forbidden_rational_coxeter_ratio(ratio: Fraction) -> bool:
    # If cos(pi/m)^2 is rational then it is one of these five values.
    return 0 < ratio < 1 and ratio not in {
        Fraction(1, 4),
        Fraction(1, 2),
        Fraction(3, 4),
    }


def _discover_strict_bad_pair(
    normals: np.ndarray,
    base: np.ndarray,
    scaled_dual: np.ndarray,
) -> tuple[int, int]:
    """Find an obstruction; final certification is object-integer exact."""

    squares = np.einsum("ij,jk,ik->i", normals, scaled_dual, normals)
    gram_rows = normals @ scaled_dual
    strict: list[int] = []
    for index, normal in enumerate(normals):
        square = int(squares[index])
        if square >= 0:
            continue
        product_with_base = int(base @ scaled_dual @ normal)
        witness = (-square) * base + product_with_base * normal
        witness_square = int(witness @ scaled_dual @ witness)
        pairings = gram_rows @ witness
        if witness_square <= 0 or pairings[index] != 0:
            continue
        if not np.all(np.delete(pairings, index) > 0):
            continue
        for other in strict:
            product = int(normal @ scaled_dual @ normals[other])
            ratio = Fraction(
                product * product,
                int(squares[index] * squares[other]),
            )
            if _is_forbidden_rational_coxeter_ratio(ratio):
                return other, index
        strict.append(index)
    raise ArithmeticError("no strict forbidden pair was found")


def _object_dot(left: np.ndarray, form: np.ndarray, right: np.ndarray) -> int:
    left_object = left.astype(object, copy=False)
    form_object = form.astype(object, copy=False)
    right_object = right.astype(object, copy=False)
    return int(left_object @ form_object @ right_object)


def _strict_facet_certificate(
    index: int,
    cycles: tuple[tuple[int, ...], ...],
    normals: np.ndarray,
    base: np.ndarray,
    scaled_dual: np.ndarray,
) -> dict[str, Any]:
    normal = normals[index]
    square = _object_dot(normal, scaled_dual, normal)
    product_with_base = _object_dot(base, scaled_dual, normal)
    witness = (-square) * base.astype(object) + product_with_base * normal.astype(object)
    witness = np.asarray(witness, dtype=object)
    witness_square = _object_dot(witness, scaled_dual, witness)
    pairings = np.asarray(
        normals.astype(object) @ scaled_dual.astype(object) @ witness,
        dtype=object,
    )
    if pairings[index] != 0:
        raise ArithmeticError("strict witness is not on its wall")
    other_pairings = np.delete(pairings, index)
    if witness_square <= 0 or any(int(value) <= 0 for value in other_pairings):
        raise ArithmeticError("strict witness is not in the relative interior")
    return {
        "orbit_index_zero_based": index,
        "orbit_index_one_based": index + 1,
        "source_wall_indices_one_based": [position + 1 for position in cycles[index]],
        "restricted_normal_covector": [int(value) for value in normal],
        "scaled_square": square,
        "base_product": product_with_base,
        "witness_covector": [int(value) for value in witness],
        "scaled_witness_square": witness_square,
        "minimum_other_wall_pairing": min(map(int, other_pairings)),
    }


def _bad_pair_certificate(
    cycles: tuple[tuple[int, ...], ...],
    normals: np.ndarray,
    base: np.ndarray,
    scaled_dual: np.ndarray,
) -> dict[str, Any]:
    left_index, right_index = _discover_strict_bad_pair(
        normals, base, scaled_dual
    )
    left = _strict_facet_certificate(
        left_index, cycles, normals, base, scaled_dual
    )
    right = _strict_facet_certificate(
        right_index, cycles, normals, base, scaled_dual
    )
    product = _object_dot(normals[left_index], scaled_dual, normals[right_index])
    ratio = Fraction(
        product * product,
        left["scaled_square"] * right["scaled_square"],
    )
    if not _is_forbidden_rational_coxeter_ratio(ratio):
        raise ArithmeticError("the selected pair is not a forbidden Coxeter angle")
    return {
        "left": left,
        "right": right,
        "scaled_product": product,
        "normalized_squared_product": str(ratio),
        "rational_coxeter_values_between_zero_and_one": ["1/4", "1/2", "3/4"],
        "decision": "not_coxeter",
    }


def _matrix_checks(
    matrix: np.ndarray,
    gram: np.ndarray,
    base: np.ndarray,
) -> None:
    identity = np.eye(RANK, dtype=np.int64)
    if not np.array_equal(matrix @ matrix, identity):
        raise ArithmeticError("class representative is not an involution")
    if int(np.trace(matrix)) != 6:
        raise ArithmeticError("class representative does not have fixed rank 14")
    if not np.array_equal(matrix @ gram @ matrix.T, gram):
        raise ArithmeticError("class representative is not a lattice isometry")
    if not np.array_equal(base @ matrix.T, base):
        raise ArithmeticError("class representative does not fix the Weyl covector")


def _selector_payload(
    selector: tuple[int, int], t_factor_count: int
) -> dict[str, int]:
    return {
        "s_index_zero_based": selector[0],
        "t_index_zero_based": selector[1],
        "s_index_one_based": selector[0] + 1,
        "t_index_one_based": selector[1] + 1,
        "flat_index_zero_based": selector[0] * t_factor_count + selector[1],
        "flat_index_one_based": selector[0] * t_factor_count + selector[1] + 1,
    }


def build_artifact(
    basic_path: Path = BASIC_SOURCE,
    chamber_path: Path = CHAMBER_SOURCE,
) -> dict[str, Any]:
    source = parse_source(basic_path, chamber_path)
    gram = source["gram"]
    scaled_dual = source["scaled_dual"]
    gram_sympy = source["gram_sympy"]
    if gram_sympy.det() != -25:
        raise ArithmeticError("GramSX does not have determinant -25")

    chamber_payloads: list[dict[str, Any]] = []
    total_classes = 0
    for chamber in source["chambers"]:
        index = chamber["index"]
        high = _full_high_rank_census(chamber)
        wall_blocks = chamber["wall_blocks"]
        block_sizes = tuple(len(block) for block in wall_blocks)
        group_order = len(chamber["s_factors"]) * len(chamber["t_factors"])
        high_selector_set = set(high["involution_selectors"])
        representatives = CLASS_REPRESENTATIVES[index]
        classes: list[dict[str, Any]] = []
        class_sizes: list[int] = []
        invariant_signatures: set[tuple[int, ...]] = set()
        for class_number, selector in enumerate(representatives, start=1):
            if selector not in high_selector_set:
                raise ArithmeticError("class representative was omitted by high-rank census")
            matrix = _matrix_at(chamber, selector)
            _matrix_checks(matrix, gram, chamber["weyl_covector"])
            cycles, normals = _wall_cycles(chamber["walls"], matrix)
            fixed_counts = _fixed_counts_by_block(cycles, block_sizes)
            if fixed_counts in invariant_signatures:
                raise ArithmeticError("class representatives lack distinct wall signatures")
            invariant_signatures.add(fixed_counts)
            centralizer = _centralizer_size(chamber, matrix)
            if group_order % centralizer:
                raise ArithmeticError("centralizer does not divide group order")
            class_size = group_order // centralizer
            class_sizes.append(class_size)
            bad_pair = _bad_pair_certificate(
                cycles,
                normals,
                chamber["weyl_covector"],
                scaled_dual,
            )
            classes.append(
                {
                    "class_number": class_number,
                    "representative_selector": _selector_payload(
                        selector, len(chamber["t_factors"])
                    ),
                    "matrix_trace": int(np.trace(matrix)),
                    "matrix_order": 2,
                    "fixed_rank": 14,
                    "hyperbolic_dimension": 13,
                    "centralizer_order": centralizer,
                    "conjugacy_class_size": class_size,
                    "fixed_wall_counts_by_published_block": list(fixed_counts),
                    "wall_orbit_count": len(cycles),
                    "singleton_wall_orbit_count": sum(len(cycle) == 1 for cycle in cycles),
                    "two_cycle_wall_orbit_count": sum(len(cycle) == 2 for cycle in cycles),
                    "strict_bad_pair": bad_pair,
                }
            )
        if sum(class_sizes) != len(high["involutions"]):
            raise ArithmeticError("listed conjugacy classes do not exhaust high involutions")
        total_classes += len(classes)
        chamber_payloads.append(
            {
                "chamber": f"D{index}",
                "published_wall_block_sizes": list(block_sizes),
                "published_wall_count": len(chamber["walls"]),
                "stabilizer_factor_lengths": [
                    len(chamber["s_factors"]),
                    len(chamber["t_factors"]),
                ],
                "stabilizer_order": group_order,
                "exhaustive_fixed_rank_screen": {
                    "field": f"F_{MODULUS}",
                    "logic": (
                        "rank_Q(g-I) >= rank_Fp(g-I), so 22-rank_Fp is an "
                        "exact upper bound for the rational fixed rank"
                    ),
                    "fixed_rank_upper_bound_histogram": high[
                        "modular_upper_bound_histogram"
                    ],
                    "survivor_count_above_rank_10_including_identity": high[
                        "survivor_count"
                    ],
                    "rank_14_involution_count": len(high["involutions"]),
                    "identity_count": 1,
                },
                "distinct_rank_14_pair_screen": {
                    "fixed_rank_upper_bound_histogram": high[
                        "pair_upper_bound_histogram"
                    ],
                    "maximum_common_fixed_rank": high[
                        "maximum_distinct_pair_fixed_rank"
                    ],
                },
                "high_rank_conjugacy_classes": classes,
                "decision": "all_fixed_slices_above_H9_are_non_coxeter",
            }
        )

    return {
        "schema": "ss-k3-char5-d0-d2-high-rank-no-go-v1",
        "claim": (
            "Every nontrivial subgroup-fixed slice of hyperbolic dimension "
            ">9 inside the three published period-preserving groups "
            "Aut_X(D0), Aut_X(D1), Aut_X(D2) is one of seven H13 involution "
            "casts, and every cast fails the exact Coxeter-angle gate."
        ),
        "scope_warning": (
            "This is not a theorem about the full lattice stabilizers Aut(D_i), "
            "nor about every induced chamber of the surface. The source lists "
            "Aut_X(D_i)=Aut(D_i) intersect G_X; it does not publish Aut(D_i). "
            "It also reports many further chamber classes without complete "
            "wall/stabilizer matrices."
        ),
        "source": {
            "paper": (
                "T. Katsura, S. Kondo, I. Shimada, On the supersingular K3 "
                "surface in characteristic 5 with Artin invariant 1, "
                "Michigan Math. J. 63 (2014), 803-844"
            ),
            "doi": "10.1307/mmj/1417799227",
            "arxiv": "1312.0687",
            "paper_url": (
                "https://home.hiroshima-u.ac.jp/ichiro-shimada/preprints/"
                "KKSssK3char5/KKSssK3char5.pdf"
            ),
            "data_urls": {
                "compdataB.txt": (
                    "https://home.hiroshima-u.ac.jp/ichiro-shimada/preprints/"
                    "KKSssK3char5/compdatassK3char5/compdataB.txt"
                ),
                "compdataChams.txt": (
                    "https://home.hiroshima-u.ac.jp/ichiro-shimada/preprints/"
                    "KKSssK3char5/compdatassK3char5/compdataChams.txt"
                ),
            },
            "sha256": source["hashes"],
            "lineage": {
                "surface": "supersingular K3 in characteristic 5, Artin invariant 1",
                "lattice_assignment": "GramSX",
                "chambers": ["D0", "D1", "D2"],
                "wall_assignments": ["walls[0]", "walls[1]", "walls[2]"],
                "stabilizer_factor_assignments": ["SAutXD[i]", "TAutXD[i]"],
                "matrix_scope": (
                    "the complete published period-preserving groups Aut_X(D_i), "
                    "not the unreported full lattice groups Aut(D_i); within each "
                    "published group no further element or subgroup filter is applied"
                ),
            },
        },
        "lattice": {
            "rank": RANK,
            "determinant": int(gram_sympy.det()),
            "signature": [1, 21],
            "dual_form_used": "5*GramSX^-1",
            "dual_form_integral": True,
        },
        "method": {
            "dimension_threshold": "hyperbolic dimension > 9, equivalently fixed rank > 10",
            "group_completeness": (
                "every published product SAutXD[i][s]*TAutXD[i][t] is screened; "
                "this exhausts Aut_X(D_i), not necessarily Aut(D_i)"
            ),
            "noncyclic_exclusion": (
                "every nonidentity element of a subgroup above the threshold must "
                "be rank 14; two distinct such involutions have common fixed rank <=10"
            ),
            "facet_gate": (
                "a wall is accepted only with a timelike point on it pairing "
                "strictly positively with every other restricted inequality"
            ),
            "volume_policy": (
                "Vinberg finite-volume/compactness tests are run only after the "
                "Coxeter-angle gate; all seven casts fail before that stage"
            ),
        },
        "chambers": chamber_payloads,
        "summary": {
            "high_rank_conjugacy_class_count": total_classes,
            "hyperbolic_dimensions": [13],
            "positive_examples": 0,
            "non_coxeter_casts": total_classes,
            "finite_volume_tests_run": 0,
            "decision": "published_AutX_D0_D2_groups_closed_negative_above_H9",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--basic-source", type=Path, default=BASIC_SOURCE)
    parser.add_argument("--chamber-source", type=Path, default=CHAMBER_SOURCE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    artifact = build_artifact(args.basic_source, args.chamber_source)
    if args.verify:
        archived = json.loads(args.output.read_text(encoding="utf-8"))
        if archived != artifact:
            raise SystemExit("archived certificate differs from exact replay")
        print(f"verified {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
        print(f"wrote {args.output}")
        print(json.dumps(artifact["summary"], indent=2))


if __name__ == "__main__":
    main()
