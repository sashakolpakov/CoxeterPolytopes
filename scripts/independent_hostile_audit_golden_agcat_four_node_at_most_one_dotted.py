#!/usr/bin/env python3
"""No-import hostile replay of the golden four-node AGCat tranche.

This script intentionally does not import ``golden_four_node_packets``, its
certifier, or the golden triangle packet implementation.  It reconstructs
the two-place golden edge alphabet, element squareclasses, S4 orbit
quotients, simultaneous graph-cycle gates, Gram determinants, exact inertia
signs, stream hashes, and the rank-23--30 terminal-dotted path family.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from math import isqrt
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OWNER_JSON = ROOT / "data/enumerations/golden_agcat_four_node_at_most_one_dotted.json"
OUTPUT = ROOT / "data/enumerations/golden_agcat_four_node_at_most_one_dotted_hostile_audit.json"
REPORT = ROOT / "results/golden_agcat_four_node_at_most_one_dotted_hostile_audit.md"

Pair = tuple[int, int]
Edge = tuple[int, int]
ElementClass = tuple[tuple[tuple[int, int], ...], int]

EDGES: tuple[Edge, ...] = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
EDGE_SLOT = {edge: slot for slot, edge in enumerate(EDGES)}
FINITE: tuple[Pair, ...] = ((0, 0), (1, 0), (2, 0), (1, 1), (3, 0), (2, 1))
TRIVIAL_CLASS: ElementClass = ((), 1)


def pair_add(*values: Pair) -> Pair:
    return sum(value[0] for value in values), sum(value[1] for value in values)


def pair_scale(multiplier: int, value: Pair) -> Pair:
    return multiplier * value[0], multiplier * value[1]


def pair_mul(left: Pair, right: Pair) -> Pair:
    a, b = left
    c, d = right
    return a * c + b * d, a * d + b * c + b * d


def pair_conjugate(value: Pair) -> Pair:
    a, b = value
    return a + b, -b


def pair_norm(value: Pair) -> int:
    a, b = value
    return a * a + a * b - b * b


def pair_sign(value: Pair, place: int = 0) -> int:
    if place not in (0, 1):
        raise ValueError("place must be 0 or 1")
    a, b = pair_conjugate(value) if place else value
    rational = 2 * a + b
    if not b:
        return (rational > 0) - (rational < 0)
    if not rational:
        return (b > 0) - (b < 0)
    rational_sign = (rational > 0) - (rational < 0)
    radical_sign = (b > 0) - (b < 0)
    if rational_sign == radical_sign:
        return rational_sign
    comparison = rational * rational - 5 * b * b
    return rational_sign if comparison > 0 else radical_sign


def floor_pair(value: Pair, place: int = 0) -> int:
    a, b = pair_conjugate(value) if place else value
    if not b:
        return a
    magnitude = abs(b)
    floor_magnitude_phi = (magnitude + isqrt(5 * magnitude * magnitude)) // 2
    return a + floor_magnitude_phi if b > 0 else a - floor_magnitude_phi - 1


def ceil_pair(value: Pair, place: int = 0) -> int:
    floor = floor_pair(value, place)
    _a, b = pair_conjugate(value) if place else value
    return floor if not b else floor + 1


def independent_edge_alphabet() -> tuple[Pair, ...]:
    """Rebuild the A=62 two-place alphabet without owner imports."""

    finite = set(FINITE)
    output: list[Pair] = []
    # If x=a+b phi lies in [0,3844) x [0,4), then
    # -4 < b*sqrt(5) < 3844; [-2,1922] is a proved generous interval.
    for b in range(-2, 1923):
        lower = max(ceil_pair((0, -b), 0), ceil_pair((0, -b), 1))
        upper = min(
            ceil_pair((3844, -b), 0) - 1,
            ceil_pair((4, -b), 1) - 1,
        )
        for a in range(lower, upper + 1):
            value = (a, b)
            if not (
                pair_sign(value, 0) >= 0
                and pair_sign(pair_add(value, (-3844, 0)), 0) < 0
                and pair_sign(value, 1) >= 0
                and pair_sign(pair_add(value, (-4, 0)), 1) < 0
            ):
                raise AssertionError("independent alphabet interval leaked")
            if value in finite or pair_sign(pair_add(value, (-4, 0)), 0) > 0:
                output.append(value)
    if len(output) != len(set(output)):
        raise AssertionError("independent edge alphabet has duplicates")
    return tuple(output)


def rational_prime_factors(value: int) -> tuple[tuple[int, int], ...]:
    if value <= 0:
        raise ValueError("factorization input must be positive")
    output = []
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor:
            divisor = 3 if divisor == 2 else divisor + 2
            continue
        exponent = 0
        while value % divisor == 0:
            value //= divisor
            exponent += 1
        output.append((divisor, exponent))
        divisor = 3 if divisor == 2 else divisor + 2
    if value > 1:
        output.append((value, 1))
    return tuple(output)


def roots_mod_prime(prime: int) -> tuple[int, ...]:
    return tuple(r for r in range(prime) if (r * r - r - 1) % prime == 0)


def hensel_root(residue: int, prime: int, exponent: int) -> int:
    modulus = prime
    for _ in range(1, exponent):
        quotient = (residue * residue - residue - 1) // modulus
        correction = -quotient * pow((2 * residue - 1) % prime, -1, prime)
        residue += (correction % prime) * modulus
        modulus *= prime
    if (residue * residue - residue - 1) % modulus:
        raise AssertionError("independent Hensel lift failed")
    return residue


def valuation_integer(value: int, prime: int, cap: int) -> int:
    if not value:
        return cap
    exponent = 0
    while exponent < cap and value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def divide_ramified_five(value: Pair) -> Pair | None:
    # (2+phi)(3-phi)=5.
    numerator = pair_mul(value, (3, -1))
    if numerator[0] % 5 or numerator[1] % 5:
        return None
    return numerator[0] // 5, numerator[1] // 5


def element_squareclass(value: Pair) -> ElementClass:
    if value == (0, 0) or pair_sign(value, 0) <= 0:
        raise ValueError("squareclass expects an identity-positive nonzero value")
    norm = abs(pair_norm(value))
    a, b = value
    parity: list[tuple[int, int]] = []
    for prime, norm_exponent in rational_prime_factors(norm):
        if prime == 5:
            quotient = value
            exponent = 0
            while (next_quotient := divide_ramified_five(quotient)) is not None:
                quotient = next_quotient
                exponent += 1
            if exponent != norm_exponent:
                raise AssertionError("ramified valuation/norm mismatch")
            if exponent % 2:
                parity.append((5, 0))
            continue
        roots = roots_mod_prime(prime)
        if not roots:
            exponent = min(
                valuation_integer(a, prime, norm_exponent),
                valuation_integer(b, prime, norm_exponent),
            )
            if 2 * exponent != norm_exponent:
                raise AssertionError("inert valuation/norm mismatch")
            if exponent % 2:
                parity.append((prime, 0))
            continue
        valuations = []
        for residue in roots:
            lifted = hensel_root(residue, prime, norm_exponent)
            exponent = valuation_integer(a + b * lifted, prime, norm_exponent)
            valuations.append(exponent)
            if exponent % 2:
                parity.append((prime, residue))
        if sum(valuations) != norm_exponent:
            raise AssertionError("split valuation/norm mismatch")
    # The selected sign is positive throughout.  Since O_K^*=+/-phi^Z,
    # the conjugate sign is the remaining unit-squareclass invariant.
    return tuple(sorted(parity)), pair_sign(value, 1)


def class_xor(left: ElementClass, right: ElementClass) -> ElementClass:
    return tuple(sorted(set(left[0]) ^ set(right[0]))), left[1] * right[1]


def xor_all(classes: tuple[ElementClass, ...]) -> ElementClass:
    output = TRIVIAL_CLASS
    for value in classes:
        output = class_xor(output, value)
    return output


@lru_cache(maxsize=None)
def square_root(value: Pair) -> Pair | None:
    """Solve A+B phi=(m+n phi)^2 by the independent n^2 quadratic."""

    if pair_sign(value, 0) <= 0 or pair_sign(value, 1) <= 0:
        return None
    a, b = value
    linear = 4 * a + 2 * b
    discriminant = linear * linear - 20 * b * b
    if discriminant < 0:
        return None
    delta = isqrt(discriminant)
    if delta * delta != discriminant:
        return None
    roots = set()
    for numerator in (linear - delta, linear + delta):
        if numerator < 0 or numerator % 10:
            continue
        n2 = numerator // 10
        n_abs = isqrt(n2)
        if n_abs * n_abs != n2:
            continue
        for n in ({n_abs, -n_abs} if n_abs else {0}):
            if n:
                m_num = b - n * n
                if m_num % (2 * n):
                    continue
                candidate = (m_num // (2 * n), n)
            else:
                m_abs = isqrt(a)
                if b or m_abs * m_abs != a:
                    continue
                candidate = (m_abs, 0)
            if pair_mul(candidate, candidate) == value and pair_sign(candidate, 0) > 0:
                roots.add(candidate)
    return next(iter(roots)) if len(roots) == 1 else None


def edge_slot(first: int, second: int) -> int:
    return EDGE_SLOT[tuple(sorted((first, second)))]


def permute_pattern(pattern: tuple[int, ...], permutation: tuple[int, ...]) -> tuple[int, ...]:
    output = [0] * 6
    for slot, (first, second) in enumerate(EDGES):
        output[edge_slot(permutation[first], permutation[second])] = pattern[slot]
    return tuple(output)


def connected(pattern: tuple[int, ...]) -> bool:
    reached = {0}
    for _ in range(4):
        for slot, (first, second) in enumerate(EDGES):
            if pattern[slot] and (first in reached or second in reached):
                reached.update((first, second))
    return len(reached) == 4


def finite_patterns() -> tuple[tuple[int, ...], ...]:
    group = tuple(permutations(range(4)))
    return tuple(
        pattern
        for pattern in product(range(6), repeat=6)
        if connected(pattern)
        and pattern == min(permute_pattern(pattern, permutation) for permutation in group)
    )


def one_dotted_patterns() -> tuple[tuple[int, ...], ...]:
    stabilizer = ((0, 1, 2, 3), (1, 0, 2, 3), (0, 1, 3, 2), (1, 0, 3, 2))
    output = []
    for complement in product(range(6), repeat=5):
        pattern = (6, *complement)
        if connected(pattern) and pattern == min(
            permute_pattern(pattern, permutation) for permutation in stabilizer
        ):
            output.append(pattern)
    return tuple(output)


def graph_cycles(pattern: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    cycles = []
    for vertices in combinations(range(4), 3):
        slots = tuple(
            edge_slot(vertices[index], vertices[(index + 1) % 3])
            for index in range(3)
        )
        if all(pattern[slot] for slot in slots):
            cycles.append(tuple(sorted(slots)))
    for order in permutations((1, 2, 3)):
        if order[0] > order[-1]:
            continue
        vertices = (0, *order)
        slots = tuple(
            edge_slot(vertices[index], vertices[(index + 1) % 4])
            for index in range(4)
        )
        if all(pattern[slot] for slot in slots):
            cycles.append(tuple(sorted(slots)))
    return tuple(sorted(set(cycles), key=lambda cycle: (len(cycle), cycle)))


def permutation_sign(permutation: tuple[int, ...]) -> int:
    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(len(permutation))
        for j in range(i + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def permutation_cycles(permutation: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    unseen = set(range(len(permutation)))
    output = []
    while unseen:
        start = min(unseen)
        cycle = [start]
        unseen.remove(start)
        while permutation[cycle[-1]] != start:
            cycle.append(permutation[cycle[-1]])
            unseen.remove(cycle[-1])
        output.append(tuple(cycle))
    return tuple(output)


def gram_determinant(edges: tuple[Pair, ...], vertices: tuple[int, ...]) -> Pair:
    """Independent permutation expansion of the normalized Gram determinant."""

    total = (0, 0)
    for local_permutation in permutations(range(len(vertices))):
        coefficient = permutation_sign(local_permutation)
        term = (1, 0)
        for local_cycle in permutation_cycles(local_permutation):
            if len(local_cycle) == 1:
                term = pair_scale(2, term)
                continue
            global_cycle = tuple(vertices[index] for index in local_cycle)
            cycle_edges = tuple(
                edges[edge_slot(global_cycle[index], global_cycle[(index + 1) % len(global_cycle)])]
                for index in range(len(global_cycle))
            )
            if any(value == (0, 0) for value in cycle_edges):
                term = (0, 0)
                break
            if len(local_cycle) == 2:
                term = pair_mul(term, cycle_edges[0])
            else:
                product_value = (1, 0)
                for value in cycle_edges:
                    product_value = pair_mul(product_value, value)
                root = square_root(product_value)
                if root is None:
                    raise AssertionError("cycle-compatible packet lost a cycle root")
                term = pair_mul(term, pair_scale((-1) ** len(local_cycle), root))
        total = pair_add(total, pair_scale(coefficient, term))
    return total


def topology(pattern: tuple[int, ...]) -> tuple[str, int, tuple[int, ...], int]:
    degrees = [0] * 4
    edge_count = 0
    for slot, (first, second) in enumerate(EDGES):
        if pattern[slot]:
            edge_count += 1
            degrees[first] += 1
            degrees[second] += 1
    sequence = tuple(sorted(degrees))
    names = {
        (3, (1, 1, 1, 3)): "star_K1_3",
        (3, (1, 1, 2, 2)): "path_P4",
        (4, (1, 2, 2, 3)): "triangle_with_pendant",
        (4, (2, 2, 2, 2)): "cycle_C4",
        (5, (2, 2, 3, 3)): "diamond_K4_minus_edge",
        (6, (3, 3, 3, 3)): "complete_K4",
    }
    return names[(edge_count, sequence)], edge_count, sequence, edge_count - 3


def record(edges: tuple[Pair, ...], leading: Pair, full: Pair) -> bytes:
    return (
        ";".join(f"{a},{b}" for a, b in (*edges, leading, full)) + "\n"
    ).encode()


def run_replay() -> dict[str, object]:
    alphabet = independent_edge_alphabet()
    dotted_b_major = tuple(
        value for value in alphabet if pair_sign(pair_add(value, (-4, 0)), 0) > 0
    )
    dotted_lex = tuple(sorted(dotted_b_major))
    if len(alphabet) != 6_875 or len(dotted_b_major) != 6_869:
        raise AssertionError("independent golden alphabet counts changed")

    finite_classes = tuple(
        None if value == (0, 0) else element_squareclass(value) for value in FINITE
    )
    dotted_by_class: dict[ElementClass, list[Pair]] = defaultdict(list)
    for value in dotted_b_major:
        dotted_by_class[element_squareclass(value)].append(value)
    for values in dotted_by_class.values():
        values.sort()

    finite_orbits = finite_patterns()
    dotted_orbits = one_dotted_patterns()
    strata: Counter[tuple[object, ...]] = Counter()
    topology_counts: Counter[tuple[str, str]] = Counter()
    topology_retained: Counter[tuple[str, str]] = Counter()
    counts: Counter[str] = Counter()
    digests = {name: sha256() for name in ("finite_all", "finite_kept", "dotted_all", "dotted_kept")}

    for pattern in finite_orbits:
        cycles = graph_cycles(pattern)
        if any(
            xor_all(tuple(finite_classes[pattern[slot]] for slot in cycle)) != TRIVIAL_CLASS
            for cycle in cycles
        ):
            continue
        edges = tuple(FINITE[index] for index in pattern)
        # Direct multiplication-back checks every accepted graph cycle.
        for cycle in cycles:
            value = (1, 0)
            for slot in cycle:
                value = pair_mul(value, edges[slot])
            root = square_root(value)
            if root is None or pair_mul(root, root) != value:
                raise AssertionError("finite compatible cycle failed direct square replay")
        leading = gram_determinant(edges, (0, 1, 2))
        full = gram_determinant(edges, (0, 1, 2, 3))
        signs = (pair_sign(leading, 0), pair_sign(full, 0), pair_sign(leading, 1), pair_sign(full, 1))
        sl, sf, cl, cf = signs
        conjugate_pd = cl > 0 and cf > 0
        selected_pd = sl > 0 and sf > 0
        selected_lorentzian = sf < 0
        retained = conjugate_pd and (selected_pd or selected_lorentzian)
        counts["finite_compatible"] += 1
        counts["finite_conjugate_pd"] += conjugate_pd
        counts["finite_retained"] += retained
        counts["finite_positive"] += retained and selected_pd
        counts["finite_lorentzian"] += retained and selected_lorentzian
        counts["finite_selected_leading_zero"] += sl == 0
        counts["finite_selected_full_zero"] += sf == 0
        counts["finite_conjugate_leading_zero"] += cl == 0
        counts["finite_conjugate_full_zero"] += cf == 0
        name, edge_count, degrees, cycle_rank = topology(pattern)
        key = ("finite_only", name)
        topology_counts[key] += 1
        topology_retained[key] += retained
        strata[("finite_only", name, edge_count, degrees, cycle_rank, *signs)] += 1
        encoded = record(edges, leading, full)
        digests["finite_all"].update(encoded)
        if retained:
            digests["finite_kept"].update(encoded)

    for pattern in dotted_orbits:
        cycles = graph_cycles(pattern)
        requirements: list[ElementClass] = []
        valid = True
        for cycle in cycles:
            finite_cycle = tuple(
                finite_classes[pattern[slot]] for slot in cycle if slot != 0
            )
            value_class = xor_all(finite_cycle)
            if 0 in cycle:
                requirements.append(value_class)
            elif value_class != TRIVIAL_CLASS:
                valid = False
                break
        if requirements and any(value != requirements[0] for value in requirements[1:]):
            valid = False
        choices = (
            tuple(dotted_by_class.get(requirements[0], ()))
            if valid and requirements
            else dotted_lex if valid else ()
        )
        name, edge_count, degrees, cycle_rank = topology(pattern)
        for dotted in choices:
            edges = (dotted,) + tuple(FINITE[index] for index in pattern[1:])
            for cycle in cycles:
                value = (1, 0)
                for slot in cycle:
                    value = pair_mul(value, edges[slot])
                root = square_root(value)
                if root is None or pair_mul(root, root) != value:
                    raise AssertionError("dotted compatible cycle failed direct square replay")
            leading = gram_determinant(edges, (0, 1, 2))
            full = gram_determinant(edges, (0, 1, 2, 3))
            signs = (pair_sign(leading, 0), pair_sign(full, 0), pair_sign(leading, 1), pair_sign(full, 1))
            sl, sf, cl, cf = signs
            conjugate_pd = cl > 0 and cf > 0
            selected_lorentzian = sl < 0 and sf < 0
            retained = conjugate_pd and selected_lorentzian
            counts["dotted_compatible"] += 1
            counts["dotted_conjugate_pd"] += conjugate_pd
            counts["dotted_retained"] += retained
            counts["dotted_lorentzian"] += retained and selected_lorentzian
            counts["dotted_nonlorentzian"] += not selected_lorentzian
            counts["dotted_selected_leading_zero"] += sl == 0
            counts["dotted_selected_full_zero"] += sf == 0
            counts["dotted_conjugate_leading_zero"] += cl == 0
            counts["dotted_conjugate_full_zero"] += cf == 0
            key = ("exactly_one_dotted", name)
            topology_counts[key] += 1
            topology_retained[key] += retained
            strata[("exactly_one_dotted", name, edge_count, degrees, cycle_rank, *signs)] += 1
            encoded = record(edges, leading, full)
            digests["dotted_all"].update(encoded)
            if retained:
                digests["dotted_kept"].update(encoded)

    path_rows = []
    for rank in range(23, 31):
        admissible = tuple(
            value
            for value in dotted_b_major
            if pair_sign(pair_add(pair_scale(rank - 1, value), (-2 * rank, 0)), 1) < 0
        )
        digest = sha256()
        for a, b in admissible:
            digest.update(f"{a},{b}\n".encode())
        sample = (2, 2)
        sample_det = pair_add((2 * rank, 0), pair_scale(-(rank - 1), sample))
        if sample not in admissible:
            raise AssertionError("uniform path sample disappeared")
        path_rows.append(
            {
                "native_rank": rank,
                "hyperbolic_dimension": rank - 1,
                "admissible_dotted_value_count": len(admissible),
                "admissible_dotted_stream_sha256": digest.hexdigest(),
                "sample_dotted_product": list(sample),
                "sample_full_determinant": list(sample_det),
            }
        )

    strata_rows = [
        {
            "tranche": key[0],
            "topology": key[1],
            "nonzero_edge_count": key[2],
            "degree_sequence": list(key[3]),
            "cycle_rank": key[4],
            "selected_leading_triple_sign": key[5],
            "selected_full_determinant_sign": key[6],
            "conjugate_leading_triple_sign": key[7],
            "conjugate_full_determinant_sign": key[8],
            "count": count,
        }
        for key, count in sorted(strata.items())
    ]
    topology_rows = [
        {
            "tranche": tranche,
            "topology": name,
            "cycle_compatible": count,
            "two_place_retained": topology_retained[(tranche, name)],
        }
        for (tranche, name), count in sorted(topology_counts.items())
    ]
    return {
        "alphabet": {
            "full_count": len(alphabet),
            "finite_count": len(FINITE),
            "dotted_count": len(dotted_b_major),
            "dotted_b_major_sha256": sha256(
                b"".join(f"{a},{b}\n".encode() for a, b in dotted_b_major)
            ).hexdigest(),
        },
        "canonical_orbits": {
            "finite_only_connected": len(finite_orbits),
            "one_dotted_complements": len(dotted_orbits),
            "one_dotted_expanded": len(dotted_orbits) * len(dotted_b_major),
        },
        "counts": dict(sorted(counts.items())),
        "topology_counts": topology_rows,
        "strata": strata_rows,
        "stream_digests": {
            "finite_only_cycle_compatible_sha256": digests["finite_all"].hexdigest(),
            "finite_only_two_place_retained_sha256": digests["finite_kept"].hexdigest(),
            "one_dotted_cycle_compatible_sha256": digests["dotted_all"].hexdigest(),
            "one_dotted_two_place_retained_sha256": digests["dotted_kept"].hexdigest(),
        },
        "global_path_rows": path_rows,
        "mathematical_audit": {
            "cycle_squareclasses_sufficient": True,
            "reason": (
                "tree gauge fixes directed entries; every chord fundamental-cycle "
                "product is the unique acute signed golden root. Direct square and "
                "multiplication-back checks were repeated on every retained graph cycle."
            ),
            "determinant_implementation": (
                "fresh 24-permutation expansion, not the owner's closed formula"
            ),
            "conjugate_Sylvester_valid": True,
            "finite_selected_inertia_logic": (
                "all selected 2x2 finite minors are positive; full determinant <0 "
                "and interlacing force exactly one negative eigenvalue, while leading3>0 "
                "and full>0 give positive definiteness"
            ),
            "one_dotted_selected_inertia_logic": (
                "the leading dotted 2x2 minor is negative; leading3<0 and full<0 "
                "make every later LDL pivot positive, hence index one"
            ),
            "path_cocompact_host_logic": (
                "the integral Cartan gauge has symmetrizer diag(c,1,...,1); conjugate "
                "positive definiteness makes the K-form anisotropic, so the full "
                "arithmetic orthogonal lattice is uniform. Reflectivity remains open."
            ),
            "local_scope_only": True,
        },
    }


def compare_owner(replay: dict[str, object], owner: dict[str, object]) -> None:
    counts = replay["counts"]
    expected = {
        "finite_compatible": owner["cycle_compatibility"]["finite_only_compatible_orbits"],
        "dotted_compatible": owner["cycle_compatibility"]["one_dotted_compatible_orbits"],
        "finite_conjugate_pd": owner["two_place_gram_gate"]["finite_only_conjugate_positive"],
        "dotted_conjugate_pd": owner["two_place_gram_gate"]["one_dotted_conjugate_positive"],
        "finite_retained": owner["two_place_gram_gate"]["finite_only_retained"],
        "dotted_retained": owner["two_place_gram_gate"]["one_dotted_retained"],
        "finite_positive": owner["two_place_gram_gate"]["finite_only_retained_selected_positive"],
        "finite_lorentzian": owner["two_place_gram_gate"]["finite_only_retained_selected_lorentzian"],
        "dotted_lorentzian": owner["two_place_gram_gate"]["one_dotted_retained_selected_lorentzian"],
    }
    for key, value in expected.items():
        if counts[key] != value:
            raise AssertionError(f"owner count mismatch {key}: {counts[key]} != {value}")
    degeneracy = owner["two_place_gram_gate"]["zero_determinant_counts"]
    split = {
        "finite_selected_leading_zero": degeneracy["finite_only_selected_leading_triple"],
        "finite_selected_full_zero": degeneracy["finite_only_selected_full_four_node"],
        "finite_conjugate_leading_zero": degeneracy["finite_only_conjugate_leading_triple"],
        "finite_conjugate_full_zero": degeneracy["finite_only_conjugate_full_four_node"],
        "dotted_selected_leading_zero": degeneracy["one_dotted_selected_leading_triple"],
        "dotted_selected_full_zero": degeneracy["one_dotted_selected_full_four_node"],
        "dotted_conjugate_leading_zero": degeneracy["one_dotted_conjugate_leading_triple"],
        "dotted_conjugate_full_zero": degeneracy["one_dotted_conjugate_full_four_node"],
    }
    for key, value in split.items():
        if counts[key] != value:
            raise AssertionError(f"owner degeneracy mismatch {key}")
    if replay["topology_counts"] != owner["topology_counts"]:
        raise AssertionError("owner topology aggregation differs")
    if replay["strata"] != owner["strata"]:
        raise AssertionError("owner exact sign strata differ")
    if replay["stream_digests"] != {
        key: owner["stream_digests"][key] for key in replay["stream_digests"]
    }:
        raise AssertionError("owner stream digests differ")
    if replay["global_path_rows"] != owner["global_one_dotted_label3_path_family"]["rows"]:
        raise AssertionError("owner global path rows differ")


def build_payload() -> dict[str, object]:
    owner_bytes = OWNER_JSON.read_bytes()
    owner = json.loads(owner_bytes)
    replay = run_replay()
    compare_owner(replay, owner)
    return {
        "schema": "golden-agcat-four-node-at-most-one-dotted-hostile-audit-0.1.0",
        "status": "PASS",
        "independence": {
            "owner_module_imported": False,
            "owner_certifier_imported": False,
            "triangle_helper_imported": False,
            "edge_alphabet_reconstructed": True,
            "squareclasses_reconstructed": True,
            "orbit_actions_reconstructed": True,
            "determinant_recomputed_by_permutation_expansion": True,
        },
        "owner_artifact": {
            "path": str(OWNER_JSON.relative_to(ROOT)),
            "sha256": sha256(owner_bytes).hexdigest(),
        },
        "replay": replay,
        "defect_history": {
            "found": True,
            "severity": "terminology_and_schema_only; retained sets unaffected",
            "original_problem": (
                "the owner called the union of 21 zero leading-three determinants and "
                "3 zero full-four determinants '24 degenerate packets'"
            ),
            "repair_verified": (
                "the frozen owner schema now reports the disjoint 21/3 split at both "
                "places; the one-dotted split is 0/0"
            ),
        },
        "verdict": {
            "four_node_local_tranche": "INDEPENDENT HOSTILE PASS",
            "global_path_family": "INDEPENDENT HOSTILE PASS",
            "new_compact_chamber": False,
            "golden_H22_H29_nonexistence": "OPEN",
        },
    }


def encoded_payload() -> bytes:
    return (json.dumps(build_payload(), indent=2, sort_keys=True) + "\n").encode()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    encoded = encoded_payload()
    if args.verify:
        if not OUTPUT.exists() or OUTPUT.read_bytes() != encoded:
            raise SystemExit("FAIL: hostile four-node audit differs from frozen replay")
        status = "PASS"
    else:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_bytes(encoded)
        status = "WROTE"
    payload = json.loads(encoded)
    print(
        json.dumps(
            {
                "status": status,
                "output": str(OUTPUT),
                "sha256": sha256(encoded).hexdigest(),
                "finite_retained": payload["replay"]["counts"]["finite_retained"],
                "one_dotted_retained": payload["replay"]["counts"]["dotted_retained"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
