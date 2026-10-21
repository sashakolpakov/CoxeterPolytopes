#!/usr/bin/env python3
"""Independent exact replay of the sparse 726-wall silver diagnostic.

This file deliberately does not import either silver Vinberg volume verifier
or either production 726-wall diagnose/certify script.  It rebuilds the Gram
matrix from the literal roots, uses the exact finite Coxeter classification,
and excludes affine boundary diagrams by determinant witnesses in two finite
quadratic fields with exact-inertia fallback.  It then compares the resulting
literal sets with the frozen production transcript.  Floating-point
eigenvalues are collected only as diagnostics; no decision depends on them.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
import math
from pathlib import Path
import sys
from typing import Mapping, Sequence

from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.exact.silver import (
    SILVER_ZERO,
    SQRT2,
    Silver,
    silver_inertia,
)
from coxeter_polytopes.search.silver_k3_vinberg import (
    silver_k3_alexeev_field_gram,
)


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/candidates/silver_k3_norm2_vinberg_extension_4_plus_3sqrt2.json"
DIAGNOSTIC = ROOT / "data/candidates/silver_k3_norm2_vinberg_726_global_diagnostic.json"
ARCHIVE = ROOT / "data/candidates/silver_k3_norm2_vinberg_726_independent_replay.json"
MAX_AFFINE_ORDER = 11


def _note(message: str) -> None:
    print(f"[independent silver-726] {message}", file=sys.stderr, flush=True)


def _canonical_digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _file_digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _decode_roots(document: Mapping[str, object]) -> tuple[tuple[Silver, ...], ...]:
    rows = document.get("ordered_prefix_roots")
    if not isinstance(rows, list):
        raise ValueError("input has no literal ordered_prefix_roots")
    roots: list[tuple[Silver, ...]] = []
    for index, row in enumerate(rows):
        if not isinstance(row, list) or len(row) != 10:
            raise ValueError(f"root {index} is not rank ten")
        if not all(isinstance(entry, Mapping) for entry in row):
            raise ValueError(f"root {index} has a malformed field entry")
        roots.append(tuple(Silver.from_json(entry) for entry in row))
    if len(roots) != 726 or document.get("wall_count") != 726:
        raise ValueError("the literal prefix is not the frozen 726-wall prefix")
    return tuple(roots)


def _build_gram(
    roots: Sequence[Sequence[Silver]],
) -> tuple[tuple[Silver, ...], ...]:
    lattice = silver_k3_alexeev_field_gram()
    images = tuple(
        tuple(
            sum(
                (lattice[row][column] * root[column] for column in range(10)),
                SILVER_ZERO,
            )
            for row in range(10)
        )
        for root in roots
    )
    work = [[SILVER_ZERO for _ in roots] for _ in roots]
    for left, image in enumerate(images):
        for right in range(left + 1):
            value = sum(
                (roots[right][index] * image[index] for index in range(10)),
                SILVER_ZERO,
            )
            work[left][right] = value
            work[right][left] = value
        if (left + 1) % 100 == 0 or left + 1 == len(roots):
            _note(f"Gram rows {left + 1}/{len(roots)}")
    return tuple(tuple(row) for row in work)


def _jsonable_gram(matrix: Sequence[Sequence[Silver]]) -> list[list[dict[str, str]]]:
    return [[entry.to_json() for entry in row] for row in matrix]


def _subgram(
    matrix: Sequence[Sequence[Silver]], vertices: Sequence[int]
) -> tuple[tuple[Silver, ...], ...]:
    return tuple(tuple(matrix[left][right] for right in vertices) for left in vertices)


def _least_bit(bits: int) -> tuple[int, int]:
    bit = bits & -bits
    return bit, bit.bit_length() - 1


def _finite_34(
    matrix: Sequence[Sequence[Silver]], vertices: tuple[int, ...]
) -> bool:
    """Exact finite-type classification for a diagram with labels 3 and 4."""

    neighbours = {vertex: set() for vertex in vertices}
    four_edges: set[frozenset[int]] = set()
    for offset, left in enumerate(vertices):
        for right in vertices[:offset]:
            value = matrix[left][right]
            if value == 0:
                continue
            if value not in (-1, -SQRT2):
                raise AssertionError("compatible diagram left the 3/4 alphabet")
            neighbours[left].add(right)
            neighbours[right].add(left)
            if value == -SQRT2:
                four_edges.add(frozenset((left, right)))

    unseen = set(vertices)
    while unseen:
        first = min(unseen)
        component = {first}
        frontier = [first]
        unseen.remove(first)
        while frontier:
            node = frontier.pop()
            added = neighbours[node] & unseen
            unseen.difference_update(added)
            component.update(added)
            frontier.extend(added)
        degrees = {node: len(neighbours[node] & component) for node in component}
        edge_count = sum(degrees.values()) // 2
        if edge_count != len(component) - 1:
            return False
        local_fours = [edge for edge in four_edges if edge.issubset(component)]
        if not local_fours:
            if max(degrees.values(), default=0) <= 2:  # A
                continue
            branch = [node for node, degree in degrees.items() if degree == 3]
            if len(branch) != 1 or any(degree > 3 for degree in degrees.values()):
                return False
            centre = branch[0]
            arms: list[int] = []
            for neighbour in neighbours[centre] & component:
                length = 1
                previous, current = centre, neighbour
                while degrees[current] == 2:
                    following = next(
                        node for node in neighbours[current] if node != previous
                    )
                    previous, current = current, following
                    length += 1
                arms.append(length)
            shape = tuple(sorted(arms))
            if shape[:2] == (1, 1) or shape in {
                (1, 2, 2),
                (1, 2, 3),
                (1, 2, 4),
            }:  # D or E
                continue
            return False
        if len(local_fours) != 1 or max(degrees.values(), default=0) > 2:
            return False
        four = local_fours[0]
        leaves = {node for node, degree in degrees.items() if degree <= 1}
        if four & leaves or len(component) == 4:  # B/C or F4
            continue
        return False
    return True


def _field_mul(
    left: tuple[int, int], right: tuple[int, int], prime: int
) -> tuple[int, int]:
    return (
        (left[0] * right[0] + 2 * left[1] * right[1]) % prime,
        (left[0] * right[1] + left[1] * right[0]) % prime,
    )


def _field_inverse(value: tuple[int, int], prime: int) -> tuple[int, int]:
    norm = (value[0] * value[0] - 2 * value[1] * value[1]) % prime
    inverse_norm = pow(norm, -1, prime)
    return value[0] * inverse_norm % prime, -value[1] * inverse_norm % prime


def _determinant_nonzero_mod_prime(
    matrix: Sequence[Sequence[Silver]], vertices: tuple[int, ...], prime: int
) -> bool:
    """Return a rigorous nonzero determinant witness in F_p[s]/(s^2-2)."""

    work: list[list[tuple[int, int]]] = []
    for left in vertices:
        row = []
        for right in vertices:
            value = matrix[left][right]
            if not value.is_integral:
                raise AssertionError("wall Gram unexpectedly left Z[sqrt(2)]")
            row.append((int(value.a) % prime, int(value.b) % prime))
        work.append(row)
    order = len(work)
    determinant = (1, 0)
    for column in range(order):
        pivot = next(
            (row for row in range(column, order) if work[row][column] != (0, 0)),
            None,
        )
        if pivot is None:
            return False
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            determinant = (-determinant[0] % prime, -determinant[1] % prime)
        pivot_value = work[column][column]
        determinant = _field_mul(determinant, pivot_value, prime)
        pivot_inverse = _field_inverse(pivot_value, prime)
        for row in range(column + 1, order):
            if work[row][column] == (0, 0):
                continue
            factor = _field_mul(work[row][column], pivot_inverse, prime)
            for index in range(column, order):
                product = _field_mul(factor, work[column][index], prime)
                work[row][index] = (
                    (work[row][index][0] - product[0]) % prime,
                    (work[row][index][1] - product[1]) % prime,
                )
    return determinant != (0, 0)


def _connected(vertices: tuple[int, ...], nonorthogonal: Sequence[int]) -> bool:
    vertices_bits = sum(1 << vertex for vertex in vertices)
    seen = 0
    frontier = 1 << vertices[0]
    while frontier:
        bit, vertex = _least_bit(frontier)
        frontier ^= bit
        seen |= bit
        frontier |= nonorthogonal[vertex] & vertices_bits & ~seen
    return seen == vertices_bits


def _numeric_matrix(matrix: Sequence[Sequence[Silver]]):
    try:
        import numpy as np
    except ImportError:
        return None, None
    result = np.array(
        [
            [float(entry.a) + float(entry.b) * math.sqrt(2.0) for entry in row]
            for row in matrix
        ],
        dtype=float,
    )
    return np, result


def _enumerate_exact(
    matrix: tuple[tuple[Silver, ...], ...]
) -> dict[str, object]:
    """Enumerate all elliptic subsets and the affine boundary exactly."""

    size = len(matrix)
    compatible = [0] * size
    nonorthogonal = [0] * size
    affine: set[tuple[int, ...]] = set()
    invalid_finite_pairs: list[tuple[int, int, Silver]] = []
    for left in range(size):
        if matrix[left][left] != 2:
            raise ValueError("independent replay requires norm-two walls")
        for right in range(left):
            pairing = matrix[left][right]
            if pairing.sign(0) > 0:
                raise ValueError("literal wall Gram is not acute")
            if pairing:
                nonorthogonal[left] |= 1 << right
                nonorthogonal[right] |= 1 << left
            minor = Silver(4) - pairing * pairing
            sign = minor.sign(0)
            if sign > 0:
                compatible[left] |= 1 << right
                compatible[right] |= 1 << left
                if pairing not in (0, -1, -SQRT2):
                    invalid_finite_pairs.append((right, left, pairing))
            elif sign == 0 and pairing:
                pair = (right, left)
                observed = silver_inertia(_subgram(matrix, pair), embedding=0)
                if observed != Inertia(1, 0, 1):
                    raise AssertionError("zero two-node minor has unexpected inertia")
                affine.add(pair)

    np, numerical = _numeric_matrix(matrix)
    spectral_min_accepted = math.inf
    spectral_max_rejected = -math.inf
    spectral_sign_mismatches = 0
    tested = Counter({1: size})
    positive = Counter({1: size})
    edges: list[tuple[int, ...]] = []
    vertices: list[tuple[int, ...]] = []
    exact_zero_extensions = len(affine)
    modular_nonzero_witnesses = 0
    exact_inertia_fallbacks = 0

    def visit(
        prefix: tuple[int, ...],
        candidates: int,
    ) -> None:
        nonlocal spectral_min_accepted
        nonlocal spectral_max_rejected
        nonlocal spectral_sign_mismatches
        nonlocal exact_zero_extensions
        nonlocal modular_nonzero_witnesses
        nonlocal exact_inertia_fallbacks
        if len(prefix) >= MAX_AFFINE_ORDER:
            return
        remaining = candidates
        while remaining:
            bit, vertex = _least_bit(remaining)
            remaining ^= bit
            extended = (*prefix, vertex)
            order = len(extended)
            tested[order] += 1
            is_positive = _finite_34(matrix, extended)
            is_affine = False
            if not is_positive and _connected(extended, nonorthogonal):
                if any(
                    _determinant_nonzero_mod_prime(matrix, extended, prime)
                    for prime in (5, 13)
                ):
                    modular_nonzero_witnesses += 1
                else:
                    exact_inertia_fallbacks += 1
                    observed = silver_inertia(_subgram(matrix, extended), embedding=0)
                    is_affine = observed == Inertia(order - 1, 0, 1) and all(
                        silver_inertia(
                            _subgram(
                                matrix,
                                extended[:omitted] + extended[omitted + 1 :],
                            ),
                            embedding=0,
                        )
                        == Inertia(order - 1, 0, 0)
                        for omitted in range(order)
                    )

            if numerical is not None and order >= 3:
                smallest = float(
                    np.linalg.eigvalsh(numerical[np.ix_(extended, extended)])[0]
                )
                if is_positive:
                    spectral_min_accepted = min(spectral_min_accepted, smallest)
                    if smallest <= 0:
                        spectral_sign_mismatches += 1
                elif not is_affine:
                    spectral_max_rejected = max(spectral_max_rejected, smallest)
                    if smallest >= 0:
                        spectral_sign_mismatches += 1
                elif abs(smallest) > 1e-8:
                    spectral_sign_mismatches += 1

            if is_positive:
                positive[order] += 1
                if order == 8:
                    edges.append(extended)
                elif order == 9:
                    vertices.append(extended)
                if order < MAX_AFFINE_ORDER:
                    visit(extended, remaining & compatible[vertex])
                continue

            if is_affine:
                exact_zero_extensions += 1
                affine.add(extended)

    all_bits = (1 << size) - 1
    for first in range(size):
        greater = all_bits & ~((1 << (first + 1)) - 1)
        visit((first,), greater & compatible[first])
        if (first + 1) % 100 == 0 or first + 1 == size:
            _note(
                f"exact diagram roots {first + 1}/{size}; "
                f"rank8={len(edges)}, rank9={len(vertices)}, affine={len(affine)}"
            )

    spectral = {
        "method": "numpy.linalg.eigvalsh on each candidate extension of order at least three",
        "proof_role": "none; finite type and affine determinant/inertia decisions are exact",
        "available": numerical is not None,
        "smallest_accepted_minimum_eigenvalue": (
            None if spectral_min_accepted == math.inf else round(spectral_min_accepted, 12)
        ),
        "closest_rejected_minimum_eigenvalue": (
            None if spectral_max_rejected == -math.inf else round(spectral_max_rejected, 12)
        ),
        "sign_mismatch_count": spectral_sign_mismatches,
    }
    return {
        "compatible_pair_count": sum(bits.bit_count() for bits in compatible) // 2,
        "invalid_finite_angle_pair_count": len(invalid_finite_pairs),
        "tested_subsets_by_order": {str(k): tested[k] for k in sorted(tested)},
        "elliptic_subsets_by_order": {str(k): positive[k] for k in sorted(positive)},
        "elliptic_edges": tuple(edges),
        "finite_vertices": tuple(vertices),
        "connected_affine_subdiagrams": tuple(sorted(affine)),
        "exact_affine_extension_count": exact_zero_extensions,
        "modular_nonzero_determinant_witness_count": modular_nonzero_witnesses,
        "exact_inertia_affine_fallback_count": exact_inertia_fallbacks,
        "spectral_diagnostics": spectral,
    }


def _literal_transcript(document: Mapping[str, object]) -> tuple[
    tuple[tuple[int, ...], ...],
    tuple[tuple[int, ...], ...],
    tuple[tuple[int, ...], ...],
]:
    raw = document.get("literal_transcript")
    if not isinstance(raw, Mapping):
        raise ValueError("production diagnostic lacks its literal transcript")
    values = []
    for key in ("elliptic_edges", "finite_vertices", "finite_endpoint_indices_by_edge"):
        rows = raw.get(key)
        if not isinstance(rows, list):
            raise ValueError(f"production transcript {key} is malformed")
        values.append(tuple(tuple(map(int, row)) for row in rows))
    return values[0], values[1], values[2]


def build_replay() -> dict[str, object]:
    input_document = json.loads(INPUT.read_text(encoding="utf-8"))
    diagnostic = json.loads(DIAGNOSTIC.read_text(encoding="utf-8"))
    roots = _decode_roots(input_document)
    matrix = _build_gram(roots)
    rebuilt_gram_digest = _canonical_digest(_jsonable_gram(matrix))
    if rebuilt_gram_digest != input_document.get("wall_gram_sha256"):
        raise AssertionError("independently rebuilt wall Gram digest differs")

    production_edges, production_vertices, production_rows = _literal_transcript(
        diagnostic
    )
    independent = _enumerate_exact(matrix)
    edges = independent.pop("elliptic_edges")
    vertices = independent.pop("finite_vertices")
    affine = independent.pop("connected_affine_subdiagrams")
    if edges != production_edges or vertices != production_vertices:
        raise AssertionError("independent literal elliptic sets differ from production")

    _note("direct exact inertia check of every archived edge and vertex")
    exact_edge_checks = 0
    exact_vertex_checks = 0
    for name, rows, order in (
        ("edge", edges, 8),
        ("vertex", vertices, 9),
    ):
        for index, row in enumerate(rows, 1):
            if silver_inertia(_subgram(matrix, row), embedding=0) != Inertia(order, 0, 0):
                raise AssertionError(f"archived {name} is not exactly elliptic")
            if name == "edge":
                exact_edge_checks += 1
            else:
                exact_vertex_checks += 1
            if index % 5_000 == 0:
                _note(f"exact {name} inertia {index}/{len(rows)}")

    edge_index = {edge: index for index, edge in enumerate(edges)}
    incidence = [set() for _ in edges]
    for vertex in vertices:
        for omitted, node in enumerate(vertex):
            edge = vertex[:omitted] + vertex[omitted + 1 :]
            found = edge_index.get(edge)
            if found is None:
                raise AssertionError("finite vertex deletion is absent from edge census")
            incidence[found].add(node)
    rebuilt_rows = tuple(tuple(sorted(row)) for row in incidence)
    if rebuilt_rows != production_rows:
        raise AssertionError("independent finite endpoint rows differ from production")
    histogram = dict(sorted(Counter(map(len, rebuilt_rows)).items()))
    incidence_sum = sum(map(len, rebuilt_rows))
    if incidence_sum != 9 * len(vertices):
        raise AssertionError("edge-vertex incidence double count failed")

    payload: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "silver_k3_726_independent_exact_replay",
        "input_manifest": str(INPUT.relative_to(ROOT)),
        "input_file_sha256": _file_digest(INPUT),
        "input_certificate_sha256": input_document.get("certificate_sha256"),
        "production_diagnostic": str(DIAGNOSTIC.relative_to(ROOT)),
        "production_diagnostic_file_sha256": _file_digest(DIAGNOSTIC),
        "production_diagnostic_sha256": diagnostic.get("diagnostic_sha256"),
        "facet_count": len(matrix),
        "dimension": 9,
        "rebuilt_wall_gram_sha256": rebuilt_gram_digest,
        "enumeration": independent,
        "literal_set_comparison": {
            "elliptic_rank_eight_count": len(edges),
            "elliptic_rank_eight_sha256": _canonical_digest(edges),
            "elliptic_rank_eight_exact_ordered_equality": True,
            "elliptic_rank_nine_count": len(vertices),
            "elliptic_rank_nine_sha256": _canonical_digest(vertices),
            "elliptic_rank_nine_exact_ordered_equality": True,
        },
        "direct_exact_inertia_checks": {
            "elliptic_rank_eight_checked": exact_edge_checks,
            "elliptic_rank_nine_checked": exact_vertex_checks,
            "all_positive_definite": True,
        },
        "finite_endpoint_incidence": {
            "row_count": len(rebuilt_rows),
            "literal_rows_sha256": _canonical_digest(rebuilt_rows),
            "exact_ordered_equality": True,
            "endpoint_histogram": {str(k): value for k, value in histogram.items()},
            "incidence_sum": incidence_sum,
            "nine_times_finite_vertex_count": 9 * len(vertices),
        },
        "affine_census": {
            "maximum_order_checked": MAX_AFFINE_ORDER,
            "connected_affine_subdiagram_count": len(affine),
            "connected_affine_subdiagrams": [list(row) for row in affine],
            "ideal_endpoint_count": 0 if not affine else None,
            "full_zero_neighborhood_test": (
                "vacuous because the independently replayed connected-affine census is empty"
                if not affine
                else "not evaluated by this independent replay"
            ),
        },
        "finite_volume_status": {
            "finite_endpoint_histogram": {str(k): value for k, value in histogram.items()},
            "finite_volume_open_edge_count": sum(
                count for endpoints, count in histogram.items() if endpoints != 2
            ),
            "finite_volume_stopping_available": False,
            "new_finite_volume_hyperbolic_example": False,
        },
        "scope": {
            "simple_root_prefix_complete": False,
            "failure_proves_eventual_infinite_volume": False,
            "reason": (
                "Later simple walls may close a current open edge or remove it from the final chamber."
            ),
            "other_reflective_norms_exhausted": False,
            "spectral_values_are_proof_inputs": False,
        },
    }
    payload["certificate_sha256"] = _canonical_digest(payload)
    return payload


def verify_archive(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "silver_k3_726_independent_exact_replay":
        raise ValueError("wrong independent replay kind")
    payload = dict(document)
    expected = payload.pop("certificate_sha256", None)
    if not isinstance(expected, str) or _canonical_digest(payload) != expected:
        raise ValueError("independent replay certificate digest changed")
    if document.get("input_file_sha256") != _file_digest(INPUT):
        raise ValueError("independent replay input file changed")
    if document.get("production_diagnostic_file_sha256") != _file_digest(DIAGNOSTIC):
        raise ValueError("production diagnostic linked by the replay changed")

    diagnostic = json.loads(DIAGNOSTIC.read_text(encoding="utf-8"))
    edges, vertices, endpoint_rows = _literal_transcript(diagnostic)
    literal = document.get("literal_set_comparison")
    inertia_checks = document.get("direct_exact_inertia_checks")
    incidence = document.get("finite_endpoint_incidence")
    affine = document.get("affine_census")
    status = document.get("finite_volume_status")
    scope = document.get("scope")
    enumeration = document.get("enumeration")
    if not all(
        isinstance(value, Mapping)
        for value in (literal, inertia_checks, incidence, affine, status, scope, enumeration)
    ):
        raise ValueError("independent replay sections are malformed")
    if (
        literal.get("elliptic_rank_eight_count") != 25_764
        or literal.get("elliptic_rank_nine_count") != 3_604
        or literal.get("elliptic_rank_eight_sha256") != _canonical_digest(edges)
        or literal.get("elliptic_rank_nine_sha256") != _canonical_digest(vertices)
        or literal.get("elliptic_rank_eight_exact_ordered_equality") is not True
        or literal.get("elliptic_rank_nine_exact_ordered_equality") is not True
    ):
        raise ValueError("independent literal set comparison changed")
    if (
        inertia_checks.get("elliptic_rank_eight_checked") != 25_764
        or inertia_checks.get("elliptic_rank_nine_checked") != 3_604
        or inertia_checks.get("all_positive_definite") is not True
    ):
        raise ValueError("direct exact inertia transcript changed")
    histogram = dict(sorted(Counter(map(len, endpoint_rows)).items()))
    if (
        incidence.get("literal_rows_sha256") != _canonical_digest(endpoint_rows)
        or incidence.get("exact_ordered_equality") is not True
        or incidence.get("incidence_sum") != 32_436
        or incidence.get("nine_times_finite_vertex_count") != 32_436
        or {int(k): int(v) for k, v in incidence.get("endpoint_histogram", {}).items()}
        != histogram
        or histogram != {0: 2_332, 1: 14_428, 2: 9_004}
    ):
        raise ValueError("independent incidence transcript changed")
    spectral = enumeration.get("spectral_diagnostics")
    if not isinstance(spectral, Mapping):
        raise ValueError("spectral diagnostic is malformed")
    if (
        enumeration.get("elliptic_subsets_by_order", {}).get("8") != 25_764
        or enumeration.get("elliptic_subsets_by_order", {}).get("9") != 3_604
        or affine.get("maximum_order_checked") != 11
        or affine.get("connected_affine_subdiagram_count") != 0
        or affine.get("connected_affine_subdiagrams") != []
        or affine.get("ideal_endpoint_count") != 0
        or spectral.get("proof_role")
        != "none; finite type and affine determinant/inertia decisions are exact"
        or scope.get("spectral_values_are_proof_inputs") is not False
    ):
        raise ValueError("independent affine/exactness transcript changed")
    if (
        status.get("finite_volume_open_edge_count") != 16_760
        or status.get("finite_volume_stopping_available") is not False
        or status.get("new_finite_volume_hyperbolic_example") is not False
        or scope.get("simple_root_prefix_complete") is not False
        or scope.get("failure_proves_eventual_infinite_volume") is not False
    ):
        raise ValueError("independent replay scope changed")
    return {
        "elliptic_rank_eight_count": 25_764,
        "elliptic_rank_nine_count": 3_604,
        "exact_inertia_check_count": 29_368,
        "finite_incidence_sum": 32_436,
        "connected_affine_count": 0,
        "maximum_affine_order_checked": 11,
        "finite_volume_open_edge_count": 16_760,
        "new_hyperbolic_example": False,
        "eventual_infinite_volume_proved": False,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="replace the compact replay archive")
    parser.add_argument("--verify", action="store_true", help="verify the frozen compact archive")
    args = parser.parse_args(argv)
    if args.verify:
        result = verify_archive(json.loads(ARCHIVE.read_text(encoding="utf-8")))
    else:
        result = build_replay()
        if args.write:
            ARCHIVE.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
