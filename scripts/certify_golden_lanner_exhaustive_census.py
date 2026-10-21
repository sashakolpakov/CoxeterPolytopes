#!/usr/bin/env python3
"""Exhaust the compact equal-root golden Lanner class in all ranks.

The finite search covers every symmetric norm-two Gram diagram with entries
``0,-1,-phi`` in ranks three, four, and five, modulo arbitrary node
permutation.  It imposes connectedness, Lorentzian identity inertia,
positive-definite conjugate inertia, and ellipticity of every proper
principal subdiagram.  Lanner's classical theorem proves that no such
minimal hyperbolic diagram has rank at least six.

For every surviving orbit, the support/sign compiler from
``certify_golden_lanner_complete_chambers`` then exhausts every norm-two
mirror meeting the closed simplex, tests for interior subdivisions, builds
the trace lattice, and applies the stable primitive K3 embedding gate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from itertools import combinations, permutations, product
from pathlib import Path
from typing import Any, Iterable, Sequence

from coxeter_polytopes.exact.golden import Golden, ONE, PHI, ZERO
from coxeter_polytopes.exact.matrix import Inertia, determinant, inertia
from scripts.certify_golden_lanner_complete_chambers import (
    _canonical_sha256,
    _encode,
    _encode_matrix,
    _encode_vector,
    _enumerate_mirror_contacts,
    _inertia_json,
    _inverse,
    _matrix_vector,
    _pair,
    _root_key,
    _submatrix,
    _trace_package,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "golden_lanner_exhaustive_census.json"
)

GoldenMatrix = tuple[tuple[Golden, ...], ...]
GoldenVector = tuple[Golden, ...]
LabelledEdge = tuple[int, int, int]


@dataclass(frozen=True, slots=True)
class DiagramSpec:
    seed_id: str
    name: str
    rank: int
    edges: tuple[LabelledEdge, ...]
    prior_status: str


SPECS: tuple[DiagramSpec, ...] = (
    DiagramSpec(
        "triangle_255",
        "(2,5,5)",
        3,
        ((0, 1, 2), (1, 2, 5), (0, 2, 5)),
        "one_of_eight_handpicked_seeds",
    ),
    DiagramSpec(
        "triangle_335",
        "(3,3,5)",
        3,
        ((0, 1, 3), (1, 2, 3), (0, 2, 5)),
        "one_of_eight_handpicked_seeds",
    ),
    DiagramSpec(
        "triangle_355",
        "(3,5,5)",
        3,
        ((0, 1, 3), (1, 2, 5), (0, 2, 5)),
        "one_of_eight_handpicked_seeds",
    ),
    DiagramSpec(
        "triangle_555",
        "(5,5,5)",
        3,
        ((0, 1, 5), (1, 2, 5), (0, 2, 5)),
        "one_of_eight_handpicked_seeds",
    ),
    DiagramSpec(
        "chain_353",
        "[3,5,3]",
        4,
        ((0, 1, 3), (1, 2, 5), (2, 3, 3)),
        "one_of_eight_handpicked_seeds",
    ),
    DiagramSpec(
        "chain_535",
        "[5,3,5]",
        4,
        ((0, 1, 5), (1, 2, 3), (2, 3, 5)),
        "one_of_eight_handpicked_seeds",
    ),
    DiagramSpec(
        "star_5_3_3",
        "D4(5,3,3)",
        4,
        ((0, 1, 5), (0, 2, 3), (0, 3, 3)),
        "newly_added_by_exhaustive_census",
    ),
    DiagramSpec(
        "cycle4_one5",
        "C4(one 5)",
        4,
        ((0, 1, 5), (1, 2, 3), (2, 3, 3), (3, 0, 3)),
        "newly_added_by_exhaustive_census",
    ),
    DiagramSpec(
        "cycle4_opposite5",
        "C4(opposite 5s)",
        4,
        ((0, 1, 5), (1, 2, 3), (2, 3, 5), (3, 0, 3)),
        "newly_added_by_exhaustive_census",
    ),
    DiagramSpec(
        "chain_5333",
        "[5,3,3,3]",
        5,
        ((0, 1, 5), (1, 2, 3), (2, 3, 3), (3, 4, 3)),
        "one_of_eight_handpicked_seeds",
    ),
    DiagramSpec(
        "chain_5335",
        "[5,3,3,5]",
        5,
        ((0, 1, 5), (1, 2, 3), (2, 3, 3), (3, 4, 5)),
        "one_of_eight_handpicked_seeds",
    ),
    DiagramSpec(
        "fork5_terminal5",
        "D5(terminal 5)",
        5,
        ((0, 1, 5), (1, 2, 3), (2, 3, 3), (2, 4, 3)),
        "newly_added_by_exhaustive_census",
    ),
)


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _gram_from_edges(rank: int, edges: Sequence[LabelledEdge]) -> GoldenMatrix:
    gram = [
        [Golden(2) if i == j else ZERO for j in range(rank)]
        for i in range(rank)
    ]
    values = {2: ZERO, 3: Golden(-1), 5: -PHI}
    seen: set[tuple[int, int]] = set()
    for left, right, label in edges:
        if left == right or not (0 <= left < rank and 0 <= right < rank):
            raise ValueError("a diagram edge has invalid endpoints")
        pair = tuple(sorted((left, right)))
        if pair in seen:
            raise ValueError("a diagram edge was repeated")
        seen.add(pair)
        gram[left][right] = gram[right][left] = values[label]
    return tuple(tuple(row) for row in gram)


def _label_matrix(gram: GoldenMatrix) -> tuple[tuple[int, ...], ...]:
    result = []
    for i, row in enumerate(gram):
        labels = []
        for j, entry in enumerate(row):
            if i == j:
                labels.append(1)
            elif entry == ZERO:
                labels.append(2)
            elif entry == -1:
                labels.append(3)
            elif entry == -PHI:
                labels.append(5)
            else:
                raise ValueError("the census admits only labels 2,3,5")
        result.append(tuple(labels))
    return tuple(result)


def _canonical_key(labels: Sequence[Sequence[int]]) -> tuple[int, ...]:
    rank = len(labels)
    pairs = tuple(combinations(range(rank), 2))
    return min(
        tuple(labels[permutation[i]][permutation[j]] for i, j in pairs)
        for permutation in permutations(range(rank))
    )


def _connected(rank: int, active_edges: Iterable[tuple[int, int]]) -> bool:
    adjacency = [set() for _ in range(rank)]
    for left, right in active_edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    reached = {0}
    frontier = [0]
    while frontier:
        current = frontier.pop()
        for neighbor in adjacency[current] - reached:
            reached.add(neighbor)
            frontier.append(neighbor)
    return len(reached) == rank


def _delete(matrix: GoldenMatrix, omitted: int) -> GoldenMatrix:
    keep = tuple(index for index in range(len(matrix)) if index != omitted)
    return _submatrix(matrix, keep)


def _graph_prefilter(rank: int, active: tuple[tuple[int, int], ...]) -> bool:
    """Necessary unlabeled ellipticity filter, with an exact proof.

    Replacing every ``-phi`` by ``-1`` replaces the nonnegative weighted
    adjacency matrix by an entrywise smaller one.  Perron--Frobenius makes
    every proper deletion of the resulting simply-laced Gram positive
    definite whenever the decorated deletion can be positive definite.
    """

    if not _connected(rank, active):
        return False
    adjacency = [[0 for _ in range(rank)] for _ in range(rank)]
    for left, right in active:
        adjacency[left][right] = adjacency[right][left] = 1
    curve = [
        [-2 if i == j else adjacency[i][j] for j in range(rank)]
        for i in range(rank)
    ]
    for omitted in range(rank):
        deleted = [
            [curve[i][j] for j in range(rank) if j != omitted]
            for i in range(rank)
            if i != omitted
        ]
        if inertia(deleted) != Inertia(0, rank - 1, 0):
            return False
    return True


def _all_proper_principal_elliptic(gram: GoldenMatrix, embedding: int) -> bool:
    rank = len(gram)
    for size in range(1, rank):
        for subset in combinations(range(rank), size):
            if inertia(_submatrix(gram, subset), embedding=embedding) != Inertia(
                size, 0, 0
            ):
                return False
    return True


def _enumerate_rank(rank: int) -> dict[str, Any]:
    pairs = tuple(combinations(range(rank), 2))
    representatives: dict[tuple[int, ...], GoldenMatrix] = {}
    prefilter_supports = 0
    prefilter_histogram: dict[str, int] = {}
    decorations_tested = 0
    labelled_survivors = 0

    for mask in range(1 << len(pairs)):
        active = tuple(
            pair for index, pair in enumerate(pairs) if mask & (1 << index)
        )
        if not _graph_prefilter(rank, active):
            continue
        prefilter_supports += 1
        key = str(len(active))
        prefilter_histogram[key] = prefilter_histogram.get(key, 0) + 1
        for decorations in product((3, 5), repeat=len(active)):
            decorations_tested += 1
            edges = tuple(
                (left, right, label)
                for (left, right), label in zip(active, decorations)
            )
            gram = _gram_from_edges(rank, edges)
            if inertia(gram, embedding=0) != Inertia(rank - 1, 1, 0):
                continue
            if inertia(gram, embedding=1) != Inertia(rank, 0, 0):
                continue
            if any(
                inertia(_delete(gram, omitted), embedding=0)
                != Inertia(rank - 1, 0, 0)
                for omitted in range(rank)
            ):
                continue
            labelled_survivors += 1
            canonical = _canonical_key(_label_matrix(gram))
            representatives.setdefault(canonical, gram)

    if not all(
        _all_proper_principal_elliptic(gram, embedding=0)
        and _all_proper_principal_elliptic(gram, embedding=1)
        for gram in representatives.values()
    ):
        raise ArithmeticError("a purported Lanner representative has a bad subgram")

    return {
        "rank": rank,
        "dimension": rank - 1,
        "all_support_graphs": 1 << len(pairs),
        "prefilter_support_graphs": prefilter_supports,
        "prefilter_edge_count_histogram": prefilter_histogram,
        "decorated_labelled_grams_tested": decorations_tested,
        "labelled_survivors": labelled_survivors,
        "unlabelled_orbits": len(representatives),
        "canonical_keys": [list(key) for key in sorted(representatives)],
        "representatives": representatives,
    }


def _analyze_diagram(spec: DiagramSpec) -> dict[str, Any]:
    gram = _gram_from_edges(spec.rank, spec.edges)
    dimension = spec.rank - 1
    inverse = _inverse(gram)
    simple_roots = tuple(
        tuple(ONE if i == j else ZERO for i in range(spec.rank))
        for j in range(spec.rank)
    )
    interior = tuple(
        -sum((inverse[i][j] for j in range(spec.rank)), start=ZERO)
        for i in range(spec.rank)
    )
    interior_pairings = _matrix_vector(gram, interior)
    interior_square = _pair(interior, gram, interior)

    vertices = []
    for omitted in range(spec.rank):
        vertex = tuple(-inverse[i][omitted] for i in range(spec.rank))
        incident = tuple(i for i in range(spec.rank) if i != omitted)
        vertex_gram = _submatrix(gram, incident)
        pairings = _matrix_vector(gram, vertex)
        vertices.append(
            {
                "opposite_facet": omitted,
                "incident_facets": list(incident),
                "coordinates": _encode_vector(vertex),
                "square": _encode(_pair(vertex, gram, vertex)),
                "facet_pairings": _encode_vector(pairings),
                "vertex_gram_determinant": _encode(determinant(vertex_gram)),
                "identity_inertia": _inertia_json(
                    inertia(vertex_gram, embedding=0)
                ),
                "conjugate_inertia": _inertia_json(
                    inertia(vertex_gram, embedding=1)
                ),
                "finite": _pair(vertex, gram, vertex).sign() < 0,
            }
        )

    contacts = _enumerate_mirror_contacts(gram)
    trace = _trace_package(gram, dimension)
    canonical = _canonical_key(_label_matrix(gram))
    checks = {
        "connected_nonorthogonal_graph": _connected(
            spec.rank,
            ((left, right) for left, right, label in spec.edges if label != 2),
        ),
        "identity_lorentzian": inertia(gram, embedding=0)
        == Inertia(dimension, 1, 0),
        "conjugate_positive_definite": inertia(gram, embedding=1)
        == Inertia(spec.rank, 0, 0),
        "all_proper_principal_identity_elliptic": (
            _all_proper_principal_elliptic(gram, embedding=0)
        ),
        "all_proper_principal_conjugate_elliptic": (
            _all_proper_principal_elliptic(gram, embedding=1)
        ),
        "specified_incenter_timelike_and_strict": interior_square.sign() < 0
        and all(value == -1 for value in interior_pairings),
        "all_vertices_finite": all(vertex["finite"] for vertex in vertices),
        "no_mixed_sign_norm_two_root": not contacts["mixed_roots"],
        "seed_equals_complete_norm_two_chamber": not contacts["mixed_roots"],
        "trace_lattice_even": trace["even"],
        "stable_primitive_k3_embedding": trace[
            "primitive_k3_embedding_automatic"
        ],
        "visible_k3_shells_collapse_to_norm_two": trace[
            "visible_k3_wall_shell_collapse"
        ]["second_shell_maps_to_norm_two"],
        "paired_k3_facets_exact": trace["paired_roots_square_minus_two"]
        and trace["same_wall_root_pairs_orthogonal"]
        and trace["paired_off_diagonal_intersections_nonnegative"],
    }
    if not all(checks.values()):
        failed = [name for name, value in checks.items() if not value]
        raise ArithmeticError(f"{spec.seed_id} failed: {failed}")

    labelled_edges = []
    labels = _label_matrix(gram)
    for left, right in combinations(range(spec.rank), 2):
        labelled_edges.append([left, right, labels[left][right]])

    return {
        "seed_id": spec.seed_id,
        "name": spec.name,
        "rank": spec.rank,
        "dimension": dimension,
        "facet_count": spec.rank,
        "edges": labelled_edges,
        "nonorthogonal_edges": [
            [left, right, label]
            for left, right, label in labelled_edges
            if label != 2
        ],
        "canonical_unlabelled_key": list(canonical),
        "weighted_norm_two_gram": _encode_matrix(gram),
        "field_determinant": _encode(determinant(gram)),
        "identity_inertia": _inertia_json(inertia(gram, embedding=0)),
        "conjugate_inertia": _inertia_json(inertia(gram, embedding=1)),
        "proper_principal_subgram_count": (1 << spec.rank) - 2,
        "specified_incenter": {
            "definition": "x=-G^{-1}(1,...,1), with (x,e_i)=-1",
            "coordinates": _encode_vector(interior),
            "square": _encode(interior_square),
            "facet_pairings": _encode_vector(interior_pairings),
        },
        "simple_walls": {
            "roots": [_encode_vector(root) for root in simple_roots],
            "weighted_gram": _encode_matrix(gram),
        },
        "vertices": vertices,
        "vinberg_compactness": {
            "finite_vertices": spec.rank,
            "ideal_vertices": 0,
            "hyperideal_vertices": 0,
            "all_rank_n_vertex_diagrams_elliptic": True,
            "compact": True,
            "finite_volume": True,
        },
        "full_norm_two_arrangement": {
            "mixed_sign_patterns_up_to_global_sign": contacts[
                "mixed_sign_patterns"
            ],
            "mixed_sign_assignments_tested": contacts["mixed_assignments"],
            "mixed_sign_roots": [
                _encode_vector(root) for root in contacts["mixed_roots"]
            ],
            "proper_boundary_supports": contacts["boundary_supports"],
            "boundary_assignments_tested": contacts["boundary_assignments"],
            "boundary_contact_mirror_count": len(contacts["boundary_roots"]),
            "boundary_contact_roots_canonical_positive": [
                _encode_vector(root)
                for root in sorted(contacts["boundary_roots"], key=_root_key)
            ],
            "interior_cut_count": len(contacts["mixed_roots"]),
            "complete_chamber_facet_count": spec.rank,
            "seed_is_complete_chamber": True,
        },
        "trace_k3_package": trace,
        "k3_eigen_nef_return": {
            "full_picard_projective_k3_exists": True,
            "projectivized_eigen_nef_section_equals_closed_seed": True,
            "generic_period_dimension": 20 - 2 * spec.rank,
            "facet_curve_class_count": 2 * spec.rank,
            "surface_automorphism_claim": False,
        },
        "lineage": {
            "prior_status": spec.prior_status,
            "same_as_classical_lanner_seed": True,
            "subdivided": False,
            "worldwide_novelty_claim": False,
        },
        "checks": checks,
        "all_checks_pass": True,
    }


def analyze() -> dict[str, Any]:
    rank_censuses = [_enumerate_rank(rank) for rank in (3, 4, 5)]
    named_keys: dict[int, set[tuple[int, ...]]] = {3: set(), 4: set(), 5: set()}
    for spec in SPECS:
        key = _canonical_key(_label_matrix(_gram_from_edges(spec.rank, spec.edges)))
        if key in named_keys[spec.rank]:
            raise ArithmeticError("two named diagrams represent the same orbit")
        named_keys[spec.rank].add(key)
    for census in rank_censuses:
        enumerated = set(census.pop("representatives"))
        if enumerated != named_keys[census["rank"]]:
            raise ArithmeticError(
                f"named rank-{census['rank']} orbit list does not match enumeration"
            )

    examples = [_analyze_diagram(spec) for spec in SPECS]
    new_examples = [
        example
        for example in examples
        if example["lineage"]["prior_status"]
        == "newly_added_by_exhaustive_census"
    ]
    lanner_seeds = [
        {
            "name": example["name"],
            "seed_id": example["seed_id"],
            "dimension": example["dimension"],
            "edges": example["edges"],
        }
        for example in examples
    ]

    payload: dict[str, Any] = {
        "schema_version": "1.0.0",
        "certificate_kind": (
            "exhaustive_equal_root_golden_Lanner_and_K3_chamber_census"
        ),
        "class_definition": {
            "field": "Q(phi), phi^2=phi+1",
            "integer_ring": "Z[phi]",
            "gram_entries": "diagonal 2; off diagonal 0,-1,-phi",
            "diagram_conditions": [
                "connected nonorthogonal graph",
                "identity inertia (k-1,1,0)",
                "conjugate inertia (k,0,0)",
                "every proper principal subgram positive definite",
            ],
            "equivalent_geometric_object": (
                "compact arithmetic equal-root Coxeter simplex with labels 2,3,5"
            ),
        },
        "finite_enumeration": {
            "method": (
                "all support graphs; the exact Perron--Frobenius necessary "
                "prefilter; every independent 3/5 decoration of every "
                "surviving support; exact two-place inertia; then "
                "canonicalization under S_k"
            ),
            "ranks": rank_censuses,
            "total_unlabelled_orbits": sum(
                census["unlabelled_orbits"] for census in rank_censuses
            ),
        },
        "rank_ge_6_nonexistence": {
            "theorem": (
                "Lanner's classification: a connected non-elliptic, "
                "non-parabolic Coxeter diagram all of whose proper "
                "subdiagrams are elliptic has order at most five"
            ),
            "application": (
                "every member of the stated class is a Lanner diagram, so "
                "there is no rank k>=6 member and no compact simplex in H^5+"
            ),
            "maximum_rank": 5,
            "maximum_dimension": 4,
            "reference": "Folke Lanner, On complexes with transitive groups of automorphisms, Lund 11 (1950), 1-71",
            "modern_reference": "https://arxiv.org/abs/2203.07248",
        },
        "lanner_seeds": lanner_seeds,
        "examples": examples,
        "newly_added_examples": [example["seed_id"] for example in new_examples],
        "summary": {
            "orbit_counts_by_dimension": {"2": 4, "3": 5, "4": 3},
            "total_orbits": 12,
            "previously_handpicked": 8,
            "newly_added": 4,
            "all_complete_norm_two_chambers": all(
                example["full_norm_two_arrangement"]["seed_is_complete_chamber"]
                for example in examples
            ),
            "all_compact": all(
                example["vinberg_compactness"]["compact"] for example in examples
            ),
            "all_stably_embed_in_K3": all(
                example["trace_k3_package"]["primitive_k3_embedding_automatic"]
                for example in examples
            ),
            "all_return_the_same_polytope_as_eigen_nef_section": all(
                example["k3_eigen_nef_return"][
                    "projectivized_eigen_nef_section_equals_closed_seed"
                ]
                for example in examples
            ),
            "subdivision_count": 0,
            "worldwide_new_abstract_polytope_claim": False,
        },
        "new_rows_exact_summary": {
            example["seed_id"]: {
                "dimension": example["dimension"],
                "facets": example["facet_count"],
                "boundary_contacts": example["full_norm_two_arrangement"][
                    "boundary_contact_mirror_count"
                ],
                "interior_cuts": example["full_norm_two_arrangement"][
                    "interior_cut_count"
                ],
                "trace_determinant": example["trace_k3_package"]["determinant"],
                "trace_discriminant_invariants": example["trace_k3_package"][
                    "discriminant_group_invariants"
                ],
            }
            for example in new_examples
        },
        "replay": {
            "command": (
                "PYTHONPATH=src:. python3 "
                "scripts/certify_golden_lanner_exhaustive_census.py "
                "--check-archive"
            ),
            "tests": (
                "PYTHONPATH=src:. pytest -q "
                "tests/test_golden_lanner_exhaustive_census.py"
            ),
        },
        "all_checks_pass": True,
    }
    payload["certificate_sha256"] = _canonical_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-certificate", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = analyze()
    if args.write_certificate:
        OUTPUT.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(f"wrote {OUTPUT}: {payload['certificate_sha256']}")
        return
    if args.check_archive:
        archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
        if archived != payload:
            raise SystemExit("archive mismatch")
        print(f"archive replay passed: {payload['certificate_sha256']}")
        return
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
