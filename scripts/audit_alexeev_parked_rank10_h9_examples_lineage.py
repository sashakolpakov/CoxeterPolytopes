#!/usr/bin/env python3
"""Hostile replay and bounded lineage audit for five parked Alexeev H9 sections.

The script deliberately does not import the batch certifier.  It replays the
archived H-cones and rays, reruns GMP double description, checks Vinberg's
literal edge/cusp criterion, exhausts reflection-stable integral overlattices,
and compares the five labelled diagrams with pinned project and CoxIter stocks.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter, defaultdict
from fractions import Fraction
from functools import reduce
import gzip
from hashlib import sha256
import itertools
import json
import math
from pathlib import Path
import subprocess
import tarfile
from typing import Any, Iterable, Sequence
import zipfile

import networkx as nx
import sympy as sp

from scripts.certify_alexeev_remaining_stable_overlattices import (
    _overlattice_basis,
    _q_histogram,
    _smith,
    _stable_subgroups,
)
from scripts.certify_k3_s3_h10_redoubling_family import _parse_coxiter_text
from scripts.certify_shimada_disc16_d5_h10_30_vinberg_faces import DiagramEnumerator
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
BATCH = ROOT / "data" / "enumerations" / "alexeev_parked_rank10_h9_examples.json"
FACTORY = ROOT / "data" / "enumerations" / "alexeev_k3_h9_14_factory.json"
CUBES = ROOT / "data" / "enumerations" / "alexeev_h9_doubling_cubes.json"
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_parked_rank10_h9_examples_hostile_lineage.json"
)
COXITER_ZIP = Path("/private/tmp/coxiter-current.zip")
DOUBLING_TAR = Path(
    "/private/tmp/coxiter-current/CoxIter-master/graphs/testing_index2.tar.gz"
)
SK_ROOT = Path("/private/tmp/Scharlau-Kirschmer-Lorentz")
SK_LIST9 = SK_ROOT / "List9"
SK_GRAM9 = SK_ROOT / "Gram9"

SCHEMA = "coxeter.alexeev_parked_rank10_h9.hostile_lineage.v1"
BATCH_SHA256 = "e181b8c021fe0f37d23b53b33e4e3900559502d0f8c982a4a9e9cbcf0866e430"
FACTORY_SHA256 = "b1ee6fd769527ae6dab01a2bb58099e7c7d1776d0417b6a8d710316697e8a775"
CUBES_SHA256 = "4919728a5c5ff8adebd80b250cfedfb562a040c8bd1e72e551b398b1f13e7e4c"
COXITER_ZIP_SHA256 = "ae04648b291c1e5f3d0e47d666c846dc805b244469a20486bca59674b4d171c3"
DOUBLING_TAR_SHA256 = "3330ca41fb3e047123ace718f08542bc1dd0019b35922d4fa836a217647e7036"
SK_LIST9_SHA256 = "b7f6982a5c03b1ff07a5a9f383517438a926a499d27fb33ea89a8a69d7fb66c3"
SK_GRAM9_SHA256 = "772187d28a067c2d10f2fe00d7cee07a063a6c3d79be5cc8400c2c0ca8c96c36"
IMHOF_PATH = "CoxIter-master/graphs/doubly truncated orthoscheme/9-imhof85_03.coxiter"
IMHOF_SHA256 = "8f61b20c7e8ccc658073a14d020cec94f67fcb6a8b62ffa7a7dc836cbf5b7171"
IMHOF_PERMUTATION_ONE_BASED = [6, 11, 10, 9, 8, 1, 2, 3, 4, 5, 7, 12]

EXPECTED: dict[str, dict[str, Any]] = {
    "A14_8_0_C01": {
        "facets": 12,
        "finite": 44,
        "ideal": 5,
        "edges": 260,
        "affine": 10,
        "stable_orders": [1, 4],
        "stable_even_orders": [1, 4],
        "maximal_index": 4,
        "two_elementary": [10, 6, 0],
        "sk_operation": "scaled_dual_M_star_2",
        "sk_entry": 23,
    },
    "A14_8_1_C02": {
        "facets": 16,
        "finite": 165,
        "ideal": 23,
        "edges": 1112,
        "affine": 55,
        "stable_orders": [1, 2, 4, 8],
        "stable_even_orders": [1, 2, 4, 8],
        "maximal_index": 8,
        "two_elementary": [10, 4, 0],
        "sk_operation": "sign_reversal",
        "sk_entry": 23,
    },
    "A14_10_1_C01": {
        "facets": 13,
        "finite": 70,
        "ideal": 7,
        "edges": 393,
        "affine": 14,
        "stable_orders": [1, 2, 4],
        "stable_even_orders": [1, 2, 4],
        "maximal_index": 4,
        "two_elementary": [10, 8, 0],
        "sk_operation": "scaled_dual_M_star_2",
        "sk_entry": 7,
    },
    "A14_10_1_C02": {
        "facets": 15,
        "finite": 123,
        "ideal": 15,
        "edges": 778,
        "affine": 35,
        "stable_orders": [1, 2, 4, 8],
        "stable_even_orders": [1, 2, 4, 8],
        "maximal_index": 8,
        "two_elementary": [10, 6, 0],
        "sk_operation": "scaled_dual_M_star_2",
        "sk_entry": 23,
    },
    "A16_14_1_C01": {
        "facets": 15,
        "finite": 100,
        "ideal": 9,
        "edges": 588,
        "affine": 23,
        "stable_orders": [1, 16, 32],
        "stable_even_orders": [1, 16, 32],
        "maximal_index": 32,
        "two_elementary": [10, 4, 0],
        "sk_operation": "sign_reversal",
        "sk_entry": 23,
    },
}


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _matrix(payload: Sequence[Sequence[str | int]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(value) for value in row] for row in payload])


def _pair(gram: sp.Matrix, left: sp.Matrix, right: sp.Matrix) -> sp.Rational:
    return sp.Rational((left.T * gram * right)[0])


def _primitive(values: Iterable[sp.Rational | int]) -> tuple[int, ...]:
    entries = [sp.Rational(value) for value in values]
    denominator = int(sp.ilcm(*(entry.q for entry in entries)))
    integers = [int(entry * denominator) for entry in entries]
    divisor = reduce(math.gcd, (abs(value) for value in integers if value))
    return tuple(value // divisor for value in integers)


def _reflection(gram: sp.Matrix, root: sp.Matrix) -> sp.Matrix:
    return sp.simplify(
        sp.eye(gram.rows)
        - sp.Rational(2, _pair(gram, root, root)) * root * (root.T * gram)
    )


def _token(gram: sp.Matrix, left: int, right: int) -> str:
    if left == right:
        return "1"
    ratio = sp.factor(
        gram[left, right] ** 2 / (gram[left, left] * gram[right, right])
    )
    known = {
        sp.Rational(0): "2",
        sp.Rational(1, 4): "3",
        sp.Rational(1, 2): "4",
        sp.Rational(3, 4): "6",
        sp.Rational(1): "parallel",
    }
    if ratio in known:
        return known[ratio]
    if ratio > 1:
        return f"dotted:{ratio}"
    raise ArithmeticError(f"unexpected normalized product {ratio}")


def _relations(gram: sp.Matrix) -> tuple[tuple[str, ...], ...]:
    return tuple(
        tuple(_token(gram, left, right) for right in range(gram.rows))
        for left in range(gram.rows)
    )


def _histogram(relations: Sequence[Sequence[str]]) -> dict[str, int]:
    return dict(
        sorted(
            Counter(
                relations[left][right]
                for left in range(len(relations))
                for right in range(left)
            ).items()
        )
    )


def _graph(relations: Sequence[Sequence[str]]) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(range(len(relations)))
    for left in range(len(relations)):
        for right in range(left):
            if relations[left][right] != "2":
                graph.add_edge(left, right, relation=relations[left][right])
    return graph


EDGE_MATCH = nx.algorithms.isomorphism.categorical_edge_match("relation", None)


def _permutation(
    source: Sequence[Sequence[str]], target: Sequence[Sequence[str]]
) -> list[int] | None:
    if len(source) != len(target):
        return None
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        _graph(source), _graph(target), edge_match=EDGE_MATCH
    )
    mapping = next(matcher.isomorphisms_iter(), None)
    return None if mapping is None else [mapping[index] for index in range(len(source))]


def _corpus_relations(matrix: Sequence[Sequence[str]]) -> tuple[tuple[str, ...], ...]:
    replacements = {
        "0": "parallel",
        "infinity": "parallel",
        "∞": "parallel",
        "1": "dotted:unspecified",
    }
    return tuple(
        tuple(
            "1" if left == right else replacements.get(str(value), str(value))
            for right, value in enumerate(row)
        )
        for left, row in enumerate(matrix)
    )


def _parse_inequality(path: Path) -> sp.Matrix:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()]
    begin = lines.index("begin")
    rows, columns, kind = lines[begin + 1].split()
    if kind != "rational":
        raise ArithmeticError("non-rational cdd input")
    payload = [
        [sp.Rational(value) for value in lines[begin + 2 + index].split()]
        for index in range(int(rows))
    ]
    if any(row[0] != 0 or len(row) != int(columns) for row in payload):
        raise ArithmeticError("unexpected inhomogeneous cdd row")
    return sp.Matrix([row[1:] for row in payload])


def _cdd_rays(path: Path, executable: str) -> tuple[tuple[int, ...], ...]:
    completed = subprocess.run(
        [executable, "--rep"],
        input=path.read_bytes(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    lines = completed.stdout.decode("utf-8", errors="strict").splitlines()
    marker = lines.index("V-representation")
    begin = lines.index("begin", marker)
    count, columns, kind = lines[begin + 1].split()
    if kind != "rational":
        raise ArithmeticError("non-rational cdd output")
    result = []
    for index in range(int(count)):
        row = [sp.Rational(value) for value in lines[begin + 2 + index].split()]
        if len(row) != int(columns) or row[0] != 0:
            raise ArithmeticError("the homogeneous cone returned a non-ray generator")
        result.append(_primitive(row[1:]))
    return tuple(result)


def _ray_archive(path: Path) -> tuple[dict[str, Any], tuple[list[Any], ...]]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        lines = [json.loads(line) for line in handle if line.strip()]
    return lines[0], tuple(lines[1:])


def _case_geometry(case: dict[str, Any], executable: str) -> dict[str, Any]:
    identifier = case["id"]
    expected = EXPECTED[identifier]
    section = case["section"]
    gram = _matrix(section["saturated_fixed_lattice"]["gram"])
    interior = sp.Matrix(
        [sp.Rational(value) for value in section["strict_fixed_interior"]["fixed_lattice_coordinates"]]
    )
    active = tuple(
        row
        for row in section["restricted_mirrors"]["records"]
        if row["hyperbolic_status"] == "active_facet"
    )
    raw_roots = tuple(
        sp.Matrix(row["primitive_integral_same_direction_coordinates"])
        for row in active
    )
    roots = tuple(
        root if _pair(gram, interior, root) > 0 else -root for root in raw_roots
    )
    facet_gram = sp.Matrix(
        [[_pair(gram, left, right) for right in roots] for left in roots]
    )
    archived_facet_gram = _matrix(section["restricted_mirrors"]["facet_gram"])
    inequality_path = ROOT / case["artifacts"]["inequality_path"]
    ray_path = ROOT / case["artifacts"]["ray_archive_path"]
    if not (
        _sha(inequality_path) == case["artifacts"]["inequality_sha256"]
        and _sha(ray_path) == case["artifacts"]["ray_archive_sha256"]
    ):
        raise ArithmeticError(f"{identifier} artifact hash changed")
    inequality = _parse_inequality(inequality_path)
    if not (
        gram.shape == (10, 10)
        and exact_inertia(gram) == (1, 9, 0)
        and len(roots) == expected["facets"]
        and facet_gram == archived_facet_gram
        and inequality == sp.Matrix.vstack(*(root.T * gram for root in roots))
        and inequality.rank() == 10
        and all(value > 0 for value in inequality * interior)
    ):
        raise ArithmeticError(f"{identifier} lattice/cone reconstruction failed")

    reflections = tuple(_reflection(gram, root) for root in roots)
    if not all(
        all(value.is_Integer for value in reflection)
        and abs(reflection.det()) == 1
        and reflection.T * gram * reflection == gram
        for reflection in reflections
    ):
        raise ArithmeticError(f"{identifier} has a nonintegral facet reflection")
    witness_map = {
        int(row["orbit_index"]): sp.Matrix(
            [sp.Rational(value) for value in row["fixed_lattice_coordinates"]]
        )
        for row in section["exact_reduction"]["strict_facet_witnesses"]
    }
    witnesses = tuple(witness_map[int(row["orbit_index"])] for row in active)
    if not all(
        _pair(gram, witnesses[target], roots[target]) == 0
        and _pair(gram, witnesses[target], witnesses[target]) > 0
        and all(
            _pair(gram, witnesses[target], roots[index]) > 0
            for index in range(len(roots))
            if index != target
        )
        for target in range(len(roots))
    ):
        raise ArithmeticError(f"{identifier} lacks a strict timelike facet witness")

    header, ray_records = _ray_archive(ray_path)
    archived_rays = tuple(tuple(map(int, row[0])) for row in ray_records)
    fresh_rays = _cdd_rays(inequality_path, executable)
    if not (
        header["case_id"] == identifier
        and header["ray_count"] == len(ray_records)
        and len(set(archived_rays)) == len(archived_rays)
        and set(fresh_rays) == set(archived_rays)
    ):
        raise ArithmeticError(f"{identifier} fresh GMP rays differ from archive")

    masks = []
    kinds = []
    pairings = []
    squares = []
    for vector_payload, active_positions, square_payload, kind in ray_records:
        vector = sp.Matrix(vector_payload)
        slacks = inequality * vector
        actual_active = [index + 1 for index, value in enumerate(slacks) if value == 0]
        square = _pair(gram, vector, vector)
        expected_kind = "finite" if square > 0 else "ideal" if square == 0 else "hyperideal"
        if not (
            all(value >= 0 for value in slacks)
            and actual_active == active_positions
            and str(square) == str(square_payload)
            and kind == expected_kind
            and sp.Matrix([inequality[index - 1, :] for index in active_positions]).rank() == 9
        ):
            raise ArithmeticError(f"{identifier} ray incidence/type replay failed")
        masks.append(sum(1 << (index - 1) for index in active_positions))
        kinds.append(kind)
        pairings.append(_pair(gram, interior, vector))
        squares.append(square)
    if not pairings or any(value <= 0 for value in pairings):
        raise ArithmeticError(f"{identifier} ray is on the wrong time sheet")

    positive_gram = -facet_gram
    enumerator = DiagramEnumerator(positive_gram, dimension=9)
    finite_masks = {mask for mask, kind in zip(masks, kinds, strict=True) if kind == "finite"}
    ideal_masks = {mask for mask, kind in zip(masks, kinds, strict=True) if kind == "ideal"}
    finite_edges = enumerator.finite_subdiagrams(8)
    finite_vertices = enumerator.finite_subdiagrams(9)
    if set(finite_vertices) != finite_masks:
        raise ArithmeticError(f"{identifier} elliptic endpoints differ from GMP")
    edge_endpoints: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for right in range(len(masks)):
        for left in range(right):
            intersection = masks[left] & masks[right]
            if intersection.bit_count() >= 8:
                edge_endpoints[intersection].append((left, right))
    if not (
        set(edge_endpoints) == set(finite_edges)
        and all(mask.bit_count() == 8 for mask in edge_endpoints)
        and all(len(pairs) == 1 for pairs in edge_endpoints.values())
    ):
        raise ArithmeticError(f"{identifier} elliptic edge endpoint gate failed")
    affine = enumerator.connected_affine_subdiagrams()
    neighborhoods = {}
    for component_mask in affine:
        component = tuple(
            index for index in range(len(roots)) if component_mask & (1 << index)
        )
        zero = tuple(
            facet
            for facet in range(len(roots))
            if facet not in component
            and all(positive_gram[facet, node] == 0 for node in component)
        )
        neighborhood = component_mask | sum(1 << index for index in zero)
        if neighborhood not in ideal_masks:
            raise ArithmeticError(f"{identifier} affine neighborhood is not a cusp")
        neighborhoods[component_mask] = neighborhood
    for cusp in ideal_masks:
        components = enumerator._components(cusp, all_nonorthogonal=True)
        if not (
            sum(component.bit_count() - 1 for component in components) == 8
            and all(component in affine for component in components)
            and all(neighborhoods[component] == cusp for component in components)
        ):
            raise ArithmeticError(f"{identifier} cusp is not a full parabolic neighborhood")

    exact = case["exact_cone"]["exact_double_description"]
    vinberg = case["literal_vinberg"]
    if not (
        Counter(kinds) == {"finite": expected["finite"], "ideal": expected["ideal"]}
        and len(finite_edges) == expected["edges"]
        and len(affine) == expected["affine"]
        and exact["hyperideal_vertex_count"] == 0
        and case["exact_cone"]["decision"]["classification"] == "finite_volume_noncompact"
        and case["time_orientation"]["all_extreme_rays_same_time_cone"] is True
        and str(min(pairings)) == case["time_orientation"]["minimum_strict_interior_pairing"]
        and str(max(pairings)) == case["time_orientation"]["maximum_strict_interior_pairing"]
        and vinberg["literal_vinberg_passes"] is True
        and vinberg["elliptic_edge_count"] == len(finite_edges)
        and vinberg["connected_affine_component_count"] == len(affine)
        and vinberg["bad_edge_count"] == vinberg["bad_parabolic_neighborhood_count"] == 0
    ):
        raise ArithmeticError(f"{identifier} archived volume summary differs from replay")
    return {
        "id": identifier,
        "fixed_lattice_determinant": int(gram.det()),
        "fixed_lattice_smith_invariants": _smith(gram),
        "facet_count": len(roots),
        "integral_reflection_count": len(reflections),
        "strict_timelike_facet_witness_count": len(witnesses),
        "fresh_GMP_extreme_ray_count": len(fresh_rays),
        "finite_vertex_count": len(finite_masks),
        "ideal_vertex_count": len(ideal_masks),
        "hyperideal_vertex_count": 0,
        "elliptic_edge_count": len(finite_edges),
        "connected_affine_component_count": len(affine),
        "time_pairing_range": [str(min(pairings)), str(max(pairings))],
        "relation_histogram": _histogram(_relations(facet_gram)),
        "all_exact_geometry_checks_pass": True,
        "_gram": gram,
        "_roots": roots,
        "_relations": _relations(facet_gram),
    }


def _overlattice(case: dict[str, Any], geometry: dict[str, Any], sk_grams: list[list[int]]) -> dict[str, Any]:
    identifier = case["id"]
    expected = EXPECTED[identifier]
    gram: sp.Matrix = geometry["_gram"]
    roots: tuple[sp.Matrix, ...] = geometry["_roots"]
    reflections = tuple(_reflection(gram, root) for root in roots)
    audit = _stable_subgroups(gram, reflections)
    groups = audit.pop("groups_internal")
    left_inverse = audit.pop("left_smith_inverse")
    if not (
        audit["stable_integral_subgroup_orders"] == expected["stable_orders"]
        and audit["stable_even_subgroup_orders"] == expected["stable_even_orders"]
        and len(audit["maximal_stable_subgroup_positions_zero_based"]) == 1
    ):
        raise ArithmeticError(f"{identifier} stable-overlattice census changed")
    position = audit["maximal_stable_subgroup_positions_zero_based"][0]
    group = groups[position]
    cosets = audit["stable_subgroups"][position]["cosets_in_smith_coordinates"]
    basis = _overlattice_basis(gram, left_inverse, cosets)
    maximal = sp.simplify(basis.T * gram * basis)
    actions = tuple(sp.simplify(basis.inv() * reflection * basis) for reflection in reflections)
    rank, a, delta = expected["two_elementary"]
    if not (
        len(group) == expected["maximal_index"]
        and abs(sp.Rational(1, basis.det())) == len(group)
        and maximal.det() == gram.det() / len(group) ** 2
        and all(maximal[index, index] % 2 == 0 for index in range(maximal.rows))
        and abs(int(maximal.det())) == 2**a
        and _q_histogram(maximal).keys() <= {"0", "1"}
        and all(
            all(value.is_Integer for value in action)
            and abs(action.det()) == 1
            and action.T * maximal * action == maximal
            for action in actions
        )
    ):
        raise ArithmeticError(f"{identifier} maximal stable lattice failed replay")

    if expected["sk_operation"] == "sign_reversal":
        normalized = -maximal
    else:
        normalized = -sp.simplify(2 * maximal.inv())
    target = sp.Matrix(10, 10, sk_grams[expected["sk_entry"] - 1])
    if not (
        all(value.is_Integer for value in normalized)
        and normalized.det() == target.det()
        and _smith(normalized) == _smith(target)
        and exact_inertia(normalized) == exact_inertia(target) == (9, 1, 0)
        and _q_histogram(normalized) == _q_histogram(target)
    ):
        raise ArithmeticError(f"{identifier} Scharlau--Kirschmer class invariant mismatch")
    sk_metadata = {
        7: {
            "genus_symbol": "II_{9,1}(2^2)",
            "full_chamber_facets": 10,
            "full_chamber_finite_vertices": 8,
            "full_chamber_cusps": 2,
        },
        23: {
            "genus_symbol": "II_{9,1}(2^4)",
            "full_chamber_facets": 11,
            "full_chamber_finite_vertices": 24,
            "full_chamber_cusps": 2,
        },
    }[expected["sk_entry"]]
    return {
        "stable_integral_subgroup_orders": audit["stable_integral_subgroup_orders"],
        "stable_even_subgroup_orders": audit["stable_even_subgroup_orders"],
        "unique_maximal_stable_overlattice": True,
        "maximal_index_over_fixed_lattice": len(group),
        "maximal_determinant": int(maximal.det()),
        "maximal_smith_invariants": _smith(maximal),
        "maximal_even_two_elementary_invariants_rank_a_delta": [rank, a, delta],
        "maximal_discriminant_q_mod_2_histogram": _q_histogram(maximal),
        "normalization_to_SK": expected["sk_operation"],
        "scharlau_kirschmer_List9_entry_one_based": expected["sk_entry"],
        **sk_metadata,
        "classification_method": (
            "exact signature, Smith group and discriminant quadratic form; "
            "uniqueness of the indefinite even 2-elementary class"
        ),
        "new_arithmetic_class": False,
    }


def _coxiter_lineage(targets: dict[str, tuple[tuple[str, ...], ...]]) -> dict[str, Any]:
    if _sha(COXITER_ZIP) != COXITER_ZIP_SHA256 or _sha(DOUBLING_TAR) != DOUBLING_TAR_SHA256:
        raise ArithmeticError("a pinned CoxIter archive changed")
    target_facets = {len(relations) for relations in targets.values()}
    ordinary_histogram: Counter[int] = Counter()
    doubling_histogram: Counter[int] = Counter()
    records: list[dict[str, Any]] = []
    with zipfile.ZipFile(COXITER_ZIP) as archive:
        ordinary_names = sorted(name for name in archive.namelist() if name.endswith(".coxiter"))
        for name in ordinary_names:
            raw = archive.read(name)
            facets, dimension, matrix = _parse_coxiter_text(raw.decode("utf-8", errors="replace"))
            if dimension == 9:
                ordinary_histogram[facets] += 1
                if facets in target_facets:
                    records.append(
                        {
                            "lane": "ordinary",
                            "path": name,
                            "sha256": sha256(raw).hexdigest(),
                            "facets": facets,
                            "relations": _corpus_relations(matrix),
                        }
                    )
    doubling_count = 0
    with tarfile.open(DOUBLING_TAR, "r:gz") as archive:
        for member in archive:
            if not member.isfile() or not member.name.endswith(".coxiter"):
                continue
            doubling_count += 1
            handle = archive.extractfile(member)
            if handle is None:
                raise ArithmeticError(f"cannot read {member.name}")
            raw = handle.read()
            facets, dimension, matrix = _parse_coxiter_text(raw.decode("utf-8", errors="replace"))
            if dimension == 9:
                doubling_histogram[facets] += 1
                if facets in target_facets:
                    records.append(
                        {
                            "lane": "index_two_doubling",
                            "path": member.name,
                            "sha256": sha256(raw).hexdigest(),
                            "facets": facets,
                            "relations": _corpus_relations(matrix),
                        }
                    )
    if not (
        len(ordinary_names) == 1214
        and doubling_count == 2959
        and dict(sorted(ordinary_histogram.items())) == {10: 38, 11: 10, 12: 13}
        and dict(sorted(doubling_histogram.items())) == {12: 8, 13: 8, 15: 8, 19: 8, 27: 8, 43: 8}
    ):
        raise ArithmeticError("pinned H9 CoxIter header inventory changed")
    comparisons = []
    total_matches = []
    for identifier, relations in targets.items():
        same = [record for record in records if record["facets"] == len(relations)]
        matches = []
        histogram_collisions = []
        for record in same:
            permutation = _permutation(relations, record["relations"])
            if permutation is not None:
                matches.append(
                    {
                        "lane": record["lane"],
                        "path": record["path"],
                        "sha256": record["sha256"],
                        "candidate_to_corpus_permutation_one_based": [value + 1 for value in permutation],
                    }
                )
            elif _histogram(relations) == _histogram(record["relations"]):
                histogram_collisions.append(record["path"])
        total_matches.extend((identifier, match) for match in matches)
        comparisons.append(
            {
                "id": identifier,
                "facet_count": len(relations),
                "same_header_record_count": len(same),
                "relation_histogram_collision_paths": histogram_collisions,
                "exact_labelled_graph_matches": matches,
            }
        )
    if not (
        len(total_matches) == 1
        and total_matches[0][0] == "A14_8_0_C01"
        and total_matches[0][1]["path"] == IMHOF_PATH
        and total_matches[0][1]["sha256"] == IMHOF_SHA256
        and total_matches[0][1]["candidate_to_corpus_permutation_one_based"]
        == IMHOF_PERMUTATION_ONE_BASED
    ):
        raise ArithmeticError("the unique Im Hof collision changed")
    return {
        "scope": "bounded exact corpus comparison, not a worldwide priority theorem",
        "ordinary_archive_sha256": COXITER_ZIP_SHA256,
        "doubling_archive_sha256": DOUBLING_TAR_SHA256,
        "ordinary_file_count": len(ordinary_names),
        "doubling_file_count": doubling_count,
        "ordinary_H9_header_histogram": {str(key): value for key, value in sorted(ordinary_histogram.items())},
        "doubling_H9_header_histogram": {str(key): value for key, value in sorted(doubling_histogram.items())},
        "comparisons": comparisons,
        "unique_exact_collision": {
            "candidate": "A14_8_0_C01",
            "published_stock": "Im Hof doubly truncated orthoscheme 9-imhof85_03",
            "path": IMHOF_PATH,
            "candidate_to_ImHof_permutation_one_based": IMHOF_PERMUTATION_ONE_BASED,
        },
        "other_four_exact_match_count": 0,
    }


def _project_lineage(targets: dict[str, tuple[tuple[str, ...], ...]]) -> dict[str, Any]:
    if _sha(FACTORY) != FACTORY_SHA256 or _sha(CUBES) != CUBES_SHA256:
        raise ArithmeticError("a pinned prior H9 project certificate changed")
    factory = json.loads(FACTORY.read_text(encoding="utf-8"))
    cubes = json.loads(CUBES.read_text(encoding="utf-8"))
    prior = []
    for construction in factory["constructions"]:
        prior.append(
            (
                construction["id"],
                _relations(_matrix(construction["section"]["restricted_mirrors"]["facet_gram"])),
            )
        )
    for record in cubes["certified_types"]:
        prior.append((record["id"], _relations(_matrix(record["facet_gram"]))))
    matches = {
        identifier: [prior_id for prior_id, relation in prior if _permutation(target, relation) is not None]
        for identifier, target in targets.items()
    }
    if any(matches.values()):
        raise ArithmeticError("a parked H9 diagram duplicates the pinned project H9 factory/cubes")
    return {
        "factory_sha256": FACTORY_SHA256,
        "doubling_cubes_sha256": CUBES_SHA256,
        "raw_factory_marking_count": len(factory["constructions"]),
        "doubling_cube_representative_count": len(cubes["certified_types"]),
        "exact_matches": matches,
        "total_exact_match_count": 0,
    }


def audit(executable: str = "cddexec_gmp") -> dict[str, Any]:
    pinned = {
        BATCH: BATCH_SHA256,
        SK_LIST9: SK_LIST9_SHA256,
        SK_GRAM9: SK_GRAM9_SHA256,
    }
    if any(_sha(path) != expected for path, expected in pinned.items()):
        raise ArithmeticError("a pinned batch or Scharlau--Kirschmer input changed")
    batch = json.loads(BATCH.read_text(encoding="utf-8"))
    cases = {case["id"]: case for case in batch["cases"]}
    if set(cases) != set(EXPECTED):
        raise ArithmeticError("the five-case denominator changed")
    sk_grams = ast.literal_eval(SK_GRAM9.read_text(encoding="utf-8"))

    geometry_internal = [_case_geometry(cases[identifier], executable) for identifier in EXPECTED]
    targets = {record["id"]: record["_relations"] for record in geometry_internal}
    within = []
    for left_index, left_id in enumerate(EXPECTED):
        for right_id in tuple(EXPECTED)[left_index + 1 :]:
            permutation = _permutation(targets[left_id], targets[right_id])
            within.append(
                {
                    "source": left_id,
                    "target": right_id,
                    "same_labelled_relation_graph": permutation is not None,
                    "source_to_target_permutation_one_based": (
                        None if permutation is None else [value + 1 for value in permutation]
                    ),
                }
            )
    if any(record["same_labelled_relation_graph"] for record in within):
        raise ArithmeticError("two parked H9 diagrams are isomorphic")

    overlattices = {
        identifier: _overlattice(
            cases[identifier],
            next(record for record in geometry_internal if record["id"] == identifier),
            sk_grams,
        )
        for identifier in EXPECTED
    }
    coxiter = _coxiter_lineage(targets)
    project = _project_lineage(targets)
    public_geometry = [
        {key: value for key, value in record.items() if not key.startswith("_")}
        for record in geometry_internal
    ]
    for record in public_geometry:
        record["reflection_stable_overlattice_lineage"] = overlattices[record["id"]]

    return {
        "schema": SCHEMA,
        "pinned_inputs": {
            "batch_path": str(BATCH.relative_to(ROOT)),
            "batch_sha256": BATCH_SHA256,
            "Scharlau_Kirschmer_List9_sha256": SK_LIST9_SHA256,
            "Scharlau_Kirschmer_Gram9_sha256": SK_GRAM9_SHA256,
        },
        "independence": {
            "imports_batch_certifier": False,
            "calls_batch_build": False,
            "fresh_GMP_double_description": True,
            "independent_ray_incidence_and_time_sheet_replay": True,
            "exact_Vinberg_edge_and_full_parabolic_neighborhood_replay": True,
            "exhaustive_reflection_stable_discriminant_subgroup_census": True,
        },
        "scope": {
            "five_finite_volume_noncompact_H9_polytopes_replayed": True,
            "one_literal_diagram_known_as_ImHof": True,
            "four_literal_diagrams_bounded_unmatched": True,
            "all_five_arithmetic_classes_known": True,
            "worldwide_novelty_claimed": False,
            "surface_automorphism_or_nef_claimed": False,
        },
        "cases": public_geometry,
        "within_batch_comparisons": within,
        "prior_project_H9_comparison": project,
        "coxiter_bounded_comparison": coxiter,
        "named_family_comparison": {
            "Vinberg_P9_and_hyperbolic_H9_simplices": (
                "10 facets, so no same-header candidate"
            ),
            "Tumarkin_H9_n_plus_2_pyramids": (
                "11 facets, so no same-header candidate; the pinned CoxIter ordinary scan includes all ten named H9/11 records"
            ),
            "Tumarkin_n_plus_3_scope": (
                "the 12-facet candidate is already identified exactly with Im Hof; Tumarkin's n+3 result is not treated as a complete H9 literal census"
            ),
            "Bredon_ADEG_P_star": (
                "14 facets; moreover ADEG requires mutually intersecting facets and labels only 2,3,6, while every candidate has a label 4"
            ),
            "Allcock_2006_H9_F_and_Q": (
                "Figure 4.1 has 11-wall F and 14-wall Q, hence no same header"
            ),
            "Allcock_2006_Vinberg_minus2_and_D7D9_H9_stocks": (
                "both have 12 facets; the minus-two stock is a pinned ordinary CoxIter exact nonmatch, and D7D9 has a dotted bond whereas A14_8_0_C01 has none"
            ),
            "historical_sources": [
                "H. C. Im Hof, Napier cycles and hyperbolic Coxeter groups",
                "P. Tumarkin, Hyperbolic Coxeter n-polytopes with n+2 and n+3 facets",
                "D. Allcock, Infinitely many hyperbolic Coxeter groups through dimension 19 (2006)",
                "N. Bredon, On ADEG-polyhedra in hyperbolic spaces (2025)",
            ],
        },
        "decision": {
            "geometry": "five_exact_finite_volume_noncompact_H9_polytopes_replay",
            "literal_diagram_lineage": (
                "one_known_ImHof_and_four_bounded_unmatched_worldwide_priority_open"
            ),
            "arithmetic_lineage": "zero_new_arithmetic_classes",
        },
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cdd", default="cddexec_gmp")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = audit(args.cdd)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if OUTPUT.read_text(encoding="utf-8") != rendered:
            raise ArithmeticError("hostile H9 lineage archive changed")
    else:
        OUTPUT.write_text(rendered, encoding="utf-8")
    print(json.dumps(payload["decision"], sort_keys=True))
    print(OUTPUT)


if __name__ == "__main__":
    main()
