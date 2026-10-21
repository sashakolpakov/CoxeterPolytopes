#!/usr/bin/env python3
"""Certify the full-NS Type-VII Coble nef chamber in characteristic five.

This is a surface-first positive certificate.  Kondō's one-boundary Type-VII
Coble surface supplies twenty complete Coble--Mukai (-2)-roots.  The exact
fifty-one ideal rays of that root chamber lift bijectively to all (-1)-curves
by ``E=f-b``, where ``b=B/2=-K``.  Thus the full Neron--Severi wall list is

    b; 20 CM roots; 51 cusp lifts,

and has 72 primitive walls in ``NS(S)=E10 + <-1> = I_(1,10)``.

The volume check is deliberately literal.  GMP cddlib supplies the complete
exact ray incidence and adjacency lattice.  A separate compiled exhaustive
clique classifier enumerates all elliptic rank-nine and rank-ten diagrams.
Finally every connected affine diagram is enumerated intrinsically and its
full zero-neighbourhood is compared with the exact ideal vertices.  This is
Vinberg's edge/endpoint and parabolic-neighbourhood criterion, not an
inference from "noncompact" or from pairwise Coxeter compatibility.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import gzip
import hashlib
from itertools import combinations, product
import json
import math
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile
from typing import Any, Iterable, Sequence

import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form
from sympy.polys.domains import ZZ

from scripts.audit_jacobian_kummer_h12_parallel_cusps import _determinant
from scripts.certify_coble_one_boundary_full_ns_h10 import (
    _calibration_walls,
    _load_templates,
    _primitive_positive,
    _root_span_and_reflections,
)
from scripts.certify_jacobian_kummer_h10_reflection_index import (
    _euler_from_connected_components,
    _relation_matrix,
)
from scripts.certify_shimada_disc16_d5_h10_lineage import (
    COXITER_COMMIT,
    _connected_elliptic_masks,
    _parse_coxiter_text,
)


ROOT = Path(__file__).resolve().parents[1]
SUMMARY_OUTPUT = (
    ROOT
    / "data/enumerations/coble_type_vii_char5_one_boundary_full_ns_h10.json"
)
RAY_OUTPUT = (
    ROOT
    / "data/enumerations/coble_type_vii_char5_one_boundary_full_ns_h10_rays.jsonl.gz"
)
EDGE_OUTPUT = (
    ROOT
    / "data/enumerations/coble_type_vii_char5_one_boundary_full_ns_h10_edges.jsonl.gz"
)
DIAGRAM_SOURCE = (
    ROOT
    / "scripts/certify_coble_type_vii_char5_one_boundary_diagrams.cpp"
)
ENRIQUES_SUMMARY = (
    ROOT / "data/known_examples/enriques_finite_automorphism_nef_chambers.json"
)
ENRIQUES_RAYS = (
    ROOT / "data/enumerations/enriques_finite_automorphism_nef_chamber_rays.jsonl.gz"
)
DEFAULT_COXITER = Path("/private/tmp/coxiter-current/CoxIter-master")

SCHEMA = "coxeter.coble_type_vii_char5_one_boundary_full_ns_h10.v1"
DIMENSION = 10
FACET_COUNT = 72
EXPECTED_RAY_HASH = "532f9c6ab31fc8bfe5223333d16e5a6b6c99bb52ba1bb17cb8ede7694b9b2e05"
P10_EULER = Fraction(-1, 11_147_673_600)
EXPECTED_EULER = Fraction(-527, 7_560)
EXPECTED_P10_INDEX = 777_093_120


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _matrix_rows(matrix: sp.Matrix) -> list[list[int]]:
    if any(value.is_Integer is not True for value in matrix):
        raise ArithmeticError("an expected integral matrix changed")
    return [[int(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)]


def _cm_overlattice_bridge(
    e10: sp.Matrix, cm_root_matrix: sp.Matrix
) -> dict[str, Any]:
    """Prove intrinsically that the Type-VII root Gram has one E10 closure."""

    root_span_basis = hermite_normal_form(cm_root_matrix)
    if root_span_basis.shape != (10, 10):
        raise ArithmeticError("the Type-VII roots did not yield a rank-ten HNF basis")
    root_lattice_gram = root_span_basis.T * e10 * root_span_basis
    smith = smith_normal_form(root_lattice_gram, domain=ZZ)
    smith_diagonal = [abs(int(smith[index, index])) for index in range(10)]
    if root_lattice_gram.det() != -4 or smith_diagonal != [1] * 8 + [2, 2]:
        raise ArithmeticError("the abstract Type-VII root lattice discriminant changed")

    discriminant_classes = []
    even_nonzero = []
    for bits in product(range(2), repeat=10):
        half_vector = sp.Matrix([sp.Rational(bit, 2) for bit in bits])
        if not all(value.is_Integer is True for value in root_lattice_gram * half_vector):
            continue
        square = sp.Rational((half_vector.T * root_lattice_gram * half_vector)[0])
        record = {
            "twice_representative_in_HNF_root_basis": list(bits),
            "q_value": str(square),
            "zero": not any(bits),
            "even_isotropic_mod_2Z": square.q == 1 and int(square) % 2 == 0,
        }
        discriminant_classes.append(record)
        if not record["zero"] and record["even_isotropic_mod_2Z"]:
            even_nonzero.append((bits, half_vector))
    if len(discriminant_classes) != 4 or len(even_nonzero) != 1:
        raise ArithmeticError("the unique even Type-VII discriminant glue changed")
    nonzero_q_values = {
        record["q_value"] for record in discriminant_classes if not record["zero"]
    }
    if nonzero_q_values != {"-5/2", "-2", "-7/2"}:
        raise ArithmeticError("the Type-VII discriminant q-values changed")

    bits, half_vector = even_nonzero[0]
    pivot = next(index for index, bit in enumerate(bits) if bit)
    extension = sp.eye(10)
    extension[:, pivot] = half_vector
    overbasis = root_span_basis * extension
    overgram = overbasis.T * e10 * overbasis
    if not (
        all(value.is_Integer is True for value in overbasis)
        and abs(int(overbasis.det())) == 1
        and overgram.det() == -1
        and all(int(overgram[index, index]) % 2 == 0 for index in range(10))
    ):
        raise ArithmeticError("the unique even glue did not produce an E10 overclosure")
    return {
        "root_lattice_HNF_basis_in_E10": _matrix_rows(root_span_basis),
        "root_lattice_gram": _matrix_rows(root_lattice_gram),
        "root_lattice_determinant": int(root_lattice_gram.det()),
        "root_lattice_discriminant_smith_diagonal": smith_diagonal,
        "discriminant_group": "(Z/2)^2",
        "discriminant_classes": discriminant_classes,
        "nonzero_q_value_set": sorted(nonzero_q_values),
        "unique_nonzero_even_isotropic_class_count": 1,
        "unique_even_overlattice_basis_in_E10": _matrix_rows(overbasis),
        "unique_even_overlattice_basis_determinant": int(overbasis.det()),
        "unique_even_overlattice_gram": _matrix_rows(overgram),
        "unique_even_overlattice_is_even_unimodular_E10": True,
        "transport_conclusion": "every entrywise Gamma_VII root-Gram isometry extends across the unique even index-two overclosure to an E10 isometry",
        "marking_uniqueness_not_claimed": True,
    }


def _mask(indices: Iterable[int]) -> int:
    return sum(1 << index for index in indices)


def _indices(mask: int) -> tuple[int, ...]:
    result = []
    while mask:
        bit = mask & -mask
        result.append(bit.bit_length() - 1)
        mask ^= bit
    return tuple(result)


def _quadratic(vector: Sequence[int], gram: Sequence[Sequence[int]]) -> int:
    return sum(
        vector[left] * gram[left][right] * vector[right]
        for left in range(len(vector))
        for right in range(len(vector))
    )


def _bilinear(
    left_vector: Sequence[int],
    gram: Sequence[Sequence[int]],
    right_vector: Sequence[int],
) -> int:
    return sum(
        left_vector[left] * gram[left][right] * right_vector[right]
        for left in range(len(left_vector))
        for right in range(len(right_vector))
    )


def _cdd_input(inequalities: sp.Matrix) -> str:
    rows = _matrix_rows(inequalities)
    lines = [
        "H-representation",
        "begin",
        f"{len(rows)} {inequalities.cols + 1} rational",
    ]
    lines.extend("0 " + " ".join(map(str, row)) for row in rows)
    lines.append("end")
    return "\n".join(lines) + "\n"


def _parse_index_block(
    lines: Sequence[str], heading: str
) -> tuple[tuple[int, ...], ...]:
    try:
        heading_position = next(
            index for index, line in enumerate(lines) if line.strip() == heading
        )
        begin = next(
            index
            for index in range(heading_position + 1, len(lines))
            if lines[index].strip() == "begin"
        )
    except StopIteration as error:
        available = [
            line.strip()
            for line in lines
            if "incidence" in line.lower() or "adjacency" in line.lower()
        ]
        raise ArithmeticError(f"cddlib omitted {heading!r}; saw {available}") from error
    row_count, universe = map(int, lines[begin + 1].split()[:2])
    records: list[tuple[int, ...] | None] = [None] * row_count
    for offset in range(row_count):
        line = lines[begin + 2 + offset]
        left, separator, right = line.partition(":")
        if not separator:
            raise ArithmeticError(f"malformed cddlib {heading} row")
        row_number, signed_count = map(int, left.split()[:2])
        listed = {int(token) - 1 for token in right.split()}
        if any(index < 0 or index >= universe for index in listed):
            raise ArithmeticError(f"out-of-range cddlib {heading} entry")
        values = set(range(universe)) - listed if signed_count < 0 else listed
        if len(values) != abs(signed_count):
            raise ArithmeticError(f"cddlib {heading} complement count failed")
        if not 1 <= row_number <= row_count or records[row_number - 1] is not None:
            raise ArithmeticError(f"cddlib {heading} row numbering failed")
        records[row_number - 1] = tuple(sorted(values))
    if lines[begin + 2 + row_count].strip() != "end" or any(
        record is None for record in records
    ):
        raise ArithmeticError(f"cddlib {heading} block did not terminate cleanly")
    return tuple(record for record in records if record is not None)


def _run_cdd_repall(
    inequalities: sp.Matrix, executable: str
) -> tuple[list[tuple[int, ...]], tuple[tuple[int, ...], ...], tuple[tuple[int, ...], ...]]:
    completed = subprocess.run(
        [executable, "--repall"],
        input=_cdd_input(inequalities),
        text=True,
        capture_output=True,
        check=False,
        timeout=600,
    )
    if completed.returncode:
        raise RuntimeError(completed.stderr[-4000:] or completed.stdout[-4000:])
    lines = completed.stdout.splitlines()
    try:
        representation = lines.index("V-representation")
        begin = next(
            index
            for index in range(representation + 1, len(lines))
            if lines[index].strip() == "begin"
        )
    except (ValueError, StopIteration) as error:
        raise ArithmeticError("cddlib omitted its V-representation") from error
    count, width, kind = lines[begin + 1].split()[:3]
    if kind != "rational" or int(width) != inequalities.cols + 1:
        raise ArithmeticError("cddlib returned an unexpected generator format")
    rays = []
    for line in lines[begin + 2 : begin + 2 + int(count)]:
        values = [Fraction(token) for token in line.split()]
        if len(values) != int(width) or values[0] != 0:
            raise ArithmeticError("the homogeneous cone acquired a vertex row")
        rays.append(_primitive_positive(values[1:]))
    if lines[begin + 2 + int(count)].strip() != "end":
        raise ArithmeticError("cddlib V-representation did not terminate cleanly")
    incidence = _parse_index_block(lines, "Vertex incidence")
    adjacency = _parse_index_block(lines, "Vertex adjacency")
    if len(rays) != len(incidence) or len(rays) != len(adjacency):
        raise ArithmeticError("cddlib ray/incidence/adjacency counts disagree")
    if len(rays) != len(set(rays)):
        raise ArithmeticError("cddlib returned duplicate oriented rays")
    return rays, incidence, adjacency


def _relation_payload(positive_gram: Sequence[Sequence[int]]) -> tuple[list[list[str]], list[list[int]]]:
    relations = _relation_matrix(positive_gram)
    numeric = [[1 if left == right else 0 for right in range(FACET_COUNT)] for left in range(FACET_COUNT)]
    supported = {"2": 2, "3": 3, "4": 4, "6": 6}
    for left in range(FACET_COUNT):
        for right in range(left):
            label = str(relations[left][right])
            value = supported.get(label, 0)
            numeric[left][right] = numeric[right][left] = value
            # The shared Solomon helper writes 0 for a parallel pair and 1
            # for an ultraparallel pair.
            if label not in {*supported, "0", "1"}:
                raise ArithmeticError(f"unsupported Coxeter relation {label}")
    return relations, numeric


def _compile_diagram_enumerator(compiler: str | None = None) -> tuple[tempfile.TemporaryDirectory[str], Path, str]:
    selected = compiler or shutil.which("c++")
    if selected is None:
        raise RuntimeError("a C++17 compiler is required for the exhaustive diagram pass")
    temporary = tempfile.TemporaryDirectory(prefix="coble-vii-diagrams-")
    executable = Path(temporary.name) / "enumerate"
    completed = subprocess.run(
        [selected, "-O3", "-std=c++17", str(DIAGRAM_SOURCE), "-o", str(executable)],
        text=True,
        capture_output=True,
        check=False,
        timeout=120,
    )
    if completed.returncode:
        temporary.cleanup()
        raise RuntimeError(completed.stderr[-4000:] or completed.stdout[-4000:])
    return temporary, executable, selected


def _finite_diagrams(
    numeric_relations: Sequence[Sequence[int]],
    executable: Path,
) -> tuple[dict[int, set[int]], dict[int, int]]:
    text = (
        f"{FACET_COUNT} {DIMENSION - 1} {DIMENSION}\n"
        + "\n".join(" ".join(map(str, row)) for row in numeric_relations)
        + "\n"
    )
    completed = subprocess.run(
        [str(executable)],
        input=text,
        text=True,
        capture_output=True,
        check=False,
        timeout=600,
    )
    if completed.returncode:
        raise RuntimeError(completed.stderr[-4000:] or completed.stdout[-4000:])
    finite = {DIMENSION - 1: set(), DIMENSION: set()}
    clique_counts: dict[int, int] = {}
    finite_summary: dict[int, int] = {}
    for line in completed.stdout.splitlines():
        tokens = line.split()
        if not tokens:
            continue
        if tokens[0] == "SUMMARY":
            order, cliques, count = map(int, tokens[1:])
            clique_counts[order] = cliques
            finite_summary[order] = count
            continue
        order, low, high = map(int, tokens)
        if order not in finite:
            raise ArithmeticError("diagram enumerator returned an unexpected order")
        finite[order].add(low | (high << 64))
    if any(len(finite[order]) != finite_summary.get(order) for order in finite):
        raise ArithmeticError("diagram enumerator output count or uniqueness failed")
    return finite, clique_counts


def _components(mask: int, nonorthogonal: Sequence[int]) -> tuple[int, ...]:
    remaining = mask
    result = []
    while remaining:
        start = remaining & -remaining
        remaining ^= start
        component = 0
        frontier = start
        while frontier:
            bit = frontier & -frontier
            frontier ^= bit
            node = bit.bit_length() - 1
            component |= bit
            fresh = nonorthogonal[node] & remaining
            remaining ^= fresh
            frontier |= fresh
        result.append(component)
    return tuple(result)


def _finite_connected_type(mask: int, numeric_relations: Sequence[Sequence[int]]) -> str | None:
    nodes = _indices(mask)
    order = len(nodes)
    if order == 1:
        return "A1"
    local = {node: [] for node in nodes}
    edges = []
    for position, left in enumerate(nodes):
        for right in nodes[position + 1 :]:
            label = numeric_relations[left][right]
            if label in {3, 4, 6}:
                local[left].append(right)
                local[right].append(left)
                edges.append((left, right, label))
    if len(edges) != order - 1:
        return None
    degrees = {node: len(local[node]) for node in nodes}
    six = [(left, right) for left, right, label in edges if label == 6]
    four = [(left, right) for left, right, label in edges if label == 4]
    if six:
        return "G2" if order == 2 and len(six) == 1 else None
    if max(degrees.values()) <= 2:
        if not four:
            return f"A{order}"
        if len(four) != 1:
            return None
        left, right = four[0]
        if degrees[left] == 1 or degrees[right] == 1:
            return f"B{order}"
        return "F4" if order == 4 else None
    if four or max(degrees.values()) > 3:
        return None
    branch = [node for node in nodes if degrees[node] == 3]
    if len(branch) != 1:
        return None
    branch_node = branch[0]
    arms = []
    for neighbour in local[branch_node]:
        previous, current, length = branch_node, neighbour, 1
        while degrees[current] == 2:
            following = next(node for node in local[current] if node != previous)
            previous, current = current, following
            length += 1
        arms.append(length)
    arm_tuple = tuple(sorted(arms))
    if arm_tuple[:2] == (1, 1):
        return f"D{order}"
    return {
        (1, 2, 2): "E6",
        (1, 2, 3): "E7",
        (1, 2, 4): "E8",
    }.get(arm_tuple)


def _finite_signature(
    mask: int,
    numeric_relations: Sequence[Sequence[int]],
    finite_nonorthogonal: Sequence[int],
    cache: dict[int, str | None],
) -> tuple[str, ...] | None:
    types = []
    for component in _components(mask, finite_nonorthogonal):
        if component not in cache:
            cache[component] = _finite_connected_type(component, numeric_relations)
        kind = cache[component]
        if kind is None:
            return None
        types.append(kind)
    return tuple(sorted(types))


def _connected_affine_masks(
    positive_gram: Sequence[Sequence[int]],
    connected_finite: set[int],
) -> tuple[dict[int, dict[str, int]], dict[int, int]]:
    size = len(positive_gram)
    compatible = []
    nonorthogonal = []
    for left in range(size):
        compatible_mask = 0
        nonorthogonal_mask = 0
        for right in range(size):
            if left == right:
                continue
            pairing = positive_gram[left][right]
            if pairing:
                nonorthogonal_mask |= 1 << right
            if pairing == 0 or (
                positive_gram[left][left] * positive_gram[right][right]
                > pairing * pairing
            ):
                compatible_mask |= 1 << right
        compatible.append(compatible_mask)
        nonorthogonal.append(nonorthogonal_mask)

    affine: dict[int, dict[str, int]] = {}
    for right in range(size):
        for left in range(right):
            if (
                positive_gram[left][left] * positive_gram[right][right]
                == positive_gram[left][right] * positive_gram[left][right]
            ):
                mask = (1 << left) | (1 << right)
                affine[mask] = {"order": 2, "coxeter_rank": 1}

    candidate_counts = {2: len(affine)}
    finite_by_order: dict[int, list[int]] = defaultdict(list)
    for mask in connected_finite:
        finite_by_order[mask.bit_count()].append(mask)
    universe = (1 << size) - 1
    for order in range(3, DIMENSION + 3):
        candidates = set()
        for base in finite_by_order.get(order - 1, []):
            common = universe
            touching = 0
            for node in _indices(base):
                common &= compatible[node]
                touching |= nonorthogonal[node]
            additions = common & touching & ~base
            while additions:
                bit = additions & -additions
                additions ^= bit
                candidates.add(base | bit)
        candidate_counts[order] = len(candidates)
        for mask in candidates:
            nodes = _indices(mask)
            matrix = [[positive_gram[left][right] for right in nodes] for left in nodes]
            # The generating deletion is positive definite of order-1.  Thus
            # det=0 is already an exact PSD-corank-one (affine) certificate.
            if _determinant(matrix) == 0:
                affine[mask] = {"order": order, "coxeter_rank": order - 1}
    return affine, candidate_counts


def _face_and_vinberg_audit(
    rays: Sequence[tuple[int, ...]],
    incidence: Sequence[tuple[int, ...]],
    adjacency: Sequence[tuple[int, ...]],
    inequalities: sp.Matrix,
    basis_gram: sp.Matrix,
    positive_gram: list[list[int]],
    numeric_relations: list[list[int]],
    finite_diagrams: dict[int, set[int]],
    connected_finite: set[int],
) -> tuple[dict[str, Any], list[dict[str, Any]], list[list[Any]]]:
    inequality_rows = _matrix_rows(inequalities)
    basis_rows = _matrix_rows(basis_gram)
    finite_nonorthogonal = [
        _mask(
            right
            for right in range(FACET_COUNT)
            if numeric_relations[left][right] in {3, 4, 6}
        )
        for left in range(FACET_COUNT)
    ]
    all_nonorthogonal = [
        _mask(
            right
            for right in range(FACET_COUNT)
            if right != left and positive_gram[left][right] != 0
        )
        for left in range(FACET_COUNT)
    ]

    kinds = []
    squares = []
    active_masks = []
    face_bases: list[list[list[int]]] = [[] for _ in range(FACET_COUNT)]
    ray_records: list[dict[str, Any]] = []
    for ray_index, (ray, active) in enumerate(zip(rays, incidence, strict=True)):
        slacks = [
            sum(row[column] * ray[column] for column in range(len(ray)))
            for row in inequality_rows
        ]
        if any(value < 0 for value in slacks):
            raise ArithmeticError("an exact cdd ray violates a wall inequality")
        zeros = tuple(index for index, value in enumerate(slacks) if value == 0)
        if zeros != active:
            raise ArithmeticError("cdd incidence disagrees with exact zero slacks")
        square = _quadratic(ray, basis_rows)
        kind = "finite" if square > 0 else "ideal" if square == 0 else "hyperideal"
        mask = _mask(active)
        kinds.append(kind)
        squares.append(square)
        active_masks.append(mask)
        for facet in active:
            if len(face_bases[facet]) == DIMENSION:
                continue
            trial = [*face_bases[facet], list(ray)]
            if sp.Matrix(trial).rank() > len(face_bases[facet]):
                face_bases[facet].append(list(ray))
        ray_records.append(
            {
                "ray_index": ray_index + 1,
                "primitive_coordinate_basis_vector": list(ray),
                "active_facets": [index + 1 for index in active],
                "square": square,
                "kind": kind,
            }
        )

    kind_histogram = Counter(kinds)
    incidence_histogram = Counter(len(active) for active in incidence)
    incidence_by_kind = {
        kind: {
            str(order): count
            for order, count in sorted(
                Counter(
                    len(incidence[index])
                    for index, candidate in enumerate(kinds)
                    if candidate == kind
                ).items()
            )
        }
        for kind in ("finite", "ideal", "hyperideal")
    }

    edge_pairs: dict[int, tuple[int, int]] = {}
    adjacency_pair_count = 0
    for left, neighbours in enumerate(adjacency):
        if left in neighbours:
            raise ArithmeticError("cdd vertex adjacency contains a loop")
        for right in neighbours:
            if left not in adjacency[right]:
                raise ArithmeticError("cdd vertex adjacency is not symmetric")
            if right <= left:
                continue
            adjacency_pair_count += 1
            support = active_masks[left] & active_masks[right]
            if support.bit_count() != DIMENSION - 1 or support in edge_pairs:
                raise ArithmeticError("cdd edge support is not a unique rank-nine set")
            edge_pairs[support] = (left, right)

    finite_rank_9 = finite_diagrams[DIMENSION - 1]
    finite_rank_10 = finite_diagrams[DIMENSION]
    finite_endpoint_masks = {
        active_masks[index] for index, kind in enumerate(kinds) if kind == "finite"
    }
    ideal_endpoint_masks = {
        active_masks[index] for index, kind in enumerate(kinds) if kind == "ideal"
    }
    if len(finite_endpoint_masks) != kind_histogram["finite"]:
        raise ArithmeticError("finite cdd endpoint incidences are not unique")
    if len(ideal_endpoint_masks) != kind_histogram["ideal"]:
        raise ArithmeticError("ideal cdd endpoint incidences are not unique")
    if set(edge_pairs) != finite_rank_9:
        raise ArithmeticError("elliptic rank-nine diagrams do not equal exact edges")
    if finite_endpoint_masks != finite_rank_10:
        raise ArithmeticError("elliptic rank-ten diagrams do not equal finite endpoints")

    diagram_extensions: dict[int, list[int]] = defaultdict(list)
    for vertex, active in enumerate(active_masks):
        active_indices = _indices(active)
        for subset in combinations(active_indices, DIMENSION - 1):
            edge = _mask(subset)
            if edge in finite_rank_9:
                diagram_extensions[edge].append(vertex)
    if set(diagram_extensions) != finite_rank_9:
        raise ArithmeticError("an elliptic edge has no endpoint extension")

    endpoint_histogram: Counter[str] = Counter()
    edge_type_histogram: Counter[str] = Counter()
    signature_cache: dict[int, str | None] = {}
    edge_records: list[list[Any]] = []
    for edge in sorted(finite_rank_9):
        extensions = tuple(sorted(set(diagram_extensions[edge])))
        if len(extensions) != 2 or extensions != edge_pairs[edge]:
            raise ArithmeticError("an elliptic rank-nine edge does not have exactly two endpoints")
        endpoint_code = "".join(sorted(kinds[index][0].upper() for index in extensions))
        endpoint_histogram[endpoint_code] += 1
        signature = _finite_signature(
            edge, numeric_relations, finite_nonorthogonal, signature_cache
        )
        if signature is None:
            raise ArithmeticError("a purported elliptic edge has nonfinite type")
        signature_text = " + ".join(signature)
        edge_type_histogram[signature_text] += 1
        edge_records.append(
            [
                [index + 1 for index in _indices(edge)],
                [index + 1 for index in extensions],
                [kinds[index] for index in extensions],
                list(signature),
            ]
        )

    affine, affine_candidate_counts = _connected_affine_masks(
        positive_gram, connected_finite
    )
    ideal_by_mask = {
        active_masks[index]: index
        for index, kind in enumerate(kinds)
        if kind == "ideal"
    }
    affine_neighbourhoods = {}
    for component, record in affine.items():
        nodes = _indices(component)
        zero_neighbours = _mask(
            facet
            for facet in range(FACET_COUNT)
            if facet not in nodes
            and all(positive_gram[facet][node] == 0 for node in nodes)
        )
        neighbourhood = component | zero_neighbours
        if neighbourhood not in ideal_by_mask:
            raise ArithmeticError("an affine component has a non-cusp full zero-neighbourhood")
        affine_neighbourhoods[component] = neighbourhood
        record["full_zero_neighbourhood_ideal_ray_index"] = ideal_by_mask[neighbourhood] + 1

    cusp_component_orders: Counter[tuple[int, ...]] = Counter()
    cusp_component_counts: Counter[int] = Counter()
    for cusp in ideal_endpoint_masks:
        components = _components(cusp, all_nonorthogonal)
        if not (
            sum(component.bit_count() - 1 for component in components)
            == DIMENSION - 1
            and all(component in affine for component in components)
            and all(affine_neighbourhoods[component] == cusp for component in components)
        ):
            raise ArithmeticError("an ideal endpoint is not a full rank-nine parabolic neighbourhood")
        orders = tuple(sorted(component.bit_count() for component in components))
        cusp_component_orders[orders] += 1
        cusp_component_counts[len(components)] += 1

    interior = tuple(sum(ray[column] for ray in rays) for column in range(len(rays[0])))
    divisor = math.gcd(*(abs(value) for value in interior))
    primitive_interior = tuple(value // divisor for value in interior)
    interior_slacks = [
        sum(row[column] * primitive_interior[column] for column in range(len(interior)))
        for row in inequality_rows
    ]
    interior_square = _quadratic(primitive_interior, basis_rows)
    if min(interior_slacks) <= 0 or interior_square <= 0:
        raise ArithmeticError("the sum of exact rays did not give a strict timelike interior point")
    time_orientation_pairings = [
        _bilinear(primitive_interior, basis_rows, ray) for ray in rays
    ]
    if min(time_orientation_pairings) <= 0:
        raise ArithmeticError("the causal cdd rays do not lie on one Lorentz time-sheet")

    cone = {
        "extreme_ray_count": len(rays),
        "finite_vertex_count": kind_histogram["finite"],
        "ideal_vertex_count": kind_histogram["ideal"],
        "hyperideal_vertex_count": kind_histogram["hyperideal"],
        "vertex_incidence_size_histogram": {
            str(order): count for order, count in sorted(incidence_histogram.items())
        },
        "vertex_incidence_size_histogram_by_kind": incidence_by_kind,
        "all_facet_face_ranks": [len(rows) for rows in face_bases],
        "all_facets_irredundant": all(len(rows) == DIMENSION for rows in face_bases),
        "primitive_ray_list_sha256": hashlib.sha256(
            json.dumps(sorted(rays), separators=(",", ":")).encode()
        ).hexdigest(),
        "strict_interior_coordinate_basis_vector": list(primitive_interior),
        "strict_interior_square": interior_square,
        "strict_interior_minimum_wall_slack": min(interior_slacks),
        "minimum_strict_interior_to_extreme_ray_pairing": min(
            time_orientation_pairings
        ),
        "all_extreme_rays_on_one_lorentz_time_sheet": True,
    }
    literal = {
        "criterion": "Vinberg elliptic-edge endpoint and full affine zero-neighbourhood criterion",
        "exact_face_engine": "cddexec_gmp --repall over Q",
        "all_pair_compatible_rank_9_cliques": 15_619_250,
        "all_pair_compatible_rank_10_cliques": 24_466_316,
        "elliptic_rank_9_edge_diagram_count": len(finite_rank_9),
        "elliptic_rank_10_finite_endpoint_diagram_count": len(finite_rank_10),
        "exact_geometric_edge_count": adjacency_pair_count,
        "edge_endpoint_count_histogram": {"2": len(edge_pairs)},
        "edge_endpoint_type_histogram": dict(sorted(endpoint_histogram.items())),
        "edge_finite_type_histogram": dict(sorted(edge_type_histogram.items())),
        "connected_affine_component_count": len(affine),
        "connected_affine_order_histogram": {
            str(order): count
            for order, count in sorted(
                Counter(record["order"] for record in affine.values()).items()
            )
        },
        "connected_affine_candidate_count_by_order": {
            str(order): count for order, count in sorted(affine_candidate_counts.items())
        },
        "distinct_full_parabolic_neighbourhood_count": len(set(affine_neighbourhoods.values())),
        "cusp_affine_component_count_histogram": {
            str(count): value for count, value in sorted(cusp_component_counts.items())
        },
        "cusp_component_order_histogram": {
            "+".join(map(str, orders)): count
            for orders, count in sorted(cusp_component_orders.items())
        },
        "bad_edge_count": 0,
        "bad_parabolic_neighbourhood_count": 0,
        "every_elliptic_rank_9_diagram_is_exactly_one_geometric_edge": True,
        "every_edge_has_exactly_two_finite_or_ideal_endpoints": True,
        "every_connected_affine_full_zero_neighbourhood_is_a_rank_9_cusp": True,
        "every_ideal_vertex_is_componentwise_affine_of_total_rank_9": True,
    }
    return cone | {"literal_vinberg": literal}, ray_records, edge_records


def _gzip_jsonl(header: dict[str, Any], records: Iterable[Any]) -> bytes:
    lines = [json.dumps(header, sort_keys=True, separators=(",", ":"))]
    lines.extend(json.dumps(record, separators=(",", ":")) for record in records)
    return gzip.compress(("\n".join(lines) + "\n").encode(), compresslevel=9, mtime=0)


def _bounded_coxiter_header_audit(corpus: Path | None) -> dict[str, Any]:
    if corpus is None or not corpus.is_dir():
        return {
            "available": False,
            "status": "not replayed; worldwide literal priority remains open",
        }
    ordinary_count = 0
    ordinary_matches = []
    for path in sorted(corpus.rglob("*.coxiter")):
        facets, dimension, _ = _parse_coxiter_text(path.read_text(encoding="utf-8"))
        ordinary_count += 1
        if (facets, dimension) == (FACET_COUNT, DIMENSION):
            ordinary_matches.append(path.relative_to(corpus).as_posix())
    archive = corpus / "graphs/testing_index2.tar.gz"
    doubling_count = 0
    doubling_matches = []
    with tarfile.open(archive, "r:gz") as handle:
        for member in handle:
            if not member.isfile() or not member.name.endswith(".coxiter"):
                continue
            extracted = handle.extractfile(member)
            if extracted is None:
                raise ArithmeticError("failed to read a CoxIter archive member")
            facets, dimension, _ = _parse_coxiter_text(
                extracted.read().decode("utf-8")
            )
            doubling_count += 1
            if (facets, dimension) == (FACET_COUNT, DIMENSION):
                doubling_matches.append(member.name)
    return {
        "available": True,
        "repository": "https://github.com/rgugl/CoxIter",
        "commit": COXITER_COMMIT,
        "ordinary_file_count": ordinary_count,
        "ordinary_H10_72_header_matches": ordinary_matches,
        "doubling_archive_sha256": _sha256(archive),
        "doubling_file_count": doubling_count,
        "doubling_H10_72_header_matches": doubling_matches,
        "bounded_exclusion_passes": not ordinary_matches and not doubling_matches,
        "scope": "header exclusion in this pinned corpus, not worldwide priority",
    }


def build_certificate(
    cdd_executable: str | None = None,
    compiler: str | None = None,
    coxiter_corpus: Path | None = DEFAULT_COXITER,
) -> tuple[dict[str, Any], bytes, bytes]:
    cdd = cdd_executable or shutil.which("cddexec_gmp")
    if cdd is None:
        raise RuntimeError("cddexec_gmp is required for the exact face replay")
    chamber_data, ray_data = _load_templates()
    lattice, walls, construction = _calibration_walls(
        "VII", chamber_data, ray_data
    )
    if construction["cm_root_count"] != 20 or construction["cm_cusp_count"] != 51:
        raise ArithmeticError("the pinned Type-VII CM input changed")
    wall_matrix = sp.Matrix.hstack(*walls).T
    negative_gram = wall_matrix * lattice * wall_matrix.T
    type_vii_template = next(
        record
        for record in chamber_data["types"]
        if record["nikulin_kondo_type"] == "VII"
    )
    if negative_gram[1:21, 1:21] != sp.Matrix(
        type_vii_template["intersection_matrix"]
    ):
        raise ArithmeticError("the Type-VII Coble/20C root Gram bridge changed")
    cm_root_matrix = sp.Matrix.hstack(*walls[1:21])[:10, :]
    cm_smith = smith_normal_form(cm_root_matrix, domain=ZZ)
    cm_smith_diagonal = [abs(int(cm_smith[index, index])) for index in range(10)]
    if cm_smith_diagonal != [1] * 9 + [2]:
        raise ArithmeticError("the Type-VII root embedding no longer has index two in E10")
    cm_overlattice_bridge = _cm_overlattice_bridge(lattice[:10, :10], cm_root_matrix)
    pivots = list(negative_gram.rref()[1])
    if len(pivots) != DIMENSION + 1 or len(walls) != FACET_COUNT:
        raise ArithmeticError("the Type-VII full wall rank/count changed")
    inequalities = negative_gram[:, pivots]
    basis_gram = negative_gram.extract(pivots, pivots)
    rays, incidence, adjacency = _run_cdd_repall(inequalities, cdd)

    positive_gram = [
        [-int(negative_gram[left, right]) for right in range(FACET_COUNT)]
        for left in range(FACET_COUNT)
    ]
    relations, numeric_relations = _relation_payload(positive_gram)
    connected_finite, connected_search = _connected_elliptic_masks(positive_gram)
    solomon = _euler_from_connected_components(relations, connected_finite)
    euler = Fraction(solomon["rational_euler_characteristic"])
    if euler != EXPECTED_EULER or euler / P10_EULER != EXPECTED_P10_INDEX:
        raise ArithmeticError("the Solomon/P10 index calculation changed")

    temporary, diagram_executable, compiler_used = _compile_diagram_enumerator(compiler)
    try:
        finite_diagrams, clique_counts = _finite_diagrams(
            numeric_relations, diagram_executable
        )
    finally:
        temporary.cleanup()
    expected_cliques = {9: 15_619_250, 10: 24_466_316}
    if clique_counts != expected_cliques:
        raise ArithmeticError("the exhaustive compatible-clique census changed")
    if not (
        len(finite_diagrams[9]) == 143_565
        and len(finite_diagrams[10]) == 18_761
        and solomon["finite_special_subset_count_by_rank"][9:] == [143_565, 18_761]
    ):
        raise ArithmeticError("independent finite-diagram enumerations disagree")

    face_payload, ray_records, edge_records = _face_and_vinberg_audit(
        rays,
        incidence,
        adjacency,
        inequalities,
        basis_gram,
        positive_gram,
        numeric_relations,
        finite_diagrams,
        connected_finite,
    )
    cone = {key: value for key, value in face_payload.items() if key != "literal_vinberg"}
    literal = face_payload["literal_vinberg"]
    arithmetic = _root_span_and_reflections(lattice, walls)
    # Keep the in-memory certificate identical to its JSON round trip.
    arithmetic["norm_histogram"] = {
        str(norm): count for norm, count in arithmetic["norm_histogram"].items()
    }
    if not (
        cone["finite_vertex_count"] == 18_761
        and cone["ideal_vertex_count"] == 1_731
        and cone["hyperideal_vertex_count"] == 0
        and cone["all_facets_irredundant"]
        and cone["primitive_ray_list_sha256"] == EXPECTED_RAY_HASH
        and literal["exact_geometric_edge_count"] == 143_565
        and literal["edge_endpoint_type_histogram"]
        == {"FF": 60_250, "FI": 67_110, "II": 16_205}
        and literal["connected_affine_component_count"] == 4_943
        and literal["distinct_full_parabolic_neighbourhood_count"] == 1_731
        and arithmetic["roots_span_whole_lattice"]
        and arithmetic["all_reflections_integral_involutive_isometries"]
        and arithmetic["all_pairs_coxeter"]
    ):
        raise ArithmeticError("the Type-VII full-NS chamber census changed")

    ray_header = {
        "schema": "coxeter.coble_type_vii_char5_one_boundary_full_ns_h10_rays.v1",
        "record_format": {
            "ray_index": "one based",
            "primitive_coordinate_basis_vector": "coordinates in the 11 pivot wall normals",
            "active_facets": "one-based complete incidence set",
            "square": "exact square in the pivot-wall Gram",
            "kind": "finite, ideal, or hyperideal",
        },
        "coordinate_basis_facets": [index + 1 for index in pivots],
        "ray_count": len(ray_records),
    }
    edge_header = {
        "schema": "coxeter.coble_type_vii_char5_one_boundary_full_ns_h10_edges.v1",
        "record_format": "[edge_facets, endpoint_ray_indices, endpoint_kinds, finite_components]",
        "indexing": "one based",
        "edge_count": len(edge_records),
    }
    ray_archive = _gzip_jsonl(ray_header, ray_records)
    edge_archive = _gzip_jsonl(edge_header, edge_records)

    bounded_lineage = _bounded_coxiter_header_audit(coxiter_corpus)
    payload = {
        "schema": SCHEMA,
        "status": "proved exact complete-nef finite-volume noncompact arithmetic Coxeter chamber",
        "direction": "SURFACE_TO_POLYTOPE",
        "inverse_carrier_score": {
            "level": "N/A/L0",
            "reason": "the construction starts from Kondō's surface and makes no polytope-to-surface implication",
        },
        "surface": {
            "classification": "Kondo Type VII strong Coble surface",
            "base_field": "algebraically closed, characteristic 5",
            "uniqueness": "unique isomorphism class in Kondo Theorem 1.1/Table 2",
            "boundary_components": 1,
            "automorphism_group": "S5",
            "effective_irreducible_CM_roots": 20,
            "construction": "Kondo Example 4.15 and Proposition 4.16",
        },
        "algebro_geometric_completeness": {
            "integral_CM_marking": "Kondo Proposition 3.1 gives CM(S)=E10; Example 4.15 gives the same Gamma_VII root graph and Proposition 4.16 a compatible (A9+A1,Z/2) glue invariant, supplying an E10 marking to the Shimada 20C index-two embedding; uniqueness of the marking is not claimed",
            "CM_root_span_smith_diagonal": cm_smith_diagonal,
            "CM_root_span_index_in_E10": 2,
            "entrywise_Coble_Gamma_VII_equals_Shimada_20C_Gram": True,
            "intrinsic_unique_even_E10_overlattice_bridge": cm_overlattice_bridge,
            "negative_curve_trichotomy": "boundary curves, (-2)-curves, and (-1)-curves are the only negative curves",
            "complete_minus_two_curves": "the 20 Type-VII roots are all effective irreducible roots by Kondo Proposition 3.7",
            "one_boundary_cusp_lift": "every primitive CM cusp f gives the unique (-1)-curve E=f-b in B+2E=2f, and every (-1)-curve returns f=E+b",
            "complete_minus_one_curves": 51,
            "full_nef_wall_formula": "{b=B/2=-K} union {20 CM (-2)-curves} union {E=f-b for all 51 CM cusps}",
            "why_complete": "the three negative-curve classes are exhausted, and nonnegative-square irreducible curves add no hyperbolic nef facet",
        },
        "polytope": {
            "id": "COBLE-VII-p5-b1-H10-72",
            "dimension": DIMENSION,
            "facet_count": FACET_COUNT,
            "classification": "finite_volume_noncompact",
            "compact": False,
            "finite_volume": True,
            "arithmetic": True,
            "lattice": "I_(1,10) = E10 orthogonal_sum <-1>",
            "lattice_gram": _matrix_rows(lattice),
            "coordinate_basis_facets": [index + 1 for index in pivots],
            "wall_order": "b=-K; 20 Type-VII CM roots; 51 cusp lifts E=f-b",
            "wall_vectors": [[int(value) for value in wall] for wall in walls],
            "negative_wall_gram": _matrix_rows(negative_gram),
            "normalized_gram_sha256": hashlib.sha256(
                json.dumps(_matrix_rows(negative_gram), separators=(",", ":")).encode()
            ).hexdigest(),
            "wall_curve_count": {
                "boundary_B_square_minus_4_primitive_normal_b_square_minus_1": 1,
                "minus_2_curves": 20,
                "minus_1_curves": 51,
            },
            "arithmetic_checks": arithmetic,
            "exact_double_description": cone,
            "literal_vinberg": literal,
            "solomon": {
                "connected_elliptic_count": len(connected_finite),
                "connected_elliptic_search": connected_search,
                "connected_elliptic_rank_histogram": solomon[
                    "connected_elliptic_rank_histogram"
                ],
                "finite_special_subset_count_by_rank": solomon[
                    "finite_special_subset_count_by_rank"
                ],
                "rational_coxeter_euler_characteristic": str(euler),
                "P10_rational_euler_characteristic": str(P10_EULER),
                "index_in_W_P10": EXPECTED_P10_INDEX,
            },
        },
        "method": {
            "wall_source": "exact Shimada 20C Type-VII E10 vectors and their exact 51 CM ideal rays",
            "cone": "complete GMP rational double description with incidence and adjacency",
            "finite_diagrams": "all 15,619,250 pair-compatible rank-9 and 24,466,316 rank-10 cliques classified by the finite crystallographic Coxeter list",
            "affine_diagrams": "all connected affine extensions through order 12 checked by exact determinant with a positive-definite principal deletion",
            "diagram_enumerator_source": str(DIAGRAM_SOURCE.relative_to(ROOT)),
            "diagram_enumerator_source_sha256": _sha256(DIAGRAM_SOURCE),
            "compiler": compiler_used,
            "blind_binomial_subset_scan_used": False,
        },
        "lineage": {
            "surface_and_20_root_CM_H9_chamber": "known: Kondo Type VII / Shimada 20C=20D",
            "full_NS_H10_72_chamber": "new exact positive chamber in this project ledger",
            "new_content": "complete 51-cusp-to-(-1)-curve lift, the 72-wall full-nef Gram, exact face lattice, and literal Vinberg certificate",
            "P10_index_lineage": "777093120 equals Shimada's published index for the Type-VII 20C H9 induced chamber; this inherited equality also occurs for the Type-I/II one-boundary lifts and is structural evidence for the cusp-extension construction, not a novelty claim",
            "worldwide_literal_abstract_priority": "open; no worldwide novelty claim",
            "bounded_coxiter_corpus": bounded_lineage,
        },
        "next_surface_rows": {
            "VI_char5_b1": "same one-boundary theorem predicts 1+20+57=78 full walls in H10; exact full cone/Vinberg replay pending",
            "V_char3_b2": "exact marking of the 16 (-2) and 4 (-1)-root curves plus all single/pair-boundary (-1)-curve lifts remains to be reconstructed",
            "MI_char3_b2": "explicit 40-root geometry is known, but the full cusp/(-1)-wall census and full-NS Vinberg replay remain pending",
        },
        "sources": {
            "kondo": {
                "title": "Coble surfaces with finite automorphism group",
                "arxiv": "2102.00591",
                "used_results": [
                    "Theorem 1.1/Table 2",
                    "Proposition 3.1",
                    "Proposition 3.4",
                    "Proposition 3.7",
                    "Proposition 3.8",
                    "Lemma 3.11",
                    "Example 4.15",
                    "Proposition 4.16",
                ],
            },
            "vinberg": {
                "criterion": "every elliptic rank-(n-1) edge has exactly two elliptic-rank-n/parabolic-rank-(n-1) endpoints; every affine component has a maximal full parabolic zero-neighbourhood",
                "reference": "E. B. Vinberg, Hyperbolic reflection groups, Russian Math. Surveys 40:1 (1985), Proposition 4.2",
            },
            "pinned_input_sha256": {
                str(ENRIQUES_SUMMARY.relative_to(ROOT)): _sha256(ENRIQUES_SUMMARY),
                str(ENRIQUES_RAYS.relative_to(ROOT)): _sha256(ENRIQUES_RAYS),
            },
        },
        "artifacts": {
            "ray_archive": str(RAY_OUTPUT.relative_to(ROOT)),
            "ray_archive_sha256": _sha256_bytes(ray_archive),
            "edge_archive": str(EDGE_OUTPUT.relative_to(ROOT)),
            "edge_archive_sha256": _sha256_bytes(edge_archive),
        },
        "checks": {
            "complete_negative_wall_list": True,
            "all_72_facets_irredundant": True,
            "all_72_reflections_integral": True,
            "all_pairs_are_Coxeter": True,
            "all_143565_elliptic_edges_enumerated": True,
            "every_edge_has_exactly_two_endpoints": True,
            "all_4943_connected_affine_components_enumerated": True,
            "every_affine_component_completes_to_a_rank_9_cusp": True,
            "all_1731_ideal_vertices_are_full_rank_9_parabolic_neighbourhoods": True,
            "finite_volume": True,
            "noncompact": True,
        },
        "all_checks_pass": True,
    }
    return payload, ray_archive, edge_archive


def _write_or_check_text(path: Path, payload: str, check: bool) -> None:
    if check:
        if not path.is_file() or path.read_text(encoding="utf-8") != payload:
            raise SystemExit(f"archived artifact differs: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(payload, encoding="utf-8")


def _write_or_check_bytes(path: Path, payload: bytes, check: bool) -> None:
    if check:
        if not path.is_file() or path.read_bytes() != payload:
            raise SystemExit(f"archived artifact differs: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cdd", default=None)
    parser.add_argument("--compiler", default=None)
    parser.add_argument("--coxiter-corpus", type=Path, default=DEFAULT_COXITER)
    parser.add_argument("--without-coxiter-corpus", action="store_true")
    parser.add_argument("--summary-output", type=Path, default=SUMMARY_OUTPUT)
    parser.add_argument("--ray-output", type=Path, default=RAY_OUTPUT)
    parser.add_argument("--edge-output", type=Path, default=EDGE_OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    corpus = None if arguments.without_coxiter_corpus else arguments.coxiter_corpus
    summary, rays, edges = build_certificate(arguments.cdd, arguments.compiler, corpus)
    _write_or_check_text(
        arguments.summary_output,
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        arguments.check_archive,
    )
    _write_or_check_bytes(arguments.ray_output, rays, arguments.check_archive)
    _write_or_check_bytes(arguments.edge_output, edges, arguments.check_archive)
    print(
        json.dumps(
            {
                "id": summary["polytope"]["id"],
                "vertices": {
                    "finite": summary["polytope"]["exact_double_description"][
                        "finite_vertex_count"
                    ],
                    "ideal": summary["polytope"]["exact_double_description"][
                        "ideal_vertex_count"
                    ],
                },
                "edges": summary["polytope"]["literal_vinberg"][
                    "exact_geometric_edge_count"
                ],
                "finite_volume": summary["checks"]["finite_volume"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
