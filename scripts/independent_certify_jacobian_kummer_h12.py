"""Independent exact cone/Vinberg audit of the Jacobian--Kummer H12 fold.

The production fold artifact supplies an exact 172 by 172 wall Gram matrix.
This script first reconstructs that matrix from Shimada's published D16 data
and replays every facet/redundancy certificate.  It then presents the fixed
13-dimensional chamber cone to the GMP build of cddlib and checks every
primitive extreme ray by integer arithmetic.

The compressed ray archive stores one JSON record per projective vertex:
``[primitive_ray_in_the_13_normal_basis, active_facets]``.  Facets are numbered
from one.  The archive is enough to replay the finite/ideal classification and
the complete edge-endpoint incidence audit without rerunning cddlib.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import gzip
import itertools
import json
import math
from pathlib import Path
import shutil
import subprocess
from typing import Any, Iterable

import numpy as np
import sympy as sp

from scripts.search_apery_fermi_chamber_folds import _orbits
from scripts.search_jacobian_kummer_chamber_folds import (
    _parse_source,
    _wall_permutation,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
INTERIM = (
    ROOT
    / "data"
    / "enumerations"
    / "jacobian_kummer_h12_interim_facet_gram.json"
)
AFFINE_AUDIT = (
    ROOT
    / "data"
    / "enumerations"
    / "jacobian_kummer_h12_affine_neighborhood_audit.json"
)
INEQUALITY_OUTPUT = (
    ROOT / "data" / "enumerations" / "jacobian_kummer_h12_exact_cone.ine"
)
RAY_OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "jacobian_kummer_h12_exact_cone_rays.jsonl.gz"
)
SUMMARY_OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "jacobian_kummer_h12_exact_cone_certificate.json"
)


def _fraction_text(value: int | Fraction | sp.Rational) -> str:
    value = Fraction(value)
    return str(value.numerator) if value.denominator == 1 else str(value)


def load_ray_archive(
    path: Path = RAY_OUTPUT,
) -> tuple[dict[str, Any], tuple[tuple[tuple[int, ...], tuple[int, ...]], ...]]:
    """Load the compact exact-ray archive for a cddlib-free replay.

    Returned active-facet indices remain one-based, as documented by the
    archive header.  The full archive is about 119,000 short records, so this
    convenience loader is appropriate for tests and independent audits.
    """

    with gzip.open(path, "rt", encoding="utf-8") as stream:
        header = json.loads(next(stream))
        records = tuple(
            (tuple(record[0]), tuple(record[1]))
            for record in (json.loads(line) for line in stream)
        )
    if len(records) != header["ray_count"]:
        raise ArithmeticError("the compressed exact-ray archive is truncated")
    return header, records


def _sympy_matrix(payload: list[list[str]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(value) for value in row] for row in payload])


def _primitive_integer_ray(tokens: Iterable[str]) -> tuple[int, ...]:
    fractions = tuple(Fraction(token) for token in tokens)
    denominator = math.lcm(*(value.denominator for value in fractions))
    integers = tuple(
        value.numerator * (denominator // value.denominator)
        for value in fractions
    )
    divisor = math.gcd(*(abs(value) for value in integers))
    if divisor == 0:
        raise ArithmeticError("cddlib emitted the zero ray")
    return tuple(value // divisor for value in integers)


def _matrix_strings(matrix: sp.Matrix) -> list[list[str]]:
    return [
        [_fraction_text(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _rebuild_exact_cone() -> dict[str, Any]:
    interim = json.loads(INTERIM.read_text(encoding="utf-8"))
    source = _parse_source()
    gram: sp.Matrix = source["gram"]
    ample: sp.Matrix = source["ample"]
    representative = _sympy_matrix(
        interim["group_class"]["representative_matrix"]
    )
    identity = sp.eye(gram.rows)
    if not (
        representative**2 == identity
        and representative * gram * representative.T == gram
        and ample * representative == ample
    ):
        raise ArithmeticError("the archived H12 representative failed replay")
    fixed_rank = gram.rows - (representative.T - identity).rank()
    if fixed_rank != 13:
        raise ArithmeticError("the archived representative does not fix rank 13")

    permutation = _wall_permutation(
        representative,
        gram,
        ample,
        source["wall_vectors"],
    )
    subgroup = frozenset((tuple(range(len(permutation))), permutation))
    wall_orbits = _orbits(subgroup, len(permutation))
    walls = tuple(
        sp.Matrix(1, gram.rows, vector) for vector in source["wall_vectors"]
    )
    orbit_normals = tuple(
        sum((walls[index] for index in orbit), sp.zeros(1, gram.rows))
        for orbit in wall_orbits
    )
    if len(orbit_normals) != 180:
        raise ArithmeticError("the H12 representative does not have 180 wall orbits")
    orbit_gram = sp.Matrix(
        [
            [(left * gram * right.T)[0] for right in orbit_normals]
            for left in orbit_normals
        ]
    )
    if exact_inertia(orbit_gram) != (1, 12, 167):
        raise ArithmeticError("the orbit-normal span is not the claimed H12 space")

    facet_indices = tuple(index - 1 for index in interim["facet_orbit_indices"])
    facet_set = set(facet_indices)
    redundancy = {
        record["orbit_index"] - 1: record
        for record in interim["redundancy_certificates"]
    }
    if not (
        len(facet_indices) == 172
        and len(redundancy) == 8
        and facet_set.isdisjoint(redundancy)
        and facet_set | set(redundancy) == set(range(180))
    ):
        raise ArithmeticError("the facet/redundancy partition is not exhaustive")

    for target, record in redundancy.items():
        recovered = sum(
            (
                sp.Rational(term["coefficient"])
                * orbit_normals[term["orbit_index"] - 1]
                for term in record["terms"]
            ),
            sp.zeros(1, gram.rows),
        )
        if recovered != orbit_normals[target]:
            raise ArithmeticError("an exact conic redundancy failed replay")

    for target in facet_indices:
        normal = orbit_normals[target]
        norm = sp.Rational(orbit_gram[target, target])
        ample_pairing = sp.Rational((ample * gram * normal.T)[0])
        witness = ample - ample_pairing / norm * normal
        pairings = tuple(
            sp.Rational((witness * gram * other.T)[0])
            for other in orbit_normals
        )
        if not (
            norm < 0
            and (witness * gram * witness.T)[0] > 0
            and pairings[target] == 0
            and all(
                pairing > 0
                for index, pairing in enumerate(pairings)
                if index != target
            )
        ):
            raise ArithmeticError("a retained inequality is not a strict facet")

    raw_facet_gram = orbit_gram.extract(facet_indices, facet_indices)
    if _matrix_strings(raw_facet_gram) != interim["facet_gram"]:
        raise ArithmeticError("the independently rebuilt facet Gram changed")
    positive_facet_gram = -raw_facet_gram
    if _matrix_strings(positive_facet_gram) != interim["positive_wall_gram"]:
        raise ArithmeticError("the positive-wall Gram changed")
    if exact_inertia(positive_facet_gram) != (12, 1, 159):
        raise ArithmeticError("the 172-wall Gram does not have H12 inertia")

    _, pivot_columns = positive_facet_gram.rref()
    basis = tuple(pivot_columns)
    if len(basis) != 13:
        raise ArithmeticError("the facet normals do not span dimension 13")
    basis_normals = sp.Matrix.vstack(*(orbit_normals[facet_indices[i]] for i in basis))
    ample_coordinates_solution = sp.linsolve(
        (basis_normals.T, ample.T)
    )
    ample_coordinates = tuple(next(iter(ample_coordinates_solution)))
    if any(not coordinate.is_Rational for coordinate in ample_coordinates):
        raise ArithmeticError("the ample class has no rational normal coordinates")

    # If x = sum_j c_j n_j in the selected normal basis, then the i-th
    # chamber inequality is <x,n_i> = sum_j c_j <n_j,n_i> >= 0.
    inequality_matrix = raw_facet_gram[:, list(basis)]
    ample_slacks = inequality_matrix * sp.Matrix(ample_coordinates)
    if any(slack <= 0 for slack in ample_slacks):
        raise ArithmeticError("the ample class is not a strict cone interior point")
    if inequality_matrix.rank() != 13:
        raise ArithmeticError("the chamber cone has nonzero lineality")

    lines = [
        "H-representation",
        "begin",
        f"{inequality_matrix.rows} {inequality_matrix.cols + 1} rational",
    ]
    lines.extend(
        "0 "
        + " ".join(
            _fraction_text(inequality_matrix[row, column])
            for column in range(inequality_matrix.cols)
        )
        for row in range(inequality_matrix.rows)
    )
    lines.append("end")
    inequality_text = "\n".join(lines) + "\n"

    return {
        "interim": interim,
        "positive_facet_gram": positive_facet_gram,
        "raw_basis_gram": raw_facet_gram.extract(basis, basis),
        "inequality_matrix": inequality_matrix,
        "inequality_text": inequality_text,
        "basis": basis,
        "ample_coordinates": ample_coordinates,
        "minimum_ample_slack": min(ample_slacks),
    }


def _run_cdd(inequality_text: str, executable: str) -> str:
    completed = subprocess.run(
        [executable, "--rep"],
        input=inequality_text,
        text=True,
        capture_output=True,
        timeout=1800,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "cddexec_gmp failed:\n"
            + completed.stdout[-2000:]
            + "\n"
            + completed.stderr[-2000:]
        )
    return completed.stdout


def _parse_cdd_rays(text: str, dimension: int) -> tuple[tuple[int, ...], ...]:
    lines = text.splitlines()
    try:
        begin = lines.index("begin")
    except ValueError as error:
        raise ArithmeticError("cddlib output has no begin marker") from error
    header = lines[begin + 1].split()
    if len(header) != 3 or header[2] != "rational":
        raise ArithmeticError("cddlib did not use exact rational arithmetic")
    count, width = map(int, header[:2])
    if width != dimension + 1:
        raise ArithmeticError("cddlib output has the wrong ray dimension")
    ray_lines = lines[begin + 2 : begin + 2 + count]
    if len(ray_lines) != count or lines[begin + 2 + count].strip() != "end":
        raise ArithmeticError("cddlib output is truncated")
    rays: list[tuple[int, ...]] = []
    for line in ray_lines:
        tokens = line.split()
        if len(tokens) != width or Fraction(tokens[0]) != 0:
            raise ArithmeticError("the homogeneous cone output contains a vertex row")
        rays.append(_primitive_integer_ray(tokens[1:]))
    if len(set(rays)) != len(rays):
        raise ArithmeticError("cddlib emitted duplicate primitive rays")
    return tuple(rays)


def _connected_components(
    facets: tuple[int, ...], positive_gram: np.ndarray
) -> tuple[tuple[int, ...], ...]:
    unseen = set(facets)
    components: list[tuple[int, ...]] = []
    while unseen:
        start = min(unseen)
        unseen.remove(start)
        component = {start}
        frontier = [start]
        while frontier:
            node = frontier.pop()
            new = {other for other in unseen if positive_gram[node, other] != 0}
            unseen.difference_update(new)
            component.update(new)
            frontier.extend(new)
        components.append(tuple(sorted(component)))
    return tuple(sorted(components))


def _analyze_rays(
    cone: dict[str, Any], rays: tuple[tuple[int, ...], ...]
) -> tuple[dict[str, Any], bytes]:
    positive_gram = np.array(
        cone["positive_facet_gram"].tolist(), dtype=np.int64
    )
    inequality_matrix = np.array(
        cone["inequality_matrix"].tolist(), dtype=np.int64
    )
    raw_basis_gram = np.array(cone["raw_basis_gram"].tolist(), dtype=np.int64)
    ray_matrix = np.array(rays, dtype=np.int64)
    if ray_matrix.shape != (119_488, 13):
        raise ArithmeticError("the exact H12 extreme-ray count changed")

    coordinate_bound = int(np.max(np.abs(ray_matrix)))
    arithmetic_bound = (
        ray_matrix.shape[1] ** 2
        * coordinate_bound**2
        * max(1, int(np.max(np.abs(raw_basis_gram))))
    )
    if arithmetic_bound >= np.iinfo(np.int64).max:
        raise OverflowError("the exact integer replay exceeds int64")

    norms = np.einsum(
        "bi,ij,bj->b", ray_matrix, raw_basis_gram, ray_matrix
    )
    if np.any(norms < 0):
        raise ArithmeticError("a hyperideal/spacelike extreme ray was found")

    incidences: list[tuple[int, ...]] = []
    active_histogram: Counter[int] = Counter()
    finite_active_histogram: Counter[int] = Counter()
    ideal_active_histogram: Counter[int] = Counter()
    for start in range(0, len(rays), 2048):
        slacks = ray_matrix[start : start + 2048] @ inequality_matrix.T
        if np.any(slacks < 0):
            raise ArithmeticError("an exact cddlib ray violates a cone inequality")
        active_counts = np.sum(slacks == 0, axis=1)
        for local, row in enumerate(slacks):
            incidence = tuple(int(value) for value in np.flatnonzero(row == 0))
            incidences.append(incidence)
            active_histogram[len(incidence)] += 1
            target = (
                ideal_active_histogram
                if norms[start + local] == 0
                else finite_active_histogram
            )
            target[len(incidence)] += 1
    if len(set(incidences)) != len(incidences):
        raise ArithmeticError("two projective vertices have the same active facets")
    if finite_active_histogram != Counter({12: 113_274}):
        raise ArithmeticError("the finite vertices are not all simple")

    ideal_incidence_set = {
        incidence
        for incidence, norm in zip(incidences, norms, strict=True)
        if norm == 0
    }
    affine = json.loads(AFFINE_AUDIT.read_text(encoding="utf-8"))
    affine_neighborhoods: set[tuple[int, ...]] = set()
    for component_one_based in affine["connected_affine_components"]:
        component = tuple(index - 1 for index in component_one_based)
        zero_nodes = tuple(
            vertex
            for vertex in range(positive_gram.shape[0])
            if vertex not in component
            and all(positive_gram[vertex, node] == 0 for node in component)
        )
        affine_neighborhoods.add(tuple(sorted((*component, *zero_nodes))))
    if affine_neighborhoods != ideal_incidence_set:
        raise ArithmeticError(
            "the exact cone cusps and exhaustive affine neighborhoods disagree"
        )
    if affine["invalid_full_zero_neighborhood_count"] != 0:
        raise ArithmeticError("the affine audit contains an invalid neighborhood")

    # Encode endpoint multiplicities in one byte pair: low byte finite,
    # high bytes ideal.  An edge is the 172-bit mask of its eleven facets.
    edges: dict[int, int] = {}
    finite_endpoint_incidences = 0
    ideal_endpoint_incidences = 0
    cusp_component_count_histogram: Counter[int] = Counter()
    cusp_component_order_pattern_histogram: Counter[tuple[int, ...]] = Counter()
    for incidence, norm in zip(incidences, norms, strict=True):
        mask = sum(1 << facet for facet in incidence)
        if norm > 0:
            finite_endpoint_incidences += len(incidence)
            for omitted in incidence:
                edge = mask ^ (1 << omitted)
                edges[edge] = edges.get(edge, 0) + 1
            continue

        components = _connected_components(incidence, positive_gram)
        if len(incidence) != 11 + len(components):
            raise ArithmeticError("a cusp is not parabolic of rank eleven")
        cusp_component_count_histogram[len(components)] += 1
        cusp_component_order_pattern_histogram[
            tuple(sorted(len(component) for component in components))
        ] += 1
        local_edges = math.prod(len(component) for component in components)
        ideal_endpoint_incidences += local_edges
        for omitted_nodes in itertools.product(*components):
            edge = mask
            for omitted in omitted_nodes:
                edge ^= 1 << omitted
            if edge.bit_count() != 11:
                raise ArithmeticError("a cusp deletion did not produce rank eleven")
            edges[edge] = edges.get(edge, 0) + 256

    endpoint_histogram: Counter[int] = Counter()
    endpoint_kind_histogram: Counter[tuple[int, int]] = Counter()
    for edge, encoded in edges.items():
        if edge.bit_count() != 11:
            raise ArithmeticError("an edge key does not have eleven facets")
        finite_count = encoded & 255
        ideal_count = encoded >> 8
        endpoint_histogram[finite_count + ideal_count] += 1
        endpoint_kind_histogram[(finite_count, ideal_count)] += 1
    if endpoint_histogram != Counter({2: 910_316}):
        raise ArithmeticError("a rank-eleven edge does not have two endpoints")

    records = sorted(
        zip(rays, incidences, strict=True), key=lambda item: (item[1], item[0])
    )
    archive_header = {
        "schema": "jacobian_kummer_h12_exact_cone_rays_v1",
        "coordinate_basis_facet_indices": [index + 1 for index in cone["basis"]],
        "facet_indexing": "one_based",
        "record_format": "[primitive_integer_ray, active_facet_indices]",
        "ray_count": len(records),
    }
    lines = [json.dumps(archive_header, sort_keys=True, separators=(",", ":"))]
    lines.extend(
        json.dumps(
            [list(ray), [index + 1 for index in incidence]],
            separators=(",", ":"),
        )
        for ray, incidence in records
    )
    ray_archive = gzip.compress(
        ("\n".join(lines) + "\n").encode("utf-8"),
        compresslevel=9,
        mtime=0,
    )

    norm_histogram = Counter(int(value) for value in norms)
    summary = {
        "schema": "jacobian_kummer_h12_exact_cone_certificate_v1",
        "source": {
            "surface": "general Jacobian Kummer K3 surface Y16",
            "published_chamber": "D16",
            "published_data": "D16WallRecs in 15nodalcompdata.txt",
            "paper_doi": "10.1007/s12215-019-00464-7",
            "fold_group": "trace-9 involution in O(S16,D16)",
        },
        "cone": {
            "hyperbolic_dimension": 12,
            "linear_dimension": 13,
            "facet_count": positive_gram.shape[0],
            "coordinate_basis_facet_indices": [
                index + 1 for index in cone["basis"]
            ],
            "inequality_rank": int(sp.Matrix(inequality_matrix.tolist()).rank()),
            "lineality_dimension": 0,
            "strict_interior_witness": "published ample16",
            "minimum_strict_ample_pairing": _fraction_text(
                cone["minimum_ample_slack"]
            ),
            "full_dimensional": True,
            "pointed": True,
        },
        "exact_double_description": {
            "engine": "cddlib GMP rational double-description",
            "command": "cddexec_gmp --rep",
            "extreme_ray_count": len(rays),
            "primitive_ray_coordinate_absolute_maximum": coordinate_bound,
            "lorentz_square_histogram": {
                str(key): value for key, value in sorted(norm_histogram.items())
            },
            "finite_vertex_count": int(np.sum(norms > 0)),
            "ideal_vertex_count": int(np.sum(norms == 0)),
            "hyperideal_vertex_count": int(np.sum(norms < 0)),
            "active_facet_count_histogram": {
                str(key): value for key, value in sorted(active_histogram.items())
            },
            "finite_vertex_active_facet_count_histogram": {
                str(key): value
                for key, value in sorted(finite_active_histogram.items())
            },
            "ideal_vertex_active_facet_count_histogram": {
                str(key): value
                for key, value in sorted(ideal_active_histogram.items())
            },
        },
        "affine_full_neighborhood_cross_check": {
            "connected_affine_component_count": len(
                affine["connected_affine_components"]
            ),
            "distinct_full_neighborhood_count": len(affine_neighborhoods),
            "invalid_full_neighborhood_count": affine[
                "invalid_full_zero_neighborhood_count"
            ],
            "equals_exact_lightlike_ray_incidence_sets": True,
            "components_per_cusp_histogram": {
                str(key): value
                for key, value in sorted(cusp_component_count_histogram.items())
            },
            "component_order_pattern_histogram": {
                "+".join(map(str, key)): value
                for key, value in sorted(
                    cusp_component_order_pattern_histogram.items()
                )
            },
        },
        "vinberg_edge_endpoint_audit": {
            "elliptic_rank_11_edge_count": len(edges),
            "endpoint_count_histogram": {
                str(key): value for key, value in sorted(endpoint_histogram.items())
            },
            "finite_finite_edge_count": endpoint_kind_histogram[(2, 0)],
            "finite_ideal_edge_count": endpoint_kind_histogram[(1, 1)],
            "ideal_ideal_edge_count": endpoint_kind_histogram[(0, 2)],
            "finite_endpoint_incidence_count": finite_endpoint_incidences,
            "ideal_endpoint_incidence_count": ideal_endpoint_incidences,
            "bad_edge_count": sum(
                value for key, value in endpoint_histogram.items() if key != 2
            ),
        },
        "decision": {
            "finite_volume": True,
            "compact": False,
            "classification": "finite_volume_noncompact",
            "projective_cone_reason": (
                "the full-dimensional pointed chamber cone has only timelike "
                "or lightlike extreme rays and no hyperideal extreme ray"
            ),
            "vinberg_reason": (
                "every exhaustive affine full neighborhood is rank-11 "
                "parabolic and every elliptic rank-11 edge has two endpoints"
            ),
        },
        "checks": {
            "published_lattice_and_wall_action_replayed": True,
            "all_redundancies_replayed_exactly": True,
            "all_172_facets_have_strict_timelike_witnesses": True,
            "gmp_output_is_rational": True,
            "all_rays_satisfy_all_172_inequalities_exactly": True,
            "all_119488_ray_incidence_sets_are_distinct": True,
            "no_hyperideal_extreme_ray": True,
            "all_6214_cusps_equal_exhaustive_affine_neighborhoods": True,
            "all_910316_edges_have_exactly_two_endpoints": True,
        },
    }
    return summary, ray_archive


def _write_or_check_text(path: Path, text: str, check: bool) -> None:
    if check:
        if not path.is_file() or path.read_text(encoding="utf-8") != text:
            raise SystemExit(f"archived artifact differs: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _write_or_check_bytes(path: Path, payload: bytes, check: bool) -> None:
    if check:
        if not path.is_file() or path.read_bytes() != payload:
            raise SystemExit(f"archived artifact differs: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cdd-executable",
        default=shutil.which("cddexec_gmp") or "cddexec_gmp",
    )
    parser.add_argument(
        "--ext-input",
        type=Path,
        help="reuse an exact cddexec_gmp V-representation instead of rerunning cdd",
    )
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--inequality-output", type=Path, default=INEQUALITY_OUTPUT)
    parser.add_argument("--ray-output", type=Path, default=RAY_OUTPUT)
    parser.add_argument("--summary-output", type=Path, default=SUMMARY_OUTPUT)
    arguments = parser.parse_args()

    cone = _rebuild_exact_cone()
    if arguments.ext_input is not None:
        cdd_output = arguments.ext_input.read_text(encoding="utf-8")
    else:
        cdd_output = _run_cdd(cone["inequality_text"], arguments.cdd_executable)
    rays = _parse_cdd_rays(cdd_output, 13)
    summary, ray_archive = _analyze_rays(cone, rays)

    _write_or_check_text(
        arguments.inequality_output,
        cone["inequality_text"],
        arguments.check_archive,
    )
    _write_or_check_bytes(arguments.ray_output, ray_archive, arguments.check_archive)
    summary_text = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    _write_or_check_text(
        arguments.summary_output,
        summary_text,
        arguments.check_archive,
    )
    print(json.dumps(summary["decision"] | summary["vinberg_edge_endpoint_audit"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
