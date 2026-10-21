"""Independent GMP exact-cone certificate for Shimada disc16 D5 H14.

The cone is reconstructed from the official disc16 atlas through the
intrinsic D5 rank-15 stabilizer selector in
``independent_screen_shimada_disc16_d5_h14``.  Every facet/redundancy witness
is replayed before the 60 exact inequalities are sent to cddlib's GMP
rational double-description engine.  Finite volume is decided only from the
complete extreme-ray census, never from pair compatibility alone.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import gzip
import json
import math
from pathlib import Path
import shutil
import subprocess
from typing import Any, Iterable

import numpy as np
import sympy as sp

from scripts.independent_screen_shimada_disc16_d5_h14 import (
    OUTPUT as SCREEN_OUTPUT,
    ROOT,
    SOURCE,
    build_screen,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


INEQUALITY_OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_d5_h14_independent_exact_cone.ine"
)
RAY_OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_d5_h14_independent_exact_cone_rays.jsonl.gz"
)
SUMMARY_OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_d5_h14_independent_exact_cone_certificate.json"
)
DIMENSION = 14
LINEAR_DIMENSION = 15


def _fraction_text(value: int | Fraction | sp.Rational) -> str:
    value = Fraction(value)
    return str(value.numerator) if value.denominator == 1 else str(value)


def _sympy_matrix(payload: list[list[str]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(value) for value in row] for row in payload])


def _matrix_strings(matrix: sp.Matrix) -> list[list[str]]:
    return [
        [_fraction_text(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _primitive_integer_ray(tokens: Iterable[str]) -> tuple[int, ...]:
    values = tuple(Fraction(token) for token in tokens)
    denominator = math.lcm(*(value.denominator for value in values))
    integers = tuple(
        value.numerator * (denominator // value.denominator) for value in values
    )
    divisor = math.gcd(*(abs(value) for value in integers))
    if divisor == 0:
        raise ArithmeticError("cddlib emitted the zero ray")
    return tuple(value // divisor for value in integers)


def load_ray_archive(
    path: Path = RAY_OUTPUT,
) -> tuple[dict[str, Any], tuple[tuple[tuple[int, ...], tuple[int, ...]], ...]]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        header = json.loads(next(stream))
        records = tuple(
            (tuple(record[0]), tuple(record[1]))
            for record in (json.loads(line) for line in stream)
        )
    if len(records) != header["ray_count"]:
        raise ArithmeticError("the compressed ray archive is truncated")
    return header, records


def _rebuild_exact_cone(path: Path = SOURCE) -> dict[str, Any]:
    screen = build_screen(path)
    if not (
        screen["status"]["all_facet_pairs_are_coxeter"]
        and screen["representative"]["fixed_rank"] == LINEAR_DIMENSION
        and screen["facet_count"] == 60
        and screen["redundant_orbit_count"] == 40
    ):
        raise ArithmeticError("the independent D5 screen is not the expected H14 candidate")

    dual_gram = _sympy_matrix(screen["ambient_lattice"]["dual_gram"])
    walls = tuple(
        sp.Matrix(1, 20, wall)
        for wall in screen["chamber"]["published_primitive_dual_walls"]
    )
    orbits = tuple(
        tuple(index - 1 for index in orbit) for orbit in screen["wall_orbits"]
    )
    orbit_normals = tuple(
        sum((walls[index] for index in orbit), sp.zeros(1, 20))
        for orbit in orbits
    )
    stacked = sp.Matrix.vstack(*orbit_normals)
    orbit_gram = stacked * dual_gram * stacked.T
    if exact_inertia(orbit_gram) != (1, 14, 85):
        raise ArithmeticError("the independently rebuilt orbit Gram is not H14")

    facet_indices = tuple(index - 1 for index in screen["facet_orbit_indices"])
    redundancies = {
        record["orbit_index"] - 1: record
        for record in screen["redundancy_certificates"]
    }
    if not (
        len(facet_indices) == 60
        and len(redundancies) == 40
        and set(facet_indices).isdisjoint(redundancies)
        and set(facet_indices) | set(redundancies) == set(range(100))
    ):
        raise ArithmeticError("the facet/redundancy partition is not exhaustive")
    for target, record in redundancies.items():
        recovered = sum(
            (
                sp.Rational(term["coefficient"])
                * orbit_normals[term["orbit_index"] - 1]
                for term in record["terms"]
            ),
            sp.zeros(1, 20),
        )
        if recovered != orbit_normals[target]:
            raise ArithmeticError("a frozen conic redundancy failed exact replay")

    ample = sp.Matrix(1, 20, screen["chamber"]["ample_dual_coordinates"])
    for target in facet_indices:
        normal = orbit_normals[target]
        norm = sp.Rational(orbit_gram[target, target])
        pairing = sp.Rational((ample * dual_gram * normal.T)[0])
        witness = ample - pairing / norm * normal
        all_pairings = tuple(
            sp.Rational((witness * dual_gram * other.T)[0])
            for other in orbit_normals
        )
        if not (
            norm < 0
            and (witness * dual_gram * witness.T)[0] > 0
            and all_pairings[target] == 0
            and all(
                value > 0
                for index, value in enumerate(all_pairings)
                if index != target
            )
        ):
            raise ArithmeticError("a frozen facet witness failed exact replay")

    raw_facet_gram = orbit_gram.extract(facet_indices, facet_indices)
    if _matrix_strings(raw_facet_gram) != screen["facet_gram"]:
        raise ArithmeticError("the independently rebuilt 60-wall Gram changed")
    positive_facet_gram = -raw_facet_gram
    if exact_inertia(positive_facet_gram) != (14, 1, 45):
        raise ArithmeticError("the positive 60-wall Gram has wrong H14 inertia")

    _, pivot_columns = positive_facet_gram.rref()
    basis = tuple(pivot_columns)
    if len(basis) != LINEAR_DIMENSION:
        raise ArithmeticError("the 60 facet normals do not span rank 15")
    basis_normals = sp.Matrix.vstack(
        *(orbit_normals[facet_indices[index]] for index in basis)
    )
    solution = sp.linsolve((basis_normals.T, ample.T))
    ample_coordinates = tuple(next(iter(solution)))
    if any(not coordinate.is_Rational for coordinate in ample_coordinates):
        raise ArithmeticError("the ample class has no rational facet-basis coordinates")

    inequality_matrix = raw_facet_gram[:, list(basis)]
    ample_slacks = inequality_matrix * sp.Matrix(ample_coordinates)
    if inequality_matrix.rank() != LINEAR_DIMENSION or any(
        slack <= 0 for slack in ample_slacks
    ):
        raise ArithmeticError("the exact cone is not pointed with strict ample interior")
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
    return {
        "screen": screen,
        "positive_facet_gram": positive_facet_gram,
        "raw_basis_gram": raw_facet_gram.extract(basis, basis),
        "inequality_matrix": inequality_matrix,
        "inequality_text": "\n".join(lines) + "\n",
        "basis": basis,
        "minimum_ample_slack": min(ample_slacks),
    }


def _run_cdd(inequality_text: str, executable: str) -> str:
    completed = subprocess.run(
        [executable, "--rep"],
        input=inequality_text,
        text=True,
        capture_output=True,
        timeout=3600,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "cddexec_gmp failed:\n"
            + completed.stdout[-4000:]
            + "\n"
            + completed.stderr[-4000:]
        )
    return completed.stdout


def _parse_cdd_rays(text: str) -> tuple[tuple[int, ...], ...]:
    lines = text.splitlines()
    try:
        begin = lines.index("begin")
    except ValueError as error:
        raise ArithmeticError("cddlib output has no begin marker") from error
    header = lines[begin + 1].split()
    if len(header) != 3 or header[2] != "rational":
        raise ArithmeticError("cddlib did not emit rational generators")
    count, width = map(int, header[:2])
    if width != LINEAR_DIMENSION + 1:
        raise ArithmeticError("cddlib output has wrong generator dimension")
    ray_lines = lines[begin + 2 : begin + 2 + count]
    if len(ray_lines) != count or lines[begin + 2 + count].strip() != "end":
        raise ArithmeticError("cddlib output is truncated")
    rays = []
    for line in ray_lines:
        tokens = line.split()
        if len(tokens) != width or Fraction(tokens[0]) != 0:
            raise ArithmeticError("the homogeneous cone output contains a vertex row")
        rays.append(_primitive_integer_ray(tokens[1:]))
    if len(set(rays)) != len(rays):
        raise ArithmeticError("cddlib emitted duplicate primitive rays")
    return tuple(rays)


def _scaled_integer_matrix(matrix: sp.Matrix) -> tuple[np.ndarray, int]:
    scale = math.lcm(*(int(sp.Rational(value).q) for value in matrix))
    integral = np.array(
        [
            [int(scale * sp.Rational(value)) for value in row]
            for row in matrix.tolist()
        ],
        dtype=np.int64,
    )
    return integral, scale


def _analyze_rays(
    cone: dict[str, Any], rays: tuple[tuple[int, ...], ...]
) -> tuple[dict[str, Any], bytes]:
    inequality, inequality_scale = _scaled_integer_matrix(cone["inequality_matrix"])
    basis_gram, basis_gram_scale = _scaled_integer_matrix(cone["raw_basis_gram"])
    ray_matrix = np.array(rays, dtype=np.int64)
    if ray_matrix.ndim != 2 or ray_matrix.shape[1] != LINEAR_DIMENSION:
        raise ArithmeticError("the exact ray matrix has wrong shape")
    coordinate_bound = int(np.max(np.abs(ray_matrix)))
    norm_bound = (
        LINEAR_DIMENSION**2
        * coordinate_bound**2
        * max(1, int(np.max(np.abs(basis_gram))))
    )
    slack_bound = (
        LINEAR_DIMENSION
        * coordinate_bound
        * max(1, int(np.max(np.abs(inequality))))
    )
    if max(norm_bound, slack_bound) >= np.iinfo(np.int64).max:
        raise OverflowError("the exact integer ray replay exceeds int64")

    norms = np.einsum("bi,ij,bj->b", ray_matrix, basis_gram, ray_matrix)
    incidences: list[tuple[int, ...]] = []
    active_histogram: Counter[int] = Counter()
    finite_active_histogram: Counter[int] = Counter()
    ideal_active_histogram: Counter[int] = Counter()
    hyperideal_active_histogram: Counter[int] = Counter()
    for start in range(0, len(rays), 2048):
        slacks = ray_matrix[start : start + 2048] @ inequality.T
        if np.any(slacks < 0):
            raise ArithmeticError("an exact cddlib ray violates a cone inequality")
        for local, row in enumerate(slacks):
            incidence = tuple(int(value) for value in np.flatnonzero(row == 0))
            if len(incidence) < DIMENSION - 1:
                raise ArithmeticError("an extreme ray has too few active facets")
            incidences.append(incidence)
            active_histogram[len(incidence)] += 1
            norm = norms[start + local]
            target = (
                ideal_active_histogram
                if norm == 0
                else finite_active_histogram
                if norm > 0
                else hyperideal_active_histogram
            )
            target[len(incidence)] += 1
    if len(set(incidences)) != len(incidences):
        raise ArithmeticError("two exact rays have the same complete incidence set")

    finite_count = int(np.sum(norms > 0))
    ideal_count = int(np.sum(norms == 0))
    hyperideal_count = int(np.sum(norms < 0))
    finite_volume = bool(rays) and hyperideal_count == 0
    compact = finite_volume and ideal_count == 0

    records = sorted(zip(rays, incidences, strict=True), key=lambda item: (item[1], item[0]))
    header = {
        "schema": "shimada_disc16_d5_h14_independent_exact_cone_rays_v1",
        "coordinate_basis_facet_indices": [index + 1 for index in cone["basis"]],
        "facet_indexing": "one_based",
        "record_format": "[primitive_integer_ray, active_facet_indices]",
        "ray_count": len(records),
    }
    lines = [json.dumps(header, sort_keys=True, separators=(",", ":"))]
    lines.extend(
        json.dumps(
            [list(ray), [index + 1 for index in incidence]],
            separators=(",", ":"),
        )
        for ray, incidence in records
    )
    archive = gzip.compress(
        ("\n".join(lines) + "\n").encode("utf-8"),
        compresslevel=9,
        mtime=0,
    )

    norm_histogram = Counter(int(value) for value in norms)
    summary = {
        "schema": "shimada_disc16_d5_h14_independent_exact_cone_certificate_v1",
        "source": {
            "archive_url": cone["screen"]["source"]["archive_url"],
            "member": cone["screen"]["source"]["member"],
            "sha256": cone["screen"]["source"]["sha256"],
            "chamber": "D5",
            "representative_selector": "unique D5 stabilizer element of fixed rank 15",
        },
        "cone": {
            "hyperbolic_dimension": DIMENSION,
            "linear_dimension": LINEAR_DIMENSION,
            "facet_count": 60,
            "coordinate_basis_facet_indices": [index + 1 for index in cone["basis"]],
            "inequality_rank": int(cone["inequality_matrix"].rank()),
            "lineality_dimension": 0,
            "strict_interior_witness": "published chamberdata[5][weyl] projected to S dual",
            "minimum_strict_ample_pairing": _fraction_text(cone["minimum_ample_slack"]),
            "full_dimensional": True,
            "pointed": True,
        },
        "exact_double_description": {
            "engine": "cddlib GMP rational double-description",
            "command": "cddexec_gmp --rep",
            "extreme_ray_count": len(rays),
            "primitive_ray_coordinate_absolute_maximum": coordinate_bound,
            "inequality_integral_scale": inequality_scale,
            "basis_gram_integral_scale": basis_gram_scale,
            "lorentz_square_histogram_in_scaled_basis_gram": {
                str(key): value for key, value in sorted(norm_histogram.items())
            },
            "finite_vertex_count": finite_count,
            "ideal_vertex_count": ideal_count,
            "hyperideal_vertex_count": hyperideal_count,
            "active_facet_count_histogram": {
                str(key): value for key, value in sorted(active_histogram.items())
            },
            "finite_vertex_active_facet_count_histogram": {
                str(key): value for key, value in sorted(finite_active_histogram.items())
            },
            "ideal_vertex_active_facet_count_histogram": {
                str(key): value for key, value in sorted(ideal_active_histogram.items())
            },
            "hyperideal_vertex_active_facet_count_histogram": {
                str(key): value for key, value in sorted(hyperideal_active_histogram.items())
            },
        },
        "decision": {
            "finite_volume": finite_volume,
            "compact": compact if finite_volume else None,
            "classification": (
                "compact"
                if compact
                else "finite_volume_noncompact"
                if finite_volume
                else "infinite_volume"
            ),
            "criterion": (
                "the full-dimensional pointed projective chamber has no "
                "spacelike/hyperideal extreme ray"
            ),
        },
        "checks": {
            "primary_disc16_D5_reconstructed_independently": True,
            "all_40_conic_redundancies_replayed_exactly": True,
            "all_60_strict_facet_witnesses_replayed_exactly": True,
            "gmp_output_is_rational": True,
            "all_rays_satisfy_all_60_inequalities_exactly": True,
            "all_ray_incidence_sets_are_distinct": True,
            "no_hyperideal_extreme_ray": hyperideal_count == 0,
        },
    }
    return summary, archive


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
    parser.add_argument("disc16", type=Path, nargs="?", default=SOURCE)
    parser.add_argument("--cdd-executable", default=shutil.which("cddexec_gmp") or "cddexec_gmp")
    parser.add_argument("--ext-input", type=Path)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--inequality-output", type=Path, default=INEQUALITY_OUTPUT)
    parser.add_argument("--ray-output", type=Path, default=RAY_OUTPUT)
    parser.add_argument("--summary-output", type=Path, default=SUMMARY_OUTPUT)
    arguments = parser.parse_args()

    cone = _rebuild_exact_cone(arguments.disc16)
    cdd_output = (
        arguments.ext_input.read_text(encoding="utf-8")
        if arguments.ext_input is not None
        else _run_cdd(cone["inequality_text"], arguments.cdd_executable)
    )
    rays = _parse_cdd_rays(cdd_output)
    summary, archive = _analyze_rays(cone, rays)
    _write_or_check_text(arguments.inequality_output, cone["inequality_text"], arguments.check_archive)
    _write_or_check_bytes(arguments.ray_output, archive, arguments.check_archive)
    _write_or_check_text(
        arguments.summary_output,
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        arguments.check_archive,
    )
    print(json.dumps(summary["decision"] | summary["exact_double_description"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
