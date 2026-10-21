"""Exact cddlib cone audit for the two Fermat-quartic V4 Coxeter folds."""

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
import tempfile
from typing import Any, Iterable
import sys

import numpy as np
import sympy as sp

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_ROOT / "src"))

from scripts.screen_ss_fermat_quartic_char3_v4_folds import (
    OUTPUT as SCREEN_OUTPUT,
    ROOT,
)
from scripts.search_ss_fermat_quartic_char3_folds import SOURCE, parse_source
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


def _fraction(value: str | int | Fraction | sp.Rational) -> Fraction:
    return Fraction(value)


def _fraction_text(value: str | int | Fraction | sp.Rational) -> str:
    value = _fraction(value)
    return str(value.numerator) if value.denominator == 1 else str(value)


def _primitive_integer_ray(tokens: Iterable[str]) -> tuple[int, ...]:
    values = tuple(Fraction(token) for token in tokens)
    denominator = math.lcm(*(value.denominator for value in values))
    integers = tuple(
        value.numerator * (denominator // value.denominator) for value in values
    )
    divisor = math.gcd(*(abs(value) for value in integers))
    if divisor == 0:
        raise ArithmeticError("cddlib emitted the zero ray")
    integers = tuple(value // divisor for value in integers)
    for value in integers:
        if value:
            return integers if value > 0 else tuple(-entry for entry in integers)
    raise AssertionError("unreachable zero-ray branch")


def _matrix_integral_scale(matrix: sp.Matrix) -> tuple[list[list[int]], int]:
    scale = math.lcm(*(int(sp.Rational(value).q) for value in matrix))
    return (
        [
            [int(scale * sp.Rational(value)) for value in matrix.row(row)]
            for row in range(matrix.rows)
        ],
        scale,
    )


def rebuild_cone(dimension: int, screen_path: Path = SCREEN_OUTPUT) -> dict[str, Any]:
    document = json.loads(screen_path.read_text(encoding="utf-8"))
    candidates = [
        cast
        for cast in document["casts"]
        if int(cast["subgroup"]["hyperbolic_dimension"]) == dimension
    ]
    if len(candidates) != 1:
        raise ArithmeticError(f"screen does not contain one H{dimension} cast")
    cast = candidates[0]
    reduction = cast["global_facet_reduction"]
    gate = cast["strict_facet_gate"]
    if not (
        reduction["exact_complete"]
        and reduction["unresolved_after_strict_facet_cone_count"] == 0
        and gate["bad_pair_count"] == 0
    ):
        raise ArithmeticError(f"H{dimension} is not a complete Coxeter survivor")

    source = parse_source(SOURCE)
    gram: sp.Matrix = source["gram"]
    inverse = gram.inv()
    ample: sp.Matrix = source["ample"]
    normals = tuple(
        sp.Matrix(1, gram.rows, row)
        for row in cast["walls"]["primitive_restricted_normal_covectors"]
    )
    normal_matrix = sp.Matrix.vstack(*normals)
    facet_indices = tuple(
        int(certificate["restricted_ray_index"]) - 1
        for certificate in gate["certificates"]
    )
    redundancies = {
        int(certificate["restricted_ray_index"]) - 1: certificate
        for certificate in reduction["redundancy_certificates"]
    }
    if not (
        len(facet_indices) == int(reduction["irredundant_facet_count"])
        and len(redundancies) == int(reduction["redundant_inequality_count"])
        and set(facet_indices).isdisjoint(redundancies)
        and set(facet_indices) | set(redundancies) == set(range(len(normals)))
    ):
        raise ArithmeticError("facet/redundancy partition is not exhaustive")
    for target, certificate in redundancies.items():
        recovered = sum(
            (
                sp.Rational(term["coefficient"])
                * normals[int(term["restricted_ray_index"]) - 1]
                for term in certificate["terms"]
            ),
            sp.zeros(1, gram.rows),
        )
        if recovered != normals[target]:
            raise ArithmeticError("a conic redundancy certificate failed replay")

    facets = sp.Matrix.vstack(*(normals[index] for index in facet_indices))
    if facets.rank() != dimension + 1:
        raise ArithmeticError("facet normals do not span the fixed Lorentzian space")
    _, basis_positions = facets.T.rref()
    if len(basis_positions) != dimension + 1:
        raise ArithmeticError("failed to select a facet-normal coordinate basis")
    basis_normals = facets[list(basis_positions), :]
    basis_gram = basis_normals * inverse * basis_normals.T
    if exact_inertia(basis_gram) != (1, dimension, 0):
        raise ArithmeticError("facet coordinate basis has wrong Lorentzian inertia")
    inequality_matrix = facets * inverse * basis_normals.T
    if inequality_matrix.rank() != dimension + 1:
        raise ArithmeticError("inequality matrix is not full rank")

    ample_covector = ample * gram
    solution = sp.linsolve((basis_normals.T, ample_covector.T))
    ample_coordinates = tuple(next(iter(solution)))
    if any(not value.is_Rational for value in ample_coordinates):
        raise ArithmeticError("ample covector has no rational basis coordinates")
    ample_slacks = inequality_matrix * sp.Matrix(ample_coordinates)
    direct_slacks = facets * ample.T
    if ample_slacks != direct_slacks or any(value <= 0 for value in ample_slacks):
        raise ArithmeticError("ample class is not a strict cone interior point")

    for position, target in enumerate(facet_indices):
        normal = normals[target]
        norm = sp.Rational((normal * inverse * normal.T)[0])
        pairing = sp.Rational((ample_covector * inverse * normal.T)[0])
        witness = ample_covector - pairing / norm * normal
        witness_pairings = witness * inverse * normal_matrix.T
        if not (
            norm < 0
            and (witness * inverse * witness.T)[0] > 0
            and witness_pairings[target] == 0
            and all(
                value > 0
                for index, value in enumerate(witness_pairings)
                if index != target
            )
        ):
            raise ArithmeticError(f"strict facet witness {position + 1} failed replay")

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
        "dimension": dimension,
        "cast": cast,
        "facet_indices": facet_indices,
        "basis_positions": tuple(basis_positions),
        "basis_gram": basis_gram,
        "inequality_matrix": inequality_matrix,
        "inequality_text": "\n".join(lines) + "\n",
        "minimum_ample_slack": min(ample_slacks),
    }


def run_cdd(inequality_text: str, executable: str, timeout: int) -> str:
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as output:
        process = subprocess.Popen(
            [executable, "--rep"],
            stdin=subprocess.PIPE,
            stdout=output,
            stderr=subprocess.PIPE,
            text=True,
        )
        try:
            _, stderr = process.communicate(inequality_text, timeout=timeout)
        except subprocess.TimeoutExpired:
            process.kill()
            process.communicate()
            raise TimeoutError(f"cddlib exceeded {timeout} seconds") from None
        if process.returncode != 0:
            raise RuntimeError(f"cddexec_gmp failed: {stderr[-4000:]}")
        output.seek(0)
        return output.read()


def parse_cdd_rays(text: str, linear_dimension: int) -> tuple[tuple[int, ...], ...]:
    lines = text.splitlines()
    try:
        begin = lines.index("begin")
    except ValueError as error:
        raise ArithmeticError("cddlib output has no begin marker") from error
    header = lines[begin + 1].split()
    if len(header) != 3 or header[2] != "rational":
        raise ArithmeticError("cddlib output is not rational")
    count, width = map(int, header[:2])
    if width != linear_dimension + 1:
        raise ArithmeticError("cddlib output has wrong generator width")
    ray_lines = lines[begin + 2 : begin + 2 + count]
    if len(ray_lines) != count or lines[begin + 2 + count].strip() != "end":
        raise ArithmeticError("cddlib output is truncated")
    rays = []
    for line in ray_lines:
        tokens = line.split()
        if len(tokens) != width or Fraction(tokens[0]) != 0:
            raise ArithmeticError("homogeneous cone output contains a vertex row")
        rays.append(_primitive_integer_ray(tokens[1:]))
    if len(set(rays)) != len(rays):
        raise ArithmeticError("cddlib emitted duplicate projective rays")
    return tuple(rays)


def analyze_rays(
    cone: dict[str, Any],
    rays: tuple[tuple[int, ...], ...],
    *,
    schema_prefix: str = "ss_fermat_quartic_char3_v4",
) -> tuple[dict[str, Any], bytes]:
    dimension = int(cone["dimension"])
    linear_dimension = dimension + 1
    inequality, inequality_scale = _matrix_integral_scale(cone["inequality_matrix"])
    basis_gram, gram_scale = _matrix_integral_scale(cone["basis_gram"])
    active_histogram: Counter[int] = Counter()
    finite_active: Counter[int] = Counter()
    ideal_active: Counter[int] = Counter()
    hyperideal_active: Counter[int] = Counter()
    norm_histogram: Counter[int] = Counter()
    records: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
    coordinate_maximum = 0
    finite_count = ideal_count = hyperideal_count = 0
    for ray in rays:
        coordinate_maximum = max(coordinate_maximum, *(abs(value) for value in ray))
        slacks = tuple(
            sum(row[column] * ray[column] for column in range(linear_dimension))
            for row in inequality
        )
        if any(value < 0 for value in slacks):
            # Projective normalization may have flipped the cdd orientation.
            flipped = tuple(-value for value in ray)
            flipped_slacks = tuple(-value for value in slacks)
            if any(value < 0 for value in flipped_slacks):
                raise ArithmeticError("an exact cdd ray violates cone inequalities")
            ray = flipped
            slacks = flipped_slacks
        incidence = tuple(index for index, value in enumerate(slacks) if value == 0)
        if len(incidence) < dimension:
            raise ArithmeticError("an extreme ray has too few active facets")
        norm = sum(
            ray[left] * basis_gram[left][right] * ray[right]
            for left in range(linear_dimension)
            for right in range(linear_dimension)
        )
        norm_histogram[norm] += 1
        active_histogram[len(incidence)] += 1
        if norm > 0:
            finite_count += 1
            finite_active[len(incidence)] += 1
        elif norm == 0:
            ideal_count += 1
            ideal_active[len(incidence)] += 1
        else:
            hyperideal_count += 1
            hyperideal_active[len(incidence)] += 1
        records.append((ray, incidence))
    if len({incidence for _, incidence in records}) != len(records):
        raise ArithmeticError("two cdd rays have the same complete incidence set")

    records.sort(key=lambda item: (item[1], item[0]))
    header = {
        "schema": f"{schema_prefix}_h{dimension}_exact_rays_v1",
        "ray_count": len(records),
        "coordinate_basis_facet_positions": [
            index + 1 for index in cone["basis_positions"]
        ],
        "record_format": "[primitive_integer_ray, active_facet_positions]",
    }
    archive_lines = [json.dumps(header, sort_keys=True, separators=(",", ":"))]
    archive_lines.extend(
        json.dumps(
            [list(ray), [index + 1 for index in incidence]],
            separators=(",", ":"),
        )
        for ray, incidence in records
    )
    archive = gzip.compress(
        ("\n".join(archive_lines) + "\n").encode(), compresslevel=9, mtime=0
    )
    finite_volume = bool(rays) and hyperideal_count == 0
    compact = finite_volume and ideal_count == 0
    summary = {
        "schema": f"{schema_prefix}_h{dimension}_exact_cone_certificate_v1",
        "cone": {
            "hyperbolic_dimension": dimension,
            "linear_dimension": linear_dimension,
            "facet_count": int(cone["inequality_matrix"].rows),
            "inequality_rank": int(cone["inequality_matrix"].rank()),
            "coordinate_basis_facet_positions": [
                index + 1 for index in cone["basis_positions"]
            ],
            "minimum_strict_ample_slack": _fraction_text(cone["minimum_ample_slack"]),
            "full_dimensional": True,
            "pointed": True,
        },
        "exact_double_description": {
            "engine": "cddlib GMP rational double-description",
            "command": "cddexec_gmp --rep",
            "extreme_ray_count": len(rays),
            "finite_vertex_count": finite_count,
            "ideal_vertex_count": ideal_count,
            "hyperideal_vertex_count": hyperideal_count,
            "primitive_ray_coordinate_absolute_maximum": coordinate_maximum,
            "inequality_integral_scale": inequality_scale,
            "basis_gram_integral_scale": gram_scale,
            "lorentz_square_histogram_in_scaled_basis_gram": {
                str(key): value for key, value in sorted(norm_histogram.items())
            },
            "active_facet_count_histogram": {
                str(key): value for key, value in sorted(active_histogram.items())
            },
            "finite_vertex_active_facet_count_histogram": {
                str(key): value for key, value in sorted(finite_active.items())
            },
            "ideal_vertex_active_facet_count_histogram": {
                str(key): value for key, value in sorted(ideal_active.items())
            },
            "hyperideal_vertex_active_facet_count_histogram": {
                str(key): value for key, value in sorted(hyperideal_active.items())
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
            "criterion": "complete pointed chamber has no spacelike/hyperideal extreme ray",
        },
        "checks": {
            "all_global_conic_redundancies_replayed_exactly": True,
            "all_strict_facet_witnesses_replayed_exactly": True,
            "gmp_output_is_rational": True,
            "all_rays_satisfy_all_inequalities_exactly": True,
            "all_ray_incidence_sets_are_distinct": True,
            "no_hyperideal_extreme_ray": hyperideal_count == 0,
        },
    }
    return summary, archive


def paths_for_dimension(dimension: int) -> tuple[Path, Path, Path]:
    stem = f"ss_fermat_quartic_char3_v4_h{dimension}_exact_cone"
    root = ROOT / "data" / "enumerations"
    return root / f"{stem}.ine", root / f"{stem}_rays.jsonl.gz", root / f"{stem}_certificate.json"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dimension", type=int, choices=(10, 11), required=True)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--timeout", type=int, default=3600)
    parser.add_argument(
        "--cdd-executable", default=shutil.which("cddexec_gmp") or "cddexec_gmp"
    )
    args = parser.parse_args()
    cone = rebuild_cone(args.dimension)
    inequality_path, ray_path, summary_path = paths_for_dimension(args.dimension)
    inequality_path.parent.mkdir(parents=True, exist_ok=True)
    inequality_path.write_text(cone["inequality_text"], encoding="utf-8")
    print(
        f"prepared H{args.dimension}: {cone['inequality_matrix'].rows} facets, "
        f"rank {cone['inequality_matrix'].rank()} -> {inequality_path}",
        flush=True,
    )
    if args.prepare_only:
        return
    cdd_text = run_cdd(cone["inequality_text"], args.cdd_executable, args.timeout)
    rays = parse_cdd_rays(cdd_text, args.dimension + 1)
    print(f"cdd emitted {len(rays)} H{args.dimension} rays", flush=True)
    summary, archive = analyze_rays(cone, rays)
    ray_path.write_bytes(archive)
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary["decision"], indent=2, sort_keys=True), flush=True)
    print(
        json.dumps(summary["exact_double_description"], indent=2, sort_keys=True),
        flush=True,
    )
    print(f"wrote {summary_path} and {ray_path}")


if __name__ == "__main__":
    main()
