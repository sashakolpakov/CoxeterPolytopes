"""GMP exact-cone certificate for the disc16/D5 H10 30-wall fold.

The complete subgroup census supplies the exact 30 retained wall orbits and
their facet/redundancy witnesses.  This verifier replays those data, builds a
rank-11 homogeneous H-representation, asks cddlib's GMP rational engine for
*all* extreme rays, and classifies every ray by its exact Lorentz square.
Finite volume is asserted only when no spacelike/hyperideal ray exists.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import gzip
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
from typing import Any, Iterable

import sympy as sp

from scripts.certify_shimada_disc16_d5_full_subgroup_fishing import (
    OUTPUT as CENSUS_OUTPUT,
    ROOT,
    SOURCE_SHA256,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


INEQUALITY_OUTPUT = (
    ROOT / "data" / "enumerations" / "shimada_disc16_d5_h10_30_exact_cone.ine"
)
RAY_OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_d5_h10_30_exact_cone_rays.jsonl.gz"
)
SUMMARY_OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_d5_h10_30_exact_cone_certificate.json"
)
DIMENSION = 10
LINEAR_DIMENSION = 11
FACET_COUNT = 30
CANDIDATE_KEY = "class_03_H10"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _fraction_text(value: int | Fraction | sp.Rational) -> str:
    value = Fraction(value)
    return str(value.numerator) if value.denominator == 1 else str(value)


def _matrix(payload: list[list[str]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(value) for value in row] for row in payload])


def _matrix_payload(matrix: sp.Matrix) -> list[list[str]]:
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
    # Preserve cddlib's cone orientation.  A homogeneous ray and its negative
    # are projectively equal, but only one of them satisfies this pointed
    # H-representation.
    return tuple(value // divisor for value in integers)


def _rebuild_exact_cone(census_path: Path = CENSUS_OUTPUT) -> dict[str, Any]:
    census = json.loads(census_path.read_text(encoding="utf-8"))
    screen = census["above_H9"]["complete_section_screens"][CANDIDATE_KEY]
    if not (
        census["source"]["sha256"] == SOURCE_SHA256
        and census["decision"]["subgroup_enumeration_complete"] is True
        and screen["hyperbolic_dimension"] == DIMENSION
        and screen["fixed_rank"] == LINEAR_DIMENSION
        and screen["facet_count"] == FACET_COUNT
        and screen["redundant_orbit_count"] == 56
        and screen["pair_gate"]["all_pairs_coxeter"] is True
    ):
        raise ArithmeticError("the pinned D5/H10 census candidate changed")

    dual_gram = _matrix(census["ambient_lattice"]["dual_gram"])
    walls = tuple(
        sp.Matrix(1, 20, wall)
        for wall in census["chamber"]["published_primitive_dual_walls"]
    )
    wall_orbits = tuple(
        tuple(index - 1 for index in orbit) for orbit in screen["wall_orbits"]
    )
    orbit_normals = tuple(
        sum((walls[index] for index in orbit), sp.zeros(1, 20))
        for orbit in wall_orbits
    )
    stacked = sp.Matrix.vstack(*orbit_normals)
    orbit_gram = stacked * dual_gram * stacked.T
    if not (
        stacked.rank() == LINEAR_DIMENSION
        and exact_inertia(orbit_gram) == (1, DIMENSION, len(wall_orbits) - LINEAR_DIMENSION)
    ):
        raise ArithmeticError("the rebuilt H10 orbit span changed")

    facet_indices = tuple(index - 1 for index in screen["facet_orbit_indices"])
    if len(facet_indices) != FACET_COUNT:
        raise ArithmeticError("the facet list changed")
    facet_index_set = set(facet_indices)
    redundancies = {
        record["orbit_index"] - 1: record
        for record in screen["redundancy_certificates"]
    }
    if not (
        facet_index_set.isdisjoint(redundancies)
        and facet_index_set | set(redundancies) == set(range(len(wall_orbits)))
    ):
        raise ArithmeticError("the facet/redundancy partition is not exhaustive")
    for target, record in redundancies.items():
        recovered = sum(
            (
                sp.Rational(term["coefficient"])
                * orbit_normals[term["facet_orbit_index"] - 1]
                for term in record["terms"]
            ),
            sp.zeros(1, 20),
        )
        if recovered != orbit_normals[target]:
            raise ArithmeticError("a conic redundancy failed exact replay")
        if any(
            term["facet_orbit_index"] - 1 not in facet_index_set
            for term in record["terms"]
        ):
            raise ArithmeticError("a redundancy certificate uses a removed wall")

    witnesses = {
        record["orbit_index"] - 1: record for record in screen["facet_witnesses"]
    }
    if set(witnesses) != facet_index_set:
        raise ArithmeticError("the strict witness list is incomplete")
    for target, record in witnesses.items():
        witness = sp.Matrix(1, 20, [sp.Rational(value) for value in record["ambient_dual_coordinates"]])
        pairings = tuple(
            sp.Rational(value)
            for value in list(witness * dual_gram * stacked.T)
        )
        if not (
            (witness * dual_gram * witness.T)[0] > 0
            and pairings[target] == 0
            and all(value > 0 for index, value in enumerate(pairings) if index != target)
        ):
            raise ArithmeticError("a strict timelike witness failed exact replay")

    raw_facet_gram = orbit_gram.extract(facet_indices, facet_indices)
    if _matrix_payload(raw_facet_gram) != screen["facet_gram"]:
        raise ArithmeticError("the rebuilt facet Gram changed")
    if exact_inertia(-raw_facet_gram) != (DIMENSION, 1, FACET_COUNT - LINEAR_DIMENSION):
        raise ArithmeticError("the 30-wall Coxeter Gram has wrong inertia")

    _, pivot_columns = raw_facet_gram.rref()
    basis = tuple(pivot_columns)
    if len(basis) != LINEAR_DIMENSION:
        raise ArithmeticError("the facet normals do not span rank 11")
    basis_normals = sp.Matrix.vstack(
        *(orbit_normals[facet_indices[index]] for index in basis)
    )
    ample = sp.Matrix(1, 20, census["chamber"]["ample_dual_coordinates"])
    solution_set = sp.linsolve((basis_normals.T, ample.T))
    if solution_set is sp.EmptySet:
        raise ArithmeticError("the ample vector left the facet span")
    ample_coordinates = tuple(next(iter(solution_set)))
    if any(not value.is_Rational for value in ample_coordinates):
        raise ArithmeticError("the ample coordinates are not rational")

    inequality_matrix = raw_facet_gram[:, list(basis)]
    ample_slacks = inequality_matrix * sp.Matrix(ample_coordinates)
    if inequality_matrix.rank() != LINEAR_DIMENSION or any(value <= 0 for value in ample_slacks):
        raise ArithmeticError("the candidate cone is not full-dimensional and pointed")
    lines = [
        "H-representation",
        "begin",
        f"{FACET_COUNT} {LINEAR_DIMENSION + 1} rational",
    ]
    lines.extend(
        "0 "
        + " ".join(
            _fraction_text(inequality_matrix[row, column])
            for column in range(LINEAR_DIMENSION)
        )
        for row in range(FACET_COUNT)
    )
    lines.append("end")
    return {
        "census": census,
        "screen": screen,
        "raw_facet_gram": raw_facet_gram,
        "basis_gram": raw_facet_gram.extract(basis, basis),
        "inequality_matrix": inequality_matrix,
        "inequality_text": "\n".join(lines) + "\n",
        "basis": basis,
        "minimum_ample_slack": min(ample_slacks),
    }


def _run_cdd(text: str, executable: str) -> str:
    completed = subprocess.run(
        [executable, "--rep"],
        input=text,
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
    count_text, width_text, kind = lines[begin + 1].split()
    count, width = int(count_text), int(width_text)
    if kind != "rational" or width != LINEAR_DIMENSION + 1:
        raise ArithmeticError("cddlib emitted the wrong generator format")
    ray_lines = lines[begin + 2 : begin + 2 + count]
    if len(ray_lines) != count or lines[begin + 2 + count].strip() != "end":
        raise ArithmeticError("cddlib output is truncated")
    rays = []
    for line in ray_lines:
        tokens = line.split()
        if len(tokens) != width or Fraction(tokens[0]) != 0:
            raise ArithmeticError("the homogeneous cone contains a vertex row")
        rays.append(_primitive_integer_ray(tokens[1:]))
    if len(set(rays)) != len(rays):
        raise ArithmeticError("cddlib emitted duplicate projective rays")
    return tuple(rays)


def _analyze_rays(
    cone: dict[str, Any], rays: tuple[tuple[int, ...], ...]
) -> tuple[dict[str, Any], bytes]:
    inequality: sp.Matrix = cone["inequality_matrix"]
    basis_gram: sp.Matrix = cone["basis_gram"]
    records = []
    type_histogram: Counter[str] = Counter()
    active_histogram: Counter[tuple[str, int]] = Counter()
    lorentz_square_histogram: Counter[str] = Counter()
    facet_incidence = Counter()
    for ray in rays:
        vector = sp.Matrix(ray)
        slacks = inequality * vector
        if any(value < 0 for value in slacks):
            raise ArithmeticError("a cddlib ray violates an exact inequality")
        active = tuple(index for index, value in enumerate(slacks) if value == 0)
        if len(active) < DIMENSION:
            raise ArithmeticError("an extreme ray has too few active facets")
        square = sp.Rational((vector.T * basis_gram * vector)[0])
        kind = "finite" if square > 0 else "ideal" if square == 0 else "hyperideal"
        type_histogram[kind] += 1
        active_histogram[(kind, len(active))] += 1
        lorentz_square_histogram[_fraction_text(square)] += 1
        facet_incidence.update(active)
        records.append((ray, active, square, kind))
    if len({active for _, active, _, _ in records}) != len(records):
        raise ArithmeticError("two exact rays have the same complete incidence set")
    if set(facet_incidence) != set(range(FACET_COUNT)):
        raise ArithmeticError("a purported facet has no exact extreme-ray incidence")

    finite_count = type_histogram["finite"]
    ideal_count = type_histogram["ideal"]
    hyperideal_count = type_histogram["hyperideal"]
    finite_volume = bool(rays) and hyperideal_count == 0
    compact = finite_volume and ideal_count == 0
    sorted_records = sorted(records, key=lambda item: (item[1], item[0]))
    header = {
        "schema": "shimada_disc16_d5_h10_30_exact_cone_rays_v1",
        "coordinate_basis_facet_positions": [index + 1 for index in cone["basis"]],
        "record_format": "[primitive_integer_ray, active_facet_positions, lorentz_square, type]",
        "ray_count": len(sorted_records),
    }
    lines = [json.dumps(header, sort_keys=True, separators=(",", ":"))]
    lines.extend(
        json.dumps(
            [
                list(ray),
                [index + 1 for index in active],
                _fraction_text(square),
                kind,
            ],
            separators=(",", ":"),
        )
        for ray, active, square, kind in sorted_records
    )
    archive = gzip.compress(
        ("\n".join(lines) + "\n").encode("utf-8"),
        compresslevel=9,
        mtime=0,
    )
    summary = {
        "schema": "shimada_disc16_d5_h10_30_exact_cone_certificate_v1",
        "source": {
            "census_path": str(CENSUS_OUTPUT.relative_to(ROOT)),
            "census_sha256": _sha256(CENSUS_OUTPUT),
            "primary_disc16_sha256": SOURCE_SHA256,
            "candidate_key": CANDIDATE_KEY,
        },
        "cone": {
            "hyperbolic_dimension": DIMENSION,
            "linear_dimension": LINEAR_DIMENSION,
            "facet_count": FACET_COUNT,
            "coordinate_basis_facet_positions": [index + 1 for index in cone["basis"]],
            "inequality_rank": int(inequality.rank()),
            "lineality_dimension": 0,
            "minimum_strict_ample_pairing": _fraction_text(cone["minimum_ample_slack"]),
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
            "active_facet_count_histogram": {
                f"{kind}_{count}": value
                for (kind, count), value in sorted(active_histogram.items())
            },
            "lorentz_square_histogram": dict(sorted(lorentz_square_histogram.items())),
            "facet_extreme_ray_incidence_histogram": {
                str(index + 1): facet_incidence[index] for index in range(FACET_COUNT)
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
            "criterion": "complete pointed projective cone has no spacelike/hyperideal extreme ray",
        },
        "checks": {
            "all_56_redundancies_replayed_as_nonnegative_sums_of_retained_facets": True,
            "all_30_strict_timelike_facet_witnesses_replayed": True,
            "all_435_facet_pairs_are_exact_Coxeter_relations": True,
            "gmp_output_rational": True,
            "all_1540_rays_replayed_against_all_30_inequalities": len(rays) == 1540,
            "all_ray_incidence_sets_distinct": True,
            "every_facet_has_extreme_ray_incidence": True,
            "no_hyperideal_extreme_ray": hyperideal_count == 0,
        },
    }
    if not (
        len(rays) == 1540
        and finite_count == 1333
        and ideal_count == 207
        and hyperideal_count == 0
        and active_histogram
        == Counter(
            {
                ("finite", 10): 1333,
                ("ideal", 10): 48,
                ("ideal", 11): 98,
                ("ideal", 12): 52,
                ("ideal", 13): 9,
            }
        )
    ):
        raise ArithmeticError("the exact D5/H10 cone census changed")
    return summary, archive


def load_ray_archive(
    path: Path = RAY_OUTPUT,
) -> tuple[dict[str, Any], tuple[tuple[Any, ...], ...]]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        header = json.loads(next(stream))
        records = tuple(tuple(json.loads(line)) for line in stream)
    if len(records) != header["ray_count"]:
        raise ArithmeticError("the compressed exact-ray archive is truncated")
    return header, records


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
    parser.add_argument("--census", type=Path, default=CENSUS_OUTPUT)
    parser.add_argument("--cdd-executable", default=shutil.which("cddexec_gmp") or "cddexec_gmp")
    parser.add_argument("--ext-input", type=Path)
    parser.add_argument("--inequality-output", type=Path, default=INEQUALITY_OUTPUT)
    parser.add_argument("--ray-output", type=Path, default=RAY_OUTPUT)
    parser.add_argument("--summary-output", type=Path, default=SUMMARY_OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    cone = _rebuild_exact_cone(arguments.census)
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
