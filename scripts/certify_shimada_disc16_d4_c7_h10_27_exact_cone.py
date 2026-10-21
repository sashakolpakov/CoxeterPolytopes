"""Exact GMP cone certificate for the archived D0004/C007 H10 fold.

The archived-stabilizer screen leaves 27 genuine walls in a Lorentzian
rank-11 fixed lattice.  This verifier replays the saturated fixed-lattice
data, all 35 conic redundancies, all 27 strict timelike facet witnesses, and
the complete Coxeter pair gate.  It then asks cddlib's GMP rational engine
for every extreme ray of the homogeneous chamber cone and classifies every
ray by its exact Lorentz square.

No volume conclusion is inferred from the pair diagram alone: finite volume
is certified only after the complete pointed cone has no spacelike
(hyperideal) extreme ray.
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

from scripts.screen_shimada_disc16_atlas_archived_high_sections import ROOT
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


CASE_INPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_atlas_archived_high_sections"
    / "D0004_C007.json"
)
INEQUALITY_OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_d4_c7_h10_27_exact_cone.ine"
)
RAY_OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_d4_c7_h10_27_exact_cone_rays.jsonl.gz"
)
SUMMARY_OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_d4_c7_h10_27_exact_cone_certificate.json"
)
CASE_ID = "D0004_C007"
DIMENSION = 10
LINEAR_DIMENSION = 11
FACET_COUNT = 27
ORBIT_COUNT = 62
REDUNDANCY_COUNT = 35


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
    # Keep the cone orientation chosen by the inequalities.  Only one sign
    # of a projective ray belongs to this pointed cone.
    return tuple(value // divisor for value in integers)


def _rebuild_exact_cone(case_path: Path = CASE_INPUT) -> dict[str, Any]:
    case = json.loads(case_path.read_text(encoding="utf-8"))
    reduction = case["exact_reduction"]
    gate = case["coxeter_gate"]
    if not (
        case["id"] == CASE_ID
        and case["case"]["chamber_index"] == 4
        and case["case"]["class_index"] == 7
        and case["case"]["subgroup_order"] == 4
        and case["folded_walls"]["wall_orbit_count"] == ORBIT_COUNT
        and reduction["hyperbolic_facet_count"] == FACET_COUNT
        and reduction["conic_redundant_count"] == REDUNDANCY_COUNT
        and reduction["unresolved_conic_facet_count"] == 0
        and reduction["causally_inactive_conic_facet_count"] == 0
        and gate["all_pairs_coxeter"] is True
        and gate["non_coxeter_pair_count"] == 0
        and case["decision"]["classification"]
        == "coxeter_candidate_requires_volume_audit"
    ):
        raise ArithmeticError("the pinned D0004/C007 candidate changed")

    fixed_gram = _matrix(case["saturated_fixed_lattice"]["gram"])
    ample = sp.Matrix(
        [
            sp.Rational(value)
            for value in case["saturated_fixed_lattice"][
                "ample_fixed_coordinates"
            ]
        ]
    )
    if not (
        fixed_gram.shape == (LINEAR_DIMENSION, LINEAR_DIMENSION)
        and exact_inertia(fixed_gram) == (1, DIMENSION, 0)
        and fixed_gram.det() == case["saturated_fixed_lattice"]["determinant"]
        and (ample.T * fixed_gram * ample)[0]
        == sp.Rational(case["saturated_fixed_lattice"]["ample_square"])
        > 0
    ):
        raise ArithmeticError("the saturated Lorentz lattice failed exact replay")

    orbit_records = {
        int(record["orbit_index"]): record
        for record in case["folded_walls"]["orbit_records"]
    }
    if set(orbit_records) != set(range(1, ORBIT_COUNT + 1)):
        raise ArithmeticError("the folded wall-orbit list is incomplete")
    orbit_normals: dict[int, sp.Matrix] = {}
    orbit_pairings: dict[int, sp.Matrix] = {}
    for index, record in orbit_records.items():
        normal = sp.Matrix(
            [sp.Rational(value) for value in record["fixed_normal_coordinates"]]
        )
        pairing = sp.Matrix(
            [sp.Rational(value) for value in record["fixed_pairing_column"]]
        )
        if not (
            normal.rows == LINEAR_DIMENSION
            and pairing == fixed_gram * normal
            and (normal.T * fixed_gram * normal)[0]
            == sp.Rational(record["square"])
        ):
            raise ArithmeticError(f"folded orbit {index} failed exact replay")
        orbit_normals[index] = normal
        orbit_pairings[index] = pairing

    facet_indices = tuple(int(index) for index in reduction["hyperbolic_facet_orbit_indices"])
    conic_facets = tuple(int(index) for index in reduction["conic_facet_orbit_indices"])
    if not (facet_indices == conic_facets and len(facet_indices) == FACET_COUNT):
        raise ArithmeticError("the active hyperbolic facet list changed")
    facet_set = set(facet_indices)

    redundancies = {
        int(record["orbit_index"]): record
        for record in reduction["conic_redundancy_certificates"]
    }
    if not (
        len(redundancies) == REDUNDANCY_COUNT
        and facet_set.isdisjoint(redundancies)
        and facet_set | set(redundancies) == set(orbit_records)
    ):
        raise ArithmeticError("the facet/redundancy partition is not exhaustive")
    for target, certificate in redundancies.items():
        recovered = sum(
            (
                sp.Rational(term["coefficient"])
                * orbit_normals[int(term["conic_facet_orbit_index"])]
                for term in certificate["terms"]
            ),
            sp.zeros(LINEAR_DIMENSION, 1),
        )
        if recovered != orbit_normals[target]:
            raise ArithmeticError(
                f"conic redundancy for orbit {target} failed exact replay"
            )
        if any(
            int(term["conic_facet_orbit_index"]) not in facet_set
            or sp.Rational(term["coefficient"]) < 0
            for term in certificate["terms"]
        ):
            raise ArithmeticError("a redundancy uses a removed wall or negative weight")

    witnesses = {
        int(record["orbit_index"]): record
        for record in reduction["strict_timelike_facet_witnesses"]
    }
    if set(witnesses) != facet_set:
        raise ArithmeticError("the strict timelike witness list is incomplete")
    for target, certificate in witnesses.items():
        witness = sp.Matrix(
            [sp.Rational(value) for value in certificate["fixed_lattice_coordinates"]]
        )
        pairings = {
            index: (witness.T * orbit_pairings[index])[0]
            for index in orbit_records
        }
        if not (
            (witness.T * fixed_gram * witness)[0]
            == sp.Rational(certificate["witness_square"])
            > 0
            and pairings[target] == 0
            and all(value > 0 for index, value in pairings.items() if index != target)
            and min(value for index, value in pairings.items() if index != target)
            == sp.Rational(certificate["minimum_strict_pairing"])
        ):
            raise ArithmeticError(
                f"strict timelike witness for orbit {target} failed exact replay"
            )

    normals = sp.Matrix.hstack(*(orbit_normals[index] for index in facet_indices))
    inequality = sp.Matrix.vstack(
        *(orbit_pairings[index].T for index in facet_indices)
    )
    raw_facet_gram = normals.T * fixed_gram * normals
    if not (
        inequality.shape == (FACET_COUNT, LINEAR_DIMENSION)
        and inequality.rank() == LINEAR_DIMENSION
        and normals.rank() == LINEAR_DIMENSION
        and _matrix_payload(raw_facet_gram) == gate["facet_gram"]
        and _matrix_payload(-raw_facet_gram) == gate["positive_wall_gram"]
        and exact_inertia(-raw_facet_gram)
        == (DIMENSION, 1, FACET_COUNT - LINEAR_DIMENSION)
    ):
        raise ArithmeticError("the 27-wall Coxeter Gram failed exact replay")
    ample_slacks = inequality * ample
    if any(value <= 0 for value in ample_slacks):
        raise ArithmeticError("the ample class is not strictly inside the cone")

    lines = [
        "H-representation",
        "begin",
        f"{FACET_COUNT} {LINEAR_DIMENSION + 1} rational",
    ]
    lines.extend(
        "0 "
        + " ".join(
            _fraction_text(inequality[row, column])
            for column in range(LINEAR_DIMENSION)
        )
        for row in range(FACET_COUNT)
    )
    lines.append("end")
    return {
        "case": case,
        "fixed_gram": fixed_gram,
        "inequality_matrix": inequality,
        "inequality_text": "\n".join(lines) + "\n",
        "facet_indices": facet_indices,
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
    cone: dict[str, Any], rays: tuple[tuple[int, ...], ...], case_path: Path
) -> tuple[dict[str, Any], bytes]:
    inequality: sp.Matrix = cone["inequality_matrix"]
    fixed_gram: sp.Matrix = cone["fixed_gram"]
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
        square = sp.Rational((vector.T * fixed_gram * vector)[0])
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
        "schema": "shimada_disc16_d4_c7_h10_27_exact_cone_rays_v1",
        "coordinate_system": "saturated fixed-lattice basis in D0004_C007",
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
        "schema": "shimada_disc16_d4_c7_h10_27_exact_cone_certificate_v1",
        "source": {
            "case_path": str(case_path.relative_to(ROOT)),
            "case_sha256": _sha256(case_path),
            "case_id": CASE_ID,
            "scope": "subgroup of Shimada archived chamber stabilizer Aut(D)",
        },
        "cone": {
            "hyperbolic_dimension": DIMENSION,
            "linear_dimension": LINEAR_DIMENSION,
            "facet_count": FACET_COUNT,
            "facet_source_orbit_indices": list(cone["facet_indices"]),
            "inequality_rank": int(inequality.rank()),
            "lineality_dimension": 0,
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
            "finite_vertex_count": finite_count,
            "ideal_vertex_count": ideal_count,
            "hyperideal_vertex_count": hyperideal_count,
            "active_facet_count_histogram": {
                f"{kind}_{count}": value
                for (kind, count), value in sorted(active_histogram.items())
            },
            "lorentz_square_histogram": dict(sorted(lorentz_square_histogram.items())),
            "facet_extreme_ray_incidence_histogram": {
                str(index + 1): facet_incidence[index]
                for index in range(FACET_COUNT)
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
            "all_35_redundancies_replayed_as_nonnegative_sums_of_retained_facets": True,
            "all_27_strict_timelike_facet_witnesses_replayed": True,
            "all_351_facet_pairs_are_exact_Coxeter_relations": True,
            "gmp_output_rational": True,
            "all_962_rays_replayed_against_all_27_inequalities": len(rays) == 962,
            "all_ray_incidence_sets_distinct": True,
            "every_facet_has_extreme_ray_incidence": True,
            "no_hyperideal_extreme_ray": hyperideal_count == 0,
        },
    }
    expected_histogram = Counter(
        {
            ("finite", 10): 895,
            ("ideal", 11): 6,
            ("ideal", 12): 58,
            ("ideal", 13): 3,
        }
    )
    if not (
        len(rays) == 962
        and finite_count == 895
        and ideal_count == 67
        and hyperideal_count == 0
        and active_histogram == expected_histogram
    ):
        raise ArithmeticError("the exact D0004/C007 H10 cone census changed")
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
    parser.add_argument("--case", type=Path, default=CASE_INPUT)
    parser.add_argument(
        "--cdd-executable", default=shutil.which("cddexec_gmp") or "cddexec_gmp"
    )
    parser.add_argument("--ext-input", type=Path)
    parser.add_argument("--inequality-output", type=Path, default=INEQUALITY_OUTPUT)
    parser.add_argument("--ray-output", type=Path, default=RAY_OUTPUT)
    parser.add_argument("--summary-output", type=Path, default=SUMMARY_OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    cone = _rebuild_exact_cone(arguments.case)
    cdd_output = (
        arguments.ext_input.read_text(encoding="utf-8")
        if arguments.ext_input is not None
        else _run_cdd(cone["inequality_text"], arguments.cdd_executable)
    )
    rays = _parse_cdd_rays(cdd_output)
    summary, archive = _analyze_rays(cone, rays, arguments.case)
    _write_or_check_text(
        arguments.inequality_output, cone["inequality_text"], arguments.check_archive
    )
    _write_or_check_bytes(arguments.ray_output, archive, arguments.check_archive)
    _write_or_check_text(
        arguments.summary_output,
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        arguments.check_archive,
    )
    print(
        json.dumps(
            summary["decision"] | summary["exact_double_description"],
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
