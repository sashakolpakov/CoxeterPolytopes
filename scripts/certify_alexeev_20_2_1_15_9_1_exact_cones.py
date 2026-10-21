"""Independent GMP double-description volume audit of three Alexeev sections.

For every exact Coxeter survivor in the ``(20,2,1)`` and ``(15,9,1)``
fixed-section census, replay the saturated fixed lattice, restricted mirror
directions, conic redundancies, strict facet witnesses, and Coxeter pair
gate.  Then enumerate every extreme ray of the complete homogeneous chamber
cone with cddlib's GMP rational engine and classify it by exact Lorentz
square.  Finite volume is asserted only when no hyperideal ray exists.
"""

from __future__ import annotations

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

from scripts.screen_alexeev_20_2_1_15_9_1_fixed_sections import (
    OUTPUT as CENSUS,
    ROOT,
    _matrix,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_20_2_1_15_9_1_exact_cone_census.json"
)
RAY_DIRECTORY = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_20_2_1_15_9_1_exact_cone_rays"
)
INEQUALITY_DIRECTORY = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_20_2_1_15_9_1_exact_cone_inequalities"
)
SCHEMA = "coxeter.alexeev_20_2_1_15_9_1_exact_cone_census.v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _fraction_text(value: int | Fraction | sp.Rational) -> str:
    value = Fraction(value)
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
    return tuple(value // divisor for value in integers)


def _replay_case(case: dict[str, Any]) -> dict[str, Any]:
    if not (
        case["decision"]["classification"]
        == "coxeter_candidate_requires_volume_audit"
        and case["coxeter_gate"]["all_pairs_coxeter"] is True
        and case["coxeter_gate"]["non_coxeter_pair_count"] == 0
        and not case["exact_reduction"]["unresolved_orbit_indices"]
    ):
        raise ArithmeticError(f"case {case['id']} is not a resolved Coxeter survivor")
    fixed_gram = _matrix(case["saturated_fixed_lattice"]["gram"])
    rank = int(case["case"]["fixed_rank"])
    dimension = int(case["case"]["hyperbolic_dimension"])
    if not (
        fixed_gram.shape == (rank, rank)
        and rank == dimension + 1
        and exact_inertia(fixed_gram) == (1, dimension, 0)
        and int(fixed_gram.det()) == int(case["saturated_fixed_lattice"]["determinant"])
    ):
        raise ArithmeticError("fixed Lorentz lattice changed")

    records = {
        int(record["orbit_index"]): record
        for record in case["restricted_mirrors"]["records"]
    }
    orbit_pairings = {}
    orbit_normals = {}
    for orbit_index, record in records.items():
        pairing = sp.Matrix([sp.Rational(value) for value in record["fixed_pairing_column"]])
        normal = sp.Matrix(
            [sp.Rational(value) for value in record["rational_fixed_normal_coordinates"]]
        )
        if not (
            pairing == fixed_gram * normal
            and (normal.T * fixed_gram * normal)[0]
            == sp.Rational(record["square"])
        ):
            raise ArithmeticError(f"restricted mirror {orbit_index} failed replay")
        orbit_pairings[orbit_index] = pairing
        orbit_normals[orbit_index] = normal

    active = tuple(map(int, case["exact_reduction"]["active_facet_orbit_indices"]))
    conic_facets = tuple(map(int, case["exact_reduction"]["conic_facet_orbit_indices"]))
    if active != conic_facets:
        raise ArithmeticError("a survivor has a causally inactive conic facet")
    active_set = set(active)
    redundancies = {
        int(record["orbit_index"]): record
        for record in case["exact_reduction"]["redundancy_certificates"]
    }
    if not (
        active_set.isdisjoint(redundancies)
        and active_set | set(redundancies) == set(records)
    ):
        raise ArithmeticError("facet/redundancy partition is incomplete")
    for target, certificate in redundancies.items():
        recovered = sum(
            (
                sp.Rational(term["coefficient"])
                * orbit_pairings[int(term["facet_orbit_index"])]
                for term in certificate["terms"]
            ),
            sp.zeros(rank, 1),
        )
        if recovered != orbit_pairings[target]:
            raise ArithmeticError("a conic redundancy failed exact replay")

    inequality = sp.Matrix.vstack(*(orbit_pairings[index].T for index in active))
    normals = sp.Matrix.hstack(*(orbit_normals[index] for index in active))
    facet_gram = normals.T * fixed_gram * normals
    if not (
        inequality.rank() == rank
        and normals.rank() == rank
        and [[str(value) for value in row] for row in facet_gram.tolist()]
        == case["restricted_mirrors"]["facet_gram"]
    ):
        raise ArithmeticError("active Coxeter Gram failed replay")

    ample = sp.Matrix(
        [sp.Rational(value) for value in case["strict_fixed_interior"]["fixed_lattice_coordinates"]]
    )
    ample_slacks = inequality * ample
    if not (
        (ample.T * fixed_gram * ample)[0]
        == sp.Rational(case["strict_fixed_interior"]["averaged_fixed_square"])
        > 0
        and all(value > 0 for value in ample_slacks)
    ):
        raise ArithmeticError("strict fixed interior failed cone replay")

    witnesses = {
        int(record["orbit_index"]): record
        for record in case["exact_reduction"]["strict_facet_witnesses"]
    }
    if set(witnesses) != active_set:
        raise ArithmeticError("strict facet witness list is incomplete")
    all_pairing = {
        index: orbit_pairings[index] for index in records
    }
    for target, certificate in witnesses.items():
        witness = sp.Matrix(
            [sp.Rational(value) for value in certificate["fixed_lattice_coordinates"]]
        )
        pairings = {
            index: (witness.T * pairing)[0]
            for index, pairing in all_pairing.items()
        }
        if not (
            (witness.T * fixed_gram * witness)[0]
            == sp.Rational(certificate["witness_square"])
            > 0
            and pairings[target] == 0
            and all(value > 0 for index, value in pairings.items() if index != target)
        ):
            raise ArithmeticError("a strict hyperbolic facet witness failed replay")

    lines = [
        "H-representation",
        "begin",
        f"{len(active)} {rank + 1} rational",
    ]
    lines.extend(
        "0 "
        + " ".join(
            _fraction_text(inequality[row, column]) for column in range(rank)
        )
        for row in range(len(active))
    )
    lines.append("end")
    return {
        "fixed_gram": fixed_gram,
        "inequality": inequality,
        "active_orbit_indices": active,
        "dimension": dimension,
        "rank": rank,
        "inequality_text": "\n".join(lines) + "\n",
        "minimum_strict_slack": min(ample_slacks),
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


def _parse_cdd_rays(text: str, rank: int) -> tuple[tuple[int, ...], ...]:
    lines = text.splitlines()
    try:
        begin = lines.index("begin")
    except ValueError as error:
        raise ArithmeticError("cddlib output has no begin marker") from error
    count_text, width_text, kind = lines[begin + 1].split()
    count, width = int(count_text), int(width_text)
    if kind != "rational" or width != rank + 1:
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
        raise ArithmeticError("cddlib emitted duplicate rays")
    return tuple(rays)


def _analyze_case(
    case: dict[str, Any], cone: dict[str, Any], rays: tuple[tuple[int, ...], ...]
) -> tuple[dict[str, Any], bytes]:
    inequality: sp.Matrix = cone["inequality"]
    gram: sp.Matrix = cone["fixed_gram"]
    dimension = int(cone["dimension"])
    records = []
    type_histogram: Counter[str] = Counter()
    active_histogram: Counter[tuple[str, int]] = Counter()
    square_histogram: Counter[str] = Counter()
    facet_incidence = Counter()
    for ray in rays:
        vector = sp.Matrix(ray)
        slacks = inequality * vector
        if any(value < 0 for value in slacks):
            raise ArithmeticError("an extreme ray violates a chamber inequality")
        active = tuple(index for index, value in enumerate(slacks) if value == 0)
        if len(active) < dimension:
            raise ArithmeticError("an extreme ray has too few active facets")
        square = sp.Rational((vector.T * gram * vector)[0])
        kind = "finite" if square > 0 else "ideal" if square == 0 else "hyperideal"
        type_histogram[kind] += 1
        active_histogram[(kind, len(active))] += 1
        square_histogram[_fraction_text(square)] += 1
        facet_incidence.update(active)
        records.append((ray, active, square, kind))
    if len({active for _, active, _, _ in records}) != len(records):
        raise ArithmeticError("two rays have the same full incidence set")
    if set(facet_incidence) != set(range(inequality.rows)):
        raise ArithmeticError("a claimed facet has no extreme-ray incidence")

    finite_count = type_histogram["finite"]
    ideal_count = type_histogram["ideal"]
    hyperideal_count = type_histogram["hyperideal"]
    finite_volume = bool(rays) and hyperideal_count == 0
    compact = finite_volume and ideal_count == 0
    sorted_records = sorted(records, key=lambda item: (item[1], item[0]))
    header = {
        "schema": "coxeter.alexeev_fixed_section_exact_cone_rays.v1",
        "case_id": case["id"],
        "coordinate_system": "saturated fixed-lattice basis",
        "record_format": "[primitive_integer_ray,active_facet_positions,lorentz_square,type]",
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
        "id": case["id"],
        "ambient_key": case["ambient_key"],
        "hyperbolic_dimension": dimension,
        "facet_count": inequality.rows,
        "facet_source_orbit_indices": list(cone["active_orbit_indices"]),
        "fixed_lattice_determinant": case["saturated_fixed_lattice"]["determinant"],
        "fixed_lattice_smith_diagonal_absolute": case["saturated_fixed_lattice"][
            "smith_diagonal_absolute"
        ],
        "cone": {
            "linear_dimension": cone["rank"],
            "inequality_rank": int(inequality.rank()),
            "full_dimensional": True,
            "pointed": True,
            "minimum_strict_interior_slack": _fraction_text(
                cone["minimum_strict_slack"]
            ),
        },
        "exact_double_description": {
            "engine": "cddlib GMP rational double-description",
            "extreme_ray_count": len(rays),
            "finite_vertex_count": finite_count,
            "ideal_vertex_count": ideal_count,
            "hyperideal_vertex_count": hyperideal_count,
            "active_facet_count_histogram": {
                f"{kind}_{count}": value
                for (kind, count), value in sorted(active_histogram.items())
            },
            "lorentz_square_histogram": dict(sorted(square_histogram.items())),
            "facet_extreme_ray_incidence_histogram": {
                str(index + 1): facet_incidence[index]
                for index in range(inequality.rows)
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
            "criterion": "the complete pointed chamber cone has no hyperideal extreme ray",
        },
        "checks": {
            "all_redundancies_replayed": True,
            "all_strict_facet_witnesses_replayed": True,
            "all_facet_pairs_exactly_Coxeter": True,
            "gmp_output_rational": True,
            "all_rays_replayed_against_all_inequalities": True,
            "all_ray_incidence_sets_distinct": True,
            "every_facet_has_extreme_ray_incidence": True,
            "no_hyperideal_extreme_ray": hyperideal_count == 0,
        },
    }
    return summary, archive


def build(executable: str) -> dict[str, Any]:
    census_hash = _sha256(CENSUS)
    census = json.loads(CENSUS.read_text())
    candidates = census["survivors_requiring_volume_audit"]
    if len(candidates) != 3:
        raise ArithmeticError("expected exactly three angle-gate survivors")
    summaries = []
    RAY_DIRECTORY.mkdir(parents=True, exist_ok=True)
    INEQUALITY_DIRECTORY.mkdir(parents=True, exist_ok=True)
    for case in census["cases"]:
        if case["id"] not in candidates:
            continue
        cone = _replay_case(case)
        inequality_path = INEQUALITY_DIRECTORY / f"{case['id']}.ine"
        inequality_path.write_text(cone["inequality_text"])
        output = _run_cdd(cone["inequality_text"], executable)
        rays = _parse_cdd_rays(output, int(cone["rank"]))
        summary, archive = _analyze_case(case, cone, rays)
        ray_path = RAY_DIRECTORY / f"{case['id']}.jsonl.gz"
        ray_path.write_bytes(archive)
        summary["artifacts"] = {
            "inequality_path": str(inequality_path.relative_to(ROOT)),
            "inequality_sha256": _sha256(inequality_path),
            "ray_archive_path": str(ray_path.relative_to(ROOT)),
            "ray_archive_sha256": _sha256(ray_path),
        }
        summaries.append(summary)
        print(
            json.dumps(
                {
                    "id": case["id"],
                    **summary["decision"],
                    "rays": summary["exact_double_description"]["extreme_ray_count"],
                    "finite": summary["exact_double_description"]["finite_vertex_count"],
                    "ideal": summary["exact_double_description"]["ideal_vertex_count"],
                    "hyperideal": summary["exact_double_description"]["hyperideal_vertex_count"],
                }
            ),
            flush=True,
        )
    classification_histogram = Counter(row["decision"]["classification"] for row in summaries)
    return {
        "schema": SCHEMA,
        "source": {
            "fixed_section_census": str(CENSUS.relative_to(ROOT)),
            "fixed_section_census_sha256": census_hash,
        },
        "scope": {
            "case_count": len(summaries),
            "method": "independent complete GMP rational double description of each pointed homogeneous cone",
            "literal_Vinberg_audit_also_required_for_positive_ledger": True,
        },
        "summary": {
            "case_count": len(summaries),
            "classification_histogram": dict(sorted(classification_histogram.items())),
            "finite_volume_case_count": sum(row["decision"]["finite_volume"] for row in summaries),
            "compact_case_count": sum(row["decision"]["compact"] is True for row in summaries),
            "total_extreme_ray_count": sum(
                row["exact_double_description"]["extreme_ray_count"] for row in summaries
            ),
            "total_hyperideal_ray_count": sum(
                row["exact_double_description"]["hyperideal_vertex_count"] for row in summaries
            ),
        },
        "cases": summaries,
    }


def main() -> None:
    executable = shutil.which("cddexec_gmp") or "cddexec_gmp"
    payload = build(executable)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["summary"], sort_keys=True))
    print(OUTPUT)


if __name__ == "__main__":
    main()
