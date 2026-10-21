"""Hostile raw-data audit of the seven finite-Aut Enriques nef chambers.

This verifier is intentionally independent of
``certify_enriques_finite_aut_nef_chambers.py`` and of that certifier's JSON
and compressed ray outputs.  Its only chamber input is Shimada's pinned raw
GAP array ``L10L26compdata2.txt``.  It uses a record-splitting parser, chooses
its own lexicographically greedy wall bases, invokes cddlib/GMP directly, and
then calls only the repository's generic exact Vinberg implementation.

The published classification supplies the geometric assertion that the
listed nodal curves are the complete nef walls for the seven Enriques types.
This script adversarially checks every arithmetic and hyperbolic consequence
of those complete lists; it does not reprove the algebro-geometric
classification from the arrays alone.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
from typing import Any, Iterable

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

from coxeter_polytopes.exact.matrix import Inertia, inertia
from coxeter_polytopes.verify.vinberg_volume import (
    verify_vinberg_volume_criterion,
)


ROOT = Path(__file__).resolve().parents[1]
RAW_DATA = (
    ROOT
    / "data"
    / "external"
    / "shimada_l10l26note"
    / "L10L26compdata2.txt"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "enriques_finite_aut_nef_chambers_hostile_audit.json"
)

RAW_DATA_SHA256 = "218c8188e467bde29369cdeda3771d767c3f70105d8fb56def8d698b90e75a6c"
DIMENSION = 9
LATTICE_RANK = 10

TYPE_TO_RECORD = {
    "I": "12A",
    "II": "12B",
    "III": "20B",
    "IV": "20F",
    "V": "20A",
    "VI": "20E",
    "VII": "20C",
}

ALL_RAW_RECORDS = {
    "12A",
    "12B",
    "20A",
    "20B",
    "20C",
    "20D",
    "20E",
    "20F",
    "40A",
    "40B",
    "40C",
    "40D",
    "40E",
    "96A",
    "96B",
    "96C",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _extract_integer_array(record_text: str, field: str) -> list[Any]:
    """Extract one balanced GAP integer array from an isolated record."""

    matches = tuple(re.finditer(rf"\b{re.escape(field)}\s*:=", record_text))
    if len(matches) != 1:
        raise ArithmeticError(
            f"isolated GAP record has {len(matches)} occurrences of {field}"
        )
    opening = record_text.find("[", matches[0].end())
    if opening < 0:
        raise ArithmeticError(f"missing array after {field}")
    level = 0
    for position in range(opening, len(record_text)):
        character = record_text[position]
        if character == "[":
            level += 1
        elif character == "]":
            level -= 1
            if level == 0:
                parsed = ast.literal_eval(record_text[opening : position + 1])
                if not isinstance(parsed, list):
                    raise ArithmeticError(f"{field} is not a list")
                return parsed
    raise ArithmeticError(f"truncated balanced array for {field}")


def parse_gap_records(path: Path = RAW_DATA) -> dict[str, dict[str, list[Any]]]:
    """Parse the raw file by first isolating its sixteen ``rec`` blocks.

    This deliberately differs from the production parser, which aligns global
    field-occurrence arrays by index.
    """

    text = path.read_text(encoding="utf-8")
    chunks = re.split(r"\)\s*,\s*rec\(", text)
    if len(chunks) != 16:
        raise ArithmeticError(f"expected 16 GAP record blocks, found {len(chunks)}")

    records: dict[str, dict[str, list[Any]]] = {}
    for chunk in chunks:
        names = re.findall(r'\bname\s*:=\s*"([^"]+)"', chunk)
        if len(names) != 1:
            raise ArithmeticError(f"isolated record has name census {names}")
        name = names[0]
        if name in records:
            raise ArithmeticError(f"duplicate GAP record {name}")
        records[name] = {
            field: _extract_integer_array(chunk, field)
            for field in ("GramL10", "configuration", "walls")
        }
    if set(records) != ALL_RAW_RECORDS:
        raise ArithmeticError(
            f"raw record names changed: {sorted(set(records) ^ ALL_RAW_RECORDS)}"
        )
    return records


def _integer_rows(matrix: sp.Matrix) -> list[list[int]]:
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _fraction_text(value: int | Fraction | sp.Rational) -> str:
    fraction = Fraction(value)
    if fraction.denominator == 1:
        return str(fraction.numerator)
    return f"{fraction.numerator}/{fraction.denominator}"


def _gcd(values: Iterable[int]) -> int:
    result = 0
    for value in values:
        result = math.gcd(result, abs(int(value)))
    return result


def _primitive_integer_vector(tokens: Iterable[str]) -> tuple[int, ...]:
    fractions = tuple(Fraction(token) for token in tokens)
    denominator = math.lcm(*(value.denominator for value in fractions))
    integers = tuple(
        value.numerator * (denominator // value.denominator)
        for value in fractions
    )
    divisor = _gcd(integers)
    if divisor == 0:
        raise ArithmeticError("cddlib returned a zero generator")
    return tuple(value // divisor for value in integers)


def _primitive_ambient(vector: sp.Matrix) -> tuple[int, ...]:
    values = tuple(int(value) for value in vector)
    divisor = _gcd(values)
    if divisor == 0:
        raise ArithmeticError("zero ambient ray")
    return tuple(value // divisor for value in values)


def _greedy_wall_basis(walls: sp.Matrix) -> tuple[int, ...]:
    """Choose the first rank-increasing wall rows, independent of RREF pivots."""

    selected: list[int] = []
    rank = 0
    for index in range(walls.rows):
        candidate = walls[selected + [index], :]
        candidate_rank = int(candidate.rank())
        if candidate_rank > rank:
            selected.append(index)
            rank = candidate_rank
            if rank == LATTICE_RANK:
                break
    if len(selected) != LATTICE_RANK:
        raise ArithmeticError("wall rows do not contain a rank-ten basis")
    return tuple(selected)


def _cdd_input(inequalities: sp.Matrix) -> str:
    lines = [
        "H-representation",
        "begin",
        f"{inequalities.rows} {inequalities.cols + 1} rational",
    ]
    for row in range(inequalities.rows):
        coefficients = " ".join(
            _fraction_text(inequalities[row, column])
            for column in range(inequalities.cols)
        )
        lines.append(f"0 {coefficients}")
    lines.append("end")
    return "\n".join(lines) + "\n"


def _run_cdd(inequalities: sp.Matrix, executable: str) -> tuple[tuple[int, ...], ...]:
    completed = subprocess.run(
        [executable, "--rep"],
        input=_cdd_input(inequalities),
        text=True,
        capture_output=True,
        timeout=600,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "cddexec_gmp failed:\n"
            + completed.stdout[-3000:]
            + "\n"
            + completed.stderr[-3000:]
        )
    lines = completed.stdout.splitlines()
    try:
        begin = lines.index("begin")
    except ValueError as error:
        raise ArithmeticError("cddlib output has no begin marker") from error
    header = lines[begin + 1].split()
    if len(header) != 3:
        raise ArithmeticError(f"bad cddlib generator header {header}")
    count, width = int(header[0]), int(header[1])
    if header[2] != "rational" or width != inequalities.cols + 1:
        raise ArithmeticError("cddlib generator width/type changed")
    generator_lines = lines[begin + 2 : begin + 2 + count]
    if len(generator_lines) != count:
        raise ArithmeticError("truncated cddlib generator block")
    rays: list[tuple[int, ...]] = []
    for line in generator_lines:
        tokens = line.split()
        if len(tokens) != width or Fraction(tokens[0]) != 0:
            raise ArithmeticError("nonhomogeneous cddlib generator")
        rays.append(_primitive_integer_vector(tokens[1:]))
    if len(set(rays)) != len(rays):
        raise ArithmeticError("duplicate oriented ray from cddlib")
    return tuple(sorted(rays))


def _reflection_gate(root: sp.Matrix, lattice_gram: sp.Matrix) -> bool:
    reflection = sp.eye(LATTICE_RANK) + root * (root.T * lattice_gram)
    return bool(
        all(value.q == 1 for value in reflection)
        and reflection.T * lattice_gram * reflection == lattice_gram
        and reflection * root == -root
        and reflection * reflection == sp.eye(LATTICE_RANK)
    )


def _wall_smith_data(walls: sp.Matrix) -> tuple[list[int], int]:
    smith = smith_normal_form(walls, domain=sp.ZZ)
    diagonal = [
        abs(int(smith[index, index])) for index in range(LATTICE_RANK)
    ]
    if not all(diagonal):
        raise ArithmeticError("zero Smith invariant in the rank-ten wall span")
    return diagonal, math.prod(diagonal)


def _pair_histogram(configuration: sp.Matrix) -> dict[str, int]:
    histogram = Counter(
        int(configuration[left, right])
        for left in range(configuration.rows)
        for right in range(left)
    )
    return {str(value): int(histogram[value]) for value in sorted(histogram)}


def _cone_audit(
    configuration: sp.Matrix,
    walls: sp.Matrix,
    lattice_gram: sp.Matrix,
    *,
    cdd_executable: str,
) -> dict[str, Any]:
    basis_indices = _greedy_wall_basis(walls)
    basis_walls = walls[list(basis_indices), :]
    basis_gram = basis_walls * lattice_gram * basis_walls.T
    if basis_gram.det() == 0:
        raise ArithmeticError("greedy wall basis has singular restricted Gram")
    inequalities = walls * lattice_gram * basis_walls.T
    if inequalities != configuration[:, list(basis_indices)]:
        raise ArithmeticError("ambient and configuration inequalities disagree")
    if inequalities.rank() != LATTICE_RANK:
        raise ArithmeticError("inequality matrix is not full-dimensional")

    rays = _run_cdd(inequalities, cdd_executable)
    ray_records: list[dict[str, Any]] = []
    kind_histogram: Counter[str] = Counter()
    activity_histogram: Counter[tuple[str, int]] = Counter()
    incident_rays: list[list[tuple[int, ...]]] = [
        [] for _ in range(walls.rows)
    ]
    incidence_sets: set[tuple[int, ...]] = set()

    for coordinates in rays:
        coordinate_vector = sp.Matrix(coordinates)
        ambient_vector = basis_walls.T * coordinate_vector
        slacks = walls * lattice_gram * ambient_vector
        if slacks != inequalities * coordinate_vector:
            raise ArithmeticError("two exact slack computations disagree")
        if any(value < 0 for value in slacks):
            raise ArithmeticError("cdd ray violates a chamber inequality")
        active = tuple(index for index, value in enumerate(slacks) if value == 0)
        if len(active) < DIMENSION:
            raise ArithmeticError("extreme ray lies on fewer than nine walls")
        if active in incidence_sets:
            raise ArithmeticError("two rays have the same complete incidence set")
        incidence_sets.add(active)

        square_from_basis = (coordinate_vector.T * basis_gram * coordinate_vector)[0]
        square_from_ambient = (ambient_vector.T * lattice_gram * ambient_vector)[0]
        if square_from_basis != square_from_ambient:
            raise ArithmeticError("ray square disagrees between coordinate systems")
        square = sp.Rational(square_from_ambient)
        kind = "finite" if square > 0 else "ideal" if square == 0 else "hyperideal"
        kind_histogram[kind] += 1
        activity_histogram[(kind, len(active))] += 1
        primitive_ambient = _primitive_ambient(ambient_vector)
        for facet in active:
            incident_rays[facet].append(primitive_ambient)
        ray_records.append(
            {
                "basis_coordinates": list(coordinates),
                "primitive_L10_vector": list(primitive_ambient),
                "active_facets": [index + 1 for index in active],
                "square": _fraction_text(square),
                "kind": kind,
            }
        )

    face_ranks = [int(sp.Matrix(rows).rank()) for rows in incident_rays]
    if any(rank != DIMENSION for rank in face_ranks):
        raise ArithmeticError(
            f"a proposed facet is redundant or has wrong face rank: {face_ranks}"
        )

    interior_coordinates = sp.Matrix(
        [sum(ray[column] for ray in rays) for column in range(LATTICE_RANK)]
    )
    interior_ambient = basis_walls.T * interior_coordinates
    interior_slacks = walls * lattice_gram * interior_ambient
    interior_square = sp.Rational(
        (interior_ambient.T * lattice_gram * interior_ambient)[0]
    )
    if not (all(value > 0 for value in interior_slacks) and interior_square > 0):
        raise ArithmeticError("sum of exact rays is not a strict positive witness")

    time_pairings = []
    for record in ray_records:
        ray = sp.Matrix(record["primitive_L10_vector"])
        time_pairings.append(
            sp.Rational((ray.T * lattice_gram * interior_ambient)[0])
        )
    if not all(value > 0 for value in time_pairings):
        raise ArithmeticError("causal rays do not occupy one future time sheet")

    ray_records.sort(
        key=lambda record: (record["active_facets"], record["primitive_L10_vector"])
    )
    return {
        "engine": "direct cddexec_gmp rational H-to-V conversion",
        "coordinate_basis_wall_numbers": [index + 1 for index in basis_indices],
        "extreme_ray_count": len(ray_records),
        "finite_vertex_count": int(kind_histogram["finite"]),
        "ideal_vertex_count": int(kind_histogram["ideal"]),
        "hyperideal_ray_count": int(kind_histogram["hyperideal"]),
        "active_facet_count_histogram": {
            f"{kind}_{count}": int(number)
            for (kind, count), number in sorted(activity_histogram.items())
        },
        "facet_face_ranks": face_ranks,
        "all_facets_irredundant": all(rank == DIMENSION for rank in face_ranks),
        "strict_interior_L10_vector": [int(value) for value in interior_ambient],
        "minimum_strict_wall_pairing": _fraction_text(min(interior_slacks)),
        "strict_interior_square": _fraction_text(interior_square),
        "all_rays_in_one_positive_time_sheet": True,
        "minimum_ray_interior_pairing": _fraction_text(min(time_pairings)),
        "finite_volume_from_exact_rays": kind_histogram["hyperideal"] == 0,
        "compact_from_exact_rays": (
            kind_histogram["hyperideal"] == 0 and kind_histogram["ideal"] == 0
        ),
        "rays": ray_records,
    }


def _vinberg_audit(configuration: sp.Matrix) -> dict[str, Any]:
    normal_gram = [[-int(value) for value in row] for row in configuration.tolist()]
    verification = verify_vinberg_volume_criterion(
        normal_gram,
        DIMENSION,
        simple_roots_complete=True,
    )
    endpoint_histogram: Counter[str] = Counter()
    for audit in verification.edge_audits:
        endpoint_histogram[
            f"F{len(audit.finite_endpoints)}I{len(audit.ideal_endpoints)}"
        ] += 1
    invalid_parabolic = sum(
        not audit.is_rank_dimension_minus_one_parabolic
        for audit in verification.parabolic_neighborhood_audits
    )
    return {
        "method": (
            "literal elliptic-rank-8 edge endpoints plus complete affine "
            "zero-neighborhoods"
        ),
        "simple_roots_complete_external_classification_gate": (
            verification.simple_roots_complete
        ),
        "decision": verification.decision,
        "finite_volume": verification.finite_volume,
        "compact": verification.compact,
        "finite_vertex_count": verification.finite_vertex_count,
        "ideal_vertex_count": verification.ideal_vertex_count,
        "elliptic_rank_8_edge_count": verification.elliptic_edge_count,
        "edge_endpoint_histogram": dict(sorted(endpoint_histogram.items())),
        "bad_finite_volume_edge_count": len(verification.finite_volume_bad_edges),
        "bad_compact_edge_count": len(verification.compact_bad_edges),
        "affine_zero_neighborhood_audit_count": len(
            verification.parabolic_neighborhood_audits
        ),
        "invalid_affine_zero_neighborhood_count": int(invalid_parabolic),
    }


def _analyze_type(
    type_name: str,
    record_name: str,
    raw_record: dict[str, list[Any]],
    lattice_gram: sp.Matrix,
    *,
    cdd_executable: str,
) -> dict[str, Any]:
    configuration = sp.Matrix(raw_record["configuration"])
    walls = sp.Matrix(raw_record["walls"])
    facets = walls.rows
    if configuration.shape != (facets, facets) or walls.cols != LATTICE_RANK:
        raise ArithmeticError(f"shape mismatch for {type_name}/{record_name}")
    if configuration != configuration.T:
        raise ArithmeticError(f"nonsymmetric configuration for {type_name}")
    if configuration != walls * lattice_gram * walls.T:
        raise ArithmeticError(f"W*GramL10*W^T mismatch for {type_name}")
    if any(configuration[index, index] != -2 for index in range(facets)):
        raise ArithmeticError(f"non-root diagonal for {type_name}")
    if walls.rank() != LATTICE_RANK:
        raise ArithmeticError(f"walls fail to span rank ten for {type_name}")

    primitive = [
        _gcd(int(value) for value in walls.row(index)) == 1
        for index in range(facets)
    ]
    if not all(primitive):
        raise ArithmeticError(f"imprimitive wall in {type_name}")
    if not all(
        _reflection_gate(sp.Matrix(walls.row(index)).T, lattice_gram)
        for index in range(facets)
    ):
        raise ArithmeticError(f"nonintegral reflection in {type_name}")

    off_diagonal = {
        int(configuration[left, right])
        for left in range(facets)
        for right in range(left)
    }
    if not off_diagonal <= {0, 1, 2}:
        raise ArithmeticError(
            f"non-Coxeter intersection entry in {type_name}: {off_diagonal}"
        )
    normal_inertia = inertia(_integer_rows(-configuration))
    expected_inertia = Inertia(DIMENSION, 1, facets - LATTICE_RANK)
    if normal_inertia != expected_inertia:
        raise ArithmeticError(
            f"normal Gram inertia {normal_inertia} != {expected_inertia}"
        )

    smith_diagonal, wall_index = _wall_smith_data(walls)
    pair_histogram = _pair_histogram(configuration)
    cone = _cone_audit(
        configuration,
        walls,
        lattice_gram,
        cdd_executable=cdd_executable,
    )
    vinberg = _vinberg_audit(configuration)
    if not (
        cone["finite_vertex_count"] == vinberg["finite_vertex_count"]
        and cone["ideal_vertex_count"] == vinberg["ideal_vertex_count"]
        and cone["hyperideal_ray_count"] == 0
        and cone["finite_volume_from_exact_rays"] is True
        and cone["compact_from_exact_rays"] is False
        and vinberg["decision"] == "finite_volume_noncompact"
        and vinberg["finite_volume"] is True
        and vinberg["compact"] is False
        and vinberg["bad_finite_volume_edge_count"] == 0
        and vinberg["invalid_affine_zero_neighborhood_count"] == 0
    ):
        raise ArithmeticError(f"cone/Vinberg disagreement for {type_name}")

    return {
        "nikulin_kondo_type": type_name,
        "raw_shimada_record": record_name,
        "facet_count": facets,
        "raw_reconstruction": {
            "configuration_equals_W_GramL10_Wt": True,
            "wall_rank": int(walls.rank()),
            "all_walls_primitive": True,
            "all_minus_two_reflections_integral": True,
            "wall_span_smith_diagonal": smith_diagonal,
            "wall_span_index_in_L10": wall_index,
        },
        "gram_inertia": {
            "positive": normal_inertia.positive,
            "negative": normal_inertia.negative,
            "zero": normal_inertia.zero,
        },
        "coxeter_pair_histogram": {
            "orthogonal_m2": int(pair_histogram.get("0", 0)),
            "angle_pi_over_3_m3": int(pair_histogram.get("1", 0)),
            "parallel_infinity": int(pair_histogram.get("2", 0)),
        },
        "exact_cone": cone,
        "literal_vinberg": vinberg,
        "classification": "finite_volume_noncompact_H9_Coxeter_polytope",
    }


def certify(
    raw_path: Path = RAW_DATA,
    *,
    cdd_executable: str = "cddexec_gmp",
) -> dict[str, Any]:
    if _sha256(raw_path) != RAW_DATA_SHA256:
        raise ArithmeticError("pinned Shimada raw-data SHA-256 changed")
    if shutil.which(cdd_executable) is None:
        raise RuntimeError(f"missing exact cone engine {cdd_executable}")
    records = parse_gap_records(raw_path)

    lattice_gram = sp.Matrix(records["12A"]["GramL10"])
    if not (
        lattice_gram.shape == (LATTICE_RANK, LATTICE_RANK)
        and lattice_gram == lattice_gram.T
        and lattice_gram.det() == -1
        and inertia(_integer_rows(lattice_gram)) == Inertia(1, 9, 0)
        and all(
            sp.Matrix(record["GramL10"]) == lattice_gram
            for record in records.values()
        )
    ):
        raise ArithmeticError("raw files do not share the even unimodular L10 Gram")
    if any(lattice_gram[index, index] % 2 for index in range(LATTICE_RANK)):
        raise ArithmeticError("raw L10 Gram is not even")

    types = [
        _analyze_type(
            type_name,
            record_name,
            records[record_name],
            lattice_gram,
            cdd_executable=cdd_executable,
        )
        for type_name, record_name in TYPE_TO_RECORD.items()
    ]

    configuration_20c = records["20C"]["configuration"]
    configuration_20d = records["20D"]["configuration"]
    if configuration_20c != configuration_20d:
        raise ArithmeticError("published 20C and 20D configurations differ")
    walls_20d = sp.Matrix(records["20D"]["walls"])
    configuration_20d_matrix = sp.Matrix(configuration_20d)
    if not (
        configuration_20d_matrix == walls_20d * lattice_gram * walls_20d.T
        and walls_20d.rank() == LATTICE_RANK
        and all(
            _gcd(int(value) for value in walls_20d.row(index)) == 1
            for index in range(walls_20d.rows)
        )
        and inertia(_integer_rows(-configuration_20d_matrix))
        == Inertia(DIMENSION, 1, walls_20d.rows - LATTICE_RANK)
    ):
        raise ArithmeticError("alternate 20D raw wall realization fails replay")
    smith_20d, index_20d = _wall_smith_data(walls_20d)

    signatures = {
        record["nikulin_kondo_type"]: (
            record["facet_count"],
            record["coxeter_pair_histogram"]["orthogonal_m2"],
            record["coxeter_pair_histogram"]["angle_pi_over_3_m3"],
            record["coxeter_pair_histogram"]["parallel_infinity"],
        )
        for record in types
    }
    if len(set(signatures.values())) != 7:
        raise ArithmeticError("pair histograms fail to separate the seven images")

    return {
        "schema": "enriques_finite_aut_nef_chambers_hostile_audit_v1",
        "independence_boundary": {
            "sole_chamber_input": str(raw_path.relative_to(ROOT)),
            "production_certifier_imported": False,
            "production_summary_used_as_input": False,
            "production_ray_archive_used_as_input": False,
            "raw_record_parser": "independent record split plus balanced local arrays",
            "cone_basis": "independent lexicographically greedy wall-row basis",
            "generic_exact_modules_used": [
                "coxeter_polytopes.exact.matrix.inertia",
                "coxeter_polytopes.verify.vinberg_volume",
                "cddexec_gmp",
            ],
        },
        "raw_source": {
            "path": str(raw_path.relative_to(ROOT)),
            "sha256": RAW_DATA_SHA256,
            "record_count": len(records),
            "record_names": sorted(records),
        },
        "classification_input": {
            "type_to_raw_record": TYPE_TO_RECORD,
            "complete_nef_wall_lists": (
                "external Kondo-Nikulin-Martin-Shimada classification input"
            ),
        },
        "lattice": {
            "name": "L10 = U + E8(-1) = Num(Y)",
            "gram": _integer_rows(lattice_gram),
            "rank": LATTICE_RANK,
            "determinant": int(lattice_gram.det()),
            "even": True,
            "unimodular": True,
            "inertia": {"positive": 1, "negative": 9, "zero": 0},
        },
        "types": types,
        "cross_type_checks": {
            "20C_configuration_equals_20D_literally": True,
            "20C_wall_vectors_equal_20D_literally": (
                records["20C"]["walls"] == records["20D"]["walls"]
            ),
            "20D_duplicate_raw_replay": {
                "configuration_equals_W_GramL10_Wt": True,
                "wall_rank": int(walls_20d.rank()),
                "all_walls_primitive": True,
                "normal_gram_has_H9_inertia": True,
                "wall_span_smith_diagonal": smith_20d,
                "wall_span_index_in_L10": index_20d,
            },
            "pair_histogram_signatures": {
                type_name: list(signature)
                for type_name, signature in signatures.items()
            },
            "all_seven_pair_histograms_distinct": True,
            "facet_count_multiset": sorted(record["facet_count"] for record in types),
            "no_fourteen_facet_image": all(
                record["facet_count"] != 14 for record in types
            ),
        },
        "global_checks": {
            "all_W_GramL10_Wt_replays_pass": True,
            "all_walls_primitive_integral_minus_two_reflection_roots": True,
            "all_normal_grams_have_H9_inertia": True,
            "all_exact_cones_full_dimensional_and_irredundant": True,
            "all_exact_rays_share_one_positive_time_sheet": True,
            "all_exact_cones_have_no_hyperideal_ray": True,
            "all_literal_vinberg_tests_pass": True,
            "all_finite_volume": True,
            "all_noncompact": True,
        },
        "conclusion": (
            "Conditional on the published completeness of the seven nodal-curve "
            "lists, the raw arrays independently yield exactly seven distinct "
            "finite-volume noncompact full-nef Coxeter chambers in H9: two with "
            "12 facets and five with 20 facets."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", type=Path, default=RAW_DATA)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--cdd", default="cddexec_gmp")
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()

    result = certify(arguments.raw, cdd_executable=arguments.cdd)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file():
            raise SystemExit(f"missing hostile archive: {arguments.output}")
        if arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"hostile archive differs: {arguments.output}")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")

    print("archive_check=" + ("PASS" if arguments.check_archive else "WRITTEN"))
    print("all_checks_pass=true")
    print("types=I,II,III,IV,V,VI,VII")
    print("output_sha256=" + hashlib.sha256(rendered.encode("utf-8")).hexdigest())


if __name__ == "__main__":
    main()
