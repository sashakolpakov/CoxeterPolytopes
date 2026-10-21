"""Certify the seven finite-automorphism Enriques nef Coxeter chambers.

The mathematical input is deliberately split into two parts.

* Kondō--Nikulin--Martin classify complex Enriques surfaces with finite
  automorphism group by seven *complete* configurations of smooth rational
  curves.
* Shimada's official ``L10L26compdata2.txt`` gives exact wall vectors in the
  even unimodular Enriques lattice ``L10`` and their intersection matrices.

This script parses that raw GAP file without evaluating GAP code, verifies the
wall vectors against the lattice, and then applies two independent exact
hyperbolic-volume tests:

1. cddlib/GMP enumerates every extreme ray of the homogeneous nef cone; and
2. the literal Vinberg edge-endpoint/full-affine-neighbourhood criterion is
   applied to the complete Coxeter diagram.

The output is a finite seven-image transfer certificate.  It does not claim
that the underlying surfaces or reflection configurations are new.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
from fractions import Fraction
import gzip
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
    / "known_examples"
    / "enriques_finite_automorphism_nef_chambers.json"
)
RAY_OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "enriques_finite_automorphism_nef_chamber_rays.jsonl.gz"
)

DIMENSION = 9
LATTICE_RANK = 10
RAW_DATA_SHA256 = "218c8188e467bde29369cdeda3771d767c3f70105d8fb56def8d698b90e75a6c"

TYPE_SPECS: tuple[dict[str, Any], ...] = (
    {
        "type": "I",
        "record": "12A",
        "facets": 12,
        "moduli": "one-dimensional irreducible family",
        "surface_automorphism_group": "D4 (order 8)",
        "finite_vertices": 14,
        "ideal_vertices": 9,
        "elliptic_edges": 114,
    },
    {
        "type": "II",
        "record": "12B",
        "facets": 12,
        "moduli": "one-dimensional irreducible family",
        "surface_automorphism_group": "S4 (order 24)",
        "finite_vertices": 38,
        "ideal_vertices": 17,
        "elliptic_edges": 270,
    },
    {
        "type": "III",
        "record": "20B",
        "facets": 20,
        "moduli": "unique complex surface",
        "surface_automorphism_group": "(Z/4 x (Z/2)^2) semidirect D4",
        "finite_vertices": 200,
        "ideal_vertices": 44,
        "elliptic_edges": 1301,
    },
    {
        "type": "IV",
        "record": "20F",
        "facets": 20,
        "moduli": "unique complex surface",
        "surface_automorphism_group": "(Z/2)^4 semidirect (Z/5 semidirect Z/4)",
        "finite_vertices": 720,
        "ideal_vertices": 87,
        "elliptic_edges": 4405,
    },
    {
        "type": "V",
        "record": "20A",
        "facets": 20,
        "moduli": "unique complex surface",
        "surface_automorphism_group": "S4 x Z/2 (order 48)",
        "finite_vertices": 192,
        "ideal_vertices": 29,
        "elliptic_edges": 1182,
    },
    {
        "type": "VI",
        "record": "20E",
        "facets": 20,
        "moduli": "unique complex surface",
        "surface_automorphism_group": "S5 (order 120)",
        "finite_vertices": 600,
        "ideal_vertices": 57,
        "elliptic_edges": 3420,
    },
    {
        "type": "VII",
        "record": "20C",
        "alternate_records": ["20D"],
        "facets": 20,
        "moduli": "unique complex surface",
        "surface_automorphism_group": "S5 (order 120)",
        "finite_vertices": 450,
        "ideal_vertices": 51,
        "elliptic_edges": 2490,
    },
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _balanced_list(text: str, marker_position: int) -> tuple[list[Any], int]:
    """Parse the first GAP/Python-compatible integer list after a marker."""

    start = text.index("[", marker_position)
    depth = 0
    for end in range(start, len(text)):
        if text[end] == "[":
            depth += 1
        elif text[end] == "]":
            depth -= 1
            if depth == 0:
                value = ast.literal_eval(text[start : end + 1])
                if not isinstance(value, list):
                    raise ArithmeticError("a GAP matrix field is not a list")
                return value, end + 1
    raise ArithmeticError("an array in the raw GAP file is truncated")


def parse_raw_gap_records(path: Path = RAW_DATA) -> dict[str, dict[str, Any]]:
    """Extract only integer matrix/list fields from Shimada's GAP records."""

    text = path.read_text(encoding="utf-8")
    field_patterns = {
        field: tuple(re.finditer(rf"\b{field}\s*:=", text))
        for field in ("GramL10", "configuration", "walls")
    }
    names = tuple(
        match.group(1)
        for match in re.finditer(r'\bname\s*:=\s*"([^"]+)"', text)
    )
    counts = {field: len(matches) for field, matches in field_patterns.items()}
    if set(counts.values()) != {16} or len(names) != 16:
        raise ArithmeticError(
            f"the raw GAP record census changed: fields={counts}, names={len(names)}"
        )

    records: dict[str, dict[str, Any]] = {}
    for index, name in enumerate(names):
        record: dict[str, Any] = {}
        positions = []
        for field, matches in field_patterns.items():
            value, end = _balanced_list(text, matches[index].end())
            record[field] = value
            positions.append((matches[index].start(), end))
        name_position = text.index(f'name := "{name}"')
        if not (
            positions[0][0] < positions[1][0] < name_position < positions[2][0]
        ):
            raise ArithmeticError(f"raw GAP field order changed in record {name}")
        if name in records:
            raise ArithmeticError(f"duplicate raw GAP record {name}")
        records[name] = record
    return records


def _fraction_text(value: int | Fraction | sp.Rational) -> str:
    fraction = Fraction(value)
    if fraction.denominator == 1:
        return str(fraction.numerator)
    return f"{fraction.numerator}/{fraction.denominator}"


def _integer_rows(matrix: sp.Matrix) -> list[list[int]]:
    """Convert a SymPy integral matrix to plain JSON/exact-library integers."""

    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
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
        raise ArithmeticError("cddlib emitted the zero vector")
    return tuple(value // divisor for value in integers)


def _cdd_h_representation(inequalities: sp.Matrix) -> str:
    if inequalities.cols != LATTICE_RANK:
        raise ArithmeticError("the cdd inequality matrix has wrong width")
    lines = [
        "H-representation",
        "begin",
        f"{inequalities.rows} {LATTICE_RANK + 1} rational",
    ]
    lines.extend(
        "0 "
        + " ".join(
            _fraction_text(inequalities[row, column])
            for column in range(inequalities.cols)
        )
        for row in range(inequalities.rows)
    )
    lines.append("end")
    return "\n".join(lines) + "\n"


def _run_cdd(text: str, executable: str) -> str:
    completed = subprocess.run(
        [executable, "--rep"],
        input=text,
        text=True,
        capture_output=True,
        timeout=600,
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
    if kind != "rational" or width != LATTICE_RANK + 1:
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
        raise ArithmeticError("cddlib emitted duplicate oriented rays")
    return tuple(rays)


def _wall_sublattice_index(walls: sp.Matrix) -> int:
    smith = smith_normal_form(walls, domain=sp.ZZ)
    diagonal = tuple(abs(int(smith[index, index])) for index in range(LATTICE_RANK))
    if not all(diagonal):
        raise ArithmeticError("the displayed walls fail to span rank ten")
    return math.prod(diagonal)


def _off_diagonal_histogram(configuration: sp.Matrix) -> dict[str, int]:
    histogram = Counter(
        int(configuration[left, right])
        for left in range(configuration.rows)
        for right in range(left)
    )
    return {str(value): histogram[value] for value in sorted(histogram)}


def _endpoint_type_histogram(verification: Any) -> dict[str, int]:
    histogram: Counter[str] = Counter()
    for audit in verification.edge_audits:
        finite_count = len(audit.finite_endpoints)
        ideal_count = len(audit.ideal_endpoints)
        histogram["F" * finite_count + "I" * ideal_count] += 1
    return dict(sorted(histogram.items()))


def _analyze_exact_cone(
    configuration: sp.Matrix,
    *,
    cdd_executable: str,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    pivot_facets = tuple(int(index) for index in configuration.rref()[1])
    if len(pivot_facets) != LATTICE_RANK:
        raise ArithmeticError("the wall configuration does not span rank ten")
    inequalities = configuration[:, list(pivot_facets)]
    basis_gram = configuration.extract(pivot_facets, pivot_facets)
    if inequalities.rank() != LATTICE_RANK:
        raise ArithmeticError("the cone inequality matrix has wrong rank")

    cdd_text = _cdd_h_representation(inequalities)
    rays = _parse_cdd_rays(_run_cdd(cdd_text, cdd_executable))
    records: list[dict[str, Any]] = []
    kind_histogram: Counter[str] = Counter()
    active_histogram: Counter[tuple[str, int]] = Counter()
    facet_active_rays: list[list[tuple[int, ...]]] = [
        [] for _ in range(configuration.rows)
    ]
    incidence_sets: set[tuple[int, ...]] = set()

    for ray in rays:
        vector = sp.Matrix(ray)
        slacks = inequalities * vector
        if any(slack < 0 for slack in slacks):
            raise ArithmeticError("an exact cdd ray violates a nef inequality")
        active = tuple(index for index, slack in enumerate(slacks) if slack == 0)
        if len(active) < DIMENSION:
            raise ArithmeticError("an extreme projective ray has too few active walls")
        square = sp.Rational((vector.T * basis_gram * vector)[0])
        kind = "finite" if square > 0 else "ideal" if square == 0 else "hyperideal"
        kind_histogram[kind] += 1
        active_histogram[(kind, len(active))] += 1
        if active in incidence_sets:
            raise ArithmeticError("two cone rays have the same full wall incidence")
        incidence_sets.add(active)
        for facet in active:
            facet_active_rays[facet].append(ray)
        records.append(
            {
                "ray": list(ray),
                "active_facets": [index + 1 for index in active],
                "lorentz_square": _fraction_text(square),
                "kind": kind,
            }
        )

    facet_face_ranks = []
    for active_rays in facet_active_rays:
        if not active_rays:
            raise ArithmeticError("a proposed nef wall has no incident cone ray")
        face_rank = int(sp.Matrix(active_rays).rank())
        if face_rank != DIMENSION:
            raise ArithmeticError(
                "a proposed wall does not cut a codimension-one cone face"
            )
        facet_face_ranks.append(face_rank)

    interior = sp.Matrix(
        [sum(ray[column] for ray in rays) for column in range(LATTICE_RANK)]
    )
    interior_slacks = inequalities * interior
    interior_square = sp.Rational((interior.T * basis_gram * interior)[0])
    if not (all(value > 0 for value in interior_slacks) and interior_square > 0):
        raise ArithmeticError("the sum of exact cone rays is not a strict ample witness")
    time_pairings = tuple(
        sp.Rational((sp.Matrix(ray).T * basis_gram * interior)[0]) for ray in rays
    )
    if not all(value > 0 for value in time_pairings):
        raise ArithmeticError("the exact cone rays do not lie in one future time sheet")

    records.sort(key=lambda record: (record["active_facets"], record["ray"]))
    summary = {
        "engine": "cddlib GMP rational double-description",
        "coordinate_basis_facet_positions": [index + 1 for index in pivot_facets],
        "extreme_ray_count": len(rays),
        "finite_vertex_count": kind_histogram["finite"],
        "ideal_vertex_count": kind_histogram["ideal"],
        "hyperideal_vertex_count": kind_histogram["hyperideal"],
        "active_facet_count_histogram": {
            f"{kind}_{count}": value
            for (kind, count), value in sorted(active_histogram.items())
        },
        "all_facet_face_ranks": facet_face_ranks,
        "strict_ample_witness_in_wall_basis": [int(value) for value in interior],
        "minimum_ample_pairing": _fraction_text(min(interior_slacks)),
        "ample_witness_square": _fraction_text(interior_square),
        "all_causal_rays_in_ample_time_sheet": True,
        "minimum_ray_ample_lorentz_pairing": _fraction_text(min(time_pairings)),
        "finite_volume": kind_histogram["hyperideal"] == 0,
        "compact": kind_histogram["hyperideal"] == 0
        and kind_histogram["ideal"] == 0,
        "criterion": (
            "the complete pointed rational projective cone has no spacelike/"
            "hyperideal extreme ray"
        ),
    }
    return summary, records


def _vinberg_payload(configuration: sp.Matrix) -> dict[str, Any]:
    wall_gram = [[-int(value) for value in row] for row in configuration.tolist()]
    verification = verify_vinberg_volume_criterion(
        wall_gram,
        DIMENSION,
        simple_roots_complete=True,
    )
    return {
        "method": (
            "literal Vinberg elliptic-edge endpoint and full affine "
            "zero-neighborhood criterion"
        ),
        "simple_roots_complete": verification.simple_roots_complete,
        "decision": verification.decision,
        "finite_volume": verification.finite_volume,
        "compact": verification.compact,
        "elliptic_rank_8_edge_count": verification.elliptic_edge_count,
        "finite_vertex_count": verification.finite_vertex_count,
        "ideal_vertex_count": verification.ideal_vertex_count,
        "connected_affine_component_audit_count": len(
            verification.parabolic_neighborhood_audits
        ),
        "invalid_full_zero_neighborhood_count": sum(
            not audit.is_rank_dimension_minus_one_parabolic
            for audit in verification.parabolic_neighborhood_audits
        ),
        "endpoint_type_histogram": _endpoint_type_histogram(verification),
        "bad_finite_volume_edge_count": len(verification.finite_volume_bad_edges),
        "bad_compact_edge_count": len(verification.compact_bad_edges),
        "every_elliptic_edge_has_exactly_two_finite_or_ideal_endpoints": not verification.finite_volume_bad_edges,
        "every_full_affine_zero_neighborhood_is_parabolic_rank_8": all(
            audit.is_rank_dimension_minus_one_parabolic
            for audit in verification.parabolic_neighborhood_audits
        ),
    }


def certify(
    raw_path: Path = RAW_DATA,
    *,
    cdd_executable: str = "cddexec_gmp",
) -> tuple[dict[str, Any], bytes]:
    if _sha256(raw_path) != RAW_DATA_SHA256:
        raise ArithmeticError("the pinned Shimada GAP data changed")
    if shutil.which(cdd_executable) is None:
        raise RuntimeError(f"required exact double-description engine missing: {cdd_executable}")
    records = parse_raw_gap_records(raw_path)
    expected_names = {
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
    if set(records) != expected_names:
        raise ArithmeticError("the sixteen Shimada configuration names changed")

    common_gram = sp.Matrix(records["12A"]["GramL10"])
    if not (
        common_gram.rows == common_gram.cols == LATTICE_RANK
        and common_gram == common_gram.T
        and common_gram.det() == -1
        and inertia(_integer_rows(common_gram)) == Inertia(1, 9, 0)
        and all(sp.Matrix(record["GramL10"]) == common_gram for record in records.values())
    ):
        raise ArithmeticError("the raw Enriques lattice is not the pinned L10")

    selected_records: list[dict[str, Any]] = []
    ray_lines = [
        json.dumps(
            {
                "schema": "enriques_finite_aut_nef_chamber_rays_v1",
                "record_format": {
                    "type": "Nikulin-Kondo type",
                    "basis": "one-based independent wall positions",
                    "ray": "primitive integer coordinates in that wall basis",
                    "active_facets": "one-based complete incidence set",
                    "lorentz_square": "exact integer/rational",
                    "kind": "finite, ideal, or hyperideal",
                },
            },
            sort_keys=True,
            separators=(",", ":"),
        )
    ]
    distinguishing_histograms: dict[str, dict[str, int]] = {}

    for specification in TYPE_SPECS:
        type_name = specification["type"]
        record_name = specification["record"]
        raw_record = records[record_name]
        configuration = sp.Matrix(raw_record["configuration"])
        walls = sp.Matrix(raw_record["walls"])
        facet_count = int(specification["facets"])
        if not (
            configuration.rows == configuration.cols == facet_count
            and walls.rows == facet_count
            and walls.cols == LATTICE_RANK
            and configuration == configuration.T
            and all(configuration[index, index] == -2 for index in range(facet_count))
            and set(configuration) <= {-2, 0, 1, 2}
            and configuration == walls * common_gram * walls.T
            and walls.rank() == LATTICE_RANK
            and all(math.gcd(*(abs(int(value)) for value in walls.row(index))) == 1 for index in range(facet_count))
        ):
            raise ArithmeticError(f"raw wall-vector replay failed for type {type_name}")

        normal_gram = -configuration
        expected_inertia = Inertia(DIMENSION, 1, facet_count - LATTICE_RANK)
        if inertia(_integer_rows(normal_gram)) != expected_inertia:
            raise ArithmeticError(f"wrong hyperbolic Gram inertia for type {type_name}")

        pair_histogram = _off_diagonal_histogram(configuration)
        distinguishing_histograms[type_name] = pair_histogram
        cone, ray_records = _analyze_exact_cone(
            configuration,
            cdd_executable=cdd_executable,
        )
        vinberg = _vinberg_payload(configuration)
        if not (
            cone["finite_vertex_count"] == specification["finite_vertices"]
            == vinberg["finite_vertex_count"]
            and cone["ideal_vertex_count"] == specification["ideal_vertices"]
            == vinberg["ideal_vertex_count"]
            and vinberg["elliptic_rank_8_edge_count"] == specification["elliptic_edges"]
            and cone["hyperideal_vertex_count"] == 0
            and cone["finite_volume"] is True
            and cone["compact"] is False
            and vinberg["decision"] == "finite_volume_noncompact"
            and vinberg["bad_finite_volume_edge_count"] == 0
            and vinberg["invalid_full_zero_neighborhood_count"] == 0
        ):
            raise ArithmeticError(f"finite-volume audit changed for type {type_name}")

        for ray in ray_records:
            ray_lines.append(
                json.dumps(
                    {
                        "type": type_name,
                        "shimada_record": record_name,
                        "basis": cone["coordinate_basis_facet_positions"],
                        **ray,
                    },
                    sort_keys=True,
                    separators=(",", ":"),
                )
            )

        selected_records.append(
            {
                "nikulin_kondo_type": type_name,
                "shimada_record": record_name,
                "alternate_shimada_records": specification.get("alternate_records", []),
                "surface_moduli": specification["moduli"],
                "surface_automorphism_group": specification[
                    "surface_automorphism_group"
                ],
                "hyperbolic_dimension": DIMENSION,
                "facet_count": facet_count,
                "intersection_matrix": _integer_rows(configuration),
                "wall_vectors_in_L10_basis": _integer_rows(walls),
                "normalized_positive_normal_gram": _integer_rows(normal_gram),
                "wall_vector_checks": {
                    "all_primitive": True,
                    "all_square_minus_two": True,
                    "span_rank": int(walls.rank()),
                    "wall_sublattice_index_in_L10": _wall_sublattice_index(walls),
                    "configuration_rebuilt_as_W_GramL10_Wt": True,
                },
                "coxeter_pair_counts": {
                    "orthogonal_m2": int(pair_histogram.get("0", 0)),
                    "angle_pi_over_3_m3": int(pair_histogram.get("1", 0)),
                    "parallel_at_infinity": int(pair_histogram.get("2", 0)),
                    "ultraparallel": 0,
                },
                "normal_gram_inertia": {
                    "positive": DIMENSION,
                    "negative": 1,
                    "zero": facet_count - LATTICE_RANK,
                },
                "exact_cone": cone,
                "literal_vinberg": vinberg,
                "classification": "finite_volume_noncompact_H9_Coxeter_polytope",
                "carrier_strength": "complete full-nef Enriques chamber",
            }
        )

    # Shimada proves 20C and 20D isomorphic; in the pinned data their
    # configuration matrices are literally equal in the published ordering.
    if records["20C"]["configuration"] != records["20D"]["configuration"]:
        raise ArithmeticError("the two published Type-VII configurations diverged")

    # No graph-isomorphism black box is needed: the weighted pair-count triples
    # already separate all seven complete Gram configurations.
    signatures = {
        type_name: tuple(sorted(histogram.items()))
        for type_name, histogram in distinguishing_histograms.items()
    }
    if len(set(signatures.values())) != 7:
        raise ArithmeticError("the weighted pair counts no longer separate seven types")

    rays_bytes = gzip.compress(
        ("\n".join(ray_lines) + "\n").encode("utf-8"),
        compresslevel=9,
        mtime=0,
    )
    result = {
        "schema": "enriques_finite_automorphism_nef_chambers_v1",
        "scope": {
            "base_field": "complex numbers",
            "surface_class": "Enriques surfaces with finite automorphism group",
            "polytope_class": (
                "complete full-nef chambers cut out by every smooth rational "
                "curve in Num(Y)"
            ),
            "worldwide_novelty_claim": False,
        },
        "primary_sources": {
            "kondo_1986": {
                "title": "Enriques surfaces with finite automorphism groups",
                "doi": "10.4099/math1924.12.191",
                "role": (
                    "complex seven-type classification; 12/20 complete nodal-curve "
                    "counts; specialized Vinberg criterion and no-hidden-curve proposition"
                ),
            },
            "martin_2019": {
                "title": "Enriques surfaces with finite automorphism group in positive characteristic",
                "doi": "10.14231/AG-2019-027",
                "arxiv": "1703.08419",
                "role": "modern iff classification statement and moduli table",
            },
            "shimada_configuration_note": {
                "title": "A note on configurations of (-2)-vectors on Enriques surfaces",
                "role": "exact dictionary from records 12A--20F to types I--VII",
            },
            "shimada_raw_gap_data": {
                "path": str(raw_path.relative_to(ROOT)),
                "sha256": RAW_DATA_SHA256,
                "record_count": len(records),
                "role": "exact L10 Gram, wall vectors, and intersection matrices",
            },
        },
        "lattice": {
            "name": "L10 = U + E8(-1) = Num(Y)",
            "gram": _integer_rows(common_gram),
            "rank": LATTICE_RANK,
            "determinant": int(common_gram.det()),
            "even": all(common_gram[index, index] % 2 == 0 for index in range(LATTICE_RANK)),
            "unimodular": abs(int(common_gram.det())) == 1,
            "inertia": {"positive": 1, "negative": 9, "zero": 0},
        },
        "surface_to_polytope_map": {
            "notation": "Phi_Enr",
            "domain": (
                "numerically marked complex Enriques surfaces with finite "
                "automorphism group"
            ),
            "map": "forget the surface and retain the congruence class of its full nef chamber",
            "image_cardinality": 7,
            "surjectivity_reason": (
                "Kondo-Martin prove existence of every one of the seven complete "
                "nodal configurations; the exact weighted Gram signatures below "
                "show that their chamber images do not collapse"
            ),
            "universality_scope": (
                "exact for this surface-free finite list of seven normalized Gram "
                "types; not a claim that every H9 Coxeter polytope has an Enriques carrier"
            ),
        },
        "types": selected_records,
        "pairwise_distinction": {
            "method": (
                "facet count plus the exact histogram of off-diagonal root "
                "intersections 0,1,2"
            ),
            "weighted_pair_histograms": distinguishing_histograms,
            "all_seven_distinct": True,
            "20C_and_20D_configuration_matrices_literally_equal": True,
        },
        "global_checks": {
            "all_seven_exact_wall_lists_rebuilt_from_raw_vectors": True,
            "all_seven_complete_nef_chambers": True,
            "all_seven_exact_cdd_cones_have_no_hyperideal_ray": True,
            "all_seven_literal_Vinberg_audits_pass": True,
            "all_seven_finite_volume": True,
            "all_seven_noncompact": True,
            "facet_count_multiset": [12, 12, 20, 20, 20, 20, 20],
            "no_full_nef_finite_Aut_Enriques_H9_14_member": True,
        },
        "conclusion": (
            "Phi_Enr has exactly seven full-nef H9 Coxeter-polytope images: "
            "two 12-facet and five 20-facet types. Every image has finite volume "
            "and at least one cusp, hence none is compact."
        ),
    }
    return result, rays_bytes


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", type=Path, default=RAW_DATA)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--rays-output", type=Path, default=RAY_OUTPUT)
    parser.add_argument("--cdd", default="cddexec_gmp")
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()

    result, rays_bytes = certify(arguments.raw, cdd_executable=arguments.cdd)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file() or arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archived summary differs: {arguments.output}")
        if not arguments.rays_output.is_file() or arguments.rays_output.read_bytes() != rays_bytes:
            raise SystemExit(f"archived ray census differs: {arguments.rays_output}")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.rays_output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
        arguments.rays_output.write_bytes(rays_bytes)

    print(
        json.dumps(
            {
                "summary": str(arguments.output),
                "summary_sha256": hashlib.sha256(rendered.encode("utf-8")).hexdigest(),
                "rays": str(arguments.rays_output),
                "rays_sha256": hashlib.sha256(rays_bytes).hexdigest(),
                "type_counts": [
                    {
                        "type": record["nikulin_kondo_type"],
                        "facets": record["facet_count"],
                        "finite_vertices": record["exact_cone"]["finite_vertex_count"],
                        "cusps": record["exact_cone"]["ideal_vertex_count"],
                        "edges": record["literal_vinberg"]["elliptic_rank_8_edge_count"],
                    }
                    for record in result["types"]
                ],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
