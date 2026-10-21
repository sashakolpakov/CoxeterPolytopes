#!/usr/bin/env python3
"""Certify uniform full-Picard K3 carriers for three H10 chambers.

The three remaining machine-Gram rows have carrier lattices

    S_a = <2a> + <-2>^10,  a=1,2,3.

They lie outside Nikulin's easy stable range.  This verifier instead builds
an explicit discriminant anti-isometry to

    T_a = <2>^2 + <-2a> + <-2>^8,

forms the graph-glued even unimodular lattice of signature (3,19), and
checks the resulting primitive embedding and complement matrices exactly.
It also replays every oriented wall/reflection and both finite-volume lanes.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
from typing import Any

import sympy as sp

from scripts.certify_three_h11_full_picard_k3_carriers import (
    _carrier_lattice,
    _matrix_payload,
    _payload_sha256,
    _smith,
    _wall_system,
)
from scripts.search_ss_fermat_quartic_char3_folds import parse_source
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
D5 = ROOT / "data/enumerations/shimada_disc16_d5_h10_30_lineage.json"
D5_CONE = (
    ROOT / "data/enumerations/shimada_disc16_d5_h10_30_exact_cone_certificate.json"
)
D5_VINBERG = (
    ROOT / "data/enumerations/shimada_disc16_d5_h10_30_vinberg_face_certificate.json"
)
SS2 = ROOT / "data/enumerations/borcherds_char2_2d_h10_47_certificate.json"
SS2_LINEAGE = ROOT / "data/enumerations/borcherds_char2_2d_h10_47_lineage.json"
F3 = ROOT / "data/enumerations/ss_fermat_quartic_char3_survivor_fixed_lattices.json"
F3_CONE = (
    ROOT / "data/enumerations/ss_fermat_quartic_char3_c2_h10_exact_cone_certificate.json"
)
F3_VINBERG = (
    ROOT / "data/enumerations/ss_fermat_quartic_char3_c2_h10_vinberg_certificate.json"
)
FERMAT_SOURCE = (
    ROOT
    / "data/external/kondo_shimada_fermat_quartic_char3/compdataAutFQChar3.txt"
)
OUTPUT = (
    ROOT / "data/enumerations/three_h10_uniform_full_picard_k3_carriers.json"
)

SCHEMA = "coxeter.three_h10.uniform_full_picard_complex_k3_carriers.v1"
EXPECTED_SHA256 = {
    "d5": "1bcdb45dfbf15b21f7b22f72c475f5946f3c06cba6732ac4995bf7b8e434a7e2",
    "d5_cone": "1160d8e5e4425a20a8a12bc9df6c6dd7fbd5031c5ebe8965b8d4357c82dad083",
    "d5_vinberg": "4e0fc3d33da785537527df982d422312c0596ca8e93ced73f812934133255f8b",
    "ss2": "bb40f8428485168ed49b5cc26c53df4a3d8ce27a1a1abf4632c42ac3f43febfb",
    "ss2_lineage": "165fef61e196bbd08d0a42b3c4b4db861879ceee3154be78aa4d6da38aff56d4",
    "f3": "431ddb789fcb30e843efe870a67f9ce276ccd75c4fb65ea4f01f51afdc836124",
    "f3_cone": "defc07c5a2228a2582cadfb9729e1122d61d7fac81f28e8218397162e0050045",
    "f3_vinberg": "57d025394f0f1b78e7d65f6380d9bd18f0982e9ff4daf2914785d98071af33bb",
    "fermat_source": "48d82456889f9fc378e02a6ac7a83ee16f24914a4d636b837bfeeac489a07baf",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _carrier_gram(a: int) -> sp.Matrix:
    return sp.diag(2 * a, *([-2] * 10))


def _complement_gram(a: int) -> sp.Matrix:
    return sp.diag(2, 2, -2 * a, *([-2] * 8))


def _four_sign_flip() -> sp.Matrix:
    return sp.ones(4) - sp.eye(4)


def _two_primary_anti_isometry() -> sp.Matrix:
    return sp.diag(sp.eye(2), _four_sign_flip(), _four_sign_flip())


def _mod_two(value: sp.Rational) -> sp.Rational:
    value = sp.Rational(value)
    return value - 2 * sp.floor(value / 2)


def _explicit_k3_gluing(a: int) -> dict[str, Any]:
    carrier = _carrier_gram(a)
    complement = _complement_gram(a)
    block = _two_primary_anti_isometry()

    anti_isometry_checks = []
    for cyclic in range(2 * a):
        for mask in range(1 << 10):
            source_bits = sp.Matrix([(mask >> index) & 1 for index in range(10)])
            target_bits = (block * source_bits).applyfunc(lambda value: int(value) % 2)
            source_q = sp.Rational(cyclic * cyclic, 2 * a) - sp.Rational(
                sum(source_bits), 2
            )
            target_q = -sp.Rational(cyclic * cyclic, 2 * a) + sp.Rational(
                int(target_bits[0]) + int(target_bits[1]), 2
            ) - sp.Rational(sum(target_bits[2:, 0]), 2)
            anti_isometry_checks.append(_mod_two(source_q + target_q) == 0)

    # Columns 0..10 are graph-glue generators; columns 11..21 are the
    # standard T_a basis.  These columns also recover every standard S_a
    # basis vector, hence form a basis of the full graph-glued overlattice.
    glue_basis = sp.zeros(22, 22)
    glue_basis[0, 0] = sp.Rational(1, 2 * a)
    glue_basis[11 + 2, 0] = sp.Rational(1, 2 * a)
    target_order_two_rows = [0, 1] + list(range(3, 11))
    for source_index in range(10):
        column = 1 + source_index
        glue_basis[1 + source_index, column] = sp.Rational(1, 2)
        for target_index, target_row in enumerate(target_order_two_rows):
            if int(block[target_index, source_index]) % 2:
                glue_basis[11 + target_row, column] = sp.Rational(1, 2)
    for target_index in range(11):
        glue_basis[11 + target_index, 11 + target_index] = 1

    direct_sum = sp.diag(carrier, complement)
    ambient = sp.simplify(glue_basis.T * direct_sum * glue_basis)
    standard_to_glued = sp.simplify(glue_basis.inv())
    embedding = standard_to_glued[:, :11]
    orthogonal_complement = standard_to_glued[:, 11:]
    graph_order = (2 * a) * 2**10
    combined = embedding.row_join(orthogonal_complement)
    checks = {
        "two_primary_matrix_invertible_mod_2": int(block.det()) % 2 == 1,
        "all_discriminant_classes_checked": len(anti_isometry_checks)
        == graph_order,
        "finite_quadratic_forms_are_opposites": all(anti_isometry_checks),
        "glue_basis_has_expected_index": sp.Rational(abs(glue_basis.det()))
        == sp.Rational(1, graph_order),
        "ambient_integral": all(value.is_Integer for value in ambient),
        "ambient_even": all(int(ambient[index, index]) % 2 == 0 for index in range(22)),
        "ambient_unimodular": abs(int(ambient.det())) == 1,
        "ambient_signature_3_19": exact_inertia(ambient) == (3, 19, 0),
        "carrier_embedding_integral": all(value.is_Integer for value in embedding),
        "carrier_embedding_primitive": _smith(embedding) == [1] * 11,
        "carrier_embedding_isometric": embedding.T * ambient * embedding
        == carrier,
        "complement_embedding_integral": all(
            value.is_Integer for value in orthogonal_complement
        ),
        "complement_embedding_primitive": _smith(orthogonal_complement)
        == [1] * 11,
        "complement_embedding_isometric": (
            orthogonal_complement.T * ambient * orthogonal_complement
            == complement
        ),
        "carrier_and_complement_orthogonal": (
            embedding.T * ambient * orthogonal_complement == sp.zeros(11)
        ),
        "combined_index_is_discriminant_order": abs(int(combined.det()))
        == graph_order,
    }
    return {
        "criterion": "explicit full-discriminant graph gluing",
        "ambient_identification": (
            "the unique even unimodular lattice of signature (3,19), hence Lambda_K3"
        ),
        "carrier": f"S_{a}=< {2*a} >+<-2>^10",
        "orthogonal_complement": f"T_{a}=<2>^2+<-{2*a}>+<-2>^8",
        "carrier_gram": _matrix_payload(carrier),
        "orthogonal_complement_gram": _matrix_payload(complement),
        "two_primary_anti_isometry_matrix_mod_2": _matrix_payload(block),
        "cyclic_generator_map": (
            f"1 in Z/{2*a} maps to 1 in Z/{2*a}; quadratic values "
            f"1/{2*a} and -1/{2*a}"
        ),
        "graph_order": graph_order,
        "graph_glue_basis_columns_in_S_plus_T_Q_space": _matrix_payload(glue_basis),
        "glued_even_unimodular_gram": _matrix_payload(ambient),
        "carrier_embedding_columns_in_glued_basis": _matrix_payload(embedding),
        "complement_embedding_columns_in_glued_basis": _matrix_payload(
            orthogonal_complement
        ),
        "combined_carrier_complement_index": abs(int(combined.det())),
        "checked_discriminant_class_count": len(anti_isometry_checks),
        "checks": checks,
    }


def _content(vector: sp.Matrix) -> int:
    return math.gcd(*(abs(int(value)) for value in vector))


def _ss2_exact_cone(gram: sp.Matrix, walls: dict[str, Any]) -> dict[str, Any]:
    oriented = [sp.Matrix(row["oriented_coordinates"]) for row in walls["wall_records"]]
    inequality = sp.Matrix.vstack(*[(root.T * gram) for root in oriented])
    lines = ["H-representation", "begin", f"{len(oriented)} 12 rational"]
    lines.extend(
        "0 " + " ".join(str(int(value)) for value in inequality.row(row))
        for row in range(inequality.rows)
    )
    lines.append("end")
    executable = shutil.which("cddexec_gmp") or "cddexec_gmp"
    completed = subprocess.run(
        [executable, "--rep"],
        input="\n".join(lines) + "\n",
        text=True,
        capture_output=True,
        timeout=600,
        check=False,
    )
    if completed.returncode:
        raise RuntimeError(completed.stdout[-4000:] + completed.stderr[-4000:])
    output_lines = completed.stdout.splitlines()
    begin = output_lines.index("begin")
    count_text, width_text, kind = output_lines[begin + 1].split()
    count, width = int(count_text), int(width_text)
    if width != 12 or kind != "rational":
        raise ArithmeticError("cddlib emitted the wrong SS2 generator format")
    ray_lines = output_lines[begin + 2 : begin + 2 + count]
    rays: list[tuple[int, ...]] = []
    for line in ray_lines:
        tokens = line.split()
        if len(tokens) != width or Fraction(tokens[0]) != 0:
            raise ArithmeticError("SS2 homogeneous cone emitted a vertex row")
        fractions = [Fraction(token) for token in tokens[1:]]
        denominator = math.lcm(*(value.denominator for value in fractions))
        values = [value.numerator * (denominator // value.denominator) for value in fractions]
        divisor = math.gcd(*(abs(value) for value in values))
        primitive = tuple(value // divisor for value in values)
        rays.append(primitive)
    if len(set(rays)) != len(rays):
        raise ArithmeticError("cddlib emitted duplicate SS2 projective rays")

    kinds: Counter[str] = Counter()
    active: Counter[str] = Counter()
    for ray in rays:
        vector = sp.Matrix(ray)
        slacks = inequality * vector
        if any(value < 0 for value in slacks):
            raise ArithmeticError("SS2 exact ray violates an inequality")
        active_count = sum(value == 0 for value in slacks)
        if active_count < 10:
            raise ArithmeticError("SS2 exact ray has too few active facets")
        square = int((vector.T * gram * vector)[0])
        ray_kind = "finite" if square > 0 else "ideal" if square == 0 else "hyperideal"
        kinds[ray_kind] += 1
        active[f"{ray_kind}_{active_count}"] += 1
    return {
        "method": "fresh cddlib GMP rational double description with exact replay",
        "extreme_ray_count": len(rays),
        "finite_vertex_count": kinds["finite"],
        "ideal_vertex_count": kinds["ideal"],
        "hyperideal_vertex_count": kinds["hyperideal"],
        "active_facet_count_histogram": dict(sorted(active.items())),
        "primitive_ray_matrix_sha256": _payload_sha256(rays),
        "checks": {
            "inequality_rank_11": inequality.rank() == 11,
            "all_projective_rays_distinct": len(set(rays)) == len(rays),
            "no_hyperideal_extreme_ray": kinds["hyperideal"] == 0,
        },
    }


def _volume_payloads(
    d5_cone: dict[str, Any],
    d5_vinberg: dict[str, Any],
    ss2: dict[str, Any],
    ss2_cone: dict[str, Any],
    f3_cone: dict[str, Any],
    f3_vinberg: dict[str, Any],
) -> dict[str, dict[str, Any]]:
    return {
        "D16-D5-H10-30": {
            "literal_Vinberg": {
                "classification": d5_vinberg["decision"]["classification"],
                "elliptic_rank9_edges": d5_vinberg["literal_endpoint_audit"]["edge_count"],
                "finite_vertices": d5_vinberg["literal_endpoint_audit"]["finite_vertex_count"],
                "cusps": d5_vinberg["literal_endpoint_audit"]["ideal_vertex_count"],
                "bad_edges": d5_vinberg["literal_endpoint_audit"]["bad_edge_count"],
                "bad_parabolic_neighborhoods": d5_vinberg["literal_endpoint_audit"]["bad_parabolic_neighborhood_count"],
            },
            "exact_GMP_cone": {
                "classification": d5_cone["decision"]["classification"],
                "extreme_rays": d5_cone["exact_double_description"]["extreme_ray_count"],
                "finite_vertices": d5_cone["exact_double_description"]["finite_vertex_count"],
                "cusps": d5_cone["exact_double_description"]["ideal_vertex_count"],
                "hyperideal_rays": d5_cone["exact_double_description"]["hyperideal_vertex_count"],
            },
        },
        "SS2-2D-H10-47": {
            "literal_Vinberg": {
                "classification": ss2["vinberg"]["decision"],
                "elliptic_rank9_edges": ss2["vinberg"]["elliptic_edge_count"],
                "finite_vertices": ss2["vinberg"]["finite_vertex_count"],
                "cusps": ss2["vinberg"]["ideal_vertex_count"],
                "bad_edges": ss2["vinberg"]["bad_edge_count"],
                "bad_parabolic_neighborhoods": ss2["vinberg"]["bad_parabolic_neighborhood_count"],
            },
            "exact_GMP_cone": {
                "classification": "finite_volume_noncompact",
                "extreme_rays": ss2_cone["extreme_ray_count"],
                "finite_vertices": ss2_cone["finite_vertex_count"],
                "cusps": ss2_cone["ideal_vertex_count"],
                "hyperideal_rays": ss2_cone["hyperideal_vertex_count"],
                "fresh_replay": ss2_cone,
            },
        },
        "F3-C2-H10-196": {
            "literal_Vinberg": {
                "classification": f3_vinberg["vinberg_decision"],
                "elliptic_rank9_edges": f3_vinberg["elliptic_edge_count"],
                "finite_vertices": f3_vinberg["finite_vertex_count"],
                "cusps": f3_vinberg["ideal_vertex_count"],
                "bad_edges": f3_vinberg["bad_edge_count"],
                "bad_parabolic_neighborhoods": f3_vinberg["invalid_full_zero_neighborhood_count"],
            },
            "exact_GMP_cone": {
                "classification": "finite_volume_noncompact",
                "extreme_rays": f3_cone["ray_count"],
                "finite_vertices": f3_cone["finite_ray_count"],
                "cusps": f3_cone["ideal_ray_count"],
                "hyperideal_rays": f3_cone["hyperideal_ray_count"],
            },
        },
    }


def build_certificate() -> dict[str, Any]:
    paths = {
        "d5": D5,
        "d5_cone": D5_CONE,
        "d5_vinberg": D5_VINBERG,
        "ss2": SS2,
        "ss2_lineage": SS2_LINEAGE,
        "f3": F3,
        "f3_cone": F3_CONE,
        "f3_vinberg": F3_VINBERG,
        "fermat_source": FERMAT_SOURCE,
    }
    for key, path in paths.items():
        if _sha256(path) != EXPECTED_SHA256[key]:
            raise ArithmeticError(f"pinned {key} archive changed")

    d5 = json.loads(D5.read_text(encoding="utf-8"))
    d5_cone = json.loads(D5_CONE.read_text(encoding="utf-8"))
    d5_vinberg = json.loads(D5_VINBERG.read_text(encoding="utf-8"))
    ss2 = json.loads(SS2.read_text(encoding="utf-8"))
    ss2_lineage = json.loads(SS2_LINEAGE.read_text(encoding="utf-8"))
    f3 = json.loads(F3.read_text(encoding="utf-8"))
    f3_cone = json.loads(F3_CONE.read_text(encoding="utf-8"))
    f3_vinberg = json.loads(F3_VINBERG.read_text(encoding="utf-8"))

    d5_records = d5["primitive_mirrors"]["records"]
    d5_roots = [row["primitive_I_1_10_coordinates"] for row in d5_records]
    d5_source_signs = [
        1 if sp.Rational(row["presented_equals_scale_times_primitive"]) > 0 else -1
        for row in d5_records
    ]
    d5_oriented = [
        sign * sp.Matrix(root) for sign, root in zip(d5_source_signs, d5_roots)
    ]
    d5_sum = sum(d5_oriented, sp.zeros(11, 1))
    d5_interior = d5_sum / _content(d5_sum)
    d5_walls = _wall_system(
        _carrier_gram(1),
        d5_roots,
        d5_interior,
        {"-2,2": 18, "-4,2": 12},
    )

    ss2_change = sp.Matrix(ss2_lineage["lattice_isometry"]["basis_change_columns"])
    ss2_roots_columns = [
        sp.simplify(ss2_change.inv() * sp.Matrix(row["primitive_fixed_coordinates"]))
        for row in ss2["restricted_chamber"]["primitive_facet_lattice_records"]
    ]
    if not all(value.is_Integer for root in ss2_roots_columns for value in root):
        raise ArithmeticError("SS2 roots did not transport integrally")
    ss2_roots = [[int(value) for value in root] for root in ss2_roots_columns]
    ss2_interior = sum(ss2_roots_columns, sp.zeros(11, 1))
    ss2_walls = _wall_system(
        _carrier_gram(2),
        ss2_roots,
        ss2_interior,
        {"-2,2": 5, "-4,2": 16, "-4,4": 26},
    )
    ss2_fresh_cone = _ss2_exact_cone(_carrier_gram(2), ss2_walls)

    f3_case = next(row for row in f3["cases"] if row["id"] == "char3-C2-H10")
    f3_fixed_basis = sp.Matrix(f3_case["basis_rows_in_published_NS_basis"])
    f3_isometry = sp.Matrix(f3["c2_h10_identification"]["explicit_isometry_matrix"])
    published_ample = parse_source(FERMAT_SOURCE)["ample"].T
    fixed_solution = sp.linsolve((f3_fixed_basis.T, published_ample))
    f3_fixed_interior = sp.Matrix(next(iter(fixed_solution)))
    f3_interior = sp.simplify(f3_isometry.inv() * f3_fixed_interior)
    f3_roots = [
        row["mcleod_lattice_coordinates"]
        for row in f3["c2_h10_196_mirrors"]["mirrors"]
    ]
    f3_walls = _wall_system(
        _carrier_gram(3),
        f3_roots,
        f3_interior,
        {"-2,2": 80, "-4,2": 36, "-12,6": 80},
    )

    volume = _volume_payloads(
        d5_cone, d5_vinberg, ss2, ss2_fresh_cone, f3_cone, f3_vinberg
    )
    case_specs = [
        ("D16-D5-H10-30", 1, d5_walls, d5),
        ("SS2-2D-H10-47", 2, ss2_walls, ss2),
        ("F3-C2-H10-196", 3, f3_walls, f3),
    ]
    chambers = []
    embeddings = {}
    for chamber_id, a, walls, upstream in case_specs:
        embedding = _explicit_k3_gluing(a)
        embeddings[chamber_id] = embedding
        carrier = _carrier_lattice(
            f"<{2*a}>+<-2>^10",
            _carrier_gram(a),
            embedding,
        )
        chamber = {
            "id": chamber_id,
            "dimension": 10,
            "facet_count": len(walls["wall_records"]),
            "carrier": carrier,
            "certified_reflection_stable_overclosures": [
                {
                    "index": 1,
                    "carrier": "original even wall lattice is already reflection-stable",
                }
            ],
            "oriented_wall_system": walls,
            "finite_volume": volume[chamber_id],
            "status": "new nine-dimensional very-general full-Picard complex-K3 carrier",
        }
        if chamber_id == "D16-D5-H10-30":
            chamber["interior_provenance"] = {
                "formula": "primitive reduction of the signed sum of all 30 source-oriented roots",
                "source_orientation_signs_match": d5_walls[
                    "orientation_sign_histogram"
                ]
                == {
                    "-1": d5_source_signs.count(-1),
                    "1": d5_source_signs.count(1),
                },
            }
        elif chamber_id == "SS2-2D-H10-47":
            chamber["interior_provenance"] = {
                "formula": "sum of all 47 transported source-oriented primitive roots"
            }
        else:
            chamber["interior_provenance"] = {
                "formula": "exact transport of the published Fermat ample class",
                "fixed_coordinates": [int(value) for value in f3_fixed_interior],
                "published_ample_recovered": (
                    f3_fixed_basis.T * f3_isometry * f3_interior
                    == published_ample
                ),
            }
        chambers.append(chamber)

    checks = {
        "all_upstream_lattice_certificates_green": (
            d5["all_checks_pass"]
            and ss2["all_checks_pass"]
            and ss2_lineage["all_checks_pass"]
            and f3["all_checks_pass"]
        ),
        "all_three_explicit_K3_gluings_green": all(
            all(embedding["checks"].values()) for embedding in embeddings.values()
        ),
        "all_273_walls_replayed": sum(
            len(chamber["oriented_wall_system"]["wall_records"])
            for chamber in chambers
        )
        == 273,
        "all_wall_sets_generate_full_Picard": all(
            chamber["oriented_wall_system"]["wall_span"]["rank"] == 11
            and chamber["oriented_wall_system"]["wall_span"][
                "index_in_carrier_lattice"
            ]
            == 1
            for chamber in chambers
        ),
        "all_reflections_integral_on_full_Picard": all(
            all(record["checks"].values())
            for chamber in chambers
            for record in chamber["oriented_wall_system"]["wall_records"]
        ),
        "all_literal_Vinberg_gates_green": all(
            chamber["finite_volume"]["literal_Vinberg"]["classification"]
            == "finite_volume_noncompact"
            and chamber["finite_volume"]["literal_Vinberg"]["bad_edges"] == 0
            and chamber["finite_volume"]["literal_Vinberg"][
                "bad_parabolic_neighborhoods"
            ]
            == 0
            for chamber in chambers
        ),
        "all_exact_GMP_gates_green": all(
            chamber["finite_volume"]["exact_GMP_cone"]["classification"]
            == "finite_volume_noncompact"
            and chamber["finite_volume"]["exact_GMP_cone"]["hyperideal_rays"]
            == 0
            for chamber in chambers
        ),
        "both_volume_lanes_agree_for_all_three": all(
            chamber["finite_volume"]["literal_Vinberg"]["finite_vertices"]
            == chamber["finite_volume"]["exact_GMP_cone"]["finite_vertices"]
            and chamber["finite_volume"]["literal_Vinberg"]["cusps"]
            == chamber["finite_volume"]["exact_GMP_cone"]["cusps"]
            for chamber in chambers
        ),
        "all_interiors_have_exact_provenance": all(
            chamber["interior_provenance"].get(
                "source_orientation_signs_match",
                chamber["interior_provenance"].get("published_ample_recovered", True),
            )
            for chamber in chambers
        ),
    }
    return {
        "schema": SCHEMA,
        "status_date": "2026-08-07",
        "sources": {
            key: {"path": str(path.relative_to(ROOT)), "sha256": EXPECTED_SHA256[key]}
            for key, path in paths.items()
        },
        "uniform_theorem": (
            "For a=1,2,3 the lattices S_a=<2a>+<-2>^10 embed primitively "
            "in Lambda_K3 by the displayed full-discriminant graph gluing. "
            "The D16-D5 H10/30, SS2-2D H10/47, and F3-C2 H10/196 "
            "finite-volume noncompact chambers therefore have "
            "nine-dimensional very-general complex-K3 families with "
            "NS(X)=S_a carrying every wall and reflection."
        ),
        "uniform_gluing_lemma": {
            "carrier_family": "S_a=<2a>+<-2>^10",
            "complement_family": "T_a=<2>^2+<-2a>+<-2>^8",
            "two_primary_identity": (
                "q_plus^10 is isometric to q_plus^2 + q_minus^8 via two "
                "copies of J_4-I_4, leaving two coordinates fixed"
            ),
            "embeddings": embeddings,
        },
        "chambers": chambers,
        "scope": {
            "positive_gain": "three known cofinite H10 chambers gain very-general full-Picard complex-K3 carriers",
            "machine_ingested_residue": (
                "these three exhaust only the presently audited repository "
                "queue of certified finite-volume H>=10 chambers with "
                "machine-ingested full Gram/wall data and no prior full-Picard carrier"
            ),
            "not_all_ledger_rows": True,
            "not_world_literature_exhaustion": True,
            "not_new_abstract_polytopes": True,
            "not_full_nef": "mixed square-minus-four and square-minus-twelve walls occur",
            "not_surface_automorphisms": "Picard reflections are not asserted to extend to K3 automorphisms",
            "worldwide_priority_claim": "none",
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate()
    if not payload["all_checks_pass"]:
        raise ArithmeticError("three-H10 full-Picard carrier certificate failed")
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.exists() or args.output.read_text(encoding="utf-8") != rendered:
            raise ArithmeticError("archived three-H10 carrier changed")
        print(f"verified {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(f"wrote {args.output}")
    print(f"payload_sha256={_payload_sha256(payload)}")


if __name__ == "__main__":
    main()
