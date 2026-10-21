#!/usr/bin/env python3
"""Certify three characteristic-p H11 chambers on complex full-Picard K3s.

The inputs are already complete cofinite Coxeter chambers.  This verifier
does the previously missing characteristic-zero carrier step:

* F3-C2-H11/206: use the full primitive K3 embedding criterion, exhibited
  here by an explicit embedding of U+E8(-2)+<-2>^2;
* F3-V4-H11/328: first pass to the certified minimal even reflection-stable
  index-four closure, then apply Nikulin's primitive-embedding stable range;
* SS2-3B-H11/22: apply the same stable range directly to U+A2(-1)^5.

All oriented primitive wall vectors are replayed, every reflection is checked
on the full carrier lattice, and the wall-span index is computed exactly.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
from itertools import product
import json
import math
from pathlib import Path
from typing import Any

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form

from scripts.certify_jacobian_kummer_h10_class0007_fixed_lattice import (
    _discriminant_basis,
)
from scripts.search_ss_fermat_quartic_char3_folds import parse_source
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
FERMAT_SOURCE = (
    ROOT
    / "data/external/kondo_shimada_fermat_quartic_char3/compdataAutFQChar3.txt"
)
F3_FIXED = (
    ROOT
    / "data/enumerations/ss_fermat_quartic_char3_survivor_fixed_lattices.json"
)
F3_C2 = ROOT / "data/enumerations/ss_fermat_quartic_char3_c2_h11_lineage.json"
F3_C2_VINBERG = (
    ROOT
    / "data/enumerations/ss_fermat_quartic_char3_c2_h11_vinberg_certificate.json"
)
F3_V4 = ROOT / "data/enumerations/ss_fermat_quartic_char3_v4_h11_lineage.json"
F3_V4_CONE = (
    ROOT
    / "data/enumerations/ss_fermat_quartic_char3_v4_h11_exact_cone_certificate.json"
)
F3_MIRRORS = (
    ROOT / "data/enumerations/ss_fermat_quartic_char3_h11_mirror_lattices.json"
)
SS2_3B = ROOT / "data/enumerations/borcherds_char2_3b_h11_22_certificate.json"
SS2_3B_LINEAGE = (
    ROOT / "data/enumerations/borcherds_char2_3b_h11_22_lineage.json"
)
OUTPUT = ROOT / "data/enumerations/three_h11_full_picard_k3_carriers.json"

SCHEMA = "coxeter.three_h11.full_picard_complex_k3_carriers.v1"
EXPECTED_SHA256 = {
    "fermat_source": "48d82456889f9fc378e02a6ac7a83ee16f24914a4d636b837bfeeac489a07baf",
    "f3_fixed": "431ddb789fcb30e843efe870a67f9ce276ccd75c4fb65ea4f01f51afdc836124",
    "f3_c2": "e1b20c853d3f2f254900c0dd46cdf187093292fb4356e1d6e9de23ca839faf89",
    "f3_c2_vinberg": "22f2104680afc9f20f9a2a6fd38f7486037f9890b57cd881c1fda6c96027c9bb",
    "f3_v4": "b796d3c368a25867ba7e7c7f5e6d64f7245b4731faadd90d324d72358a564f28",
    "f3_v4_cone": "283fa156fb466224b2548b20f8ccf6eafdac591faa56d6602f7f8d8a67f12e8b",
    "f3_mirrors": "28abf2e31b1af3c0ff101509303ccccb7aff94f81c68f873a1dcab1161a9ee07",
    "ss2_3b": "cf1766077eddc74cfd159633c12fb4828fb8b957d08e2edb881dfca89bf7e044",
    "ss2_3b_lineage": "01525acf079e9a8a76e11593f84bba55fe423af46cb0ff0faf1bf2ee26d7d0bd",
}

C2_INTERIOR = sp.Matrix([-4, -3, 1, 3, 4, 6, 5, 3, 3, 1, -1, 1])
V4_INTERIOR = sp.Matrix([0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0])


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _payload_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _matrix_payload(matrix: sp.Matrix) -> list[list[int | str]]:
    return [
        [
            int(value) if sp.Rational(value).q == 1 else str(sp.Rational(value))
            for value in matrix.row(row)
        ]
        for row in range(matrix.rows)
    ]


def _smith(matrix: sp.Matrix) -> list[int]:
    diagonal = smith_normal_form(matrix, domain=ZZ)
    return [
        abs(int(diagonal[index, index]))
        for index in range(min(diagonal.shape))
        if diagonal[index, index]
    ]


def _content(vector: sp.Matrix) -> int:
    return math.gcd(*(abs(int(value)) for value in vector))


def _divisibility(gram: sp.Matrix, vector: sp.Matrix) -> int:
    return _content(gram * vector)


def _reflection(gram: sp.Matrix, root: sp.Matrix) -> sp.Matrix:
    square = int((root.T * gram * root)[0])
    return sp.simplify(
        sp.eye(gram.rows)
        - sp.Rational(2, square) * root * (root.T * gram)
    )


def _discriminant_data(gram: sp.Matrix) -> dict[str, Any]:
    orders, _, form = _discriminant_basis(gram)
    histogram: Counter[str] = Counter()
    for coefficients in product(*(range(order) for order in orders)):
        vector = sp.Matrix(coefficients)
        value = sp.Rational((vector.T * form * vector)[0])
        value -= 2 * sp.floor(value / 2)
        histogram[str(value)] += 1
    return {
        "generator_orders": orders,
        "length": len(orders),
        "reduced_finite_quadratic_Gram_mod_2Z": _matrix_payload(form),
        "q_value_histogram_mod_2Z": dict(sorted(histogram.items())),
    }


def _minimal_even_scaled_dual(gram: sp.Matrix) -> dict[str, Any]:
    for scale in range(1, 257):
        dual = sp.simplify(scale * gram.inv())
        if all(value.is_Integer for value in dual) and all(
            int(dual[index, index]) % 2 == 0 for index in range(gram.rows)
        ):
            smith = _smith(dual)
            return {
                "minimal_scale": scale,
                "gram": _matrix_payload(dual),
                "determinant_absolute": abs(int(dual.det())),
                "smith_invariants": smith,
                "discriminant_length": sum(value > 1 for value in smith),
            }
    raise ArithmeticError("minimal even scaled dual exceeded search bound")


def _wall_system(
    gram: sp.Matrix,
    raw_roots: list[list[int]],
    interior: sp.Matrix,
    archived_types: dict[str, int],
) -> dict[str, Any]:
    rows: list[list[int]] = []
    records: list[dict[str, Any]] = []
    types: Counter[tuple[int, int]] = Counter()
    signs: Counter[int] = Counter()
    for facet, raw in enumerate(raw_roots, start=1):
        root = sp.Matrix(raw)
        raw_pairing = int((interior.T * gram * root)[0])
        if raw_pairing == 0:
            raise ArithmeticError("interior lies on a facet")
        sign = 1 if raw_pairing > 0 else -1
        oriented = sign * root
        square = int((oriented.T * gram * oriented)[0])
        divisibility = _divisibility(gram, oriented)
        pairing = int((interior.T * gram * oriented)[0])
        reflection = _reflection(gram, oriented)
        checks = {
            "primitive": _content(oriented) == 1,
            "negative_square": square < 0,
            "strict_positive_orientation": pairing > 0,
            "reflection_integral": all(value.is_Integer for value in reflection),
            "reflection_involutive": reflection**2 == sp.eye(gram.rows),
            "reflection_Gram_isometry": reflection.T * gram * reflection == gram,
            "reflection_negates_root": reflection * oriented == -oriented,
        }
        if not all(checks.values()):
            raise ArithmeticError(f"facet {facet} failed")
        signs[sign] += 1
        types[(square, divisibility)] += 1
        rows.append([int(value) for value in oriented])
        records.append(
            {
                "facet": facet,
                "orientation_sign": sign,
                "oriented_coordinates": [int(value) for value in oriented],
                "square": square,
                "divisibility": divisibility,
                "interior_pairing": pairing,
                "checks": checks,
            }
        )
    root_matrix = sp.Matrix(rows)
    type_payload = {
        f"{square},{divisibility}": count
        for (square, divisibility), count in sorted(types.items())
    }
    if type_payload != archived_types:
        raise ArithmeticError("primitive mirror histogram changed")
    smith = _smith(root_matrix)
    rank = root_matrix.rank()
    return {
        "strict_interior_coordinates": [int(value) for value in interior],
        "strict_interior_square": int((interior.T * gram * interior)[0]),
        "minimum_wall_pairing": min(row["interior_pairing"] for row in records),
        "orientation_sign_histogram": {
            str(sign): count for sign, count in sorted(signs.items())
        },
        "primitive_mirror_type_histogram": type_payload,
        "wall_records": records,
        "wall_span": {
            "rank": rank,
            "smith_nonzero_invariants": smith,
            "index_in_carrier_lattice": (
                math.prod(smith) if rank == gram.rows else None
            ),
            "oriented_root_matrix_sha256": _payload_sha256(rows),
        },
    }


def _explicit_c2_k3_embedding(gram: sp.Matrix) -> dict[str, Any]:
    U = sp.Matrix([[0, 1], [1, 0]])
    E8m = sp.simplify(gram[2:10, 2:10] / 2)
    ambient = sp.diag(U, U, U, E8m, E8m)
    embedding = sp.zeros(22, 12)
    embedding[0, 0] = embedding[1, 1] = 1
    for index in range(8):
        embedding[6 + index, 2 + index] = 1
        embedding[14 + index, 2 + index] = 1
    embedding[2, 10] = 1
    embedding[3, 10] = -1
    embedding[4, 11] = 1
    embedding[5, 11] = -1

    complement = sp.zeros(22, 10)
    for index in range(8):
        complement[6 + index, index] = 1
        complement[14 + index, index] = -1
    complement[2, 8] = complement[3, 8] = 1
    complement[4, 9] = complement[5, 9] = 1
    complement_gram = sp.diag(gram[2:10, 2:10], sp.Matrix([[2]]), sp.Matrix([[2]]))
    combined = embedding.row_join(complement)
    checks = {
        "ambient_K3_lattice": (
            ambient.det() == -1 and exact_inertia(ambient) == (3, 19, 0)
        ),
        "embedding_isometric": embedding.T * ambient * embedding == gram,
        "embedding_primitive": _smith(embedding) == [1] * 12,
        "complement_isometric": complement.T * ambient * complement == complement_gram,
        "orthogonal": embedding.T * ambient * complement == sp.zeros(12, 10),
        "gluing_index": abs(int(combined.det())) == abs(int(gram.det())) == 1024,
    }
    return {
        "criterion": "explicit full primitive embedding",
        "ambient": "Lambda_K3=U^3+E8(-1)^2",
        "embedding_matrix_columns": _matrix_payload(embedding),
        "description": (
            "U as a summand, E8(-2) diagonally in E8(-1)^2, and the two "
            "<-2> generators as e-f in the remaining two U summands"
        ),
        "orthogonal_complement": "E8(-2)+<2>^2",
        "orthogonal_complement_gram": _matrix_payload(complement_gram),
        "orthogonal_complement_matrix_columns": _matrix_payload(complement),
        "combined_index": abs(int(combined.det())),
        "checks": checks,
    }


def _stable_range_embedding(gram: sp.Matrix) -> dict[str, Any]:
    smith = _smith(gram)
    length = sum(value > 1 for value in smith)
    rank = gram.rows
    complement_rank = 22 - rank
    signature = exact_inertia(gram)
    checks = {
        "even": all(gram[index, index] % 2 == 0 for index in range(rank)),
        "hyperbolic_signature": signature == (1, rank - 1, 0),
        "both_complement_signs_available": 3 - 1 > 0 and 19 - (rank - 1) > 0,
        "Nikulin_stable_range": complement_rank >= length + 2,
    }
    return {
        "criterion": "Nikulin primitive-embedding stable range",
        "ambient": "Lambda_K3=U^3+E8(-1)^2",
        "rank": rank,
        "discriminant_length": length,
        "complement_rank": complement_rank,
        "inequality": f"{complement_rank} >= {length}+2",
        "complement_signature": [2, 20 - rank],
        "primitive_embedding_exists": all(checks.values()),
        "checks": checks,
    }


def _carrier_lattice(
    name: str,
    gram: sp.Matrix,
    embedding: dict[str, Any],
) -> dict[str, Any]:
    smith = _smith(gram)
    return {
        "name": name,
        "gram": _matrix_payload(gram),
        "rank": gram.rows,
        "signature": list(exact_inertia(gram)),
        "determinant_absolute": abs(int(gram.det())),
        "smith_invariants": smith,
        "discriminant": _discriminant_data(gram),
        "primitive_K3_embedding": embedding,
        "very_general_period_consequence": {
            "moduli_dimension": 20 - gram.rows,
            "full_Picard_lattice": f"NS(X)={name}",
            "reason": (
                "period surjectivity plus avoidance of countably many extra "
                "algebraic hyperplanes"
            ),
        },
        "minimal_even_scaled_dual": _minimal_even_scaled_dual(gram),
    }


def build_certificate() -> dict[str, Any]:
    paths = {
        "fermat_source": FERMAT_SOURCE,
        "f3_fixed": F3_FIXED,
        "f3_c2": F3_C2,
        "f3_c2_vinberg": F3_C2_VINBERG,
        "f3_v4": F3_V4,
        "f3_v4_cone": F3_V4_CONE,
        "f3_mirrors": F3_MIRRORS,
        "ss2_3b": SS2_3B,
        "ss2_3b_lineage": SS2_3B_LINEAGE,
    }
    for key, path in paths.items():
        if _sha256(path) != EXPECTED_SHA256[key]:
            raise ArithmeticError(f"pinned {key} archive changed")

    fermat_source = parse_source(FERMAT_SOURCE)
    f3_fixed = json.loads(F3_FIXED.read_text(encoding="utf-8"))
    c2 = json.loads(F3_C2.read_text(encoding="utf-8"))
    c2_vinberg = json.loads(F3_C2_VINBERG.read_text(encoding="utf-8"))
    v4 = json.loads(F3_V4.read_text(encoding="utf-8"))
    v4_cone = json.loads(F3_V4_CONE.read_text(encoding="utf-8"))
    mirrors = json.loads(F3_MIRRORS.read_text(encoding="utf-8"))
    ss2 = json.loads(SS2_3B.read_text(encoding="utf-8"))
    ss2_lineage = json.loads(SS2_3B_LINEAGE.read_text(encoding="utf-8"))
    fixed_cases = {row["id"]: row for row in f3_fixed["cases"]}
    published_ample = fermat_source["ample"]

    c2_gram = sp.Matrix(c2["fixed_lattice_identification"]["target_gram"])
    c2_fixed_basis = sp.Matrix(
        fixed_cases["char3-C2-H11"]["basis_rows_in_published_NS_basis"]
    )
    c2_isometry = sp.Matrix(
        c2["fixed_lattice_identification"]["explicit_isometry_P"]
    )
    c2_fixed_interior = c2_isometry * C2_INTERIOR
    c2_published_interior = c2_fixed_basis.T * c2_fixed_interior
    c2_ample_transport = {
        "standard_carrier_coordinates": [int(value) for value in C2_INTERIOR],
        "fixed_lattice_coordinates": [int(value) for value in c2_fixed_interior],
        "equals_published_Fermat_ample_class": (
            c2_published_interior == published_ample.T
        ),
    }
    c2_walls = _wall_system(
        c2_gram,
        [row["standard_lattice_coordinates"] for row in c2["induced_206_wall_reflection_subgroup"]["mirrors"]],
        C2_INTERIOR,
        {"-2,1": 12, "-2,2": 164, "-4,2": 30},
    )
    c2_embedding = _explicit_c2_k3_embedding(c2_gram)
    c2_lattice = _carrier_lattice(
        "U+E8(-2)+<-2>^2", c2_gram, c2_embedding
    )

    stable = v4["minimal_reflection_stable_lattice"]
    v4_gram = sp.Matrix(stable["gram"])
    v4_walls = _wall_system(
        v4_gram,
        [row["stable_lattice_coordinates"] for row in v4["induced_328_wall_reflection_subgroup"]["mirrors_in_stable_and_sk_bases"]],
        V4_INTERIOR,
        {"-2,1": 28, "-4,2": 24, "-6,3": 276},
    )
    v4_lattice = _carrier_lattice(
        (
            "negative of Scharlau--Kirschmer List11/Gram11 entry 20 "
            "(signature (1,11), determinant 48)"
        ),
        v4_gram,
        _stable_range_embedding(v4_gram),
    )
    closure = mirrors["v4_h11"]["minimal_reflection_stable_overlattice"]
    closure_basis = sp.Matrix(closure["basis_columns_in_fixed_lattice_rational_coordinates"])
    v4_fixed_basis = sp.Matrix(
        fixed_cases["char3-V4-H11"]["basis_rows_in_published_NS_basis"]
    )
    v4_fixed_interior = closure_basis * V4_INTERIOR
    v4_published_interior = v4_fixed_basis.T * v4_fixed_interior
    v4_ample_transport = {
        "stable_carrier_coordinates": [int(value) for value in V4_INTERIOR],
        "fixed_lattice_coordinates": [int(value) for value in v4_fixed_interior],
        "equals_published_Fermat_ample_class": (
            v4_published_interior == published_ample.T
        ),
    }
    original_v4_gram = sp.simplify(
        closure_basis.T.inv() * v4_gram * closure_basis.inv()
    )

    ss2_gram = sp.Matrix(ss2["fixed_lattice"]["gram"])
    ss2_raw_roots = [
        row["primitive_fixed_coordinates"]
        for row in ss2["restricted_chamber"]["primitive_facet_lattice_records"]
    ]
    ss2_root_sum = sum(
        (sp.Matrix(root) for root in ss2_raw_roots), sp.zeros(12, 1)
    )
    if any(value % 3 for value in ss2_root_sum):
        raise ArithmeticError("SS2-3B root sum is not divisible by three")
    ss2_interior = ss2_root_sum / 3
    ss2_walls = _wall_system(
        ss2_gram,
        ss2_raw_roots,
        ss2_interior,
        {"-2,1": 12, "-6,3": 10},
    )
    ss2_lattice = _carrier_lattice(
        "U+A2(-1)^5", ss2_gram, _stable_range_embedding(ss2_gram)
    )

    chamber_c2 = {
        "id": "F3-C2-H11-206",
        "dimension": 11,
        "facet_count": 206,
        "original_wall_lattice": c2_lattice,
        "certified_reflection_stable_overclosures": [
            {"index": 1, "carrier": "original lattice already reflection-stable"}
        ],
        "oriented_wall_system": c2_walls,
        "published_ample_transport": c2_ample_transport,
        "volume": {
            "proof": "literal Vinberg edge-endpoint and full zero-neighborhood criterion",
            "classification": c2_vinberg["vinberg_decision"],
            "finite_vertices": c2_vinberg["finite_vertex_count"],
            "cusps": c2_vinberg["ideal_vertex_count"],
            "elliptic_rank10_edges": c2_vinberg["elliptic_edge_count"],
            "bad_edges": c2_vinberg["bad_edge_count"],
            "bad_parabolic_neighborhoods": c2_vinberg[
                "invalid_full_zero_neighborhood_count"
            ],
        },
        "status": "new characteristic-zero very-general full-Picard K3 carrier",
    }
    chamber_v4 = {
        "id": "F3-V4-H11-328",
        "dimension": 11,
        "facet_count": 328,
        "original_wall_lattice": {
            "name": "saturated characteristic-three divisor-fixed lattice",
            "gram": _matrix_payload(original_v4_gram),
            "determinant_absolute": abs(int(original_v4_gram.det())),
            "smith_invariants": _smith(original_v4_gram),
            "reflection_failure": "36 of 328 reflections are nonintegral",
            "minimal_even_scaled_dual": _minimal_even_scaled_dual(original_v4_gram),
        },
        "certified_reflection_stable_overclosures": [
            {
                "index": 1,
                "even": True,
                "determinant_absolute": 768,
                "all_reflections_integral": False,
            },
            {
                "index": 4,
                "unique_minimal_for_all_328_reflections": True,
                "basis_columns_in_original_Q_space": _matrix_payload(closure_basis),
                "carrier_lattice": v4_lattice,
                "all_reflections_integral": True,
            },
        ],
        "oriented_wall_system_on_closure": v4_walls,
        "published_ample_transport": v4_ample_transport,
        "volume": {
            "proof": "complete exact GMP double description",
            "classification": v4_cone["decision"]["classification"],
            "finite_vertices": v4_cone["exact_double_description"]["finite_vertex_count"],
            "cusps": v4_cone["exact_double_description"]["ideal_vertex_count"],
            "hyperideal_vertices": v4_cone["exact_double_description"]["hyperideal_vertex_count"],
        },
        "status": "new characteristic-zero very-general full-Picard K3 carrier after the mandatory stable closure",
    }
    chamber_ss2 = {
        "id": "SS2-3B-H11-22",
        "dimension": 11,
        "facet_count": 22,
        "original_wall_lattice": ss2_lattice,
        "certified_reflection_stable_overclosures": [
            {"index": 1, "carrier": "original lattice already reflection-stable"}
        ],
        "oriented_wall_system": ss2_walls,
        "strict_interior_derivation": {
            "formula": "one third of the sum of all 22 archived primitive roots",
            "coordinates": [int(value) for value in ss2_interior],
            "square": int((ss2_interior.T * ss2_gram * ss2_interior)[0]),
        },
        "standard_model_identification": {
            "name": ss2_lineage["census_match"]["identified_lattice"],
            "standard_model_gram": ss2_lineage["lattice"]["standard_model_gram"],
            "class_number_one_source_semantics": ss2_lineage["census_match"][
                "class_number_one_source_semantics"
            ],
        },
        "volume": {
            "proof": ss2["vinberg"]["criterion"],
            "classification": ss2["vinberg"]["decision"],
            "finite_vertices": ss2["vinberg"]["finite_vertex_count"],
            "cusps": ss2["vinberg"]["ideal_vertex_count"],
            "elliptic_rank10_edges": ss2["vinberg"]["elliptic_edge_count"],
            "bad_edges": ss2["vinberg"]["bad_edge_count"],
            "bad_parabolic_neighborhoods": ss2["vinberg"]["bad_parabolic_neighborhood_count"],
        },
        "status": "new characteristic-zero very-general full-Picard K3 carrier",
    }

    checks = {
        "all_upstream_certificates_green": (
            c2["all_checks_pass"]
            and v4["all_checks_pass"]
            and mirrors["all_checks_pass"]
            and ss2["all_checks_pass"]
            and f3_fixed["all_checks_pass"]
            and ss2_lineage["all_checks_pass"]
        ),
        "C2_interior_is_published_Fermat_ample": c2_ample_transport[
            "equals_published_Fermat_ample_class"
        ],
        "V4_interior_is_published_Fermat_ample": v4_ample_transport[
            "equals_published_Fermat_ample_class"
        ],
        "SS2_interior_derived_from_complete_root_set": (
            ss2_walls["strict_interior_square"] == 14
            and ss2_walls["minimum_wall_pairing"] == 1
        ),
        "SS2_standard_model_lineage_replayed": (
            ss2_gram == sp.Matrix(ss2_lineage["lattice"]["candidate_gram"])
            and ss2_lineage["checks"]["entry_27_lattice_genus_match"]
            and ss2_lineage["lattice"]["checks"][
                "finite_discriminant_forms_match"
            ]
            and ss2_lineage["census_match"]["class_number_one_source_semantics"]
        ),
        "C2_explicit_primitive_K3_embedding": all(c2_embedding["checks"].values()),
        "V4_stable_range_K3_embedding": v4_lattice["primitive_K3_embedding"]["primitive_embedding_exists"],
        "SS2_3B_stable_range_K3_embedding": ss2_lattice["primitive_K3_embedding"]["primitive_embedding_exists"],
        "all_three_wall_sets_generate_full_carrier": (
            c2_walls["wall_span"]["rank"] == 12
            and c2_walls["wall_span"]["index_in_carrier_lattice"] == 1
            and v4_walls["wall_span"]["rank"] == 12
            and v4_walls["wall_span"]["index_in_carrier_lattice"] == 1
            and ss2_walls["wall_span"]["rank"] == 12
            and ss2_walls["wall_span"]["index_in_carrier_lattice"] == 1
        ),
        "all_three_volume_gates_green": (
            chamber_c2["volume"]["bad_edges"] == 0
            and chamber_c2["volume"]["bad_parabolic_neighborhoods"] == 0
            and chamber_v4["volume"]["hyperideal_vertices"] == 0
            and chamber_ss2["volume"]["bad_edges"] == 0
            and chamber_ss2["volume"]["bad_parabolic_neighborhoods"] == 0
        ),
    }
    return {
        "schema": SCHEMA,
        "status_date": "2026-08-07",
        "sources": {
            key: {"path": str(path.relative_to(ROOT)), "sha256": EXPECTED_SHA256[key]}
            for key, path in paths.items()
        },
        "theorem": (
            "The finite-volume noncompact F3-C2 H11/206, F3-V4 H11/328, "
            "and SS2-3B H11/22 chambers all have characteristic-zero "
            "eight-dimensional very-general complex-K3 families whose full "
            "Picard lattices carry every oriented wall and reflection.  The "
            "328-wall case requires, and uses, its canonical even index-four "
            "reflection-stable closure."
        ),
        "chambers": [chamber_c2, chamber_v4, chamber_ss2],
        "scope": {
            "positive_gain": "three known cofinite H11 chambers gain very-general characteristic-zero full-Picard K3 carriers",
            "not_new_abstract_polytopes": True,
            "not_full_nef": "mixed square-minus-four and square-minus-six walls prevent a complete K3 nef-wall interpretation",
            "not_surface_automorphisms": "integral Picard reflections are not asserted to extend to K3 automorphisms",
            "all_wall_sets_generate_full_Picard": True,
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
        raise ArithmeticError("three-H11 K3 carrier certificate failed")
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.exists() or args.output.read_text(encoding="utf-8") != rendered:
            raise ArithmeticError("archived three-H11 carrier changed")
        print(f"verified {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(f"wrote {args.output}")
    print(f"payload_sha256={_payload_sha256(payload)}")


if __name__ == "__main__":
    main()
