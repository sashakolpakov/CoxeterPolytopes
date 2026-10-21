#!/usr/bin/env python3
"""Certify a U+D8+A1 K3 scaled-dual carrier of A18 C04 H10/28.

The A18 C04 fixed lattice has a unique reflection-stable odd index-two
closure M.  Its doubled dual N=M^vee(2) is even of determinant eight.  An
explicit determinant-one isometry identifies N with Roulleau's actual K3
Picard stock U+D8(-1)+A1(-1).  Primitive covectors transport all 28 oriented
halfspaces, generate the whole target Picard lattice, and give integral
reflections.  The chamber is already known; this certifies a second,
different K3 divisor-lattice carrier rather than a new abstract polytope.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form

from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
CONSTRUCTION = (
    ROOT / "data/enumerations/alexeev_remaining_high_fixed_section_construction.json"
)
STABLE = (
    ROOT
    / "data/enumerations/alexeev_remaining_high_fixed_section_stable_overlattices.json"
)
VINBERG = (
    ROOT / "data/enumerations/alexeev_remaining_high_fixed_section_vinberg.json"
)
K3_SOURCE = (
    ROOT / "data/enumerations/roulleau_rank11_direct_graph_denominator.json"
)
OUTPUT = (
    ROOT / "data/enumerations/a18_c04_ud8a1_scaled_dual_carrier.json"
)

SCHEMA = "coxeter.a18_c04.ud8a1_scaled_dual_divisor_carrier.v1"
EXPECTED_SHA256 = {
    "construction": "c339b56a6d0b68c0939dab85ba6390c7429a30e56eeafec822ac0d5373ec0e91",
    "stable": "882204e5457e96a3930399201dd061bd6bc35839497693a384d0b931f10924a0",
    "vinberg": "d5e56b83f66182a0c52758ee2be8d09217c66f0e4ffe4bfc82bc73b867f2ab0c",
    "k3_source": "d82e31de7dd9d1f567cba9cc9a98588168628714521842c87943a9c76ec20d43",
}

# Coordinates in N=M^vee(2) are sent to Roulleau's archived A1,...,A11
# basis of U+D8(-1)+A1(-1).  It was found by splitting a primitive U in
# both lattices and applying PARI/GP qfisom to their positive-definite
# rank-nine complements.  The permanent replay verifies determinant one and
# X^T G_K3 X=G_N directly, so it has no external-program dependency.
DOUBLED_DUAL_TO_UD8A1 = sp.Matrix(
    [
        [-695, 570, 656, -1209, 1063, 1338, 1663, 1528, 1343, 1857, -510],
        [-811, 664, 766, -1410, 1241, 1562, 1942, 1785, 1569, 2168, -596],
        [-769, 630, 726, -1337, 1176, 1481, 1841, 1692, 1487, 2055, -565],
        [-695, 570, 656, -1208, 1062, 1337, 1662, 1527, 1342, 1856, -510],
        [-621, 510, 586, -1080, 949, 1195, 1486, 1365, 1200, 1659, -456],
        [-515, 424, 486, -896, 787, 991, 1232, 1132, 995, 1376, -378],
        [-552, 454, 520, -961, 844, 1063, 1321, 1213, 1066, 1475, -405],
        [-510, 420, 481, -889, 780, 982, 1220, 1120, 984, 1363, -374],
        [-499, 411, 471, -870, 763, 961, 1194, 1096, 963, 1334, -366],
        [-621, 510, 586, -1081, 950, 1196, 1486, 1365, 1200, 1660, -456],
        [-190, 154, 180, -329, 291, 366, 455, 419, 369, 508, -140],
    ]
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _payload(matrix: sp.Matrix) -> list[list[int | str]]:
    return [
        [
            int(value) if sp.Rational(value).q == 1 else str(sp.Rational(value))
            for value in matrix.row(row)
        ]
        for row in range(matrix.rows)
    ]


def _payload_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _smith(matrix: sp.Matrix) -> list[int]:
    normal = smith_normal_form(matrix, domain=ZZ)
    return [abs(int(normal[index, index])) for index in range(matrix.rows)]


def _content(vector: sp.Matrix) -> int:
    return math.gcd(*(abs(int(value)) for value in vector))


def _reflection(gram: sp.Matrix, root: sp.Matrix) -> sp.Matrix:
    square = (root.T * gram * root)[0]
    return sp.simplify(
        sp.eye(gram.rows) - sp.Rational(2, square) * root * root.T * gram
    )


def _histogram_payload(counter: Counter[tuple[int, ...]]) -> dict[str, int]:
    return {
        ",".join(str(value) for value in key): count
        for key, count in sorted(counter.items())
    }


def build_certificate() -> dict[str, Any]:
    paths = {
        "construction": CONSTRUCTION,
        "stable": STABLE,
        "vinberg": VINBERG,
        "k3_source": K3_SOURCE,
    }
    for key, path in paths.items():
        if _sha256(path) != EXPECTED_SHA256[key]:
            raise ArithmeticError(f"pinned {key} archive changed")

    construction = json.loads(CONSTRUCTION.read_text(encoding="utf-8"))
    stable = json.loads(STABLE.read_text(encoding="utf-8"))
    vinberg = json.loads(VINBERG.read_text(encoding="utf-8"))
    k3_source = json.loads(K3_SOURCE.read_text(encoding="utf-8"))

    case = next(row for row in construction["cases"] if row["id"] == "A18_4_0_C04")
    stable_case = next(row for row in stable["cases"] if row["id"] == "A18_4_0_C04")
    vinberg_case = next(row for row in vinberg["cases"] if row["id"] == "A18_4_0_C04")
    k3_case = next(row for row in k3_source["cases"] if row["id"] == "R11_UD8A1")

    gram_l = sp.Matrix(case["saturated_fixed_lattice"]["gram"])
    closure = stable_case["unique_maximal_reflection_stable_overlattice"]
    closure_basis = sp.Matrix(closure["basis_columns_in_original_fixed_lattice"])
    gram_m = sp.simplify(closure_basis.T * gram_l * closure_basis)
    gram_n = sp.simplify(2 * gram_m.inv())
    gram_k3 = sp.Matrix(k3_case["source"]["basis_gram"])
    transport = DOUBLED_DUAL_TO_UD8A1

    lattice_checks = {
        "all_four_upstream_certificates_green": (
            construction["all_checks_pass"]
            and stable["all_checks_pass"]
            and vinberg["all_checks_pass"]
            and k3_source["all_checks_pass"]
            and all(case["checks"].values())
            and k3_case["all_checks_pass"]
        ),
        "source_fixed_lattice_rank11_det1024_even": (
            exact_inertia(gram_l) == (1, 10, 0)
            and int(gram_l.det()) == 1024
            and _smith(gram_l) == [1] * 2 + [2] * 8 + [4]
            and all(int(gram_l[index, index]) % 2 == 0 for index in range(11))
        ),
        "unique_stable_odd_index_two_closure": (
            closure["unique"] is True
            and closure["index_over_original_lattice"] == 2
            and abs(closure_basis.det()) == sp.Rational(1, 2)
            and gram_m == sp.Matrix(closure["gram"])
            and gram_m.det() == 256
            and exact_inertia(gram_m) == (1, 10, 0)
            and _smith(gram_m) == [1] * 3 + [2] * 8
            and any(int(gram_m[index, index]) % 2 for index in range(11))
        ),
        "doubled_dual_integral_even_determinant_eight": (
            all(value.is_Integer for value in gram_n)
            and gram_n.det() == 8
            and exact_inertia(gram_n) == (1, 10, 0)
            and _smith(gram_n) == [1] * 8 + [2] * 3
            and all(int(gram_n[index, index]) % 2 == 0 for index in range(11))
        ),
        "target_is_actual_roulleau_K3_picard_stock": (
            k3_case["source"]["lattice"] == "U+D8(-1)+A1(-1)"
            and k3_case["source"]["picard_rank"] == 11
            and all(k3_case["source"]["checks"].values())
            and k3_case["published_full_nef_chamber"]["geometry"][
                "finite_volume_certified"
            ]
            and gram_k3.det() == 8
            and exact_inertia(gram_k3) == (1, 10, 0)
            and _smith(gram_k3) == [1] * 8 + [2] * 3
        ),
        "explicit_unimodular_isometry": (
            all(value.is_Integer for value in transport)
            and transport.det() == 1
            and transport.T * gram_k3 * transport == gram_n
        ),
    }
    if not all(lattice_checks.values()):
        raise ArithmeticError(
            f"A18 C04 scaled-dual bridge failed: "
            f"{[key for key, value in lattice_checks.items() if not value]}"
        )

    active_orbits = set(case["exact_reduction"]["active_facet_orbit_indices"])
    active_records = [
        row
        for row in case["restricted_mirrors"]["records"]
        if row["orbit_index"] in active_orbits
    ]
    if len(active_records) != 28:
        raise ArithmeticError("active C04 facet list changed")

    inverse_closure_basis = closure_basis.inv()
    records: list[dict[str, Any]] = []
    roots_m: list[sp.Matrix] = []
    roots_n: list[sp.Matrix] = []
    roots_k3: list[sp.Matrix] = []
    closure_types: Counter[tuple[int, int]] = Counter()
    target_types: Counter[tuple[int, int]] = Counter()
    source_to_target_types: Counter[tuple[int, ...]] = Counter()

    for index, source in enumerate(active_records):
        root_l = sp.Matrix(source["primitive_integral_same_direction_coordinates"])
        orientation = int(sp.Rational(source["presentation_scale_from_primitive"]))
        if orientation not in (-1, 1):
            raise ArithmeticError("unexpected non-unit C04 presentation scale")
        raw_root_m = inverse_closure_basis * root_l
        if any(not value.is_Integer for value in raw_root_m):
            raise ArithmeticError(f"facet {index} is absent from the stable closure")
        closure_content = _content(raw_root_m)
        # The orbit-sum presentation, hence the actual halfspace, is oriented
        # by `orientation`; primitive reduction uses a positive content.
        root_m = orientation * (raw_root_m / closure_content).applyfunc(int)
        div_m = _content(gram_m * root_m)
        root_n = (gram_m * root_m / div_m).applyfunc(int)
        root_k3 = transport * root_n

        square_l = int((root_l.T * gram_l * root_l)[0])
        div_l = _content(gram_l * root_l)
        square_m = int((root_m.T * gram_m * root_m)[0])
        square_n = int((root_n.T * gram_n * root_n)[0])
        div_n = _content(gram_n * root_n)
        reflection_m = _reflection(gram_m, root_m)
        reflection_n = _reflection(gram_n, root_n)
        reflection_k3 = _reflection(gram_k3, root_k3)
        factor_from_m = sp.Rational(2, div_m)
        factor_from_l = sp.Rational(2, div_m * closure_content)

        checks = {
            "source_square_archived": square_l == int(source["square"]),
            "source_pairing_column_and_orientation": (
                gram_l * (orientation * root_l)
                == sp.Matrix(source["fixed_pairing_column"])
            ),
            "primitive_in_M_N_and_K3": (
                _content(root_m) == _content(root_n) == _content(root_k3) == 1
            ),
            "scaled_dual_square_law": square_n * div_m**2 == 2 * square_m,
            "exact_positive_halfspace_transport": (
                factor_from_m > 0
                and factor_from_l > 0
                and root_n.T * gram_n * gram_m
                == factor_from_m * root_m.T * gram_m
            ),
            "closure_reflection_integral": all(
                value.is_Integer for value in reflection_m
            ),
            "contragredient_reflection": reflection_n == reflection_m.T,
            "dual_and_target_reflections_integral": (
                all(value.is_Integer for value in reflection_n)
                and all(value.is_Integer for value in reflection_k3)
            ),
            "target_reflection_intertwines": (
                reflection_k3 * transport == transport * reflection_n
            ),
            "target_reflection_unimodular_involution_isometry": (
                reflection_k3.det() == -1
                and reflection_k3**2 == sp.eye(11)
                and reflection_k3.T * gram_k3 * reflection_k3 == gram_k3
            ),
        }
        if not all(checks.values()):
            raise ArithmeticError(
                f"facet {index} failed: "
                f"{[key for key, value in checks.items() if not value]}"
            )

        closure_types[(square_m, div_m)] += 1
        target_types[(square_n, div_n)] += 1
        source_to_target_types[
            (square_l, div_l, closure_content, square_m, div_m, square_n, div_n)
        ] += 1
        roots_m.append(root_m)
        roots_n.append(root_n)
        roots_k3.append(root_k3)
        records.append(
            {
                "facet_zero_based": index,
                "source_orbit_index_one_based": source["orbit_index"],
                "source_L_coordinates_unoriented": [int(value) for value in root_l],
                "source_presentation_orientation": orientation,
                "source_L_square": square_l,
                "source_L_divisibility": div_l,
                "coordinates_in_M_before_primitive_reduction_unoriented": [
                    int(value) for value in raw_root_m
                ],
                "M_coordinate_content": closure_content,
                "primitive_oriented_M_coordinates": [int(value) for value in root_m],
                "primitive_M_square": square_m,
                "primitive_M_divisibility": div_m,
                "primitive_oriented_N_coordinates": [int(value) for value in root_n],
                "primitive_oriented_K3_coordinates": [int(value) for value in root_k3],
                "primitive_K3_square": square_n,
                "primitive_K3_divisibility": div_n,
                "positive_halfspace_factor_from_M": str(factor_from_m),
                "positive_halfspace_factor_from_original_L": str(factor_from_l),
                "checks": checks,
            }
        )

    # Pairwise Gram entries differ by positive row/column scalings, so the
    # normalized Coxeter relations and the complete chamber are unchanged.
    pairwise_transport = all(
        (roots_n[left].T * gram_n * roots_n[right])[0]
        * _content(gram_m * roots_m[left])
        * _content(gram_m * roots_m[right])
        == 2 * (roots_m[left].T * gram_m * roots_m[right])[0]
        for left in range(28)
        for right in range(28)
    )
    roots_n_matrix = sp.Matrix.hstack(*roots_n).applyfunc(int)
    roots_k3_matrix = sp.Matrix.hstack(*roots_k3).applyfunc(int)
    hnf_n = hermite_normal_form(roots_n_matrix)
    hnf_k3 = hermite_normal_form(roots_k3_matrix)
    carrier = sum(roots_n, sp.zeros(11, 1))
    carrier_k3 = transport * carrier
    carrier_pairings = [
        int((carrier.T * gram_n * root)[0]) for root in roots_n
    ]
    generation_checks = {
        "all_28_active_facets_transport": len(records) == 28,
        "exact_pairwise_positive_Gram_rescaling": pairwise_transport,
        "primitive_N_roots_generate_full_lattice": hnf_n == sp.eye(11),
        "primitive_K3_roots_generate_full_Picard_lattice": hnf_k3 == sp.eye(11),
        "sum_is_strict_timelike_interior": (
            int((carrier.T * gram_n * carrier)[0]) == 508
            and min(carrier_pairings) == 3
            and all(value > 0 for value in carrier_pairings)
            and carrier_k3.T * gram_k3 * carrier_k3
            == carrier.T * gram_n * carrier
        ),
        "closure_type_histogram": closure_types
        == Counter({(-4, 2): 16, (-1, 1): 8, (-2, 1): 4}),
        "target_type_histogram": target_types
        == Counter({(-2, 1): 16, (-2, 2): 8, (-4, 2): 4}),
    }
    if not all(generation_checks.values()):
        raise ArithmeticError(
            f"generation/interior audit failed: "
            f"{[key for key, value in generation_checks.items() if not value]}"
        )

    vinberg_checks = {
        "literal_case_green": vinberg_case["all_checks_pass"] is True,
        "finite_volume_noncompact": (
            vinberg_case["finite_volume"] is True
            and vinberg_case["compact"] is False
            and vinberg_case["classification"]
            == "finite_volume_noncompact_Coxeter_polytope"
        ),
        "every_elliptic_edge_has_two_endpoints": (
            vinberg_case["bad_edge_count"] == 0
            and vinberg_case["elliptic_edge_count"] == 6956
        ),
        "all_parabolic_neighborhoods_good": (
            vinberg_case["bad_parabolic_neighborhood_count"] == 0
            and vinberg_case["connected_affine_component_count"] == 204
        ),
        "finite_and_ideal_vertices": (
            vinberg_case["finite_vertex_count"] == 1087
            and vinberg_case["ideal_vertex_count"] == 74
        ),
    }
    if not all(vinberg_checks.values()):
        raise ArithmeticError(f"inherited Vinberg audit failed: {vinberg_checks}")

    checks = {
        "lattice_bridge": all(lattice_checks.values()),
        "all_facet_records": all(all(row["checks"].values()) for row in records),
        "full_generation_and_interior": all(generation_checks.values()),
        "literal_Vinberg_finite_volume": all(vinberg_checks.values()),
    }
    return {
        "schema": SCHEMA,
        "status_date": "2026-08-07",
        "id": "A18-C04-H10-28-UD8A1-SCALED-DUAL-CARRIER",
        "upstream": {
            key: {
                "path": str(path.relative_to(ROOT)),
                "sha256": EXPECTED_SHA256[key],
            }
            for key, path in paths.items()
        },
        "theorem": (
            "The finite-volume noncompact A18 C04 H10/28 chamber has a "
            "second K3 divisor-lattice carrier.  After its unique "
            "reflection-stable odd index-two closure M, the doubled dual "
            "M^vee(2) is explicitly the actual Roulleau Picard lattice "
            "U+D8(-1)+A1(-1); all 28 positive halfspaces and reflections "
            "transport and the primitive mirrors generate the full Picard lattice."
        ),
        "stable_closure": {
            "operation": "unique reflection-stable odd index-two closure before duality",
            "basis_columns_in_original_fixed_Q_space": _payload(closure_basis),
            "gram": _payload(gram_m),
            "determinant": int(gram_m.det()),
            "smith_invariants": _smith(gram_m),
            "parity": "odd",
            "primitive_mirror_types_square_divisibility": _histogram_payload(
                closure_types
            ),
        },
        "scaled_dual_bridge": {
            "operation": "N=M^vee(2), with Gram 2*Gram(M)^(-1)",
            "gram_N": _payload(gram_n),
            "N_determinant": int(gram_n.det()),
            "N_smith_invariants": _smith(gram_n),
            "N_parity": "even",
            "target": {
                "variety": "Roulleau's complex K3 surface stock",
                "divisor_lattice": "NS(X)=U+D8(-1)+A1(-1)",
                "archived_basis": "Roulleau curve classes A1,...,A11",
                "archived_basis_gram": _payload(gram_k3),
                "published_complete_nef_wall_count": k3_case["source"]["wall_count"],
            },
            "N_to_K3_isometry": _payload(transport),
            "isometry_determinant": int(transport.det()),
            "primitive_target_mirror_types_square_divisibility": (
                _histogram_payload(target_types)
            ),
            "source_to_closure_to_target_type_histogram": (
                _histogram_payload(source_to_target_types)
            ),
            "positive_halfspace_identity": (
                "For oriented primitive q in M, p=Gram(M)q/div_M(q) and "
                "Phi(x)=Gram(M)x satisfy (Phi(x),p)_N="
                "(2/div_M(q))*(x,q)_M.  The archived presentation normal "
                "is a positive multiple of q after its recorded sign is applied."
            ),
            "checks": lattice_checks,
        },
        "facet_records": records,
        "full_lattice_generation": {
            "N_root_matrix_HNF": _payload(hnf_n),
            "K3_root_matrix_HNF": _payload(hnf_k3),
            "generation_index": 1,
            "strict_interior_N_coordinates": [int(value) for value in carrier],
            "strict_interior_K3_coordinates": [int(value) for value in carrier_k3],
            "strict_interior_square": int((carrier.T * gram_n * carrier)[0]),
            "minimum_wall_pairing": min(carrier_pairings),
            "checks": generation_checks,
        },
        "inherited_literal_Vinberg_audit": {
            "criterion": vinberg_case["criterion"],
            "dimension": vinberg_case["dimension"],
            "facet_count": vinberg_case["facet_count"],
            "finite_vertex_count": vinberg_case["finite_vertex_count"],
            "ideal_vertex_count": vinberg_case["ideal_vertex_count"],
            "elliptic_edge_count": vinberg_case["elliptic_edge_count"],
            "connected_affine_component_count": vinberg_case[
                "connected_affine_component_count"
            ],
            "bad_edge_count": vinberg_case["bad_edge_count"],
            "bad_parabolic_neighborhood_count": vinberg_case[
                "bad_parabolic_neighborhood_count"
            ],
            "classification": vinberg_case["classification"],
            "checks": vinberg_checks,
        },
        "lineage_and_scope": {
            "existing_geometric_carrier": (
                "a special Picard-rank-19 K3 with symplectic Nikulin "
                "involution and fixed rank-11 lattice of determinant 1024"
            ),
            "new_positive_gain": (
                "a different full-rank carrier on the actual determinant-eight "
                "Roulleau K3 Picard lattice, with coordinate-level scaled-dual "
                "wall transport"
            ),
            "known_abstract_polytope": True,
            "not_a_new_abstract_polytope": True,
            "not_the_complete_K3_nef_chamber": (
                "the imported arrangement has 28 walls including four "
                "primitive square-minus-four mirrors; the archived complete "
                "Roulleau nef chamber has 14 curve walls"
            ),
            "not_an_automorphism_fixed_chamber_claim": True,
            "minus_two_divisors_not_asserted_irreducible": True,
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
        raise ArithmeticError("A18 C04 scaled-dual carrier certificate failed")
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.exists():
            raise ArithmeticError(f"missing archive: {args.output}")
        if args.output.read_text(encoding="utf-8") != serialized:
            raise ArithmeticError("archived A18 C04 scaled-dual payload changed")
        print(f"verified {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized, encoding="utf-8")
        print(f"wrote {args.output}")
    print(f"payload_sha256={_payload_sha256(payload)}")


if __name__ == "__main__":
    main()
