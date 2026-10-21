#!/usr/bin/env python3
"""Certify the two-tier K3 carrier of Alexeev's A15_11 H10/21 chamber.

Tier A is a very general K3 surface Y with
NS(Y)=U+E8(-2)+<-2>.  All 21 facet reflections act integrally on its full
Neron--Severi lattice.  Tier B is Shimada's named singular K3 X_2: the same
lattice occurs primitively as NS(X_2)^iota for a rational involution, but the
ten square-minus-four reflections do not extend to the rank-20 NS(X_2).
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import shutil
from typing import Any

import sympy as sp

from scripts.certify_alexeev_20_2_1_15_9_1_exact_cones import (
    _parse_cdd_rays,
    _replay_case,
    _run_cdd,
)
from scripts.certify_alexeev_remaining_stable_overlattices import (
    _matrix,
    _overlattice_basis,
    _q_histogram,
    _reflection,
    _smith,
    _stable_subgroups,
)
from scripts.certify_roulleau_r14_c2t_pure_order4_obstruction import (
    saturated_row_kernel,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
CONSTRUCTION = ROOT / "data/enumerations/alexeev_remaining_high_fixed_section_construction.json"
STABLE = ROOT / "data/enumerations/alexeev_remaining_high_fixed_section_stable_overlattices.json"
VINBERG = ROOT / "data/enumerations/alexeev_remaining_high_fixed_section_vinberg.json"
LINEAGE = ROOT / "data/enumerations/alexeev_remaining_high_fixed_section_lineage.json"
X2 = ROOT / "data/enumerations/shimada_aut3singk3_x2_rational_involution_h10_37.json"
X2_STOCK = ROOT / "data/enumerations/shimada_aut3singk3_x2_478_wall_stock.json"
OUTPUT = ROOT / "data/enumerations/a15_11_h10_21_x2_k3_carrier.json"
SCHEMA = "coxeter.a15_11.h10_21.x2_k3_carrier.v1"

EXPECTED_SHA256 = {
    "construction": "c339b56a6d0b68c0939dab85ba6390c7429a30e56eeafec822ac0d5373ec0e91",
    "stable": "882204e5457e96a3930399201dd061bd6bc35839497693a384d0b931f10924a0",
    "vinberg": "d5e56b83f66182a0c52758ee2be8d09217c66f0e4ffe4bfc82bc73b867f2ab0c",
    "lineage": "dea8024f9e6e58abd0483e339774e8023b1bd511a6ac208d7db6e5e709187173",
    "x2": "12914344403b7d3d6b69a2866c4150f0a0420a64914e4bc8584fa372799599a7",
    "x2_stock": "10d1c767d21d7194c6ccb603c5dff29520201954b4a07fd3c5f214d452ea2f74",
}

# Columns are standard-model coordinates expressed in the internal index-four
# overclosure basis: P^T G_N P = T.
STANDARD_TO_OVERCLOSURE = sp.Matrix([
    [-13,-13,14,1,0,-1,6,0,-17,-2,8],[-12,-12,11,0,1,0,5,0,-15,-3,8],
    [-1,-1,0,0,0,0,1,0,-1,-1,1],[-9,-9,9,0,0,0,4,0,-11,-2,6],
    [7,7,-7,0,-1,1,-4,1,8,2,-4],[-7,-7,8,1,0,-1,4,-1,-9,-1,4],
    [19,19,-20,0,0,0,-8,0,24,4,-12],[26,25,-26,0,-2,2,-12,0,32,6,-16],
    [6,6,-7,-1,0,1,-3,0,8,1,-3],[-8,-8,8,0,1,-1,4,0,-10,-2,5],
    [2,2,-2,0,0,0,-1,0,3,0,-1],
])

# Direct inclusion M -> T=U+E8(-2)+<-2>; J^T T J=G_M and det J=4.
M_TO_STANDARD = sp.Matrix([
    [-2,0,0,0,2,2,0,2,0,4,0],[-2,0,0,0,2,2,2,-2,0,0,0],
    [-3,0,-1,2,2,2,1,0,-1,2,0],[-3,1,-1,3,2,2,1,0,-2,1,0],
    [-5,1,-2,4,3,3,2,0,-2,3,-1],[-7,1,-2,6,4,4,2,0,-3,3,-2],
    [-5,0,-1,5,3,3,2,0,-2,3,-2],[-3,0,-1,4,2,1,1,0,-2,2,-1],
    [-1,0,-1,3,0,0,0,0,-1,0,0],[0,0,-1,2,0,0,0,0,0,0,-1],
    [0,0,0,0,1,1,0,0,1,1,1],
])

# Columns are T coordinates expressed in the actual archived fixed basis of
# F=NS(X_2)^iota: Q^T G_F Q=T.
STANDARD_TO_X2_FIXED = sp.Matrix([
    [7,154,21,0,0,0,14,0,-35,0,-47],[-14,-308,-40,0,0,0,-28,-1,69,1,95],
    [-25,-550,-71,0,-1,0,-50,0,123,0,170],[-10,-220,-28,0,-1,1,-21,0,49,1,68],
    [-15,-330,-42,0,0,0,-31,0,74,1,102],[23,506,67,0,0,0,46,0,-113,1,-156],
    [25,550,72,1,1,-1,50,0,-123,1,-170],[-4,-88,-12,0,0,0,-8,0,20,0,27],
    [5,110,16,0,0,0,9,0,-25,3,-34],[-10,-220,-29,0,0,0,-20,0,49,0,68],
    [-1,-23,-3,0,0,0,-2,0,5,0,7],
])


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _payload(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)]


def _payload_hash(value: Any) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def _root_reflection(gram: sp.Matrix, root: sp.Matrix) -> sp.Matrix:
    norm = (root.T * gram * root)[0]
    return sp.simplify(sp.eye(gram.rows) - sp.Rational(2, norm) * root * root.T * gram)


def _content(vector: sp.Matrix) -> int:
    return int(sp.gcd_list([abs(int(value)) for value in vector]))


def _active_roots(case: dict[str, Any]) -> list[tuple[dict[str, Any], sp.Matrix]]:
    return [
        (record, sp.Matrix(record["primitive_integral_same_direction_coordinates"]))
        for record in case["restricted_mirrors"]["records"]
        if record["hyperbolic_status"] == "active_facet"
    ]


def _cdd_replay(case: dict[str, Any]) -> dict[str, Any]:
    executable = shutil.which("cddexec_gmp")
    if executable is None:
        raise RuntimeError("--replay-cdd requires cddexec_gmp")
    cone = _replay_case(case)
    output = _run_cdd(cone["inequality_text"], executable)
    rays = _parse_cdd_rays(output, cone["rank"])
    gram = cone["fixed_gram"]
    signs = Counter()
    ideal_incidence = Counter()
    for ray in rays:
        vector = sp.Matrix(ray)
        square = int((vector.T * gram * vector)[0])
        signs["finite" if square > 0 else "ideal" if square == 0 else "hyperideal"] += 1
        if square == 0:
            ideal_incidence[str(sum(value == 0 for value in cone["inequality"] * vector))] += 1
    result = {
        "extreme_ray_count": len(rays),
        "causal_type_histogram": dict(signs),
        "ideal_facet_incidence_histogram": dict(sorted(ideal_incidence.items())),
    }
    expected = {
        "extreme_ray_count": 657,
        "causal_type_histogram": {"finite": 607, "ideal": 50},
        "ideal_facet_incidence_histogram": {"11": 18, "12": 32},
    }
    if result != expected:
        raise ArithmeticError(f"independent cdd replay changed: {result}")
    return result


def build_certificate(*, replay_cdd: bool = False) -> dict[str, Any]:
    paths = {
        "construction": CONSTRUCTION, "stable": STABLE, "vinberg": VINBERG,
        "lineage": LINEAGE, "x2": X2, "x2_stock": X2_STOCK,
    }
    for key, path in paths.items():
        if _sha256(path) != EXPECTED_SHA256[key]:
            raise ArithmeticError(f"pinned {key} archive changed")

    construction = json.loads(CONSTRUCTION.read_text())
    stable = json.loads(STABLE.read_text())
    vinberg = json.loads(VINBERG.read_text())
    lineage = json.loads(LINEAGE.read_text())
    x2 = json.loads(X2.read_text())
    stock = json.loads(X2_STOCK.read_text())
    case = next(row for row in construction["cases"] if row["id"] == "A15_11_1_C02")
    stable_case = next(row for row in stable["cases"] if row["id"] == "A15_11_1_C02")
    vinberg_case = next(row for row in vinberg["cases"] if row["id"] == "A15_11_1_C02")
    lineage_case = next(row for row in lineage["cases"] if row["id"] == "A15_11_1_C02")
    upstream_green = {
        "construction": construction["all_checks_pass"],
        "stable": stable["all_checks_pass"],
        "vinberg": vinberg["all_checks_pass"],
        "lineage": lineage["all_checks_pass"],
        "x2": x2["all_checks_pass"],
        "x2_stock": stock["all_checks_pass"],
    }
    if not all(upstream_green.values()):
        raise ArithmeticError("an upstream archive is not green")

    gram_m = _matrix(case["saturated_fixed_lattice"]["gram"])
    active = _active_roots(case)
    roots_m = [root for _, root in active]
    reflections_m = [_reflection(gram_m, root) for root in roots_m]
    closure = _stable_subgroups(gram_m, reflections_m)
    even_four = [row for row in closure["stable_subgroups"] if row["even_overlattice"] and row["order_and_overlattice_index"] == 4]
    if len(even_four) != 1:
        raise ArithmeticError("unique even index-four overclosure changed")
    basis_n = _overlattice_basis(gram_m, closure["left_smith_inverse"], even_four[0]["cosets_in_smith_coordinates"])
    gram_n = sp.simplify(basis_n.T * gram_m * basis_n)

    fixed = x2["saturated_fixed_divisor_lattice"]
    gram_f = sp.Matrix(fixed["gram"])
    standard = sp.Matrix(fixed["standard_model_gram"])
    fixed_rows = sp.Matrix(fixed["basis_rows_in_published_NS_basis"])
    gram_x2 = sp.Matrix(stock["surface"]["published_neron_severi_gram"])
    involution = sp.Matrix(stock["ample_stabilizer"]["all_integral_20_by_20_matrices_in_supplement_order"][x2["rational_involution_fold"]["representative_zero_based"]])
    ample_x2 = sp.Matrix(stock["surface"]["published_ample_vector_a2"])
    kernel_rows = saturated_row_kernel(involution - sp.eye(20))
    kernel_change_t = kernel_rows.T.gauss_jordan_solve(fixed_rows.T)[0]
    kernel_change = kernel_change_t.T

    n_to_f = sp.simplify(STANDARD_TO_X2_FIXED * STANDARD_TO_OVERCLOSURE.inv())
    inclusion_f = sp.simplify(STANDARD_TO_X2_FIXED * M_TO_STANDARD)
    inclusion_x2 = fixed_rows.T * inclusion_f
    ample_f = fixed_rows.T.gauss_jordan_solve(ample_x2)[0]
    strict_m = sp.Matrix(
        [sp.Rational(value) for value in case["strict_fixed_interior"]["fixed_lattice_coordinates"]]
    )
    strict_f = inclusion_f * strict_m

    bridge_checks = {
        "M_rank_signature_determinant": gram_m.rows == 11 and exact_inertia(gram_m) == (1,10,0) and gram_m.det() == 8192,
        "stable_integral_orders_1_4": closure["stable_integral_subgroup_orders"] == [1,4],
        "stable_even_orders_1_4": closure["stable_even_subgroup_orders"] == [1,4],
        "unique_even_index_four": len(even_four) == 1 and abs(basis_n.det()) == sp.Rational(1,4),
        "N_even_det512_smith": gram_n.det() == 512 and all(gram_n[i,i] % 2 == 0 for i in range(11)) and _smith(gram_n) == [1,1]+[2]*9,
        "N_q_invariants": _q_histogram(gram_n) == {"0":136,"1":120,"1/2":120,"3/2":136},
        "P_congruence_to_exact_standard": STANDARD_TO_OVERCLOSURE.det() == 1 and STANDARD_TO_OVERCLOSURE.T * gram_n * STANDARD_TO_OVERCLOSURE == standard,
        "J_direct_index_four_inclusion": abs(M_TO_STANDARD.det()) == 4 and M_TO_STANDARD.T * standard * M_TO_STANDARD == gram_m,
        "J_equals_composed_map": M_TO_STANDARD == STANDARD_TO_OVERCLOSURE.inv() * basis_n.inv(),
        "Q_actual_X2_fixed_to_standard": STANDARD_TO_X2_FIXED.det() == 1 and STANDARD_TO_X2_FIXED.T * gram_f * STANDARD_TO_X2_FIXED == standard,
        "explicit_N_to_actual_F": abs(n_to_f.det()) == 1 and n_to_f.T * gram_f * n_to_f == gram_n,
        "explicit_M_index_four_in_actual_F": abs(inclusion_f.det()) == 4 and inclusion_f.T * gram_f * inclusion_f == gram_m,
        "actual_fixed_gram_inside_full_NS": fixed_rows * gram_x2 * fixed_rows.T == gram_f,
        "fixed_kernel_saturated_and_exact": all(value.is_Integer for value in kernel_change) and abs(kernel_change.det()) == 1 and fixed_rows == kernel_change * kernel_rows,
        "full_NS_inclusion_gram": inclusion_x2.T * gram_x2 * inclusion_x2 == gram_m,
        "published_ample_is_fixed_and_in_F": involution.T * ample_x2 == ample_x2 and fixed_rows.T * ample_f == ample_x2 and (ample_f.T * gram_f * ample_f)[0] == 12,
        "transported_strict_point_square": (strict_f.T * gram_f * strict_f)[0] == sp.Rational(159, 2),
    }
    if not all(bridge_checks.values()):
        raise ArithmeticError([key for key,value in bridge_checks.items() if not value])

    wall_records = []
    for index, ((source_record, root_m), reflection_m) in enumerate(zip(active, reflections_m)):
        presentation_f = inclusion_f * root_m
        content = _content(presentation_f)
        primitive_f = presentation_f / content
        primitive_x2 = fixed_rows.T * primitive_f
        norm = int((primitive_f.T * gram_f * primitive_f)[0])
        div_f = _content(gram_f * primitive_f)
        div_x2 = _content(gram_x2 * primitive_x2)
        reflection_f = _root_reflection(gram_f, primitive_f)
        reflection_x2 = _root_reflection(gram_x2, primitive_x2)
        scale = sp.Rational(source_record["presentation_scale_from_primitive"])
        oriented_f = sp.simplify(scale * primitive_f)
        strict_pairing = sp.simplify((strict_f.T * gram_f * oriented_f)[0])
        checks = {
            "primitive_in_F": _content(primitive_f) == 1,
            "primitive_in_full_NS": _content(primitive_x2) == 1,
            "same_hyperplane_after_normalization": presentation_f == content * primitive_f,
            "F_reflection_integral": all(value.is_Integer for value in reflection_f),
            "F_reflection_unimodular_involution": reflection_f.det() == -1 and reflection_f**2 == sp.eye(11),
            "F_reflection_preserves_gram": reflection_f.T * gram_f * reflection_f == gram_f,
            "F_reflection_intertwines_M": reflection_f * inclusion_f == inclusion_f * reflection_m,
            "full_NS_integrality_matches_divisibility": all(value.is_Integer for value in reflection_x2) == (2*div_x2 % abs(norm) == 0),
            "fixed_by_rational_involution": involution.T * primitive_x2 == primitive_x2,
            "oriented_strict_pairing_positive": bool(strict_pairing > 0),
        }
        if not all(checks.values()):
            raise ArithmeticError(f"wall {index+1} failed")
        wall_records.append({
            "facet_one_based": index+1,
            "source_orbit_index": source_record["orbit_index"],
            "source_orientation_scale": str(scale),
            "oriented_strict_pairing": str(strict_pairing),
            "primitive_oriented_F_coordinates": [str(value) for value in oriented_f],
            "primitive_F_coordinates": [str(value) for value in primitive_f],
            "primitive_full_NS_coordinates": [str(value) for value in primitive_x2],
            "presentation_content_in_F": content,
            "square": norm, "divisibility_in_F": div_f, "divisibility_in_full_NS": div_x2,
            "reflection_integral_on_F": True,
            "reflection_integral_on_full_NS": all(value.is_Integer for value in reflection_x2),
            "checks": checks,
        })

    type_f = Counter((row["square"],row["divisibility_in_F"]) for row in wall_records)
    type_full = Counter((row["square"],row["divisibility_in_F"],row["divisibility_in_full_NS"],row["reflection_integral_on_full_NS"]) for row in wall_records)
    expected_f = {(-2,1):4,(-2,2):7,(-4,2):10}
    expected_full = {(-2,1,1,True):4,(-2,2,1,True):7,(-4,2,1,False):10}
    strict_pairings = [sp.Rational(row["oriented_strict_pairing"]) for row in wall_records]
    expected_pairings = list(map(sp.Rational, [1,1,11,1,1,1,1,1,2,5,1,1,11,14,7,7,6,11,13,15,15]))
    if dict(type_f) != expected_f or dict(type_full) != expected_full or strict_pairings != expected_pairings:
        raise ArithmeticError("primitive wall split changed")

    vinberg_checks = {
        "classification": vinberg_case["classification"] == "finite_volume_noncompact_Coxeter_polytope",
        "dimension_facets": (vinberg_case["dimension"],vinberg_case["facet_count"]) == (10,21),
        "edge_endpoint_criterion": vinberg_case["checks"]["every_elliptic_edge_has_exactly_two_endpoints"],
        "affine_neighborhood_criterion": vinberg_case["checks"]["all_affine_zero_neighborhoods_are_full_rank_parabolic"],
        "zero_bad_audits": vinberg_case["bad_edge_count"] == vinberg_case["bad_parabolic_neighborhood_count"] == 0,
        "noncompact": vinberg_case["ideal_vertex_count"] == 50,
    }
    if not all(vinberg_checks.values()):
        raise ArithmeticError("Vinberg handoff changed")

    expected_cdd = {"extreme_ray_count":657,"causal_type_histogram":{"finite":607,"ideal":50},"ideal_facet_incidence_histogram":{"11":18,"12":32}}
    if replay_cdd:
        _cdd_replay(case)

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "theorem": "The A15_11 C02 H10/21 finite-volume chamber has a very-general K3 carrier with full NS=U+E8(-2)+<-2>; the same lattice also occurs as the actual invariant divisor lattice NS(X_2)^iota on Shimada's named X_2.",
        "reusable_carrier_criterion": {
            "statement": "If an even Lorentzian chamber lattice M has a finite-index even reflection-stable overclosure N, N embeds primitively in the K3 lattice, and the chamber passes literal finite-volume tests, then a very general N-polarized projective K3 has NS=N and carries the chamber as a full-rank divisor-reflection subsystem. Weyl re-marking puts a strict chamber point in the ample cone and K3 Riemann--Roch orients every primitive square-minus-two facet effectively.",
            "full_nef_gate": "Promotion from divisor-reflection subsystem to full nef chamber is separate and requires the ambient-root one-sign criterion; it is not proved here.",
        },
        "scope": {
            "generic_K3_Y": "all 21 reflections are integral on full NS(Y)=F, but they are not asserted to extend to Hodge isometries or surface automorphisms; extra roots may refine the chamber",
            "named_X2": "F=NS(X_2)^iota is a primitive invariant divisor sublattice of rank-20 NS(X_2); only the 11 square-minus-two reflections extend integrally to full NS(X_2), while all ten square-minus-four reflections fail",
            "not_fixed_nef_section": "the 21-wall chamber is not Shimada's 37-wall restriction of the complete 478-wall induced nef chamber",
            "novelty": "new project-level K3/divisor realization of a known arithmetic chamber; no new abstract polytope, lattice class, or maximal reflection group is claimed",
        },
        "sources": {key:{"path":str(path.relative_to(ROOT)),"sha256":EXPECTED_SHA256[key],"upstream_all_checks_pass":upstream_green[key]} for key,path in paths.items()},
        "lattice_bridge": {
            "M":{"rank":11,"signature":[1,10,0],"determinant":8192,"smith":_smith(gram_m),"gram":_payload(gram_m)},
            "N":{"isometry_class":"U+E8(-2)+<-2>","determinant":512,"smith":_smith(gram_n),"basis_columns_in_M_tensor_Q":_payload(basis_n),"gram":_payload(gram_n)},
            "standard_model_gram":_payload(standard),
            "standard_to_N_P":_payload(STANDARD_TO_OVERCLOSURE),
            "M_to_standard_J":_payload(M_TO_STANDARD),
            "standard_to_actual_X2_fixed_Q":_payload(STANDARD_TO_X2_FIXED),
            "N_to_actual_X2_fixed":_payload(n_to_f),
            "M_to_actual_X2_fixed_index_four":_payload(inclusion_f),
            "M_to_rank20_X2_NS":_payload(inclusion_x2),
            "checks":bridge_checks,
        },
        "generic_F_polarized_K3_carrier": {
            "existence": "F occurs primitively in the K3 lattice through its saturated occurrence in NS(X_2); outside the countable Noether--Lefschetz hyperplanes, the F-polarized period has NS exactly F",
            "full_NS":"NS(Y)=F=U+E8(-2)+<-2>",
            "all_21_reflections_integral_on_full_NS":True,
            "effective_minus_two_facets":11,
            "transported_strict_point_coordinates_in_F":[str(value) for value in strict_f],
            "transported_strict_point_square":"159/2",
            "oriented_wall_pairings":[str(value) for value in strict_pairings],
            "minimum_oriented_wall_pairing":"1",
            "effectivity_scope":"after Weyl re-marking and sign orientation, Riemann--Roch makes the 11 primitive square-minus-two classes effective; irreducibility is not asserted",
        },
        "named_X2_occurrence": {
            "surface":"Shimada X_2","full_NS_rank":20,"full_NS_determinant":-15,
            "involution":"rational involution, supplement representative 596 (zero based)",
            "fixed_lattice_rank":11,"fixed_lattice_determinant":512,
            "published_invariant_ample_coordinates_in_F":[str(value) for value in ample_f],
            "all_21_divisor_vectors_fixed_by_involution":True,
            "full_NS_reflection_integral_count":11,"full_NS_reflection_nonintegral_count":10,
        },
        "walls": {
            "count":21,
            "primitive_F_type_histogram":{f"{key[0]},{key[1]}":value for key,value in sorted(type_f.items())},
            "full_NS_extension_histogram":{f"{key[0]},{key[1]},{key[2]},{key[3]}":value for key,value in sorted(type_full.items())},
            "all_21_integral_on_F":True,"records":wall_records,
        },
        "finite_volume": {
            "literal_vinberg":{"classification":vinberg_case["classification"],"finite_vertices":607,"ideal_vertices":50,"elliptic_edges":4159,"connected_affine_components":132,"endpoint_type_histogram":vinberg_case["edge_endpoint_type_histogram"],"checks":vinberg_checks},
            "independent_GMP_cdd_double_description":{**expected_cdd,"hyperideal_ray_count":0,"replay":"run this certifier with --replay-cdd"},
        },
        "arithmetic_lineage": {
            "full_reflection_chamber":"Scharlau--Kirschmer List10 entry 6, 13 facets",
            "index_in_full_reflection_group":lineage_case["containing_full_reflection_chamber"]["section_to_full_index"],
            "coxeter_euler_characteristic":lineage_case["solomon"]["rational_euler_characteristic"],
            "bounded_bridge_audit":"Alexeev records the source chamber, Shimada records X_2 and its 478-wall chamber, and the prior project ledger recorded this A15 row without a K3 carrier; the exact composed bridge appears unrecorded in that bounded comparison",
            "priority_scope":"bounded evidence only; no worldwide-priority claim",
        },
        "all_checks_pass":True,
    }
    payload["certificate_payload_sha256"] = _payload_hash(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--replay-cdd", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate(replay_cdd=args.replay_cdd)
    if args.check_archive:
        if certificate != json.loads(OUTPUT.read_text()):
            raise ArithmeticError("archived A15_11/X2 K3 carrier changed")
        print(f"archive check passed: {OUTPUT}")
    else:
        OUTPUT.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
        print(OUTPUT)


if __name__ == "__main__":
    main()
