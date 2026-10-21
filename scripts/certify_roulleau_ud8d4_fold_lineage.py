#!/usr/bin/env python3
"""Close the integral/reflection lineage of the five U+D8+D4 folds.

The construction certificate proves that Roulleau's twenty-wall K3 nef stock
has five distinct cofinite fixed-section diagrams.  Here we exhaust every
integral overlattice invariant under their displayed facet reflections.  The
maximal stable lattices are then identified with the appropriate
Scharlau--Kirschmer/Vinberg full reflection lattices, keeping incomparable
overclosures and abstract CoxIter matches separate.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
from typing import Any, Sequence

import sympy as sp

from scripts.certify_alexeev_remaining_stable_overlattices import (
    _canonical_sha256,
    _discriminant_data,
    _overlattice_basis,
    _payload,
    _q_histogram,
    _reflection,
    _smith,
    _stable_subgroups,
)
from scripts.certify_jacobian_kummer_h10_maximal_reflection_lineage import (
    _parse_coxiter,
    _permutation as _relation_permutation,
)
from scripts.certify_jacobian_kummer_h10_reflection_index import (
    _connected_masks_direct,
    _euler_from_connected_components,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
CONSTRUCTION = ROOT / "data/enumerations/roulleau_ud8d4_all_subgroup_folds.json"
OUTPUT = ROOT / "data/enumerations/roulleau_ud8d4_fold_lineage.json"
SK_ROOT = Path("/private/tmp/Scharlau-Kirschmer-Lorentz")
COXITER_ROOT = Path("/private/tmp/coxiter-current/CoxIter-master")
SCHEMA = "coxeter.roulleau_ud8d4_fold_lineage.v1"


EXPECTED: dict[str, dict[str, Any]] = {
    "R14_UD8D4_TYPE_1": {
        "representative_class": 4,
        "stable_orders": [1, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4],
        "stable_even_orders": [1, 2, 2, 2],
        "maximal_count": 3,
        "maximal_index": 4,
        "maximal_determinant": 1,
        "maximal_parities": ["odd", "odd", "odd"],
        "classification": "I_(1,12), the unique odd unimodular Lorentzian lattice",
        "sk_list": "List12",
        "sk_entries": [1, 1, 1],
        "sk_symbols": ["I_{12,1}()"] * 3,
        "full_facets": 14,
        "full_finite_vertices": 35,
        "full_cusps": 2,
        "full_coxiter": "graphs/group units qf/12-1,1,1,1,1,1,1,1,1,1,1,1,1.coxiter",
        "full_euler": "691/382588157952000",
        "candidate_euler": "21421/1471492915200",
        "index": 8060,
        "index_status": "exact Solomon Euler-characteristic ratio, identical for all three isometric maximal overclosures",
    },
    "R14_UD8D4_TYPE_2": {
        "representative_class": 6,
        "stable_orders": [1, 2, 2, 2, 4],
        "stable_even_orders": [1, 2],
        "maximal_count": 1,
        "maximal_index": 4,
        "maximal_determinant": -1,
        "maximal_parities": ["odd"],
        "classification": "I_(1,11), the unique odd unimodular Lorentzian lattice",
        "sk_list": "List11",
        "sk_entries": [1],
        "sk_symbols": ["I_{11,1}()"],
        "full_facets": 13,
        "full_finite_vertices": 26,
        "full_cusps": 2,
        "full_coxiter": "graphs/group units qf/11-1,1,1,1,1,1,1,1,1,1,1,1.coxiter",
        "full_euler": "0",
        "candidate_euler": "0",
        "index": None,
        "index_status": "finite containment proved; numerical index is not recovered from the zero Euler characteristics in odd hyperbolic dimension",
    },
    "R14_UD8D4_TYPE_3": {
        "representative_class": 7,
        "stable_orders": [1, 2, 2, 2],
        "stable_even_orders": [1, 2, 2],
        "maximal_count": 3,
        "maximal_index": 2,
        "maximal_determinant": -3,
        "maximal_parities": ["even", "odd", "even"],
        "classification": (
            "two incomparable even maxima are U+E8(-1)+A2(-1), SK List11 entry 4; "
            "the incomparable odd maximum is SK List11 entry 3"
        ),
        "sk_list": "List11",
        "sk_entries": [4, 3, 4],
        "sk_symbols": [
            "II_{11,1}(3^{-1})",
            "I_{11,1}(3^{-1})",
            "II_{11,1}(3^{-1})",
        ],
        "full_facets": [13, 15, 13],
        "full_finite_vertices": [27, 76, 27],
        "full_cusps": [1, 5, 1],
        "full_coxiter": None,
        "full_euler": "0",
        "candidate_euler": "0",
        "index": None,
        "index_status": "finite containment in each of three maximal reflection lattices; numerical indices are not recovered from zero Euler characteristics",
    },
    "R14_UD8D4_TYPE_4": {
        "representative_class": 8,
        "stable_orders": [1, 2],
        "stable_even_orders": [1, 2],
        "maximal_count": 1,
        "maximal_index": 2,
        "maximal_determinant": 2,
        "maximal_parities": ["even"],
        "classification": "U+E8(-1)+A1(-1), the unique even 2-elementary (r,a,delta)=(11,1,1) lattice",
        "sk_list": "List10",
        "sk_entries": [3],
        "sk_symbols": ["II_{10,1}([2^{1}]_{1})"],
        "full_facets": 12,
        "full_finite_vertices": 18,
        "full_cusps": 1,
        "full_coxiter": "graphs/10-thierryP_1-10.coxiter",
        "full_euler": "-1/183936614400",
        "candidate_euler": "-31/5748019200",
        "index": 992,
        "index_status": "exact Solomon Euler-characteristic ratio",
    },
    "R14_UD8D4_TYPE_5": {
        "representative_class": 19,
        "stable_orders": [1, 2],
        "stable_even_orders": [1],
        "maximal_count": 1,
        "maximal_index": 2,
        "maximal_determinant": 3,
        "maximal_parities": ["odd"],
        "classification": "I_(10,1)(3^-1), distinguished from the other determinant-three genus by its discriminant bilinear value 1/3",
        "sk_list": "List10",
        "sk_entries": [4],
        "sk_symbols": ["I_{10,1}(3^{-1})"],
        "full_facets": 14,
        "full_finite_vertices": 56,
        "full_cusps": 3,
        "full_coxiter": "graphs/group units qf/10-3,1,1,1,1,1,1,1,1,1,1.coxiter",
        "full_euler": "-341/16721510400",
        "candidate_euler": "-341/16721510400",
        "index": 1,
        "index_status": "exact normalized relation-graph equality with the complete determinant-three chamber",
    },
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _matrix(payload: Sequence[Sequence[int | str]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(value) for value in row] for row in payload])


def _bilinear_histogram_mod_one(gram: sp.Matrix) -> dict[str, int]:
    moduli, _, discriminant_gram = _discriminant_data(gram)
    histogram: Counter[str] = Counter()
    for element in itertools.product(*(range(modulus) for modulus in moduli)):
        vector = sp.Matrix(element)
        value = sp.factor((vector.T * discriminant_gram * vector)[0]) % 1
        histogram[str(value)] += 1
    return dict(sorted(histogram.items()))


def _sk_gram(list_name: str, entry: int) -> sp.Matrix:
    dimension = int(list_name.removeprefix("List"))
    rows = ast.literal_eval((SK_ROOT / f"Gram{dimension}").read_text())
    return -sp.Matrix(dimension + 1, dimension + 1, rows[entry - 1])


def _sk_reference(
    list_name: str, entry: int, symbol: str, maximal_gram: sp.Matrix
) -> dict[str, Any]:
    list_path = SK_ROOT / list_name
    gram_path = SK_ROOT / list_name.replace("List", "Gram")
    line = list_path.read_text().splitlines()[entry]
    sk_gram = _sk_gram(list_name, entry)
    checks = {
        "list_symbol": symbol in line,
        "rank": sk_gram.rows == maximal_gram.rows,
        "signature": exact_inertia(sk_gram) == exact_inertia(maximal_gram),
        "determinant": sk_gram.det() == maximal_gram.det(),
        "smith": _smith(sk_gram) == _smith(maximal_gram),
        "parity": all(sk_gram[i, i] % 2 == 0 for i in range(sk_gram.rows))
        == all(maximal_gram[i, i] % 2 == 0 for i in range(maximal_gram.rows)),
        "discriminant_bilinear_histogram_mod_one": (
            _bilinear_histogram_mod_one(sk_gram)
            == _bilinear_histogram_mod_one(maximal_gram)
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"SK lattice invariants failed for {list_name}/{entry}: {checks}")
    return {
        "list": list_name,
        "entry_one_based": entry,
        "genus_symbol": symbol,
        "list_row": line.strip(),
        "list_sha256": _sha256(list_path),
        "gram_sha256": _sha256(gram_path),
        "negative_SK_gram": _payload(sk_gram),
        "classification_method": (
            "signature, parity, determinant/Smith and discriminant bilinear form; "
            "uniqueness of the indicated indefinite unimodular, 2-elementary, "
            "or rank-at-least-length-plus-two genus"
        ),
        "checks": checks,
    }


def _coxiter_reference(path_text: str) -> dict[str, Any]:
    path = COXITER_ROOT / path_text
    facets, dimension, relations = _parse_coxiter(path)
    solomon = _euler_from_connected_components(
        relations, _connected_masks_direct(relations)
    )
    return {
        "path": path_text,
        "sha256": _sha256(path),
        "dimension": dimension,
        "facet_count": facets,
        "normalized_relation_matrix": relations,
        "solomon": solomon,
    }


def _one_type(
    type_record: dict[str, Any], high_by_class: dict[int, dict[str, Any]]
) -> dict[str, Any]:
    identifier = type_record["type_id"]
    expected = EXPECTED[identifier]
    high = high_by_class[expected["representative_class"]]
    gram = _matrix(high["fixed_lattice"]["gram"])
    roots = [
        sp.Matrix(row)
        for row in high["fixed_lattice"]["orbit_root_coordinate_rows"]
    ]
    reflections = [_reflection(gram, root) for root in roots]
    audit = _stable_subgroups(gram, reflections)
    groups = audit.pop("groups_internal")
    left_inverse = audit.pop("left_smith_inverse")
    if not (
        audit["stable_integral_subgroup_orders"] == expected["stable_orders"]
        and audit["stable_even_subgroup_orders"] == expected["stable_even_orders"]
        and len(audit["maximal_stable_subgroup_positions_zero_based"])
        == expected["maximal_count"]
    ):
        raise ArithmeticError(f"stable-overlattice denominator changed for {identifier}")

    maxima = []
    for ordinal, position in enumerate(
        audit["maximal_stable_subgroup_positions_zero_based"]
    ):
        subgroup = groups[position]
        cosets = audit["stable_subgroups"][position][
            "cosets_in_smith_coordinates"
        ]
        basis = _overlattice_basis(gram, left_inverse, cosets)
        maximal_gram = sp.simplify(basis.T * gram * basis)
        actions = [sp.simplify(basis.inv() * reflection * basis) for reflection in reflections]
        maximal_roots = [sp.simplify(basis.inv() * root) for root in roots]
        parity = (
            "even"
            if all(maximal_gram[index, index] % 2 == 0 for index in range(maximal_gram.rows))
            else "odd"
        )
        sk_entry = expected["sk_entries"][ordinal]
        sk_symbol = expected["sk_symbols"][ordinal]
        checks = {
            "index": abs(sp.Rational(1, basis.det())) == expected["maximal_index"],
            "determinant_square_law": (
                maximal_gram.det()
                == gram.det() / expected["maximal_index"] ** 2
            ),
            "expected_determinant": (
                maximal_gram.det() == expected["maximal_determinant"]
            ),
            "expected_parity": parity == expected["maximal_parities"][ordinal],
            "roots_integral": all(value.is_Integer for root in maximal_roots for value in root),
            "actions_integral": all(value.is_Integer for action in actions for value in action),
            "actions_unimodular": all(abs(int(action.det())) == 1 for action in actions),
            "actions_preserve_gram": all(
                action.T * maximal_gram * action == maximal_gram for action in actions
            ),
        }
        if not all(checks.values()):
            raise ArithmeticError(f"maximal stable lattice failed for {identifier}: {checks}")
        maxima.append(
            {
                "stable_subgroup_position_zero_based": position,
                "index_over_fixed_lattice": len(subgroup),
                "basis_columns_in_fixed_lattice": _payload(basis),
                "gram": _payload(maximal_gram),
                "rank": maximal_gram.rows,
                "inertia": list(exact_inertia(maximal_gram)),
                "determinant": int(maximal_gram.det()),
                "smith_invariants": _smith(maximal_gram),
                "parity": parity,
                "discriminant_bilinear_histogram_mod_one": (
                    _bilinear_histogram_mod_one(maximal_gram)
                ),
                "discriminant_q_histogram_mod_two": (
                    _q_histogram(maximal_gram) if parity == "even" else None
                ),
                "candidate_root_columns": [_payload(root) for root in maximal_roots],
                "facet_reflection_actions": [_payload(action) for action in actions],
                "facet_reflection_actions_sha256": _canonical_sha256(
                    [_payload(action) for action in actions]
                ),
                "scharlau_kirschmer": _sk_reference(
                    expected["sk_list"], sk_entry, sk_symbol, maximal_gram
                ),
                "checks": checks,
            }
        )

    coxiter = (
        _coxiter_reference(expected["full_coxiter"])
        if expected["full_coxiter"] is not None
        else None
    )
    full_facets = expected["full_facets"]
    full_finite = expected["full_finite_vertices"]
    full_cusps = expected["full_cusps"]
    if coxiter is not None:
        if not (
            coxiter["facet_count"] == full_facets
            and coxiter["dimension"] == type_record["hyperbolic_dimension"]
            and coxiter["solomon"]["rational_euler_characteristic"]
            == expected["full_euler"]
            and coxiter["solomon"]["finite_special_subset_count_by_rank"][-1]
            == full_finite
        ):
            raise ArithmeticError(f"full CoxIter reference changed for {identifier}")

    candidate_euler = Fraction(expected["candidate_euler"])
    full_euler = Fraction(expected["full_euler"])
    if expected["index"] is None:
        index_check = candidate_euler == full_euler == 0
    else:
        index_check = candidate_euler / full_euler == expected["index"]
    exact_full_match = None
    if expected["index"] == 1:
        if coxiter is None:
            raise ArithmeticError("index-one type lacks a CoxIter reference")
        exact_full_match = _relation_permutation(
            type_record["normalized_relation_matrix"],
            coxiter["normalized_relation_matrix"],
        )
        index_check = index_check and exact_full_match is not None
    if not index_check:
        raise ArithmeticError(f"reflection subgroup index failed for {identifier}")

    return {
        "type_id": identifier,
        "hyperbolic_dimension": type_record["hyperbolic_dimension"],
        "facet_count": type_record["facet_count"],
        "source_subgroup_classes_one_based": type_record["class_indices"],
        "original_fixed_lattice": {
            "gram": _payload(gram),
            "rank": gram.rows,
            "inertia": list(exact_inertia(gram)),
            "determinant": int(gram.det()),
            "smith_invariants": _smith(gram),
            "parity": "even",
            "discriminant_bilinear_histogram_mod_one": (
                _bilinear_histogram_mod_one(gram)
            ),
            "discriminant_q_histogram_mod_two": _q_histogram(gram),
        },
        "stable_overlattice_audit": audit,
        "maximal_stable_overlattices": maxima,
        "maximal_overlattice_count": len(maxima),
        "maximal_lattice_classification": expected["classification"],
        "containing_full_reflection_chambers": {
            "scharlau_kirschmer_entries_one_based": expected["sk_entries"],
            "full_facet_counts": full_facets,
            "full_finite_vertex_counts": full_finite,
            "full_cusp_counts": full_cusps,
            "coxiter_reference": coxiter,
        },
        "reflection_subgroup_index": expected["index"],
        "reflection_subgroup_index_status": expected["index_status"],
        "candidate_coxeter_euler_characteristic": expected["candidate_euler"],
        "full_coxeter_euler_characteristic": expected["full_euler"],
        "candidate_to_full_relation_permutation_zero_based": exact_full_match,
        "surface_action_status": (
            "diagram/lattice fixed section only; the general K3 has Aut(X)=C2xC2, "
            "strictly smaller than the C2xS5 weighted-diagram group, and no "
            "specific subgroup lift is asserted"
        ),
        "all_checks_pass": True,
    }


def build_certificate() -> dict[str, Any]:
    source = json.loads(CONSTRUCTION.read_text())
    if not source["all_checks_pass"]:
        raise ArithmeticError("construction certificate is not certified")
    high_by_class = {
        record["class_index"]: record for record in source["high_fixed_sections"]
    }
    records = [
        _one_type(type_record, high_by_class)
        for type_record in source["abstract_diagram_types"]
    ]
    if [record["type_id"] for record in records] != list(EXPECTED):
        raise ArithmeticError("five-type order changed")
    checks = {
        "five_types": len(records) == 5,
        "all_stable_overlattice_denominators_closed": all(
            record["all_checks_pass"] for record in records
        ),
        "exact_index_types": {
            record["type_id"]: record["reflection_subgroup_index"]
            for record in records
            if record["reflection_subgroup_index"] is not None
        }
        == {
            "R14_UD8D4_TYPE_1": 8060,
            "R14_UD8D4_TYPE_4": 992,
            "R14_UD8D4_TYPE_5": 1,
        },
        "odd_dimension_index_open_types": [
            record["type_id"]
            for record in records
            if record["reflection_subgroup_index"] is None
        ]
        == ["R14_UD8D4_TYPE_2", "R14_UD8D4_TYPE_3"],
    }
    if not all(checks.values()):
        raise ArithmeticError(f"global lineage checks failed: {checks}")
    return {
        "schema": SCHEMA,
        "artifact_id": "roulleau_ud8d4_fold_lineage",
        "source": {
            "construction_certificate": str(CONSTRUCTION.relative_to(ROOT)),
            "construction_sha256": _sha256(CONSTRUCTION),
            "scharlau_kirschmer_root": str(SK_ROOT),
            "coxiter_root": str(COXITER_ROOT),
        },
        "types": records,
        "summary": {
            "type_count": 5,
            "unique_maximal_stable_overlattice_type_count": 3,
            "multiple_incomparable_maxima_type_count": 2,
            "exact_numerical_index_count": 3,
            "finite_index_value_open_count": 2,
            "indices": {
                record["type_id"]: record["reflection_subgroup_index"]
                for record in records
            },
            "abstract_lineage": (
                "TYPE_5 is an exact known full chamber; TYPE_1--TYPE_4 are "
                "proper finite-index reflection subgroups of known maximal "
                "arithmetic chambers, while their literal diagrams remain "
                "unmatched in the pinned CoxIter corpus"
            ),
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    certificate = build_certificate()
    if args.check_archive:
        if certificate != json.loads(OUTPUT.read_text()):
            raise SystemExit("recomputed certificate differs from archive")
    elif args.output:
        args.output.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    else:
        print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
