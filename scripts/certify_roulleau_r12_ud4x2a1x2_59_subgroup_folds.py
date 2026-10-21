#!/usr/bin/env python3
"""Certify every high fixed section of Roulleau's 59-wall rank-12 stock.

The complete nef chamber belongs to
``NS(X)=U+D4(-1)^2+A1(-1)^2``.  Roulleau prints an integral basis, the
first nineteen curves, and closed formulae for the remaining forty curves.
The companion scout reconstructs those data exactly.  This certificate
enumerates all subgroup conjugacy classes of the full weighted wall group,
applies the exact Coxeter and Vinberg gates, separates the central p67 diagram
action from both general-surface deck actions, and computes the full stable-
overlattice and reflection-group lineage of every survivor.  Semantic action
IDs are primary: exploratory GAP classes C3/C9 are the same two classes called
C4/C12 by the pinned small-generator presentation used here.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Any, Sequence

import sympy as sp

from scripts.certify_borcherds_char2_high_fixed_sections_vinberg import (
    _parse_output,
)
from scripts.certify_jacobian_kummer_h10_reflection_index import (
    _connected_masks_by_nonright_components,
    _connected_masks_direct,
    _euler_from_connected_components,
    _relation_matrix,
)
from scripts.certify_roulleau_non2elementary_involution_folds import (
    _scan_coxiter,
)
from scripts.certify_shimada_disc16_d5_h10_lineage import (
    _p10_roots,
    _parse_sk_list10_first,
)
from scripts.scout_roulleau_r12_ud4x2a1x2_59 import coordinate_model
from scripts.scout_roulleau_rank12_graph_folds import (
    _fixed_lattice_data,
    _lift_actions,
)
from scripts.search_k3_high_rank_involution_folds import (
    _coxeter_audit,
    _fold_reflection_audit,
    _group_checks,
    _matrix_payload,
)
from scripts.search_k3_nef_subgroup_folds import (
    Permutation,
    Subgroup,
    _compose,
    _generated_subgroup,
    _identity,
    _orbit_gram,
    _orbits,
    _order,
    _pairwise_disjoint,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/roulleau_r12_ud4x2a1x2_59_all_subgroup_folds.json"
SCOUT = ROOT / "data/enumerations/roulleau_r12_ud4x2a1x2_59_scout.json"
INVENTORY = ROOT / "data/enumerations/roulleau_rank11_14_source_stock_inventory.json"
SURFACE_ACTIONS = ROOT / "data/enumerations/roulleau_r12_ud4x2a1x2_surface_actions.json"
SOURCE_TEX = Path("/private/tmp/roulleau_source/epiga_vol6_Roulleau.tex")
CPP_SOURCE = ROOT / "scripts/certify_borcherds_char2_high_fixed_sections_vinberg.cpp"
SK_LIST10 = Path("/private/tmp/Scharlau-Kirschmer-Lorentz/List10")

SCHEMA = "coxeter.roulleau_r12_ud4x2a1x2_59.all_subgroup_folds.v1"
RANK = 12
WALL_COUNT = 59
SOURCE_TEX_SHA256 = "e6eeaaa4a3f375f209762718e8ef9a361b0a95356cc8be2e9f99a6713abf5df3"
SCOUT_SHA256 = "47acc3aeb02c0b95669caaffcb1ec46dece31dd99261c09b554ec7c6b289e23b"
CPP_SHA256 = "45174f5461d28affc1fdb4cb55047b7a1581037d307e32800651badb3c2df045"
SK_LIST10_SHA256 = "9e627c4e1303015fd96acdcf6ccbeb264ebadba49abfcd761705f54d211b2e21"
SURFACE_ACTIONS_SHA256 = "c5ce379b73f544d1204d7bcdffff03545d376478a6a7a7be9e653d92f555d154"


def _zero_based(rows: Sequence[Sequence[int]]) -> tuple[Permutation, ...]:
    return tuple(tuple(value - 1 for value in row) for row in rows)


# GAP's SmallGeneratingSet for the complete order-2880 wall group.  Pinning
# these avoids the much slower generic Python greedy generator search.
SMALL_GENERATORS = _zero_based(
    (
        (13, 2, 1, 3, 4, 7, 6, 8, 9, 10, 11, 12, 5, 19, 15, 14, 16, 17, 18,
         36, 33, 38, 39, 47, 42, 41, 45, 44, 40, 26, 23, 28, 29, 57, 52, 51,
         55, 54, 50, 56, 53, 58, 59, 27, 22, 21, 25, 24, 20, 46, 43, 48, 49,
         37, 32, 31, 35, 34, 30),
        (17, 16, 19, 14, 15, 7, 6, 8, 10, 9, 12, 11, 18, 4, 3, 13, 1, 2, 5,
         30, 49, 34, 31, 48, 35, 32, 43, 46, 37, 20, 59, 24, 21, 58, 25, 22,
         53, 56, 27, 50, 29, 54, 51, 28, 55, 52, 23, 26, 57, 40, 39, 44, 41,
         38, 45, 42, 33, 36, 47),
        (18, 19, 15, 16, 17, 7, 6, 8, 10, 9, 12, 11, 14, 5, 13, 2, 3, 4, 1,
         43, 37, 35, 34, 40, 41, 42, 39, 38, 36, 53, 27, 25, 24, 50, 51, 52,
         29, 28, 26, 23, 57, 55, 54, 20, 21, 22, 59, 58, 56, 33, 47, 45, 44,
         30, 31, 32, 49, 48, 46),
    )
)


# The central p67 formula/diagram symmetry; the surface-action audit proves
# that it is neither deck action on the general K3.
P67: Permutation = tuple(
    (
        0, 1, 2, 3, 4, 6, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
        29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 19, 20, 21, 22, 23, 24, 25,
        26, 27, 28, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 39, 40, 41, 42,
        43, 44, 45, 46, 47, 48,
    )
)


EXPECTED = {
    4: {
        "id": "R12-D4x2A1x2-p67-H10-38",
        "order": 2,
        "class_size": 1,
        "facets": 38,
        "finite_counts": [38, 367, 1820, 5654, 11934, 17572, 17850, 11900, 4518, 526],
        "cusps": 141,
        "affine": 307,
        "endpoints": {"0,2": 920, "1,1": 1936, "2,0": 1662},
        "relation_sha256": "1cf3bda953ca40f8093c435dce74cd71b5692e34316fe3f5f4f9ed83e3a0d30a",
        "connected_elliptic": 465,
        "euler": "-527/2580480",
        "index": 2276640,
        "stable_count": 122,
        "maximal_count": 15,
        "stable_action_sha256": "b6cf11bbb80d140ba78b864de2acb7c3c03af992de8fc38dffee1f24d5dbdd06",
        "q_histogram": {"0": 16, "1": 16, "3/4": 12, "7/4": 20},
        "surface_status": "diagram_and_integral_NS_isometry_general_surface_excluded_special_period_open",
        "exploratory_all_element_GAP_class_index": 3,
    },
    12: {
        "id": "R12-D4x2A1x2-S6transposition-H10-45",
        "order": 2,
        "class_size": 15,
        "facets": 45,
        "finite_counts": [45, 507, 2793, 9201, 19804, 28971, 28982, 19136, 7332, 1048],
        "cusps": 93,
        "affine": 302,
        "endpoints": {"0,2": 578, "1,1": 3028, "2,0": 3726},
        "relation_sha256": "632685c1d62444d2a9d949977cb31acb48dddbd1d4a9264954d896d1e9be0558",
        "connected_elliptic": 332,
        "euler": "-527/2903040",
        "index": 2023680,
        "stable_count": 38,
        "maximal_count": 3,
        "stable_action_sha256": "7660b69da8201bddea3053b5e317f62a9ba8c40cee81653d70d1a0393bcb1107",
        "q_histogram": {
            "0": 8, "1": 8, "1/2": 8, "1/4": 8, "3/2": 8,
            "3/4": 4, "5/4": 8, "7/4": 12,
        },
        "surface_status": "diagram_and_integral_NS_isometry_general_surface_excluded_special_period_open",
        "exploratory_all_element_GAP_class_index": 9,
    },
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _canonical_sha256(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _matrix(payload: Sequence[Sequence[str | int]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(value) for value in row] for row in payload])


def _gap_permutation(permutation: Permutation) -> str:
    return "PermList(" + str([value + 1 for value in permutation]) + ")"


def _gap_census(group: tuple[Permutation, ...]) -> dict[str, Any]:
    generators = ",".join(_gap_permutation(row) for row in SMALL_GENERATORS)
    program = f"""
SizeScreen([1000000,1000000]);;
G:=Group([{generators}]);;
cc:=ConjugacyClassesSubgroups(G);;
Print("GROUP|",GAPInfo.Version,"|",Size(G),"|",Length(cc),"|",
      StructureDescription(G),"\\n");;
i:=0;;
for c in cc do
  i:=i+1;; h:=Representative(c);;
  Print("CLASS|",i,"|",Size(h),"|",Size(c),"|",StructureDescription(h),"|",
        List(GeneratorsOfGroup(h),x->ListPerm(x,{WALL_COUNT})),"\\n");;
od;
QUIT;
"""
    completed = subprocess.run(
        ["gap", "-q"], input=program, text=True, capture_output=True,
        timeout=300, check=False,
    )
    if completed.returncode != 0 or completed.stderr.strip():
        raise RuntimeError(
            "GAP subgroup census failed:\n"
            + completed.stdout[-4000:] + "\n" + completed.stderr[-4000:]
        )
    group_set = set(group)
    header = None
    rows = []
    for line in completed.stdout.splitlines():
        if line.startswith("GROUP|"):
            _, version, order, count, structure = line.split("|")
            header = {
                "gap_version": version,
                "order": int(order),
                "subgroup_conjugacy_class_count": int(count),
                "structure": structure,
            }
        elif line.startswith("CLASS|"):
            _, index, order, class_size, structure, raw = line.split("|", 5)
            representative_generators = tuple(
                tuple(int(value) - 1 for value in row)
                for row in ast.literal_eval(raw)
            )
            if not set(representative_generators).issubset(group_set):
                raise ArithmeticError("GAP returned a generator outside the wall group")
            rows.append(
                {
                    "class_index": int(index),
                    "order": int(order),
                    "conjugacy_class_size": int(class_size),
                    "structure": structure,
                    "representative_generators": representative_generators,
                }
            )
    if header is None or not (
        header == {
            "gap_version": "4.15.1",
            "order": 2880,
            "subgroup_conjugacy_class_count": 914,
            "structure": "C2 x C2 x S6",
        }
        and len(rows) == 914
    ):
        raise ArithmeticError(f"GAP denominator changed: {header}")
    return {"header": header, "rows": rows}


def _action(permutation: Permutation, wall_rows: sp.Matrix) -> sp.Matrix:
    return wall_rows[[permutation[index] for index in range(RANK)], :]


def _fixed_rank(
    generators: Sequence[Permutation], wall_rows: sp.Matrix
) -> int:
    equations = sp.Matrix.vstack(
        *(_action(generator, wall_rows).T - sp.eye(RANK) for generator in generators)
    )
    return RANK - equations.rank()


def _inverse(permutation: Permutation) -> Permutation:
    result = [0] * len(permutation)
    for source, target in enumerate(permutation):
        result[target] = source
    return tuple(result)


def _conjugate(subgroup: Subgroup, action: Permutation) -> Subgroup:
    inverse = _inverse(action)
    return frozenset(
        _compose(_compose(action, element), inverse) for element in subgroup
    )


def _cycles(permutation: Permutation, labels: Sequence[str]) -> list[list[str]]:
    seen = set()
    result = []
    for start in range(len(permutation)):
        if start in seen or permutation[start] == start:
            continue
        cycle = []
        current = start
        while current not in seen:
            seen.add(current)
            cycle.append(labels[current])
            current = permutation[current]
        result.append(cycle)
    return result


def _source_audit(labels: Sequence[str]) -> dict[str, Any]:
    inventory = json.loads(INVENTORY.read_text())
    stock = next(row for row in inventory["stocks"] if row["key"] == "R12_UD4x2A1x2")
    source = SOURCE_TEX.read_text()
    section = "\n".join(source.splitlines()[8803:9033])
    cycles = _cycles(P67, labels)
    expected_pairs = (
        [["A6", "A7"]]
        + [[f"B1_{i}{j}", f"B2_{i}{j}"] for i, j in
           ((1,2),(1,3),(1,4),(1,5),(2,3),(2,4),(2,5),(3,4),(3,5),(4,5))]
        + [[f"C1_{i}{j}", f"C2_{i}{j}"] for i, j in
           ((1,2),(1,3),(1,4),(1,5),(2,3),(2,4),(2,5),(3,4),(3,5),(4,5))]
    )
    checks = {
        "source_tex_pinned": _sha256(SOURCE_TEX) == SOURCE_TEX_SHA256,
        "inventory_certified": inventory["all_checks_pass"] is True,
        "inventory_stock": (
            stock["picard_rank"] == 12 and stock["wall_count"] == 59
            and stock["general_surface_automorphism_group"] == "C2xC2"
            and stock["roulleau_ns_action_tag"] == "tn"
            and stock["source_completeness_method"]
            == "induced_chamber_completion_inside_II_1_17"
        ),
        "section_has_59_curves": "contains $59$ $(-2)$-curves" in section,
        "section_has_integral_basis": "A_{1},\\dots,A_{12}$ generate" in section,
        "section_has_both_models": (
            r"\subsubsection{First involution}" in section
            and r"\subsubsection{Second involution}" in section
            and r"\varphi'\colon X\to\mathbf{F}_{4}" in section
        ),
        "second_model_is_hyperelliptic": "is hyperelliptic since $D_{8}F'=2$" in section,
        "source_action_criterion_present": (
            "singularities  of type" in source
            and r"\mathbf{d_{4}},\,\mathbf{d_{6}},\,\mathbf{d_{8}}" in source
            and "the action on $\\cu$-curves" in source
            and "is trivial" in source
        ),
        "closed_BC_formulae_present": all(
            token in section
            for token in (r"B_{ij}^{(1)}", r"B_{ij}^{(2)}", r"C_{ij}^{(1)}", r"C_{ij}^{(2)}")
        ),
        "p67_has_exact_printed_superscript_cycles": cycles == expected_pairs,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"Roulleau source/geometric audit changed: {checks}")
    return {
        "tex_path": str(SOURCE_TEX),
        "tex_sha256": SOURCE_TEX_SHA256,
        "section_lines_one_based": [8804, 9033],
        "table_lines_one_based": [10926, 10926],
        "completeness_method": stock["source_completeness_method"],
        "general_surface_automorphism_group": "C2xC2",
        "neron_severi_action_tag": "tn",
        "p67_moved_wall_pairs": cycles,
        "geometric_interpretation": (
            "The p67 formula symmetry swaps A6/A7 and every B1/B2 and C1/C2 "
            "pair, but it is not either deck action on the general surface.  The "
            "second D8 model has only d4/a1 branch singularities and is NS-trivial "
            "by Roulleau's criterion; the first D2 model is the nontrivial NS action."
        ),
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def _strict_facet_witnesses(
    subgroup: Subgroup,
    actions: dict[Permutation, sp.Matrix],
    fixed_lattice: dict[str, Any],
) -> dict[str, Any]:
    fixed_basis = _matrix(fixed_lattice["basis_rows_in_source_basis"])
    fixed_gram = _matrix(fixed_lattice["gram"])
    roots = _matrix(fixed_lattice["orbit_root_coordinate_columns"])
    d40 = sp.Matrix([[3, 3, 3, 0, 0, 2, 2, 6, 6, 3, 1, 7]])
    invariant = sum((d40 * actions[element] for element in subgroup), sp.zeros(1, RANK))
    ample = sp.Matrix(tuple(next(iter(sp.linsolve((fixed_basis.T, invariant.T))))))
    rows = []
    for index in range(roots.cols):
        root = roots[:, index]
        norm = (root.T * fixed_gram * root)[0]
        pairing = (ample.T * fixed_gram * root)[0]
        witness = ample - pairing / norm * root
        pairings = list(witness.T * fixed_gram * roots)
        square = (witness.T * fixed_gram * witness)[0]
        other = [value for position, value in enumerate(pairings) if position != index]
        rows.append(
            {
                "facet_one_based": index + 1,
                "square": str(square),
                "target_pairing": str(pairings[index]),
                "minimum_other_pairing": str(min(other)),
                "strict_timelike_facet_witness": (
                    square > 0 and pairings[index] == 0 and all(value > 0 for value in other)
                ),
            }
        )
    return {
        "method": "orthogonal projections of the subgroup-averaged ample D40 class",
        "records": rows,
        "all_orbit_inequalities_are_genuine_facets": all(
            row["strict_timelike_facet_witness"] for row in rows
        ),
    }


def _vinberg(executable: Path, positive_gram: sp.Matrix, expected: dict[str, Any]) -> dict[str, Any]:
    lines = [f"10 {positive_gram.rows}"]
    lines.extend(
        " ".join(str(int(positive_gram[left, right])) for right in range(positive_gram.cols))
        for left in range(positive_gram.rows)
    )
    completed = subprocess.run(
        [str(executable)], input="\n".join(lines) + "\n", text=True,
        capture_output=True, timeout=1800, check=True,
    )
    result = _parse_output(completed.stdout)
    counts = result["finite_special_subset_counts_by_rank"]
    endpoints = result["edge_endpoint_type_histogram"]
    checks = {
        "dimension_and_facets": result["dimension"] == 10 and result["facet_count"] == expected["facets"],
        "finite_special_census": (
            counts == expected["finite_counts"]
            and result["all_finite_subset_search_node_count"] == 1 + sum(counts)
        ),
        "cusps_and_affine_components": (
            result["ideal_vertex_count"] == expected["cusps"]
            and result["connected_affine_component_count"] == expected["affine"]
        ),
        "literal_vinberg_failures_zero": (
            result["bad_edge_count"] == 0
            and result["bad_parabolic_neighborhood_count"] == 0
        ),
        "endpoint_histogram": endpoints == expected["endpoints"],
        "two_endpoints_per_edge": (
            sum((int(key[0]) + int(key[2])) * count for key, count in endpoints.items())
            == 2 * counts[-2]
        ),
        "finite_vertex_edge_incidence": (
            2 * endpoints.get("2,0", 0) + endpoints.get("1,1", 0)
            == 10 * counts[-1]
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"literal Vinberg audit changed: {checks}")
    return result | {
        "elliptic_edge_count": counts[-2],
        "finite_vertex_count": counts[-1],
        "finite_volume": True,
        "compact": False,
        "classification": "finite_volume_noncompact_hyperbolic_Coxeter_polytope",
        "criterion": (
            "every elliptic codimension-one edge has two finite/ideal endpoints, "
            "and every maximal affine component has a valid full zero-neighborhood"
        ),
        "checks": checks,
        "all_checks_pass": True,
    }


def _p10_lineage() -> dict[str, Any]:
    roots = _p10_roots()
    odd_unimodular = sp.diag(1, *([-1] * 10))
    positive_gram = -roots.T * odd_unimodular * roots
    payload = [[str(positive_gram[i, j]) for j in range(12)] for i in range(12)]
    relations = _relation_matrix(payload)
    direct = _connected_masks_direct(relations)
    componentwise = _connected_masks_by_nonright_components(relations)
    solomon = _euler_from_connected_components(relations, componentwise)
    sk = _parse_sk_list10_first(SK_LIST10)
    checks = {
        "sk_source_pinned": _sha256(SK_LIST10) == SK_LIST10_SHA256,
        "sk_entry_one": sk == {
            "determinant": 1,
            "genus_symbol": "I_{10,1}()",
            "finite_face_counts_facets_to_vertices": [12,65,210,450,672,714,540,285,99,17],
            "cusp_count": 2,
        },
        "componentwise_helper_equals_direct": componentwise == direct,
        "p10_special_counts": solomon["finite_special_subset_count_by_rank"]
        == [1,12,65,210,450,672,714,540,285,99,17],
        "p10_euler": solomon["rational_euler_characteristic"] == "-1/11147673600",
    }
    if not all(checks.values()):
        raise ArithmeticError(f"P10 denominator changed: {checks}")
    return {
        "name": "Vinberg P10",
        "lattice": "I_(1,10), odd unimodular",
        "scharlau_kirschmer": sk,
        "simple_root_columns": _matrix_payload(roots),
        "positive_normal_gram": _matrix_payload(positive_gram),
        "solomon": solomon,
        "checks": checks,
        "all_checks_pass": True,
    }


def _root_type_histogram(fixed_lattice: dict[str, Any]) -> dict[str, int]:
    return dict(
        sorted(
            Counter(
                f"{row['norm']},{row['divisibility']},{len(row['source_walls_one_based'])}"
                for row in fixed_lattice["reflection_checks"]
            ).items()
        )
    )


def build_certificate() -> dict[str, Any]:
    if not (
        _sha256(SCOUT) == SCOUT_SHA256
        and _sha256(CPP_SOURCE) == CPP_SHA256
        and _sha256(SURFACE_ACTIONS) == SURFACE_ACTIONS_SHA256
    ):
        raise ArithmeticError("a pinned construction source changed")
    scout = json.loads(SCOUT.read_text())
    surface_actions = json.loads(SURFACE_ACTIONS.read_text())
    if not (
        surface_actions["all_checks_pass"] is True
        and surface_actions["checks"]["p67_is_neither_geometric_NS_action"] is True
        and surface_actions["first_double_plane_deck_action"]["fixed_rank"] == 4
        and surface_actions["second_F4_deck_action"]["fixed_rank"] == 12
        and surface_actions["second_F4_deck_action"]["status"]
        == "actual_general_surface_automorphism_trivial_on_NS"
    ):
        raise ArithmeticError("independent surface-action certificate changed")
    basis_gram, wall_records, wall_gram = coordinate_model()
    wall_rows = sp.Matrix.vstack(*(row["coordinates"] for row in wall_records))
    labels = [str(row["name"]) for row in wall_records]
    group = tuple(
        tuple(row) for row in scout["weighted_group"]["elements_zero_based"]
    )
    source = _source_audit(labels)
    reconstruction_checks = {
        "scout_checks": all(scout["coordinate_checks"].values()),
        "basis_determinant_and_inertia": (
            int(basis_gram.det()) == -64 and exact_inertia(basis_gram) == (1, 11, 0)
        ),
        "wall_rank_and_inertia": (
            wall_gram.rank() == 12 and exact_inertia(wall_gram) == (1, 11, 47)
        ),
        "all_norm_minus_two": all(wall_gram[index, index] == -2 for index in range(59)),
        "scout_matrices_replay": (
            _matrix(scout["basis_gram"]) == basis_gram
            and _matrix(scout["wall_gram"]) == wall_gram
        ),
        "source_formulae_replay": all(
            wall_rows[index, column] == sp.Rational(scout["walls"][index]["coordinates"][column])
            for index in range(59) for column in range(12)
        ),
    }
    if not all(reconstruction_checks.values()):
        raise ArithmeticError(f"coordinate reconstruction changed: {reconstruction_checks}")

    group_set = set(group)
    generator_closure_check = all(
        _compose(element, generator) in group_set
        for element in group for generator in SMALL_GENERATORS
    )
    if not (
        len(group) == len(group_set) == 2880
        and all(scout["weighted_group"]["checks"].values())
        and set(SMALL_GENERATORS).issubset(group_set)
        and generator_closure_check
    ):
        raise ArithmeticError("weighted wall group archive changed")
    census = _gap_census(group)

    center = tuple(
        element for element in group
        if all(_compose(element, generator) == _compose(generator, element)
               for generator in SMALL_GENERATORS)
    )
    action_generators = set(SMALL_GENERATORS) | set(center)
    for row in census["rows"]:
        if row["class_index"] in EXPECTED:
            action_generators.update(row["representative_generators"])
    actions = _lift_actions(tuple(action_generators), wall_rows, basis_gram)
    center_ranks = sorted(
        RANK - (actions[element] - sp.eye(RANK)).rank() for element in center
    )
    label_position = {label: index for index, label in enumerate(labels)}
    d2_deck = next(
        element for element in center
        if element[label_position["B1_12"]] == label_position["C1_12"]
        and element[label_position["B2_12"]] == label_position["C2_12"]
    )
    d2_deck_fixed_rank = RANK - (actions[d2_deck] - sp.eye(RANK)).rank()
    group_audit = {
        "order": len(group),
        "structure": census["header"]["structure"],
        "small_generators_zero_based": [list(row) for row in SMALL_GENERATORS],
        "element_order_histogram": {
            str(key): value for key, value in sorted(Counter(_order(row) for row in group).items())
        },
        "center_order": len(center),
        "center_fixed_ranks": center_ranks,
        "p67_is_central": P67 in center,
        "p67_fixed_rank": RANK - (actions[P67] - sp.eye(RANK)).rank(),
        "first_D2_deck_action_zero_based": list(d2_deck),
        "first_D2_deck_fixed_rank": d2_deck_fixed_rank,
        "second_D8_deck_NS_action": "identity",
        "general_surface_NS_action_image_order": 2,
        "checks": {
            "scout_group_checks": all(scout["weighted_group"]["checks"].values()),
            "closure_under_small_generators": generator_closure_check,
            "gap_structure": census["header"]["structure"] == "C2 x C2 x S6",
            "center_is_V4": len(center) == 4 and center_ranks == [4, 5, 11, 12],
            "p67_exact": P67 in center and RANK - (actions[P67] - sp.eye(RANK)).rank() == 11,
            "deck_actions_distinguished": (
                d2_deck_fixed_rank == 4
                and d2_deck != P67
                and _identity(WALL_COUNT) in center
            ),
            "small_generators_lift_integrally": len(actions) >= len(SMALL_GENERATORS),
        },
    }
    group_audit["all_checks_pass"] = all(group_audit["checks"].values())
    if not group_audit["all_checks_pass"]:
        raise ArithmeticError(f"group recognition changed: {group_audit}")

    public_classes = []
    high_internal = []
    for row in census["rows"]:
        generators = row["representative_generators"]
        fixed_rank = _fixed_rank(generators, wall_rows)
        proper_high = row["order"] > 1 and fixed_rank >= 11
        public = {
            **{key: value for key, value in row.items() if key != "representative_generators"},
            "fixed_rank": fixed_rank,
            "hyperbolic_dimension": fixed_rank - 1,
            "proper_fixed_section_above_H9": proper_high,
            "representative_generators_zero_based": [list(value) for value in generators],
        }
        public_classes.append(public)
        if proper_high:
            subgroup = _generated_subgroup(generators)
            if len(subgroup) != row["order"]:
                raise ArithmeticError("high GAP subgroup failed direct generation")
            conjugates = {_conjugate(subgroup, action) for action in group}
            if len(conjugates) != row["conjugacy_class_size"]:
                raise ArithmeticError("high GAP conjugacy class failed direct replay")
            orbits = _orbits(subgroup, WALL_COUNT)
            orbit_gram = _orbit_gram(wall_gram, orbits)
            if orbit_gram.rank() != fixed_rank:
                raise ArithmeticError("action and orbit fixed ranks disagree")
            high_internal.append((public, subgroup, generators, orbits, orbit_gram))

    actual_high = {
        row[0]["class_index"]: {
            "order": row[0]["order"],
            "class_size": row[0]["conjugacy_class_size"],
            "fixed_rank": row[0]["fixed_rank"],
        }
        for row in high_internal
    }
    expected_high = {
        index: {"order": expected["order"], "class_size": expected["class_size"], "fixed_rank": 11}
        for index, expected in EXPECTED.items()
    }
    if actual_high != expected_high:
        raise ArithmeticError(f"high-rank subgroup partition changed: {actual_high}")

    compiler = shutil.which("clang++") or shutil.which("g++")
    if compiler is None:
        raise RuntimeError("a C++17 compiler is required")
    p10 = _p10_lineage()
    p10_euler = Fraction(p10["solomon"]["rational_euler_characteristic"])
    cases = []
    with tempfile.TemporaryDirectory(prefix="roulleau-r12-59-vinberg-") as raw:
        executable = Path(raw) / "vinberg"
        subprocess.run(
            [compiler, "-O3", "-std=c++17", "-Wall", "-Wextra", "-Werror",
             str(CPP_SOURCE), "-o", str(executable)],
            check=True, capture_output=True, text=True,
        )
        for public, subgroup, generators, orbits, orbit_gram in high_internal:
            class_index = public["class_index"]
            expected = EXPECTED[class_index]
            cycles_one = [[value + 1 for value in orbit] for orbit in orbits]
            coxeter = _coxeter_audit(
                orbit_gram, list(range(orbit_gram.rows)), cycles_one
            )
            reflection = _fold_reflection_audit(
                wall_gram, orbit_gram, list(range(orbit_gram.rows)), cycles_one
            )
            exact_gate = {
                "all_orbits_pairwise_disjoint": all(
                    _pairwise_disjoint(wall_gram, orbit) for orbit in orbits
                ),
                "all_orbit_normals_spacelike": all(
                    orbit_gram[index, index] < 0 for index in range(orbit_gram.rows)
                ),
                "fixed_space_lorentzian": exact_inertia(orbit_gram)
                == (1, 10, orbit_gram.rows - 11),
                "all_pairs_exact_coxeter": coxeter["all_pairs_are_exact_coxeter_relations"],
                "all_folded_reflections_integral": reflection["all_folded_reflections_are_integral"],
            }
            if not all(exact_gate.values()):
                raise ArithmeticError(f"class {class_index} failed exact Coxeter gate: {exact_gate}")

            fixed_lattice = _fixed_lattice_data(
                generators, actions, basis_gram, wall_rows, orbits, orbit_gram
            )
            witnesses = _strict_facet_witnesses(subgroup, actions, fixed_lattice)
            positive_gram = -orbit_gram
            relations = _relation_matrix(
                [[str(positive_gram[i, j]) for j in range(positive_gram.cols)]
                 for i in range(positive_gram.rows)]
            )
            connected = _connected_masks_by_nonright_components(relations)
            solomon = _euler_from_connected_components(relations, connected)
            section_euler = Fraction(solomon["rational_euler_characteristic"])
            index = section_euler / p10_euler
            vinberg = _vinberg(executable, positive_gram, expected)
            coxiter = _scan_coxiter(relations, 10, orbit_gram.rows)
            stable = fixed_lattice["stable_overlattice_audit"]
            maxima = fixed_lattice["maximal_stable_overlattices"]
            relation_sha = _canonical_sha256(relations)
            geometric = {
                "status": expected["surface_status"],
                "is_general_surface_automorphism_action": False,
                "p67_exact_permutation_match": generators == (P67,),
                "central_in_full_wall_group": all(
                    generator in center for generator in generators
                ),
                "general_surface_NS_action_image": {
                    "order": 2,
                    "identity_from_second_D8_deck": True,
                    "nontrivial_generator": "first D2 deck action",
                    "nontrivial_generator_fixed_rank": d2_deck_fixed_rank,
                },
                "p67_excluded_from_general_surface_image": (
                    P67 != d2_deck and P67 != _identity(WALL_COUNT)
                ),
                "selected_action_excluded_from_general_surface_image": all(
                    generator not in {_identity(WALL_COUNT), d2_deck}
                    for generator in generators
                    if generator != _identity(WALL_COUNT)
                ),
                "qualification": (
                    "The fixed-nef section is an exact integral diagram/NS-lattice "
                    "construction.  It is not induced by an automorphism of the "
                    "general K3; realization at a special period remains open."
                ),
            }
            checks = {
                "coxeter_gate": all(exact_gate.values()),
                "genuine_facets": witnesses["all_orbit_inequalities_are_genuine_facets"],
                "finite_volume": vinberg["all_checks_pass"],
                "fixed_lattice": (
                    fixed_lattice["rank"] == 11
                    and fixed_lattice["determinant"] == 64
                    and fixed_lattice["smith_invariants"] == [1] * 6 + [2] * 4 + [4]
                    and fixed_lattice["even"]
                    and fixed_lattice["discriminant_q_histogram_mod_two"] == expected["q_histogram"]
                    and fixed_lattice["all_reflections_integral"]
                ),
                "stable_denominator": (
                    len(stable["stable_subgroups"]) == expected["stable_count"]
                    and stable["induced_discriminant_actions_sha256"]
                    == expected["stable_action_sha256"]
                    and len(maxima) == expected["maximal_count"]
                    and all(
                        maximum["index_over_fixed_lattice"] == 8
                        and maximum["determinant"] == 1
                        and maximum["smith_invariants"] == [1] * 11
                        and maximum["parity"] == "odd"
                        and maximum["all_reflections_integral"]
                        and maximum["all_reflections_preserve_gram"]
                        for maximum in maxima
                    )
                ),
                "solomon": (
                    len(connected) == expected["connected_elliptic"]
                    and solomon["finite_special_subset_count_by_rank"]
                    == [1] + expected["finite_counts"]
                    and solomon["rational_euler_characteristic"] == expected["euler"]
                    and index == expected["index"]
                ),
                "relation_digest": relation_sha == expected["relation_sha256"],
                "bounded_corpora_nonmatch": (
                    coxiter["ordinary_header_candidate_count"] == 0
                    and coxiter["doubling_header_candidate_count"] == 0
                    and not coxiter["ordinary_exact_matches"]
                    and not coxiter["doubling_exact_matches"]
                ),
                "geometric_status": (
                    not geometric["is_general_surface_automorphism_action"]
                    and geometric["p67_excluded_from_general_surface_image"]
                    and geometric["selected_action_excluded_from_general_surface_image"]
                    and (
                        (class_index == 4 and geometric["p67_exact_permutation_match"]
                         and geometric["central_in_full_wall_group"])
                        or (class_index == 12 and not geometric["central_in_full_wall_group"])
                    )
                ),
            }
            if not all(checks.values()):
                raise ArithmeticError(
                    f"class {class_index} certificate changed: {checks}; "
                    f"actual_stable_count={len(stable['stable_subgroups'])}; "
                    f"actual_stable_sha={stable['induced_discriminant_actions_sha256']}; "
                    f"actual_maxima={len(maxima)}; actual_relation_sha={relation_sha}"
                )
            cases.append(
                {
                    "id": expected["id"],
                    "class_index": class_index,
                    "exploratory_all_element_GAP_class_index": expected[
                        "exploratory_all_element_GAP_class_index"
                    ],
                    "subgroup": {
                        "order": len(subgroup),
                        "conjugacy_class_size": public["conjugacy_class_size"],
                        "representative_generators_zero_based": [list(row) for row in generators],
                        "wall_orbits_one_based": cycles_one,
                    },
                    "geometric_realization": geometric,
                    "exact_coxeter_gate": exact_gate,
                    "coxeter_audit": coxeter,
                    "reflection_audit": reflection,
                    "strict_facet_witnesses": witnesses,
                    "negative_orbit_gram": _matrix_payload(orbit_gram),
                    "positive_normal_gram": _matrix_payload(positive_gram),
                    "normalized_relation_matrix": relations,
                    "relation_matrix_sha256": relation_sha,
                    "relation_histogram": dict(
                        sorted(Counter(
                            relations[left][right]
                            for left in range(len(relations))
                            for right in range(left + 1, len(relations))
                        ).items())
                    ),
                    "fixed_lattice": fixed_lattice,
                    "root_type_histogram": _root_type_histogram(fixed_lattice),
                    "vinberg_finite_volume_audit": vinberg,
                    "solomon": solomon,
                    "maximal_reflection_lineage": {
                        "containing_lattice": "I_(1,10), odd unimodular",
                        "containing_full_chamber": "Vinberg P10 / Scharlau-Kirschmer List10 entry 1",
                        "maximal_stable_embedding_count": len(maxima),
                        "index_in_each_P10_reflection_group": int(index),
                        "index_formula": "chi(section chamber) / chi(P10)",
                        "section_euler_characteristic": str(section_euler),
                        "P10_euler_characteristic": str(p10_euler),
                        "arithmetic_status": "known odd-unimodular arithmetic class",
                    },
                    "bounded_lineage": {
                        "coxiter_scan": coxiter,
                        "ordinary_H10_facet_sizes": [12,13,14,15,16,19,20,28,44],
                        "doubling_H10_facet_sizes": [15,16,18,22,30,46],
                        "prior_project_same_dimension_facet_count": 0,
                        "literal_priority_status": (
                            "unmatched in the pinned CoxIter snapshots and pre-catch project ledger; "
                            "worldwide literal priority remains open"
                        ),
                    },
                    "checks": checks,
                    "all_checks_pass": True,
                }
            )

    cases.sort(key=lambda row: row["class_index"])
    top_checks = {
        "source": source["all_checks_pass"],
        "reconstruction": all(reconstruction_checks.values()),
        "group": group_audit["all_checks_pass"],
        "complete_subgroup_denominator": len(public_classes) == 914,
        "exact_high_partition": [row["class_index"] for row in cases] == [4, 12],
        "two_positive_examples": len(cases) == 2 and all(row["all_checks_pass"] for row in cases),
        "p10_denominator": p10["all_checks_pass"],
    }
    if not all(top_checks.values()):
        raise ArithmeticError(f"top-level R12/59 audit changed: {top_checks}")
    return {
        "schema": SCHEMA,
        "artifact_id": "roulleau_r12_ud4x2a1x2_59_all_subgroup_folds",
        "source": source,
        "independent_surface_action_certificate": {
            "path": str(SURFACE_ACTIONS.relative_to(ROOT)),
            "sha256": SURFACE_ACTIONS_SHA256,
            "artifact_id": surface_actions["artifact_id"],
            "all_checks_pass": surface_actions["all_checks_pass"],
        },
        "coordinate_reconstruction": {
            "lattice": "U+D4(-1)^2+A1(-1)^2",
            "rank": 12,
            "wall_count": 59,
            "basis_gram": _matrix_payload(basis_gram),
            "wall_coordinate_rows": _matrix_payload(wall_rows),
            "wall_gram": _matrix_payload(wall_gram),
            "checks": reconstruction_checks,
        },
        "weighted_wall_group": group_audit,
        "subgroup_census": {
            "gap": census["header"],
            "class_count": len(public_classes),
            "proper_high_class_indices": [4, 12],
            "exploratory_all_element_GAP_aliases": {"4": 3, "12": 9},
            "classes": public_classes,
        },
        "cases": cases,
        "P10_denominator": p10,
        "toolchain": {
            "cpp_source": str(CPP_SOURCE.relative_to(ROOT)),
            "cpp_sha256": CPP_SHA256,
            "scout_artifact": str(SCOUT.relative_to(ROOT)),
            "scout_sha256": SCOUT_SHA256,
        },
        "theorem": (
            "The complete 59-wall nef chamber has exactly two proper subgroup-fixed "
            "sections above H9.  Both are finite-volume noncompact H10 Coxeter "
            "polytopes, with 38 and 45 facets.  Both are exact integral divisor-"
            "lattice sections; neither action is induced by an automorphism of the "
            "general K3 surface, and special-period realization remains open."
        ),
        "checks": top_checks,
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate()
    if args.check_archive:
        if payload != json.loads(args.output.read_text()):
            raise SystemExit("recomputed R12/59 certificate differs from archive")
    else:
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
