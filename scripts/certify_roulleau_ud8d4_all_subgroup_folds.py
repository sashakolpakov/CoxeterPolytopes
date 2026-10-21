#!/usr/bin/env python3
"""Certify every high-dimensional diagram fold of Roulleau's U+D8+D4 stock.

Roulleau prints the complete weighted intersection graph of the twenty
(-2)-curves bounding the nef chamber of a K3 surface with
``NS(X)=U+D8(-1)+D4(-1)``.  This certificate transcribes that graph, recovers
the integral Neron--Severi lattice from A1,...,A14, enumerates the complete
weighted-graph group and every subgroup conjugacy class, and certifies every
proper fixed section of rank at least eleven.

The fixed-section construction uses diagram/lattice isometries.  Roulleau's
general surface has automorphism group C2 x C2, whereas the complete weighted
diagram has group C2 x S5.  Consequently no subgroup is silently promoted to
an algebraic surface action.  This distinction does not affect the exact
hyperbolic-polytope construction.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
from fractions import Fraction
import hashlib
import json
from math import gcd
from pathlib import Path
import subprocess
from typing import Any, Iterable, Sequence

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

from scripts.certify_jacobian_kummer_h10_reflection_index import (
    _connected_masks_direct,
    _euler_from_connected_components,
    _relation_matrix,
)
from scripts.certify_jacobian_kummer_h10_maximal_reflection_lineage import (
    _permutation as _relation_permutation,
)
from scripts.certify_roulleau_non2elementary_involution_folds import _scan_coxiter
from scripts.certify_ss_fermat_char3_survivor_fixed_lattices import (
    _saturated_fixed_basis,
)
from scripts.search_k3_high_rank_involution_folds import (
    _canonical_vinberg_audit_digest,
    _coxeter_audit,
    _geometry_certificate,
    _group_checks,
    _json_digest,
    _matrix_payload,
    _weighted_graph_automorphisms,
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
OUTPUT = ROOT / "data/enumerations/roulleau_ud8d4_all_subgroup_folds.json"
SOURCE_INVENTORY = (
    ROOT / "data/enumerations/roulleau_rank11_14_source_stock_inventory.json"
)
SOURCE_TEX = Path("/private/tmp/roulleau_source/epiga_vol6_Roulleau.tex")
SOURCE_TEX_SHA256 = (
    "e6eeaaa4a3f375f209762718e8ef9a361b0a95356cc8be2e9f99a6713abf5df3"
)
SCHEMA = "coxeter.roulleau_ud8d4_all_subgroup_folds.v1"
SOURCE_STOCK = "R14_UD8D4"
RANK = 14
WALL_COUNT = 20
MIN_FIXED_RANK = 11


# (left label, right label, intersection).  The long horizontal chain closes
# at A1.  The four remaining arms are the four parallel length-three paths
# between A10 and A1 printed in the TikZ graph.
SOURCE_EDGES: tuple[tuple[int, int, int], ...] = (
    (1, 2, 1),
    (2, 3, 1),
    (3, 4, 1),
    (4, 6, 1),
    (6, 7, 1),
    (7, 8, 1),
    (8, 9, 1),
    (9, 10, 1),
    (10, 13, 1),
    (13, 18, 2),
    (18, 1, 1),
    (6, 5, 1),
    (10, 11, 1),
    (11, 16, 2),
    (16, 1, 1),
    (10, 12, 1),
    (12, 17, 2),
    (17, 1, 1),
    (10, 14, 1),
    (14, 19, 2),
    (19, 1, 1),
    (10, 15, 1),
    (15, 20, 2),
    (20, 1, 1),
)


EXPECTED_HIGH_CLASSES: dict[int, dict[str, Any]] = {
    4: {"order": 2, "class_size": 10, "fixed_rank": 13},
    6: {"order": 2, "class_size": 15, "fixed_rank": 12},
    7: {"order": 3, "class_size": 10, "fixed_rank": 12},
    8: {"order": 4, "class_size": 5, "fixed_rank": 11},
    12: {"order": 4, "class_size": 15, "fixed_rank": 12},
    14: {"order": 4, "class_size": 15, "fixed_rank": 11},
    19: {"order": 6, "class_size": 10, "fixed_rank": 11},
    22: {"order": 6, "class_size": 10, "fixed_rank": 12},
    23: {"order": 6, "class_size": 10, "fixed_rank": 11},
    29: {"order": 8, "class_size": 15, "fixed_rank": 11},
    35: {"order": 12, "class_size": 5, "fixed_rank": 11},
    41: {"order": 12, "class_size": 10, "fixed_rank": 11},
    48: {"order": 24, "class_size": 5, "fixed_rank": 11},
}


EXPECTED_TYPES: tuple[dict[str, Any], ...] = (
    {
        "classes": [4],
        "dimension": 12,
        "facets": 18,
        "edges": 1628,
        "finite_vertices": 212,
        "cusps": 21,
        "relations": {"0": 4, "2": 132, "3": 15, "4": 2},
    },
    {
        "classes": [6, 12],
        "dimension": 11,
        "facets": 16,
        "edges": 750,
        "finite_vertices": 112,
        "cusps": 9,
        "relations": {"0": 3, "2": 102, "3": 11, "4": 4},
    },
    {
        "classes": [7, 22],
        "dimension": 11,
        "facets": 16,
        "edges": 630,
        "finite_vertices": 94,
        "cusps": 7,
        "relations": {"0": 3, "2": 102, "3": 13, "6": 2},
    },
    {
        "classes": [8, 14, 29, 35, 48],
        "dimension": 10,
        "facets": 14,
        "edges": 254,
        "finite_vertices": 44,
        "cusps": 3,
        "relations": {"0": 4, "2": 76, "3": 11},
    },
    {
        "classes": [19, 23, 41],
        "dimension": 10,
        "facets": 14,
        "edges": 320,
        "finite_vertices": 56,
        "cusps": 3,
        "relations": {"0": 2, "2": 76, "3": 9, "4": 2, "6": 2},
    },
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _source_gram() -> sp.Matrix:
    gram = -2 * sp.eye(WALL_COUNT)
    seen = set()
    for left, right, product in SOURCE_EDGES:
        pair = tuple(sorted((left, right)))
        if pair in seen or left == right:
            raise ArithmeticError("source edge list is not simple")
        seen.add(pair)
        gram[left - 1, right - 1] = gram[right - 1, left - 1] = product
    return gram


def _ambient_reconstruction(
    wall_gram: sp.Matrix,
) -> tuple[sp.Matrix, sp.Matrix]:
    basis_gram = wall_gram[:RANK, :RANK]
    wall_rows = wall_gram[:, :RANK] * basis_gram.inv()
    if not (
        int(basis_gram.det()) == -16
        and exact_inertia(basis_gram) == (1, 13, 0)
        and wall_rows[:RANK, :] == sp.eye(RANK)
        and all(value.is_Integer for value in wall_rows)
        and wall_rows * basis_gram * wall_rows.T == wall_gram
    ):
        raise ArithmeticError("A1,...,A14 did not reconstruct the stated NS lattice")
    return basis_gram, wall_rows


def _greedy_generators(group: Sequence[Permutation]) -> tuple[Permutation, ...]:
    generated: Subgroup = frozenset({_identity(len(group[0]))})
    generators: list[Permutation] = []
    for element in group:
        if element in generated:
            continue
        generators.append(element)
        generated = _generated_subgroup(generators)
        if len(generated) == len(group):
            break
    if generated != frozenset(group):
        raise ArithmeticError("greedy generators do not generate the diagram group")
    return tuple(generators)


def _gap_permutation(permutation: Permutation) -> str:
    return "PermList(" + str([value + 1 for value in permutation]) + ")"


def _gap_subgroup_census(
    group: tuple[Permutation, ...], generators: tuple[Permutation, ...]
) -> dict[str, Any]:
    degree = len(group[0])
    gap_generators = ",".join(_gap_permutation(generator) for generator in generators)
    program = f"""
SizeScreen([1000000,1000000]);;
G:=Group([{gap_generators}]);;
cc:=ConjugacyClassesSubgroups(G);;
gid:=IdGroup(G);;
Print("GROUP|",GAPInfo.Version,"|",Size(G),"|",Length(cc),"|",
      StructureDescription(G),"|",gid[1],"|",gid[2],"\\n");;
i:=0;;
for c in cc do
  i:=i+1;;
  h:=Representative(c);;
  Print("CLASS|",i,"|",Size(h),"|",Size(c),"|",StructureDescription(h),"|",
        List(GeneratorsOfGroup(h),x->ListPerm(x,{degree})),"\\n");;
od;
QUIT;
"""
    completed = subprocess.run(
        ["gap", "-q"],
        input=program,
        text=True,
        capture_output=True,
        timeout=300,
        check=False,
    )
    if completed.returncode != 0 or completed.stderr.strip():
        raise RuntimeError(
            "GAP subgroup census failed:\n"
            + completed.stdout[-4000:]
            + "\n"
            + completed.stderr[-4000:]
        )
    header = None
    rows: list[dict[str, Any]] = []
    for line in completed.stdout.splitlines():
        if line.startswith("GROUP|"):
            _, version, order, count, structure, id_order, id_index = line.split("|")
            header = {
                "gap_version": version,
                "order": int(order),
                "subgroup_class_count": int(count),
                "structure": structure,
                "id_group": [int(id_order), int(id_index)],
            }
        elif line.startswith("CLASS|"):
            _, index, order, class_size, structure, generator_text = line.split("|", 5)
            representative_generators = tuple(
                tuple(int(value) - 1 for value in row)
                for row in ast.literal_eval(generator_text)
            )
            subgroup = (
                _generated_subgroup(representative_generators)
                if representative_generators
                else frozenset({_identity(degree)})
            )
            if len(subgroup) != int(order) or not subgroup.issubset(group):
                raise ArithmeticError("GAP returned an invalid subgroup representative")
            conjugates = {
                frozenset(
                    _compose(_compose(action, element), _inverse(action))
                    for element in subgroup
                )
                for action in group
            }
            if len(conjugates) != int(class_size):
                raise ArithmeticError("GAP subgroup class size failed direct replay")
            rows.append(
                {
                    "class_index": int(index),
                    "order": int(order),
                    "conjugacy_class_size": int(class_size),
                    "structure": structure,
                    "representative_generators": representative_generators,
                    "representative_subgroup": subgroup,
                }
            )
    if header is None or not (
        header["order"] == len(group)
        and header["subgroup_class_count"] == len(rows)
    ):
        raise ArithmeticError("GAP subgroup census header mismatch")
    return {"header": header, "rows": rows}


def _inverse(permutation: Permutation) -> Permutation:
    result = [0] * len(permutation)
    for source, target in enumerate(permutation):
        result[target] = source
    return tuple(result)


def _lift_actions(
    group: Sequence[Permutation], wall_rows: sp.Matrix, basis_gram: sp.Matrix
) -> dict[Permutation, sp.Matrix]:
    actions = {}
    for permutation in group:
        action = wall_rows[[permutation[index] for index in range(RANK)], :]
        expected = wall_rows[[permutation[index] for index in range(WALL_COUNT)], :]
        if not (
            all(value.is_Integer for value in action)
            and abs(int(action.det())) == 1
            and wall_rows * action == expected
            and action * basis_gram * action.T == basis_gram
        ):
            raise ArithmeticError("a diagram symmetry failed its integral NS lift")
        actions[permutation] = action
    return actions


def _divisibility(gram: sp.Matrix, root: sp.Matrix) -> int:
    value = 0
    for product in root * gram:
        value = gcd(value, abs(int(product)))
    return value


def _fixed_lattice_and_reflections(
    generators: Sequence[Permutation],
    actions: dict[Permutation, sp.Matrix],
    subgroup: Subgroup,
    cycles_zero_based: Sequence[Sequence[int]],
    basis_gram: sp.Matrix,
    wall_rows: sp.Matrix,
    wall_gram: sp.Matrix,
) -> tuple[dict[str, Any], dict[str, Any], sp.Matrix]:
    generator_actions = tuple(actions[generator] for generator in generators)
    fixed_basis, smith_kernel = _saturated_fixed_basis(generator_actions, RANK)
    fixed_gram = fixed_basis * basis_gram * fixed_basis.T
    fixed_rank = fixed_gram.rows
    orbit_roots = []
    reflection_checks = []
    for orbit_index, orbit in enumerate(cycles_zero_based, start=1):
        normal = sum(
            (wall_rows.row(index) for index in orbit),
            sp.zeros(1, RANK),
        )
        solution = sp.linsolve((fixed_basis.T, normal.T))
        if solution is sp.EmptySet:
            raise ArithmeticError("an orbit root is absent from the fixed lattice")
        root = sp.Matrix(tuple(next(iter(solution)))).T
        if not (root * fixed_basis == normal and all(value.is_Integer for value in root)):
            raise ArithmeticError("an orbit root is not integral in the fixed lattice")
        norm = int((root * fixed_gram * root.T)[0])
        reflection = (
            sp.eye(fixed_rank)
            - sp.Rational(2, norm) * fixed_gram * root.T * root
        )
        disjoint = _pairwise_disjoint(wall_gram, tuple(orbit))
        check = {
            "orbit_index": orbit_index,
            "source_walls_one_based": [index + 1 for index in orbit],
            "fixed_basis_coordinates": [int(value) for value in root],
            "orbit_size": len(orbit),
            "pairwise_disjoint_source_roots": disjoint,
            "norm": norm,
            "divisibility": _divisibility(fixed_gram, root),
            "reflection_integral": all(value.is_Integer for value in reflection),
            "reflection_is_involutive_isometry": (
                reflection * reflection == sp.eye(fixed_rank)
                and reflection * fixed_gram * reflection.T == fixed_gram
            ),
        }
        check["integral_reflection_on_fixed_lattice"] = all(
            (
                disjoint,
                norm == -2 * len(orbit),
                check["reflection_integral"],
                check["reflection_is_involutive_isometry"],
            )
        )
        reflection_checks.append(check)
        orbit_roots.append(root)
    root_rows = sp.Matrix.vstack(*orbit_roots)
    orbit_gram = root_rows * fixed_gram * root_rows.T
    direct_orbit_gram = _orbit_gram(wall_gram, tuple(tuple(o) for o in cycles_zero_based))
    if orbit_gram != direct_orbit_gram:
        raise ArithmeticError("fixed-lattice and direct orbit Grams disagree")
    smith = smith_normal_form(fixed_gram, domain=ZZ)
    lattice = {
        "basis_rows_in_A1_to_A14_basis": _matrix_payload(fixed_basis),
        "gram": _matrix_payload(fixed_gram),
        "rank": fixed_rank,
        "inertia": list(exact_inertia(fixed_gram)),
        "determinant": int(fixed_gram.det()),
        "smith_diagonal_absolute": [
            abs(int(smith[index, index])) for index in range(fixed_rank)
        ],
        "even": all(int(fixed_gram[index, index]) % 2 == 0 for index in range(fixed_rank)),
        "smith_kernel_audit": smith_kernel,
        "orbit_root_coordinate_rows": _matrix_payload(root_rows),
    }
    reflections = {
        "all_folded_reflections_are_integral": all(
            check["integral_reflection_on_fixed_lattice"]
            for check in reflection_checks
        ),
        "checks": reflection_checks,
    }
    return lattice, reflections, orbit_gram


def _ns_action_audit(
    generators: Sequence[Permutation],
    actions: dict[Permutation, sp.Matrix],
    subgroup: Subgroup,
    fixed_rank: int,
) -> dict[str, Any]:
    checks = {
        "every_generator_is_integral": all(
            all(value.is_Integer for value in actions[generator])
            for generator in generators
        ),
        "every_generator_is_unimodular": all(
            abs(int(actions[generator].det())) == 1 for generator in generators
        ),
        "generator_actions_are_faithful": len(
            {
                tuple(int(value) for value in actions[element])
                for element in subgroup
            }
        )
        == len(subgroup),
    }
    return {
        "basis_wall_labels": list(range(1, RANK + 1)),
        "representative_generator_permutations_one_based": [
            [value + 1 for value in generator] for generator in generators
        ],
        "representative_generator_action_matrices": [
            _matrix_payload(actions[generator]) for generator in generators
        ],
        "fixed_rank": fixed_rank,
        "checks": checks,
        "all_integral_lattice_checks_pass": all(checks.values()),
        "surface_realization_status": (
            "diagram_and_NS_lattice_isometry; general-surface Aut(X)=C2xC2 is "
            "strictly smaller than the C2xS5 wall-diagram group, so no surface "
            "automorphism realization is asserted"
        ),
    }


def _relation_histogram(relations: list[list[str]]) -> dict[str, int]:
    return dict(
        sorted(
            Counter(
                relations[left][right]
                for left in range(len(relations))
                for right in range(left + 1, len(relations))
            ).items()
        )
    )


def _cluster_relations(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    clusters: list[dict[str, Any]] = []
    for record in records:
        relations = record["normalized_relation_matrix"]
        match = None
        match_permutation = None
        for cluster in clusters:
            permutation = _relation_permutation(
                cluster["normalized_relation_matrix"], relations
            )
            if permutation is not None:
                match = cluster
                match_permutation = permutation
                break
        if match is None:
            cluster = {
                "type_id": f"R14_UD8D4_TYPE_{len(clusters) + 1}",
                "representative_class_index": record["class_index"],
                "class_indices": [record["class_index"]],
                "normalized_relation_matrix": relations,
                "transports_from_representative": {
                    str(record["class_index"]): list(range(len(relations)))
                },
            }
            clusters.append(cluster)
        else:
            match["class_indices"].append(record["class_index"])
            match["transports_from_representative"][str(record["class_index"])] = (
                match_permutation
            )
    return clusters


def _source_checks() -> dict[str, Any]:
    if _sha256(SOURCE_TEX) != SOURCE_TEX_SHA256:
        raise ArithmeticError("pinned Roulleau TeX hash changed")
    inventory = json.loads(SOURCE_INVENTORY.read_text())
    stock = next(stock for stock in inventory["stocks"] if stock["key"] == SOURCE_STOCK)
    lines = SOURCE_TEX.read_text().splitlines()
    excerpt = "\n".join(lines[9668:9827])
    checks = {
        "source_inventory_is_certified": inventory["all_checks_pass"] is True,
        "inventory_stock_rank": stock["picard_rank"] == RANK,
        "inventory_stock_wall_count": stock["wall_count"] == WALL_COUNT,
        "section_lattice_token": r"U\oplus\mathbf{D}_{8}\oplus\mathbf{D}_{4}" in excerpt,
        "section_wall_count": "contains $20$ $(-2)$-curves" in excerpt,
        "basis_statement": "$A_{1},\\dots,A_{14}$ generate" in excerpt,
        "four_thick_arms": excerpt.count(r"\draw [very thick]") == 5,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"Roulleau source transcription checks failed: {checks}")
    return checks


def build_certificate() -> dict[str, Any]:
    source_checks = _source_checks()
    wall_gram = _source_gram()
    basis_gram, wall_rows = _ambient_reconstruction(wall_gram)
    group, search_nodes = _weighted_graph_automorphisms(wall_gram)
    group_checks = _group_checks(group)
    if not (len(group) == 240 and all(group_checks.values())):
        raise ArithmeticError("unexpected complete weighted-diagram group")
    actions = _lift_actions(group, wall_rows, basis_gram)
    generators = _greedy_generators(group)
    census = _gap_subgroup_census(group, generators)
    print(
        json.dumps(
            {
                "stage": "subgroup_census",
                "group_order": len(group),
                "subgroup_classes": len(census["rows"]),
            }
        ),
        flush=True,
    )
    if not (
        census["header"]["structure"] == "C2 x S5"
        and census["header"]["id_group"] == [240, 189]
        and census["header"]["subgroup_class_count"] == 57
    ):
        raise ArithmeticError(f"unexpected GAP group description: {census['header']}")

    class_records = []
    high_internal = []
    for row in census["rows"]:
        subgroup: Subgroup = row.pop("representative_subgroup")
        representative_generators: tuple[Permutation, ...] = row.pop(
            "representative_generators"
        )
        cycles_zero = _orbits(subgroup, WALL_COUNT)
        orbit_gram = _orbit_gram(wall_gram, cycles_zero)
        fixed_rank = orbit_gram.rank()
        element_histogram = {
            str(order): count
            for order, count in sorted(
                Counter(_order(element) for element in subgroup).items()
            )
        }
        proper_high = len(subgroup) > 1 and fixed_rank >= MIN_FIXED_RANK
        public = {
            **row,
            "fixed_rank": fixed_rank,
            "hyperbolic_dimension": fixed_rank - 1,
            "proper_fixed_section_above_H9": proper_high,
            "element_order_histogram": element_histogram,
            "representative_generators_zero_based": [
                list(generator) for generator in representative_generators
            ],
            "wall_orbits_one_based": [
                [value + 1 for value in orbit] for orbit in cycles_zero
            ] if proper_high else None,
        }
        class_records.append(public)
        if proper_high:
            lattice, reflection_audit, exact_orbit_gram = (
                _fixed_lattice_and_reflections(
                    representative_generators,
                    actions,
                    subgroup,
                    cycles_zero,
                    basis_gram,
                    wall_rows,
                    wall_gram,
                )
            )
            cycles_one = [[value + 1 for value in orbit] for orbit in cycles_zero]
            facets = list(range(exact_orbit_gram.rows))
            coxeter = _coxeter_audit(exact_orbit_gram, facets, cycles_one)
            positive_gram = [
                [-int(exact_orbit_gram[left, right]) for right in facets]
                for left in facets
            ]
            relations = _relation_matrix(positive_gram)
            ns_action = _ns_action_audit(
                representative_generators, actions, subgroup, fixed_rank
            )
            if not (
                exact_inertia(exact_orbit_gram)
                == (1, fixed_rank - 1, exact_orbit_gram.rows - fixed_rank)
                and coxeter["all_pairs_are_exact_coxeter_relations"]
                and reflection_audit["all_folded_reflections_are_integral"]
                and ns_action["all_integral_lattice_checks_pass"]
            ):
                raise ArithmeticError("a high subgroup failed the exact pre-volume gate")
            high_internal.append(
                {
                    **public,
                    "subgroup_elements_zero_based": [
                        list(element) for element in sorted(subgroup)
                    ],
                    "fixed_lattice": lattice,
                    "reflection_audit": reflection_audit,
                    "orbit_gram": _matrix_payload(exact_orbit_gram),
                    "normalized_relation_matrix": relations,
                    "relation_histogram": _relation_histogram(relations),
                    "coxeter": coxeter,
                    "ns_action": ns_action,
                    "_subgroup": subgroup,
                    "_generators": representative_generators,
                    "_cycles_one": cycles_one,
                    "_orbit_gram_matrix": exact_orbit_gram,
                }
            )

    actual_high = {
        record["class_index"]: {
            "order": record["order"],
            "class_size": record["conjugacy_class_size"],
            "fixed_rank": record["fixed_rank"],
        }
        for record in high_internal
    }
    if actual_high != EXPECTED_HIGH_CLASSES:
        raise ArithmeticError(f"high-subgroup denominator changed: {actual_high}")
    print(
        json.dumps(
            {"stage": "high_fixed_sections", "class_count": len(high_internal)}
        ),
        flush=True,
    )

    clusters = _cluster_relations(high_internal)
    actual_clusters = [cluster["class_indices"] for cluster in clusters]
    expected_clusters = [row["classes"] for row in EXPECTED_TYPES]
    if actual_clusters != expected_clusters:
        raise ArithmeticError(
            f"normalized relation clustering changed: {actual_clusters}"
        )

    type_records = []
    class_to_type: dict[int, str] = {}
    for cluster, expected in zip(clusters, EXPECTED_TYPES):
        representative = next(
            record
            for record in high_internal
            if record["class_index"] == cluster["representative_class_index"]
        )
        print(
            json.dumps(
                {
                    "stage": "vinberg_type",
                    "type_id": cluster["type_id"],
                    "classes": cluster["class_indices"],
                    "dimension": expected["dimension"],
                    "facets": expected["facets"],
                }
            ),
            flush=True,
        )
        geometry = _geometry_certificate(
            RANK,
            representative["_cycles_one"],
            representative["_orbit_gram_matrix"],
            representative["coxeter"],
            representative["reflection_audit"],
            representative["ns_action"],
        )
        checks = {
            "full_geometry_certificate": geometry["all_checks_pass"],
            "finite_volume_noncompact": (
                geometry["finite_volume_certified"]
                and geometry["compact"] is False
            ),
            "dimension": geometry["hyperbolic_dimension"] == expected["dimension"],
            "facets": geometry["facet_count"] == expected["facets"],
            "elliptic_edges": geometry["vinberg_volume_criterion"]["elliptic_edge_count"]
            == expected["edges"],
            "finite_vertices": geometry["finite_vertex_count"]
            == expected["finite_vertices"],
            "cusps": geometry["ideal_vertex_count"] == expected["cusps"],
            "relations": representative["relation_histogram"] == expected["relations"],
        }
        if not all(checks.values()):
            raise ArithmeticError(f"type geometry checks failed: {checks}")
        connected = _connected_masks_direct(cluster["normalized_relation_matrix"])
        solomon = _euler_from_connected_components(
            cluster["normalized_relation_matrix"], connected
        )
        coxiter = _scan_coxiter(
            cluster["normalized_relation_matrix"],
            expected["dimension"],
            expected["facets"],
        )
        type_record = {
            **cluster,
            "hyperbolic_dimension": expected["dimension"],
            "facet_count": expected["facets"],
            "volume_type": "finite_volume_noncompact",
            "geometry": geometry,
            "solomon": solomon,
            "lineage": {
                "coxiter": coxiter,
                "exact_CoxIter_match_count": (
                    len(coxiter["ordinary_exact_matches"])
                    + len(coxiter["doubling_exact_matches"])
                ),
                "worldwide_priority_status": (
                    "known_exact_CoxIter_diagram" if (
                        coxiter["ordinary_exact_matches"]
                        or coxiter["doubling_exact_matches"]
                    ) else "no_exact_match_in_pinned_CoxIter_corpus; broader_literature_priority_open"
                ),
            },
            "checks": checks,
            "all_checks_pass": all(checks.values()),
        }
        type_records.append(type_record)
        for class_index in cluster["class_indices"]:
            class_to_type[class_index] = cluster["type_id"]

    high_public = []
    for record in high_internal:
        public = {key: value for key, value in record.items() if not key.startswith("_")}
        public["geometry_type_id"] = class_to_type[record["class_index"]]
        public["status"] = "certified_finite_volume_noncompact_Coxeter_polytope"
        high_public.append(public)

    all_checks = {
        "source_transcription": all(source_checks.values()),
        "full_wall_gram": (
            wall_gram.rank() == RANK
            and exact_inertia(wall_gram) == (1, 13, 6)
        ),
        "integral_NS_reconstruction": int(basis_gram.det()) == -16,
        "diagram_group": (
            len(group) == 240
            and census["header"]["structure"] == "C2 x S5"
            and census["header"]["id_group"] == [240, 189]
        ),
        "all_57_subgroup_classes": len(class_records) == 57,
        "all_13_high_classes": len(high_public) == 13,
        "exactly_five_diagram_types": len(type_records) == 5,
        "all_types_finite_volume_noncompact": all(
            record["all_checks_pass"] for record in type_records
        ),
    }
    if not all(all_checks.values()):
        raise ArithmeticError(f"global checks failed: {all_checks}")

    return {
        "schema": SCHEMA,
        "artifact_id": "roulleau_ud8d4_all_subgroup_folds",
        "source": {
            "stock_key": SOURCE_STOCK,
            "lattice": "U+D8(-1)+D4(-1)",
            "picard_rank": RANK,
            "complete_minus_two_curve_count": WALL_COUNT,
            "roulleau_tex_section_lines": [9669, 9827],
            "roulleau_tex_sha256": SOURCE_TEX_SHA256,
            "source_inventory": str(SOURCE_INVENTORY.relative_to(ROOT)),
            "source_inventory_sha256": _sha256(SOURCE_INVENTORY),
            "source_edges_one_based": [list(edge) for edge in SOURCE_EDGES],
            "wall_gram": _matrix_payload(wall_gram),
            "basis_wall_labels": list(range(1, RANK + 1)),
            "basis_gram": _matrix_payload(basis_gram),
            "wall_coordinate_rows": _matrix_payload(wall_rows),
            "completeness": "Roulleau_published_complete_curve_configuration",
            "checks": source_checks,
        },
        "weighted_diagram_group": {
            **census["header"],
            "exhaustive_backtracking_node_count": search_nodes,
            "element_order_histogram": {
                str(order): count
                for order, count in sorted(
                    Counter(_order(element) for element in group).items()
                )
            },
            "greedy_generators_zero_based": [list(generator) for generator in generators],
            "group_checks": group_checks,
            "surface_automorphism_warning": (
                "Roulleau records Aut(X)=C2xC2 for a general surface, while the "
                "weighted wall diagram has C2xS5; the full diagram group is used "
                "only as an exact divisor-lattice folding group"
            ),
        },
        "subgroup_census": {
            "class_count_including_identity_and_full_group": len(class_records),
            "proper_fixed_rank_at_least_11_count": len(high_public),
            "classes": class_records,
        },
        "high_fixed_sections": high_public,
        "abstract_diagram_types": type_records,
        "summary": {
            "proper_high_subgroup_class_count": len(high_public),
            "distinct_normalized_diagram_count": len(type_records),
            "finite_volume_noncompact_type_count": len(type_records),
            "dimension_facet_pairs": sorted(
                {
                    (record["hyperbolic_dimension"], record["facet_count"])
                    for record in type_records
                }
            ),
            "new_divisor_lattice_fixed_section_realizations": len(type_records),
            "new_abstract_polytope_priority": (
                "resolved_per_type_against_pinned_CoxIter; any unmatched type still "
                "requires broader literature search"
            ),
        },
        "checks": all_checks,
        "all_checks_pass": all(all_checks.values()),
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
