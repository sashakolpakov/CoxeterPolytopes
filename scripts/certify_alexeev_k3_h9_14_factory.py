#!/usr/bin/env python3
"""Freeze five marked K3-divisor-section rows yielding three ``H9/14`` diagrams."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from typing import Any

import networkx as nx
import sympy as sp

from scripts.certify_alexeev_20_2_1_15_9_1_exact_cones import (
    _analyze_case,
    _parse_cdd_rays,
    _replay_case,
    _run_cdd,
)
from scripts.certify_alexeev_high_dim_reflection_examples import _extract_matrix
from scripts.explore_alexeev_14_10_1_rank10_c1 import (
    DENOMINATOR as SMALL_DENOMINATOR,
    _literal_vinberg,
    _time_orientation,
    build as build_section,
)
from scripts.inventory_alexeev_20_2_1_15_9_1_wall_symmetries import (
    SOURCE,
    SOURCE_SHA256,
    _ambient_reconstruction,
    _lift_actions,
)
from scripts.screen_alexeev_20_2_1_15_9_1_fixed_sections import _matrix
from scripts.search_apery_fermi_chamber_folds import (
    _compose,
    _generated_subgroup,
    _identity,
)


ROOT = Path(__file__).resolve().parents[1]
REMAINING_DENOMINATOR = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_remaining_high_fixed_section_construction.json"
)
OUTPUT = ROOT / "data" / "enumerations" / "alexeev_k3_h9_14_factory.json"
INEQUALITY_DIRECTORY = (
    ROOT / "data" / "enumerations" / "alexeev_k3_h9_14_factory_inequalities"
)
RAY_DIRECTORY = (
    ROOT / "data" / "enumerations" / "alexeev_k3_h9_14_factory_rays"
)
SCHEMA = "coxeter.alexeev_k3_h9_14_factory.v1"


SPECS: tuple[dict[str, Any], ...] = (
    {
        "key": (14, 8, 1),
        "class_index": 1,
        "denominator": SMALL_DENOMINATOR,
        "raw": 16,
        "redundant": 2,
        "determinant": -1024,
        "relations": {"2": 74, "3": 9, "4": 8},
        "finite": 112,
        "ideal": 15,
        "edges": 707,
        "affine": 32,
        "time_min": "21",
        "time_max": "226",
        "source_norms": {"-4": 10, "-2": 14},
        "active_source_norms": {"-4": 8, "-2": 6},
        "folded_facet_norms": {"-8": 2, "-4": 10, "-2": 2},
        "type": "T1",
    },
    {
        "key": (18, 4, 0),
        "class_index": 1,
        "denominator": REMAINING_DENOMINATOR,
        "raw": 26,
        "redundant": 12,
        "determinant": -256,
        "relations": {"2": 75, "3": 8, "4": 8},
        "finite": 90,
        "ideal": 20,
        "edges": 683,
        "affine": 42,
        "time_min": "22",
        "time_max": "136",
        "source_norms": {"-4": 24, "-2": 24},
        "active_source_norms": {"-4": 0, "-2": 14},
        "folded_facet_norms": {"-4": 10, "-2": 4},
        "type": "T2",
    },
    {
        "key": (18, 4, 0),
        "class_index": 7,
        "denominator": REMAINING_DENOMINATOR,
        "raw": 22,
        "redundant": 8,
        "determinant": -81,
        "relations": {"2": 75, "3": 12, "6": 4},
        "finite": 128,
        "ideal": 6,
        "edges": 671,
        "affine": 14,
        "time_min": "26",
        "time_max": "276",
        "source_norms": {"-4": 24, "-2": 24},
        "active_source_norms": {"-4": 0, "-2": 14},
        "folded_facet_norms": {"-6": 5, "-2": 9},
        "type": "T3",
    },
    {
        "key": (18, 4, 0),
        "class_index": 14,
        "denominator": REMAINING_DENOMINATOR,
        "raw": 20,
        "redundant": 6,
        "determinant": -256,
        "relations": {"2": 75, "3": 8, "4": 8},
        "finite": 90,
        "ideal": 20,
        "edges": 683,
        "affine": 42,
        "time_min": "22",
        "time_max": "129",
        "source_norms": {"-4": 24, "-2": 24},
        "active_source_norms": {"-4": 0, "-2": 14},
        "folded_facet_norms": {"-4": 10, "-2": 4},
        "type": "T2",
    },
    {
        "key": (18, 4, 0),
        "class_index": 35,
        "denominator": REMAINING_DENOMINATOR,
        "raw": 18,
        "redundant": 4,
        "determinant": -81,
        "relations": {"2": 75, "3": 12, "6": 4},
        "finite": 128,
        "ideal": 6,
        "edges": 671,
        "affine": 14,
        "time_min": "26",
        "time_max": "269",
        "source_norms": {"-4": 24, "-2": 24},
        "active_source_norms": {"-4": 0, "-2": 14},
        "folded_facet_norms": {"-6": 5, "-2": 9},
        "type": "T3",
    },
)


def _digest(payload: Any) -> str:
    return sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _file_digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _relation_graph(section: dict[str, Any]) -> nx.Graph:
    gram = _matrix(section["restricted_mirrors"]["positive_wall_gram"])
    graph = nx.Graph()
    graph.add_nodes_from(range(gram.rows))
    for left in range(gram.rows):
        for right in range(left + 1, gram.rows):
            ratio = sp.cancel(
                gram[left, right] ** 2 / (gram[left, left] * gram[right, right])
            )
            if ratio:
                graph.add_edge(left, right, ratio=str(ratio))
    return graph


def _isomorphism(left: nx.Graph, right: nx.Graph) -> dict[str, int] | None:
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        left,
        right,
        edge_match=nx.algorithms.isomorphism.categorical_edge_match("ratio", None),
    )
    mapping = next(matcher.isomorphisms_iter(), None)
    return (
        {str(index + 1): int(image) + 1 for index, image in sorted(mapping.items())}
        if mapping is not None
        else None
    )


def _inverse(permutation: tuple[int, ...]) -> tuple[int, ...]:
    inverse = [0] * len(permutation)
    for index, image in enumerate(permutation):
        inverse[image] = index
    return tuple(inverse)


def _ambient_section_orbit_audit(
    bases: dict[str, sp.Matrix],
) -> dict[str, Any]:
    """Prove that the two duplicate A18 diagram rows are congruent sections.

    C14 and C35 are larger marked subgroups than C01 and C07.  In each case
    the larger subgroup has the same fixed space as one of its elements, and
    that element is conjugate in the complete weighted-diagram group to an
    element defining the smaller row.  Since the complete source chamber is
    invariant under that group, this certifies equality up to ambient action
    of the resulting section cones, not merely graph isomorphism.
    """

    denominator = json.loads(REMAINING_DENOMINATOR.read_text(encoding="utf-8"))
    example = next(
        row
        for row in denominator["inventory"]["examples"]
        if tuple(row["key"]) == (18, 4, 0)
    )
    group = tuple(
        tuple(map(int, permutation))
        for permutation in example["weighted_diagram_group"]["elements"]
    )
    if len(group) != 1152 or frozenset(group) != _generated_subgroup(
        tuple(
            tuple(map(int, permutation))
            for permutation in example["weighted_diagram_group"]["generators"]
        )
    ):
        raise ArithmeticError("the complete A18 weighted-diagram group changed")
    classes = {
        int(row["class_index"]): row for row in example["subgroup_classes"]
    }
    subgroups = {
        index: _generated_subgroup(
            tuple(
                tuple(map(int, permutation))
                for permutation in classes[index]["representative_generators"]
            )
        )
        for index in (1, 7, 14, 35)
    }

    wall_gram = _extract_matrix(SOURCE.read_text(encoding="utf-8"), (18, 4, 0))
    ambient_gram, wall_rows = _ambient_reconstruction(wall_gram, 18)
    needed = tuple(sorted(set().union(*subgroups.values())))
    actions = _lift_actions(needed, wall_rows, ambient_gram)
    identity = _identity(wall_gram.rows)

    def full_fixed_elements(class_index: int, section_id: str) -> tuple[tuple[int, ...], ...]:
        basis = bases[section_id]
        candidates = []
        for element in subgroups[class_index]:
            if element == identity:
                continue
            action = actions[element]
            constraint = action - sp.eye(18)
            if constraint.rank() == 8 and basis * constraint == sp.zeros(10, 18):
                candidates.append(element)
        return tuple(sorted(candidates))

    pairs = (
        (1, "A18_4_0_C01", 14, "A18_4_0_C14", 1, 1, 192),
        (7, "A18_4_0_C07", 35, "A18_4_0_C35", 2, 2, 72),
    )
    audits: dict[str, Any] = {}
    for (
        base_index,
        base_id,
        secondary_index,
        secondary_id,
        expected_base_count,
        expected_secondary_count,
        expected_conjugators,
    ) in pairs:
        base_elements = full_fixed_elements(base_index, base_id)
        secondary_elements = full_fixed_elements(secondary_index, secondary_id)
        if not (
            len(base_elements) == expected_base_count
            and len(secondary_elements) == expected_secondary_count
        ):
            raise ArithmeticError("the full-fixed-element census changed")
        base = base_elements[0]
        secondary = secondary_elements[0]
        conjugators = tuple(
            element
            for element in group
            if _compose(_compose(element, base), _inverse(element)) == secondary
        )
        if len(conjugators) != expected_conjugators:
            raise ArithmeticError("the ambient section-conjugacy count changed")
        audits[f"{base_id}__{secondary_id}"] = {
            "base_subgroup_class_index": base_index,
            "secondary_subgroup_class_index": secondary_index,
            "base_full_fixed_element_count": len(base_elements),
            "secondary_full_fixed_element_count": len(secondary_elements),
            "selected_base_element": list(base),
            "selected_secondary_element": list(secondary),
            "conjugator_count": len(conjugators),
            "first_conjugator": list(conjugators[0]),
            "fixed_spaces_conjugate_under_complete_weighted_diagram_group": True,
            "source_section_cones_congruent": True,
        }
    return {
        "complete_weighted_diagram_group_order": len(group),
        "marked_row_count": 5,
        "ambient_section_cone_orbit_count": 3,
        "duplicate_pair_audits": audits,
    }
def build(
    executable: str | None = None,
) -> tuple[dict[str, Any], dict[Path, bytes]]:
    if _file_digest(SOURCE) != SOURCE_SHA256:
        raise ArithmeticError("Alexeev ancillary source hash changed")
    cdd = executable or shutil.which("cddexec_gmp") or "cddexec_gmp"
    records: list[dict[str, Any]] = []
    artifacts: dict[Path, bytes] = {}
    graphs: dict[str, nx.Graph] = {}
    bases: dict[str, sp.Matrix] = {}
    for spec in SPECS:
        section = build_section(
            spec["key"], spec["class_index"], spec["denominator"]
        )
        source_gram = _extract_matrix(
            SOURCE.read_text(encoding="utf-8"), spec["key"]
        )
        source_norms = {
            str(norm): sum(
                source_gram[index, index] == norm
                for index in range(source_gram.rows)
            )
            for norm in (-4, -2)
        }
        active_source_norms = {"-4": 0, "-2": 0}
        for mirror in section["restricted_mirrors"]["records"]:
            if mirror["hyperbolic_status"] != "active_facet":
                continue
            orbit_norms = {
                int(source_gram[int(index) - 1, int(index) - 1])
                for index in mirror["source_wall_indices"]
            }
            if len(orbit_norms) != 1:
                raise ArithmeticError("a retained wall orbit mixes source norms")
            active_source_norms[str(orbit_norms.pop())] += 1
        facet_gram = _matrix(section["restricted_mirrors"]["facet_gram"])
        folded_facet_norms = {
            str(norm): sum(
                facet_gram[index, index] == norm
                for index in range(facet_gram.rows)
            )
            for norm in sorted({int(facet_gram[index, index]) for index in range(facet_gram.rows)})
        }
        reduction = section["exact_reduction"]
        if not (
            section["case"]["fixed_rank"] == 10
            and section["case"]["hyperbolic_dimension"] == 9
            and section["case"]["wall_orbit_count"] == spec["raw"]
            and reduction["active_hyperbolic_facet_count"] == 14
            and reduction["conic_redundant_count"] == spec["redundant"]
            and reduction["unresolved_orbit_indices"] == []
            and section["saturated_fixed_lattice"]["determinant"]
            == spec["determinant"]
            and section["coxeter_gate"]["coxeter_edge_counts"]
            == spec["relations"]
            and section["coxeter_gate"]["all_pairs_coxeter"] is True
            and section["restricted_mirrors"]
            ["all_active_reflections_integral_on_saturated_fixed_lattice"]
            is True
            and source_norms == spec["source_norms"]
            and active_source_norms == spec["active_source_norms"]
            and folded_facet_norms == spec["folded_facet_norms"]
        ):
            raise ArithmeticError(f"section gate changed for {section['id']}")
        cone = _replay_case(section)
        cdd_output = _run_cdd(cone["inequality_text"], cdd)
        rays = _parse_cdd_rays(cdd_output, int(cone["rank"]))
        volume, ray_archive = _analyze_case(section, cone, rays)
        time_orientation = _time_orientation(section, cone, rays)
        vinberg = _literal_vinberg(section, cone, rays)
        exact = volume["exact_double_description"]
        if not (
            volume["decision"]["classification"] == "finite_volume_noncompact"
            and exact["finite_vertex_count"] == spec["finite"]
            and exact["ideal_vertex_count"] == spec["ideal"]
            and exact["hyperideal_vertex_count"] == 0
            and vinberg["literal_vinberg_passes"] is True
            and vinberg["elliptic_edge_count"] == spec["edges"]
            and vinberg["connected_affine_component_count"] == spec["affine"]
            and vinberg["bad_edge_count"] == 0
            and vinberg["bad_parabolic_neighborhood_count"] == 0
            and time_orientation["all_extreme_rays_same_time_cone"] is True
            and time_orientation["minimum_strict_interior_pairing"]
            == spec["time_min"]
            and time_orientation["maximum_strict_interior_pairing"]
            == spec["time_max"]
        ):
            raise ArithmeticError(f"volume gate changed for {section['id']}")
        inequality_path = INEQUALITY_DIRECTORY / f"{section['id']}.ine"
        ray_path = RAY_DIRECTORY / f"{section['id']}.jsonl.gz"
        inequality = cone["inequality_text"].encode()
        artifacts[inequality_path] = inequality
        artifacts[ray_path] = ray_archive
        graphs[section["id"]] = _relation_graph(section)
        bases[section["id"]] = _matrix(
            section["saturated_fixed_lattice"]["basis_primal_rows"]
        )
        records.append(
            {
                "id": section["id"],
                "diagram_type": spec["type"],
                "section": section,
                "exact_cone": volume,
                "time_orientation": time_orientation,
                "literal_vinberg": vinberg,
                "source_wall_audit": {
                    "source_wall_norm_histogram": source_norms,
                    "active_facet_source_norm_histogram": active_source_norms,
                    "folded_primitive_facet_norm_histogram": folded_facet_norms,
                },
                "artifacts": {
                    "inequality_path": str(inequality_path.relative_to(ROOT)),
                    "inequality_sha256": sha256(inequality).hexdigest(),
                    "ray_archive_path": str(ray_path.relative_to(ROOT)),
                    "ray_archive_sha256": sha256(ray_archive).hexdigest(),
                    "facet_gram_sha256": _digest(
                        section["restricted_mirrors"]["facet_gram"]
                    ),
                },
            }
        )

    by_id = {record["id"]: record for record in records}
    identifiers = tuple(by_id)
    pairwise = {}
    for right_index, right in enumerate(identifiers):
        for left in identifiers[:right_index]:
            pairwise[f"{left}__{right}"] = _isomorphism(
                graphs[left], graphs[right]
            )
    expected_pairs = {
        "A18_4_0_C01__A18_4_0_C14",
        "A18_4_0_C07__A18_4_0_C35",
    }
    if {key for key, value in pairwise.items() if value is not None} != expected_pairs:
        raise ArithmeticError("three-type relation-diagram partition changed")
    combined_fixed_ranks = {
        "A18_4_0_C01__A18_4_0_C14": int(
            sp.Matrix.vstack(bases["A18_4_0_C01"], bases["A18_4_0_C14"]).rank()
        ),
        "A18_4_0_C07__A18_4_0_C35": int(
            sp.Matrix.vstack(bases["A18_4_0_C07"], bases["A18_4_0_C35"]).rank()
        ),
    }
    if any(rank == 10 for rank in combined_fixed_ranks.values()):
        raise ArithmeticError("paired marked fixed subspaces unexpectedly coincide")
    ambient_section_orbits = _ambient_section_orbit_audit(bases)
    if ambient_section_orbits["ambient_section_cone_orbit_count"] != 3:
        raise ArithmeticError("the ambient section-cone orbit count changed")

    payload = {
        "schema": SCHEMA,
        "source": {
            "alexeev_ancillary_path": str(SOURCE),
            "alexeev_ancillary_sha256": SOURCE_SHA256,
            "denominators": {
                str(path.relative_to(ROOT)): _file_digest(path)
                for path in (SMALL_DENOMINATOR, REMAINING_DENOMINATOR)
            },
        },
        "scope": {
            "abstract_coxeter_polytopes_certified": True,
            "finite_volume_and_noncompactness_certified": True,
            "integral_fixed_lattice_reflections_certified": True,
            "weak_full_picard_k3_divisor_carriers_certified": True,
            "complete_nef_section_carriers_certified": False,
            "mixed_source_walls_include_square_minus_four": True,
            "surface_automorphism_or_quotient_nef_claimed": False,
            "worldwide_novelty_claimed": False,
        },
        "summary": {
            "raw_marked_subgroup_row_count": 5,
            "ambient_weighted_diagram_section_cone_orbit_count": 3,
            "distinct_abstract_relation_diagram_count": 3,
            "finite_volume_noncompact_count": 5,
            "compact_count": 0,
            "total_bad_edge_count": 0,
            "total_bad_parabolic_neighborhood_count": 0,
        },
        "diagram_types": {
            "T1": ["A14_8_1_C01"],
            "T2": ["A18_4_0_C01", "A18_4_0_C14"],
            "T3": ["A18_4_0_C07", "A18_4_0_C35"],
        },
        "pairwise_relation_diagram_isomorphisms": pairwise,
        "paired_combined_fixed_subspace_ranks": combined_fixed_ranks,
        "ambient_section_cone_orbit_audit": ambient_section_orbits,
        "lineage": {
            "bounded_audit_cutoff": "2026-08-08",
            "T1": {
                "abstract_diagram_status": (
                    "apparently_new_after_bounded_audit_worldwide_priority_open"
                ),
                "arithmetic_status": "rationally_odd_unimodular_integral_overgroup_open",
            },
            "T2": {
                "abstract_diagram_status": (
                    "apparently_new_after_bounded_audit_worldwide_priority_open"
                ),
                "arithmetic_status": (
                    "known_II_9_1_2_power_2_class_literal_14_facet_chamber_unmatched"
                ),
            },
            "T3": {
                "abstract_diagram_status": (
                    "known_Bredon_P_star_and_Scharlau_Kirschmer_II_9_1_3_power_4"
                ),
                "arithmetic_status": "known_full_reflection_chamber",
            },
            "worldwide_new_diagram_count_certified": 0,
        },
        "constructions": records,
        "decision": {
            "classification": "three_k3_divisor_section_orbits_three_H9_14_finite_volume_noncompact_diagrams",
            "new_project_ledger_diagram_count": 3,
            "worldwide_novelty_status": "mixed_T1_T2_open_T3_known",
            "worldwide_new_diagram_count_certified": 0,
        },
    }
    return payload, artifacts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--cdd", default=None)
    args = parser.parse_args()
    payload, artifacts = build(args.cdd)
    rendered = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()
    if args.check_archive:
        if OUTPUT.read_bytes() != rendered:
            raise ArithmeticError("factory certificate archive changed")
        for path, content in artifacts.items():
            if path.read_bytes() != content:
                raise ArithmeticError(f"factory artifact changed: {path}")
    else:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_bytes(rendered)
        for path, content in artifacts.items():
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(content)
    print(json.dumps(payload["summary"], sort_keys=True))
    print(OUTPUT)


if __name__ == "__main__":
    main()
