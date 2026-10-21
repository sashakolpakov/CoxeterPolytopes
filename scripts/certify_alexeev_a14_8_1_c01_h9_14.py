#!/usr/bin/env python3
"""Freeze the ``A14_8_1_C01`` finite-volume ``H9/14`` Coxeter polytope.

The source is Alexeev's complete 24-wall chamber for the 2-elementary K3
lattice ``(r,a,delta)=(14,8,1)``.  A diagram involution has saturated fixed
rank ten.  Its sixteen wall-orbit inequalities reduce exactly to fourteen
facets.  This certificate replays the fixed lattice, Coxeter and integral
reflection gates, GMP extreme-ray census, and the literal Vinberg
edge/affine-neighbourhood criterion.  A geometric K3 automorphism lift and
worldwide novelty are deliberately outside its scope.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from typing import Any

from scripts.certify_alexeev_high_dim_reflection_examples import _extract_matrix
from scripts.certify_alexeev_20_2_1_15_9_1_exact_cones import (
    _analyze_case,
    _parse_cdd_rays,
    _replay_case,
    _run_cdd,
)
from scripts.explore_alexeev_14_10_1_rank10_c1 import (
    DENOMINATOR,
    _literal_vinberg,
    _time_orientation,
    build as build_section,
)
from scripts.inventory_alexeev_20_2_1_15_9_1_wall_symmetries import (
    SOURCE,
    SOURCE_SHA256,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_a14_8_1_c01_h9_14.json"
)
INEQUALITY = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_a14_8_1_c01_h9_14.ine"
)
RAYS = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_a14_8_1_c01_h9_14_rays.jsonl.gz"
)
SCHEMA = "coxeter.alexeev_a14_8_1_c01_h9_14.v1"
KEY = (14, 8, 1)
CLASS_INDEX = 1
EXPECTED_FACET_GRAM_SHA256 = (
    "55611ee3e651320e9e7f6016ebf83724ad637e75332e23f111c7c3b30b19e7a5"
)


def _sha256_bytes(payload: bytes) -> str:
    return sha256(payload).hexdigest()


def _sha256_path(path: Path) -> str:
    return _sha256_bytes(path.read_bytes())


def _payload_sha256(payload: Any) -> str:
    return _sha256_bytes(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    )


def build(executable: str | None = None) -> tuple[dict[str, Any], bytes, bytes]:
    if _sha256_path(SOURCE) != SOURCE_SHA256:
        raise ArithmeticError("Alexeev ancillary source hash changed")
    section = build_section(KEY, CLASS_INDEX)
    source_gram = _extract_matrix(SOURCE.read_text(encoding="utf-8"), KEY)
    source_diagonal_histogram = {
        str(value): sum(source_gram[index, index] == value for index in range(source_gram.rows))
        for value in (-4, -2)
    }
    active_source_norm_histogram = {"-4": 0, "-2": 0}
    for record in section["restricted_mirrors"]["records"]:
        if record["hyperbolic_status"] != "active_facet":
            continue
        source_indices = [int(index) - 1 for index in record["source_wall_indices"]]
        source_norms = {int(source_gram[index, index]) for index in source_indices}
        if len(source_norms) != 1:
            raise ArithmeticError("a wall orbit mixes source norms")
        active_source_norm_histogram[str(source_norms.pop())] += 1
    facet_hash = _payload_sha256(section["restricted_mirrors"]["facet_gram"])
    if not (
        section["id"] == "A14_8_1_C01"
        and section["case"]["fixed_rank"] == 10
        and section["case"]["hyperbolic_dimension"] == 9
        and section["case"]["source_wall_count"] == 24
        and section["case"]["wall_orbit_count"] == 16
        and section["exact_reduction"]["active_hyperbolic_facet_count"] == 14
        and section["exact_reduction"]["conic_redundant_count"] == 2
        and section["exact_reduction"]["unresolved_orbit_indices"] == []
        and section["coxeter_gate"]["all_pairs_coxeter"] is True
        and section["coxeter_gate"]["non_coxeter_pair_count"] == 0
        and section["restricted_mirrors"]
        ["all_active_reflections_integral_on_saturated_fixed_lattice"]
        is True
        and section["saturated_fixed_lattice"]["determinant"] == -1024
        and facet_hash == EXPECTED_FACET_GRAM_SHA256
        and source_diagonal_histogram == {"-4": 10, "-2": 14}
        and active_source_norm_histogram == {"-4": 8, "-2": 6}
    ):
        raise ArithmeticError("A14_8_1_C01 section gate changed")

    cone = _replay_case(section)
    cdd = executable or shutil.which("cddexec_gmp") or "cddexec_gmp"
    cdd_output = _run_cdd(cone["inequality_text"], cdd)
    rays = _parse_cdd_rays(cdd_output, int(cone["rank"]))
    volume, ray_archive = _analyze_case(section, cone, rays)
    time_orientation = _time_orientation(section, cone, rays)
    vinberg = _literal_vinberg(section, cone, rays)
    exact = volume["exact_double_description"]
    if not (
        volume["decision"]["classification"] == "finite_volume_noncompact"
        and exact["extreme_ray_count"] == 127
        and exact["finite_vertex_count"] == 112
        and exact["ideal_vertex_count"] == 15
        and exact["hyperideal_vertex_count"] == 0
        and vinberg["literal_vinberg_passes"] is True
        and vinberg["elliptic_edge_count"] == 707
        and vinberg["finite_endpoint_count"] == 112
        and vinberg["ideal_endpoint_count"] == 15
        and vinberg["connected_affine_component_count"] == 32
        and vinberg["bad_edge_count"] == 0
        and vinberg["bad_parabolic_neighborhood_count"] == 0
        and time_orientation["all_extreme_rays_same_time_cone"] is True
        and time_orientation["minimum_strict_interior_pairing"] == "21"
        and time_orientation["maximum_strict_interior_pairing"] == "226"
    ):
        raise ArithmeticError("A14_8_1_C01 volume/Vinberg gate changed")

    inequality_bytes = cone["inequality_text"].encode()
    payload = {
        "schema": SCHEMA,
        "source": {
            "alexeev_ancillary_path": str(SOURCE),
            "alexeev_ancillary_sha256": SOURCE_SHA256,
            "subgroup_denominator_path": str(DENOMINATOR.relative_to(ROOT)),
            "subgroup_denominator_sha256": _sha256_path(DENOMINATOR),
            "ambient_k3_lattice_invariants": list(KEY),
            "weighted_diagram_subgroup_class_index": CLASS_INDEX,
        },
        "scope": {
            "abstract_coxeter_polytope_certified": True,
            "finite_volume_and_compactness_certified": True,
            "integral_fixed_lattice_reflections_certified": True,
            "geometric_k3_action_lift_certified": False,
            "weak_full_picard_k3_divisor_carrier_certified": True,
            "complete_nef_section_certified": False,
            "source_wall_norm_histogram": source_diagonal_histogram,
            "active_facet_source_norm_histogram": active_source_norm_histogram,
            "worldwide_novelty_certified": False,
            "finite_volume_noncompact_only": True,
        },
        "section": section,
        "exact_cone": volume,
        "time_orientation": time_orientation,
        "literal_vinberg": vinberg,
        "artifacts": {
            "inequality_path": str(INEQUALITY.relative_to(ROOT)),
            "inequality_sha256": _sha256_bytes(inequality_bytes),
            "ray_archive_path": str(RAYS.relative_to(ROOT)),
            "ray_archive_sha256": _sha256_bytes(ray_archive),
            "facet_gram_sha256": facet_hash,
        },
        "decision": {
            "classification": "finite_volume_noncompact_hyperbolic_coxeter_polytope",
            "hyperbolic_dimension": 9,
            "facet_count": 14,
            "finite_vertex_count": 112,
            "ideal_vertex_count": 15,
            "hyperideal_vertex_count": 0,
            "literal_vinberg_bad_edge_count": 0,
            "literal_vinberg_bad_parabolic_neighborhood_count": 0,
            "new_project_ledger_item": True,
            "worldwide_novelty_status": (
                "apparently_new_after_bounded_audit_worldwide_priority_open"
            ),
        },
    }
    return payload, inequality_bytes, ray_archive


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--cdd", default=None)
    args = parser.parse_args()
    payload, inequality, rays = build(args.cdd)
    rendered = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()
    if args.check_archive:
        if not (
            OUTPUT.read_bytes() == rendered
            and INEQUALITY.read_bytes() == inequality
            and RAYS.read_bytes() == rays
        ):
            raise ArithmeticError("A14_8_1_C01 frozen archive changed")
    else:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_bytes(rendered)
        INEQUALITY.write_bytes(inequality)
        RAYS.write_bytes(rays)
    print(json.dumps(payload["decision"], sort_keys=True))
    print(OUTPUT)


if __name__ == "__main__":
    main()
