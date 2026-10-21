#!/usr/bin/env python3
"""Certify the exact McLeod-group lineage and index of F3-C2-H10-196.

The fixed-lattice certificate supplies an explicit unimodular isometry from
the half-scaled divisor lattice to ``<3> + <-1>^10`` and verifies that all
196 facet mirrors are integral reflecting roots.  This script combines that
inclusion with complete Solomon Euler sums for the 196-wall chamber and the
14-wall McLeod chamber in the pinned CoxIter corpus.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any

from scripts.certify_jacobian_kummer_h10_maximal_reflection_lineage import (
    _parse_coxiter,
)
from scripts.certify_jacobian_kummer_h10_reflection_index import (
    _connected_masks_direct,
    _euler_from_connected_components,
    _factorization,
    _relation_matrix,
)


ROOT = Path(__file__).resolve().parents[1]
SCREEN = (
    ROOT
    / "data"
    / "enumerations"
    / "ss_fermat_quartic_char3_c2_h10_coxeter_screen.json"
)
COMPONENTS = (
    ROOT
    / "data"
    / "enumerations"
    / "ss_fermat_quartic_char3_c2_h10_component_audit.json"
)
VINBERG = (
    ROOT
    / "data"
    / "enumerations"
    / "ss_fermat_quartic_char3_c2_h10_vinberg_certificate.json"
)
CONE = (
    ROOT
    / "data"
    / "enumerations"
    / "ss_fermat_quartic_char3_c2_h10_exact_cone_certificate.json"
)
FIXED = (
    ROOT
    / "data"
    / "enumerations"
    / "ss_fermat_quartic_char3_survivor_fixed_lattices.json"
)
CORPUS = Path("/private/tmp/coxiter-current/CoxIter-master")
MCLEOD_RELATIVE = Path(
    "graphs/group units qf/10-3,1,1,1,1,1,1,1,1,1,1.coxiter"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "ss_fermat_quartic_char3_c2_h10_mcleod_index.json"
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _fraction(value: str) -> Fraction:
    return Fraction(value)


def build_certificate(corpus: Path = CORPUS) -> dict[str, Any]:
    screen = json.loads(SCREEN.read_text(encoding="utf-8"))
    components = json.loads(COMPONENTS.read_text(encoding="utf-8"))
    vinberg = json.loads(VINBERG.read_text(encoding="utf-8"))
    cone = json.loads(CONE.read_text(encoding="utf-8"))
    fixed = json.loads(FIXED.read_text(encoding="utf-8"))

    mirrors = fixed["c2_h10_196_mirrors"]
    upstream = {
        "screen_has_196_coxeter_facets": (
            screen["facet_count"] == 196
            and screen["status"]["all_facet_pairs_are_coxeter"]
        ),
        "component_audit_has_no_invalid_cusp_neighborhood": (
            components["facet_count"] == 196
            and components["invalid_full_zero_neighborhood_count"] == 0
        ),
        "literal_vinberg_is_finite_volume_noncompact": (
            vinberg["finite_volume"]
            and not vinberg["compact"]
            and vinberg["bad_edge_count"] == 0
            and vinberg["unrecognized_endpoint_edge_count"] == 0
        ),
        "independent_cone_has_no_hyperideal_rays": (
            cone["finite_volume_certified"]
            and cone["hyperideal_ray_count"] == 0
            and cone["literal_vinberg_endpoint_counts_match"]
        ),
        "fixed_lattice_certificate_passes": fixed["all_checks_pass"],
        "all_196_mirrors_are_integral_mcleod_roots": (
            mirrors["facet_count"] == 196
            and mirrors["all_reflections_integral_involutive_isometries"]
        ),
    }
    if not all(upstream.values()):
        raise ArithmeticError(f"an upstream exact gate failed: {upstream}")

    subgroup_relations = _relation_matrix(screen["positive_wall_gram"])
    subgroup_connected = {
        sum(1 << (vertex - 1) for vertex in component)
        for component in components["connected_elliptic_components"]
    }
    subgroup = _euler_from_connected_components(
        subgroup_relations, subgroup_connected
    )

    mcleod_path = corpus / MCLEOD_RELATIVE
    facets, dimension, base_relations = _parse_coxiter(mcleod_path)
    base_connected = _connected_masks_direct(base_relations)
    base = _euler_from_connected_components(base_relations, base_connected)

    subgroup_euler = _fraction(subgroup["rational_euler_characteristic"])
    base_euler = _fraction(base["rational_euler_characteristic"])
    index = subgroup_euler / base_euler
    volume_coefficient = -subgroup_euler * Fraction(32, 945)
    reconstructed_f_vector = [
        subgroup["finite_special_subset_count_by_rank"][10]
        + vinberg["ideal_vertex_count"],
        *reversed(subgroup["finite_special_subset_count_by_rank"][1:10]),
        1,
    ]
    expected_f_vector = [
        223_384,
        1_427_720,
        3_435_180,
        4_355_400,
        3_290_232,
        1_543_260,
        447_360,
        76_260,
        6_680,
        196,
        1,
    ]
    checks = {
        **upstream,
        "mcleod_corpus_chamber_has_14_facets_in_h10": (
            facets == 14 and dimension == 10
        ),
        "subgroup_connected_count_matches_archive": (
            len(subgroup_connected)
            == len(components["connected_elliptic_components"])
            == 38_266
        ),
        "subgroup_rank_10_count_matches_finite_vertices": (
            subgroup["finite_special_subset_count_by_rank"][10]
            == vinberg["finite_vertex_count"]
            == 212_824
        ),
        "subgroup_rank_9_count_matches_vinberg_edges": (
            subgroup["finite_special_subset_count_by_rank"][9]
            == vinberg["elliptic_edge_count"]
            == 1_427_720
        ),
        "subgroup_euler_is_expected": subgroup_euler
        == Fraction(-5_797, 3_780),
        "mcleod_base_euler_is_expected": base_euler
        == Fraction(-341, 16_721_510_400),
        "euler_ratio_is_positive_integer": index.denominator == 1 and index > 0,
        "index_is_expected": index == 75_202_560,
        "volume_coefficient_is_expected": volume_coefficient
        == Fraction(46_376, 893_025),
        "full_face_vector_reconstructs": reconstructed_f_vector
        == expected_f_vector,
        "mirror_type_histogram_is_expected": mirrors[
            "half_scaled_square_divisibility_histogram"
        ]
        == {"-6,3": 80, "-2,1": 36, "-1,1": 80},
    }
    if not all(checks.values()):
        raise ArithmeticError(f"McLeod-index checks failed: {checks}")

    return {
        "schema": "ss_fermat_quartic_char3_c2_h10_mcleod_index_v1",
        "sources": {
            "coxeter_screen": str(SCREEN.relative_to(ROOT)),
            "coxeter_screen_sha256": _sha256(SCREEN),
            "connected_component_audit": str(COMPONENTS.relative_to(ROOT)),
            "connected_component_audit_sha256": _sha256(COMPONENTS),
            "literal_vinberg_certificate": str(VINBERG.relative_to(ROOT)),
            "literal_vinberg_certificate_sha256": _sha256(VINBERG),
            "exact_cone_certificate": str(CONE.relative_to(ROOT)),
            "exact_cone_certificate_sha256": _sha256(CONE),
            "fixed_lattice_and_mirror_certificate": str(FIXED.relative_to(ROOT)),
            "fixed_lattice_and_mirror_certificate_sha256": _sha256(FIXED),
            "mcleod_primary_source": "https://arxiv.org/abs/1007.2299",
            "coxiter_snapshot_relative_path": MCLEOD_RELATIVE.as_posix(),
            "coxiter_snapshot_sha256": _sha256(mcleod_path),
        },
        "method": "Solomon rational Euler sum over every finite special subgroup",
        "subgroup_196_wall_chamber": subgroup,
        "mcleod_f3_h10_14_wall_chamber": base,
        "integral_lineage": {
            "half_scaled_fixed_lattice": "<3> + <-1>^10",
            "explicit_unimodular_isometry": True,
            "primitive_root_type_histogram": mirrors[
                "half_scaled_square_divisibility_histogram"
            ],
            "all_196_reflections_integral": True,
            "reflection_subgroup_inclusion": "W_196 < W(f_3^10)",
        },
        "index": {
            "formula": "chi(W_196) / chi(W_f3_H10)",
            "value": index.numerator,
            "prime_factorization": _factorization(index.numerator),
        },
        "hyperbolic_covolume": {
            "gauss_bonnet_formula": "vol = -(32*pi^5/945)*chi in H10",
            "value": "pi^5 * 46376/893025",
        },
        "full_face_vector": reconstructed_f_vector,
        "lineage_decision": {
            "maximal_reflection_group": "McLeod f_3^10",
            "maximal_chamber_facet_count": 14,
            "subgroup_is_proper": True,
            "subgroup_index_is_exact": True,
            "literal_196_wall_subgroup_priority_remains_open": True,
            "new_lattice_or_commensurability_class": False,
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, default=CORPUS)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    result = build_certificate(arguments.corpus)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if (
            not arguments.output.is_file()
            or arguments.output.read_text(encoding="utf-8") != rendered
        ):
            raise SystemExit(f"archived artifact differs: {arguments.output}")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
    print(arguments.output)


if __name__ == "__main__":
    main()
