#!/usr/bin/env python3
"""Certify a full-Picard complex-K3 carrier for the JK H12/172 chamber.

The 172-facet chamber, its saturated rank-13 fixed divisor lattice, and two
independent finite-volume proofs are already archived.  This verifier closes
the remaining carrier step.  It reconstructs all primitive facet normals in
the saturated fixed-lattice basis, orients them by the published ample class,
checks every reflection on the full lattice, computes the wall-span index,
and applies Nikulin's primitive-embedding stable range to the exact
discriminant form.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

from scripts.certify_jacobian_kummer_h12_fixed_lattice import (
    _fixed_basis,
    _primitive_integer_column,
)
from scripts.certify_three_h11_full_picard_k3_carriers import (
    _carrier_lattice,
    _payload_sha256,
    _stable_range_embedding,
    _wall_system,
)
from scripts.search_jacobian_kummer_chamber_folds import _parse_source


ROOT = Path(__file__).resolve().parents[1]
PUBLISHED_DATA = (
    ROOT
    / "data/external/dolgachev_shimada_15_nodal_2019/15nodalcompdata.txt"
)
CHAMBER = ROOT / "data/enumerations/jacobian_kummer_h12_interim_facet_gram.json"
FIXED_LATTICE = (
    ROOT / "data/enumerations/jacobian_kummer_h12_fixed_lattice_certificate.json"
)
VINBERG = (
    ROOT / "data/enumerations/jacobian_kummer_h12_vinberg_volume_certificate.json"
)
EXACT_CONE = (
    ROOT / "data/enumerations/jacobian_kummer_h12_exact_cone_certificate.json"
)
OUTPUT = (
    ROOT / "data/enumerations/jacobian_kummer_h12_full_picard_k3_carrier.json"
)

SCHEMA = "coxeter.jacobian_kummer_h12.full_picard_complex_k3_carrier.v1"
EXPECTED_SHA256 = {
    "published_data": "413d7363faf4574595e9c28329737eb4ad8010e0aebae53afb76762ce9f8966a",
    "chamber": "53c39642b81bbd7f5be172408bfecb857bd4ed57467c91bdbaed63f9e734fb92",
    "fixed_lattice": "208059cf885b8a3e56d57270af2595dde06b66b478c870af5b08b93b9bea2614",
    "vinberg": "c3474e33120cab541205225dd8df31beea7ae8006cd1667c4fd66408f9b931f1",
    "exact_cone": "84d9660b5990d89fb4a78b73180fa815843ee7af46d4d8c56be73f51e37a7c0f",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _oriented_wall_data(
    chamber: dict[str, Any],
    source: dict[str, Any],
    fixed_gram: sp.Matrix,
    basis: sp.Matrix,
) -> tuple[dict[str, Any], sp.Matrix]:
    ambient_gram: sp.Matrix = source["gram"]
    ample: sp.Matrix = source["ample"]
    ample_coordinates = sp.simplify(
        fixed_gram.inv() * basis * ambient_gram * ample.T
    )
    if not all(value.is_Integer for value in ample_coordinates):
        raise ArithmeticError("published ample class is not integral in fixed lattice")

    published_walls = tuple(
        sp.Matrix(1, ambient_gram.rows, vector)
        for vector in source["wall_vectors"]
    )
    roots: list[list[int]] = []
    for wall_orbit in chamber["facet_source_wall_orbits"]:
        orbit_sum = sum(
            (published_walls[int(index) - 1] for index in wall_orbit),
            sp.zeros(1, ambient_gram.rows),
        )
        coordinates = sp.simplify(
            fixed_gram.inv() * basis * ambient_gram * orbit_sum.T
        )
        primitive = _primitive_integer_column(coordinates)
        roots.append([int(value) for value in primitive])

    walls = _wall_system(
        fixed_gram,
        roots,
        ample_coordinates,
        {"-2,1": 16, "-4,2": 36, "-4,4": 96, "-8,4": 24},
    )
    return walls, ample_coordinates


def build_certificate() -> dict[str, Any]:
    paths = {
        "published_data": PUBLISHED_DATA,
        "chamber": CHAMBER,
        "fixed_lattice": FIXED_LATTICE,
        "vinberg": VINBERG,
        "exact_cone": EXACT_CONE,
    }
    for key, path in paths.items():
        if _sha256(path) != EXPECTED_SHA256[key]:
            raise ArithmeticError(f"pinned {key} input changed")

    chamber = json.loads(CHAMBER.read_text(encoding="utf-8"))
    fixed = json.loads(FIXED_LATTICE.read_text(encoding="utf-8"))
    vinberg = json.loads(VINBERG.read_text(encoding="utf-8"))
    cone = json.loads(EXACT_CONE.read_text(encoding="utf-8"))
    source = _parse_source(PUBLISHED_DATA)

    involution = sp.Matrix(chamber["group_class"]["representative_matrix"])
    basis = _fixed_basis(involution)
    fixed_gram = sp.simplify(basis * source["gram"] * basis.T)
    if fixed_gram != sp.Matrix(fixed["fixed_lattice"]["gram"]):
        raise ArithmeticError("reconstructed fixed Gram differs from archive")

    embedding = _stable_range_embedding(fixed_gram)
    carrier = _carrier_lattice(
        "U(2)+D5(-1)^2+<-4>", fixed_gram, embedding
    )
    carrier["canonical_discriminant_form"] = fixed["fixed_lattice"][
        "canonical_discriminant_form"
    ]
    carrier["source_occurrence"] = {
        "surface": "general genus-2 Jacobian--Kummer K3 surface Y16",
        "role": "saturated fixed divisor sublattice inside a rank-16 Picard lattice",
        "new_gain": (
            "a separate seven-dimensional characteristic-zero K3 family "
            "with this rank-13 carrier as the full Picard lattice"
        ),
    }

    walls, ample_coordinates = _oriented_wall_data(
        chamber, source, fixed_gram, basis
    )
    literal_vinberg = {
        "method": vinberg["method"],
        "classification": vinberg["vinberg_decision"],
        "elliptic_rank11_edges": vinberg["elliptic_rank_11_edge_count"],
        "endpoint_count_histogram": vinberg["endpoint_count_histogram"],
        "finite_vertices": vinberg["finite_vertex_count"],
        "cusps": vinberg["ideal_vertex_count"],
        "bad_edges": vinberg["bad_edge_count"],
        "unrecognized_endpoint_edges": vinberg["unrecognized_endpoint_edge_count"],
        "invalid_full_zero_neighborhoods": vinberg[
            "invalid_full_zero_neighborhood_count"
        ],
    }
    exact = cone["exact_double_description"]
    exact_cone = {
        "method": exact["engine"],
        "classification": cone["decision"]["classification"],
        "extreme_rays": exact["extreme_ray_count"],
        "finite_vertices": exact["finite_vertex_count"],
        "cusps": exact["ideal_vertex_count"],
        "hyperideal_rays": exact["hyperideal_vertex_count"],
        "all_exact_replay_checks": cone["checks"],
    }

    checks = {
        "upstream_metric_certificate_green": all(chamber["checks"].values()),
        "upstream_fixed_lattice_certificate_green": (
            fixed["all_checks_pass"] and all(fixed["checks"].values())
        ),
        "fixed_gram_replayed_exactly": fixed_gram
        == sp.Matrix(fixed["fixed_lattice"]["gram"]),
        "published_ample_replayed_in_full_carrier": (
            [int(value) for value in ample_coordinates]
            == walls["strict_interior_coordinates"]
        ),
        "Nikulin_stable_range_primitive_embedding": (
            embedding["primitive_embedding_exists"]
            and embedding["inequality"] == "9 >= 5+2"
        ),
        "all_172_walls_replayed": len(walls["wall_records"]) == 172,
        "all_172_reflections_integral_on_full_carrier": all(
            all(record["checks"].values()) for record in walls["wall_records"]
        ),
        "walls_generate_full_Picard_lattice": (
            walls["wall_span"]["rank"] == 13
            and walls["wall_span"]["index_in_carrier_lattice"] == 1
        ),
        "literal_Vinberg_gate_green": (
            literal_vinberg["classification"] == "finite_volume_noncompact"
            and literal_vinberg["bad_edges"] == 0
            and literal_vinberg["unrecognized_endpoint_edges"] == 0
            and literal_vinberg["invalid_full_zero_neighborhoods"] == 0
        ),
        "exact_GMP_gate_green": (
            exact_cone["classification"] == "finite_volume_noncompact"
            and exact_cone["hyperideal_rays"] == 0
            and all(exact_cone["all_exact_replay_checks"].values())
        ),
        "independent_volume_counts_agree": (
            literal_vinberg["finite_vertices"] == exact_cone["finite_vertices"]
            and literal_vinberg["cusps"] == exact_cone["cusps"]
        ),
    }

    payload = {
        "schema": SCHEMA,
        "status_date": "2026-08-07",
        "sources": {
            key: {
                "path": str(path.relative_to(ROOT)),
                "sha256": EXPECTED_SHA256[key],
            }
            for key, path in paths.items()
        },
        "theorem": (
            "The certified finite-volume noncompact JK-H12/172 chamber has "
            "a seven-dimensional very-general complex-projective K3 family "
            "with full Picard lattice U(2)+D5(-1)^2+<-4>; all 172 oriented "
            "wall normals generate that Picard lattice and all 172 "
            "reflections act integrally on it."
        ),
        "chamber": {
            "id": "JK-H12-172",
            "dimension": 12,
            "facet_count": 172,
            "carrier": carrier,
            "certified_reflection_stable_overclosures": [
                {
                    "index": 1,
                    "carrier": "original saturated fixed lattice is already reflection-stable",
                }
            ],
            "oriented_wall_system": walls,
            "finite_volume": {
                "literal_Vinberg": literal_vinberg,
                "independent_exact_GMP_cone": exact_cone,
            },
            "status": "new very-general characteristic-zero full-Picard K3 carrier",
        },
        "scope": {
            "positive_gain": (
                "the known chamber's fixed divisor lattice is promoted from "
                "a sublattice occurrence on Y16 to the full Picard lattice "
                "of a seven-dimensional K3 family"
            ),
            "not_new_abstract_polytope": True,
            "not_first_algebraic_occurrence": True,
            "not_full_nef": (
                "square-minus-four and square-minus-eight walls prevent a "
                "complete K3 nef-wall interpretation"
            ),
            "not_surface_automorphisms": (
                "integral Picard reflections are not asserted to extend to "
                "surface automorphisms or Hodge isometries"
            ),
            "worldwide_priority_claim": "none",
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate()
    if not payload["all_checks_pass"]:
        raise ArithmeticError("JK-H12 full-Picard K3 carrier certificate failed")
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.exists() or args.output.read_text(encoding="utf-8") != rendered:
            raise ArithmeticError("archived JK-H12 full-Picard carrier changed")
        print(f"verified {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(f"wrote {args.output}")
    print(f"payload_sha256={_payload_sha256(payload)}")


if __name__ == "__main__":
    main()
