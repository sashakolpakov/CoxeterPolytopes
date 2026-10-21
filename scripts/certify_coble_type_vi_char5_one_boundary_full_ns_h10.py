#!/usr/bin/env python3
"""Certify the characteristic-five one-boundary Type-VI full-nef chamber.

Kondo's Type-VI Coble surface supplies the complete twenty-root 20E
Coble--Mukai chamber.  Its 57 primitive cusps lift bijectively to all
(-1)-curves by E=f-b, where b=B/2=-K.  Together with the boundary this gives
78 primitive full-Neron--Severi walls in I_(1,10).

The finite-volume gate is exact: GMP cddlib supplies the complete rational
ray/incidence/adjacency data; a compiled exhaustive classifier checks every
rank-nine and rank-ten pair-compatible Coxeter mask; and all connected and
disconnected affine cusp neighbourhoods are checked literally.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import shutil
from typing import Any, Iterator

import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form
from sympy.polys.domains import ZZ

from scripts.certify_coble_one_boundary_full_ns_h10 import (
    _calibration_walls,
    _load_templates,
    _root_span_and_reflections,
)
from scripts.certify_jacobian_kummer_h10_reflection_index import (
    _euler_from_connected_components,
)
from scripts.certify_shimada_disc16_d5_h10_lineage import (
    _connected_elliptic_masks,
)
import scripts.certify_coble_type_vii_char5_one_boundary_full_ns_h10 as engine


ROOT = Path(__file__).resolve().parents[1]
SUMMARY_OUTPUT = (
    ROOT / "data/enumerations/coble_type_vi_char5_one_boundary_full_ns_h10.json"
)
RAY_OUTPUT = (
    ROOT
    / "data/enumerations/coble_type_vi_char5_one_boundary_full_ns_h10_rays.jsonl.gz"
)
EDGE_OUTPUT = (
    ROOT
    / "data/enumerations/coble_type_vi_char5_one_boundary_full_ns_h10_edges.jsonl.gz"
)
ENRIQUES_SUMMARY = (
    ROOT / "data/known_examples/enriques_finite_automorphism_nef_chambers.json"
)
ENRIQUES_RAYS = (
    ROOT / "data/enumerations/enriques_finite_automorphism_nef_chamber_rays.jsonl.gz"
)
DEFAULT_COXITER = Path("/private/tmp/coxiter-current/CoxIter-master")

SCHEMA = "coxeter.coble_type_vi_char5_one_boundary_full_ns_h10.v1"
DIMENSION = 10
FACET_COUNT = 78
EXPECTED_RAY_HASH = "58a0fbb0fccc636d5544e6aaa5c6e665a67f36eb1be01d2a4c06bcaffcddd903"
EXPECTED_CLIQUES = {9: 26_769_085, 10: 40_619_993}
EXPECTED_FINITE_DIAGRAMS = {9: 166_320, 10: 21_661}
P10_EULER = Fraction(-1, 11_147_673_600)
EXPECTED_EULER = Fraction(-527, 6_480)
EXPECTED_P10_INDEX = 906_608_640
H10_VOLUME_COEFFICIENT = Fraction(1_054, 382_725)


@contextmanager
def _type_vi_engine_context() -> Iterator[None]:
    """Parameterize the generic Type-VII face helpers for 78 walls."""

    previous_facet_count = engine.FACET_COUNT
    previous_dimension = engine.DIMENSION
    engine.FACET_COUNT = FACET_COUNT
    engine.DIMENSION = DIMENSION
    try:
        yield
    finally:
        engine.FACET_COUNT = previous_facet_count
        engine.DIMENSION = previous_dimension


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _matrix_rows(matrix: sp.Matrix) -> list[list[int]]:
    if any(value.is_Integer is not True for value in matrix):
        raise ArithmeticError("an expected integral matrix changed")
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _type_vi_integral_marking_bridge(
    e10: sp.Matrix, cm_root_matrix: sp.Matrix
) -> dict[str, Any]:
    """Prove that the labelled Gamma_VI roots themselves form E10."""

    smith = smith_normal_form(cm_root_matrix, domain=ZZ)
    smith_diagonal = [abs(int(smith[index, index])) for index in range(10)]
    hnf_basis = hermite_normal_form(cm_root_matrix)
    root_gram = hnf_basis.T * e10 * hnf_basis
    root_gram_smith = smith_normal_form(root_gram, domain=ZZ)
    root_gram_smith_diagonal = [
        abs(int(root_gram_smith[index, index])) for index in range(10)
    ]
    if not (
        cm_root_matrix.rank() == 10
        and smith_diagonal == [1] * 10
        and hnf_basis.shape == (10, 10)
        and abs(int(hnf_basis.det())) == 1
        and root_gram.det() == -1
        and root_gram_smith_diagonal == [1] * 10
    ):
        raise ArithmeticError("the intrinsic Type-VI E10 marking bridge changed")
    return {
        "root_coordinate_span_smith_diagonal": smith_diagonal,
        "root_coordinate_span_is_all_E10": True,
        "HNF_basis_in_E10": _matrix_rows(hnf_basis),
        "HNF_basis_determinant": int(hnf_basis.det()),
        "root_lattice_gram": _matrix_rows(root_gram),
        "root_lattice_determinant": int(root_gram.det()),
        "root_lattice_gram_smith_diagonal": root_gram_smith_diagonal,
        "transport_conclusion": (
            "the entrywise Kondo Gamma_VI to Shimada 20E root-Gram isometry "
            "already spans and therefore extends to an integral E10 isometry"
        ),
        "graph_only_match_used": False,
    }


def build_certificate(
    cdd_executable: str | None = None,
    compiler: str | None = None,
    coxiter_corpus: Path | None = DEFAULT_COXITER,
) -> tuple[dict[str, Any], bytes, bytes]:
    cdd = cdd_executable or shutil.which("cddexec_gmp")
    if cdd is None:
        raise RuntimeError("cddexec_gmp is required for the exact face replay")

    chamber_data, ray_data = _load_templates()
    lattice, walls, construction = _calibration_walls(
        "VI", chamber_data, ray_data
    )
    if construction != {
        "cm_root_count": 20,
        "cm_cusp_count": 57,
        "wall_order": "b=-K; CM (-2)-roots; cusp lifts E=f-b",
    }:
        raise ArithmeticError("the pinned Type-VI CM seed changed")

    wall_matrix = sp.Matrix.hstack(*walls).T
    negative_gram = wall_matrix * lattice * wall_matrix.T
    type_vi_template = next(
        record
        for record in chamber_data["types"]
        if record["nikulin_kondo_type"] == "VI"
    )
    if negative_gram[1:21, 1:21] != sp.Matrix(
        type_vi_template["intersection_matrix"]
    ):
        raise ArithmeticError("the Type-VI Coble/20E root Gram bridge changed")
    cm_root_matrix = sp.Matrix.hstack(*walls[1:21])[:10, :]
    marking_bridge = _type_vi_integral_marking_bridge(
        lattice[:10, :10], cm_root_matrix
    )

    pivots = list(negative_gram.rref()[1])
    if len(pivots) != DIMENSION + 1 or len(walls) != FACET_COUNT:
        raise ArithmeticError("the Type-VI full wall rank/count changed")
    inequalities = negative_gram[:, pivots]
    basis_gram = negative_gram.extract(pivots, pivots)

    rays, incidence, adjacency = engine._run_cdd_repall(inequalities, cdd)
    positive_gram = [
        [-int(negative_gram[left, right]) for right in range(FACET_COUNT)]
        for left in range(FACET_COUNT)
    ]

    with _type_vi_engine_context():
        relations, numeric_relations = engine._relation_payload(positive_gram)
        connected_finite, connected_search = _connected_elliptic_masks(
            positive_gram
        )
        solomon = _euler_from_connected_components(relations, connected_finite)
        euler = Fraction(solomon["rational_euler_characteristic"])
        if (
            euler != EXPECTED_EULER
            or euler / P10_EULER != EXPECTED_P10_INDEX
        ):
            raise ArithmeticError("the Solomon/P10 index calculation changed")

        temporary, diagram_executable, compiler_used = (
            engine._compile_diagram_enumerator(compiler)
        )
        try:
            finite_diagrams, clique_counts = engine._finite_diagrams(
                numeric_relations, diagram_executable
            )
        finally:
            temporary.cleanup()
        if clique_counts != EXPECTED_CLIQUES or {
            order: len(masks) for order, masks in finite_diagrams.items()
        } != EXPECTED_FINITE_DIAGRAMS:
            raise ArithmeticError("the exhaustive finite-diagram census changed")
        if solomon["finite_special_subset_count_by_rank"][9:] != [
            EXPECTED_FINITE_DIAGRAMS[9],
            EXPECTED_FINITE_DIAGRAMS[10],
        ]:
            raise ArithmeticError("the two finite-diagram engines disagree")

        face_payload, ray_records, edge_records = engine._face_and_vinberg_audit(
            rays,
            incidence,
            adjacency,
            inequalities,
            basis_gram,
            positive_gram,
            numeric_relations,
            finite_diagrams,
            connected_finite,
        )
        bounded_lineage = engine._bounded_coxiter_header_audit(coxiter_corpus)
        if bounded_lineage.get("available"):
            # The reusable Type-VII scanner correctly sees the temporarily
            # parameterized 78-wall target, but its display keys contain the
            # old target's literal label.  Normalize only the labels here.
            bounded_lineage["ordinary_H10_78_header_matches"] = (
                bounded_lineage.pop("ordinary_H10_72_header_matches")
            )
            bounded_lineage["doubling_H10_78_header_matches"] = (
                bounded_lineage.pop("doubling_H10_72_header_matches")
            )

    cone = {
        key: value
        for key, value in face_payload.items()
        if key != "literal_vinberg"
    }
    literal = face_payload["literal_vinberg"]
    literal["all_pair_compatible_rank_9_cliques"] = clique_counts[9]
    literal["all_pair_compatible_rank_10_cliques"] = clique_counts[10]
    arithmetic = _root_span_and_reflections(lattice, walls)
    arithmetic["norm_histogram"] = {
        str(norm): count
        for norm, count in arithmetic["norm_histogram"].items()
    }

    if not (
        cone["extreme_ray_count"] == 23_458
        and cone["finite_vertex_count"] == 21_661
        and cone["ideal_vertex_count"] == 1_797
        and cone["hyperideal_vertex_count"] == 0
        and cone["all_facets_irredundant"]
        and cone["all_extreme_rays_on_one_lorentz_time_sheet"]
        and cone["primitive_ray_list_sha256"] == EXPECTED_RAY_HASH
        and literal["exact_geometric_edge_count"] == 166_320
        and literal["connected_affine_component_count"] == 5_551
        and literal["distinct_full_parabolic_neighbourhood_count"] == 1_797
        and literal["bad_edge_count"] == 0
        and literal["bad_parabolic_neighbourhood_count"] == 0
        and arithmetic["roots_span_whole_lattice"]
        and arithmetic["all_reflections_integral_involutive_isometries"]
        and arithmetic["all_pairs_coxeter"]
    ):
        raise ArithmeticError("the Type-VI full-NS chamber census changed")

    ray_header = {
        "schema": "coxeter.coble_type_vi_char5_one_boundary_full_ns_h10_rays.v1",
        "record_format": {
            "ray_index": "one based",
            "primitive_coordinate_basis_vector": (
                "coordinates in the 11 pivot wall normals"
            ),
            "active_facets": "one-based complete incidence set",
            "square": "exact square in the pivot-wall Gram",
            "kind": "finite, ideal, or hyperideal",
        },
        "coordinate_basis_facets": [index + 1 for index in pivots],
        "ray_count": len(ray_records),
    }
    edge_header = {
        "schema": "coxeter.coble_type_vi_char5_one_boundary_full_ns_h10_edges.v1",
        "record_format": (
            "[edge_facets, endpoint_ray_indices, endpoint_kinds, finite_components]"
        ),
        "indexing": "one based",
        "edge_count": len(edge_records),
    }
    ray_archive = engine._gzip_jsonl(ray_header, ray_records)
    edge_archive = engine._gzip_jsonl(edge_header, edge_records)

    payload = {
        "schema": SCHEMA,
        "status": (
            "proved exact complete-nef finite-volume noncompact arithmetic "
            "Coxeter chamber"
        ),
        "direction": "SURFACE_TO_POLYTOPE",
        "inverse_carrier_score": {
            "level": "N/A/L0",
            "reason": (
                "the construction starts from Kondo's Type-VI surface and "
                "makes no polytope-to-surface implication"
            ),
        },
        "surface": {
            "classification": "Kondo Type VI strong Coble surface",
            "base_field": "algebraically closed, characteristic 5",
            "uniqueness": "unique isomorphism class in Kondo Theorem 1.1/Table 2",
            "boundary_components": 1,
            "automorphism_group": "S5",
            "R_invariant": "(E6+A4,{0})",
            "effective_irreducible_CM_roots": 20,
            "construction": "Kondo Example 4.12 and Proposition 4.13",
        },
        "algebro_geometric_completeness": {
            "intrinsic_full_E10_marking_bridge": marking_bridge,
            "entrywise_Coble_Gamma_VI_equals_Shimada_20E_Gram": True,
            "negative_curve_trichotomy": (
                "boundary curves, (-2)-curves, and (-1)-curves are the only "
                "negative curves"
            ),
            "complete_minus_two_curves": (
                "the 20 Type-VI roots are all effective irreducible roots by "
                "Kondo Proposition 3.7"
            ),
            "one_boundary_cusp_lift": (
                "every primitive CM cusp f gives the unique (-1)-curve E=f-b "
                "in B+2E=2f, and every (-1)-curve returns f=E+b"
            ),
            "complete_minus_one_curves": 57,
            "full_nef_wall_formula": (
                "{b=B/2=-K} union {20 CM (-2)-curves} union "
                "{E=f-b for all 57 CM cusps}"
            ),
            "why_complete": (
                "the three negative-curve classes are exhausted, and "
                "nonnegative-square irreducible curves add no hyperbolic nef facet"
            ),
        },
        "polytope": {
            "id": "COBLE-VI-p5-b1-H10-78",
            "dimension": DIMENSION,
            "facet_count": FACET_COUNT,
            "classification": "finite_volume_noncompact",
            "compact": False,
            "finite_volume": True,
            "arithmetic": True,
            "lattice": "I_(1,10) = E10 orthogonal_sum <-1>",
            "lattice_gram": _matrix_rows(lattice),
            "coordinate_basis_facets": [index + 1 for index in pivots],
            "wall_order": "b=-K; 20 Type-VI CM roots; 57 cusp lifts E=f-b",
            "wall_vectors": [[int(value) for value in wall] for wall in walls],
            "negative_wall_gram": _matrix_rows(negative_gram),
            "normalized_gram_sha256": hashlib.sha256(
                json.dumps(
                    _matrix_rows(negative_gram), separators=(",", ":")
                ).encode()
            ).hexdigest(),
            "wall_curve_count": {
                "boundary_B_square_minus_4_primitive_normal_b_square_minus_1": 1,
                "minus_2_curves": 20,
                "minus_1_curves": 57,
            },
            "arithmetic_checks": arithmetic,
            "exact_double_description": cone,
            "literal_vinberg": literal,
            "solomon": {
                "connected_elliptic_count": len(connected_finite),
                "connected_elliptic_search": connected_search,
                "connected_elliptic_rank_histogram": solomon[
                    "connected_elliptic_rank_histogram"
                ],
                "finite_special_subset_count_by_rank": solomon[
                    "finite_special_subset_count_by_rank"
                ],
                "rational_coxeter_euler_characteristic": str(EXPECTED_EULER),
                "P10_rational_euler_characteristic": str(P10_EULER),
                "index_in_W_P10": EXPECTED_P10_INDEX,
                "H10_Gauss_Bonnet_volume": (
                    f"{H10_VOLUME_COEFFICIENT.numerator}*pi^5/"
                    f"{H10_VOLUME_COEFFICIENT.denominator}"
                ),
            },
        },
        "method": {
            "wall_source": (
                "exact Shimada 20E Type-VI E10 vectors and their exact 57 CM "
                "ideal rays, transported by the full root-span marking theorem"
            ),
            "cone": (
                "complete GMP rational double description with incidence, "
                "adjacency, and a common Lorentz time-sheet check"
            ),
            "finite_diagrams": (
                "all 26,769,085 pair-compatible rank-9 and 40,619,993 rank-10 "
                "masks classified componentwise by the finite crystallographic list"
            ),
            "affine_diagrams": (
                "all connected affine extensions through order 12 checked by "
                "exact determinant with a positive-definite principal deletion"
            ),
            "diagram_enumerator_source": str(
                engine.DIAGRAM_SOURCE.relative_to(ROOT)
            ),
            "diagram_enumerator_source_sha256": _sha256(engine.DIAGRAM_SOURCE),
            "compiler": compiler_used,
            "blind_binomial_subset_scan_used": False,
        },
        "lineage": {
            "surface_and_20_root_CM_H9_chamber": (
                "known: Kondo Type VI / Shimada 20E"
            ),
            "H9_seed_cusp_count": (
                "57 is known seed data and is not a novelty claim"
            ),
            "full_NS_H10_78_chamber": (
                "new exact positive chamber in this project ledger"
            ),
            "new_content": (
                "the characteristic-five one-boundary bijection of all 57 cusps "
                "with all 57 (-1)-curves, complete 78-wall full-nef Gram, exact "
                "face lattice, and literal Vinberg certificate"
            ),
            "P10_index_lineage": (
                "906608640 equals Shimada's published index for the Type-VI 20E "
                "H9 induced chamber; inherited equality is structural evidence, "
                "not novelty evidence"
            ),
            "worldwide_literal_abstract_priority": (
                "open; no worldwide novelty claim"
            ),
            "bounded_coxiter_corpus": bounded_lineage,
        },
        "sources": {
            "kondo": {
                "title": "Coble surfaces with finite automorphism group",
                "arxiv": "2102.00591",
                "used_results": [
                    "Theorem 1.1/Table 2",
                    "Proposition 3.1",
                    "Proposition 3.4",
                    "Proposition 3.7",
                    "Proposition 3.8",
                    "Lemma 3.11",
                    "Example 4.12",
                    "Proposition 4.13",
                ],
            },
            "shimada": {
                "record": "20E",
                "published_index": EXPECTED_P10_INDEX,
                "official_page": (
                    "https://www.math.sci.hiroshima-u.ac.jp/shimada/"
                    "K3andEnriques.html"
                ),
            },
            "vinberg": {
                "criterion": (
                    "every elliptic rank-(n-1) edge has exactly two finite/ideal "
                    "endpoints and every affine component has a maximal full "
                    "rank-(n-1) parabolic zero-neighbourhood"
                ),
                "reference": (
                    "E. B. Vinberg, Hyperbolic reflection groups, Russian Math. "
                    "Surveys 40:1 (1985), Proposition 4.2"
                ),
            },
            "pinned_input_sha256": {
                str(ENRIQUES_SUMMARY.relative_to(ROOT)): _sha256(
                    ENRIQUES_SUMMARY
                ),
                str(ENRIQUES_RAYS.relative_to(ROOT)): _sha256(ENRIQUES_RAYS),
            },
        },
        "artifacts": {
            "ray_archive": str(RAY_OUTPUT.relative_to(ROOT)),
            "ray_archive_sha256": hashlib.sha256(ray_archive).hexdigest(),
            "edge_archive": str(EDGE_OUTPUT.relative_to(ROOT)),
            "edge_archive_sha256": hashlib.sha256(edge_archive).hexdigest(),
        },
        "checks": {
            "complete_negative_wall_list": True,
            "all_78_facets_irredundant": True,
            "all_78_reflections_integral": True,
            "all_pairs_are_Coxeter": True,
            "all_166320_elliptic_edges_enumerated": True,
            "every_edge_has_exactly_two_endpoints": True,
            "all_5551_connected_affine_components_enumerated": True,
            "every_affine_component_completes_to_a_rank_9_cusp": True,
            "all_1797_ideal_vertices_are_full_rank_9_parabolic_neighbourhoods": True,
            "finite_volume": True,
            "noncompact": True,
        },
        "all_checks_pass": True,
    }
    return payload, ray_archive, edge_archive


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cdd", default=None)
    parser.add_argument("--compiler", default=None)
    parser.add_argument("--coxiter-corpus", type=Path, default=DEFAULT_COXITER)
    parser.add_argument("--without-coxiter-corpus", action="store_true")
    parser.add_argument("--summary-output", type=Path, default=SUMMARY_OUTPUT)
    parser.add_argument("--ray-output", type=Path, default=RAY_OUTPUT)
    parser.add_argument("--edge-output", type=Path, default=EDGE_OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    corpus = None if arguments.without_coxiter_corpus else arguments.coxiter_corpus
    summary, rays, edges = build_certificate(
        arguments.cdd, arguments.compiler, corpus
    )
    engine._write_or_check_text(
        arguments.summary_output,
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        arguments.check_archive,
    )
    engine._write_or_check_bytes(
        arguments.ray_output, rays, arguments.check_archive
    )
    engine._write_or_check_bytes(
        arguments.edge_output, edges, arguments.check_archive
    )
    print(
        json.dumps(
            {
                "id": summary["polytope"]["id"],
                "vertices": {
                    "finite": summary["polytope"]["exact_double_description"][
                        "finite_vertex_count"
                    ],
                    "ideal": summary["polytope"]["exact_double_description"][
                        "ideal_vertex_count"
                    ],
                },
                "edges": summary["polytope"]["literal_vinberg"][
                    "exact_geometric_edge_count"
                ],
                "finite_volume": summary["checks"]["finite_volume"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
