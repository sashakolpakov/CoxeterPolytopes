#!/usr/bin/env python3
"""Certify full-domain discrete Coxeter labels for H9/14 records 30 and 42.

This certificate is deliberately local/combinatorial.  It checks every
rank-nine vertex diagram and every higher minimal nonface, but assigns no
dotted Gram entry and makes no global inertia or rank claim.
"""

from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
import json
from math import comb
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import sympy as sp

from coxeter_polytopes.combinatorics.facet_completion import (
    FULL_VERTEX_MASK,
    explicit_product_facet_selection,
    generated_face_table,
    minimal_nonfaces_from_table,
)
from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
    FINITE_LABELS,
    audit_h9_14_coxeter_labels,
    is_finite_coxeter_labeling,
    lanner_isomorphism_codes,
    lanner_placement_codes,
    minimal_nonface_component_partition,
    minimum_cross_component_nonright_edges,
    pair,
    representative_lanner_labels,
)


ROOT = Path(__file__).resolve().parents[1]
REFINEMENT = ROOT / "data/enumerations/h9_14_facet_refinement.json"
LABEL_MODULE = (
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_coxeter_labels.py"
)
OUTPUT = ROOT / "data/enumerations/h9_14_product_coxeter_labels.json"
EXPECTED_CERTIFICATE_SHA256 = (
    "743bb6f42588a35dc59293d19e7d27e762ce7ae392497504bc652cd026c6167a"
)


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_sha256(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _certificate_digest(manifest: Mapping[str, object]) -> str:
    critical = dict(manifest)
    critical.pop("certificate_sha256", None)
    return _canonical_sha256(critical)


def _vertices(mask: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(14) if mask & (1 << vertex))


def _label_map(order: int, code: Sequence[int]) -> dict[tuple[int, int], int]:
    return {
        pair(left, right): int(label)
        for (left, right), label in zip(
            combinations(range(order), 2), code, strict=True
        )
    }


def _lanner_gram(order: int, code: Sequence[int]) -> sp.Matrix:
    gram = 2 * sp.eye(order)
    for (left, right), label in zip(
        combinations(range(order), 2), code, strict=True
    ):
        gram[left, right] = gram[right, left] = -2 * sp.cos(sp.pi / label)
    return gram


@lru_cache(maxsize=None)
def _exhaustive_lanner_placement_codes(
    order: int,
) -> tuple[tuple[tuple[int, ...], ...], dict[str, int]]:
    """Generate the Lanner table from all ``4 choose pair`` label words.

    We first restrict one fixed node deletion to the complete finite Coxeter
    list; this is an exhaustive pruning, not an ansatz.  Every remaining
    extension is tested on all one-node deletions.  Interlacing then leaves
    only the exact determinant alternatives: negative is Lanner and zero is
    affine.
    """

    if order not in (4, 5):
        raise ValueError("only the H9/14 Lanner orders four and five are used")
    pairs = tuple(combinations(range(order), 2))
    prefix_pairs = tuple(combinations(range(order - 1), 2))
    finite_prefixes: list[dict[tuple[int, int], int]] = []
    for prefix_code in product(FINITE_LABELS, repeat=len(prefix_pairs)):
        prefix = dict(zip(prefix_pairs, prefix_code, strict=True))
        if is_finite_coxeter_labeling(range(order - 1), prefix):
            finite_prefixes.append(prefix)

    deletion_sets = tuple(
        tuple(vertex for vertex in range(order) if vertex != deleted)
        for deleted in range(order)
    )
    lanner_codes = []
    minimally_nonfinite = 0
    affine = 0
    for prefix in finite_prefixes:
        for extension in product(FINITE_LABELS, repeat=order - 1):
            labels = dict(prefix)
            labels.update(
                {(vertex, order - 1): label for vertex, label in enumerate(extension)}
            )
            if not all(
                is_finite_coxeter_labeling(deletion, labels)
                for deletion in deletion_sets
            ):
                continue
            if is_finite_coxeter_labeling(range(order), labels):
                continue
            minimally_nonfinite += 1
            code = tuple(labels[edge] for edge in pairs)
            determinant = sp.expand(_lanner_gram(order, code).det())
            if determinant.is_negative is True:
                lanner_codes.append(code)
            elif determinant.is_zero is True:
                affine += 1
            else:
                raise AssertionError(
                    "a minimally nonfinite word was neither Lanner nor affine"
                )
    return tuple(sorted(lanner_codes)), {
        "all_label_words": len(FINITE_LABELS) ** len(pairs),
        "finite_fixed_deletion_prefixes": len(finite_prefixes),
        "minimally_nonfinite_words": minimally_nonfinite,
        "affine_words": affine,
        "lanner_words": len(lanner_codes),
    }


def _lanner_catalogue_record(order: int) -> dict[str, object]:
    types = lanner_isomorphism_codes(order)
    exhaustive_codes, exhaustive_counts = _exhaustive_lanner_placement_codes(order)
    if exhaustive_codes != lanner_placement_codes(order):
        raise AssertionError(
            f"order-{order} hard-coded Lanner orbit table is not exhaustive"
        )
    determinant_records = []
    for index, code in enumerate(types):
        labels = _label_map(order, code)
        deletion_finite = [
            is_finite_coxeter_labeling(
                (vertex for vertex in range(order) if vertex != deleted), labels
            )
            for deleted in range(order)
        ]
        gram = _lanner_gram(order, code)
        determinant = sp.factor(
            gram.det(), extension=[sp.sqrt(2), sp.sqrt(5)]
        )
        if determinant.is_negative is not True or not all(deletion_finite):
            raise AssertionError(f"order-{order} Lanner type {index} failed replay")
        determinant_records.append(
            {
                "type_index": index,
                "canonical_pair_code": list(code),
                "full_gram_determinant": sp.sstr(determinant),
                "full_gram_determinant_is_exactly_negative": True,
                "every_one_node_deletion_is_finite": True,
            }
        )
    return {
        "order": order,
        "isomorphism_type_count": len(types),
        "labelled_placement_count": len(lanner_placement_codes(order)),
        "complete_generation_audit": {
            **exhaustive_counts,
            "generated_codes_equal_orbit_table": True,
            "method": (
                "Exhaust all label words in {2,3,4,5}; prune only by one "
                "finite node deletion; test every node deletion with the "
                "complete finite Coxeter classifier; separate Lanner from "
                "affine by the exact Gram determinant."
            ),
        },
        "types": determinant_records,
        "inertia_argument": (
            "Every one-node principal submatrix is positive definite by the "
            "finite Coxeter classification. Interlacing and the negative, "
            "nonzero full determinant give inertia (order-1,1,0)."
        ),
    }


def _bracelet_code(colors: Sequence[int]) -> tuple[int, ...]:
    word = tuple(colors)
    rotations = [word[index:] + word[:index] for index in range(len(word))]
    reverse = tuple(reversed(word))
    rotations.extend(
        reverse[index:] + reverse[:index] for index in range(len(word))
    )
    return min(rotations)


def _five_cycle_bracelet_count() -> int:
    bracelets = {
        _bracelet_code(colors) for colors in product(FINITE_LABELS, repeat=5)
    }
    if len(bracelets) != 136:
        raise AssertionError("the four-color length-five bracelet count changed")
    return len(bracelets)


def _record_data(index: int) -> tuple[str, tuple[int, ...]]:
    archive = json.loads(REFINEMENT.read_text(encoding="utf-8"))
    record = archive["records"][index]
    if record["index"] != index or record["disposition"] != "explicit_product_witness":
        raise ValueError(f"record {index} is no longer an explicit product witness")
    graph6 = str(record["graph6"])
    selected = tuple(int(value) for value in record["selected_facet_complement_masks"])
    if explicit_product_facet_selection(graph6) != selected:
        raise ValueError(f"record {index} no longer matches its product construction")
    return graph6, selected


def _minimal_nonfaces(selected: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    masks_by_order = minimal_nonfaces_from_table(generated_face_table(selected))
    return tuple(
        _vertices(mask)
        for order in sorted(masks_by_order)
        for mask in masks_by_order[order]
    )


def _meeting_pairs(minimal_nonfaces: Iterable[Sequence[int]]) -> set[tuple[int, int]]:
    dotted = {pair(*face) for face in minimal_nonfaces if len(face) == 2}
    return set(combinations(range(14), 2)) - dotted


def _audit_record(index: int) -> dict[str, object]:
    graph6, selected = _record_data(index)
    vertices = tuple(_vertices(FULL_VERTEX_MASK ^ complement) for complement in selected)
    nonfaces = _minimal_nonfaces(selected)
    nonface_components = minimal_nonface_component_partition(nonfaces)
    minimum_cross_edges = minimum_cross_component_nonright_edges(nonfaces)
    by_order: dict[int, list[tuple[int, ...]]] = {}
    for face in nonfaces:
        by_order.setdefault(len(face), []).append(face)
    labels = {edge: 2 for edge in _meeting_pairs(nonfaces)}

    if index == 30:
        if Counter(map(len, nonfaces)) != {2: 3, 4: 2}:
            raise AssertionError("record 30 product factorization changed")
        type_indices = (0, 8)
        lanner_recipe = []
        for face, type_index in zip(by_order[4], type_indices, strict=True):
            labels.update(representative_lanner_labels(face, type_index))
            lanner_recipe.append(
                {"facets": list(face), "order": 4, "type_index": type_index}
            )
        cycle_recipe: list[dict[str, object]] = []
        orbit_parameterization = {
            "description": (
                "Choose an unordered pair, with repetition, of the nine "
                "order-four Lanner isomorphism types; make every cross-factor "
                "meeting pair right-angled."
            ),
            "raw_labelled_placements": 108**2,
            "isomorphism_orbits": comb(9 + 2 - 1, 2),
            "calculation": "binomial(9+2-1,2)=45",
            "nerve_automorphism_group_order": 55_296,
            "symmetry_quotient": (
                "The full minimal-nonface-hypergraph automorphism group is "
                "(S2 wr S3) x (S4 wr S2), of order 55296.  Internal S4 "
                "actions reduce each four-block to its Lanner isomorphism "
                "type and the outer S2 makes the pair of types unordered."
            ),
            "completeness_within_ansatz": (
                "Every orthogonal-join assignment is obtained uniquely up "
                "to nerve symmetry by this unordered pair of Lanner types; "
                "all remaining meeting labels are fixed to two."
            ),
            "global_gram_obstruction": {
                "status": "unsat",
                "proof": (
                    "The normalized Gram matrix is the orthogonal direct sum "
                    "of two Lanner blocks of inertia (3,1,0) and three dotted "
                    "two-facet blocks of inertia (1,1,0).  Hence every "
                    "legitimate dotted choice has inertia (9,5,0), not "
                    "(9,1,4)."
                ),
                "cross_factor_nonright_edge_is_necessary": True,
            },
        }
        expected_vertices = 128
        if len(nonface_components) != 5 or minimum_cross_edges != 4:
            raise AssertionError("record 30 minimal-nonface components changed")
    elif index == 42:
        if Counter(map(len, nonfaces)) != {2: 5, 4: 1, 5: 1}:
            raise AssertionError("record 42 product factorization changed")
        lanner_recipe = []
        for order, type_index in ((4, 7), (5, 4)):
            face = by_order[order][0]
            labels.update(representative_lanner_labels(face, type_index))
            lanner_recipe.append(
                {"facets": list(face), "order": order, "type_index": type_index}
            )
        polygon_facets = tuple(sorted({vertex for edge in by_order[2] for vertex in edge}))
        polygon_meeting_edges = tuple(
            edge
            for edge in combinations(polygon_facets, 2)
            if edge in labels
        )
        if len(polygon_facets) != 5 or len(polygon_meeting_edges) != 5:
            raise AssertionError("record 42 polygon factor changed")
        cycle_colors = (2, 3, 4, 5, 3)
        cycle_recipe = []
        for edge, label in zip(polygon_meeting_edges, cycle_colors, strict=True):
            labels[edge] = label
            cycle_recipe.append({"facets": list(edge), "label": label})
        bracelet_count = _five_cycle_bracelet_count()
        orbit_parameterization = {
            "description": (
                "Choose one of nine order-four Lanner types, one of five "
                "order-five types, and a four-coloring of the five polygon "
                "meeting edges modulo the dihedral group D5; make every "
                "other cross-factor meeting pair right-angled."
            ),
            "raw_labelled_placements": 108 * 420 * 4**5,
            "four_color_length_five_bracelets": bracelet_count,
            "bracelet_burnside_calculation": (
                "(4^5 + 4*4 + 5*4^3)/10 = 136"
            ),
            "isomorphism_orbits": 9 * 5 * bracelet_count,
            "calculation": "9*5*136=6120",
            "nerve_automorphism_group_order": 28_800,
            "symmetry_quotient": (
                "The full minimal-nonface-hypergraph automorphism group is "
                "D5 x S4 x S5, of order 28800.  The simplex-block actions "
                "leave only the 9 and 5 Lanner types; D5 leaves the 136 "
                "four-color bracelets on the polygon meeting-edge cycle."
            ),
            "completeness_within_ansatz": (
                "Every orthogonal-join assignment is obtained uniquely up "
                "to nerve symmetry by one order-four type, one order-five "
                "type, and one D5 bracelet; all remaining meeting labels are "
                "fixed to two."
            ),
            "global_gram_obstruction": {
                "status": "unsat",
                "proof": (
                    "The order-four and order-five Lanner principal blocks "
                    "are mutually orthogonal.  Each contributes one negative "
                    "direction, so the full Gram matrix has negative index at "
                    "least two for every choice of polygon dotted entries."
                ),
                "cross_factor_nonright_edge_is_necessary": True,
            },
        }
        expected_vertices = 100
        if len(nonface_components) != 3 or minimum_cross_edges != 2:
            raise AssertionError("record 42 minimal-nonface components changed")
    else:
        raise ValueError("this product certificate covers records 30 and 42")

    audit = audit_h9_14_coxeter_labels(vertices, nonfaces, labels)
    if len(vertices) != expected_vertices or not audit.passes_discrete_gate:
        raise AssertionError(f"record {index} representative failed the label gate")
    if orbit_parameterization["isomorphism_orbits"] not in (45, 6120):
        raise AssertionError("orthogonal-join orbit count changed")

    nonright = [
        [left, right, label]
        for (left, right), label in sorted(labels.items())
        if label > 2
    ]
    return {
        "record_index": index,
        "graph6": graph6,
        "ordinary_polytope": (
            "Delta_1^3 x Delta_3^2"
            if index == 30
            else "C_5 x Delta_3 x Delta_4"
        ),
        "rank_nine_vertex_count": len(vertices),
        "all_rank_nine_vertices_checked": True,
        "minimal_nonfaces": [list(face) for face in nonfaces],
        "minimal_nonface_order_histogram": {
            str(order): len(faces) for order, faces in sorted(by_order.items())
        },
        "representative_recipe": {
            "default_on_every_meeting_pair": 2,
            "lanner_blocks": lanner_recipe,
            "polygon_meeting_edge_overrides": cycle_recipe,
            "nonright_labelled_pairs": nonright,
        },
        "representative_audit": {
            "facet_count": audit.facet_count,
            "vertex_count": audit.vertex_count,
            "meeting_pair_count": audit.meeting_pair_count,
            "dotted_pair_count": audit.dotted_pair_count,
            "finite_label_histogram": {
                str(label): count for label, count in audit.finite_label_histogram
            },
            "bad_vertex_indices": list(audit.bad_vertex_indices),
            "bad_lanner_nonfaces": [list(face) for face in audit.bad_lanner_nonfaces],
            "missing_or_extra_label_pairs": [
                list(edge) for edge in audit.missing_or_extra_label_pairs
            ],
            "passes_discrete_gate": audit.passes_discrete_gate,
        },
        "orthogonal_join_subfamily": orbit_parameterization,
        "unrestricted_forced_cross_edge_gate": {
            "minimal_nonface_components": [
                list(component) for component in nonface_components
            ],
            "component_count": len(nonface_components),
            "minimum_nonright_meeting_edges_between_components": (
                minimum_cross_edges
            ),
            "proof": (
                "Every component contains an indefinite Lanner or dotted "
                "principal block.  With fewer than c-1 cross-component "
                "non-right edges, the nonzero Gram graph leaves at least two "
                "such blocks in distinct orthogonal components, forcing "
                "negative index at least two."
            ),
            "scope": "necessary for every global labeling, not only this ansatz",
        },
        "scope": {
            "discrete_label_status": "sat",
            "full_admissible_label_census": "not claimed",
            "dotted_entries_assigned": False,
            "global_gram_rank_and_inertia": (
                "unsat for the orthogonal-join subfamily; unrestricted "
                "cross-factor labels remain untested"
            ),
            "hyperbolic_coxeter_example": False,
        },
    }


def build() -> dict[str, object]:
    catalogues = [_lanner_catalogue_record(order) for order in (4, 5)]
    records = [_audit_record(index) for index in (30, 42)]
    if [record["orthogonal_join_subfamily"]["isomorphism_orbits"] for record in records] != [45, 6120]:
        raise AssertionError("product orbit census changed")
    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_full_domain_discrete_coxeter_label_gate",
        "claim": (
            "The ordinary-polytopal H9/14 records 30 and 42 pass the complete "
            "discrete finite/Lanner Coxeter-label gate over {2,3,4,5}.  "
            "Explicit orthogonal-join subfamilies have respectively 45 and "
            "6120 isomorphism orbits, and exact inertia arguments exclude "
            "every member of both subfamilies from the global Lorentzian gate. "
            "Every unrestricted candidate needs at least four and two "
            "cross-component non-right edges respectively."
        ),
        "finite_label_domain": {
            "labels": list(FINITE_LABELS),
            "reason_no_label_at_least_six": (
                "Burcroff, Lemma 10.10: a compact hyperbolic Coxeter "
                "9-polytope with 14 facets has no multi-multiple edge.  In "
                "Burcroff's terminology this excludes finite labels m>=6; "
                "the remaining finite labels are exactly 2,3,4,5."
            ),
            "reference": (
                "A. Burcroff, Compact hyperbolic Coxeter d-polytopes with "
                "d+4 facets and related dimension bounds, arXiv:2201.03437, "
                "Lemma 10.10"
            ),
            "golden_subalphabet_only": False,
        },
        "classification_basis": {
            "finite_vertex_diagrams": (
                "Complete irreducible finite Coxeter classification "
                "A,B,D,E,F,H,I2, applied componentwise without floating point."
            ),
            "lanner_diagrams": (
                "Complete classical Lanner classification in orders four and five."
            ),
            "catalogues": catalogues,
        },
        "source": {
            "facet_refinement_path": str(REFINEMENT.relative_to(ROOT)),
            "facet_refinement_sha256": _file_sha256(REFINEMENT),
            "label_verifier_path": str(LABEL_MODULE.relative_to(ROOT)),
            "label_verifier_sha256": _file_sha256(LABEL_MODULE),
        },
        "records": records,
        "global_scope": {
            "proved": (
                "Existence and exact symmetry-reduced counts of the displayed "
                "orthogonal-join discrete label families, followed by a "
                "global Gram obstruction for every member of those families; "
                "the minimal-nonface component theorem gives unrestricted "
                "cross-edge lower bounds four and two."
            ),
            "not_proved": [
                "A census of all discrete labelings of either nerve.",
                "Existence of compatible dotted Gram entries.",
                "A symmetric 14x14 Gram matrix of rank 10 and inertia (9,1,4).",
                "A new compact hyperbolic Coxeter polytope.",
            ],
        },
    }
    manifest["certificate_sha256"] = _certificate_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    if manifest.get("kind") != "h9_14_full_domain_discrete_coxeter_label_gate":
        raise ValueError("unexpected product label certificate kind")
    rebuilt = build()
    if manifest != rebuilt:
        raise ValueError("product label certificate did not reproduce")
    digest = _certificate_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("product label certificate digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != (
        EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("product label certificate is not frozen")
    return {
        "record_indices": [record["record_index"] for record in rebuilt["records"]],
        "discrete_label_statuses": [
            record["scope"]["discrete_label_status"] for record in rebuilt["records"]
        ],
        "orthogonal_join_isomorphism_orbits": [
            record["orthogonal_join_subfamily"]["isomorphism_orbits"]
            for record in rebuilt["records"]
        ],
        "certificate_sha256": digest,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--print-certificate",
        action="store_true",
        help="print a freshly rebuilt certificate instead of replaying the archive",
    )
    args = parser.parse_args(argv)
    if args.print_certificate:
        manifest = build()
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    print(json.dumps(verify(manifest), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
