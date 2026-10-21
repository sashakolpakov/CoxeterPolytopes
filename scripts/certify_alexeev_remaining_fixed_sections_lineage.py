#!/usr/bin/env python3
"""Exact bounded lineage audit for six remaining Alexeev fixed sections.

The certificate separates five questions which are easy to conflate:

* rational similarity of the fixed lattices;
* equality between the six normalized signed facet Grams;
* equality with the one same-header chamber in the frozen prior ledger;
* literal Coxeter-relation matches in two pinned CoxIter snapshots; and
* equality with any of Alexeev's 18 explicit *full* chamber matrices.

It also records Solomon's exact Coxeter Euler characteristic.  All corpus
comparisons are bounded.  Allcock's 2006 infinite-family theorem independently
rules out a first-in-dimension interpretation in H10, H11, or H13.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path
import tarfile
from typing import Any, Sequence

import networkx as nx
import sympy as sp

from scripts.certify_alexeev_18_stock_denominator import EXPECTED_ALL_KEYS
from scripts.certify_alexeev_20_2_1_15_9_1_lineage_bounds import (
    _rational_comparisons,
    _scan_coxiter,
)
from scripts.certify_alexeev_high_dim_reflection_examples import _extract_matrix
from scripts.certify_jacobian_kummer_h10_maximal_reflection_lineage import (
    _parse_coxiter,
    _permutation as _relation_permutation,
)
from scripts.certify_jacobian_kummer_h10_reflection_index import (
    _euler_from_connected_components,
    _relation_matrix,
)
from scripts.certify_k3_rank18_19_all_subgroup_completion import (
    _parse_coxiter_text,
)
from scripts.certify_shimada_disc16_d5_h10_lineage import (
    _connected_elliptic_masks,
)
from scripts.inventory_alexeev_20_2_1_15_9_1_wall_symmetries import (
    SOURCE as ALEXEEV_SOURCE,
    SOURCE_SHA256 as ALEXEEV_SOURCE_SHA256,
)


ROOT = Path(__file__).resolve().parents[1]
CONSTRUCTION = (
    ROOT / "data/enumerations/alexeev_remaining_high_fixed_section_construction.json"
)
VINBERG = ROOT / "data/enumerations/alexeev_remaining_high_fixed_section_vinberg.json"
STABLE_OVERLATTICES = (
    ROOT
    / "data/enumerations/alexeev_remaining_high_fixed_section_stable_overlattices.json"
)
D4_C007 = (
    ROOT
    / "data/enumerations/shimada_disc16_atlas_archived_high_sections/D0004_C007.json"
)
OUTPUT = (
    ROOT / "data/enumerations/alexeev_remaining_high_fixed_section_lineage.json"
)
COXITER_ROOT = Path("/private/tmp/coxiter-current/CoxIter-master")
DOUBLING_ARCHIVE = COXITER_ROOT / "graphs/testing_index2.tar.gz"
SCHEMA = "coxeter.alexeev_remaining_high_fixed_section_lineage.v1"


EXPECTED_EULER = {
    "A18_4_0_C03": "0",
    "A18_4_0_C04": "-527/1451520",
    "A18_4_0_C06": "0",
    "A18_4_0_C18": "-527/23224320",
    "A15_11_1_C02": "-527/4976640",
    "A15_7_1_C01": "-527/1658880",
}

# This is the positive ledger immediately before the present six-case batch.
# Header inequality is an exact obstruction to a square normalized-Gram
# permutation.  D16-D4-C7-H10-27 is the sole same-header record and is replayed
# entrywise below.
FROZEN_PRIOR_CHAMBERS = (
    ("F3-C2-H10-196", 10, 196),
    ("F3-V4-H10-123", 10, 123),
    ("A15-C2-H10-34", 10, 34),
    ("SS2-V4-H10-34", 10, 34),
    ("SS2-2D-H10-47", 10, 47),
    ("K3-S3-C2-H10-13", 10, 13),
    ("D16-D5-H10-30", 10, 30),
    ("D16-D4-C7-H10-27", 10, 27),
    ("JK-H10-0007", 10, 84),
    ("F3-C2-H11-206", 11, 206),
    ("A20-C2b-H11-20", 11, 20),
    ("F3-V4-H11-328", 11, 328),
    ("SS2-3B-H11-22", 11, 22),
    ("K3-D4-H12", 12, 17),
    ("JK-H12", 12, 172),
    ("A20-C2a-H12-19", 12, 19),
    ("SS2-2A-H13-50", 13, 50),
    ("SS2-2C-H14-98", 14, 98),
    ("D16-D5-H14", 14, 60),
    ("K3-D6-H14", 14, 17),
    ("X3-H19", 19, 36),
    ("BOR-D4-H21-210", 21, 210),
)

EXPECTED_DUPLICATE_PERMUTATION = [
    3, 1, 9, 2, 4, 7, 12, 17, 5, 8, 13, 18, 11, 16, 10, 21, 14, 6,
    15, 23, 19, 20, 22, 24, 25, 26, 27,
]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _matrix(payload: Sequence[Sequence[str | int]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(value) for value in row] for row in payload])


def _normalized_label(gram: sp.Matrix, left: int, right: int) -> str:
    if left == right:
        return "1"
    pairing = sp.Rational(gram[left, right])
    if pairing == 0:
        return "0"
    ratio = sp.factor(
        pairing**2
        / (sp.Rational(gram[left, left]) * sp.Rational(gram[right, right]))
    )
    sign = "-" if pairing < 0 else "+"
    return f"{sign}sqrt({ratio})"


def _normalized_graph(gram: sp.Matrix) -> nx.Graph:
    graph = nx.Graph()
    for vertex in range(gram.rows):
        signature = Counter(
            _normalized_label(gram, vertex, other)
            for other in range(gram.rows)
            if other != vertex
        )
        graph.add_node(vertex, signature=tuple(sorted(signature.items())))
    for left in range(gram.rows):
        for right in range(left + 1, gram.rows):
            graph.add_edge(
                left,
                right,
                normalized_label=_normalized_label(gram, left, right),
            )
    return graph


def _normalized_permutation(
    source: sp.Matrix, target: sp.Matrix
) -> list[int] | None:
    if source.shape != target.shape:
        return None
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        _normalized_graph(source),
        _normalized_graph(target),
        node_match=nx.algorithms.isomorphism.categorical_node_match(
            "signature", None
        ),
        edge_match=nx.algorithms.isomorphism.categorical_edge_match(
            "normalized_label", None
        ),
    )
    mapping = next(matcher.isomorphisms_iter(), None)
    if mapping is None:
        return None
    permutation = [int(mapping[index]) for index in range(source.rows)]
    if any(
        _normalized_label(source, left, right)
        != _normalized_label(target, permutation[left], permutation[right])
        for left in range(source.rows)
        for right in range(source.rows)
    ):
        raise ArithmeticError("normalized-Gram permutation failed replay")
    return permutation


def _solomon(gram_payload: Sequence[Sequence[str | int]]) -> dict[str, Any]:
    relations = _relation_matrix(gram_payload)
    connected, _ = _connected_elliptic_masks(gram_payload)
    result = _euler_from_connected_components(relations, connected)
    value = Fraction(result["rational_euler_characteristic"])
    return {
        "rational_euler_characteristic": str(value),
        "connected_elliptic_special_subdiagram_count": len(connected),
    }


def _alexeev_full_headers() -> dict[str, Any]:
    if _sha256(ALEXEEV_SOURCE) != ALEXEEV_SOURCE_SHA256:
        raise ArithmeticError("the pinned Alexeev ancillary source changed")
    source = ALEXEEV_SOURCE.read_text(encoding="utf-8")
    rows = []
    for key in sorted(EXPECTED_ALL_KEYS):
        gram = _extract_matrix(source, key)
        rows.append(
            {
                "key": list(key),
                "hyperbolic_dimension": key[0] - 1,
                "facet_count": gram.rows,
            }
        )
    return {
        "source_path": str(ALEXEEV_SOURCE),
        "source_sha256": ALEXEEV_SOURCE_SHA256,
        "explicit_full_matrix_count": len(rows),
        "full_chamber_headers": rows,
    }


def _coxiter_exact_candidates(
    cases: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    targets = {
        (
            int(case["case"]["hyperbolic_dimension"]),
            int(case["exact_reduction"]["active_hyperbolic_facet_count"]),
        )
        for case in cases.values()
    }
    audit = _scan_coxiter(targets)
    case_relations = {
        identifier: _relation_matrix(
            case["restricted_mirrors"]["positive_wall_gram"]
        )
        for identifier, case in cases.items()
    }
    parsed: dict[tuple[int, int], list[dict[str, Any]]] = {
        target: [] for target in targets
    }

    for target in sorted(targets):
        dimension, facets = target
        key = f"H{dimension}_facets_{facets}"
        for relative in audit["targets"][key]["ordinary_header_matches"]:
            path = COXITER_ROOT / relative
            parsed_facets, parsed_dimension, relations = _parse_coxiter(path)
            if (parsed_dimension, parsed_facets) != target:
                raise ArithmeticError("ordinary CoxIter header changed")
            parsed[target].append(
                {
                    "lane": "ordinary",
                    "path": relative,
                    "sha256": _sha256(path),
                    "case_relation_permutations_zero_based": {
                        identifier: _relation_permutation(
                            case_relations[identifier], relations
                        )
                        for identifier, case in cases.items()
                        if (
                            case["case"]["hyperbolic_dimension"],
                            case["exact_reduction"][
                                "active_hyperbolic_facet_count"
                            ],
                        )
                        == target
                    },
                }
            )

    with tarfile.open(DOUBLING_ARCHIVE, "r:gz") as archive:
        for target in sorted(targets):
            dimension, facets = target
            key = f"H{dimension}_facets_{facets}"
            for member_name in audit["targets"][key][
                "doubling_header_matches"
            ]:
                member = archive.getmember(member_name)
                handle = archive.extractfile(member)
                if handle is None:
                    raise ArithmeticError(f"cannot read {member_name}")
                raw = handle.read()
                parsed_facets, parsed_dimension, relations = _parse_coxiter_text(
                    raw.decode("utf-8", errors="replace")
                )
                if (parsed_dimension, parsed_facets) != target:
                    raise ArithmeticError("doubling CoxIter header changed")
                parsed[target].append(
                    {
                        "lane": "index_two_doubling",
                        "path": member_name,
                        "sha256": hashlib.sha256(raw).hexdigest(),
                        "case_relation_permutations_zero_based": {
                            identifier: _relation_permutation(
                                case_relations[identifier], relations
                            )
                            for identifier, case in cases.items()
                            if (
                                case["case"]["hyperbolic_dimension"],
                                case["exact_reduction"][
                                    "active_hyperbolic_facet_count"
                                ],
                            )
                            == target
                        },
                    }
                )

    return {
        **audit,
        "exact_relation_graph_candidate_records": {
            f"H{dimension}_facets_{facets}": parsed[(dimension, facets)]
            for dimension, facets in sorted(parsed)
        },
        "no_exact_relation_graph_match": all(
            permutation is None
            for records in parsed.values()
            for record in records
            for permutation in record[
                "case_relation_permutations_zero_based"
            ].values()
        ),
    }


def build() -> dict[str, Any]:
    construction = json.loads(CONSTRUCTION.read_text(encoding="utf-8"))
    vinberg = json.loads(VINBERG.read_text(encoding="utf-8"))
    stable = json.loads(STABLE_OVERLATTICES.read_text(encoding="utf-8"))
    if not (
        construction["all_checks_pass"]
        and vinberg["all_checks_pass"]
        and stable["all_checks_pass"]
    ):
        raise ArithmeticError("an upstream construction/volume check failed")
    cases = {case["id"]: case for case in construction["cases"]}
    stable_cases = {case["id"]: case for case in stable["cases"]}
    if set(cases) != set(EXPECTED_EULER):
        raise ArithmeticError("the six-case lineage set changed")
    if set(stable_cases) != set(cases):
        raise ArithmeticError("the stable-overlattice case set changed")

    full = _alexeev_full_headers()
    prior_by_header: dict[tuple[int, int], list[str]] = {}
    for identifier, dimension, facets in FROZEN_PRIOR_CHAMBERS:
        prior_by_header.setdefault((dimension, facets), []).append(identifier)

    records = []
    for identifier, case in cases.items():
        dimension = int(case["case"]["hyperbolic_dimension"])
        facets = int(case["exact_reduction"]["active_hyperbolic_facet_count"])
        gram = _matrix(case["restricted_mirrors"]["positive_wall_gram"])
        lattice = _matrix(case["saturated_fixed_lattice"]["gram"])
        solomon = _solomon(
            case["restricted_mirrors"]["positive_wall_gram"]
        )
        if solomon["rational_euler_characteristic"] != EXPECTED_EULER[identifier]:
            raise ArithmeticError(f"Solomon value changed for {identifier}")
        stable_case = stable_cases[identifier]
        maximal = stable_case["unique_maximal_reflection_stable_overlattice"]
        maximal_lineage = stable_case["maximal_reflection_lineage"]

        full_header_matches = [
            row
            for row in full["full_chamber_headers"]
            if (row["hyperbolic_dimension"], row["facet_count"])
            == (dimension, facets)
        ]
        records.append(
            {
                "id": identifier,
                "hyperbolic_dimension": dimension,
                "facet_count": facets,
                "positive_wall_gram_sha256": hashlib.sha256(
                    json.dumps(
                        case["restricted_mirrors"]["positive_wall_gram"],
                        sort_keys=True,
                        separators=(",", ":"),
                    ).encode()
                ).hexdigest(),
                "fixed_lattice_determinant": int(lattice.det()),
                "fixed_lattice_smith_diagonal_absolute": case[
                    "saturated_fixed_lattice"
                ]["smith_diagonal_absolute"],
                "rational_similarity_comparisons": _rational_comparisons(lattice),
                "solomon": solomon,
                "reflection_stable_overlattice": {
                    "stable_integral_subgroup_orders": stable_case[
                        "exhaustive_discriminant_module_audit"
                    ]["stable_integral_subgroup_orders"],
                    "stable_even_subgroup_orders": stable_case[
                        "exhaustive_discriminant_module_audit"
                    ]["stable_even_subgroup_orders"],
                    "unique_maximal": maximal["unique"],
                    "index_over_fixed_lattice": maximal[
                        "index_over_original_lattice"
                    ],
                    "determinant": maximal["determinant"],
                    "parity": maximal["parity"],
                    "smith_invariants": maximal["smith_invariants"],
                    "scope": (
                        "maximal among integral overlattices preserved by every "
                        "displayed facet reflection"
                    ),
                },
                "containing_full_reflection_chamber": {
                    "scharlau_kirschmer_list": maximal_lineage["list"],
                    "entry_one_based": maximal_lineage["entry_one_based"],
                    "genus_symbol": maximal_lineage["genus_symbol"],
                    "facet_count": maximal_lineage["full_chamber_facets"],
                    "finite_vertex_count": maximal_lineage[
                        "full_chamber_finite_vertices"
                    ],
                    "cusp_count": maximal_lineage["full_chamber_cusps"],
                    "coxeter_euler_characteristic": maximal_lineage[
                        "full_coxeter_euler_characteristic"
                    ],
                    "section_to_full_index": maximal_lineage[
                        "section_to_full_index"
                    ],
                    "index_status": maximal_lineage["index_status"],
                },
                "full_Alexeev_same_header_records": full_header_matches,
                "no_full_Alexeev_normalized_gram_match_by_exact_header_obstruction": (
                    not full_header_matches
                ),
                "frozen_prior_ledger_same_header_ids": prior_by_header.get(
                    (dimension, facets), []
                ),
            }
        )

    mutual = []
    identifiers = list(cases)
    for offset, left_id in enumerate(identifiers):
        left = cases[left_id]
        left_header = (
            left["case"]["hyperbolic_dimension"],
            left["exact_reduction"]["active_hyperbolic_facet_count"],
        )
        for right_id in identifiers[offset + 1 :]:
            right = cases[right_id]
            right_header = (
                right["case"]["hyperbolic_dimension"],
                right["exact_reduction"]["active_hyperbolic_facet_count"],
            )
            if left_header != right_header:
                continue
            permutation = _normalized_permutation(
                _matrix(left["restricted_mirrors"]["positive_wall_gram"]),
                _matrix(right["restricted_mirrors"]["positive_wall_gram"]),
            )
            mutual.append(
                {
                    "source": left_id,
                    "target": right_id,
                    "source_to_target_facet_permutation_zero_based": permutation,
                    "same_signed_normalized_gram": permutation is not None,
                }
            )

    d4 = json.loads(D4_C007.read_text(encoding="utf-8"))
    duplicate_id = "A15_7_1_C01"
    duplicate_permutation = _normalized_permutation(
        _matrix(cases[duplicate_id]["restricted_mirrors"]["positive_wall_gram"]),
        _matrix(d4["coxeter_gate"]["positive_wall_gram"]),
    )
    duplicate_one_based = (
        [value + 1 for value in duplicate_permutation]
        if duplicate_permutation is not None
        else None
    )
    if duplicate_one_based != EXPECTED_DUPLICATE_PERMUTATION:
        raise ArithmeticError("the exact D0004/C007 duplicate permutation changed")

    coxiter = _coxiter_exact_candidates(cases)
    checks = {
        "six_cases": len(records) == 6,
        "all_six_vinberg_finite_volume_noncompact": vinberg["summary"] == {
            "bad_edge_count": 0,
            "bad_parabolic_neighborhood_count": 0,
            "case_count": 6,
            "finite_volume_noncompact_count": 6,
        },
        "all_six_in_odd_unimodular_rational_class": all(
            row["rational_similarity_comparisons"]["f1_odd_unimodular"]
            ["rationally_similar"]
            for row in records
        ),
        "all_six_Solomon_values_replay": all(
            row["solomon"]["rational_euler_characteristic"]
            == EXPECTED_EULER[row["id"]]
            for row in records
        ),
        "two_H10_21_chambers_are_nonisomorphic": mutual == [
            {
                "source": "A18_4_0_C18",
                "target": "A15_11_1_C02",
                "source_to_target_facet_permutation_zero_based": None,
                "same_signed_normalized_gram": False,
            }
        ],
        "exact_D0004_C007_duplicate_replays": duplicate_one_based
        == EXPECTED_DUPLICATE_PERMUTATION,
        "other_five_absent_from_frozen_prior_ledger_by_header": all(
            row["frozen_prior_ledger_same_header_ids"] == []
            for row in records
            if row["id"] != duplicate_id
        ),
        "all_six_absent_from_full_Alexeev_matrices_by_header": all(
            row[
                "no_full_Alexeev_normalized_gram_match_by_exact_header_obstruction"
            ]
            for row in records
        ),
        "pinned_CoxIter_has_no_exact_relation_graph_match": coxiter[
            "no_exact_relation_graph_match"
        ],
        "permanent_stable_overlattice_certificate_passes": (
            stable["schema"]
            == "coxeter.alexeev_remaining_high_fixed_sections.stable_overlattices.v1"
            and stable["source"]["construction_sha256"] == _sha256(CONSTRUCTION)
            and stable["summary"]["case_count"] == 6
            and stable["summary"]["unique_maximal_stable_overlattice_count"]
            == 6
        ),
        "four_exact_even_dimensional_indices": stable["summary"][
            "exact_numerical_indices"
        ]
        == {
            "A18_4_0_C04": 130560,
            "A18_4_0_C18": 144,
            "A15_11_1_C02": 6720,
            "A15_7_1_C01": 2016,
        },
        "odd_dimensional_index_values_explicitly_open": stable["summary"][
            "zero_characteristic_index_values_open"
        ]
        == ["A18_4_0_C03", "A18_4_0_C06"],
        "no_new_arithmetic_commensurability_class": stable["summary"][
            "new_arithmetic_commensurability_class_count"
        ]
        == 0,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"remaining-section lineage checks failed: {checks}")

    return {
        "schema": SCHEMA,
        "sources": {
            path.name: {
                "path": str(path.relative_to(ROOT)),
                "sha256": _sha256(path),
            }
            for path in (CONSTRUCTION, VINBERG, STABLE_OVERLATTICES, D4_C007)
        },
        "scope": {
            "proved": (
                "exact rational-class, Solomon, normalized-Gram duplicate, "
                "bounded Alexeev/prior-ledger/CoxIter comparisons, and an "
                "exhaustive reflection-stable integral-overlattice census"
            ),
            "not_proved": (
                "worldwide literal-diagram priority, a new arithmetic "
                "commensurability class, geometric lifting of the selected "
                "diagram subgroups, or numerical indices in H11 and H13"
            ),
            "historical_ceiling": {
                "source": (
                    "Daniel Allcock, Infinitely many hyperbolic Coxeter groups "
                    "through dimension 19, Geometry & Topology 10 (2006), 737-758"
                ),
                "doi": "https://doi.org/10.2140/gt.2006.10.737",
                "consequence": (
                    "H10, H11, and H13 already contain infinitely many "
                    "finite-volume Coxeter polytopes; no first-in-dimension "
                    "claim is possible"
                ),
            },
        },
        "cases": records,
        "within_batch_same_header_comparisons": mutual,
        "prior_ledger_exact_duplicate": {
            "source": duplicate_id,
            "target": "D16-D4-C7-H10-27 = D0004/C007",
            "source_to_target_facet_permutation_one_based": duplicate_one_based,
            "same_signed_normalized_gram": True,
            "target_certificate_path": str(D4_C007.relative_to(ROOT)),
            "target_certificate_sha256": _sha256(D4_C007),
        },
        "full_Alexeev_matrix_comparison": full,
        "frozen_prior_positive_ledger_headers": [
            {"id": identifier, "dimension": dimension, "facets": facets}
            for identifier, dimension, facets in FROZEN_PRIOR_CHAMBERS
        ],
        "coxiter_bounded_comparison": coxiter,
        "reflection_stable_overlattice_status": {
            "permanent_exhaustive_certificate": str(
                STABLE_OVERLATTICES.relative_to(ROOT)
            ),
            "sha256": _sha256(STABLE_OVERLATTICES),
            "schema": stable["schema"],
            "scope": stable["method"]["scope_warning"],
            "summary": stable["summary"],
            "status": "closed for all six displayed reflection groups",
        },
        "checks": checks,
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.is_file() or args.output.read_text(encoding="utf-8") != rendered:
            raise ArithmeticError(f"archived lineage artifact differs: {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(json.dumps(payload["checks"], sort_keys=True))
    print(args.output)


if __name__ == "__main__":
    main()
