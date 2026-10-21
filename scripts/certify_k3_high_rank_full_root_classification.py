#!/usr/bin/env python3
"""Certify the high-rank full-root K3 finite-volume classification transfer.

This is a synthesis certificate.  It combines the pinned Nikulin--Vinberg
classification table as reproduced by Roulleau, the rankwise count in
Artebani--Hausen--Laface, the project's exact stock inventory, the nineteen
literal Vinberg endpoint/full-neighbourhood replays, and the oriented
arithmetic period/change-of-marking classification.  The ambient
K3Root[11,20] class does not assume finite reflection index.  The certificate
identifies its finite-volume sublocus but does not claim strict-M3c closure,
explicit presentations of the infinite arithmetic groups, or independently
certified geometric family monodromy.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
INVENTORY_PATH = (
    ROOT / "data/enumerations/roulleau_rank11_14_source_stock_inventory.json"
)
HIGH_GRAM_PATH = ROOT / "data/known_examples/k3_high_rank_nef_chambers.json"
HIGH_VINBERG_PATH = (
    ROOT / "data/known_examples/k3_high_rank_nef_vinberg_volume.json"
)
LOW_VINBERG_PATH = (
    ROOT / "data/enumerations/k3_rank11_14_literal_vinberg_volume.json"
)
ORIENTED_PERIOD_PATH = (
    ROOT / "data/enumerations/k3_high_rank_oriented_monodromy.json"
)
OUTPUT_PATH = (
    ROOT / "data/enumerations/k3_high_rank_full_root_classification.json"
)

ROULLEAU_TEX = Path("/private/tmp/roulleau_source/epiga_vol6_Roulleau.tex")
AHL_PDF = Path("/private/tmp/ahl_k3.pdf")
AHL_TEXT = Path("/private/tmp/ahl_k3.txt")

PINNED_SOURCE_HASHES = {
    "roulleau_tex": "e6eeaaa4a3f375f209762718e8ef9a361b0a95356cc8be2e9f99a6713abf5df3",
    "ahl_pdf": "44553cad30080e7aa4feb32b13147d5f689261bbee0adbd95d53869f483f6f47",
    "ahl_text": "3278ddaf1215222b60206d012575c6b6fb08486128a344168453314af439a91c",
}

# AHL Theorem 2.12 gives these counts for the polyhedral-effective-cone,
# equivalently finite-volume complete-root, locus rank by rank: four for each
# of ranks 11--12, three for each of ranks 13--14, one for each rank 15--19,
# and zero for rank 20.  They do not count every lattice in ambient K3Root.
EXPECTED_RANK_COUNTS = {
    11: 4,
    12: 4,
    13: 3,
    14: 3,
    15: 1,
    16: 1,
    17: 1,
    18: 1,
    19: 1,
    20: 0,
}

# key, rho, lattice, |disc|, length of the discriminant group, facets
EXPECTED_CASES = [
    ("R11_UE8A1", 11, "U+E8(-1)+A1(-1)", 2, 1, 12),
    ("R11_UD8A1", 11, "U+D8(-1)+A1(-1)", 8, 3, 14),
    ("R11_UD4D4A1", 11, "U+D4(-1)^2+A1(-1)", 32, 5, 22),
    ("R11_UD4A1x5", 11, "U+D4(-1)+A1(-1)^5", 128, 7, 90),
    ("R12_UE8A1x2", 12, "U+E8(-1)+A1(-1)^2", 4, 2, 14),
    ("R12_UD8A1x2", 12, "U+D8(-1)+A1(-1)^2", 16, 4, 19),
    ("R12_UD4x2A1x2", 12, "U+D4(-1)^2+A1(-1)^2", 64, 6, 59),
    ("R12_UA2E8", 12, "U+A2(-1)+E8(-1)", 3, 1, 13),
    ("R13_UE8A1x3", 13, "U+E8(-1)+A1(-1)^3", 8, 3, 17),
    ("R13_UD8A1x3", 13, "U+D8(-1)+A1(-1)^3", 32, 5, 39),
    ("R13_UE8A3", 13, "U+E8(-1)+A3(-1)", 4, 1, 14),
    ("R14_UE8D4", 14, "U+E8(-1)+D4(-1)", 4, 2, 15),
    ("R14_UD8D4", 14, "U+D8(-1)+D4(-1)", 16, 4, 20),
    ("R14_UE8A1x4", 14, "U+E8(-1)+A1(-1)^4", 16, 4, 27),
    ("R15_UE8D4A1", 15, "U+E8(-1)+D4(-1)+A1(-1)", 8, 3, 21),
    ("R16_UE8D6", 16, "U+E8(-1)+D6(-1)", 4, 2, 19),
    ("R17_UE8E7", 17, "U+E8(-1)+E7(-1)", 2, 1, 19),
    ("R18_UE8E8", 18, "U+E8(-1)^2", 1, 0, 19),
    ("R19_UE8E8A1", 19, "U+E8(-1)^2+A1(-1)", 2, 1, 24),
]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _compact_lattice_name(name: str) -> str:
    return name.replace(" ", "")


def check_pinned_sources(
    roulleau_tex: Path = ROULLEAU_TEX,
    ahl_pdf: Path = AHL_PDF,
    ahl_text: Path = AHL_TEXT,
) -> dict[str, Any]:
    """Check the exact local source copies and line-local witnesses."""

    paths = {
        "roulleau_tex": roulleau_tex,
        "ahl_pdf": ahl_pdf,
        "ahl_text": ahl_text,
    }
    for key, path in paths.items():
        if not path.exists():
            raise FileNotFoundError(path)
        if _sha256(path) != PINNED_SOURCE_HASHES[key]:
            raise AssertionError(f"source hash mismatch: {key}")

    inventory = _json(INVENTORY_PATH)
    tex_lines = roulleau_tex.read_text().splitlines()
    introduction = "\n".join(tex_lines[196:218])
    assert "list of $118$ N\\'eron--Severi" in introduction
    moduli_excerpt = "\n".join(tex_lines[1124:1229])
    assert "The $118$ moduli spaces" in moduli_excerpt
    assert "There are $28$ lattices $L$ of rank at least $10$" in moduli_excerpt

    for stock in inventory["stocks"]:
        start, end = stock["source"]["table_lines"]
        excerpt = "\n".join(tex_lines[start - 1 : end])
        assert stock["tex_table_lattice_token"] in excerpt, stock["key"]
        assert f"${stock['wall_count']}$" in excerpt, stock["key"]

    ahl = ahl_text.read_text(errors="replace")
    assert "Theorem 2.12" in ahl
    assert "11−12" in ahl and "15−19" in ahl

    return {
        "checked": True,
        "all_checks_pass": True,
        "sha256": {key: _sha256(path) for key, path in paths.items()},
        "roulleau_witnesses": {
            "classification_118_lines": [203, 212],
            "irreducible_moduli_proposition_lines": [1129, 1227],
            "rank_at_least_10_count_lines": [1209, 1227],
            "rank_at_least_11_table_lines": [10906, 10966],
        },
        "ahl_witness": {
            "result": "Theorem 2.12",
            "finite_volume_rank_counts": {
                str(rank): count for rank, count in EXPECTED_RANK_COUNTS.items()
            },
        },
    }


def _validate_inputs() -> tuple[
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
]:
    inventory = _json(INVENTORY_PATH)
    high_gram = _json(HIGH_GRAM_PATH)
    high_vinberg = _json(HIGH_VINBERG_PATH)
    low_vinberg = _json(LOW_VINBERG_PATH)
    oriented_period = _json(ORIENTED_PERIOD_PATH)

    assert inventory["all_checks_pass"] is True
    assert inventory["source_check"] == {"checked": True, "all_checks_pass": True}
    assert inventory["denominator"]["rank_at_least_11_stock_count"] == 19
    assert inventory["denominator"]["already_exactly_machine_represented_count"] == 19
    assert len(inventory["stocks"]) == 19

    expected = {item[0]: item for item in EXPECTED_CASES}
    actual_keys = {stock["key"] for stock in inventory["stocks"]}
    assert actual_keys == set(expected)
    for stock in inventory["stocks"]:
        key, rho, lattice, disc, _length, facets = expected[stock["key"]]
        assert stock["key"] == key
        assert stock["picard_rank"] == rho
        assert stock["hyperbolic_dimension"] == rho - 1
        assert _compact_lattice_name(stock["lattice"]) == lattice
        assert stock["absolute_discriminant"] == disc
        assert stock["wall_count"] == facets
        assert stock["exact_wall_gram_extractable"] is True

    rank_counts = {
        rho: sum(stock["picard_rank"] == rho for stock in inventory["stocks"])
        for rho in EXPECTED_RANK_COUNTS
    }
    assert rank_counts == EXPECTED_RANK_COUNTS

    assert high_gram["artifact_status"] == "source_audited_gram_verified"
    assert {int(case["rho"]) for case in high_gram["cases"]} == set(range(15, 20))
    for case in high_gram["cases"]:
        rho = int(case["rho"])
        stock = next(s for s in inventory["stocks"] if s["picard_rank"] == rho)
        assert case["wall_count"] == stock["wall_count"]
        assert case["absolute_discriminant"] == stock["absolute_discriminant"]
        assert case["expected_rank"] == rho
        assert case["expected_signature"] == [1, rho - 1, case["wall_count"] - rho]

    literal_cases = high_vinberg["cases"]
    assert {int(case["rho"]) for case in literal_cases} == set(range(15, 20))
    assert all(case["all_checks_pass"] for case in literal_cases)
    assert all(case["decision"] == "finite_volume_noncompact" for case in literal_cases)
    assert all(case["finite_volume_bad_edge_count"] == 0 for case in literal_cases)
    assert all(
        case["bad_full_parabolic_neighborhood_count"] == 0
        for case in literal_cases
    )

    low_literal_cases = low_vinberg["cases"]
    assert low_vinberg["all_checks_pass"] is True
    assert low_vinberg["checks"]["all_literal_vinberg_audits_pass"] is True
    assert low_vinberg["checks"]["all_finite_volume_noncompact"] is True
    assert low_vinberg["totals"] == {
        "bad_edge_count": 0,
        "bad_parabolic_neighborhood_count": 0,
        "case_count": 14,
        "connected_affine_component_count": 4177,
        "elliptic_edge_count": 106630,
        "facet_count": 375,
        "finite_vertex_count": 12854,
        "ideal_vertex_count": 1119,
    }
    transcription = low_vinberg["rank11_22_source_transcription_diagnostic"]
    assert transcription["tex_lines"] == [8206, 8318]
    assert transcription["expected_lattice_absolute_determinant"] == 32
    assert transcription["literal_basis_graph_absolute_determinant"] == 128
    assert transcription["printed_dependent_curve_squares"]["A13"] == -14
    assert transcription["obstruction_is_terminal"] is False
    assert low_vinberg["rank11_22_reconstruction"]["all_checks_pass"] is True
    assert {case["key"] for case in low_literal_cases} == {
        key for key, rho, *_rest in EXPECTED_CASES if rho <= 14
    }
    for case in low_literal_cases:
        key, rho, lattice, _disc, _length, facets = expected[case["key"]]
        assert case["key"] == key
        assert case["rho"] == rho
        assert case["hyperbolic_dimension"] == rho - 1
        assert case["lattice"] == lattice
        assert case["facet_count"] == facets
        assert case["all_checks_pass"] is True
        assert case["simple_roots_complete"] is True
        assert case["decision"] == "finite_volume_noncompact"
        assert case["bad_edge_count"] == 0
        assert case["bad_parabolic_neighborhood_count"] == 0
        assert all(case["gram_checks"].values())
        if key == "R11_UD4D4A1":
            assert case["matrix_provenance"]["kind"] == (
                "independent_standard_lattice_norm_two_vinberg_reconstruction"
            )

    assert oriented_period["all_checks_pass"] is True
    assert oriented_period["scope"]["case_count"] == 19
    assert oriented_period["literal_vinberg_status"]["all_19_closed"] is True
    assert (
        oriented_period["literal_vinberg_status"]["remaining_literal_V_gap_count"]
        == 0
    )
    assert oriented_period["strict_M3c_residual"]["status"] == "open"
    assert {row["key"] for row in oriented_period["rows"]} == set(expected)
    for row in oriented_period["rows"]:
        key, rho, lattice, disc, length, facets = expected[row["key"]]
        assert row["key"] == key
        assert row["picard_rank"] == rho
        assert row["lattice"] == lattice
        assert row["absolute_discriminant"] == disc
        assert row["discriminant_length"] == length
        assert row["facet_count"] == facets
        assert row["primitive_embedding_orbits"] == {
            "O_plus_no_split_reason": (
                "compose a wrong-orientation transporter with the pointwise "
                "complement-U positive reflection"
            ),
            "under_O_plus_Lambda_K3": 1,
            "under_full_O_Lambda_K3": 1,
        }
        arithmetic = row["oriented_period_change_of_marking_data"]
        assert arithmetic["pointwise_N_polarized_T_period_image"].startswith(
            "tilde O^+(T)"
        )
        assert arithmetic["unmarked_ambient_chamber_normalizer"] == "Gamma^+_{j,C}"
        assert arithmetic["unmarked_chamber_normalizer_T_period_image"] == "O^+(T)"
        assert arithmetic["geometric_monodromy_identification"].startswith(
            "conditional_only"
        )
    return inventory, high_gram, high_vinberg, low_vinberg, oriented_period


def build_artifact(source_check: dict[str, Any] | None = None) -> dict[str, Any]:
    (
        inventory,
        high_gram,
        high_vinberg,
        low_vinberg,
        oriented_period,
    ) = _validate_inputs()
    expected = {item[0]: item for item in EXPECTED_CASES}
    low_literal_by_key = {case["key"]: case for case in low_vinberg["cases"]}
    high_key_by_rank = {rho: key for key, rho, *_rest in EXPECTED_CASES if rho >= 15}
    high_literal_by_key = {
        high_key_by_rank[int(case["rho"])]: case for case in high_vinberg["cases"]
    }
    literal_by_key = low_literal_by_key | high_literal_by_key
    high_gram_by_rank = {int(case["rho"]): case for case in high_gram["cases"]}
    oriented_by_key = {row["key"]: row for row in oriented_period["rows"]}

    cases: list[dict[str, Any]] = []
    for stock in inventory["stocks"]:
        key, rho, lattice, disc, discriminant_length, facets = expected[stock["key"]]
        complement_rank = 22 - rho
        stable_unique = discriminant_length + 2 <= complement_rank
        assert stable_unique
        literal = literal_by_key[key]
        high_case = high_gram_by_rank.get(rho)
        oriented = oriented_by_key[key]
        if rho <= 14:
            literal_summary = {
                "available": True,
                "decision": literal["decision"],
                "elliptic_edge_count": literal["elliptic_edge_count"],
                "finite_vertex_count": literal["finite_vertex_count"],
                "ideal_vertex_count": literal["ideal_vertex_count"],
                "connected_affine_audit_count": literal[
                    "connected_affine_component_count"
                ],
                "finite_volume_bad_edge_count": literal["bad_edge_count"],
                "bad_full_parabolic_neighborhood_count": literal[
                    "bad_parabolic_neighborhood_count"
                ],
                "literal_audit_sha256": literal[
                    "literal_verifier_stdout_sha256"
                ],
                "source_package": str(LOW_VINBERG_PATH.relative_to(ROOT)),
            }
            gram_summary = {
                "available": True,
                "expected_rank": rho,
                "expected_signature": [1, rho - 1, facets - rho],
                "all_gram_checks_pass": all(literal["gram_checks"].values()),
                "matrix_provenance": literal["matrix_provenance"],
                "source_package": str(LOW_VINBERG_PATH.relative_to(ROOT)),
            }
        else:
            literal_summary = {
                "available": True,
                "decision": literal["decision"],
                "elliptic_edge_count": literal["elliptic_edge_count"],
                "finite_vertex_count": literal["finite_vertex_count"],
                "ideal_vertex_count": literal["ideal_vertex_count"],
                "connected_affine_audit_count": literal[
                    "connected_affine_audit_count"
                ],
                "finite_volume_bad_edge_count": literal[
                    "finite_volume_bad_edge_count"
                ],
                "bad_full_parabolic_neighborhood_count": literal[
                    "bad_full_parabolic_neighborhood_count"
                ],
                "literal_audit_sha256": literal["canonical_full_audit_sha256"],
                "source_package": str(HIGH_VINBERG_PATH.relative_to(ROOT)),
            }
            assert high_case is not None
            gram_summary = {
                "available": True,
                "expected_rank": high_case["expected_rank"],
                "expected_signature": high_case["expected_signature"],
                "source_package": str(HIGH_GRAM_PATH.relative_to(ROOT)),
            }
        cases.append(
            {
                "key": key,
                "picard_rank": rho,
                "hyperbolic_dimension": rho - 1,
                "lattice": lattice,
                "absolute_discriminant": disc,
                "discriminant_group_length": discriminant_length,
                "k3_complement_rank": complement_rank,
                "nikulin_stable_unique_embedding_inequality": {
                    "left": discriminant_length + 2,
                    "right": complement_rank,
                    "passes": stable_unique,
                },
                "primitive_embedding_orbit": (
                    "one orbit under both full O(Lambda_K3) and "
                    "O^+(Lambda_K3)"
                ),
                "facet_count": facets,
                "pairwise_distinction_key": [rho - 1, facets],
                "complete_root_chamber": True,
                "complete_root_evidence": stock["source_completeness_method"],
                "source_section_lines": stock["source"]["section_lines"],
                "source_table_lines": stock["source"]["table_lines"],
                "general_surface_automorphism_group": stock[
                    "general_surface_automorphism_group"
                ],
                "exact_picard_moduli_dimension": 20 - rho,
                "coarse_moduli_status": (
                    "Roulleau Proposition The-118-moduli-irred proves the "
                    "L-polarized coarse moduli irreducible; its exact-Picard "
                    "Noether-Lefschetz complement is nonempty and has the stated "
                    "dimension"
                ),
                "contains_U": True,
                "noncompactness_proof": (
                    "U supplies a rational isotropic class; Weyl transport gives a "
                    "primitive nef isotropic ray and hence an ideal cusp"
                ),
                "finite_volume_proof": (
                    "classification_output: the Nikulin--Vinberg finite-"
                    "automorphism classification puts this row in the finite-index "
                    "W_2(N) locus; arithmetic covolume then gives finite volume"
                ),
                "literal_vinberg_status": literal_summary,
                "exact_gram_replay_in_this_synthesis": gram_summary,
                "oriented_arithmetic_period_data": {
                    "primitive_embedding_orbits_under_full_O": 1,
                    "primitive_embedding_orbits_under_O_plus": 1,
                    "positive_cone_component_orbits_under_O_plus": oriented[
                        "positive_cone_components"
                    ]["orbits_under_O_plus_N"],
                    "complete_chamber_orbits_in_one_component_under_O_plus": oriented[
                        "complete_root_chambers_in_one_positive_component"
                    ]["orbits_under_O_plus_N"],
                    "discriminant_group_invariants": oriented[
                        "discriminant_group_invariants"
                    ],
                    "discriminant_form_orthogonal_group_order": oriented[
                        "discriminant_form_orthogonal_group_order"
                    ],
                    "pointwise_N_polarized_T_period_image": oriented[
                        "oriented_period_change_of_marking_data"
                    ]["pointwise_N_polarized_T_period_image"],
                    "unmarked_ambient_chamber_normalizer": "Gamma^+_{j,C}",
                    "unmarked_chamber_normalizer_T_period_image": "O^+(T)",
                    "geometric_monodromy_identification": oriented[
                        "oriented_period_change_of_marking_data"
                    ]["geometric_monodromy_identification"],
                },
                "novelty": "published_known_full_nef_Coxeter_polytope",
            }
        )

    cases.sort(key=lambda case: (case["picard_rank"], case["facet_count"]))
    distinction_keys = [tuple(case["pairwise_distinction_key"]) for case in cases]
    assert len(set(distinction_keys)) == 19
    assert sum(case["facet_count"] for case in cases) == 477
    assert [case["exact_picard_moduli_dimension"] for case in cases] == [
        9, 9, 9, 9,
        8, 8, 8, 8,
        7, 7, 7,
        6, 6, 6,
        5, 4, 3, 2, 1,
    ]
    vinberg_totals = {
        "facet_count": sum(case["facet_count"] for case in cases),
        "elliptic_edge_count": sum(
            case["literal_vinberg_status"]["elliptic_edge_count"]
            for case in cases
        ),
        "finite_vertex_count": sum(
            case["literal_vinberg_status"]["finite_vertex_count"]
            for case in cases
        ),
        "ideal_vertex_count": sum(
            case["literal_vinberg_status"]["ideal_vertex_count"]
            for case in cases
        ),
        "connected_affine_audit_count": sum(
            case["literal_vinberg_status"]["connected_affine_audit_count"]
            for case in cases
        ),
        "finite_volume_bad_edge_count": sum(
            case["literal_vinberg_status"]["finite_volume_bad_edge_count"]
            for case in cases
        ),
        "bad_full_parabolic_neighborhood_count": sum(
            case["literal_vinberg_status"][
                "bad_full_parabolic_neighborhood_count"
            ]
            for case in cases
        ),
    }
    assert vinberg_totals == {
        "facet_count": 477,
        "elliptic_edge_count": 117845,
        "finite_vertex_count": 14098,
        "ideal_vertex_count": 1156,
        "connected_affine_audit_count": 4255,
        "finite_volume_bad_edge_count": 0,
        "bad_full_parabolic_neighborhood_count": 0,
    }

    return {
        "schema": "coxeter.k3_high_rank_full_root_classification.v3",
        "artifact_id": "k3_high_rank_full_root_classification",
        "all_checks_pass": True,
        "surface_free_class": {
            "name": "K3Root[11,20]",
            "definition": [
                "N is an even lattice of signature (1,rho-1), 11 <= rho <= 20",
                "N embeds primitively in Lambda_K3",
                "C is a component of the complement of every square-minus-two root hyperplane",
            ],
            "surface_language_used_in_definition": False,
            "finite_volume_or_finite_index_used_in_definition": False,
        },
        "finite_volume_sublocus": {
            "name": "K3RootFV[11,20]",
            "definition": "members of K3Root[11,20] whose complete root chamber has finite hyperbolic volume",
            "finite_index_equivalence": (
                "finite volume is equivalent to [O^+(N):W_2(N)] < infinity; "
                "this is a classification output, not an ambient hypothesis"
            ),
        },
        "classification": {
            "classified_locus": "finite_volume_sublocus_only",
            "ambient_class_exhausted": False,
            "finite_volume_rank_counts": {
                str(rank): count for rank, count in EXPECTED_RANK_COUNTS.items()
            },
            "rank_11_through_19_finite_volume_type_count": 19,
            "rank_20_finite_volume_type_count": 0,
            "finite_volume_hyperbolic_dimension_range_realized": [10, 18],
            "maximum_finite_volume_hyperbolic_dimension": 18,
            "total_facets_across_types": 477,
            "one_unmarked_chamber_type_per_lattice_reason": (
                "W_2(N) acts transitively on complete root chambers"
            ),
            "pairwise_abstract_polytope_distinctness": (
                "the pair (hyperbolic dimension, facet count) is different for "
                "every one of the 19 cases"
            ),
            "worldwide_novelty": False,
            "novelty_scope": (
                "new project synthesis and unmarked classification/finiteness "
                "packaging of the published Nikulin--Vinberg--Roulleau "
                "classification; no new abstract polytope or strict-M3c closure "
                "is claimed"
            ),
            "surface_side_equivalences": [
                "finite automorphism group",
                "rational-polyhedral effective cone",
                "finitely generated Cox ring",
            ],
            "literal_vinberg_audit_totals": vinberg_totals,
        },
        "transfer_theorem": {
            "statement": (
                "The finite-volume unmarked integral locus of K3Root[11,20], "
                "denoted K3RootFV[11,20], consists exactly of the 19 listed "
                "cases; each is the complete nef chamber of a nonempty exact-"
                "Picard locus in an irreducible lattice-polarized K3 moduli "
                "family, and there is no finite-volume rank-20 case."
            ),
            "U": {
                "status": "closed",
                "proof": (
                    "primitive K3 embedding plus period surjectivity produces an "
                    "exact-Picard K3 surface for every ambient member; Nikulin--"
                    "Vinberg exhaust the finite-volume/finite-index locus"
                ),
            },
            "R": {
                "status": "closed_for_integrally_marked_N_and_complete_chamber",
                "proof": (
                    "Weyl transitivity changes the NS marking so that extraction "
                    "recovers the prescribed labelled lattice and chamber"
                ),
                "warning": (
                    "forgetting the integral marking is many-to-one on period points"
                ),
            },
            "C": {
                "status": "closed",
                "proof": (
                    "the K3 Weyl-chamber theorem identifies the nef cone with the "
                    "complete component cut by all square-minus-two roots; Roulleau's "
                    "lists are complete curve configurations"
                ),
            },
            "V": {
                "status": "closed_literally_and_structurally_for_all_19",
                "proof": (
                    "the surface classification puts each listed row in the finite-"
                    "automorphism/polyhedral-effective-cone locus, hence W_2(N) "
                    "has finite index in arithmetic O^+(N) and its fundamental "
                    "chamber has finite hyperbolic covolume; independently, exact "
                    "Vinberg/CoxIter endpoint and full affine-neighbourhood audits "
                    "pass for all nineteen complete wall systems; finite index was "
                    "not assumed in the ambient class"
                ),
                "literal_vinberg_replay_count": 19,
                "structural_only_at_synthesis_level_count": 0,
                "noncompact_all_19": True,
                "literal_audit_totals": vinberg_totals,
            },
            "strict_m3c": {
                "status": "open",
                "reasons": [
                    "the ambient class is not the bounded raw IC0 Coxeter specification required by strict M3c",
                    "no explicit generators-and-relations presentations are supplied for the 19 infinite arithmetic period groups",
                    "geometric family monodromy is not independently certified beyond the conditional standard moduli-stack/deck-group identification",
                ],
                "completed_scope": (
                    "unmarked full-root K3 finite-volume classification, all-19 "
                    "literal Vinberg V, and the oriented arithmetic embedding/"
                    "chamber/period-image layer"
                ),
            },
        },
        "oriented_arithmetic_period_classification": {
            "status": "closed_for_all_19",
            "claim_level": oriented_period["scope"]["claim_level"],
            "primitive_embedding_orbits_under_O_plus_per_case": 1,
            "complete_chamber_orbits_in_a_fixed_component_under_O_plus_per_case": 1,
            "pointwise_N_polarized_T_period_image": "tilde O^+(T)",
            "pointwise_discriminant_image": "identity",
            "unmarked_ambient_chamber_normalizer": "Gamma^+_{j,C}",
            "unmarked_chamber_normalizer_T_period_image": "O^+(T)",
            "normalizer_warning": (
                "O^+(T) is the T/period image of Gamma^+_{j,C}, not an "
                "identification Gamma^+_{j,C}=O^+(T)"
            ),
            "geometric_monodromy_status": oriented_period["uniform_theorem"][
                "geometric_monodromy_status"
            ],
            "finite_discriminant_quotients": oriented_period[
                "finite_discriminant_quotients"
            ],
            "finite_quadratic_form_provenance": oriented_period[
                "finite_quadratic_form_provenance"
            ],
        },
        "moduli": {
            "irreducible_family_count": 19,
            "irreducibility_scope": (
                "Roulleau proves the full coarse N-polarized moduli space "
                "irreducible; the exact-Picard locus is asserted only nonempty "
                "and of the same dimension, not a separate finite-type "
                "irreducible algebraic moduli space"
            ),
            "dimensions_by_rank": {str(rho): 20 - rho for rho in range(11, 20)},
            "primitive_embedding_orbits_under_full_k3_lattice_isometry": 19,
            "primitive_embedding_orbits_under_oriented_k3_lattice_isometry": 19,
            "orientation_and_monodromy_boundary": (
                "the oriented arithmetic period/change-of-marking images and "
                "chamber orbits are classified; their identification with "
                "geometric family monodromy remains conditional on the standard "
                "moduli-stack/deck-group interpretation"
            ),
        },
        "source_artifacts": {
            "inventory": {
                "path": str(INVENTORY_PATH.relative_to(ROOT)),
                "sha256": _sha256(INVENTORY_PATH),
            },
            "rank15_19_exact_grams": {
                "path": str(HIGH_GRAM_PATH.relative_to(ROOT)),
                "sha256": _sha256(HIGH_GRAM_PATH),
            },
            "rank15_19_literal_vinberg": {
                "path": str(HIGH_VINBERG_PATH.relative_to(ROOT)),
                "sha256": _sha256(HIGH_VINBERG_PATH),
            },
            "rank11_14_literal_vinberg": {
                "path": str(LOW_VINBERG_PATH.relative_to(ROOT)),
                "sha256": _sha256(LOW_VINBERG_PATH),
            },
            "oriented_arithmetic_period_groups": {
                "path": str(ORIENTED_PERIOD_PATH.relative_to(ROOT)),
                "sha256": _sha256(ORIENTED_PERIOD_PATH),
            },
            "pinned_external_source_hashes": PINNED_SOURCE_HASHES,
        },
        "source_check": source_check or {"checked": False},
        "cases": cases,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    parser.add_argument("--check-sources", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()

    if args.check_archive:
        archived = _json(args.output)
        artifact = build_artifact(archived["source_check"])
        if artifact != archived:
            raise SystemExit("recomputed artifact differs from archive")
        print("archive replay: PASS")
        return

    source_check = check_pinned_sources() if args.check_sources else None
    artifact = build_artifact(source_check)
    if not args.no_write:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(json.dumps(artifact["classification"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
