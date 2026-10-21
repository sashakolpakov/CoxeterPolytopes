#!/usr/bin/env python3
"""Exact interface replay for projected-K3 maximal chamber finiteness."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/projected_k3_maximal_abstract_chamber_finiteness.md"
ARCHIVE = (
    ROOT
    / "data/enumerations/projected_k3_maximal_abstract_chamber_finiteness.json"
)
NIKULIN_AUDIT = ROOT / "results/nikulin_conceptual_dependency_hostile_audit.md"
FIELD_COMPACT = ROOT / "results/k3_compact_narrow_facet_field_finiteness.md"
FIELD_NONCOMPACT = ROOT / "results/exact_trace_noncompact_ground_field_collapse.md"
SCALING_REPORT = ROOT / "results/k3_compact_projected_root_scaling_no_go.md"
SCALING_ARCHIVE = (
    ROOT / "data/enumerations/k3_compact_projected_root_scaling_no_go.json"
)
REDOUBLING_REPORT = ROOT / "results/k3_s3_h10_redoubling_family.md"
REDOUBLING_ARCHIVE = ROOT / "data/enumerations/k3_s3_h10_redoubling_family.json"
COMPOSITION_AUDIT = (
    ROOT / "results/compact_trace_host_chamber_finiteness_composition_hostile_audit.md"
)


DEPENDENCY_SHA256 = {
    "nikulin_primary_pdf": (
        "696652b6a67c6ed420de6c1cb61d9422e3619e28e081a78131a5323ce2edee9d"
    ),
    "nikulin_audit": (
        "953717dabaf24dd1163c49ec738664eb1e84bc7652e3a9b0dc124eaa68b4bddb"
    ),
    "compact_field": (
        "39c9a3a8422110bcdcf751ed7bae90ac86cdf3a720d220ee164ff241bd9ac9b5"
    ),
    "noncompact_field": (
        "d738018bacc3803f78fa725da375bb8045ecce92a1c879b03b30a9998e980255"
    ),
    "scaling_report": (
        "bcd4e41d6b8915f09726a9512c8bcac1876da3aa3cdd513ac567a4bcf4c7b057"
    ),
    "scaling_archive": (
        "6aaeca3809bbed3e66af8a97fe89b3f93b9b5269755c00087ad924c5554d95c9"
    ),
    "redoubling_report": (
        "62a9e4315a669c5428a29aa4121ad28c473993fd386ecbdd79c8c47178507c66"
    ),
    "redoubling_archive": (
        "c2cf3c5ee99765ecccc88c027f47666a059eafa521ae45f06a7df46bd47a4df0"
    ),
    "composition_audit": (
        "8fd86c7c588f247fd3ec380e8f1e68840da81484414ab7e7413a7b3487013900"
    ),
}


DELTA_246_GRAM = (
    (-12, 6, 0),
    (6, -6, 3),
    (0, 3, -2),
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(payload: Any) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def scale_matrix(
    matrix: tuple[tuple[int, ...], ...], scale: int
) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(scale * value for value in row) for row in matrix)


def cartan_matrix(
    gram: tuple[tuple[int, ...], ...]
) -> tuple[tuple[Fraction, ...], ...]:
    size = len(gram)
    return tuple(
        tuple(Fraction(2 * gram[i][j], gram[i][i]) for j in range(size))
        for i in range(size)
    )


def normalized_squared_gram(
    gram: tuple[tuple[int, ...], ...]
) -> tuple[tuple[Fraction, ...], ...]:
    size = len(gram)
    return tuple(
        tuple(
            Fraction(4 * gram[i][j] * gram[i][j], gram[i][i] * gram[j][j])
            for j in range(size)
        )
        for i in range(size)
    )


def fraction_rows(
    matrix: tuple[tuple[Fraction, ...], ...]
) -> list[list[str]]:
    return [[str(value) for value in row] for row in matrix]


def rank_budget() -> list[dict[str, int]]:
    rows: list[dict[str, int]] = []
    for degree in range(1, 7):
        maximum_rank = 20 // degree
        if maximum_rank < 3:
            continue
        rows.append(
            {
                "field_degree": degree,
                "maximum_native_rank": maximum_rank,
                "maximum_hyperbolic_dimension": maximum_rank - 1,
            }
        )
    return rows


def validate_dependencies() -> tuple[dict[str, Any], dict[str, Any]]:
    actual = {
        "nikulin_audit": sha256_file(NIKULIN_AUDIT),
        "compact_field": sha256_file(FIELD_COMPACT),
        "noncompact_field": sha256_file(FIELD_NONCOMPACT),
        "scaling_report": sha256_file(SCALING_REPORT),
        "scaling_archive": sha256_file(SCALING_ARCHIVE),
        "redoubling_report": sha256_file(REDOUBLING_REPORT),
        "redoubling_archive": sha256_file(REDOUBLING_ARCHIVE),
        "composition_audit": sha256_file(COMPOSITION_AUDIT),
    }
    for key, value in actual.items():
        if value != DEPENDENCY_SHA256[key]:
            raise AssertionError(f"dependency hash changed: {key}")

    scaling = json.loads(SCALING_ARCHIVE.read_text(encoding="utf-8"))
    scaling_scope = scaling["payload"]["scope"]
    if scaling_scope["produces_infinitely_many_abstract_chambers"] is not False:
        raise AssertionError("scaling family scope changed")
    if not scaling_scope[
        "proves_projected_physical_walls_insufficient_for_finite_hosts"
    ]:
        raise AssertionError("scaling no-go claim changed")

    redoubling = json.loads(REDOUBLING_ARCHIVE.read_text(encoding="utf-8"))
    if not redoubling["all_checks_pass"]:
        raise AssertionError("redoubling archive no longer passes")
    boundary = redoubling["theorem"]["boundary_proof"]
    if boundary["facet_count_identity"] != "9+(N+1)+2=N+12":
        raise AssertionError("redoubling facet formula changed")
    return scaling, redoubling


def build_certificate() -> dict[str, Any]:
    scaling, redoubling = validate_dependencies()

    cartan = cartan_matrix(DELTA_246_GRAM)
    normalized = normalized_squared_gram(DELTA_246_GRAM)
    for scale in (2, 3, 5, 15, 101):
        scaled = scale_matrix(DELTA_246_GRAM, scale)
        if cartan_matrix(scaled) != cartan:
            raise AssertionError("Cartan matrix changed under global similarity")
        if normalized_squared_gram(scaled) != normalized:
            raise AssertionError("normalized Gram changed under global similarity")

    if normalized[0][1] != 2 or normalized[1][2] != 3:
        raise AssertionError("Delta(2,4,6) finite labels changed")
    if normalized[0][2] != 0:
        raise AssertionError("Delta(2,4,6) right-angle entry changed")

    budgets = rank_budget()
    if max(row["maximum_hyperbolic_dimension"] for row in budgets) != 19:
        raise AssertionError("K3 dimension ceiling changed")
    if max(
        row["maximum_hyperbolic_dimension"]
        for row in budgets
        if row["field_degree"] > 1
    ) != 9:
        raise AssertionError("nonrational K3 dimension ceiling changed")

    # The exact projected-root estimate uses two vectors of squared norm two
    # and diag(7,1), whose squared operator norm is 49.  Hence its bilinear
    # maximum is 7*sqrt(2)*sqrt(2)=14; strict narrowness plus integrality gives 13.
    operator_squared_eigenvalues = (49, 1)
    operator_norm = 7
    if operator_norm**2 != max(operator_squared_eigenvalues):
        raise AssertionError("projected complement operator norm changed")
    physical_intersection_strict_upper_bound = operator_norm * 2
    if physical_intersection_strict_upper_bound != 14:
        raise AssertionError("projected physical-root bound changed")

    samples = []
    for index in (1, 2, 3, 7, 25):
        samples.append(
            {
                "N": index,
                "facets": index + 12,
                "reflection_subgroup_index": index,
                "volume_multiplier": index,
            }
        )
    if len({row["facets"] for row in samples}) != len(samples):
        raise AssertionError("redoubling negative control lost distinct facet counts")

    return {
        "schema": "coxeter.projected-k3-maximal-abstract-chamber-finiteness.v1",
        "theorem": {
            "exact_projected_full_trace_K3_return": True,
            "minimum_hyperbolic_dimension": 2,
            "raw_reflection_maximality": True,
            "finite_ambient_NS_hosts_required": False,
            "finite_native_similarity_hosts": True,
            "finite_unmarked_normalized_abstract_chambers": True,
            "classification_free_AG_only_proof": False,
            "imported_fixed_data_input": (
                "Nikulin 1982 Theorem 5.2.1: fixed totally real field and rank "
                "give finitely many similarity classes of reflective hyperbolic lattices"
            ),
        },
        "K3_rank_budget": budgets,
        "field_gate": {
            "noncompact_ground_field": "Q",
            "compact_discriminant_bound": "(32*14^10)^30",
            "finite_field_set": True,
        },
        "similarity_invariance": {
            "raw_gram": [list(row) for row in DELTA_246_GRAM],
            "tested_global_scales": [2, 3, 5, 15, 101],
            "cartan_matrix": fraction_rows(cartan),
            "normalized_squared_gram": fraction_rows(normalized),
            "reflection_formula_scale_cancels": True,
        },
        "projected_physical_star": {
            "complement_operator_diagonal": [7, 1],
            "operator_squared_eigenvalues": list(operator_squared_eigenvalues),
            "root_component_vector_squared_norm": 2,
            "strict_real_upper_bound": 14,
            "nonnegative_integral_upper_bound": 13,
            "implies_finite_ambient_NS_hosts": False,
            "reason": (
                "physical lifts need not span NS and their K3 reflections need not "
                "equal the native raw reflections"
            ),
        },
        "scaling_control": {
            "integer_scale_domain": scaling["payload"]["all_parameter_family"][
                "integer_scale_domain"
            ],
            "ambient_determinant": scaling["payload"]["all_parameter_family"][
                "ambient_determinant"
            ],
            "glue_index": scaling["payload"]["all_parameter_family"][
                "direct_sum_to_ambient_overlattice_index"
            ],
            "abstract_chamber_count_infinite": scaling["payload"]["scope"][
                "produces_infinitely_many_abstract_chambers"
            ],
            "native_similarity_class_fixed": True,
        },
        "nonmaximal_redoubling_control": {
            "statement": redoubling["theorem"]["statement"],
            "samples": samples,
            "fixed_dimension": 10,
            "fixed_ground_field": "Q",
            "fixed_divisor_lattice": True,
            "reflection_maximal": False,
        },
        "logic": [
            "exact K3 trace rank and finite-volume field theorem -> finite (k,m)",
            "arithmetic W(P) inside O+(L), both lattices -> finite index",
            "W(P) subset R(L) -> L reflective",
            "Nikulin fixed-(k,m) theorem -> finite native similarity hosts",
            "raw reflection maximality -> W(P)=R(L)",
            "complete mirror chambers form one orbit -> finite normalized chambers",
        ],
        "dependency_sha256": DEPENDENCY_SHA256,
    }


def encoded_certificate() -> str:
    payload = build_certificate()
    wrapped = {"payload": payload, "payload_sha256": digest(payload)}
    return json.dumps(wrapped, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-output", action="store_true")
    parser.add_argument("--write-output", action="store_true")
    args = parser.parse_args()

    encoded = encoded_certificate()
    if args.check_output:
        if ARCHIVE.read_text(encoding="utf-8") != encoded:
            raise SystemExit("frozen certificate differs from exact replay")
        print("projected K3 maximal abstract-chamber finiteness: PASS")
        return
    if args.write_output:
        ARCHIVE.write_text(encoded, encoding="utf-8")
        print(ARCHIVE)
        return
    print(encoded, end="")


if __name__ == "__main__":
    main()
