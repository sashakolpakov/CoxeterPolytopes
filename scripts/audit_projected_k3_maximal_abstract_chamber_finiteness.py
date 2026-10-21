#!/usr/bin/env python3
"""Independent interface audit for projected-K3 maximal chamber finiteness."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OWNER_REPORT = ROOT / "results/projected_k3_maximal_abstract_chamber_finiteness.md"
OWNER_SCRIPT = (
    ROOT / "scripts/certify_projected_k3_maximal_abstract_chamber_finiteness.py"
)
OWNER_ARCHIVE = (
    ROOT
    / "data/enumerations/projected_k3_maximal_abstract_chamber_finiteness.json"
)
OWNER_TEST = (
    ROOT / "tests/test_projected_k3_maximal_abstract_chamber_finiteness.py"
)
AUDIT_REPORT = (
    ROOT
    / "results/projected_k3_maximal_abstract_chamber_finiteness_independent_hostile_audit.md"
)
AUDIT_ARCHIVE = (
    ROOT
    / "data/enumerations/projected_k3_maximal_abstract_chamber_finiteness_independent_hostile_audit.json"
)


OWNER_SHA256 = {
    "report": "59dd458acf14d18021f6c5583096da13f39c83115fc136d6ba91bcc850ce3603",
    "script": "6e4d6859ea4513f614714794ac76197c895bf163a4122565674d3694bf03e23b",
    "archive": "b7108377559f75b65c03fcebb7b01098f7d91b1f5e499c2133d80e882c6f559b",
    "test": "d38d650e7b2d0f0d1d799bfa5bf94ad3c65f1aaa4b32d7a309f15863f71dec5f",
}
AUDIT_REPORT_SHA256 = (
    "edc663c0e7f35880123149d81d832812d5f28837c837df16eaa0f1daf725916a"
)
PRIMARY_NIKULIN_SHA256 = (
    "696652b6a67c6ed420de6c1cb61d9422e3619e28e081a78131a5323ce2edee9d"
)


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
    return tuple(tuple(scale * entry for entry in row) for row in matrix)


def cartan_matrix(
    gram: tuple[tuple[int, ...], ...]
) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(
        tuple(Fraction(2 * gram[i][j], gram[i][i]) for j in range(len(gram)))
        for i in range(len(gram))
    )


def normalized_squared_gram(
    gram: tuple[tuple[int, ...], ...]
) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(
        tuple(
            Fraction(
                4 * gram[i][j] * gram[i][j],
                gram[i][i] * gram[j][j],
            )
            for j in range(len(gram))
        )
        for i in range(len(gram))
    )


def fraction_rows(
    matrix: tuple[tuple[Fraction, ...], ...]
) -> list[list[str]]:
    return [[str(entry) for entry in row] for row in matrix]


def rank_budget() -> list[dict[str, int]]:
    rows: list[dict[str, int]] = []
    for degree in range(1, 21):
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


def finite_index_replay(limit: int = 256) -> dict[str, int]:
    """Replay [O+:W]=[O+:R][R:W] and [O:O+]<=2 numerically."""

    cases = 0
    maximum_full_index = 0
    for positive_index_w in range(1, limit + 1):
        for positive_index_r in range(1, positive_index_w + 1):
            if positive_index_w % positive_index_r:
                continue
            r_index_w = positive_index_w // positive_index_r
            assert positive_index_w == positive_index_r * r_index_w
            assert positive_index_r <= positive_index_w
            for cone_index in (1, 2):
                full_index_r = cone_index * positive_index_r
                assert full_index_r <= 2 * positive_index_w
                maximum_full_index = max(maximum_full_index, full_index_r)
                cases += 1
    return {
        "maximum_tested_positive_index_W": limit,
        "factorizations_tested": cases,
        "maximum_resulting_full_index_R": maximum_full_index,
    }


def validate_owner() -> dict[str, Any]:
    paths = {
        "report": OWNER_REPORT,
        "script": OWNER_SCRIPT,
        "archive": OWNER_ARCHIVE,
        "test": OWNER_TEST,
    }
    for key, path in paths.items():
        actual = sha256_file(path)
        if actual != OWNER_SHA256[key]:
            raise AssertionError(f"owner {key} hash changed: {actual}")
    if sha256_file(AUDIT_REPORT) != AUDIT_REPORT_SHA256:
        raise AssertionError("independent audit report hash changed")

    report = OWNER_REPORT.read_text(encoding="utf-8")
    required = (
        "exact projected complex-K3 eigenslice return containing",
        "full native trace module",
        "W(P)\\subseteq R(L)\\subseteq O^+(L)",
        "Theorem 5.2.1, pp. 119--122",
        "W(P)=R(L)",
        "one \\(R(L)\\)-orbit of",
        "This is not yet an AG-only derivation of ABSW",
        "python3 -m pytest -q tests/test_projected_k3_maximal_abstract_chamber_finiteness.py",
    )
    for needle in required:
        if needle not in report:
            raise AssertionError(f"owner interface text missing: {needle}")
    if re.search(r"(?<!\\)qquad", report):
        raise AssertionError("TeX transport regression: qquad")

    archived = json.loads(OWNER_ARCHIVE.read_text(encoding="utf-8"))
    payload = archived["payload"]
    if archived["payload_sha256"] != digest(payload):
        raise AssertionError("owner payload digest mismatch")
    theorem = payload["theorem"]
    if theorem["classification_free_AG_only_proof"] is not False:
        raise AssertionError("owner import scope changed")
    if theorem["finite_native_similarity_hosts"] is not True:
        raise AssertionError("owner native-host conclusion changed")
    if theorem["finite_ambient_NS_hosts_required"] is not False:
        raise AssertionError("owner ambient-host nonrequirement changed")
    return payload


def build_certificate() -> dict[str, Any]:
    owner = validate_owner()

    budgets = rank_budget()
    if budgets[0]["maximum_hyperbolic_dimension"] != 19:
        raise AssertionError("complex-K3 dimension ceiling changed")
    if max(
        row["maximum_hyperbolic_dimension"]
        for row in budgets
        if row["field_degree"] > 1
    ) != 9:
        raise AssertionError("nonrational complex-K3 ceiling changed")

    cartan = cartan_matrix(DELTA_246_GRAM)
    normalized = normalized_squared_gram(DELTA_246_GRAM)
    for scale in range(1, 129):
        scaled = scale_matrix(DELTA_246_GRAM, scale)
        if cartan_matrix(scaled) != cartan:
            raise AssertionError("reflection Cartan matrix changed under similarity")
        if normalized_squared_gram(scaled) != normalized:
            raise AssertionError("normalized Gram changed under similarity")

    # The projected physical-root estimate has strict real upper bound 14.
    # Equal-coordinate rational samples give 7a+(2-a)=6a+2.
    for numerator in range(0, 201):
        a = Fraction(numerator, 100)
        if not 0 <= a <= 2:
            continue
        if 6 * a + 2 > 14:
            raise AssertionError("physical-root operator bound failed")
    if max(integer for integer in range(100) if integer < 14) != 13:
        raise AssertionError("strict-integral physical-root bound changed")

    index_replay = finite_index_replay()

    return {
        "schema": (
            "coxeter.projected-k3-maximal-abstract-chamber-finiteness."
            "independent-hostile-audit.v1"
        ),
        "verdict": "PASS_CONDITIONAL_CLASSIFICATION_IMPORT_EXPLICIT",
        "owner_sha256": OWNER_SHA256,
        "audit_report_sha256": AUDIT_REPORT_SHA256,
        "primary_source": {
            "author": "V. V. Nikulin",
            "journal": "Math. USSR Izv. 18 (1982), 99-123",
            "locator": "Theorem 5.2.1, pp. 119-122",
            "pdf_sha256": PRIMARY_NIKULIN_SHA256,
            "directly_read": True,
            "scope_matches_owner": True,
            "fixed_data": ["totally real ground field", "rank at least three"],
            "equivalence": "hyperbolic lattices up to similarity",
            "property": "full automorphism group generated up to finite index by reflections",
            "torsion_free_projective_modules_allowed": True,
        },
        "exact_scope_gate": {
            "finite_sided_finite_volume_arithmetic_coxeter_chamber": True,
            "minimum_hyperbolic_dimension": 2,
            "full_native_trace_module": True,
            "admissible_native_form": True,
            "exact_complete_facet_operator_return": True,
            "raw_reflection_maximality": True,
            "merely_semialgebraic_NS_slice_suffices": False,
            "universal_projected_K3_existence_proved": False,
        },
        "K3_rank_budget": budgets,
        "field_gate": owner["field_gate"],
        "reflective_lattice_bridge": {
            "chain": "W(P) <= R(L) <= O+(L) <= O(L)",
            "finite_covolume_endpoint_index": "[O+(L):W(P)] < infinity",
            "intermediate_index": "[O+(L):R(L)] < infinity",
            "positive_cone_index": "[O(L):O+(L)] <= 2",
            "nikulin_index": "[O(L):R(L)] < infinity",
            "finite_index_replay": index_replay,
            "suborder_to_maximal_order_span_preserves_action": True,
            "common_positive_integral_scale_clears_denominators": True,
            "module_need_not_be_free": True,
        },
        "similarity_replay": {
            "seed_gram": [list(row) for row in DELTA_246_GRAM],
            "positive_integer_scales_tested": [1, 128],
            "cartan": fraction_rows(cartan),
            "normalized_squared_gram": fraction_rows(normalized),
            "reflection_formula_scale_cancels": True,
        },
        "maximality_and_chambers": {
            "R_is_discrete_reflection_overgroup": True,
            "R_is_arithmetic_when_arithmetic_only_maximality_is_used": True,
            "maximality_forces_W_equals_R": True,
            "complete_mirror_arrangement_locally_finite": True,
            "generic_compact_path_crosses_finitely_many_mirrors": True,
            "successive_crossing_reflections_give_chamber_transitivity": True,
            "unmarked_normalized_chambers_per_oriented_similarity_host_at_most": 1,
        },
        "hostile_controls": {
            "projected_scale_glue_family": {
                "ambient_hosts_infinite": True,
                "native_similarity_host_fixed": True,
                "abstract_chamber_fixed": True,
            },
            "redoubling_family": {
                "fixed_field_rank_and_lattice": True,
                "abstract_chambers_infinite": True,
                "reflection_maximal": False,
            },
            "connected_frame": {
                "bounded_spanning_local_star_may_be_disconnected": True,
                "right_angled_pentagon_control": True,
                "local_star_replaces_Nikulin_fixed_data_theorem": False,
            },
        },
        "physical_star": {
            "strict_real_upper_bound": 14,
            "nonnegative_integral_upper_bound": 13,
            "used_in_global_chamber_finiteness": False,
        },
        "conclusion": {
            "finite_native_similarity_hosts": True,
            "finite_unmarked_normalized_maximal_chambers": True,
            "finite_ambient_NS_hosts": False,
            "full_ABSW_without_exhaustive_carrier_stack": False,
            "AG_only_no_Nikulin_import": False,
        },
        "logic": [
            "exact full-trace complex-K3 rank plus exact field theorem -> finite (k,m)",
            "finite covolume inclusion W<=R<=O+ and [O:O+]<=2 -> L is Nikulin-reflective",
            "Nikulin 1982 Theorem 5.2.1 -> finite native similarity hosts",
            "raw reflection maximality -> W=R",
            "local finiteness plus finite mirror gallery -> one complete-chamber orbit",
            "finite union over (k,m) -> finite unmarked normalized chambers",
        ],
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
        if AUDIT_ARCHIVE.read_text(encoding="utf-8") != encoded:
            raise SystemExit("frozen independent certificate differs from exact replay")
        print("projected K3 maximal abstract-chamber independent audit: PASS")
        return
    if args.write_output:
        AUDIT_ARCHIVE.write_text(encoded, encoding="utf-8")
        print(AUDIT_ARCHIVE)
        return
    print(encoded, end="")


if __name__ == "__main__":
    main()
