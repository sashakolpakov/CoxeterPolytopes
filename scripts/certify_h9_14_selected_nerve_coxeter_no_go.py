#!/usr/bin/env python3
"""Freeze the re-audited compact no-go for seven displayed H9/14 nerves."""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter
from typing import Mapping, Sequence

import sympy as sp

from coxeter_polytopes.combinatorics.h9_14_record32_coxeter_labels import (
    audit_record32_coxeter_labels,
)
from coxeter_polytopes.combinatorics.h9_14_selected_nerve_coxeter_no_go import (
    SelectedNerveAudit,
    audit_selected_nerve_coxeter_no_gos,
)
from scripts.certify_h9_14_record35_coxeter_labels import (
    verify as verify_record35,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/h9_14_selected_nerve_coxeter_no_go.json"
SOURCE_PATHS = (
    "src/coxeter_polytopes/combinatorics/h9_14_selected_nerve_coxeter_no_go.py",
    "src/coxeter_polytopes/combinatorics/h9_14_disjoint_lanner_gate.py",
    "src/coxeter_polytopes/combinatorics/h9_14_record32_coxeter_labels.py",
    "src/coxeter_polytopes/combinatorics/h9_14_coxeter_labels.py",
    "src/coxeter_polytopes/combinatorics/h9_14_record32.py",
    "src/coxeter_polytopes/combinatorics/h9_14_record43.py",
    "src/coxeter_polytopes/combinatorics/h9_14_record44.py",
    "src/coxeter_polytopes/combinatorics/h9_14_record47.py",
    "scripts/independent_h9_14_selected_nerve_coxeter_no_go_audit.py",
)
INPUT_PATHS = (
    "data/enumerations/h9_14_facet_completion.json",
    "data/enumerations/h9_14_multiwedge_polytopality.json",
    "data/enumerations/h9_14_record32_polytopality.json",
    "data/enumerations/h9_14_record32_coxeter_labels.json",
    "data/enumerations/h9_14_record35_coxeter_labels.json",
)
EXPECTED_CERTIFICATE_SHA256 = (
    "c2dfcbfc03dc5dbae34938ab32f3a1c1111fd729d7844fd76d04dfc6cfa44de8"
)


def _file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _certificate_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    return sha256(
        json.dumps(
            payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode("utf-8")
    ).hexdigest()


def _path_records(paths: Sequence[str]) -> list[dict[str, str]]:
    return [
        {"path": relative, "sha256": _file_sha256(ROOT / relative)}
        for relative in paths
    ]


def _gate_summary(gate: object) -> dict[str, object]:
    passing_by_pair = gate.passing_count_by_cross_pair
    return {
        "left_support": list(gate.left_support),
        "right_support": list(gate.right_support),
        "left_state_count": gate.left_state_count,
        "right_state_count": gate.right_state_count,
        "sparse_state_count": gate.sparse_state_count,
        "cross_pair_count": len(gate.cross_pairs),
        "actual_restriction_size_histogram": {
            str(size): count
            for size, count in gate.unique_actual_restriction_size_histogram
        },
        "single_edge_label": 3,
        "single_edge_test_count": gate.label_three_test_count,
        "single_edge_passing_count": len(gate.passing_trials),
        "passing_count_per_cross_pair_histogram": {
            str(count): multiplicity
            for count, multiplicity in sorted(
                Counter(value for _edge, value in passing_by_pair).items()
            )
        },
        "first_bad_restriction_size_histogram": {
            str(size): count
            for size, count in gate.first_bad_restriction_size_histogram
        },
    }


def _record_summary(audit: SelectedNerveAudit) -> dict[str, object]:
    return {
        "record_index": audit.record_index,
        "rank_nine_vertex_count": audit.vertex_count,
        "selected_nerve_only": True,
        "support_gates": [_gate_summary(gate) for gate in audit.support_gates],
        "stronger_single_edge_survivors": {
            str(label): count
            for label, count in audit.stronger_single_edge_survivors
        },
        "raw_two_edge_trial_count": audit.raw_double_edge_trial_count,
        "monotonicity_residual_two_edge_trial_count": (
            audit.monotonicity_residual_double_trial_count
        ),
        "two_edge_passing_count": audit.double_edge_passing_count,
        "exceptional_class_histogram": {
            name: count for name, count in audit.exceptional_class_histogram
        },
        "exact_index_two_witnesses": [
            witness.to_record() for witness in audit.exact_index_two_witnesses
        ],
        "normalizing_support_symmetry_orders": list(
            audit.normalizing_support_symmetry_orders
        ),
        "unrestricted_compact_coxeter_no_go": (
            audit.unrestricted_compact_coxeter_no_go
        ),
    }


def _exact_radical_sign_proofs() -> dict[str, object]:
    r2, r5, r10 = sp.sqrt(2), sp.sqrt(5), sp.sqrt(10)
    bounds = {
        "sqrt(2)": (sp.Rational(7, 5), sp.Rational(3, 2), 2),
        "sqrt(5)": (sp.Rational(11, 5), sp.Rational(9, 4), 5),
        "sqrt(10)": (sp.Rational(79, 25), sp.Rational(19, 6), 10),
    }
    if not all(lower**2 < radicand < upper**2 for lower, upper, radicand in bounds.values()):
        raise AssertionError("an exact rational radical enclosure changed")
    decisive_rational_bounds = {
        "record42_Delta8_upper": 8 - 4 * sp.Rational(11, 5),
        "record42_Delta9_lower": (
            -11
            - 8 * sp.Rational(3, 2)
            + 5 * sp.Rational(11, 5)
            + 4 * sp.Rational(79, 25)
        ),
        "A5_A_Delta4_lower": 6 - 2 * sp.Rational(9, 4),
        "A5_A_Delta5_upper": 4 - 2 * sp.Rational(11, 5),
        "A5_A_Delta6_numerator_upper": 9 - 5 * sp.Rational(11, 5),
        "A5_A_Delta9_upper": 26 - 12 * sp.Rational(11, 5),
        "A5_A_Delta10_numerator_lower": 25 - 11 * sp.Rational(9, 4),
        "A5_B_Delta3_lower": 5 - sp.Rational(9, 4),
        "A5_B_Delta4_numerator_lower": 7 - 3 * sp.Rational(9, 4),
        "A5_B_Delta9_upper": 1 - sp.Rational(11, 5),
        "A5_B_Delta10_lower": (
            -23
            - 16 * sp.Rational(3, 2)
            + 11 * sp.Rational(11, 5)
            + 8 * sp.Rational(79, 25)
        )
        / 2,
    }
    expected_bound_signs = {
        "record42_Delta8_upper": "-",
        "record42_Delta9_lower": "+",
        "A5_A_Delta4_lower": "+",
        "A5_A_Delta5_upper": "-",
        "A5_A_Delta6_numerator_upper": "-",
        "A5_A_Delta9_upper": "-",
        "A5_A_Delta10_numerator_lower": "+",
        "A5_B_Delta3_lower": "+",
        "A5_B_Delta4_numerator_lower": "+",
        "A5_B_Delta9_upper": "-",
        "A5_B_Delta10_lower": "+",
    }
    for name, value in decisive_rational_bounds.items():
        observed = "+" if value > 0 else "-" if value < 0 else "0"
        if observed != expected_bound_signs[name]:
            raise AssertionError(f"the exact rational sign bound {name} changed")
    identities = {
        "record42_last_square_difference": sp.simplify(
            (5 * r5 + 4 * r10) ** 2 - (11 + 8 * r2) ** 2
        ),
        "A5_class_B_last_square_difference": sp.simplify(
            (11 * r5 + 8 * r10) ** 2 - (23 + 16 * r2) ** 2
        ),
    }
    if identities != {
        "record42_last_square_difference": 36 + 24 * r2,
        "A5_class_B_last_square_difference": 204 + 144 * r2,
    }:
        raise AssertionError("the exact radical square identities changed")
    return {
        "rational_radical_enclosures": {
            name: {
                "lower": str(lower),
                "upper": str(upper),
                "proof": f"{lower}^2 < {radicand} < {upper}^2",
            }
            for name, (lower, upper, radicand) in bounds.items()
        },
        "decisive_rational_bounds": {
            name: sp.sstr(value) for name, value in decisive_rational_bounds.items()
        },
        "record42": [
            "8-4*sqrt(5)<0 because 5>2^2.",
            (
                "The final minor is positive: both sides are positive and "
                "(5*sqrt(5)+4*sqrt(10))^2-(11+8*sqrt(2))^2="
                "36+24*sqrt(2)>0."
            ),
            (
                "Including Delta_0=1, the signs +++++++-+ have two changes; "
                "exact LDL/Jacobi therefore gives inertia (7,2,0)."
            ),
        ],
        "A5_class_A": [
            (
                "3>sqrt(5), sqrt(5)>2, 5*sqrt(5)>9, "
                "12*sqrt(5)>26, and 25>11*sqrt(5); these follow by "
                "squaring positive rationals (9>5, 5>4, 125>81, "
                "720>676, 625>605)."
            ),
            (
                "Thus the ten nonzero leading-minor signs are ++++-----+; "
                "with Delta_0=1 they give exact inertia (8,2,0)."
            ),
        ],
        "A5_class_B": [
            "5-sqrt(5)>0 and 7-3*sqrt(5)>0 since 25>5 and 49>45.",
            (
                "The final minor is positive: both sides are positive and "
                "(11*sqrt(5)+8*sqrt(10))^2-(23+16*sqrt(2))^2="
                "204+144*sqrt(2)>0."
            ),
            (
                "Together with sqrt(5)>2, 5*sqrt(5)>9, and sqrt(5)>1, "
                "the signs are ++++-----+, hence inertia (8,2,0)."
            ),
        ],
        "machine_checked_exact_expressions": {
            key: sp.sstr(value) for key, value in identities.items()
        },
    }


def generate() -> dict[str, object]:
    support_audits = audit_selected_nerve_coxeter_no_gos()
    if tuple(audit.record_index for audit in support_audits) != (
        30,
        32,
        42,
        43,
        44,
        47,
    ):
        raise AssertionError("the six selected-support audits changed")
    if not all(audit.unrestricted_compact_coxeter_no_go for audit in support_audits):
        raise AssertionError("an unrestricted selected-support theorem disappeared")
    record32 = audit_record32_coxeter_labels()
    if not record32.verified:
        raise AssertionError("the selected record-32 certificate failed")
    record35_document = json.loads(
        (ROOT / INPUT_PATHS[-1]).read_text(encoding="utf-8")
    )
    record35_summary = verify_record35(record35_document)
    if record35_summary["compact_coxeter_realization_exists"] is not False:
        raise AssertionError("the selected record-35 no-go disappeared")

    document: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_seven_selected_nerve_compact_coxeter_no_go_reaudited",
        "claim": (
            "Each of the seven explicitly displayed ordinary-polytopal H9/14 "
            "nerves numbered 30, 32, 35, 42, 43, 44, and 47 is excluded as "
            "the nerve of a compact hyperbolic Coxeter polytope. A hostile "
            "scope re-audit explicitly verifies the support symmetries and "
            "proves that labels outside each decisive principal block cannot "
            "alter its index-two obstruction."
        ),
        "theorem": {
            "dimension": 9,
            "facet_count": 14,
            "selected_nerve_indices": [30, 32, 35, 42, 43, 44, 47],
            "selected_nerve_count": 7,
            "selected_nerve_only": True,
            "compact_coxeter_realization_exists": False,
            "scope": (
                "Exactly the seven frozen/displayed nerve completions. Record "
                "indices name graph shadows, so this certificate does not by "
                "itself exclude another completion of the same graph shadow."
            ),
            "nonclaims": [
                "No graph shadow is excluded merely because one displayed completion is excluded.",
                "No finite-volume noncompact chamber is excluded.",
                "No compact H9 Coxeter polytope with another nerve is excluded.",
                "No new Coxeter polytope is constructed.",
            ],
        },
        "scope_reaudit": {
            "status": "hostile_false_alarm_rejected_by_principal_block_proof",
            "reaudited_indices": [30, 42, 43, 44, 47],
            "normalizing_support_symmetry_orders": {
                "30": [24, 24],
                "42": [24, 120],
                "43": [120],
                "44": [120],
                "47": [120],
            },
            "principal_block_lemma": (
                "Entries incident to nodes outside a principal submatrix do "
                "not change that submatrix or its inertia. Hence an outside "
                "dotted or Lanner component cannot mediate away two negative "
                "directions already present in the selected block."
            ),
            "bad_vertex_lemma": (
                "A failed actual-vertex restriction is a principal subdiagram "
                "of the required elliptic vertex diagram. Adding or strengthening "
                "non-right entries cannot restore ellipticity by nonnegative "
                "adjacency Perron-Frobenius monotonicity."
            ),
            "multiple_support_gate_lemma": (
                "For records 43,44,47, fix any one audited residual A5 support. "
                "The zero/one/at-least-two case split exhausts all cross entries "
                "inside that ten-node block. Edges to nodes outside the fixed "
                "support do not alter the block; an edge associated with an "
                "overlapping support whose endpoints lie inside it is already "
                "one of the audited block entries."
            ),
        },
        "vinberg_tests_used": {
            "minimal_nonfaces": (
                "For a compact chamber, every missing face supports a Lanner "
                "diagram; the complete order-four and order-five placement "
                "tables are exhausted exactly."
            ),
            "actual_vertices": (
                "Every displayed rank-nine vertex must have a finite (elliptic) "
                "Coxeter diagram; the complete actual vertex lists, not only "
                "abstract support intersections, are tested."
            ),
            "global_gram": (
                "A 14-wall compact chamber in H9 needs Gram inertia (9,1,4). "
                "A principal block with two negative eigenvalues contradicts "
                "this by Cauchy interlacing."
            ),
            "monotone_minimal_label_reduction": (
                "In normalized convention G=I-A, A_ij=cos(pi/m_ij)>=0. "
                "The smallest non-right label m=3 minimizes A_ij. Increasing "
                "a label or adding any non-right edge increases A entrywise, "
                "so Perron-Frobenius monotonicity cannot repair a nonelliptic "
                "actual-vertex principal diagram."
            ),
        },
        "case_exhaustion": {
            "no_cross_edge_inside_decisive_block": (
                "The two Lanner supports are orthogonal and their principal "
                "Gram direct sum has negative index two."
            ),
            "locally_forbidden_cross_edge": (
                "The minimal m=3 trial is already nonelliptic in an actual "
                "vertex; stronger or additional entries cannot repair it."
            ),
            "exceptional_cross_edges": (
                "All individually admissible edges, stronger labels, and "
                "monotonicity-residual pairs are exhausted. Every remaining "
                "one-edge block has exact negative index two."
            ),
            "outside_entries": (
                "They cannot alter the decisive principal submatrix, so they "
                "are irrelevant to Cauchy interlacing."
            ),
        },
        "exact_radical_sign_proofs": _exact_radical_sign_proofs(),
        "six_nerve_census": [
            _record_summary(audit) for audit in support_audits
        ],
        "record32_detailed_census": record32.to_record(),
        "record35_chained_certificate": {
            "selected_nerve_only": True,
            "certificate_sha256": record35_document["certificate_sha256"],
            "summary": record35_summary,
        },
        "sources": {
            "code": _path_records(SOURCE_PATHS),
            "frozen_inputs": _path_records(INPUT_PATHS),
        },
    }
    document["certificate_sha256"] = _certificate_digest(document)
    return document


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != (
        "h9_14_seven_selected_nerve_compact_coxeter_no_go_reaudited"
    ):
        raise ValueError("unexpected selected-nerve certificate kind")
    rebuilt = generate()
    if document != rebuilt:
        raise ValueError("selected-nerve certificate did not replay exactly")
    digest = _certificate_digest(document)
    if document.get("certificate_sha256") != digest:
        raise ValueError("selected-nerve certificate digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 not in ("TO_BE_FROZEN", digest):
        raise ValueError("selected-nerve certificate is not the frozen artifact")
    theorem = rebuilt["theorem"]
    return {
        "selected_nerve_indices": theorem["selected_nerve_indices"],
        "selected_nerve_count": theorem["selected_nerve_count"],
        "compact_coxeter_realization_exists": theorem[
            "compact_coxeter_realization_exists"
        ],
        "scope": theorem["scope"],
        "reaudited_indices": rebuilt["scope_reaudit"]["reaudited_indices"],
        "certificate_sha256": digest,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args(argv)
    started = perf_counter()
    if args.verify:
        document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    else:
        document = generate()
    if args.write:
        OUTPUT.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    summary = verify(document)
    summary["runtime_seconds"] = round(perf_counter() - started, 6)
    summary["written"] = args.write
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
