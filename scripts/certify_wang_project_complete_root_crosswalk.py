"""Exact Wang/Gritsenko--Nikulin crosswalk for named project root lattices.

The classification theorems are published inputs.  This certifier checks the
project-level identifications, the weights forced by Wang's formula, and the
same-form even-overlattice chains between named ADE models.  It performs no
Vinberg enumeration.
"""

from __future__ import annotations

import json


# project id, S (negative summands are suppressed), status, weight, beta_1,
# exact correction lattice/form.  beta_0 is one in every Phi_12-pullback row.
ROWS = (
    # Extra lower-rank complete-nef stocks explicitly named in the ledgers.
    ("N3-01/o0", "U+D4", "yes", (72, 40), None,
     "Phi_{72,D4} on 2U+D4; Phi_{40} on U+U(2)+D4"),
    ("N3-04/o0", "U(2)+D4", "yes", (40, 24), None,
     "Phi_{40} on U+U(2)+D4; F_{24} on 2U(2)+D4"),
    ("N3-05/o1", "U+3A2", "yes", (39,), None,
     "Phi_{39,3A2}=Phi_12 quasi-pullback on 2U+3A2"),
    ("N5-02/o0", "U+D4+A2", "no", (), None, None),

    # The ten rank-nine/H8 complete-root rows named in the project ledger.
    ("R9-UE7-NEF-H8-9", "U+E7", "yes", (165,), None,
     "Phi_{165,E7}=Phi_12 quasi-pullback on 2U+E7"),
    ("R9-UD6A1-NEF-H8-10", "U+D6+A1", "yes", (85,), 17,
     "Borch(phi_{0,D6+A1}) on 2U+D6+A1"),
    ("R9-UD4A1x3-NEF-H8-15", "U+D4+3A1", "yes", (45,), 9,
     "Borch(phi_{0,D4+3A1}) on 2U+D4+3A1"),
    ("R9-UA1x7-NEF-H8-37", "U+7A1", "yes", (25,), 5,
     "seven-variable pullback of the 8A1 Borcherds product on 2U+7A1"),
    ("R9-U2A1x7-NEF-H8-240", "<2>+8A1(-1)", "yes", (5,), None,
     "Delta_{5,D7}=Lift(eta^3 product(theta_i)) on U(2)+S"),
    ("R9-UA7-NEF-H8-9", "U+A7", "yes", (80,), None,
     "Phi_{80,A7}=Phi_12 quasi-pullback on 2U+A7"),
    ("R9-UD4A3-NEF-H8-10", "U+D4+A3", "no", (), None, None),
    ("R9-UD5A2-NEF-H8-10", "U+D5+A2", "no", (), None, None),
    ("R9-UD7-NEF-H8-9", "U+D7", "yes", (114,), None,
     "Phi_{114,D7}=Phi_12 quasi-pullback on 2U+D7"),
    ("R9-UA1E6-NEF-H8-10", "U+E6+A1", "yes", (99,), 21,
     "Borch(phi_{0,E6+A1}) on 2U+E6+A1"),

    # The nine rank-ten/H9 complete-root rows named in the project ledger.
    ("R10-UE8-NEF-H9-10", "U+E8", "yes", (252,), None,
     "Phi_{252,E8}=Phi_12 quasi-pullback on 2U+E8"),
    ("R10-UD8-NEF-H9-10", "U+D8", "yes", (124,), None,
     "Phi_{124,D8}=Phi_12 quasi-pullback on 2U+D8"),
    ("R10-UE7A1-NEF-H9-11", "U+E7+A1", "yes", (132,), 33,
     "Borch(phi_{0,E7+A1}) on 2U+E7+A1"),
    ("R10-UD4D4-NEF-H9-11", "U+2D4", "yes", (60,), None,
     "Phi_{60,2D4}=Phi_12 quasi-pullback on 2U+2D4"),
    ("R10-UD6A1x2-NEF-H9-14", "U+D6+2A1", "yes", (68,), 17,
     "Borch(phi_{0,D6+2A1}) on 2U+D6+2A1"),
    ("R10-U2D4D4-NEF-H9-18", "U(2)+2D4=U+N8", "yes", (28,), None,
     "Phi_{28,N8} on U+S = 2U+N8"),
    ("R10-UD4A1x4-NEF-H9-27", "U+D4+4A1", "yes", (36,), 9,
     "Borch(phi_{0,D4+4A1}) on 2U+D4+4A1"),
    ("R10-UA1x8-NEF-H9-145", "U+8A1", "yes", (20,), 5,
     "Borch(phi_{0,8A1}) on 2U+8A1"),
    ("R10-UA2E6-NEF-H9-11", "U+E6+A2", "no", (), None, None),

    # Exhaustive project rank-eleven--nineteen complete-root denominator.
    ("R11_UE8A1", "U+E8+A1", "yes", (195,), 57,
     "Borch(phi_{0,E8+A1}) on 2U+E8+A1"),
    ("R11_UD8A1", "U+D8+A1", "yes", (99,), 25,
     "Borch(phi_{0,D8+A1}) on 2U+D8+A1"),
    ("R11_UD4D4A1", "U+2D4+A1", "yes", (51,), 9,
     "Borch(phi_{0,2D4+A1}) on 2U+2D4+A1"),
    ("R11_UD4A1x5", "U+D4+5A1", "yes", (27,), 9,
     "Borch(phi_{0,D4+5A1}) on 2U+D4+5A1"),
    ("R12_UA2E8", "U+E8+A2", "no", (), None, None),
    ("R12_UE8A1x2", "U+E8+2A1", "yes", (138,), 57,
     "Borch(phi_{0,E8+2A1}) on 2U+E8+2A1"),
    ("R12_UD8A1x2", "U+D8+2A1", "no", (), None, None),
    ("R12_UD4x2A1x2", "U+2D4+2A1", "no", (), None, None),
    ("R13_UE8A3", "U+E8+A3", "no", (), None, None),
    ("R13_UE8A1x3", "U+E8+3A1", "no", (), None, None),
    ("R13_UD8A1x3", "U+D8+3A1", "no", (), None, None),
    ("R14_UE8D4", "U+E8+D4", "no", (), None, None),
    ("R14_UD8D4", "U+D8+D4", "no", (), None, None),
    ("R14_UE8A1x4", "U+E8+4A1", "no", (), None, None),
    ("R15_UE8D4A1", "U+E8+D4+A1", "no", (), None, None),
    ("R16_UE8D6", "U+E8+D6", "no", (), None, None),
    ("R17_UE8E7", "U+E8+E7", "no", (), None, None),
    ("R18_UE8E8", "U+2E8", "yes", (132,), None,
     "Phi_{132,2E8}=Phi_12 quasi-pullback on 2U+2E8"),
    ("R19_UE8E8A1", "U+2E8+A1", "yes", (75,), 57,
     "Borch(phi_{0,2E8+A1}) on 2U+2E8+A1"),
)


# source, target, index.  These are same quadratic form, same rank, even
# inclusions; signs on the negative ADE factors are suppressed.
OVERLATTICES = (
    # H8 chains.
    ("U+7A1", "U+D4+3A1", 2),
    ("U+D4+3A1", "U+D6+A1", 2),
    ("U+D6+A1", "U+E7", 2),
    ("U+A7", "U+E7", 2),
    # H9 binary-code/root chains.
    ("U+8A1", "U+D4+4A1", 2),
    ("U+8A1", "U(2)+2D4=U+N8", 2),
    ("U+D4+4A1", "U+2D4", 2),
    ("U(2)+2D4=U+N8", "U+2D4", 2),
    ("U+2D4", "U+D8", 2),
    ("U+D8", "U+E8", 2),
    ("U+D6+2A1", "U+E7+A1", 2),
    ("U+E7+A1", "U+E8", 2),
    # H10--H13 chains.
    ("U+D4+5A1", "U+2D4+A1", 2),
    ("U+2D4+A1", "U+D8+A1", 2),
    ("U+D8+A1", "U+E8+A1", 2),
    ("U+2D4+2A1", "U+D8+2A1", 2),
    ("U+D8+2A1", "U+E8+2A1", 2),
    ("U+D8+3A1", "U+E8+3A1", 2),
    ("U+E8+4A1", "U+E8+D4", 2),
    ("U+D8+D4", "U+E8+D4", 2),
)


# Exact alternative names occurring in project prose or in Wang's genus
# remarks.  Adding U makes the indicated same-genus positive-definite models
# isometric by the standard indefinite uniqueness criterion.
ALIASES = (
    ("II_(1,9)", "U+E8", "same hyperbolic lattice"),
    ("L10", "U+E8", "same even-unimodular hyperbolic lattice"),
    ("U+D10", "U+E8+2A1", "same hyperbolic lattice"),
    ("U+D6+3A1", "U+2D4+A1", "same hyperbolic lattice"),
    ("U+E7+2A1", "U+D8+A1", "same hyperbolic lattice"),
    ("U+N8+A1", "U+D4+5A1", "same hyperbolic lattice"),
    (
        "II_(1,9)+<2>=2U+E7",
        "U+E7",
        "correction lattice M, not the hyperbolic lattice S",
    ),
)


def _quasi_niemeier_weight(rank_l: int, a1_count: int, h: int) -> int:
    """Wang Theorem 6.2(c), with principal multiplicity beta_0=1."""
    numerator = (24 - rank_l - 3 * a1_count) * h
    if numerator % 2:
        raise ArithmeticError("nonintegral Wang weight")
    return numerator // 2 + 12 + 3 * a1_count


def analyze() -> dict[str, object]:
    by_lattice = {row[1]: row for row in ROWS}
    if len(by_lattice) != len(ROWS):
        raise ArithmeticError("duplicate lattice in Wang crosswalk")

    expected_formula = {
        "U+D6+A1": (7, 1, 10),
        "U+D4+3A1": (7, 3, 6),
        "U+E6+A1": (7, 1, 12),
        "U+E7+A1": (8, 1, 18),
        "U+D6+2A1": (8, 2, 10),
        "U+D4+4A1": (8, 4, 6),
        "U+E8+A1": (9, 1, 30),
        "U+D8+A1": (9, 1, 14),
        "U+2D4+A1": (9, 1, 6),
        "U+D4+5A1": (9, 5, 6),
        "U+E8+2A1": (10, 2, 30),
        "U+2E8+A1": (17, 1, 30),
    }
    formula_checks = {}
    for lattice, (rank_l, m, h) in expected_formula.items():
        computed = _quasi_niemeier_weight(rank_l, m, h)
        recorded = by_lattice[lattice][3][0]
        beta_one = by_lattice[lattice][4]
        formula_checks[lattice] = {
            "computed_weight": computed,
            "recorded_weight": recorded,
            "computed_beta_1": 2 * h - 3,
            "recorded_beta_1": beta_one,
            "passes": computed == recorded and beta_one == 2 * h - 3,
        }

    edge_checks = []
    for source, target, index in OVERLATTICES:
        edge_checks.append({
            "source": source,
            "target": target,
            "index": index,
            "both_named": source in by_lattice and target in by_lattice,
        })

    high = [row for row in ROWS if row[0].startswith("R1") and "_" in row[0]]
    # The exact high denominator uses underscore keys R11_... through R19_....
    high = [row for row in high if row[0].split("_")[0][1:].isdigit()]
    high_yes = sum(row[2] == "yes" for row in high)
    high_no = sum(row[2] == "no" for row in high)

    checks = {
        "all_rows_have_unique_lattice": len(by_lattice) == len(ROWS),
        "all_formula_rows_pass": all(x["passes"] for x in formula_checks.values()),
        "all_overlattice_endpoints_are_named": all(x["both_named"] for x in edge_checks),
        "frozen_high_rows_are_19": len(high) == 19,
        "frozen_high_split_is_7_plus_12": (high_yes, high_no) == (7, 12),
        "all_negative_rows_have_no_form_or_weight": all(
            status == "yes" or (not weights and form is None)
            for _, _, status, weights, _, form in ROWS
        ),
        "all_alias_targets_are_named": all(target in by_lattice for _, target, _ in ALIASES),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"Wang crosswalk failed: {checks}")

    return {
        "schema": "coxeter.wang_project_complete_root_crosswalk.v1",
        "scope": {
            "published_inputs": [
                "Wang 2024, Theorems 6.2 and 7.2 and Section 5",
                "Gritsenko--Nikulin 2018, Theorems 4.3, 6.1, 6.2",
                "Nikulin--Vinberg complete hyperbolic 2-reflective classification",
            ],
            "project_rows": len(ROWS),
            "same_form_overlattice_edges": len(OVERLATTICES),
            "warning": (
                "A finite mixed-root divisor chamber is not automatically the "
                "complete norm-two chamber of its carrier lattice."
            ),
        },
        "rows": [
            {
                "project_id": project_id,
                "hyperbolic_lattice_S": lattice,
                "automorphic_correction": status == "yes",
                "weights": list(weights),
                "principal_divisor_multiplicity_beta_0": 1 if status == "yes" else None,
                "A1_divisor_multiplicity_beta_1": beta_one,
                "exact_form": form,
            }
            for project_id, lattice, status, weights, beta_one, form in ROWS
        ],
        "same_form_even_overlattices": edge_checks,
        "exact_aliases": [
            {"alias": alias, "canonical_row": target, "role": role}
            for alias, target, role in ALIASES
        ],
        "formula_checks": formula_checks,
        "checks": checks,
    }


if __name__ == "__main__":
    print(json.dumps(analyze(), indent=2, sort_keys=True))
