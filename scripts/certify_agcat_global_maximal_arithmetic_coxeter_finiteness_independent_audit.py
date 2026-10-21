#!/usr/bin/env python3
"""Independent interface certificate for the terminal AGCat composition.

The substantive arithmetic calculations are frozen in the owner and hostile
certificates listed below.  This script checks the *composition*: exact
category/return data, exhaustive dimensional field branches, the independent
dimension-to-rank closure, maximality scopes, dependency firewalls, pins, and
local links.  It deliberately does not pretend that string checks prove the
mathematical lemmas audited in the companion report.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
OWNER = ROOT / "results/agcat_global_maximal_arithmetic_coxeter_finiteness.md"

FROZEN_INPUTS = {
    "composite_owner": (
        OWNER,
        "441deb5a379fc59c6b4c1e1e86212f5f84db4e603a4b37259200db43a4508dfd",
    ),
    "h2_owner": (
        ROOT / "results/arithmetic_coxeter_polygon_orbifold_curve_ag_finiteness.md",
        "68ab8d85002106fde4493540ee6c6409247211cdda5364a3181000ae0135b659",
    ),
    "h2_audit": (
        ROOT / "results/arithmetic_coxeter_polygon_orbifold_curve_ag_finiteness_hostile_audit.md",
        "7be0ea8da1eb6dcc97ced1b68f40a47c9eb3ae382b47ab054e77735df4afcfc9",
    ),
    "n_ge_6_owner": (
        ROOT / "results/global_minimality14_integral_arrangement_field_funnel.md",
        "05acee48bd7a6c03e92febc32434cc38bcd3581c0e9044dafec3f07f9dbab678",
    ),
    "n_ge_6_audit": (
        ROOT / "results/global_minimality14_integral_arrangement_field_funnel_hostile_audit.md",
        "e5e2754bf886d922ba50173e8afef67ff7eceffc62b706bb3c4dea4a309cac69",
    ),
    "low_dimension_owner": (
        ROOT / "results/low_dimensional_minimality14_field_descent_to_h3.md",
        "8157c0779705f5f265e75126d10361441bf2e2b504e522d917cbda848b75f019",
    ),
    "low_dimension_audit": (
        ROOT / "results/low_dimensional_minimality14_field_descent_to_h3_independent_hostile_audit.md",
        "74c32b41b593419951ecb6db0c43c3364134cc033e77038175156bffb385b5da",
    ),
    "compact_h3_owner": (
        ROOT / "results/compact_h3_general_fork_arrangement_field_finiteness.md",
        "4a971c64460a7d234114ecb12b793e613b080945850098da37b903d2c50876d5",
    ),
    "compact_h3_audit": (
        ROOT / "results/compact_h3_general_fork_arrangement_field_finiteness_independent_hostile_audit.md",
        "7bded0c2dc46c1652e06dd15ab691685736db4b327c6fb8a218ca84066cef996",
    ),
    "native_host_owner": (
        ROOT / "results/ag_field_funnel_bounded_frame_native_host_closure.md",
        "7ba392dd833358379ea92d1385909a52b2dc65417a9cb4ccb856170118df02a3",
    ),
    "native_host_audit": (
        ROOT / "results/ag_field_funnel_bounded_frame_native_host_closure_independent_hostile_audit.md",
        "0b05127f62a89b16060d9cbaacde93009a56ac056e08bb7e16e0b8f79322c3b1",
    ),
}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def field_branch(dimension: int, compact: bool) -> str:
    """Return the theorem component covering one n >= 2 field input."""
    if dimension < 2:
        raise ValueError("the global theorem is scoped to hyperbolic n >= 2")
    if not compact:
        return "rational_cusp"
    if dimension == 2:
        return "real_dm_orbicurve"
    if dimension == 3:
        return "compact_h3_general_fork"
    if dimension in (4, 5):
        return "codimension_descent"
    return "minimality14_n_ge_6"


def determinant3(matrix: tuple[tuple[Fraction, ...], ...]) -> Fraction:
    a, b, c = matrix[0]
    d, e, f = matrix[1]
    g, h, i = matrix[2]
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def fork_identity(u: Fraction, a: Fraction, b: Fraction) -> dict[str, object]:
    gram = (
        (Fraction(-2), u, a),
        (u, Fraction(-2), b),
        (a, b, Fraction(-2)),
    )
    beta_from_determinant = determinant3(gram) / 2
    # With c_ij=-G_ij in the displayed negative-square convention,
    # z=c12*c23*c31=-uab.
    t12, t13, t23, z = u * u, a * a, b * b, -u * a * b
    beta_from_cycles = t12 + t13 + t23 - z - 4
    return {
        "determinant_half": str(beta_from_determinant),
        "cycle_coordinate": str(beta_from_cycles),
        "equal": beta_from_determinant == beta_from_cycles,
    }


def local_markdown_links(text: str) -> list[str]:
    links = re.findall(r"(?<!!)\[[^\]]+\]\(([^)]+)\)", text)
    return [
        link.split("#", 1)[0]
        for link in links
        if link and not link.startswith(("http://", "https://", "#"))
    ]


def build_certificate() -> dict[str, object]:
    text = OWNER.read_text(encoding="utf-8")
    actual_hashes = {name: sha256_file(path) for name, (path, _) in FROZEN_INPUTS.items()}
    expected_hashes = {name: expected for name, (_, expected) in FROZEN_INPUTS.items()}
    pins_pass = actual_hashes == expected_hashes

    required_fragments = {
        "all_chamber_source_groupoid": r"\mathsf{Cox}^{\mathrm{mark}}_{\mathrm{arith}}",
        "displayed_base_field": "defined over a displayed totally real base field",
        "target_morphisms": "Morphisms are field-isomorphism-compatible projective isomorphisms",
        "exact_return": r"R_{\mathrm{arr}}F_{\mathrm{arr}}(P)\cong P",
        "structural_vinberg_reduction": "Vinberg's structural reflection lemma",
        "global_field_scope": "In hyperbolic dimensions\n> \\(n\\ge2\\)",
        "h3_primitive_coordinate": r"k=\mathbf Q(\beta)",
        "noncompact_rational": r"Noncompact chambers have \(k=\mathbf Q\)",
        "admissible_native_hosts": "integral admissible arithmetic hosts",
        "dimension_ceiling": r"n\le995",
        "rank_range": r"\{3,\ldots,996\}",
        "arithmetic_maximality": "The same conclusion holds for maximality among arithmetic",
        "ag_hybrid_attribution": "AG-mediated**, not AG-only",
        "direct_ag_degree_open": "is still not known and is not claimed",
        "fixed_data_nonimport": "does not invoke that theorem's finiteness",
    }
    fragment_checks = {name: fragment in text for name, fragment in required_fragments.items()}

    branch_checks = {
        "noncompact_all_dimensions": all(
            field_branch(n, compact=False) == "rational_cusp" for n in range(2, 1201)
        ),
        "compact_n2": field_branch(2, True) == "real_dm_orbicurve",
        "compact_n3": field_branch(3, True) == "compact_h3_general_fork",
        "compact_n4_n5": all(
            field_branch(n, True) == "codimension_descent" for n in (4, 5)
        ),
        "compact_n_ge_6": all(
            field_branch(n, True) == "minimality14_n_ge_6" for n in range(6, 1201)
        ),
    }

    links = local_markdown_links(text)
    missing_links = sorted(
        link for link in links if not (OWNER.parent / link).resolve().exists()
    )

    tex_checks = {
        "display_delimiters_balanced": text.count(r"\[") == text.count(r"\]"),
        "inline_delimiters_balanced": text.count(r"\(") == text.count(r"\)"),
    }
    dimension_quantifier = {
        "hyperbolic_dimensions": [2, 995],
        "native_ranks": [3, 996],
        "rank_count": 996 - 3 + 1,
    }
    fork_checks = [
        fork_identity(Fraction(5, 2), Fraction(1), Fraction(3, 2)),
        fork_identity(Fraction(3), Fraction(4, 3), Fraction(5, 4)),
    ]

    all_checks_pass = (
        pins_pass
        and all(fragment_checks.values())
        and all(branch_checks.values())
        and not missing_links
        and all(tex_checks.values())
        and dimension_quantifier["rank_count"] == 994
        and all(check["equal"] for check in fork_checks)
    )
    return {
        "schema": "agcat-global-maximal-arithmetic-coxeter-independent-audit-v1",
        "all_checks_pass": all_checks_pass,
        "verdict": "PASS" if all_checks_pass else "FAIL",
        "pins": {
            "passed": pins_pass,
            "actual": actual_hashes,
            "expected": expected_hashes,
        },
        "categorical_and_scope_interface": fragment_checks,
        "field_branch_partition": branch_checks,
        "dimension_quantifier": dimension_quantifier,
        "fork_cycle_identity_samples": fork_checks,
        "local_links": {"checked": links, "missing": missing_links},
        "tex": tex_checks,
    }


def main() -> None:
    certificate = build_certificate()
    print(json.dumps(certificate, indent=2, sort_keys=True))
    if not certificate["all_checks_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
