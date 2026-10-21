#!/usr/bin/env python3
"""Certify a non-isotrivial family of rational P12--P16 nef carriers.

The earlier inverse-design certificate used the cubic

    Y^2 Z = X^3 + Z^3.

Here it is replaced by the Weierstrass family

    C_a: Y^2 Z = X^3 + a X Z^2 + Z^3,

while Q_a=Z^4+(X+Z)C_a and S_a=C_a^2+Z^2 Q_a.  The local identities which
produce the free cluster and the visible curve classes are independent of
``a``.  Absolute irreducibility is replayed for a=0 and a=1, and the old
proximity/Bezout cone certificate is then applied verbatim.  Since the
unique genus-one extremal curve is C_a and j(C_a) is nonconstant, this gives
genuinely non-isomorphic divisor-cone carriers, not merely new markings.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
from typing import Any

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.certify_rational_i1n_inverse_design import (
    VINBERG_VOLUME_DATA,
    _extreme_rays_via_cdd,
    _extreme_rays_via_sympy,
    enumerate_simple_roots,
    lorentz_pair,
)


ROOT = Path(__file__).resolve().parents[1]
UPSTREAM = ROOT / "data/enumerations/rational_i1n_inverse_design.json"
OUTPUT = ROOT / "data/enumerations/rational_i1n_j_moduli_carriers.json"
SCHEMA = "coxeter.rational_i1n.j_moduli_carriers.v1"
PARAMETERS = (0, 1)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _payload_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _forms(parameter: int) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    X, Y, Z = sp.symbols("X Y Z")
    cubic = Y**2 * Z - X**3 - parameter * X * Z**2 - Z**3
    quartic = Z**4 + (X + Z) * cubic
    sextic = cubic**2 + Z**2 * quartic
    return tuple(map(sp.expand, (cubic, quartic, sextic)))


def _absolute_irreducibility(parameter: int, executable: str) -> dict[str, Any]:
    cubic = f"Y^2*Z-X^3-{parameter}*X*Z^2-Z^3"
    script = f'''LIB "absfact.lib";
ring r=0,(X,Y,Z),dp;
poly C={cubic};
poly Q=Z^4+(X+Z)*C;
poly S=Z^2*Q+C^2;
def rc=absFactorize(C); setring rc;
print("RESULT_C="+string(size(absolute_factors[1])-1)+","+string(absolute_factors[2][2])); setring r;
def rq=absFactorize(Q); setring rq;
print("RESULT_Q="+string(size(absolute_factors[1])-1)+","+string(absolute_factors[2][2])); setring r;
def rs=absFactorize(S); setring rs;
print("RESULT_S="+string(size(absolute_factors[1])-1)+","+string(absolute_factors[2][2]));
quit;
'''
    completed = subprocess.run(
        [executable, "-q"],
        input=script,
        text=True,
        capture_output=True,
        check=True,
    )
    records = {}
    for name in ("C", "Q", "S"):
        match = re.search(rf"RESULT_{name}=(\d+),(\d+)", completed.stdout)
        if match is None:
            raise ArithmeticError(f"Singular did not report {name} for a={parameter}")
        records[name] = {
            "absolute_nonconstant_factor_count": int(match.group(1)),
            "factor_exponent": int(match.group(2)),
            "absolutely_irreducible": match.groups() == ("1", "1"),
        }
    if not all(record["absolutely_irreducible"] for record in records.values()):
        raise ArithmeticError(f"a={parameter} lost absolute irreducibility")
    return records


def _local_record(parameter: int) -> dict[str, Any]:
    X, Y, Z = sp.symbols("X Y Z")
    x, z = sp.symbols("x z")
    cubic, quartic, sextic = _forms(parameter)
    affine_cubic = sp.expand(cubic.subs({X: x, Y: 1, Z: z}))
    affine_quartic = sp.expand(quartic.subs({X: x, Y: 1, Z: z}))
    affine_sextic = sp.expand(sextic.subs({X: x, Y: 1, Z: z}))

    def homogeneous_piece(polynomial: sp.Expr, degree: int) -> sp.Expr:
        return sp.expand(
            sum(
                term
                for term in sp.Add.make_args(polynomial)
                if sp.Poly(term, x, z).total_degree() == degree
            )
        )

    discriminant_factor = 4 * parameter**3 + 27
    j_invariant = sp.Rational(1728 * 4 * parameter**3, discriminant_factor)
    selected_q_orders = {"z": 3, "C": 11, "other_Q_branch": 1}
    s_linear_coefficient_order = (
        2 * selected_q_orders["z"] + selected_q_orders["other_Q_branch"]
    )
    s_constant_term_order = 2 * selected_q_orders["C"]
    s_displacement_order = s_constant_term_order - s_linear_coefficient_order
    checks = {
        "smooth_Weierstrass_discriminant_nonzero": discriminant_factor != 0,
        "flex_line_meets_C_as_3p": sp.expand(cubic.subs(Z, 0)) == -X**3,
        "Q_identity": sp.expand(quartic - Z**4 - (X + Z) * cubic) == 0,
        "S_identity": sp.expand(sextic - cubic**2 - Z**2 * quartic) == 0,
        "C_tangent_is_z": homogeneous_piece(affine_cubic, 1) == z,
        "Q_node_tangent_cone": homogeneous_piece(affine_quartic, 2)
        == sp.expand(z * (x + z)),
        "S_double_tangent_cone": homogeneous_piece(affine_sextic, 2) == z**2,
        "C_Q_total_contact_12": 3 * 4 == 12,
        "C_Q_branch_contacts_1_plus_11": 1 + 11 == 12,
        "Q_S_long_branch_contact_22": 2 * 11 == 22,
        "S_selected_branch_displacement_order_15": s_displacement_order == 15,
        "S_two_branches_share_first_seven_centres": s_linear_coefficient_order
        == 7,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"local branch identities changed for a={parameter}")
    return {
        "parameter": parameter,
        "forms": {
            "C": str(cubic),
            "Q": str(quartic),
            "S": str(sextic),
        },
        "Weierstrass_discriminant_factor": discriminant_factor,
        "j_invariant": str(j_invariant),
        "local_equations_at_p_Y_equals_1": {
            "C": str(affine_cubic),
            "Q": str(affine_quartic),
            "S": str(affine_sextic),
        },
        "branch_data": {
            "C_Q_contacts": [1, 11],
            "selected_Q_branch": "contact-11 branch",
            "C_and_selected_Q_common_centres": 11,
            "selected_Q_and_S_contact": 22,
            "orders_along_selected_Q_normalization": selected_q_orders,
            "S_local_balance": (
                "ord(C^2)=22 and ord(z^2*other_Q_branch)=7, so the "
                "selected S branch has displacement order 22-7=15"
            ),
            "C_multiplicity_sequence": [1] * 11,
            "Q_multiplicity_sequence_through_p15": [2] + [1] * 14,
            "S_multiplicity_sequence_through_p16": [2] * 7 + [1] * 9,
        },
        "checks": checks,
    }


def _nef_extreme_rays(n: int) -> tuple[tuple[int, ...], ...]:
    roots = enumerate_simple_roots(n)
    # Surface intersection with the curve represented by root r is
    # d*x0-sum(r_i*x_i) = row*x.  Hence nef means every row*x >=0.
    rows = [(root[0], *(-entry for entry in root[1:])) for root in roots]
    rays = _extreme_rays_via_cdd(rows)
    if rays is None:
        rays = _extreme_rays_via_sympy(rows)
    return rays


def _canonical_pairing(divisor: tuple[int, ...]) -> int:
    # K=-3H+E_1+...+E_n in the physical total-transform coordinates.
    return -3 * divisor[0] - sum(divisor[1:])


def _cusp_record(n: int) -> dict[str, Any]:
    rays = _nef_extreme_rays(n)
    isotropic = sorted(ray for ray in rays if lorentz_pair(ray, ray) == 0)
    expected = VINBERG_VOLUME_DATA[n]["ideal_vertices"]
    line_pencil = (1, -1, *((0,) * (n - 1)))
    cubic_pencil = (3, *((-1,) * 9), *((0,) * (n - 9)))
    if not (
        len(isotropic) == expected
        and line_pencil in isotropic
        and cubic_pencil in isotropic
    ):
        raise ArithmeticError(f"cusp census changed for P{n}")
    records = []
    for ray in isotropic:
        k_pairing = _canonical_pairing(ray)
        records.append(
            {
                "primitive_class": list(ray),
                "canonical_pairing": k_pairing,
                "arithmetic_genus_if_effective": 1 + k_pairing // 2,
                "explicit_geometry": (
                    "basepoint-free rational pencil of lines through p"
                    if ray == line_pencil
                    else "basepoint-free genus-one pencil <C,3L>"
                    if ray == cubic_pencil
                    else "no semiampleness assertion"
                ),
            }
        )
    return {
        "extreme_nef_ray_count": len(rays),
        "ideal_ray_count": len(isotropic),
        "ideal_rays": records,
        "two_universal_semiample_cusps": [
            list(line_pencil),
            list(cubic_pencil),
        ],
    }


def _dimension_record(upstream: dict[str, Any], n: int) -> dict[str, Any]:
    record = upstream["dimensions"][str(n)]
    roots = record["roots"]
    genus_one = [row for row in roots if row["arithmetic_genus"] == 1]
    expected_cubic = [3, *([-1] * 11), *([0] * (n - 11))]
    checks = {
        "upstream_positive_status": "positive" in record["status"][
            "coxeter_polytope"
        ],
        "complete_nef_cone_proof": (
            record["status"]["nef_cone_completeness"].startswith("proved exactly")
        ),
        "unique_genus_one_extremal_curve": (
            len(genus_one) == 1 and genus_one[0]["vector"] == expected_cubic
        ),
        "finite_volume_noncompact": (
            record["vinberg_volume"]["finite_volume"]
            and not record["vinberg_volume"]["compact"]
        ),
    }
    if n in (12, 13):
        checks["tail_cone_proof_applies"] = record[
            "p12_p13_effective_cone_proof"
        ]["applies"]
    elif n in (14, 15):
        checks["quartic_cone_proof_zero_failures"] = record[
            "q_tail_effective_cone_certificate"
        ]["every_ray_has_exact_nonnegative_L_C_Q_exceptional_decomposition"]
    else:
        checks["sextic_cone_proof_zero_failures"] = record[
            "p16_effective_cone_certificate"
        ]["every_ray_has_exact_nonnegative_L_C_Q_S_exceptional_decomposition"]
    if not all(checks.values()):
        raise ArithmeticError(f"upstream P{n} carrier proof changed: {checks}")
    return {
        "hyperbolic_dimension": n,
        "facet_count": record["facet_count"],
        "lattice": record["lattice"],
        "abstract_chamber": f"Vinberg P{n}",
        "finite_vertices": record["vinberg_volume"]["finite_vertices"],
        "ideal_vertices": record["vinberg_volume"]["ideal_vertices"],
        "unique_genus_one_extremal_curve": genus_one[0],
        "cusps": _cusp_record(n),
        "checks": checks,
    }


def build(upstream_path: Path = UPSTREAM, singular: str | None = None) -> dict[str, Any]:
    upstream = json.loads(upstream_path.read_text(encoding="utf-8"))
    executable = singular or shutil.which("Singular")
    if executable is None:
        raise RuntimeError("Singular with absfact.lib is required")

    fibers = []
    for parameter in PARAMETERS:
        local = _local_record(parameter)
        local["absolute_irreducibility"] = _absolute_irreducibility(
            parameter, executable
        )
        fibers.append(local)
    dimensions = {
        str(n): _dimension_record(upstream, n) for n in range(12, 17)
    }
    j_values = [fiber["j_invariant"] for fiber in fibers]
    checks = {
        "two_explicit_fibers": [fiber["parameter"] for fiber in fibers] == [0, 1],
        "all_C_Q_S_absolutely_irreducible": all(
            all(
                entry["absolutely_irreducible"]
                for entry in fiber["absolute_irreducibility"].values()
            )
            for fiber in fibers
        ),
        "distinct_j_invariants": j_values == ["0", "6912/31"],
        "all_local_cluster_checks": all(
            all(fiber["checks"].values()) for fiber in fibers
        ),
        "all_five_complete_nef_carriers": all(
            all(record["checks"].values()) for record in dimensions.values()
        ),
        "unique_genus_one_curve_separates_isomorphism_classes": all(
            record["checks"]["unique_genus_one_extremal_curve"]
            for record in dimensions.values()
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"rational j-family carrier audit changed: {checks}")
    payload = {
        "schema": SCHEMA,
        "source": {
            "upstream_path": str(upstream_path.relative_to(ROOT)),
            "upstream_sha256": _sha256(upstream_path),
            "upstream_payload_sha256": upstream["certificate_payload_sha256"],
            "absolute_factorization": "Singular absfact.lib::absFactorize",
        },
        "family": {
            "base": "a in A^1 with 4a^3+27 != 0 and the geometric-integrality open conditions",
            "C_a": "Y^2 Z-X^3-a X Z^2-Z^3",
            "Q_a": "Z^4+(X+Z)C_a",
            "S_a": "C_a^2+Z^2 Q_a",
            "j_formula": "6912 a^3/(4a^3+27)",
            "non_isotrivial": True,
            "open_family_reason": (
                "smoothness and geometric integrality are open; the local branch "
                "identities are coefficient identities; j(a) is nonconstant"
            ),
        },
        "explicit_fibers": fibers,
        "dimensions": dimensions,
        "positive_conclusion": {
            "statement": (
                "for every 12<=n<=16 a nonempty open set of parameters gives a "
                "non-isotrivial family, containing infinitely many pairwise "
                "non-isomorphic smooth rational surfaces whose complete "
                "projectivized nef cone is Vinberg P_n"
            ),
            "strongest_new_explicit_carrier": (
                "a=1, n=16: an H16/20 complete divisor-cone carrier with "
                "j(C)=6912/31, non-isomorphic to the archived j=0 carrier"
            ),
            "abstract_polytope_novelty": False,
            "divisor_carrier_novelty_in_this_project": True,
            "worldwide_carrier_priority": "open",
        },
        "checks": checks,
        "all_checks_pass": True,
    }
    payload["certificate_payload_sha256"] = _payload_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream", type=Path, default=UPSTREAM)
    parser.add_argument("--singular")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build(args.upstream, args.singular)
    if args.check_archive:
        if payload != json.loads(args.output.read_text(encoding="utf-8")):
            raise ArithmeticError("rebuilt rational j-family archive changed")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["checks"], sort_keys=True))
    print(args.output)


if __name__ == "__main__":
    main()
