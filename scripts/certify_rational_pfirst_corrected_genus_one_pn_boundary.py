#!/usr/bin/env python3
"""Certify the rank-blind corrected genus-one packet and scoped P_n bound."""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/rational_pfirst_corrected_genus_one_pn_boundary.json"
SOURCE = ROOT / "data/enumerations/rational_i1n_inverse_design.json"
SOURCE_SHA256 = "351a0a8525209313de13f9dd3752320bd074daf0985da704fc634684c75efc18"
P17_MMP = ROOT / "data/enumerations/rational_p17_mmp.json"
P17_MMP_SHA256 = "2c35b61bb5357111a0f74823c98d192ca2eaef86d5a939071261eeb371be060b"
P18_P19_MMP = ROOT / "data/enumerations/rational_p18_p19_mmp.json"
P18_P19_MMP_SHA256 = "b2bf17c357a605ff0f39f8645628079f3fb8a808255ad5e184fb1da99995751d"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dot(x: Sequence[int], y: Sequence[int]) -> int:
    return x[0] * y[0] - sum(a * b for a, b in zip(x[1:], y[1:]))


def reflection_integral(root: Sequence[int]) -> bool:
    square = dot(root, root)
    for j in range(len(root)):
        e = [0] * len(root)
        e[j] = 1
        coefficient = Fraction(-2 * dot(e, root), square)
        if any((Fraction(e[i]) + coefficient * root[i]).denominator != 1 for i in range(len(root))):
            return False
    return True


def check_mmp_archives() -> tuple[dict[str, object], dict[str, object]]:
    """Pin the exhaustive marking-classification inputs used for n=17,18,19."""

    for path, expected in (
        (P17_MMP, P17_MMP_SHA256),
        (P18_P19_MMP, P18_P19_MMP_SHA256),
    ):
        if sha256(path) != expected:
            raise ArithmeticError(f"frozen surface-MMP archive changed: {path.name}")
    p17 = json.loads(P17_MMP.read_text())
    p18_p19 = json.loads(P18_P19_MMP.read_text())
    if p17["schema"] != "rational_p17_mmp_v1":
        raise ArithmeticError("P17 MMP schema changed")
    if p18_p19["schema"] != "rational_p18_p19_mmp_v1":
        raise ArithmeticError("P18/P19 MMP schema changed")
    return p17, p18_p19


def build() -> dict[str, object]:
    if sha256(SOURCE) != SOURCE_SHA256:
        raise ArithmeticError("frozen rational I1n inverse-design archive changed")
    source = json.loads(SOURCE.read_text())
    p17_mmp, p18_p19_mmp = check_mmp_archives()
    rows = []
    abel_quartics: dict[str, list[int]] = {}
    for n in range(12, 20):
        case = source["dimensions"][str(n)]
        genus_one = [root for root in case["roots"] if root.get("arithmetic_genus") == 1]
        if len(genus_one) != 1:
            raise ArithmeticError(f"P{n} no longer has one displayed genus-one root")
        root = genus_one[0]
        c = root["vector"]
        expected_c = [3] + [-1] * 11 + [0] * (n - 11)
        k = [-3] + [1] * n
        z = [c[i] + k[i] for i in range(n + 1)]
        expected_z = [0] * 12 + [1] * (n - 11)
        if c != expected_c or z != expected_z:
            raise ArithmeticError(f"P{n} corrected anticanonical decomposition changed")
        if not (
            dot(c, c) == -2
            and dot(k, c) == 2
            and dot(c, z) == 0
            and dot(z, z) == 11 - n
            and reflection_integral(c)
        ):
            raise ArithmeticError(f"P{n} genus-one lattice identities changed")
        status_text = json.dumps(case["status"], sort_keys=True)
        carrier = n <= 16
        if carrier and "positive rational-surface nef realization" not in status_text:
            raise ArithmeticError(f"P{n} positive carrier status changed")
        if not carrier and "impossible on every smooth projective rational surface" not in status_text:
            raise ArithmeticError(f"P{n} AG no-go status changed")
        obstruction = case["flex_cubic_obstruction"]
        if not carrier and not (
            obstruction["applies"]
            and obstruction["marking_independent_mmp_upgrade"] is not None
            and "Abel" in " ".join(obstruction["proof"])
        ):
            raise ArithmeticError(f"P{n} MMP/Abel obstruction changed")
        if not carrier:
            a_label = case["named_root_labels"]["A"]
            a = case["roots"][a_label - 1]["vector"]
            expected_a_length = min(n, 18)
            expected_a = [4] + [-1] * expected_a_length + [0] * (n - expected_a_length)
            k_dot_a = dot(k, a)
            a_genus = 1 + (dot(a, a) + k_dot_a) // 2
            if not (
                a == expected_a
                and dot(c, a) == 1
                and reflection_integral(a)
                and a_genus == 3
            ):
                raise ArithmeticError(f"P{n} Abel quartic class changed")
            abel_quartics[str(n)] = a
            dedicated = p17_mmp if n == 17 else p18_p19_mmp["dimensions"][str(n)]
            upgrade = obstruction["marking_independent_mmp_upgrade"]
            expected_count = {17: 8, 18: 24, 19: 120}[n]
            rank_two = dedicated["rank_two_endpoint_types"]
            expected_rank_two = (
                {"F1_nonminimal": expected_count, "F2_minimal": expected_count}
                if n == 17
                else {"F1": expected_count, "F2": expected_count}
            )
            if not (
                rank_two == expected_rank_two
                and dedicated["rank_one_P2_endpoint_count"] == expected_count
                and dedicated["canonical_class_count"] == expected_count
                and dedicated["weighted_diagram_automorphism_order"] == expected_count
                and dedicated["all_diagram_automorphisms_lift_to_integral_unimodular_isometries"]
                and dedicated["p2_and_f2_reconstructions_agree"]
                and dedicated["canonical_classes_equal_standard_diagram_orbit"]
                and upgrade["canonical_class_count"] == expected_count
                and upgrade["rank_one_P2_endpoints"] == expected_count
                and upgrade["weighted_diagram_automorphism_order"] == expected_count
            ):
                raise ArithmeticError(f"P{n} dedicated MMP classification changed")
        rows.append(
            {
                "n": n,
                "rho": n + 1,
                "facet_count": case["facet_count"],
                "C": c,
                "K": k,
                "Z_equals_C_plus_K": z,
                "C_square": dot(c, c),
                "K_dot_C": dot(k, c),
                "arithmetic_genus": root["arithmetic_genus"],
                "C_reflection_integral": True,
                "C_dot_Z": dot(c, z),
                "Z_square": dot(z, z),
                "complex_rational_full_nef_carrier": carrier,
                "status": case["status"],
                "marking_independent_MMP_Abel_obstruction": (
                    None if carrier else obstruction["marking_independent_mmp_upgrade"]
                ),
            }
        )

    invariant_tuples = {
        (row["C_square"], row["K_dot_C"], row["arithmetic_genus"], row["C_reflection_integral"])
        for row in rows
    }
    if invariant_tuples != {(-2, 2, 1, True)}:
        raise ArithmeticError("local genus-one packet ceased to be rank-blind")
    if [row["n"] for row in rows if row["complex_rational_full_nef_carrier"]] != [12, 13, 14, 15, 16]:
        raise ArithmeticError("positive corrected-genus-one range changed")
    if [row["n"] for row in rows if not row["complex_rational_full_nef_carrier"]] != [17, 18, 19]:
        raise ArithmeticError("negative corrected-genus-one range changed")

    return {
        "schema": "rational-pfirst-corrected-genus-one-pn-boundary-v1",
        "status": "PASS_RANK_BLIND_LOCAL_PACKET_AND_SCOPED_PN_N_LE_16",
        "all_checks_pass": True,
        "local_NS_theorem": {
            "formula": "C=-K+Z, with C^2=-2, K.C=2, p_a(C)=1, C.Z=0, Z^2=11-n",
            "rank_blind_data": ["C^2=-2", "K.C=2", "p_a(C)=1", "reflection in C is integral"],
            "dimensions_replayed": list(range(12, 20)),
            "conclusion": "reflection integrality plus adjunction of C alone gives no n-bound once C != -K",
        },
        "scoped_Pn_AG_theorem": {
            "class": "the exact marked classical P_n facet systems for 12<=n<=19",
            "complex_rational_full_nef_carrier_iff": "n<=16",
            "positive_dimensions": [12, 13, 14, 15, 16],
            "negative_dimensions": [17, 18, 19],
            "positive_method": "explicit cluster plus exact proximity/Bezout Mori-cone decomposition",
            "negative_method": "facet-contraction MMP recovers the free chain; Pic^1 Abel injectivity contradicts the required quartic wall",
            "obstructing_quartic_classes": abel_quartics,
            "not_from": ["bare reflection maximality", "Coxeter finite volume", "the local genus-one norm/adjunction packet"],
        },
        "rows": rows,
        "property_4_native_invariant": (
            "the oriented facet-contraction/proximity design together with a nonempty global curve-incidence "
            "certificate (or its Picard-group obstruction), not merely the NS class and reflection"
        ),
        "scope": {
            "base_field_for_n17_n19_no_go": "characteristic zero",
            "not_claimed": [
                "a universal n<=16 theorem for all reflection-maximal chambers with a genus-one facet",
                "classification of abstract P_n chambers beyond the declared finite family",
                "nonexistence of P17, P18, or P19 as abstract hyperbolic polytopes",
                "a bound from local NS/canonical data alone",
            ],
        },
        "upstream_MMP_archives": [
            {
                "path": str(P17_MMP.relative_to(ROOT)),
                "sha256": P17_MMP_SHA256,
                "certificate_payload_sha256": p17_mmp["certificate_payload_sha256"],
                "dimensions": [17],
            },
            {
                "path": str(P18_P19_MMP.relative_to(ROOT)),
                "sha256": P18_P19_MMP_SHA256,
                "certificate_payload_sha256": p18_p19_mmp["certificate_payload_sha256"],
                "dimensions": [18, 19],
            },
        ],
        "source_archive": str(SOURCE.relative_to(ROOT)),
        "source_archive_sha256": SOURCE_SHA256,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    payload = build()
    if args.check_archive:
        if payload != json.loads(ARCHIVE.read_text()):
            raise SystemExit("live corrected-genus-one certificate differs from archive")
        print("archived corrected-genus-one P_n certificate matches")
    elif args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(payload["status"])
        print("positive", payload["scoped_Pn_AG_theorem"]["positive_dimensions"])
        print("negative", payload["scoped_Pn_AG_theorem"]["negative_dimensions"])


if __name__ == "__main__":
    main()
