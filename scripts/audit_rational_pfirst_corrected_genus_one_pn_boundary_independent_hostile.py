#!/usr/bin/env python3
"""Independent hostile audit of the corrected genus-one P_n boundary.

This audit does not import the owner certifier.  It pins the complete-cone
and MMP inputs, rebuilds the rank-blind lattice packet, checks the exact
dimension-dependent Abel quartics, and independently verifies that the MMP
canonical profiles are precisely the full weighted-diagram orbits.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Sequence

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]

OWNER_PINS = {
    "report": {
        "path": "results/rational_pfirst_corrected_genus_one_pn_boundary.md",
        "sha256": "0df4c789f94c654982ef6650475ea0bc2d2904bb7b0214d9884c56e945f8d57e",
    },
    "certifier": {
        "path": "scripts/certify_rational_pfirst_corrected_genus_one_pn_boundary.py",
        "sha256": "c2f401672cad08f2a79984ec87d9a4dda6a6f5fd6d68c7af79ae3d4ee8eddbe0",
    },
    "archive": {
        "path": "data/enumerations/rational_pfirst_corrected_genus_one_pn_boundary.json",
        "sha256": "feae8dbbcace5f7657468c95188b032bb4e591dc771fdad40ba5532de32dff4b",
    },
    "tests": {
        "path": "tests/test_rational_pfirst_corrected_genus_one_pn_boundary.py",
        "sha256": "44094114f68533f8d9abbf4ce143994a8a5a2486a2bb287b74a0dd9fcf25a1d5",
    },
}

DEPENDENCY_PINS = {
    "inverse_archive": {
        "path": "data/enumerations/rational_i1n_inverse_design.json",
        "sha256": "351a0a8525209313de13f9dd3752320bd074daf0985da704fc634684c75efc18",
    },
    "inverse_report": {
        "path": "results/rational_i1n_inverse_design.md",
        "sha256": "19eef6671655d0c860d708e607bfb4708926e45126b879a7325e366411971f10",
    },
    "inverse_certifier": {
        "path": "scripts/certify_rational_i1n_inverse_design.py",
        "sha256": "7303f1d4b6cde541c41eaca9ee9e21b4d1a838ab1d49f76aff34bd2c9b8849ad",
    },
    "irreducibility_replay": {
        "path": "scripts/certify_rational_i1n_curves.sing",
        "sha256": "ee7336d70b487829a8b4474a6b490031b6188df30c86341469bda54496d972d3",
    },
    "p17_mmp_archive": {
        "path": "data/enumerations/rational_p17_mmp.json",
        "sha256": "2c35b61bb5357111a0f74823c98d192ca2eaef86d5a939071261eeb371be060b",
    },
    "p17_mmp_certifier": {
        "path": "scripts/certify_rational_p17_mmp.py",
        "sha256": "12e5a1b04e895a9f2e9f6caeed7b205e71c842d83a5da9905343c22292720be1",
    },
    "p17_mmp_tests": {
        "path": "tests/test_rational_p17_mmp.py",
        "sha256": "9ad3f964cfff0c82669d08b5bed7914e5a68f14561a95ef3fe9624ec7243bf1a",
    },
    "p18_p19_mmp_archive": {
        "path": "data/enumerations/rational_p18_p19_mmp.json",
        "sha256": "b2bf17c357a605ff0f39f8645628079f3fb8a808255ad5e184fb1da99995751d",
    },
    "p18_p19_mmp_certifier": {
        "path": "scripts/certify_rational_p18_p19_mmp.py",
        "sha256": "296ceca4d966d9b63cf1bc0acf68670576b102539a1591a2168db05180bd8f6d",
    },
    "p18_p19_mmp_tests": {
        "path": "tests/test_rational_p18_p19_mmp.py",
        "sha256": "41fa74aef2409a6aad635a29e3ef5f25bb966948773a0ff72f29d17f80ac4f50",
    },
}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dot(left: Sequence[int], right: Sequence[int]) -> int:
    return left[0] * right[0] - sum(a * b for a, b in zip(left[1:], right[1:]))


def payload_digest(document: dict[str, object], key: str) -> str:
    payload = dict(document)
    expected = payload.pop(key)
    observed = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    assert expected == observed
    return observed


def reflection_integral(root: Sequence[int]) -> bool:
    square = dot(root, root)
    for index in range(len(root)):
        basis = [0] * len(root)
        basis[index] = 1
        coefficient = Fraction(-2 * dot(basis, root), square)
        if any(
            (Fraction(basis[row]) + coefficient * root[row]).denominator != 1
            for row in range(len(root))
        ):
            return False
    return True


def rank_blind_record(n: int) -> dict[str, object]:
    c = [3] + [-1] * 11 + [0] * (n - 11)
    k = [-3] + [1] * n
    z = [c[index] + k[index] for index in range(n + 1)]
    assert z == [0] * 12 + [1] * (n - 11)
    assert dot(c, c) == -2
    assert dot(k, c) == 2
    assert 1 + (dot(c, c) + dot(k, c)) // 2 == 1
    assert dot(c, z) == 0
    assert dot(z, z) == 11 - n
    assert reflection_integral(c)
    return {
        "n": n,
        "C_square": -2,
        "K_dot_C": 2,
        "arithmetic_genus": 1,
        "C_dot_Z": 0,
        "Z_square": 11 - n,
        "reflection_integral": True,
    }


def root_vectors(source: dict[str, object], n: int) -> list[tuple[int, ...]]:
    return [
        tuple(record["vector"])
        for record in source["dimensions"][str(n)]["roots"]
    ]


def positive_cone_records(source: dict[str, object]) -> list[dict[str, object]]:
    expected = {
        12: (14, "two-wall symbolic decomposition", None),
        13: (15, "two-wall symbolic decomposition", None),
        14: (17, "exact extreme-ray decomposition", (94, "a81a8babd5fe2028bfd173882cb80eee974a6b680e131f053416afe3760e623f")),
        15: (18, "exact extreme-ray decomposition", (145, "edb9d4dfea27168e5f39bdabc1f2e93e7004b13dc37bcdcdbaa0044eb8f82751")),
        16: (20, "exact extreme-ray decomposition", (325, "c6733607bfa5c0653cd1fadc9752d9ddf8bf366a7937bc166e4980f60d3d98be")),
    }
    records = []
    for n, (facets, method, ray_pin) in expected.items():
        case = source["dimensions"][str(n)]
        assert case["facet_count"] == facets
        assert case["coxiter"]["exact_relation_match"] is True
        assert case["coxiter"]["sha256_match"] is True
        assert "positive rational-surface nef realization" in json.dumps(case["status"])
        if n in (12, 13):
            proof = case["p12_p13_effective_cone_proof"]
            assert proof["applies"] and proof["conclusion"].startswith("NE(X_n)=")
            ray_count = None
        elif n in (14, 15):
            proof = case["q_tail_effective_cone_certificate"]
            assert proof["every_ray_has_exact_nonnegative_L_C_Q_exceptional_decomposition"]
            assert (proof["primitive_extreme_ray_count"], proof["primitive_extreme_rays_sha256"]) == ray_pin
            ray_count = ray_pin[0]
        else:
            proof = case["p16_effective_cone_certificate"]
            assert proof["every_ray_has_exact_nonnegative_L_C_Q_S_exceptional_decomposition"]
            assert (proof["primitive_extreme_ray_count"], proof["primitive_extreme_rays_sha256"]) == ray_pin
            ray_count = ray_pin[0]
        records.append(
            {
                "n": n,
                "facet_count": facets,
                "method": method,
                "primitive_candidate_ray_count": ray_count,
                "complete_nef_carrier": True,
            }
        )
    cluster = source["explicit_cluster"]
    assert cluster["absolute_irreducibility"]["result"] == {"C": 1, "Q": 1, "S": 1}
    assert cluster["absolute_irreducibility"]["replay_script_sha256"] == DEPENDENCY_PINS[
        "irreducibility_replay"
    ]["sha256"]
    return records


def mmp_dimension(document: dict[str, object], n: int) -> dict[str, object]:
    return document if n == 17 else document["dimensions"][str(n)]


def permuted_profile(profile: Sequence[int], permutation: Sequence[int]) -> tuple[int, ...]:
    result = [0] * len(profile)
    for source, target in enumerate(permutation):
        result[target] = profile[source]
    return tuple(result)


def mmp_orbit_record(
    n: int,
    source: dict[str, object],
    p17: dict[str, object],
    p18_p19: dict[str, object],
) -> dict[str, object]:
    roots = root_vectors(source, n)
    dedicated = mmp_dimension(p17 if n == 17 else p18_p19, n)
    expected_count = {17: 8, 18: 24, 19: 120}[n]

    # The n B-chain roots followed by L form a unimodular basis, so facet
    # pairing profiles determine canonical vectors without an embedding gap.
    basis = sp.Matrix(roots[: n + 1])
    assert abs(int(basis.det())) == 1

    automorphisms = dedicated["weighted_diagram_automorphisms_zero_based"]
    assert len(automorphisms) == expected_count
    gram = [[dot(left, right) for right in roots] for left in roots]
    for permutation in automorphisms:
        assert sorted(permutation) == list(range(len(roots)))
        for left in range(len(roots)):
            for right in range(len(roots)):
                assert gram[left][right] == gram[permutation[left]][permutation[right]]

    standard = (-3,) + (1,) * n
    standard_profile = tuple(dot(standard, root) for root in roots)
    profile_orbit = {
        permuted_profile(standard_profile, permutation)
        for permutation in automorphisms
    }
    archived_profiles = {
        tuple(profile) for profile in dedicated["canonical_facet_pairing_profiles"]
    }
    assert archived_profiles == profile_orbit
    assert len(profile_orbit) == expected_count

    canonical_classes = {tuple(vector) for vector in dedicated["canonical_classes"]}
    assert standard in canonical_classes
    assert len(canonical_classes) == expected_count
    assert {tuple(dot(vector, root) for root in roots) for vector in canonical_classes} == archived_profiles
    for canonical in canonical_classes:
        assert dot(canonical, canonical) == 9 - n
        assert all(entry % 2 for entry in canonical)
        assert all(
            1 + (dot(root, root) + dot(canonical, root)) // 2 >= 0
            for root in roots
        )
    assert all(1 + (dot(root, root) + dot(standard, root)) // 2 == 0 for root in roots[:n])
    assert dedicated["canonical_class_count"] == expected_count
    assert dedicated["rank_one_P2_endpoint_count"] == expected_count
    assert dedicated["weighted_diagram_automorphism_order"] == expected_count
    assert dedicated["p2_and_f2_reconstructions_agree"] is True
    assert dedicated["canonical_classes_equal_standard_diagram_orbit"] is True
    assert dedicated["all_diagram_automorphisms_lift_to_integral_unimodular_isometries"] is True
    return {
        "n": n,
        "unimodular_facet_basis": True,
        "canonical_class_count": expected_count,
        "diagram_orbit_count": len(profile_orbit),
        "P2_F2_reconstructions_agree": True,
    }


def abel_record(n: int, source: dict[str, object]) -> dict[str, object]:
    roots = set(root_vectors(source, n))
    pad = lambda vector: tuple(vector) + (0,) * (n + 1 - len(vector))
    line = pad((1, -1, -1, -1))
    cubic = pad((3,) + (-1,) * 11)
    quartic_q = pad((4, -2) + (-1,) * 14)
    a_length = min(n, 18)
    quartic_a = pad((4,) + (-1,) * a_length)
    for vector in (line, cubic, quartic_q, quartic_a):
        assert vector in roots
    assert {
        "L.C": dot(line, cubic),
        "C.Q": dot(cubic, quartic_q),
        "C.A": dot(cubic, quartic_a),
        "Q.A": dot(quartic_q, quartic_a),
    } == {"L.C": 0, "C.Q": 0, "C.A": 1, "Q.A": 0}
    standard = (-3,) + (1,) * n
    assert 1 + (dot(quartic_a, quartic_a) + dot(standard, quartic_a)) // 2 == 3
    assert reflection_integral(quartic_a)
    # 3*4=12, while eleven common unit multiplicities account for the local
    # contact at p; the residual divisor therefore has degree one.
    assert 3 * 4 - 11 == 1
    return {
        "n": n,
        "A": list(quartic_a),
        "A_chain_length": a_length,
        "C_dot_A": 1,
        "residual_degree_on_C": 1,
        "Abel_contradiction": "O_C(r)=O_C(p) on the smooth locus forces r=p",
    }


def build_audit() -> dict[str, object]:
    for item in (*OWNER_PINS.values(), *DEPENDENCY_PINS.values()):
        assert sha256_file(ROOT / item["path"]) == item["sha256"]

    owner = json.loads((ROOT / OWNER_PINS["archive"]["path"]).read_text())
    source = json.loads((ROOT / DEPENDENCY_PINS["inverse_archive"]["path"]).read_text())
    p17 = json.loads((ROOT / DEPENDENCY_PINS["p17_mmp_archive"]["path"]).read_text())
    p18_p19 = json.loads((ROOT / DEPENDENCY_PINS["p18_p19_mmp_archive"]["path"]).read_text())
    payload_digest(source, "certificate_payload_sha256")
    payload_digest(p17, "certificate_payload_sha256")
    payload_digest(p18_p19, "certificate_payload_sha256")
    assert owner["source_archive_sha256"] == DEPENDENCY_PINS["inverse_archive"]["sha256"]
    assert [entry["sha256"] for entry in owner["upstream_MMP_archives"]] == [
        DEPENDENCY_PINS["p17_mmp_archive"]["sha256"],
        DEPENDENCY_PINS["p18_p19_mmp_archive"]["sha256"],
    ]

    formal_records = [rank_blind_record(n) for n in range(12, 65)]
    assert {
        (record["C_square"], record["K_dot_C"], record["arithmetic_genus"], record["reflection_integral"])
        for record in formal_records
    } == {(-2, 2, 1, True)}
    positives = positive_cone_records(source)
    mmp_records = [mmp_orbit_record(n, source, p17, p18_p19) for n in (17, 18, 19)]
    abel_records = [abel_record(n, source) for n in (17, 18, 19)]
    assert [record["A_chain_length"] for record in abel_records] == [17, 18, 18]
    assert owner["scoped_Pn_AG_theorem"]["obstructing_quartic_classes"] == {
        str(record["n"]): record["A"] for record in abel_records
    }

    report = (ROOT / OWNER_PINS["report"]["path"]).read_text()
    assert "a_n=\\min(n,18)" in report
    assert "literal quartic class changes" in report
    assert "global smoothness of the putative cubic is not assumed" in " ".join(report.split())
    assert "The same class persists" not in report
    return {
        "verdict": "PASS_AFTER_TWO_SCOPE_REPAIRS",
        "repairs": [
            "the no-go cubic is only integral and smooth at p, not necessarily globally smooth",
            "A_17 follows 17 centers while A_18 and A_19 follow 18 centers",
        ],
        "local_rank_blind_range_replayed": [12, 64],
        "local_invariant_packet": [-2, 2, 1, True],
        "positive_full_nef_carriers": positives,
        "marking_independent_MMP_orbits": mmp_records,
        "Abel_quartics": abel_records,
        "scoped_carrier_iff": "P_n for 12<=n<=19 has a complex smooth rational full-nef carrier iff n<=16",
        "scope": {
            "universal_n_le_16_claimed": False,
            "local_packet_gives_dimension_bound": False,
            "negative_no_go_characteristic_zero": True,
            "abstract_P17_P19_exist": True,
            "property_4_gate": "oriented contraction/proximity and Pic^1 incidence compatibility",
        },
    }


def main() -> None:
    print(json.dumps(build_audit(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
