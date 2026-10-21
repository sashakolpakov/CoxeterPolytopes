#!/usr/bin/env python3
"""Independent hostile replay of the enriched AG/ABSW transfer.

The owner certifier is intentionally not imported.  This audit recomputes
the degree, generator, grid, conjugate, and discriminant constants and pins
both the branch theorems and their relevant hostile audits.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from math import floor, log2
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]

OWNER_PINS = {
    "report": {
        "path": "results/enriched_canonical_chamber_ag_absw_transfer.md",
        "sha256": "8f6956eb5fe0ccd54d427a2bef199be1f5d5ff4d3236e7551ab385b1aa4e03c4",
    },
    "certifier": {
        "path": "scripts/certify_enriched_canonical_chamber_ag_absw_transfer.py",
        "sha256": "237148b6d61ab9fea3d8758ea80d13202133453c4b2f057fc22dd0d82b592284",
    },
    "archive": {
        "path": "data/enumerations/enriched_canonical_chamber_ag_absw_transfer.json",
        "sha256": "a6bb38844f89b3f448f05634b997e08e2f74654b1b74cf3ca7bcf85e62dfd11f",
    },
    "tests": {
        "path": "tests/test_enriched_canonical_chamber_ag_absw_transfer.py",
        "sha256": "c11696b78de36bac4207a8aa190db1c3e33950b99ef76d89c6b1f6446dbbc65b",
    },
}

DEPENDENCY_PINS = {
    "canonical_transfer": {
        "path": "results/canonical_chamber_rigidity_absw_transfer_theorem.md",
        "sha256": "850d2b5832c3f7bd9fa85c39a6e84449c8f1fb5ede6c98de931f966218da3767",
    },
    "canonical_transfer_audit": {
        "path": "results/canonical_chamber_rigidity_absw_transfer_theorem_hostile_audit.md",
        "sha256": "808106d7d71a0df68e8ee50c56c7e3b7ae77ec451a71e8227bbb5f75af384fe5",
    },
    "direct_rational": {
        "path": "results/rational_direct_full_ns_n11_dimension_theorem.md",
        "sha256": "a5dbd6e2fbbd7cd55025f5af324de6b778a196cf4f72cf06903a1dce2e4ae27f",
    },
    "direct_rational_audit": {
        "path": "results/rational_direct_full_ns_n11_dimension_theorem_hostile_audit.md",
        "sha256": "9293c924bca806212ccac9a776ff1f477d9f9009c9ba1ce328dd81b1bc193193",
    },
    "weak_log": {
        "path": "results/pfirst_weak_log_complete_carrier_transfer.md",
        "sha256": "bff338b5e2f4dde2c918742a853656661c2ac074684938031ebbc95df8fd14b2",
    },
    "weak_log_audit": {
        "path": "results/pfirst_weak_log_complete_carrier_transfer_independent_hostile_audit.md",
        "sha256": "ce7cc977b2172418a7d6667c4a4a6e8ac8c9a3f3ec150d151ee91f0f7109cb3f",
    },
    "trace_descent": {
        "path": "results/equivariant_minimalization_trace_descent_theorem.md",
        "sha256": "ff8a89ae37a24ad2a10515f5acd931b15a280ed8f7fafb77e1db0c3d6b61ae5d",
    },
    "trace_descent_audit": {
        "path": "results/equivariant_minimalization_trace_descent_hostile_audit.md",
        "sha256": "fa76984c0fe8c82a7616fd3b49e2401ed7e2351b7e37be8db796a7bb00ec647c",
    },
    "ordinary_property4": {
        "path": "results/rational_pfirst_ordinary_dimension_bound_and_exact_carrier.md",
        "sha256": "8ea1ede6cc1717ede57dc902e459ceb831682f9f6cd6328c70fc7657452920af",
    },
    "ordinary_property4_independent_audit": {
        "path": "results/rational_pfirst_ordinary_dimension_bound_and_exact_carrier_independent_hostile_audit.md",
        "sha256": "83a43a4f497fb2dd0cbb2fa8875e037a0ec7aab68cd99083a34b4c29b616f007",
    },
    "chordless_height": {
        "path": "results/k3_chordless_cycle_height_field_bound.md",
        "sha256": "3372cddd560471830fc0a9cbaac02560ce4a682781db1cdc877b8cd108ea471c",
    },
    "chordless_height_independent_audit": {
        "path": "results/k3_chordless_cycle_height_field_bound_independent_hostile_audit.md",
        "sha256": "eb9e8b40ce646e7c7806eb157fbcd6a80393fd9d52a6bb2fb06061fb3fcd3914",
    },
    "rank4_escape": {
        "path": "results/rank4_k3_exact_maximal_eigen_interval_infinite_fields.md",
        "sha256": "9351875964fbc73da81b55ec3ff68db54af3dbe1ec90272555cff1938ab13c95",
    },
    "rank4_escape_audit": {
        "path": "results/rank4_k3_exact_maximal_eigen_interval_infinite_fields_hostile_scope_audit.md",
        "sha256": "3204a7bcfc580e24a0a015dd46f74e11d284a31ab33f42acdb6fa70c27220264",
    },
}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_digest(payload: Any) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def verify_pins() -> dict[str, str]:
    observed: dict[str, str] = {}
    for label, item in {**OWNER_PINS, **DEPENDENCY_PINS}.items():
        digest = sha256_file(ROOT / item["path"])
        assert digest == item["sha256"], (label, item["sha256"], digest)
        observed[label] = digest
    return observed


def field_constant_record(n: int) -> dict[str, int | str | bool]:
    assert 2 <= n <= 25
    trace_ceiling = 20 // (n + 1)
    if trace_ceiling < 2:
        return {
            "n": n,
            "full_trace_K3_degree_ceiling": trace_ceiling,
            "actual_field_degree_bound": 1,
            "field": "Q",
            "nonrational_possible": False,
        }
    degree = trace_ceiling
    length = n + 3
    selected = floor(log2(degree))
    grid = degree * (degree - 1) // 2 + 1
    small = 2 ** (length + 1)
    factor = selected * grid
    exponent = degree * (degree - 1)
    assert 2**selected <= degree < 2 ** (selected + 1)
    assert exponent == 2 * (degree * (degree - 1) // 2)
    return {
        "n": n,
        "full_trace_K3_degree_ceiling": degree,
        "actual_field_degree_bound": degree,
        "chordless_length_bound": length,
        "selected_generator_count": selected,
        "primitive_grid_side": grid,
        "small_coordinate_bound": small,
        "primitive_bound_factor": factor,
        "discriminant_exponent": exponent,
        "nonrational_possible": True,
    }


def uniform_constant_certificate() -> dict[str, object]:
    records = [field_constant_record(n) for n in range(2, 26)]
    nonrational = [record for record in records if record["nonrational_possible"]]
    assert [record["n"] for record in nonrational] == list(range(2, 10))
    assert max(int(record["primitive_bound_factor"]) for record in nonrational) == 32
    assert max(int(record["small_coordinate_bound"]) for record in nonrational) == 8192
    assert max(int(record["discriminant_exponent"]) for record in nonrational) == 30
    n2 = nonrational[0]
    assert n2 == {
        "n": 2,
        "full_trace_K3_degree_ceiling": 6,
        "actual_field_degree_bound": 6,
        "chordless_length_bound": 5,
        "selected_generator_count": 2,
        "primitive_grid_side": 16,
        "small_coordinate_bound": 64,
        "primitive_bound_factor": 32,
        "discriminant_exponent": 30,
        "nonrational_possible": True,
    }

    # Exact domination of every dimension-specific base by the coarse base
    # for arbitrary x=C/s0>=0.  Break points suffice because both sides are
    # piecewise linear maxima.
    hostile_ratios = [
        Fraction(0),
        Fraction(1),
        Fraction(64),
        Fraction(128),
        Fraction(8192),
        Fraction(8193),
        Fraction(10**6),
    ]
    for record in nonrational:
        factor = int(record["primitive_bound_factor"])
        small = int(record["small_coordinate_bound"])
        exponent = int(record["discriminant_exponent"])
        for ratio in hostile_ratios:
            local_base = factor * max(Fraction(small), ratio)
            coarse_base = 32 * max(Fraction(8192), ratio)
            assert local_base <= coarse_base
            assert exponent <= 30

    return {
        "records": records,
        "nonrational_dimension_range": [2, 9],
        "maximum_primitive_factor": 32,
        "maximum_small_coordinate_bound": 8192,
        "maximum_discriminant_exponent": 30,
        "uniform_formula": "(32*max(8192,C/s0))^30",
        "n2_formula": "(32*max(64,C/s0))^30",
    }


def primitive_grid_certificate() -> list[dict[str, int]]:
    records = []
    for degree in range(2, 7):
        selected = floor(log2(degree))
        bad_hyperplanes = degree * (degree - 1) // 2
        side = bad_hyperplanes + 1
        total_grid = side**selected
        union_bound = bad_hyperplanes * side ** (selected - 1)
        assert union_bound < total_grid
        records.append(
            {
                "degree_ceiling": degree,
                "selected_generators": selected,
                "embedding_pair_hyperplanes": bad_hyperplanes,
                "grid_side": side,
                "grid_points": total_grid,
                "union_bound_for_bad_points": union_bound,
                "vandermonde_exponent": degree * (degree - 1),
            }
        )
    return records


def dominant_height_certificate() -> dict[str, object]:
    # q(H,H)=p-sum(a_i), with p>0 and each non-distinguished summand -a_i.
    # q(H,tau_theta H)-theta_0 q(H,H)
    #     =sum_i (theta_0-theta_i)a_i >=0.
    samples = 0
    for positive in (10, 25, 100):
        for negatives in ((1,), (1, 2), (2, 3, 5)):
            if positive <= sum(negatives):
                continue
            for distinguished in (3, 11, 29):
                other = tuple(Fraction(distinguished, j + 2) for j in range(len(negatives)))
                mixed = distinguished * positive - sum(
                    theta * amount for theta, amount in zip(other, negatives)
                )
                square = positive - sum(negatives)
                gap = mixed - distinguished * square
                expected = sum(
                    (distinguished - theta) * amount
                    for theta, amount in zip(other, negatives)
                )
                assert gap == expected >= 0
                samples += 1
    return {
        "identity": (
            "q(H,tau_theta H)-theta0*q(H,H)="
            "sum_(sigma!=0)(theta0-theta_sigma)*(-q_sigma(H_sigma))"
        ),
        "exact_samples": samples,
        "dominance_and_negative_definite_other_summands_imply_height": True,
    }


def verify_owner_archive() -> dict[str, object]:
    envelope = json.loads(
        (ROOT / OWNER_PINS["archive"]["path"]).read_text(encoding="utf-8")
    )
    payload = envelope["payload"]
    assert envelope["payload_sha256"] == canonical_digest(payload)
    assert payload["report_sha256"] == OWNER_PINS["report"]["sha256"]
    assert payload["dependency_pins_match"] is True
    rows = payload["field_rows_C256_s2"]
    assert [row["dimension"] for row in rows] == list(range(2, 26))
    for index, independent in enumerate(
        field_constant_record(n) for n in range(2, 26)
    ):
        owner = rows[index]
        assert owner["full_trace_K3_degree_ceiling"] == independent[
            "full_trace_K3_degree_ceiling"
        ]
        assert owner["field_degree_bound"] == independent[
            "actual_field_degree_bound"
        ]
        assert owner["nonrational_field_possible"] == independent[
            "nonrational_possible"
        ]
        if independent["nonrational_possible"]:
            assert owner["chordless_cycle_length_bound"] == independent[
                "chordless_length_bound"
            ]
            assert owner["selected_generator_count_bound"] == independent[
                "selected_generator_count"
            ]
            assert owner["primitive_combination_coefficient_bound"] == independent[
                "primitive_grid_side"
            ]
            assert owner["small_affine_cycle_conjugate_bound"] == independent[
                "small_coordinate_bound"
            ]
            assert owner["field_discriminant_exponent"] == independent[
                "discriminant_exponent"
            ]
        else:
            assert owner["ground_field"] == "Q"
    return {
        "payload_digest_verified": True,
        "report_pin_verified": True,
        "dimensions_replayed": [2, 25],
        "n20_trace_ceiling": rows[18]["full_trace_K3_degree_ceiling"],
        "n20_actual_field_degree_bound": rows[18]["field_degree_bound"],
    }


def scope_certificate() -> dict[str, bool]:
    report = " ".join(
        (ROOT / OWNER_PINS["report"]["path"]).read_text(encoding="utf-8").split()
    )
    required = {
        "conditional_not_raw_compiler": (
            "This is not an assertion that every maximal arithmetic Coxeter chamber passes the gates."
        ),
        "full_picard_complex_K3": "exact full-Picard complex K3 trace carrier",
        "H_in_full_trace": r"big-and-nef class \(H\in T\)",
        "trace_sign": r"\(-q_\sigma\), is positive definite",
        "complex_weak_log": "smooth projective complex rational weak-log carrier",
        "D_below_two": r"If \(D_n<2\), then \(k=\mathbf Q\)",
        "signed_positive_negative": "difference of positive and negative effective parts",
        "both_parts_bounded": "Hilbert/Chow degree of **both** parts",
        "denominators_bounded": "bound every rational denominator",
        "fields_not_hosts": "bounded field discriminant does not by itself bound",
        "n1_counterexample_scope": "once hyperbolic dimension one is allowed",
        "n2_automatic_height_open": (
            r"does not prove that Gate 5 is logically necessary in the present \(n\ge2\) scope"
        ),
        "minus_four_open": r"twice-primitive \((-4)\) lane",
    }
    for label, phrase in required.items():
        assert phrase in report, (label, phrase)
    return {label: True for label in required}


def build_audit() -> dict[str, object]:
    pins = verify_pins()
    constants = uniform_constant_certificate()
    archive = verify_owner_archive()
    scope = scope_certificate()
    return {
        "schema": "enriched_canonical_chamber_ag_absw_transfer/independent_hostile/v1",
        "verdict": "PASS_CONDITIONAL",
        "pins_verified": pins,
        "dimension_transfer": {
            "direct_rational_dimension": 11,
            "weak_log_complex_rational_dimension": 25,
            "complex_full_trace_rank_ceiling": 20,
            "nonrational_dimension_ceiling": 9,
            "global_conditional_dimension_ceiling": 25,
            "ordinary_exact_property4_dimensions": [3, 9],
        },
        "full_picard_height": dominant_height_certificate(),
        "field_constants": constants,
        "primitive_grid": primitive_grid_certificate(),
        "owner_archive": archive,
        "scope_checks": scope,
        "logical_separation": {
            "finite_field_set_proved_under_uniform_height": True,
            "finite_fields_imply_finite_hosts": False,
            "finite_hosts_plus_relative_saturation_give_finite_chambers": True,
            "gate5_automatic_from_exact_maximality_in_n_ge_2": False,
            "raw_essential_surjectivity": False,
        },
    }


def main() -> None:
    print(json.dumps(build_audit(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
