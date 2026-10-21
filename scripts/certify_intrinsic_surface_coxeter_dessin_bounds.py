#!/usr/bin/env python3
"""Exact arithmetic certificate for the intrinsic-surface angle bounds.

This certifies the elementary part of
results/marked_facet_degree_mostow_dessin_bridge.md.  The geometric
Mostow--Prasad step is a proof theorem, not a numerical experiment.
"""

from __future__ import annotations

import argparse
import json
from typing import Any


def euler_phi(n: int) -> int:
    if n < 1:
        raise ValueError("n must be positive")
    result = n
    p = 2
    remaining = n
    while p * p <= remaining:
        if remaining % p == 0:
            while remaining % p == 0:
                remaining //= p
            result -= result // p
        p += 1
    if remaining > 1:
        result -= result // remaining
    return result


def finite_angle_alphabet(field_degree: int) -> list[int]:
    """Return {2} union {m>=3: phi(m)/2 <= field_degree}."""
    if field_degree < 1:
        raise ValueError("field degree must be positive")
    cutoff = 8 * field_degree * field_degree
    return [2] + [
        m for m in range(3, cutoff + 1)
        if euler_phi(m) <= 2 * field_degree
    ]


def verify_totient_cutoff(limit: int) -> None:
    """Check 2*phi(m)^2 >= m on a declared finite audit range."""
    for m in range(1, limit + 1):
        phi = euler_phi(m)
        if 2 * phi * phi < m:
            raise AssertionError((m, phi))


def build_certificate() -> dict[str, Any]:
    verify_totient_cutoff(100_000)
    k3_rank_bound = 20
    h9_dmax = k3_rank_bound // 10
    h10_dmax = k3_rank_bound // 11
    return {
        "theorem": "intrinsic_surface_coxeter_dessin_bounds",
        "totient_identity_consequence": "2*phi(m)^2 >= m",
        "totient_audit_limit": 100_000,
        "angle_alphabets": {
            str(d): finite_angle_alphabet(d) for d in range(1, 6)
        },
        "k3": {
            "picard_rank_bound": k3_rank_bound,
            "H9_max_field_degree": h9_dmax,
            "H9_compact_forces_field_degree": 2,
            "H9_finite_angle_alphabet": finite_angle_alphabet(2),
            "H10_max_field_degree_from_NS": h10_dmax,
            "H10_plus_compact_full_trace_possible": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pretty", action="store_true")
    args = parser.parse_args()
    print(json.dumps(
        build_certificate(),
        indent=2 if args.pretty else None,
        sort_keys=True,
    ))


if __name__ == "__main__":
    main()
