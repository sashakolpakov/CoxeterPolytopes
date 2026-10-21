#!/usr/bin/env python3
"""Exact arithmetic audit for the polytope-first golden K3 theorem.

The script checks the eight compact golden Lanner seeds used to prove that
the inverse theorem has a non-singleton positive class.  It verifies full and
vertex-link signatures at both real places, trace-lattice parity/signature,
Smith invariants, Nikulin's stable-range inequality, period dimensions, and
the paired trace-root arithmetic.  It also pins the Bugaenko H8 survivor and
the silver H9 failure to their existing frozen certificates.

Period surjectivity, Nikulin's primitive-embedding theorem, the K3 Weyl
theorem, and local finiteness of the arithmetic mirror arrangement are cited
mathematical inputs; this script does not pretend to prove them numerically.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Iterable

from sympy import Matrix as SympyMatrix
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

from coxeter_polytopes.exact.golden import Golden, PHI, ZERO
from coxeter_polytopes.exact.matrix import determinant, inertia
from coxeter_polytopes.exact.trace_lattice import (
    INV_SQRT5,
    integral_gram,
    is_even_integral,
    trace_gram,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/polytope_first_golden_k3_leverage.json"
BUGAENKO_POLYTOPE = ROOT / "data/known_examples/bugaenko_compact_h8.json"
BUGAENKO_TRACE = ROOT / "data/trace_lattices/bugaenko_l8_trace.json"
BUGAENKO_K3 = ROOT / "data/trace_lattices/bugaenko_k3_eigen_nef.json"
SILVER_K3 = ROOT / "data/candidates/silver_k3_carrier.json"
SILVER_STABILIZER = (
    ROOT / "data/candidates/silver_k3_norm2_infinite_stabilizer.json"
)

Edge = tuple[int, int, int]


SEEDS: tuple[tuple[str, int, tuple[Edge, ...], Golden], ...] = (
    (
        "triangle_2_5_5",
        2,
        ((0, 1, 2), (1, 2, 5), (0, 2, 5)),
        Golden(4, -4),
    ),
    (
        "triangle_3_3_5",
        2,
        ((0, 1, 3), (1, 2, 3), (0, 2, 5)),
        Golden(2, -4),
    ),
    (
        "triangle_3_5_5",
        2,
        ((0, 1, 3), (1, 2, 5), (0, 2, 5)),
        Golden(0, -6),
    ),
    (
        "triangle_5_5_5",
        2,
        ((0, 1, 5), (1, 2, 5), (0, 2, 5)),
        Golden(0, -10),
    ),
    (
        "chain_3_5_3",
        3,
        ((0, 1, 3), (1, 2, 5), (2, 3, 3)),
        Golden(5, -4),
    ),
    (
        "chain_5_3_5",
        3,
        ((0, 1, 5), (1, 2, 3), (2, 3, 5)),
        Golden(6, -5),
    ),
    (
        "chain_5_3_3_3",
        4,
        ((0, 1, 5), (1, 2, 3), (2, 3, 3), (3, 4, 3)),
        Golden(6, -4),
    ),
    (
        "chain_5_3_3_5",
        4,
        ((0, 1, 5), (1, 2, 3), (2, 3, 3), (3, 4, 5)),
        Golden(8, -6),
    ),
)


EXPECTED_TRACE_SMITH: dict[str, tuple[int, ...]] = {
    "triangle_2_5_5": (1, 1, 1, 1, 4, 4),
    "triangle_3_3_5": (1, 1, 1, 1, 2, 10),
    "triangle_3_5_5": (1, 1, 1, 1, 6, 6),
    "triangle_5_5_5": (1, 1, 1, 1, 10, 10),
    "chain_3_5_3": (1, 1, 1, 1, 1, 1, 1, 11),
    "chain_5_3_5": (1, 1, 1, 1, 1, 1, 1, 19),
    "chain_5_3_3_3": (1, 1, 1, 1, 1, 1, 1, 1, 2, 2),
    "chain_5_3_3_5": (1, 1, 1, 1, 1, 1, 1, 1, 2, 10),
}


def _canonical_sha256(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return sha256(encoded).hexdigest()


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _weight(label: int) -> Golden:
    if label == 2:
        return ZERO
    if label == 3:
        return Golden(1)
    if label == 5:
        return PHI
    raise ValueError(f"unsupported equal-root label {label}")


def field_gram(dimension: int, edges: Iterable[Edge]) -> list[list[Golden]]:
    rank = dimension + 1
    gram = [[ZERO for _ in range(rank)] for _ in range(rank)]
    for index in range(rank):
        gram[index][index] = Golden(2)
    for left, right, label in edges:
        value = -_weight(label)
        gram[left][right] = value
        gram[right][left] = value
    return gram


def _principal_delete(
    gram: list[list[Golden]], deleted: int
) -> list[list[Golden]]:
    keep = [index for index in range(len(gram)) if index != deleted]
    return [[gram[i][j] for j in keep] for i in keep]


def _smith_diagonal(gram: list[list[int]]) -> tuple[int, ...]:
    smith = smith_normal_form(SympyMatrix(gram), domain=ZZ)
    return tuple(abs(int(smith[index, index])) for index in range(smith.rows))


def _encode_inertia(value: object) -> list[int]:
    return [value.positive, value.negative, value.zero]  # type: ignore[attr-defined]


def _analyze_seed(
    name: str,
    dimension: int,
    edges: tuple[Edge, ...],
    expected_determinant: Golden,
) -> dict[str, Any]:
    gram = field_gram(dimension, edges)
    rank = dimension + 1
    field_determinant = determinant(gram)
    identity_inertia = inertia(gram, embedding=0)
    conjugate_inertia = inertia(gram, embedding=1)
    vertex_links = [inertia(_principal_delete(gram, index)) for index in range(rank)]

    traced = trace_gram(gram, -PHI * INV_SQRT5)
    trace_integer = integral_gram(traced)
    trace_inertia = inertia(trace_integer)
    trace_smith = _smith_diagonal(trace_integer)
    discriminant_length = sum(value > 1 for value in trace_smith)
    trace_rank = 2 * rank
    complement_rank = 22 - trace_rank
    period_dimension = 20 - trace_rank

    unit = PHI - 1
    root_square = ((-PHI * INV_SQRT5) * Golden(2)).trace()
    partner_square = (
        (-PHI * INV_SQRT5) * Golden(2) * unit * unit
    ).trace()
    root_partner_pairing = (
        (-PHI * INV_SQRT5) * Golden(2) * unit
    ).trace()

    checks = {
        "field_determinant_exact": field_determinant == expected_determinant,
        "identity_lorentzian": _encode_inertia(identity_inertia)
        == [dimension, 1, 0],
        "conjugate_positive_definite": _encode_inertia(conjugate_inertia)
        == [rank, 0, 0],
        "all_vertex_links_finite": all(
            _encode_inertia(link) == [dimension, 0, 0] for link in vertex_links
        ),
        "field_gram_even_integral": all(
            gram[index][index].a.denominator == 1
            and gram[index][index].b.denominator == 1
            and int(gram[index][index].a) % 2 == 0
            and int(gram[index][index].b) % 2 == 0
            for index in range(rank)
        )
        and all(entry.is_integral for row in gram for entry in row),
        "trace_even_integral": is_even_integral(traced),
        "trace_signature": _encode_inertia(trace_inertia)
        == [1, 2 * dimension + 1, 0],
        "trace_smith_exact": trace_smith == EXPECTED_TRACE_SMITH[name],
        "nikulin_stable_range": complement_rank >= discriminant_length + 2,
        "worst_case_automatic_range": complement_rank >= trace_rank + 2,
        "paired_trace_roots_are_orthogonal_minus_two": (
            root_square == -2
            and partner_square == -2
            and root_partner_pairing == 0
        ),
    }
    if not all(checks.values()):
        failed = [key for key, value in checks.items() if not value]
        raise ArithmeticError(f"golden Lanner audit failed for {name}: {failed}")

    return {
        "name": name,
        "dimension": dimension,
        "rank_over_Z_phi": rank,
        "edges": [list(edge) for edge in edges],
        "field_determinant": field_determinant.to_json(),
        "conjugate_field_determinant": field_determinant.conjugate().to_json(),
        "identity_inertia": _encode_inertia(identity_inertia),
        "conjugate_inertia": _encode_inertia(conjugate_inertia),
        "vertex_link_inertias": [_encode_inertia(link) for link in vertex_links],
        "trace_rank": trace_rank,
        "trace_determinant": int(SympyMatrix(trace_integer).det()),
        "trace_smith": list(trace_smith),
        "discriminant_length": discriminant_length,
        "k3_complement_rank": complement_rank,
        "k3_period_dimension": period_dimension,
        "checks": checks,
    }


def build() -> dict[str, Any]:
    seeds = [_analyze_seed(*seed) for seed in SEEDS]
    bugaenko_polytope = _read(BUGAENKO_POLYTOPE)
    bugaenko_trace = _read(BUGAENKO_TRACE)
    bugaenko_k3 = _read(BUGAENKO_K3)
    silver_k3 = _read(SILVER_K3)
    silver_stabilizer = _read(SILVER_STABILIZER)

    silver_model = silver_k3["carrier_models"][0]
    global_checks = {
        "eight_compact_lanner_seeds": len(seeds) == 8,
        "dimensions_2_3_4_represented": sorted({row["dimension"] for row in seeds})
        == [2, 3, 4],
        "automatic_k3_range_all_seeds": all(
            row["checks"]["worst_case_automatic_range"] for row in seeds
        ),
        "period_dimensions_14_12_10": sorted(
            {row["k3_period_dimension"] for row in seeds}, reverse=True
        )
        == [14, 12, 10],
        "bugaenko_stable_boundary": (
            bugaenko_polytope["dimension"] == 8
            and bugaenko_trace["rank"] == 18
            and bugaenko_trace["determinant"] == "-4"
            and bugaenko_k3["trace_lattice_primitive_k3_embedding_verified"]
            is True
            and bugaenko_k3["vinberg_completion_verified"] is True
            and bugaenko_k3["eigen_nef_realization"]
            ["restricted_chamber_equals_bugaenko_chamber"]
            is True
        ),
        "silver_rank20_carrier_but_infinite_volume": (
            silver_model["trace_inertia"] == [1, 19, 0]
            and silver_k3["k3_embedding"]["primitive_embedding_verified"] is True
            and silver_stabilizer["verified"] is True
            and silver_stabilizer["arithmetic_conclusion"]
            ["complete_chamber_hyperbolic_volume"]
            == "infinite"
        ),
        "quadratic_k3_full_trace_ceiling_H9": 2 * (9 + 1) == 20
        and 2 * (10 + 1) > 20,
    }
    if not all(global_checks.values()):
        failed = [key for key, value in global_checks.items() if not value]
        raise ArithmeticError(f"global golden K3 leverage audit failed: {failed}")

    payload: dict[str, Any] = {
        "schema_version": "1.0.0",
        "certificate_kind": "polytope_first_golden_k3_leverage",
        "theorem_scope": {
            "input": (
                "fixed compact even norm-two Z[phi] Coxeter realization; "
                "not an abstract diagram with unenumerated normal scalings"
            ),
            "automatic_result": (
                "every root-complete package in H2-H4 has a full-Picard K3 "
                "eigen-nef realization; compact Lanner seeds make the class "
                "nonempty and non-singleton"
            ),
            "stable_result": (
                "ell(A_N)<=18-2n is sufficient through H8; Bugaenko realizes "
                "equality at n=8"
            ),
            "nonexistence_scope": (
                "O-stable integral full-trace golden quadratic K3 carriers, "
                "not arithmetic Coxeter polytopes themselves, are impossible "
                "from H10 onward"
            ),
            "embedding_orbit_scope": (
                "Nikulin local discriminant forms decide existence/complement "
                "genera; O(Lambda_K3)-orbit enumeration additionally requires "
                "complement isometry classes and discriminant anti-isometry "
                "double cosets"
            ),
            "trace_even_wall_scope": (
                "for a general trace-even overclosure, visible roots range "
                "over the finite five-shell list and equality with P must be "
                "tested explicitly; identification with the norm-two field "
                "arrangement is restricted to field-even packages"
            ),
            "compactness_bound_scope": (
                "the two embedded field norms are the finitely many positive "
                "constants sigma(t) from the five-shell list, not uniformly 2"
            ),
            "external_theorems": [
                "Nikulin stable primitive-embedding theorem",
                "surjectivity of the marked K3 period map",
                "K3 Weyl-chamber transitivity and Riemann-Roch",
                "local finiteness of fixed-shell arithmetic mirrors",
            ],
        },
        "lanner_seeds": seeds,
        "paired_root_arithmetic": {
            "root_trace_square": -2,
            "unit_partner": "phi-1",
            "unit_partner_trace_square": -2,
            "trace_pairing": 0,
            "geometric_scope": (
                "for a generic restricted facet in an even root-complete "
                "package, shell uniqueness leaves a local A1+A1 root system; "
                "Weyl framing makes its two simple roots irreducible curves"
            ),
        },
        "high_dimensional_boundary": {
            "bugaenko_H8": {
                "trace_rank": 18,
                "discriminant_length": 2,
                "stable_bound_right_side": 2,
                "period_dimension": 2,
                "complete_chamber": "compact_11_facets",
            },
            "silver_H9": {
                "trace_rank": 20,
                "period_dimension": 0,
                "primitive_k3_carrier": True,
                "complete_chamber_volume": "infinite",
            },
            "full_trace_quadratic_k3_ceiling": 9,
        },
        "inputs": {
            str(path.relative_to(ROOT)): _file_sha256(path)
            for path in (
                BUGAENKO_POLYTOPE,
                BUGAENKO_TRACE,
                BUGAENKO_K3,
                SILVER_K3,
                SILVER_STABILIZER,
            )
        },
        "checks": global_checks,
        "verified": True,
    }
    payload["certificate_sha256"] = _canonical_sha256(payload)
    return payload


def verify(archived: dict[str, Any]) -> dict[str, Any]:
    fresh = build()
    if archived != fresh:
        raise AssertionError("the archived polytope-first golden K3 audit drifted")
    return {
        "certificate_sha256": fresh["certificate_sha256"],
        "seed_count": len(fresh["lanner_seeds"]),
        "dimensions": sorted(
            {row["dimension"] for row in fresh["lanner_seeds"]}
        ),
        "full_trace_quadratic_k3_ceiling": fresh["high_dimensional_boundary"]
        ["full_trace_quadratic_k3_ceiling"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--check-archive",
        action="store_true",
        help="compare the exact replay with the frozen JSON certificate",
    )
    args = parser.parse_args()
    payload = build()
    if args.check_archive:
        archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
        summary = verify(archived)
        print(
            "archive replay passed: "
            f"{summary['certificate_sha256']} "
            f"({summary['seed_count']} seeds)"
        )
        return
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
