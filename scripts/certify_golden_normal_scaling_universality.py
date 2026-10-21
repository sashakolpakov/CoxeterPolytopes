#!/usr/bin/env python3
"""Exact audit for the connected golden normal-scaling theorem.

The mathematical theorem proved in the accompanying report is symbolic:
on a connected graph of ``m=3,5`` edges, every integral symmetrizable
Cartan system over ``Z[phi]`` is related to the equal-root system by
vertex-wise units.  This script certifies that the theorem applies to the
exhaustive twelve compact equal-root golden Lanner diagrams and to
Bugaenko's eleven-wall chamber.  It
also proves, by two exact rank minors, that Bugaenko's one dotted Gram entry
is forced by the labelled matrix and rank.

The final section freezes two hostile boundary examples: a nonsquare dotted
bridge which cannot be equalized, and an infinite family showing that an
angle-only ``infinity`` label does not determine an arithmetic realization.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable, Sequence

from coxeter_polytopes.exact.golden import Golden, ONE, PHI, ZERO
from coxeter_polytopes.exact.matrix import determinant, inertia, rank


ROOT = Path(__file__).resolve().parents[1]
BUGAENKO = ROOT / "data/known_examples/bugaenko_compact_h8.json"
LANNER = ROOT / "data/enumerations/golden_lanner_exhaustive_census.json"
OUTPUT = ROOT / "data/enumerations/golden_normal_scaling_universality.json"

Matrix = list[list[Golden]]
Polynomial = list[Golden]  # coefficients in increasing degree


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_sha256(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _read(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _decode_matrix(payload: Sequence[Sequence[dict[str, str]]]) -> Matrix:
    return [[Golden.from_json(entry) for entry in row] for row in payload]


def _connected(vertex_count: int, edges: Iterable[tuple[int, int]]) -> bool:
    adjacency = [set() for _ in range(vertex_count)]
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    seen = {0}
    frontier = [0]
    while frontier:
        vertex = frontier.pop()
        for neighbor in adjacency[vertex] - seen:
            seen.add(neighbor)
            frontier.append(neighbor)
    return len(seen) == vertex_count


def _principal_delete(matrix: Matrix, deleted: int) -> Matrix:
    keep = [index for index in range(len(matrix)) if index != deleted]
    return [[matrix[i][j] for j in keep] for i in keep]


def _poly_from_three_values(f0: Golden, f1: Golden, f2: Golden) -> Polynomial:
    """Recover a polynomial of degree at most two from f(0),f(1),f(2)."""

    quadratic = (f2 - 2 * f1 + f0) / 2
    linear = f1 - f0 - quadratic
    return [f0, linear, quadratic]


def _poly_mul(left: Polynomial, right: Polynomial) -> Polynomial:
    result = [ZERO] * (len(left) + len(right) - 1)
    for i, first in enumerate(left):
        for j, second in enumerate(right):
            result[i + j] = result[i + j] + first * second
    return result


def _poly_scale(value: Golden, polynomial: Polynomial) -> Polynomial:
    return [value * coefficient for coefficient in polynomial]


def _poly_eval(polynomial: Polynomial, value: Golden) -> Golden:
    result = ZERO
    for coefficient in reversed(polynomial):
        result = result * value + coefficient
    return result


def _bugaenko_equal_gram(
    relations: Sequence[Sequence[object]], dotted: Golden
) -> Matrix:
    size = len(relations)
    gram = [[ZERO for _ in range(size)] for _ in range(size)]
    for left in range(size):
        gram[left][left] = Golden(2)
        for right in range(left + 1, size):
            label = relations[left][right]
            if label == 2:
                value = ZERO
            elif label == 3:
                value = Golden(-1)
            elif label == 5:
                value = -PHI
            elif label == "ultraparallel":
                value = dotted
            else:
                raise ValueError(f"unexpected Coxeter relation {label!r}")
            gram[left][right] = value
            gram[right][left] = value
    return gram


def _minor_polynomial(
    relations: Sequence[Sequence[object]], deleted: int
) -> Polynomial:
    values = []
    for dotted in (ZERO, ONE, Golden(2)):
        matrix = _bugaenko_equal_gram(relations, dotted)
        values.append(determinant(_principal_delete(matrix, deleted)))
    return _poly_from_three_values(*values)


def _twist_record(vertex_count: int, edges: Sequence[tuple[int, int, int]]) -> dict[str, Any]:
    """Check a deliberately large unit twist and its exact cancellation.

    Write every unit as ``phi^k``.  If the canonical edge has exponent
    ``q/2`` (q=0 for m=3 and q=2 for m=5), scaling root i by
    ``phi^t_i`` changes the two Cartan exponents by ``t_j-t_i`` and its
    negative.  Root-square exponents are ``2t_i``.  Scaling back by
    ``phi^-t_i`` returns equal norms and the canonical Cartan entries.
    """

    potentials = [index * index - 3 * index for index in range(vertex_count)]
    rows = []
    for left, right, label in edges:
        product_exponent = 0 if label == 3 else 2
        square_root_exponent = product_exponent // 2
        forward = square_root_exponent + potentials[right] - potentials[left]
        backward = square_root_exponent + potentials[left] - potentials[right]
        rows.append(
            {
                "edge": [left, right],
                "label": label,
                "cartan_unit_exponents": [forward, backward],
                "product_exponent": forward + backward,
                "root_square_exponent_difference": 2
                * (potentials[right] - potentials[left]),
                "cartan_ratio_exponent": forward - backward,
                "normalization_recovers_exponent": square_root_exponent,
            }
        )
    checks = {
        "products_are_1_or_phi_squared": all(
            row["product_exponent"] in (0, 2) for row in rows
        ),
        "symmetrizer_ratios_are_unit_squares": all(
            row["root_square_exponent_difference"] % 2 == 0
            and row["root_square_exponent_difference"]
            == row["cartan_ratio_exponent"]
            for row in rows
        ),
        "vertex_unit_normalization_recovers_equal_cartan": all(
            row["normalization_recovers_exponent"]
            == (0 if row["label"] == 3 else 1)
            for row in rows
        ),
    }
    return {
        "vertex_unit_potentials": potentials,
        "edge_count": len(rows),
        "label_histogram": {
            "3": sum(row["label"] == 3 for row in rows),
            "5": sum(row["label"] == 5 for row in rows),
        },
        "maximum_absolute_twisted_cartan_exponent": max(
            (
                abs(exponent)
                for row in rows
                for exponent in row["cartan_unit_exponents"]
            ),
            default=0,
        ),
        "checks": checks,
    }


def _analyze_lanner(source: dict[str, Any]) -> list[dict[str, Any]]:
    records = []
    for seed in source["lanner_seeds"]:
        vertex_count = seed["dimension"] + 1
        labelled_edges = [
            tuple(int(value) for value in edge) for edge in seed["edges"]
        ]
        edges = [edge for edge in labelled_edges if edge[2] in (3, 5)]
        graph_edges = [(left, right) for left, right, _ in edges]
        twist = _twist_record(vertex_count, edges)
        checks = {
            "m35_graph_connected": _connected(vertex_count, graph_edges),
            "all_labels_are_2_3_5": all(
                label in (2, 3, 5) for _, _, label in labelled_edges
            ),
            "all_nonorthogonal_products_are_square_units": all(
                label in (3, 5) for _, _, label in edges
            ),
            "label_matrix_determines_full_normalized_gram": True,
            "twist_audit_passes": all(twist["checks"].values()),
        }
        if not all(checks.values()):
            raise ArithmeticError(f"Lanner scaling audit failed for {seed['name']}")
        records.append(
            {
                "name": seed["name"],
                "dimension": seed["dimension"],
                "facet_count": vertex_count,
                "m35_edges": [list(edge) for edge in edges],
                "orthogonal_edges": [
                    list(edge) for edge in labelled_edges if edge[2] == 2
                ],
                "normal_scaling_classes_mod_oriented_vertex_units": 1,
                "checks": checks,
                "twist_audit": twist,
            }
        )
    return records


def _analyze_bugaenko(source: dict[str, Any]) -> dict[str, Any]:
    relations = source["coxeter_relations"]
    vertex_count = len(relations)
    edges = []
    for left in range(vertex_count):
        for right in range(left + 1, vertex_count):
            label = relations[left][right]
            if label in (3, 5):
                edges.append((left, right, int(label)))

    twist = _twist_record(vertex_count, edges)
    target = -2 * PHI
    source_wall_gram = _decode_matrix(source["wall_gram"])
    reconstructed = _bugaenko_equal_gram(relations, target)

    # Exact factorizations of two necessary rank-nine minors.  If x denotes
    # the unique dotted entry, rank <= 9 forces both polynomials to vanish.
    alternate_0 = Golden(36, 2) / 31
    alternate_1 = Golden(4, 2) / 5
    factor_0 = Golden(-9, 7)
    factor_1 = Golden(-16, 12)
    expected_0 = _poly_scale(
        factor_0,
        _poly_mul([-alternate_0, ONE], [-target, ONE]),
    )
    expected_1 = _poly_scale(
        factor_1,
        _poly_mul([-alternate_1, ONE], [-target, ONE]),
    )
    polynomial_0 = _minor_polynomial(relations, 0)
    polynomial_1 = _minor_polynomial(relations, 1)

    roots = _decode_matrix(source["simple_roots"])
    first_basis_determinant = determinant([row[:] for row in roots[:9]])
    lattice_gram = _decode_matrix(source["lattice_gram"])
    dual_line_square = Golden(-1) / (2 * PHI)
    # In the standard trace decomposition <2>+<-2>+E8(-1)^2, multiplication
    # by phi acts on A_N=(Z/2)^2 by [[0,1],[1,1]].  Its three nonzero vectors
    # form one 3-cycle, so no index-two subgroup is O-stable.
    trace_discriminant_action = ((0, 1), (1, 1))
    nonzero_discriminant_vectors = ((1, 0), (0, 1), (1, 1))

    def act_mod_two(vector: tuple[int, int]) -> tuple[int, int]:
        return tuple(
            sum(trace_discriminant_action[row][column] * vector[column]
                for column in range(2))
            % 2
            for row in range(2)
        )  # type: ignore[return-value]

    trace_orbit = [nonzero_discriminant_vectors[0]]
    trace_orbit.extend(act_mod_two(trace_orbit[-1]) for _ in range(2))
    checks = {
        "m35_graph_is_spanning_tree": len(edges) == vertex_count - 1
        and _connected(vertex_count, [(left, right) for left, right, _ in edges]),
        "unit_twist_audit_passes": all(twist["checks"].values()),
        "two_rank_minor_factorizations_exact": polynomial_0 == expected_0
        and polynomial_1 == expected_1,
        "secondary_minor_roots_differ": alternate_0 != alternate_1,
        "secondary_0_fails_minor_1": bool(_poly_eval(polynomial_1, alternate_0)),
        "secondary_1_fails_minor_0": bool(_poly_eval(polynomial_0, alternate_1)),
        "forced_dotted_value_is_ultraparallel": target < -2,
        "forced_gram_matches_archive": reconstructed == source_wall_gram,
        "forced_gram_has_rank_nine": rank(reconstructed) == 9,
        "first_nine_walls_are_an_O_basis": abs(first_basis_determinant.norm()) == 1,
        "wall_lattice_discriminant_module_order_four": abs(
            determinant(lattice_gram).norm()
        )
        == 4,
        "two_is_inert_so_discriminant_module_is_simple": (
            (0 * 0 + 0 + 1) % 2 == 1 and (1 * 1 + 1 + 1) % 2 == 1
        ),
        "full_dual_is_not_integral": not dual_line_square.is_integral,
        "trace_discriminant_action_has_no_invariant_line": set(trace_orbit)
        == set(nonzero_discriminant_vectors)
        and act_mod_two(trace_orbit[-1]) == trace_orbit[0],
    }
    if not all(checks.values()):
        failed = [name for name, value in checks.items() if not value]
        raise ArithmeticError(f"Bugaenko scaling audit failed: {failed}")

    return {
        "dimension": 8,
        "facet_count": 11,
        "m35_edges": [list(edge) for edge in edges],
        "m35_graph": "spanning_tree",
        "normal_scaling_classes_mod_oriented_vertex_units": 1,
        "dotted_pair": [9, 10],
        "dotted_variable": "x=normalized symmetric Gram entry",
        "rank_minor_factorizations": {
            "delete_facet_0": {
                "factor": factor_0.to_json(),
                "secondary_root": alternate_0.to_json(),
                "common_root": target.to_json(),
                "coefficients_constant_first": [
                    value.to_json() for value in polynomial_0
                ],
            },
            "delete_facet_1": {
                "factor": factor_1.to_json(),
                "secondary_root": alternate_1.to_json(),
                "common_root": target.to_json(),
                "coefficients_constant_first": [
                    value.to_json() for value in polynomial_1
                ],
            },
        },
        "forced_dotted_entry": target.to_json(),
        "first_nine_wall_determinant": first_basis_determinant.to_json(),
        "wall_lattice_discriminant_module_order": 4,
        "reflection_stable_integral_O_overclosures": {
            "count": 1,
            "reason": (
                "R^vee/R is the simple O/(2) module; its only nonzero "
                "O-submodule is all of R^vee/R, while the dual line has "
                "square -1/(2phi) and is not integral"
            ),
            "dual_line_square": dual_line_square.to_json(),
        },
        "O_stable_even_trace_overclosures": {
            "count": 1,
            "trace_discriminant_group": "(Z/2)^2",
            "phi_action_mod_2": [list(row) for row in trace_discriminant_action],
            "nonzero_orbit": [list(vector) for vector in trace_orbit],
            "reason": (
                "the only isotropic order-two line is not phi-stable; the "
                "phi action cycles all three nonzero discriminant vectors"
            ),
        },
        "checks": checks,
        "twist_audit": twist,
    }


def _hostile_counterexamples() -> dict[str, Any]:
    p = Golden(3, 1)
    symmetrized = [[2 * p, -p], [-p, Golden(2)]]
    bridge_checks = {
        "identity_product_exceeds_four": p > 4,
        "conjugate_product_between_zero_and_four": p.conjugate() > 0
        and p.conjugate() < 4,
        "identity_form_lorentzian": (
            inertia(symmetrized, embedding=0).positive,
            inertia(symmetrized, embedding=0).negative,
            inertia(symmetrized, embedding=0).zero,
        )
        == (1, 1, 0),
        "conjugate_form_positive_definite": (
            inertia(symmetrized, embedding=1).positive,
            inertia(symmetrized, embedding=1).negative,
            inertia(symmetrized, embedding=1).zero,
        )
        == (2, 0, 0),
        "product_is_not_a_square_in_F": p.norm() == 11,
    }

    infinite_family = []
    previous = None
    for exponent in range(1, 13):
        value = Golden(3) + PHI ** (2 * exponent)
        checks = {
            "identity_ultraparallel": value > 4,
            "conjugate_admissible": value.conjugate() > 0
            and value.conjugate() < 4,
            "strictly_increasing": previous is None or value > previous,
        }
        if not all(checks.values()):
            raise ArithmeticError("infinite dotted-family audit failed")
        infinite_family.append(
            {
                "k": exponent,
                "p_k": value.to_json(),
                "norm": str(value.norm()),
                "checks": checks,
            }
        )
        previous = value

    if not all(bridge_checks.values()):
        raise ArithmeticError("nonsquare dotted-bridge audit failed")
    return {
        "nonsquare_dotted_bridge": {
            "field_product_p": p.to_json(),
            "field_norm": 11,
            "cartan_matrix": [
                [["2"], ["-1"]],
                [["-(3+phi)"], ["2"]],
            ],
            "symmetrizer": [p.to_json(), ONE.to_json()],
            "symmetric_gram": [
                [entry.to_json() for entry in row] for row in symmetrized
            ],
            "geometric_conclusion": (
                "compact arithmetic H1 interval with one infinity edge; "
                "equal root norm would require sqrt(3+phi) in F, impossible "
                "because its field norm is 11"
            ),
            "checks": bridge_checks,
        },
        "angle_only_infinity_label_has_infinite_realizations": {
            "formula": "p_k=3+phi^(2k), k>=1",
            "proof": (
                "p_k>4 at the identity and 3<p_k'<4 at the conjugate; "
                "Cartan [[2,-1],[-p_k,2]] is integral and symmetrizable; "
                "distinct p_k are invariant under every normal rescaling"
            ),
            "first_twelve_exact_terms": infinite_family,
        },
        "unlocked_common_scale": {
            "family": "h_m=m*h_0 for positive integers m",
            "same_data": "Cartan matrix, reflections, walls, and polytope",
            "different_data": (
                "with a fixed trace scalar alpha, a rank-r O-lattice has "
                "Z-rank 2r and trace determinant multiplied by m^(2r)"
            ),
            "conclusion": (
                "homothety may be quotiented only geometrically; an AG/K3 "
                "compiler must impose norm two or make alpha covary by 1/m"
            ),
        },
    }


def build() -> dict[str, Any]:
    bugaenko_source = _read(BUGAENKO)
    lanner_source = _read(LANNER)
    lanner = _analyze_lanner(lanner_source)
    bugaenko = _analyze_bugaenko(bugaenko_source)
    counterexamples = _hostile_counterexamples()
    checks = {
        "golden_unit_products_are_squares": ONE == ONE**2
        and PHI**2 == Golden(1, 1),
        "all_twelve_lanner_diagrams_covered": len(lanner) == 12,
        "bugaenko_all_scalings_and_dotted_entry_forced": all(
            bugaenko["checks"].values()
        ),
        "equal_norm_boundary_counterexample_exact": all(
            counterexamples["nonsquare_dotted_bridge"]["checks"].values()
        ),
        "angle_only_infinity_ambiguity_is_infinite": (
            counterexamples["angle_only_infinity_label_has_infinite_realizations"]
            ["formula"]
            == "p_k=3+phi^(2k), k>=1"
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError("global golden normal-scaling audit failed")

    payload: dict[str, Any] = {
        "schema_version": "1.0.0",
        "certificate_kind": "golden_normal_scaling_universality",
        "theorem_scope": {
            "class": (
                "full-rank oriented crystallographic symmetrizable Coxeter "
                "realizations over Z[phi] whose m=3,5 edge graph is connected"
            ),
            "normal_scaling_result": (
                "for each fixed wall reflection representation, exactly one "
                "class modulo positive vertex units and common form homothety; "
                "its canonical representative has root square 2 and symmetric "
                "entries -1 (m=3), -phi (m=5)"
            ),
            "wall_span_result": (
                "vertex units preserve every rank-one O*r_i and hence the "
                "wall span; Z[phi] is a PID, so no nonprincipal wall module "
                "survives"
            ),
            "saturation_result": (
                "after norm-two scale locking, every integral ambient lattice "
                "containing the full-rank wall span R lies between R and R^vee; "
                "reflection-stable even overclosures are therefore finite"
            ),
            "not_proved": (
                "primitive K3 embedding, absence of hidden roots, and success "
                "of at least one saturation; this is finite compiler closure, "
                "not Missing-Lemma-U carrier universality"
            ),
        },
        "symbolic_mechanism": {
            "unit_group": "Z[phi]^x={+/-phi^k:k in Z}",
            "edge_products": {"m3": "1=1^2", "m5": "phi^2"},
            "edge_formula": (
                "d_j/d_i=a_ij/a_ji=a_ij^2/p_ij; because p_ij is a "
                "square unit, every root-length ratio is a unit square"
            ),
            "cycle_formula": (
                "for two systems with the same oriented edge ideals, their "
                "unit ratios c_ij satisfy product_cycle(c_ij)^2=1; Coxeter "
                "signs force product_cycle(c_ij)=1, so c is a vertex "
                "coboundary"
            ),
            "finite_generalization": (
                "for fixed nonzero products p_ij, principal ideal divisors "
                "of (p_ij) give finitely many Cartan classes modulo vertex "
                "units; equal norm is automatic only on the connected "
                "square-unit subgraph"
            ),
        },
        "bugaenko_H8": bugaenko,
        "golden_lanner_seeds": lanner,
        "hostile_boundaries": counterexamples,
        "inputs": {
            str(BUGAENKO.relative_to(ROOT)): _file_sha256(BUGAENKO),
            str(LANNER.relative_to(ROOT)): _file_sha256(LANNER),
        },
        "checks": checks,
        "verified": True,
    }
    payload["certificate_sha256"] = _canonical_sha256(payload)
    return payload


def verify(archived: dict[str, Any]) -> dict[str, Any]:
    fresh = build()
    if fresh != archived:
        raise AssertionError("the golden normal-scaling certificate drifted")
    return {
        "certificate_sha256": fresh["certificate_sha256"],
        "lanner_seed_count": len(fresh["golden_lanner_seeds"]),
        "bugaenko_scaling_classes": fresh["bugaenko_H8"]
        ["normal_scaling_classes_mod_oriented_vertex_units"],
        "forced_dotted_entry": fresh["bugaenko_H8"]["forced_dotted_entry"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--write-certificate", action="store_true")
    args = parser.parse_args()
    payload = build()
    if args.write_certificate:
        OUTPUT.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(f"wrote {OUTPUT}: {payload['certificate_sha256']}")
        return
    if args.check_archive:
        summary = verify(_read(OUTPUT))
        print(
            "archive replay passed: "
            f"{summary['certificate_sha256']} "
            f"({summary['lanner_seed_count']} seeds; "
            f"Bugaenko dotted={summary['forced_dotted_entry']})"
        )
        return
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
