#!/usr/bin/env python3
"""Independent hostile replay of the JK-H12 maximal doubling-wall stocks.

This audit does not import the production certifier or read its output.  It
starts from the frozen 172-by-172 facet Gram matrix, classifies every pair in
integer arithmetic, and enumerates maximal cliques by a maximal-independent-
set recurrence on the conflict graph.  The recurrence is different from the
production Bron--Kerbosch implementation:

* maximal independent sets containing v are v plus those of G-N[v];
* those omitting v are the maximal independent sets of G-v that contain a
  neighbour of v (and therefore dominate v).

The carrier check separately reads the full-Picard wall coordinates, verifies
that their relation matrix agrees with the raw Gram, and recomputes the Smith
cores and integral reflection action.  Frozen Vinberg and exact-cone archives
are used only for the already independent finite-volume seed gate.
"""

from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form


ROOT = Path(__file__).resolve().parents[1]
RAW_GRAM = (
    ROOT / "data/enumerations/jacobian_kummer_h12_interim_facet_gram.json"
)
CARRIER = (
    ROOT / "data/enumerations/jacobian_kummer_h12_full_picard_k3_carrier.json"
)
VINBERG = (
    ROOT / "data/enumerations/jacobian_kummer_h12_vinberg_volume_certificate.json"
)
EXACT_CONE = (
    ROOT / "data/enumerations/jacobian_kummer_h12_exact_cone_certificate.json"
)
OUTPUT = (
    ROOT
    / "data/enumerations/"
    "jacobian_kummer_h12_maximal_disjoint_doubling_family_hostile_audit.json"
)

SCHEMA = "coxeter.jacobian_kummer.h12_maximal_disjoint_hostile_audit.v1"
STATUS_DATE = "2026-08-08"
PINNED_SHA256 = {
    "raw_172_facet_Gram": (
        "53c39642b81bbd7f5be172408bfecb857bd4ed57467c91bdbaed63f9e734fb92"
    ),
    "full_Picard_carrier": (
        "f2c57b9fb0423560c9e238541b0ac1a81ca141917c864929501ebdfa7f294ff2"
    ),
    "literal_Vinberg_seed_gate": (
        "c3474e33120cab541205225dd8df31beea7ae8006cd1667c4fd66408f9b931f1"
    ),
    "exact_cone_seed_gate": (
        "84d9660b5990d89fb4a78b73180fa815843ee7af46d4d8c56be73f51e37a7c0f"
    ),
}

DISJOINT_CONFLICT_MATCHING_ONE_BASED = (
    (31, 162),
    (32, 169),
    (33, 170),
    (36, 130),
    (39, 141),
    (42, 145),
    (48, 161),
    (51, 159),
    (56, 164),
    (58, 172),
    (59, 171),
    (60, 168),
)
ULTRAPARALLEL_CONFLICT_MATCHING_ONE_BASED = (
    (31, 172),
    (32, 169),
    (33, 170),
    (36, 159),
    (39, 161),
    (42, 164),
    (48, 168),
    (51, 167),
    (56, 171),
    (58, 166),
    (59, 165),
    (60, 163),
    (61, 112),
    (62, 162),
    (63, 158),
    (64, 160),
    (65, 150),
    (66, 154),
    (67, 156),
    (68, 155),
    (69, 153),
    (70, 152),
    (71, 147),
    (72, 149),
    (73, 157),
    (74, 151),
    (75, 148),
    (76, 146),
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_sha256(payload: Any) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def _content(values: Iterable[int]) -> int:
    nonzero = tuple(abs(int(value)) for value in values if int(value))
    return math.gcd(*nonzero) if nonzero else 0


def _pair(gram: Sequence[Sequence[int]], left: Sequence[int], right: Sequence[int]) -> int:
    return sum(
        int(left[row]) * int(gram[row][column]) * int(right[column])
        for row in range(len(gram))
        for column in range(len(gram))
    )


def _relation(pairing: int, left_square: int, right_square: int) -> str:
    if pairing < 0 or left_square >= 0 or right_square >= 0:
        raise ArithmeticError("invalid oriented hyperbolic Gram entry")
    square = pairing * pairing
    product = left_square * right_square
    if pairing == 0:
        return "m2"
    if 4 * square == product:
        return "m3"
    if 2 * square == product:
        return "m4"
    if 4 * square == 3 * product:
        return "m6"
    if square == product:
        return "parallel"
    if square > product:
        return "ultraparallel"
    raise ArithmeticError(
        f"unrecognized facet relation: p={pairing}, squares={left_square},{right_square}"
    )


def _load_inputs() -> dict[str, Any]:
    paths = {
        "raw_172_facet_Gram": RAW_GRAM,
        "full_Picard_carrier": CARRIER,
        "literal_Vinberg_seed_gate": VINBERG,
        "exact_cone_seed_gate": EXACT_CONE,
    }
    hashes = {name: _sha256(path) for name, path in paths.items()}
    if hashes != PINNED_SHA256:
        raise ArithmeticError(f"a hostile-audit input changed: {hashes}")
    raw = json.loads(RAW_GRAM.read_text(encoding="utf-8"))
    gram = tuple(tuple(int(value) for value in row) for row in raw["facet_gram"])
    if len(gram) != 172 or any(len(row) != 172 for row in gram):
        raise ArithmeticError("raw facet Gram is not 172 by 172")
    if any(gram[i][j] != gram[j][i] for i in range(172) for j in range(172)):
        raise ArithmeticError("raw facet Gram is not symmetric")
    return {
        "hashes": hashes,
        "raw_gram": gram,
        "carrier": json.loads(CARRIER.read_text(encoding="utf-8")),
        "vinberg": json.loads(VINBERG.read_text(encoding="utf-8")),
        "cone": json.loads(EXACT_CONE.read_text(encoding="utf-8")),
    }


def _raw_relations(gram: Sequence[Sequence[int]]) -> dict[str, Any]:
    relations: list[tuple[str, ...]] = []
    histogram: Counter[str] = Counter()
    for left in range(172):
        row: list[str] = []
        for right in range(172):
            if left == right:
                row.append("m1")
                continue
            label = _relation(gram[left][right], gram[left][left], gram[right][right])
            row.append(label)
            if right < left:
                histogram[label] += 1
        relations.append(tuple(row))
    doubling = tuple(
        wall for wall in range(172) if "m3" not in relations[wall]
    )
    expected = tuple(
        [30, 31, 32, 35, 38, 41, 47, 50, 55] + list(range(57, 172))
    )
    stock_histogram = Counter(
        relations[left][right]
        for position, left in enumerate(doubling)
        for right in doubling[:position]
    )
    checks = {
        "all_14706_pairs_classified": sum(histogram.values()) == 14_706,
        "whole_relation_census_exact": histogram
        == Counter(
            {
                "m2": 2_820,
                "m3": 144,
                "m4": 800,
                "parallel": 1_750,
                "ultraparallel": 9_192,
            }
        ),
        "only_m3_is_a_finite_odd_label": not (
            set(histogram) - {"m2", "m3", "m4", "m6", "parallel", "ultraparallel"}
        )
        and histogram.get("m6", 0) == 0,
        "exact_124_doubling_labels": doubling == expected,
        "doubling_stock_census_exact": stock_histogram
        == Counter({"m2": 264, "parallel": 1_138, "ultraparallel": 6_224}),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"raw relation replay failed: {checks}")
    return {
        "relations": tuple(relations),
        "doubling": doubling,
        "whole_histogram": dict(sorted(histogram.items())),
        "stock_histogram": dict(sorted(stock_histogram.items())),
        "checks": checks,
    }


def _maximal_independent_sets(conflicts: Sequence[int]) -> tuple[int, ...]:
    """Enumerate all maximal independent sets of a bitset graph exactly."""

    @lru_cache(maxsize=None)
    def visit(vertices: int) -> tuple[int, ...]:
        if vertices == 0:
            return (0,)
        candidates: list[int] = []
        remaining = vertices
        while remaining:
            bit = remaining & -remaining
            candidates.append(bit.bit_length() - 1)
            remaining ^= bit
        vertex = max(
            candidates,
            key=lambda item: ((conflicts[item] & vertices).bit_count(), -item),
        )
        bit = 1 << vertex
        neighbours = conflicts[vertex] & vertices
        containing = tuple(
            independent | bit
            for independent in visit(vertices & ~bit & ~neighbours)
        )
        omitting = tuple(
            independent
            for independent in visit(vertices & ~bit)
            if independent & neighbours
        )
        result = containing + omitting
        if len(result) != len(set(result)):
            raise ArithmeticError("maximal-independent-set recurrence duplicated output")
        return result

    return tuple(sorted(visit((1 << len(conflicts)) - 1), key=lambda item: (item.bit_count(), item)))


def _graph_audit(
    raw: dict[str, Any],
    *,
    name: str,
    allowed: Callable[[str], bool],
    expected_histogram: Counter[int],
    expected_degrees: Counter[int],
    expected_maximum_one_based: tuple[int, ...],
    conflict_matching_one_based: tuple[tuple[int, int], ...],
) -> dict[str, Any]:
    stock: tuple[int, ...] = raw["doubling"]
    position = {wall: index for index, wall in enumerate(stock)}
    relations: tuple[tuple[str, ...], ...] = raw["relations"]
    conflicts = [0] * len(stock)
    edge_count = 0
    for left_position, left in enumerate(stock):
        for right_position, right in enumerate(stock[:left_position]):
            if allowed(relations[left][right]):
                edge_count += 1
            else:
                conflicts[left_position] |= 1 << right_position
                conflicts[right_position] |= 1 << left_position

    cliques_as_bits = _maximal_independent_sets(conflicts)
    cliques = tuple(
        tuple(stock[index] for index in range(len(stock)) if bits & (1 << index))
        for bits in cliques_as_bits
    )
    size_histogram = Counter(map(len, cliques))
    maximum_size = max(size_histogram)
    maxima = tuple(clique for clique in cliques if len(clique) == maximum_size)
    expected_maximum = tuple(wall - 1 for wall in expected_maximum_one_based)
    degrees = Counter(len(stock) - 1 - conflicts[index].bit_count() for index in range(len(stock)))

    matching = tuple((left - 1, right - 1) for left, right in conflict_matching_one_based)
    matching_vertices = tuple(wall for pair in matching for wall in pair)
    matching_is_conflicting = all(
        conflicts[position[left]] & (1 << position[right]) for left, right in matching
    )
    matching_bound = len(stock) - len(matching)
    selected_histogram = Counter(
        relations[left][right]
        for left_position, left in enumerate(expected_maximum)
        for right in expected_maximum[:left_position]
    )
    one_based_cliques = [[wall + 1 for wall in clique] for clique in cliques]
    checks = {
        "maximal_clique_census_exact": size_histogram == expected_histogram,
        "degree_census_exact": degrees == expected_degrees,
        "unique_maximum_exact": maxima == (expected_maximum,),
        "every_output_is_a_clique": all(
            allowed(relations[left][right])
            for clique in cliques
            for left_position, left in enumerate(clique)
            for right in clique[:left_position]
        ),
        "every_output_is_maximal": all(
            all(
                any(not allowed(relations[other][wall]) for wall in clique)
                for other in stock
                if other not in clique
            )
            for clique in cliques
        ),
        "matching_vertices_distinct": len(set(matching_vertices)) == 2 * len(matching),
        "matching_edges_are_conflicts": matching_is_conflicting,
        "matching_upper_bound_closes": matching_bound == maximum_size,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"{name} hostile graph replay failed: {checks}")
    return {
        "relation": name,
        "vertex_count": len(stock),
        "edge_count": edge_count,
        "degree_histogram": {str(key): value for key, value in sorted(degrees.items())},
        "maximal_clique_count": len(cliques),
        "maximal_clique_size_histogram": {
            str(key): value for key, value in sorted(size_histogram.items())
        },
        "all_maximal_cliques_sha256": _json_sha256(one_based_cliques),
        "maximum_size": maximum_size,
        "unique_maximum_one_based": list(expected_maximum_one_based),
        "unique_maximum_relation_histogram": dict(sorted(selected_histogram.items())),
        "conflict_matching_one_based": [list(pair) for pair in conflict_matching_one_based],
        "conflict_matching_size": len(matching),
        "matching_upper_bound": matching_bound,
        "checks": checks,
    }


def _smith_diagonal(columns: Sequence[Sequence[int]]) -> list[int]:
    matrix = sp.Matrix.hstack(*(sp.Matrix(column) for column in columns))
    smith = smith_normal_form(matrix)
    return [
        abs(int(smith[index, index]))
        for index in range(min(smith.rows, smith.cols))
        if smith[index, index]
    ]


def _reflection_check(gram: sp.Matrix, root: sp.Matrix) -> bool:
    square = int((root.T * gram * root)[0])
    covector = root.T * gram
    reflection = sp.eye(gram.rows) - sp.Rational(2, square) * root * covector
    return bool(
        all(sp.Rational(value).q == 1 for value in reflection)
        and reflection**2 == sp.eye(gram.rows)
        and reflection.T * gram * reflection == gram
        and reflection * root == -root
    )


def _carrier_audit(
    inputs: dict[str, Any], raw: dict[str, Any], disjoint: dict[str, Any], ultra: dict[str, Any]
) -> dict[str, Any]:
    carrier_payload = inputs["carrier"]
    carrier = carrier_payload["chamber"]["carrier"]
    wall_system = carrier_payload["chamber"]["oriented_wall_system"]
    records = wall_system["wall_records"]
    gram = sp.Matrix(carrier["gram"])
    roots = tuple(tuple(int(value) for value in row["oriented_coordinates"]) for row in records)
    raw_relations: tuple[tuple[str, ...], ...] = raw["relations"]
    carrier_relations: list[list[str]] = [["m1"] * 172 for _ in range(172)]
    squares = tuple(_pair(carrier["gram"], root, root) for root in roots)
    for left in range(172):
        for right in range(left):
            label = _relation(
                _pair(carrier["gram"], roots[left], roots[right]),
                squares[left],
                squares[right],
            )
            carrier_relations[left][right] = carrier_relations[right][left] = label

    doubling: tuple[int, ...] = raw["doubling"]
    disjoint_core = tuple(range(60))
    ultra_core = tuple(range(76))
    disjoint_smith = _smith_diagonal([roots[index] for index in disjoint_core])
    ultra_smith = _smith_diagonal([roots[index] for index in ultra_core])
    checks = {
        "carrier_is_rank13_signature_1_12": carrier["rank"] == 13
        and carrier["signature"] == [1, 12, 0],
        "carrier_is_declared_full_Picard_family": carrier["name"]
        == "U(2)+D5(-1)^2+<-4>"
        and carrier["very_general_period_consequence"]["moduli_dimension"] == 7,
        "all_172_wall_records_present_in_order": len(records) == 172
        and [record["facet"] for record in records] == list(range(1, 173)),
        "all_carrier_roots_primitive_negative": all(
            _content(root) == 1 and square < 0 for root, square in zip(roots, squares)
        ),
        "carrier_and_raw_relation_matrices_agree": all(
            carrier_relations[left][right] == raw_relations[left][right]
            for left in range(172)
            for right in range(172)
        ),
        "all_124_doubling_reflections_integral_on_Picard": all(
            _reflection_check(gram, sp.Matrix(roots[index])) for index in doubling
        ),
        "112_stock_complement_is_W1_through_W60": set(range(172))
        - {wall - 1 for wall in disjoint["unique_maximum_one_based"]}
        == set(disjoint_core),
        "96_stock_complement_is_W1_through_W76": set(range(172))
        - {wall - 1 for wall in ultra["unique_maximum_one_based"]}
        == set(ultra_core),
        "W1_through_W60_has_Smith_1_power_13": disjoint_smith == [1] * 13,
        "W1_through_W76_has_Smith_1_power_13": ultra_smith == [1] * 13,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"hostile carrier replay failed: {checks}")
    return {
        "lattice": carrier["name"],
        "rank": carrier["rank"],
        "signature": carrier["signature"],
        "moduli_dimension": carrier["very_general_period_consequence"]["moduli_dimension"],
        "disjoint_112_exterior_core": {
            "walls_one_based": list(range(1, 61)),
            "rank": 13,
            "Smith_diagonal": disjoint_smith,
            "index": 1,
        },
        "ultraparallel_96_exterior_core": {
            "walls_one_based": list(range(1, 77)),
            "rank": 13,
            "Smith_diagonal": ultra_smith,
            "index": 1,
        },
        "inheritance_scope": (
            "weak full-Picard divisor-reflection carrier: unglued identity walls "
            "remain exterior, while every translated wall is a primitive root "
            "under an integral Picard isometry; no nef, effectivity, Hodge, or "
            "surface-automorphism statement is inferred"
        ),
        "checks": checks,
    }


def _volume_and_allcock_audit(inputs: dict[str, Any], ultra: dict[str, Any]) -> dict[str, Any]:
    vinberg = inputs["vinberg"]
    cone = inputs["cone"]
    seed_checks = {
        "literal_Vinberg_gate_finite_volume_noncompact": vinberg["vinberg_decision"]
        == "finite_volume_noncompact"
        and vinberg["finite_volume"]
        and not vinberg["compact"]
        and vinberg["bad_edge_count"] == 0
        and vinberg["unrecognized_endpoint_edge_count"] == 0
        and vinberg["invalid_full_zero_neighborhood_count"] == 0
        and vinberg["ideal_vertex_count"] == 6_214,
        "exact_cone_gate_has_no_hyperideal_ray": cone["decision"]["classification"]
        == "finite_volume_noncompact"
        and cone["decision"]["finite_volume"]
        and not cone["decision"]["compact"]
        and cone["exact_double_description"]["hyperideal_vertex_count"] == 0
        and all(cone["checks"].values()),
        "three_ultraparallel_doubling_walls_available_directly": ultra["maximum_size"]
        >= 3,
    }
    if not all(seed_checks.values()):
        raise ArithmeticError(f"hostile Allcock/volume gate failed: {seed_checks}")
    return {
        "primary_source": {
            "author": "Daniel Allcock",
            "title": "Infinitely many hyperbolic Coxeter groups through dimension 19",
            "journal": "Geometry & Topology 10 (2006), 737-758",
            "doi": "10.2140/gt.2006.10.737",
            "primary_pdf": "https://web.ma.utexas.edu/users/allcock/research/doubling.pdf",
            "section_4_setup": (
                "k pairwise disjoint doubling-wall reflections generate the "
                "free product of k copies of C2; its Cayley graph is the "
                "k-regular tree"
            ),
            "Theorem_4_1": (
                "the union Q_T indexed by any subtree T of that Cayley tree "
                "is a Coxeter polyhedron"
            ),
            "Corollary_4_2": (
                "an I-vertex subtree gives a reflection subgroup of index I"
            ),
            "Theorem_4_3": (
                "if Q has finite volume and at least three pairwise disjoint "
                "doubling walls, the number N(I) of index-I reflection "
                "subgroups up to conjugacy by hyperbolic isometries is bounded "
                "below exponentially in I; the proof constructs nonisometric Q_T"
            ),
        },
        "verified_application": {
            "dimension": 12,
            "seed": "JK-H12-172",
            "seed_status": "finite_volume_noncompact",
            "generator_counts": [112, 96],
            "finite_subtree_consequence": (
                "every finite subtree gives a finite-sided finite-volume "
                "noncompact H12 Coxeter polytope"
            ),
            "index_and_volume": (
                "[W(Q):W(Q_T)]=|V(T)| and Vol(Q_T)=|V(T)|*Vol(Q)"
            ),
            "noncompact_reason": (
                "Q_T contains an isometric unbounded copy of the noncompact seed Q"
            ),
            "exponential_scope": (
                "Allcock 4.3 applies inside any chosen three of the 96 directly "
                "ultraparallel walls; it does not assert that every pair of "
                "nonisomorphic finite subtrees gives nonisometric polytopes"
            ),
        },
        "seed_gate": {
            "finite_vertices": vinberg["finite_vertex_count"],
            "cusps": vinberg["ideal_vertex_count"],
            "hyperideal_rays": cone["exact_double_description"]["hyperideal_vertex_count"],
            "checks": seed_checks,
        },
    }


def build() -> dict[str, Any]:
    inputs = _load_inputs()
    raw = _raw_relations(inputs["raw_gram"])
    disjoint = _graph_audit(
        raw,
        name="parallel_or_ultraparallel",
        allowed=lambda label: label in {"parallel", "ultraparallel"},
        expected_histogram=Counter({36: 9, 55: 36, 74: 42, 93: 12, 112: 1}),
        expected_degrees=Counter({99: 12, 120: 16, 121: 96}),
        expected_maximum_one_based=tuple(range(61, 173)),
        conflict_matching_one_based=DISJOINT_CONFLICT_MATCHING_ONE_BASED,
    )
    ultra = _graph_audit(
        raw,
        name="strictly_ultraparallel",
        allowed=lambda label: label == "ultraparallel",
        expected_histogram=Counter({32: 32, 57: 12, 61: 16, 96: 1}),
        expected_degrees=Counter({56: 12, 64: 16, 112: 96}),
        expected_maximum_one_based=tuple(range(77, 173)),
        conflict_matching_one_based=ULTRAPARALLEL_CONFLICT_MATCHING_ONE_BASED,
    )
    carrier = _carrier_audit(inputs, raw, disjoint, ultra)
    transfer = _volume_and_allcock_audit(inputs, ultra)
    checks = {
        "raw_relation_replay_green": all(raw["checks"].values()),
        "disjoint_graph_replay_green": all(disjoint["checks"].values()),
        "ultraparallel_graph_replay_green": all(ultra["checks"].values()),
        "carrier_and_Smith_replay_green": all(carrier["checks"].values()),
        "Allcock_and_volume_transfer_green": all(
            transfer["seed_gate"]["checks"].values()
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"hostile theorem audit failed: {checks}")
    return {
        "schema": SCHEMA,
        "status_date": STATUS_DATE,
        "verdict": "PASS_FREEZE",
        "all_checks_pass": True,
        "independence": {
            "production_certifier_imported": False,
            "production_certificate_read": False,
            "raw_172_facet_Gram_parsed_directly": True,
            "graph_method": (
                "memoized maximal-independent-set recurrence on conflict graph; "
                "not production Bron--Kerbosch"
            ),
            "carrier_coordinates_parsed_directly": True,
            "carrier_relation_matrix_cross_checked_against_raw_Gram": True,
        },
        "pinned_inputs": inputs["hashes"],
        "raw_relation_replay": {
            "facet_count": 172,
            "doubling_wall_count": len(raw["doubling"]),
            "doubling_walls_one_based": [wall + 1 for wall in raw["doubling"]],
            "whole_relation_histogram": raw["whole_histogram"],
            "doubling_stock_relation_histogram": raw["stock_histogram"],
            "checks": raw["checks"],
        },
        "maximal_disjoint_graph": disjoint,
        "maximal_ultraparallel_graph": ultra,
        "full_Picard_inheritance": carrier,
        "Allcock_transfer": transfer,
        "claim_boundary": {
            "maxima_scoped_to_exact_124_seed_doubling_walls": True,
            "arbitrary_descendant_roots_classified": False,
            "closed_arbitrary_tree_facet_formula_claimed": False,
            "worldwide_priority_claimed": False,
            "new_dimension_claimed": False,
            "new_commensurability_class_claimed": False,
            "full_nef_or_effectivity_claimed": False,
            "surface_automorphism_claimed": False,
        },
        "checks": checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    payload = build()
    if arguments.check_archive:
        archived = json.loads(arguments.output.read_text(encoding="utf-8"))
        if payload != archived:
            raise SystemExit("hostile audit archive differs from exact rebuild")
        print(f"hostile archive replay passed: {arguments.output}")
        return
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(arguments.output)
    print(f"sha256={_sha256(arguments.output)}")


if __name__ == "__main__":
    main()
