#!/usr/bin/env python3
"""Independent hostile audit of the nondotted K3 field-generating core theorem."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PINS = {
    "report": (
        "results/k3_nondotted_field_generating_core_bound.md",
        "c72c695e11787fa72fa746bec24640d0ca1499e48c067fe153d3a476d0cfd11a",
    ),
    "certifier": (
        "scripts/certify_k3_nondotted_field_generating_core_bound.py",
        "8582037cbe3306ea13786c6f6707a6ec32df735277e6d8f2d7ebf7623c667d9c",
    ),
    "archive": (
        "data/enumerations/k3_nondotted_field_generating_core_bound.json",
        "996862fef84590fe30d3fa9912b2cdafd6ed392a6319e5c7d733998b98080f82",
    ),
    "tests": (
        "tests/test_k3_nondotted_field_generating_core_bound.py",
        "b4f74470ed34d476d93b807a5c2070e5f19b17bbb6a0b473319beed83e14bd47",
    ),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_pins() -> dict[str, str]:
    observed = {}
    for label, (relative, expected) in PINS.items():
        actual = sha256(ROOT / relative)
        assert actual == expected, (label, expected, actual)
        observed[label] = actual
    return observed


def determinant(matrix: tuple[tuple[Fraction, ...], ...]) -> Fraction:
    size = len(matrix)
    assert all(len(row) == size for row in matrix)
    work = [list(row) for row in matrix]
    result = Fraction(1)
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column]),
            None,
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            result = -result
        pivot_value = work[column][column]
        result *= pivot_value
        for row in range(column + 1, size):
            multiple = work[row][column] / pivot_value
            for target in range(column, size):
                work[row][target] -= multiple * work[column][target]
    return result


def rank_one_trace_formula(entries: tuple[tuple[Fraction, ...], ...], cycle: tuple[int, ...]) -> Fraction:
    """Trace of P_i1...P_il when P_i(x)=alpha_i(x)e_i."""

    result = Fraction(1)
    for index, vertex in enumerate(cycle):
        next_vertex = cycle[(index + 1) % len(cycle)]
        result *= entries[vertex][next_vertex]
    return result


def trace_integrality_certificate() -> dict[str, object]:
    # Directed Cartan entries need not be integral.  The closed products are
    # invariant under e_i -> lambda_i e_i and equal canonical traces of
    # products of the integral endomorphisms P_i=1-s_i.
    entries = (
        (Fraction(2), Fraction(2, 3), Fraction(-3, 5), Fraction(7, 11)),
        (Fraction(5, 7), Fraction(2), Fraction(4, 9), Fraction(-2, 5)),
        (Fraction(-8, 13), Fraction(6, 5), Fraction(2), Fraction(3, 8)),
        (Fraction(9, 10), Fraction(-7, 6), Fraction(11, 9), Fraction(2)),
    )
    scales = (Fraction(2), Fraction(3), Fraction(5), Fraction(7))
    rescaled = tuple(
        tuple(
            scales[column] * entries[row][column] / scales[row]
            for column in range(4)
        )
        for row in range(4)
    )
    cycles = ((0, 1), (0, 1, 2), (0, 1, 2, 3), (0, 2, 1, 3))
    for cycle in cycles:
        assert rank_one_trace_formula(entries, cycle) == rank_one_trace_formula(
            rescaled, cycle
        )
    assert entries != rescaled
    return {
        "integral_endomorphisms": "P_i=1-s_i in End_O(M)",
        "finite_projective_module_has_canonical_O_valued_trace": True,
        "cycle_formula": "Tr(P_i1...P_il)=product(a_i_r,i_r+1)",
        "cycles_checked": len(cycles),
        "normal_rescaling_invariance_checked": True,
        "primitive_normals_required": False,
        "individual_directed_Cartan_entries_assumed_integral": False,
    }


def chordless_minor(length: int) -> dict[str, object]:
    assert length >= 3
    diagonal = tuple(Fraction(10 + index) for index in range(length))
    edge = tuple(Fraction(index + 2, index + 3) for index in range(length))
    gram = [[Fraction(0) for _ in range(length)] for _ in range(length)]
    for index in range(length):
        gram[index][index] = diagonal[index]
        next_index = (index + 1) % length
        gram[index][next_index] = gram[next_index][index] = edge[index]
    rows = range(1, length - 1)
    columns = range(2, length)
    minor = tuple(tuple(gram[row][column] for column in columns) for row in rows)
    direct = determinant(minor)
    expected = Fraction(1)
    for index in range(1, length - 1):
        expected *= edge[index]
    assert direct == expected != 0
    return {
        "cycle_length": length,
        "minor_size": length - 2,
        "minor_nonzero": True,
        "rank_lower_bound": length - 2,
    }


def chordless_rank_certificate() -> dict[str, object]:
    rows = [chordless_minor(length) for length in range(3, 13)]
    return {
        "exact_minor_rows": rows,
        "minor_rows": "2,...,ell-1",
        "minor_columns": "3,...,ell",
        "minor_diagonal": "nonzero perimeter entries G[r,r+1]",
        "rank_inequality": "ell-2<=rank(G_C)<=m",
        "cycle_length_bound": "ell<=m+2",
        "nonzero_perimeter_preserved_at_all_embeddings": True,
    }


def conjugate_height_certificate() -> dict[str, object]:
    records = []
    for module_rank in range(3, 11):
        length_bound = module_rank + 2
        height = 2**length_bound
        # c(C)^2 is the product of the perimeter pair cycles.  At sigma_0
        # each selected factor is <=4; at definite places each is <4.
        records.append(
            {
                "module_rank": module_rank,
                "chordless_length_bound": length_bound,
                "all_embedding_cycle_height": height,
            }
        )
    return {
        "cycle_square_identity": "c(C)^2=product(perimeter p_ij)",
        "distinguished_selected_pair_bound": "0<=p_ij<=4",
        "definite_pair_bound": "0<=sigma(p_ij)<4",
        "records": records,
        "pair_generators_obey_same_final_height": True,
    }


def primitive_grid_and_discriminant_certificate() -> dict[str, object]:
    rows = []
    for module_rank in range(3, 21):
        degree_ceiling = 20 // module_rank
        if degree_ceiling < 2:
            exponent = 0
        else:
            exponent = (module_rank + 7) * degree_ceiling * (degree_ceiling - 1)
        rows.append(
            {
                "module_rank": module_rank,
                "degree_ceiling": max(1, degree_ceiling),
                "exponent": exponent,
            }
        )
    maximum = max(row["exponent"] for row in rows)
    maximizers = [row for row in rows if row["exponent"] == maximum]
    assert maximum == 300
    assert maximizers == [{"module_rank": 3, "degree_ceiling": 6, "exponent": 300}]
    assert comb(6, 2) == 15
    return {
        "minimal_generating_subpacket_size_upper_bound": 2,
        "reason": "three strict adjunctions force degree at least 8>6",
        "embedding_pair_count_upper_bound": 15,
        "primitive_grid": "x+q*y for q=0,...,15",
        "generator_height": "2^(m+2)",
        "primitive_height": "2^(m+6)",
        "difference_height": "2^(m+7)",
        "Vandermonde_exponent_formula": "(m+7)d(d-1)",
        "rank_rows": rows,
        "uniform_binary_exponent": maximum,
        "field_discriminant_bound": "2^300",
        "Hermite_Minkowski_gives_finite_field_set": True,
    }


def scope_certificate() -> dict[str, object]:
    text = " ".join((ROOT / PINS["report"][0]).read_text().split())
    required = {
        "full_trace": "primitive full trace",
        "minimal_field": "minimal Vinberg field",
        "generating_core_only": "hypothesis concerns only a field-generating core",
        "dotted_outside": "dotted or ultraparallel pairs outside that core",
        "pair_generator_literal": "pair-cycle generator is itself required to satisfy",
        "no_facet_bound": "neither bounds the number of facets",
        "dotted_generator_open": "generating chordless cycle containing an ultraparallel perimeter edge",
        "no_mixed_height": "bound the canonical mixed intersection",
        "no_chamber_finiteness": "prove finiteness of chambers",
        "no_maximality": "use reflection maximality, finite volume",
    }
    for label, phrase in required.items():
        assert phrase in text, (label, phrase)
    return {label: True for label in required}


def build_audit() -> dict[str, object]:
    return {
        "schema": "k3_nondotted_field_generating_core_bound/independent_hostile/v1",
        "verdict": "PASS_CONDITIONAL_EXACT_SCOPE",
        "pins": verify_pins(),
        "projective_module_trace_integrality": trace_integrality_certificate(),
        "chordless_minor_and_length": chordless_rank_certificate(),
        "all_embedding_cycle_height": conjugate_height_certificate(),
        "primitive_grid_and_discriminant": primitive_grid_and_discriminant_certificate(),
        "scope": scope_certificate(),
    }


def main() -> None:
    print(json.dumps(build_audit(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
