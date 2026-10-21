"""Independent hostile checks for the raw adjunction/log Stage-0 compiler.

These tests deliberately avoid the compiler's minimum-basis and HNF
membership helpers.  They recompute the finite adjunction search from
alternative wall bases and the adapted-module gate from its finite residue
group modulo the original Neron--Severi lattice.
"""

from fractions import Fraction
import itertools
import json
import math
import random

import sympy as sp

from scripts.certify_raw_adjunction_log_stage0 import (
    FIXTURES,
    solve_raw_input,
)


def _fixtures():
    return json.loads(FIXTURES.read_text(encoding="utf-8"))["cases"]


def _dot(g, left, right):
    return (left.T * g * right)[0]


def _local_rows(q):
    rows = ((1, q - 2),)
    return rows + (((2, 2 * q - 1),) if q % 2 else ())


def _alternative_minimum_weight_bases(payload, count=8):
    """Choose several minimum-ambiguity bases with independent tie orders."""

    g = sp.Matrix(payload["picard_gram"])
    walls = [sp.Matrix(row) for row in payload["wall_vectors"]]
    q = [-int(_dot(g, wall, wall)) for wall in walls]
    rank = g.rows
    bases = set()
    for seed in range(count * 4):
        order = list(range(len(walls)))
        random.Random(0x5A17E0 + seed).shuffle(order)
        order.sort(key=lambda index: q[index] % 2)
        chosen = []
        columns = sp.zeros(rank, 0)
        for index in order:
            trial = columns.row_join(walls[index])
            if trial.rank() > columns.rank():
                chosen.append(index)
                columns = trial
            if len(chosen) == rank:
                break
        assert len(chosen) == rank
        bases.add(tuple(chosen))
        if len(bases) == count:
            break
    return sorted(bases)


def _enumerate_from_specified_basis(payload, basis):
    """Reimplement the adjunction enumeration without compiler helpers."""

    g = sp.Matrix(payload["picard_gram"])
    walls = [sp.Matrix(row) for row in payload["wall_vectors"]]
    q = [-int(_dot(g, wall, wall)) for wall in walls]
    covectors = sp.Matrix([list(walls[index].T * g) for index in basis])
    assert covectors.det() != 0
    results = set()
    for choices in itertools.product(*[_local_rows(q[index]) for index in basis]):
        canonical = covectors.inv() * sp.Matrix([pairing for _, pairing in choices])
        if any(sp.denom(entry) != 1 for entry in canonical):
            continue
        scales = []
        for wall_q, wall in zip(q, walls):
            pairing = int(_dot(g, canonical, wall))
            matches = [scale for scale, target in _local_rows(wall_q) if target == pairing]
            if len(matches) != 1:
                break
            scales.append(matches[0])
        else:
            # Characteristic means K.e_i == e_i^2 (mod 2) on the coordinate
            # basis; no compiler parity helper is used here.
            if any(
                (int((sp.eye(g.rows)[:, i].T * g * canonical)[0]) - int(g[i, i])) % 2
                for i in range(g.rows)
            ):
                continue
            if _dot(g, canonical, canonical) != 10 - g.rows:
                continue
            results.add(
                (
                    tuple(int(entry) for entry in canonical),
                    tuple(scales),
                )
            )
    return results


def _compiler_candidate_set(certificate):
    return {
        (tuple(candidate["canonical_class"]), tuple(candidate["physical_scales"]))
        for candidate in certificate["canonical_candidates"]
    }


def _candidate_signature(certificate):
    return {
        (
            tuple(candidate["canonical_class"]),
            tuple(candidate["physical_scales"]),
            candidate["branch"],
            candidate["accepted"],
        )
        for candidate in certificate["canonical_candidates"]
    }


def _reindex_signature(certificate, permutation):
    rows = set()
    for candidate in certificate["canonical_candidates"]:
        old_scales = [None] * len(permutation)
        for new_index, old_index in enumerate(permutation):
            old_scales[old_index] = candidate["physical_scales"][new_index]
        rows.add(
            (
                tuple(candidate["canonical_class"]),
                tuple(old_scales),
                candidate["branch"],
                candidate["accepted"],
            )
        )
    return rows


def _fraction(value):
    if isinstance(value, str):
        return Fraction(value)
    return Fraction(value)


def _fractional_part(vector):
    return tuple(_fraction(entry) % 1 for entry in vector)


def _residue_module(g, physical, boundary, ramification):
    """Return M/L as residues, without HNF or compiler module code."""

    generators = [
        [Fraction(int(entry), ramification[position]) for entry in physical[index]]
        for position, index in enumerate(boundary)
    ]
    residues = set()
    for coefficients in itertools.product(*(range(m) for m in ramification)):
        value = [Fraction(0) for _ in range(g.rows)]
        for coefficient, generator in zip(coefficients, generators):
            for i, entry in enumerate(generator):
                value[i] += coefficient * entry
        residues.add(_fractional_part(value))
    return generators, residues


def _reflect_fraction_vector(g, root, vector):
    root = [Fraction(int(entry)) for entry in root]
    pairing = sum(
        vector[i] * int(g[i, j]) * root[j]
        for i in range(g.rows)
        for j in range(g.cols)
    )
    square = sum(
        root[i] * int(g[i, j]) * root[j]
        for i in range(g.rows)
        for j in range(g.cols)
    )
    return [vector[i] - Fraction(2) * pairing * root[i] / square for i in range(g.rows)]


def _independent_module_checks(payload, candidate):
    g = sp.Matrix(payload["picard_gram"])
    walls = payload["wall_vectors"]
    physical = [
        [scale * entry for entry in wall]
        for scale, wall in zip(candidate["physical_scales"], walls)
    ]
    boundary = candidate["boundary_indices_zero_based_inferred"]
    ramification = candidate["ramification_indices_inferred"]
    boundary_generators, residues = _residue_module(
        g, physical, boundary, ramification
    )
    lattice_generators = [
        [Fraction(int(i == j)) for i in range(g.rows)] for j in range(g.rows)
    ] + boundary_generators
    all_preserve = True
    ns_integral = True
    for root in walls:
        for generator in lattice_generators:
            image = _reflect_fraction_vector(g, root, generator)
            all_preserve &= _fractional_part(image) in residues
        for generator in lattice_generators[: g.rows]:
            image = _reflect_fraction_vector(g, root, generator)
            ns_integral &= all(entry.denominator == 1 for entry in image)

    ell = math.lcm(*ramification) if ramification else 1
    scaled_form_integral = True
    for left in lattice_generators:
        for right in lattice_generators:
            pairing = sum(
                left[i] * int(g[i, j]) * right[j]
                for i in range(g.rows)
                for j in range(g.cols)
            )
            scaled_form_integral &= (ell * pairing).denominator == 1
    return all_preserve, ns_integral, scaled_form_integral


def test_every_frozen_case_is_exhaustive_from_alternative_wall_bases():
    for case in _fixtures():
        payload = case["raw_input"]
        certificate = solve_raw_input(payload)
        expected = _compiler_candidate_set(certificate)
        bases = _alternative_minimum_weight_bases(payload)
        assert bases
        for basis in bases:
            assert _enumerate_from_specified_basis(payload, basis) == expected


def test_full_candidate_data_are_invariant_under_hostile_wall_permutations():
    for case_number, case in enumerate(_fixtures()):
        payload = case["raw_input"]
        expected = _candidate_signature(solve_raw_input(payload))
        wall_count = len(payload["wall_vectors"])
        permutations = [
            list(reversed(range(wall_count))),
            list(range(1, wall_count)) + [0],
            random.Random(0xBAD5EED + case_number).sample(range(wall_count), wall_count),
        ]
        for permutation in permutations:
            permuted = {
                "schema": payload["schema"],
                "picard_gram": payload["picard_gram"],
                "wall_vectors": [payload["wall_vectors"][i] for i in permutation],
            }
            assert _reindex_signature(solve_raw_input(permuted), permutation) == expected


def test_half_lattice_reflection_gate_matches_direct_residue_group_calculation():
    for case in _fixtures():
        certificate = solve_raw_input(case["raw_input"])
        candidates = list(certificate["accepted_candidates"])
        if case["case_id"] == "blowup_F4_maximal_triangle_negative_control":
            candidates = certificate["canonical_candidates"]
        for candidate in candidates:
            direct = _independent_module_checks(case["raw_input"], candidate)
            module = candidate["reflection_module"]
            assert direct == (
                module["all_wall_reflections_preserve_inferred_adapted_module"],
                candidate["all_reflections_integral_on_NS"],
                module["ell_scaled_form_is_integral"],
            )


def test_primitive_and_twice_primitive_minus_four_modules_do_not_collapse():
    by_id = {case["case_id"]: case for case in _fixtures()}
    h10 = solve_raw_input(by_id["coble_type_ii_h10_30"]["raw_input"])[
        "accepted_candidates"
    ][0]
    assert h10["physical_scales"][h10["boundary_indices_zero_based_inferred"][0]] == 2
    assert h10["reflection_module"]["common_coordinate_denominator"] == 1
    assert h10["all_reflections_integral_on_NS"]

    for case_id in (
        "coble_type_i_h11_26",
        "coble_type_vi_char3_h14_60",
        "coble_type_mii_char3_h17_48",
    ):
        candidate = solve_raw_input(by_id[case_id]["raw_input"])["accepted_candidates"][0]
        assert all(
            candidate["physical_scales"][index] == 1
            for index in candidate["boundary_indices_zero_based_inferred"]
        )
        assert candidate["reflection_module"]["common_coordinate_denominator"] == 2
        assert not candidate["all_reflections_integral_on_NS"]
        assert candidate["reflection_module"][
            "all_wall_reflections_preserve_inferred_adapted_module"
        ]


def test_F4_canonical_class_is_unique_by_coordinate_box_brute_force():
    case = next(
        case
        for case in _fixtures()
        if case["case_id"] == "blowup_F4_maximal_triangle_negative_control"
    )
    payload = case["raw_input"]
    g = sp.Matrix(payload["picard_gram"])
    walls = [sp.Matrix(row) for row in payload["wall_vectors"]]
    q = [-int(_dot(g, wall, wall)) for wall in walls]
    # The three walls are a Q-basis.  Solving all four local pairing rows
    # first proves that the displayed box contains every possible K, so the
    # subsequent coordinate search is exhaustive rather than heuristic.
    covectors = sp.Matrix([list(wall.T * g) for wall in walls])
    local_solutions = [
        covectors.inv() * sp.Matrix(pairings)
        for pairings in itertools.product(
            *[[pairing for _, pairing in _local_rows(wall_q)] for wall_q in q]
        )
    ]
    assert all(abs(entry) <= 12 for canonical in local_solutions for entry in canonical)
    found = []
    for coordinates in itertools.product(range(-12, 13), repeat=3):
        canonical = sp.Matrix(coordinates)
        if _dot(g, canonical, canonical) != 7:
            continue
        if any(
            (int((sp.eye(3)[:, i].T * g * canonical)[0]) - int(g[i, i])) % 2
            for i in range(3)
        ):
            continue
        if all(
            int(_dot(g, canonical, wall)) in {pairing for _, pairing in _local_rows(wall_q)}
            for wall_q, wall in zip(q, walls)
        ):
            found.append(coordinates)
    assert found == [(-2, -6, 1)]


def test_all_six_calibration_counts_and_alphabets_are_recomputed_exactly():
    expected = {
        "coble_type_ii_h10_30": (11, 30, 2, 2, 1, {"-4": 1, "-2": 12, "-1": 17}),
        "coble_type_i_h11_26": (12, 26, 1, 1, 1, {"-4": 2, "-2": 10, "-1": 14}),
        "coble_type_vi_char3_h14_60": (15, 60, 1, 1, 1, {"-4": 5, "-2": 10, "-1": 45}),
        "coble_type_mii_char3_h17_48": (18, 48, 1, 1, 1, {"-4": 8, "-2": 24, "-1": 16}),
        "toric_h6_9": (7, 9, 2, 2, 1, {"-2": 6, "-1": 3}),
        "blowup_F4_maximal_triangle_negative_control": (3, 3, 4, 1, 0, None),
    }
    for case in _fixtures():
        certificate = solve_raw_input(case["raw_input"])
        rank, walls, assignments, candidates, accepted, histogram = expected[
            case["case_id"]
        ]
        assert certificate["validation"]["rank"] == rank
        assert certificate["validation"]["wall_count"] == walls
        assert certificate["finite_enumeration"]["basis_pairing_assignments_checked"] == assignments
        assert len(certificate["canonical_candidates"]) == candidates
        assert certificate["accepted_candidate_count"] == accepted
        if histogram is not None:
            assert certificate["accepted_candidates"][0]["physical_square_histogram"] == histogram
        else:
            rejected = certificate["canonical_candidates"][0]
            assert rejected["log_residual_K_plus_Delta"] == ["-3/2", -6, 1]
            assert not rejected["log_equation"]


def test_orientation_primitivity_and_ray_distinctness_are_real_input_gates():
    case = next(
        case
        for case in _fixtures()
        if case["case_id"] == "blowup_F4_maximal_triangle_negative_control"
    )
    raw = case["raw_input"]

    flipped = json.loads(json.dumps(raw))
    flipped["wall_vectors"][0] = [-entry for entry in flipped["wall_vectors"][0]]
    assert solve_raw_input(flipped)["status"] == "INVALID_INPUT"

    nonprimitive = json.loads(json.dumps(raw))
    nonprimitive["wall_vectors"][0] = [2 * entry for entry in nonprimitive["wall_vectors"][0]]
    invalid = solve_raw_input(nonprimitive)
    assert invalid["status"] == "INVALID_INPUT"
    assert not invalid["validation"]["every_wall_primitive"]

    opposite_duplicate = json.loads(json.dumps(raw))
    opposite_duplicate["wall_vectors"].append(
        [-entry for entry in opposite_duplicate["wall_vectors"][0]]
    )
    invalid = solve_raw_input(opposite_duplicate)
    assert invalid["status"] == "INVALID_INPUT"
    assert not invalid["validation"]["distinct_unoriented_wall_rays"]
