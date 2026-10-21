"""Independent hostile checks for PB-Atlas numerical Stage A.

The level traversal below deliberately tests every remaining wall at every
orbit representative; it does not use the core's stabilizer-orbit shortcut.
The endpoint checks consume only frozen arrays and emitted programs and
reconstruct every lattice identity independently.
"""

from functools import lru_cache
import copy
import json
import math
from pathlib import Path

import pytest
import sympy as sp

from scripts.pb_atlas_contraction_core import (
    enumerate_frame_orbits,
    parse_physical_branch,
    validate_automorphisms,
    wall_gram_automorphisms,
)


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "data/fixtures/pb_atlas_numerical_stage_a_inputs.json"
ARCHIVE = ROOT / "data/enumerations/pb_atlas_numerical_contraction_stage_a.json"


@lru_cache(maxsize=1)
def _fixtures():
    return json.loads(FIXTURES.read_text(encoding="utf-8"))


@lru_cache(maxsize=1)
def _archive():
    return json.loads(ARCHIVE.read_text(encoding="utf-8"))


def _mask_indices(mask):
    result = []
    while mask:
        bit = mask & -mask
        result.append(bit.bit_length() - 1)
        mask -= bit
    return tuple(result)


def _transform_mask(mask, mapped_bits):
    result = 0
    while mask:
        bit = mask & -mask
        result |= mapped_bits[bit.bit_length() - 1]
        mask -= bit
    return result


def _independent_all_candidate_quotient(payload, permutations):
    """Traverse every candidate, quotienting children only after testing."""

    branch = parse_physical_branch(payload)
    group = [tuple(map(int, permutation)) for permutation in permutations]
    identity = tuple(range(branch.wall_count))
    assert identity in group
    mapped_bits = [tuple(1 << permutation[i] for i in identity) for permutation in group]
    canonical_cache = {0: 0}
    orbit_sizes = {0: 1}

    def canonical(mask):
        if mask in canonical_cache:
            return canonical_cache[mask]
        members = {
            _transform_mask(mask, permutation_bits)
            for permutation_bits in mapped_bits
        }
        representative = min(members)
        for member in members:
            canonical_cache[member] = representative
        orbit_sizes[representative] = len(members)
        return representative

    wall_gram = branch.wall_gram
    canonical_pairings = branch.canonical_pairings
    current = {0}
    orbit_counts = [1]
    labelled_counts = [1]
    checked_eligible_transitions = 0
    for depth in range(1, branch.rank):
        following = set()
        for mask in current:
            indices = _mask_indices(mask)
            if indices:
                inverse = wall_gram.extract(indices, indices).inv()
                assert all(entry.q == 1 for entry in inverse)
            else:
                inverse = sp.zeros(0)
            for candidate in branch.allowed_indices:
                if mask & (1 << candidate):
                    continue
                vector = sp.Matrix(
                    [wall_gram[index, candidate] for index in indices]
                )
                product = inverse * vector
                current_square = wall_gram[candidate, candidate] - (
                    vector.T * product if indices else sp.zeros(1)
                )[0]
                current_K_pairing = canonical_pairings[candidate] - sum(
                    canonical_pairings[index] * product[position]
                    for position, index in enumerate(indices)
                )
                if current_square != -1 or current_K_pairing != -1:
                    continue

                # Independent determinant/negative-definiteness audit of the
                # Schur transition, without the core's inverse-update helper.
                child_indices = (*indices, candidate)
                child_gram = wall_gram.extract(child_indices, child_indices)
                assert child_gram.det() == (-1) ** depth
                assert all(
                    value > 0
                    for value in (
                        child_gram[:size, :size].det()
                        * (-1) ** size
                        for size in range(1, depth + 1)
                    )
                )
                checked_eligible_transitions += 1
                following.add(canonical(mask | (1 << candidate)))
        current = following
        orbit_counts.append(len(current))
        labelled_counts.append(sum(orbit_sizes[mask] for mask in current))
    return orbit_counts, labelled_counts, checked_eligible_transitions


def _compose(left, right):
    return tuple(left[right[index]] for index in range(len(left)))


def _inverse_permutation(permutation):
    inverse = [0] * len(permutation)
    for source, target in enumerate(permutation):
        inverse[target] = source
    return tuple(inverse)


def _assert_group_and_marked_action(payload, permutations):
    branch = parse_physical_branch(payload)
    group = {tuple(map(int, row)) for row in permutations}
    identity = tuple(range(branch.wall_count))
    assert identity in group
    assert len(group) == len(permutations)
    # Inverses plus closure under a greedily found generating set provide an
    # independent group audit without a quadratic multiplication table.
    assert all(_inverse_permutation(row) in group for row in group)
    generated = {identity}
    generators = []
    for candidate in sorted(group):
        if candidate in generated:
            continue
        generators.append(candidate)
        generated = {identity}
        frontier = [identity]
        while frontier:
            current = frontier.pop()
            for generator in generators:
                product = _compose(generator, current)
                assert product in group
                if product not in generated:
                    generated.add(product)
                    frontier.append(product)
    assert generated == group

    boundary_weights = branch.boundary_weight_by_index
    for permutation in group:
        assert tuple(sorted(permutation)) == identity
        for left in range(branch.wall_count):
            assert boundary_weights[left] == boundary_weights[permutation[left]]
            assert branch.canonical_pairings[left] == branch.canonical_pairings[
                permutation[left]
            ]
            for right in range(branch.wall_count):
                assert branch.wall_gram[left, right] == branch.wall_gram[
                    permutation[left], permutation[right]
                ]


def _program_matrix_rows_as_columns(rows):
    return sp.Matrix(rows).T


def _primitive_pair(values):
    divisor = math.gcd(*(abs(int(value)) for value in values))
    assert divisor
    return tuple(int(value) // divisor for value in values)


def _independent_nef_rays(model_gram, endpoint_rows):
    inequalities = [
        tuple(int(value) for value in (sp.Matrix(row).T * model_gram))
        for row in endpoint_rows
    ]
    rays = set()
    for left, right in inequalities:
        if left == right == 0:
            continue
        for candidate in ((right, -left), (-right, left)):
            if all(a * candidate[0] + b * candidate[1] >= 0 for a, b in inequalities):
                rays.add(_primitive_pair(candidate))
    return rays


def _verify_endpoint_program(payload, record):
    program = record["program"]
    g = sp.Matrix(payload["picard_gram"])
    walls = [sp.Matrix(row) for row in payload["wall_vectors"]]
    wall_matrix = sp.Matrix.hstack(*walls)
    canonical = sp.Matrix(payload["canonical_class"])
    selected = tuple(program["contracted_wall_indices_zero_based"])
    order = tuple(program["sequential_contraction_order_zero_based"])
    assert set(selected) == set(order)

    strict_matrix = sp.Matrix.hstack(*(walls[index] for index in selected))
    strict_gram = strict_matrix.T * g * strict_matrix
    assert strict_gram == sp.Matrix(program["strict_contracted_wall_gram"])

    independently_recovered_exceptionals = []
    prior_strict = []
    for wall_index in order:
        curve = walls[wall_index]
        if prior_strict:
            prior = sp.Matrix.hstack(*prior_strict)
            prior_gram = prior.T * g * prior
            exceptional = curve - prior * prior_gram.inv() * prior.T * g * curve
        else:
            exceptional = curve
        assert all(entry.q == 1 for entry in exceptional)
        assert (exceptional.T * g * exceptional)[0] == -1
        assert (canonical.T * g * exceptional)[0] == -1
        assert all((old.T * g * exceptional)[0] == 0 for old in independently_recovered_exceptionals)
        independently_recovered_exceptionals.append(exceptional)
        prior_strict.append(curve)
    assert [list(map(int, vector)) for vector in independently_recovered_exceptionals] == program[
        "total_exceptional_basis_in_input_picard_coordinates"
    ]

    model = _program_matrix_rows_as_columns(
        program["model_basis_in_input_picard_coordinates"]
    )
    basis = sp.Matrix.hstack(model, *independently_recovered_exceptionals)
    assert basis.shape == (g.rows, g.rows)
    assert abs(int(basis.det())) == 1
    actual_basis_gram = basis.T * g * basis
    assert actual_basis_gram == sp.Matrix(program["endpoint_blowup_basis_gram"])
    coordinates = basis.inv() * wall_matrix
    assert coordinates.T == sp.Matrix(
        program["all_wall_coordinates_in_endpoint_blowup_basis"]
    )
    model_rank = model.cols
    endpoint_rows = [
        row[:model_rank]
        for row in program["all_wall_coordinates_in_endpoint_blowup_basis"]
    ]
    assert endpoint_rows == program["all_wall_endpoint_image_coordinates"]
    assert all(endpoint_rows[index] == [0] * model_rank for index in selected)
    assert list(map(int, basis.inv() * canonical)) == program[
        "canonical_coordinates_in_endpoint_blowup_basis"
    ]

    if program["endpoint_kind"] == "P2":
        assert program["hirzebruch_e"] is None
        assert model.T * g * model == sp.Matrix([[1]])
        assert program["canonical_coordinates_in_endpoint_blowup_basis"][0] == -3
        assert all(row[0] >= 0 for row in endpoint_rows)
        assert program["nef_rays_in_model_basis"] == [[1]]
    else:
        e = program["hirzebruch_e"]
        model_gram = sp.Matrix([[-e, 1], [1, 0]])
        assert model.T * g * model == model_gram
        assert program["canonical_coordinates_in_endpoint_blowup_basis"][:2] == [
            -2,
            -(e + 2),
        ]
        emitted_rays = {tuple(row) for row in program["nef_rays_in_model_basis"]}
        assert emitted_rays == {(0, 1), (1, e)}
        assert emitted_rays == _independent_nef_rays(model_gram, endpoint_rows)

    expected_boundary = []
    for index, weight in zip(
        payload["boundary_indices_zero_based"], payload["boundary_weights"]
    ):
        expected_boundary.append(
            {
                "wall_index_zero_based": index,
                "anti_bicanonical_weight": weight,
                "endpoint_blowup_coordinates": list(map(int, coordinates[:, index])),
                "endpoint_image_coordinates": list(map(int, coordinates[:model_rank, index])),
            }
        )
    assert program["boundary"] == expected_boundary


def test_all_four_level_censuses_survive_an_independent_no_stabilizer_traversal():
    for fixture, frozen in zip(_fixtures()["cases"], _archive()["cases"]):
        permutations = frozen["certificate"]["symmetry_quotient"][
            "automorphism_permutations_zero_based"
        ]
        _assert_group_and_marked_action(fixture["input"], permutations)
        orbit_counts, labelled_counts, transitions = _independent_all_candidate_quotient(
            fixture["input"], permutations
        )
        census = frozen["certificate"]["contraction_level_census"]
        assert orbit_counts == [row["orbit_count"] for row in census]
        assert labelled_counts == [row["labelled_frame_count"] for row in census]
        assert labelled_counts == fixture["expected"]["labelled_frame_counts_by_depth"]
        assert transitions > 0


def test_every_frozen_terminal_orbit_and_endpoint_packet_reconstructs_exactly():
    for fixture, frozen in zip(_fixtures()["cases"], _archive()["cases"]):
        payload = fixture["input"]
        certificate = frozen["certificate"]
        permutations = [
            tuple(row)
            for row in certificate["symmetry_quotient"][
                "automorphism_permutations_zero_based"
            ]
        ]
        histogram = {}
        for family in ("hirzebruch_endpoint_orbits", "P2_endpoint_orbits"):
            for record in certificate[family]:
                representative = int(record["representative_mask"], 16)
                independent_members = {
                    _transform_mask(
                        representative,
                        tuple(1 << permutation[i] for i in range(len(permutation))),
                    )
                    for permutation in permutations
                }
                assert independent_members == {
                    int(encoded, 16) for encoded in record["member_masks"]
                }
                assert len(independent_members) == record["orbit_size"]
                _verify_endpoint_program(payload, record)
                program = record["program"]
                key = (
                    "P2"
                    if program["endpoint_kind"] == "P2"
                    else f"F_{program['hirzebruch_e']}"
                )
                histogram[key] = histogram.get(key, 0) + record["orbit_size"]
        assert histogram.get("P2", 0) == certificate["P2_labelled_count"]
        assert {
            key: count for key, count in histogram.items() if key != "P2"
        } == certificate["hirzebruch_labelled_histogram"]


def test_nonsaturated_physical_subset_has_same_labelled_census_with_or_without_symmetry():
    payload = copy.deepcopy(_fixtures()["cases"][0]["input"])
    subset = (0, 4, 5, 6, 10, 18, 19, 22, 26, 27, 28)
    payload["wall_vectors"] = [payload["wall_vectors"][index] for index in subset]
    payload["boundary_indices_zero_based"] = [0]
    branch = parse_physical_branch(payload)
    assert branch.wall_sublattice_index == 2
    safe_group = wall_gram_automorphisms(branch)
    identity = (tuple(range(branch.wall_count)),)
    safe = enumerate_frame_orbits(branch, automorphisms=safe_group)
    unquotiented = enumerate_frame_orbits(branch, automorphisms=identity)
    assert [row["labelled_frame_count"] for row in safe["level_census"]] == [
        row["labelled_frame_count"] for row in unquotiented["level_census"]
    ]


def test_archived_nonsaturated_swap_really_has_a_half_integral_lift():
    hostile = _archive()["nonsaturated_symmetry_hostile_fixture"]
    gram = sp.Matrix(hostile["synthetic_picard_gram"])
    wall_matrix = sp.Matrix(hostile["synthetic_wall_vectors"]).T
    wall_gram = wall_matrix.T * gram * wall_matrix
    assert wall_gram == sp.Matrix(hostile["synthetic_wall_gram"])
    permutation = hostile["unsafe_marked_Gram_swap"]
    assert all(
        wall_gram[left, right] == wall_gram[permutation[left], permutation[right]]
        for left in range(wall_matrix.cols)
        for right in range(wall_matrix.cols)
    )
    pivots = tuple(int(index) for index in wall_matrix.rref()[1])
    source = wall_matrix[:, pivots]
    target = sp.Matrix.hstack(*(wall_matrix[:, permutation[index]] for index in pivots))
    lift = target * source.inv()
    assert any(entry.q == 2 for entry in lift)


def test_unclosed_or_duplicate_automorphism_lists_are_rejected_before_orbit_use():
    fixture = _fixtures()["cases"][0]
    branch = parse_physical_branch(fixture["input"])
    full_group = [
        tuple(row)
        for row in _archive()["cases"][0]["certificate"]["symmetry_quotient"][
            "automorphism_permutations_zero_based"
        ]
    ]
    identity = tuple(range(branch.wall_count))
    element = next(row for row in full_group if _compose(row, row) != identity)
    with pytest.raises(ValueError, match="not closed under composition"):
        validate_automorphisms(branch, (identity, element))
    with pytest.raises(ValueError, match="contains duplicates"):
        validate_automorphisms(branch, (identity, identity))
