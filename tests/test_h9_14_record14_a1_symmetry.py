from itertools import product
import json
from pathlib import Path

import z3
from ortools.sat.python import cp_model
import pytest

from scripts.audit_h9_14_record14_ridge_capacity import audit as capacity_audit
from scripts.audit_h9_14_record14_ridge_capacity import require_capacity
from scripts.independent_h9_14_record14_a0_audit import (
    _build_a0_data,
    _source_actions as _independent_a0_actions,
)
from scripts.independent_h9_14_record14_a1_audit import (
    _add_involution_lex_leader,
    _build_a1_data,
    _source_actions,
)
from scripts.explore_h9_14_record14_a0_symmetry_z3 import (
    build_a0_problem,
    source_generators,
)
from scripts.explore_h9_14_record14_a1_symmetry_z3 import (
    _add_lex_leader,
    _move,
    stabilizer_generators,
)
from scripts.explore_h9_14_record14_fixed_orbit_z3 import build_problem


ROOT = Path(__file__).resolve().parents[1]


def test_involution_pair_reduction_is_exact_lex_order() -> None:
    permutation = (4, 3, 2, 1, 0, 5)
    context = z3.Context()
    values = [z3.Bool(f"x_{index}", ctx=context) for index in range(6)]
    solver = z3.Solver(ctx=context)
    pair_count = _add_lex_leader(solver, values, permutation, "test", context)
    assert pair_count == 2
    for bits in product((False, True), repeat=6):
        solver.push()
        solver.add([value == bit for value, bit in zip(values, bits, strict=True)])
        encoded = solver.check() == z3.sat
        literal = bits <= tuple(bits[permutation[index]] for index in range(6))
        solver.pop()
        assert encoded == literal


def test_independent_cp_sat_involution_encoding_is_exact() -> None:
    permutation = (4, 3, 2, 1, 0, 5)
    model = cp_model.CpModel()
    values = [model.new_bool_var(f"x_{index}") for index in range(6)]
    assert _add_involution_lex_leader(model, values, permutation, "test") == 2

    class Collector(cp_model.CpSolverSolutionCallback):
        def __init__(self) -> None:
            super().__init__()
            self.answers: set[tuple[bool, ...]] = set()

        def on_solution_callback(self) -> None:
            self.answers.add(tuple(bool(self.value(value)) for value in values))

    collector = Collector()
    solver = cp_model.CpSolver()
    solver.parameters.enumerate_all_solutions = True
    solver.parameters.num_search_workers = 1
    assert solver.solve(model, collector) == cp_model.OPTIMAL
    expected = {
        bits
        for bits in product((False, True), repeat=6)
        if bits <= tuple(bits[permutation[index]] for index in range(6))
    }
    assert collector.answers == expected


def test_all_a1_generator_actions_preserve_candidate_sets() -> None:
    for orbit_index in range(5):
        problem = build_problem(1, orbit_index)
        independent = _build_a1_data(orbit_index)
        assert independent["fours"] == problem["fours"]
        assert independent["forced"] == problem["shadow"]
        assert independent["addition_count"] == problem["addition_count"]
        assert independent["lower"] == problem["lower"]
        assert independent["upper"] == problem["upper"]
        assert independent["candidates"] == problem["candidates"]
        assert independent["base"] == problem["base"]
        assert independent["small"] == problem["small"]
        assert independent["facets"] == problem["fives"]
        four = problem["fours"][0]
        assert _source_actions(four) == stabilizer_generators(four)
        candidates = tuple(problem["candidates"])
        candidate_set = set(candidates)
        candidate_index = {face: index for index, face in enumerate(candidates)}
        for action in stabilizer_generators(four):
            assert _move(four, action) == four
            permutation = tuple(
                candidate_index[_move(face, action)] for face in candidates
            )
            assert sorted(permutation) == list(range(len(candidates)))
            assert all(
                permutation[permutation[index]] == index
                for index in range(len(candidates))
            )
            assert {_move(face, action) for face in candidates} == candidate_set


def test_a0_problem_and_full_source_generators() -> None:
    problem = build_a0_problem()
    independent = _build_a0_data()
    assert problem["addition_count"] == 8
    assert problem["scheduled_facets"] == 110
    assert len(problem["candidates"]) == 1002
    assert len(problem["small"]) == 101
    assert len(problem["fives"]) == 310
    assert problem["lower"] == (0,) * 14
    assert problem["upper"] == (4, 4, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 1, 0)
    assert independent["lower"] == problem["lower"]
    assert independent["upper"] == problem["upper"]
    assert independent["candidates"] == problem["candidates"]
    assert independent["base"] == problem["base"]
    assert independent["small"] == problem["small"]
    assert independent["facets"] == problem["fives"]
    candidates = tuple(problem["candidates"])
    candidate_index = {face: index for index, face in enumerate(candidates)}
    actions = source_generators()
    assert _independent_a0_actions() == actions
    assert len(actions) == 9
    for action in actions:
        permutation = tuple(
            candidate_index[_move(face, action)] for face in candidates
        )
        assert sorted(permutation) == list(range(len(candidates)))
        assert all(
            permutation[permutation[index]] == index
            for index in range(len(candidates))
        )


def test_record14_rejects_undersized_reused_ridge_bitset() -> None:
    require_capacity(784, 13)
    with pytest.raises(
        ValueError,
        match="candidate count 994 exceeds CandidateWord capacity 832",
    ):
        require_capacity(994, 13)
    require_capacity(994, 16)
    document = capacity_audit()
    assert document["shared_record18_capacity"] == 832
    assert document["isolated_record14_capacity"] == 1024
    assert [row["rejected_by_13_word_backend"] for row in document["rows"]] == [
        True,
        True,
        False,
        True,
        True,
    ]
    frozen = json.loads(
        (
            ROOT
            / "data/enumerations/h9_14_record14_ridge_capacity_hazard.json"
        ).read_text(encoding="utf-8")
    )
    assert frozen == document
