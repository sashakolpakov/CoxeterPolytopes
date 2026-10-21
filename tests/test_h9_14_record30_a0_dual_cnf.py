from __future__ import annotations

from itertools import product
import json

from scripts.certify_h9_14_record30_a0_dual_exhaustion import (
    ARCHIVE,
    verify,
)

from scripts.independent_h9_14_record30_a0_dual_cnf import (
    CNF,
    DEFAULT_CENSUS,
    SEED,
    normalize_targets,
    semantic_target_audit,
    seed_stabilizer_generators,
    transform,
)


def _has_extension(formula: CNF, primary: tuple[bool, ...]) -> bool:
    auxiliary_count = formula.variable_count - len(primary)
    for auxiliary in product((False, True), repeat=auxiliary_count):
        assignment = (False,) + primary + auxiliary
        if all(
            any(
                assignment[abs(literal)] == (literal > 0)
                for literal in clause
            )
            for clause in formula.clauses
        ):
            return True
    return False


def test_exact_unary_counter_is_logically_exact() -> None:
    for size in range(1, 5):
        for target in range(size + 1):
            formula = CNF()
            variables = tuple(formula.variable() for _ in range(size))
            formula.exactly(variables, target)
            for primary in product((False, True), repeat=size):
                assert _has_extension(formula, primary) == (
                    sum(primary) == target
                )


def test_at_most_unary_counter_is_logically_exact() -> None:
    for size in range(1, 5):
        for bound in range(size + 1):
            formula = CNF()
            variables = tuple(formula.variable() for _ in range(size))
            formula.at_most(variables, bound)
            for primary in product((False, True), repeat=size):
                assert _has_extension(formula, primary) == (
                    sum(primary) <= bound
                )


def test_lex_encoding_matches_boolean_lex_order() -> None:
    size = 3
    formula = CNF()
    left = tuple(formula.variable() for _ in range(size))
    right = tuple(formula.variable() for _ in range(size))
    formula.lex_greater_equal(left, right)
    for primary in product((False, True), repeat=2 * size):
        assert _has_extension(formula, primary) == (
            primary[:size] >= primary[size:]
        )


def test_seed_stabilizer_generators_generate_order_8640() -> None:
    generators = seed_stabilizer_generators()
    identity = tuple(range(14))
    reached = {identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            composed = tuple(generator[current[index]] for index in range(14))
            if composed not in reached:
                reached.add(composed)
                frontier.append(composed)
    assert len(reached) == 8640
    assert all(transform(SEED, action) == SEED for action in reached)


def test_all_eight_input_targets_satisfy_dual_formula_semantics() -> None:
    document = json.loads(DEFAULT_CENSUS.read_text(encoding="utf-8"))
    targets = normalize_targets(document)
    assert len(targets) == 8
    audits = [semantic_target_audit(target) for target in targets]
    assert len({row["A5_family_sha256"] for row in audits}) == 8
    assert all(row["blocker_count"] == 132 for row in audits)
    assert all(row["ridge_degree_histogram"] == {"2": 594} for row in audits)
    assert all(row["surviving_small_transversal_count"] == 0 for row in audits)


def test_frozen_dual_cnf_drat_theorem_and_scope() -> None:
    document = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert verify(document) == {
        "target_orbit_count": 8,
        "blocked_seed_form_count": 63960,
        "cnf_variable_count": 226358,
        "cnf_clause_count": 984935,
        "drat_verified": True,
        "positive_witness_count": 8,
        "compact_label_survivor_count": 0,
        "finite_volume_noncompact_statement": False,
    }
