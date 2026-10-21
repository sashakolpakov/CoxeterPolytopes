from __future__ import annotations

from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path

from scripts.independent_h9_14_record30_a0_dual_cnf import CNF
from scripts.hostile_h9_14_record30_a0_dual_cnf_orbits import (
    DEFAULT_CENSUS,
    SEED,
    family_orbit,
    full_group_generators,
    independent,
    normalize_targets,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "scripts/hostile_h9_14_record30_a0_dual_cnf_orbits.py"
ARCHIVE = (
    ROOT / "data/enumerations/h9_14_record30_a0_dual_orbit_hostile.json"
)
EXPECTED_SOURCE_SHA256 = (
    "f49d28b445adffa22b326582dc8b4bc7d6018e3699c093d2140dc21ca81ca3b7"
)
EXPECTED_ARCHIVE_SHA256 = (
    "b33237d52eb8206cad7c8c6c5bc95e43c7ec54dac5367030a0c55fc9803122b5"
)


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _satisfies(formula: CNF, values: tuple[bool, ...]) -> bool:
    assignment = (False,) + values
    return all(
        any(assignment[abs(literal)] == (literal > 0) for literal in clause)
        for clause in formula.clauses
    )


def _has_extension(formula: CNF, primary: tuple[bool, ...]) -> bool:
    auxiliary_count = formula.variable_count - len(primary)
    return any(
        _satisfies(formula, primary + auxiliary)
        for auxiliary in product((False, True), repeat=auxiliary_count)
    )


def test_frozen_no_import_orbit_audit_and_source() -> None:
    assert _digest(SOURCE) == EXPECTED_SOURCE_SHA256
    assert _digest(ARCHIVE) == EXPECTED_ARCHIVE_SHA256
    document = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    payload = dict(document)
    declared_summary = payload.pop("summary_sha256")
    assert sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest() == declared_summary
    assert document["input_target_orbit_count"] == 8
    assert document["blocked_seed_form_count"] == 63_960
    assert document["blocked_seed_families_sha256"] == (
        "591e61e6cb33cbd7fb039cb6313336fa022c0ce33ad15143a81cc94778f5a25e"
    )
    assert document["blocking_clauses_sha256"] == (
        "76c2b85f93cdc562af26bb9cee0841b872a0b94b24ed86d3926556c0d52d56d2"
    )
    assert document["cnf_sha256"] == (
        "531eedff4e8ba07f3750351d3416f9834abf6fbfd6e0d338fedbb7cc463f511c"
    )
    assert document["cnf_tail_is_exactly_the_independent_blocking_clause_set"] is True
    assert document["all_target_orbits_pairwise_disjoint"] is True
    assert document["every_seed_count_matches_orbit_incidence_double_count"] is True
    assert document["every_orbit_has_a_seed_stabilizer_lex_witness"] is True


def test_all_frozen_orbit_stabilizer_double_counts_are_exact() -> None:
    document = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    expected_seed_counts = (1200, 14760, 14640, 3600, 3720, 14880, 3720, 7440)
    audits = document["target_orbit_audits"]
    assert tuple(item["direct_seed_form_count"] for item in audits) == expected_seed_counts
    for item in audits:
        orbit = item["full_family_orbit_size"]
        stabilizer = item["family_stabilizer_order"]
        type_three = item["type_three_blocker_count"]
        seed_forms = item["direct_seed_form_count"]
        assert orbit * stabilizer == document["full_group_order"]
        assert orbit * type_three == (
            document["type_three_complement_orbit_size"] * seed_forms
        )
        assert item["orbit_incidence_seed_form_count"] == seed_forms


def test_first_target_orbit_is_recomputed_from_elementary_generators() -> None:
    census = json.loads(DEFAULT_CENSUS.read_text(encoding="utf-8"))
    target = normalize_targets(census)[0]
    orbit = family_orbit(target, full_group_generators())
    seed_forms = {
        family for family in orbit if all(SEED & row for row in family)
    }
    assert len(orbit) == 2_240
    assert len(seed_forms) == 1_200


def test_blocker_biconditional_pattern_is_exact_by_truth_table() -> None:
    formula = CNF()
    facet = formula.variable()
    disjoint_rows = tuple(formula.variable() for _ in range(4))
    for selected in disjoint_rows:
        formula.clause(-facet, -selected)
    formula.clause(facet, *disjoint_rows)
    for values in product((False, True), repeat=5):
        expected_facet = not any(values[1:])
        assert _satisfies(formula, values) == (values[0] == expected_facet)


def test_ridge_zero_or_two_pattern_is_exact_by_truth_table() -> None:
    for size in range(1, 7):
        formula = CNF()
        members = tuple(formula.variable() for _ in range(size))
        for triple in combinations(members, 3):
            formula.clause(*(-item for item in triple))
        for item in members:
            formula.clause(
                -item, *(other for other in members if other != item)
            )
        for values in product((False, True), repeat=size):
            assert _satisfies(formula, values) == (sum(values) in (0, 2))


def test_missing_face_indicator_pattern_is_exact_by_truth_table() -> None:
    formula = CNF()
    selected_rows = tuple(formula.variable() for _ in range(4))
    missing = formula.variable()
    for selected in selected_rows:
        formula.clause(-selected, missing)
    formula.clause(-missing, *selected_rows)
    for values in product((False, True), repeat=5):
        assert _satisfies(formula, values) == (values[-1] == any(values[:-1]))


def test_lex_encoder_remains_exact_when_right_is_a_permutation_of_left() -> None:
    for size in range(1, 5):
        for permutation in permutations(range(size)):
            formula = CNF()
            left = tuple(formula.variable() for _ in range(size))
            right = tuple(left[index] for index in permutation)
            formula.lex_greater_equal(left, right)
            for values in product((False, True), repeat=size):
                moved = tuple(values[index] for index in permutation)
                assert _has_extension(formula, values) == (values >= moved)


def test_optional_face_vector_counts_accept_every_known_target() -> None:
    census = json.loads(DEFAULT_CENSUS.read_text(encoding="utf-8"))
    for family in normalize_targets(census):
        missing_counts = []
        for size in (6, 7):
            missing = 0
            for chosen in combinations(range(14), size):
                support = sum(1 << vertex for vertex in chosen)
                if independent(support) and any(
                    support & row == row for row in family
                ):
                    missing += 1
            missing_counts.append(missing)
        assert tuple(missing_counts) == (56, 160)
