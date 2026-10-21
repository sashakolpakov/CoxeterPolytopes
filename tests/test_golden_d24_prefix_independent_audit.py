from pathlib import Path

from scripts.audit_golden_d24_prefix_independent import build_audit


def test_independent_golden_d24_and_d24_plus_prefix_counts_and_membership():
    payload = build_audit()
    conclusion = payload["conclusion"]
    assert conclusion["root_lattice_candidate_counts"] == [1, 3, 4]
    assert conclusion["root_lattice_accepted_counts"] == [1, 2, 1]
    assert conclusion["d24_plus_candidate_counts"] == [1, 3, 4]
    assert conclusion["d24_plus_accepted_counts"] == [1, 2, 1]
    assert conclusion["spinor_candidates_through_phi_squared"] == 0
    assert conclusion["accepted_root_count"] == 4
    assert conclusion["accepted_roots_all_in_d24_root_coset"]
    assert conclusion["all_accepted_norm_two"]
    assert conclusion["all_accepted_pairings_acute"]

    accepted = payload["accepted_roots_canonical_order"]
    assert [len(root["nonzero_support"]) for root in accepted] == [2, 4, 1, 6]
    assert all(root["membership_pass"] and root["d24_root_member"] for root in accepted)
    assert all(root["d24_plus_spinor_residue"] is None for root in accepted)


def test_independent_golden_d24_pairing_matrix_and_rejection_witnesses():
    payload = build_audit()
    matrix = [
        [entry["display"] for entry in row]
        for row in payload["accepted_root_pairing_matrix"]
    ]
    assert matrix == [
        ["2", "0", "0", "0"],
        ["0", "2", "-2*phi", "0"],
        ["0", "-2*phi", "2", "-2 - 2*phi"],
        ["0", "0", "-2 - 2*phi", "2"],
    ]
    levels = payload["runs"]["D24-plus"]["levels"]
    rejected = [
        candidate
        for level in levels
        for candidate in level["candidates"]
        if not candidate["accepted"]
    ]
    assert len(rejected) == 4
    assert all(candidate["positive_pairing_witnesses"] for candidate in rejected)


def test_independent_audit_has_no_production_search_dependency():
    source = Path("scripts/audit_golden_d24_prefix_independent.py").read_text(
        encoding="utf-8"
    )
    assert "coxeter_polytopes.search.golden_dn" not in source
    assert "coxeter_polytopes.exact.golden" not in source

