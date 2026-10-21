from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "scripts/certify_gram_belyi_m5_degree29_five_black_ww_partition33_one_critical_towers.py"
)
ARCHIVE = (
    ROOT
    / "data/enumerations/gram_belyi_m5_degree29_five_black_ww_partition33_one_critical_towers.json"
)
REPORT = (
    ROOT
    / "results/gram_belyi_m5_degree29_five_black_ww_partition33_cs_full_sc_tranche.md"
)


def load_module():
    sys.path.insert(0, str(ROOT / "scripts"))
    spec = importlib.util.spec_from_file_location("partition33_one_critical", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def load_archive():
    return json.loads(ARCHIVE.read_text())


def test_exact_scoped_row_counts_and_order():
    module = load_module()
    archive = load_archive()
    records = module.records(False)
    cs = [row for row in archive["rows"] if row["orientation"] == "CS"]
    sc = [row for row in archive["rows"] if row["orientation"] == "SC"]
    assert len(records) == 1280
    assert len(cs) == 1280
    assert len(sc) == 41
    assert [tuple(row[key] for key in ("a_lambda", "b", "c", "d", "e")) for row in cs] == records
    assert all(17 <= row["a_lambda"] <= 20 for row in sc)


def test_every_row_has_exact_successive_tower_certificate():
    archive = load_archive()
    expected_stdout = "DERIVATIVE_BASE\n1\nCRITICAL_VALUE_BASE\n0\n696\nUNIT"
    for row in archive["rows"]:
        assert row["derivative_base_dimension"] == 1
        assert row["critical_value_base_dimension"] == 0
        assert row["critical_value_quotient_dimension"] == 696
        assert row["successive_simple_endpoint_standard_basis"] == [1]
        assert row["singular_stdout"] == expected_stdout
        assert row["uniform_partition_infinity_lemma_applied"] is True


def test_orientations_are_formulated_separately():
    module = load_module()
    _sqrt5, lam = module.lambda_mod(101)
    _context, derivative = module.residual_derivative(101, lam, (1, 1, 1, 1, 25))
    cs = module.equation_strings((1, 1, 1, 1, 25), derivative, "CS")
    sc = module.equation_strings((1, 1, 1, 1, 25), derivative, "SC")
    assert cs != sc
    assert "29*r^2" in cs[2] and "+58*r" in cs[3]
    assert "29*r^2" in sc[0] and "58*(1+r)" in sc[3]
    assert archive_reflection_flag() is False


def archive_reflection_flag():
    return load_archive()["reflection_symmetry_used"]


def test_exact_ledger_and_open_boundary():
    archive = load_archive()
    assert archive["new_marked_passports_eliminated"] == 1321
    assert archive["partition_marked_passports_remaining"] == 2519
    assert archive["partition_remaining_endpoint_status_counts"] == {
        "SC": 1239,
        "SS": 1280,
    }
    assert archive["five_black_ww_marked_passports_remaining"] == 27067
    assert archive["five_black_marked_passports_remaining"] == 128467
    assert archive["global_marked_passports_remaining"] == 42500755
    assert archive["global_degree_29_status"] == "OPEN"


def test_claim_flags_do_not_overreach():
    archive = load_archive()
    assert archive["complete_CS_no_go_claimed"] is True
    assert archive["complete_SC_tranche_no_go_claimed"] is True
    for key in (
        "complete_SC_no_go_claimed",
        "complete_partition_no_go_claimed",
        "complete_five_black_ww_no_go_claimed",
        "complete_degree_29_no_go_claimed",
        "degree_29_witness_found",
        "broad_plane_tree_search_used",
    ):
        assert archive[key] is False


def test_archive_is_canonical_and_report_scoped():
    module = load_module()
    archive = load_archive()
    assert ARCHIVE.read_bytes() == module.canonical_bytes(archive)
    assert sum(archive["prime_histogram"].values()) == 1321
    text = REPORT.read_text()
    assert "1,239" in text and "2,519" in text and "42,500,755" in text
    assert "remain open" in text.lower()
    assert "not stored, counted, or promoted" in text
