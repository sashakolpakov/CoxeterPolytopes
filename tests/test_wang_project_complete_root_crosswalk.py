import json
from pathlib import Path

from scripts.certify_wang_project_complete_root_crosswalk import analyze


ROOT = Path(__file__).resolve().parents[1]


def _rows():
    return {row["hyperbolic_lattice_S"]: row for row in analyze()["rows"]}


def test_high_rank_wang_split_and_exact_weights() -> None:
    report = analyze()
    assert report["checks"]["frozen_high_rows_are_19"]
    assert report["checks"]["frozen_high_split_is_7_plus_12"]
    rows = _rows()
    assert {
        lattice: row["weights"]
        for lattice, row in rows.items()
        if lattice in {
            "U+E8+A1",
            "U+D8+A1",
            "U+2D4+A1",
            "U+D4+5A1",
            "U+E8+2A1",
            "U+2E8",
            "U+2E8+A1",
        }
    } == {
        "U+E8+A1": [195],
        "U+D8+A1": [99],
        "U+2D4+A1": [51],
        "U+D4+5A1": [27],
        "U+E8+2A1": [138],
        "U+2E8": [132],
        "U+2E8+A1": [75],
    }


def test_special_nonsplit_corrections_are_not_misidentified() -> None:
    rows = _rows()
    assert rows["<2>+8A1(-1)"]["weights"] == [5]
    assert "U(2)+S" in rows["<2>+8A1(-1)"]["exact_form"]
    assert rows["U(2)+2D4=U+N8"]["weights"] == [28]
    assert "2U+N8" in rows["U(2)+2D4=U+N8"]["exact_form"]
    assert rows["U(2)+D4"]["weights"] == [40, 24]


def test_same_form_overlattice_chains_and_nonmonotonicity() -> None:
    report = analyze()
    edges = {
        (edge["source"], edge["target"], edge["index"])
        for edge in report["same_form_even_overlattices"]
    }
    assert {
        ("U+D4+5A1", "U+2D4+A1", 2),
        ("U+2D4+A1", "U+D8+A1", 2),
        ("U+D8+A1", "U+E8+A1", 2),
        ("U+2D4+2A1", "U+D8+2A1", 2),
        ("U+D8+2A1", "U+E8+2A1", 2),
    } <= edges
    rows = _rows()
    assert not rows["U+2D4+2A1"]["automorphic_correction"]
    assert not rows["U+D8+2A1"]["automorphic_correction"]
    assert rows["U+E8+2A1"]["automorphic_correction"]


def test_wang_weight_formula_replay() -> None:
    report = analyze()
    assert report["checks"]["all_formula_rows_pass"]
    assert report["formula_checks"]["U+E8+A1"]["computed_weight"] == 195
    assert report["formula_checks"]["U+D4+5A1"]["computed_weight"] == 27


def test_exact_project_aliases_land_on_the_weighted_rows() -> None:
    report = analyze()
    aliases = {
        item["alias"]: item["canonical_row"] for item in report["exact_aliases"]
    }
    assert aliases["U+D10"] == "U+E8+2A1"
    assert aliases["U+D6+3A1"] == "U+2D4+A1"
    assert aliases["U+E7+2A1"] == "U+D8+A1"
    assert aliases["U+N8+A1"] == "U+D4+5A1"
    rows = _rows()
    assert rows[aliases["U+D10"]]["weights"] == [138]
    assert rows[aliases["U+D6+3A1"]]["weights"] == [51]
    assert rows[aliases["U+E7+2A1"]]["weights"] == [99]
    assert rows[aliases["U+N8+A1"]]["weights"] == [27]


def test_frozen_archive_is_exact() -> None:
    path = ROOT / "data/enumerations/wang_project_complete_root_crosswalk.json"
    assert json.loads(path.read_text()) == analyze()
