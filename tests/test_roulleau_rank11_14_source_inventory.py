import json
from pathlib import Path

from scripts.certify_roulleau_rank11_14_source_inventory import (
    CASES,
    SOURCE_HASHES,
    build_artifact,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT / "data/enumerations/roulleau_rank11_14_source_stock_inventory.json"
)


def test_frozen_artifact_rebuilds_exactly():
    frozen = json.loads(ARTIFACT.read_text())
    assert build_artifact(frozen["source_check"]) == frozen


def test_denominator_is_finite_and_exact():
    frozen = json.loads(ARTIFACT.read_text())
    denominator = frozen["denominator"]
    assert len(CASES) == 19
    assert denominator["rank11_14_stock_count"] == 14
    assert denominator["rank11_14_total_wall_count"] == 375
    assert denominator["rank11_14_direct_complete_graph_count"] == 9
    assert denominator["rank11_14_basis_graph_plus_formula_count"] == 5
    assert denominator["rank11_14_induced_chamber_completion_count"] == 3
    assert denominator["already_exactly_machine_represented_count"] == 19
    assert denominator["new_transcription_queue_count"] == 0
    assert denominator["new_transcription_queue_keys"] == []


def test_every_stock_has_exact_source_and_extraction_metadata():
    frozen = json.loads(ARTIFACT.read_text())
    keys = set()
    for stock in frozen["stocks"]:
        assert stock["key"] not in keys
        keys.add(stock["key"])
        assert stock["picard_rank"] >= 11
        assert stock["hyperbolic_dimension"] == stock["picard_rank"] - 1
        assert stock["wall_count"] > stock["hyperbolic_dimension"]
        assert stock["exact_wall_gram_extractable"] is True
        start, end = stock["source"]["section_lines"]
        assert 1 <= start <= end
        tstart, tend = stock["source"]["table_lines"]
        assert 10906 <= tstart <= tend <= 10966


def test_pinned_source_hashes_and_positive_status():
    frozen = json.loads(ARTIFACT.read_text())
    assert frozen["source"]["sha256"] == SOURCE_HASHES
    assert frozen["source_check"] == {"checked": True, "all_checks_pass": True}
    status = frozen["mathematical_status"]
    assert status["full_nef_chambers_are_positive_examples"] is True
    assert status["volume_type"] == "finite-volume_noncompact"
