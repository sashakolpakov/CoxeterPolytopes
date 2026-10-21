import json
from pathlib import Path

from scripts.certify_roulleau_ue8d4_all_subgroup_folds import (
    _ambient_reconstruction,
    _source_gram,
    build_certificate,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/roulleau_ue8d4_all_subgroup_folds.json"


def test_source_transcription_and_ns_basis() -> None:
    gram = _source_gram()
    basis_gram, wall_rows = _ambient_reconstruction(gram)
    assert gram.rank() == 14
    assert exact_inertia(gram) == (1, 13, 1)
    assert basis_gram.det() == -4
    assert wall_rows * basis_gram * wall_rows.T == gram


def test_archive_rebuilds_exactly() -> None:
    assert build_certificate() == json.loads(ARTIFACT.read_text())


def test_all_subgroups_give_two_known_index_one_chambers() -> None:
    payload = json.loads(ARTIFACT.read_text())
    assert payload["all_checks_pass"] is True
    assert payload["weighted_diagram_group"]["structure"] == "S3"
    assert payload["subgroup_census"]["conjugacy_class_count"] == 4
    assert payload["subgroup_census"]["proper_high_class_count"] == 3
    assert [
        (record["geometry"]["hyperbolic_dimension"], record["geometry"]["facet_count"])
        for record in payload["types"]
    ] == [(12, 14), (11, 13)]
    assert all(
        record["geometry"]["finite_volume_certified"]
        and record["geometry"]["compact"] is False
        and record["maximal_reflection_lineage"]["reflection_subgroup_index"] == 1
        and record["ns_action"]["surface_realization_status"].startswith(
            "not_realized_on_the_general_K3"
        )
        for record in payload["types"]
    )

