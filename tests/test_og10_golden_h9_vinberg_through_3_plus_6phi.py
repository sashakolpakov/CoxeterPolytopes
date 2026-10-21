import json
from pathlib import Path

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.search.bugaenko_vinberg import _pair
from coxeter_polytopes.search.og10_golden_h9_vinberg import h9_field_gram
from coxeter_polytopes.verify.chamber import verify_hyperideal_extreme_ray
from scripts.certify_og10_golden_h9_vinberg_through_3_plus_6phi import (
    _digest,
    build,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT / "data/candidates/og10_golden_h9_vinberg_through_3_plus_6phi.json"
)


def _roots(payload):
    return tuple(
        tuple(Golden.from_json(entry) for entry in root)
        for root in payload["simple_roots"]
    )


def test_direct_h9_vinberg_through_3_plus_6phi_is_reproducible() -> None:
    frozen = json.loads(ARTIFACT.read_text())
    result = build()
    assert result == frozen
    assert result["verified"] is True
    assert all(result["checks"].values())

    assert result["completed_through_control_coordinate"] == Golden(3, 6).to_json()
    assert result["next_admissible_control_coordinate"] == Golden(4, 6).to_json()
    assert result["facet_count"] == 71
    assert result["positive_form_conditional_diagonal_coefficients"] == [
        "1",
        "4/3",
        "3/2",
        "8/5",
        "5/3",
        "12/7",
        "7/4",
        "16",
        "4",
    ]

    calibration = result["calibration"]
    assert calibration["ordered_root_match"] is True
    assert calibration["ambient_fixed_checkpoint_root_count"] == 61
    assert calibration["replayed_root_count"] == 61
    assert (
        calibration["ambient_fixed_checkpoint_sha256"]
        == calibration["replayed_roots_sha256"]
        == "b83e9e3018debdccc15a9934e2ab735e9882f643f1f04ab155bb52f5a881eabc"
    )
    assert [level["root_count"] for level in calibration["levels"]] == [
        2,
        2,
        2,
        0,
        4,
        3,
        2,
        7,
        0,
        6,
        7,
        1,
        13,
        3,
    ]

    level = result["new_level"]
    assert level["control_coordinate"] == Golden(3, 6).to_json()
    assert level["root_count"] == 10
    assert level["visited_node_count"] == 156_899
    assert level["visited_by_depth"] == [
        1,
        49,
        784,
        2793,
        4802,
        8134,
        14504,
        27832,
        65170,
        32830,
    ]
    assert level["roots_sha256"] == (
        "cb22d2c82dd89c8b423a5a5e2cc78b31488b836c3d12aa642d9dd7011582498f"
    )
    assert level["old_hyperideal_ray_pairing_signs"] == {"negative": 10}

    cross_check = result["pell_fold_g_cross_check"]
    assert cross_check["images_of_61_outside_61_count"] == 12
    assert cross_check["next_level_forced_root_indices"] == [1, 4, 5, 8]
    assert cross_check["source_root_indices_in_61"] == [37, 51, 59, 60]

    roots = _roots(result)
    gram = h9_field_gram()
    assert len(roots) == len(set(roots)) == 71
    assert all(_pair(root, gram, root) == 2 for root in roots)
    assert _digest(result["simple_roots"]) == result["roots_sha256"] == (
        "3637dbdf1e42813fca18190fc02e07c831b4715fe6ad5fee32f09db5d34a8c39"
    )

    archived_witness = result["hyperideal_witness"]
    witness = verify_hyperideal_extreme_ray(
        roots,
        gram,
        archived_witness["active_facets"],
        dimension=9,
    )
    assert witness.ray == tuple(
        Golden.from_json(entry) for entry in archived_witness["ray"]
    )
    assert witness.norm == Golden(210, -128) / 169
    assert witness.norm.sign() > 0
    assert witness.ray == (
        Golden(1),
        2 * PHI,
        4 * PHI,
        6 * PHI,
        Golden(0, 64) / 13,
        Golden(0, 48) / 13,
        Golden(0, 32) / 13,
        Golden(0, 16) / 13,
        Golden(0, 40) / 13,
        PHI / 13,
    )
    assert (
        _digest([value.to_json() for value in witness.wall_pairings])
        == archived_witness["wall_pairings_sha256"]
        == "c051e5e266ba7a9bf39dbe725b892c5885a70f0c5055a7a335bd933d270e2816"
    )
