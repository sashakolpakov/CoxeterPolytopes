import importlib.util
import json
import sys
from hashlib import sha256
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/independent_hostile_audit_silver_p7_d24_plus_h500.py"
ARTIFACT = (
    ROOT
    / "data/enumerations/silver_p7_d24_plus_h500_independent_hostile_audit.json"
)


def _load_audit_module():
    spec = importlib.util.spec_from_file_location("silver_h500_hostile_audit", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_frozen_hostile_certificate_records_the_full_replay() -> None:
    payload = json.loads(ARTIFACT.read_text())
    assert payload["verdict"] == "PASS"
    assert payload["control_completeness"]["complete_control_count"] == 64
    assert payload["control_completeness"]["extension_control_count"] == 16
    assert payload["dense_packet"] == {
        "height": {"a": "181", "b": "128"},
        "raw_reflective": 4173,
        "imprimitive_sqrt2": 83,
        "primitive_candidates": 4090,
        "simultaneous_survivors": 4003,
        "same_height_unordered_pairs": 8_010_003,
        "same_height_positive_pairs": 0,
        "ordering_leakage": False,
    }
    matrix = payload["integer_matrix_audit"]
    assert matrix["matrix_entries"] == 17_587_030
    assert matrix["maximum_proved_intermediate_absolute_bound"] < matrix["int64_max"]
    assert matrix["maximum_abs_pair_coefficient"] == 4632
    assert matrix["scalar_zero_replays"] == 78_258
    assert payload["old_h300_ray_cut"]["positive"] == 74
    ray = payload["replacement_hyperideal_extreme_ray"]
    assert (ray["negative_pairings"], ray["zero_pairings"], ray["positive_pairings"]) == (
        4350,
        24,
        0,
    )
    assert ray["active_rank"] == 24
    assert ray["norm"] == {"a": "693/5", "b": "-98"}
    assert ray["limiting_wall_index"] == 4361
    assert payload["conclusion"]["h500_prefix_finite_volume"] is False
    assert payload["conclusion"]["eventual_chamber_status"] == "unresolved"


def test_clean_room_control_and_geometry_replay_without_production_helpers() -> None:
    audit = _load_audit_module()
    levels, certificate = audit.enumerate_control_levels(audit.H500)
    assert len(levels) == 64
    assert [item["energy_ceiling"] for item in certificate["boxes"]] == [
        137,
        274,
        604,
        1208,
    ]

    # Rebuild the archived wall stream with only the standalone K2 parser and
    # directly replay the final ray pairings and active rank.
    h300 = json.loads(audit.H300_PATH.read_text())
    h500 = json.loads(audit.H500_PATH.read_text())
    assert sha256(audit.H300_PATH.read_bytes()).hexdigest() == audit.H300_SHA256
    walls = list(audit.negative_d24_simple_walls())
    for packet in h300["enumeration"]["packets"]:
        walls.extend(
            audit.vector(item["vector_orthonormal"])
            for item in packet["accepted_roots"]
        )
    for packet in h500["enumeration"]["extension_packets"]:
        walls.extend(
            audit.vector(item["vector_orthonormal"])
            for item in packet["accepted_roots"]
        )
    assert len(walls) == 4374
    c = audit.K2.make(4, "-12/5")
    ray = (audit.ONE, *(c for _ in range(5)), *(audit.ZERO for _ in range(19)))
    pairings = tuple(audit.lorentz(ray, wall) for wall in walls)
    assert [item.sign() for item in pairings].count(1) == 0
    active = tuple(index for index, item in enumerate(pairings) if not item)
    assert active == (*range(4), *range(5, 24), 4361)
    rank, _pivots = audit.row_rank(
        tuple(audit.functional(walls[index]) for index in active), 25
    )
    assert rank == 24
    assert audit.lorentz(ray, ray) == audit.K2.make("693/5", -98)


def test_first_dense_shell_is_reconstructed_from_all_four_glue_residues() -> None:
    audit = _load_audit_module()
    shell = audit.enumerate_d24_plus_shell(audit.K2.make(182, 128))
    assert len(shell) == 83
    counts = {label: sum(item_label == label for _vector, item_label in shell) for label in {
        "D24",
        "spinor-1",
        "spinor-1*sqrt(2)",
        "spinor-1 + 1*sqrt(2)",
    }}
    assert counts == {
        "D24": 62,
        "spinor-1": 12,
        "spinor-1*sqrt(2)": 4,
        "spinor-1 + 1*sqrt(2)": 5,
    }
