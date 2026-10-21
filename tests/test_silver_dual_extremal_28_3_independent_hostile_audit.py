import importlib.util
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "scripts/independent_hostile_audit_silver_dual_extremal_28_3_first_prefix.py"
)
ARTIFACT = (
    ROOT
    / "data/enumerations/"
    "silver_dual_extremal_28_3_first_prefix_independent_hostile_audit.json"
)
REPORT = (
    ROOT
    / "results/"
    "silver_dual_extremal_28_3_first_prefix_independent_hostile_audit.md"
)


def load_audit():
    spec = importlib.util.spec_from_file_location("dual_extremal_28_3_hostile", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def payload() -> dict[str, object]:
    return json.loads(ARTIFACT.read_text())


def test_hostile_replay_is_owner_import_free_and_regenerates() -> None:
    source = SCRIPT.read_text()
    assert "from coxeter_polytopes.search.silver_dual_extremal_28_3 import" not in source
    assert "import coxeter_polytopes.search.silver_dual_extremal_28_3" not in source
    subprocess.run(
        [sys.executable, str(SCRIPT), "--verify"],
        cwd=ROOT,
        check=True,
    )


def test_literal_gram_and_independent_s3_replay() -> None:
    data = payload()
    assert data["verdict"] == "PASS"
    gram = data["literal_gram"]
    assert gram == {
        "determinant": 3**26,
        "dimension": 28,
        "even_integral_symmetric": True,
        "minimum": 8,
        "minimum_shell_signed_count": 1_260,
        "normalized_stream_sha256": (
            "d61f9a887f3b31db364b6c38c9b0e59b615c4cd178201ae3e8d0036a87629c54"
        ),
        "rank_mod_three": 2,
        "smith_invariants": [1, 1] + [3] * 26,
    }
    s3 = data["pairing_radical_S3"]
    assert s3["preimage_index"] == 9
    assert s3["minimum"] == 12
    assert s3["minimum_shell_signed_count"] == 65_520
    assert s3["gram_pairings_divisible_by_three"] is True


def test_standalone_kernel_and_complete_local_alphabet() -> None:
    audit = load_audit()
    gram, digest, _source_digest = audit.parse_literal_gram()
    assert digest == audit.EXPECTED_GRAM_STREAM_SHA256
    assert audit.rank_mod_prime(gram, 3) == 2
    s3_basis = audit.pairing_radical_basis_mod_three(gram)
    assert abs(int(s3_basis.det())) == 9
    assert all(int(value) % 3 == 0 for value in gram * s3_basis)
    assert [str(gate.norm) for gate in audit.NORM_GATES] == [
        "1",
        "3",
        "3 + 1*sqrt(2)",
        "9 + 3*sqrt(2)",
        "2",
        "6",
        "6 + 2*sqrt(2)",
        "18 + 6*sqrt(2)",
    ]
    assert all(not audit.height_zero_triples(gate) for gate in audit.NORM_GATES)


def test_complete_chronology_endpoint_types_and_augmented_maxima() -> None:
    data = payload()
    chronology = data["chronology"]
    assert chronology["control_count"] == 51
    assert chronology["distinct_height_count"] == 35
    assert chronology["full_control_stream_sha256"] == (
        "ff6d4b175d8c9268ec565e4ccf06e979342b6cf55dc325c91bd6ba8574bdcedd"
    )
    assert chronology["earlier_nonempty_necessary_controls"] == 0
    assert chronology["endpoint_component_count"] == 2
    assert chronology["endpoint_type_count"] == 6
    rows = [row for component in chronology["endpoint"] for row in component["types"]]
    assert [row["gram"] for row in rows] == [
        [138, 72, 99],
        [300, 132, 198],
        [288, 138, 198],
        [276, 144, 198],
        [264, 150, 198],
        [252, 156, 198],
    ]
    assert [row["maximum"] for row in rows] == [
        "24 + 18*sqrt(2)",
        "45 + 30*sqrt(2)",
        "44 + 31*sqrt(2)",
        "44 + 32*sqrt(2)",
        "40 + 31*sqrt(2)",
        "33 + 27*sqrt(2)",
    ]
    assert [row["retained_extensions"] for row in rows] == [
        169, 361, 481, 529, 417, 199,
    ]
    assert [row["augmented_determinant"] for row in rows] == [
        432, 468, 528, 480, 504, 108,
    ]


def test_explicit_d_over_two_witness_is_strict_and_scope_is_bounded() -> None:
    data = payload()
    witness = data["strict_hyperideal_witness"]
    assert witness["explicit_direction"].startswith("d=e0")
    assert witness["norm"] == "3 - 2*sqrt(2)"
    assert witness["norm_exact"] == {"a": "3", "b": "-2"}
    assert witness["all_six_margins_strict"] is True
    assert witness["bounded_prefix_finite_volume"] is False
    assert data["production_comparison"] == {
        "all_requested_claims_match": True,
        "defects_found": [],
    }
    assert data["conclusion"] == {
        "bounded_prefix_has_open_hyperideal_patch": True,
        "compact_or_finite_volume_example": False,
        "eventual_chamber": "UNRESOLVED",
    }
    assert "no later-wall or termination claim" in data["scope"]
    assert "No production defect was found" in REPORT.read_text()
