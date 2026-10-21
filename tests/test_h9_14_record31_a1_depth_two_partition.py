import json
from pathlib import Path

from scripts.certify_h9_14_record31_a1_depth_two_partition import verify
from scripts.generate_h9_14_record31_a1_depth_two_class_engine import (
    generate_depth_two_class_engine,
)


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = (
    ROOT
    / "data/enumerations/"
    "h9_14_record31_a1_depth_two_partition_equivalence.json"
)


def test_depth_two_partition_equivalence_archive() -> None:
    result = verify(json.loads(ARCHIVE.read_text()))
    assert result["root_fragment_count"] == 4
    assert result["individual_depth_two_fragment_count"] == 22
    assert result["canonical_target_union_count"] == 1


def test_depth_two_partition_records_skipped_codes_before_return() -> None:
    source = generate_depth_two_class_engine()
    assert "audit->depth_two_class_codes.insert(canonical_code)" in source
    insertion = source.index("audit->depth_two_class_codes.insert(canonical_code)")
    selection = source.index("depth_two_class_id < gDepthTwoClassBegin")
    assert insertion < selection
    assert "--depth-two-classes ROOT BEGIN END" in source
