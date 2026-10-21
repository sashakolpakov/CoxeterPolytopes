from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from scripts.generate_h9_14_record31_a1_addition_engine import generate
from scripts.generate_h9_14_record31_a1_compact_key_engine import generate_compact
from scripts.generate_h9_14_record31_a1_inline_key_engine import generate_inline
from scripts.verify_h9_14_record31_a1_inline_key_equivalence import (
    REFERENCE_DEPTH_NODES,
    REFERENCE_SEEN_COUNTS,
    OUTPUT,
    validate_transcripts,
)


ROOT = Path(__file__).resolve().parents[1]
ORBIT0 = (
    ROOT
    / "data"
    / "enumerations"
    / "h9_14_record31_a1_orbit0_corrected.json"
)


def test_inline_generator_changes_only_quotient_memo_storage_path() -> None:
    dense = generate()
    compact = generate_compact()
    inline = generate_inline()

    assert len({dense, compact, inline}) == 3
    assert "struct InlineQuotientKey" in inline
    assert "std::array<std::uint64_t, 3> words{}" in inline
    assert "code.size() > 22" in inline
    assert "row_count * kN + 7" in inline
    assert "static_assert(sizeof(InlineQuotientKey) == 3" in inline
    assert "InlineQuotientKeyHash" in inline
    assert "InlineQuotientKeyEqual" in inline
    assert "left.words == right.words" in inline
    assert "matrix_quotient_seen" in inline
    assert "const InlineQuotientKey canonical(CanonicalPartialCode" in inline

    # Only the quotient memo adopts the fixed-width key.  The independent
    # matrix quotient and target-orbit canonical maps deliberately keep their
    # pre-existing string representations.
    assert "std::array<std::unordered_set<std::string>" in inline
    assert "std::map<std::string, std::pair<std::vector<Mask>" in inline
    assert "std::string CanonicalFamilyCode" in inline
    assert "std::string CanonicalMatrixCode" in inline
    assert "compact canonical key failed dense-graph round trip" in inline


def test_zero_padded_reference_packing_keeps_row_count_and_all_bytes() -> None:
    def pack(code: bytes) -> tuple[int, int, int]:
        assert 0 < len(code) <= 22
        words = [0, 0, 0]
        for index, value in enumerate(code):
            words[index // 8] |= value << (8 * (index % 8))
        return tuple(words)  # type: ignore[return-value]

    five_rows = bytes([5]) + bytes(9)
    six_rows = bytes([6]) + bytes(10)
    assert pack(five_rows) != pack(six_rows)

    prefix = bytes([12]) + bytes(20)
    differing_last_byte = prefix + b"\x01"
    zero_last_byte = prefix + b"\x00"
    assert len(differing_last_byte) == len(zero_last_byte) == 22
    assert pack(differing_last_byte) != pack(zero_last_byte)
    assert pack(zero_last_byte)[2] == 0


def test_triple_transcript_validator_is_exact_without_running_a_replay() -> None:
    document = json.loads(ORBIT0.read_text(encoding="utf-8"))
    row = validate_transcripts(document, document, document)
    assert row["partition_depth_node_counts"] == REFERENCE_DEPTH_NODES
    assert row["quotient_seen_counts"] == REFERENCE_SEEN_COUNTS

    changed = copy.deepcopy(document)
    changed["orbits"][0]["quotient_seen_counts"][-1] += 1
    with pytest.raises(ArithmeticError, match="transcripts differ"):
        validate_transcripts(document, document, changed)


def test_archived_triple_equivalence_certificate() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert document["dense_compact_inline_transcripts_identical"] is True
    assert document["compact_and_inline_dense_round_trip_enabled"] is True
    assert document["inline_key_bytes"] == 24
    assert document["maximum_compact_code_bytes"] == 22
    assert document["inline_source_sha256"] == (
        "ec577df47018f1d03d48e8304c56a2541ef02180144937600569ec0a9a7e68c7"
    )
