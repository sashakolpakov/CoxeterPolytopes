#!/usr/bin/env python3
"""Audit the candidate-bitset capacity hazard in Record 14's ridge engine."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import re

from scripts.explore_h9_14_record14_fixed_orbit_z3 import build_problem


ROOT = Path(__file__).resolve().parents[1]
SHARED_RECORD18 = ROOT / "scripts/explore_h9_14_record18_ridge_aware_additions.cpp"
ISOLATED_RECORD14 = ROOT / "scripts/explore_h9_14_record14_ridge_aware_engine.cpp"
DECLARATION = re.compile(
    r"using CandidateWord = std::array<std::uint64_t, (\d+)>;"
)
ISOLATED_HEADER = (
    "// Record-14 isolated fork of the ridge-aware engine.\n"
    "// The only structural change from the shared Record-18 engine is the\n"
    "// checked 16-word candidate bitset (capacity 1024 instead of 832).\n"
    "// Do not replace this with the 13-word engine: Record 14 has 994 candidates.\n\n"
)
CAPACITY_GUARD = (
    "  constexpr std::size_t kCandidateWordCapacity =\n"
    "      64 * std::tuple_size_v<CandidateWord>;\n"
    "  if (g_candidates.size() > kCandidateWordCapacity) {\n"
    "    throw std::runtime_error(\n"
    "        \"candidate count \" + std::to_string(g_candidates.size()) +\n"
    "        \" exceeds ridge CandidateWord capacity \" +\n"
    "        std::to_string(kCandidateWordCapacity));\n"
    "  }\n"
)


def candidate_word_count(source: Path) -> int:
    matches = DECLARATION.findall(source.read_text(encoding="utf-8"))
    if len(matches) != 1:
        raise AssertionError(f"CandidateWord declaration changed in {source}")
    return int(matches[0])


def require_capacity(candidate_count: int, word_count: int) -> None:
    capacity = 64 * word_count
    if candidate_count > capacity:
        raise ValueError(
            f"candidate count {candidate_count} exceeds CandidateWord capacity {capacity}"
        )


def audit() -> dict[str, object]:
    shared_words = candidate_word_count(SHARED_RECORD18)
    isolated_words = candidate_word_count(ISOLATED_RECORD14)
    if shared_words != 13 or isolated_words != 16:
        raise AssertionError("ridge-engine word capacities changed")
    rows = []
    for orbit_index in range(5):
        candidate_count = len(build_problem(1, orbit_index)["candidates"])
        rejected_by_pinned = False
        try:
            require_capacity(candidate_count, shared_words)
        except ValueError:
            rejected_by_pinned = True
        require_capacity(candidate_count, isolated_words)
        rows.append(
            {
                "source_positive_A4_orbit_index": orbit_index,
                "candidate_count": candidate_count,
                "rejected_by_13_word_backend": rejected_by_pinned,
                "admitted_by_checked_16_word_backend": True,
            }
        )
    expected = [994, 994, 784, 994, 994]
    if [row["candidate_count"] for row in rows] != expected:
        raise AssertionError("Record-14 a=1 candidate counts changed")
    if [row["rejected_by_13_word_backend"] for row in rows] != [
        True,
        True,
        False,
        True,
        True,
    ]:
        raise AssertionError("13-word rejection pattern changed")
    isolated_text = ISOLATED_RECORD14.read_text(encoding="utf-8")
    if "if (g_candidates.size() > kCandidateWordCapacity)" not in isolated_text:
        raise AssertionError("isolated engine lost its runtime capacity assertion")
    if not isolated_text.startswith(ISOLATED_HEADER):
        raise AssertionError("isolated engine provenance header changed")
    normalized = isolated_text.removeprefix(ISOLATED_HEADER)
    normalized = normalized.replace(
        "using CandidateWord = std::array<std::uint64_t, 16>;",
        "using CandidateWord = std::array<std::uint64_t, 13>;",
        1,
    )
    normalized = normalized.replace(CAPACITY_GUARD, "", 1)
    if normalized != SHARED_RECORD18.read_text(encoding="utf-8"):
        raise AssertionError(
            "isolated engine differs from the shared engine beyond its "
            "provenance header, 16-word type, and runtime capacity guard"
        )
    return {
        "kind": "h9_14_record14_ridge_candidate_word_capacity_hazard",
        "shared_record18_word_count": shared_words,
        "shared_record18_capacity": 64 * shared_words,
        "isolated_record14_word_count": isolated_words,
        "isolated_record14_capacity": 64 * isolated_words,
        "isolated_fork_diff_scope_verified": True,
        "rows": rows,
        "scope": (
            "Implementation-safety certificate only. A rejected backend run "
            "has no mathematical meaning for the Coxeter completion problem."
        ),
        "source_sha256": {
            str(SHARED_RECORD18.relative_to(ROOT)): sha256(
                SHARED_RECORD18.read_bytes()
            ).hexdigest(),
            str(ISOLATED_RECORD14.relative_to(ROOT)): sha256(
                ISOLATED_RECORD14.read_bytes()
            ).hexdigest(),
        },
    }


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2, sort_keys=True))
