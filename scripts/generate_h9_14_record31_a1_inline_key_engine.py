#!/usr/bin/env python3
"""Generate a fixed-width quotient-key variant of the compact a=1 engine.

The compact Record-31 engine serializes a canonical colored-bipartite graph
into a string of at most 22 bytes.  ``std::unordered_set<std::string>`` still
stores a node-sized string object for every quotient state and leaves its
allocation policy to the standard library.  This generator changes only the
production quotient memo key to three inline ``uint64_t`` words.

The words are value-initialized to zero before the compact bytes are packed.
The compact code's first byte is its canonical row count, and that row count
determines the unpadded byte length.  Consequently zero padding cannot merge
two valid compact codes.  Equality compares all three words exactly; hash
collisions therefore affect performance only, never quotient semantics.

Matrix-quotient keys and all target/family canonical keys deliberately remain
``std::string``.  The live production source and binary are not modified by
this generator.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from scripts.generate_h9_14_record31_a1_addition_engine import replace_exact
from scripts.generate_h9_14_record31_a1_compact_key_engine import generate_compact


WORDS_ALIAS = "using Words = std::vector<std::uint64_t>;\n"

INLINE_KEY_TYPES = r"""using Words = std::vector<std::uint64_t>;

struct InlineQuotientKey {
  std::array<std::uint64_t, 3> words{};

  explicit InlineQuotientKey(const std::string& code) {
    if (code.empty() || code.size() > 22) {
      throw std::runtime_error("compact quotient key length changed");
    }
    const std::size_t row_count =
        static_cast<unsigned char>(code.front());
    const std::size_t expected_size = 1 + (row_count * kN + 7) / 8;
    if (code.size() != expected_size) {
      throw std::runtime_error("compact quotient row-count prefix changed");
    }
    for (std::size_t index = 0; index < code.size(); ++index) {
      words[index / 8] |=
          std::uint64_t{static_cast<unsigned char>(code[index])}
          << (8 * (index % 8));
    }
  }
};

static_assert(sizeof(InlineQuotientKey) == 3 * sizeof(std::uint64_t));

struct InlineQuotientKeyHash {
  std::size_t operator()(const InlineQuotientKey& key) const noexcept {
    std::uint64_t state = 0x243f6a8885a308d3ULL;
    for (std::uint64_t word : key.words) {
      word += 0x9e3779b97f4a7c15ULL;
      word = (word ^ (word >> 30)) * 0xbf58476d1ce4e5b9ULL;
      word = (word ^ (word >> 27)) * 0x94d049bb133111ebULL;
      word ^= word >> 31;
      state ^= word + 0x9e3779b97f4a7c15ULL +
               (state << 6) + (state >> 2);
    }
    return static_cast<std::size_t>(state);
  }
};

struct InlineQuotientKeyEqual {
  bool operator()(const InlineQuotientKey& left,
                  const InlineQuotientKey& right) const noexcept {
    return left.words == right.words;
  }
};
"""

STRING_QUOTIENT_SET = r"""  std::array<std::unordered_set<std::string>, kMaxAdditions + 1>
      quotient_seen;
"""

INLINE_QUOTIENT_SETS = r"""  std::array<
      std::unordered_set<InlineQuotientKey, InlineQuotientKeyHash,
                         InlineQuotientKeyEqual>,
      kMaxAdditions + 1>
      quotient_seen;
  // The independent matrix quotient has a different canonical code.  Keep
  // that mode byte-for-byte unchanged and separate from the inline nauty key.
  std::array<std::unordered_set<std::string>, kMaxAdditions + 1>
      matrix_quotient_seen;
"""

PROGRESS_SIZE = r"""      std::cerr << audit->quotient_seen[depth].size() << ',';
"""

INLINE_PROGRESS_SIZE = r"""      std::cerr
          << (audit->matrix_quotient
                  ? audit->matrix_quotient_seen[depth].size()
                  : audit->quotient_seen[depth].size())
          << ',';
"""

STRING_INSERT = r"""  const std::string canonical = audit->matrix_quotient
      ? CanonicalMatrixCode(*audit, state)
      : CanonicalPartialCode(*audit, state);
  if (!audit->quotient_seen[depth].insert(canonical).second) {
    ++audit->quotient_isomorphic_prunes;
    return;
  }
"""

INLINE_INSERT = r"""  bool inserted = false;
  if (audit->matrix_quotient) {
    const std::string canonical = CanonicalMatrixCode(*audit, state);
    inserted = audit->matrix_quotient_seen[depth].insert(canonical).second;
  } else {
    const InlineQuotientKey canonical(CanonicalPartialCode(*audit, state));
    inserted = audit->quotient_seen[depth].insert(canonical).second;
  }
  if (!inserted) {
    ++audit->quotient_isomorphic_prunes;
    return;
  }
"""

STRING_FINALIZE = r"""    audit.quotient_seen_final_counts[depth] = audit.quotient_seen[depth].size();
    audit.quotient_seen[depth].clear();
    audit.quotient_seen[depth].rehash(0);
"""

INLINE_FINALIZE = r"""    audit.quotient_seen_final_counts[depth] =
        audit.matrix_quotient
            ? audit.matrix_quotient_seen[depth].size()
            : audit.quotient_seen[depth].size();
    audit.quotient_seen[depth].clear();
    audit.quotient_seen[depth].rehash(0);
    audit.matrix_quotient_seen[depth].clear();
    audit.matrix_quotient_seen[depth].rehash(0);
"""


def generate_inline() -> str:
    """Return the compact engine with exact 24-byte quotient memo keys."""

    text = generate_compact()
    text = replace_exact(text, WORDS_ALIAS, INLINE_KEY_TYPES)
    text = replace_exact(text, STRING_QUOTIENT_SET, INLINE_QUOTIENT_SETS)
    text = replace_exact(text, PROGRESS_SIZE, INLINE_PROGRESS_SIZE)
    text = replace_exact(text, STRING_INSERT, INLINE_INSERT)
    text = replace_exact(text, STRING_FINALIZE, INLINE_FINALIZE)
    return text


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    arguments = parser.parse_args()
    arguments.output.write_text(generate_inline(), encoding="utf-8")


if __name__ == "__main__":
    main()
