#!/usr/bin/env python3
"""Generate the inline a=1 engine with a deterministic root-class range.

The ordinary quotient discovers canonical depth-one states in a deterministic
one-based order while it scans the 196 raw root branches.  The option

    --quotient --root-classes BEGIN END

retains only canonical root classes in that inclusive range.  A separate
root-code set still records skipped classes, so later raw representatives of a
skipped class cannot be mistaken for a new selected class.

Fragments deliberately restart the deeper quotient memo.  Their node counts,
pure-target multiplicities, and histograms can therefore overcount states
reachable from two root fragments and must not be added to an interrupted
full-run transcript.  The union of canonical ridge-valid target orbits is
exhaustive.  In the current source-orbit-4 application the completed first
fragment has no target, so a zero-target complementary 4--6 fragment would be
a complete incidence no-go.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from scripts.generate_h9_14_record31_a1_addition_engine import replace_exact
from scripts.generate_h9_14_record31_a1_inline_key_engine import generate_inline


CONSTANTS = """constexpr int kMaxAdditions = 7;
"""

RANGE_GLOBALS = """constexpr int kMaxAdditions = 7;

int gRootClassBegin = 1;
int gRootClassEnd = 1000000;
bool gRootClassRestricted = false;
"""

AUDIT_FIELDS = """  bool fixed_local_failure = false;
  bool matrix_quotient = false;
  std::uint64_t nodes = 0;
"""

AUDIT_FIELDS_WITH_ROOT = """  bool fixed_local_failure = false;
  bool matrix_quotient = false;
  std::unordered_set<std::string> root_class_codes;
  int root_class_discovered = 0;
  int root_class_selected = 0;
  int root_class_skipped = 0;
  std::uint64_t nodes = 0;
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

RANGED_INSERT = r"""  const std::string canonical_code = audit->matrix_quotient
      ? CanonicalMatrixCode(*audit, state)
      : CanonicalPartialCode(*audit, state);
  if (depth == 1) {
    if (!audit->root_class_codes.insert(canonical_code).second) {
      ++audit->quotient_isomorphic_prunes;
      return;
    }
    const int root_class_id = ++audit->root_class_discovered;
    if (root_class_id < gRootClassBegin || root_class_id > gRootClassEnd) {
      ++audit->root_class_skipped;
      return;
    }
    ++audit->root_class_selected;
  }

  bool inserted = false;
  if (audit->matrix_quotient) {
    inserted =
        audit->matrix_quotient_seen[depth].insert(canonical_code).second;
  } else {
    const InlineQuotientKey canonical(canonical_code);
    inserted = audit->quotient_seen[depth].insert(canonical).second;
  }
  if (!inserted) {
    ++audit->quotient_isomorphic_prunes;
    return;
  }
"""

ARGUMENT_PARSER = r"""    bool quotient_mode = false;
    bool matrix_quotient = false;
    if (argc == 2 && std::string(argv[1]) == "--quotient") {
      quotient_mode = true;
    } else if (argc == 2 && std::string(argv[1]) == "--matrix-quotient") {
      quotient_mode = true;
      matrix_quotient = true;
    } else if (argc != 1) {
      throw std::runtime_error("usage: explore_record31_a1 [--quotient]");
    }
"""

RANGED_ARGUMENT_PARSER = r"""    bool quotient_mode = false;
    bool matrix_quotient = false;
    if (argc == 2 && std::string(argv[1]) == "--quotient") {
      quotient_mode = true;
    } else if (argc == 2 && std::string(argv[1]) == "--matrix-quotient") {
      quotient_mode = true;
      matrix_quotient = true;
    } else if (argc == 5 &&
               (std::string(argv[1]) == "--quotient" ||
                std::string(argv[1]) == "--matrix-quotient") &&
               std::string(argv[2]) == "--root-classes") {
      quotient_mode = true;
      matrix_quotient = std::string(argv[1]) == "--matrix-quotient";
      gRootClassBegin = std::stoi(argv[3]);
      gRootClassEnd = std::stoi(argv[4]);
      if (gRootClassBegin < 1 || gRootClassEnd < gRootClassBegin) {
        throw std::runtime_error("invalid inclusive root-class range");
      }
      gRootClassRestricted = true;
    } else if (argc != 1) {
      throw std::runtime_error(
          "usage: explore_record31_a1 [--quotient|--matrix-quotient] "
          "[--root-classes BEGIN END]");
    }
"""

ROW_OUTPUT = r"""      std::cout << ']'
                << ",\"no_small_transversal_count\":" << audit.no_small
"""

RANGED_ROW_OUTPUT = r"""      std::cout << ']'
                << ",\"root_class_range_restricted\":"
                << (gRootClassRestricted ? "true" : "false")
                << ",\"requested_root_class_begin\":" << gRootClassBegin
                << ",\"requested_root_class_end\":" << gRootClassEnd
                << ",\"root_class_discovered_count\":"
                << audit.root_class_discovered
                << ",\"root_class_selected_count\":"
                << audit.root_class_selected
                << ",\"root_class_skipped_count\":"
                << audit.root_class_skipped
                << ",\"no_small_transversal_count\":" << audit.no_small
"""

SCOPE_OUTPUT = r"""              << ",\"scope\":\"Exact necessary compact incidence census for Record 31 a=1; no finite-volume inference.\"}\n";
"""

RANGED_SCOPE_OUTPUT = r"""              << ",\"scope\":\"Exact necessary compact incidence census for Record 31 a=1; a restricted root-class run is an exhaustive overcover of that fragment, has nonadditive internal counts, and makes no finite-volume inference.\"}\n";
"""


def generate_root_class_engine() -> str:
    """Return the inline engine with an inclusive root-class selector."""

    text = generate_inline()
    text = replace_exact(text, CONSTANTS, RANGE_GLOBALS)
    text = replace_exact(text, AUDIT_FIELDS, AUDIT_FIELDS_WITH_ROOT)
    text = replace_exact(text, INLINE_INSERT, RANGED_INSERT)
    text = replace_exact(text, ARGUMENT_PARSER, RANGED_ARGUMENT_PARSER)
    text = replace_exact(text, ROW_OUTPUT, RANGED_ROW_OUTPUT)
    text = replace_exact(text, SCOPE_OUTPUT, RANGED_SCOPE_OUTPUT)
    return text


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    arguments = parser.parse_args()
    arguments.output.write_text(generate_root_class_engine(), encoding="utf-8")


if __name__ == "__main__":
    main()
