#!/usr/bin/env python3
"""Generate the inline Record-31 ``a=1`` engine with a depth-two partition.

The existing root-class engine can isolate one canonical depth-one class, but
some such fragments still require too much quotient-memory.  This generator
adds the mode

    --quotient --depth-two-classes ROOT BEGIN END

(and the analogous ``--matrix-quotient`` mode).  It first retains exactly the
one-based canonical depth-one class ``ROOT`` and then retains the inclusive
range ``BEGIN..END`` of canonical depth-two classes reached below that root.
Both the root and depth-two code sets record skipped classes before returning,
so a later raw representative cannot be rediscovered under a new class id.

Every fragment restarts the deeper quotient memo.  Therefore node counts and
internal target multiplicities are not additive across fragments.  The union
of the canonical ridge-target family keys is nevertheless exhaustive; an
empty union is an exact incidence no-go for the partitioned source.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from scripts.generate_h9_14_record31_a1_addition_engine import replace_exact
from scripts.generate_h9_14_record31_a1_root_class_engine import (
    generate_root_class_engine,
)


ROOT_GLOBALS = """int gRootClassBegin = 1;
int gRootClassEnd = 1000000;
bool gRootClassRestricted = false;
"""

DEPTH_TWO_GLOBALS = """int gRootClassBegin = 1;
int gRootClassEnd = 1000000;
bool gRootClassRestricted = false;
int gDepthTwoClassBegin = 1;
int gDepthTwoClassEnd = 1000000;
bool gDepthTwoClassRestricted = false;
"""

ROOT_AUDIT_FIELDS = """  std::unordered_set<std::string> root_class_codes;
  int root_class_discovered = 0;
  int root_class_selected = 0;
  int root_class_skipped = 0;
  std::uint64_t nodes = 0;
"""

DEPTH_TWO_AUDIT_FIELDS = """  std::unordered_set<std::string> root_class_codes;
  int root_class_discovered = 0;
  int root_class_selected = 0;
  int root_class_skipped = 0;
  std::unordered_set<std::string> depth_two_class_codes;
  int depth_two_class_discovered = 0;
  int depth_two_class_selected = 0;
  int depth_two_class_skipped = 0;
  std::uint64_t nodes = 0;
"""

ROOT_FILTER = r"""  if (depth == 1) {
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
"""

DEPTH_TWO_FILTER = r"""  if (depth == 1) {
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
  if (depth == 2 && gDepthTwoClassRestricted) {
    if (!audit->depth_two_class_codes.insert(canonical_code).second) {
      ++audit->quotient_isomorphic_prunes;
      return;
    }
    const int depth_two_class_id = ++audit->depth_two_class_discovered;
    if (depth_two_class_id < gDepthTwoClassBegin ||
        depth_two_class_id > gDepthTwoClassEnd) {
      ++audit->depth_two_class_skipped;
      return;
    }
    ++audit->depth_two_class_selected;
  }

  bool inserted = false;
"""

ROOT_PARSER_TAIL = r"""    } else if (argc != 1) {
      throw std::runtime_error(
          "usage: explore_record31_a1 [--quotient|--matrix-quotient] "
          "[--root-classes BEGIN END]");
    }
"""

DEPTH_TWO_PARSER_TAIL = r"""    } else if (argc == 6 &&
               (std::string(argv[1]) == "--quotient" ||
                std::string(argv[1]) == "--matrix-quotient") &&
               std::string(argv[2]) == "--depth-two-classes") {
      quotient_mode = true;
      matrix_quotient = std::string(argv[1]) == "--matrix-quotient";
      gRootClassBegin = std::stoi(argv[3]);
      gRootClassEnd = gRootClassBegin;
      gDepthTwoClassBegin = std::stoi(argv[4]);
      gDepthTwoClassEnd = std::stoi(argv[5]);
      if (gRootClassBegin < 1 || gDepthTwoClassBegin < 1 ||
          gDepthTwoClassEnd < gDepthTwoClassBegin) {
        throw std::runtime_error("invalid depth-two class partition");
      }
      gRootClassRestricted = true;
      gDepthTwoClassRestricted = true;
    } else if (argc != 1) {
      throw std::runtime_error(
          "usage: explore_record31_a1 [--quotient|--matrix-quotient] "
          "[--root-classes BEGIN END|--depth-two-classes ROOT BEGIN END]");
    }
"""

ROOT_OUTPUT = r"""                << ",\"root_class_skipped_count\":"
                << audit.root_class_skipped
                << ",\"no_small_transversal_count\":" << audit.no_small
"""

DEPTH_TWO_OUTPUT = r"""                << ",\"root_class_skipped_count\":"
                << audit.root_class_skipped
                << ",\"depth_two_class_range_restricted\":"
                << (gDepthTwoClassRestricted ? "true" : "false")
                << ",\"requested_depth_two_class_begin\":"
                << gDepthTwoClassBegin
                << ",\"requested_depth_two_class_end\":"
                << gDepthTwoClassEnd
                << ",\"depth_two_class_discovered_count\":"
                << audit.depth_two_class_discovered
                << ",\"depth_two_class_selected_count\":"
                << audit.depth_two_class_selected
                << ",\"depth_two_class_skipped_count\":"
                << audit.depth_two_class_skipped
                << ",\"no_small_transversal_count\":" << audit.no_small
"""

ROOT_SCOPE = (
    'Exact necessary compact incidence census for Record 31 a=1; a restricted '
    'root-class run is an exhaustive overcover of that fragment, has '
    'nonadditive internal counts, and makes no finite-volume inference.'
)

DEPTH_TWO_SCOPE = (
    'Exact necessary compact incidence census for Record 31 a=1; a restricted '
    'depth-two class run is an exhaustive overcover of that fragment, has '
    'nonadditive internal counts, and makes no finite-volume inference.'
)


def generate_depth_two_class_engine() -> str:
    """Return the inline engine with root and depth-two class selectors."""

    text = generate_root_class_engine()
    text = replace_exact(text, ROOT_GLOBALS, DEPTH_TWO_GLOBALS)
    text = replace_exact(text, ROOT_AUDIT_FIELDS, DEPTH_TWO_AUDIT_FIELDS)
    text = replace_exact(text, ROOT_FILTER, DEPTH_TWO_FILTER)
    text = replace_exact(text, ROOT_PARSER_TAIL, DEPTH_TWO_PARSER_TAIL)
    text = replace_exact(text, ROOT_OUTPUT, DEPTH_TWO_OUTPUT)
    text = replace_exact(text, ROOT_SCOPE, DEPTH_TWO_SCOPE)
    return text


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    arguments = parser.parse_args()
    arguments.output.write_text(generate_depth_two_class_engine(), encoding="utf-8")


if __name__ == "__main__":
    main()
