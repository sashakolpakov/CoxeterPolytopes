#!/usr/bin/env python3
"""Generate the Record-31 compact a=1 C++ census from the frozen a=2 engine.

The search logic is deliberately shared with the fully replayed a=2 engine.
Every schedule- and row-colour-dependent edit is an exact textual assertion:
generation aborts if the frozen template changes at any expected site.
"""

from __future__ import annotations

import argparse
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "scripts/explore_h9_14_record31_a2_additions.cpp"


def replace_exact(text: str, old: str, new: str, count: int = 1) -> str:
    actual = text.count(old)
    if actual != count:
        raise ValueError(
            f"a=1 engine template mismatch: expected {count} occurrence(s) "
            f"of {old!r}, found {actual}"
        )
    return text.replace(old, new)


def generate() -> str:
    text = TEMPLATE.read_text(encoding="utf-8")

    # Schedule constants and the number of positive-A4 rows.
    for old, new in (
        ("constexpr int kExpectedOrbitCount = 49;", "constexpr int kExpectedOrbitCount = 6;"),
        (
            "constexpr int kScheduledFives = 21;",
            "constexpr int kScheduledFours = 1;\nconstexpr int kScheduledFives = 15;",
        ),
        ("constexpr int kScheduledFacets = 118;", "constexpr int kScheduledFacets = 120;"),
        ("constexpr int kMaxAdditions = 6;", "constexpr int kMaxAdditions = 7;"),
    ):
        text = replace_exact(text, old, new)
    text = replace_exact(
        text,
        "//   orbit_index a4_mask_0 a4_mask_1",
        "//   orbit_index a4_mask_0",
    )

    # The a=1 source family has one A4 rather than two.  Keep A4 and A5 rows
    # in different nauty/matrix colour cells exactly as in the a=2 proof.
    text = replace_exact(text, "std::array<Mask, 2>", "std::array<Mask, 1>", 5)

    independent_block = """bool Independent(Mask value) {
  for (const Mask edge : kEdges) {
    if ((value & edge) == edge) return false;
  }
  return true;
}
"""
    derived_ridges = independent_block + """
int DerivedScheduledRidgeCount() {
  // Derive f_7 from the graph/A4 prefix and Dehn--Sommerville symmetry.
  // The leading entry is f_{-1}; a missing A4 removes one independent
  // four-set from f_3.
  constexpr int dimension = 9;
  std::array<long long, 5> face_prefix{};
  face_prefix[0] = 1;
  for (Mask face = 1; face <= kFull; ++face) {
    const int order = std::popcount(face);
    if (order <= 4 && Independent(face)) ++face_prefix[order];
  }
  face_prefix[4] -= kScheduledFours;

  std::array<long long, dimension + 1> h{};
  for (int index = 0; index <= 4; ++index) {
    for (int order = 0; order <= index; ++order) {
      h[index] += (((index - order) & 1) ? -1 : 1) *
          Binomial(dimension - order, dimension - index) *
          face_prefix[order];
    }
    h[dimension - index] = h[index];
  }
  long long facets = 0;
  for (const long long entry : h) facets += entry;
  if (facets != kScheduledFacets) {
    throw std::runtime_error("a=1 Dehn--Sommerville facet count changed");
  }
  long long ridges = 0;
  for (int index = 0; index <= 8; ++index) {
    ridges += Binomial(dimension - index, dimension - 8) * h[index];
  }
  return static_cast<int>(ridges);
}
"""
    text = replace_exact(text, independent_block, derived_ridges)
    replacements = (
        ("family.size() < 2 || family.size() > 8", "family.size() < 1 || family.size() > 8", 1),
        ("const int five_start = four_start + 2;", "const int five_start = four_start + 1;", 1),
        ("static_cast<int>(family.size()) - 2", "static_cast<int>(family.size()) - 1", 1),
        ("row < 2 ? four_start + row : five_start + row - 2", "row < 1 ? four_start + row : five_start + row - 1", 1),
        ("family.begin(), family.begin() + 2", "family.begin(), family.begin() + 1", 2),
        ("family.begin() + 2, family.end()", "family.begin() + 1, family.end()", 2),
        ("EqualPatternOrders(rows, 0, 2)", "EqualPatternOrders(rows, 0, 1)", 1),
        ("EqualPatternOrders(rows, 2, rows.size())", "EqualPatternOrders(rows, 1, rows.size())", 1),
        ("audit.additions < 3", "audit.additions < 5", 1),
    )
    for old, new, count in replacements:
        text = replace_exact(text, old, new, count)

    old_input = """    std::vector<std::pair<int, std::array<Mask, 1>>> input;
    int previous = -1;
    for (int row = 0; row < count; ++row) {
      int index = -1;
      unsigned left = 0, right = 0;
      if (!(std::cin >> index >> left >> right) || index <= previous ||
          index >= kExpectedOrbitCount || left > kFull || right > kFull ||
          std::popcount(left) != 4 || std::popcount(right) != 4) {
        throw std::runtime_error("malformed positive-A4 representative input");
      }
      previous = index;
      input.push_back({index, {static_cast<Mask>(left), static_cast<Mask>(right)}});
    }
"""
    new_input = """    std::vector<std::pair<int, std::array<Mask, 1>>> input;
    int previous = -1;
    for (int row = 0; row < count; ++row) {
      int index = -1;
      unsigned only = 0;
      if (!(std::cin >> index >> only) || index <= previous ||
          index >= kExpectedOrbitCount || only > kFull ||
          std::popcount(only) != 4) {
        throw std::runtime_error("malformed positive-A4 representative input");
      }
      previous = index;
      input.push_back({index, {static_cast<Mask>(only)}});
    }
"""
    text = replace_exact(text, old_input, new_input)

    old_print = """      PrintMask(audit.fours[0]);
      std::cout << ',';
      PrintMask(audit.fours[1]);
      std::cout << "]"
"""
    new_print = """      PrintMask(audit.fours[0]);
      std::cout << "]"
"""
    text = replace_exact(text, old_print, new_print)

    text = replace_exact(
        text,
        "  if (degree_histogram != std::map<int, int>{{2, 531}}) return;",
        "  static const int derived_scheduled_ridges = DerivedScheduledRidgeCount();\n"
        "  if (degree_histogram !=\n"
        "      std::map<int, int>{{2, derived_scheduled_ridges}}) return;",
    )
    target_insert = """  auto [iterator, inserted] = audit->ridge_valid_orbits.try_emplace(
      code, std::make_pair(family, std::uint64_t{0}));
  ++iterator->second.second;
"""
    target_insert_with_log = """  auto [iterator, inserted] = audit->ridge_valid_orbits.try_emplace(
      code, std::make_pair(family, std::uint64_t{0}));
  ++iterator->second.second;
  if (inserted) {
    std::cerr << "TARGET orbit=" << audit->index;
    for (const Mask face : family) std::cerr << ' ' << face;
    std::cerr << '\\n';
  }
"""
    text = replace_exact(text, target_insert, target_insert_with_log)

    for old, new, count in (
        ("expected between one and 49 representatives", "expected between one and 6 representatives", 1),
        ("record31_a2", "record31_a1", 2),
        ("pure_118_facet", "pure_120_facet", 2),
        ("Record-31 a=2", "Record-31 a=1", 2),
        ("Record 31 a=2", "Record 31 a=1", 1),
        ("compact a=2 schedule", "compact a=1 schedule", 1),
    ):
        text = replace_exact(text, old, new, count)

    root_branch = """  if (audit->matrix_quotient) std::reverse(branches.begin(), branches.end());
  for (const int candidate : branches) {
"""
    root_branch_with_log = """  if (audit->matrix_quotient) std::reverse(branches.begin(), branches.end());
  if (depth == 0) {
    std::cerr << "orbit " << audit->index
              << " root-branch-count=" << branches.size() << '\\n';
  }
  for (const int candidate : branches) {
"""
    text = replace_exact(text, root_branch, root_branch_with_log)

    return text


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.write_text(generate(), encoding="utf-8")


if __name__ == "__main__":
    main()
