// Exact addition gate for the exceptional a=4 orbit in the record-18 scout.
//
// The four fixed missing A4s have forced A5 shadow 28 and require four
// additional minimal A5s.  This program derives all candidates, all base
// transversals of order at most four, and exhaustively searches the exact
// four-row cover problem subject to the strict vertex-link incidence caps.
// It is a compact-only necessary-incidence calculation, not a theorem.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <unordered_set>
#include <vector>

namespace {

constexpr int kN = 14;
constexpr std::uint16_t kW = (1u << 12) - 1;
constexpr std::array<std::uint16_t, 5> kEdges = {
    (1u << 0) | (1u << 12),
    (1u << 1) | (1u << 12),
    (1u << 2) | (1u << 13),
    (1u << 3) | (1u << 13),
    (1u << 12) | (1u << 13),
};
constexpr std::array<std::uint16_t, 4> kFours = {
    (1u << 0) | (1u << 2) | (1u << 10) | (1u << 11),
    (1u << 1) | (1u << 2) | (1u << 10) | (1u << 11),
    (1u << 0) | (1u << 3) | (1u << 10) | (1u << 11),
    (1u << 1) | (1u << 3) | (1u << 10) | (1u << 11),
};
constexpr std::array<int, kN> kAdditionLower = {
    0,0,0,0,0,0,0,0,0,0,1,1,0,0};
constexpr std::array<int, kN> kAdditionUpper = {
    2,2,2,2,2,2,2,2,2,2,2,2,0,0};

std::vector<std::uint16_t> g_shadow;
std::vector<std::uint16_t> g_candidates;
std::vector<std::uint16_t> g_small_transversals;
std::vector<std::uint64_t> g_coverage;
std::vector<std::vector<int>> g_coverers;
std::array<std::unordered_set<std::uint64_t>, 5> g_seen;
std::set<std::array<int, 4>> g_solutions;
std::uint64_t g_nodes = 0;
std::uint64_t g_dead_uncoverable = 0;

bool HitsAll(std::uint16_t value, const std::vector<std::uint16_t>& edges) {
  return std::ranges::all_of(edges, [value](auto edge) { return value & edge; });
}

std::uint64_t ChosenKey(std::vector<int> chosen) {
  std::sort(chosen.begin(), chosen.end());
  std::uint64_t answer = chosen.size();
  for (const int index : chosen) answer = (answer << 10) | index;
  return answer;
}

bool Contains(const std::vector<int>& chosen, int index) {
  return std::find(chosen.begin(), chosen.end(), index) != chosen.end();
}

void Search(std::vector<int> chosen, std::uint64_t covered,
            std::uint16_t once, std::uint16_t twice) {
  ++g_nodes;
  const auto state_key = ChosenKey(chosen);
  if (!g_seen[chosen.size()].insert(state_key).second) return;
  const std::uint64_t full = (std::uint64_t{1} << g_small_transversals.size()) - 1;
  if (chosen.size() == 4) {
    bool lower = true;
    for (int vertex = 0; vertex < kN; ++vertex) {
      lower = lower && (kAdditionLower[vertex] == 0 ||
                        ((once | twice) & (1u << vertex)));
    }
    if (covered == full && lower) {
      std::sort(chosen.begin(), chosen.end());
      g_solutions.insert({chosen[0], chosen[1], chosen[2], chosen[3]});
    }
    return;
  }

  std::vector<int> best;
  std::uint64_t feasible_coverage = covered;
  bool has_uncovered = false;
  for (int constraint = 0;
       constraint < static_cast<int>(g_small_transversals.size()); ++constraint) {
    if ((covered >> constraint) & 1) continue;
    has_uncovered = true;
    std::vector<int> feasible;
    for (const int candidate : g_coverers[constraint]) {
      if (Contains(chosen, candidate) || (g_candidates[candidate] & twice)) continue;
      feasible.push_back(candidate);
      feasible_coverage |= g_coverage[candidate];
    }
    if (feasible.empty()) {
      ++g_dead_uncoverable;
      return;
    }
    if (best.empty() || feasible.size() < best.size()) best = std::move(feasible);
  }
  if ((feasible_coverage & full) != full) {
    ++g_dead_uncoverable;
    return;
  }
  if (!has_uncovered) {
    for (int candidate = 0; candidate < static_cast<int>(g_candidates.size()); ++candidate) {
      if (!Contains(chosen, candidate) && !(g_candidates[candidate] & twice)) {
        best.push_back(candidate);
      }
    }
  }

  for (const int candidate : best) {
    const auto face = g_candidates[candidate];
    const auto next_twice = static_cast<std::uint16_t>(twice | (once & face));
    const auto next_once = static_cast<std::uint16_t>((once ^ face) & ~next_twice);
    auto next_chosen = chosen;
    next_chosen.push_back(candidate);
    Search(std::move(next_chosen), covered | g_coverage[candidate],
           next_once, next_twice);
  }
}

std::vector<std::uint16_t> MinimalTransversals(
    const std::vector<std::uint16_t>& clutter) {
  std::vector<std::uint16_t> answers;
  for (std::uint16_t value = 1; value < (1u << kN); ++value) {
    if (!HitsAll(value, clutter)) continue;
    bool minimal = true;
    for (int vertex = 0; vertex < kN; ++vertex) {
      if ((value >> vertex) & 1) {
        minimal = minimal && !HitsAll(value ^ (1u << vertex), clutter);
      }
    }
    if (minimal) answers.push_back(value);
  }
  return answers;
}

}  // namespace

int main() {
  try {
    for (std::uint16_t value = 0; value <= kW; ++value) {
      if (std::popcount(value) != 5) continue;
      if (std::ranges::any_of(kFours, [value](auto four) {
            return (value & four) == four;
          })) {
        g_shadow.push_back(value);
      } else {
        g_candidates.push_back(value);
      }
    }
    if (g_shadow.size() != 28 || g_candidates.size() != 764) {
      throw std::runtime_error("record-18 exceptional a4 shadow/candidates changed");
    }
    for (int vertex = 0; vertex < kN; ++vertex) {
      const int expected_upper = vertex < 12 ? 2 : 0;
      if (kAdditionUpper[vertex] != expected_upper) {
        throw std::runtime_error("record-18 exceptional upper interval changed");
      }
    }
    if (std::ranges::any_of(g_candidates, [](auto face) {
          return face & ((1u << 12) | (1u << 13));
        })) {
      throw std::runtime_error("record-18 exceptional candidate escaped W");
    }

    std::vector<std::uint16_t> base(kEdges.begin(), kEdges.end());
    base.insert(base.end(), kFours.begin(), kFours.end());
    base.insert(base.end(), g_shadow.begin(), g_shadow.end());
    for (std::uint16_t value = 1; value < (1u << kN); ++value) {
      if (std::popcount(value) <= 4 && HitsAll(value, base)) {
        g_small_transversals.push_back(value);
      }
    }
    if (g_small_transversals.size() != 47) {
      throw std::runtime_error("record-18 exceptional a4 small-transversal set changed");
    }
    g_coverers.assign(g_small_transversals.size(), {});
    for (int candidate = 0; candidate < static_cast<int>(g_candidates.size()); ++candidate) {
      std::uint64_t word = 0;
      for (int constraint = 0;
           constraint < static_cast<int>(g_small_transversals.size()); ++constraint) {
        if (!(g_candidates[candidate] & g_small_transversals[constraint])) {
          word |= std::uint64_t{1} << constraint;
          g_coverers[constraint].push_back(candidate);
        }
      }
      g_coverage.push_back(word);
    }

    Search({}, 0, 0, 0);

    std::map<int, int> facet_histogram;
    for (const auto& solution : g_solutions) {
      auto clutter = base;
      for (const int index : solution) clutter.push_back(g_candidates[index]);
      ++facet_histogram[MinimalTransversals(clutter).size()];
    }
    const int target_count = facet_histogram[102];
    std::cout << "{\n"
              << "  \"kind\":\"h9_14_record18_a4_exception_addition_scout\",\n"
              << "  \"fixed_A4s\":[[0,2,10,11],[1,2,10,11],[0,3,10,11],[1,3,10,11]],\n"
              << "  \"forced_A5_shadow_count\":" << g_shadow.size() << ",\n"
              << "  \"additional_A5_count\":4,\n"
              << "  \"candidate_count\":" << g_candidates.size() << ",\n"
              << "  \"addition_lower\":[0,0,0,0,0,0,0,0,0,0,1,1,0,0],\n"
              << "  \"addition_upper\":[2,2,2,2,2,2,2,2,2,2,2,2,0,0],\n"
              << "  \"raw_unordered_addition_count\":14084657251,\n"
              << "  \"base_small_transversal_count\":"
              << g_small_transversals.size() << ",\n"
              << "  \"search_node_count\":" << g_nodes << ",\n"
              << "  \"memo_state_counts\":[";
    for (int depth = 0; depth <= 4; ++depth) {
      if (depth) std::cout << ',';
      std::cout << g_seen[depth].size();
    }
    std::cout << "],\n  \"dead_uncoverable_count\":" << g_dead_uncoverable
              << ",\n  \"no_small_transversal_addition_count\":"
              << g_solutions.size() << ",\n  \"facet_count_histogram\":{";
    bool first = true;
    for (const auto& [facets, count] : facet_histogram) {
      if (!first) std::cout << ',';
      std::cout << '"' << facets << "\":" << count;
      first = false;
    }
    std::cout << "},\n  \"scheduled_102_facet_target_count\":" << target_count
              << ",\n  \"scope\":\"Necessary incidence addition gate only; no whole-shadow or finite-volume conclusion.\"\n"
              << "}\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
