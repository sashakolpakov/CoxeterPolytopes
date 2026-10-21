// Generic exact addition-cover driver for a fixed positive-A4 record-18 orbit.
//
// Usage: program a A4mask ... A4mask
// where there are exactly a decimal (or 0x...) masks.  The driver derives the
// forced A5 shadow, exact local addition intervals, all candidates, every base
// transversal of order at most four, the cap-feasible subset count, and every
// addition family eliminating the small transversals.  Surviving families are
// fully blocker-dualized.  This is a necessary compact-incidence scout only.

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

constexpr int kN = 14;
constexpr std::uint16_t kFull = (1u << kN) - 1;

#ifndef H9_14_ADDITION_EDGE_INITIALIZER
#define H9_14_ADDITION_EDGE_INITIALIZER                                      \
  {                                                                          \
    (1u << 0) | (1u << 12), (1u << 1) | (1u << 12),                         \
        (1u << 2) | (1u << 13), (1u << 3) | (1u << 13),                     \
        (1u << 12) | (1u << 13)                                             \
  }
#endif
#ifndef H9_14_ADDITION_BASE_G3_INITIALIZER
#define H9_14_ADDITION_BASE_G3_INITIALIZER \
  { 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 1, 1 }
#endif
#ifndef H9_14_ADDITION_BASE_G4_INITIALIZER
#define H9_14_ADDITION_BASE_G4_INITIALIZER \
  { 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 1, 1 }
#endif
#ifndef H9_14_ADDITION_KIND
#define H9_14_ADDITION_KIND "h9_14_record18_fixed_positive_addition_cover"
#endif

constexpr std::array<std::uint16_t, 5> kEdges = {
    H9_14_ADDITION_EDGE_INITIALIZER};
constexpr std::array<int, kN> kBaseG3 =
    H9_14_ADDITION_BASE_G3_INITIALIZER;
constexpr std::array<int, kN> kBaseG4 =
    H9_14_ADDITION_BASE_G4_INITIALIZER;

using Word = std::vector<std::uint64_t>;

int g_a = 0;
int g_addition_count = 0;
int g_scheduled_facets = 0;
std::vector<std::uint16_t> g_fours;
std::vector<std::uint16_t> g_shadow;
std::vector<std::uint16_t> g_candidates;
std::vector<std::uint16_t> g_base;
std::vector<std::uint16_t> g_small;
std::array<int, kN> g_lower{};
std::array<int, kN> g_upper{};
std::vector<Word> g_coverage;
std::vector<std::vector<int>> g_coverers;
std::vector<std::unordered_set<std::uint64_t>> g_seen;
std::set<std::vector<int>> g_solutions;
std::uint64_t g_nodes = 0;
std::uint64_t g_dead_uncoverable = 0;
std::uint64_t g_dead_lower = 0;

std::uint64_t Binomial(int n, int k) {
  if (k < 0 || k > n) return 0;
  k = std::min(k, n - k);
  std::uint64_t answer = 1;
  for (int i = 1; i <= k; ++i) answer = answer * (n - k + i) / i;
  return answer;
}

int MacaulayRaise(int value, int degree) {
  if (value < 0 || degree < 1) throw std::runtime_error("bad Macaulay input");
  int remaining = value;
  int next_top = -1;
  int raised = 0;
  for (int lower = degree; lower >= 1; --lower) {
    int top = lower - 1;
    while ((next_top < 0 || top + 1 < next_top) &&
           Binomial(top + 1, lower) <= static_cast<std::uint64_t>(remaining)) {
      ++top;
    }
    remaining -= static_cast<int>(Binomial(top, lower));
    raised += static_cast<int>(Binomial(top + 1, lower + 1));
    next_top = top;
  }
  if (remaining != 0) throw std::runtime_error("inexact Macaulay expansion");
  return raised;
}

bool Independent(std::uint16_t value) {
  return std::ranges::none_of(kEdges, [value](auto edge) {
    return (value & edge) == edge;
  });
}

bool HitsAll(std::uint16_t value, const std::vector<std::uint16_t>& clutter) {
  return std::ranges::all_of(clutter, [value](auto edge) { return value & edge; });
}

bool WordFull(const Word& value, const Word& full) {
  for (int block = 0; block < static_cast<int>(full.size()); ++block) {
    if ((value[block] & full[block]) != full[block]) return false;
  }
  return true;
}

void WordOr(Word& target, const Word& source) {
  for (int block = 0; block < static_cast<int>(target.size()); ++block) {
    target[block] |= source[block];
  }
}

bool Contains(const std::vector<int>& chosen, int candidate) {
  return std::find(chosen.begin(), chosen.end(), candidate) != chosen.end();
}

std::uint64_t ChosenKey(std::vector<int> chosen) {
  std::sort(chosen.begin(), chosen.end());
  std::uint64_t answer = chosen.size();
  for (const int index : chosen) {
    if (index >= 1024) throw std::runtime_error("candidate index escaped ten bits");
    answer = (answer << 10) | index;
  }
  return answer;
}

bool CapFeasible(const std::array<unsigned char, kN>& degrees,
                 std::uint16_t face) {
  for (int vertex = 0; vertex < kN; ++vertex) {
    if (((face >> vertex) & 1) && degrees[vertex] >= g_upper[vertex]) return false;
  }
  return true;
}

std::array<unsigned char, kN> AddDegree(
    std::array<unsigned char, kN> degrees, std::uint16_t face) {
  for (int vertex = 0; vertex < kN; ++vertex) {
    degrees[vertex] += (face >> vertex) & 1;
  }
  return degrees;
}

void Search(std::vector<int> chosen, Word covered,
            std::array<unsigned char, kN> degrees) {
  ++g_nodes;
  const auto state_key = ChosenKey(chosen);
  if (!g_seen[chosen.size()].insert(state_key).second) return;
  const int remaining = g_addition_count - static_cast<int>(chosen.size());
  for (int vertex = 0; vertex < kN; ++vertex) {
    if (degrees[vertex] + remaining < g_lower[vertex]) {
      ++g_dead_lower;
      return;
    }
  }
  Word full((g_small.size() + 63) / 64, ~std::uint64_t{0});
  if (!g_small.empty() && g_small.size() % 64) {
    full.back() = (std::uint64_t{1} << (g_small.size() % 64)) - 1;
  }
  if (chosen.size() == static_cast<std::size_t>(g_addition_count)) {
    bool lower = true;
    for (int vertex = 0; vertex < kN; ++vertex) {
      lower = lower && degrees[vertex] >= g_lower[vertex];
    }
    if (WordFull(covered, full) && lower) {
      std::sort(chosen.begin(), chosen.end());
      g_solutions.insert(std::move(chosen));
    }
    return;
  }

  std::vector<int> best;
  Word feasible_coverage = covered;
  bool has_uncovered = false;
  for (int constraint = 0; constraint < static_cast<int>(g_small.size()); ++constraint) {
    if ((covered[constraint / 64] >> (constraint % 64)) & 1) continue;
    has_uncovered = true;
    std::vector<int> feasible;
    for (const int candidate : g_coverers[constraint]) {
      if (Contains(chosen, candidate) ||
          !CapFeasible(degrees, g_candidates[candidate])) continue;
      feasible.push_back(candidate);
      WordOr(feasible_coverage, g_coverage[candidate]);
    }
    if (feasible.empty()) {
      ++g_dead_uncoverable;
      return;
    }
    if (best.empty() || feasible.size() < best.size()) best = std::move(feasible);
  }
  if (!WordFull(feasible_coverage, full)) {
    ++g_dead_uncoverable;
    return;
  }
  if (!has_uncovered) {
    for (int candidate = 0; candidate < static_cast<int>(g_candidates.size()); ++candidate) {
      if (!Contains(chosen, candidate) &&
          CapFeasible(degrees, g_candidates[candidate])) best.push_back(candidate);
    }
  }
  for (const int candidate : best) {
    auto next = covered;
    WordOr(next, g_coverage[candidate]);
    auto next_chosen = chosen;
    next_chosen.push_back(candidate);
    Search(std::move(next_chosen), std::move(next),
           AddDegree(degrees, g_candidates[candidate]));
  }
}

std::uint64_t EncodeDegrees(const std::array<unsigned char, kN>& degrees) {
  std::uint64_t answer = 0;
  for (int vertex = 0; vertex < kN; ++vertex) {
    answer |= std::uint64_t{degrees[vertex]} << (3 * vertex);
  }
  return answer;
}

std::array<unsigned char, kN> DecodeDegrees(std::uint64_t word) {
  std::array<unsigned char, kN> answer{};
  for (int vertex = 0; vertex < kN; ++vertex) {
    answer[vertex] = (word >> (3 * vertex)) & 7;
  }
  return answer;
}

std::uint64_t CountCapFeasible() {
  std::vector<std::unordered_map<std::uint64_t, std::uint64_t>> dp(
      g_addition_count + 1);
  dp[0][0] = 1;
  int processed = 0;
  for (const auto face : g_candidates) {
    for (int depth = std::min(processed, g_addition_count - 1); depth >= 0; --depth) {
      for (const auto& [word, multiplicity] : dp[depth]) {
        const auto degrees = DecodeDegrees(word);
        if (!CapFeasible(degrees, face)) continue;
        const auto next = EncodeDegrees(AddDegree(degrees, face));
        dp[depth + 1][next] += multiplicity;
      }
    }
    ++processed;
  }
  std::uint64_t answer = 0;
  for (const auto& [word, multiplicity] : dp[g_addition_count]) {
    const auto degrees = DecodeDegrees(word);
    bool lower = true;
    for (int vertex = 0; vertex < kN; ++vertex) {
      lower = lower && degrees[vertex] >= g_lower[vertex];
    }
    if (lower) answer += multiplicity;
  }
  return answer;
}

std::vector<std::uint16_t> MinimalTransversals(
    const std::vector<std::uint16_t>& clutter) {
  std::vector<std::uint16_t> answers;
  for (std::uint16_t value = 1; value <= kFull; ++value) {
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

void PrintMask(std::uint16_t value) {
  std::cout << '[';
  bool first = true;
  for (int vertex = 0; vertex < kN; ++vertex) {
    if (!((value >> vertex) & 1)) continue;
    if (!first) std::cout << ',';
    std::cout << vertex;
    first = false;
  }
  std::cout << ']';
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 2) throw std::runtime_error("usage: program a A4mask...");
    g_a = std::stoi(argv[1]);
    const bool fast_scheduled_gate = argc == g_a + 3 &&
        std::string(argv[argc - 1]) == "--fast-scheduled-gate";
    if (g_a < 1 || g_a > 4 ||
        (argc != g_a + 2 && !fast_scheduled_gate)) {
      throw std::runtime_error("a must lie in 1..4 and match the A4 mask count");
    }
    for (int index = 0; index < g_a; ++index) {
      const auto value = static_cast<std::uint16_t>(
          std::stoul(argv[index + 2], nullptr, 0));
      if (std::popcount(value) != 4 || !Independent(value)) {
        throw std::runtime_error("invalid independent A4 mask");
      }
      g_fours.push_back(value);
    }
    if (std::set<std::uint16_t>(g_fours.begin(), g_fours.end()).size() !=
        g_fours.size()) throw std::runtime_error("repeated A4 mask");

    std::vector<std::uint16_t> fives;
    for (std::uint16_t value = 0; value <= kFull; ++value) {
      if (std::popcount(value) == 5 && Independent(value)) fives.push_back(value);
    }
    for (const auto five : fives) {
      if (std::ranges::any_of(g_fours, [five](auto four) {
            return (five & four) == four;
          })) g_shadow.push_back(five);
    }
    g_addition_count = 8 + 6 * g_a - static_cast<int>(g_shadow.size());
    g_scheduled_facets = 110 - 2 * g_a;
    if (g_addition_count < 0 || g_addition_count > 6) {
      throw std::runtime_error("addition count escaped the supported range");
    }

    for (int vertex = 0; vertex < kN; ++vertex) {
      int a4_incidence = 0;
      for (const auto four : g_fours) a4_incidence += (four >> vertex) & 1;
      int forced = 0;
      for (const auto five : g_shadow) forced += (five >> vertex) & 1;
      const int previous = kBaseG3[vertex] - a4_incidence;
      const int before = kBaseG4[vertex] + 6 * a4_incidence - forced;
      if (previous < 0) throw std::runtime_error("negative previous local g3");
      g_lower[vertex] = std::max(0, before - MacaulayRaise(previous, 3));
      g_upper[vertex] = before - 1;
      if (g_lower[vertex] > g_upper[vertex]) {
        throw std::runtime_error(
            "empty addition interval at vertex " + std::to_string(vertex) +
            ": previous=" + std::to_string(previous) +
            " before=" + std::to_string(before) +
            " lower=" + std::to_string(g_lower[vertex]) +
            " upper=" + std::to_string(g_upper[vertex]));
      }
    }
    for (const auto five : fives) {
      if (std::find(g_shadow.begin(), g_shadow.end(), five) != g_shadow.end()) continue;
      bool admitted = true;
      for (int vertex = 0; vertex < kN; ++vertex) {
        if (((five >> vertex) & 1) && g_upper[vertex] < 1) admitted = false;
      }
      if (admitted) g_candidates.push_back(five);
    }

    g_base.assign(kEdges.begin(), kEdges.end());
    g_base.insert(g_base.end(), g_fours.begin(), g_fours.end());
    g_base.insert(g_base.end(), g_shadow.begin(), g_shadow.end());
    for (std::uint16_t value = 1; value <= kFull; ++value) {
      if (std::popcount(value) <= 4 && HitsAll(value, g_base)) g_small.push_back(value);
    }
    const int blocks = (g_small.size() + 63) / 64;
    g_coverers.assign(g_small.size(), {});
    for (int candidate = 0; candidate < static_cast<int>(g_candidates.size()); ++candidate) {
      Word word(blocks);
      for (int constraint = 0; constraint < static_cast<int>(g_small.size()); ++constraint) {
        if (!(g_candidates[candidate] & g_small[constraint])) {
          word[constraint / 64] |= std::uint64_t{1} << (constraint % 64);
          g_coverers[constraint].push_back(candidate);
        }
      }
      g_coverage.push_back(std::move(word));
    }

    const auto cap_start = std::chrono::steady_clock::now();
    const auto cap_count = CountCapFeasible();
    const auto cap_end = std::chrono::steady_clock::now();
    g_seen.resize(g_addition_count + 1);
    const auto search_start = std::chrono::steady_clock::now();
    Search({}, Word(blocks), {});
    const auto search_end = std::chrono::steady_clock::now();

    std::map<int, int> facet_histogram;
    std::map<std::string, int> blocker_size_histogram;
    std::map<int, int> order_five_histogram;
    int target_count = 0;
    std::vector<std::vector<std::uint16_t>> targets;
    std::vector<std::uint16_t> base_order_five;
    if (fast_scheduled_gate) {
      for (std::uint16_t value = 1; value <= kFull; ++value) {
        if (std::popcount(value) == 5 && HitsAll(value, g_base)) {
          base_order_five.push_back(value);
        }
      }
    }
    const auto blocker_start = std::chrono::steady_clock::now();
    for (const auto& solution : g_solutions) {
      auto clutter = g_base;
      std::vector<std::uint16_t> additions;
      for (const int index : solution) {
        clutter.push_back(g_candidates[index]);
        additions.push_back(g_candidates[index]);
      }
      if (fast_scheduled_gate) {
        const int order_five_count = std::ranges::count_if(
            base_order_five, [&additions](auto candidate) {
              return std::ranges::all_of(additions, [candidate](auto face) {
                return candidate & face;
              });
            });
        ++order_five_histogram[order_five_count];
        if (order_five_count != g_scheduled_facets) continue;
      }
      const auto blocker = MinimalTransversals(clutter);
      ++facet_histogram[blocker.size()];
      std::map<int, int> sizes;
      for (const auto value : blocker) ++sizes[std::popcount(value)];
      std::string size_key;
      for (const auto& [size, count] : sizes) {
        if (!size_key.empty()) size_key += ',';
        size_key += std::to_string(size) + ':' + std::to_string(count);
      }
      ++blocker_size_histogram[size_key];
      if (sizes == std::map<int, int>{{5, g_scheduled_facets}}) {
        ++target_count;
        targets.push_back(std::move(additions));
      }
    }
    const auto blocker_end = std::chrono::steady_clock::now();
    const auto seconds = [](auto start, auto end) {
      return std::chrono::duration<double>(end - start).count();
    };

    std::cout << "{\n  \"kind\":\"" H9_14_ADDITION_KIND "\",\n"
              << "  \"a\":" << g_a << ",\n  \"representative_A4_masks\":[";
    for (int index = 0; index < g_a; ++index) {
      if (index) std::cout << ',';
      std::cout << g_fours[index];
    }
    std::cout << "],\n  \"representative_A4s\":[";
    for (int index = 0; index < g_a; ++index) {
      if (index) std::cout << ',';
      PrintMask(g_fours[index]);
    }
    std::cout << "],\n  \"forced_A5_shadow_count\":" << g_shadow.size()
              << ",\n  \"additional_A5_count\":" << g_addition_count
              << ",\n  \"scheduled_facet_count\":" << g_scheduled_facets
              << ",\n  \"addition_lower\":[";
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (vertex) std::cout << ',';
      std::cout << g_lower[vertex];
    }
    std::cout << "],\n  \"addition_upper\":[";
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (vertex) std::cout << ',';
      std::cout << g_upper[vertex];
    }
    std::cout << "],\n  \"candidate_count\":" << g_candidates.size()
              << ",\n  \"raw_unordered_addition_count\":"
              << Binomial(g_candidates.size(), g_addition_count)
              << ",\n  \"vertex_interval_feasible_addition_count\":" << cap_count
              << ",\n  \"base_small_transversal_count\":" << g_small.size()
              << ",\n  \"search_node_count\":" << g_nodes
              << ",\n  \"memo_state_counts\":[";
    for (int depth = 0; depth <= g_addition_count; ++depth) {
      if (depth) std::cout << ',';
      std::cout << g_seen[depth].size();
    }
    std::cout << "],\n  \"dead_uncoverable_count\":" << g_dead_uncoverable
              << ",\n  \"dead_lower_count\":" << g_dead_lower
              << ",\n  \"no_small_transversal_addition_count\":" << g_solutions.size()
              << ",\n  \"cap_DP_seconds\":" << seconds(cap_start, cap_end)
              << ",\n  \"cover_search_seconds\":" << seconds(search_start, search_end)
              << ",\n  \"blocker_gate_seconds\":" << seconds(blocker_start, blocker_end)
              << ",\n  \"fast_scheduled_gate\":"
              << (fast_scheduled_gate ? "true" : "false")
              << ",\n  \"order_five_transversal_count_histogram\":{";
    bool first = true;
    for (const auto& [facets, count] : order_five_histogram) {
      if (!first) std::cout << ',';
      std::cout << '"' << facets << "\":" << count;
      first = false;
    }
    std::cout << "},\n  \"facet_count_histogram\":{";
    first = true;
    for (const auto& [facets, count] : facet_histogram) {
      if (!first) std::cout << ',';
      std::cout << '"' << facets << "\":" << count;
      first = false;
    }
    std::cout << "},\n  \"blocker_size_histogram\":{";
    first = true;
    for (const auto& [sizes, count] : blocker_size_histogram) {
      if (!first) std::cout << ',';
      std::cout << '"' << sizes << "\":" << count;
      first = false;
    }
    std::cout << "},\n  \"scheduled_pure_target_count\":" << target_count
              << ",\n  \"scheduled_target_additions\":[";
    for (int target = 0; target < static_cast<int>(targets.size()); ++target) {
      if (target) std::cout << ',';
      std::cout << '[';
      for (int row = 0; row < static_cast<int>(targets[target].size()); ++row) {
        if (row) std::cout << ',';
        PrintMask(targets[target][row]);
      }
      std::cout << ']';
    }
    std::cout << "],\n  \"scope\":\"Necessary compact incidence gate only; no topology, Coxeter, or finite-volume conclusion.\"\n}\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
  return 0;
}
