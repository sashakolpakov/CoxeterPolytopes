// No-project generic addition audit for a fixed positive-A4 record-18 orbit.
//
// Usage: audit a A4mask ... A4mask
//
// The program independently derives the forced A5 shadow, local Macaulay
// intervals, candidate additions, cap-feasible count, small transversals, and
// complete blocker.  Its search differs from the production driver: candidates
// are ordered by decreasing transversal coverage, and an earliest-coverer
// pivot rule makes the recursion isomorphism-free without a chosen-set memo.
// The final blocker is generated through the three minimal vertex covers of
// the missing-edge double star rather than by scanning all masks uniformly.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace {

constexpr int kN = 14;
constexpr int kWSize = 12;
using Mask = std::uint16_t;
using Word = std::vector<std::uint64_t>;

constexpr std::array<Mask, 5> kEdges = {
    (1u << 0) | (1u << 12),
    (1u << 1) | (1u << 12),
    (1u << 2) | (1u << 13),
    (1u << 3) | (1u << 13),
    (1u << 12) | (1u << 13),
};
constexpr std::array<Mask, 3> kGraphCoverCores = {
    (1u << 12) | (1u << 13),
    (1u << 0) | (1u << 1) | (1u << 13),
    (1u << 2) | (1u << 3) | (1u << 12),
};

struct Candidate {
  Mask face;
  Word coverage;
  int coverage_count;
};

int g_a = 0;
int g_addition_count = 0;
int g_target_facets = 0;
bool g_ridge_first = false;
std::vector<Mask> g_fours;
std::vector<Mask> g_shadow;
std::vector<Mask> g_base;
std::vector<Mask> g_small;
std::vector<Candidate> g_candidates;
std::vector<Word> g_coverers;
Word g_full_coverage;
std::array<int, kN> g_lower{};
std::array<int, kN> g_upper{};
std::set<std::vector<Mask>> g_solutions;
std::uint64_t g_nodes = 0;
std::uint64_t g_accepted_paths = 0;
std::uint64_t g_dead_lower = 0;
std::uint64_t g_dead_uncoverable = 0;

std::uint64_t Binomial(int n, int k) {
  if (k < 0 || k > n) return 0;
  k = std::min(k, n - k);
  std::uint64_t value = 1;
  for (int i = 1; i <= k; ++i) {
    value = value * static_cast<std::uint64_t>(n - k + i) /
            static_cast<std::uint64_t>(i);
  }
  return value;
}

int MacaulayRaise(int value, int degree) {
  int remainder = value;
  int next_top = -1;
  int raised = 0;
  for (int lower = degree; lower >= 1; --lower) {
    int top = lower - 1;
    while ((next_top < 0 || top + 1 < next_top) &&
           Binomial(top + 1, lower) <=
               static_cast<std::uint64_t>(remainder)) {
      ++top;
    }
    remainder -= static_cast<int>(Binomial(top, lower));
    raised += static_cast<int>(Binomial(top + 1, lower + 1));
    next_top = top;
  }
  if (remainder != 0) throw std::runtime_error("Macaulay expansion failed");
  return raised;
}

bool Independent(Mask value) {
  return std::none_of(kEdges.begin(), kEdges.end(),
                      [value](Mask edge) { return (value & edge) == edge; });
}

std::vector<Mask> IndependentMasks(int cardinality) {
  std::vector<Mask> answer;
  for (Mask value = 0; value < (Mask{1} << kN); ++value) {
    if (std::popcount(value) == cardinality && Independent(value)) {
      answer.push_back(value);
    }
  }
  return answer;
}

bool HitsAll(Mask value, const std::vector<Mask>& clutter) {
  return std::all_of(clutter.begin(), clutter.end(),
                     [value](Mask edge) { return (value & edge) != 0; });
}

void SetBit(Word* word, int index) {
  (*word)[index / 64] |= std::uint64_t{1} << (index % 64);
}

void ClearBit(Word* word, int index) {
  (*word)[index / 64] &= ~(std::uint64_t{1} << (index % 64));
}

bool TestBit(const Word& word, int index) {
  return (word[index / 64] >> (index % 64)) & 1u;
}

bool ContainsAll(const Word& value, const Word& target) {
  for (int block = 0; block < static_cast<int>(target.size()); ++block) {
    if ((value[block] & target[block]) != target[block]) return false;
  }
  return true;
}

Word Intersect(const Word& left, const Word& right) {
  Word answer(left.size());
  for (int block = 0; block < static_cast<int>(left.size()); ++block) {
    answer[block] = left[block] & right[block];
  }
  return answer;
}

void UnionInto(Word* target, const Word& source) {
  for (int block = 0; block < static_cast<int>(target->size()); ++block) {
    (*target)[block] |= source[block];
  }
}

int CountBits(const Word& word) {
  int answer = 0;
  for (const std::uint64_t block : word) answer += std::popcount(block);
  return answer;
}

std::vector<int> BitIndices(const Word& word) {
  std::vector<int> answer;
  for (int block = 0; block < static_cast<int>(word.size()); ++block) {
    std::uint64_t value = word[block];
    while (value) {
      const int offset = std::countr_zero(value);
      answer.push_back(64 * block + offset);
      value &= value - 1;
    }
  }
  return answer;
}

bool CapFeasible(const std::array<std::uint8_t, kN>& degrees, Mask face) {
  for (int vertex = 0; vertex < kN; ++vertex) {
    if (((face >> vertex) & 1u) && degrees[vertex] >= g_upper[vertex]) {
      return false;
    }
  }
  return true;
}

std::array<std::uint8_t, kN> AddFace(
    std::array<std::uint8_t, kN> degrees, Mask face) {
  for (int vertex = 0; vertex < kN; ++vertex) {
    degrees[vertex] += (face >> vertex) & 1u;
  }
  return degrees;
}

Word CompatibleSubset(const Word& source,
                      const std::array<std::uint8_t, kN>& degrees) {
  Word answer(source.size());
  for (const int index : BitIndices(source)) {
    if (index < static_cast<int>(g_candidates.size()) &&
        CapFeasible(degrees, g_candidates[index].face)) {
      SetBit(&answer, index);
    }
  }
  return answer;
}

void Search(std::vector<int>* chosen, const Word& covered, const Word& allowed,
            const std::array<std::uint8_t, kN>& degrees) {
  ++g_nodes;
  const int remaining = g_addition_count - static_cast<int>(chosen->size());
  for (int vertex = 0; vertex < kN; ++vertex) {
    if (static_cast<int>(degrees[vertex]) + remaining < g_lower[vertex]) {
      ++g_dead_lower;
      return;
    }
  }
  if (remaining == 0) {
    if (!ContainsAll(covered, g_full_coverage)) return;
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (degrees[vertex] < g_lower[vertex] ||
          degrees[vertex] > g_upper[vertex]) {
        return;
      }
    }
    ++g_accepted_paths;
    std::vector<Mask> solution;
    for (const int index : *chosen) solution.push_back(g_candidates[index].face);
    std::sort(solution.begin(), solution.end());
    g_solutions.insert(std::move(solution));
    return;
  }

  int pivot = -1;
  int pivot_count = static_cast<int>(g_candidates.size()) + 1;
  Word choices(allowed.size());
  for (int requirement = 0; requirement < static_cast<int>(g_small.size());
       ++requirement) {
    if (TestBit(covered, requirement)) continue;
    const Word raw = Intersect(allowed, g_coverers[requirement]);
    const Word compatible = CompatibleSubset(raw, degrees);
    const int count = CountBits(compatible);
    if (count == 0) {
      ++g_dead_uncoverable;
      return;
    }
    if (count < pivot_count) {
      pivot = requirement;
      pivot_count = count;
      choices = compatible;
    }
  }

  // If every small transversal is already destroyed, fill the scheduled row
  // count in increasing candidate-order; otherwise use the earliest candidate
  // destroying the selected pivot.  Both cases remove earlier choices, so an
  // unordered final family has a unique recursion path.
  if (pivot < 0) choices = CompatibleSubset(allowed, degrees);
  if (CountBits(choices) == 0) {
    ++g_dead_uncoverable;
    return;
  }
  Word forbidden_earlier(allowed.size());
  for (const int index : BitIndices(choices)) {
    Word next_allowed = allowed;
    for (int block = 0; block < static_cast<int>(next_allowed.size()); ++block) {
      next_allowed[block] &= ~forbidden_earlier[block];
    }
    ClearBit(&next_allowed, index);
    Word next_covered = covered;
    UnionInto(&next_covered, g_candidates[index].coverage);
    chosen->push_back(index);
    Search(chosen, next_covered, next_allowed,
           AddFace(degrees, g_candidates[index].face));
    chosen->pop_back();
    SetBit(&forbidden_earlier, index);
  }
}

std::uint64_t EncodeDegrees(const std::array<std::uint8_t, kN>& degrees) {
  std::uint64_t word = 0;
  for (int vertex = 0; vertex < kN; ++vertex) {
    word |= std::uint64_t{degrees[vertex]} << (3 * vertex);
  }
  return word;
}

std::array<std::uint8_t, kN> DecodeDegrees(std::uint64_t word) {
  std::array<std::uint8_t, kN> answer{};
  for (int vertex = 0; vertex < kN; ++vertex) {
    answer[vertex] = static_cast<std::uint8_t>((word >> (3 * vertex)) & 7u);
  }
  return answer;
}

std::uint64_t CountCapFeasible() {
  std::vector<std::unordered_map<std::uint64_t, std::uint64_t>> states(
      g_addition_count + 1);
  states[0][0] = 1;
  int processed = 0;
  for (const Candidate& candidate : g_candidates) {
    for (int depth = std::min(processed, g_addition_count - 1); depth >= 0;
         --depth) {
      for (const auto& [word, multiplicity] : states[depth]) {
        const auto degrees = DecodeDegrees(word);
        if (!CapFeasible(degrees, candidate.face)) continue;
        const auto next = EncodeDegrees(AddFace(degrees, candidate.face));
        states[depth + 1][next] += multiplicity;
      }
    }
    ++processed;
  }
  std::uint64_t answer = 0;
  for (const auto& [word, multiplicity] : states[g_addition_count]) {
    const auto degrees = DecodeDegrees(word);
    bool lower = true;
    for (int vertex = 0; vertex < kN; ++vertex) {
      lower = lower && degrees[vertex] >= g_lower[vertex];
    }
    if (lower) answer += multiplicity;
  }
  return answer;
}

std::set<Mask> MinimalTransversalsByGraphCovers(
    const std::vector<Mask>& clutter) {
  std::set<Mask> answer;
  for (const Mask core : kGraphCoverCores) {
    for (Mask supplement = 0; supplement < (Mask{1} << kWSize); ++supplement) {
      const Mask value = core | supplement;
      if (!HitsAll(value, clutter)) continue;
      bool minimal = true;
      for (int vertex = 0; vertex < kN; ++vertex) {
        if ((value >> vertex) & 1u) {
          minimal = minimal &&
                    !HitsAll(value ^ (Mask{1} << vertex), clutter);
        }
      }
      if (minimal) answer.insert(value);
    }
  }
  return answer;
}

std::array<int, kN> BaseLinkG(int which) {
  std::array<std::vector<Mask>, 6> independent;
  for (int order = 2; order <= 5; ++order) {
    independent[order] = IndependentMasks(order);
  }
  std::array<int, kN> answer{};
  for (int vertex = 0; vertex < kN; ++vertex) {
    std::array<int, 6> count{};
    for (int order = 2; order <= 5; ++order) {
      for (const Mask face : independent[order]) {
        count[order] += (face >> vertex) & 1u;
      }
    }
    const int f0 = count[2];
    const int f1 = count[3];
    const int f2 = count[4];
    const int f3 = count[5];
    answer[vertex] = which == 3
                         ? f2 - 7 * f1 + 28 * f0 - 84
                         : f3 - 6 * f2 + 21 * f1 - 56 * f0 + 126;
  }
  return answer;
}

void Initialize() {
  const std::vector<Mask> fives = IndependentMasks(5);
  for (const Mask five : fives) {
    if (std::any_of(g_fours.begin(), g_fours.end(),
                    [five](Mask four) { return (five & four) == four; })) {
      g_shadow.push_back(five);
    }
  }
  g_addition_count = 8 + 6 * g_a - static_cast<int>(g_shadow.size());
  g_target_facets = 110 - 2 * g_a;
  if (g_addition_count < 0 || g_addition_count > 6) {
    throw std::runtime_error("addition count escaped 0..6");
  }

  const auto base_g3 = BaseLinkG(3);
  const auto base_g4 = BaseLinkG(4);
  for (int vertex = 0; vertex < kN; ++vertex) {
    int q = 0;
    int shadow_incidence = 0;
    for (const Mask four : g_fours) q += (four >> vertex) & 1u;
    for (const Mask five : g_shadow) {
      shadow_incidence += (five >> vertex) & 1u;
    }
    const int previous = base_g3[vertex] - q;
    const int before = base_g4[vertex] + 6 * q - shadow_incidence;
    g_lower[vertex] = std::max(0, before - MacaulayRaise(previous, 3));
    g_upper[vertex] = before - 1;
    if (g_lower[vertex] > g_upper[vertex]) {
      throw std::runtime_error("empty addition interval");
    }
  }

  const std::set<Mask> shadow_set(g_shadow.begin(), g_shadow.end());
  std::vector<Mask> candidate_faces;
  for (const Mask five : fives) {
    if (shadow_set.contains(five)) continue;
    bool admitted = true;
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (((five >> vertex) & 1u) && g_upper[vertex] < 1) admitted = false;
    }
    if (admitted) candidate_faces.push_back(five);
  }

  g_base.assign(kEdges.begin(), kEdges.end());
  g_base.insert(g_base.end(), g_fours.begin(), g_fours.end());
  g_base.insert(g_base.end(), g_shadow.begin(), g_shadow.end());
  for (int order = 1; order <= 4; ++order) {
    for (Mask value = 1; value < (Mask{1} << kN); ++value) {
      if (std::popcount(value) == order && HitsAll(value, g_base)) {
        g_small.push_back(value);
      }
    }
  }
  const int coverage_blocks = (static_cast<int>(g_small.size()) + 63) / 64;
  for (const Mask face : candidate_faces) {
    Word coverage(coverage_blocks);
    int count = 0;
    for (int requirement = 0; requirement < static_cast<int>(g_small.size());
         ++requirement) {
      if ((face & g_small[requirement]) == 0) {
        SetBit(&coverage, requirement);
        ++count;
      }
    }
    g_candidates.push_back({face, std::move(coverage), count});
  }
  std::sort(g_candidates.begin(), g_candidates.end(),
            [](const Candidate& left, const Candidate& right) {
              return std::tuple(left.coverage_count, left.face) >
                     std::tuple(right.coverage_count, right.face);
            });

  const int candidate_blocks =
      (static_cast<int>(g_candidates.size()) + 63) / 64;
  g_coverers.assign(g_small.size(), Word(candidate_blocks));
  for (int candidate = 0; candidate < static_cast<int>(g_candidates.size());
       ++candidate) {
    for (int requirement = 0; requirement < static_cast<int>(g_small.size());
         ++requirement) {
      if (TestBit(g_candidates[candidate].coverage, requirement)) {
        SetBit(&g_coverers[requirement], candidate);
      }
    }
  }
  g_full_coverage.assign(coverage_blocks, ~std::uint64_t{0});
  if (!g_small.empty() && g_small.size() % 64 != 0) {
    g_full_coverage.back() =
        (std::uint64_t{1} << (g_small.size() % 64)) - 1;
  }
}

void PrintIntArray(const std::array<int, kN>& values) {
  std::cout << '[';
  for (int vertex = 0; vertex < kN; ++vertex) {
    if (vertex) std::cout << ',';
    std::cout << values[vertex];
  }
  std::cout << ']';
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 2) throw std::runtime_error("usage: audit a A4mask...");
    int argument = 1;
    if (std::string(argv[argument]) == "--ridge-first") {
      g_ridge_first = true;
      ++argument;
    }
    if (argument >= argc) throw std::runtime_error("missing a argument");
    g_a = std::stoi(argv[argument]);
    if (g_a < 1 || g_a > 4 || argc != argument + g_a + 1) {
      throw std::runtime_error("a must lie in 1..4 and match A4 argument count");
    }
    for (int index = 0; index < g_a; ++index) {
      const Mask value =
          static_cast<Mask>(std::stoul(argv[argument + index + 1], nullptr, 0));
      if (std::popcount(value) != 4 || !Independent(value)) {
        throw std::runtime_error("invalid A4 mask");
      }
      g_fours.push_back(value);
    }
    if (std::set<Mask>(g_fours.begin(), g_fours.end()).size() !=
        g_fours.size()) {
      throw std::runtime_error("repeated A4 mask");
    }
    Initialize();
    const std::uint64_t cap_count = CountCapFeasible();

    const int candidate_blocks =
        (static_cast<int>(g_candidates.size()) + 63) / 64;
    Word allowed(candidate_blocks, ~std::uint64_t{0});
    if (!g_candidates.empty() && g_candidates.size() % 64 != 0) {
      allowed.back() =
          (std::uint64_t{1} << (g_candidates.size() % 64)) - 1;
    }
    std::vector<int> chosen;
    Search(&chosen, Word(g_full_coverage.size()), allowed, {});
    if (g_accepted_paths != g_solutions.size()) {
      throw std::runtime_error("earliest-coverer recursion produced duplicates");
    }

    std::map<int, int> facet_histogram;
    std::map<std::string, int> blocker_size_histogram;
    std::map<int, int> size_five_histogram;
    std::map<std::string, int> ridge_histogram;
    int scheduled_five_count = 0;
    int ridge_valid_count = 0;
    int target_count = 0;
    std::vector<std::vector<Mask>> pure_targets;
    for (const std::vector<Mask>& additions : g_solutions) {
      std::vector<Mask> clutter = g_base;
      clutter.insert(clutter.end(), additions.begin(), additions.end());
      if (g_ridge_first) {
        std::vector<Mask> facets;
        for (Mask value = 0; value < (Mask{1} << kN); ++value) {
          if (std::popcount(value) == 5 && HitsAll(value, clutter)) {
            facets.push_back(value);
          }
        }
        ++size_five_histogram[static_cast<int>(facets.size())];
        if (static_cast<int>(facets.size()) != g_target_facets) continue;
        ++scheduled_five_count;
        std::map<int, int> degrees;
        bool ridge_valid = true;
        for (Mask six = 0; six < (Mask{1} << kN); ++six) {
          if (std::popcount(six) != 6) continue;
          int degree = 0;
          for (const Mask facet : facets) degree += (facet & six) == facet;
          if (degree > 0) {
            ++degrees[degree];
            ridge_valid = ridge_valid && degree == 2;
          }
        }
        std::string ridge_key;
        for (const auto& [degree, count] : degrees) {
          if (!ridge_key.empty()) ridge_key += ',';
          ridge_key += std::to_string(degree) + ':' + std::to_string(count);
        }
        ++ridge_histogram[ridge_key];
        if (!ridge_valid) continue;
        ++ridge_valid_count;
      }
      const auto blocker = MinimalTransversalsByGraphCovers(clutter);
      ++facet_histogram[static_cast<int>(blocker.size())];
      std::map<int, int> sizes;
      for (const Mask value : blocker) ++sizes[std::popcount(value)];
      std::string key;
      for (const auto& [size, count] : sizes) {
        if (!key.empty()) key += ',';
        key += std::to_string(size) + ':' + std::to_string(count);
      }
      ++blocker_size_histogram[key];
      if (sizes == std::map<int, int>{{5, g_target_facets}}) {
        ++target_count;
        pure_targets.push_back(additions);
      }
    }

    std::cout << "{\n"
              << "  \"kind\":\"independent_h9_14_record18_fixed_addition_cover\",\n"
              << "  \"a\":" << g_a << ",\n"
              << "  \"representative_A4_masks\":[";
    for (int index = 0; index < g_a; ++index) {
      if (index) std::cout << ',';
      std::cout << g_fours[index];
    }
    std::cout << "],\n  \"forced_A5_shadow_count\":" << g_shadow.size()
              << ",\n  \"additional_A5_count\":" << g_addition_count
              << ",\n  \"scheduled_facet_count\":" << g_target_facets
              << ",\n  \"addition_lower\":";
    PrintIntArray(g_lower);
    std::cout << ",\n  \"addition_upper\":";
    PrintIntArray(g_upper);
    std::cout << ",\n  \"candidate_count\":" << g_candidates.size()
              << ",\n  \"raw_unordered_addition_count\":"
              << Binomial(static_cast<int>(g_candidates.size()),
                          g_addition_count)
              << ",\n  \"vertex_interval_feasible_addition_count\":"
              << cap_count
              << ",\n  \"base_small_transversal_count\":"
              << g_small.size()
              << ",\n  \"earliest_pivot_search_node_count\":" << g_nodes
              << ",\n  \"dead_lower_count\":" << g_dead_lower
              << ",\n  \"dead_uncoverable_count\":" << g_dead_uncoverable
              << ",\n  \"no_small_transversal_addition_count\":"
              << g_solutions.size() << ",\n  \"facet_count_histogram\":{";
    bool first = true;
    for (const auto& [facets, count] : facet_histogram) {
      if (!first) std::cout << ',';
      std::cout << '\"' << facets << "\":" << count;
      first = false;
    }
    std::cout << "},\n  \"blocker_size_histogram\":{";
    first = true;
    for (const auto& [key, count] : blocker_size_histogram) {
      if (!first) std::cout << ',';
      std::cout << '\"' << key << "\":" << count;
      first = false;
    }
    std::cout << '}';
    if (g_ridge_first) {
      std::cout << ",\n  \"size_five_transversal_count_histogram\":{";
      first = true;
      for (const auto& [facets, count] : size_five_histogram) {
        if (!first) std::cout << ',';
        std::cout << '\"' << facets << "\":" << count;
        first = false;
      }
      std::cout << "},\n  \"scheduled_five_transversal_count\":"
                << scheduled_five_count
                << ",\n  \"scheduled_five_ridge_histogram\":{";
      first = true;
      for (const auto& [key, count] : ridge_histogram) {
        if (!first) std::cout << ',';
        std::cout << '\"' << key << "\":" << count;
        first = false;
      }
      std::cout << "},\n  \"ridge_valid_scheduled_five_count\":"
                << ridge_valid_count;
    }
    std::cout << ",\n  \"scheduled_pure_target_count\":" << target_count
              << ",\n  \"scheduled_pure_target_addition_masks\":[";
    for (int target = 0; target < static_cast<int>(pure_targets.size());
         ++target) {
      if (target) std::cout << ',';
      std::cout << '[';
      for (int row = 0; row < static_cast<int>(pure_targets[target].size());
           ++row) {
        if (row) std::cout << ',';
        std::cout << pure_targets[target][row];
      }
      std::cout << ']';
    }
    std::cout << ']'
              << ",\n  \"compact_only\":true,\n"
              << "  \"whole_shadow_compact_coxeter_theorem\":false,\n"
              << "  \"finite_volume_noncompact_chamber_excluded\":false\n"
              << "}\n";
  } catch (const std::exception& error) {
    std::cerr << "audit failed: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
