// Independent exceptional-a=4 addition audit for H9/14 record 18.
//
// The input is the intrinsic size-28 positive-A4 orbit representative.  The
// search is a transversal-cover recursion: an added A5 destroys a current
// small transversal precisely when it is disjoint from it.  At each node we
// pivot on an uncovered transversal and distinguish the earliest selected
// candidate that destroys it.  Candidate order is decreasing coverage size
// and then decreasing mask, intentionally different from lexicographic
// combination enumeration.  No project source or production archive is used.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace {

constexpr int kN = 14;
constexpr int kCandidateVertices = 12;
constexpr int kCandidateWordCount = 12;  // ceil(764/64)
using Mask = std::uint16_t;
using CandidateBits = std::array<std::uint64_t, kCandidateWordCount>;

constexpr std::array<Mask, 5> kMissingEdges = {
    (1u << 0) | (1u << 12),
    (1u << 1) | (1u << 12),
    (1u << 2) | (1u << 13),
    (1u << 3) | (1u << 13),
    (1u << 12) | (1u << 13),
};

constexpr std::array<Mask, 4> kMissingFours = {
    (1u << 0) | (1u << 2) | (1u << 10) | (1u << 11),
    (1u << 1) | (1u << 2) | (1u << 10) | (1u << 11),
    (1u << 0) | (1u << 3) | (1u << 10) | (1u << 11),
    (1u << 1) | (1u << 3) | (1u << 10) | (1u << 11),
};

struct Candidate {
  Mask mask;
  std::uint64_t destroyed_transversals;
};

std::vector<Mask> g_small_transversals;
std::vector<Candidate> g_candidates;
std::vector<CandidateBits> g_coverers;
std::array<int, kCandidateVertices> g_lower{};
std::array<int, kCandidateVertices> g_upper{};
std::set<std::array<Mask, 4>> g_solutions;
std::uint64_t g_nodes = 0;
std::uint64_t g_early_cover_nodes = 0;

bool HitsEvery(Mask value, const std::vector<Mask>& clutter) {
  return std::all_of(clutter.begin(), clutter.end(),
                     [value](Mask edge) { return (value & edge) != 0; });
}

bool ContainsAnyFour(Mask value) {
  return std::any_of(kMissingFours.begin(), kMissingFours.end(),
                     [value](Mask four) { return (value & four) == four; });
}

std::vector<Mask> MasksOfSize(int vertex_count, int cardinality) {
  std::vector<Mask> answer;
  for (Mask value = 0; value < (Mask{1} << vertex_count); ++value) {
    if (std::popcount(value) == cardinality) answer.push_back(value);
  }
  return answer;
}

void SetBit(CandidateBits* bits, int index) {
  (*bits)[index / 64] |= std::uint64_t{1} << (index % 64);
}

void ClearBit(CandidateBits* bits, int index) {
  (*bits)[index / 64] &= ~(std::uint64_t{1} << (index % 64));
}

CandidateBits Intersect(const CandidateBits& left, const CandidateBits& right) {
  CandidateBits answer{};
  for (int word = 0; word < kCandidateWordCount; ++word) {
    answer[word] = left[word] & right[word];
  }
  return answer;
}

int CountBits(const CandidateBits& bits) {
  int answer = 0;
  for (const std::uint64_t word : bits) answer += std::popcount(word);
  return answer;
}

std::vector<int> BitIndices(const CandidateBits& bits) {
  std::vector<int> answer;
  for (int word = 0; word < kCandidateWordCount; ++word) {
    std::uint64_t value = bits[word];
    while (value) {
      const int offset = std::countr_zero(value);
      answer.push_back(64 * word + offset);
      value &= value - 1;
    }
  }
  return answer;
}

bool Compatible(const Candidate& candidate,
                const std::array<std::uint8_t, kCandidateVertices>& incidence) {
  for (int vertex = 0; vertex < kCandidateVertices; ++vertex) {
    if (static_cast<int>(incidence[vertex]) +
            static_cast<int>((candidate.mask >> vertex) & 1u) >
        g_upper[vertex]) {
      return false;
    }
  }
  return true;
}

CandidateBits CompatibleSubset(
    const CandidateBits& bits,
    const std::array<std::uint8_t, kCandidateVertices>& incidence) {
  CandidateBits answer{};
  for (const int index : BitIndices(bits)) {
    if (index < static_cast<int>(g_candidates.size()) &&
        Compatible(g_candidates[index], incidence)) {
      SetBit(&answer, index);
    }
  }
  return answer;
}

void Search(std::vector<int>* chosen, std::uint64_t covered,
            const CandidateBits& allowed,
            const std::array<std::uint8_t, kCandidateVertices>& incidence) {
  ++g_nodes;
  const std::uint64_t all =
      (std::uint64_t{1} << g_small_transversals.size()) - 1;
  if (chosen->size() == 4) {
    if (covered != all) return;
    for (int vertex = 0; vertex < kCandidateVertices; ++vertex) {
      if (incidence[vertex] < g_lower[vertex] ||
          incidence[vertex] > g_upper[vertex]) {
        return;
      }
    }
    std::array<Mask, 4> solution{};
    for (int i = 0; i < 4; ++i) solution[i] = g_candidates[(*chosen)[i]].mask;
    std::sort(solution.begin(), solution.end());
    g_solutions.insert(solution);
    return;
  }
  if (covered == all) {
    ++g_early_cover_nodes;
    return;
  }

  // Every still-uncovered transversal must be destroyed by a future A5.
  // Choose the one with the smallest currently compatible candidate set.
  int pivot = -1;
  int pivot_count = static_cast<int>(g_candidates.size()) + 1;
  CandidateBits pivot_candidates{};
  const std::uint64_t uncovered = all & ~covered;
  for (int requirement = 0;
       requirement < static_cast<int>(g_small_transversals.size());
       ++requirement) {
    if (!((uncovered >> requirement) & 1u)) continue;
    const CandidateBits raw = Intersect(allowed, g_coverers[requirement]);
    const CandidateBits compatible = CompatibleSubset(raw, incidence);
    const int count = CountBits(compatible);
    if (count == 0) return;
    if (count < pivot_count) {
      pivot = requirement;
      pivot_count = count;
      pivot_candidates = compatible;
    }
  }
  if (pivot < 0) throw std::runtime_error("uncovered pivot was not found");

  // In each solution extending this state, distinguish its earliest candidate
  // (in the declared nonlexicographic order) which destroys the pivot.  Earlier
  // pivot-coverers are removed from later branches, making enumeration unique.
  CandidateBits forbidden_earlier{};
  for (const int index : BitIndices(pivot_candidates)) {
    CandidateBits next_allowed = allowed;
    for (int word = 0; word < kCandidateWordCount; ++word) {
      next_allowed[word] &= ~forbidden_earlier[word];
    }
    ClearBit(&next_allowed, index);
    auto next_incidence = incidence;
    for (int vertex = 0; vertex < kCandidateVertices; ++vertex) {
      next_incidence[vertex] += (g_candidates[index].mask >> vertex) & 1u;
    }
    chosen->push_back(index);
    Search(chosen, covered | g_candidates[index].destroyed_transversals,
           next_allowed, next_incidence);
    chosen->pop_back();
    SetBit(&forbidden_earlier, index);
  }
}

void Initialize() {
  std::vector<Mask> base(kMissingEdges.begin(), kMissingEdges.end());
  base.insert(base.end(), kMissingFours.begin(), kMissingFours.end());
  for (int size = 0; size <= 4; ++size) {
    for (const Mask candidate : MasksOfSize(kN, size)) {
      if (HitsEvery(candidate, base)) g_small_transversals.push_back(candidate);
    }
  }
  if (g_small_transversals.size() != 47) {
    throw std::runtime_error("exceptional base does not have 47 small transversals");
  }

  int shadow_count = 0;
  for (const Mask five : MasksOfSize(kCandidateVertices, 5)) {
    if (ContainsAnyFour(five)) {
      ++shadow_count;
      continue;
    }
    std::uint64_t coverage = 0;
    for (int index = 0; index < static_cast<int>(g_small_transversals.size());
         ++index) {
      if ((five & g_small_transversals[index]) == 0) {
        coverage |= std::uint64_t{1} << index;
      }
    }
    g_candidates.push_back({five, coverage});
  }
  if (shadow_count != 28 || g_candidates.size() != 764) {
    throw std::runtime_error("exceptional shadow or candidate count changed");
  }

  std::sort(g_candidates.begin(), g_candidates.end(),
            [](const Candidate& left, const Candidate& right) {
              return std::tuple(std::popcount(left.destroyed_transversals),
                                left.mask) >
                     std::tuple(std::popcount(right.destroyed_transversals),
                                right.mask);
            });
  g_coverers.assign(g_small_transversals.size(), CandidateBits{});
  for (int index = 0; index < static_cast<int>(g_candidates.size()); ++index) {
    for (int requirement = 0;
         requirement < static_cast<int>(g_small_transversals.size());
         ++requirement) {
      if ((g_candidates[index].destroyed_transversals >> requirement) & 1u) {
        SetBit(&g_coverers[requirement], index);
      }
    }
  }

  g_lower.fill(0);
  g_upper.fill(2);
  g_lower[10] = 1;
  g_lower[11] = 1;
}

std::map<int, int> FacetHistogram() {
  std::map<int, int> histogram;
  for (const auto& additions : g_solutions) {
    std::vector<Mask> clutter(kMissingEdges.begin(), kMissingEdges.end());
    clutter.insert(clutter.end(), kMissingFours.begin(), kMissingFours.end());
    clutter.insert(clutter.end(), additions.begin(), additions.end());
    int facets = 0;
    for (Mask candidate = 1; candidate < (Mask{1} << kN); ++candidate) {
      if (!HitsEvery(candidate, clutter)) continue;
      bool minimal = true;
      for (int vertex = 0; vertex < kN; ++vertex) {
        if ((candidate >> vertex) & 1u) {
          minimal = minimal &&
                    !HitsEvery(candidate ^ (Mask{1} << vertex), clutter);
        }
      }
      facets += minimal;
    }
    ++histogram[facets];
  }
  histogram[102] += 0;
  return histogram;
}

std::uint64_t Binomial(int n, int k) {
  k = std::min(k, n - k);
  std::uint64_t answer = 1;
  for (int i = 1; i <= k; ++i) {
    answer = answer * static_cast<std::uint64_t>(n - k + i) /
             static_cast<std::uint64_t>(i);
  }
  return answer;
}

}  // namespace

int main() {
  try {
    Initialize();
    CandidateBits allowed{};
    for (int index = 0; index < static_cast<int>(g_candidates.size()); ++index) {
      SetBit(&allowed, index);
    }
    std::vector<int> chosen;
    std::array<std::uint8_t, kCandidateVertices> incidence{};
    Search(&chosen, 0, allowed, incidence);
    if (g_early_cover_nodes != 0) {
      throw std::runtime_error("a transversal cover completed before four A5s");
    }
    if (g_solutions.size() != 3960) {
      throw std::runtime_error("exceptional no-small-transversal count changed");
    }
    const std::map<int, int> histogram = FacetHistogram();
    const std::map<int, int> expected = {
        {82, 160}, {83, 1920}, {85, 240}, {87, 120},
        {88, 900}, {92, 260},  {94, 360}, {102, 0},
    };
    if (histogram != expected) {
      throw std::runtime_error("exceptional facet histogram changed");
    }

    std::cout << "{\n"
              << "  \"kind\":\"independent_h9_14_record18_a4_exception_audit\",\n"
              << "  \"missing_four_representative\":[[0,2,10,11],[1,2,10,11],"
                 "[0,3,10,11],[1,3,10,11]],\n"
              << "  \"forced_A5_shadow_count\":28,\n"
              << "  \"addition_intervals\":{\"zero_to_two\":10,"
                 "\"one_to_two\":2,\"zero_to_zero\":2},\n"
              << "  \"additional_A5_count\":4,\n"
              << "  \"candidate_count\":" << g_candidates.size() << ",\n"
              << "  \"raw_unordered_addition_count\":"
              << Binomial(static_cast<int>(g_candidates.size()), 4) << ",\n"
              << "  \"base_small_transversal_count\":"
              << g_small_transversals.size() << ",\n"
              << "  \"pivot_search_nodes\":" << g_nodes << ",\n"
              << "  \"no_small_transversal_addition_count\":"
              << g_solutions.size() << ",\n"
              << "  \"facet_count_histogram\":{";
    bool first = true;
    for (const auto& [facets, count] : histogram) {
      if (!first) std::cout << ',';
      std::cout << '\"' << facets << "\":" << count;
      first = false;
    }
    std::cout << "},\n"
              << "  \"scheduled_102_facet_target_count\":0,\n"
              << "  \"compact_only\":true,\n"
              << "  \"whole_shadow_compact_coxeter_theorem\":false,\n"
              << "  \"finite_volume_noncompact_chamber_excluded\":false\n"
              << "}\n";
  } catch (const std::exception& error) {
    std::cerr << "audit failed: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
