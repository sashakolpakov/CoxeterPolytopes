// Exact pivoted three-A5 addition census for the three residual Record-31
// compact a=4 orbit representatives.  The other seventeen representatives
// are exhausted by explore_h9_14_record31_a4_additions.py.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

namespace {

constexpr int kN = 14;
using Mask = std::uint16_t;
constexpr Mask kFull = (Mask{1} << kN) - 1;
constexpr std::array<Mask, 4> kEdges = {
    (1u << 0) | (1u << 11),
    (1u << 1) | (1u << 12),
    (1u << 0) | (1u << 13),
    (1u << 1) | (1u << 13),
};
constexpr std::array<std::array<Mask, 4>, 3> kRepresentatives = {{
    {{3840, 5888, 7680, 9984}},   // production orbit 2
    {{3840, 7680, 11520, 13824}}, // production orbit 18
    {{6528, 7680, 10624, 13824}}, // production orbit 19
}};
constexpr std::array<int, 3> kProductionIndices = {2, 18, 19};
constexpr std::array<int, kN> kBaseG3 = {2,2,7,7,7,7,7,7,7,7,7,5,5,4};
constexpr std::array<int, kN> kBaseG4 = {2,2,8,8,8,8,8,8,8,8,8,6,6,5};

struct Audit {
  int production_index = -1;
  std::array<Mask, 4> fours{};
  std::vector<Mask> shadow;
  std::vector<Mask> candidates;
  std::vector<Mask> small;
  std::vector<std::uint64_t> coverage;
  std::vector<std::vector<int>> coverers;
  std::array<int, kN> lower{};
  std::array<int, kN> upper{};
  std::array<std::set<std::uint64_t>, 4> seen;
  std::set<std::array<Mask, 3>> solutions;
  std::uint64_t nodes = 0;
  std::uint64_t dead_uncoverable = 0;
  std::map<int, std::uint64_t> facet_histogram;
};

long long Binomial(int n, int k) {
  if (k < 0 || k > n) return 0;
  k = std::min(k, n - k);
  long long answer = 1;
  for (int i = 1; i <= k; ++i) answer = answer * (n - k + i) / i;
  return answer;
}

int MacaulayRaise(int value, int degree) {
  int remaining = value;
  int next_top = -1;
  int raised = 0;
  for (int lower = degree; lower >= 1; --lower) {
    int top = lower - 1;
    while ((next_top < 0 || top + 1 < next_top) &&
           Binomial(top + 1, lower) <= remaining) ++top;
    remaining -= Binomial(top, lower);
    raised += Binomial(top + 1, lower + 1);
    next_top = top;
  }
  if (remaining != 0) throw std::runtime_error("Macaulay expansion failed");
  return raised;
}

bool Independent(Mask value) {
  return std::ranges::none_of(kEdges, [value](Mask edge) {
    return (value & edge) == edge;
  });
}

bool HitsAll(Mask value, const std::vector<Mask>& clutter) {
  return std::ranges::all_of(clutter, [value](Mask face) {
    return (value & face) != 0;
  });
}

bool ContainsFour(Mask five, const std::array<Mask, 4>& fours) {
  return std::ranges::any_of(fours, [five](Mask four) {
    return (five & four) == four;
  });
}

std::vector<Mask> IndependentFives() {
  std::vector<Mask> answer;
  for (Mask value = 0; value <= kFull; ++value)
    if (std::popcount(value) == 5 && Independent(value)) answer.push_back(value);
  if (answer.size() != 1296) throw std::runtime_error("A5 ground set changed");
  return answer;
}

std::uint64_t ChosenKey(std::vector<int> chosen) {
  std::sort(chosen.begin(), chosen.end());
  std::uint64_t key = chosen.size();
  for (const int index : chosen) key = (key << 10) | index;
  return key;
}

bool Compatible(const Audit& audit, Mask candidate,
                const std::array<std::uint8_t, kN>& incidence) {
  for (int vertex = 0; vertex < kN; ++vertex) {
    if (static_cast<int>(incidence[vertex]) +
            static_cast<int>((candidate >> vertex) & 1u) >
        audit.upper[vertex]) return false;
  }
  return true;
}

void Search(Audit* audit, std::vector<int> chosen, std::uint64_t covered,
            std::array<std::uint8_t, kN> incidence) {
  ++audit->nodes;
  if (!audit->seen[chosen.size()].insert(ChosenKey(chosen)).second) return;
  const std::uint64_t all =
      (std::uint64_t{1} << audit->small.size()) - 1;
  if (chosen.size() == 3) {
    if (covered != all) return;
    for (int vertex = 0; vertex < kN; ++vertex)
      if (incidence[vertex] < audit->lower[vertex] ||
          incidence[vertex] > audit->upper[vertex]) return;
    std::array<Mask, 3> solution = {
        audit->candidates[chosen[0]], audit->candidates[chosen[1]],
        audit->candidates[chosen[2]]};
    std::sort(solution.begin(), solution.end());
    audit->solutions.insert(solution);
    return;
  }

  std::vector<int> best;
  std::uint64_t feasible_union = covered;
  bool has_uncovered = false;
  for (int requirement = 0;
       requirement < static_cast<int>(audit->small.size()); ++requirement) {
    if ((covered >> requirement) & 1u) continue;
    has_uncovered = true;
    std::vector<int> feasible;
    for (const int candidate : audit->coverers[requirement]) {
      if (std::find(chosen.begin(), chosen.end(), candidate) != chosen.end())
        continue;
      if (!Compatible(*audit, audit->candidates[candidate], incidence)) continue;
      feasible.push_back(candidate);
      feasible_union |= audit->coverage[candidate];
    }
    if (feasible.empty()) {
      ++audit->dead_uncoverable;
      return;
    }
    if (best.empty() || feasible.size() < best.size()) best = std::move(feasible);
  }
  if ((feasible_union & all) != all) {
    ++audit->dead_uncoverable;
    return;
  }
  if (!has_uncovered) {
    for (int candidate = 0;
         candidate < static_cast<int>(audit->candidates.size()); ++candidate) {
      if (std::find(chosen.begin(), chosen.end(), candidate) == chosen.end() &&
          Compatible(*audit, audit->candidates[candidate], incidence))
        best.push_back(candidate);
    }
  }
  for (const int candidate : best) {
    auto next_incidence = incidence;
    for (int vertex = 0; vertex < kN; ++vertex)
      next_incidence[vertex] += (audit->candidates[candidate] >> vertex) & 1u;
    auto next = chosen;
    next.push_back(candidate);
    Search(audit, std::move(next), covered | audit->coverage[candidate],
           next_incidence);
  }
}

std::vector<Mask> MinimalTransversals(const std::vector<Mask>& clutter) {
  std::vector<Mask> answers;
  for (Mask value = 1; value <= kFull; ++value) {
    if (!HitsAll(value, clutter)) continue;
    bool minimal = true;
    for (int vertex = 0; vertex < kN; ++vertex)
      if ((value >> vertex) & 1u)
        minimal = minimal && !HitsAll(value ^ (Mask{1} << vertex), clutter);
    if (minimal) answers.push_back(value);
  }
  return answers;
}

Audit Run(int slot, const std::vector<Mask>& fives) {
  Audit audit;
  audit.production_index = kProductionIndices[slot];
  audit.fours = kRepresentatives[slot];
  for (const Mask five : fives)
    if (ContainsFour(five, audit.fours)) audit.shadow.push_back(five);
  if (audit.shadow.size() != 30) throw std::runtime_error("shadow size changed");

  for (int vertex = 0; vertex < kN; ++vertex) {
    int a4_incidence = 0;
    for (const Mask four : audit.fours) a4_incidence += (four >> vertex) & 1u;
    int forced_incidence = 0;
    for (const Mask five : audit.shadow) forced_incidence += (five >> vertex) & 1u;
    const int previous = kBaseG3[vertex] - a4_incidence;
    const int before = kBaseG4[vertex] + 6 * a4_incidence - forced_incidence;
    if (previous < 1 || before < 1) throw std::runtime_error("strict gate changed");
    audit.lower[vertex] = std::max(0, before - MacaulayRaise(previous, 3));
    audit.upper[vertex] = before - 1;
  }
  int lower_sum = 0, upper_sum = 0;
  for (int v = 0; v < kN; ++v) {
    lower_sum += audit.lower[v];
    upper_sum += audit.upper[v];
  }
  if (!(lower_sum <= 15 && 15 <= upper_sum))
    throw std::runtime_error("aggregate gate changed");

  for (const Mask five : fives) {
    if (ContainsFour(five, audit.fours)) continue;
    bool allowed = true;
    for (int vertex = 0; vertex < kN; ++vertex)
      if (((five >> vertex) & 1u) && audit.upper[vertex] < 1) allowed = false;
    if (allowed) audit.candidates.push_back(five);
  }

  std::vector<Mask> base(kEdges.begin(), kEdges.end());
  base.insert(base.end(), audit.fours.begin(), audit.fours.end());
  base.insert(base.end(), audit.shadow.begin(), audit.shadow.end());
  for (Mask value = 1; value <= kFull; ++value)
    if (std::popcount(value) <= 4 && HitsAll(value, base))
      audit.small.push_back(value);
  if (audit.small.size() > 63) throw std::runtime_error("small gate escaped word");

  audit.coverers.assign(audit.small.size(), {});
  for (int candidate = 0;
       candidate < static_cast<int>(audit.candidates.size()); ++candidate) {
    std::uint64_t coverage = 0;
    for (int requirement = 0;
         requirement < static_cast<int>(audit.small.size()); ++requirement) {
      if (!(audit.candidates[candidate] & audit.small[requirement])) {
        coverage |= std::uint64_t{1} << requirement;
        audit.coverers[requirement].push_back(candidate);
      }
    }
    audit.coverage.push_back(coverage);
  }
  Search(&audit, {}, 0, {});
  for (const auto& solution : audit.solutions) {
    auto clutter = base;
    clutter.insert(clutter.end(), solution.begin(), solution.end());
    const auto blocker = MinimalTransversals(clutter);
    ++audit.facet_histogram[blocker.size()];
  }
  return audit;
}

}  // namespace

int main() {
  try {
    const auto fives = IndependentFives();
    std::vector<Audit> audits;
    for (int slot = 0; slot < 3; ++slot) {
      audits.push_back(Run(slot, fives));
      std::cerr << "orbit " << audits.back().production_index
                << " solutions=" << audits.back().solutions.size() << '\n';
    }
    if (audits[0].candidates.size() != 764 || audits[0].nodes != 9297 ||
        audits[0].solutions.size() != 0 || !audits[0].facet_histogram.empty() ||
        audits[1].candidates.size() != 462 || audits[1].nodes != 57289 ||
        audits[1].solutions.size() != 270 ||
        audits[1].facet_histogram != std::map<int,std::uint64_t>{{105,150},{109,120}} ||
        audits[2].candidates.size() != 462 || audits[2].nodes != 26503 ||
        audits[2].solutions.size() != 1121 ||
        audits[2].facet_histogram != std::map<int,std::uint64_t>{
            {94,160},{95,120},{96,340},{99,80},{100,260},{102,150},
            {104,10},{115,1}}) {
      throw std::runtime_error("Record-31 a=4 add3 census changed");
    }
    std::cout << "{\n  \"kind\":\"h9_14_record31_a4_add3_pivot_scout\",\n"
              << "  \"compact_only\":true,\n"
              << "  \"representative_audits\":[";
    for (int slot = 0; slot < 3; ++slot) {
      if (slot) std::cout << ',';
      const auto& audit = audits[slot];
      std::cout << "{\"orbit_representative_index\":" << audit.production_index
                << ",\"forced_five_shadow_count\":" << audit.shadow.size()
                << ",\"additional_minimal_A5_count\":3"
                << ",\"individual_candidate_count\":" << audit.candidates.size()
                << ",\"raw_unordered_addition_count\":"
                << Binomial(audit.candidates.size(), 3)
                << ",\"addition_lower\":[";
      for (int v = 0; v < kN; ++v) {
        if (v) std::cout << ',';
        std::cout << audit.lower[v];
      }
      std::cout << "],\"addition_upper\":[";
      for (int v = 0; v < kN; ++v) {
        if (v) std::cout << ',';
        std::cout << audit.upper[v];
      }
      std::cout << "],\"base_small_transversal_count\":" << audit.small.size()
                << ",\"search_node_count\":" << audit.nodes
                << ",\"memo_state_counts\":[";
      for (int depth = 0; depth <= 3; ++depth) {
        if (depth) std::cout << ',';
        std::cout << audit.seen[depth].size();
      }
      std::cout << "],\"dead_uncoverable_count\":" << audit.dead_uncoverable
                << ",\"no_small_transversal_addition_count\":"
                << audit.solutions.size() << ",\"facet_count_histogram\":{";
      bool first = true;
      for (const auto& [facets, count] : audit.facet_histogram) {
        if (!first) std::cout << ',';
        std::cout << '\"' << facets << "\":" << count;
        first = false;
      }
      const auto target = audit.facet_histogram.find(114);
      std::cout << "},\"scheduled_114_facet_target_count\":"
                << (target == audit.facet_histogram.end() ? 0 : target->second)
                << '}';
    }
    std::cout << "],\n  \"scope\":\"Exact literal pivoted coverage of the three Record-31 a=4 representatives requiring three A5 additions. Compact-only; no finite-volume inference.\"\n}\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
