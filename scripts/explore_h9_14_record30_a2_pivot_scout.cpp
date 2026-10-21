// Exact pivoted compact-incidence scout for all 21 positive-A4 orbit
// representatives in H9/14 Record 30's a=2 schedule.
//
// Input on stdin:
//   21
//   orbit_index a4_mask_0 a4_mask_1
//
// The schedule requires 20 minimal A5 nonfaces and exactly 128 facets.
// The search reconstructs the local Macaulay intervals, pivots on a smallest
// uncovered transversal of order at most four, memoizes unordered addition
// sets, and audits order-five and (when needed) complete blockers literally.

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

using Mask = std::uint16_t;
constexpr int kN = 14;
constexpr Mask kFull = (Mask{1} << kN) - 1;
constexpr int kExpectedOrbitCount = 21;
constexpr int kScheduledFives = 20;
constexpr int kScheduledFacets = 128;
constexpr std::array<Mask, 3> kEdges = {
    (1u << 0) | (1u << 11),
    (1u << 1) | (1u << 12),
    (1u << 2) | (1u << 13),
};
constexpr std::array<int, kN> kBaseG3 = {
    4, 4, 4, 8, 8, 8, 8, 8, 8, 8, 8, 4, 4, 4,
};
constexpr std::array<int, kN> kBaseG4 = kBaseG3;

struct Audit {
  int index = -1;
  std::array<Mask, 2> fours{};
  std::vector<Mask> shadow;
  std::vector<Mask> candidates;
  std::vector<Mask> small;
  std::vector<std::uint64_t> coverage;
  std::vector<std::vector<int>> coverers;
  std::vector<Mask> graph_five_covers;
  std::array<int, kN> lower{};
  std::array<int, kN> upper{};
  int additions = 0;
  std::array<std::set<std::uint64_t>, 6> seen;
  std::uint64_t nodes = 0;
  std::uint64_t dead_uncoverable = 0;
  std::uint64_t no_small = 0;
  std::map<int, std::uint64_t> order_five_histogram;
  std::map<int, std::uint64_t> scheduled_blocker_count_histogram;
  std::uint64_t pure_targets = 0;
};

long long Binomial(int n, int k) {
  if (k < 0 || k > n) return 0;
  k = std::min(k, n - k);
  long long answer = 1;
  for (int index = 1; index <= k; ++index) {
    answer = answer * (n - k + index) / index;
  }
  return answer;
}

int MacaulayRaise(int value, int degree) {
  int remaining = value;
  int next_top = -1;
  int answer = 0;
  for (int lower = degree; lower >= 1; --lower) {
    int top = lower - 1;
    while ((next_top < 0 || top + 1 < next_top) &&
           Binomial(top + 1, lower) <= remaining) {
      ++top;
    }
    remaining -= static_cast<int>(Binomial(top, lower));
    answer += static_cast<int>(Binomial(top + 1, lower + 1));
    next_top = top;
  }
  if (remaining != 0) throw std::runtime_error("Macaulay expansion failed");
  return answer;
}

bool Independent(Mask mask) {
  for (const Mask edge : kEdges) {
    if ((mask & edge) == edge) return false;
  }
  return true;
}

bool ContainsFour(Mask five, const std::array<Mask, 2>& fours) {
  for (const Mask four : fours) {
    if ((five & four) == four) return true;
  }
  return false;
}

bool HitsAll(Mask mask, const std::vector<Mask>& clutter) {
  for (const Mask face : clutter) {
    if ((mask & face) == 0) return false;
  }
  return true;
}

std::vector<Mask> IndependentFives() {
  std::vector<Mask> answer;
  for (Mask mask = 0; mask <= kFull; ++mask) {
    if (std::popcount(mask) == 5 && Independent(mask)) answer.push_back(mask);
  }
  if (answer.size() != 1372) {
    throw std::runtime_error("Record-30 independent-five count changed");
  }
  return answer;
}

std::vector<Mask> GraphFiveCovers() {
  std::vector<Mask> answer;
  for (Mask mask = 0; mask <= kFull; ++mask) {
    if (std::popcount(mask) != 5) continue;
    bool cover = true;
    for (const Mask edge : kEdges) cover = cover && ((mask & edge) != 0);
    if (cover) answer.push_back(mask);
  }
  return answer;
}

std::uint64_t ChosenKey(std::vector<int> chosen) {
  std::sort(chosen.begin(), chosen.end());
  std::uint64_t key = chosen.size();
  for (const int index : chosen) {
    if (index >= 2048) throw std::runtime_error("candidate key overflow");
    key = (key << 11) | static_cast<unsigned>(index);
  }
  return key;
}

bool Compatible(const Audit& audit, Mask face,
                const std::array<std::uint8_t, kN>& incidence) {
  for (int vertex = 0; vertex < kN; ++vertex) {
    if (static_cast<int>(incidence[vertex]) +
            static_cast<int>((face >> vertex) & 1u) >
        audit.upper[vertex]) {
      return false;
    }
  }
  return true;
}

std::vector<Mask> FullMinimalTransversals(const std::vector<Mask>& clutter) {
  std::vector<Mask> answer;
  for (Mask mask = 1; mask <= kFull; ++mask) {
    if (!HitsAll(mask, clutter)) continue;
    bool minimal = true;
    for (int vertex = 0; vertex < kN; ++vertex) {
      if ((mask >> vertex) & 1u) {
        minimal = minimal &&
                  !HitsAll(mask ^ (Mask{1} << vertex), clutter);
      }
    }
    if (minimal) answer.push_back(mask);
  }
  return answer;
}

void Finish(Audit* audit, const std::vector<int>& chosen,
            const std::vector<Mask>& base) {
  ++audit->no_small;
  int order_five = 0;
  for (const Mask cover : audit->graph_five_covers) {
    bool hits = true;
    for (const Mask four : audit->fours) hits = hits && ((cover & four) != 0);
    for (const int index : chosen) {
      hits = hits && ((cover & audit->candidates[index]) != 0);
    }
    order_five += hits;
  }
  ++audit->order_five_histogram[order_five];
  if (order_five != kScheduledFacets) return;

  std::vector<Mask> clutter = base;
  for (const int index : chosen) clutter.push_back(audit->candidates[index]);
  const auto blockers = FullMinimalTransversals(clutter);
  ++audit->scheduled_blocker_count_histogram[blockers.size()];
  bool pure = blockers.size() == kScheduledFacets;
  for (const Mask blocker : blockers) {
    pure = pure && std::popcount(blocker) == 5;
  }
  audit->pure_targets += pure;
}

void Search(Audit* audit, std::vector<int> chosen, std::uint64_t covered,
            std::array<std::uint8_t, kN> incidence,
            const std::vector<Mask>& base) {
  ++audit->nodes;
  if (!audit->seen[chosen.size()].insert(ChosenKey(chosen)).second) return;

  const std::uint64_t full = audit->small.empty()
                                 ? 0
                                 : (std::uint64_t{1} << audit->small.size()) - 1;
  if (static_cast<int>(chosen.size()) == audit->additions) {
    if (covered != full) return;
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (incidence[vertex] < audit->lower[vertex] ||
          incidence[vertex] > audit->upper[vertex]) {
        return;
      }
    }
    Finish(audit, chosen, base);
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
      if (std::ranges::find(chosen, candidate) != chosen.end()) continue;
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
  if ((feasible_union & full) != full) {
    ++audit->dead_uncoverable;
    return;
  }
  if (!has_uncovered) {
    for (int candidate = 0;
         candidate < static_cast<int>(audit->candidates.size()); ++candidate) {
      if (std::ranges::find(chosen, candidate) == chosen.end() &&
          Compatible(*audit, audit->candidates[candidate], incidence)) {
        best.push_back(candidate);
      }
    }
  }
  for (const int candidate : best) {
    auto next_incidence = incidence;
    for (int vertex = 0; vertex < kN; ++vertex) {
      next_incidence[vertex] +=
          (audit->candidates[candidate] >> vertex) & 1u;
    }
    auto next = chosen;
    next.push_back(candidate);
    Search(audit, std::move(next), covered | audit->coverage[candidate],
           next_incidence, base);
  }
}

Audit Run(int index, const std::array<Mask, 2>& fours,
          const std::vector<Mask>& fives,
          const std::vector<Mask>& graph_five_covers) {
  Audit audit;
  audit.index = index;
  audit.fours = fours;
  audit.graph_five_covers = graph_five_covers;
  for (const Mask five : fives) {
    if (ContainsFour(five, fours)) audit.shadow.push_back(five);
  }
  audit.additions = kScheduledFives - audit.shadow.size();
  if (audit.additions < 0 || audit.additions > 5) {
    throw std::runtime_error("addition count changed");
  }

  for (int vertex = 0; vertex < kN; ++vertex) {
    int a4_incidence = 0;
    int forced_incidence = 0;
    for (const Mask four : fours) a4_incidence += (four >> vertex) & 1u;
    for (const Mask five : audit.shadow) {
      forced_incidence += (five >> vertex) & 1u;
    }
    const int previous = kBaseG3[vertex] - a4_incidence;
    const int before =
        kBaseG4[vertex] + 6 * a4_incidence - forced_incidence;
    audit.lower[vertex] =
        std::max(0, before - MacaulayRaise(previous, 3));
    audit.upper[vertex] = before - 1;
    if (audit.upper[vertex] < 0) {
      throw std::runtime_error("positive orbit escaped local cap");
    }
  }

  for (const Mask five : fives) {
    if (ContainsFour(five, fours)) continue;
    bool allowed = true;
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (((five >> vertex) & 1u) && audit.upper[vertex] == 0) {
        allowed = false;
      }
    }
    if (allowed) audit.candidates.push_back(five);
  }

  std::vector<Mask> base(kEdges.begin(), kEdges.end());
  base.insert(base.end(), fours.begin(), fours.end());
  base.insert(base.end(), audit.shadow.begin(), audit.shadow.end());
  for (Mask mask = 1; mask <= kFull; ++mask) {
    if (std::popcount(mask) <= 4 && HitsAll(mask, base)) {
      audit.small.push_back(mask);
    }
  }
  if (audit.small.size() > 63) {
    throw std::runtime_error("small-transversal gate escaped uint64");
  }
  audit.coverers.assign(audit.small.size(), {});
  for (int candidate = 0;
       candidate < static_cast<int>(audit.candidates.size()); ++candidate) {
    std::uint64_t coverage = 0;
    for (int row = 0; row < static_cast<int>(audit.small.size()); ++row) {
      if ((audit.candidates[candidate] & audit.small[row]) == 0) {
        coverage |= std::uint64_t{1} << row;
        audit.coverers[row].push_back(candidate);
      }
    }
    audit.coverage.push_back(coverage);
  }
  Search(&audit, {}, 0, {}, base);
  return audit;
}

void PrintHistogram(const std::map<int, std::uint64_t>& histogram) {
  std::cout << '{';
  bool first = true;
  for (const auto& [key, value] : histogram) {
    if (!first) std::cout << ',';
    std::cout << '\"' << key << "\":" << value;
    first = false;
  }
  std::cout << '}';
}

}  // namespace

int main() {
  try {
    int count = 0;
    if (!(std::cin >> count) || count <= 0 || count > kExpectedOrbitCount) {
      throw std::runtime_error("expected between one and 21 representatives");
    }
    std::vector<std::pair<int, std::array<Mask, 2>>> input;
    int previous_index = -1;
    for (int row = 0; row < count; ++row) {
      int index = -1;
      std::array<Mask, 2> fours{};
      unsigned raw = 0;
      if (!(std::cin >> index) || index <= previous_index ||
          index >= kExpectedOrbitCount) {
        throw std::runtime_error("malformed positive-orbit index");
      }
      previous_index = index;
      for (Mask& four : fours) {
        if (!(std::cin >> raw) || raw > kFull || std::popcount(raw) != 4) {
          throw std::runtime_error("malformed positive A4 mask");
        }
        four = static_cast<Mask>(raw);
      }
      input.push_back({index, fours});
    }

    const auto fives = IndependentFives();
    const auto graph_five_covers = GraphFiveCovers();
    std::vector<Audit> audits;
    std::uint64_t total_targets = 0;
    for (const auto& [index, fours] : input) {
      audits.push_back(Run(index, fours, fives, graph_five_covers));
      total_targets += audits.back().pure_targets;
      std::cerr << "orbit " << index << " add=" << audits.back().additions
                << " no-small=" << audits.back().no_small
                << " targets=" << audits.back().pure_targets << '\n';
    }

    std::cout << "{\"kind\":\"h9_14_record30_a2_pivot_scout\",";
    std::cout << "\"compact_only\":true,\"orbit_count\":" << audits.size()
              << ",\"orbits\":[";
    for (int row = 0; row < count; ++row) {
      if (row) std::cout << ',';
      const auto& audit = audits[row];
      std::cout << "{\"orbit_index\":" << audit.index
                << ",\"forced_A5_shadow_count\":" << audit.shadow.size()
                << ",\"additional_A5_count\":" << audit.additions
                << ",\"candidate_count\":" << audit.candidates.size()
                << ",\"raw_unordered_addition_count\":"
                << Binomial(audit.candidates.size(), audit.additions)
                << ",\"base_small_transversal_count\":" << audit.small.size()
                << ",\"search_node_count\":" << audit.nodes
                << ",\"memo_state_counts\":[";
      for (int depth = 0; depth <= audit.additions; ++depth) {
        if (depth) std::cout << ',';
        std::cout << audit.seen[depth].size();
      }
      std::cout << "],\"dead_uncoverable_count\":" << audit.dead_uncoverable
                << ",\"no_small_transversal_family_count\":" << audit.no_small
                << ",\"order_five_transversal_count_histogram\":";
      PrintHistogram(audit.order_five_histogram);
      std::cout << ",\"scheduled_order_five_blocker_count_histogram\":";
      PrintHistogram(audit.scheduled_blocker_count_histogram);
      std::cout << ",\"scheduled_pure_target_count\":" << audit.pure_targets
                << '}';
    }
    std::cout << "],\"scheduled_pure_target_count_across_orbits\":"
              << total_targets
              << ",\"scope\":\"Exact pivoted compact Record-30 a=2 incidence scout; no finite-volume inference.\"}\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
