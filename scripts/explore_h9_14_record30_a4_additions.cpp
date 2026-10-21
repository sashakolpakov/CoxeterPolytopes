// Literal compact-incidence census for all seven positive-A4 orbit
// representatives in H9/14 Record 30's a=4 schedule.
//
// Input on stdin:
//   7
//   orbit_index a4_mask_0 a4_mask_1 a4_mask_2 a4_mask_3
//
// This production route traverses ordinary increasing pairs or triples of
// additional A5 masks.  No SAT/SMT solver or symmetry quotient is used after
// the frozen seven-orbit input.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

namespace {

using Mask = std::uint16_t;
constexpr int kN = 14;
constexpr Mask kFull = (Mask{1} << kN) - 1;
constexpr int kOrbitCount = 7;
constexpr int kA4Count = 4;
constexpr int kScheduledA5Count = 32;
constexpr int kScheduledFacetCount = 124;
constexpr std::array<Mask, 3> kEdges = {
    (1u << 0) | (1u << 11),
    (1u << 1) | (1u << 12),
    (1u << 2) | (1u << 13),
};
constexpr std::array<int, kN> kBaseG3 = {
    4, 4, 4, 8, 8, 8, 8, 8, 8, 8, 8, 4, 4, 4,
};
constexpr std::array<int, kN> kBaseG4 = kBaseG3;

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

bool HitsAll(Mask mask, const std::vector<Mask>& clutter) {
  return std::ranges::all_of(clutter, [mask](Mask face) {
    return (mask & face) != 0;
  });
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

struct OrbitAudit {
  int index = -1;
  std::array<Mask, kA4Count> fours{};
  int forced_shadow_count = 0;
  int additions = 0;
  int candidate_count = 0;
  long long raw_unordered_count = 0;
  int small_count = 0;
  long long interval_feasible_count = 0;
  long long no_small_count = 0;
  std::map<int, long long> order_five_histogram;
  std::map<int, long long> scheduled_blocker_count_histogram;
  long long pure_target_count = 0;
};

OrbitAudit RunOrbit(int index, const std::array<Mask, kA4Count>& fours,
                    const std::vector<Mask>& fives) {
  OrbitAudit audit;
  audit.index = index;
  audit.fours = fours;

  std::vector<Mask> shadow;
  for (const Mask five : fives) {
    if (std::ranges::any_of(fours, [five](Mask four) {
          return (five & four) == four;
        })) {
      shadow.push_back(five);
    }
  }
  audit.forced_shadow_count = shadow.size();
  audit.additions = kScheduledA5Count - audit.forced_shadow_count;
  if (audit.additions < 2 || audit.additions > 3) {
    throw std::runtime_error("Record-30 a=4 addition count changed");
  }

  std::array<int, kN> lower{};
  std::array<int, kN> upper{};
  for (int vertex = 0; vertex < kN; ++vertex) {
    int a4_incidence = 0;
    int forced_incidence = 0;
    for (const Mask four : fours) a4_incidence += (four >> vertex) & 1u;
    for (const Mask five : shadow) forced_incidence += (five >> vertex) & 1u;
    const int previous = kBaseG3[vertex] - a4_incidence;
    const int before =
        kBaseG4[vertex] + 6 * a4_incidence - forced_incidence;
    lower[vertex] = std::max(0, before - MacaulayRaise(previous, 3));
    upper[vertex] = before - 1;
    if (upper[vertex] < 0) {
      throw std::runtime_error("positive orbit escaped local cap");
    }
  }

  std::vector<Mask> candidates;
  for (const Mask five : fives) {
    if (std::ranges::find(shadow, five) != shadow.end()) continue;
    bool allowed = true;
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (((five >> vertex) & 1u) && upper[vertex] == 0) allowed = false;
    }
    if (allowed) candidates.push_back(five);
  }
  audit.candidate_count = candidates.size();
  audit.raw_unordered_count = Binomial(candidates.size(), audit.additions);

  std::vector<Mask> base(kEdges.begin(), kEdges.end());
  base.insert(base.end(), fours.begin(), fours.end());
  base.insert(base.end(), shadow.begin(), shadow.end());
  std::vector<Mask> small;
  std::vector<Mask> base_order_five;
  for (Mask mask = 1; mask <= kFull; ++mask) {
    if (std::popcount(mask) <= 4 && HitsAll(mask, base)) small.push_back(mask);
    if (std::popcount(mask) == 5 && HitsAll(mask, base)) {
      base_order_five.push_back(mask);
    }
  }
  audit.small_count = small.size();
  if (small.size() > 63) {
    throw std::runtime_error("small-transversal gate escaped uint64");
  }
  const std::uint64_t full = small.empty()
                                 ? 0
                                 : (std::uint64_t{1} << small.size()) - 1;
  std::vector<std::uint64_t> coverage(candidates.size());
  for (int candidate = 0; candidate < static_cast<int>(candidates.size());
       ++candidate) {
    for (int row = 0; row < static_cast<int>(small.size()); ++row) {
      if ((candidates[candidate] & small[row]) == 0) {
        coverage[candidate] |= std::uint64_t{1} << row;
      }
    }
  }

  auto finish = [&](const std::vector<int>& chosen) {
    std::array<int, kN> incidence{};
    std::uint64_t covered = 0;
    for (const int candidate : chosen) {
      covered |= coverage[candidate];
      for (int vertex = 0; vertex < kN; ++vertex) {
        incidence[vertex] += (candidates[candidate] >> vertex) & 1u;
      }
    }
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (incidence[vertex] < lower[vertex] ||
          incidence[vertex] > upper[vertex]) {
        return;
      }
    }
    ++audit.interval_feasible_count;
    if (covered != full) return;
    ++audit.no_small_count;

    int order_five = 0;
    for (const Mask transversal : base_order_five) {
      bool hits = true;
      for (const int candidate : chosen) {
        hits = hits && ((transversal & candidates[candidate]) != 0);
      }
      order_five += hits;
    }
    ++audit.order_five_histogram[order_five];
    if (order_five != kScheduledFacetCount) return;

    auto clutter = base;
    for (const int candidate : chosen) clutter.push_back(candidates[candidate]);
    const auto blockers = FullMinimalTransversals(clutter);
    ++audit.scheduled_blocker_count_histogram[blockers.size()];
    bool pure = blockers.size() == kScheduledFacetCount;
    for (const Mask blocker : blockers) pure = pure && std::popcount(blocker) == 5;
    audit.pure_target_count += pure;
  };

  if (audit.additions == 2) {
    for (int first = 0; first < static_cast<int>(candidates.size()); ++first) {
      for (int second = first + 1;
           second < static_cast<int>(candidates.size()); ++second) {
        finish({first, second});
      }
    }
  } else {
    for (int first = 0; first < static_cast<int>(candidates.size()); ++first) {
      for (int second = first + 1;
           second < static_cast<int>(candidates.size()); ++second) {
        for (int third = second + 1;
             third < static_cast<int>(candidates.size()); ++third) {
          finish({first, second, third});
        }
      }
    }
  }
  return audit;
}

void PrintHistogram(const std::map<int, long long>& histogram) {
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
    if (!(std::cin >> count) || count != kOrbitCount) {
      throw std::runtime_error("expected seven positive-orbit representatives");
    }
    std::vector<std::pair<int, std::array<Mask, kA4Count>>> input;
    for (int row = 0; row < count; ++row) {
      int index = -1;
      std::array<Mask, kA4Count> fours{};
      unsigned raw = 0;
      if (!(std::cin >> index) || index != row) {
        throw std::runtime_error("malformed positive-orbit index");
      }
      for (Mask& four : fours) {
        if (!(std::cin >> raw) || raw > kFull || std::popcount(raw) != 4) {
          throw std::runtime_error("malformed positive A4 mask");
        }
        four = static_cast<Mask>(raw);
      }
      input.push_back({index, fours});
    }

    const auto fives = IndependentFives();
    std::vector<OrbitAudit> audits;
    long long total_targets = 0;
    for (const auto& [index, fours] : input) {
      audits.push_back(RunOrbit(index, fours, fives));
      total_targets += audits.back().pure_target_count;
      std::cerr << "orbit " << index << " add=" << audits.back().additions
                << " no-small=" << audits.back().no_small_count
                << " targets=" << audits.back().pure_target_count << '\n';
    }

    std::cout << "{\"kind\":\"h9_14_record30_a4_literal_addition_audit\",";
    std::cout << "\"compact_only\":true,\"orbit_count\":" << audits.size()
              << ",\"orbits\":[";
    for (int row = 0; row < count; ++row) {
      if (row) std::cout << ',';
      const auto& audit = audits[row];
      std::cout << "{\"orbit_index\":" << audit.index << ",\"A4_masks\":[";
      for (int index = 0; index < kA4Count; ++index) {
        if (index) std::cout << ',';
        std::cout << audit.fours[index];
      }
      std::cout << "],\"forced_A5_shadow_count\":"
                << audit.forced_shadow_count
                << ",\"additional_A5_count\":" << audit.additions
                << ",\"candidate_count\":" << audit.candidate_count
                << ",\"raw_unordered_addition_count\":"
                << audit.raw_unordered_count
                << ",\"base_small_transversal_count\":" << audit.small_count
                << ",\"interval_feasible_family_count\":"
                << audit.interval_feasible_count
                << ",\"no_small_transversal_family_count\":"
                << audit.no_small_count
                << ",\"order_five_transversal_count_histogram\":";
      PrintHistogram(audit.order_five_histogram);
      std::cout << ",\"scheduled_order_five_blocker_count_histogram\":";
      PrintHistogram(audit.scheduled_blocker_count_histogram);
      std::cout << ",\"scheduled_pure_target_count\":"
                << audit.pure_target_count << '}';
    }
    std::cout << "],\"scheduled_pure_target_count_across_orbits\":"
              << total_targets
              << ",\"scope\":\"Literal compact Record-30 a=4 incidence audit; no finite-volume inference.\"}\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
