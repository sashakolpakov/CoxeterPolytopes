// Exact addition census for every positive-A4 orbit in the compact
// H9/14 record-30 a=3 schedule.
//
// Input on stdin:
//   orbit_count
//   orbit_index a4_mask_0 a4_mask_1 a4_mask_2
// for each orbit.  Masks use the fourteen low bits.
//
// The program deliberately enumerates a superset of genuine polytope
// completions.  It imposes the exact vertex-link incidence intervals and
// eliminates every transversal of order at most four.  Every surviving
// addition family is then audited by literal traversal of all 2^14 masks.
// Consequently, zero pure families with the scheduled 126 facets is a safe
// compact-incidence exclusion without using topology or Coxeter labels.

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
constexpr int kA4Count = 3;
constexpr int kScheduledA5Count = 26;
constexpr int kScheduledFacetCount = 126;
constexpr std::array<std::uint16_t, 3> kMissingEdges = {
    (1u << 0) | (1u << 11),
    (1u << 1) | (1u << 12),
    (1u << 2) | (1u << 13),
};
constexpr std::array<int, kN> kBaseG3 = {
    4, 4, 4, 8, 8, 8, 8, 8, 8, 8, 8, 4, 4, 4,
};
constexpr std::array<int, kN> kBaseG4 = kBaseG3;

using Family = std::vector<std::uint16_t>;

long long Binomial(int n, int k) {
  if (k < 0 || k > n) return 0;
  k = std::min(k, n - k);
  long long value = 1;
  for (int index = 1; index <= k; ++index) {
    value = value * (n - k + index) / index;
  }
  return value;
}

int MacaulayRaise(int value, int degree) {
  int remaining = value;
  int next_top = -1;
  int raised = 0;
  for (int lower = degree; lower >= 1; --lower) {
    int top = lower - 1;
    while ((next_top < 0 || top + 1 < next_top) &&
           Binomial(top + 1, lower) <= remaining) {
      ++top;
    }
    remaining -= static_cast<int>(Binomial(top, lower));
    raised += static_cast<int>(Binomial(top + 1, lower + 1));
    next_top = top;
  }
  if (remaining != 0) throw std::runtime_error("inexact Macaulay expansion");
  return raised;
}

bool Independent(std::uint16_t mask) {
  for (const auto edge : kMissingEdges) {
    if ((mask & edge) == edge) return false;
  }
  return true;
}

std::vector<std::uint16_t> IndependentFives() {
  std::vector<std::uint16_t> answer;
  for (std::uint16_t mask = 0; mask < (1u << kN); ++mask) {
    if (std::popcount(mask) == 5 && Independent(mask)) answer.push_back(mask);
  }
  if (answer.size() != 1372) {
    throw std::runtime_error("record-30 independent-five count changed");
  }
  return answer;
}

bool HitsAll(std::uint16_t mask, const Family& clutter) {
  return std::ranges::all_of(clutter, [mask](std::uint16_t face) {
    return (mask & face) != 0;
  });
}

bool MinimalTransversal(std::uint16_t mask, const Family& clutter) {
  if (!HitsAll(mask, clutter)) return false;
  for (int vertex = 0; vertex < kN; ++vertex) {
    const auto bit = static_cast<std::uint16_t>(1u << vertex);
    if ((mask & bit) && HitsAll(mask ^ bit, clutter)) return false;
  }
  return true;
}

struct SolutionAudit {
  std::map<int, int> blocker_size_histogram;
  int blocker_count = 0;

  auto operator<=>(const SolutionAudit&) const = default;
};

SolutionAudit AuditSolution(const Family& clutter) {
  SolutionAudit audit;
  for (std::uint16_t mask = 1; mask < (1u << kN); ++mask) {
    if (MinimalTransversal(mask, clutter)) {
      ++audit.blocker_count;
      ++audit.blocker_size_histogram[std::popcount(mask)];
    }
  }
  return audit;
}

struct OrbitResult {
  int orbit_index = -1;
  std::array<std::uint16_t, kA4Count> fours{};
  int forced_shadow_count = 0;
  int addition_count = 0;
  int candidate_count = 0;
  long long raw_addition_count = 0;
  int small_transversal_count = 0;
  std::set<Family> no_small_families;
  std::map<int, int> order_five_transversal_histogram;
  std::map<SolutionAudit, int> scheduled_order_five_audit_histogram;
  std::vector<Family> scheduled_pure_targets;
};

class AdditionEnumerator {
 public:
  AdditionEnumerator(std::vector<std::uint16_t> candidates,
                     std::vector<std::uint16_t> small_transversals,
                     std::array<int, kN> lower,
                     std::array<int, kN> upper,
                     int target_size)
      : candidates_(std::move(candidates)),
        small_(std::move(small_transversals)),
        lower_(lower),
        upper_(upper),
        target_size_(target_size) {
    if (small_.size() > 63) {
      throw std::runtime_error("small-transversal word no longer fits uint64");
    }
    full_cover_ = small_.empty() ? 0 : (std::uint64_t{1} << small_.size()) - 1;
    cover_.resize(candidates_.size());
    for (int index = 0; index < static_cast<int>(candidates_.size()); ++index) {
      for (int row = 0; row < static_cast<int>(small_.size()); ++row) {
        if ((candidates_[index] & small_[row]) == 0) {
          cover_[index] |= std::uint64_t{1} << row;
        }
      }
    }
  }

  std::set<Family> Run() {
    if (target_size_ == 4) return RunFourByPairs();
    Family chosen;
    std::vector<bool> used(candidates_.size());
    std::array<int, kN> incidence{};
    Search(chosen, used, incidence, 0);
    return solutions_;
  }

 private:
  struct PairRow {
    int left;
    int right;
    std::uint64_t covered;
  };

  bool CanAdd(int index, const std::array<int, kN>& incidence) const {
    const auto candidate = candidates_[index];
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (((candidate >> vertex) & 1) && incidence[vertex] >= upper_[vertex]) {
        return false;
      }
    }
    return true;
  }

  bool LowerCanStillBeMet(const std::array<int, kN>& incidence,
                          int slots_left) const {
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (incidence[vertex] + slots_left < lower_[vertex]) return false;
    }
    return true;
  }

  bool PairRowsCompatible(const PairRow& left, const PairRow& right) const {
    for (int vertex = 0; vertex < kN; ++vertex) {
      const int incidence =
          ((candidates_[left.left] >> vertex) & 1) +
          ((candidates_[left.right] >> vertex) & 1) +
          ((candidates_[right.left] >> vertex) & 1) +
          ((candidates_[right.right] >> vertex) & 1);
      if (incidence < lower_[vertex] || incidence > upper_[vertex]) return false;
    }
    return true;
  }

  std::set<Family> RunFourByPairs() const {
    std::vector<PairRow> pairs;
    for (int left = 0; left < static_cast<int>(candidates_.size()); ++left) {
      std::array<int, kN> first{};
      for (int vertex = 0; vertex < kN; ++vertex) {
        first[vertex] = (candidates_[left] >> vertex) & 1;
      }
      for (int right = left + 1; right < static_cast<int>(candidates_.size()); ++right) {
        if (!CanAdd(right, first)) continue;
        pairs.push_back({left, right, cover_[left] | cover_[right]});
      }
    }

    const int rows = static_cast<int>(small_.size());
    std::vector<std::vector<int>> single(rows);
    std::vector<std::vector<int>> joint(rows * rows);
    for (int index = 0; index < static_cast<int>(pairs.size()); ++index) {
      std::vector<int> bits;
      for (int bit = 0; bit < rows; ++bit) {
        if ((pairs[index].covered >> bit) & 1) {
          single[bit].push_back(index);
          bits.push_back(bit);
        }
      }
      for (int a = 0; a < static_cast<int>(bits.size()); ++a) {
        for (int b = a + 1; b < static_cast<int>(bits.size()); ++b) {
          joint[bits[a] * rows + bits[b]].push_back(index);
        }
      }
    }

    std::set<Family> answer;
    for (const auto& first : pairs) {
      const std::uint64_t missing = full_cover_ & ~first.covered;
      std::vector<int> bits;
      for (int bit = 0; bit < rows; ++bit) {
        if ((missing >> bit) & 1) bits.push_back(bit);
      }
      const std::vector<int>* query = nullptr;
      if (bits.empty()) {
        // This branch is not expected for record 30, but retaining it keeps
        // the pair decomposition lossless.
      } else if (bits.size() == 1) {
        query = &single[bits[0]];
      } else {
        for (int a = 0; a < static_cast<int>(bits.size()); ++a) {
          for (int b = a + 1; b < static_cast<int>(bits.size()); ++b) {
            const auto& candidate = joint[bits[a] * rows + bits[b]];
            if (query == nullptr || candidate.size() < query->size()) query = &candidate;
          }
        }
      }

      auto consume = [&](const PairRow& second) {
        // Every sorted four-set has the unique split (first two | last two).
        if (first.right >= second.left ||
            (first.covered | second.covered) != full_cover_ ||
            !PairRowsCompatible(first, second)) return;
        Family family = {
            candidates_[first.left], candidates_[first.right],
            candidates_[second.left], candidates_[second.right],
        };
        std::sort(family.begin(), family.end());
        answer.insert(std::move(family));
      };
      if (query == nullptr) {
        for (const auto& second : pairs) consume(second);
      } else {
        for (const int index : *query) consume(pairs[index]);
      }
    }
    return answer;
  }

  void InsertSolution(const Family& chosen,
                      const std::array<int, kN>& incidence,
                      std::uint64_t covered) {
    if (static_cast<int>(chosen.size()) != target_size_ ||
        covered != full_cover_) {
      return;
    }
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (incidence[vertex] < lower_[vertex] ||
          incidence[vertex] > upper_[vertex]) return;
    }
    Family canonical = chosen;
    std::sort(canonical.begin(), canonical.end());
    solutions_.insert(std::move(canonical));
  }

  void FillAfterCovered(Family& chosen, std::vector<bool>& used,
                        std::array<int, kN>& incidence,
                        std::uint64_t covered, int start) {
    if (static_cast<int>(chosen.size()) == target_size_) {
      InsertSolution(chosen, incidence, covered);
      return;
    }
    const int slots_left = target_size_ - static_cast<int>(chosen.size());
    if (static_cast<int>(candidates_.size()) - start < slots_left ||
        !LowerCanStillBeMet(incidence, slots_left)) return;
    for (int index = start; index < static_cast<int>(candidates_.size()); ++index) {
      if (used[index] || !CanAdd(index, incidence)) continue;
      used[index] = true;
      chosen.push_back(candidates_[index]);
      for (int vertex = 0; vertex < kN; ++vertex) {
        incidence[vertex] += (candidates_[index] >> vertex) & 1;
      }
      FillAfterCovered(chosen, used, incidence, covered, index + 1);
      for (int vertex = 0; vertex < kN; ++vertex) {
        incidence[vertex] -= (candidates_[index] >> vertex) & 1;
      }
      chosen.pop_back();
      used[index] = false;
    }
  }

  void Search(Family& chosen, std::vector<bool>& used,
              std::array<int, kN>& incidence, std::uint64_t covered) {
    if (static_cast<int>(chosen.size()) == target_size_) {
      InsertSolution(chosen, incidence, covered);
      return;
    }
    const int slots_left = target_size_ - static_cast<int>(chosen.size());
    if (!LowerCanStillBeMet(incidence, slots_left)) return;
    if (covered == full_cover_) {
      FillAfterCovered(chosen, used, incidence, covered, 0);
      return;
    }

    int pivot = -1;
    int best_count = static_cast<int>(candidates_.size()) + 1;
    for (int row = 0; row < static_cast<int>(small_.size()); ++row) {
      if ((covered >> row) & 1) continue;
      int count = 0;
      for (int index = 0; index < static_cast<int>(candidates_.size()); ++index) {
        if (!used[index] && ((cover_[index] >> row) & 1) &&
            CanAdd(index, incidence)) ++count;
      }
      if (count < best_count) {
        best_count = count;
        pivot = row;
      }
    }
    if (pivot < 0 || best_count == 0) return;

    for (int index = 0; index < static_cast<int>(candidates_.size()); ++index) {
      if (used[index] || !((cover_[index] >> pivot) & 1) ||
          !CanAdd(index, incidence)) continue;
      used[index] = true;
      chosen.push_back(candidates_[index]);
      for (int vertex = 0; vertex < kN; ++vertex) {
        incidence[vertex] += (candidates_[index] >> vertex) & 1;
      }
      Search(chosen, used, incidence, covered | cover_[index]);
      for (int vertex = 0; vertex < kN; ++vertex) {
        incidence[vertex] -= (candidates_[index] >> vertex) & 1;
      }
      chosen.pop_back();
      used[index] = false;
    }
  }

  std::vector<std::uint16_t> candidates_;
  std::vector<std::uint16_t> small_;
  std::array<int, kN> lower_{};
  std::array<int, kN> upper_{};
  int target_size_ = 0;
  std::uint64_t full_cover_ = 0;
  std::vector<std::uint64_t> cover_;
  std::set<Family> solutions_;
};

OrbitResult AuditOrbit(int orbit_index,
                       std::array<std::uint16_t, kA4Count> fours,
                       const std::vector<std::uint16_t>& fives) {
  OrbitResult result;
  result.orbit_index = orbit_index;
  result.fours = fours;

  Family shadow;
  for (const auto five : fives) {
    if (std::ranges::any_of(fours, [five](std::uint16_t four) {
          return (five & four) == four;
        })) shadow.push_back(five);
  }
  result.forced_shadow_count = static_cast<int>(shadow.size());
  result.addition_count = kScheduledA5Count - result.forced_shadow_count;
  if (result.addition_count < 0 || result.addition_count > 4) {
    throw std::runtime_error("record-30 a3 addition count left 0..4");
  }

  std::array<int, kN> lower{};
  std::array<int, kN> upper{};
  for (int vertex = 0; vertex < kN; ++vertex) {
    int a4_incidence = 0;
    for (const auto four : fours) a4_incidence += (four >> vertex) & 1;
    int forced_incidence = 0;
    for (const auto five : shadow) forced_incidence += (five >> vertex) & 1;
    const int previous = kBaseG3[vertex] - a4_incidence;
    const int before = kBaseG4[vertex] + 6 * a4_incidence - forced_incidence;
    lower[vertex] = std::max(0, before - MacaulayRaise(previous, 3));
    upper[vertex] = before - 1;
    if (upper[vertex] < 0) {
      throw std::runtime_error("positive scout passed a negative vertex cap");
    }
  }

  Family candidates;
  for (const auto five : fives) {
    if (std::ranges::find(shadow, five) != shadow.end()) continue;
    bool possible = true;
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (((five >> vertex) & 1) && upper[vertex] == 0) possible = false;
    }
    if (possible) candidates.push_back(five);
  }
  result.candidate_count = static_cast<int>(candidates.size());
  result.raw_addition_count = Binomial(result.candidate_count,
                                       result.addition_count);

  Family base(kMissingEdges.begin(), kMissingEdges.end());
  base.insert(base.end(), fours.begin(), fours.end());
  base.insert(base.end(), shadow.begin(), shadow.end());
  std::vector<std::uint16_t> small;
  for (std::uint16_t mask = 1; mask < (1u << kN); ++mask) {
    if (std::popcount(mask) <= 4 && HitsAll(mask, base)) small.push_back(mask);
  }
  result.small_transversal_count = static_cast<int>(small.size());

  // Once all transversals of order at most four have been eliminated, every
  // order-five transversal is automatically minimal.  Precompute the
  // order-five transversals of the fixed part of the clutter so that the
  // scheduled facet-count gate can be applied before the much more expensive
  // literal all-subsets blocker audit.
  Family base_order_five_transversals;
  for (std::uint16_t mask = 1; mask < (1u << kN); ++mask) {
    if (std::popcount(mask) == 5 && HitsAll(mask, base)) {
      base_order_five_transversals.push_back(mask);
    }
  }

  if (result.addition_count == 0) {
    bool interval_ok = std::ranges::all_of(lower, [](int value) {
      return value == 0;
    });
    if (interval_ok && small.empty()) result.no_small_families.insert({});
  } else {
    AdditionEnumerator enumerator(candidates, small, lower, upper,
                                  result.addition_count);
    result.no_small_families = enumerator.Run();
  }

  for (const auto& additions : result.no_small_families) {
    int order_five_count = 0;
    for (const auto transversal : base_order_five_transversals) {
      if (std::ranges::all_of(additions, [transversal](std::uint16_t face) {
            return (transversal & face) != 0;
          })) {
        ++order_five_count;
      }
    }
    ++result.order_five_transversal_histogram[order_five_count];
    if (order_five_count != kScheduledFacetCount) continue;

    Family clutter = base;
    clutter.insert(clutter.end(), additions.begin(), additions.end());
    const auto audit = AuditSolution(clutter);
    ++result.scheduled_order_five_audit_histogram[audit];
    if (audit.blocker_size_histogram == std::map<int, int>{{5, kScheduledFacetCount}}) {
      result.scheduled_pure_targets.push_back(additions);
    }
  }
  return result;
}

void PrintFamily(const Family& family) {
  std::cout << '[';
  for (int index = 0; index < static_cast<int>(family.size()); ++index) {
    if (index) std::cout << ',';
    std::cout << family[index];
  }
  std::cout << ']';
}

void PrintResult(const OrbitResult& result) {
  std::cout << "{\"orbit_index\":" << result.orbit_index << ",\"A4_masks\":[";
  for (int index = 0; index < kA4Count; ++index) {
    if (index) std::cout << ',';
    std::cout << result.fours[index];
  }
  std::cout << "],\"forced_A5_shadow_count\":" << result.forced_shadow_count
            << ",\"additional_A5_count\":" << result.addition_count
            << ",\"candidate_count\":" << result.candidate_count
            << ",\"raw_unordered_addition_count\":" << result.raw_addition_count
            << ",\"base_small_transversal_count\":"
            << result.small_transversal_count
            << ",\"no_small_transversal_family_count\":"
            << result.no_small_families.size()
            << ",\"order_five_transversal_count_histogram\":{";
  bool first = true;
  for (const auto& [facets, count] : result.order_five_transversal_histogram) {
    if (!first) std::cout << ',';
    std::cout << '\"' << facets << "\":" << count;
    first = false;
  }
  std::cout << "},\"scheduled_order_five_blocker_audit_histogram\":[";
  first = true;
  for (const auto& [audit, count] : result.scheduled_order_five_audit_histogram) {
    if (!first) std::cout << ',';
    std::cout << "{\"blocker_count\":" << audit.blocker_count
              << ",\"blocker_size_histogram\":{";
    bool first_size = true;
    for (const auto& [size, number] : audit.blocker_size_histogram) {
      if (!first_size) std::cout << ',';
      std::cout << '\"' << size << "\":" << number;
      first_size = false;
    }
    std::cout << "},\"family_count\":" << count << '}';
    first = false;
  }
  std::cout << "],\"scheduled_pure_target_count\":"
            << result.scheduled_pure_targets.size()
            << ",\"scheduled_pure_target_additions\":[";
  first = true;
  for (const auto& family : result.scheduled_pure_targets) {
    if (!first) std::cout << ',';
    PrintFamily(family);
    first = false;
  }
  std::cout << "]}";
}

}  // namespace

int main() {
  try {
    int orbit_count = 0;
    if (!(std::cin >> orbit_count) || orbit_count <= 0) {
      throw std::runtime_error("missing positive-orbit input");
    }
    const auto fives = IndependentFives();
    std::vector<OrbitResult> results;
    for (int row = 0; row < orbit_count; ++row) {
      int orbit_index = -1;
      std::array<std::uint16_t, kA4Count> fours{};
      unsigned raw = 0;
      if (!(std::cin >> orbit_index)) throw std::runtime_error("short orbit input");
      for (auto& four : fours) {
        if (!(std::cin >> raw) || raw >= (1u << kN) || std::popcount(raw) != 4) {
          throw std::runtime_error("invalid A4 input mask");
        }
        four = static_cast<std::uint16_t>(raw);
      }
      auto result = AuditOrbit(orbit_index, fours, fives);
      std::cerr << "orbit " << orbit_index
                << ": add=" << result.addition_count
                << " no-small=" << result.no_small_families.size()
                << " targets=" << result.scheduled_pure_targets.size() << '\n';
      results.push_back(std::move(result));
    }
    std::cout << "{\"kind\":\"h9_14_record30_a3_addition_census\","
              << "\"compact_only\":true,\"orbit_count\":" << results.size()
              << ",\"orbits\":[";
    for (int index = 0; index < static_cast<int>(results.size()); ++index) {
      if (index) std::cout << ',';
      PrintResult(results[index]);
    }
    std::cout << "],\"scope\":\"Exact necessary compact-incidence census; "
                 "no finite-volume noncompact conclusion.\"}\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
