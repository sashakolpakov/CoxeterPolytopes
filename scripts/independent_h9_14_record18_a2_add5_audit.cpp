// Independent exact audit of the high-addition record-18 branches:
// all fourteen a=2 representatives and the three a=1 representatives.
//
// This is deliberately not the production ridge-aware DFS and does not use
// a SAT/SMT backend.  It quotients every partial five-row incidence matrix by
// row permutations, the relevant A4 stabilizer on fixed columns, and all
// permutations of the residual isolate columns.  Mandatory small-cover,
// facet-deletion, and bad-ridge reducer sets then drive an exact recursion.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

using Mask = std::uint16_t;
using FacetWord = std::array<std::uint64_t, 4>;
using CandidateWord = std::array<std::uint64_t, 13>;
using StateKey = unsigned __int128;

struct StateKeyHash {
  std::size_t operator()(StateKey value) const {
    const auto low = static_cast<std::uint64_t>(value);
    const auto high = static_cast<std::uint64_t>(value >> 64);
    return std::hash<std::uint64_t>{}(
        low ^ (high + 0x9e3779b97f4a7c15ULL + (low << 6) + (low >> 2)));
  }
};

constexpr int kN = 14;
constexpr int kMaxAdditionCount = 6;
constexpr Mask kFull = (Mask{1} << kN) - 1;
constexpr std::array<Mask, 5> kEdges = {
    (Mask{1} << 0) | (Mask{1} << 12),
    (Mask{1} << 1) | (Mask{1} << 12),
    (Mask{1} << 2) | (Mask{1} << 13),
    (Mask{1} << 3) | (Mask{1} << 13),
    (Mask{1} << 12) | (Mask{1} << 13),
};
constexpr std::array<int, kN> kBaseG3 = {
    5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 1, 1};
constexpr std::array<int, kN> kBaseG4 = {
    6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 1, 1};

struct Candidate {
  Mask row = 0;
  std::uint64_t small_coverage = 0;
  FacetWord facet_keep{};
};

struct Transcript {
  std::uint64_t nodes = 0;
  std::uint64_t symmetry_duplicates = 0;
  std::uint64_t dead_cap = 0;
  std::uint64_t dead_small = 0;
  std::uint64_t dead_facet = 0;
  std::uint64_t dead_ridge = 0;
  std::uint64_t scheduled_ridge_terminals = 0;
  std::uint64_t literal_pure_targets = 0;
  std::array<std::uint64_t, kMaxAdditionCount + 1> distinct_states{};
};

std::uint64_t Binomial(int n, int k) {
  if (k < 0 || k > n) return 0;
  k = std::min(k, n - k);
  std::uint64_t answer = 1;
  for (int index = 1; index <= k; ++index) {
    answer = answer * (n - k + index) / index;
  }
  return answer;
}

int MacaulayRaise(int value, int degree) {
  int remaining = value;
  int next_top = -1;
  int raised = 0;
  for (int lower = degree; lower >= 1; --lower) {
    int top = lower - 1;
    while ((next_top < 0 || top + 1 < next_top) &&
           Binomial(top + 1, lower) <=
               static_cast<std::uint64_t>(remaining)) {
      ++top;
    }
    remaining -= static_cast<int>(Binomial(top, lower));
    raised += static_cast<int>(Binomial(top + 1, lower + 1));
    next_top = top;
  }
  if (remaining != 0) throw std::runtime_error("Macaulay expansion failed");
  return raised;
}

bool Independent(Mask value) {
  return std::all_of(kEdges.begin(), kEdges.end(), [value](Mask edge) {
    return (value & edge) != edge;
  });
}

bool HitsAll(Mask value, const std::vector<Mask>& clutter) {
  return std::all_of(clutter.begin(), clutter.end(), [value](Mask row) {
    return (value & row) != 0;
  });
}

std::vector<Mask> MinimalTransversals(const std::vector<Mask>& clutter) {
  std::vector<Mask> answer;
  for (Mask value = 1; value <= kFull; ++value) {
    if (!HitsAll(value, clutter)) continue;
    bool minimal = true;
    for (int vertex = 0; vertex < kN; ++vertex) {
      if ((value >> vertex) & 1u) {
        minimal = minimal && !HitsAll(value ^ (Mask{1} << vertex), clutter);
      }
    }
    if (minimal) answer.push_back(value);
  }
  return answer;
}

template <std::size_t Blocks>
int Count(const std::array<std::uint64_t, Blocks>& word) {
  int answer = 0;
  for (const auto block : word) answer += std::popcount(block);
  return answer;
}

template <std::size_t Blocks>
bool Test(const std::array<std::uint64_t, Blocks>& word, int index) {
  return (word[index / 64] >> (index % 64)) & 1u;
}

template <std::size_t Blocks>
void Set(std::array<std::uint64_t, Blocks>* word, int index) {
  (*word)[index / 64] |= std::uint64_t{1} << (index % 64);
}

template <std::size_t Blocks>
std::array<std::uint64_t, Blocks> Intersection(
    std::array<std::uint64_t, Blocks> left,
    const std::array<std::uint64_t, Blocks>& right) {
  for (int block = 0; block < static_cast<int>(Blocks); ++block) {
    left[block] &= right[block];
  }
  return left;
}

template <std::size_t Blocks>
void UnionInto(std::array<std::uint64_t, Blocks>* target,
               const std::array<std::uint64_t, Blocks>& source) {
  for (int block = 0; block < static_cast<int>(Blocks); ++block) {
    (*target)[block] |= source[block];
  }
}

template <std::size_t Blocks>
std::vector<int> Indices(const std::array<std::uint64_t, Blocks>& word) {
  std::vector<int> answer;
  for (int block = 0; block < static_cast<int>(Blocks); ++block) {
    std::uint64_t value = word[block];
    while (value) {
      const int offset = std::countr_zero(value);
      answer.push_back(64 * block + offset);
      value &= value - 1;
    }
  }
  return answer;
}

class Audit {
 public:
  explicit Audit(std::vector<Mask> fours) : fours_(std::move(fours)) {
    std::sort(fours_.begin(), fours_.end());
    if (fours_ == std::vector<Mask>{53, 54}) {
      fixed_count_ = 6;
      isolate_start_ = 6;
      actions_ = GeneratedActions({Swap(0, 1), Swap(4, 5)});
    } else if (fours_ == std::vector<Mask>{23, 27}) {
      fixed_count_ = 5;
      isolate_start_ = 5;
      actions_ = GeneratedActions({Swap(0, 1), Swap(2, 3)});
    } else if (fours_ == std::vector<Mask>{53}) {
      fixed_count_ = 6;
      isolate_start_ = 6;
      actions_ = GeneratedActions({ExchangeLeafPairs(), Swap(4, 5)});
    } else if (fours_ == std::vector<Mask>{23}) {
      fixed_count_ = 5;
      isolate_start_ = 5;
      actions_ = GeneratedActions({Swap(0, 1)});
    } else if (fours_ == std::vector<Mask>{15}) {
      fixed_count_ = 4;
      isolate_start_ = 4;
      actions_ = GeneratedActions(
          {Swap(0, 1), Swap(2, 3), ExchangeLeafPairs()});
    } else if (fours_.size() == 2) {
      // The remaining a=2 representatives have shadow size 16 and require
      // four added A5 rows.  Quotient only by the manifestly free isolate
      // columns here: using the identity subgroup of the A4 stabilizer is
      // slower but cannot identify inequivalent states.
      int largest_used_isolate = 3;
      for (const Mask four : fours_) {
        for (int vertex = 4; vertex < 12; ++vertex) {
          if ((four >> vertex) & 1u) largest_used_isolate = vertex;
        }
      }
      fixed_count_ = largest_used_isolate + 1;
      isolate_start_ = fixed_count_;
      actions_ = {Identity()};
    } else {
      throw std::runtime_error("unexpected high-addition representative");
    }
    Initialize();
  }

  Transcript Run() {
    FacetWord initial{};
    for (int facet = 0; facet < static_cast<int>(facets_.size()); ++facet) {
      Set(&initial, facet);
    }
    Search({}, 0, {}, initial);
    return transcript_;
  }

  int CandidateCount() const { return candidates_.size(); }
  int SmallCount() const { return small_.size(); }
  int FacetCount() const { return facets_.size(); }
  int RidgeCount() const { return ridges_.size(); }
  int FixedCount() const { return fixed_count_; }
  int IsolateCount() const { return 12 - isolate_start_; }
  int AdditionCount() const { return addition_count_; }
  int ScheduledFacets() const { return scheduled_facets_; }
  int ShadowCount() const { return shadow_.size(); }
  int ActionCount() const { return actions_.size(); }
  const std::vector<Mask>& Fours() const { return fours_; }
  const std::vector<std::vector<Mask>>& TargetFamilies() const {
    return target_families_;
  }
  const std::vector<std::uint64_t>& TargetOrbitSizes() const {
    return target_orbit_sizes_;
  }
  const std::array<int, kN>& Upper() const { return upper_; }
  const std::array<int, kN>& Lower() const { return lower_; }

 private:
  using Move = std::array<int, 12>;

  static Move Identity() {
    Move answer{};
    for (int vertex = 0; vertex < 12; ++vertex) answer[vertex] = vertex;
    return answer;
  }

  static Move Swap(int first, int second) {
    Move answer = Identity();
    std::swap(answer[first], answer[second]);
    return answer;
  }

  static Move ExchangeLeafPairs() {
    Move answer = Identity();
    answer[0] = 2;
    answer[2] = 0;
    answer[1] = 3;
    answer[3] = 1;
    return answer;
  }

  static Move Compose(const Move& left, const Move& right) {
    Move answer{};
    for (int vertex = 0; vertex < 12; ++vertex) {
      answer[vertex] = left[right[vertex]];
    }
    return answer;
  }

  static std::vector<Move> GeneratedActions(
      const std::vector<Move>& generators) {
    std::set<Move> found{Identity()};
    std::vector<Move> frontier{Identity()};
    while (!frontier.empty()) {
      const Move current = frontier.back();
      frontier.pop_back();
      for (const Move& generator : generators) {
        for (const Move composite : {Compose(current, generator),
                                     Compose(generator, current)}) {
          if (found.insert(composite).second) frontier.push_back(composite);
        }
      }
    }
    return {found.begin(), found.end()};
  }

  void Initialize() {
    if (fours_.empty() || fours_.size() > 2 ||
        std::any_of(fours_.begin(), fours_.end(), [](Mask value) {
          return std::popcount(value) != 4 || !Independent(value);
        })) {
      throw std::runtime_error("invalid A4 input");
    }
    std::vector<Mask> independent_fives;
    for (Mask value = 0; value <= kFull; ++value) {
      if (std::popcount(value) == 5 && Independent(value)) {
        independent_fives.push_back(value);
      }
    }
    for (const Mask value : independent_fives) {
      if (std::any_of(fours_.begin(), fours_.end(), [value](Mask four) {
            return (value & four) == four;
          })) {
        shadow_.push_back(value);
      }
    }
    const int a = fours_.size();
    addition_count_ = 8 + 6 * a - static_cast<int>(shadow_.size());
    scheduled_facets_ = 110 - 2 * a;
    const bool valid_schedule =
        (a == 1 && shadow_.size() == 8 && addition_count_ == 6) ||
        (a == 2 &&
         ((shadow_.size() == 15 && addition_count_ == 5) ||
          (shadow_.size() == 16 && addition_count_ == 4)));
    if (!valid_schedule) {
      throw std::runtime_error("forced A5 shadow/addition schedule changed");
    }
    for (int vertex = 0; vertex < kN; ++vertex) {
      int a4_incidence = 0;
      for (const Mask four : fours_) a4_incidence += (four >> vertex) & 1u;
      int forced = 0;
      for (const Mask value : shadow_) forced += (value >> vertex) & 1u;
      const int previous = kBaseG3[vertex] - a4_incidence;
      const int before = kBaseG4[vertex] + 6 * a4_incidence - forced;
      lower_[vertex] = std::max(0, before - MacaulayRaise(previous, 3));
      upper_[vertex] = before - 1;
    }
    if (std::any_of(lower_.begin(), lower_.end(), [](int value) {
          return value != 0;
        })) {
      throw std::runtime_error("high-addition lower interval unexpectedly nonzero");
    }

    base_.assign(kEdges.begin(), kEdges.end());
    base_.insert(base_.end(), fours_.begin(), fours_.end());
    base_.insert(base_.end(), shadow_.begin(), shadow_.end());
    for (Mask value = 1; value <= kFull; ++value) {
      if (std::popcount(value) <= 4 && HitsAll(value, base_)) {
        small_.push_back(value);
      }
      if (std::popcount(value) == 5 && HitsAll(value, base_)) {
        facets_.push_back(value);
      }
    }
    if (small_.size() > 64 || facets_.size() > 256) {
      throw std::runtime_error("fixed audit word capacity exceeded");
    }

    std::set<Mask> shadow_set(shadow_.begin(), shadow_.end());
    for (const Mask row : independent_fives) {
      if (shadow_set.contains(row)) continue;
      bool admitted = true;
      for (int vertex = 0; vertex < kN; ++vertex) {
        if (((row >> vertex) & 1u) && upper_[vertex] < 1) admitted = false;
      }
      if (!admitted) continue;
      Candidate candidate;
      candidate.row = row;
      for (int constraint = 0; constraint < static_cast<int>(small_.size());
           ++constraint) {
        if ((row & small_[constraint]) == 0) {
          candidate.small_coverage |= std::uint64_t{1} << constraint;
        }
      }
      for (int facet = 0; facet < static_cast<int>(facets_.size()); ++facet) {
        if (row & facets_[facet]) Set(&candidate.facet_keep, facet);
      }
      candidates_.push_back(candidate);
    }
    if (candidates_.size() > 13 * 64) {
      throw std::runtime_error("candidate audit word capacity exceeded");
    }
    for (int constraint = 0; constraint < static_cast<int>(small_.size());
         ++constraint) {
      CandidateWord coverers{};
      for (int candidate = 0;
           candidate < static_cast<int>(candidates_.size()); ++candidate) {
        if ((candidates_[candidate].small_coverage >> constraint) & 1u) {
          Set(&coverers, candidate);
        }
      }
      small_coverers_.push_back(coverers);
    }
    facet_killers_.resize(facets_.size());
    for (int facet = 0; facet < static_cast<int>(facets_.size()); ++facet) {
      for (int candidate = 0;
           candidate < static_cast<int>(candidates_.size()); ++candidate) {
        if (!Test(candidates_[candidate].facet_keep, facet)) {
          Set(&facet_killers_[facet], candidate);
        }
      }
    }
    for (Mask six = 0; six <= kFull; ++six) {
      if (std::popcount(six) != 6) continue;
      std::vector<int> members;
      for (int facet = 0; facet < static_cast<int>(facets_.size()); ++facet) {
        if ((facets_[facet] & six) == facets_[facet]) {
          members.push_back(facet);
        }
      }
      if (!members.empty()) ridges_.push_back(std::move(members));
    }
    seen_.resize(addition_count_ + 1);
    full_small_ = small_.size() == 64
                      ? ~std::uint64_t{0}
                      : (std::uint64_t{1} << small_.size()) - 1;
  }

  StateKey CanonicalKey(const std::vector<int>& selected) const {
    const int depth = selected.size();
    if (depth == 0) return 0;
    StateKey best = 0;
    bool first = true;
    for (const Move& move : actions_) {
      std::vector<std::pair<int, int>> fixed_rows;
      for (int row = 0; row < depth; ++row) {
        fixed_rows.push_back(
            {MovedFixedMask(candidates_[selected[row]].row, move), row});
      }
      std::sort(fixed_rows.begin(), fixed_rows.end());
      std::vector<std::vector<int>> groups;
      for (const auto& [fixed, row] : fixed_rows) {
        if (groups.empty() ||
            MovedFixedMask(candidates_[selected[groups.back().front()]].row,
                           move) != fixed) {
          groups.push_back({});
        }
        groups.back().push_back(row);
      }
      std::vector<std::vector<int>> orders(1);
      for (auto group : groups) {
        std::sort(group.begin(), group.end());
        std::vector<std::vector<int>> next;
        do {
          for (const auto& prefix : orders) {
            auto combined = prefix;
            combined.insert(combined.end(), group.begin(), group.end());
            next.push_back(std::move(combined));
          }
        } while (std::next_permutation(group.begin(), group.end()));
        orders = std::move(next);
      }
      for (const auto& order : orders) {
        StateKey key = 0;
        for (int position = 0; position < depth; ++position) {
          const Mask row = candidates_[selected[order[position]]].row;
          key = (key << fixed_count_) | MovedFixedMask(row, move);
        }
        std::vector<int> columns;
        for (int isolate = isolate_start_; isolate < 12; ++isolate) {
          int column = 0;
          for (int position = 0; position < depth; ++position) {
            const Mask row = candidates_[selected[order[position]]].row;
            column |= int((row >> isolate) & 1u) << position;
          }
          columns.push_back(column);
        }
        std::sort(columns.begin(), columns.end());
        for (const int column : columns) key = (key << depth) | column;
        if (first || key < best) {
          best = key;
          first = false;
        }
      }
    }
    return best;
  }

  bool Contains(const std::vector<int>& selected, int candidate) const {
    return std::find(selected.begin(), selected.end(), candidate) !=
           selected.end();
  }

  bool CapFeasible(const std::array<std::uint8_t, kN>& degrees,
                   Mask row) const {
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (((row >> vertex) & 1u) && degrees[vertex] >= upper_[vertex]) {
        return false;
      }
    }
    return true;
  }

  std::array<std::uint8_t, kN> Add(
      std::array<std::uint8_t, kN> degrees, Mask row) const {
    for (int vertex = 0; vertex < kN; ++vertex) {
      degrees[vertex] += (row >> vertex) & 1u;
    }
    return degrees;
  }

  CandidateWord Feasible(
      const std::vector<int>& selected,
      const std::array<std::uint8_t, kN>& degrees) const {
    CandidateWord answer{};
    for (int candidate = 0;
         candidate < static_cast<int>(candidates_.size()); ++candidate) {
      if (!Contains(selected, candidate) &&
          CapFeasible(degrees, candidates_[candidate].row)) {
        Set(&answer, candidate);
      }
    }
    return answer;
  }

  int RidgeDegree(const FacetWord& survivors,
                  const std::vector<int>& ridge) const {
    int degree = 0;
    for (const int facet : ridge) degree += Test(survivors, facet);
    return degree;
  }

  CandidateWord Reducers(const FacetWord& survivors,
                         const std::vector<int>& facets) const {
    CandidateWord answer{};
    for (const int facet : facets) {
      if (Test(survivors, facet)) UnionInto(&answer, facet_killers_[facet]);
    }
    return answer;
  }

  void KeepSmaller(CandidateWord* best, bool* has_best,
                   const CandidateWord& constraint) const {
    const int count = Count(constraint);
    if (count == 0) return;
    if (!*has_best || count < Count(*best)) {
      *best = constraint;
      *has_best = true;
    }
  }

  bool DeletionUpperBound(const CandidateWord& feasible,
                          const FacetWord& survivors, int remaining) const {
    const int excess = Count(survivors) - scheduled_facets_;
    if (excess <= 0) return true;
    std::vector<int> kills;
    for (const int candidate : Indices(feasible)) {
      const int kept = Count(Intersection(survivors,
                                          candidates_[candidate].facet_keep));
      kills.push_back(Count(survivors) - kept);
    }
    std::sort(kills.begin(), kills.end(), std::greater<>());
    int optimistic = 0;
    for (int index = 0;
         index < remaining && index < static_cast<int>(kills.size()); ++index) {
      optimistic += kills[index];
    }
    return optimistic >= excess;
  }

  static std::uint64_t Factorial(int value) {
    std::uint64_t answer = 1;
    for (int factor = 2; factor <= value; ++factor) answer *= factor;
    return answer;
  }

  int MovedFixedMask(Mask row, const Move& move) const {
    int answer = 0;
    for (int source = 0; source < fixed_count_; ++source) {
      if ((row >> source) & 1u) answer |= 1 << move[source];
    }
    return answer;
  }

  std::vector<int> IsolateColumns(const std::vector<int>& selected,
                                  const std::vector<int>& order) const {
    std::vector<int> columns;
    for (int isolate = isolate_start_; isolate < 12; ++isolate) {
      int column = 0;
      for (int position = 0; position < static_cast<int>(order.size());
           ++position) {
        const Mask row = candidates_[selected[order[position]]].row;
        column |= int((row >> isolate) & 1u) << position;
      }
      columns.push_back(column);
    }
    std::sort(columns.begin(), columns.end());
    return columns;
  }

  std::uint64_t OrbitSize(const std::vector<int>& selected) const {
    const int depth = selected.size();
    std::vector<int> identity(depth);
    for (int index = 0; index < depth; ++index) identity[index] = index;
    const auto target_columns = IsolateColumns(selected, identity);
    std::map<int, int> column_multiplicity;
    for (const int column : target_columns) ++column_multiplicity[column];
    std::uint64_t isolate_matching_permutations = 1;
    for (const auto& [column, multiplicity] : column_multiplicity) {
      static_cast<void>(column);
      isolate_matching_permutations *= Factorial(multiplicity);
    }

    std::uint64_t stabilizer = 0;
    std::vector<int> order = identity;
    do {
      if (IsolateColumns(selected, order) != target_columns) continue;
      for (const Move& move : actions_) {
        bool fixed_match = true;
        for (int position = 0; position < depth; ++position) {
          const Mask source = candidates_[selected[order[position]]].row;
          const Mask target = candidates_[selected[position]].row;
          const int target_fixed = target & ((1 << fixed_count_) - 1);
          if (MovedFixedMask(source, move) != target_fixed) {
            fixed_match = false;
            break;
          }
        }
        if (fixed_match) stabilizer += isolate_matching_permutations;
      }
    } while (std::next_permutation(order.begin(), order.end()));
    if (stabilizer == 0) throw std::runtime_error("zero family stabilizer");
    const std::uint64_t group_order =
        actions_.size() * Factorial(12 - isolate_start_);
    if (group_order % stabilizer != 0) {
      throw std::runtime_error("nonintegral target orbit size");
    }
    return group_order / stabilizer;
  }

  bool LiteralAudit(const std::vector<int>& selected) {
    std::vector<Mask> clutter = base_;
    for (const int index : selected) clutter.push_back(candidates_[index].row);
    const auto blocker = MinimalTransversals(clutter);
    std::map<int, int> histogram;
    for (const Mask value : blocker) ++histogram[std::popcount(value)];
    if (histogram == std::map<int, int>{{5, scheduled_facets_}}) {
      ++transcript_.literal_pure_targets;
      return true;
    }
    return false;
  }

  void Search(std::vector<int> selected, std::uint64_t covered,
              std::array<std::uint8_t, kN> degrees,
              const FacetWord& survivors) {
    ++transcript_.nodes;
    const int depth = selected.size();
    const auto key = CanonicalKey(selected);
    if (!seen_[depth].insert(key).second) {
      ++transcript_.symmetry_duplicates;
      return;
    }
    ++transcript_.distinct_states[depth];
    const int remaining = addition_count_ - depth;
    const int facet_count = Count(survivors);
    if (facet_count < scheduled_facets_) {
      ++transcript_.dead_facet;
      return;
    }
    if (depth == addition_count_) {
      if (covered != full_small_ || facet_count != scheduled_facets_) return;
      for (const auto& ridge : ridges_) {
        const int degree = RidgeDegree(survivors, ridge);
        if (degree != 0 && degree != 2) return;
      }
      ++transcript_.scheduled_ridge_terminals;
      if (LiteralAudit(selected)) {
        std::vector<Mask> family;
        for (const int index : selected) family.push_back(candidates_[index].row);
        std::sort(family.begin(), family.end());
        target_orbit_sizes_.push_back(OrbitSize(selected));
        target_families_.push_back(std::move(family));
      }
      return;
    }

    const CandidateWord feasible = Feasible(selected, degrees);
    if (Count(feasible) < remaining) {
      ++transcript_.dead_cap;
      return;
    }
    if (!DeletionUpperBound(feasible, survivors, remaining)) {
      ++transcript_.dead_facet;
      return;
    }

    CandidateWord best{};
    bool has_best = false;
    std::uint64_t possible_coverage = covered;
    for (int constraint = 0; constraint < static_cast<int>(small_.size());
         ++constraint) {
      if ((covered >> constraint) & 1u) continue;
      const auto choices = Intersection(feasible, small_coverers_[constraint]);
      if (Count(choices) == 0) {
        ++transcript_.dead_small;
        return;
      }
      KeepSmaller(&best, &has_best, choices);
    }
    for (const int candidate : Indices(feasible)) {
      possible_coverage |= candidates_[candidate].small_coverage;
    }
    if (possible_coverage != full_small_) {
      ++transcript_.dead_small;
      return;
    }

    for (const auto& ridge : ridges_) {
      const int degree = RidgeDegree(survivors, ridge);
      if (degree == 0 || degree == 2) continue;
      const auto choices = Intersection(feasible, Reducers(survivors, ridge));
      if (Count(choices) == 0) {
        ++transcript_.dead_ridge;
        return;
      }
      KeepSmaller(&best, &has_best, choices);
    }
    if (facet_count > scheduled_facets_) {
      std::vector<int> live_facets;
      for (int facet = 0; facet < static_cast<int>(facets_.size()); ++facet) {
        if (Test(survivors, facet)) live_facets.push_back(facet);
      }
      const auto choices =
          Intersection(feasible, Reducers(survivors, live_facets));
      if (Count(choices) == 0) {
        ++transcript_.dead_facet;
        return;
      }
      KeepSmaller(&best, &has_best, choices);
    }
    if (!has_best) best = feasible;

    for (const int candidate : Indices(best)) {
      auto next_selected = selected;
      next_selected.push_back(candidate);
      Search(std::move(next_selected),
             covered | candidates_[candidate].small_coverage,
             Add(degrees, candidates_[candidate].row),
             Intersection(survivors, candidates_[candidate].facet_keep));
    }
  }

  std::vector<Mask> fours_;
  int fixed_count_ = 0;
  int isolate_start_ = 0;
  int addition_count_ = 0;
  int scheduled_facets_ = 0;
  std::vector<Move> actions_;
  std::array<int, kN> lower_{};
  std::array<int, kN> upper_{};
  std::vector<Mask> shadow_;
  std::vector<Mask> base_;
  std::vector<Mask> small_;
  std::vector<Mask> facets_;
  std::vector<Candidate> candidates_;
  std::vector<CandidateWord> small_coverers_;
  std::vector<CandidateWord> facet_killers_;
  std::vector<std::vector<int>> ridges_;
  std::uint64_t full_small_ = 0;
  std::vector<std::unordered_set<StateKey, StateKeyHash>> seen_;
  std::vector<std::vector<Mask>> target_families_;
  std::vector<std::uint64_t> target_orbit_sizes_;
  Transcript transcript_;
};

void PrintIntArray(const std::array<int, kN>& values) {
  std::cout << '[';
  for (int index = 0; index < kN; ++index) {
    if (index) std::cout << ',';
    std::cout << values[index];
  }
  std::cout << ']';
}

void PrintMaskVertices(Mask value) {
  std::cout << '[';
  bool first = true;
  for (int vertex = 0; vertex < kN; ++vertex) {
    if (!((value >> vertex) & 1u)) continue;
    if (!first) std::cout << ',';
    std::cout << vertex;
    first = false;
  }
  std::cout << ']';
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 2 && argc != 3) {
      throw std::runtime_error("usage: program A4mask [A4mask]");
    }
    std::vector<Mask> fours;
    for (int index = 1; index < argc; ++index) {
      fours.push_back(static_cast<Mask>(std::stoul(argv[index], nullptr, 0)));
    }
    Audit audit(fours);
    const Transcript result = audit.Run();
    fours = audit.Fours();
    std::uint64_t symmetry_order = audit.ActionCount();
    for (int value = 2; value <= audit.IsolateCount(); ++value) {
      symmetry_order *= value;
    }
    std::cout << "{\n  \"kind\":\"independent_h9_14_record18_high_addition_audit\",\n"
              << "  \"a\":" << fours.size()
              << ",\n  \"representative_A4_masks\":[";
    for (int index = 0; index < static_cast<int>(fours.size()); ++index) {
      if (index) std::cout << ',';
      std::cout << fours[index];
    }
    std::cout << "],\n  \"forced_A5_shadow_count\":" << audit.ShadowCount()
              << ",\n  \"additional_A5_count\":" << audit.AdditionCount()
              << ",\n  \"scheduled_facet_count\":" << audit.ScheduledFacets()
              << ",\n  \"candidate_count\":" << audit.CandidateCount()
              << ",\n  \"base_small_transversal_count\":"
              << audit.SmallCount()
              << ",\n  \"base_order_five_transversal_count\":"
              << audit.FacetCount()
              << ",\n  \"nonempty_ridge_constraint_count\":"
              << audit.RidgeCount()
              << ",\n  \"incidence_symmetry_order\":" << symmetry_order
              << ",\n  \"addition_lower\":";
    PrintIntArray(audit.Lower());
    std::cout << ",\n  \"addition_upper\":";
    PrintIntArray(audit.Upper());
    std::cout << ",\n  \"search_node_count\":" << result.nodes
              << ",\n  \"distinct_state_counts\":[";
    for (int depth = 0; depth <= audit.AdditionCount(); ++depth) {
      if (depth) std::cout << ',';
      std::cout << result.distinct_states[depth];
    }
    std::cout << "],\n  \"symmetry_duplicate_count\":"
              << result.symmetry_duplicates
              << ",\n  \"scheduled_ridge_terminal_count\":"
              << result.scheduled_ridge_terminals
              << ",\n  \"literal_pure_target_count\":"
              << result.literal_pure_targets
              << ",\n  \"canonical_target_count\":"
              << audit.TargetFamilies().size()
              << ",\n  \"canonical_target_orbit_sizes\":[";
    std::uint64_t orbit_weight_sum = 0;
    for (int index = 0;
         index < static_cast<int>(audit.TargetOrbitSizes().size()); ++index) {
      if (index) std::cout << ',';
      std::cout << audit.TargetOrbitSizes()[index];
      orbit_weight_sum += audit.TargetOrbitSizes()[index];
    }
    std::cout << "],\n  \"literal_target_orbit_weight_sum\":"
              << orbit_weight_sum
              << ",\n  \"canonical_target_additions\":[";
    for (int family = 0;
         family < static_cast<int>(audit.TargetFamilies().size()); ++family) {
      if (family) std::cout << ',';
      std::cout << '[';
      for (int row = 0;
           row < static_cast<int>(audit.TargetFamilies()[family].size()); ++row) {
        if (row) std::cout << ',';
        PrintMaskVertices(audit.TargetFamilies()[family][row]);
      }
      std::cout << ']';
    }
    std::cout << ']'
              << ",\n  \"compact_incidence_unsat\":"
              << (audit.TargetFamilies().empty() ? "true" : "false")
              << ",\n  \"compact_only\":true,\n"
              << "  \"finite_volume_noncompact_checked\":false,\n"
              << "  \"scope\":\"Exact no-SAT symmetry-quotiented incidence/ridge audit; no ideal-vertex conclusion.\"\n"
              << "}\n";
  } catch (const std::exception& error) {
    std::cerr << "audit failed: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
