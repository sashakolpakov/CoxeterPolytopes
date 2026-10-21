// Standalone audit of the positive-missing-four branches for H9/14 record 18.
//
// This program deliberately uses no project headers or generated search data.
// It starts from the five missing edges, derives the local g/Macaulay bounds,
// and quotients unordered A4 families by viewing them as binary incidence
// matrices.  The four noncentral active columns are acted on by the order-eight
// automorphism group of the double star; the other eight columns are stored as
// multiplicities of their row patterns, which is exactly quotienting by S8.
// Row permutations are handled explicitly.  This avoids enumerating the more
// than 2.4 billion literal four-row families.

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

constexpr int kN = 14;
constexpr int kLeafCount = 4;
constexpr int kIsolateFirst = 4;
constexpr int kIsolateCount = 8;
constexpr std::array<int, 6> kActive = {0, 1, 2, 3, 12, 13};
constexpr std::array<std::uint16_t, 5> kMissingEdges = {
    (1u << 0) | (1u << 12),
    (1u << 1) | (1u << 12),
    (1u << 2) | (1u << 13),
    (1u << 3) | (1u << 13),
    (1u << 12) | (1u << 13),
};

using Mask = std::uint16_t;
using ActiveRows = std::vector<int>;
using Permutation = std::vector<int>;
using ColumnCounts = std::array<std::uint8_t, 16>;
using Key = std::vector<int>;

struct LocalBound {
  Mask support;
  int lower;
  int upper;
};

struct ActiveAction {
  std::array<int, kN> image{};
};

struct ActiveMatch {
  int action_index;
  Permutation old_row_at_new_row;
};

struct GateResult {
  bool strict;
  bool aggregate;
  int shadow_size;
  int additional_fives;
  int lower_sum;
  int upper_sum;
};

struct OrbitRecord {
  Key key;
  std::uint64_t literal_size;
  GateResult gate;
};

struct BranchTotals {
  int a = 0;
  std::uint64_t all = 0;
  std::uint64_t interval = 0;
  std::uint64_t strict = 0;
  std::uint64_t viable = 0;
  std::uint64_t quotient_all = 0;
  std::uint64_t quotient_interval = 0;
  std::uint64_t quotient_strict = 0;
  std::vector<OrbitRecord> viable_orbits;
};

std::vector<LocalBound> g_bounds;
std::vector<ActiveAction> g_actions;
std::vector<Mask> g_independent_fives;
std::array<int, kN> g_base_g3{};
std::array<int, kN> g_base_g4{};

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

std::uint64_t Factorial(int n) {
  std::uint64_t value = 1;
  for (int i = 2; i <= n; ++i) value *= static_cast<std::uint64_t>(i);
  return value;
}

int MacaulayRaise(int value, int degree) {
  if (value < 0 || degree < 1) {
    throw std::runtime_error("invalid Macaulay input");
  }
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

bool IsIndependent(Mask value) {
  for (const Mask edge : kMissingEdges) {
    if ((value & edge) == edge) return false;
  }
  return true;
}

std::vector<Mask> IndependentMasks(int cardinality) {
  std::vector<Mask> answer;
  for (Mask value = 0; value < (Mask{1} << kN); ++value) {
    if (std::popcount(value) == cardinality && IsIndependent(value)) {
      answer.push_back(value);
    }
  }
  return answer;
}

Mask MoveMask(Mask value, const ActiveAction& action) {
  Mask moved = 0;
  for (int vertex = 0; vertex < kN; ++vertex) {
    if ((value >> vertex) & 1u) moved |= Mask{1} << action.image[vertex];
  }
  return moved;
}

void BuildActiveGroup() {
  std::array<int, 6> images = kActive;
  const std::set<Mask> target(kMissingEdges.begin(), kMissingEdges.end());
  do {
    ActiveAction action;
    std::iota(action.image.begin(), action.image.end(), 0);
    for (int i = 0; i < 6; ++i) action.image[kActive[i]] = images[i];
    std::set<Mask> moved;
    for (const Mask edge : kMissingEdges) moved.insert(MoveMask(edge, action));
    if (moved == target) g_actions.push_back(action);
  } while (std::next_permutation(images.begin(), images.end()));
  if (g_actions.size() != 8) {
    throw std::runtime_error("active double-star automorphism order is not eight");
  }
}

std::pair<int, int> StrictInterval(int previous, int beta,
                                   int candidate_count, int degree) {
  const int lower = std::max(0, beta - MacaulayRaise(previous, degree));
  const int upper = std::min(candidate_count, beta - 1);
  if (lower > upper) throw std::runtime_error("empty strict local interval");
  return {lower, upper};
}

void BuildLocalBounds() {
  for (int vertex = 0; vertex < kN; ++vertex) {
    std::vector<int> link;
    for (int other = 0; other < kN; ++other) {
      const Mask pair = (Mask{1} << vertex) | (Mask{1} << other);
      if (other != vertex && IsIndependent(pair)) link.push_back(other);
    }
    int link_edges = 0;
    int link_triangles = 0;
    for (int i = 0; i < static_cast<int>(link.size()); ++i) {
      for (int j = i + 1; j < static_cast<int>(link.size()); ++j) {
        const Mask pair = (Mask{1} << link[i]) | (Mask{1} << link[j]);
        link_edges += IsIndependent(pair);
        for (int k = j + 1; k < static_cast<int>(link.size()); ++k) {
          const Mask triple = pair | (Mask{1} << link[k]);
          link_triangles += IsIndependent(triple);
        }
      }
    }
    const int vertices = static_cast<int>(link.size());
    const int g1 = vertices - 9;
    const int g2 = link_edges - 8 * vertices + 36;
    const int beta =
        link_triangles - 7 * link_edges + 28 * vertices - 84;
    if (g1 <= 0 || g2 < 1 || g2 > MacaulayRaise(g1, 1)) {
      throw std::runtime_error("vertex-link prerequisite failed");
    }
    const auto [lower, upper] =
        StrictInterval(g2, beta, link_triangles, 2);
    g_bounds.push_back({static_cast<Mask>(Mask{1} << vertex), lower, upper});
  }

  for (int left = 0; left < kN; ++left) {
    for (int right = left + 1; right < kN; ++right) {
      const Mask support = (Mask{1} << left) | (Mask{1} << right);
      if (!IsIndependent(support)) continue;
      std::vector<int> link;
      for (int other = 0; other < kN; ++other) {
        if (other == left || other == right) continue;
        if (IsIndependent(support | (Mask{1} << other))) link.push_back(other);
      }
      int link_edges = 0;
      for (int i = 0; i < static_cast<int>(link.size()); ++i) {
        for (int j = i + 1; j < static_cast<int>(link.size()); ++j) {
          link_edges += IsIndependent((Mask{1} << link[i]) |
                                      (Mask{1} << link[j]));
        }
      }
      const int vertices = static_cast<int>(link.size());
      const int g1 = vertices - 8;
      const int beta = link_edges - 7 * vertices + 28;
      if (g1 <= 0) throw std::runtime_error("edge-link prerequisite failed");
      const auto [lower, upper] = StrictInterval(g1, beta, link_edges, 1);
      g_bounds.push_back({support, lower, upper});
    }
  }

  for (int first = 0; first < kN; ++first) {
    for (int second = first + 1; second < kN; ++second) {
      for (int third = second + 1; third < kN; ++third) {
        const Mask support = (Mask{1} << first) | (Mask{1} << second) |
                             (Mask{1} << third);
        if (!IsIndependent(support)) continue;
        int candidates = 0;
        for (int other = 0; other < kN; ++other) {
          if ((support >> other) & 1u) continue;
          candidates += IsIndependent(support | (Mask{1} << other));
        }
        g_bounds.push_back({support, 0, candidates - 8});
      }
    }
  }

  if (g_bounds.size() != 410 ||
      std::any_of(g_bounds.begin(), g_bounds.end(),
                  [](const LocalBound& bound) { return bound.lower != 0; })) {
    throw std::runtime_error("derived local-bound table changed");
  }
}

bool PassesIntervals(const std::vector<Mask>& family) {
  for (const LocalBound& bound : g_bounds) {
    int incidence = 0;
    for (const Mask row : family) {
      incidence += (row & bound.support) == bound.support;
    }
    if (incidence < bound.lower || incidence > bound.upper) return false;
  }
  return true;
}

void BuildFiveSetsAndLinkBaselines() {
  g_independent_fives = IndependentMasks(5);
  std::array<std::vector<Mask>, 6> independent;
  for (int size = 2; size <= 5; ++size) {
    independent[size] = IndependentMasks(size);
  }
  for (int vertex = 0; vertex < kN; ++vertex) {
    std::array<int, 6> count{};
    for (int size = 2; size <= 5; ++size) {
      for (const Mask face : independent[size]) count[size] += (face >> vertex) & 1u;
    }
    const int f0 = count[2];
    const int f1 = count[3];
    const int f2 = count[4];
    const int f3 = count[5];
    g_base_g3[vertex] = f2 - 7 * f1 + 28 * f0 - 84;
    g_base_g4[vertex] =
        f3 - 6 * f2 + 21 * f1 - 56 * f0 + 126;
  }
  const std::array<int, kN> expected_g3 =
      {5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 1, 1};
  const std::array<int, kN> expected_g4 =
      {6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 1, 1};
  if (g_independent_fives.size() != 1212 || g_base_g3 != expected_g3 ||
      g_base_g4 != expected_g4) {
    throw std::runtime_error("independent-face baseline changed");
  }
}

GateResult EvaluateShadowGate(const std::vector<Mask>& family) {
  std::array<int, kN> row_incidence{};
  std::array<int, kN> forced_incidence{};
  for (const Mask row : family) {
    for (int vertex = 0; vertex < kN; ++vertex) {
      row_incidence[vertex] += (row >> vertex) & 1u;
    }
  }
  int shadow_size = 0;
  for (const Mask five : g_independent_fives) {
    bool forced = false;
    for (const Mask row : family) {
      if ((five & row) == row) {
        forced = true;
        break;
      }
    }
    if (!forced) continue;
    ++shadow_size;
    for (int vertex = 0; vertex < kN; ++vertex) {
      forced_incidence[vertex] += (five >> vertex) & 1u;
    }
  }
  const int total_fives = 8 + 6 * static_cast<int>(family.size());
  const int additional = total_fives - shadow_size;
  bool strict = additional >= 0;
  int lower_sum = 0;
  int upper_sum = 0;
  for (int vertex = 0; vertex < kN; ++vertex) {
    const int previous = g_base_g3[vertex] - row_incidence[vertex];
    const int before = g_base_g4[vertex] + 6 * row_incidence[vertex] -
                       forced_incidence[vertex];
    strict = strict && previous >= 1 && before >= 1;
    if (previous >= 0) {
      lower_sum += std::max(0, before - MacaulayRaise(previous, 3));
    }
    upper_sum += before - 1;
  }
  const bool aggregate = strict && lower_sum <= 5 * additional &&
                         5 * additional <= upper_sum;
  return {strict, aggregate, shadow_size, additional, lower_sum, upper_sum};
}

int MoveLeafMask(int mask, const ActiveAction& action) {
  int moved = 0;
  for (int leaf = 0; leaf < kLeafCount; ++leaf) {
    if ((mask >> leaf) & 1) moved |= 1 << action.image[leaf];
  }
  return moved;
}

void EnumeratePermutations(int a, std::vector<Permutation>* output) {
  Permutation value(a);
  std::iota(value.begin(), value.end(), 0);
  do {
    output->push_back(value);
  } while (std::next_permutation(value.begin(), value.end()));
}

bool IsCanonicalActiveRows(const ActiveRows& rows) {
  for (const ActiveAction& action : g_actions) {
    ActiveRows moved;
    moved.reserve(rows.size());
    for (const int row : rows) moved.push_back(MoveLeafMask(row, action));
    std::sort(moved.begin(), moved.end());
    if (moved < rows) return false;
  }
  return true;
}

std::vector<ActiveMatch> ActiveStabilizer(const ActiveRows& rows,
                                          const std::vector<Permutation>& perms) {
  std::vector<ActiveMatch> answer;
  for (int action_index = 0; action_index < static_cast<int>(g_actions.size());
       ++action_index) {
    for (const Permutation& perm : perms) {
      bool matches = true;
      for (int new_row = 0; new_row < static_cast<int>(rows.size()); ++new_row) {
        if (MoveLeafMask(rows[perm[new_row]], g_actions[action_index]) !=
            rows[new_row]) {
          matches = false;
          break;
        }
      }
      if (matches) answer.push_back({action_index, perm});
    }
  }
  if (answer.empty()) throw std::runtime_error("active stabilizer is empty");
  return answer;
}

int MoveColumnPattern(int pattern, const Permutation& old_row_at_new_row) {
  int moved = 0;
  for (int new_row = 0; new_row < static_cast<int>(old_row_at_new_row.size());
       ++new_row) {
    moved |= ((pattern >> old_row_at_new_row[new_row]) & 1) << new_row;
  }
  return moved;
}

ColumnCounts MoveColumns(const ColumnCounts& counts, int pattern_count,
                         const Permutation& perm) {
  ColumnCounts moved{};
  for (int pattern = 0; pattern < pattern_count; ++pattern) {
    moved[MoveColumnPattern(pattern, perm)] += counts[pattern];
  }
  return moved;
}

std::vector<int> ExpandedColumns(const ColumnCounts& counts, int pattern_count) {
  std::vector<int> answer;
  answer.reserve(kIsolateCount);
  for (int pattern = 0; pattern < pattern_count; ++pattern) {
    for (int repeat = 0; repeat < counts[pattern]; ++repeat) {
      answer.push_back(pattern);
    }
  }
  if (answer.size() != kIsolateCount) {
    throw std::runtime_error("isolate-column multiplicity does not sum to eight");
  }
  return answer;
}

bool IsCanonicalColumns(const ColumnCounts& counts, int pattern_count,
                        const std::vector<ActiveMatch>& stabilizer) {
  const std::vector<int> current = ExpandedColumns(counts, pattern_count);
  for (const ActiveMatch& match : stabilizer) {
    const ColumnCounts moved = MoveColumns(
        counts, pattern_count, match.old_row_at_new_row);
    if (ExpandedColumns(moved, pattern_count) < current) return false;
  }
  return true;
}

bool RowsAreDistinct(const ActiveRows& active, const ColumnCounts& counts,
                     int pattern_count) {
  const int a = static_cast<int>(active.size());
  for (int left = 0; left < a; ++left) {
    for (int right = left + 1; right < a; ++right) {
      if (active[left] != active[right]) continue;
      bool distinguished = false;
      for (int pattern = 0; pattern < pattern_count; ++pattern) {
        if (counts[pattern] == 0) continue;
        if (((pattern >> left) & 1) != ((pattern >> right) & 1)) {
          distinguished = true;
          break;
        }
      }
      if (!distinguished) return false;
    }
  }
  return true;
}

std::vector<Mask> ReconstructFamily(const ActiveRows& active,
                                    const ColumnCounts& counts,
                                    int pattern_count) {
  const int a = static_cast<int>(active.size());
  std::vector<Mask> family(a);
  for (int row = 0; row < a; ++row) family[row] = active[row];
  int isolate = kIsolateFirst;
  for (int pattern = 0; pattern < pattern_count; ++pattern) {
    for (int repeat = 0; repeat < counts[pattern]; ++repeat) {
      for (int row = 0; row < a; ++row) {
        if ((pattern >> row) & 1) family[row] |= Mask{1} << isolate;
      }
      ++isolate;
    }
  }
  if (isolate != kIsolateFirst + kIsolateCount) {
    throw std::runtime_error("wrong reconstructed isolate count");
  }
  for (const Mask row : family) {
    if (std::popcount(row) != 4) throw std::runtime_error("row is not an A4");
  }
  return family;
}

std::uint64_t LiteralOrbitSize(const ColumnCounts& counts, int pattern_count,
                               const std::vector<ActiveMatch>& stabilizer) {
  std::uint64_t internal = 1;
  for (int pattern = 0; pattern < pattern_count; ++pattern) {
    internal *= Factorial(counts[pattern]);
  }
  std::uint64_t stabilizer_size = 0;
  for (const ActiveMatch& match : stabilizer) {
    if (MoveColumns(counts, pattern_count, match.old_row_at_new_row) == counts) {
      stabilizer_size += internal;
    }
  }
  const std::uint64_t group_size = g_actions.size() * Factorial(kIsolateCount);
  if (stabilizer_size == 0 || group_size % stabilizer_size != 0) {
    throw std::runtime_error("nonintegral literal orbit size");
  }
  return group_size / stabilizer_size;
}

Key MakeKey(const ActiveRows& active, const ColumnCounts& counts,
            int pattern_count) {
  Key key(active.begin(), active.end());
  const std::vector<int> columns = ExpandedColumns(counts, pattern_count);
  key.insert(key.end(), columns.begin(), columns.end());
  return key;
}

using ColumnBuckets = std::map<int, std::vector<ColumnCounts>>;

int RowSumCode(const std::vector<int>& sums) {
  int code = 0;
  int scale = 1;
  for (const int sum : sums) {
    code += scale * sum;
    scale *= 5;
  }
  return code;
}

void GenerateColumnCountsRecursive(int a, int pattern, int remaining,
                                   std::vector<int>* sums,
                                   ColumnCounts* counts,
                                   ColumnBuckets* buckets) {
  const int pattern_count = 1 << a;
  if (pattern == pattern_count - 1) {
    bool valid = true;
    for (int row = 0; row < a; ++row) {
      if ((pattern >> row) & 1) (*sums)[row] += remaining;
      if ((*sums)[row] > 4) valid = false;
    }
    if (valid) {
      (*counts)[pattern] = static_cast<std::uint8_t>(remaining);
      (*buckets)[RowSumCode(*sums)].push_back(*counts);
      (*counts)[pattern] = 0;
    }
    for (int row = 0; row < a; ++row) {
      if ((pattern >> row) & 1) (*sums)[row] -= remaining;
    }
    return;
  }
  for (int multiplicity = 0; multiplicity <= remaining; ++multiplicity) {
    bool valid = true;
    for (int row = 0; row < a; ++row) {
      if ((pattern >> row) & 1) (*sums)[row] += multiplicity;
      if ((*sums)[row] > 4) valid = false;
    }
    if (valid) {
      (*counts)[pattern] = static_cast<std::uint8_t>(multiplicity);
      GenerateColumnCountsRecursive(a, pattern + 1, remaining - multiplicity,
                                    sums, counts, buckets);
      (*counts)[pattern] = 0;
    }
    for (int row = 0; row < a; ++row) {
      if ((pattern >> row) & 1) (*sums)[row] -= multiplicity;
    }
  }
}

ColumnBuckets GenerateColumnBuckets(int a) {
  ColumnBuckets buckets;
  std::vector<int> sums(a, 0);
  ColumnCounts counts{};
  GenerateColumnCountsRecursive(a, 0, kIsolateCount, &sums, &counts, &buckets);
  return buckets;
}

void GenerateActiveRowsRecursive(int a, int position, int minimum,
                                 ActiveRows* rows,
                                 std::vector<ActiveRows>* output) {
  if (position == a) {
    output->push_back(*rows);
    return;
  }
  for (int mask = minimum; mask < (1 << kLeafCount); ++mask) {
    (*rows)[position] = mask;
    GenerateActiveRowsRecursive(a, position + 1, mask, rows, output);
  }
}

BranchTotals RunBranch(int a) {
  const auto start = std::chrono::steady_clock::now();
  const int pattern_count = 1 << a;
  std::vector<Permutation> permutations;
  EnumeratePermutations(a, &permutations);
  const ColumnBuckets buckets = GenerateColumnBuckets(a);
  std::vector<ActiveRows> active_multisets;
  ActiveRows scratch(a);
  GenerateActiveRowsRecursive(a, 0, 0, &scratch, &active_multisets);

  BranchTotals totals;
  totals.a = a;
  std::set<Key> seen_keys;
  for (const ActiveRows& active : active_multisets) {
    if (!IsCanonicalActiveRows(active)) continue;
    const std::vector<ActiveMatch> stabilizer =
        ActiveStabilizer(active, permutations);
    std::vector<int> demands;
    demands.reserve(a);
    for (const int mask : active) {
      demands.push_back(4 - std::popcount(static_cast<unsigned>(mask)));
    }
    const auto bucket_it = buckets.find(RowSumCode(demands));
    if (bucket_it == buckets.end()) continue;
    for (const ColumnCounts& counts : bucket_it->second) {
      if (!RowsAreDistinct(active, counts, pattern_count)) continue;
      if (!IsCanonicalColumns(counts, pattern_count, stabilizer)) continue;
      const Key key = MakeKey(active, counts, pattern_count);
      if (!seen_keys.insert(key).second) {
        throw std::runtime_error("duplicate canonical matrix key");
      }
      const std::uint64_t orbit_size =
          LiteralOrbitSize(counts, pattern_count, stabilizer);
      totals.all += orbit_size;
      ++totals.quotient_all;
      const std::vector<Mask> family =
          ReconstructFamily(active, counts, pattern_count);
      if (!PassesIntervals(family)) continue;
      totals.interval += orbit_size;
      ++totals.quotient_interval;
      const GateResult gate = EvaluateShadowGate(family);
      if (!gate.strict) continue;
      totals.strict += orbit_size;
      ++totals.quotient_strict;
      if (!gate.aggregate) continue;
      totals.viable += orbit_size;
      totals.viable_orbits.push_back({key, orbit_size, gate});
    }
  }

  const std::uint64_t expected_all = Binomial(495, a);
  if (totals.all != expected_all) {
    throw std::runtime_error("orbit decomposition did not recover C(495,a)");
  }
  std::sort(totals.viable_orbits.begin(), totals.viable_orbits.end(),
            [](const OrbitRecord& left, const OrbitRecord& right) {
              return left.key < right.key;
            });
  const auto elapsed = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - start).count();
  std::cerr << "a=" << a << " completed in " << elapsed << " seconds\n";
  return totals;
}

void PrintKey(const Key& key) {
  std::cout << '[';
  for (int i = 0; i < static_cast<int>(key.size()); ++i) {
    if (i) std::cout << ',';
    std::cout << key[i];
  }
  std::cout << ']';
}

void PrintBranch(const BranchTotals& branch) {
  std::map<std::uint64_t, int> viable_size_histogram;
  for (const OrbitRecord& orbit : branch.viable_orbits) {
    ++viable_size_histogram[orbit.literal_size];
  }
  std::cout << "    {\"a\":" << branch.a
            << ",\"raw_literal_families\":" << branch.all
            << ",\"interval_literal_families\":" << branch.interval
            << ",\"strict_shadow_literal_families\":" << branch.strict
            << ",\"viable_literal_families\":" << branch.viable
            << ",\"raw_orbits\":" << branch.quotient_all
            << ",\"interval_orbits\":" << branch.quotient_interval
            << ",\"strict_shadow_orbits\":" << branch.quotient_strict
            << ",\"viable_orbits\":" << branch.viable_orbits.size()
            << ",\"viable_orbit_size_histogram\":{";
  bool first = true;
  for (const auto& [size, count] : viable_size_histogram) {
    if (!first) std::cout << ',';
    std::cout << '\"' << size << "\":" << count;
    first = false;
  }
  std::cout << "},\"viable_orbit_keys\":[";
  first = true;
  for (const OrbitRecord& orbit : branch.viable_orbits) {
    if (!first) std::cout << ',';
    std::cout << "{\"key\":";
    PrintKey(orbit.key);
    std::cout << ",\"literal_size\":" << orbit.literal_size
              << ",\"shadow_size\":" << orbit.gate.shadow_size
              << ",\"additional_fives\":" << orbit.gate.additional_fives
              << ",\"macaulay_sum_interval\":[" << orbit.gate.lower_sum
              << ',' << orbit.gate.upper_sum << "]}";
    first = false;
  }
  std::cout << "]}";
}

}  // namespace

int main() {
  try {
    BuildActiveGroup();
    BuildLocalBounds();
    BuildFiveSetsAndLinkBaselines();

    const std::vector<Mask> independent_fours = IndependentMasks(4);
    std::vector<Mask> interval_fours;
    for (const Mask face : independent_fours) {
      if (PassesIntervals({face})) interval_fours.push_back(face);
    }
    if (independent_fours.size() != 735 || interval_fours.size() != 495 ||
        std::any_of(interval_fours.begin(), interval_fours.end(),
                    [](Mask face) { return face & ((Mask{1} << 12) | (Mask{1} << 13)); })) {
      throw std::runtime_error("record-18 A4 ground-set reduction changed");
    }

    std::vector<BranchTotals> branches;
    for (int a = 1; a <= 4; ++a) branches.push_back(RunBranch(a));

    std::cout << "{\n"
              << "  \"kind\":\"independent_h9_14_record18_positive_orbit_audit\",\n"
              << "  \"active_missing_edge_graph\":\"double_star_S2_wr_S2\",\n"
              << "  \"active_automorphism_order\":" << g_actions.size() << ",\n"
              << "  \"isolate_count\":" << kIsolateCount << ",\n"
              << "  \"local_bound_count\":" << g_bounds.size() << ",\n"
              << "  \"independent_A4_count\":" << independent_fours.size() << ",\n"
              << "  \"interval_A4_count\":" << interval_fours.size() << ",\n"
              << "  \"branches\":[\n";
    for (int i = 0; i < static_cast<int>(branches.size()); ++i) {
      PrintBranch(branches[i]);
      std::cout << (i + 1 == static_cast<int>(branches.size()) ? "\n" : ",\n");
    }
    std::cout << "  ]\n}\n";
  } catch (const std::exception& error) {
    std::cerr << "audit failed: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
