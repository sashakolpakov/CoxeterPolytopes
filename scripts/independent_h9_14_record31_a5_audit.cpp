// Independent direct quotient audit of the H9/14 record-31 a=5 branch.
//
// This verifier deliberately does not use the production extension search.
// It writes each unordered five-row A4 family as five active P5 row masks
// together with a multiset of nine binary isolate columns.  Nondecreasing
// active rows and nondecreasing columns quotient row order and S9 directly;
// repeated active rows may create harmless duplicate paths.  Thus every
// literal family occurs at least once, while no conclusion depends on orbit
// weights.  Every terminal matrix is checked from scratch against all local
// intervals and the forced-A5 strict/Macaulay gate.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>

namespace {

constexpr int kN = 14;
constexpr int kRows = 5;
constexpr int kScheduledFives = 39;
constexpr std::array<int, 5> kActive = {0, 1, 11, 12, 13};
constexpr std::array<int, 9> kIsolates = {2, 3, 4, 5, 6, 7, 8, 9, 10};
constexpr std::array<std::uint16_t, 4> kEdges = {
    (1u << 0) | (1u << 11),
    (1u << 1) | (1u << 12),
    (1u << 0) | (1u << 13),
    (1u << 1) | (1u << 13),
};

struct Constraint {
  std::uint16_t mask;
  int upper;
};

struct GateKey {
  bool strict;
  bool aggregate;
  int shadow;
  int additional;
  auto operator<=>(const GateKey&) const = default;
};

using Family = std::array<std::uint16_t, kRows>;
using Word = std::vector<std::uint64_t>;

std::vector<Constraint> constraints;
std::vector<std::uint16_t> fours;
std::vector<std::uint16_t> fives;
std::vector<Word> shadows;
std::array<Word, kN> vertex_five_words;
std::array<int, kN> base_g3{};
std::array<int, kN> base_g4{};
std::array<int, 1 << kN> four_index{};
int five_word_count = 0;

std::uint64_t active_multisets = 0;
std::uint64_t active_multisets_passing = 0;
std::uint64_t terminal_column_paths = 0;
std::uint64_t distinct_row_paths = 0;
std::uint64_t ground_row_paths = 0;
std::uint64_t interval_paths = 0;
std::uint64_t shadow_feasible_paths = 0;
std::uint64_t viable_paths = 0;
std::map<int, std::uint64_t> shadow_histogram;
std::map<GateKey, std::uint64_t> gate_histogram;

long long Binomial(int n, int k) {
  if (k < 0 || k > n) return 0;
  k = std::min(k, n - k);
  long long value = 1;
  for (int i = 1; i <= k; ++i) value = value * (n - k + i) / i;
  return value;
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
    remaining -= Binomial(top, lower);
    answer += Binomial(top + 1, lower + 1);
    next_top = top;
  }
  if (remaining != 0) throw std::runtime_error("Macaulay expansion failed");
  return answer;
}

bool MissingEdge(int a, int b) {
  const auto edge = static_cast<std::uint16_t>(
      (std::uint16_t{1} << a) | (std::uint16_t{1} << b));
  return std::find(kEdges.begin(), kEdges.end(), edge) != kEdges.end();
}

bool Independent(std::uint16_t mask) {
  for (const auto edge : kEdges) if ((mask & edge) == edge) return false;
  return true;
}

std::vector<std::uint16_t> Masks(int size) {
  std::vector<std::uint16_t> result;
  for (std::uint16_t mask = 0; mask < (1u << kN); ++mask) {
    if (std::popcount(mask) == size && Independent(mask)) result.push_back(mask);
  }
  return result;
}

std::pair<int, int> Interval(int previous, int beta,
                             int candidates, int degree) {
  const int lower = std::max(0, beta - MacaulayRaise(previous, degree));
  const int upper = std::min(candidates, beta - 1);
  return {lower, upper};
}

void BuildConstraints() {
  for (int v = 0; v < kN; ++v) {
    std::vector<int> link;
    for (int x = 0; x < kN; ++x) {
      if (x != v && !MissingEdge(v, x)) link.push_back(x);
    }
    int inherited = 0;
    for (int i = 0; i < static_cast<int>(link.size()); ++i)
      for (int j = i + 1; j < static_cast<int>(link.size()); ++j)
        inherited += MissingEdge(link[i], link[j]);
    const int n = link.size();
    const int f1 = Binomial(n, 2) - inherited;
    int f2 = 0;
    for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j)
      for (int k = j + 1; k < n; ++k) {
        const auto mask = static_cast<std::uint16_t>(
            (std::uint16_t{1} << link[i]) |
            (std::uint16_t{1} << link[j]) |
            (std::uint16_t{1} << link[k]));
        f2 += Independent(mask);
      }
    const int g1 = n - 9;
    const int g2 = f1 - 8 * n + 36;
    const int beta = f2 - 7 * f1 + 28 * n - 84;
    const auto [lower, upper] = Interval(g2, beta, f2, 2);
    if (g1 <= 0 || lower != 0) throw std::runtime_error("vertex interval changed");
    constraints.push_back({static_cast<std::uint16_t>(1u << v), upper});
  }
  for (int a = 0; a < kN; ++a) for (int b = a + 1; b < kN; ++b) {
    if (MissingEdge(a, b)) continue;
    std::vector<int> link;
    for (int x = 0; x < kN; ++x) {
      if (x != a && x != b && !MissingEdge(x, a) && !MissingEdge(x, b))
        link.push_back(x);
    }
    int f1 = 0;
    for (int i = 0; i < static_cast<int>(link.size()); ++i)
      for (int j = i + 1; j < static_cast<int>(link.size()); ++j)
        f1 += !MissingEdge(link[i], link[j]);
    const int n = link.size();
    const int g1 = n - 8;
    const int beta = f1 - 7 * n + 28;
    const auto [lower, upper] = Interval(g1, beta, f1, 1);
    if (lower != 0) throw std::runtime_error("edge interval changed");
    constraints.push_back({static_cast<std::uint16_t>(
        (std::uint16_t{1} << a) | (std::uint16_t{1} << b)), upper});
  }
  for (int a = 0; a < kN; ++a) for (int b = a + 1; b < kN; ++b)
    for (int c = b + 1; c < kN; ++c) {
      const auto triple = static_cast<std::uint16_t>(
          (std::uint16_t{1} << a) | (std::uint16_t{1} << b) |
          (std::uint16_t{1} << c));
      if (!Independent(triple)) continue;
      int candidates = 0;
      for (int x = 0; x < kN; ++x) {
        if (x == a || x == b || x == c) continue;
        candidates += !MissingEdge(x, a) && !MissingEdge(x, b) &&
                      !MissingEdge(x, c);
      }
      constraints.push_back({triple, candidates - 8});
    }
  std::map<int, int> histogram;
  for (const auto row : constraints) ++histogram[row.upper];
  if (constraints.size() != 420 ||
      histogram != std::map<int, int>{{0,30},{1,156},{2,84},{3,103},
                                      {4,2},{5,36},{6,9}}) {
    throw std::runtime_error("independent constraint census changed");
  }
}

bool IntervalFeasible(const Family& family) {
  for (const auto constraint : constraints) {
    int count = 0;
    for (const auto face : family) {
      count += (face & constraint.mask) == constraint.mask;
    }
    if (count > constraint.upper) return false;
  }
  return true;
}

void BuildGroundData() {
  const auto all_fours = Masks(4);
  for (const auto face : all_fours) {
    bool feasible = true;
    for (const auto constraint : constraints) {
      if ((face & constraint.mask) == constraint.mask && constraint.upper == 0) {
        feasible = false;
        break;
      }
    }
    if (feasible) fours.push_back(face);
  }
  fives = Masks(5);
  if (all_fours.size() != 771 || fours.size() != 663 || fives.size() != 1296)
    throw std::runtime_error("independent ground sets changed");
  four_index.fill(-1);
  for (int i = 0; i < static_cast<int>(fours.size()); ++i) four_index[fours[i]] = i;
  five_word_count = (fives.size() + 63) / 64;
  shadows.assign(fours.size(), Word(five_word_count));
  for (auto& word : vertex_five_words) word.assign(five_word_count, 0);
  for (int j = 0; j < static_cast<int>(fives.size()); ++j) {
    for (int i = 0; i < static_cast<int>(fours.size()); ++i) {
      if ((fives[j] & fours[i]) == fours[i])
        shadows[i][j / 64] |= std::uint64_t{1} << (j % 64);
    }
    for (int v = 0; v < kN; ++v) if ((fives[j] >> v) & 1)
      vertex_five_words[v][j / 64] |= std::uint64_t{1} << (j % 64);
  }
  std::array<std::vector<std::uint16_t>, 6> independent;
  for (int size = 2; size <= 5; ++size) independent[size] = Masks(size);
  for (int v = 0; v < kN; ++v) {
    std::array<int, 6> count{};
    for (int size = 2; size <= 5; ++size)
      for (const auto face : independent[size]) count[size] += (face >> v) & 1;
    const int f0 = count[2], f1 = count[3], f2 = count[4], f3 = count[5];
    base_g3[v] = f2 - 7 * f1 + 28 * f0 - 84;
    base_g4[v] = f3 - 6 * f2 + 21 * f1 - 56 * f0 + 126;
  }
  if (base_g3 != std::array<int,kN>{2,2,7,7,7,7,7,7,7,7,7,5,5,4} ||
      base_g4 != std::array<int,kN>{2,2,8,8,8,8,8,8,8,8,8,6,6,5})
    throw std::runtime_error("independent local g data changed");
}

GateKey CheckGate(const Family& family) {
  Word shadow(five_word_count);
  std::array<int, kRows> indices{};
  for (int row = 0; row < kRows; ++row) {
    indices[row] = four_index[family[row]];
    for (int block = 0; block < five_word_count; ++block)
      shadow[block] |= shadows[indices[row]][block];
  }
  int shadow_count = 0;
  for (const auto block : shadow) shadow_count += std::popcount(block);
  ++shadow_histogram[shadow_count];
  const int additional = kScheduledFives - shadow_count;
  bool strict = additional >= 0;
  int lower_sum = 0;
  int upper_sum = 0;
  for (int v = 0; v < kN; ++v) {
    int incidence = 0;
    for (const auto face : family) incidence += (face >> v) & 1;
    int forced = 0;
    for (int block = 0; block < five_word_count; ++block)
      forced += std::popcount(shadow[block] & vertex_five_words[v][block]);
    const int previous = base_g3[v] - incidence;
    const int before = base_g4[v] + 6 * incidence - forced;
    strict = strict && previous >= 1 && before >= 1;
    const int lower = previous >= 0
        ? std::max(0, before - MacaulayRaise(previous, 3)) : 0;
    const int upper = before - 1;
    lower_sum += lower;
    upper_sum += upper;
    if (upper < 0) strict = false;
  }
  const bool aggregate = strict && lower_sum <= 5 * additional &&
                         5 * additional <= upper_sum;
  return {strict, aggregate, shadow_count, additional};
}

bool ActiveConstraintsHold(const std::array<std::uint16_t, kRows>& active) {
  const auto active_union = std::accumulate(
      kActive.begin(), kActive.end(), std::uint16_t{0},
      [](std::uint16_t value, int v) {
        return value | static_cast<std::uint16_t>(1u << v);
      });
  for (const auto constraint : constraints) {
    if (constraint.mask & ~active_union) continue;
    int count = 0;
    for (const auto row : active)
      count += (row & constraint.mask) == constraint.mask;
    if (count > constraint.upper) return false;
  }
  return true;
}

void VisitTerminal(const std::array<std::uint16_t, kRows>& active,
                   const std::vector<int>& columns) {
  ++terminal_column_paths;
  Family family = active;
  for (int isolate = 0; isolate < static_cast<int>(columns.size()); ++isolate)
    for (int row = 0; row < kRows; ++row)
      if ((columns[isolate] >> row) & 1)
        family[row] |= std::uint16_t{1} << kIsolates[isolate];
  if (std::set<std::uint16_t>(family.begin(), family.end()).size() != kRows) return;
  ++distinct_row_paths;
  for (const auto face : family) if (four_index[face] < 0) return;
  ++ground_row_paths;
  if (!IntervalFeasible(family)) return;
  ++interval_paths;
  const auto gate = CheckGate(family);
  if (gate.shadow <= kScheduledFives) ++shadow_feasible_paths;
  ++gate_histogram[gate];
  if (gate.strict && gate.aggregate) ++viable_paths;
}

void EnumerateColumns(const std::array<std::uint16_t, kRows>& active,
                      int pattern, int left,
                      std::array<int, kRows> remaining,
                      std::vector<int>& columns) {
  if (pattern == (1 << kRows)) {
    if (left == 0 && std::ranges::all_of(remaining,
          [](int value) { return value == 0; })) VisitTerminal(active, columns);
    return;
  }
  int maximum = left;
  for (int row = 0; row < kRows; ++row)
    if ((pattern >> row) & 1) maximum = std::min(maximum, remaining[row]);
  for (int multiplicity = 0; multiplicity <= maximum; ++multiplicity) {
    auto next = remaining;
    const int next_left = left - multiplicity;
    bool feasible = true;
    for (int row = 0; row < kRows; ++row) {
      next[row] -= multiplicity * ((pattern >> row) & 1);
      feasible = feasible && next[row] >= 0 && next[row] <= next_left;
    }
    if (!feasible) continue;
    for (int copy = 0; copy < multiplicity; ++copy) columns.push_back(pattern);
    EnumerateColumns(active, pattern + 1, next_left, next, columns);
    for (int copy = 0; copy < multiplicity; ++copy) columns.pop_back();
  }
}

void EnumerateActive(const std::vector<std::uint16_t>& active_masks,
                     int row, int start,
                     std::array<std::uint16_t, kRows>& active) {
  if (row == kRows) {
    ++active_multisets;
    if (!ActiveConstraintsHold(active)) return;
    ++active_multisets_passing;
    std::array<int, kRows> remaining{};
    for (int i = 0; i < kRows; ++i) remaining[i] = 4 - std::popcount(active[i]);
    std::vector<int> columns;
    EnumerateColumns(active, 0, kIsolates.size(), remaining, columns);
    return;
  }
  for (int index = start; index < static_cast<int>(active_masks.size()); ++index) {
    active[row] = active_masks[index];
    EnumerateActive(active_masks, row + 1, index, active);
  }
}

void PrintMap(const std::map<int, std::uint64_t>& values) {
  bool first = true;
  for (const auto& [key, value] : values) {
    if (!first) std::cout << ',';
    std::cout << '\"' << key << "\":" << value;
    first = false;
  }
}

}  // namespace

int main() {
  try {
    BuildConstraints();
    BuildGroundData();
    const auto active_vertex_mask = std::accumulate(
        kActive.begin(), kActive.end(), std::uint16_t{0},
        [](std::uint16_t value, int v) {
          return value | static_cast<std::uint16_t>(1u << v);
        });
    std::vector<std::uint16_t> active_masks;
    for (std::uint16_t mask = 0; mask < (1u << kN); ++mask)
      if ((mask & ~active_vertex_mask) == 0 && Independent(mask))
        active_masks.push_back(mask);
    if (active_masks.size() != 13) throw std::runtime_error("active masks changed");
    std::array<std::uint16_t, kRows> active{};
    EnumerateActive(active_masks, 0, 0, active);
    if (active_multisets != 6188 || active_multisets_passing != 1000 ||
        terminal_column_paths != 15889463 || distinct_row_paths != 15128612 ||
        ground_row_paths != 15128612 || interval_paths != 12797343 ||
        shadow_feasible_paths != 31078 || viable_paths != 0) {
      throw std::runtime_error("independent record-31 a=5 census changed");
    }
    std::cout << "{\n"
              << "  \"kind\":\"h9_14_record31_a5_independent_direct_matrix_audit\",\n"
              << "  \"compact_only\":true,\n"
              << "  \"active_mask_count\":" << active_masks.size() << ",\n"
              << "  \"active_multisets\":" << active_multisets << ",\n"
              << "  \"active_multisets_passing_active_caps\":" << active_multisets_passing << ",\n"
              << "  \"terminal_column_paths\":" << terminal_column_paths << ",\n"
              << "  \"distinct_row_paths\":" << distinct_row_paths << ",\n"
              << "  \"ground_row_paths\":" << ground_row_paths << ",\n"
              << "  \"interval_paths\":" << interval_paths << ",\n"
              << "  \"shadow_feasible_paths\":" << shadow_feasible_paths << ",\n"
              << "  \"viable_paths\":" << viable_paths << ",\n"
              << "  \"shadow_histogram\":{";
    PrintMap(shadow_histogram);
    std::cout << "},\n  \"gate_histogram\":[";
    bool first = true;
    for (const auto& [gate, count] : gate_histogram) {
      if (!first) std::cout << ',';
      std::cout << "{\"strict\":" << (gate.strict ? "true" : "false")
                << ",\"aggregate\":" << (gate.aggregate ? "true" : "false")
                << ",\"shadow\":" << gate.shadow
                << ",\"additional\":" << gate.additional
                << ",\"path_count\":" << count << '}';
      first = false;
    }
    std::cout << "],\n"
              << "  \"a5_excluded_by_necessary_gate\":"
              << (viable_paths == 0 ? "true" : "false") << "\n}\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
