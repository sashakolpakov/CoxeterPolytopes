// Lossless positive-A4 scout for the H9/14 record-31 graph shadow.
//
// This program derives the complete Dehn--Sommerville schedule and every
// vertex/edge/triple link interval directly from the five missing edges.  It
// enumerates all interval-feasible A4 families for every positive schedule
// surviving the summed strict-link bound, applies the forced-A5-shadow
// strict-link/Macaulay necessary gate, and quotients the survivors by the
// complete graph automorphism group.  It is a compact-only incidence scout:
// additional minimal A5s, topology, Coxeter labels, and finite-volume
// noncompact chambers are not decided here.

#include <algorithm>
#include <array>
#include <bit>
#include <compare>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace {

constexpr int kN = 14;
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
  int lower;
  int upper;
};

struct VertexLinkRow {
  int vertex;
  int link_vertices;
  int inherited_missing_edges;
  int candidate_triangles;
  int g1;
  int g2;
  int g3_without_missing_fours;
  int lower;
  int upper;
};

struct EdgeLinkRow {
  int left;
  int right;
  int link_vertices;
  int candidate_edges;
  int g1;
  int g2_without_missing_fours;
  int lower;
  int upper;
};

struct TripleLinkRow {
  int first;
  int second;
  int third;
  int candidate_vertices;
  int lower;
  int upper;
};

struct Schedule {
  int missing_fours;
  int missing_fives;
  int global_g4;
  std::array<int, 9> f_vector;
  std::array<int, 10> h_vector;
  int summed_vertex_link_g4;
  bool strict_vertex_sum_feasible;
};

struct Gate {
  bool strict;
  bool aggregate;
  int shadow;
  int additional;
  auto operator<=>(const Gate&) const = default;
};

using Family = std::vector<std::uint16_t>;
using Word = std::vector<std::uint64_t>;

std::vector<Constraint> g_constraints;
std::vector<VertexLinkRow> g_vertex_link_rows;
std::vector<EdgeLinkRow> g_edge_link_rows;
std::vector<TripleLinkRow> g_triple_link_rows;
std::vector<std::uint16_t> g_fours;
std::vector<std::uint16_t> g_fives;
std::vector<Word> g_shadow_words;
std::array<Word, kN> g_vertex_words;
std::array<int, kN> g_base_g3{};
std::array<int, kN> g_base_g4{};
std::vector<std::array<int, kN>> g_automorphisms;
std::vector<Schedule> g_schedules;
int g_five_word_count = 0;

long long Binomial(int n, int k) {
  if (k < 0 || k > n) return 0;
  k = std::min(k, n - k);
  long long answer = 1;
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

bool MissingEdge(int left, int right) {
  const auto value = static_cast<std::uint16_t>(
      (std::uint16_t{1} << left) | (std::uint16_t{1} << right));
  return std::find(kEdges.begin(), kEdges.end(), value) != kEdges.end();
}

bool Independent(std::uint16_t value) {
  for (const auto edge : kEdges) {
    if ((value & edge) == edge) return false;
  }
  return true;
}

std::vector<std::uint16_t> MasksOfSize(int size) {
  std::vector<std::uint16_t> answer;
  for (std::uint16_t value = 0; value < (1u << kN); ++value) {
    if (std::popcount(value) == size && Independent(value)) answer.push_back(value);
  }
  return answer;
}

std::pair<int, int> StrictInterval(int previous, int beta,
                                   int candidate_count, int degree) {
  const int lower = std::max(0, beta - MacaulayRaise(previous, degree));
  const int upper = std::min(candidate_count, beta - 1);
  if (lower > upper) throw std::runtime_error("empty local interval");
  return {lower, upper};
}

void InitializeLocalConstraints() {
  // Vertex links.
  for (int vertex = 0; vertex < kN; ++vertex) {
    std::vector<int> link;
    for (int other = 0; other < kN; ++other) {
      if (other != vertex && !MissingEdge(vertex, other)) link.push_back(other);
    }
    const int n = static_cast<int>(link.size());
    int inherited_missing = 0;
    for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) {
      inherited_missing += MissingEdge(link[i], link[j]);
    }
    const int f1 = static_cast<int>(Binomial(n, 2)) - inherited_missing;
    int f2 = 0;
    for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j)
      for (int k = j + 1; k < n; ++k) {
        const auto value = static_cast<std::uint16_t>(
            (std::uint16_t{1} << link[i]) |
            (std::uint16_t{1} << link[j]) |
            (std::uint16_t{1} << link[k]));
        f2 += Independent(value);
      }
    const int g1 = n - 9;
    const int g2 = f1 - 8 * n + 36;
    const int beta = f2 - 7 * f1 + 28 * n - 84;
    if (!(g1 > 0 && g2 >= 1 && g2 <= MacaulayRaise(g1, 1))) {
      throw std::runtime_error("vertex-link prerequisite changed");
    }
    const auto interval = StrictInterval(g2, beta, f2, 2);
    g_constraints.push_back({static_cast<std::uint16_t>(1u << vertex),
                             interval.first, interval.second});
    g_vertex_link_rows.push_back({vertex, n, inherited_missing, f2, g1, g2,
                                  beta, interval.first, interval.second});
  }

  // Edge links for every meeting pair.
  for (int left = 0; left < kN; ++left) for (int right = left + 1; right < kN; ++right) {
    if (MissingEdge(left, right)) continue;
    std::vector<int> link;
    for (int other = 0; other < kN; ++other) {
      if (other != left && other != right &&
          !MissingEdge(other, left) && !MissingEdge(other, right)) {
        link.push_back(other);
      }
    }
    const int n = static_cast<int>(link.size());
    int f1 = 0;
    for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) {
      f1 += !MissingEdge(link[i], link[j]);
    }
    const int g1 = n - 8;
    const int beta = f1 - 7 * n + 28;
    if (g1 <= 0) throw std::runtime_error("edge-link prerequisite changed");
    const auto interval = StrictInterval(g1, beta, f1, 1);
    g_constraints.push_back({static_cast<std::uint16_t>(
                                 (std::uint16_t{1} << left) |
                                 (std::uint16_t{1} << right)),
                             interval.first, interval.second});
    g_edge_link_rows.push_back({left, right, n, f1, g1, beta,
                                interval.first, interval.second});
  }

  // Triple links.  The rank-six strict g1 bound is simply the cap below.
  for (int a = 0; a < kN; ++a) for (int b = a + 1; b < kN; ++b)
    for (int c = b + 1; c < kN; ++c) {
      const auto triple = static_cast<std::uint16_t>(
          (std::uint16_t{1} << a) |
          (std::uint16_t{1} << b) |
          (std::uint16_t{1} << c));
      if (!Independent(triple)) continue;
      int candidates = 0;
      for (int other = 0; other < kN; ++other) {
        if (other == a || other == b || other == c) continue;
        candidates += !MissingEdge(other, a) && !MissingEdge(other, b) &&
                      !MissingEdge(other, c);
      }
      g_constraints.push_back({triple, 0, candidates - 8});
      g_triple_link_rows.push_back({a, b, c, candidates, 0, candidates - 8});
    }
  std::map<int, int> cap_histogram;
  for (const auto& row : g_constraints) {
    if (row.lower != 0) throw std::runtime_error("positive local lower bound appeared");
    ++cap_histogram[row.upper];
  }
  if (g_constraints.size() != 420 ||
      cap_histogram != std::map<int, int>{{0, 30}, {1, 156}, {2, 84},
                                          {3, 103}, {4, 2}, {5, 36},
                                          {6, 9}}) {
    throw std::runtime_error("record-31 local constraint table changed");
  }
  if (g_vertex_link_rows.size() != 14 || g_edge_link_rows.size() != 87 ||
      g_triple_link_rows.size() != 319) {
    throw std::runtime_error("record-31 exact link table sizes changed");
  }
}

bool IntervalFeasible(const Family& family) {
  for (const auto& constraint : g_constraints) {
    int incidence = 0;
    for (const auto face : family) {
      incidence += (face & constraint.mask) == constraint.mask;
    }
    if (incidence < constraint.lower || incidence > constraint.upper) return false;
  }
  return true;
}

int WordCount(const Word& word) {
  int answer = 0;
  for (const auto block : word) answer += std::popcount(block);
  return answer;
}

Word ShadowUnion(const std::vector<int>& indices) {
  Word answer(g_five_word_count);
  for (const int index : indices) {
    for (int block = 0; block < g_five_word_count; ++block) {
      answer[block] |= g_shadow_words[index][block];
    }
  }
  return answer;
}

int WordIntersectionCount(const Word& left, const Word& right) {
  int answer = 0;
  for (int block = 0; block < g_five_word_count; ++block) {
    answer += std::popcount(left[block] & right[block]);
  }
  return answer;
}

Gate ShadowGate(const std::vector<int>& indices, int total_fives) {
  const Word shadow = ShadowUnion(indices);
  const int shadow_count = WordCount(shadow);
  const int additional = total_fives - shadow_count;
  bool strict = additional >= 0;
  int lower_sum = 0;
  int upper_sum = 0;
  for (int vertex = 0; vertex < kN; ++vertex) {
    int a4_incidence = 0;
    for (const int index : indices) {
      a4_incidence += (g_fours[index] >> vertex) & 1;
    }
    const int previous = g_base_g3[vertex] - a4_incidence;
    const int forced = WordIntersectionCount(shadow, g_vertex_words[vertex]);
    const int before = g_base_g4[vertex] + 6 * a4_incidence - forced;
    strict = strict && previous >= 1 && before >= 1;
    const int lower = previous >= 0
                          ? std::max(0, before - MacaulayRaise(previous, 3))
                          : 0;
    const int upper = before - 1;
    lower_sum += lower;
    upper_sum += upper;
    if (upper < 0) strict = false;
  }
  const bool aggregate = strict && lower_sum <= 5 * additional &&
                         5 * additional <= upper_sum;
  return {strict, aggregate, shadow_count, additional};
}

std::uint16_t Transform(std::uint16_t value, const std::array<int, kN>& action) {
  std::uint16_t answer = 0;
  for (int vertex = 0; vertex < kN; ++vertex) {
    if ((value >> vertex) & 1) answer |= std::uint16_t{1} << action[vertex];
  }
  return answer;
}

void InitializeAutomorphisms() {
  auto images = kActive;
  const std::set<std::uint16_t> expected(kEdges.begin(), kEdges.end());
  do {
    std::array<int, kN> action{};
    std::iota(action.begin(), action.end(), 0);
    for (int index = 0; index < static_cast<int>(kActive.size()); ++index) action[kActive[index]] = images[index];
    std::set<std::uint16_t> moved;
    for (const auto edge : kEdges) moved.insert(Transform(edge, action));
    if (moved == expected) g_automorphisms.push_back(action);
  } while (std::next_permutation(images.begin(), images.end()));
  if (g_automorphisms.size() != 2) {
    throw std::runtime_error("record-31 active automorphism order changed");
  }
}

std::vector<int> OrbitKey(const Family& family) {
  const std::uint16_t active_mask = std::accumulate(
      kActive.begin(), kActive.end(), std::uint16_t{0},
      [](std::uint16_t value, int vertex) {
        return value | (std::uint16_t{1} << vertex);
      });
  std::vector<int> best;
  bool first = true;
  for (const auto& action : g_automorphisms) {
    Family moved;
    for (const auto face : family) moved.push_back(Transform(face, action));
    std::vector<int> order(family.size());
    std::iota(order.begin(), order.end(), 0);
    do {
      std::vector<int> key;
      for (const int row : order) key.push_back(moved[row] & active_mask);
      std::vector<int> columns;
      for (const int isolate : kIsolates) {
        int column = 0;
        for (int row = 0; row < static_cast<int>(order.size()); ++row) {
          column |= ((moved[order[row]] >> isolate) & 1) << row;
        }
        columns.push_back(column);
      }
      std::sort(columns.begin(), columns.end());
      key.insert(key.end(), columns.begin(), columns.end());
      if (first || key < best) {
        best = key;
        first = false;
      }
    } while (std::next_permutation(order.begin(), order.end()));
  }
  return best;
}

void InitializeA4AndShadowData() {
  const auto all_fours = MasksOfSize(4);
  for (const auto face : all_fours) {
    if (IntervalFeasible({face})) g_fours.push_back(face);
  }
  if (all_fours.size() != 771 || g_fours.size() != 663) {
    throw std::runtime_error("record-31 A4 ground set changed");
  }
  g_fives = MasksOfSize(5);
  if (g_fives.size() != 1296) {
    throw std::runtime_error("record-31 A5 ground set changed");
  }
  g_five_word_count = (static_cast<int>(g_fives.size()) + 63) / 64;
  g_shadow_words.assign(g_fours.size(), Word(g_five_word_count));
  for (auto& word : g_vertex_words) word.assign(g_five_word_count, 0);
  for (int index = 0; index < static_cast<int>(g_fives.size()); ++index) {
    const auto five = g_fives[index];
    for (int row = 0; row < static_cast<int>(g_fours.size()); ++row) {
      if ((five & g_fours[row]) == g_fours[row]) {
        g_shadow_words[row][index / 64] |= std::uint64_t{1} << (index % 64);
      }
    }
    for (int vertex = 0; vertex < kN; ++vertex) {
      if ((five >> vertex) & 1) {
        g_vertex_words[vertex][index / 64] |= std::uint64_t{1} << (index % 64);
      }
    }
  }

  std::array<std::vector<std::uint16_t>, 6> independent;
  for (int size = 2; size <= 5; ++size) independent[size] = MasksOfSize(size);
  for (int vertex = 0; vertex < kN; ++vertex) {
    std::array<int, 6> count{};
    for (int size = 2; size <= 5; ++size) {
      for (const auto face : independent[size]) count[size] += (face >> vertex) & 1;
    }
    const int f0 = count[2], f1 = count[3], f2 = count[4], f3 = count[5];
    g_base_g3[vertex] = f2 - 7 * f1 + 28 * f0 - 84;
    g_base_g4[vertex] = f3 - 6 * f2 + 21 * f1 - 56 * f0 + 126;
  }
  if (g_base_g3 != std::array<int, kN>{2,2,7,7,7,7,7,7,7,7,7,5,5,4} ||
      g_base_g4 != std::array<int, kN>{2,2,8,8,8,8,8,8,8,8,8,6,6,5}) {
    throw std::runtime_error("record-31 base local g-vectors changed");
  }
}

void DeriveSchedules() {
  constexpr int dimension = 9;
  const int independent_two = static_cast<int>(MasksOfSize(2).size());
  const int independent_three = static_cast<int>(MasksOfSize(3).size());
  const int independent_four = static_cast<int>(MasksOfSize(4).size());
  const int independent_five = static_cast<int>(MasksOfSize(5).size());
  const int h1 = kN - dimension;
  const int h2 = independent_two - 8 * kN + 36;
  const int h3 = independent_three - 7 * independent_two + 28 * kN - 84;
  const int h4_without_missing_fours =
      independent_four - 6 * independent_three + 21 * independent_two -
      56 * kN + 126;

  for (int missing_fours = 0; missing_fours <= independent_four;
       ++missing_fours) {
    const int h4 = h4_without_missing_fours - missing_fours;
    const int global_g4 = h4 - h3;
    // The strict GLBT/Lanner argument rules out g4=0, so this is the complete
    // Dehn--Sommerville schedule interval relevant to a compact Coxeter nerve.
    if (global_g4 < 1) break;
    const std::array<int, 10> h = {
        1, h1, h2, h3, h4, h4, h3, h2, h1, 1};
    std::array<int, 9> f{};
    for (int size = 1; size <= dimension; ++size) {
      for (int index = 0; index <= size; ++index) {
        f[size - 1] += static_cast<int>(
            Binomial(dimension - index, dimension - size) * h[index]);
      }
    }
    const int missing_fives = independent_five - f[4];
    const int summed_vertex_link_g4 =
        5 * f[4] - 24 * f[3] + 63 * f[2] - 112 * f[1] + 126 * kN;
    g_schedules.push_back({missing_fours, missing_fives, global_g4, f, h,
                           summed_vertex_link_g4,
                           summed_vertex_link_g4 >= kN});
  }
  if (independent_two != 87 || independent_three != 319 ||
      independent_four != 771 || independent_five != 1296 ||
      g_schedules.size() != 8) {
    throw std::runtime_error("record-31 graph-derived schedule data changed");
  }
}

const Schedule& ScheduleFor(int missing_fours) {
  const auto found = std::find_if(
      g_schedules.begin(), g_schedules.end(),
      [missing_fours](const Schedule& row) {
        return row.missing_fours == missing_fours;
      });
  if (found == g_schedules.end()) {
    throw std::runtime_error("record-31 schedule lookup failed");
  }
  return *found;
}

struct BranchResult {
  int a;
  std::uint64_t raw;
  std::uint64_t interval_count;
  std::map<Gate, std::uint64_t> gate_histogram;
  std::vector<Family> viable;
  std::map<std::vector<int>, std::vector<Family>> orbits;
};

int Factorial(int value);

void BuildOrbits(BranchResult& result) {
  for (const auto& family : result.viable) {
    result.orbits[OrbitKey(family)].push_back(family);
  }
  const int full_group_order =
      static_cast<int>(g_automorphisms.size()) *
      Factorial(static_cast<int>(kIsolates.size()));
  std::uint64_t recovered = 0;
  for (const auto& [key, orbit] : result.orbits) {
    (void)key;
    if (orbit.empty() || full_group_order % static_cast<int>(orbit.size()) != 0) {
      throw std::runtime_error("record-31 direct orbit size is invalid");
    }
    recovered += orbit.size();
  }
  if (recovered != result.viable.size()) {
    throw std::runtime_error("record-31 direct orbit weights do not sum");
  }
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

void PrintBranch(const BranchResult& result) {
  std::map<int, int> orbit_sizes;
  for (const auto& [key, orbit] : result.orbits) ++orbit_sizes[orbit.size()];
  std::cout << "    {\"a\":" << result.a
            << ",\"raw\":" << result.raw
            << ",\"interval_feasible\":" << result.interval_count
            << ",\"strict_aggregate_viable\":" << result.viable.size()
            << ",\"orbit_count\":" << result.orbits.size()
            << ",\"orbit_size_histogram\":{";
  bool first = true;
  for (const auto& [size, count] : orbit_sizes) {
    if (!first) std::cout << ',';
    std::cout << '"' << size << "\":" << count;
    first = false;
  }
  std::cout << "},\"orbits\":[";
  first = true;
  for (const auto& [key, orbit] : result.orbits) {
    if (!first) std::cout << ',';
    std::cout << "{\"key\":[";
    for (int index = 0; index < static_cast<int>(key.size()); ++index) {
      if (index) std::cout << ',';
      std::cout << key[index];
    }
    std::cout << "],\"size\":" << orbit.size() << ",\"representative\":[";
    const auto representative = *std::min_element(orbit.begin(), orbit.end());
    for (int row = 0; row < static_cast<int>(representative.size()); ++row) {
      if (row) std::cout << ',';
      PrintMask(representative[row]);
    }
    std::cout << "]}";
    first = false;
  }
  std::cout << "],\"gate_histogram\":[";
  first = true;
  for (const auto& [gate, count] : result.gate_histogram) {
    if (!first) std::cout << ',';
    std::cout << "{\"strict\":" << (gate.strict ? "true" : "false")
              << ",\"aggregate\":" << (gate.aggregate ? "true" : "false")
              << ",\"shadow\":" << gate.shadow
              << ",\"additional\":" << gate.additional
              << ",\"count\":" << count << '}';
    first = false;
  }
  std::cout << "]}";
}

Word Signature(std::uint16_t face, const std::vector<Constraint>& rows) {
  const int blocks = (static_cast<int>(rows.size()) + 63) / 64;
  Word answer(blocks);
  for (int index = 0; index < static_cast<int>(rows.size()); ++index) {
    if ((face & rows[index].mask) == rows[index].mask) {
      answer[index / 64] |= std::uint64_t{1} << (index % 64);
    }
  }
  return answer;
}

bool TripleCapFeasible(int i, int j, int k,
                       const std::vector<Word>& cap_two_signatures) {
  for (int block = 0; block < static_cast<int>(cap_two_signatures[i].size()); ++block) {
    if (cap_two_signatures[i][block] & cap_two_signatures[j][block] &
        cap_two_signatures[k][block]) return false;
  }
  return true;
}

bool QuadCapFeasible(int i, int j, int k, int l,
                     const std::vector<Word>& cap_three_signatures) {
  for (int block = 0; block < static_cast<int>(cap_three_signatures[i].size()); ++block) {
    if (cap_three_signatures[i][block] & cap_three_signatures[j][block] &
        cap_three_signatures[k][block] & cap_three_signatures[l][block]) return false;
  }
  return true;
}

int Factorial(int value) {
  int answer = 1;
  for (int factor = 2; factor <= value; ++factor) answer *= factor;
  return answer;
}

int TransformColumn(int column, const std::vector<int>& order) {
  int answer = 0;
  for (int row = 0; row < static_cast<int>(order.size()); ++row) {
    answer |= ((column >> order[row]) & 1) << row;
  }
  return answer;
}

std::vector<int> EqualActiveColumnKey(const std::vector<std::uint16_t>& active,
                                      const std::vector<int>& columns) {
  std::vector<int> order(active.size());
  std::iota(order.begin(), order.end(), 0);
  std::vector<int> best;
  bool first = true;
  do {
    bool preserves = true;
    for (int row = 0; row < static_cast<int>(active.size()); ++row) {
      preserves = preserves && active[order[row]] == active[row];
    }
    if (!preserves) continue;
    std::vector<int> moved;
    for (const int column : columns) moved.push_back(TransformColumn(column, order));
    std::sort(moved.begin(), moved.end());
    if (first || moved < best) {
      best = std::move(moved);
      first = false;
    }
  } while (std::next_permutation(order.begin(), order.end()));
  if (first) throw std::runtime_error("equal-active column key is empty");
  return best;
}

Family FamilyFromColumns(const std::vector<std::uint16_t>& active,
                         const std::vector<int>& columns) {
  Family family(active.begin(), active.end());
  for (int isolate = 0; isolate < static_cast<int>(columns.size()); ++isolate) {
    for (int row = 0; row < static_cast<int>(active.size()); ++row) {
      if ((columns[isolate] >> row) & 1) {
        family[row] |= std::uint16_t{1} << kIsolates[isolate];
      }
    }
  }
  return family;
}

std::uint64_t IsolateOrbitSize(const std::vector<std::uint16_t>& active,
                               const std::vector<int>& columns) {
  std::map<int, int> multiplicities;
  for (const int column : columns) ++multiplicities[column];
  int column_factor = 1;
  for (const auto& [column, count] : multiplicities) column_factor *= Factorial(count);

  std::vector<int> order(active.size());
  std::iota(order.begin(), order.end(), 0);
  int compatible_row_orders = 0;
  do {
    bool preserves = true;
    for (int row = 0; row < static_cast<int>(active.size()); ++row) {
      preserves = preserves && active[order[row]] == active[row];
    }
    if (!preserves) continue;
    std::vector<int> moved;
    for (const int column : columns) moved.push_back(TransformColumn(column, order));
    std::sort(moved.begin(), moved.end());
    compatible_row_orders += moved == columns;
  } while (std::next_permutation(order.begin(), order.end()));
  const int stabilizer = column_factor * compatible_row_orders;
  if (stabilizer <= 0 || Factorial(kIsolates.size()) % stabilizer != 0) {
    throw std::runtime_error("invalid record-31 isolate-column stabilizer");
  }
  return Factorial(kIsolates.size()) / stabilizer;
}

struct QuotientBranchResult {
  int a = 4;
  std::uint64_t raw_labelled_family_count = 0;
  std::uint64_t active_row_multiset_count = 0;
  std::uint64_t column_multiset_path_count = 0;
  std::uint64_t equal_active_duplicate_path_count = 0;
  std::uint64_t repeated_full_row_path_count = 0;
  std::uint64_t distinct_matrix_orbit_count = 0;
  std::uint64_t ground_isolate_orbit_count = 0;
  std::uint64_t ground_labelled_family_count = 0;
  std::uint64_t interval_isolate_orbit_count = 0;
  std::uint64_t interval_labelled_family_count = 0;
  std::uint64_t viable_isolate_orbit_count = 0;
  std::uint64_t viable_labelled_family_count = 0;
  std::map<Gate, std::uint64_t> labelled_gate_histogram;
  std::map<std::vector<int>, Family> interval_representatives;
  std::map<std::vector<int>, std::uint64_t> full_orbit_sizes;
  std::map<std::vector<int>, Family> representatives;
};

QuotientBranchResult EnumerateA4Quotient(
    const std::vector<Word>& neighbors,
    const std::vector<Word>& cap_two_signatures,
    const std::vector<Word>& cap_three_signatures) {
  QuotientBranchResult result;
  result.raw_labelled_family_count = Binomial(static_cast<int>(g_fours.size()), 4);

  std::vector<std::uint16_t> active_masks;
  const std::uint16_t active_vertex_mask = std::accumulate(
      kActive.begin(), kActive.end(), std::uint16_t{0},
      [](std::uint16_t value, int vertex) {
        return value | (std::uint16_t{1} << vertex);
      });
  for (std::uint16_t value = 0; value < (1u << kN); ++value) {
    if ((value & ~active_vertex_mask) == 0 && Independent(value)) {
      active_masks.push_back(value);
    }
  }
  if (active_masks.size() != 13) {
    throw std::runtime_error("record-31 independent active masks changed");
  }

  std::array<int, 1 << kN> four_index;
  four_index.fill(-1);
  for (int index = 0; index < static_cast<int>(g_fours.size()); ++index) {
    four_index[g_fours[index]] = index;
  }

  for (int a0 = 0; a0 < static_cast<int>(active_masks.size()); ++a0)
    for (int a1 = a0; a1 < static_cast<int>(active_masks.size()); ++a1)
      for (int a2 = a1; a2 < static_cast<int>(active_masks.size()); ++a2)
        for (int a3 = a2; a3 < static_cast<int>(active_masks.size()); ++a3) {
          ++result.active_row_multiset_count;
          const std::vector<std::uint16_t> active = {
              active_masks[a0], active_masks[a1], active_masks[a2], active_masks[a3]};
          std::array<int, 4> sums{};
          for (int row = 0; row < 4; ++row) {
            sums[row] = 4 - std::popcount(active[row]);
          }
          std::set<std::vector<int>> seen_equal_active_keys;
          std::vector<int> columns;

          auto visit = [&](auto&& self, int pattern, int columns_left,
                           std::array<int, 4> remaining) -> void {
            if (pattern == 16) {
              if (columns_left != 0 || std::ranges::any_of(remaining,
                    [](int value) { return value != 0; })) return;
              ++result.column_multiset_path_count;
              const auto column_key = EqualActiveColumnKey(active, columns);
              if (!seen_equal_active_keys.insert(column_key).second) {
                ++result.equal_active_duplicate_path_count;
                return;
              }
              const Family family = FamilyFromColumns(active, column_key);
              if (std::set<std::uint16_t>(family.begin(), family.end()).size() != 4) {
                ++result.repeated_full_row_path_count;
                return;
              }
              ++result.distinct_matrix_orbit_count;
              const auto weight = IsolateOrbitSize(active, column_key);
              std::array<int, 4> indices{};
              for (int row = 0; row < 4; ++row) {
                indices[row] = four_index[family[row]];
                if (indices[row] < 0) return;
              }
              ++result.ground_isolate_orbit_count;
              result.ground_labelled_family_count += weight;
              for (int left = 0; left < 4; ++left) for (int right = left + 1; right < 4; ++right) {
                if (!(neighbors[indices[left]][indices[right] / 64] &
                      (std::uint64_t{1} << (indices[right] % 64)))) return;
              }
              if (!TripleCapFeasible(indices[0], indices[1], indices[2], cap_two_signatures) ||
                  !TripleCapFeasible(indices[0], indices[1], indices[3], cap_two_signatures) ||
                  !TripleCapFeasible(indices[0], indices[2], indices[3], cap_two_signatures) ||
                  !TripleCapFeasible(indices[1], indices[2], indices[3], cap_two_signatures) ||
                  !QuadCapFeasible(indices[0], indices[1], indices[2], indices[3],
                                   cap_three_signatures)) return;
              ++result.interval_isolate_orbit_count;
              result.interval_labelled_family_count += weight;
              result.interval_representatives.try_emplace(
                  OrbitKey(family), family);
              const auto gate = ShadowGate(
                  {indices[0], indices[1], indices[2], indices[3]},
                  ScheduleFor(4).missing_fives);
              result.labelled_gate_histogram[gate] += weight;
              if (!(gate.strict && gate.aggregate)) return;
              ++result.viable_isolate_orbit_count;
              result.viable_labelled_family_count += weight;
              const auto full_key = OrbitKey(family);
              result.full_orbit_sizes[full_key] += weight;
              result.representatives.try_emplace(full_key, family);
              return;
            }
            int maximum = columns_left;
            for (int row = 0; row < 4; ++row) {
              if ((pattern >> row) & 1) maximum = std::min(maximum, remaining[row]);
            }
            for (int multiplicity = 0; multiplicity <= maximum; ++multiplicity) {
              const int next_left = columns_left - multiplicity;
              auto next = remaining;
              bool feasible = true;
              for (int row = 0; row < 4; ++row) {
                next[row] -= multiplicity * ((pattern >> row) & 1);
                feasible = feasible && next[row] >= 0 && next[row] <= next_left;
              }
              if (!feasible) continue;
              for (int copy = 0; copy < multiplicity; ++copy) columns.push_back(pattern);
              self(self, pattern + 1, next_left, next);
              for (int copy = 0; copy < multiplicity; ++copy) columns.pop_back();
            }
          };
          visit(visit, 0, kIsolates.size(), sums);
        }
  if (result.ground_labelled_family_count != result.raw_labelled_family_count) {
    throw std::runtime_error("record-31 a=4 quotient does not recover C(579,4)");
  }
  if (result.viable_labelled_family_count != [&] {
        std::uint64_t total = 0;
        for (const auto& [key, size] : result.full_orbit_sizes) total += size;
        return total;
      }()) {
    throw std::runtime_error("record-31 a=4 full-orbit weights do not sum");
  }
  const std::uint64_t full_group_order = 2 * Factorial(kIsolates.size());
  for (const auto& [key, size] : result.full_orbit_sizes) {
    if (size == 0 || full_group_order % size != 0) {
      throw std::runtime_error("record-31 a=4 full orbit size is invalid");
    }
  }
  if (result.interval_representatives.size() != 62234) {
    throw std::runtime_error("record-31 a=4 full interval quotient changed");
  }
  return result;
}

void PrintQuotientBranch(const QuotientBranchResult& result) {
  std::map<std::uint64_t, int> orbit_size_histogram;
  for (const auto& [key, size] : result.full_orbit_sizes) ++orbit_size_histogram[size];
  std::cout << "    {\"a\":4"
            << ",\"raw\":" << result.raw_labelled_family_count
            << ",\"active_row_multisets\":" << result.active_row_multiset_count
            << ",\"column_paths\":" << result.column_multiset_path_count
            << ",\"equal_active_duplicate_paths\":" << result.equal_active_duplicate_path_count
            << ",\"repeated_full_row_paths\":" << result.repeated_full_row_path_count
            << ",\"distinct_matrix_orbits\":" << result.distinct_matrix_orbit_count
            << ",\"ground_isolate_column_orbits\":" << result.ground_isolate_orbit_count
            << ",\"ground_labelled_family_count\":" << result.ground_labelled_family_count
            << ",\"interval_isolate_column_orbits\":" << result.interval_isolate_orbit_count
            << ",\"interval_full_graph_orbits\":" << result.interval_representatives.size()
            << ",\"interval_feasible\":" << result.interval_labelled_family_count
            << ",\"viable_isolate_column_orbits\":" << result.viable_isolate_orbit_count
            << ",\"strict_aggregate_viable\":" << result.viable_labelled_family_count
            << ",\"orbit_count\":" << result.full_orbit_sizes.size()
            << ",\"orbit_size_histogram\":{";
  bool first = true;
  for (const auto& [size, count] : orbit_size_histogram) {
    if (!first) std::cout << ',';
    std::cout << '"' << size << "\":" << count;
    first = false;
  }
  std::cout << "},\"orbits\":[";
  first = true;
  for (const auto& [key, size] : result.full_orbit_sizes) {
    if (!first) std::cout << ',';
    std::cout << "{\"key\":[";
    for (int index = 0; index < static_cast<int>(key.size()); ++index) {
      if (index) std::cout << ',';
      std::cout << key[index];
    }
    std::cout << "],\"size\":" << size << ",\"representative\":[";
    const auto& representative = result.representatives.at(key);
    for (int row = 0; row < static_cast<int>(representative.size()); ++row) {
      if (row) std::cout << ',';
      PrintMask(representative[row]);
    }
    std::cout << "]}";
    first = false;
  }
  std::cout << "],\"gate_histogram\":[";
  first = true;
  for (const auto& [gate, count] : result.labelled_gate_histogram) {
    if (!first) std::cout << ',';
    std::cout << "{\"strict\":" << (gate.strict ? "true" : "false")
              << ",\"aggregate\":" << (gate.aggregate ? "true" : "false")
              << ",\"shadow\":" << gate.shadow
              << ",\"additional\":" << gate.additional
              << ",\"count\":" << count << '}';
    first = false;
  }
  std::cout << "]}";
}

struct A5ExtensionResult {
  std::uint64_t raw_labelled_family_count = 0;
  std::uint64_t source_a4_full_orbit_count = 0;
  std::uint64_t extension_path_count = 0;
  std::uint64_t interval_extension_path_count = 0;
  std::uint64_t shadow_feasible_extension_path_count = 0;
  std::map<int, std::uint64_t> shadow_count_histogram;
  std::map<Gate, std::uint64_t> path_gate_histogram;
  std::map<std::vector<int>, Family> viable_full_orbits;
};

A5ExtensionResult EnumerateA5Extensions(
    const QuotientBranchResult& a4_quotient) {
  A5ExtensionResult result;
  result.raw_labelled_family_count = Binomial(g_fours.size(), 5);
  result.source_a4_full_orbit_count =
      a4_quotient.interval_representatives.size();

  const int count = static_cast<int>(g_fours.size());
  const int face_word_count = (count + 63) / 64;
  std::vector<Word> containing_words(
      g_constraints.size(), Word(face_word_count));
  for (int constraint = 0;
       constraint < static_cast<int>(g_constraints.size()); ++constraint) {
    for (int face = 0; face < count; ++face) {
      if ((g_fours[face] & g_constraints[constraint].mask) ==
          g_constraints[constraint].mask) {
        containing_words[constraint][face / 64] |=
            std::uint64_t{1} << (face % 64);
      }
    }
  }

  std::array<int, 1 << kN> four_index;
  four_index.fill(-1);
  for (int index = 0; index < count; ++index) {
    four_index[g_fours[index]] = index;
  }
  const int total_fives = ScheduleFor(5).missing_fives;

  for (const auto& [source_key, source_family] :
       a4_quotient.interval_representatives) {
    (void)source_key;
    if (source_family.size() != 4 || !IntervalFeasible(source_family)) {
      throw std::runtime_error("record-31 a=5 source quotient is invalid");
    }
    std::array<int, 4> source_indices{};
    std::vector<int> incidence(g_constraints.size(), 0);
    Word forbidden(face_word_count);
    Word source_shadow(g_five_word_count);
    for (int row = 0; row < 4; ++row) {
      source_indices[row] = four_index[source_family[row]];
      if (source_indices[row] < 0) {
        throw std::runtime_error("record-31 a=5 source row is absent");
      }
      for (int block = 0; block < g_five_word_count; ++block) {
        source_shadow[block] |= g_shadow_words[source_indices[row]][block];
      }
    }
    for (int constraint = 0;
         constraint < static_cast<int>(g_constraints.size()); ++constraint) {
      for (const auto face : source_family) {
        incidence[constraint] +=
            (face & g_constraints[constraint].mask) ==
            g_constraints[constraint].mask;
      }
      if (incidence[constraint] == g_constraints[constraint].upper) {
        for (int block = 0; block < face_word_count; ++block) {
          forbidden[block] |= containing_words[constraint][block];
        }
      }
    }

    for (int candidate = 0; candidate < count; ++candidate) {
      if (std::find(source_indices.begin(), source_indices.end(), candidate) !=
          source_indices.end()) {
        continue;
      }
      ++result.extension_path_count;
      if (forbidden[candidate / 64] &
          (std::uint64_t{1} << (candidate % 64))) {
        continue;
      }
      ++result.interval_extension_path_count;
      int shadow_count = 0;
      for (int block = 0; block < g_five_word_count; ++block) {
        shadow_count += std::popcount(
            source_shadow[block] | g_shadow_words[candidate][block]);
      }
      ++result.shadow_count_histogram[shadow_count];
      if (shadow_count > total_fives) continue;
      ++result.shadow_feasible_extension_path_count;
      const std::vector<int> indices = {
          source_indices[0], source_indices[1], source_indices[2],
          source_indices[3], candidate};
      const auto gate = ShadowGate(indices, total_fives);
      ++result.path_gate_histogram[gate];
      if (!(gate.strict && gate.aggregate)) continue;
      Family family = source_family;
      family.push_back(g_fours[candidate]);
      std::sort(family.begin(), family.end());
      result.viable_full_orbits.try_emplace(OrbitKey(family), family);
    }
  }
  return result;
}

void PrintA5ExtensionBranch(const A5ExtensionResult& result) {
  std::cout << "    {\"a\":5"
            << ",\"enumeration_method\":\"extend every full-graph orbit of interval-feasible a=4 families by every interval-feasible A4\""
            << ",\"raw\":" << result.raw_labelled_family_count
            << ",\"source_a4_full_orbits\":" << result.source_a4_full_orbit_count
            << ",\"extension_paths\":" << result.extension_path_count
            << ",\"interval_extension_paths\":" << result.interval_extension_path_count
            << ",\"shadow_feasible_extension_paths\":"
            << result.shadow_feasible_extension_path_count
            << ",\"orbit_count\":" << result.viable_full_orbits.size()
            << ",\"orbits\":[";
  bool first = true;
  for (const auto& [key, family] : result.viable_full_orbits) {
    if (!first) std::cout << ',';
    std::cout << "{\"key\":[";
    for (int index = 0; index < static_cast<int>(key.size()); ++index) {
      if (index) std::cout << ',';
      std::cout << key[index];
    }
    std::cout << "],\"representative\":[";
    for (int row = 0; row < static_cast<int>(family.size()); ++row) {
      if (row) std::cout << ',';
      PrintMask(family[row]);
    }
    std::cout << "]}";
    first = false;
  }
  std::cout << "],\"shadow_count_histogram\":{";
  first = true;
  for (const auto& [shadow, count] : result.shadow_count_histogram) {
    if (!first) std::cout << ',';
    std::cout << '\"' << shadow << "\":" << count;
    first = false;
  }
  std::cout << "},\"gate_histogram\":[";
  first = true;
  for (const auto& [gate, count] : result.path_gate_histogram) {
    if (!first) std::cout << ',';
    std::cout << "{\"strict\":" << (gate.strict ? "true" : "false")
              << ",\"aggregate\":" << (gate.aggregate ? "true" : "false")
              << ",\"shadow\":" << gate.shadow
              << ",\"additional\":" << gate.additional
              << ",\"path_count\":" << count << '}';
    first = false;
  }
  std::cout << "]}";
}

template <std::size_t Size>
void PrintArray(const std::array<int, Size>& values) {
  std::cout << '[';
  for (std::size_t index = 0; index < Size; ++index) {
    if (index) std::cout << ',';
    std::cout << values[index];
  }
  std::cout << ']';
}

void PrintSchedules() {
  std::cout << "  \"schedules\":[";
  for (std::size_t index = 0; index < g_schedules.size(); ++index) {
    if (index) std::cout << ',';
    const auto& row = g_schedules[index];
    std::cout << "{\"a\":" << row.missing_fours
              << ",\"missing_independent_fives\":" << row.missing_fives
              << ",\"global_g4\":" << row.global_g4
              << ",\"f_vector\":";
    PrintArray(row.f_vector);
    std::cout << ",\"h_vector\":";
    PrintArray(row.h_vector);
    std::cout << ",\"summed_vertex_link_g4\":"
              << row.summed_vertex_link_g4
              << ",\"strict_vertex_sum_feasible\":"
              << (row.strict_vertex_sum_feasible ? "true" : "false") << '}';
  }
  std::cout << "],\n";
}

void PrintLocalLinkRows() {
  std::cout << "  \"local_link_constraints\":{\n    \"vertices\":[";
  for (std::size_t index = 0; index < g_vertex_link_rows.size(); ++index) {
    if (index) std::cout << ',';
    const auto& row = g_vertex_link_rows[index];
    std::cout << "{\"vertex\":" << row.vertex
              << ",\"link_vertex_count\":" << row.link_vertices
              << ",\"inherited_missing_edge_count\":"
              << row.inherited_missing_edges
              << ",\"candidate_link_triangle_count\":"
              << row.candidate_triangles << ",\"g1\":" << row.g1
              << ",\"g2\":" << row.g2
              << ",\"g3_without_missing_fours\":"
              << row.g3_without_missing_fours
              << ",\"missing_four_incidence_interval\":["
              << row.lower << ',' << row.upper << "]}";
  }
  std::cout << "],\n    \"edges\":[";
  for (std::size_t index = 0; index < g_edge_link_rows.size(); ++index) {
    if (index) std::cout << ',';
    const auto& row = g_edge_link_rows[index];
    std::cout << "{\"edge\":[" << row.left << ',' << row.right
              << "],\"link_vertex_count\":" << row.link_vertices
              << ",\"candidate_link_edge_count\":" << row.candidate_edges
              << ",\"g1\":" << row.g1
              << ",\"g2_without_missing_fours\":"
              << row.g2_without_missing_fours
              << ",\"missing_four_incidence_interval\":["
              << row.lower << ',' << row.upper << "]}";
  }
  std::cout << "],\n    \"triples\":[";
  for (std::size_t index = 0; index < g_triple_link_rows.size(); ++index) {
    if (index) std::cout << ',';
    const auto& row = g_triple_link_rows[index];
    std::cout << "{\"triple\":[" << row.first << ',' << row.second << ','
              << row.third << "],\"candidate_link_vertex_count\":"
              << row.candidate_vertices
              << ",\"missing_four_incidence_interval\":["
              << row.lower << ',' << row.upper << "]}";
  }
  std::map<int, int> histogram;
  for (const auto& row : g_constraints) ++histogram[row.upper];
  std::cout << "],\n    \"upper_endpoint_histogram\":{";
  bool first = true;
  for (const auto& [upper, count] : histogram) {
    if (!first) std::cout << ',';
    std::cout << '\"' << upper << "\":" << count;
    first = false;
  }
  std::cout << "}\n  },\n";
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const int max_a = argc == 1 ? 5 : std::stoi(argv[1]);
    if (max_a < 1 || max_a > 5) throw std::runtime_error("max_a must lie in 1..5");
    InitializeLocalConstraints();
    InitializeAutomorphisms();
    InitializeA4AndShadowData();
    DeriveSchedules();

    const int count = static_cast<int>(g_fours.size());
    const int neighbor_words = (count + 63) / 64;
    std::vector<Word> neighbors(count, Word(neighbor_words));
    std::vector<BranchResult> results;

    BranchResult a1{1, static_cast<std::uint64_t>(count), 0, {}, {}, {}};
    for (int i = 0; i < count; ++i) {
      ++a1.interval_count;
      const auto gate = ShadowGate({i}, ScheduleFor(1).missing_fives);
      ++a1.gate_histogram[gate];
      if (gate.strict && gate.aggregate) a1.viable.push_back({g_fours[i]});
    }
    BuildOrbits(a1);
    results.push_back(std::move(a1));
    std::cerr << "a=1 done: viable=" << results.back().viable.size()
              << " orbits=" << results.back().orbits.size() << '\n';

    if (max_a >= 2) {
      BranchResult a2{2, static_cast<std::uint64_t>(Binomial(count, 2)), 0, {}, {}, {}};
      for (int i = 0; i < count; ++i) for (int j = i + 1; j < count; ++j) {
        if (!IntervalFeasible({g_fours[i], g_fours[j]})) continue;
        ++a2.interval_count;
        neighbors[i][j / 64] |= std::uint64_t{1} << (j % 64);
        neighbors[j][i / 64] |= std::uint64_t{1} << (i % 64);
        const auto gate = ShadowGate({i, j}, ScheduleFor(2).missing_fives);
        ++a2.gate_histogram[gate];
        if (gate.strict && gate.aggregate) a2.viable.push_back({g_fours[i], g_fours[j]});
      }
      BuildOrbits(a2);
      results.push_back(std::move(a2));
      std::cerr << "a=2 done: viable=" << results.back().viable.size()
                << " orbits=" << results.back().orbits.size() << '\n';
    } else {
      // Pair compatibility is not needed for an a=1-only run.
    }

    std::vector<Constraint> cap_two;
    std::vector<Constraint> cap_three;
    for (const auto& row : g_constraints) {
      if (row.upper == 2) cap_two.push_back(row);
      if (row.upper == 3) cap_three.push_back(row);
    }
    if (cap_two.size() != 84 || cap_three.size() != 103) {
      throw std::runtime_error("record-31 cap signature sizes changed");
    }
    std::vector<Word> cap_two_signatures;
    std::vector<Word> cap_three_signatures;
    for (const auto face : g_fours) {
      cap_two_signatures.push_back(Signature(face, cap_two));
      cap_three_signatures.push_back(Signature(face, cap_three));
    }

    if (max_a >= 3) {
      BranchResult a3{3, static_cast<std::uint64_t>(Binomial(count, 3)), 0, {}, {}, {}};
      for (int i = 0; i < count; ++i) for (int j = i + 1; j < count; ++j) {
        if (!(neighbors[i][j / 64] & (std::uint64_t{1} << (j % 64)))) continue;
        for (int k = j + 1; k < count; ++k) {
          if (!(neighbors[i][k / 64] & (std::uint64_t{1} << (k % 64))) ||
              !(neighbors[j][k / 64] & (std::uint64_t{1} << (k % 64))) ||
              !TripleCapFeasible(i, j, k, cap_two_signatures)) continue;
          ++a3.interval_count;
          const auto gate = ShadowGate({i, j, k}, ScheduleFor(3).missing_fives);
          ++a3.gate_histogram[gate];
          if (gate.strict && gate.aggregate) {
            a3.viable.push_back({g_fours[i], g_fours[j], g_fours[k]});
          }
        }
      }
      BuildOrbits(a3);
      results.push_back(std::move(a3));
      std::cerr << "a=3 done: interval=" << results.back().interval_count
                << " viable=" << results.back().viable.size()
                << " orbits=" << results.back().orbits.size() << '\n';
    }

    QuotientBranchResult a4_quotient;
    const bool has_a4 = max_a >= 4;
    if (has_a4) {
      a4_quotient = EnumerateA4Quotient(
          neighbors, cap_two_signatures, cap_three_signatures);
      std::cerr << "a=4 done: interval="
                << a4_quotient.interval_labelled_family_count
                << " viable=" << a4_quotient.viable_labelled_family_count
                << " orbits=" << a4_quotient.full_orbit_sizes.size() << '\n';
    }

    A5ExtensionResult a5_extension;
    const bool has_a5 = max_a >= 5;
    if (has_a5) {
      a5_extension = EnumerateA5Extensions(a4_quotient);
      if (a5_extension.raw_labelled_family_count != 1051528580817ULL ||
          a5_extension.source_a4_full_orbit_count != 62234 ||
          a5_extension.extension_path_count != 41012206 ||
          a5_extension.interval_extension_path_count != 34615149 ||
          a5_extension.shadow_feasible_extension_path_count != 158736 ||
          !a5_extension.viable_full_orbits.empty()) {
        throw std::runtime_error("record-31 a=5 exact census changed");
      }
      std::cerr << "a=5 done: interval paths="
                << a5_extension.interval_extension_path_count
                << " shadow-feasible="
                << a5_extension.shadow_feasible_extension_path_count
                << " orbits=" << a5_extension.viable_full_orbits.size()
                << '\n';
    }

    std::cout << "{\n  \"kind\":\"h9_14_record31_positive_A4_scout\",\n"
              << "  \"graph6\":\"M?????????O?O?o??\",\n"
              << "  \"missing_edges\":[[0,11],[1,12],[0,13],[1,13]],\n"
              << "  \"theorem_certificate\":false,\n"
              << "  \"independent_replay_complete\":false,\n"
              << "  \"whole_shadow_compact_coxeter_theorem\":false,\n"
              << "  \"compact_only\":true,\n"
              << "  \"finite_volume_noncompact_chamber_excluded\":false,\n"
              << "  \"collision_free_symmetry_quotient\":true,\n"
              << "  \"quotient_key\":\"active row masks plus the sorted complete binary isolate-column patterns, minimized over every row permutation and active-graph automorphism\",\n"
              << "  \"all_positive_strict_schedules_enumerated\":"
              << (max_a == 5 ? "true" : "false") << ",\n"
              << "  \"independent_set_counts\":{\"2\":87,\"3\":319,"
                 "\"4\":771,\"5\":1296},\n";
    PrintSchedules();
    std::cout << "  \"strict_compact_surviving_a\":[";
    bool first_schedule = true;
    for (const auto& schedule : g_schedules) {
      if (!schedule.strict_vertex_sum_feasible) continue;
      if (!first_schedule) std::cout << ',';
      std::cout << schedule.missing_fours;
      first_schedule = false;
    }
    std::cout << "],\n  \"strict_compact_excluded_a\":[";
    first_schedule = true;
    for (const auto& schedule : g_schedules) {
      if (schedule.strict_vertex_sum_feasible) continue;
      if (!first_schedule) std::cout << ',';
      std::cout << schedule.missing_fours;
      first_schedule = false;
    }
    std::cout << "],\n"
              << "  \"independent_A4_count\":" << MasksOfSize(4).size() << ",\n"
              << "  \"interval_A4_count\":" << g_fours.size() << ",\n"
              << "  \"active_automorphism_order\":2,\n"
              << "  \"isolate_count\":9,\n";
    PrintLocalLinkRows();
    std::cout << "  \"branches\":[\n";
    for (int index = 0; index < static_cast<int>(results.size()); ++index) {
      PrintBranch(results[index]);
      std::cout << (index + 1 == static_cast<int>(results.size()) && !has_a4
                        ? "\n" : ",\n");
    }
    if (has_a4) {
      PrintQuotientBranch(a4_quotient);
      std::cout << (has_a5 ? ",\n" : "\n");
    }
    if (has_a5) {
      PrintA5ExtensionBranch(a5_extension);
      std::cout << '\n';
    }
    std::cout << "  ],\n"
              << "  \"scope\":\"Compact-only necessary incidence reduction: exact graph-derived schedules, local intervals, forced A5 shadow, and strict/aggregate link gates. No topology, Coxeter-label, hyperbolic-realizability, or finite-volume conclusion.\"\n"
              << "}\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
