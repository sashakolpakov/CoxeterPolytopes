// Independent positive-A4 orbit replay for H9/14 record 32.
//
// No project header or generated producer table is used.  The program derives
// the vertex-, edge-, and triple-figure intervals directly from the four
// missing edges, enumerates the 602 admissible A4s, exhausts the a=1,2,3
// interval families, applies the forced-A5-shadow/Macaulay gate, and
// canonicalizes the viable families under Aut(P4 + K2) x S8 x S_a.

#include <algorithm>
#include <array>
#include <bit>
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
constexpr std::array<int, 6> kActive = {0, 1, 2, 11, 12, 13};
constexpr std::array<int, 8> kIsolates = {3, 4, 5, 6, 7, 8, 9, 10};
constexpr std::array<std::uint16_t, 4> kEdges = {
    (1u << 0) | (1u << 11),
    (1u << 1) | (1u << 12),
    (1u << 0) | (1u << 13),
    (1u << 2) | (1u << 13),
};

struct Constraint {
  std::uint16_t mask;
  int lower;
  int upper;
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
std::vector<std::uint16_t> g_fours;
std::vector<std::uint16_t> g_fives;
std::vector<Word> g_shadow_words;
std::array<Word, kN> g_vertex_words;
std::array<int, kN> g_base_g3{};
std::array<int, kN> g_base_g4{};
std::vector<std::array<int, kN>> g_automorphisms;
int g_word_count = 0;

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
  const auto mask = (std::uint16_t{1} << left) |
                    (std::uint16_t{1} << right);
  return std::find(kEdges.begin(), kEdges.end(), mask) != kEdges.end();
}

bool Independent(std::uint16_t value) {
  for (const auto edge : kEdges) if ((value & edge) == edge) return false;
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
  // Vertex figures.
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
        const auto mask = (std::uint16_t{1} << link[i]) |
                          (std::uint16_t{1} << link[j]) |
                          (std::uint16_t{1} << link[k]);
        f2 += Independent(mask);
      }
    const int g1 = n - 9;
    const int g2 = f1 - 8 * n + 36;
    const int beta = f2 - 7 * f1 + 28 * n - 84;
    if (!(g1 > 0 && g2 >= 1 && g2 <= MacaulayRaise(g1, 1))) {
      throw std::runtime_error("vertex prerequisite changed");
    }
    const auto interval = StrictInterval(g2, beta, f2, 2);
    g_constraints.push_back({static_cast<std::uint16_t>(std::uint16_t{1} << vertex),
                             interval.first, interval.second});
  }

  // Edge figures for every meeting pair.
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
    if (g1 <= 0) throw std::runtime_error("edge prerequisite changed");
    const auto interval = StrictInterval(g1, beta, f1, 1);
    g_constraints.push_back({static_cast<std::uint16_t>(
                                 (std::uint16_t{1} << left) |
                                 (std::uint16_t{1} << right)),
                             interval.first, interval.second});
  }

  // Triple figures.
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
    }
  if (g_constraints.size() != 419 ||
      std::any_of(g_constraints.begin(), g_constraints.end(),
                  [](const auto& row) { return row.lower != 0; })) {
    throw std::runtime_error("record-32 local constraint system changed");
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

Word WordUnion(const std::vector<int>& indices) {
  Word answer(g_word_count);
  for (const int index : indices) {
    for (int block = 0; block < g_word_count; ++block) {
      answer[block] |= g_shadow_words[index][block];
    }
  }
  return answer;
}

int WordIntersectionCount(const Word& left, const Word& right) {
  int answer = 0;
  for (int block = 0; block < g_word_count; ++block) {
    answer += std::popcount(left[block] & right[block]);
  }
  return answer;
}

Gate ShadowGate(const std::vector<int>& indices, int total_fives) {
  const Word shadow = WordUnion(indices);
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

void InitializeAutomorphisms() {
  std::array<int, 6> images = kActive;
  const std::set<std::uint16_t> expected(kEdges.begin(), kEdges.end());
  do {
    std::array<int, kN> action{};
    std::iota(action.begin(), action.end(), 0);
    for (int index = 0; index < 6; ++index) action[kActive[index]] = images[index];
    std::set<std::uint16_t> moved;
    for (const auto edge : kEdges) moved.insert(Transform(edge, action));
    if (moved == expected) g_automorphisms.push_back(action);
  } while (std::next_permutation(images.begin(), images.end()));
  if (g_automorphisms.size() != 4) throw std::runtime_error("active automorphisms changed");
}

void InitializeA4AndShadowData() {
  const auto all_fours = MasksOfSize(4);
  for (const auto face : all_fours) {
    if (IntervalFeasible({face})) g_fours.push_back(face);
  }
  if (all_fours.size() != 762 || g_fours.size() != 602) {
    throw std::runtime_error("record-32 A4 ground set changed");
  }
  g_fives = MasksOfSize(5);
  g_word_count = (static_cast<int>(g_fives.size()) + 63) / 64;
  g_shadow_words.assign(g_fours.size(), Word(g_word_count));
  for (auto& word : g_vertex_words) word.assign(g_word_count, 0);
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
}

struct BranchResult {
  int a;
  std::uint64_t raw;
  std::uint64_t interval_count;
  std::map<Gate, std::uint64_t> gate_histogram;
  std::vector<Family> viable;
  std::map<std::vector<int>, std::vector<Family>> orbits;
};

void BuildOrbits(BranchResult& result) {
  for (const auto& family : result.viable) result.orbits[OrbitKey(family)].push_back(family);
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
  std::cout << "},\"orbit_keys_and_sizes\":[";
  first = true;
  for (const auto& [key, orbit] : result.orbits) {
    if (!first) std::cout << ',';
    std::cout << "{\"key\":[";
    for (int index = 0; index < static_cast<int>(key.size()); ++index) {
      if (index) std::cout << ',';
      std::cout << key[index];
    }
    std::cout << "],\"size\":" << orbit.size() << '}';
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

}  // namespace

int main() {
  try {
    InitializeLocalConstraints();
    InitializeAutomorphisms();
    InitializeA4AndShadowData();

    const int count = static_cast<int>(g_fours.size());
    const int neighbor_words = (count + 63) / 64;
    std::vector<Word> neighbors(count, Word(neighbor_words));
    std::vector<BranchResult> results;

    // a=1.
    BranchResult a1{1, static_cast<std::uint64_t>(count), 0, {}, {}, {}};
    for (int i = 0; i < count; ++i) {
      ++a1.interval_count;
      const auto gate = ShadowGate({i}, 12);
      ++a1.gate_histogram[gate];
      if (gate.strict && gate.aggregate) a1.viable.push_back({g_fours[i]});
    }
    BuildOrbits(a1);
    results.push_back(std::move(a1));

    // a=2 and the exact pair-neighbor graph used by a=3.
    BranchResult a2{2, static_cast<std::uint64_t>(Binomial(count, 2)), 0, {}, {}, {}};
    for (int i = 0; i < count; ++i) for (int j = i + 1; j < count; ++j) {
      if (!IntervalFeasible({g_fours[i], g_fours[j]})) continue;
      ++a2.interval_count;
      neighbors[i][j / 64] |= std::uint64_t{1} << (j % 64);
      neighbors[j][i / 64] |= std::uint64_t{1} << (i % 64);
      const auto gate = ShadowGate({i, j}, 18);
      ++a2.gate_histogram[gate];
      if (gate.strict && gate.aggregate) a2.viable.push_back({g_fours[i], g_fours[j]});
    }
    BuildOrbits(a2);
    results.push_back(std::move(a2));

    // Cap-two signatures.  Individual cap-zero and pair cap-one constraints
    // have already been imposed by the A4 ground set and neighbor graph.
    std::vector<Constraint> cap_two;
    for (const auto& row : g_constraints) if (row.upper == 2) cap_two.push_back(row);
    if (cap_two.size() != 130) throw std::runtime_error("cap-two table changed");
    const int cap_words = (static_cast<int>(cap_two.size()) + 63) / 64;
    std::vector<Word> signatures(count, Word(cap_words));
    for (int i = 0; i < count; ++i) {
      for (int c = 0; c < static_cast<int>(cap_two.size()); ++c) {
        if ((g_fours[i] & cap_two[c].mask) == cap_two[c].mask) {
          signatures[i][c / 64] |= std::uint64_t{1} << (c % 64);
        }
      }
    }

    BranchResult a3{3, static_cast<std::uint64_t>(Binomial(count, 3)), 0, {}, {}, {}};
    for (int i = 0; i < count; ++i) for (int j = i + 1; j < count; ++j) {
      if (!(neighbors[i][j / 64] & (std::uint64_t{1} << (j % 64)))) continue;
      Word common(cap_words);
      for (int block = 0; block < cap_words; ++block) {
        common[block] = signatures[i][block] & signatures[j][block];
      }
      for (int k = j + 1; k < count; ++k) {
        if (!(neighbors[i][k / 64] & (std::uint64_t{1} << (k % 64))) ||
            !(neighbors[j][k / 64] & (std::uint64_t{1} << (k % 64)))) continue;
        bool violates = false;
        for (int block = 0; block < cap_words; ++block) {
          if (common[block] & signatures[k][block]) {
            violates = true;
            break;
          }
        }
        if (violates) continue;
        ++a3.interval_count;
        const auto gate = ShadowGate({i, j, k}, 24);
        ++a3.gate_histogram[gate];
        if (gate.strict && gate.aggregate) {
          a3.viable.push_back({g_fours[i], g_fours[j], g_fours[k]});
        }
      }
    }
    BuildOrbits(a3);
    results.push_back(std::move(a3));

    if (results[0].interval_count != 602 || results[0].viable.size() != 322 ||
        results[0].orbits.size() != 4 ||
        results[1].interval_count != 172039 || results[1].viable.size() != 2156 ||
        results[1].orbits.size() != 6 ||
        results[2].interval_count != 30954448 || results[2].viable.size() != 56 ||
        results[2].orbits.size() != 1) {
      throw std::runtime_error("record-32 independent positive orbit census changed");
    }

    std::cout << "{\n  \"kind\":\"independent_h9_14_record32_positive_orbits\",\n"
              << "  \"independent_A4_count\":762,\n"
              << "  \"interval_A4_count\":602,\n  \"branches\":[\n";
    for (int index = 0; index < 3; ++index) {
      PrintBranch(results[index]);
      std::cout << (index == 2 ? "\n" : ",\n");
    }
    std::cout << "  ]\n}\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
