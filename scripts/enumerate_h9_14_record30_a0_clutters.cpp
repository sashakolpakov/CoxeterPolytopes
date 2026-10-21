// Exact orbit census of record-30 a=0 A5 clutters.
//
// The three prescribed missing edges form 3K2.  Any further minimal
// nonface is an independent five-set.  We enumerate unordered families of
// eight such A5s, quotient partial families by Aut(3K2) x S8 with nauty,
// and maintain their five-transversal blocker.  The terminal gates are:
// no transversal of size three/four, exactly 132 size-five transversals,
// the target f-vector, and a connected degree-two ridge incidence table.

#include <algorithm>
#include <array>
#include <bit>
#include <bitset>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include <nauty.h>

namespace {
constexpr int N = 14;
constexpr int ROWS = 8;
constexpr int TARGET = 132;
constexpr std::uint16_t FULL = (1u << N) - 1;
constexpr std::array<int, 6> ACTIVE = {0, 1, 2, 11, 12, 13};
constexpr std::array<std::uint16_t, 3> EDGES = {
    (1u << 0) | (1u << 11),
    (1u << 1) | (1u << 12),
    (1u << 2) | (1u << 13),
};
using Mask = std::uint16_t;
using FacetBits = std::bitset<326>;
struct SmallBits {
  std::array<std::uint64_t, 2> words{};
  void set(int index) { words[index / 64] |= std::uint64_t{1} << (index % 64); }
  bool operator[](int index) const { return words[index / 64] & (std::uint64_t{1} << (index % 64)); }
  bool any() const { return words[0] || words[1]; }
  SmallBits& operator&=(const SmallBits& other) {
    words[0] &= other.words[0];
    words[1] &= other.words[1];
    return *this;
  }
};

SmallBits operator&(const SmallBits& left, const SmallBits& right) {
  return {{left.words[0] & right.words[0], left.words[1] & right.words[1]}};
}

std::vector<Mask> g_rows;
std::array<int, 1 << N> g_row_index;
std::vector<Mask> g_facets;
std::array<int, 1 << N> g_facet_index;
std::vector<FacetBits> g_row_facets;
std::vector<std::vector<Mask>> g_row_local_subsets;
std::vector<Mask> g_small_transversals;
std::vector<SmallBits> g_row_small_hits;
std::vector<std::vector<int>> g_rows_disjoint_small;
std::array<std::vector<int>, 1 << N> g_rows_within;
std::vector<std::vector<int>> g_ridge_members;
std::array<std::unordered_set<std::string>, ROWS + 1> g_partial_orbits;
std::array<std::uint8_t, 1 << N> g_local_incidence{};
std::map<std::string, std::vector<Mask>> g_targets;
const auto g_started = std::chrono::steady_clock::now();

struct Census {
  std::uint64_t nodes = 0;
  std::uint64_t isomorphic_prunes = 0;
  std::uint64_t local_cap_rejects = 0;
  std::uint64_t blocker_below_132 = 0;
  std::uint64_t terminal = 0;
  std::uint64_t terminal_small_transversal = 0;
  std::uint64_t terminal_wrong_blocker_count = 0;
  std::uint64_t terminal_ridge = 0;
  std::uint64_t terminal_face_table = 0;
  std::uint64_t target_paths = 0;
  std::array<std::uint64_t, ROWS + 1> orbit_counts{};
} g_census;

bool HitsEdges(Mask value) {
  return std::ranges::all_of(EDGES, [value](Mask edge) { return value & edge; });
}

bool Independent(Mask value) {
  return std::ranges::none_of(EDGES, [value](Mask edge) { return (value & edge) == edge; });
}

Mask SmallUnion(const SmallBits& selected) {
  Mask result = 0;
  for (int word = 0; word < 2; ++word) {
    std::uint64_t bits = selected.words[word];
    while (bits) {
      const int bit = std::countr_zero(bits);
      const int index = 64 * word + bit;
      if (index < static_cast<int>(g_small_transversals.size())) {
        result |= g_small_transversals[index];
      }
      bits &= bits - 1;
    }
  }
  return result;
}

int SmallPivot(const SmallBits& survivors) {
  int pivot = -1;
  std::size_t best = g_rows.size() + 1;
  for (int index = 0; index < static_cast<int>(g_small_transversals.size()); ++index) {
    if (survivors[index] && g_rows_disjoint_small[index].size() < best) {
      best = g_rows_disjoint_small[index].size();
      pivot = index;
    }
  }
  return pivot;
}

int LocalCap(Mask subset) {
  const int size = std::popcount(subset);
  int active = 0;
  for (int vertex : ACTIVE) active += (subset >> vertex) & 1u;
  if (size == 1) return active ? 3 : 7;
  if (size == 2) return std::array<int, 3>{7, 3, 1}[active];
  return std::array<int, 4>{6, 3, 1, 0}[active];
}

void Initialize() {
  g_facet_index.fill(-1);
  g_row_index.fill(-1);
  for (Mask value = 0; value < (1u << N); ++value) {
    if (std::popcount(value) == 5 && Independent(value)) {
      const int row = g_rows.size();
      g_rows.push_back(value);
      g_row_index[value] = row;
      std::vector<Mask> local;
      Mask subset = value;
      while (subset) {
        if (std::popcount(subset) <= 3) local.push_back(subset);
        subset = (subset - 1) & value;
      }
      if (local.size() != 25) throw std::runtime_error("A5 local-subset count changed");
      g_row_local_subsets.push_back(std::move(local));
      (void)row;
    }
    if (std::popcount(value) == 5 && HitsEdges(value)) {
      g_facet_index[value] = g_facets.size();
      g_facets.push_back(value);
    }
    if ((std::popcount(value) == 3 || std::popcount(value) == 4) && HitsEdges(value)) {
      g_small_transversals.push_back(value);
    }
  }
  if (g_rows.size() != 1372 || g_facets.size() != 326 ||
      g_small_transversals.size() > 512) {
    throw std::runtime_error("record-30 universe count changed");
  }
  g_row_facets.resize(g_rows.size());
  g_row_small_hits.resize(g_rows.size());
  g_rows_disjoint_small.resize(g_small_transversals.size());
  for (int row = 0; row < static_cast<int>(g_rows.size()); ++row) {
    for (int index = 0; index < static_cast<int>(g_facets.size()); ++index) {
      if (g_rows[row] & g_facets[index]) g_row_facets[row].set(index);
    }
    for (int index = 0; index < static_cast<int>(g_small_transversals.size()); ++index) {
      if (g_rows[row] & g_small_transversals[index]) {
        g_row_small_hits[row].set(index);
      } else {
        g_rows_disjoint_small[index].push_back(row);
      }
    }
    const Mask remaining = FULL ^ g_rows[row];
    Mask extra = remaining;
    while (true) {
      g_rows_within[g_rows[row] | extra].push_back(row);
      if (!extra) break;
      extra = (extra - 1) & remaining;
    }
  }
  std::map<Mask, int> ridge_index;
  for (int index = 0; index < static_cast<int>(g_facets.size()); ++index) {
    for (int vertex = 0; vertex < N; ++vertex) if (!(g_facets[index] & (Mask{1} << vertex))) {
      const Mask six = g_facets[index] | (Mask{1} << vertex);
      auto [where, inserted] = ridge_index.emplace(six, ridge_index.size());
      if (inserted) g_ridge_members.emplace_back();
      g_ridge_members[where->second].push_back(index);
    }
  }
}

std::string CanonicalCode(const std::vector<int>& family) {
  const int row_count = family.size();
  const int node_count = N + 3 + row_count;
  const int words = SETWORDSNEEDED(node_count);
  std::vector<graph> input(static_cast<std::size_t>(words) * node_count);
  std::vector<graph> canonical(static_cast<std::size_t>(words) * node_count);
  EMPTYGRAPH(input.data(), words, node_count);
  for (int edge = 0; edge < 3; ++edge) {
    ADDONEEDGE(input.data(), N + edge, edge, words);
    ADDONEEDGE(input.data(), N + edge, 11 + edge, words);
  }
  for (int position = 0; position < row_count; ++position) {
    const Mask value = g_rows[family[position]];
    for (int vertex = 0; vertex < N; ++vertex) if (value & (Mask{1} << vertex)) {
      ADDONEEDGE(input.data(), N + 3 + position, vertex, words);
    }
  }
  std::vector<int> lab(node_count), ptn(node_count), orbits(node_count);
  int position = 0;
  auto add_cell = [&](const std::vector<int>& cell) {
    if (cell.empty()) return;
    for (int value : cell) {
      lab[position] = value;
      ptn[position++] = 1;
    }
    ptn[position - 1] = 0;
  };
  add_cell({0, 1, 2, 11, 12, 13});
  add_cell({3, 4, 5, 6, 7, 8, 9, 10});
  add_cell({14, 15, 16});
  std::vector<int> row_nodes(row_count);
  std::iota(row_nodes.begin(), row_nodes.end(), N + 3);
  add_cell(row_nodes);
  DEFAULTOPTIONS_GRAPH(options);
  statsblk stats{};
  options.defaultptn = FALSE;
  options.getcanon = TRUE;
  densenauty(input.data(), lab.data(), ptn.data(), orbits.data(), &options,
             &stats, words, node_count, canonical.data());
  if (stats.errstatus != 0) throw std::runtime_error("nauty canonicalization failed");
  std::vector<Mask> code;
  code.reserve(3 + row_count);
  for (int node = N; node < N + 3; ++node) {
    Mask value = 0;
    const set* neighbors = GRAPHROW(canonical.data(), node, words);
    for (int vertex = 0; vertex < N; ++vertex) if (ISELEMENT(neighbors, vertex)) {
      value |= Mask{1} << vertex;
    }
    code.push_back(value);
  }
  std::sort(code.begin(), code.end());
  const std::size_t row_start = code.size();
  for (int node = N + 3; node < node_count; ++node) {
    Mask value = 0;
    const set* neighbors = GRAPHROW(canonical.data(), node, words);
    for (int vertex = 0; vertex < N; ++vertex) if (ISELEMENT(neighbors, vertex)) {
      value |= Mask{1} << vertex;
    }
    code.push_back(value);
  }
  std::sort(code.begin() + row_start, code.end());
  std::string key(2 * code.size(), '\0');
  for (std::size_t index = 0; index < code.size(); ++index) {
    key[2 * index] = static_cast<char>(code[index] & 0xffu);
    key[2 * index + 1] = static_cast<char>(code[index] >> 8);
  }
  return key;
}

bool LocalAdmissible(int row) {
  for (Mask subset : g_row_local_subsets[row]) {
    if (g_local_incidence[subset] >= LocalCap(subset)) return false;
  }
  return true;
}

void ChangeLocal(int row, int delta) {
  for (Mask subset : g_row_local_subsets[row]) {
    g_local_incidence[subset] = static_cast<std::uint8_t>(g_local_incidence[subset] + delta);
  }
}

bool ClosedConnected(const FacetBits& selected) {
  std::array<std::uint8_t, 326> degrees{};
  std::vector<std::array<int, 2>> edges;
  for (const auto& members : g_ridge_members) {
    std::array<int, 2> chosen{};
    int count = 0;
    for (int member : members) if (selected[member]) {
      if (count < 2) chosen[count] = member;
      ++count;
    }
    if (!count) continue;
    if (count != 2) return false;
    ++degrees[chosen[0]];
    ++degrees[chosen[1]];
    edges.push_back(chosen);
  }
  int root = -1;
  for (int index = 0; index < 326; ++index) if (selected[index]) {
    if (degrees[index] != 9) return false;
    if (root < 0) root = index;
  }
  FacetBits reached;
  reached.set(root);
  bool changed = true;
  while (changed) {
    changed = false;
    for (auto [a, b] : edges) {
      if (reached[a] != reached[b]) {
        reached.set(a);
        reached.set(b);
        changed = true;
      }
    }
  }
  return (reached & selected).count() == selected.count();
}

bool FaceTableGate(const std::vector<int>& family, const FacetBits& facets) {
  std::array<bool, 1 << N> face{};
  for (int index = 0; index < 326; ++index) if (facets[index]) {
    const Mask facet = FULL ^ g_facets[index];
    Mask subset = facet;
    while (true) {
      face[subset] = true;
      if (!subset) break;
      subset = (subset - 1) & facet;
    }
  }
  constexpr std::array<int, 9> expected = {14, 88, 328, 806, 1364, 1596, 1248, 594, 132};
  for (int size = 1; size <= 9; ++size) {
    int count = 0;
    for (Mask value = 0; value < (1u << N); ++value) {
      count += face[value] && std::popcount(value) == size;
    }
    if (count != expected[size - 1]) return false;
  }
  std::set<Mask> expected_nonfaces(EDGES.begin(), EDGES.end());
  for (int row : family) expected_nonfaces.insert(g_rows[row]);
  std::set<Mask> actual_nonfaces;
  for (Mask value = 1; value < (1u << N); ++value) if (!face[value]) {
    bool minimal = true;
    for (int vertex = 0; vertex < N; ++vertex) if (value & (Mask{1} << vertex)) {
      minimal = minimal && face[value ^ (Mask{1} << vertex)];
    }
    if (minimal) actual_nonfaces.insert(value);
  }
  return actual_nonfaces == expected_nonfaces;
}

void Search(std::vector<int>& family, const FacetBits& possible,
            const SmallBits& small_survivors) {
  ++g_census.nodes;
  const int depth = family.size();
  if (depth == ROWS) {
    ++g_census.terminal;
    if (small_survivors.any()) {
      ++g_census.terminal_small_transversal;
      return;
    }
    if (possible.count() != TARGET) {
      ++g_census.terminal_wrong_blocker_count;
      return;
    }
    if (!ClosedConnected(possible)) {
      ++g_census.terminal_ridge;
      return;
    }
    if (!FaceTableGate(family, possible)) {
      ++g_census.terminal_face_table;
      return;
    }
    const std::string code = CanonicalCode(family);
    if (g_targets.contains(code)) return;
    ++g_census.target_paths;
    g_targets.emplace(code, std::vector<Mask>{});
    auto& target = g_targets.at(code);
    for (int row : family) target.push_back(g_rows[row]);
    std::sort(target.begin(), target.end());
    return;
  }
  if (depth + 1 < ROWS) {
    const std::string code = CanonicalCode(family);
    if (!g_partial_orbits[depth].insert(code).second) {
      ++g_census.isomorphic_prunes;
      return;
    }
  }
  // At depth seven it is cheaper and memory-free to enumerate the handful
  // of possible final rows directly.  Duplicates are harmless because every
  // accepted terminal family is canonically deduplicated below.
  ++g_census.orbit_counts[depth];
  if (g_census.nodes % 100000 == 0) {
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - g_started).count();
    std::cerr << "nodes=" << g_census.nodes << " depth=" << depth
              << " orbits=" << g_census.orbit_counts[depth]
              << " elapsed=" << seconds << '\n';
  }
  std::array<bool, 1372> present{};
  for (int row : family) present[row] = true;
  const int small_pivot = SmallPivot(small_survivors);
  if (depth + 2 == ROWS) {
    // Complete the last two rows directly from this quotiented six-row
    // parent.  After choosing the seventh row, the eighth must be disjoint
    // from the union of every surviving small transversal.
    auto try_penultimate = [&](int penultimate) {
      if (present[penultimate] || !LocalAdmissible(penultimate)) return;
      const FacetBits after_penultimate = possible & g_row_facets[penultimate];
      if (after_penultimate.count() < TARGET) return;
      const SmallBits surviving = small_survivors & g_row_small_hits[penultimate];
      const Mask forbidden_vertices = SmallUnion(surviving);
      const Mask allowed = FULL ^ forbidden_vertices;
      if (std::popcount(allowed) < 5) return;
      family.push_back(penultimate);
      ChangeLocal(penultimate, +1);
      auto try_final = [&](int final_row) {
        if (final_row == penultimate || present[final_row] ||
            (g_rows[final_row] & forbidden_vertices) || !LocalAdmissible(final_row)) return;
        const FacetBits completion = after_penultimate & g_row_facets[final_row];
        if (completion.count() != TARGET) return;
        family.push_back(final_row);
        ChangeLocal(final_row, +1);
        SmallBits none;
        Search(family, completion, none);
        ChangeLocal(final_row, -1);
        family.pop_back();
      };
      for (int final_row : g_rows_within[allowed]) try_final(final_row);
      ChangeLocal(penultimate, -1);
      family.pop_back();
    };
    if (small_pivot >= 0) {
      for (int penultimate : g_rows_disjoint_small[small_pivot]) try_penultimate(penultimate);
    } else {
      for (int penultimate = 0; penultimate < static_cast<int>(g_rows.size()); ++penultimate) {
        try_penultimate(penultimate);
      }
    }
    return;
  }
  Mask final_forbidden_vertices = 0;
  if (depth + 1 == ROWS) {
    final_forbidden_vertices = SmallUnion(small_survivors);
    if (std::popcount(static_cast<Mask>(FULL ^ final_forbidden_vertices)) < 5) return;
  }
  auto try_row = [&](int row) {
    if (present[row]) return;
    if (depth + 1 == ROWS && (g_rows[row] & final_forbidden_vertices)) return;
    if (!LocalAdmissible(row)) {
      ++g_census.local_cap_rejects;
      return;
    }
    const SmallBits next_small = small_survivors & g_row_small_hits[row];
    if (depth + 1 == ROWS && next_small.any()) {
      ++g_census.terminal;
      ++g_census.terminal_small_transversal;
      return;
    }
    if (depth + 2 == ROWS) {
      const Mask one_row_forbidden = SmallUnion(next_small);
      if (std::popcount(static_cast<Mask>(FULL ^ one_row_forbidden)) < 5) return;
    }
    const FacetBits next_possible = possible & g_row_facets[row];
    if (next_possible.count() < TARGET) {
      ++g_census.blocker_below_132;
      return;
    }
    if (depth + 1 == ROWS && next_possible.count() != TARGET) {
      ++g_census.terminal;
      ++g_census.terminal_wrong_blocker_count;
      return;
    }
    family.push_back(row);
    ChangeLocal(row, +1);
    Search(family, next_possible, next_small);
    ChangeLocal(row, -1);
    family.pop_back();
  };
  if (depth + 1 == ROWS &&
      std::popcount(static_cast<Mask>(FULL ^ final_forbidden_vertices)) <= 10) {
    const Mask allowed = FULL ^ final_forbidden_vertices;
    Mask subset = allowed;
    while (subset) {
      if (std::popcount(subset) == 5 && g_row_index[subset] >= 0) try_row(g_row_index[subset]);
      subset = (subset - 1) & allowed;
    }
  } else if (small_pivot >= 0) {
    for (int row : g_rows_disjoint_small[small_pivot]) try_row(row);
  } else {
    for (int row = 0; row < static_cast<int>(g_rows.size()); ++row) try_row(row);
  }
}

void PrintMask(Mask value) {
  std::cout << '[';
  bool first = true;
  for (int vertex = 0; vertex < N; ++vertex) if (value & (Mask{1} << vertex)) {
    if (!first) std::cout << ',';
    std::cout << vertex;
    first = false;
  }
  std::cout << ']';
}
}  // namespace

int main(int argc, char** argv) {
  try {
    int seed_active_count = -1;
    if (argc == 2) {
      const std::string argument = argv[1];
      const std::string prefix = "--seed-active-count=";
      if (!argument.starts_with(prefix)) throw std::runtime_error("unknown argument");
      seed_active_count = std::stoi(argument.substr(prefix.size()));
      if (seed_active_count < 0 || seed_active_count > 2) {
        throw std::runtime_error("seed active count must be 0, 1, or 2");
      }
    } else if (argc != 1) {
      throw std::runtime_error("expected at most one seed argument");
    }
    Initialize();
    FacetBits possible;
    for (int index = 0; index < static_cast<int>(g_facets.size()); ++index) possible.set(index);
    SmallBits small;
    for (int index = 0; index < static_cast<int>(g_small_transversals.size()); ++index) small.set(index);
    std::vector<int> family;
    if (seed_active_count >= 0) {
      constexpr std::array<Mask, 3> seeds = {
          (1u << 3) | (1u << 4) | (1u << 5) | (1u << 6) | (1u << 7),
          (1u << 0) | (1u << 3) | (1u << 4) | (1u << 5) | (1u << 6),
          (1u << 0) | (1u << 1) | (1u << 3) | (1u << 4) | (1u << 5),
      };
      const int row = g_row_index[seeds[seed_active_count]];
      if (row < 0) throw std::runtime_error("canonical seed absent");
      family.push_back(row);
      ChangeLocal(row, +1);
      possible &= g_row_facets[row];
      small &= g_row_small_hits[row];
    }
    Search(family, possible, small);
    std::cout << "{\n  \"kind\":\"h9_14_record30_a0_clutter_orbit_census\",\n"
              << "  \"A5_universe_count\":" << g_rows.size() << ",\n"
              << "  \"five_cover_count\":" << g_facets.size() << ",\n"
              << "  \"small_transversal_count\":" << g_small_transversals.size() << ",\n"
              << "  \"seed_active_count\":" << seed_active_count << ",\n"
              << "  \"search_nodes\":" << g_census.nodes << ",\n"
              << "  \"isomorphic_prunes\":" << g_census.isomorphic_prunes << ",\n"
              << "  \"local_cap_rejects\":" << g_census.local_cap_rejects << ",\n"
              << "  \"blocker_below_132\":" << g_census.blocker_below_132 << ",\n"
              << "  \"terminal\":" << g_census.terminal << ",\n"
              << "  \"terminal_small_transversal\":" << g_census.terminal_small_transversal << ",\n"
              << "  \"terminal_wrong_blocker_count\":" << g_census.terminal_wrong_blocker_count << ",\n"
              << "  \"terminal_ridge\":" << g_census.terminal_ridge << ",\n"
              << "  \"terminal_face_table\":" << g_census.terminal_face_table << ",\n"
              << "  \"target_paths\":" << g_census.target_paths << ",\n"
              << "  \"target_orbit_count\":" << g_targets.size() << ",\n"
              << "  \"partial_orbit_counts\":[";
    for (int depth = 0; depth <= ROWS; ++depth) {
      if (depth) std::cout << ',';
      std::cout << g_census.orbit_counts[depth];
    }
    std::cout << "],\n  \"target_orbit_representatives\":[";
    bool first_family = true;
    for (const auto& [code, family_masks] : g_targets) {
      (void)code;
      if (!first_family) std::cout << ',';
      std::cout << '[';
      for (int row = 0; row < ROWS; ++row) {
        if (row) std::cout << ',';
        PrintMask(family_masks[row]);
      }
      std::cout << ']';
      first_family = false;
    }
    std::cout << "],\n  \"scope\":\"Compact a=0 incidence census under the frozen local-cap table; no finite-volume claim.\"\n}\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
