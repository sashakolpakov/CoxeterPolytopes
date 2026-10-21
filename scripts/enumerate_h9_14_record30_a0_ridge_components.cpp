// Solver-independent connected ridge-component census for record-30 a=0.
//
// Candidate facet complements are the 326 five-covers of 3K2.  In a pure
// simple nine-polytope boundary, every selected complement C lies in nine
// six-set ridge groups C+v, and each such group contains exactly two selected
// complements.  Starting with the canonical type-(3 active + 2 isolate)
// cover {0,1,2,3,4}, this recursion chooses the unique partner in the most
// constrained open ridge.  Every connected target is reached; 132 selected
// complements are then dualized literally and filtered to the a=0 schedule.

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
constexpr int TARGET = 132;
constexpr std::uint16_t FULL = (1u << N) - 1;
constexpr std::array<int, 6> ACTIVE = {0, 1, 2, 11, 12, 13};
constexpr std::array<int, 8> ISOLATE = {3, 4, 5, 6, 7, 8, 9, 10};
constexpr std::array<std::uint16_t, 3> EDGES = {
    (1u << 0) | (1u << 11),
    (1u << 1) | (1u << 12),
    (1u << 2) | (1u << 13),
};
using Mask = std::uint16_t;
using Bits = std::bitset<326>;
using Family = std::array<Mask, 8>;
struct FaceCoverage {
  std::bitset<1 << N> present;
  std::array<int, 10> counts{};
};

std::vector<Mask> g_candidates;
std::array<int, 1 << N> g_index;
std::vector<Mask> g_group_masks;
std::vector<std::vector<int>> g_group_members;
std::vector<std::array<int, 9>> g_block_groups;
std::vector<std::vector<Mask>> g_high_faces;
std::vector<std::array<int, N>> g_automorphisms;
std::set<Family> g_canonical_targets;
std::array<std::unordered_set<std::string>, TARGET + 1> g_partial_orbits;
const auto g_started = std::chrono::steady_clock::now();

struct Census {
  std::uint64_t nodes = 0;
  std::uint64_t contradictions = 0;
  std::uint64_t premature_closed = 0;
  std::uint64_t size_132_closed = 0;
  std::uint64_t wrong_schedule = 0;
  std::uint64_t target_paths = 0;
  std::uint64_t isomorphic_prunes = 0;
  std::uint64_t forced_additions = 0;
  std::uint64_t forced_completion_attempts = 0;
  std::uint64_t forced_completion_rejects = 0;
  std::map<int, std::uint64_t> closed_size_histogram;
} g_census;

bool HitsEdges(Mask value) {
  return std::ranges::all_of(EDGES, [value](Mask edge) { return value & edge; });
}

bool Independent(Mask value) {
  return std::ranges::none_of(EDGES, [value](Mask edge) { return (value & edge) == edge; });
}

void Initialize() {
  g_index.fill(-1);
  for (Mask value = 0; value < (1u << N); ++value) {
    if (std::popcount(value) == 5 && HitsEdges(value)) {
      g_index[value] = static_cast<int>(g_candidates.size());
      g_candidates.push_back(value);
    }
  }
  if (g_candidates.size() != 326) throw std::runtime_error("candidate cover count changed");

  std::map<Mask, int> group_index;
  for (int index = 0; index < static_cast<int>(g_candidates.size()); ++index) {
    const Mask complement = g_candidates[index];
    std::array<int, 9> groups{};
    int position = 0;
    for (int vertex = 0; vertex < N; ++vertex) if (!(complement & (Mask{1} << vertex))) {
      const Mask six = complement | (Mask{1} << vertex);
      auto [where, inserted] = group_index.emplace(six, group_index.size());
      if (inserted) g_group_masks.push_back(six);
      groups[position++] = where->second;
    }
    if (position != 9) throw std::runtime_error("candidate does not have nine ridge groups");
    g_block_groups.push_back(groups);
    std::vector<Mask> high_faces;
    const Mask facet = FULL ^ complement;
    Mask subset = facet;
    while (true) {
      if (std::popcount(subset) >= 5) high_faces.push_back(subset);
      if (!subset) break;
      subset = (subset - 1) & facet;
    }
    if (high_faces.size() != 256) throw std::runtime_error("high-face table changed");
    g_high_faces.push_back(std::move(high_faces));
  }
  g_group_members.resize(g_group_masks.size());
  for (int group = 0; group < static_cast<int>(g_group_masks.size()); ++group) {
    const Mask six = g_group_masks[group];
    for (int vertex = 0; vertex < N; ++vertex) if (six & (Mask{1} << vertex)) {
      const int index = g_index[six ^ (Mask{1} << vertex)];
      if (index >= 0) g_group_members[group].push_back(index);
    }
  }

  std::array<int, 6> images = {0, 1, 2, 3, 4, 5};
  const std::set<std::pair<int, int>> expected = {{0, 3}, {1, 4}, {2, 5}};
  do {
    std::set<std::pair<int, int>> moved;
    for (auto [a, b] : expected) moved.insert(std::minmax(images[a], images[b]));
    if (moved != expected) continue;
    std::array<int, N> action{};
    std::iota(action.begin(), action.end(), 0);
    for (int j = 0; j < 6; ++j) action[ACTIVE[j]] = ACTIVE[images[j]];
    g_automorphisms.push_back(action);
  } while (std::next_permutation(images.begin(), images.end()));
  if (g_automorphisms.size() != 48) throw std::runtime_error("active automorphism count changed");
}

Mask Transform(Mask value, const std::array<int, N>& action) {
  Mask out = 0;
  for (int v = 0; v < N; ++v) if (value & (Mask{1} << v)) out |= Mask{1} << action[v];
  return out;
}

Family CanonicalFamily(Family input) {
  Family best{};
  bool first = true;
  std::array<std::uint8_t, 8> original_columns{};
  for (int j = 0; j < 8; ++j) for (int row = 0; row < 8; ++row) {
    original_columns[j] |= ((input[row] >> ISOLATE[j]) & 1u) << row;
  }
  for (const auto& action : g_automorphisms) {
    std::array<Mask, 8> active_parts{};
    for (int row = 0; row < 8; ++row) active_parts[row] = Transform(input[row], action) &
        ((Mask{1} << 0) | (Mask{1} << 1) | (Mask{1} << 2) |
         (Mask{1} << 11) | (Mask{1} << 12) | (Mask{1} << 13));
    std::array<int, 8> order = {0, 1, 2, 3, 4, 5, 6, 7};
    do {
      std::array<std::uint8_t, 8> columns{};
      for (int j = 0; j < 8; ++j) for (int row = 0; row < 8; ++row) {
        columns[j] |= ((original_columns[j] >> order[row]) & 1u) << row;
      }
      std::sort(columns.begin(), columns.end());
      Family value{};
      for (int row = 0; row < 8; ++row) value[row] = active_parts[order[row]];
      for (int j = 0; j < 8; ++j) for (int row = 0; row < 8; ++row) {
        value[row] |= ((columns[j] >> row) & 1u) << ISOLATE[j];
      }
      if (first || value < best) { best = value; first = false; }
    } while (std::next_permutation(order.begin(), order.end()));
  }
  return best;
}

int LocalCap(Mask subset) {
  const int size = std::popcount(subset);
  int active = 0;
  for (int v : ACTIVE) active += (subset >> v) & 1u;
  if (size == 1) return active ? 3 : 7;
  if (size == 2) return std::array<int, 3>{7, 3, 1}[active];
  return std::array<int, 4>{6, 3, 1, 0}[active];
}

void AuditTarget(const Bits& selected) {
  std::array<bool, 1 << N> face{};
  for (int index = 0; index < 326; ++index) if (selected[index]) {
    const Mask facet = FULL ^ g_candidates[index];
    Mask subset = facet;
    while (true) {
      face[subset] = true;
      if (!subset) break;
      subset = (subset - 1) & facet;
    }
  }
  std::vector<Mask> nonfaces;
  for (Mask value = 1; value < (1u << N); ++value) {
    if (face[value]) continue;
    bool minimal = true;
    for (int v = 0; v < N; ++v) if (value & (Mask{1} << v)) minimal &= face[value ^ (Mask{1} << v)];
    if (minimal) nonfaces.push_back(value);
  }
  const std::vector<Mask> edges(EDGES.begin(), EDGES.end());
  Family fives{};
  int five_count = 0;
  bool correct = true;
  for (Mask value : nonfaces) {
    if (std::popcount(value) == 2) {
      correct = correct && std::find(edges.begin(), edges.end(), value) != edges.end();
    } else if (std::popcount(value) == 5 && five_count < 8) {
      fives[five_count++] = value;
    } else {
      correct = false;
    }
  }
  correct = correct && nonfaces.size() == 11 && five_count == 8;
  if (correct) {
    for (int size = 1; size <= 9; ++size) {
      int count = 0;
      for (Mask value = 0; value < (1u << N); ++value) count += face[value] && std::popcount(value) == size;
      const std::array<int, 9> expected = {14,88,328,806,1364,1596,1248,594,132};
      correct = correct && count == expected[size - 1];
    }
  }
  if (correct) {
    for (int size = 1; size <= 3; ++size) for (Mask subset = 1; subset < (1u << N); ++subset) {
      if (std::popcount(subset) != size || !Independent(subset)) continue;
      int incidence = 0;
      for (Mask five : fives) incidence += (five & subset) == subset;
      correct = correct && incidence <= LocalCap(subset);
    }
  }
  if (!correct) { ++g_census.wrong_schedule; return; }
  ++g_census.target_paths;
  std::sort(fives.begin(), fives.end());
  g_canonical_targets.insert(CanonicalFamily(fives));
}

bool ClosedConnected(const Bits& selected) {
  std::array<std::uint8_t, 326> degree{};
  for (const auto& members : g_group_members) {
    std::vector<int> chosen;
    for (int member : members) if (selected[member]) chosen.push_back(member);
    if (chosen.empty()) continue;
    if (chosen.size() != 2) return false;
    ++degree[chosen[0]];
    ++degree[chosen[1]];
  }
  int root = -1;
  for (int index = 0; index < 326; ++index) if (selected[index]) {
    if (degree[index] != 9) return false;
    if (root < 0) root = index;
  }
  if (root < 0) return false;
  Bits reached;
  reached.set(root);
  std::vector<int> queue = {root};
  for (std::size_t position = 0; position < queue.size(); ++position) {
    const int current = queue[position];
    for (int group : g_block_groups[current]) {
      for (int member : g_group_members[group]) {
        if (selected[member] && !reached[member]) {
          reached.set(member);
          queue.push_back(member);
        }
      }
    }
  }
  return (reached & selected).count() == selected.count();
}

// Once 1364 independent five-faces have appeared, the remaining eight are
// the complete A5 schedule: adding a facet covering any of them would exceed
// the target f-vector.  Those eight sets determine their blocker uniquely, so
// audit that sole possible completion immediately.  Returning true means the
// current branch is completely decided, whether accepted or rejected.
bool TryForcedCompletion(const Bits& selected, const FaceCoverage& coverage) {
  if (coverage.counts[5] != 1364) return false;
  ++g_census.forced_completion_attempts;
  Family missing{};
  int missing_count = 0;
  for (Mask value = 0; value < (1u << N); ++value) {
    if (std::popcount(value) != 5 || !Independent(value) || coverage.present[value]) continue;
    if (missing_count < 8) missing[missing_count] = value;
    ++missing_count;
  }
  if (missing_count != 8) {
    ++g_census.forced_completion_rejects;
    return true;
  }
  // A complement of size <5 would give a facet of dimension >8 in the dual
  // boundary, contradicting the fixed nine-dimensional schedule.
  for (Mask value = 0; value < (1u << N); ++value) {
    const int size = std::popcount(value);
    if (size < 3 || size > 4 || !HitsEdges(value)) continue;
    bool hits = true;
    for (Mask five : missing) hits = hits && (value & five);
    if (hits) {
      ++g_census.forced_completion_rejects;
      return true;
    }
  }
  Bits completion;
  for (int index = 0; index < 326; ++index) {
    bool hits = true;
    for (Mask five : missing) hits = hits && (g_candidates[index] & five);
    if (hits) completion.set(index);
  }
  if (completion.count() != TARGET || (selected & ~completion).any() ||
      !ClosedConnected(completion)) {
    ++g_census.forced_completion_rejects;
    return true;
  }
  AuditTarget(completion);
  return true;
}

// A collision-free canonical colored incidence graph for a partial family.
// Its base vertices are six active vertices, eight isolates, and three pair
// nodes encoding 3K2; one further color class contains the selected covers.
// Thus graph isomorphism is exactly Aut(3K2) x S8 on the original vertices.
std::string CanonicalPartial(const Bits& selected) {
  const int block_count = selected.count();
  const int node_count = N + 3 + block_count;
  const int words = SETWORDSNEEDED(node_count);
  std::vector<graph> input(static_cast<std::size_t>(words) * node_count);
  std::vector<graph> canonical(static_cast<std::size_t>(words) * node_count);
  EMPTYGRAPH(input.data(), words, node_count);
  for (int edge = 0; edge < 3; ++edge) {
    ADDONEEDGE(input.data(), N + edge, edge, words);
    ADDONEEDGE(input.data(), N + edge, 11 + edge, words);
  }
  int node = N + 3;
  for (int index = 0; index < 326; ++index) if (selected[index]) {
    for (int vertex = 0; vertex < N; ++vertex) {
      if (g_candidates[index] & (Mask{1} << vertex)) {
        ADDONEEDGE(input.data(), node, vertex, words);
      }
    }
    ++node;
  }
  std::vector<int> lab(node_count), ptn(node_count), orbits(node_count);
  int position = 0;
  auto add_cell = [&](const std::vector<int>& cell) {
    for (int value : cell) {
      lab[position] = value;
      ptn[position] = 1;
      ++position;
    }
    ptn[position - 1] = 0;
  };
  add_cell({0, 1, 2, 11, 12, 13});
  add_cell({3, 4, 5, 6, 7, 8, 9, 10});
  add_cell({14, 15, 16});
  std::vector<int> blocks(block_count);
  std::iota(blocks.begin(), blocks.end(), N + 3);
  add_cell(blocks);
  if (position != node_count) throw std::runtime_error("nauty partition size changed");
  DEFAULTOPTIONS_GRAPH(options);
  statsblk stats{};
  options.defaultptn = FALSE;
  options.getcanon = TRUE;
  densenauty(input.data(), lab.data(), ptn.data(), orbits.data(), &options,
             &stats, words, node_count, canonical.data());
  if (stats.errstatus != 0) throw std::runtime_error("nauty canonicalization failed");
  // Store only the canonical matching and canonical cover masks.  This is a
  // lossless encoding of the colored incidence graph and is about six times
  // smaller than its dense adjacency matrix.
  std::vector<Mask> code;
  code.reserve(3 + block_count);
  for (int pair_node = N; pair_node < N + 3; ++pair_node) {
    Mask value = 0;
    const set* row = GRAPHROW(canonical.data(), pair_node, words);
    for (int vertex = 0; vertex < N; ++vertex) {
      if (ISELEMENT(row, vertex)) value |= Mask{1} << vertex;
    }
    code.push_back(value);
  }
  std::sort(code.begin(), code.end());
  const std::size_t cover_start = code.size();
  for (int block_node = N + 3; block_node < node_count; ++block_node) {
    Mask value = 0;
    const set* row = GRAPHROW(canonical.data(), block_node, words);
    for (int vertex = 0; vertex < N; ++vertex) {
      if (ISELEMENT(row, vertex)) value |= Mask{1} << vertex;
    }
    code.push_back(value);
  }
  std::sort(code.begin() + cover_start, code.end());
  std::string key(2 * code.size(), '\0');
  for (std::size_t index = 0; index < code.size(); ++index) {
    key[2 * index] = static_cast<char>(code[index] & 0xffu);
    key[2 * index + 1] = static_cast<char>(code[index] >> 8);
  }
  return key;
}

// Saturated ridge groups determine every forbidden cover, so the complete
// logical state is the selected family alone.  Unit ridge partners are added
// to a fixed point before branching.
bool AddFacet(int index, Bits& selected, FaceCoverage& coverage) {
  if (selected[index]) return true;
  selected.set(index);
  static constexpr std::array<int, 10> expected = {
      0, 14, 88, 328, 806, 1364, 1596, 1248, 594, 132};
  for (Mask face : g_high_faces[index]) if (!coverage.present[face]) {
    coverage.present.set(face);
    const int size = std::popcount(face);
    if (++coverage.counts[size] > expected[size]) return false;
  }
  return true;
}

bool Propagate(Bits& selected, FaceCoverage& coverage, int& pivot_group,
               std::vector<int>& pivot_options) {
  while (true) {
    if (selected.count() > TARGET) return false;
    std::vector<std::uint8_t> counts(g_group_members.size());
    for (int index = 0; index < 326; ++index) if (selected[index]) {
      for (int group : g_block_groups[index]) {
        if (++counts[group] > 2) return false;
      }
    }
    Bits forbidden;
    for (int group = 0; group < static_cast<int>(counts.size()); ++group) {
      if (counts[group] == 2) {
        for (int member : g_group_members[group]) if (!selected[member]) forbidden.set(member);
      }
    }
    if ((selected & forbidden).any()) return false;
    pivot_group = -1;
    pivot_options.clear();
    int best = 100;
    int forced = -1;
    for (int group = 0; group < static_cast<int>(counts.size()); ++group) {
      if (counts[group] != 1) continue;
      std::vector<int> available;
      for (int member : g_group_members[group]) {
        if (!selected[member] && !forbidden[member]) available.push_back(member);
      }
      if (available.empty()) return false;
      if (available.size() == 1) {
        forced = available.front();
        break;
      }
      if (static_cast<int>(available.size()) < best) {
        best = available.size();
        pivot_group = group;
        pivot_options = std::move(available);
      }
    }
    if (forced < 0) return true;
    if (!AddFacet(forced, selected, coverage)) return false;
    ++g_census.forced_additions;
  }
}

void Search(Bits selected, FaceCoverage coverage) {
  ++g_census.nodes;
  int pivot = -1;
  std::vector<int> options;
  if (!Propagate(selected, coverage, pivot, options)) {
    ++g_census.contradictions;
    return;
  }
  if (selected.count() > TARGET) { ++g_census.contradictions; return; }
  if (TryForcedCompletion(selected, coverage)) return;
  auto& seen = g_partial_orbits[selected.count()];
  if (!seen.insert(CanonicalPartial(selected)).second) {
    ++g_census.isomorphic_prunes;
    return;
  }
  if (g_census.nodes % 100000 == 0) {
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - g_started).count();
    std::cerr << "nodes=" << g_census.nodes << " size=" << selected.count()
              << " iso=" << g_census.isomorphic_prunes
              << " elapsed=" << seconds << '\n';
  }
  if (pivot < 0) {
    ++g_census.closed_size_histogram[selected.count()];
    if (selected.count() != TARGET) { ++g_census.premature_closed; return; }
    ++g_census.size_132_closed;
    AuditTarget(selected);
    return;
  }
  // The pivot group currently has exactly one selected member.  Each branch
  // chooses its unique partner and excludes all alternatives immediately.
  for (int partner : options) {
    Bits next_selected = selected;
    FaceCoverage next_coverage = coverage;
    if (!AddFacet(partner, next_selected, next_coverage)) {
      ++g_census.contradictions;
      continue;
    }
    Search(next_selected, std::move(next_coverage));
  }
}

void PrintFamily(const Family& family) {
  std::cout << '[';
  for (int row = 0; row < 8; ++row) {
    if (row) std::cout << ',';
    std::cout << '[';
    bool first = true;
    for (int v = 0; v < N; ++v) if (family[row] & (Mask{1} << v)) {
      if (!first) std::cout << ',';
      std::cout << v;
      first = false;
    }
    std::cout << ']';
  }
  std::cout << ']';
}
}  // namespace

int main() {
  try {
    Initialize();
    const Mask seed = (1u << 0) | (1u << 1) | (1u << 2) | (1u << 3) | (1u << 4);
    const int seed_index = g_index[seed];
    if (seed_index < 0) throw std::runtime_error("canonical seed absent");
    Bits selected;
    FaceCoverage coverage;
    if (!AddFacet(seed_index, selected, coverage)) throw std::runtime_error("seed face overflow");
    Search(selected, coverage);
    std::cout << "{\n  \"kind\":\"h9_14_record30_a0_connected_ridge_census\",\n"
              << "  \"canonical_seed\":[0,1,2,3,4],\n"
              << "  \"candidate_complement_count\":" << g_candidates.size() << ",\n"
              << "  \"ridge_group_count\":" << g_group_masks.size() << ",\n"
              << "  \"search_nodes\":" << g_census.nodes << ",\n"
              << "  \"contradictions\":" << g_census.contradictions << ",\n"
              << "  \"premature_closed\":" << g_census.premature_closed << ",\n"
              << "  \"size_132_closed\":" << g_census.size_132_closed << ",\n"
              << "  \"wrong_schedule\":" << g_census.wrong_schedule << ",\n"
              << "  \"target_paths\":" << g_census.target_paths << ",\n"
              << "  \"isomorphic_prunes\":" << g_census.isomorphic_prunes << ",\n"
              << "  \"forced_additions\":" << g_census.forced_additions << ",\n"
              << "  \"forced_completion_attempts\":" << g_census.forced_completion_attempts << ",\n"
              << "  \"forced_completion_rejects\":" << g_census.forced_completion_rejects << ",\n"
              << "  \"target_orbit_count\":" << g_canonical_targets.size() << ",\n"
              << "  \"closed_size_histogram\":{";
    bool first = true;
    for (auto [size, count] : g_census.closed_size_histogram) {
      if (!first) std::cout << ',';
      std::cout << '\"' << size << "\":" << count;
      first = false;
    }
    std::cout << "},\n  \"target_orbit_representatives\":[";
    first = true;
    for (const auto& family : g_canonical_targets) {
      if (!first) std::cout << ',';
      PrintFamily(family);
      first = false;
    }
    std::cout << "],\n  \"scope\":\"Connected compact a=0 incidence census only; no finite-volume noncompact claim.\"\n}\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
