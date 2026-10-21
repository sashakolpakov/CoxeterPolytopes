// Symmetry-quotiented exact incidence scout for Record 30, a=2 orbit 14.
//
// The fixed missing structure is 3K2 together with
//   {0,1,3,4}, {0,2,3,4}.
// The compact schedule requires five additional minimal A5 supports and 128
// facets.  A colored nauty incidence graph quotients every partial family by
// the full stabilizer of the fixed structure.  Unlike the older literal
// searches, the recursion also carries the exact surviving order-five facet
// covers and prunes below 128 immediately.
//
// This is a discovery program.  Its output still needs an independently
// implemented replay before it can support a theorem claim.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include <nauty.h>

namespace {

constexpr int kN = 14;
constexpr int kMaximumAdditions = 5;
constexpr int kScheduledFives = 20;
constexpr int kScheduledFacets = 128;
using Mask = std::uint16_t;
constexpr Mask kFull = (Mask{1} << kN) - 1;
constexpr std::array<Mask, 3> kEdges = {
    (1u << 0) | (1u << 11),
    (1u << 1) | (1u << 12),
    (1u << 2) | (1u << 13),
};
std::array<Mask, 2> g_fours = {
    (1u << 0) | (1u << 1) | (1u << 3) | (1u << 4),
    (1u << 0) | (1u << 2) | (1u << 3) | (1u << 4),
};
constexpr std::array<int, kN> kBaseG3 = {
    4, 4, 4, 8, 8, 8, 8, 8, 8, 8, 8, 4, 4, 4,
};
constexpr std::array<int, kN> kBaseG4 = kBaseG3;

struct Candidate {
  Mask face = 0;
  std::uint64_t destroys_small = 0;
  std::vector<std::uint64_t> hits_covers;
};

struct State {
  std::vector<int> chosen;
  std::uint64_t destroyed_small = 0;
  std::array<std::uint8_t, kN> incidence{};
  std::vector<std::uint64_t> surviving_covers;
};

std::vector<Mask> g_shadow;
std::vector<Mask> g_small;
std::vector<Mask> g_covers;
std::vector<Candidate> g_candidates;
std::vector<std::vector<int>> g_destroyers;
std::array<int, kN> g_upper{};
std::vector<std::unordered_set<std::string>> g_seen;
int g_orbit_index = 14;
int g_additions = 5;
std::uint64_t g_stabilizer_order = 0;
std::uint64_t g_nodes = 0;
std::uint64_t g_facet_prunes = 0;
std::uint64_t g_uncoverable_prunes = 0;
std::uint64_t g_scheduled_orbits = 0;
std::uint64_t g_pure_target_orbits = 0;
std::map<int, std::uint64_t> g_blocker_histogram;
std::vector<std::vector<Mask>> g_target_representatives;

long long Choose(int n, int k) {
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
           Choose(top + 1, lower) <= remaining) {
      ++top;
    }
    remaining -= static_cast<int>(Choose(top, lower));
    answer += static_cast<int>(Choose(top + 1, lower + 1));
    next_top = top;
  }
  if (remaining != 0) throw std::runtime_error("Macaulay expansion failed");
  return answer;
}

bool Independent(Mask value) {
  for (Mask edge : kEdges) {
    if ((value & edge) == edge) return false;
  }
  return true;
}

bool ContainsFour(Mask value) {
  for (Mask four : g_fours) {
    if ((value & four) == four) return true;
  }
  return false;
}

bool HitsAll(Mask transversal, const std::vector<Mask>& clutter) {
  for (Mask face : clutter) {
    if (!(transversal & face)) return false;
  }
  return true;
}

int BitCount(const std::vector<std::uint64_t>& words) {
  int answer = 0;
  for (std::uint64_t word : words) answer += std::popcount(word);
  return answer;
}

std::vector<std::uint64_t> Meet(
    const std::vector<std::uint64_t>& left,
    const std::vector<std::uint64_t>& right) {
  std::vector<std::uint64_t> answer(left.size());
  for (int word = 0; word < static_cast<int>(left.size()); ++word) {
    answer[word] = left[word] & right[word];
  }
  return answer;
}

bool CapCompatible(const State& state, Mask face) {
  for (int vertex = 0; vertex < kN; ++vertex) {
    if (static_cast<int>(state.incidence[vertex]) +
            static_cast<int>((face >> vertex) & 1u) >
        g_upper[vertex]) {
      return false;
    }
  }
  return true;
}

std::string CanonicalCode(const std::vector<int>& chosen) {
  const int row_count = static_cast<int>(chosen.size());
  const int edge_start = kN;
  const int four_start = edge_start + static_cast<int>(kEdges.size());
  const int row_start = four_start + static_cast<int>(g_fours.size());
  const int node_count = row_start + row_count;
  const int words = SETWORDSNEEDED(node_count);
  std::vector<graph> input(static_cast<std::size_t>(words) * node_count);
  std::vector<graph> canonical(static_cast<std::size_t>(words) * node_count);
  EMPTYGRAPH(input.data(), words, node_count);
  for (int edge = 0; edge < static_cast<int>(kEdges.size()); ++edge) {
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (kEdges[edge] & (Mask{1} << vertex)) {
        ADDONEEDGE(input.data(), edge_start + edge, vertex, words);
      }
    }
  }
  for (int four = 0; four < static_cast<int>(g_fours.size()); ++four) {
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (g_fours[four] & (Mask{1} << vertex)) {
        ADDONEEDGE(input.data(), four_start + four, vertex, words);
      }
    }
  }
  for (int row = 0; row < row_count; ++row) {
    const Mask face = g_candidates[chosen[row]].face;
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (face & (Mask{1} << vertex)) {
        ADDONEEDGE(input.data(), row_start + row, vertex, words);
      }
    }
  }

  std::vector<int> lab(node_count), ptn(node_count), orbits(node_count);
  int position = 0;
  auto add_cell = [&](int first, int last) {
    if (first == last) return;
    for (int value = first; value < last; ++value) {
      lab[position] = value;
      ptn[position++] = 1;
    }
    ptn[position - 1] = 0;
  };
  add_cell(0, kN);
  add_cell(edge_start, four_start);
  add_cell(four_start, row_start);
  add_cell(row_start, node_count);
  if (position != node_count) {
    throw std::runtime_error("nauty color partition changed");
  }
  DEFAULTOPTIONS_GRAPH(options);
  statsblk stats{};
  options.defaultptn = FALSE;
  options.getcanon = TRUE;
  densenauty(input.data(), lab.data(), ptn.data(), orbits.data(), &options,
             &stats, words, node_count, canonical.data());
  if (stats.errstatus != 0) {
    throw std::runtime_error("nauty canonicalization failed");
  }
  if (chosen.empty()) {
    const double order = stats.grpsize1 * std::pow(10.0, stats.grpsize2);
    g_stabilizer_order = static_cast<std::uint64_t>(std::llround(order));
  }
  return std::string(
      reinterpret_cast<const char*>(canonical.data()),
      canonical.size() * sizeof(graph));
}

std::vector<Mask> LiteralMinimalBlocker(const std::vector<Mask>& clutter) {
  std::vector<Mask> answer;
  for (Mask candidate = 1; candidate <= kFull; ++candidate) {
    if (!HitsAll(candidate, clutter)) continue;
    bool minimal = true;
    for (int vertex = 0; vertex < kN; ++vertex) {
      if ((candidate >> vertex) & 1u &&
          HitsAll(candidate ^ (Mask{1} << vertex), clutter)) {
        minimal = false;
      }
    }
    if (minimal) answer.push_back(candidate);
  }
  return answer;
}

void Finish(const State& state) {
  ++g_scheduled_orbits;
  std::vector<Mask> clutter(kEdges.begin(), kEdges.end());
  clutter.insert(clutter.end(), g_fours.begin(), g_fours.end());
  for (int index : state.chosen) clutter.push_back(g_candidates[index].face);
  const auto blocker = LiteralMinimalBlocker(clutter);
  ++g_blocker_histogram[blocker.size()];
  bool pure = blocker.size() == kScheduledFacets;
  for (Mask transversal : blocker) pure = pure && std::popcount(transversal) == 5;
  if (!pure) return;
  ++g_pure_target_orbits;
  std::vector<Mask> representative;
  for (int index : state.chosen) representative.push_back(g_candidates[index].face);
  std::sort(representative.begin(), representative.end());
  g_target_representatives.push_back(std::move(representative));
  std::cerr << "pure target orbit " << g_pure_target_orbits
            << " after " << g_nodes << " nodes\n";
}

void Search(const State& state) {
  ++g_nodes;
  if (g_nodes % 1000000 == 0) {
    std::cerr << "nodes=" << g_nodes << " depth=" << state.chosen.size()
              << " seen=";
    for (const auto& level : g_seen) std::cerr << level.size() << ',';
    std::cerr << " targets=" << g_pure_target_orbits << '\n';
  }
  const int depth = static_cast<int>(state.chosen.size());
  if (!g_seen[depth].insert(CanonicalCode(state.chosen)).second) return;
  const int surviving_count = BitCount(state.surviving_covers);
  if (surviving_count < kScheduledFacets) {
    ++g_facet_prunes;
    return;
  }
  const std::uint64_t full =
      g_small.empty() ? 0 : (std::uint64_t{1} << g_small.size()) - 1;
  if (depth == g_additions) {
    if (state.destroyed_small == full &&
        surviving_count == kScheduledFacets) {
      Finish(state);
    }
    return;
  }

  std::vector<int> best;
  std::uint64_t feasible_union = state.destroyed_small;
  bool has_uncovered = false;
  for (int requirement = 0;
       requirement < static_cast<int>(g_small.size()); ++requirement) {
    if ((state.destroyed_small >> requirement) & 1u) continue;
    has_uncovered = true;
    std::vector<int> feasible;
    for (int candidate : g_destroyers[requirement]) {
      if (std::find(state.chosen.begin(), state.chosen.end(), candidate) !=
              state.chosen.end() ||
          !CapCompatible(state, g_candidates[candidate].face)) {
        continue;
      }
      const auto surviving = Meet(
          state.surviving_covers, g_candidates[candidate].hits_covers);
      if (BitCount(surviving) < kScheduledFacets) continue;
      feasible.push_back(candidate);
      feasible_union |= g_candidates[candidate].destroys_small;
    }
    if (feasible.empty()) {
      ++g_uncoverable_prunes;
      return;
    }
    if (best.empty() || feasible.size() < best.size()) best = std::move(feasible);
  }
  if ((feasible_union & full) != full) {
    ++g_uncoverable_prunes;
    return;
  }
  if (!has_uncovered) {
    for (int candidate = 0;
         candidate < static_cast<int>(g_candidates.size()); ++candidate) {
      if (std::find(state.chosen.begin(), state.chosen.end(), candidate) !=
              state.chosen.end() ||
          !CapCompatible(state, g_candidates[candidate].face)) {
        continue;
      }
      const auto surviving = Meet(
          state.surviving_covers, g_candidates[candidate].hits_covers);
      if (BitCount(surviving) >= kScheduledFacets) best.push_back(candidate);
    }
  }
  for (int candidate : best) {
    State next = state;
    next.chosen.push_back(candidate);
    next.destroyed_small |= g_candidates[candidate].destroys_small;
    for (int vertex = 0; vertex < kN; ++vertex) {
      next.incidence[vertex] +=
          (g_candidates[candidate].face >> vertex) & 1u;
    }
    next.surviving_covers = Meet(
        state.surviving_covers, g_candidates[candidate].hits_covers);
    Search(next);
  }
}

void Initialize() {
  std::vector<Mask> independent_fives;
  for (Mask face = 0; face <= kFull; ++face) {
    if (std::popcount(face) == 5 && Independent(face)) {
      independent_fives.push_back(face);
    }
  }
  if (independent_fives.size() != 1372) {
    throw std::runtime_error("independent-five ground changed");
  }
  for (Mask face : independent_fives) {
    if (ContainsFour(face)) g_shadow.push_back(face);
  }
  g_additions = kScheduledFives - static_cast<int>(g_shadow.size());
  if (g_additions < 0 || g_additions > kMaximumAdditions) {
    throw std::runtime_error("a=2 addition schedule escaped the declared bound");
  }
  g_seen.resize(g_additions + 1);

  for (int vertex = 0; vertex < kN; ++vertex) {
    int four_incidence = 0;
    int shadow_incidence = 0;
    for (Mask four : g_fours) four_incidence += (four >> vertex) & 1u;
    for (Mask face : g_shadow) shadow_incidence += (face >> vertex) & 1u;
    const int previous = kBaseG3[vertex] - four_incidence;
    const int before =
        kBaseG4[vertex] + 6 * four_incidence - shadow_incidence;
    const int lower = std::max(0, before - MacaulayRaise(previous, 3));
    if (lower != 0) throw std::runtime_error("a=2 lower cap changed");
    g_upper[vertex] = before - 1;
  }

  std::vector<Mask> base(kEdges.begin(), kEdges.end());
  base.insert(base.end(), g_fours.begin(), g_fours.end());
  base.insert(base.end(), g_shadow.begin(), g_shadow.end());
  for (int order = 1; order <= 4; ++order) {
    for (Mask transversal = 1; transversal <= kFull; ++transversal) {
      if (std::popcount(transversal) == order && HitsAll(transversal, base)) {
        g_small.push_back(transversal);
      }
    }
  }
  if (g_small.size() > 63) {
    throw std::runtime_error("a=2 small-transversal gate escaped one word");
  }

  for (Mask cover = 0; cover <= kFull; ++cover) {
    if (std::popcount(cover) != 5) continue;
    bool hits = true;
    for (Mask edge : kEdges) hits = hits && (cover & edge);
    for (Mask four : g_fours) hits = hits && (cover & four);
    if (hits) g_covers.push_back(cover);
  }
  std::vector<Mask> candidate_faces;
  for (Mask face : independent_fives) {
    if (ContainsFour(face)) continue;
    bool allowed = true;
    for (int vertex = 0; vertex < kN; ++vertex) {
      if ((face >> vertex & 1u) && g_upper[vertex] < 1) allowed = false;
    }
    if (allowed) candidate_faces.push_back(face);
  }
  const int cover_words = (g_covers.size() + 63) / 64;
  for (Mask face : candidate_faces) {
    Candidate candidate;
    candidate.face = face;
    for (int requirement = 0;
         requirement < static_cast<int>(g_small.size()); ++requirement) {
      if (!(face & g_small[requirement])) {
        candidate.destroys_small |= std::uint64_t{1} << requirement;
      }
    }
    candidate.hits_covers.assign(cover_words, 0);
    for (int cover = 0; cover < static_cast<int>(g_covers.size()); ++cover) {
      if (face & g_covers[cover]) {
        candidate.hits_covers[cover / 64] |=
            std::uint64_t{1} << (cover % 64);
      }
    }
    g_candidates.push_back(std::move(candidate));
  }
  g_destroyers.assign(g_small.size(), {});
  for (int candidate = 0;
       candidate < static_cast<int>(g_candidates.size()); ++candidate) {
    for (int requirement = 0;
         requirement < static_cast<int>(g_small.size()); ++requirement) {
      if (g_candidates[candidate].destroys_small >> requirement & 1u) {
        g_destroyers[requirement].push_back(candidate);
      }
    }
  }
}

void PrintFamily(const std::vector<Mask>& family) {
  std::cout << '[';
  for (int row = 0; row < static_cast<int>(family.size()); ++row) {
    if (row) std::cout << ',';
    std::cout << family[row];
  }
  std::cout << ']';
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 1 && argc != 4) {
      throw std::runtime_error(
          "usage: scout [positive_A4_orbit_index first_A4_mask second_A4_mask]");
    }
    if (argc == 4) {
      g_orbit_index = std::stoi(argv[1]);
      if (g_orbit_index < 0 || g_orbit_index >= 21) {
        throw std::runtime_error("positive-A4 orbit index is out of range");
      }
      for (int row = 0; row < 2; ++row) {
        const unsigned long raw = std::stoul(argv[row + 2]);
        if (raw > kFull || std::popcount(static_cast<unsigned>(raw)) != 4) {
          throw std::runtime_error("each A4 mask must be a 14-bit four-set");
        }
        g_fours[row] = static_cast<Mask>(raw);
        if (!Independent(g_fours[row])) {
          throw std::runtime_error("each A4 support must avoid the missing graph");
        }
      }
      if (g_fours[0] == g_fours[1]) {
        throw std::runtime_error("the two A4 supports must be distinct");
      }
      std::sort(g_fours.begin(), g_fours.end());
    }
    Initialize();
    State initial;
    initial.surviving_covers.assign((g_covers.size() + 63) / 64, 0);
    for (int cover = 0; cover < static_cast<int>(g_covers.size()); ++cover) {
      initial.surviving_covers[cover / 64] |=
          std::uint64_t{1} << (cover % 64);
    }
    Search(initial);
    if (!g_stabilizer_order) {
      throw std::runtime_error("fixed-structure stabilizer was not computed");
    }
    std::cout << "{\"kind\":\"h9_14_record30_a2_nauty_quotient_scout\","
              << "\"compact_only\":true,\"positive_A4_orbit_index\":"
              << g_orbit_index << ",\"fixed_A4_masks\":[" << g_fours[0]
              << ',' << g_fours[1] << "],\"stabilizer_order\":"
              << g_stabilizer_order << ','
              << "\"forced_shadow_count\":" << g_shadow.size()
              << ",\"additional_A5_count\":" << g_additions
              << ",\"candidate_count\":"
              << g_candidates.size()
              << ",\"raw_unordered_addition_count\":"
              << Choose(g_candidates.size(), g_additions)
              << ",\"base_small_transversal_count\":" << g_small.size()
              << ",\"initial_order_five_cover_count\":" << g_covers.size()
              << ",\"search_node_count\":" << g_nodes
              << ",\"canonical_state_counts\":[";
    for (int depth = 0; depth <= g_additions; ++depth) {
      if (depth) std::cout << ',';
      std::cout << g_seen[depth].size();
    }
    std::cout << "],\"facet_lower_bound_prune_count\":" << g_facet_prunes
              << ",\"uncoverable_prune_count\":" << g_uncoverable_prunes
              << ",\"scheduled_incidence_orbit_count\":" << g_scheduled_orbits
              << ",\"scheduled_blocker_count_histogram\":{";
    bool first = true;
    for (const auto& [size, count] : g_blocker_histogram) {
      if (!first) std::cout << ',';
      std::cout << '\"' << size << "\":" << count;
      first = false;
    }
    std::cout << "},\"pure_128_facet_target_orbit_count\":"
              << g_pure_target_orbits
              << ",\"pure_target_representatives\":[";
    for (int index = 0;
         index < static_cast<int>(g_target_representatives.size()); ++index) {
      if (index) std::cout << ',';
      PrintFamily(g_target_representatives[index]);
    }
    std::cout << "],\"scope\":\"Exact nauty-quotiented compact incidence "
                 "scout for one Record-30 a=2 positive-A4 orbit; independent replay and "
                 "Coxeter labels remain separate. No finite-volume inference.\"}\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
