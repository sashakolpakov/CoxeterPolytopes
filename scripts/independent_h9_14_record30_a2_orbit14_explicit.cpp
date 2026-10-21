// Standalone explicit-group replay for Record 30, compact a=2 orbit 14.
//
// This program intentionally does not use nauty, Z3, Python project imports,
// or a frozen target list.  It reconstructs the graph automorphisms from the
// three disjoint missing edges and eight isolates, filters them by the
// unordered pair of missing A4 supports, and obtains the exact stabilizer.
// It then exhausts the five additional minimal-A5 families using an explicit
// canonical form for C2 x C2 x S6, replays every blocker literally, and
// expands the final quotient representatives through all stabilizer elements.
//
// Coxeter labels are deliberately left to a separately implemented replay.
// There is no finite-volume noncompact inference.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

constexpr int kN = 14;
constexpr int kAdditionCount = 5;
constexpr int kScheduledFacetCount = 128;
using Mask = std::uint16_t;
constexpr Mask kFull = (Mask{1} << kN) - 1;
constexpr std::array<Mask, 3> kEdges = {
    (Mask{1} << 0) | (Mask{1} << 11),
    (Mask{1} << 1) | (Mask{1} << 12),
    (Mask{1} << 2) | (Mask{1} << 13),
};
constexpr std::array<Mask, 2> kFours = {
    (Mask{1} << 0) | (Mask{1} << 1) | (Mask{1} << 3) | (Mask{1} << 4),
    (Mask{1} << 0) | (Mask{1} << 2) | (Mask{1} << 3) | (Mask{1} << 4),
};
constexpr std::array<int, kN> kBaseG3 = {
    4, 4, 4, 8, 8, 8, 8, 8, 8, 8, 8, 4, 4, 4,
};
constexpr std::array<int, kN> kBaseG4 = kBaseG3;

using Permutation = std::array<std::uint8_t, kN>;

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
std::array<std::unordered_set<std::string>, kAdditionCount + 1> g_seen;
std::vector<Permutation> g_stabilizer;
std::vector<std::vector<Mask>> g_targets;
std::uint64_t g_nodes = 0;
std::uint64_t g_facet_prunes = 0;
std::uint64_t g_uncoverable_prunes = 0;
std::uint64_t g_scheduled_orbits = 0;
std::map<int, std::uint64_t> g_blocker_histogram;

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
  for (Mask four : kFours) {
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

Mask PermutedMask(Mask value, const Permutation& permutation) {
  Mask answer = 0;
  for (int vertex = 0; vertex < kN; ++vertex) {
    if ((value >> vertex) & 1u) answer |= Mask{1} << permutation[vertex];
  }
  return answer;
}

std::string PermutationCode(const Permutation& permutation) {
  return std::string(
      reinterpret_cast<const char*>(permutation.data()), permutation.size());
}

bool SameMaskSet(const std::array<Mask, 2>& left,
                 const std::array<Mask, 2>& right) {
  auto first = left;
  auto second = right;
  std::sort(first.begin(), first.end());
  std::sort(second.begin(), second.end());
  return first == second;
}

std::vector<Permutation> EnumerateStabilizer() {
  // Every automorphism of 3K2 + 8K1 is a permutation of the three edge
  // components, an independent flip of each component, and a permutation of
  // the eight isolates.  Filtering this complete 48*8! list by the unordered
  // A4-support pair is therefore an exact stabilizer computation.
  constexpr std::array<std::array<int, 2>, 3> endpoints = {
      std::array<int, 2>{0, 11},
      std::array<int, 2>{1, 12},
      std::array<int, 2>{2, 13},
  };
  std::array<int, 3> component_image = {0, 1, 2};
  std::array<int, 8> isolate_image = {3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<Permutation> answer;
  std::unordered_set<std::string> distinct;
  do {
    for (int flip_bits = 0; flip_bits < 8; ++flip_bits) {
      do {
        Permutation permutation{};
        for (int component = 0; component < 3; ++component) {
          const int target = component_image[component];
          const int flip = (flip_bits >> component) & 1;
          permutation[endpoints[component][0]] = endpoints[target][flip];
          permutation[endpoints[component][1]] = endpoints[target][1 - flip];
        }
        for (int index = 0; index < 8; ++index) {
          permutation[3 + index] = isolate_image[index];
        }
        const std::array<Mask, 2> image = {
            PermutedMask(kFours[0], permutation),
            PermutedMask(kFours[1], permutation),
        };
        if (SameMaskSet(image, kFours)) {
          if (!distinct.insert(PermutationCode(permutation)).second) {
            throw std::runtime_error("duplicate graph automorphism");
          }
          answer.push_back(permutation);
        }
      } while (std::next_permutation(isolate_image.begin(),
                                     isolate_image.end()));
    }
  } while (std::next_permutation(component_image.begin(),
                                 component_image.end()));
  if (answer.size() != 2880) {
    throw std::runtime_error("explicit orbit-14 stabilizer order changed");
  }

  // Independently check the structural description used by CanonicalCode:
  // (1 2)(12 13), (3 4), and an arbitrary permutation of 5,...,10.
  std::unordered_set<std::string> structural;
  std::array<int, 6> tail = {5, 6, 7, 8, 9, 10};
  for (int swap_12 = 0; swap_12 < 2; ++swap_12) {
    for (int swap_34 = 0; swap_34 < 2; ++swap_34) {
      do {
        Permutation permutation{};
        std::iota(permutation.begin(), permutation.end(), 0);
        if (swap_12) {
          std::swap(permutation[1], permutation[2]);
          std::swap(permutation[12], permutation[13]);
        }
        if (swap_34) std::swap(permutation[3], permutation[4]);
        for (int index = 0; index < 6; ++index) {
          permutation[5 + index] = tail[index];
        }
        structural.insert(PermutationCode(permutation));
      } while (std::next_permutation(tail.begin(), tail.end()));
    }
  }
  if (structural != distinct) {
    throw std::runtime_error("structural and exhaustive stabilizers disagree");
  }
  return answer;
}

std::string CanonicalCode(const std::vector<int>& chosen) {
  // Rows are unordered.  The two indicated involutions are enumerated;
  // permutations of vertices 5,...,10 are eliminated by sorting their column
  // signatures.  Enumerating all row orders makes this a complete invariant,
  // not a hash, for C2 x C2 x S6 orbits of partial row families.
  const int depth = static_cast<int>(chosen.size());
  std::vector<int> row_order(depth);
  std::iota(row_order.begin(), row_order.end(), 0);
  std::string best;
  bool initialized = false;
  for (int swap_12 = 0; swap_12 < 2; ++swap_12) {
    for (int swap_34 = 0; swap_34 < 2; ++swap_34) {
      do {
        std::array<std::uint8_t, kN> columns{};
        for (int position = 0; position < depth; ++position) {
          Mask face = g_candidates[chosen[row_order[position]]].face;
          if (swap_12) {
            const bool one = (face >> 1) & 1u;
            const bool two = (face >> 2) & 1u;
            const bool twelve = (face >> 12) & 1u;
            const bool thirteen = (face >> 13) & 1u;
            face &= ~((Mask{1} << 1) | (Mask{1} << 2) |
                      (Mask{1} << 12) | (Mask{1} << 13));
            if (one) face |= Mask{1} << 2;
            if (two) face |= Mask{1} << 1;
            if (twelve) face |= Mask{1} << 13;
            if (thirteen) face |= Mask{1} << 12;
          }
          if (swap_34) {
            const bool three = (face >> 3) & 1u;
            const bool four = (face >> 4) & 1u;
            face &= ~((Mask{1} << 3) | (Mask{1} << 4));
            if (three) face |= Mask{1} << 4;
            if (four) face |= Mask{1} << 3;
          }
          for (int vertex = 0; vertex < kN; ++vertex) {
            if ((face >> vertex) & 1u) {
              columns[vertex] |= std::uint8_t{1} << position;
            }
          }
        }
        std::array<std::uint8_t, 6> tail = {
            columns[5], columns[6], columns[7],
            columns[8], columns[9], columns[10],
        };
        std::sort(tail.begin(), tail.end());
        const std::array<int, 8> fixed = {0, 1, 2, 3, 4, 11, 12, 13};
        std::string code;
        code.reserve(kN);
        for (int vertex : fixed) code.push_back(columns[vertex]);
        for (std::uint8_t signature : tail) code.push_back(signature);
        if (!initialized || code < best) {
          best = std::move(code);
          initialized = true;
        }
      } while (std::next_permutation(row_order.begin(), row_order.end()));
    }
  }
  if (!initialized || best.size() != kN) {
    throw std::runtime_error("explicit canonical form failed");
  }
  return best;
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

std::vector<Mask> LiteralMinimalBlocker(const std::vector<Mask>& clutter) {
  std::vector<Mask> answer;
  for (Mask candidate = 1; candidate <= kFull; ++candidate) {
    if (!HitsAll(candidate, clutter)) continue;
    bool minimal = true;
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (((candidate >> vertex) & 1u) &&
          HitsAll(candidate ^ (Mask{1} << vertex), clutter)) {
        minimal = false;
        break;
      }
    }
    if (minimal) answer.push_back(candidate);
  }
  return answer;
}

void Finish(const State& state) {
  ++g_scheduled_orbits;
  std::vector<Mask> clutter(kEdges.begin(), kEdges.end());
  clutter.insert(clutter.end(), kFours.begin(), kFours.end());
  std::vector<Mask> family;
  for (int index : state.chosen) {
    clutter.push_back(g_candidates[index].face);
    family.push_back(g_candidates[index].face);
  }
  const auto blocker = LiteralMinimalBlocker(clutter);
  ++g_blocker_histogram[blocker.size()];
  bool pure = blocker.size() == kScheduledFacetCount;
  for (Mask transversal : blocker) {
    pure = pure && std::popcount(transversal) == 5;
  }
  if (!pure) return;
  if (BitCount(state.surviving_covers) != static_cast<int>(blocker.size())) {
    throw std::runtime_error("cover intersection and blocker disagree");
  }
  std::sort(family.begin(), family.end());
  g_targets.push_back(std::move(family));
}

void Search(const State& state) {
  ++g_nodes;
  const int depth = static_cast<int>(state.chosen.size());
  if (!g_seen[depth].insert(CanonicalCode(state.chosen)).second) return;
  const int surviving_count = BitCount(state.surviving_covers);
  if (surviving_count < kScheduledFacetCount) {
    ++g_facet_prunes;
    return;
  }
  const std::uint64_t full =
      (std::uint64_t{1} << static_cast<int>(g_small.size())) - 1;
  if (depth == kAdditionCount) {
    if (state.destroyed_small == full &&
        surviving_count == kScheduledFacetCount) {
      Finish(state);
    }
    return;
  }

  std::vector<int> branch;
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
      if (BitCount(surviving) < kScheduledFacetCount) continue;
      feasible.push_back(candidate);
      feasible_union |= g_candidates[candidate].destroys_small;
    }
    if (feasible.empty()) {
      ++g_uncoverable_prunes;
      return;
    }
    if (branch.empty() || feasible.size() < branch.size()) {
      branch = std::move(feasible);
    }
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
      if (BitCount(surviving) >= kScheduledFacetCount) {
        branch.push_back(candidate);
      }
    }
  }

  for (int candidate : branch) {
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
  g_stabilizer = EnumerateStabilizer();
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
  if (g_shadow.size() != 15) {
    throw std::runtime_error("forced A5 shadow changed");
  }

  for (int vertex = 0; vertex < kN; ++vertex) {
    int four_incidence = 0;
    int shadow_incidence = 0;
    for (Mask four : kFours) four_incidence += (four >> vertex) & 1u;
    for (Mask face : g_shadow) shadow_incidence += (face >> vertex) & 1u;
    const int previous = kBaseG3[vertex] - four_incidence;
    const int before =
        kBaseG4[vertex] + 6 * four_incidence - shadow_incidence;
    const int lower = std::max(0, before - MacaulayRaise(previous, 3));
    if (lower != 0) throw std::runtime_error("local lower cap changed");
    g_upper[vertex] = before - 1;
  }

  std::vector<Mask> base(kEdges.begin(), kEdges.end());
  base.insert(base.end(), kFours.begin(), kFours.end());
  base.insert(base.end(), g_shadow.begin(), g_shadow.end());
  for (Mask transversal = 1; transversal <= kFull; ++transversal) {
    if (std::popcount(transversal) <= 4 && HitsAll(transversal, base)) {
      g_small.push_back(transversal);
    }
  }
  if (g_small.size() != 61) {
    throw std::runtime_error("small-transversal ground changed");
  }

  for (Mask cover = 0; cover <= kFull; ++cover) {
    if (std::popcount(cover) != 5) continue;
    bool hits = true;
    for (Mask edge : kEdges) hits = hits && (cover & edge);
    for (Mask four : kFours) hits = hits && (cover & four);
    if (hits) g_covers.push_back(cover);
  }
  if (g_covers.size() != 269) {
    throw std::runtime_error("initial order-five cover count changed");
  }

  for (Mask face : independent_fives) {
    if (ContainsFour(face)) continue;
    bool allowed = true;
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (((face >> vertex) & 1u) && g_upper[vertex] < 1) allowed = false;
    }
    if (!allowed) continue;
    Candidate candidate;
    candidate.face = face;
    for (int requirement = 0;
         requirement < static_cast<int>(g_small.size()); ++requirement) {
      if (!(face & g_small[requirement])) {
        candidate.destroys_small |= std::uint64_t{1} << requirement;
      }
    }
    candidate.hits_covers.assign((g_covers.size() + 63) / 64, 0);
    for (int cover = 0; cover < static_cast<int>(g_covers.size()); ++cover) {
      if (face & g_covers[cover]) {
        candidate.hits_covers[cover / 64] |=
            std::uint64_t{1} << (cover % 64);
      }
    }
    g_candidates.push_back(std::move(candidate));
  }
  if (g_candidates.size() != 966) {
    throw std::runtime_error("candidate ground changed");
  }
  g_destroyers.assign(g_small.size(), {});
  for (int candidate = 0;
       candidate < static_cast<int>(g_candidates.size()); ++candidate) {
    for (int requirement = 0;
         requirement < static_cast<int>(g_small.size()); ++requirement) {
      if ((g_candidates[candidate].destroys_small >> requirement) & 1u) {
        g_destroyers[requirement].push_back(candidate);
      }
    }
  }
}

void PrintFamily(const std::vector<Mask>& family) {
  std::cout << '[';
  for (int index = 0; index < static_cast<int>(family.size()); ++index) {
    if (index) std::cout << ',';
    std::cout << family[index];
  }
  std::cout << ']';
}

}  // namespace

int main() {
  try {
    Initialize();
    State initial;
    initial.surviving_covers.assign((g_covers.size() + 63) / 64, 0);
    for (int cover = 0; cover < static_cast<int>(g_covers.size()); ++cover) {
      initial.surviving_covers[cover / 64] |=
          std::uint64_t{1} << (cover % 64);
    }
    Search(initial);
    std::sort(g_targets.begin(), g_targets.end());
    if (g_targets.size() != g_scheduled_orbits) {
      throw std::runtime_error("a scheduled target lost blocker purity");
    }

    std::set<std::vector<Mask>> literal_union;
    std::map<int, int> orbit_histogram;
    std::vector<int> orbit_sizes;
    for (const auto& family : g_targets) {
      std::set<std::vector<Mask>> orbit;
      for (const auto& permutation : g_stabilizer) {
        std::vector<Mask> image;
        for (Mask face : family) {
          image.push_back(PermutedMask(face, permutation));
        }
        std::sort(image.begin(), image.end());
        orbit.insert(image);
      }
      orbit_sizes.push_back(static_cast<int>(orbit.size()));
      ++orbit_histogram[static_cast<int>(orbit.size())];
      for (const auto& image : orbit) {
        if (!literal_union.insert(image).second) {
          throw std::runtime_error("distinct quotient targets have overlapping orbits");
        }
      }
    }

    std::cout
        << "{\"kind\":\"h9_14_record30_a2_orbit14_explicit_group_audit\""
        << ",\"compact_only\":true"
        << ",\"graph_automorphism_candidates_checked\":" << 48 * 40320
        << ",\"stabilizer_order\":" << g_stabilizer.size()
        << ",\"stabilizer_structure\":\"C2 x C2 x S6\""
        << ",\"forced_shadow_count\":" << g_shadow.size()
        << ",\"additional_A5_count\":" << kAdditionCount
        << ",\"candidate_count\":" << g_candidates.size()
        << ",\"raw_unordered_addition_count\":"
        << Choose(static_cast<int>(g_candidates.size()), kAdditionCount)
        << ",\"base_small_transversal_count\":" << g_small.size()
        << ",\"initial_order_five_cover_count\":" << g_covers.size()
        << ",\"search_node_count\":" << g_nodes
        << ",\"canonical_state_counts\":[";
    for (int depth = 0; depth <= kAdditionCount; ++depth) {
      if (depth) std::cout << ',';
      std::cout << g_seen[depth].size();
    }
    std::cout << "]"
              << ",\"facet_lower_bound_prune_count\":" << g_facet_prunes
              << ",\"uncoverable_prune_count\":" << g_uncoverable_prunes
              << ",\"scheduled_incidence_orbit_count\":" << g_scheduled_orbits
              << ",\"scheduled_blocker_count_histogram\":{";
    bool first = true;
    for (const auto& [size, count] : g_blocker_histogram) {
      if (!first) std::cout << ',';
      std::cout << '\"' << size << "\":" << count;
      first = false;
    }
    std::cout << "}"
              << ",\"pure_128_facet_target_orbit_count\":" << g_targets.size()
              << ",\"literal_target_count\":" << literal_union.size()
              << ",\"orbit_size_histogram\":{";
    first = true;
    for (const auto& [size, count] : orbit_histogram) {
      if (!first) std::cout << ',';
      std::cout << '\"' << size << "\":" << count;
      first = false;
    }
    std::cout << "},\"target_orbit_sizes\":[";
    for (int index = 0; index < static_cast<int>(orbit_sizes.size()); ++index) {
      if (index) std::cout << ',';
      std::cout << orbit_sizes[index];
    }
    std::cout << "],\"pure_target_representatives\":[";
    for (int index = 0; index < static_cast<int>(g_targets.size()); ++index) {
      if (index) std::cout << ',';
      PrintFamily(g_targets[index]);
    }
    std::cout
        << "],\"scope\":\"Standalone exact explicit-group compact incidence replay for Record-30 a=2 positive-A4 orbit 14; no nauty, Z3, project import, Coxeter-label, or finite-volume inference.\"}\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
