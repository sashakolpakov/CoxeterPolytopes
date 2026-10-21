// Standalone exhaustive compact-incidence replay for H9/14 Record 30, a=1.
//
// This program deliberately does not use nauty, a SAT/SMT solver, or the
// production canonicalizer.  A partial family is represented by the colored
// incidence matrix whose row colors are
//
//   * the three unordered missing-edge rows,
//   * the unique fixed missing-A4 row, and
//   * the unordered added missing-A5 rows.
//
// Its canonical form enumerates the 3! edge-row orders and the admissible
// A5-row orders, then sorts the fourteen column signatures.  A5 rows having
// different one-row canonical forms cannot be exchanged by an isomorphism,
// so only permutations inside equal one-row classes are enumerated.  This is
// an invariant speedup of the full d! enumeration, not a symmetry assumption.
// Equality of the resulting arrays is therefore exactly colored-incidence
// isomorphism.  The search separately replays all literal transversal, cover,
// local-link, blocker, face-table, ridge, g/Macaulay, and Euler-link gates.
//
// The result is compact-only.  Incidence survival does not assert PL
// sphericity, convex polytopality, Coxeter realizability, or finite volume.

#include <algorithm>
#include <array>
#include <bit>
#include <bitset>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

constexpr int kN = 14;
constexpr int kScheduledMissingFives = 14;
constexpr int kScheduledFacets = 130;
constexpr int kMaximumSmall = 128;
using Mask = std::uint16_t;
using CoverBits = std::vector<std::uint64_t>;
using SmallBits = std::bitset<kMaximumSmall>;
using Permutation = std::array<std::uint8_t, kN>;
using CanonicalCode = std::array<std::uint16_t, kN>;
constexpr Mask kFull = (Mask{1} << kN) - 1;
constexpr std::array<Mask, 3> kEdges = {
    (Mask{1} << 0) | (Mask{1} << 11),
    (Mask{1} << 1) | (Mask{1} << 12),
    (Mask{1} << 2) | (Mask{1} << 13),
};
constexpr std::array<int, 6> kActive = {0, 1, 2, 11, 12, 13};
constexpr std::array<int, 8> kIsolates = {3, 4, 5, 6, 7, 8, 9, 10};
constexpr std::array<Mask, 3> kOrbitFours = {
    (Mask{1} << 3) | (Mask{1} << 4) | (Mask{1} << 5) | (Mask{1} << 6),
    (Mask{1} << 0) | (Mask{1} << 3) | (Mask{1} << 4) | (Mask{1} << 5),
    (Mask{1} << 0) | (Mask{1} << 1) | (Mask{1} << 3) | (Mask{1} << 4),
};
constexpr std::array<int, 3> kExpectedShadows = {10, 9, 8};
constexpr std::array<int, 3> kExpectedAdditions = {4, 5, 6};
constexpr std::array<std::uint64_t, 3> kExpectedStabilizers = {
    27648, 5760, 5760,
};
constexpr std::array<int, 9> kExpectedFVector = {
    14, 88, 328, 805, 1358, 1582, 1232, 585, 130,
};

struct CodeHash {
  std::size_t operator()(const CanonicalCode& code) const noexcept {
    std::uint64_t value = 1469598103934665603ULL;
    for (std::uint16_t entry : code) {
      value ^= entry & 0xffu;
      value *= 1099511628211ULL;
      value ^= entry >> 8;
      value *= 1099511628211ULL;
    }
    return static_cast<std::size_t>(value);
  }
};

struct LocalConstraint {
  Mask subset = 0;
  int lower = 0;
  int upper = 0;
  int before = 0;
  int previous = 0;
};

struct Candidate {
  Mask face = 0;
  SmallBits destroys_small;
  CoverBits hits_covers;
  std::array<int, 25> local_indices{};
  CanonicalCode singleton_code{};
};

struct State {
  std::vector<int> chosen;
  SmallBits destroyed_small;
  std::array<std::uint8_t, 430> local_incidence{};
  CoverBits surviving_covers;
};

struct LinkWitness {
  bool present = false;
  Mask subset = 0;
  int dimension = 0;
  std::vector<int> f_prefix;
  std::vector<int> h;
  std::vector<int> g;
};

int g_orbit_index = -1;
Mask g_four = 0;
int g_additions = 0;
std::vector<Mask> g_independent_fives;
std::vector<Mask> g_shadow;
std::vector<Mask> g_small;
std::vector<Mask> g_covers;
std::vector<LocalConstraint> g_local_constraints;
std::array<int, 1 << kN> g_local_index;
std::vector<Candidate> g_candidates;
std::vector<std::vector<int>> g_destroyers;
std::vector<std::unordered_set<CanonicalCode, CodeHash>> g_seen;
LinkWitness g_fixed_link_witness;

std::uint64_t g_nodes = 0;
std::uint64_t g_isomorphic_prunes = 0;
std::uint64_t g_facet_prunes = 0;
std::uint64_t g_uncoverable_prunes = 0;
std::uint64_t g_scheduled_orbits = 0;
std::uint64_t g_pure_blocker_orbits = 0;
std::uint64_t g_face_table_rejects = 0;
std::uint64_t g_ridge_rejects = 0;
std::uint64_t g_link_rejects = 0;
std::uint64_t g_full_canonical_crosschecks = 0;
std::map<int, std::uint64_t> g_blocker_histogram;
std::vector<std::vector<Mask>> g_targets;

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
  if (value < 0 || degree < 1) throw std::runtime_error("bad Macaulay input");
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

bool HitsAll(Mask transversal, const std::vector<Mask>& clutter) {
  for (Mask face : clutter) {
    if (!(transversal & face)) return false;
  }
  return true;
}

Mask Transform(Mask value, const Permutation& permutation) {
  Mask result = 0;
  for (int vertex = 0; vertex < kN; ++vertex) {
    if (value & (Mask{1} << vertex)) {
      result |= Mask{1} << permutation[vertex];
    }
  }
  return result;
}

int BitCount(const CoverBits& words) {
  int answer = 0;
  for (std::uint64_t word : words) answer += std::popcount(word);
  return answer;
}

CoverBits Meet(const CoverBits& left, const CoverBits& right) {
  if (left.size() != right.size()) throw std::runtime_error("cover width changed");
  CoverBits answer(left.size());
  for (int index = 0; index < static_cast<int>(left.size()); ++index) {
    answer[index] = left[index] & right[index];
  }
  return answer;
}

CanonicalCode OrderedIncidenceCode(
    const std::vector<Mask>& rows, const std::vector<int>& row_order,
    const std::array<int, 3>& edge_order) {
  std::array<std::uint16_t, kN> signatures{};
  for (int position = 0; position < 3; ++position) {
    const Mask edge = kEdges[edge_order[position]];
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (edge & (Mask{1} << vertex)) {
        signatures[vertex] |= std::uint16_t{1} << position;
      }
    }
  }
  for (int vertex = 0; vertex < kN; ++vertex) {
    if (g_four & (Mask{1} << vertex)) {
      signatures[vertex] |= std::uint16_t{1} << 3;
    }
  }
  for (int position = 0; position < static_cast<int>(row_order.size());
       ++position) {
    const Mask face = rows[row_order[position]];
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (face & (Mask{1} << vertex)) {
        signatures[vertex] |= std::uint16_t{1} << (4 + position);
      }
    }
  }
  std::sort(signatures.begin(), signatures.end());
  return signatures;
}

void ConsiderRowOrder(const std::vector<Mask>& rows,
                      const std::vector<int>& row_order,
                      CanonicalCode* best, bool* initialized) {
  std::array<int, 3> edge_order = {0, 1, 2};
  do {
    const CanonicalCode code = OrderedIncidenceCode(rows, row_order, edge_order);
    if (!*initialized || code < *best) {
      *best = code;
      *initialized = true;
    }
  } while (std::next_permutation(edge_order.begin(), edge_order.end()));
}

CanonicalCode FullCanonicalCode(const std::vector<Mask>& rows) {
  std::vector<int> row_order(rows.size());
  std::iota(row_order.begin(), row_order.end(), 0);
  CanonicalCode best{};
  bool initialized = false;
  do {
    ConsiderRowOrder(rows, row_order, &best, &initialized);
  } while (std::next_permutation(row_order.begin(), row_order.end()));
  if (!initialized) throw std::runtime_error("full canonicalization failed");
  return best;
}

void EnumerateClassOrders(
    const std::vector<Mask>& rows,
    const std::vector<std::vector<int>>& classes, int class_index,
    std::vector<int>* row_order, CanonicalCode* best, bool* initialized) {
  if (class_index == static_cast<int>(classes.size())) {
    ConsiderRowOrder(rows, *row_order, best, initialized);
    return;
  }
  std::vector<int> order = classes[class_index];
  std::sort(order.begin(), order.end());
  do {
    const std::size_t old_size = row_order->size();
    row_order->insert(row_order->end(), order.begin(), order.end());
    EnumerateClassOrders(rows, classes, class_index + 1, row_order, best,
                         initialized);
    row_order->resize(old_size);
  } while (std::next_permutation(order.begin(), order.end()));
}

CanonicalCode IndependentCanonicalCode(const std::vector<int>& chosen) {
  std::vector<Mask> rows;
  rows.reserve(chosen.size());
  std::map<CanonicalCode, std::vector<int>> partition;
  for (int position = 0; position < static_cast<int>(chosen.size()); ++position) {
    const Candidate& candidate = g_candidates[chosen[position]];
    rows.push_back(candidate.face);
    partition[candidate.singleton_code].push_back(position);
  }
  std::vector<std::vector<int>> classes;
  classes.reserve(partition.size());
  for (auto& [code, positions] : partition) {
    (void)code;
    classes.push_back(std::move(positions));
  }
  std::vector<int> row_order;
  row_order.reserve(rows.size());
  CanonicalCode best{};
  bool initialized = false;
  EnumerateClassOrders(rows, classes, 0, &row_order, &best, &initialized);
  if (!initialized) throw std::runtime_error("independent canonicalization failed");
  return best;
}

std::vector<int> PrefixG(Mask linked) {
  const int linked_size = std::popcount(linked);
  const int dimension = 9 - linked_size;
  const int top = 5 - linked_size;
  std::vector<int> f(top + 1);
  f[0] = 1;
  const Mask rest = kFull ^ linked;
  for (Mask extension = rest;; extension = (extension - 1) & rest) {
    const int size = std::popcount(extension);
    if (size >= 1 && size <= top) {
      const Mask face = linked | extension;
      if (Independent(face) && (face & g_four) != g_four) ++f[size];
    }
    if (!extension) break;
  }
  std::vector<int> h(top + 1);
  for (int index = 0; index <= top; ++index) {
    for (int order = 0; order <= index; ++order) {
      h[index] += (((index - order) & 1) ? -1 : 1) *
          static_cast<int>(Choose(dimension - order, dimension - index)) *
          f[order];
    }
  }
  std::vector<int> g(top + 1);
  g[0] = h[0];
  for (int index = 1; index <= top; ++index) g[index] = h[index] - h[index - 1];
  return g;
}

LinkWitness BuildFixedLinkWitness(Mask subset) {
  LinkWitness witness;
  witness.present = true;
  witness.subset = subset;
  witness.dimension = 9 - std::popcount(subset);
  const int top = 5 - std::popcount(subset);
  const Mask rest = kFull ^ subset;
  witness.f_prefix.assign(top + 1, 0);
  witness.f_prefix[0] = 1;
  for (Mask extension = rest;; extension = (extension - 1) & rest) {
    const int size = std::popcount(extension);
    if (size >= 1 && size <= top) {
      const Mask face = subset | extension;
      if (Independent(face) && (face & g_four) != g_four) {
        ++witness.f_prefix[size];
      }
    }
    if (!extension) break;
  }
  witness.h.assign(top + 1, 0);
  for (int index = 0; index <= top; ++index) {
    for (int order = 0; order <= index; ++order) {
      witness.h[index] += (((index - order) & 1) ? -1 : 1) *
          static_cast<int>(Choose(witness.dimension - order,
                                  witness.dimension - index)) *
          witness.f_prefix[order];
    }
  }
  for (int index = 1; index <= top; ++index) {
    witness.g.push_back(witness.h[index] - witness.h[index - 1]);
  }
  return witness;
}

void InitializeLocalConstraints() {
  g_local_index.fill(-1);
  for (int size = 1; size <= 3; ++size) {
    for (Mask subset = 1; subset <= kFull; ++subset) {
      if (std::popcount(subset) != size || !Independent(subset)) continue;
      const std::vector<int> g = PrefixG(subset);
      const int top = 5 - size;
      LocalConstraint constraint;
      constraint.subset = subset;
      constraint.previous = g[top - 1];
      constraint.before = g[top];
      constraint.lower = constraint.previous >= 0
          ? std::max(0, constraint.before -
                            MacaulayRaise(constraint.previous, top - 1))
          : 0;
      constraint.upper = constraint.before - 1;
      g_local_index[subset] = static_cast<int>(g_local_constraints.size());
      g_local_constraints.push_back(constraint);
    }
  }
  if (g_local_constraints.size() != 430) {
    throw std::runtime_error("local face ground changed");
  }
  for (const LocalConstraint& constraint : g_local_constraints) {
    if (constraint.upper < 0) {
      g_fixed_link_witness = BuildFixedLinkWitness(constraint.subset);
      break;
    }
  }
  if (g_orbit_index == 2) {
    const Mask expected = (Mask{1} << 0) | (Mask{1} << 1);
    if (!g_fixed_link_witness.present ||
        g_fixed_link_witness.subset != expected ||
        g_fixed_link_witness.dimension != 7 ||
        g_fixed_link_witness.f_prefix != std::vector<int>({1, 10, 43, 104}) ||
        g_fixed_link_witness.h != std::vector<int>({1, 3, 4, 4}) ||
        g_fixed_link_witness.g != std::vector<int>({2, 1, 0})) {
      throw std::runtime_error("strict two-active link witness changed");
    }
  } else if (g_fixed_link_witness.present) {
    throw std::runtime_error("viable A4 orbit acquired a strict-link failure");
  }
}

std::vector<Permutation> EnumerateStabilizer() {
  std::vector<std::array<int, 6>> active_actions;
  std::array<int, 6> active_image = {0, 1, 2, 11, 12, 13};
  const std::set<Mask> expected_edges(kEdges.begin(), kEdges.end());
  do {
    Permutation partial{};
    std::iota(partial.begin(), partial.end(), 0);
    for (int index = 0; index < 6; ++index) {
      partial[kActive[index]] = static_cast<std::uint8_t>(active_image[index]);
    }
    std::set<Mask> moved;
    for (Mask edge : kEdges) moved.insert(Transform(edge, partial));
    if (moved == expected_edges) active_actions.push_back(active_image);
  } while (std::next_permutation(active_image.begin(), active_image.end()));
  if (active_actions.size() != 48) throw std::runtime_error("active group changed");

  std::vector<Permutation> result;
  std::array<int, 8> isolate_image = kIsolates;
  for (const auto& active : active_actions) {
    isolate_image = kIsolates;
    do {
      Permutation permutation{};
      std::iota(permutation.begin(), permutation.end(), 0);
      for (int index = 0; index < 6; ++index) {
        permutation[kActive[index]] = static_cast<std::uint8_t>(active[index]);
      }
      for (int index = 0; index < 8; ++index) {
        permutation[kIsolates[index]] =
            static_cast<std::uint8_t>(isolate_image[index]);
      }
      if (Transform(g_four, permutation) == g_four) {
        result.push_back(permutation);
      }
    } while (std::next_permutation(isolate_image.begin(), isolate_image.end()));
  }
  if (result.size() != kExpectedStabilizers[g_orbit_index]) {
    throw std::runtime_error("fixed-A4 stabilizer changed");
  }
  return result;
}

std::vector<Mask> LiteralMinimalBlocker(const std::vector<Mask>& clutter) {
  std::vector<Mask> answer;
  for (Mask candidate = 1; candidate <= kFull; ++candidate) {
    if (!HitsAll(candidate, clutter)) continue;
    bool minimal = true;
    for (int vertex = 0; vertex < kN; ++vertex) {
      if ((candidate & (Mask{1} << vertex)) &&
          HitsAll(candidate ^ (Mask{1} << vertex), clutter)) {
        minimal = false;
        break;
      }
    }
    if (minimal) answer.push_back(candidate);
  }
  return answer;
}

std::array<bool, 1 << kN> BuildFaceTable(const std::vector<Mask>& complements) {
  std::array<bool, 1 << kN> table{};
  for (Mask complement : complements) {
    const Mask facet = kFull ^ complement;
    for (Mask subset = facet;; subset = (subset - 1) & facet) {
      table[subset] = true;
      if (!subset) break;
    }
  }
  return table;
}

std::vector<Mask> MinimalNonfaces(const std::array<bool, 1 << kN>& table) {
  std::vector<Mask> answer;
  for (Mask candidate = 1; candidate <= kFull; ++candidate) {
    if (table[candidate]) continue;
    bool minimal = true;
    for (int vertex = 0; vertex < kN; ++vertex) {
      if ((candidate & (Mask{1} << vertex)) &&
          !table[candidate ^ (Mask{1} << vertex)]) {
        minimal = false;
        break;
      }
    }
    if (minimal) answer.push_back(candidate);
  }
  return answer;
}

bool FaceTableGate(const std::array<bool, 1 << kN>& table,
                   std::vector<Mask> expected_nonfaces) {
  for (int size = 1; size <= 9; ++size) {
    int count = 0;
    for (Mask face = 0; face <= kFull; ++face) {
      count += table[face] && std::popcount(face) == size;
    }
    if (count != kExpectedFVector[size - 1]) return false;
  }
  std::sort(expected_nonfaces.begin(), expected_nonfaces.end());
  return MinimalNonfaces(table) == expected_nonfaces;
}

bool RidgeGate(const std::vector<Mask>& complements) {
  std::map<Mask, std::vector<int>> groups;
  for (int index = 0; index < static_cast<int>(complements.size()); ++index) {
    const Mask outside = kFull ^ complements[index];
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (outside & (Mask{1} << vertex)) {
        groups[complements[index] | (Mask{1} << vertex)].push_back(index);
      }
    }
  }
  if (groups.size() != 585) return false;
  std::vector<std::vector<int>> adjacency(complements.size());
  for (const auto& [ridge, members] : groups) {
    (void)ridge;
    if (members.size() != 2) return false;
    adjacency[members[0]].push_back(members[1]);
    adjacency[members[1]].push_back(members[0]);
  }
  for (const auto& neighbors : adjacency) {
    if (neighbors.size() != 9) return false;
  }
  std::vector<bool> reached(complements.size());
  std::queue<int> queue;
  reached[0] = true;
  queue.push(0);
  while (!queue.empty()) {
    const int current = queue.front();
    queue.pop();
    for (int other : adjacency[current]) {
      if (!reached[other]) {
        reached[other] = true;
        queue.push(other);
      }
    }
  }
  return std::ranges::all_of(reached, [](bool value) { return value; });
}

std::vector<int> LinkG(const std::array<bool, 1 << kN>& table, Mask linked) {
  const int linked_size = std::popcount(linked);
  const int dimension = 9 - linked_size;
  const int half = dimension / 2;
  std::vector<int> f(half + 1);
  f[0] = 1;
  const Mask rest = kFull ^ linked;
  for (Mask extension = rest;; extension = (extension - 1) & rest) {
    const int size = std::popcount(extension);
    if (size >= 1 && size <= half && table[linked | extension]) ++f[size];
    if (!extension) break;
  }
  std::vector<int> h(half + 1);
  for (int index = 0; index <= half; ++index) {
    for (int order = 0; order <= index; ++order) {
      h[index] += (((index - order) & 1) ? -1 : 1) *
          static_cast<int>(Choose(dimension - order, dimension - index)) *
          f[order];
    }
  }
  std::vector<int> g;
  for (int index = 1; index <= half; ++index) {
    g.push_back(h[index] - h[index - 1]);
  }
  return g;
}

bool LocalLinkGate(const std::array<bool, 1 << kN>& table) {
  int tested = 0;
  int transitions = 0;
  for (int size = 1; size <= 3; ++size) {
    for (Mask linked = 1; linked <= kFull; ++linked) {
      if (std::popcount(linked) != size || !table[linked]) continue;
      const int dimension = 9 - size;
      const int strict_max = std::min(dimension / 2, dimension - 4);
      const std::vector<int> g = LinkG(table, linked);
      for (int index = 0; index < static_cast<int>(g.size()); ++index) {
        const int order = index + 1;
        if (g[index] < (order <= strict_max ? 1 : 0)) return false;
        if (index > 0) {
          ++transitions;
          if (g[index] > MacaulayRaise(g[index - 1], index)) return false;
        }
      }
      ++tested;
    }
  }
  return tested == 430 && transitions == 874;
}

bool EulerLinkGate(const std::array<bool, 1 << kN>& table) {
  for (Mask linked = 0; linked <= kFull; ++linked) {
    if (!table[linked] || std::popcount(linked) > 8) continue;
    const Mask rest = kFull ^ linked;
    int observed = 0;
    for (Mask extension = rest; extension;
         extension = (extension - 1) & rest) {
      if (table[linked | extension]) {
        observed += (std::popcount(extension) & 1) ? 1 : -1;
      }
    }
    const int link_dimension = 8 - std::popcount(linked);
    const int expected = 1 + ((link_dimension & 1) ? -1 : 1);
    if (observed != expected) return false;
  }
  return true;
}

bool CapCompatible(const State& state, const Candidate& candidate) {
  for (int index : candidate.local_indices) {
    if (state.local_incidence[index] >= g_local_constraints[index].upper) {
      return false;
    }
  }
  return true;
}

void CrosscheckCanonicalizer(const std::vector<int>& chosen) {
  std::vector<Mask> rows;
  for (int candidate : chosen) rows.push_back(g_candidates[candidate].face);
  // The class-ordered form need not be the lexicographically least form from
  // the unrestricted d! enumeration: its class order is deliberately fixed
  // by the invariant singleton codes.  Recompute the unrestricted form here
  // as a separate terminal transcript.  EmitResult below checks both forms
  // on every distinct labelled image of every surviving target.
  (void)FullCanonicalCode(rows);
  ++g_full_canonical_crosschecks;
}

std::vector<int> CandidateIndices(const std::vector<Mask>& rows) {
  std::vector<int> answer;
  answer.reserve(rows.size());
  for (Mask face : rows) {
    const auto iterator = std::lower_bound(
        g_candidates.begin(), g_candidates.end(), face,
        [](const Candidate& candidate, Mask value) {
          return candidate.face < value;
        });
    if (iterator == g_candidates.end() || iterator->face != face) {
      throw std::runtime_error("stabilizer image escaped candidate ground");
    }
    answer.push_back(static_cast<int>(iterator - g_candidates.begin()));
  }
  return answer;
}

void Finish(const State& state) {
  ++g_scheduled_orbits;
  for (int index = 0; index < static_cast<int>(g_local_constraints.size()); ++index) {
    if (state.local_incidence[index] < g_local_constraints[index].lower ||
        state.local_incidence[index] > g_local_constraints[index].upper) {
      ++g_link_rejects;
      return;
    }
  }

  std::vector<Mask> clutter(kEdges.begin(), kEdges.end());
  clutter.push_back(g_four);
  std::vector<Mask> additions;
  for (int index : state.chosen) {
    clutter.push_back(g_candidates[index].face);
    additions.push_back(g_candidates[index].face);
  }
  std::sort(additions.begin(), additions.end());
  const std::vector<Mask> blocker = LiteralMinimalBlocker(clutter);
  ++g_blocker_histogram[static_cast<int>(blocker.size())];
  bool pure = blocker.size() == kScheduledFacets;
  for (Mask complement : blocker) pure = pure && std::popcount(complement) == 5;
  if (!pure) return;
  ++g_pure_blocker_orbits;

  std::vector<Mask> carried;
  for (int cover = 0; cover < static_cast<int>(g_covers.size()); ++cover) {
    if (state.surviving_covers[cover / 64] &
        (std::uint64_t{1} << (cover % 64))) {
      carried.push_back(g_covers[cover]);
    }
  }
  if (carried != blocker) {
    throw std::runtime_error("carried covers and literal blocker disagree");
  }
  const auto table = BuildFaceTable(blocker);
  if (!FaceTableGate(table, clutter)) {
    ++g_face_table_rejects;
    return;
  }
  if (!RidgeGate(blocker)) {
    ++g_ridge_rejects;
    return;
  }
  if (!LocalLinkGate(table) || !EulerLinkGate(table)) {
    ++g_link_rejects;
    return;
  }
  CrosscheckCanonicalizer(state.chosen);
  g_targets.push_back(std::move(additions));
  std::cerr << "orbit " << g_orbit_index << " target " << g_targets.size()
            << " after " << g_nodes << " nodes\n";
}

void Search(const State& state) {
  ++g_nodes;
  if (g_nodes % 1000000 == 0) {
    std::cerr << "orbit=" << g_orbit_index << " nodes=" << g_nodes
              << " depth=" << state.chosen.size() << " seen=";
    for (const auto& level : g_seen) std::cerr << level.size() << ',';
    std::cerr << " targets=" << g_targets.size() << '\n';
  }
  const int depth = static_cast<int>(state.chosen.size());
  const CanonicalCode code = IndependentCanonicalCode(state.chosen);
  if (!g_seen[depth].insert(code).second) {
    ++g_isomorphic_prunes;
    return;
  }
  const int surviving_count = BitCount(state.surviving_covers);
  if (surviving_count < kScheduledFacets) {
    ++g_facet_prunes;
    return;
  }
  SmallBits full;
  for (int index = 0; index < static_cast<int>(g_small.size()); ++index) {
    full.set(index);
  }
  if (depth == g_additions) {
    if (state.destroyed_small == full && surviving_count == kScheduledFacets) {
      Finish(state);
    }
    return;
  }

  std::vector<int> branch;
  SmallBits feasible_union = state.destroyed_small;
  bool has_uncovered = false;
  for (int requirement = 0; requirement < static_cast<int>(g_small.size());
       ++requirement) {
    if (state.destroyed_small[requirement]) continue;
    has_uncovered = true;
    std::vector<int> feasible;
    for (int candidate : g_destroyers[requirement]) {
      if (std::find(state.chosen.begin(), state.chosen.end(), candidate) !=
              state.chosen.end() ||
          !CapCompatible(state, g_candidates[candidate])) {
        continue;
      }
      const CoverBits surviving =
          Meet(state.surviving_covers, g_candidates[candidate].hits_covers);
      if (BitCount(surviving) < kScheduledFacets) continue;
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
    for (int candidate = 0; candidate < static_cast<int>(g_candidates.size());
         ++candidate) {
      if (std::find(state.chosen.begin(), state.chosen.end(), candidate) !=
              state.chosen.end() ||
          !CapCompatible(state, g_candidates[candidate])) {
        continue;
      }
      const CoverBits surviving =
          Meet(state.surviving_covers, g_candidates[candidate].hits_covers);
      if (BitCount(surviving) >= kScheduledFacets) branch.push_back(candidate);
    }
  }

  for (int candidate : branch) {
    State next = state;
    next.chosen.push_back(candidate);
    next.destroyed_small |= g_candidates[candidate].destroys_small;
    for (int local : g_candidates[candidate].local_indices) {
      ++next.local_incidence[local];
    }
    next.surviving_covers =
        Meet(state.surviving_covers, g_candidates[candidate].hits_covers);
    Search(next);
  }
}

void InitializeGround() {
  for (Mask face = 0; face <= kFull; ++face) {
    if (std::popcount(face) == 5 && Independent(face)) {
      g_independent_fives.push_back(face);
    }
  }
  if (g_independent_fives.size() != 1372) {
    throw std::runtime_error("independent-five ground changed");
  }
  for (Mask face : g_independent_fives) {
    if ((face & g_four) == g_four) g_shadow.push_back(face);
  }
  if (static_cast<int>(g_shadow.size()) != kExpectedShadows[g_orbit_index]) {
    throw std::runtime_error("forced A5 shadow count changed");
  }
  g_additions = kScheduledMissingFives - static_cast<int>(g_shadow.size());
  if (g_additions != kExpectedAdditions[g_orbit_index]) {
    throw std::runtime_error("additional A5 schedule changed");
  }

  std::vector<Mask> base(kEdges.begin(), kEdges.end());
  base.push_back(g_four);
  for (Mask transversal = 1; transversal <= kFull; ++transversal) {
    if (std::popcount(transversal) <= 4 && HitsAll(transversal, base)) {
      g_small.push_back(transversal);
    }
  }
  if (g_small.size() > kMaximumSmall) {
    throw std::runtime_error("small-transversal ground escaped bitset");
  }
  for (Mask cover = 0; cover <= kFull; ++cover) {
    if (std::popcount(cover) == 5 && HitsAll(cover, base)) {
      g_covers.push_back(cover);
    }
  }
  InitializeLocalConstraints();
  g_seen.resize(g_additions + 1);
  if (g_fixed_link_witness.present) return;

  const int cover_words = (static_cast<int>(g_covers.size()) + 63) / 64;
  for (Mask face : g_independent_fives) {
    if ((face & g_four) == g_four) continue;
    Candidate candidate;
    candidate.face = face;
    int local_count = 0;
    bool allowed = true;
    for (Mask subset = face; subset; subset = (subset - 1) & face) {
      if (std::popcount(subset) > 3) continue;
      const int index = g_local_index[subset];
      if (index < 0) throw std::runtime_error("candidate local subset absent");
      candidate.local_indices[local_count++] = index;
      if (g_local_constraints[index].upper < 1) allowed = false;
    }
    if (local_count != 25) throw std::runtime_error("candidate local count changed");
    if (!allowed) continue;
    for (int requirement = 0; requirement < static_cast<int>(g_small.size());
         ++requirement) {
      if (!(face & g_small[requirement])) {
        candidate.destroys_small.set(requirement);
      }
    }
    candidate.hits_covers.assign(cover_words, 0);
    for (int cover = 0; cover < static_cast<int>(g_covers.size()); ++cover) {
      if (face & g_covers[cover]) {
        candidate.hits_covers[cover / 64] |=
            std::uint64_t{1} << (cover % 64);
      }
    }
    candidate.singleton_code = FullCanonicalCode({face});
    g_candidates.push_back(std::move(candidate));
  }
  g_destroyers.assign(g_small.size(), {});
  for (int candidate = 0; candidate < static_cast<int>(g_candidates.size());
       ++candidate) {
    for (int requirement = 0; requirement < static_cast<int>(g_small.size());
         ++requirement) {
      if (g_candidates[candidate].destroys_small[requirement]) {
        g_destroyers[requirement].push_back(candidate);
      }
    }
  }
}

std::uint64_t TranscriptCode(const std::set<std::vector<Mask>>& families) {
  std::uint64_t value = 1469598103934665603ULL;
  for (const auto& family : families) {
    for (Mask face : family) {
      value ^= face & 0xffu;
      value *= 1099511628211ULL;
      value ^= face >> 8;
      value *= 1099511628211ULL;
    }
    value ^= 0xffu;
    value *= 1099511628211ULL;
  }
  return value;
}

void PrintMask(Mask value) {
  std::cout << '[';
  bool first = true;
  for (int vertex = 0; vertex < kN; ++vertex) {
    if (!(value & (Mask{1} << vertex))) continue;
    if (!first) std::cout << ',';
    std::cout << vertex;
    first = false;
  }
  std::cout << ']';
}

void PrintIntVector(const std::vector<int>& values) {
  std::cout << '[';
  for (int index = 0; index < static_cast<int>(values.size()); ++index) {
    if (index) std::cout << ',';
    std::cout << values[index];
  }
  std::cout << ']';
}

void EmitResult() {
  std::sort(g_targets.begin(), g_targets.end());
  const std::vector<Permutation> stabilizer = EnumerateStabilizer();
  std::set<std::vector<Mask>> literal_union;
  std::vector<int> orbit_sizes;
  std::vector<std::uint64_t> orbit_transcripts;
  std::map<int, int> orbit_size_histogram;
  for (const std::vector<Mask>& target : g_targets) {
    std::set<std::vector<Mask>> orbit;
    const CanonicalCode target_code = FullCanonicalCode(target);
    const CanonicalCode target_class_code =
        IndependentCanonicalCode(CandidateIndices(target));
    for (const Permutation& permutation : stabilizer) {
      std::vector<Mask> image;
      for (Mask face : target) image.push_back(Transform(face, permutation));
      std::sort(image.begin(), image.end());
      orbit.insert(std::move(image));
    }
    for (const std::vector<Mask>& image : orbit) {
      if (FullCanonicalCode(image) != target_code ||
          IndependentCanonicalCode(CandidateIndices(image)) !=
              target_class_code) {
        throw std::runtime_error("explicit stabilizer broke a canonical form");
      }
    }
    orbit_sizes.push_back(static_cast<int>(orbit.size()));
    orbit_transcripts.push_back(TranscriptCode(orbit));
    ++orbit_size_histogram[static_cast<int>(orbit.size())];
    for (const auto& image : orbit) {
      if (!literal_union.insert(image).second) {
        throw std::runtime_error("target orbits overlap literally");
      }
    }
  }

  std::cout << "{\"kind\":\"h9_14_record30_a1_explicit_incidence_exhaustion\""
            << ",\"compact_only\":true"
            << ",\"uses_nauty\":false"
            << ",\"positive_A4_orbit_index\":" << g_orbit_index
            << ",\"fixed_A4\":";
  PrintMask(g_four);
  std::cout << ",\"corrected_forced_nonminimal_A5_shadow_count\":"
            << g_shadow.size()
            << ",\"additional_minimal_A5_count\":" << g_additions
            << ",\"scheduled_missing_independent_five_count\":"
            << kScheduledMissingFives
            << ",\"scheduled_facet_count\":" << kScheduledFacets
            << ",\"canonicalizer\":\"3! edge-row orders; unique A4 row; permutations within invariant singleton-A5 classes; sorted 14-column signatures\""
            << ",\"explicit_stabilizer_order\":" << stabilizer.size()
            << ",\"base_small_transversal_count\":" << g_small.size()
            << ",\"initial_order_five_cover_count\":" << g_covers.size()
            << ",\"candidate_count\":" << g_candidates.size()
            << ",\"raw_unordered_addition_count\":"
            << Choose(static_cast<int>(g_candidates.size()), g_additions)
            << ",\"strict_fixed_link_exclusion\":"
            << (g_fixed_link_witness.present ? "true" : "false");
  if (g_fixed_link_witness.present) {
    std::cout << ",\"strict_link_witness\":{\"subset\":";
    PrintMask(g_fixed_link_witness.subset);
    std::cout << ",\"polytope_dimension\":" << g_fixed_link_witness.dimension
              << ",\"f_prefix\":";
    PrintIntVector(g_fixed_link_witness.f_prefix);
    std::cout << ",\"h\":";
    PrintIntVector(g_fixed_link_witness.h);
    std::cout << ",\"g_without_g0\":";
    PrintIntVector(g_fixed_link_witness.g);
    std::cout << ",\"strict_entry\":3,\"observed\":0,"
                 "\"equality_forced_missing_face_minimum_order\":6,"
                 "\"inherited_maximum_missing_face_order\":5}";
  }
  std::cout << ",\"search_node_count\":" << g_nodes
            << ",\"canonical_state_counts\":[";
  for (int depth = 0; depth <= g_additions; ++depth) {
    if (depth) std::cout << ',';
    std::cout << g_seen[depth].size();
  }
  std::cout << "]"
            << ",\"isomorphic_prune_count\":" << g_isomorphic_prunes
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
            << ",\"pure_blocker_orbit_count\":" << g_pure_blocker_orbits
            << ",\"face_table_reject_count\":" << g_face_table_rejects
            << ",\"ridge_reject_count\":" << g_ridge_rejects
            << ",\"necessary_link_reject_count\":" << g_link_rejects
            << ",\"terminal_full_d_factorial_canonical_recomputation_count\":"
            << g_full_canonical_crosschecks
            << ",\"necessary_incidence_target_orbit_count\":"
            << g_targets.size()
            << ",\"literal_target_count\":" << literal_union.size()
            << ",\"literal_union_FNV1a64\":" << TranscriptCode(literal_union)
            << ",\"target_orbit_size_histogram\":{";
  first = true;
  for (const auto& [size, count] : orbit_size_histogram) {
    if (!first) std::cout << ',';
    std::cout << '\"' << size << "\":" << count;
    first = false;
  }
  std::cout << "},\"target_orbit_sizes\":[";
  for (int index = 0; index < static_cast<int>(orbit_sizes.size()); ++index) {
    if (index) std::cout << ',';
    std::cout << orbit_sizes[index];
  }
  std::cout << "],\"target_orbit_FNV1a64\":[";
  for (int index = 0; index < static_cast<int>(orbit_transcripts.size()); ++index) {
    if (index) std::cout << ',';
    std::cout << orbit_transcripts[index];
  }
  std::cout << "],\"target_addition_representatives\":[";
  for (int index = 0; index < static_cast<int>(g_targets.size()); ++index) {
    if (index) std::cout << ',';
    std::cout << '[';
    for (int row = 0; row < static_cast<int>(g_targets[index].size()); ++row) {
      if (row) std::cout << ',';
      PrintMask(g_targets[index][row]);
    }
    std::cout << ']';
  }
  std::cout << "]"
            << ",\"independent_search_exhaustion_replay\":true"
            << ",\"scope\":\"Standalone exact compact Record-30 a=1 incidence exhaustion. No nauty, SAT/SMT, PL, convex-polytopal, Coxeter-label, or finite-volume inference.\"}\n";
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 2) {
      throw std::runtime_error("usage: replay positive_A4_orbit_index(0..2)");
    }
    g_orbit_index = std::stoi(argv[1]);
    if (g_orbit_index < 0 || g_orbit_index >= 3) {
      throw std::runtime_error("positive A4 orbit index must be 0, 1, or 2");
    }
    g_four = kOrbitFours[g_orbit_index];
    InitializeGround();
    if (!g_fixed_link_witness.present) {
      State initial;
      initial.surviving_covers.assign((g_covers.size() + 63) / 64, 0);
      for (int cover = 0; cover < static_cast<int>(g_covers.size()); ++cover) {
        initial.surviving_covers[cover / 64] |=
            std::uint64_t{1} << (cover % 64);
      }
      Search(initial);
    }
    EmitResult();
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
