// Standalone no-nauty replay of the three exploratory Record-30 a=1 targets.
//
// The program writes the three representatives literally, reconstructs the
// full fixed-structure stabilizers inside Aut(3K2) x S8, expands the target
// orbits to 48 + 80 + 40 = 168 labelled families, and audits every labelled
// family through the literal 2^14 blocker, face table, and ridge equations.
// Local g/Macaulay and every-face Euler-link gates are recomputed on each of
// the three representatives and transported through the explicitly verified
// vertex permutations.  It also independently reconstructs the strict
// {0,1}-link obstruction for the third A4 orbit.
//
// This validates the displayed targets and their literal expansions.  It is
// not an independent replay of the nauty search exhaustion and makes no PL,
// convex-polytopal, Coxeter-label, or finite-volume claim.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <stdexcept>
#include <vector>

namespace {

constexpr int kN = 14;
using Mask = std::uint16_t;
using Permutation = std::array<std::uint8_t, kN>;
constexpr Mask kFull = (Mask{1} << kN) - 1;
constexpr std::array<Mask, 3> kEdges = {
    (Mask{1} << 0) | (Mask{1} << 11),
    (Mask{1} << 1) | (Mask{1} << 12),
    (Mask{1} << 2) | (Mask{1} << 13),
};
constexpr std::array<int, 6> kActive = {0, 1, 2, 11, 12, 13};
constexpr std::array<int, 8> kIsolates = {3, 4, 5, 6, 7, 8, 9, 10};
constexpr std::array<int, 9> kExpectedFVector = {
    14, 88, 328, 805, 1358, 1582, 1232, 585, 130,
};

struct Target {
  Mask four;
  std::vector<Mask> additions;
  int expected_stabilizer;
  int expected_orbit_size;
};

Mask MakeMask(std::initializer_list<int> vertices) {
  Mask value = 0;
  for (int vertex : vertices) value |= Mask{1} << vertex;
  return value;
}

std::vector<Target> Targets() {
  return {
      {
          MakeMask({3, 4, 5, 6}),
          {
              MakeMask({1, 7, 8, 9, 10}),
              MakeMask({4, 7, 8, 9, 10}),
              MakeMask({7, 8, 9, 10, 11}),
              MakeMask({0, 3, 5, 6, 12}),
          },
          27648,
          48,
      },
      {
          MakeMask({0, 3, 4, 5}),
          {
              MakeMask({3, 4, 5, 7, 8}),
              MakeMask({6, 7, 8, 9, 10}),
              MakeMask({1, 6, 9, 10, 11}),
              MakeMask({6, 8, 9, 10, 11}),
              MakeMask({3, 4, 5, 7, 12}),
          },
          5760,
          80,
      },
      {
          MakeMask({0, 3, 4, 5}),
          {
              MakeMask({3, 4, 5, 7, 8}),
              MakeMask({6, 7, 8, 9, 10}),
              MakeMask({2, 6, 9, 10, 11}),
              MakeMask({6, 9, 10, 11, 12}),
              MakeMask({1, 3, 4, 5, 13}),
          },
          5760,
          40,
      },
  };
}

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
  if (remaining) throw std::runtime_error("Macaulay expansion failed");
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

std::vector<Permutation> Stabilizer(Mask four) {
  std::vector<std::array<int, 6>> active_actions;
  std::array<int, 6> active_image = {0, 1, 2, 11, 12, 13};
  const std::set<Mask> expected_edges(kEdges.begin(), kEdges.end());
  do {
    Permutation partial{};
    std::iota(partial.begin(), partial.end(), 0);
    for (int index = 0; index < 6; ++index) {
      partial[kActive[index]] = static_cast<std::uint8_t>(active_image[index]);
    }
    std::set<Mask> image;
    for (Mask edge : kEdges) image.insert(Transform(edge, partial));
    if (image == expected_edges) active_actions.push_back(active_image);
  } while (std::next_permutation(active_image.begin(), active_image.end()));
  if (active_actions.size() != 48) throw std::runtime_error("active group changed");

  std::vector<Permutation> answer;
  std::array<int, 8> isolate_image{};
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
      if (Transform(four, permutation) == four) answer.push_back(permutation);
    } while (std::next_permutation(isolate_image.begin(), isolate_image.end()));
  }
  return answer;
}

std::vector<Mask> MinimalBlocker(const std::vector<Mask>& clutter) {
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

std::array<bool, 1 << kN> FaceTable(const std::vector<Mask>& complements) {
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

bool FaceGate(const std::array<bool, 1 << kN>& table,
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
  const int size = std::popcount(linked);
  const int dimension = 9 - size;
  const int half = dimension / 2;
  std::vector<int> f(half + 1);
  f[0] = 1;
  const Mask rest = kFull ^ linked;
  for (Mask extension = rest;; extension = (extension - 1) & rest) {
    const int order = std::popcount(extension);
    if (order >= 1 && order <= half && table[linked | extension]) ++f[order];
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
      const auto g = LinkG(table, linked);
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
  int tested = 0;
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
    const int dimension = 8 - std::popcount(linked);
    const int expected = 1 + ((dimension & 1) ? -1 : 1);
    if (observed != expected) return false;
    ++tested;
  }
  return tested == 5993;
}

void AuditLiteral(Mask four, const std::vector<Mask>& additions,
                  bool expensive_links) {
  std::vector<Mask> clutter(kEdges.begin(), kEdges.end());
  clutter.push_back(four);
  clutter.insert(clutter.end(), additions.begin(), additions.end());
  const auto blocker = MinimalBlocker(clutter);
  if (blocker.size() != 130 ||
      std::ranges::any_of(blocker, [](Mask value) {
        return std::popcount(value) != 5;
      })) {
    throw std::runtime_error("literal pure blocker failed");
  }
  const auto table = FaceTable(blocker);
  if (!FaceGate(table, clutter)) throw std::runtime_error("literal face table failed");
  if (!RidgeGate(blocker)) throw std::runtime_error("literal ridge audit failed");
  if (expensive_links && (!LocalLinkGate(table) || !EulerLinkGate(table))) {
    throw std::runtime_error("representative link audit failed");
  }
}

std::vector<int> ResidualLinkG(Mask four, Mask linked) {
  const int size = std::popcount(linked);
  const int dimension = 9 - size;
  const int top = 5 - size;
  std::vector<int> f(top + 1);
  f[0] = 1;
  const Mask rest = kFull ^ linked;
  for (Mask extension = rest;; extension = (extension - 1) & rest) {
    const int order = std::popcount(extension);
    const Mask face = linked | extension;
    if (order >= 1 && order <= top && Independent(face) &&
        (face & four) != four) {
      ++f[order];
    }
    if (!extension) break;
  }
  if (f != std::vector<int>({1, 10, 43, 104})) {
    throw std::runtime_error("residual link f-prefix changed");
  }
  std::vector<int> h(top + 1);
  for (int index = 0; index <= top; ++index) {
    for (int order = 0; order <= index; ++order) {
      h[index] += (((index - order) & 1) ? -1 : 1) *
          static_cast<int>(Choose(dimension - order, dimension - index)) *
          f[order];
    }
  }
  if (h != std::vector<int>({1, 3, 4, 4})) {
    throw std::runtime_error("residual link h-vector changed");
  }
  std::vector<int> g;
  for (int index = 1; index <= top; ++index) g.push_back(h[index] - h[index - 1]);
  return g;
}

std::uint64_t TranscriptCode(const std::set<std::vector<Mask>>& families) {
  // Deterministic FNV-1a transcript fingerprint; the literal rows are audited
  // directly above, so this checksum is only a compact reproducibility aid.
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

}  // namespace

int main() {
  try {
    const auto targets = Targets();
    std::vector<int> stabilizer_orders;
    std::vector<int> orbit_sizes;
    std::vector<std::uint64_t> transcript_codes;
    int literal_total = 0;
    for (int index = 0; index < static_cast<int>(targets.size()); ++index) {
      const auto& target = targets[index];
      const auto stabilizer = Stabilizer(target.four);
      if (static_cast<int>(stabilizer.size()) != target.expected_stabilizer) {
        throw std::runtime_error("target stabilizer order changed");
      }
      stabilizer_orders.push_back(static_cast<int>(stabilizer.size()));
      std::set<std::vector<Mask>> orbit;
      for (const auto& permutation : stabilizer) {
        std::vector<Mask> image;
        for (Mask face : target.additions) image.push_back(Transform(face, permutation));
        std::sort(image.begin(), image.end());
        orbit.insert(std::move(image));
      }
      if (static_cast<int>(orbit.size()) != target.expected_orbit_size) {
        throw std::runtime_error("literal target orbit size changed");
      }
      orbit_sizes.push_back(static_cast<int>(orbit.size()));
      literal_total += static_cast<int>(orbit.size());
      transcript_codes.push_back(TranscriptCode(orbit));
      bool first = true;
      for (const auto& family : orbit) {
        AuditLiteral(target.four, family, first);
        first = false;
      }
    }
    if (literal_total != 168) throw std::runtime_error("literal total changed");

    const Mask excluded_four = MakeMask({0, 1, 3, 4});
    const Mask linked = MakeMask({0, 1});
    const auto residual_g = ResidualLinkG(excluded_four, linked);
    if (residual_g != std::vector<int>({2, 1, 0})) {
      throw std::runtime_error("strict residual link obstruction changed");
    }

    std::cout
        << "{\"kind\":\"h9_14_record30_a1_no_nauty_literal_target_replay\""
        << ",\"compact_only\":true"
        << ",\"uses_nauty\":false"
        << ",\"graph_automorphism_candidates_checked_per_A4_orbit\":"
        << 48 * 40320
        << ",\"target_representative_count\":" << targets.size()
        << ",\"stabilizer_orders\":[";
    for (int index = 0; index < static_cast<int>(stabilizer_orders.size()); ++index) {
      if (index) std::cout << ',';
      std::cout << stabilizer_orders[index];
    }
    std::cout << "],\"literal_orbit_sizes\":[";
    for (int index = 0; index < static_cast<int>(orbit_sizes.size()); ++index) {
      if (index) std::cout << ',';
      std::cout << orbit_sizes[index];
    }
    std::cout << "],\"literal_target_count\":" << literal_total
              << ",\"literal_family_FNV1a64\":[";
    for (int index = 0; index < static_cast<int>(transcript_codes.size()); ++index) {
      if (index) std::cout << ',';
      std::cout << transcript_codes[index];
    }
    std::cout
        << "],\"literal_blocker_face_ridge_audit_count\":" << literal_total
        << ",\"representative_local_g_Macaulay_audit_count\":3"
        << ",\"representative_every_face_Euler_audit_count\":3"
        << ",\"Euler_face_links_per_representative\":5993"
        << ",\"strict_two_active_A4_link_witness\":{"
           "\"subset\":[0,1],\"polytope_dimension\":7,"
           "\"f_prefix\":[1,10,43,104],\"h\":[1,3,4,4],"
           "\"g_without_g0\":[2,1,0],\"strict_g3_violation\":true,"
           "\"equality_forced_missing_face_minimum_order\":6,"
           "\"inherited_maximum_missing_face_order\":5}"
        << ",\"all_literal_targets_pass\":true"
        << ",\"independent_search_exhaustion_replay\":false"
        << ",\"scope\":\"No-nauty explicit-group replay of all 168 labelled expansions of the three displayed Record-30 a=1 targets. This validates the target orbits and independently reconstructs the strict-link exclusion, but does not replay quotient-search exhaustion. No PL, Coxeter-label, or finite-volume claim.\"}\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
