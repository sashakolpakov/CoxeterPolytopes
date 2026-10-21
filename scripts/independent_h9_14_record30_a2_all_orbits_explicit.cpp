// Standalone explicit colored-incidence replay for all 21 compact Record-30
// a=2 positive-A4 orbits.  No nauty, SAT/SMT solver, or project code is used.
//
// Input:
//   21
//   orbit_index first_A4_mask second_A4_mask
//
// The canonical form is exact: it enumerates the 3! missing-edge rows, the
// 2! A4 rows, and the d! unordered added-A5 rows, then sorts the fourteen
// facet-column signatures.  Thus equality is precisely isomorphism of the
// same colored incidence graph sent to nauty by the production search.

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
constexpr int kScheduledFives = 20;
constexpr int kScheduledFacets = 128;
using Mask = std::uint16_t;
constexpr Mask kFull = (Mask{1} << kN) - 1;
constexpr std::array<Mask, 3> kEdges = {
    (Mask{1} << 0) | (Mask{1} << 11),
    (Mask{1} << 1) | (Mask{1} << 12),
    (Mask{1} << 2) | (Mask{1} << 13),
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

struct Audit {
  int orbit_index = -1;
  std::array<Mask, 2> fours{};
  int additions = 0;
  std::uint64_t stabilizer_order = 0;
  std::vector<Permutation> stabilizer;
  std::vector<Mask> shadow;
  std::vector<Mask> small;
  std::vector<Mask> covers;
  std::vector<Candidate> candidates;
  std::vector<std::vector<int>> destroyers;
  std::array<int, kN> upper{};
  std::vector<std::unordered_set<std::string>> seen;
  std::uint64_t nodes = 0;
  std::uint64_t facet_prunes = 0;
  std::uint64_t uncoverable_prunes = 0;
  std::uint64_t scheduled_orbits = 0;
  std::map<int, std::uint64_t> blocker_histogram;
  std::vector<std::vector<Mask>> targets;
  std::vector<int> target_orbit_sizes;
  std::map<int, int> orbit_size_histogram;
  std::size_t literal_target_count = 0;
};

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

Mask PermutedMask(Mask value, const Permutation& permutation) {
  Mask answer = 0;
  for (int vertex = 0; vertex < kN; ++vertex) {
    if ((value >> vertex) & 1u) answer |= Mask{1} << permutation[vertex];
  }
  return answer;
}

bool SameMaskSet(std::array<Mask, 2> left, std::array<Mask, 2> right) {
  std::sort(left.begin(), left.end());
  std::sort(right.begin(), right.end());
  return left == right;
}

std::vector<Permutation> GraphAutomorphisms() {
  constexpr std::array<std::array<int, 2>, 3> endpoints = {
      std::array<int, 2>{0, 11},
      std::array<int, 2>{1, 12},
      std::array<int, 2>{2, 13},
  };
  std::array<int, 3> component_image = {0, 1, 2};
  std::array<int, 8> isolate_image = {3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<Permutation> answer;
  answer.reserve(48 * 40320);
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
        answer.push_back(permutation);
      } while (std::next_permutation(isolate_image.begin(),
                                     isolate_image.end()));
    }
  } while (std::next_permutation(component_image.begin(),
                                 component_image.end()));
  if (answer.size() != 1935360) {
    throw std::runtime_error("3K2 + 8K1 automorphism count changed");
  }
  return answer;
}

bool Independent(Mask value) {
  for (Mask edge : kEdges) {
    if ((value & edge) == edge) return false;
  }
  return true;
}

bool ContainsFour(Mask value, const std::array<Mask, 2>& fours) {
  for (Mask four : fours) {
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
  for (int index = 0; index < static_cast<int>(left.size()); ++index) {
    answer[index] = left[index] & right[index];
  }
  return answer;
}

std::string CanonicalCode(const Audit& audit, const std::vector<int>& chosen) {
  std::array<int, 3> edge_order = {0, 1, 2};
  std::array<int, 2> four_order = {0, 1};
  std::vector<int> row_order(chosen.size());
  std::iota(row_order.begin(), row_order.end(), 0);
  std::string best;
  bool initialized = false;
  do {
    do {
      do {
        std::array<std::uint16_t, kN> signatures{};
        for (int position = 0; position < 3; ++position) {
          for (int vertex = 0; vertex < kN; ++vertex) {
            if ((kEdges[edge_order[position]] >> vertex) & 1u) {
              signatures[vertex] |= std::uint16_t{1} << position;
            }
          }
        }
        for (int position = 0; position < 2; ++position) {
          for (int vertex = 0; vertex < kN; ++vertex) {
            if ((audit.fours[four_order[position]] >> vertex) & 1u) {
              signatures[vertex] |= std::uint16_t{1} << (3 + position);
            }
          }
        }
        for (int position = 0;
             position < static_cast<int>(row_order.size()); ++position) {
          const Mask face =
              audit.candidates[chosen[row_order[position]]].face;
          for (int vertex = 0; vertex < kN; ++vertex) {
            if ((face >> vertex) & 1u) {
              signatures[vertex] |= std::uint16_t{1} << (5 + position);
            }
          }
        }
        std::sort(signatures.begin(), signatures.end());
        std::string code;
        code.reserve(2 * kN);
        for (std::uint16_t signature : signatures) {
          code.push_back(static_cast<char>(signature & 255));
          code.push_back(static_cast<char>(signature >> 8));
        }
        if (!initialized || code < best) {
          best = std::move(code);
          initialized = true;
        }
      } while (std::next_permutation(row_order.begin(), row_order.end()));
    } while (std::next_permutation(four_order.begin(), four_order.end()));
  } while (std::next_permutation(edge_order.begin(), edge_order.end()));
  if (!initialized) throw std::runtime_error("canonical form failed");
  return best;
}

bool CapCompatible(const Audit& audit, const State& state, Mask face) {
  for (int vertex = 0; vertex < kN; ++vertex) {
    if (static_cast<int>(state.incidence[vertex]) +
            static_cast<int>((face >> vertex) & 1u) >
        audit.upper[vertex]) {
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

void Finish(Audit* audit, const State& state) {
  ++audit->scheduled_orbits;
  std::vector<Mask> clutter(kEdges.begin(), kEdges.end());
  clutter.insert(clutter.end(), audit->fours.begin(), audit->fours.end());
  std::vector<Mask> family;
  for (int index : state.chosen) {
    clutter.push_back(audit->candidates[index].face);
    family.push_back(audit->candidates[index].face);
  }
  const auto blocker = LiteralMinimalBlocker(clutter);
  ++audit->blocker_histogram[static_cast<int>(blocker.size())];
  bool pure = blocker.size() == kScheduledFacets;
  for (Mask transversal : blocker) pure = pure && std::popcount(transversal) == 5;
  if (!pure) return;
  if (BitCount(state.surviving_covers) != static_cast<int>(blocker.size())) {
    throw std::runtime_error("cover intersection and blocker disagree");
  }
  std::sort(family.begin(), family.end());
  audit->targets.push_back(std::move(family));
}

void Search(Audit* audit, const State& state) {
  ++audit->nodes;
  const int depth = static_cast<int>(state.chosen.size());
  if (!audit->seen[depth].insert(CanonicalCode(*audit, state.chosen)).second) {
    return;
  }
  const int surviving_count = BitCount(state.surviving_covers);
  if (surviving_count < kScheduledFacets) {
    ++audit->facet_prunes;
    return;
  }
  const std::uint64_t full =
      (std::uint64_t{1} << static_cast<int>(audit->small.size())) - 1;
  if (depth == audit->additions) {
    if (state.destroyed_small == full && surviving_count == kScheduledFacets) {
      Finish(audit, state);
    }
    return;
  }

  std::vector<int> branch;
  std::uint64_t feasible_union = state.destroyed_small;
  bool has_uncovered = false;
  for (int requirement = 0;
       requirement < static_cast<int>(audit->small.size()); ++requirement) {
    if ((state.destroyed_small >> requirement) & 1u) continue;
    has_uncovered = true;
    std::vector<int> feasible;
    for (int candidate : audit->destroyers[requirement]) {
      if (std::find(state.chosen.begin(), state.chosen.end(), candidate) !=
              state.chosen.end() ||
          !CapCompatible(*audit, state, audit->candidates[candidate].face)) {
        continue;
      }
      const auto surviving = Meet(
          state.surviving_covers, audit->candidates[candidate].hits_covers);
      if (BitCount(surviving) < kScheduledFacets) continue;
      feasible.push_back(candidate);
      feasible_union |= audit->candidates[candidate].destroys_small;
    }
    if (feasible.empty()) {
      ++audit->uncoverable_prunes;
      return;
    }
    if (branch.empty() || feasible.size() < branch.size()) {
      branch = std::move(feasible);
    }
  }
  if ((feasible_union & full) != full) {
    ++audit->uncoverable_prunes;
    return;
  }
  if (!has_uncovered) {
    for (int candidate = 0;
         candidate < static_cast<int>(audit->candidates.size()); ++candidate) {
      if (std::find(state.chosen.begin(), state.chosen.end(), candidate) !=
              state.chosen.end() ||
          !CapCompatible(*audit, state, audit->candidates[candidate].face)) {
        continue;
      }
      const auto surviving = Meet(
          state.surviving_covers, audit->candidates[candidate].hits_covers);
      if (BitCount(surviving) >= kScheduledFacets) branch.push_back(candidate);
    }
  }
  for (int candidate : branch) {
    State next = state;
    next.chosen.push_back(candidate);
    next.destroyed_small |= audit->candidates[candidate].destroys_small;
    for (int vertex = 0; vertex < kN; ++vertex) {
      next.incidence[vertex] +=
          (audit->candidates[candidate].face >> vertex) & 1u;
    }
    next.surviving_covers = Meet(
        state.surviving_covers, audit->candidates[candidate].hits_covers);
    Search(audit, next);
  }
}

Audit Run(
    int orbit_index, std::array<Mask, 2> fours,
    const std::vector<Mask>& independent_fives,
    const std::vector<Permutation>& graph_automorphisms) {
  Audit audit;
  audit.orbit_index = orbit_index;
  std::sort(fours.begin(), fours.end());
  audit.fours = fours;
  for (const auto& permutation : graph_automorphisms) {
    const std::array<Mask, 2> image = {
        PermutedMask(fours[0], permutation),
        PermutedMask(fours[1], permutation),
    };
    if (SameMaskSet(image, fours)) audit.stabilizer.push_back(permutation);
  }
  audit.stabilizer_order = audit.stabilizer.size();
  if (!audit.stabilizer_order) throw std::runtime_error("empty stabilizer");

  for (Mask face : independent_fives) {
    if (ContainsFour(face, fours)) audit.shadow.push_back(face);
  }
  audit.additions = kScheduledFives - static_cast<int>(audit.shadow.size());
  if (audit.additions < 0 || audit.additions > 5) {
    throw std::runtime_error("addition schedule changed");
  }
  audit.seen.resize(audit.additions + 1);
  for (int vertex = 0; vertex < kN; ++vertex) {
    int four_incidence = 0;
    int shadow_incidence = 0;
    for (Mask four : fours) four_incidence += (four >> vertex) & 1u;
    for (Mask face : audit.shadow) shadow_incidence += (face >> vertex) & 1u;
    const int previous = kBaseG3[vertex] - four_incidence;
    const int before =
        kBaseG4[vertex] + 6 * four_incidence - shadow_incidence;
    const int lower = std::max(0, before - MacaulayRaise(previous, 3));
    if (lower != 0) throw std::runtime_error("local lower cap changed");
    audit.upper[vertex] = before - 1;
  }

  std::vector<Mask> base(kEdges.begin(), kEdges.end());
  base.insert(base.end(), fours.begin(), fours.end());
  base.insert(base.end(), audit.shadow.begin(), audit.shadow.end());
  for (Mask transversal = 1; transversal <= kFull; ++transversal) {
    if (std::popcount(transversal) <= 4 && HitsAll(transversal, base)) {
      audit.small.push_back(transversal);
    }
  }
  if (audit.small.size() > 63) {
    throw std::runtime_error("small-transversal gate escaped one word");
  }
  for (Mask cover = 0; cover <= kFull; ++cover) {
    if (std::popcount(cover) != 5) continue;
    bool hits = true;
    for (Mask edge : kEdges) hits = hits && (cover & edge);
    for (Mask four : fours) hits = hits && (cover & four);
    if (hits) audit.covers.push_back(cover);
  }

  for (Mask face : independent_fives) {
    if (ContainsFour(face, fours)) continue;
    bool allowed = true;
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (((face >> vertex) & 1u) && audit.upper[vertex] < 1) allowed = false;
    }
    if (!allowed) continue;
    Candidate candidate;
    candidate.face = face;
    for (int requirement = 0;
         requirement < static_cast<int>(audit.small.size()); ++requirement) {
      if (!(face & audit.small[requirement])) {
        candidate.destroys_small |= std::uint64_t{1} << requirement;
      }
    }
    candidate.hits_covers.assign((audit.covers.size() + 63) / 64, 0);
    for (int cover = 0; cover < static_cast<int>(audit.covers.size()); ++cover) {
      if (face & audit.covers[cover]) {
        candidate.hits_covers[cover / 64] |=
            std::uint64_t{1} << (cover % 64);
      }
    }
    audit.candidates.push_back(std::move(candidate));
  }
  audit.destroyers.assign(audit.small.size(), {});
  for (int candidate = 0;
       candidate < static_cast<int>(audit.candidates.size()); ++candidate) {
    for (int requirement = 0;
         requirement < static_cast<int>(audit.small.size()); ++requirement) {
      if ((audit.candidates[candidate].destroys_small >> requirement) & 1u) {
        audit.destroyers[requirement].push_back(candidate);
      }
    }
  }

  State initial;
  initial.surviving_covers.assign((audit.covers.size() + 63) / 64, 0);
  for (int cover = 0; cover < static_cast<int>(audit.covers.size()); ++cover) {
    initial.surviving_covers[cover / 64] |=
        std::uint64_t{1} << (cover % 64);
  }
  Search(&audit, initial);
  std::sort(audit.targets.begin(), audit.targets.end());
  if (audit.targets.size() != audit.scheduled_orbits) {
    throw std::runtime_error("scheduled target lost blocker purity");
  }

  std::set<std::vector<Mask>> literal_union;
  for (const auto& family : audit.targets) {
    std::set<std::vector<Mask>> orbit;
    for (const auto& permutation : audit.stabilizer) {
      std::vector<Mask> image;
      for (Mask face : family) image.push_back(PermutedMask(face, permutation));
      std::sort(image.begin(), image.end());
      orbit.insert(image);
    }
    audit.target_orbit_sizes.push_back(static_cast<int>(orbit.size()));
    ++audit.orbit_size_histogram[static_cast<int>(orbit.size())];
    for (const auto& image : orbit) {
      if (!literal_union.insert(image).second) {
        throw std::runtime_error("target stabilizer orbits overlap");
      }
    }
  }
  audit.literal_target_count = literal_union.size();
  return audit;
}

void PrintFamily(const std::vector<Mask>& family) {
  std::cout << '[';
  for (int index = 0; index < static_cast<int>(family.size()); ++index) {
    if (index) std::cout << ',';
    std::cout << family[index];
  }
  std::cout << ']';
}

void PrintAudit(const Audit& audit) {
  std::cout << "{\"positive_A4_orbit_index\":" << audit.orbit_index
            << ",\"fixed_A4_masks\":[" << audit.fours[0] << ','
            << audit.fours[1] << ']'
            << ",\"stabilizer_order\":" << audit.stabilizer_order
            << ",\"forced_shadow_count\":" << audit.shadow.size()
            << ",\"additional_A5_count\":" << audit.additions
            << ",\"candidate_count\":" << audit.candidates.size()
            << ",\"raw_unordered_addition_count\":"
            << Choose(static_cast<int>(audit.candidates.size()), audit.additions)
            << ",\"base_small_transversal_count\":" << audit.small.size()
            << ",\"initial_order_five_cover_count\":" << audit.covers.size()
            << ",\"search_node_count\":" << audit.nodes
            << ",\"canonical_state_counts\":[";
  for (int depth = 0; depth <= audit.additions; ++depth) {
    if (depth) std::cout << ',';
    std::cout << audit.seen[depth].size();
  }
  std::cout << "]"
            << ",\"facet_lower_bound_prune_count\":" << audit.facet_prunes
            << ",\"uncoverable_prune_count\":" << audit.uncoverable_prunes
            << ",\"scheduled_incidence_orbit_count\":"
            << audit.scheduled_orbits
            << ",\"scheduled_blocker_count_histogram\":{";
  bool first = true;
  for (const auto& [size, count] : audit.blocker_histogram) {
    if (!first) std::cout << ',';
    std::cout << '\"' << size << "\":" << count;
    first = false;
  }
  std::cout << "},\"pure_128_facet_target_orbit_count\":"
            << audit.targets.size()
            << ",\"literal_target_count_with_fixed_A4_pair\":"
            << audit.literal_target_count
            << ",\"target_orbit_size_histogram\":{";
  first = true;
  for (const auto& [size, count] : audit.orbit_size_histogram) {
    if (!first) std::cout << ',';
    std::cout << '\"' << size << "\":" << count;
    first = false;
  }
  std::cout << "},\"target_orbit_sizes\":[";
  for (int index = 0;
       index < static_cast<int>(audit.target_orbit_sizes.size()); ++index) {
    if (index) std::cout << ',';
    std::cout << audit.target_orbit_sizes[index];
  }
  std::cout << "],\"pure_target_representatives\":[";
  for (int index = 0; index < static_cast<int>(audit.targets.size()); ++index) {
    if (index) std::cout << ',';
    PrintFamily(audit.targets[index]);
  }
  std::cout << "]}";
}

}  // namespace

int main() {
  try {
    int count = 0;
    if (!(std::cin >> count) || count != 21) {
      throw std::runtime_error("expected 21 positive-A4 representatives");
    }
    std::vector<std::pair<int, std::array<Mask, 2>>> input;
    for (int row = 0; row < count; ++row) {
      int index = -1;
      unsigned first = 0;
      unsigned second = 0;
      if (!(std::cin >> index >> first >> second) || index != row ||
          first > kFull || second > kFull ||
          std::popcount(first) != 4 || std::popcount(second) != 4) {
        throw std::runtime_error("malformed positive-A4 representative");
      }
      input.push_back(
          {index, {static_cast<Mask>(first), static_cast<Mask>(second)}});
    }
    std::vector<Mask> independent_fives;
    for (Mask face = 0; face <= kFull; ++face) {
      if (std::popcount(face) == 5 && Independent(face)) {
        independent_fives.push_back(face);
      }
    }
    if (independent_fives.size() != 1372) {
      throw std::runtime_error("independent-five ground changed");
    }
    const auto graph_automorphisms = GraphAutomorphisms();
    std::vector<Audit> audits;
    std::uint64_t total_targets = 0;
    for (const auto& [index, fours] : input) {
      audits.push_back(Run(index, fours, independent_fives, graph_automorphisms));
      total_targets += audits.back().targets.size();
      std::cerr << "orbit " << index << " targets="
                << audits.back().targets.size() << '\n';
    }
    std::cout
        << "{\"kind\":\"h9_14_record30_a2_all_orbits_explicit_incidence_audit\""
        << ",\"compact_only\":true"
        << ",\"graph_automorphism_count\":" << graph_automorphisms.size()
        << ",\"positive_A4_orbit_count\":" << audits.size()
        << ",\"pure_128_facet_target_orbit_count\":" << total_targets
        << ",\"orbit_audits\":[";
    for (int index = 0; index < static_cast<int>(audits.size()); ++index) {
      if (index) std::cout << ',';
      PrintAudit(audits[index]);
    }
    std::cout
        << "],\"scope\":\"Standalone exact colored-incidence and explicit-stabilizer replay of all 21 compact Record-30 a=2 positive-A4 orbits. No nauty, SAT/SMT, Coxeter-label, or finite-volume inference.\"}\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
