// Lossless pivoted A5-addition/blocker census for every positive-A4 orbit
// in H9/14 Record 31's compact a=2 schedule.
//
// Input on stdin:
//   representative_count
//   orbit_index a4_mask_0 a4_mask_1
//
// Unlike the earlier a=3 scout, this implementation supports more than 64
// simultaneous small-transversal requirements and more than 1024 candidate
// A5 rows.  Dynamic bitsets make coverage exact, and a least-destroyer
// partition visits each final addition family on exactly one branch without
// a potentially enormous chosen-set memo table.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include <nauty.h>

namespace {

using Mask = std::uint16_t;
constexpr int kN = 14;
constexpr Mask kFull = (Mask{1} << kN) - 1;
constexpr int kExpectedOrbitCount = 49;
constexpr int kScheduledFives = 21;
constexpr int kScheduledFacets = 118;
constexpr int kMaxAdditions = 6;
constexpr std::array<Mask, 4> kEdges = {
    (Mask{1} << 0) | (Mask{1} << 11),
    (Mask{1} << 1) | (Mask{1} << 12),
    (Mask{1} << 0) | (Mask{1} << 13),
    (Mask{1} << 1) | (Mask{1} << 13),
};
constexpr std::array<int, kN> kBaseG3 = {
    2, 2, 7, 7, 7, 7, 7, 7, 7, 7, 7, 5, 5, 4,
};
constexpr std::array<int, kN> kBaseG4 = {
    2, 2, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 6, 5,
};

using Words = std::vector<std::uint64_t>;

struct LocalConstraint {
  Mask subset = 0;
  int lower = 0;
  int upper = 0;
  int previous = 0;
  int before = 0;
  std::vector<int> f_prefix;
  std::vector<int> h_prefix;
  std::vector<int> g_prefix;
};

struct LinkPrefix {
  std::vector<int> f;
  std::vector<int> h;
  std::vector<int> g;
};

struct Audit {
  int index = -1;
  std::array<Mask, 2> fours{};
  std::vector<Mask> shadow;
  std::vector<Mask> candidates;
  std::vector<Mask> small;
  std::vector<Words> coverage;
  std::vector<std::vector<int>> coverers;
  std::vector<Mask> five_covers;
  Words base_five_cover_words;
  std::vector<Words> candidate_five_hit_words;
  std::vector<LocalConstraint> local_constraints;
  std::array<int, 1 << kN> local_index{};
  std::vector<std::array<int, 25>> candidate_local_indices;
  std::array<int, kN> lower{};
  std::array<int, kN> upper{};
  int additions = 0;
  bool fixed_local_failure = false;
  bool matrix_quotient = false;
  std::uint64_t nodes = 0;
  std::array<std::uint64_t, kMaxAdditions + 1> depth_nodes{};
  std::uint64_t dead_uncoverable = 0;
  std::uint64_t dead_five_count = 0;
  std::uint64_t dead_slot_bound = 0;
  std::uint64_t no_small = 0;
  std::map<int, std::uint64_t> order_five_histogram;
  std::map<int, std::uint64_t> target_blocker_size_histogram;
  std::map<std::vector<std::pair<int, int>>, std::uint64_t>
      pure_ridge_histogram;
  std::uint64_t pure_targets = 0;
  std::uint64_t ridge_valid_literal_targets = 0;
  std::map<std::string, std::pair<std::vector<Mask>, std::uint64_t>>
      ridge_valid_orbits;
  std::array<std::unordered_set<std::string>, kMaxAdditions + 1>
      quotient_seen;
  std::array<std::uint64_t, kMaxAdditions + 1> quotient_seen_final_counts{};
  std::uint64_t quotient_isomorphic_prunes = 0;
};

struct SearchState {
  std::vector<int> chosen;
  Words forbidden;
  Words covered;
  std::array<std::uint8_t, 420> local_incidence{};
  Words surviving_five_covers;
};

long long Binomial(int n, int k) {
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
           Binomial(top + 1, lower) <= remaining) {
      ++top;
    }
    remaining -= static_cast<int>(Binomial(top, lower));
    answer += static_cast<int>(Binomial(top + 1, lower + 1));
    next_top = top;
  }
  if (remaining != 0) throw std::runtime_error("Macaulay expansion failed");
  return answer;
}

bool Independent(Mask value) {
  for (const Mask edge : kEdges) {
    if ((value & edge) == edge) return false;
  }
  return true;
}

LinkPrefix PrefixData(Mask linked, const std::array<Mask, 2>& fours) {
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
      bool allowed = Independent(face);
      for (const Mask four : fours) {
        if ((face & four) == four) allowed = false;
      }
      f[size] += allowed;
    }
    if (!extension) break;
  }
  std::vector<int> h(top + 1);
  for (int index = 0; index <= top; ++index) {
    for (int order = 0; order <= index; ++order) {
      h[index] += (((index - order) & 1) ? -1 : 1) *
          static_cast<int>(Binomial(dimension - order,
                                    dimension - index)) * f[order];
    }
  }
  std::vector<int> g(top + 1);
  g[0] = h[0];
  for (int index = 1; index <= top; ++index) {
    g[index] = h[index] - h[index - 1];
  }
  return {f, h, g};
}

bool ContainsFour(Mask five, const std::array<Mask, 2>& fours) {
  for (const Mask four : fours) {
    if ((five & four) == four) return true;
  }
  return false;
}

bool HitsAll(Mask value, const std::vector<Mask>& clutter) {
  for (const Mask face : clutter) {
    if ((value & face) == 0) return false;
  }
  return true;
}

std::vector<Mask> IndependentFives() {
  std::vector<Mask> answer;
  for (Mask value = 0; value <= kFull; ++value) {
    if (std::popcount(value) == 5 && Independent(value)) answer.push_back(value);
  }
  if (answer.size() != 1296) {
    throw std::runtime_error("Record-31 independent-five count changed");
  }
  return answer;
}

std::vector<Mask> GraphFiveCovers() {
  std::vector<Mask> answer;
  for (Mask value = 0; value <= kFull; ++value) {
    if (std::popcount(value) != 5) continue;
    bool cover = true;
    for (const Mask edge : kEdges) cover = cover && ((value & edge) != 0);
    if (cover) answer.push_back(value);
  }
  return answer;
}

int WordCount(const Words& words) {
  int answer = 0;
  for (const std::uint64_t word : words) answer += std::popcount(word);
  return answer;
}

bool IsSet(const Words& words, int index) {
  return (words[index / 64] >> (index % 64)) & 1u;
}

void Set(Words* words, int index) {
  (*words)[index / 64] |= std::uint64_t{1} << (index % 64);
}

Words Meet(const Words& left, const Words& right) {
  Words answer(left.size());
  for (int word = 0; word < static_cast<int>(left.size()); ++word) {
    answer[word] = left[word] & right[word];
  }
  return answer;
}

int MeetCount(const Words& left, const Words& right) {
  int answer = 0;
  for (int word = 0; word < static_cast<int>(left.size()); ++word) {
    answer += std::popcount(left[word] & right[word]);
  }
  return answer;
}

void OrInto(Words* target, const Words& source) {
  for (int word = 0; word < static_cast<int>(target->size()); ++word) {
    (*target)[word] |= source[word];
  }
}

bool IsFull(const Words& covered, int bit_count) {
  const int whole_words = bit_count / 64;
  for (int word = 0; word < whole_words; ++word) {
    if (covered[word] != ~std::uint64_t{0}) return false;
  }
  const int remainder = bit_count % 64;
  if (remainder == 0) return true;
  const std::uint64_t tail = (std::uint64_t{1} << remainder) - 1;
  return (covered[whole_words] & tail) == tail;
}

bool Compatible(const Audit& audit, int candidate,
                const std::array<std::uint8_t, 420>& incidence) {
  for (const int local : audit.candidate_local_indices[candidate]) {
    if (incidence[local] >= audit.local_constraints[local].upper) {
      return false;
    }
  }
  return true;
}

std::vector<Mask> FullMinimalTransversals(const std::vector<Mask>& clutter) {
  std::vector<Mask> answer;
  for (Mask value = 1; value <= kFull; ++value) {
    if (!HitsAll(value, clutter)) continue;
    bool minimal = true;
    for (int vertex = 0; vertex < kN; ++vertex) {
      if ((value >> vertex) & 1u) {
        minimal = minimal && !HitsAll(value ^ (Mask{1} << vertex), clutter);
      }
    }
    if (minimal) answer.push_back(value);
  }
  return answer;
}

std::string CanonicalFamilyCode(const std::vector<Mask>& family) {
  if (family.size() < 2 || family.size() > 8) {
    throw std::runtime_error("Record-31 a=2 family shape changed");
  }
  const int edge_start = kN;
  const int four_start = edge_start + kEdges.size();
  const int five_start = four_start + 2;
  const int node_count = five_start + static_cast<int>(family.size()) - 2;
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
  for (int row = 0; row < static_cast<int>(family.size()); ++row) {
    const int node = row < 2 ? four_start + row : five_start + row - 2;
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (family[row] & (Mask{1} << vertex)) {
        ADDONEEDGE(input.data(), node, vertex, words);
      }
    }
  }
  std::vector<int> lab(node_count), ptn(node_count), orbits(node_count);
  int position = 0;
  auto add_cell = [&](int first, int last) {
    for (int value = first; value < last; ++value) {
      lab[position] = value;
      ptn[position++] = 1;
    }
    if (last > first) ptn[position - 1] = 0;
  };
  add_cell(0, kN);
  add_cell(edge_start, four_start);
  add_cell(four_start, five_start);
  add_cell(five_start, node_count);
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
    throw std::runtime_error("nauty family canonicalization failed");
  }
  return std::string(reinterpret_cast<const char*>(canonical.data()),
                     canonical.size() * sizeof(graph));
}

std::string CanonicalPartialCode(const Audit& audit,
                                 const SearchState& state) {
  std::vector<Mask> family(audit.fours.begin(), audit.fours.end());
  for (const int candidate : state.chosen) {
    family.push_back(audit.candidates[candidate]);
  }
  std::sort(family.begin(), family.begin() + 2);
  std::sort(family.begin() + 2, family.end());
  return CanonicalFamilyCode(family);
}

Mask FlipActive(Mask face) {
  Mask result = face;
  const bool zero = face & (Mask{1} << 0);
  const bool one = face & (Mask{1} << 1);
  const bool eleven = face & (Mask{1} << 11);
  const bool twelve = face & (Mask{1} << 12);
  result &= ~((Mask{1} << 0) | (Mask{1} << 1) |
              (Mask{1} << 11) | (Mask{1} << 12));
  if (zero) result |= Mask{1} << 1;
  if (one) result |= Mask{1} << 0;
  if (eleven) result |= Mask{1} << 12;
  if (twelve) result |= Mask{1} << 11;
  return result;
}

int ActivePattern(Mask face) {
  constexpr std::array<int, 5> active = {0, 1, 11, 12, 13};
  int result = 0;
  for (int index = 0; index < static_cast<int>(active.size()); ++index) {
    result |= ((face >> active[index]) & 1u) << index;
  }
  return result;
}

std::vector<std::vector<int>> EqualPatternOrders(
    const std::vector<Mask>& rows, int begin, int end) {
  std::map<int, std::vector<int>> buckets;
  for (int index = begin; index < end; ++index) {
    buckets[ActivePattern(rows[index])].push_back(index);
  }
  std::vector<std::vector<int>> orders(1);
  for (auto& [pattern, bucket] : buckets) {
    (void)pattern;
    std::sort(bucket.begin(), bucket.end());
    std::vector<std::vector<int>> next;
    do {
      for (const auto& prefix : orders) {
        auto value = prefix;
        value.insert(value.end(), bucket.begin(), bucket.end());
        next.push_back(std::move(value));
      }
    } while (std::next_permutation(bucket.begin(), bucket.end()));
    orders = std::move(next);
  }
  return orders;
}

std::string CanonicalMatrixCode(const Audit& audit,
                                const SearchState& state) {
  std::vector<Mask> original(audit.fours.begin(), audit.fours.end());
  for (const int candidate : state.chosen) {
    original.push_back(audit.candidates[candidate]);
  }
  std::string best;
  bool first = true;
  for (int flip = 0; flip < 2; ++flip) {
    std::vector<Mask> rows = original;
    if (flip) {
      for (Mask& row : rows) row = FlipActive(row);
    }
    const auto four_orders = EqualPatternOrders(rows, 0, 2);
    const auto five_orders = EqualPatternOrders(rows, 2, rows.size());
    for (const auto& four_order : four_orders) {
      for (const auto& five_order : five_orders) {
        std::vector<int> order = four_order;
        order.insert(order.end(), five_order.begin(), five_order.end());
        std::string code;
        code.reserve(order.size() + 9);
        for (const int row : order) {
          code.push_back(static_cast<char>(ActivePattern(rows[row])));
        }
        std::array<int, 9> columns{};
        for (int isolate = 0; isolate < 9; ++isolate) {
          const int vertex = isolate + 2;
          for (int position = 0; position < static_cast<int>(order.size());
               ++position) {
            columns[isolate] |=
                ((rows[order[position]] >> vertex) & 1u) << position;
          }
        }
        std::sort(columns.begin(), columns.end());
        for (const int column : columns) {
          code.push_back(static_cast<char>(column));
        }
        if (first || code < best) {
          best = std::move(code);
          first = false;
        }
      }
    }
  }
  if (first) throw std::runtime_error("matrix canonicalization produced no key");
  return best;
}

void Finish(Audit* audit, const SearchState& state,
            const std::vector<Mask>& base) {
  ++audit->no_small;
  const int order_five = WordCount(state.surviving_five_covers);
  ++audit->order_five_histogram[order_five];
  if (order_five != kScheduledFacets) return;

  std::vector<Mask> clutter = base;
  for (const int candidate : state.chosen) {
    clutter.push_back(audit->candidates[candidate]);
  }
  const auto blocker = FullMinimalTransversals(clutter);
  ++audit->target_blocker_size_histogram[blocker.size()];
  bool pure = blocker.size() == kScheduledFacets;
  for (const Mask face : blocker) pure = pure && std::popcount(face) == 5;
  if (!pure) return;

  ++audit->pure_targets;
  std::map<Mask, int> ridge_degrees;
  for (const Mask complement : blocker) {
    Mask outside = kFull ^ complement;
    while (outside) {
      const Mask bit = outside & static_cast<Mask>(-outside);
      outside ^= bit;
      ++ridge_degrees[complement | bit];
    }
  }
  std::map<int, int> degree_histogram;
  for (const auto& [ridge, degree] : ridge_degrees) {
    (void)ridge;
    ++degree_histogram[degree];
  }
  const std::vector<std::pair<int, int>> signature(
      degree_histogram.begin(), degree_histogram.end());
  ++audit->pure_ridge_histogram[signature];
  if (degree_histogram != std::map<int, int>{{2, 531}}) return;

  std::vector<Mask> family(audit->fours.begin(), audit->fours.end());
  for (const int candidate : state.chosen) {
    family.push_back(audit->candidates[candidate]);
  }
  std::sort(family.begin(), family.begin() + 2);
  std::sort(family.begin() + 2, family.end());
  ++audit->ridge_valid_literal_targets;
  const std::string code = audit->matrix_quotient
      ? CanonicalMatrixCode(*audit, state)
      : CanonicalFamilyCode(family);
  auto [iterator, inserted] = audit->ridge_valid_orbits.try_emplace(
      code, std::make_pair(family, std::uint64_t{0}));
  ++iterator->second.second;
}

SearchState AddCandidate(const Audit& audit, const SearchState& state,
                         int candidate, Words forbidden) {
  SearchState next = state;
  next.chosen.push_back(candidate);
  next.forbidden = std::move(forbidden);
  OrInto(&next.covered, audit.coverage[candidate]);
  for (const int local : audit.candidate_local_indices[candidate]) {
    ++next.local_incidence[local];
  }
  next.surviving_five_covers =
      Meet(state.surviving_five_covers,
           audit.candidate_five_hit_words[candidate]);
  return next;
}

void Search(Audit* audit, const SearchState& state,
            const std::vector<Mask>& base);

void Fill(Audit* audit, const SearchState& state, int first_candidate,
          const std::vector<Mask>& base) {
  ++audit->nodes;
  if (audit->nodes % 1000000 == 0) {
    std::cerr << "orbit " << audit->index << " progress nodes=" << audit->nodes
              << " pure=" << audit->pure_targets
              << " ridge=" << audit->ridge_valid_literal_targets
              << " q=" << audit->ridge_valid_orbits.size() << '\n';
  }
  ++audit->depth_nodes[state.chosen.size()];
  if (WordCount(state.surviving_five_covers) < kScheduledFacets) {
    ++audit->dead_five_count;
    return;
  }
  if (static_cast<int>(state.chosen.size()) == audit->additions) {
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (state.local_incidence[vertex] < audit->lower[vertex] ||
          state.local_incidence[vertex] > audit->upper[vertex]) {
        return;
      }
    }
    for (int local = 0;
         local < static_cast<int>(audit->local_constraints.size()); ++local) {
      if (state.local_incidence[local] < audit->local_constraints[local].lower ||
          state.local_incidence[local] > audit->local_constraints[local].upper) {
        return;
      }
    }
    Finish(audit, state, base);
    return;
  }
  for (int candidate = first_candidate;
       candidate < static_cast<int>(audit->candidates.size()); ++candidate) {
    if (IsSet(state.forbidden, candidate) ||
        !Compatible(*audit, candidate, state.local_incidence)) {
      continue;
    }
    const Words surviving =
        Meet(state.surviving_five_covers,
             audit->candidate_five_hit_words[candidate]);
    if (WordCount(surviving) < kScheduledFacets) continue;
    Words forbidden = state.forbidden;
    Set(&forbidden, candidate);
    auto next = AddCandidate(*audit, state, candidate, std::move(forbidden));
    next.surviving_five_covers = surviving;
    Fill(audit, next, candidate + 1, base);
  }
}

void SearchQuotient(Audit* audit, const SearchState& state,
                    const std::vector<Mask>& base) {
  ++audit->nodes;
  ++audit->depth_nodes[state.chosen.size()];
  if (audit->nodes % 1000000 == 0) {
    std::cerr << "orbit " << audit->index << " quotient-progress nodes="
              << audit->nodes << " seen=";
    for (int depth = 0; depth <= audit->additions; ++depth) {
      std::cerr << audit->quotient_seen[depth].size() << ',';
    }
    std::cerr << " q=" << audit->ridge_valid_orbits.size() << '\n';
  }
  const int depth = static_cast<int>(state.chosen.size());
  const std::string canonical = audit->matrix_quotient
      ? CanonicalMatrixCode(*audit, state)
      : CanonicalPartialCode(*audit, state);
  if (!audit->quotient_seen[depth].insert(canonical).second) {
    ++audit->quotient_isomorphic_prunes;
    return;
  }
  const int current_facets = WordCount(state.surviving_five_covers);
  if (current_facets < kScheduledFacets) {
    ++audit->dead_five_count;
    return;
  }
  if (depth == audit->additions) {
    if (IsFull(state.covered, audit->small.size())) Finish(audit, state, base);
    return;
  }

  const int remaining = audit->additions - depth;
  std::vector<std::uint8_t> available(audit->candidates.size(), 0);
  std::vector<int> losses;
  std::vector<int> new_coverage_counts;
  int uncovered_count = 0;
  for (int requirement = 0;
       requirement < static_cast<int>(audit->small.size()); ++requirement) {
    uncovered_count +=
        !((state.covered[requirement / 64] >> (requirement % 64)) & 1u);
  }
  for (int candidate = 0;
       candidate < static_cast<int>(audit->candidates.size()); ++candidate) {
    if (IsSet(state.forbidden, candidate) ||
        !Compatible(*audit, candidate, state.local_incidence)) {
      continue;
    }
    const int next_facets = MeetCount(
        state.surviving_five_covers,
        audit->candidate_five_hit_words[candidate]);
    if (next_facets < kScheduledFacets) continue;
    available[candidate] = 1;
    losses.push_back(current_facets - next_facets);
    int newly_covered = 0;
    for (int word = 0; word < static_cast<int>(state.covered.size()); ++word) {
      newly_covered += std::popcount(
          audit->coverage[candidate][word] & ~state.covered[word]);
    }
    new_coverage_counts.push_back(newly_covered);
  }
  if (losses.empty()) {
    ++audit->dead_uncoverable;
    return;
  }
  std::ranges::sort(losses, std::greater<>());
  int maximum_loss = 0;
  for (int index = 0;
       index < std::min(remaining, static_cast<int>(losses.size())); ++index) {
    maximum_loss += losses[index];
  }
  if (current_facets - maximum_loss > kScheduledFacets) {
    ++audit->dead_slot_bound;
    return;
  }
  std::ranges::sort(new_coverage_counts, std::greater<>());
  int maximum_new_coverage = 0;
  for (int index = 0;
       index < std::min(remaining,
                        static_cast<int>(new_coverage_counts.size())); ++index) {
    maximum_new_coverage += new_coverage_counts[index];
  }
  if (maximum_new_coverage < uncovered_count) {
    ++audit->dead_slot_bound;
    return;
  }

  std::vector<int> branches;
  if (IsFull(state.covered, audit->small.size())) {
    for (int candidate = 0;
         candidate < static_cast<int>(audit->candidates.size()); ++candidate) {
      if (available[candidate]) branches.push_back(candidate);
    }
  } else {
    for (int requirement = 0;
         requirement < static_cast<int>(audit->small.size()); ++requirement) {
      if ((state.covered[requirement / 64] >> (requirement % 64)) & 1u) continue;
      std::vector<int> feasible;
      for (const int candidate : audit->coverers[requirement]) {
        if (available[candidate]) feasible.push_back(candidate);
      }
      if (feasible.empty()) {
        ++audit->dead_uncoverable;
        return;
      }
      if (branches.empty() || feasible.size() < branches.size() ||
          (audit->matrix_quotient && feasible.size() == branches.size())) {
        branches = std::move(feasible);
      }
    }
  }
  if (audit->matrix_quotient) std::reverse(branches.begin(), branches.end());
  for (const int candidate : branches) {
    Words forbidden = state.forbidden;
    Set(&forbidden, candidate);
    SearchQuotient(
        audit, AddCandidate(*audit, state, candidate, std::move(forbidden)),
        base);
  }
}

void Search(Audit* audit, const SearchState& state,
            const std::vector<Mask>& base) {
  ++audit->nodes;
  if (audit->nodes % 1000000 == 0) {
    std::cerr << "orbit " << audit->index << " progress nodes=" << audit->nodes
              << " pure=" << audit->pure_targets
              << " ridge=" << audit->ridge_valid_literal_targets
              << " q=" << audit->ridge_valid_orbits.size() << '\n';
  }
  ++audit->depth_nodes[state.chosen.size()];
  if (WordCount(state.surviving_five_covers) < kScheduledFacets) {
    ++audit->dead_five_count;
    return;
  }
  if (static_cast<int>(state.chosen.size()) == audit->additions) {
    if (IsFull(state.covered, audit->small.size())) Finish(audit, state, base);
    return;
  }
  if (IsFull(state.covered, audit->small.size())) {
    Fill(audit, state, 0, base);
    return;
  }

  const int remaining =
      audit->additions - static_cast<int>(state.chosen.size());
  const int current_facets = WordCount(state.surviving_five_covers);
  std::vector<std::uint8_t> available(audit->candidates.size(), 0);
  std::vector<int> losses;
  std::vector<int> new_coverage_counts;
  losses.reserve(audit->candidates.size());
  new_coverage_counts.reserve(audit->candidates.size());
  int uncovered_count = 0;
  for (int requirement = 0;
       requirement < static_cast<int>(audit->small.size()); ++requirement) {
    uncovered_count +=
        !((state.covered[requirement / 64] >> (requirement % 64)) & 1u);
  }
  for (int candidate = 0;
       candidate < static_cast<int>(audit->candidates.size()); ++candidate) {
    if (IsSet(state.forbidden, candidate) ||
        !Compatible(*audit, candidate, state.local_incidence)) {
      continue;
    }
    const int next_facets = MeetCount(
        state.surviving_five_covers,
        audit->candidate_five_hit_words[candidate]);
    if (next_facets < kScheduledFacets) continue;
    available[candidate] = 1;
    losses.push_back(current_facets - next_facets);
    int newly_covered = 0;
    for (int word = 0; word < static_cast<int>(state.covered.size()); ++word) {
      newly_covered += std::popcount(
          audit->coverage[candidate][word] & ~state.covered[word]);
    }
    new_coverage_counts.push_back(newly_covered);
  }
  if (losses.empty()) {
    ++audit->dead_uncoverable;
    return;
  }
  std::ranges::sort(losses, std::greater<>());
  int maximum_loss = 0;
  for (int index = 0;
       index < std::min(remaining, static_cast<int>(losses.size())); ++index) {
    maximum_loss += losses[index];
  }
  if (current_facets - maximum_loss > kScheduledFacets) {
    ++audit->dead_slot_bound;
    return;
  }
  std::ranges::sort(new_coverage_counts, std::greater<>());
  int maximum_new_coverage = 0;
  for (int index = 0;
       index < std::min(remaining,
                        static_cast<int>(new_coverage_counts.size())); ++index) {
    maximum_new_coverage += new_coverage_counts[index];
  }
  if (maximum_new_coverage < uncovered_count) {
    ++audit->dead_slot_bound;
    return;
  }

  std::vector<int> best;
  Words feasible_union = state.covered;
  for (int requirement = 0;
       requirement < static_cast<int>(audit->small.size()); ++requirement) {
    if ((state.covered[requirement / 64] >> (requirement % 64)) & 1u) continue;
    std::vector<int> feasible;
    for (const int candidate : audit->coverers[requirement]) {
      if (!available[candidate]) continue;
      feasible.push_back(candidate);
      OrInto(&feasible_union, audit->coverage[candidate]);
    }
    if (feasible.empty()) {
      ++audit->dead_uncoverable;
      return;
    }
    if (best.empty() || feasible.size() < best.size()) best = std::move(feasible);
  }
  if (!IsFull(feasible_union, audit->small.size())) {
    ++audit->dead_uncoverable;
    return;
  }
  // Disjoint least-destroyer partition: branch i asserts that candidate i is
  // the first currently feasible destroyer of the chosen pivot occurring in
  // the final family, while all earlier candidates are permanently forbidden.
  Words cumulative_forbidden = state.forbidden;
  for (const int candidate : best) {
    Words branch_forbidden = cumulative_forbidden;
    Set(&branch_forbidden, candidate);
    Search(audit,
           AddCandidate(*audit, state, candidate, std::move(branch_forbidden)),
           base);
    Set(&cumulative_forbidden, candidate);
  }
}

Audit Run(int index, const std::array<Mask, 2>& fours,
          const std::vector<Mask>& fives,
          const std::vector<Mask>& five_covers, bool quotient_mode,
          bool matrix_quotient) {
  Audit audit;
  audit.index = index;
  audit.matrix_quotient = matrix_quotient;
  audit.fours = fours;
  audit.five_covers = five_covers;
  for (const Mask five : fives) {
    if (ContainsFour(five, fours)) audit.shadow.push_back(five);
  }
  audit.additions = kScheduledFives - audit.shadow.size();
  if (audit.additions < 3 || audit.additions > kMaxAdditions) {
    throw std::runtime_error("Record-31 a=2 addition count changed");
  }

  audit.local_index.fill(-1);
  for (int size = 1; size <= 3; ++size) {
    for (Mask subset = 1; subset <= kFull; ++subset) {
      if (std::popcount(subset) != size || !Independent(subset)) continue;
      const auto prefix = PrefixData(subset, fours);
      const int top = 5 - size;
      const int previous = prefix.g[top - 1];
      const int before = prefix.g[top];
      LocalConstraint constraint;
      constraint.subset = subset;
      constraint.previous = previous;
      constraint.before = before;
      constraint.f_prefix = prefix.f;
      constraint.h_prefix = prefix.h;
      constraint.g_prefix = prefix.g;
      constraint.lower = previous >= 0
          ? std::max(0, before - MacaulayRaise(previous, top - 1))
          : 0;
      constraint.upper = before - 1;
      audit.local_index[subset] = audit.local_constraints.size();
      audit.local_constraints.push_back(constraint);
      if (constraint.lower > constraint.upper) audit.fixed_local_failure = true;
    }
  }
  if (audit.local_constraints.size() != 420) {
    throw std::runtime_error("Record-31 local face ground changed");
  }

  for (int vertex = 0; vertex < kN; ++vertex) {
    int a4_incidence = 0;
    int forced_incidence = 0;
    for (const Mask four : fours) a4_incidence += (four >> vertex) & 1u;
    for (const Mask five : audit.shadow) {
      forced_incidence += (five >> vertex) & 1u;
    }
    const int previous = kBaseG3[vertex] - a4_incidence;
    const int before = kBaseG4[vertex] + 6 * a4_incidence - forced_incidence;
    if (previous < 1 || before < 1) {
      throw std::runtime_error("positive orbit escaped strict local gate");
    }
    audit.lower[vertex] =
        std::max(0, before - MacaulayRaise(previous, 3));
    audit.upper[vertex] = before - 1;
    const auto& local = audit.local_constraints[vertex];
    if (local.subset != (Mask{1} << vertex) ||
        local.previous != previous || local.before != before ||
        local.lower != audit.lower[vertex] ||
        local.upper != audit.upper[vertex]) {
      throw std::runtime_error("vertex and full local constraints disagree");
    }
  }

  for (const Mask five : fives) {
    if (ContainsFour(five, fours)) continue;
    bool allowed = true;
    std::array<int, 25> local_indices{};
    int position = 0;
    for (Mask subset = five; subset;
         subset = (subset - 1) & five) {
      const int size = std::popcount(subset);
      if (size < 1 || size > 3) continue;
      const int local = audit.local_index[subset];
      if (local < 0) throw std::runtime_error("candidate local support missing");
      local_indices[position++] = local;
      if (audit.local_constraints[local].upper < 1) allowed = false;
    }
    if (position != 25) throw std::runtime_error("A5 local support count changed");
    if (allowed) {
      audit.candidates.push_back(five);
      audit.candidate_local_indices.push_back(local_indices);
    }
  }

  std::vector<Mask> base(kEdges.begin(), kEdges.end());
  base.insert(base.end(), fours.begin(), fours.end());
  base.insert(base.end(), audit.shadow.begin(), audit.shadow.end());
  for (Mask value = 1; value <= kFull; ++value) {
    if (std::popcount(value) <= 4 && HitsAll(value, base)) {
      audit.small.push_back(value);
    }
  }
  const int coverage_words = (audit.small.size() + 63) / 64;
  audit.coverers.assign(audit.small.size(), {});
  for (int candidate = 0;
       candidate < static_cast<int>(audit.candidates.size()); ++candidate) {
    Words words(coverage_words, 0);
    for (int requirement = 0;
         requirement < static_cast<int>(audit.small.size()); ++requirement) {
      if ((audit.candidates[candidate] & audit.small[requirement]) == 0) {
        words[requirement / 64] |=
            std::uint64_t{1} << (requirement % 64);
        audit.coverers[requirement].push_back(candidate);
      }
    }
    audit.coverage.push_back(std::move(words));
  }

  const int five_words = (five_covers.size() + 63) / 64;
  audit.base_five_cover_words.assign(five_words, 0);
  for (int cover = 0; cover < static_cast<int>(five_covers.size()); ++cover) {
    bool hits = true;
    for (const Mask four : fours) hits = hits && ((five_covers[cover] & four) != 0);
    if (hits) {
      audit.base_five_cover_words[cover / 64] |=
          std::uint64_t{1} << (cover % 64);
    }
  }
  audit.candidate_five_hit_words.assign(
      audit.candidates.size(), Words(five_words, 0));
  for (int candidate = 0;
       candidate < static_cast<int>(audit.candidates.size()); ++candidate) {
    for (int cover = 0; cover < static_cast<int>(five_covers.size()); ++cover) {
      if (audit.candidates[candidate] & five_covers[cover]) {
        audit.candidate_five_hit_words[candidate][cover / 64] |=
            std::uint64_t{1} << (cover % 64);
      }
    }
  }

  SearchState initial;
  initial.forbidden.assign((audit.candidates.size() + 63) / 64, 0);
  initial.covered.assign(coverage_words, 0);
  initial.surviving_five_covers = audit.base_five_cover_words;
  if (!audit.fixed_local_failure) {
    if (quotient_mode) {
      SearchQuotient(&audit, initial, base);
    } else {
      Search(&audit, initial, base);
    }
  }
  for (int depth = 0; depth <= audit.additions; ++depth) {
    audit.quotient_seen_final_counts[depth] = audit.quotient_seen[depth].size();
    audit.quotient_seen[depth].clear();
    audit.quotient_seen[depth].rehash(0);
  }
  return audit;
}

void PrintMask(Mask value) {
  std::cout << '[';
  bool first = true;
  for (int vertex = 0; vertex < kN; ++vertex) {
    if ((value >> vertex) & 1u) {
      if (!first) std::cout << ',';
      std::cout << vertex;
      first = false;
    }
  }
  std::cout << ']';
}

void PrintHistogram(const std::map<int, std::uint64_t>& histogram) {
  std::cout << '{';
  bool first = true;
  for (const auto& [key, value] : histogram) {
    if (!first) std::cout << ',';
    std::cout << '\"' << key << "\":" << value;
    first = false;
  }
  std::cout << '}';
}

}  // namespace

int main(int argc, char** argv) {
  try {
    bool quotient_mode = false;
    bool matrix_quotient = false;
    if (argc == 2 && std::string(argv[1]) == "--quotient") {
      quotient_mode = true;
    } else if (argc == 2 && std::string(argv[1]) == "--matrix-quotient") {
      quotient_mode = true;
      matrix_quotient = true;
    } else if (argc != 1) {
      throw std::runtime_error("usage: explore_record31_a2 [--quotient]");
    }
    int count = 0;
    if (!(std::cin >> count) || count <= 0 || count > kExpectedOrbitCount) {
      throw std::runtime_error("expected between one and 49 representatives");
    }
    std::vector<std::pair<int, std::array<Mask, 2>>> input;
    int previous = -1;
    for (int row = 0; row < count; ++row) {
      int index = -1;
      unsigned left = 0, right = 0;
      if (!(std::cin >> index >> left >> right) || index <= previous ||
          index >= kExpectedOrbitCount || left > kFull || right > kFull ||
          std::popcount(left) != 4 || std::popcount(right) != 4) {
        throw std::runtime_error("malformed positive-A4 representative input");
      }
      previous = index;
      input.push_back({index, {static_cast<Mask>(left), static_cast<Mask>(right)}});
    }

    const auto fives = IndependentFives();
    const auto five_covers = GraphFiveCovers();
    std::vector<Audit> audits;
    std::uint64_t total_targets = 0;
    for (const auto& [index, fours] : input) {
      audits.push_back(Run(index, fours, fives, five_covers, quotient_mode,
                           matrix_quotient));
      total_targets += audits.back().pure_targets;
      std::cerr << "orbit " << index << " add=" << audits.back().additions
                << " nodes=" << audits.back().nodes
                << " no-small=" << audits.back().no_small
                << " pure=" << audits.back().pure_targets
                << " ridge-literal=" << audits.back().ridge_valid_literal_targets
                << " ridge-orbits=" << audits.back().ridge_valid_orbits.size()
                << '\n';
    }

    std::cout << "{\"kind\":\"h9_14_record31_a2_pivoted_addition_scout\",";
    std::cout << "\"compact_only\":true,\"enumeration_mode\":\""
              << (matrix_quotient ? "direct_binary_matrix_partial_quotient" :
                  quotient_mode ? "nauty_partial_quotient" :
                                  "literal_least_destroyer_partition")
              << "\",\"orbit_count\":" << audits.size()
              << ",\"orbits\":[";
    for (int row = 0; row < static_cast<int>(audits.size()); ++row) {
      if (row) std::cout << ',';
      const auto& audit = audits[row];
      std::cout << "{\"orbit_index\":" << audit.index
                << ",\"positive_A4_masks\":[";
      PrintMask(audit.fours[0]);
      std::cout << ',';
      PrintMask(audit.fours[1]);
      std::cout << "]"
                << ",\"forced_A5_shadow_count\":" << audit.shadow.size()
                << ",\"additional_A5_count\":" << audit.additions
                << ",\"candidate_count\":" << audit.candidates.size()
                << ",\"fixed_local_failure\":"
                << (audit.fixed_local_failure ? "true" : "false")
                << ",\"fixed_local_failure_witnesses\":[";
      bool first_local_witness = true;
      for (const auto& constraint : audit.local_constraints) {
        if (constraint.lower <= constraint.upper) continue;
        if (!first_local_witness) std::cout << ',';
        std::cout << "{\"subset\":";
        PrintMask(constraint.subset);
        std::cout << ",\"link_dimension\":"
                  << 9 - std::popcount(constraint.subset)
                  << ",\"f_prefix\":[";
        for (int entry = 0;
             entry < static_cast<int>(constraint.f_prefix.size()); ++entry) {
          if (entry) std::cout << ',';
          std::cout << constraint.f_prefix[entry];
        }
        std::cout << "],\"h_prefix\":[";
        for (int entry = 0;
             entry < static_cast<int>(constraint.h_prefix.size()); ++entry) {
          if (entry) std::cout << ',';
          std::cout << constraint.h_prefix[entry];
        }
        std::cout << "],\"g_prefix\":[";
        for (int entry = 0;
             entry < static_cast<int>(constraint.g_prefix.size()); ++entry) {
          if (entry) std::cout << ',';
          std::cout << constraint.g_prefix[entry];
        }
        std::cout << "],\"additional_A5_incidence_interval\":["
                  << constraint.lower << ',' << constraint.upper << "]} ";
        first_local_witness = false;
      }
      std::cout << ']'
                << ",\"base_small_transversal_count\":" << audit.small.size()
                << ",\"search_node_count\":" << audit.nodes
                << ",\"partition_depth_node_counts\":[";
      for (int depth = 0; depth <= audit.additions; ++depth) {
        if (depth) std::cout << ',';
        std::cout << audit.depth_nodes[depth];
      }
      std::cout << "]"
                << ",\"dead_uncoverable_count\":" << audit.dead_uncoverable
                << ",\"dead_five_count\":" << audit.dead_five_count
                << ",\"dead_slot_bound_count\":" << audit.dead_slot_bound
                << ",\"quotient_isomorphic_prune_count\":"
                << audit.quotient_isomorphic_prunes
                << ",\"quotient_seen_counts\":[";
      for (int depth = 0; depth <= audit.additions; ++depth) {
        if (depth) std::cout << ',';
        std::cout << audit.quotient_seen_final_counts[depth];
      }
      std::cout << ']'
                << ",\"no_small_transversal_count\":" << audit.no_small
                << ",\"order_five_transversal_count_histogram\":";
      PrintHistogram(audit.order_five_histogram);
      std::cout << ",\"scheduled_blocker_size_histogram\":";
      PrintHistogram(audit.target_blocker_size_histogram);
      std::cout << ",\"pure_target_ridge_histogram\":{";
      bool first_ridge_signature = true;
      for (const auto& [signature, multiplicity] :
           audit.pure_ridge_histogram) {
        if (!first_ridge_signature) std::cout << ',';
        std::cout << '\"';
        for (int entry = 0; entry < static_cast<int>(signature.size()); ++entry) {
          if (entry) std::cout << ';';
          std::cout << signature[entry].first << ':' << signature[entry].second;
        }
        std::cout << "\":" << multiplicity;
        first_ridge_signature = false;
      }
      std::cout << '}'
                << ",\"ridge_valid_pure_target_count\":"
                << audit.ridge_valid_literal_targets
                << ",\"ridge_valid_target_orbit_count\":"
                << audit.ridge_valid_orbits.size();
      std::cout << ",\"pure_118_facet_target_count\":"
                << audit.pure_targets
                << ",\"ridge_valid_target_minimal_nonfaces\":[";
      bool first_target = true;
      for (const auto& [code, value] : audit.ridge_valid_orbits) {
        (void)code;
        if (!first_target) std::cout << ',';
        std::cout << '[';
        const auto& target = value.first;
        for (int face = 0;
             face < static_cast<int>(target.size()); ++face) {
          if (face) std::cout << ',';
          PrintMask(target[face]);
        }
        std::cout << ']';
        first_target = false;
      }
      std::cout << "],\"ridge_valid_target_orbit_literal_multiplicities\":[";
      first_target = true;
      for (const auto& [code, value] : audit.ridge_valid_orbits) {
        (void)code;
        if (!first_target) std::cout << ',';
        std::cout << value.second;
        first_target = false;
      }
      std::cout << "]}";
    }
    std::cout << "],\"pure_118_facet_target_count_across_representatives\":"
              << total_targets
              << ",\"scope\":\"Exact necessary compact incidence census for Record 31 a=2; no finite-volume inference.\"}\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
