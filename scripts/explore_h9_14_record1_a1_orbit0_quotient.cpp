// Exact, interruption-safe incidence quotient for H9/14 Record 1, a=1,
// positive-A4 source orbit 0.
//
// This engine deliberately does not have a default full-search mode.  Use
// --inventory to print the 21 exact seed-row classes, --class I to exhaust
// one class, or --all-classes to request the complete census.  Target and
// checkpoint records are JSON Lines and are flushed immediately.  Only a
// final run_terminal record with whole_census_complete=true is a whole-census
// certificate; a cutoff or an interrupted log is never complete.
//
// The fixed source is
//   graph rows {0,13}, {1,13}, {2,13}; A4 row {0,3,4,5}.
// Its exact ground-set stabilizer is S2({1,2}) x S3({3,4,5}) x
// S7({6,...,12}), of order 60,480.  The 1,278 candidate A5 rows split into
// 21 orbits.  A class search distinguishes one representative seed row and
// quotients partial families only by its stabilizer.  This is essential:
// quotienting a seed-conditioned search by the full source group could prune
// an isomorphic state whose image no longer contains the fixed seed.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include <nauty.h>

namespace {

using Mask = std::uint16_t;
using Words = std::vector<std::uint64_t>;

constexpr int kN = 14;
constexpr Mask kFull = (Mask{1} << kN) - 1;
constexpr int kAdditionCount = 19;
constexpr int kScheduledFacets = 158;
constexpr int kExpectedRidges = 711;
constexpr int kExpectedCandidateCount = 1278;
constexpr int kExpectedSmallCount = 259;
constexpr int kExpectedFiveCoverCount = 634;
constexpr int kExpectedLocalCount = 433;
constexpr int kExpectedClassCount = 21;
constexpr std::uint64_t kSourceGroupOrder = 60480;
constexpr Mask kFour = (Mask{1} << 0) | (Mask{1} << 3) |
                       (Mask{1} << 4) | (Mask{1} << 5);
constexpr std::array<Mask, 3> kEdges = {
    (Mask{1} << 0) | (Mask{1} << 13),
    (Mask{1} << 1) | (Mask{1} << 13),
    (Mask{1} << 2) | (Mask{1} << 13),
};
constexpr const char* kProtocol =
    "h9_14_record1_a1_orbit0_incidence_quotient_jsonl_v1";

struct InlineKey {
  // At full depth there are 3 edge rows + 1 A4 row + 19 A5 rows.
  // The canonical bipartite incidence payload is therefore
  // 1 + ceil(23*14/8) = 42 bytes.  Six inline words are collision-free.
  std::array<std::uint64_t, 6> words{};

  bool operator==(const InlineKey& other) const noexcept {
    return words == other.words;
  }
  bool operator<(const InlineKey& other) const noexcept {
    return words < other.words;
  }
};

static_assert(sizeof(InlineKey) == 48);

struct InlineKeyHash {
  std::size_t operator()(const InlineKey& key) const noexcept {
    std::uint64_t state = 0x243f6a8885a308d3ULL;
    for (std::uint64_t word : key.words) {
      word += 0x9e3779b97f4a7c15ULL;
      word = (word ^ (word >> 30)) * 0xbf58476d1ce4e5b9ULL;
      word = (word ^ (word >> 27)) * 0x94d049bb133111ebULL;
      word ^= word >> 31;
      state ^= word + 0x9e3779b97f4a7c15ULL + (state << 6) + (state >> 2);
    }
    return static_cast<std::size_t>(state);
  }
};

struct LocalConstraint {
  Mask subset = 0;
  int lower = 0;
  int upper = 0;
  int previous = 0;
  int before = 0;
};

struct LinkPrefix {
  std::vector<int> f;
  std::vector<int> h;
  std::vector<int> g;
};

struct Problem {
  std::vector<Mask> shadow;
  std::vector<Mask> candidates;
  std::vector<Mask> small;
  std::vector<Mask> five_covers;
  std::vector<Words> coverage;
  std::vector<std::vector<int>> coverers;
  Words base_five_cover_words;
  std::vector<Words> candidate_five_hit_words;
  std::vector<LocalConstraint> local_constraints;
  std::array<int, 1 << kN> local_index{};
  std::vector<std::array<int, 25>> candidate_local_indices;
};

struct ClassInfo {
  int index = -1;
  InlineKey full_singleton_key{};
  int seed_candidate = -1;
  std::vector<int> members;
};

struct SearchState {
  std::vector<int> chosen;
  Words forbidden;
  Words covered;
  std::vector<std::uint8_t> local_incidence;
  Words surviving_five_covers;
};

struct ClassAudit {
  int class_index = -1;
  int seed_candidate = -1;
  std::uint64_t node_limit = 0;
  std::uint64_t checkpoint_nodes = 1000000;
  bool cutoff = false;
  std::uint64_t nodes = 0;
  std::array<std::uint64_t, kAdditionCount + 1> depth_nodes{};
  std::uint64_t quotient_prunes = 0;
  std::uint64_t dead_facet_count = 0;
  std::uint64_t dead_uncoverable = 0;
  std::uint64_t dead_slot_bound = 0;
  std::uint64_t terminal_covered_states = 0;
  std::uint64_t scheduled_facet_states = 0;
  std::uint64_t pure_seed_stabilizer_targets = 0;
  std::uint64_t ridge_valid_seed_stabilizer_targets = 0;
  std::array<std::unordered_set<InlineKey, InlineKeyHash>,
             kAdditionCount + 1>
      quotient_seen;
  std::array<std::uint64_t, kAdditionCount + 1> final_seen_counts{};
  std::unordered_set<InlineKey, InlineKeyHash> full_target_orbits;
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
  for (Mask edge : kEdges) {
    if ((value & edge) == edge) return false;
  }
  return true;
}

bool HitsAll(Mask value, const std::vector<Mask>& clutter) {
  for (Mask row : clutter) {
    if ((value & row) == 0) return false;
  }
  return true;
}

LinkPrefix PrefixData(Mask linked) {
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
      f[size] += Independent(face) && ((face & kFour) != kFour);
    }
    if (extension == 0) break;
  }
  std::vector<int> h(top + 1);
  for (int index = 0; index <= top; ++index) {
    for (int order = 0; order <= index; ++order) {
      h[index] += (((index - order) & 1) ? -1 : 1) *
                  static_cast<int>(Binomial(dimension - order,
                                            dimension - index)) *
                  f[order];
    }
  }
  std::vector<int> g(top + 1);
  g[0] = h[0];
  for (int index = 1; index <= top; ++index) g[index] = h[index] - h[index - 1];
  return {std::move(f), std::move(h), std::move(g)};
}

int DerivedScheduledRidgeCount() {
  constexpr int dimension = 9;
  std::array<long long, 5> face_prefix{};
  face_prefix[0] = 1;
  for (Mask face = 1; face <= kFull; ++face) {
    const int size = std::popcount(face);
    if (size <= 4 && Independent(face)) ++face_prefix[size];
  }
  --face_prefix[4];  // the fixed A4 row
  std::array<long long, dimension + 1> h{};
  for (int index = 0; index <= 4; ++index) {
    for (int order = 0; order <= index; ++order) {
      h[index] += (((index - order) & 1) ? -1 : 1) *
                  Binomial(dimension - order, dimension - index) *
                  face_prefix[order];
    }
    h[dimension - index] = h[index];
  }
  long long facets = 0;
  for (long long entry : h) facets += entry;
  if (facets != kScheduledFacets) {
    throw std::runtime_error("Record-1 a=1 facet schedule changed");
  }
  long long ridges = 0;
  for (int index = 0; index <= 8; ++index) {
    ridges += Binomial(dimension - index, dimension - 8) * h[index];
  }
  return static_cast<int>(ridges);
}

int WordCount(const Words& words) {
  int answer = 0;
  for (std::uint64_t word : words) answer += std::popcount(word);
  return answer;
}

bool IsSet(const Words& words, int index) {
  return (words[index / 64] >> (index % 64)) & 1u;
}

void Set(Words* words, int index) {
  (*words)[index / 64] |= std::uint64_t{1} << (index % 64);
}

void OrInto(Words* target, const Words& source) {
  for (int index = 0; index < static_cast<int>(target->size()); ++index) {
    (*target)[index] |= source[index];
  }
}

Words Meet(const Words& left, const Words& right) {
  Words answer(left.size());
  for (int index = 0; index < static_cast<int>(left.size()); ++index) {
    answer[index] = left[index] & right[index];
  }
  return answer;
}

int MeetCount(const Words& left, const Words& right) {
  int answer = 0;
  for (int index = 0; index < static_cast<int>(left.size()); ++index) {
    answer += std::popcount(left[index] & right[index]);
  }
  return answer;
}

bool IsFull(const Words& words, int bit_count) {
  const int whole = bit_count / 64;
  for (int index = 0; index < whole; ++index) {
    if (words[index] != ~std::uint64_t{0}) return false;
  }
  const int remainder = bit_count % 64;
  if (remainder == 0) return true;
  const std::uint64_t tail = (std::uint64_t{1} << remainder) - 1;
  return (words[whole] & tail) == tail;
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

unsigned char KeyByte(const InlineKey& key, std::size_t index) {
  return static_cast<unsigned char>((key.words[index / 8] >>
                                     (8 * (index % 8))) &
                                    0xffu);
}

std::size_t KeyByteCount(const InlineKey& key) {
  const std::size_t rows = KeyByte(key, 0);
  const std::size_t bytes = 1 + (rows * kN + 7) / 8;
  if (rows < 5 || rows > 23 || bytes > sizeof(key)) {
    throw std::runtime_error("packed canonical key row count changed");
  }
  return bytes;
}

std::string KeyHex(const InlineKey& key) {
  std::ostringstream output;
  output << std::hex << std::setfill('0');
  for (std::size_t index = 0; index < KeyByteCount(key); ++index) {
    output << std::setw(2) << static_cast<int>(KeyByte(key, index));
  }
  return output.str();
}

InlineKey CanonicalIncidenceKey(const std::vector<Mask>& additions,
                                std::optional<Mask> distinguished_seed) {
  if (additions.empty() || additions.size() > kAdditionCount) {
    throw std::runtime_error("canonical addition-family size changed");
  }
  const int edge_start = kN;
  const int four_start = edge_start + kEdges.size();
  const int addition_start = four_start + 1;
  const int node_count = addition_start + static_cast<int>(additions.size());
  const int words = SETWORDSNEEDED(node_count);
  std::vector<graph> input(static_cast<std::size_t>(words) * node_count);
  std::vector<graph> canonical(static_cast<std::size_t>(words) * node_count);
  EMPTYGRAPH(input.data(), words, node_count);
  for (int row = 0; row < static_cast<int>(kEdges.size()); ++row) {
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (kEdges[row] & (Mask{1} << vertex)) {
        ADDONEEDGE(input.data(), edge_start + row, vertex, words);
      }
    }
  }
  for (int vertex = 0; vertex < kN; ++vertex) {
    if (kFour & (Mask{1} << vertex)) {
      ADDONEEDGE(input.data(), four_start, vertex, words);
    }
  }

  std::vector<Mask> ordered;
  if (distinguished_seed.has_value()) {
    int occurrences = 0;
    for (Mask row : additions) occurrences += row == *distinguished_seed;
    if (occurrences != 1) {
      throw std::runtime_error("distinguished seed is not unique in family");
    }
    ordered.push_back(*distinguished_seed);
    for (Mask row : additions) {
      if (row != *distinguished_seed) ordered.push_back(row);
    }
  } else {
    ordered = additions;
  }
  for (int row = 0; row < static_cast<int>(ordered.size()); ++row) {
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (ordered[row] & (Mask{1} << vertex)) {
        ADDONEEDGE(input.data(), addition_start + row, vertex, words);
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
  add_cell(four_start, addition_start);
  if (distinguished_seed.has_value()) {
    add_cell(addition_start, addition_start + 1);
    add_cell(addition_start + 1, node_count);
  } else {
    add_cell(addition_start, node_count);
  }
  if (position != node_count) {
    throw std::runtime_error("canonical color partition changed");
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

  const int row_count = node_count - kN;
  const int incidence_bits = row_count * kN;
  const std::size_t byte_count = 1 + (incidence_bits + 7) / 8;
  if (byte_count > sizeof(InlineKey)) {
    throw std::runtime_error("inline canonical key capacity exceeded");
  }
  InlineKey code;
  code.words[0] = static_cast<unsigned char>(row_count);
  for (int row = 0; row < row_count; ++row) {
    const set* adjacency = GRAPHROW(canonical.data(), kN + row, words);
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (!ISELEMENT(adjacency, vertex)) continue;
      const int bit = row * kN + vertex;
      const std::size_t byte = 1 + bit / 8;
      const unsigned shift = 8 * (byte % 8) + bit % 8;
      code.words[byte / 8] |= std::uint64_t{1} << shift;
    }
  }

#ifdef RECORD1_VERIFY_PACKED_KEY
  for (int left = 0; left < node_count; ++left) {
    const set* adjacency = GRAPHROW(canonical.data(), left, words);
    for (int right = 0; right < node_count; ++right) {
      const bool dense = ISELEMENT(adjacency, right);
      bool recovered = false;
      int row = -1;
      int vertex = -1;
      if (left < kN && right >= kN) {
        vertex = left;
        row = right - kN;
      } else if (right < kN && left >= kN) {
        vertex = right;
        row = left - kN;
      }
      if (row >= 0) {
        const int bit = row * kN + vertex;
        const std::size_t byte = 1 + bit / 8;
        recovered = (KeyByte(code, byte) >> (bit % 8)) & 1u;
      }
      if (dense != recovered) {
        throw std::runtime_error("packed key failed dense adjacency replay");
      }
    }
  }
#endif
  return code;
}

std::vector<Mask> FamilyMasks(const Problem& problem,
                              const SearchState& state) {
  std::vector<Mask> answer;
  answer.reserve(state.chosen.size());
  for (int candidate : state.chosen) answer.push_back(problem.candidates[candidate]);
  return answer;
}

bool Compatible(const Problem& problem, int candidate,
                const std::vector<std::uint8_t>& incidence) {
  for (int local : problem.candidate_local_indices[candidate]) {
    if (incidence[local] >= problem.local_constraints[local].upper) return false;
  }
  return true;
}

Problem BuildProblem() {
  Problem problem;
  problem.local_index.fill(-1);
  for (int size = 1; size <= 3; ++size) {
    for (Mask subset = 1; subset <= kFull; ++subset) {
      if (std::popcount(subset) != size || !Independent(subset)) continue;
      const LinkPrefix prefix = PrefixData(subset);
      const int top = 5 - size;
      const int previous = prefix.g[top - 1];
      const int before = prefix.g[top];
      LocalConstraint constraint;
      constraint.subset = subset;
      constraint.previous = previous;
      constraint.before = before;
      constraint.lower = previous >= 0
          ? std::max(0, before - MacaulayRaise(previous, top - 1))
          : 0;
      constraint.upper = before - 1;
      if (constraint.lower > constraint.upper) {
        throw std::runtime_error("fixed Record-1 source failed local interval");
      }
      problem.local_index[subset] = problem.local_constraints.size();
      problem.local_constraints.push_back(constraint);
    }
  }
  if (problem.local_constraints.size() != kExpectedLocalCount) {
    throw std::runtime_error("Record-1 local-constraint ground changed");
  }

  std::vector<Mask> independent_fives;
  for (Mask five = 0; five <= kFull; ++five) {
    if (std::popcount(five) != 5 || !Independent(five)) continue;
    independent_fives.push_back(five);
    if ((five & kFour) == kFour) problem.shadow.push_back(five);
  }
  if (independent_fives.size() != 1497 || problem.shadow.size() != 9) {
    throw std::runtime_error("Record-1 independent-five or shadow count changed");
  }

  for (Mask five : independent_fives) {
    if ((five & kFour) == kFour) continue;
    std::array<int, 25> local_indices{};
    int position = 0;
    bool allowed = true;
    for (Mask subset = five; subset; subset = (subset - 1) & five) {
      const int size = std::popcount(subset);
      if (size < 1 || size > 3) continue;
      const int local = problem.local_index[subset];
      if (local < 0) throw std::runtime_error("candidate local support missing");
      local_indices[position++] = local;
      if (problem.local_constraints[local].upper < 1) allowed = false;
    }
    if (position != 25) throw std::runtime_error("A5 local support changed");
    if (allowed) {
      problem.candidates.push_back(five);
      problem.candidate_local_indices.push_back(local_indices);
    }
  }
  if (problem.candidates.size() != kExpectedCandidateCount) {
    throw std::runtime_error("Record-1 candidate count changed");
  }

  std::vector<Mask> base(kEdges.begin(), kEdges.end());
  base.push_back(kFour);
  base.insert(base.end(), problem.shadow.begin(), problem.shadow.end());
  for (Mask value = 1; value <= kFull; ++value) {
    if (std::popcount(value) <= 4 && HitsAll(value, base)) {
      problem.small.push_back(value);
    }
  }
  if (problem.small.size() != kExpectedSmallCount) {
    throw std::runtime_error("Record-1 small-transversal count changed");
  }
  const int coverage_words = (problem.small.size() + 63) / 64;
  problem.coverers.assign(problem.small.size(), {});
  for (int candidate = 0;
       candidate < static_cast<int>(problem.candidates.size()); ++candidate) {
    Words words(coverage_words);
    for (int requirement = 0;
         requirement < static_cast<int>(problem.small.size()); ++requirement) {
      if ((problem.candidates[candidate] & problem.small[requirement]) == 0) {
        Set(&words, requirement);
        problem.coverers[requirement].push_back(candidate);
      }
    }
    problem.coverage.push_back(std::move(words));
  }

  for (Mask cover = 0; cover <= kFull; ++cover) {
    if (std::popcount(cover) != 5) continue;
    bool hits = (cover & kFour) != 0;
    for (Mask edge : kEdges) hits = hits && ((cover & edge) != 0);
    if (hits) problem.five_covers.push_back(cover);
  }
  if (problem.five_covers.size() != kExpectedFiveCoverCount) {
    throw std::runtime_error("Record-1 five-cover count changed");
  }
  const int five_words = (problem.five_covers.size() + 63) / 64;
  problem.base_five_cover_words.assign(five_words, ~std::uint64_t{0});
  if (problem.five_covers.size() % 64) {
    problem.base_five_cover_words.back() =
        (std::uint64_t{1} << (problem.five_covers.size() % 64)) - 1;
  }
  problem.candidate_five_hit_words.assign(
      problem.candidates.size(), Words(five_words));
  for (int candidate = 0;
       candidate < static_cast<int>(problem.candidates.size()); ++candidate) {
    for (int cover = 0;
         cover < static_cast<int>(problem.five_covers.size()); ++cover) {
      if (problem.candidates[candidate] & problem.five_covers[cover]) {
        Set(&problem.candidate_five_hit_words[candidate], cover);
      }
    }
  }
  if (DerivedScheduledRidgeCount() != kExpectedRidges) {
    throw std::runtime_error("Record-1 derived ridge schedule changed");
  }
  return problem;
}

std::vector<ClassInfo> BuildClasses(const Problem& problem) {
  std::map<InlineKey, std::vector<int>> groups;
  for (int candidate = 0;
       candidate < static_cast<int>(problem.candidates.size()); ++candidate) {
    groups[CanonicalIncidenceKey({problem.candidates[candidate]}, std::nullopt)]
        .push_back(candidate);
  }
  std::vector<ClassInfo> classes;
  for (auto& [key, members] : groups) {
    ClassInfo info;
    info.full_singleton_key = key;
    info.members = std::move(members);
    info.seed_candidate = *std::min_element(
        info.members.begin(), info.members.end(), [&](int left, int right) {
          return problem.candidates[left] < problem.candidates[right];
        });
    classes.push_back(std::move(info));
  }
  std::sort(classes.begin(), classes.end(), [&](const ClassInfo& left,
                                                 const ClassInfo& right) {
    return problem.candidates[left.seed_candidate] <
           problem.candidates[right.seed_candidate];
  });
  if (classes.size() != kExpectedClassCount) {
    throw std::runtime_error("Record-1 candidate orbit count changed");
  }
  std::size_t member_sum = 0;
  for (int index = 0; index < static_cast<int>(classes.size()); ++index) {
    classes[index].index = index;
    member_sum += classes[index].members.size();
    if (kSourceGroupOrder % classes[index].members.size() != 0) {
      throw std::runtime_error("candidate orbit size does not divide source group");
    }
  }
  if (member_sum != problem.candidates.size()) {
    throw std::runtime_error("candidate classes do not partition candidates");
  }
  return classes;
}

void PrintIntegerArray(const std::array<std::uint64_t,
                                       kAdditionCount + 1>& values) {
  std::cout << '[';
  for (int index = 0; index <= kAdditionCount; ++index) {
    if (index) std::cout << ',';
    std::cout << values[index];
  }
  std::cout << ']';
}

void PrintMaskArray(const std::vector<Mask>& masks) {
  std::cout << '[';
  for (int index = 0; index < static_cast<int>(masks.size()); ++index) {
    if (index) std::cout << ',';
    std::cout << masks[index];
  }
  std::cout << ']';
}

void EmitHeader(const Problem& problem, const std::vector<ClassInfo>& classes,
                const std::string& mode) {
  std::cout << "{\"event\":\"run_header\",\"protocol\":\"" << kProtocol
            << "\",\"mode\":\"" << mode
            << "\",\"record_index\":1,\"a\":1,"
            << "\"source_orbit_index\":0,\"source_group\":\"S2xS3xS7\","
            << "\"source_group_order\":" << kSourceGroupOrder
            << ",\"candidate_count\":" << problem.candidates.size()
            << ",\"seed_class_count\":" << classes.size()
            << ",\"additional_A5_count\":" << kAdditionCount
            << ",\"scheduled_facet_count\":" << kScheduledFacets
            << ",\"small_transversal_count\":" << problem.small.size()
            << ",\"five_cover_count\":" << problem.five_covers.size()
            << ",\"local_constraint_count\":"
            << problem.local_constraints.size()
            << ",\"proof_policy\":\"only run_terminal with "
               "whole_census_complete=true certifies the whole source "
               "quotient; cutoff, checkpoint, target, and selected-class "
               "terminal events do not\"}\n"
            << std::flush;
}

void EmitInventory(const Problem& problem,
                   const std::vector<ClassInfo>& classes) {
  std::cout << "{\"event\":\"class_inventory\",\"class_count\":"
            << classes.size() << ",\"classes\":[";
  for (int index = 0; index < static_cast<int>(classes.size()); ++index) {
    if (index) std::cout << ',';
    const ClassInfo& info = classes[index];
    const Mask seed = problem.candidates[info.seed_candidate];
    std::cout << "{\"class_index\":" << info.index
              << ",\"seed_mask\":" << seed
              << ",\"candidate_orbit_size\":" << info.members.size()
              << ",\"seed_stabilizer_order\":"
              << kSourceGroupOrder / info.members.size()
              << ",\"full_singleton_key\":\""
              << KeyHex(info.full_singleton_key) << "\"}";
  }
  std::cout << "]}\n" << std::flush;
}

void EmitCheckpoint(const ClassAudit& audit) {
  std::cout << "{\"event\":\"checkpoint\",\"class_index\":"
            << audit.class_index << ",\"nodes\":" << audit.nodes
            << ",\"depth_nodes\":";
  PrintIntegerArray(audit.depth_nodes);
  std::cout << ",\"seen_counts\":[";
  for (int depth = 0; depth <= kAdditionCount; ++depth) {
    if (depth) std::cout << ',';
    std::cout << audit.quotient_seen[depth].size();
  }
  std::cout << "],\"target_orbits_streamed\":"
            << audit.full_target_orbits.size()
            << ",\"proof_complete\":false}\n" << std::flush;
}

void EmitTarget(const Problem& problem, const SearchState& state,
                const ClassAudit& audit, const InlineKey& full_key) {
  std::vector<Mask> additions = FamilyMasks(problem, state);
  std::sort(additions.begin(), additions.end());
  std::cout << "{\"event\":\"target\",\"class_index\":"
            << audit.class_index << ",\"canonical_full_key\":\""
            << KeyHex(full_key) << "\",\"source_A4_mask\":" << kFour
            << ",\"additional_A5_masks\":";
  PrintMaskArray(additions);
  std::cout << ",\"scheduled_facets\":" << kScheduledFacets
            << ",\"active_ridges\":" << kExpectedRidges
            << ",\"proof_complete\":false}\n" << std::flush;
}

SearchState AddCandidate(const Problem& problem, const SearchState& state,
                         int candidate) {
  SearchState next = state;
  next.chosen.push_back(candidate);
  Set(&next.forbidden, candidate);
  OrInto(&next.covered, problem.coverage[candidate]);
  for (int local : problem.candidate_local_indices[candidate]) {
    ++next.local_incidence[local];
  }
  next.surviving_five_covers =
      Meet(state.surviving_five_covers,
           problem.candidate_five_hit_words[candidate]);
  return next;
}

void Finish(const Problem& problem, const SearchState& state,
            ClassAudit* audit) {
  ++audit->terminal_covered_states;
  if (WordCount(state.surviving_five_covers) != kScheduledFacets) return;
  ++audit->scheduled_facet_states;
  for (int local = 0;
       local < static_cast<int>(problem.local_constraints.size()); ++local) {
    const int count = state.local_incidence[local];
    if (count < problem.local_constraints[local].lower ||
        count > problem.local_constraints[local].upper) {
      return;
    }
  }
  std::vector<Mask> clutter(kEdges.begin(), kEdges.end());
  clutter.push_back(kFour);
  clutter.insert(clutter.end(), problem.shadow.begin(), problem.shadow.end());
  const std::vector<Mask> additions = FamilyMasks(problem, state);
  clutter.insert(clutter.end(), additions.begin(), additions.end());
  const std::vector<Mask> blocker = FullMinimalTransversals(clutter);
  bool pure = blocker.size() == kScheduledFacets;
  for (Mask face : blocker) pure = pure && std::popcount(face) == 5;
  if (!pure) return;
  ++audit->pure_seed_stabilizer_targets;

  std::map<Mask, int> ridge_degrees;
  for (Mask facet : blocker) {
    Mask outside = kFull ^ facet;
    while (outside) {
      const Mask bit = outside & static_cast<Mask>(-outside);
      outside ^= bit;
      ++ridge_degrees[facet | bit];
    }
  }
  if (ridge_degrees.size() != kExpectedRidges) return;
  for (const auto& [ridge, degree] : ridge_degrees) {
    (void)ridge;
    if (degree != 2) return;
  }
  ++audit->ridge_valid_seed_stabilizer_targets;
  const InlineKey full_key =
      CanonicalIncidenceKey(additions, std::nullopt);
  if (audit->full_target_orbits.insert(full_key).second) {
    EmitTarget(problem, state, *audit, full_key);
  }
}

void Search(const Problem& problem, const ClassInfo& info,
            const SearchState& state, ClassAudit* audit) {
  if (audit->cutoff) return;
  if (audit->node_limit != 0 && audit->nodes >= audit->node_limit) {
    audit->cutoff = true;
    return;
  }
  ++audit->nodes;
  ++audit->depth_nodes[state.chosen.size()];
  if (audit->checkpoint_nodes != 0 &&
      audit->nodes % audit->checkpoint_nodes == 0) {
    EmitCheckpoint(*audit);
  }

  const int depth = state.chosen.size();
  const Mask seed = problem.candidates[info.seed_candidate];
  const InlineKey partial =
      CanonicalIncidenceKey(FamilyMasks(problem, state), seed);
  if (!audit->quotient_seen[depth].insert(partial).second) {
    ++audit->quotient_prunes;
    return;
  }

  const int current_facets = WordCount(state.surviving_five_covers);
  if (current_facets < kScheduledFacets) {
    ++audit->dead_facet_count;
    return;
  }
  if (depth == kAdditionCount) {
    if (IsFull(state.covered, problem.small.size())) Finish(problem, state, audit);
    return;
  }

  const int remaining = kAdditionCount - depth;
  std::vector<std::uint8_t> available(problem.candidates.size());
  std::vector<int> losses;
  std::vector<int> gains;
  int uncovered_count = 0;
  for (int requirement = 0;
       requirement < static_cast<int>(problem.small.size()); ++requirement) {
    uncovered_count += !IsSet(state.covered, requirement);
  }
  for (int candidate = 0;
       candidate < static_cast<int>(problem.candidates.size()); ++candidate) {
    if (IsSet(state.forbidden, candidate) ||
        !Compatible(problem, candidate, state.local_incidence)) {
      continue;
    }
    const int next_facets = MeetCount(
        state.surviving_five_covers,
        problem.candidate_five_hit_words[candidate]);
    if (next_facets < kScheduledFacets) continue;
    available[candidate] = 1;
    losses.push_back(current_facets - next_facets);
    int gain = 0;
    for (int word = 0; word < static_cast<int>(state.covered.size()); ++word) {
      gain += std::popcount(problem.coverage[candidate][word] &
                            ~state.covered[word]);
    }
    gains.push_back(gain);
  }
  if (static_cast<int>(losses.size()) < remaining) {
    ++audit->dead_uncoverable;
    return;
  }
  std::ranges::sort(losses, std::greater<>());
  int maximum_loss = 0;
  for (int index = 0; index < remaining; ++index) maximum_loss += losses[index];
  if (current_facets - maximum_loss > kScheduledFacets) {
    ++audit->dead_slot_bound;
    return;
  }
  std::ranges::sort(gains, std::greater<>());
  int maximum_gain = 0;
  for (int index = 0; index < remaining; ++index) maximum_gain += gains[index];
  if (maximum_gain < uncovered_count) {
    ++audit->dead_slot_bound;
    return;
  }

  std::vector<int> branches;
  if (uncovered_count == 0) {
    for (int candidate = 0;
         candidate < static_cast<int>(problem.candidates.size()); ++candidate) {
      if (available[candidate]) branches.push_back(candidate);
    }
  } else {
    for (int requirement = 0;
         requirement < static_cast<int>(problem.small.size()); ++requirement) {
      if (IsSet(state.covered, requirement)) continue;
      std::vector<int> feasible;
      for (int candidate : problem.coverers[requirement]) {
        if (available[candidate]) feasible.push_back(candidate);
      }
      if (feasible.empty()) {
        ++audit->dead_uncoverable;
        return;
      }
      if (branches.empty() || feasible.size() < branches.size()) {
        branches = std::move(feasible);
      }
    }
  }
  for (int candidate : branches) {
    Search(problem, info, AddCandidate(problem, state, candidate), audit);
    if (audit->cutoff) return;
  }
}

ClassAudit RunClass(const Problem& problem, const ClassInfo& info,
                    std::uint64_t node_limit,
                    std::uint64_t checkpoint_nodes) {
  ClassAudit audit;
  audit.class_index = info.index;
  audit.seed_candidate = info.seed_candidate;
  audit.node_limit = node_limit;
  audit.checkpoint_nodes = checkpoint_nodes;
  const Mask seed = problem.candidates[info.seed_candidate];
  std::cout << "{\"event\":\"class_started\",\"class_index\":"
            << info.index << ",\"seed_mask\":" << seed
            << ",\"candidate_orbit_size\":" << info.members.size()
            << ",\"seed_stabilizer_order\":"
            << kSourceGroupOrder / info.members.size()
            << ",\"proof_complete\":false}\n" << std::flush;

  SearchState initial;
  initial.forbidden.assign((problem.candidates.size() + 63) / 64, 0);
  initial.covered.assign((problem.small.size() + 63) / 64, 0);
  initial.local_incidence.assign(problem.local_constraints.size(), 0);
  initial.surviving_five_covers = problem.base_five_cover_words;
  initial = AddCandidate(problem, initial, info.seed_candidate);
  Search(problem, info, initial, &audit);

  for (int depth = 0; depth <= kAdditionCount; ++depth) {
    audit.final_seen_counts[depth] = audit.quotient_seen[depth].size();
  }
  std::cout << "{\"event\":\"class_terminal\",\"class_index\":"
            << audit.class_index << ",\"terminal_status\":\""
            << (audit.cutoff ? "cutoff_max_nodes" : "complete_exhaustive")
            << "\",\"class_proof_complete\":"
            << (audit.cutoff ? "false" : "true")
            << ",\"nodes\":" << audit.nodes << ",\"depth_nodes\":";
  PrintIntegerArray(audit.depth_nodes);
  std::cout << ",\"quotient_seen_counts\":";
  PrintIntegerArray(audit.final_seen_counts);
  std::cout << ",\"quotient_prunes\":" << audit.quotient_prunes
            << ",\"dead_facet_count\":" << audit.dead_facet_count
            << ",\"dead_uncoverable_count\":" << audit.dead_uncoverable
            << ",\"dead_slot_bound_count\":" << audit.dead_slot_bound
            << ",\"terminal_covered_states\":"
            << audit.terminal_covered_states
            << ",\"scheduled_facet_states\":"
            << audit.scheduled_facet_states
            << ",\"pure_seed_stabilizer_target_count\":"
            << audit.pure_seed_stabilizer_targets
            << ",\"ridge_valid_seed_stabilizer_target_count\":"
            << audit.ridge_valid_seed_stabilizer_targets
            << ",\"full_target_orbit_count_within_class\":"
            << audit.full_target_orbits.size() << "}\n" << std::flush;
  return audit;
}

struct Options {
  enum class Mode { kInventory, kClass, kAll, kCanonicalFull, kCanonicalSeed };
  Mode mode = Mode::kInventory;
  int class_index = -1;
  std::uint64_t max_nodes = 0;
  std::uint64_t checkpoint_nodes = 1000000;
};

std::uint64_t ParseUnsigned(const std::string& value, const char* name) {
  std::size_t used = 0;
  const unsigned long long parsed = std::stoull(value, &used);
  if (used != value.size()) throw std::runtime_error(std::string("bad ") + name);
  return parsed;
}

Options ParseOptions(int argc, char** argv) {
  Options options;
  bool mode_seen = false;
  for (int index = 1; index < argc; ++index) {
    const std::string argument = argv[index];
    if (argument == "--inventory") {
      if (mode_seen) throw std::runtime_error("choose exactly one mode");
      options.mode = Options::Mode::kInventory;
      mode_seen = true;
    } else if (argument == "--class") {
      if (mode_seen || ++index >= argc) throw std::runtime_error("bad --class");
      options.mode = Options::Mode::kClass;
      options.class_index = static_cast<int>(ParseUnsigned(argv[index], "class"));
      mode_seen = true;
    } else if (argument == "--all-classes") {
      if (mode_seen) throw std::runtime_error("choose exactly one mode");
      options.mode = Options::Mode::kAll;
      mode_seen = true;
    } else if (argument == "--canonicalize-full-stdin") {
      if (mode_seen) throw std::runtime_error("choose exactly one mode");
      options.mode = Options::Mode::kCanonicalFull;
      mode_seen = true;
    } else if (argument == "--canonicalize-seed-stdin") {
      if (mode_seen || ++index >= argc) {
        throw std::runtime_error("bad --canonicalize-seed-stdin");
      }
      options.mode = Options::Mode::kCanonicalSeed;
      options.class_index = static_cast<int>(ParseUnsigned(argv[index], "class"));
      mode_seen = true;
    } else if (argument == "--max-nodes") {
      if (++index >= argc) throw std::runtime_error("bad --max-nodes");
      options.max_nodes = ParseUnsigned(argv[index], "max-nodes");
      if (options.max_nodes == 0) throw std::runtime_error("max-nodes must be positive");
    } else if (argument == "--checkpoint-nodes") {
      if (++index >= argc) throw std::runtime_error("bad --checkpoint-nodes");
      options.checkpoint_nodes = ParseUnsigned(argv[index], "checkpoint-nodes");
    } else {
      throw std::runtime_error("unknown argument: " + argument);
    }
  }
  return options;
}

void CanonicalizeStdin(const Problem& problem,
                       const std::vector<ClassInfo>& classes,
                       std::optional<int> seed_class) {
  std::optional<Mask> seed;
  if (seed_class.has_value()) {
    if (*seed_class < 0 || *seed_class >= static_cast<int>(classes.size())) {
      throw std::runtime_error("canonical seed class is out of range");
    }
    seed = problem.candidates[classes[*seed_class].seed_candidate];
  }
  int row_index = 0;
  int count = 0;
  while (std::cin >> count) {
    if (count <= 0 || count > kAdditionCount) {
      throw std::runtime_error("canonical family count out of range");
    }
    std::vector<Mask> family;
    for (int position = 0; position < count; ++position) {
      unsigned raw = 0;
      if (!(std::cin >> raw) || raw > kFull) {
        throw std::runtime_error("malformed canonical family");
      }
      const Mask value = static_cast<Mask>(raw);
      if (!std::binary_search(problem.candidates.begin(),
                              problem.candidates.end(), value)) {
        throw std::runtime_error("canonical family contains a noncandidate");
      }
      family.push_back(value);
    }
    std::sort(family.begin(), family.end());
    if (std::adjacent_find(family.begin(), family.end()) != family.end()) {
      throw std::runtime_error("canonical family contains a duplicate");
    }
    const InlineKey key = CanonicalIncidenceKey(family, seed);
    std::cout << "{\"event\":\"canonical_family\",\"row_index\":"
              << row_index++ << ",\"family_size\":" << family.size()
              << ",\"canonical_key\":\"" << KeyHex(key) << "\"}\n";
  }
  std::cout << "{\"event\":\"canonical_terminal\",\"family_count\":"
            << row_index << ",\"terminal_status\":\"complete_eof\","
            << "\"proof_complete\":true}\n" << std::flush;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Options options = ParseOptions(argc, argv);
    const Problem problem = BuildProblem();
    const std::vector<ClassInfo> classes = BuildClasses(problem);

    if (options.mode == Options::Mode::kCanonicalFull ||
        options.mode == Options::Mode::kCanonicalSeed) {
      CanonicalizeStdin(problem, classes,
                        options.mode == Options::Mode::kCanonicalSeed
                            ? std::optional<int>(options.class_index)
                            : std::nullopt);
      return 0;
    }

    const std::string mode =
        options.mode == Options::Mode::kInventory ? "inventory" :
        options.mode == Options::Mode::kClass ? "selected_class" :
                                                "all_classes";
    EmitHeader(problem, classes, mode);
    EmitInventory(problem, classes);
    if (options.mode == Options::Mode::kInventory) {
      std::cout << "{\"event\":\"run_terminal\","
                << "\"terminal_status\":\"inventory_complete\","
                << "\"enumeration_complete\":false,"
                << "\"whole_census_complete\":false}\n" << std::flush;
      return 0;
    }

    std::vector<int> selected;
    if (options.mode == Options::Mode::kClass) {
      if (options.class_index < 0 ||
          options.class_index >= static_cast<int>(classes.size())) {
        throw std::runtime_error("class index is out of range");
      }
      selected.push_back(options.class_index);
    } else {
      for (int index = 0; index < static_cast<int>(classes.size()); ++index) {
        selected.push_back(index);
      }
    }

    bool all_selected_complete = true;
    std::set<InlineKey> all_target_orbits;
    for (int class_index : selected) {
      ClassAudit audit = RunClass(problem, classes[class_index],
                                  options.max_nodes,
                                  options.checkpoint_nodes);
      all_selected_complete = all_selected_complete && !audit.cutoff;
      all_target_orbits.insert(audit.full_target_orbits.begin(),
                               audit.full_target_orbits.end());
      if (audit.cutoff) break;
    }
    const bool whole = options.mode == Options::Mode::kAll &&
                       all_selected_complete &&
                       selected.size() == classes.size();
    std::cout << "{\"event\":\"run_terminal\",\"terminal_status\":\""
              << (!all_selected_complete ? "cutoff_max_nodes" :
                  whole ? "complete_exhaustive" :
                          "complete_selected_classes_only")
              << "\",\"enumeration_complete\":"
              << (whole ? "true" : "false")
              << ",\"whole_census_complete\":"
              << (whole ? "true" : "false")
              << ",\"selected_class_count\":" << selected.size()
              << ",\"global_target_orbit_count_in_this_run\":"
              << all_target_orbits.size() << "}\n" << std::flush;
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
