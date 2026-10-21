// Independent corrected a=0 census for H9/14 record 18.
//
// The proof input is the balanced-cut lemma for the two one-centre blocker
// sections.  Unlike the withdrawn shortcut, this program never requires the
// whole blocker B(H) to be 3-uniform.  It tracks only the visible three-
// transversals (those not containing all of L or R), applies the inhomogeneous
// six-set ridge equations, and finally constructs the literal blocker of E+H.
//
// Search is an earliest-coverer recursion with coverage-descending candidate
// order and no selected-set memo.  This is structurally separate from the
// production probe.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {

using Mask = std::uint16_t;
using Word = std::vector<std::uint64_t>;
constexpr int kN = 14;
constexpr int kWSize = 12;
constexpr Mask kW = (Mask{1} << kWSize) - 1;
constexpr Mask kL = (Mask{1} << 0) | (Mask{1} << 1);
constexpr Mask kR = (Mask{1} << 2) | (Mask{1} << 3);
constexpr Mask kI = kW ^ kL ^ kR;
constexpr std::array<Mask, 5> kEdges = {
    (Mask{1} << 0) | (Mask{1} << 12),
    (Mask{1} << 1) | (Mask{1} << 12),
    (Mask{1} << 2) | (Mask{1} << 13),
    (Mask{1} << 3) | (Mask{1} << 13),
    (Mask{1} << 12) | (Mask{1} << 13),
};

struct CutPair {
  Mask left_a;
  Mask left_b;
  Mask right_a;
  Mask right_b;
  std::array<int, 6> key;
};

struct Candidate {
  Mask row;
  std::uint64_t small_coverage;
  Word triple_keep;
  int small_count;
  int triple_reduction_count;
};

struct Survivor {
  std::array<int, 6> cut_key;
  std::vector<Mask> rows;
};

struct OrbitTranscript {
  std::array<int, 6> key{};
  int base_size = 0;
  int additions = 0;
  int candidate_count = 0;
  int small_count = 0;
  std::uint64_t nodes = 0;
  std::uint64_t no_small = 0;
  std::uint64_t visible_triple_ridge = 0;
  std::uint64_t literal_pure = 0;
  std::uint64_t ridge = 0;
};

bool HitsAll(Mask value, const std::vector<Mask>& clutter) {
  return std::all_of(clutter.begin(), clutter.end(),
                     [value](Mask row) { return (value & row) != 0; });
}

bool Visible(Mask value) {
  return (value & kL) != kL && (value & kR) != kR;
}

std::vector<Mask> MasksOfSize(int vertex_count, int cardinality) {
  std::vector<Mask> answer;
  for (Mask value = 0; value < (Mask{1} << vertex_count); ++value) {
    if (std::popcount(value) == cardinality) answer.push_back(value);
  }
  return answer;
}

std::vector<Mask> MinimalTransversals(const std::vector<Mask>& clutter) {
  std::vector<Mask> answer;
  for (Mask value = 1; value < (Mask{1} << kN); ++value) {
    if (!HitsAll(value, clutter)) continue;
    bool minimal = true;
    for (int vertex = 0; vertex < kN; ++vertex) {
      if ((value >> vertex) & 1u) {
        minimal = minimal &&
                  !HitsAll(value ^ (Mask{1} << vertex), clutter);
      }
    }
    if (minimal) answer.push_back(value);
  }
  return answer;
}

std::array<int, 6> CanonicalCutKey(Mask left, Mask right) {
  const int r_count = std::popcount(static_cast<unsigned>(left & kR));
  const int l_count = std::popcount(static_cast<unsigned>(right & kL));
  std::array<int, 4> cells{};
  for (int isolate = 4; isolate < 12; ++isolate) {
    const int code = 2 * ((left >> isolate) & 1u) +
                     ((right >> isolate) & 1u);
    ++cells[code];
  }
  bool first = true;
  std::array<int, 6> best{};
  for (int complement_left = 0; complement_left < 2; ++complement_left) {
    for (int complement_right = 0; complement_right < 2; ++complement_right) {
      std::array<int, 4> moved{};
      for (int left_bit = 0; left_bit < 2; ++left_bit) {
        for (int right_bit = 0; right_bit < 2; ++right_bit) {
          const int source = 2 * (left_bit ^ complement_left) +
                             (right_bit ^ complement_right);
          moved[2 * left_bit + right_bit] = cells[source];
        }
      }
      const int rr = complement_left ? 2 - r_count : r_count;
      const int ll = complement_right ? 2 - l_count : l_count;
      const std::array<int, 6> direct =
          {rr, ll, moved[0], moved[1], moved[2], moved[3]};
      const std::array<int, 6> swapped =
          {ll, rr, moved[0], moved[2], moved[1], moved[3]};
      for (const auto& key : {direct, swapped}) {
        if (first || key < best) {
          best = key;
          first = false;
        }
      }
    }
  }
  return best;
}

std::vector<CutPair> CutRepresentatives() {
  const Mask left_ground = kR | kI;
  const Mask right_ground = kL | kI;
  std::vector<Mask> left_shores;
  std::vector<Mask> right_shores;
  for (const Mask row : MasksOfSize(kWSize, 5)) {
    if ((row & ~left_ground) == 0 && row < (left_ground ^ row)) {
      left_shores.push_back(row);
    }
    if ((row & ~right_ground) == 0 && row < (right_ground ^ row)) {
      right_shores.push_back(row);
    }
  }
  if (left_shores.size() != 126 || right_shores.size() != 126) {
    throw std::runtime_error("balanced cut shore count changed");
  }
  std::map<std::array<int, 6>, CutPair> representatives;
  for (const Mask left : left_shores) {
    for (const Mask right : right_shores) {
      const auto key = CanonicalCutKey(left, right);
      representatives.try_emplace(
          key, CutPair{left, Mask(left_ground ^ left), right,
                       Mask(right_ground ^ right), key});
    }
  }
  std::vector<CutPair> answer;
  for (const auto& [key, pair] : representatives) answer.push_back(pair);
  if (answer.size() != 9) {
    throw std::runtime_error("balanced cut-pair orbit count is not nine");
  }
  return answer;
}

void SetBit(Word* word, int index) {
  (*word)[index / 64] |= std::uint64_t{1} << (index % 64);
}

void ClearBit(Word* word, int index) {
  (*word)[index / 64] &= ~(std::uint64_t{1} << (index % 64));
}

bool TestBit(const Word& word, int index) {
  return (word[index / 64] >> (index % 64)) & 1u;
}

Word Intersect(const Word& left, const Word& right) {
  Word answer(left.size());
  for (int block = 0; block < static_cast<int>(left.size()); ++block) {
    answer[block] = left[block] & right[block];
  }
  return answer;
}

int CountBits(const Word& word) {
  int answer = 0;
  for (const auto block : word) answer += std::popcount(block);
  return answer;
}

std::vector<int> BitIndices(const Word& word) {
  std::vector<int> answer;
  for (int block = 0; block < static_cast<int>(word.size()); ++block) {
    std::uint64_t value = word[block];
    while (value) {
      const int offset = std::countr_zero(value);
      answer.push_back(64 * block + offset);
      value &= value - 1;
    }
  }
  return answer;
}

class Search {
 public:
  Search(CutPair cuts, std::vector<Survivor>* survivors)
      : cuts_(cuts), survivors_(survivors) {
    base_ = {cuts.left_a, cuts.left_b, cuts.right_a, cuts.right_b};
    std::sort(base_.begin(), base_.end());
    base_.erase(std::unique(base_.begin(), base_.end()), base_.end());
    additions_needed_ = 8 - static_cast<int>(base_.size());
    for (const Mask row : MasksOfSize(kWSize, 5)) {
      if ((row & kL) != 0 && (row & kR) != 0) candidate_rows_.push_back(row);
    }
    if (candidate_rows_.size() != 344) {
      throw std::runtime_error("mixed-row candidate count changed");
    }
    for (int vertex = 0; vertex < kWSize; ++vertex) {
      int base_degree = 0;
      for (const Mask row : base_) base_degree += (row >> vertex) & 1u;
      caps_[vertex] = (vertex < 4 ? 5 : 6) - base_degree;
      if (caps_[vertex] < 0) throw std::runtime_error("negative residual cap");
    }

    for (int order = 1; order <= 2; ++order) {
      for (const Mask value : MasksOfSize(kWSize, order)) {
        if (Visible(value) && HitsAll(value, base_)) small_.push_back(value);
      }
    }
    if (small_.size() >= 64) {
      throw std::runtime_error("small-transversal table escaped one word");
    }
    for (const Mask triple : MasksOfSize(kWSize, 3)) {
      if (Visible(triple) && HitsAll(triple, base_)) triples_.push_back(triple);
    }
    const int triple_blocks = (static_cast<int>(triples_.size()) + 63) / 64;
    initial_triples_.assign(triple_blocks, ~std::uint64_t{0});
    if (!triples_.empty() && triples_.size() % 64 != 0) {
      initial_triples_.back() =
          (std::uint64_t{1} << (triples_.size() % 64)) - 1;
    }

    std::vector<Candidate> candidates;
    for (const Mask row : candidate_rows_) {
      std::uint64_t coverage = 0;
      for (int index = 0; index < static_cast<int>(small_.size()); ++index) {
        if ((row & small_[index]) == 0) {
          coverage |= std::uint64_t{1} << index;
        }
      }
      Word keep(triple_blocks);
      int reductions = 0;
      for (int index = 0; index < static_cast<int>(triples_.size()); ++index) {
        if ((row & triples_[index]) != 0) {
          SetBit(&keep, index);
        } else {
          ++reductions;
        }
      }
      candidates.push_back(
          {row, coverage, std::move(keep), std::popcount(coverage), reductions});
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const Candidate& left, const Candidate& right) {
                return std::tuple(left.small_count, left.triple_reduction_count,
                                  left.row) >
                       std::tuple(right.small_count, right.triple_reduction_count,
                                  right.row);
              });
    candidates_ = std::move(candidates);

    const int candidate_blocks =
        (static_cast<int>(candidates_.size()) + 63) / 64;
    small_coverers_.assign(small_.size(), Word(candidate_blocks));
    for (int candidate = 0; candidate < static_cast<int>(candidates_.size());
         ++candidate) {
      for (int constraint = 0; constraint < static_cast<int>(small_.size());
           ++constraint) {
        if ((candidates_[candidate].small_coverage >> constraint) & 1u) {
          SetBit(&small_coverers_[constraint], candidate);
        }
      }
    }
    full_small_ = small_.empty()
                      ? 0
                      : (std::uint64_t{1} << small_.size()) - 1;

    std::vector<Mask> fixed_facets;
    for (int u = 0; u < kWSize; ++u) {
      if (!((cuts_.right_a >> u) & 1u)) continue;
      for (int v = 0; v < kWSize; ++v) {
        if ((cuts_.right_b >> v) & 1u) {
          fixed_facets.push_back(kR | (Mask{1} << 12) |
                                 (Mask{1} << u) | (Mask{1} << v));
        }
      }
    }
    for (int u = 0; u < kWSize; ++u) {
      if (!((cuts_.left_a >> u) & 1u)) continue;
      for (int v = 0; v < kWSize; ++v) {
        if ((cuts_.left_b >> v) & 1u) {
          fixed_facets.push_back(kL | (Mask{1} << 13) |
                                 (Mask{1} << u) | (Mask{1} << v));
        }
      }
    }
    if (fixed_facets.size() != 50) {
      throw std::runtime_error("fixed X/Z facet count is not fifty");
    }
    for (const Mask q : MasksOfSize(kWSize, 4)) {
      const Mask six = q | (Mask{1} << 12) | (Mask{1} << 13);
      int fixed = 0;
      for (const Mask facet : fixed_facets) fixed += (facet & six) == facet;
      if (fixed > 2) cut_impossible_ = true;
      q_constraints_.push_back({q, fixed});
    }
  }

  OrbitTranscript Run() {
    OrbitTranscript transcript;
    transcript.key = cuts_.key;
    transcript.base_size = static_cast<int>(base_.size());
    transcript.additions = additions_needed_;
    transcript.candidate_count = static_cast<int>(candidates_.size());
    transcript.small_count = static_cast<int>(small_.size());
    if (!cut_impossible_) {
      const int blocks = (static_cast<int>(candidates_.size()) + 63) / 64;
      Word allowed(blocks, ~std::uint64_t{0});
      if (candidates_.size() % 64 != 0) {
        allowed.back() =
            (std::uint64_t{1} << (candidates_.size() % 64)) - 1;
      }
      Explore({}, allowed, 0, initial_triples_, {}, &transcript);
    }
    transcript.nodes = nodes_;
    return transcript;
  }

 private:
  bool CapFeasible(const std::array<std::uint8_t, kWSize>& degree,
                   Mask row) const {
    for (int vertex = 0; vertex < kWSize; ++vertex) {
      if (((row >> vertex) & 1u) && degree[vertex] >= caps_[vertex]) {
        return false;
      }
    }
    return true;
  }

  std::array<std::uint8_t, kWSize> Add(
      std::array<std::uint8_t, kWSize> degree, Mask row) const {
    for (int vertex = 0; vertex < kWSize; ++vertex) {
      degree[vertex] += (row >> vertex) & 1u;
    }
    return degree;
  }

  Word Compatible(const Word& source,
                  const std::array<std::uint8_t, kWSize>& degree) const {
    Word answer(source.size());
    for (const int index : BitIndices(source)) {
      if (index < static_cast<int>(candidates_.size()) &&
          CapFeasible(degree, candidates_[index].row)) {
        SetBit(&answer, index);
      }
    }
    return answer;
  }

  int YDegree(Mask q, const Word& current) const {
    int answer = 0;
    for (int index = 0; index < static_cast<int>(triples_.size()); ++index) {
      if (TestBit(current, index) &&
          (triples_[index] & q) == triples_[index]) {
        ++answer;
      }
    }
    return answer;
  }

  bool RidgeEquations(const Word& current) const {
    for (const auto& [q, fixed] : q_constraints_) {
      const int degree = fixed + YDegree(q, current);
      if (degree != 0 && degree != 2) return false;
    }
    return true;
  }

  Word TripleReducers(const Word& allowed, const Word& current,
                      const std::array<std::uint8_t, kWSize>& degree,
                      bool reducers) const {
    Word answer(allowed.size());
    for (const int index : BitIndices(allowed)) {
      if (index >= static_cast<int>(candidates_.size()) ||
          !CapFeasible(degree, candidates_[index].row)) {
        continue;
      }
      bool reduces = false;
      for (int block = 0; block < static_cast<int>(current.size()); ++block) {
        reduces = reduces ||
                  (current[block] & ~candidates_[index].triple_keep[block]);
      }
      if (reduces == reducers) SetBit(&answer, index);
    }
    return answer;
  }

  void Audit(const std::vector<int>& selected, OrbitTranscript* transcript) {
    std::vector<Mask> h = base_;
    for (const int index : selected) h.push_back(candidates_[index].row);
    std::sort(h.begin(), h.end());
    if (h.size() != 8 || std::adjacent_find(h.begin(), h.end()) != h.end()) {
      throw std::runtime_error("final H is not an eight-row clutter");
    }
    std::vector<Mask> clutter(kEdges.begin(), kEdges.end());
    clutter.insert(clutter.end(), h.begin(), h.end());
    const auto blocker = MinimalTransversals(clutter);
    std::map<int, int> sizes;
    for (const Mask value : blocker) ++sizes[std::popcount(value)];
    if (sizes != std::map<int, int>{{5, 110}}) return;
    ++transcript->literal_pure;
    std::map<Mask, int> ridges;
    for (const Mask facet : blocker) {
      for (int vertex = 0; vertex < kN; ++vertex) {
        if (!((facet >> vertex) & 1u)) {
          ++ridges[facet | (Mask{1} << vertex)];
        }
      }
    }
    if (!std::all_of(ridges.begin(), ridges.end(),
                     [](const auto& item) { return item.second == 2; })) {
      return;
    }
    ++transcript->ridge;
    survivors_->push_back({cuts_.key, std::move(h)});
  }

  void Explore(std::vector<int> selected, const Word& allowed,
               std::uint64_t small_covered, const Word& triples,
               const std::array<std::uint8_t, kWSize>& degree,
               OrbitTranscript* transcript) {
    ++nodes_;
    const int remaining = additions_needed_ - static_cast<int>(selected.size());
    if (remaining == 0) {
      if (small_covered != full_small_) return;
      ++transcript->no_small;
      if (CountBits(triples) != 60 || !RidgeEquations(triples)) return;
      ++transcript->visible_triple_ridge;
      Audit(selected, transcript);
      return;
    }
    if (CountBits(allowed) < remaining || CountBits(triples) < 60) return;
    if (CountBits(triples) == 60 && !RidgeEquations(triples)) return;

    Word choices(allowed.size());
    int best_count = static_cast<int>(candidates_.size()) + 1;
    for (int constraint = 0; constraint < static_cast<int>(small_.size());
         ++constraint) {
      if ((small_covered >> constraint) & 1u) continue;
      const Word possible =
          Compatible(Intersect(allowed, small_coverers_[constraint]), degree);
      const int count = CountBits(possible);
      if (count == 0) return;
      if (count < best_count) {
        best_count = count;
        choices = possible;
      }
    }
    if (CountBits(triples) > 60) {
      const Word reducers = TripleReducers(allowed, triples, degree, true);
      const int count = CountBits(reducers);
      if (count == 0) return;
      if (count < best_count) {
        best_count = count;
        choices = reducers;
      }
    } else {
      const Word preservers = TripleReducers(allowed, triples, degree, false);
      const int count = CountBits(preservers);
      if (count == 0) return;
      if (count < best_count) {
        best_count = count;
        choices = preservers;
      }
    }
    if (best_count == static_cast<int>(candidates_.size()) + 1) {
      choices = Compatible(allowed, degree);
    }

    Word forbidden_earlier(allowed.size());
    for (const int index : BitIndices(choices)) {
      Word next_allowed = allowed;
      for (int block = 0; block < static_cast<int>(allowed.size()); ++block) {
        next_allowed[block] &= ~forbidden_earlier[block];
      }
      ClearBit(&next_allowed, index);
      Word next_triples = triples;
      for (int block = 0; block < static_cast<int>(triples.size()); ++block) {
        next_triples[block] &= candidates_[index].triple_keep[block];
      }
      selected.push_back(index);
      Explore(selected, next_allowed,
              small_covered | candidates_[index].small_coverage,
              next_triples, Add(degree, candidates_[index].row), transcript);
      selected.pop_back();
      SetBit(&forbidden_earlier, index);
    }
  }

  CutPair cuts_;
  std::vector<Survivor>* survivors_;
  std::vector<Mask> base_;
  std::vector<Mask> candidate_rows_;
  std::vector<Candidate> candidates_;
  std::array<int, kWSize> caps_{};
  std::vector<Mask> small_;
  std::vector<Mask> triples_;
  Word initial_triples_;
  std::vector<Word> small_coverers_;
  std::vector<std::pair<Mask, int>> q_constraints_;
  std::uint64_t full_small_ = 0;
  int additions_needed_ = 0;
  bool cut_impossible_ = false;
  std::uint64_t nodes_ = 0;
};

void PrintIntArray(const std::array<int, 6>& values) {
  std::cout << '[';
  for (int index = 0; index < 6; ++index) {
    if (index) std::cout << ',';
    std::cout << values[index];
  }
  std::cout << ']';
}

void PrintMaskArray(const std::vector<Mask>& values) {
  std::cout << '[';
  for (int index = 0; index < static_cast<int>(values.size()); ++index) {
    if (index) std::cout << ',';
    std::cout << values[index];
  }
  std::cout << ']';
}

std::vector<std::array<int, 4>> ActiveLeafAutomorphisms() {
  std::set<std::array<int, 4>> found;
  for (int swap_left = 0; swap_left < 2; ++swap_left) {
    for (int swap_right = 0; swap_right < 2; ++swap_right) {
      for (int exchange = 0; exchange < 2; ++exchange) {
        std::array<int, 4> move{};
        for (int source = 0; source < 4; ++source) {
          const int side = source / 2;
          int offset = source % 2;
          if ((side == 0 && swap_left) || (side == 1 && swap_right)) {
            offset ^= 1;
          }
          const int target_side = side ^ exchange;
          move[source] = 2 * target_side + offset;
        }
        found.insert(move);
      }
    }
  }
  if (found.size() != 8) {
    throw std::runtime_error("active leaf automorphism order changed");
  }
  return {found.begin(), found.end()};
}

// Canonicalizes the 8-by-12 incidence matrix without graph-isomorphism
// software.  Rows may be permuted, the eight isolate columns may be permuted,
// and the four active leaf columns move under the order-eight double-star
// automorphism group.  Sorting isolate column bit-patterns quotients S_8
// exactly, so no isolate permutation loop is needed.
std::vector<int> CanonicalFamilyKey(const std::vector<Mask>& rows) {
  if (rows.size() != 8 || std::set<Mask>(rows.begin(), rows.end()).size() != 8) {
    throw std::runtime_error("canonicalizer expects eight distinct rows");
  }
  const auto actions = ActiveLeafAutomorphisms();
  std::array<int, 8> order{0, 1, 2, 3, 4, 5, 6, 7};
  std::vector<int> best;
  bool first = true;
  do {
    std::array<int, 8> isolate_columns{};
    for (int isolate = 4; isolate < 12; ++isolate) {
      int column = 0;
      for (int position = 0; position < 8; ++position) {
        column |= int((rows[order[position]] >> isolate) & 1u) << position;
      }
      isolate_columns[isolate - 4] = column;
    }
    std::sort(isolate_columns.begin(), isolate_columns.end());
    for (const auto& move : actions) {
      std::vector<int> key;
      key.reserve(16);
      for (int position = 0; position < 8; ++position) {
        int active = 0;
        for (int source = 0; source < 4; ++source) {
          if ((rows[order[position]] >> source) & 1u) {
            active |= 1 << move[source];
          }
        }
        key.push_back(active);
      }
      key.insert(key.end(), isolate_columns.begin(), isolate_columns.end());
      if (first || key < best) {
        best = std::move(key);
        first = false;
      }
    }
  } while (std::next_permutation(order.begin(), order.end()));
  return best;
}

void PrintIntVector(const std::vector<int>& values) {
  std::cout << '[';
  for (int index = 0; index < static_cast<int>(values.size()); ++index) {
    if (index) std::cout << ',';
    std::cout << values[index];
  }
  std::cout << ']';
}

std::vector<Mask> ParseRows(const std::string& text) {
  std::vector<Mask> rows;
  std::stringstream stream(text);
  std::string token;
  while (std::getline(stream, token, ',')) {
    if (token.empty()) continue;
    const long value = std::stol(token);
    if (value <= 0 || value > std::numeric_limits<Mask>::max()) {
      throw std::runtime_error("row mask outside canonicalizer range");
    }
    rows.push_back(static_cast<Mask>(value));
  }
  return rows;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc > 1) {
      if (std::string(argv[1]) != "--canonicalize") {
        throw std::runtime_error("unknown command-line mode");
      }
      std::cout << '[';
      for (int index = 2; index < argc; ++index) {
        if (index != 2) std::cout << ',';
        PrintIntVector(CanonicalFamilyKey(ParseRows(argv[index])));
      }
      std::cout << "]\n";
      return 0;
    }
    const auto cuts = CutRepresentatives();
    std::vector<Survivor> survivors;
    std::vector<OrbitTranscript> transcripts;
    for (const CutPair& pair : cuts) {
      Search search(pair, &survivors);
      transcripts.push_back(search.Run());
    }
    std::uint64_t total_nodes = 0;
    std::uint64_t total_no_small = 0;
    std::uint64_t total_triple = 0;
    std::uint64_t total_pure = 0;
    std::uint64_t total_ridge = 0;
    for (const auto& row : transcripts) {
      total_nodes += row.nodes;
      total_no_small += row.no_small;
      total_triple += row.visible_triple_ridge;
      total_pure += row.literal_pure;
      total_ridge += row.ridge;
    }
    if (survivors.size() != total_ridge) {
      throw std::runtime_error("survivor archive count does not match ridge count");
    }
    std::map<std::vector<int>, int> canonical_orbits;
    std::vector<std::vector<int>> canonical_keys;
    canonical_keys.reserve(survivors.size());
    for (const auto& survivor : survivors) {
      auto key = CanonicalFamilyKey(survivor.rows);
      ++canonical_orbits[key];
      canonical_keys.push_back(std::move(key));
    }

    std::cout << "{\n"
              << "  \"kind\":\"independent_h9_14_record18_a0_balanced_cut_audit\",\n"
              << "  \"cut_pair_orbit_count\":" << cuts.size()
              << ",\n  \"orbit_transcripts\":[";
    for (int index = 0; index < static_cast<int>(transcripts.size()); ++index) {
      if (index) std::cout << ',';
      const auto& row = transcripts[index];
      std::cout << "{\"key\":";
      PrintIntArray(row.key);
      std::cout << ",\"base_size\":" << row.base_size
                << ",\"additions\":" << row.additions
                << ",\"candidate_count\":" << row.candidate_count
                << ",\"small_count\":" << row.small_count
                << ",\"nodes\":" << row.nodes
                << ",\"no_small\":" << row.no_small
                << ",\"visible_triple_ridge\":"
                << row.visible_triple_ridge
                << ",\"literal_pure_5_110\":" << row.literal_pure
                << ",\"ridge_valid\":" << row.ridge << '}';
    }
    std::cout << "],\n  \"totals\":{\"nodes\":" << total_nodes
              << ",\"no_small\":" << total_no_small
              << ",\"visible_triple_ridge\":" << total_triple
              << ",\"literal_pure_5_110\":" << total_pure
              << ",\"ridge_valid\":" << total_ridge
              << "},\n  \"survivors\":[";
    for (int index = 0; index < static_cast<int>(survivors.size()); ++index) {
      if (index) std::cout << ',';
      std::cout << "{\"cut_key\":";
      PrintIntArray(survivors[index].cut_key);
      std::cout << ",\"missing_five_masks\":";
      PrintMaskArray(survivors[index].rows);
      std::cout << ",\"canonical_family_key\":";
      PrintIntVector(canonical_keys[index]);
      std::cout << '}';
    }
    std::cout << "],\n  \"canonical_orbit_count\":"
              << canonical_orbits.size()
              << ",\n  \"canonical_orbits\":[";
    int orbit_index = 0;
    for (const auto& [key, count] : canonical_orbits) {
      if (orbit_index++) std::cout << ',';
      std::cout << "{\"key\":";
      PrintIntVector(key);
      std::cout << ",\"multiplicity\":" << count << '}';
    }
    std::cout << "],\n  \"compact_only\":true,\n"
              << "  \"theorem_certificate\":false,\n"
              << "  \"whole_shadow_compact_coxeter_theorem\":false,\n"
              << "  \"finite_volume_noncompact_chamber_excluded\":false,\n"
              << "  \"scope\":\"Corrected visible-transversal and literal-blocker incidence audit only; topology and Coxeter labels are separate.\"\n"
              << "}\n";
  } catch (const std::exception& error) {
    std::cerr << "audit failed: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
