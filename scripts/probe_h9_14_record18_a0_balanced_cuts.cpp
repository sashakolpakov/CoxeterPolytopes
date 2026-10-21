// Exact exploratory census for the reduced a=0 branch of H9/14 record 18.
//
// The structural reduction writes the eight missing five-sets H as the
// union of the shores of two unordered balanced cuts, one on R union I and
// one on L union I, followed by four or five rows meeting both L and R.
// This program enumerates the nine cut-pair orbits, uses exact set-cover
// search to remove all order-at-most-two transversals visible in the
// both-centres graph-cover sector, and then applies the literal full-clutter
// blocker-purity and ridge-degree gates.  Transversals containing all of L or
// all of R are deliberately not rejected merely inside B(H), because their
// both-centres extensions are nonminimal.  It is deliberately a standalone
// probe and makes no Coxeter or finite-volume claim.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>
#include <unordered_set>
#include <vector>

namespace {

using Mask = std::uint16_t;
constexpr int kN = 14;
constexpr Mask kFull = (Mask{1} << kN) - 1;
constexpr Mask kW = (Mask{1} << 12) - 1;
constexpr Mask kL = (Mask{1} << 0) | (Mask{1} << 1);
constexpr Mask kR = (Mask{1} << 2) | (Mask{1} << 3);
constexpr Mask kI = kW ^ kL ^ kR;
constexpr std::array<Mask, 5> kEdges = {
    Mask((Mask{1} << 0) | (Mask{1} << 12)),
    Mask((Mask{1} << 1) | (Mask{1} << 12)),
    Mask((Mask{1} << 2) | (Mask{1} << 13)),
    Mask((Mask{1} << 3) | (Mask{1} << 13)),
    Mask((Mask{1} << 12) | (Mask{1} << 13)),
};

struct CutPair {
  Mask a;
  Mask b;
  Mask c;
  Mask d;
  std::array<int, 6> key;
};

struct Census {
  std::uint64_t nodes = 0;
  std::uint64_t cap_dead = 0;
  std::uint64_t cover_dead = 0;
  std::uint64_t no_small = 0;
  std::uint64_t triple_gate = 0;
  std::uint64_t visible_h = 0;
  std::uint64_t pure_full = 0;
  std::uint64_t ridge = 0;
  std::map<std::map<int, int>, std::uint64_t> full_blocker_histograms;
};

bool HitsAll(Mask value, const std::vector<Mask>& clutter) {
  return std::ranges::all_of(clutter,
      [value](Mask row) { return (value & row) != 0; });
}

bool VisibleWithBothCenters(Mask value) {
  return (value & kL) != kL && (value & kR) != kR;
}

std::vector<Mask> MinimalTransversals(const std::vector<Mask>& clutter,
                                      int vertex_count) {
  std::vector<Mask> answer;
  const Mask limit = (Mask{1} << vertex_count) - 1;
  for (Mask value = 1; value <= limit; ++value) {
    if (!HitsAll(value, clutter)) continue;
    bool minimal = true;
    Mask remaining = value;
    while (remaining != 0) {
      const Mask bit = Mask(remaining & Mask(-int(remaining)));
      if (HitsAll(Mask(value ^ bit), clutter)) {
        minimal = false;
        break;
      }
      remaining ^= bit;
    }
    if (minimal) answer.push_back(value);
  }
  return answer;
}

std::array<int, 6> OrbitKey(Mask a, Mask c) {
  const int r = std::popcount(unsigned(a & kR));
  const int l = std::popcount(unsigned(c & kL));
  std::array<int, 4> counts{};
  for (int vertex = 4; vertex < 12; ++vertex) {
    const int code = 2 * ((a >> vertex) & 1) + ((c >> vertex) & 1);
    ++counts[code];
  }
  std::array<int, 6> best{};
  bool initialized = false;
  for (int flip_a = 0; flip_a < 2; ++flip_a) {
    for (int flip_c = 0; flip_c < 2; ++flip_c) {
      std::array<int, 4> transformed{};
      for (int ia = 0; ia < 2; ++ia) {
        for (int ic = 0; ic < 2; ++ic) {
          const int old_code = 2 * (ia ^ flip_a) + (ic ^ flip_c);
          const int new_code = 2 * ia + ic;
          transformed[new_code] = counts[old_code];
        }
      }
      const int rr = flip_a ? 2 - r : r;
      const int ll = flip_c ? 2 - l : l;
      const std::array<int, 6> direct = {
          rr, ll, transformed[0], transformed[1],
          transformed[2], transformed[3]};
      const std::array<int, 6> swapped = {
          ll, rr, transformed[0], transformed[2],
          transformed[1], transformed[3]};
      for (const auto& candidate : {direct, swapped}) {
        if (!initialized || candidate < best) {
          best = candidate;
          initialized = true;
        }
      }
    }
  }
  return best;
}

std::vector<CutPair> CutOrbitRepresentatives() {
  std::vector<std::pair<Mask, Mask>> left;
  std::vector<std::pair<Mask, Mask>> right;
  const Mask left_ground = kR | kI;
  const Mask right_ground = kL | kI;
  for (Mask row = 0; row <= kW; ++row) {
    if (std::popcount(unsigned(row)) != 5) continue;
    if ((row & ~left_ground) == 0) {
      const Mask complement = left_ground ^ row;
      if (row < complement) left.emplace_back(row, complement);
    }
    if ((row & ~right_ground) == 0) {
      const Mask complement = right_ground ^ row;
      if (row < complement) right.emplace_back(row, complement);
    }
  }
  std::map<std::array<int, 6>, CutPair> representatives;
  for (const auto& [a, b] : left) {
    for (const auto& [c, d] : right) {
      const auto key = OrbitKey(a, c);
      representatives.try_emplace(key, CutPair{a, b, c, d, key});
    }
  }
  std::vector<CutPair> answer;
  for (const auto& [key, value] : representatives) answer.push_back(value);
  return answer;
}

class Search {
 public:
  explicit Search(CutPair cuts) : cuts_(cuts) {
    base_ = {cuts.a, cuts.b, cuts.c, cuts.d};
    std::sort(base_.begin(), base_.end());
    base_.erase(std::unique(base_.begin(), base_.end()), base_.end());
    need_ = 8 - int(base_.size());
    for (Mask row = 0; row <= kW; ++row) {
      if (std::popcount(unsigned(row)) != 5) continue;
      if ((row & kL) == 0 || (row & kR) == 0) continue;
      candidates_.push_back(row);
    }
    for (int vertex = 0; vertex < 12; ++vertex) {
      int degree = 0;
      for (Mask row : base_) degree += (row >> vertex) & 1;
      residual_caps_[vertex] = (vertex < 4 ? 5 : 6) - degree;
      if (residual_caps_[vertex] < 0) {
        throw std::runtime_error("balanced-cut base exceeds local cap");
      }
    }
    for (Mask value = 1; value <= kW; ++value) {
      if (std::popcount(unsigned(value)) <= 2 && HitsAll(value, base_) &&
          VisibleWithBothCenters(value)) {
        small_.push_back(value);
      }
    }
    if (small_.size() > 64) throw std::runtime_error("coverage escaped one word");
    covers_.resize(candidates_.size());
    triple_hits_.resize(candidates_.size());
    coverers_.resize(small_.size());
    for (int row = 0; row < int(candidates_.size()); ++row) {
      for (int constraint = 0; constraint < int(small_.size()); ++constraint) {
        if ((candidates_[row] & small_[constraint]) == 0) {
          covers_[row] |= std::uint64_t{1} << constraint;
          coverers_[constraint].push_back(row);
        }
      }
    }
    for (Mask triple = 0; triple <= kW; ++triple) {
      if (std::popcount(unsigned(triple)) == 3 && HitsAll(triple, base_) &&
          VisibleWithBothCenters(triple)) {
        base_triples_.push_back(triple);
      }
    }
    if (base_triples_.size() > 256) {
      throw std::runtime_error("triple table escaped four words");
    }
    for (int row = 0; row < int(candidates_.size()); ++row) {
      triple_hits_[row].fill(0);
      for (int index = 0; index < int(base_triples_.size()); ++index) {
        if (candidates_[row] & base_triples_[index]) {
          triple_hits_[row][index / 64] |= std::uint64_t{1} << (index % 64);
        }
      }
    }
    triple_full_.fill(~std::uint64_t{0});
    if (base_triples_.size() % 64) {
      triple_full_[(base_triples_.size() - 1) / 64] =
          (std::uint64_t{1} << (base_triples_.size() % 64)) - 1;
    }
    for (int block = (base_triples_.size() + 63) / 64; block < 4; ++block) {
      triple_full_[block] = 0;
    }
    std::vector<Mask> xz;
    for (Mask u = 0; u < 12; ++u) {
      if (!((cuts_.c >> u) & 1)) continue;
      for (Mask v = 0; v < 12; ++v) {
        if (!((cuts_.d >> v) & 1)) continue;
        xz.push_back(Mask((Mask{1} << 2) | (Mask{1} << 3) |
                          (Mask{1} << 12) | (Mask{1} << u) | (Mask{1} << v)));
      }
    }
    for (Mask u = 0; u < 12; ++u) {
      if (!((cuts_.a >> u) & 1)) continue;
      for (Mask v = 0; v < 12; ++v) {
        if (!((cuts_.b >> v) & 1)) continue;
        xz.push_back(Mask((Mask{1} << 0) | (Mask{1} << 1) |
                          (Mask{1} << 13) | (Mask{1} << u) | (Mask{1} << v)));
      }
    }
    for (Mask q = 0; q <= kW; ++q) {
      if (std::popcount(unsigned(q)) != 4) continue;
      const Mask six = Mask(q | (Mask{1} << 12) | (Mask{1} << 13));
      int fixed_degree = 0;
      for (Mask complement : xz) {
        if ((complement & six) == complement) ++fixed_degree;
      }
      if (fixed_degree > 2) {
        cut_impossible_ = true;
        continue;
      }
      q_constraints_.push_back({q, fixed_degree});
    }
    full_cover_ = small_.size() == 64
        ? ~std::uint64_t{0}
        : (std::uint64_t{1} << small_.size()) - 1;
    seen_.resize(need_ + 1);
  }

  Census Run() {
    if (!cut_impossible_) Explore({}, 0, {}, triple_full_);
    return census_;
  }

  int BaseSize() const { return base_.size(); }
  int Need() const { return need_; }
  int CandidateCount() const { return candidates_.size(); }
  int SmallCount() const { return small_.size(); }

 private:
  bool CapFeasible(const std::array<unsigned char, 12>& degrees,
                   Mask row) const {
    for (int vertex = 0; vertex < 12; ++vertex) {
      if (((row >> vertex) & 1) && degrees[vertex] >= residual_caps_[vertex]) {
        return false;
      }
    }
    return true;
  }

  std::array<unsigned char, 12> Add(
      std::array<unsigned char, 12> degrees, Mask row) const {
    for (int vertex = 0; vertex < 12; ++vertex) {
      degrees[vertex] += (row >> vertex) & 1;
    }
    return degrees;
  }

  std::uint64_t SelectedKey(std::vector<int> selected) const {
    std::sort(selected.begin(), selected.end());
    std::uint64_t key = selected.size();
    for (int index : selected) key = (key << 9) | index;
    return key;
  }

  int TripleCount(const std::array<std::uint64_t, 4>& triples) const {
    int answer = 0;
    for (auto word : triples) answer += std::popcount(word);
    return answer;
  }

  int YDegree(Mask q,
              const std::array<std::uint64_t, 4>& triples) const {
    int answer = 0;
    for (int index = 0; index < int(base_triples_.size()); ++index) {
      if (!((triples[index / 64] >> (index % 64)) & 1)) continue;
      if ((base_triples_[index] & q) == base_triples_[index]) ++answer;
    }
    return answer;
  }

  bool TripleGate(const std::array<std::uint64_t, 4>& triples,
                  bool final) const {
    if (TripleCount(triples) < 60) return false;
    for (const auto& [q, fixed_degree] : q_constraints_) {
      const int y_degree = YDegree(q, triples);
      if (!final) {
        if (fixed_degree == 1 && y_degree == 0) return false;
        continue;
      }
      const int degree = fixed_degree + y_degree;
      if (degree != 0 && degree != 2) return false;
    }
    return !final || TripleCount(triples) == 60;
  }

  void Audit(const std::vector<int>& selected,
             const std::array<std::uint64_t, 4>& triples) {
    ++census_.no_small;
    if (!TripleGate(triples, true)) return;
    ++census_.triple_gate;
    auto h = base_;
    for (int index : selected) h.push_back(candidates_[index]);
    const auto h_blocker = MinimalTransversals(h, 12);
    std::map<int, int> h_histogram;
    for (Mask value : h_blocker) {
      if (VisibleWithBothCenters(value)) {
        ++h_histogram[std::popcount(unsigned(value))];
      }
    }
    if (h_histogram != std::map<int, int>{{3, 60}}) return;
    ++census_.visible_h;

    std::vector<Mask> clutter(kEdges.begin(), kEdges.end());
    clutter.insert(clutter.end(), h.begin(), h.end());
    const auto blocker = MinimalTransversals(clutter, 14);
    std::map<int, int> histogram;
    for (Mask value : blocker) ++histogram[std::popcount(unsigned(value))];
    ++census_.full_blocker_histograms[histogram];
    if (histogram != std::map<int, int>{{5, 110}}) return;
    ++census_.pure_full;

    std::map<Mask, int> ridges;
    for (Mask complement : blocker) {
      Mask outside = kFull ^ complement;
      while (outside != 0) {
        const Mask bit = Mask(outside & Mask(-int(outside)));
        ++ridges[Mask(complement | bit)];
        outside ^= bit;
      }
    }
    if (std::ranges::all_of(ridges,
          [](const auto& item) { return item.second == 2; })) {
      ++census_.ridge;
      std::cout << "SURVIVOR key=";
      for (int value : cuts_.key) std::cout << value << ',';
      std::cout << " rows=";
      for (Mask row : h) std::cout << row << ',';
      std::cout << '\n';
    }
  }

  void Explore(std::vector<int> selected, std::uint64_t covered,
               std::array<unsigned char, 12> degrees,
               std::array<std::uint64_t, 4> triples) {
    ++census_.nodes;
    const auto state_key = SelectedKey(selected);
    if (!seen_[selected.size()].insert(state_key).second) return;
    if (int(selected.size()) == need_) {
      if (covered == full_cover_) Audit(selected, triples);
      return;
    }
    if (!TripleGate(triples, false)) return;

    std::vector<int> best;
    std::uint64_t feasible_union = covered;
    bool uncovered = false;
    for (int constraint = 0; constraint < int(small_.size()); ++constraint) {
      if ((covered >> constraint) & 1) continue;
      uncovered = true;
      std::vector<int> feasible;
      for (int candidate : coverers_[constraint]) {
        if (std::find(selected.begin(), selected.end(), candidate) != selected.end()) continue;
        if (!CapFeasible(degrees, candidates_[candidate])) continue;
        feasible.push_back(candidate);
        feasible_union |= covers_[candidate];
      }
      if (feasible.empty()) {
        ++census_.cover_dead;
        return;
      }
      if (best.empty() || feasible.size() < best.size()) best = std::move(feasible);
    }
    if (feasible_union != full_cover_) {
      ++census_.cover_dead;
      return;
    }
    auto consider_reducer_constraint = [&](Mask q, int current_degree) {
      std::vector<int> reducers;
      for (int candidate = 0; candidate < int(candidates_.size()); ++candidate) {
        if (std::find(selected.begin(), selected.end(), candidate) != selected.end()) continue;
        if (!CapFeasible(degrees, candidates_[candidate])) continue;
        auto next_triples = triples;
        for (int block = 0; block < 4; ++block) {
          next_triples[block] &= triple_hits_[candidate][block];
        }
        if (YDegree(q, next_triples) < current_degree) reducers.push_back(candidate);
      }
      if (reducers.empty()) {
        best.clear();
        return false;
      }
      if (best.empty() || reducers.size() < best.size()) best = std::move(reducers);
      return true;
    };
    for (const auto& [q, fixed_degree] : q_constraints_) {
      const int y_degree = YDegree(q, triples);
      const bool must_reduce =
          (fixed_degree == 2 && y_degree > 0) ||
          (fixed_degree == 1 && y_degree > 1) ||
          (fixed_degree == 0 && (y_degree == 1 || y_degree > 2));
      if (must_reduce) {
        if (!consider_reducer_constraint(q, y_degree)) return;
        break;
      }
    }
    if (TripleCount(triples) > 60) {
      std::vector<int> reducers;
      const int current_count = TripleCount(triples);
      for (int candidate = 0; candidate < int(candidates_.size()); ++candidate) {
        if (std::find(selected.begin(), selected.end(), candidate) != selected.end()) continue;
        if (!CapFeasible(degrees, candidates_[candidate])) continue;
        auto next_triples = triples;
        for (int block = 0; block < 4; ++block) {
          next_triples[block] &= triple_hits_[candidate][block];
        }
        if (TripleCount(next_triples) < current_count) reducers.push_back(candidate);
      }
      if (reducers.empty()) return;
      if (best.empty() || reducers.size() < best.size()) best = std::move(reducers);
    }
    if (!uncovered && best.empty()) {
      for (int candidate = 0; candidate < int(candidates_.size()); ++candidate) {
        if (std::find(selected.begin(), selected.end(), candidate) != selected.end()) continue;
        if (CapFeasible(degrees, candidates_[candidate])) best.push_back(candidate);
      }
    }
    if (best.empty()) {
      ++census_.cap_dead;
      return;
    }
    for (int candidate : best) {
      auto next = selected;
      next.push_back(candidate);
      auto next_triples = triples;
      for (int block = 0; block < 4; ++block) {
        next_triples[block] &= triple_hits_[candidate][block];
      }
      Explore(std::move(next), covered | covers_[candidate],
              Add(degrees, candidates_[candidate]), next_triples);
    }
  }

  CutPair cuts_;
  std::vector<Mask> base_;
  std::vector<Mask> candidates_;
  std::vector<Mask> small_;
  std::vector<std::uint64_t> covers_;
  std::vector<std::array<std::uint64_t, 4>> triple_hits_;
  std::vector<std::vector<int>> coverers_;
  std::array<int, 12> residual_caps_{};
  std::vector<Mask> base_triples_;
  std::array<std::uint64_t, 4> triple_full_{};
  std::vector<std::pair<Mask, int>> q_constraints_;
  bool cut_impossible_ = false;
  int need_ = 0;
  std::uint64_t full_cover_ = 0;
  std::vector<std::unordered_set<std::uint64_t>> seen_;
  Census census_;
};

}  // namespace

int main() {
  try {
    std::cout.setf(std::ios::unitbuf);
    const auto representatives = CutOrbitRepresentatives();
    if (representatives.size() != 9) {
      throw std::runtime_error("cut-pair orbit count is not nine");
    }
    std::cout << "cut_pair_orbits " << representatives.size() << '\n';
    Census total;
    int index = 0;
    for (const auto& representative : representatives) {
      Search search(representative);
      const auto result = search.Run();
      std::cout << "orbit " << index++ << " key";
      for (int value : representative.key) std::cout << ' ' << value;
      std::cout << " base " << search.BaseSize()
                << " need " << search.Need()
                << " candidates " << search.CandidateCount()
                << " small " << search.SmallCount()
                << " nodes " << result.nodes
                << " no_small " << result.no_small
                << " triple_gate " << result.triple_gate
                << " visible_h " << result.visible_h
                << " pure_full " << result.pure_full
                << " ridge " << result.ridge << '\n';
      total.nodes += result.nodes;
      total.no_small += result.no_small;
      total.triple_gate += result.triple_gate;
      total.visible_h += result.visible_h;
      total.pure_full += result.pure_full;
      total.ridge += result.ridge;
    }
    std::cout << "total nodes " << total.nodes
              << " no_small " << total.no_small
              << " triple_gate " << total.triple_gate
              << " visible_h " << total.visible_h
              << " pure_full " << total.pure_full
              << " ridge " << total.ridge << '\n';
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
