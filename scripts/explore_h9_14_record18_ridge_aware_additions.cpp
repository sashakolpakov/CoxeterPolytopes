// Lossless ridge-aware addition search for the positive-A4 branches of
// H9/14 record 18.
//
// Usage: program a A4mask ... A4mask
//
// This translation unit reuses the exact local interval and blocker helpers
// of the cover-only driver, but replaces its search by a stronger monotone
// search.  Besides eliminating every base transversal of order at most four,
// it tracks every surviving order-five transversal.  A valid completion must
// retain the scheduled number of them and every six-set must contain zero or
// two.  Both conditions are monotone under adding missing five-sets, so an
// offending ridge or excess facet supplies a lossless next-branch constraint.
// Every terminal target is still checked by literal blocker dualization.

#define main record18_cover_only_main
#include "explore_h9_14_record18_addition_cover.cpp"
#undef main

namespace {

#ifndef H9_14_ADDITION_FIRST_ISOLATE
#define H9_14_ADDITION_FIRST_ISOLATE 4
#endif
#ifndef H9_14_RIDGE_KIND
#define H9_14_RIDGE_KIND "h9_14_record18_ridge_aware_additions"
#endif

std::vector<std::uint16_t> g_base_five_transversals;
std::vector<Word> g_candidate_five_hits;
std::vector<std::vector<int>> g_six_set_members;
std::vector<std::unordered_set<std::uint64_t>> g_ridge_seen;
std::vector<std::unordered_set<std::string>> g_symmetry_seen;
std::set<std::vector<int>> g_ridge_solutions;
std::set<std::vector<std::uint16_t>> g_expanded_ridge_solutions;
std::uint64_t g_ridge_nodes = 0;
std::uint64_t g_dead_facet_lower = 0;
std::uint64_t g_dead_ridge = 0;
std::uint64_t g_dead_exact_target = 0;
bool g_symmetry_quotient = false;

using Move = std::array<int, kN>;
using CandidateWord = std::array<std::uint64_t, 13>;
std::vector<Move> g_fixed_actions;
int g_fixed_count = 0;
int g_isolate_start = 0;
std::vector<CandidateWord> g_small_coverer_words;
std::vector<CandidateWord> g_facet_killers;

int WordCount(const Word& word) {
  int answer = 0;
  for (const auto block : word) answer += std::popcount(block);
  return answer;
}

int CandidateWordCount(const CandidateWord& word) {
  int answer = 0;
  for (const auto block : word) answer += std::popcount(block);
  return answer;
}

void CandidateWordSet(CandidateWord* word, int index) {
  (*word)[index / 64] |= std::uint64_t{1} << (index % 64);
}

CandidateWord CandidateWordAnd(CandidateWord left,
                               const CandidateWord& right) {
  for (int block = 0; block < static_cast<int>(left.size()); ++block) {
    left[block] &= right[block];
  }
  return left;
}

void CandidateWordOr(CandidateWord* target, const CandidateWord& source) {
  for (int block = 0; block < static_cast<int>(target->size()); ++block) {
    (*target)[block] |= source[block];
  }
}

std::vector<int> CandidateWordIndices(const CandidateWord& word) {
  std::vector<int> answer;
  for (int block = 0; block < static_cast<int>(word.size()); ++block) {
    auto value = word[block];
    while (value) {
      const int offset = std::countr_zero(value);
      const int index = 64 * block + offset;
      if (index < static_cast<int>(g_candidates.size())) answer.push_back(index);
      value &= value - 1;
    }
  }
  return answer;
}

bool WordBit(const Word& word, int index) {
  return (word[index / 64] >> (index % 64)) & 1;
}

Word WordAnd(Word left, const Word& right) {
  for (int block = 0; block < static_cast<int>(left.size()); ++block) {
    left[block] &= right[block];
  }
  return left;
}

int RidgeDegree(const Word& survivors, const std::vector<int>& members) {
  int answer = 0;
  for (const int member : members) answer += WordBit(survivors, member);
  return answer;
}

bool SelectedContains(const std::vector<int>& selected, int candidate) {
  return std::find(selected.begin(), selected.end(), candidate) != selected.end();
}

CandidateWord FeasibleCandidateWord(
    const std::vector<int>& selected,
    const std::array<unsigned char, kN>& degrees) {
  CandidateWord answer{};
  for (int candidate = 0; candidate < static_cast<int>(g_candidates.size());
       ++candidate) {
    if (!SelectedContains(selected, candidate) &&
        CapFeasible(degrees, g_candidates[candidate])) {
      CandidateWordSet(&answer, candidate);
    }
  }
  return answer;
}

CandidateWord RidgeReducers(const Word& survivors,
                            const std::vector<int>& members) {
  CandidateWord answer{};
  for (const int facet : members) {
    if (WordBit(survivors, facet)) {
      CandidateWordOr(&answer, g_facet_killers[facet]);
    }
  }
  return answer;
}

void KeepSmallerCandidateConstraint(CandidateWord* best, bool* has_best,
                                    const CandidateWord& candidate) {
  const int size = CandidateWordCount(candidate);
  if (size != 0 && (!*has_best || size < CandidateWordCount(*best))) {
    *best = candidate;
    *has_best = true;
  }
}

Move IdentityMove() {
  Move answer{};
  for (int vertex = 0; vertex < kN; ++vertex) answer[vertex] = vertex;
  return answer;
}

Move SwapMove(int first, int second) {
  Move answer = IdentityMove();
  std::swap(answer[first], answer[second]);
  return answer;
}

#ifndef H9_14_RECORD14_RIDGE_KEY
Move ExchangeLeafPairsMove() {
  Move answer = IdentityMove();
  std::swap(answer[0], answer[2]);
  std::swap(answer[1], answer[3]);
  std::swap(answer[12], answer[13]);
  return answer;
}
#endif

Move ComposeMove(const Move& left, const Move& right) {
  Move answer{};
  for (int vertex = 0; vertex < kN; ++vertex) {
    answer[vertex] = left[right[vertex]];
  }
  return answer;
}

std::vector<Move> GeneratedMoves(const std::vector<Move>& generators) {
  std::set<Move> found{IdentityMove()};
  std::vector<Move> frontier{IdentityMove()};
  while (!frontier.empty()) {
    const auto current = frontier.back();
    frontier.pop_back();
    for (const auto& generator : generators) {
      for (const auto& composite : {ComposeMove(current, generator),
                                    ComposeMove(generator, current)}) {
        if (found.insert(composite).second) frontier.push_back(composite);
      }
    }
  }
  return {found.begin(), found.end()};
}

std::uint16_t MoveMask(std::uint16_t value, const Move& move) {
  std::uint16_t answer = 0;
  for (int source = 0; source < kN; ++source) {
    if ((value >> source) & 1u) answer |= std::uint16_t{1} << move[source];
  }
  return answer;
}

int MovedFixedMask(std::uint16_t row, const Move& move) {
#ifdef H9_14_RECORD14_RIDGE_KEY
  const auto moved = MoveMask(row, move);
  int answer = 0;
  int position = 0;
  for (const int vertex : {0, 1, 2, 12, 13}) {
    answer |= int((moved >> vertex) & 1u) << position++;
  }
  for (int vertex = H9_14_ADDITION_FIRST_ISOLATE;
       vertex < g_isolate_start; ++vertex) {
    answer |= int((moved >> vertex) & 1u) << position++;
  }
  return answer;
#else
  int answer = 0;
  for (int source = 0; source < g_fixed_count; ++source) {
    if ((row >> source) & 1u) answer |= 1 << move[source];
  }
  return answer;
#endif
}

void AppendFixedMask(std::string* key, int value) {
#ifdef H9_14_RECORD14_RIDGE_KEY
  key->push_back(static_cast<char>(value & 255));
  key->push_back(static_cast<char>((value >> 8) & 255));
#else
  key->push_back(static_cast<char>(value));
#endif
}

void EvaluateCanonicalTieOrders(
    const std::vector<std::vector<int>>& blocks, int block,
    std::vector<int>* order, const Move& move, std::string* best,
    bool* first) {
  if (block != static_cast<int>(blocks.size())) {
    auto values = blocks[block];
    std::sort(values.begin(), values.end());
    do {
      const auto old_size = order->size();
      order->insert(order->end(), values.begin(), values.end());
      EvaluateCanonicalTieOrders(blocks, block + 1, order, move, best, first);
      order->resize(old_size);
    } while (std::next_permutation(values.begin(), values.end()));
    return;
  }

  std::string key;
  key.reserve(order->size() + 12 - g_isolate_start);
  for (const int candidate : *order) {
    AppendFixedMask(&key, MovedFixedMask(g_candidates[candidate], move));
  }
  std::vector<int> columns;
  for (int isolate = g_isolate_start; isolate < 12; ++isolate) {
    int column = 0;
    for (int position = 0; position < static_cast<int>(order->size());
         ++position) {
      column |= int((g_candidates[(*order)[position]] >> isolate) & 1u)
                << position;
    }
    columns.push_back(column);
  }
  std::sort(columns.begin(), columns.end());
  for (const int column : columns) key.push_back(static_cast<char>(column));
  if (*first || key < *best) {
    *best = std::move(key);
    *first = false;
  }
}

std::string CanonicalFamilyKey(const std::vector<int>& selected) {
  if (selected.empty()) return {};
  std::string best;
  bool first = true;
  for (const auto& move : g_fixed_actions) {
    std::map<int, std::vector<int>> by_fixed_mask;
    for (const int candidate : selected) {
      by_fixed_mask[MovedFixedMask(g_candidates[candidate], move)].push_back(
          candidate);
    }
    std::vector<std::vector<int>> blocks;
    for (auto& [fixed_mask, values] : by_fixed_mask) {
      static_cast<void>(fixed_mask);
      blocks.push_back(std::move(values));
    }
    std::vector<int> order;
    EvaluateCanonicalTieOrders(blocks, 0, &order, move, &best, &first);
  }
  return best;
}

bool SeenBefore(const std::vector<int>& selected) {
  const auto depth = selected.size();
  if (g_symmetry_quotient) {
#ifndef H9_14_RECORD14_RIDGE_KEY
    // The same labelled set is often reached in several row orders.  Reject
    // those cheap duplicates before constructing the more expensive exact
    // incidence-matrix orbit key.
    if (!g_ridge_seen[depth].insert(ChosenKey(selected)).second) return true;
#endif
    return !g_symmetry_seen[depth].insert(CanonicalFamilyKey(selected)).second;
  }
  return !g_ridge_seen[depth].insert(ChosenKey(selected)).second;
}

void AuditRidgeTarget(const std::vector<int>& selected) {
  auto clutter = g_base;
  for (const int index : selected) clutter.push_back(g_candidates[index]);
  const auto blocker = MinimalTransversals(clutter);
  std::map<int, int> sizes;
  for (const auto value : blocker) ++sizes[std::popcount(value)];
  if (sizes == std::map<int, int>{{5, g_scheduled_facets}}) {
    g_ridge_solutions.insert(selected);
  }
}

void ConfigureSymmetryQuotient() {
  if (g_a < 1 || g_a > 2 || static_cast<int>(g_fours.size()) != g_a) {
    throw std::runtime_error(
        "the collision-free quotient/expansion mode supports a=1 or a=2");
  }
  int largest_used_isolate = H9_14_ADDITION_FIRST_ISOLATE - 1;
  std::set<int> used_isolates;
  for (const auto four : g_fours) {
    for (int vertex = H9_14_ADDITION_FIRST_ISOLATE; vertex < 12; ++vertex) {
      if ((four >> vertex) & 1u) {
        used_isolates.insert(vertex);
        largest_used_isolate = std::max(largest_used_isolate, vertex);
      }
    }
  }
  g_fixed_count = largest_used_isolate + 1;
  g_isolate_start = g_fixed_count;
  for (int vertex = H9_14_ADDITION_FIRST_ISOLATE;
       vertex < g_isolate_start; ++vertex) {
    if (!used_isolates.contains(vertex)) {
      throw std::runtime_error(
          "representative has a gap in its used-isolate labels");
    }
  }

#ifdef H9_14_RECORD14_RIDGE_KEY
  const auto active_actions = GeneratedMoves({SwapMove(0, 1)});
#else
  const auto active_actions = GeneratedMoves(
      {SwapMove(0, 1), SwapMove(2, 3), ExchangeLeafPairsMove()});
#endif
  std::vector<int> isolate_images(used_isolates.begin(),
                                  used_isolates.end());
  const std::set<std::uint16_t> target_fours(g_fours.begin(), g_fours.end());
  do {
    for (const auto& active : active_actions) {
      Move move = active;
      int index = 0;
      for (const int source : used_isolates) {
        move[source] = isolate_images[index++];
      }
      std::set<std::uint16_t> moved_fours;
      for (const auto four : g_fours) {
        moved_fours.insert(MoveMask(four, move));
      }
      if (moved_fours == target_fours) g_fixed_actions.push_back(move);
    }
  } while (std::next_permutation(isolate_images.begin(),
                                 isolate_images.end()));
  std::ranges::sort(g_fixed_actions);
  g_fixed_actions.erase(
      std::unique(g_fixed_actions.begin(), g_fixed_actions.end()),
      g_fixed_actions.end());
  if (g_fixed_actions.empty()) {
    throw std::runtime_error("failed to derive a nonempty A4 stabilizer");
  }
}

void ValidateOneMove(const Move& move) {
  std::set<int> images(move.begin(), move.end());
  if (images.size() != kN || *images.begin() != 0 ||
      *images.rbegin() != kN - 1) {
    throw std::runtime_error("declared symmetry is not a permutation");
  }
  const auto moved_set = [&move](const auto& values) {
    std::set<std::uint16_t> answer;
    for (const auto value : values) answer.insert(MoveMask(value, move));
    return answer;
  };
  if (moved_set(kEdges) !=
      std::set<std::uint16_t>(kEdges.begin(), kEdges.end())) {
    throw std::runtime_error("declared symmetry does not preserve the graph");
  }
  if (moved_set(g_fours) !=
      std::set<std::uint16_t>(g_fours.begin(), g_fours.end()) ||
      moved_set(g_shadow) !=
          std::set<std::uint16_t>(g_shadow.begin(), g_shadow.end()) ||
      moved_set(g_candidates) !=
          std::set<std::uint16_t>(g_candidates.begin(), g_candidates.end()) ||
      moved_set(g_small) !=
          std::set<std::uint16_t>(g_small.begin(), g_small.end()) ||
      moved_set(g_base_five_transversals) !=
          std::set<std::uint16_t>(g_base_five_transversals.begin(),
                                  g_base_five_transversals.end())) {
    throw std::runtime_error(
        "declared symmetry does not preserve a derived incidence family");
  }
  for (int source = 0; source < kN; ++source) {
    if (g_lower[source] != g_lower[move[source]] ||
        g_upper[source] != g_upper[move[source]]) {
      throw std::runtime_error(
          "declared symmetry does not preserve the local intervals");
    }
  }
}

void ValidateSymmetryQuotient() {
  for (const auto& move : g_fixed_actions) ValidateOneMove(move);
  if (12 - g_isolate_start >= 2) {
    ValidateOneMove(SwapMove(g_isolate_start, g_isolate_start + 1));
  }
}

void ExpandCanonicalTargets() {
  std::set<std::uint16_t> candidate_rows(g_candidates.begin(),
                                         g_candidates.end());
  std::vector<int> isolate_images;
  for (int vertex = g_isolate_start; vertex < 12; ++vertex) {
    isolate_images.push_back(vertex);
  }
  for (const auto& selected : g_ridge_solutions) {
    std::sort(isolate_images.begin(), isolate_images.end());
    do {
      for (const auto& fixed_move : g_fixed_actions) {
        Move move = fixed_move;
        for (int index = 0; index < static_cast<int>(isolate_images.size());
             ++index) {
          move[g_isolate_start + index] = isolate_images[index];
        }
        std::vector<std::uint16_t> family;
        for (const int candidate : selected) {
          const auto row = MoveMask(g_candidates[candidate], move);
          if (!candidate_rows.contains(row)) {
            throw std::runtime_error(
                "target expansion left the candidate row family");
          }
          family.push_back(row);
        }
        std::sort(family.begin(), family.end());
        if (std::adjacent_find(family.begin(), family.end()) != family.end()) {
          throw std::runtime_error("target expansion merged distinct rows");
        }
        g_expanded_ridge_solutions.insert(std::move(family));
      }
    } while (std::next_permutation(isolate_images.begin(),
                                   isolate_images.end()));
  }
}

std::uint64_t AuditExpandedTargets() {
  std::uint64_t answer = 0;
  for (const auto& family : g_expanded_ridge_solutions) {
    auto clutter = g_base;
    clutter.insert(clutter.end(), family.begin(), family.end());
    const auto blocker = MinimalTransversals(clutter);
    std::map<int, int> sizes;
    for (const auto value : blocker) ++sizes[std::popcount(value)];
    if (sizes != std::map<int, int>{{5, g_scheduled_facets}}) {
      throw std::runtime_error(
          "an expanded target failed literal blocker dualization");
    }
    for (std::uint16_t six = 0; six <= kFull; ++six) {
      if (std::popcount(six) != 6) continue;
      int degree = 0;
      for (const auto facet : blocker) {
        degree += (facet & six) == facet;
      }
      if (degree != 0 && degree != 2) {
        throw std::runtime_error("an expanded target failed the literal ridge law");
      }
    }
    ++answer;
    if (answer % 100 == 0) {
      std::cerr << "expanded_literal_targets " << answer << '\n';
    }
  }
  return answer;
}

void RidgeSearch(std::vector<int> selected, Word covered,
                 std::array<unsigned char, kN> degrees, Word survivors) {
  ++g_ridge_nodes;
  if (g_ridge_nodes % 10'000'000 == 0) {
    std::cerr << "ridge_nodes " << g_ridge_nodes << " depth "
              << selected.size() << '\n';
  }
  const int remaining = g_addition_count - static_cast<int>(selected.size());
  for (int vertex = 0; vertex < kN; ++vertex) {
    if (degrees[vertex] + remaining < g_lower[vertex]) {
      ++g_dead_lower;
      return;
    }
  }
  const int facet_count = WordCount(survivors);
  if (facet_count < g_scheduled_facets) {
    ++g_dead_facet_lower;
    return;
  }
  if (SeenBefore(selected)) return;

  Word full((g_small.size() + 63) / 64, ~std::uint64_t{0});
  if (!g_small.empty() && g_small.size() % 64) {
    full.back() = (std::uint64_t{1} << (g_small.size() % 64)) - 1;
  }
  if (selected.size() == static_cast<std::size_t>(g_addition_count)) {
    bool lower = true;
    for (int vertex = 0; vertex < kN; ++vertex) {
      lower = lower && degrees[vertex] >= g_lower[vertex];
    }
    if (!lower || !WordFull(covered, full) ||
        facet_count != g_scheduled_facets) {
      return;
    }
    bool ridges = true;
    for (const auto& members : g_six_set_members) {
      const int degree = RidgeDegree(survivors, members);
      if (degree != 0 && degree != 2) {
        ridges = false;
        break;
      }
    }
    if (ridges) {
      ++g_dead_exact_target;
      auto normalized = selected;
      std::sort(normalized.begin(), normalized.end());
      AuditRidgeTarget(normalized);
    }
    return;
  }

  const auto feasible_word = FeasibleCandidateWord(selected, degrees);
  const auto feasible_candidates = CandidateWordIndices(feasible_word);
  if (CandidateWordCount(feasible_word) < remaining) return;
  if (facet_count > g_scheduled_facets) {
    std::vector<int> optimistic_kills;
    optimistic_kills.reserve(feasible_candidates.size());
    for (const int candidate : feasible_candidates) {
      optimistic_kills.push_back(
          facet_count -
          WordCount(WordAnd(survivors, g_candidate_five_hits[candidate])));
    }
    std::ranges::sort(optimistic_kills, std::greater<>());
    int optimistic = 0;
    for (int index = 0;
         index < remaining && index < static_cast<int>(optimistic_kills.size());
         ++index) {
      optimistic += optimistic_kills[index];
    }
    if (facet_count - optimistic > g_scheduled_facets) {
      ++g_dead_facet_lower;
      return;
    }
  }

  CandidateWord best{};
  bool has_best = false;
  Word feasible_coverage = covered;
  bool has_uncovered = false;
  for (int constraint = 0; constraint < static_cast<int>(g_small.size());
       ++constraint) {
    if (WordBit(covered, constraint)) continue;
    has_uncovered = true;
    const auto choices = CandidateWordAnd(
        feasible_word, g_small_coverer_words[constraint]);
    if (CandidateWordCount(choices) == 0) {
      ++g_dead_uncoverable;
      return;
    }
    KeepSmallerCandidateConstraint(&best, &has_best, choices);
  }
  for (const int candidate : feasible_candidates) {
    WordOr(feasible_coverage, g_coverage[candidate]);
  }
  if (!WordFull(feasible_coverage, full)) {
    ++g_dead_uncoverable;
    return;
  }

  // A current degree one, or a degree larger than two, must decrease in any
  // final pseudomanifold.  Since additions only delete surviving facets, any
  // valid completion contains a remaining candidate that makes that decrease.
  for (const auto& members : g_six_set_members) {
    const int degree = RidgeDegree(survivors, members);
    if (degree != 1 && degree <= 2) continue;
    const auto choices = CandidateWordAnd(
        feasible_word, RidgeReducers(survivors, members));
    if (CandidateWordCount(choices) == 0) {
      ++g_dead_ridge;
      return;
    }
    KeepSmallerCandidateConstraint(&best, &has_best, choices);
  }

  if (facet_count > g_scheduled_facets) {
    CandidateWord reducers{};
    for (int facet = 0;
         facet < static_cast<int>(g_base_five_transversals.size()); ++facet) {
      if (WordBit(survivors, facet)) {
        CandidateWordOr(&reducers, g_facet_killers[facet]);
      }
    }
    const auto choices = CandidateWordAnd(feasible_word, reducers);
    if (CandidateWordCount(choices) == 0) {
      ++g_dead_facet_lower;
      return;
    }
    KeepSmallerCandidateConstraint(&best, &has_best, choices);
  }

  if (!has_uncovered && !has_best) {
    best = feasible_word;
    has_best = true;
  }
  if (!has_best) return;
  const auto best_candidates = CandidateWordIndices(best);
  if (remaining == 1) {
    for (const int candidate : best_candidates) {
      auto normalized = selected;
      normalized.push_back(candidate);
      std::sort(normalized.begin(), normalized.end());
      if (SeenBefore(normalized)) continue;
      ++g_ridge_nodes;
      auto next_covered = covered;
      WordOr(next_covered, g_coverage[candidate]);
      if (!WordFull(next_covered, full)) continue;
      const auto next_degrees = AddDegree(degrees, g_candidates[candidate]);
      bool lower = true;
      for (int vertex = 0; vertex < kN; ++vertex) {
        lower = lower && next_degrees[vertex] >= g_lower[vertex];
      }
      if (!lower) continue;
      const auto next_survivors =
          WordAnd(survivors, g_candidate_five_hits[candidate]);
      if (WordCount(next_survivors) != g_scheduled_facets) continue;
      bool ridges = true;
      for (const auto& members : g_six_set_members) {
        const int degree = RidgeDegree(next_survivors, members);
        if (degree != 0 && degree != 2) {
          ridges = false;
          break;
        }
      }
      if (!ridges) continue;
      ++g_dead_exact_target;
      AuditRidgeTarget(normalized);
    }
    return;
  }
  for (const int candidate : best_candidates) {
    auto next_selected = selected;
    next_selected.push_back(candidate);
    auto next_covered = covered;
    WordOr(next_covered, g_coverage[candidate]);
    RidgeSearch(
        std::move(next_selected), std::move(next_covered),
        AddDegree(degrees, g_candidates[candidate]),
        WordAnd(survivors, g_candidate_five_hits[candidate]));
  }
}

void InitializeRidgeTables() {
  for (std::uint16_t value = 0; value <= kFull; ++value) {
    if (std::popcount(value) == 5 && HitsAll(value, g_base)) {
      g_base_five_transversals.push_back(value);
    }
  }
  const int blocks = (g_base_five_transversals.size() + 63) / 64;
  g_candidate_five_hits.assign(g_candidates.size(), Word(blocks));
  for (int candidate = 0; candidate < static_cast<int>(g_candidates.size());
       ++candidate) {
    for (int facet = 0;
         facet < static_cast<int>(g_base_five_transversals.size()); ++facet) {
      if (g_candidates[candidate] & g_base_five_transversals[facet]) {
        g_candidate_five_hits[candidate][facet / 64] |=
            std::uint64_t{1} << (facet % 64);
      }
    }
  }
  g_small_coverer_words.assign(g_small.size(), CandidateWord{});
  for (int constraint = 0; constraint < static_cast<int>(g_small.size());
       ++constraint) {
    for (const int candidate : g_coverers[constraint]) {
      CandidateWordSet(&g_small_coverer_words[constraint], candidate);
    }
  }
  g_facet_killers.assign(g_base_five_transversals.size(), CandidateWord{});
  for (int facet = 0;
       facet < static_cast<int>(g_base_five_transversals.size()); ++facet) {
    for (int candidate = 0; candidate < static_cast<int>(g_candidates.size());
         ++candidate) {
      if (!WordBit(g_candidate_five_hits[candidate], facet)) {
        CandidateWordSet(&g_facet_killers[facet], candidate);
      }
    }
  }
  for (std::uint16_t six = 0; six <= kFull; ++six) {
    if (std::popcount(six) != 6) continue;
    std::vector<int> members;
    for (int facet = 0;
         facet < static_cast<int>(g_base_five_transversals.size()); ++facet) {
      if ((g_base_five_transversals[facet] & six) ==
          g_base_five_transversals[facet]) {
        members.push_back(facet);
      }
    }
    if (!members.empty()) g_six_set_members.push_back(std::move(members));
  }
  std::ranges::sort(g_six_set_members, [](const auto& left, const auto& right) {
    const auto reducer_count = [](const auto& members) {
      int count = 0;
      for (int candidate = 0;
           candidate < static_cast<int>(g_candidates.size()); ++candidate) {
        bool reduces = false;
        for (const int member : members) {
          if (!WordBit(g_candidate_five_hits[candidate], member)) {
            reduces = true;
            break;
          }
        }
        count += reduces;
      }
      return count;
    };
    return std::tuple(reducer_count(left), left) <
           std::tuple(reducer_count(right), right);
  });
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 2) {
      throw std::runtime_error(
          "usage: program a A4mask... [--symmetry-expand]");
    }
    g_a = std::stoi(argv[1]);
    g_symmetry_quotient =
        argc == g_a + 3 &&
        std::string(argv[argc - 1]) == "--symmetry-expand";
    if (g_a < 1 || g_a > 4 ||
        (argc != g_a + 2 && !g_symmetry_quotient)) {
      throw std::runtime_error("a must lie in 1..4 and match the A4 mask count");
    }
    for (int index = 0; index < g_a; ++index) {
      const auto value = static_cast<std::uint16_t>(
          std::stoul(argv[index + 2], nullptr, 0));
      if (std::popcount(value) != 4 || !Independent(value)) {
        throw std::runtime_error("invalid independent A4 mask");
      }
      g_fours.push_back(value);
    }
    if (std::set<std::uint16_t>(g_fours.begin(), g_fours.end()).size() !=
        g_fours.size()) {
      throw std::runtime_error("repeated A4 mask");
    }

    std::vector<std::uint16_t> fives;
    for (std::uint16_t value = 0; value <= kFull; ++value) {
      if (std::popcount(value) == 5 && Independent(value)) fives.push_back(value);
    }
    for (const auto five : fives) {
      if (std::ranges::any_of(g_fours, [five](auto four) {
            return (five & four) == four;
          })) {
        g_shadow.push_back(five);
      }
    }
    g_addition_count = 8 + 6 * g_a - static_cast<int>(g_shadow.size());
    g_scheduled_facets = 110 - 2 * g_a;
    if (g_addition_count < 0 || g_addition_count > 6) {
      throw std::runtime_error("addition count escaped supported range");
    }

    for (int vertex = 0; vertex < kN; ++vertex) {
      int a4_incidence = 0;
      for (const auto four : g_fours) a4_incidence += (four >> vertex) & 1;
      int forced = 0;
      for (const auto five : g_shadow) forced += (five >> vertex) & 1;
      const int previous = kBaseG3[vertex] - a4_incidence;
      const int before = kBaseG4[vertex] + 6 * a4_incidence - forced;
      g_lower[vertex] = std::max(0, before - MacaulayRaise(previous, 3));
      g_upper[vertex] = before - 1;
      if (g_lower[vertex] > g_upper[vertex]) {
        throw std::runtime_error("empty local addition interval");
      }
    }
    for (const auto five : fives) {
      if (std::ranges::find(g_shadow, five) != g_shadow.end()) continue;
      bool admitted = true;
      for (int vertex = 0; vertex < kN; ++vertex) {
        if (((five >> vertex) & 1) && g_upper[vertex] < 1) admitted = false;
      }
      if (admitted) g_candidates.push_back(five);
    }

    g_base.assign(kEdges.begin(), kEdges.end());
    g_base.insert(g_base.end(), g_fours.begin(), g_fours.end());
    g_base.insert(g_base.end(), g_shadow.begin(), g_shadow.end());
    for (std::uint16_t value = 1; value <= kFull; ++value) {
      if (std::popcount(value) <= 4 && HitsAll(value, g_base)) {
        g_small.push_back(value);
      }
    }
    const int cover_blocks = (g_small.size() + 63) / 64;
    g_coverers.assign(g_small.size(), {});
    for (int candidate = 0; candidate < static_cast<int>(g_candidates.size());
         ++candidate) {
      Word word(cover_blocks);
      for (int constraint = 0; constraint < static_cast<int>(g_small.size());
           ++constraint) {
        if (!(g_candidates[candidate] & g_small[constraint])) {
          word[constraint / 64] |= std::uint64_t{1} << (constraint % 64);
          g_coverers[constraint].push_back(candidate);
        }
      }
      g_coverage.push_back(std::move(word));
    }

    InitializeRidgeTables();
    if (g_symmetry_quotient) {
      ConfigureSymmetryQuotient();
      ValidateSymmetryQuotient();
    }
    Word initial((g_base_five_transversals.size() + 63) / 64,
                 ~std::uint64_t{0});
    if (!g_base_five_transversals.empty() &&
        g_base_five_transversals.size() % 64) {
      initial.back() =
          (std::uint64_t{1} << (g_base_five_transversals.size() % 64)) - 1;
    }
    g_ridge_seen.resize(g_addition_count + 1);
    g_symmetry_seen.resize(g_addition_count + 1);
    RidgeSearch({}, Word(cover_blocks), {}, std::move(initial));
    std::cerr << "ridge_search_complete nodes " << g_ridge_nodes
              << " canonical_targets " << g_ridge_solutions.size() << '\n';
    std::uint64_t expanded_literal_target_count = 0;
    if (g_symmetry_quotient) {
      ExpandCanonicalTargets();
      expanded_literal_target_count = AuditExpandedTargets();
    }

    std::cout << "{\n  \"kind\":\"" H9_14_RIDGE_KIND "\",\n"
              << "  \"a\":" << g_a << ",\n  \"representative_A4_masks\":[";
    for (int index = 0; index < g_a; ++index) {
      if (index) std::cout << ',';
      std::cout << g_fours[index];
    }
    std::cout << "],\n  \"forced_A5_shadow_count\":" << g_shadow.size()
              << ",\n  \"additional_A5_count\":" << g_addition_count
              << ",\n  \"scheduled_facet_count\":" << g_scheduled_facets
              << ",\n  \"addition_lower\":[";
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (vertex) std::cout << ',';
      std::cout << g_lower[vertex];
    }
    std::cout << "],\n  \"addition_upper\":[";
    for (int vertex = 0; vertex < kN; ++vertex) {
      if (vertex) std::cout << ',';
      std::cout << g_upper[vertex];
    }
    std::cout << ']'
              << ",\n  \"candidate_count\":" << g_candidates.size()
              << ",\n  \"base_small_transversal_count\":" << g_small.size()
              << ",\n  \"base_order_five_transversal_count\":"
              << g_base_five_transversals.size()
              << ",\n  \"nonempty_six_set_constraint_count\":"
              << g_six_set_members.size()
              << ",\n  \"search_node_count\":" << g_ridge_nodes
              << ",\n  \"collision_free_symmetry_quotient\":"
              << (g_symmetry_quotient ? "true" : "false")
              << ",\n  \"small_stabilizer_action_order\":"
              << g_fixed_actions.size()
              << ",\n  \"ridge_scheduled_candidate_count\":" << g_dead_exact_target
              << ",\n  \"canonical_literal_pure_target_count\":"
              << g_ridge_solutions.size()
              << ",\n  \"expanded_literal_pure_target_count\":"
              << expanded_literal_target_count
              << ",\n  \"canonical_literal_target_additions\":[";
    bool first_family = true;
    for (const auto& family : g_ridge_solutions) {
      if (!first_family) std::cout << ',';
      first_family = false;
      std::cout << '[';
      for (int row = 0; row < static_cast<int>(family.size()); ++row) {
        if (row) std::cout << ',';
        PrintMask(g_candidates[family[row]]);
      }
      std::cout << ']';
    }
    std::cout << "],\n  \"expanded_literal_target_additions\":[";
    first_family = true;
    for (const auto& family : g_expanded_ridge_solutions) {
      if (!first_family) std::cout << ',';
      first_family = false;
      std::cout << '[';
      for (int row = 0; row < static_cast<int>(family.size()); ++row) {
        if (row) std::cout << ',';
        PrintMask(family[row]);
      }
      std::cout << ']';
    }
    std::cout << "],\n  \"scope\":\"Exact necessary compact incidence/ridge gate; no Coxeter or finite-volume conclusion.\"\n}\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
