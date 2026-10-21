// Hostile standalone replay of the corrected record-12 paired-branch census.
//
// This implementation imports no production census code.  It reads the
// independently verified 1,308 coloured branch representatives, enumerates
// both orientations of the source cut by an MRV coloured subgraph backtrack,
// reconstructs local blockers with an incremental algorithm, and applies the
// complete coupled blocker schedule.  Passing rows are written in the same
// semantic format as the exploratory production stream for set comparison.

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

constexpr int N = 10;
constexpr int GLOBAL_N = 14;
constexpr int FULL10 = (1 << N) - 1;

using LocalRow = uint16_t;
using GlobalRow = uint16_t;
using EdgeMask = uint64_t;
using PermCode = uint64_t;

int pair_index(int first, int second) {
  if (first > second) std::swap(first, second);
  int result = 0;
  for (int prior = 0; prior < first; ++prior) result += N - 1 - prior;
  return result + second - first - 1;
}

struct State {
  std::string word;
  int cut = 0;
  int y = 0;
  long long multiplicity = 0;
  std::array<LocalRow, N> adjacency{};
  EdgeMask edges = 0;
  std::vector<LocalRow> triples;
  std::vector<LocalRow> blocker;
  std::array<int, 4> colour_profile{};
};

std::array<LocalRow, N> decode_graph6(const std::string& word) {
  if (word.empty() || static_cast<unsigned char>(word[0]) - 63 != N) {
    throw std::runtime_error("hostile pair replay expected order-ten graph6");
  }
  std::vector<int> payload;
  for (std::size_t index = 1; index < word.size(); ++index) {
    const int code = static_cast<unsigned char>(word[index]) - 63;
    if (code < 0 || code > 63) throw std::runtime_error("bad graph6 code");
    for (int shift = 5; shift >= 0; --shift) payload.push_back((code >> shift) & 1);
  }
  std::array<LocalRow, N> result{};
  int at = 0;
  for (int second = 1; second < N; ++second) {
    for (int first = 0; first < second; ++first) {
      if (payload.at(at++)) {
        result[first] |= LocalRow{1} << second;
        result[second] |= LocalRow{1} << first;
      }
    }
  }
  for (; at < static_cast<int>(payload.size()); ++at) {
    if (payload[at]) throw std::runtime_error("nonzero graph6 padding");
  }
  return result;
}

bool edge(const std::array<LocalRow, N>& adjacency, int first, int second) {
  return (adjacency[first] >> second) & 1U;
}

bool cut_edge(int cut, int first, int second) {
  return ((cut >> first) & 1) != ((cut >> second) & 1);
}

EdgeMask edge_mask(const std::array<LocalRow, N>& adjacency) {
  EdgeMask result = 0;
  for (int first = 0; first < N; ++first) {
    for (int second = first + 1; second < N; ++second) {
      if (edge(adjacency, first, second)) {
        result |= EdgeMask{1} << pair_index(first, second);
      }
    }
  }
  return result;
}

std::vector<LocalRow> branch_triples(
    const std::array<LocalRow, N>& adjacency, int cut, int y) {
  std::vector<LocalRow> result;
  for (int vertex = 0; vertex < N; ++vertex) {
    if ((y >> vertex) & 1) {
      result.push_back((LocalRow{1} << 10) | (LocalRow{1} << 11) |
                       (LocalRow{1} << vertex));
    }
  }
  for (int first = 0; first < N; ++first) {
    for (int second = first + 1; second < N; ++second) {
      const bool x = edge(adjacency, first, second);
      const bool delta = cut_edge(cut, first, second);
      if (x) {
        result.push_back((LocalRow{1} << 11) | (LocalRow{1} << first) |
                         (LocalRow{1} << second));
      }
      if (x != delta) {
        result.push_back((LocalRow{1} << 10) | (LocalRow{1} << first) |
                         (LocalRow{1} << second));
      }
    }
  }
  for (int first = 0; first < N; ++first) {
    for (int second = first + 1; second < N; ++second) {
      for (int third = second + 1; third < N; ++third) {
        const int parity = edge(adjacency, first, second) +
                           edge(adjacency, first, third) +
                           edge(adjacency, second, third);
        if (parity & 1) {
          result.push_back((LocalRow{1} << first) | (LocalRow{1} << second) |
                           (LocalRow{1} << third));
        }
      }
    }
  }
  std::sort(result.begin(), result.end());
  result.erase(std::unique(result.begin(), result.end()), result.end());
  return result;
}

template <typename Row>
std::vector<Row> minimalize(std::vector<Row> candidates) {
  std::sort(candidates.begin(), candidates.end(), [](Row first, Row second) {
    const int first_size = __builtin_popcountll(static_cast<uint64_t>(first));
    const int second_size = __builtin_popcountll(static_cast<uint64_t>(second));
    return first_size != second_size ? first_size < second_size : first < second;
  });
  candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
  std::vector<Row> result;
  for (Row row : candidates) {
    bool contains_prior = false;
    for (Row prior : result) {
      if ((prior & row) == prior) {
        contains_prior = true;
        break;
      }
    }
    if (!contains_prior) result.push_back(row);
  }
  return result;
}

std::vector<LocalRow> blocker(const std::vector<LocalRow>& clutter) {
  // In particular, the blocker of the empty clutter remains {empty}.
  std::vector<LocalRow> rows{0};
  for (LocalRow hyperedge : clutter) {
    std::vector<LocalRow> candidates;
    for (LocalRow row : rows) {
      if (row & hyperedge) {
        candidates.push_back(row);
      } else {
        for (int bits = hyperedge; bits; bits &= bits - 1) {
          const int bit = bits & -bits;
          candidates.push_back(row | static_cast<LocalRow>(bit));
        }
      }
    }
    rows = minimalize(std::move(candidates));
  }
  return rows;
}

std::array<int, 4> profile(int cut, int y) {
  std::array<int, 4> result{};
  for (int vertex = 0; vertex < N; ++vertex) {
    ++result[2 * ((cut >> vertex) & 1) + ((y >> vertex) & 1)];
  }
  return result;
}

int allowed_z_values(int x, int delta_p, int delta_q, int yi, int yj) {
  int result = 0;
  for (int z = 0; z < 2; ++z) {
    const int total = yi + yj + x + (x ^ delta_p) + z + (z ^ delta_q);
    if (total == 0 || total == 2) result |= 1 << z;
  }
  return result;
}

std::array<LocalRow, N> allowed_graph(const State& state, long long* no_forced_one) {
  std::array<LocalRow, N> result{};
  const int q = state.cut ^ state.y;
  for (int first = 0; first < N; ++first) {
    for (int second = first + 1; second < N; ++second) {
      const int allowed = allowed_z_values(
          edge(state.adjacency, first, second),
          cut_edge(state.cut, first, second), cut_edge(q, first, second),
          (state.y >> first) & 1, (state.y >> second) & 1);
      if (allowed == 0) throw std::runtime_error("locally infeasible retained state");
      if (allowed == 2) throw std::runtime_error("a forced-one z case exists");
      ++*no_forced_one;
      if (allowed & 2) {
        result[first] |= LocalRow{1} << second;
        result[second] |= LocalRow{1} << first;
      }
    }
  }
  return result;
}

struct EmbeddingEnumerator {
  const State& source;
  const std::array<LocalRow, N>& allowed;
  int source_cut;
  int target_cut;
  int target_y;
  std::array<int, N> image{};
  int used = 0;
  std::unordered_map<EdgeMask, PermCode> images;

  EmbeddingEnumerator(const State& source_value,
                      const std::array<LocalRow, N>& allowed_value,
                      int source_cut_value, int target_cut_value, int target_y_value)
      : source(source_value), allowed(allowed_value), source_cut(source_cut_value),
        target_cut(target_cut_value), target_y(target_y_value) {
    image.fill(-1);
  }

  std::vector<int> candidates(int vertex) const {
    std::vector<int> result;
    const int wanted_cut = (source_cut >> vertex) & 1;
    const int wanted_y = (source.y >> vertex) & 1;
    for (int target = 0; target < N; ++target) {
      if ((used >> target) & 1) continue;
      if (wanted_cut != ((target_cut >> target) & 1) ||
          wanted_y != ((target_y >> target) & 1)) continue;
      bool compatible = true;
      for (int prior = 0; prior < N; ++prior) {
        if (image[prior] < 0) continue;
        if (edge(source.adjacency, vertex, prior) &&
            !edge(allowed, target, image[prior])) {
          compatible = false;
          break;
        }
      }
      if (compatible) result.push_back(target);
    }
    return result;
  }

  void recurse(int mapped_count) {
    if (mapped_count == N) {
      EdgeMask transformed = 0;
      PermCode permutation = 0;
      for (int vertex = 0; vertex < N; ++vertex) {
        permutation |= PermCode(image[vertex]) << (4 * vertex);
      }
      for (int first = 0; first < N; ++first) {
        for (int second = first + 1; second < N; ++second) {
          if (edge(source.adjacency, first, second)) {
            transformed |= EdgeMask{1} << pair_index(image[first], image[second]);
          }
        }
      }
      images.emplace(transformed, permutation);
      return;
    }

    int chosen = -1;
    std::vector<int> chosen_images;
    for (int vertex = 0; vertex < N; ++vertex) {
      if (image[vertex] >= 0) continue;
      auto options = candidates(vertex);
      if (options.empty()) return;
      if (chosen < 0 || options.size() < chosen_images.size() ||
          (options.size() == chosen_images.size() &&
           __builtin_popcount(source.adjacency[vertex]) >
               __builtin_popcount(source.adjacency[chosen]))) {
        chosen = vertex;
        chosen_images = std::move(options);
      }
    }
    for (int target : chosen_images) {
      image[chosen] = target;
      used |= 1 << target;
      recurse(mapped_count + 1);
      used ^= 1 << target;
      image[chosen] = -1;
    }
  }
};

LocalRow transform_local(LocalRow row, PermCode permutation) {
  LocalRow result = row & ((LocalRow{1} << 10) | (LocalRow{1} << 11));
  for (int vertex = 0; vertex < N; ++vertex) {
    if ((row >> vertex) & 1U) {
      const int target = (permutation >> (4 * vertex)) & 15;
      result |= LocalRow{1} << target;
    }
  }
  return result;
}

std::vector<GlobalRow> coupled_blocker(
    const State& p_state, const State& q_state, PermCode permutation) {
  std::vector<GlobalRow> p_rows{GlobalRow{1} << 0, GlobalRow{1} << 1};
  for (LocalRow row : p_state.blocker) {
    GlobalRow global = (row & ((LocalRow{1} << N) - 1)) << 2;
    if ((row >> 10) & 1U) global |= GlobalRow{1} << 12;
    if ((row >> 11) & 1U) global |= GlobalRow{1} << 13;
    p_rows.push_back(global);
  }
  std::vector<GlobalRow> q_rows{GlobalRow{1} << 12, GlobalRow{1} << 13};
  for (LocalRow source_row : q_state.blocker) {
    const LocalRow row = transform_local(source_row, permutation);
    GlobalRow global = (row & ((LocalRow{1} << N) - 1)) << 2;
    if ((row >> 10) & 1U) global |= GlobalRow{1} << 0;
    if ((row >> 11) & 1U) global |= GlobalRow{1} << 1;
    q_rows.push_back(global);
  }
  std::vector<GlobalRow> candidates;
  for (GlobalRow first : p_rows) {
    for (GlobalRow second : q_rows) candidates.push_back(first | second);
  }
  return minimalize(std::move(candidates));
}

bool schedule_passes(const State& p_state, const State& q_state,
                     PermCode permutation, int facet_count, int* a_out, int* b_out) {
  const auto rows = coupled_blocker(p_state, q_state, permutation);
  std::array<int, GLOBAL_N + 1> histogram{};
  std::vector<GlobalRow> twos;
  for (GlobalRow row : rows) {
    const int size = __builtin_popcount(static_cast<unsigned>(row));
    ++histogram[size];
    if (size == 2) twos.push_back(row);
  }
  std::sort(twos.begin(), twos.end());
  const std::vector<GlobalRow> expected{
      GlobalRow((1U << 0) | (1U << 12)), GlobalRow((1U << 1) | (1U << 12)),
      GlobalRow((1U << 0) | (1U << 13)), GlobalRow((1U << 1) | (1U << 13))};
  if (twos != expected) return false;
  for (int size = 0; size <= GLOBAL_N; ++size) {
    if (histogram[size] && size != 2 && size != 4 && size != 5) return false;
  }
  const int a = histogram[4];
  if (a > 7 || facet_count != 130 - 2 * a) return false;
  *a_out = a;
  *b_out = histogram[5];
  return true;
}

}  // namespace

int main() {
  std::vector<State> states;
  std::string line;
  while (std::getline(std::cin, line)) {
    if (line.empty()) continue;
    std::istringstream input(line);
    State state;
    input >> state.word >> state.cut >> state.y >> state.multiplicity;
    if (!input) throw std::runtime_error("malformed hostile state row");
    state.adjacency = decode_graph6(state.word);
    state.edges = edge_mask(state.adjacency);
    state.triples = branch_triples(state.adjacency, state.cut, state.y);
    state.blocker = blocker(state.triples);
    state.colour_profile = profile(state.cut, state.y);
    states.push_back(std::move(state));
  }
  if (states.size() != 1308) throw std::runtime_error("expected 1308 states");

  long long type_pairs = 0;
  long long direct_profile_pairs = 0;
  long long complement_profile_pairs = 0;
  long long both_profile_pairs = 0;
  long long direct_images = 0;
  long long complement_images = 0;
  long long image_overlap = 0;
  long long embedded_graphs = 0;
  long long passing = 0;
  long long no_forced_one_checks = 0;

  for (int p_index = 0; p_index < static_cast<int>(states.size()); ++p_index) {
    const State& p_state = states[p_index];
    const auto allowed = allowed_graph(p_state, &no_forced_one_checks);
    const int target_cut = p_state.cut ^ p_state.y;
    const auto direct_target_profile = profile(target_cut, p_state.y);
    const auto complement_target_profile = profile(target_cut ^ FULL10, p_state.y);
    for (int q_index = 0; q_index < static_cast<int>(states.size()); ++q_index) {
      const State& q_state = states[q_index];
      const int facet_count = static_cast<int>(p_state.triples.size() +
          q_state.triples.size()) - __builtin_popcount(static_cast<unsigned>(p_state.y));
      if (facet_count < 116 || facet_count > 130 || (facet_count & 1)) continue;
      const bool direct = q_state.colour_profile == direct_target_profile;
      const bool complemented = q_state.colour_profile == complement_target_profile;
      if (!direct && !complemented) continue;
      ++type_pairs;
      direct_profile_pairs += direct;
      complement_profile_pairs += complemented;
      both_profile_pairs += direct && complemented;

      std::unordered_map<EdgeMask, PermCode> images;
      std::unordered_map<EdgeMask, PermCode> direct_map;
      std::unordered_map<EdgeMask, PermCode> complement_map;
      if (direct) {
        EmbeddingEnumerator search(q_state, allowed, q_state.cut, target_cut, p_state.y);
        search.recurse(0);
        direct_map = std::move(search.images);
        direct_images += direct_map.size();
        for (const auto& item : direct_map) images.emplace(item);
      }
      if (complemented) {
        EmbeddingEnumerator search(
            q_state, allowed, q_state.cut ^ FULL10, target_cut, p_state.y);
        search.recurse(0);
        complement_map = std::move(search.images);
        complement_images += complement_map.size();
        for (const auto& item : complement_map) {
          if (images.contains(item.first)) ++image_overlap;
          images.emplace(item);
        }
      }
      embedded_graphs += images.size();
      for (const auto& [z_mask, permutation] : images) {
        int a = -1;
        int b = -1;
        if (!schedule_passes(p_state, q_state, permutation, facet_count, &a, &b)) continue;
        ++passing;
        std::cout << "PAIR " << p_index << ' ' << q_index << ' '
                  << p_state.edges << ' ' << z_mask << ' ' << p_state.cut << ' '
                  << p_state.y << ' ' << facet_count << ' ' << a << ' ' << b << '\n';
      }
    }
    if ((p_index + 1) % 100 == 0) {
      std::cerr << "HOSTILE_PROGRESS " << (p_index + 1) << '/' << states.size()
                << " types=" << type_pairs << " embeddings=" << embedded_graphs
                << " passing=" << passing << '\n';
    }
  }
  std::cerr << "HOSTILE_SUMMARY states=" << states.size()
            << " type_pairs=" << type_pairs
            << " direct_profile_pairs=" << direct_profile_pairs
            << " complement_profile_pairs=" << complement_profile_pairs
            << " both_profile_pairs=" << both_profile_pairs
            << " direct_images=" << direct_images
            << " complement_images=" << complement_images
            << " image_overlap=" << image_overlap
            << " embeddings=" << embedded_graphs
            << " passing=" << passing
            << " no_forced_one_checks=" << no_forced_one_checks << '\n';
}
