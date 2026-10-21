// Exact exploratory coupling of the sound record-12 branch prefilter.
//
// Input is the coloured-orbit stream produced by
// explore_h9_14_record12_unfiltered_orbits.py after retaining precisely those
// branch states whose minimal transversals meeting a branch special have
// orders four or five.  This program enumerates every compatible relative
// embedding of the second branch, applies the full coupled blocker, and
// prints the embeddings satisfying the global missing-face sizes/counts.
// No one-sided condition is imposed on isolate-only branch transversals.

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

constexpr int N = 10;
constexpr int LOCAL_N = 12;
constexpr int GLOBAL_N = 14;

struct State {
  std::string word;
  int p = 0;
  int y = 0;
  long long multiplicity = 0;
  std::array<uint16_t, N> adjacency{};
  uint64_t edge_mask = 0;
  int triple_count = 0;
  std::array<int, 4> profile{};  // 00,01,10,11 for (p,y).
  std::vector<uint16_t> local_blocker;
};

int pair_index(int first, int second) {
  if (first > second) std::swap(first, second);
  int result = 0;
  for (int i = 0; i < first; ++i) result += 9 - i;
  result += second - first - 1;
  return result;
}

std::array<uint16_t, N> decode_graph6(const std::string& word) {
  if (word.empty() || static_cast<unsigned char>(word[0]) - 63 != N) {
    throw std::runtime_error("expected an order-ten short graph6 word");
  }
  std::vector<int> bits;
  for (std::size_t k = 1; k < word.size(); ++k) {
    const int value = static_cast<unsigned char>(word[k]) - 63;
    for (int shift = 5; shift >= 0; --shift) bits.push_back((value >> shift) & 1);
  }
  std::array<uint16_t, N> adjacency{};
  int at = 0;
  for (int second = 1; second < N; ++second) {
    for (int first = 0; first < second; ++first) {
      if (bits.at(at++)) {
        adjacency[first] |= uint16_t{1} << second;
        adjacency[second] |= uint16_t{1} << first;
      }
    }
  }
  return adjacency;
}

uint64_t edge_mask(const std::array<uint16_t, N>& adjacency) {
  uint64_t result = 0;
  for (int first = 0; first < N; ++first) {
    for (int second = first + 1; second < N; ++second) {
      if ((adjacency[first] >> second) & 1U) {
        result |= uint64_t{1} << pair_index(first, second);
      }
    }
  }
  return result;
}

std::array<uint16_t, N> adjacency_from_edge_mask(uint64_t value) {
  std::array<uint16_t, N> result{};
  for (int first = 0; first < N; ++first) {
    for (int second = first + 1; second < N; ++second) {
      if ((value >> pair_index(first, second)) & 1U) {
        result[first] |= uint16_t{1} << second;
        result[second] |= uint16_t{1} << first;
      }
    }
  }
  return result;
}

bool edge(const std::array<uint16_t, N>& adjacency, int first, int second) {
  return (adjacency[first] >> second) & 1U;
}

bool cut_edge(int cut, int first, int second) {
  return ((cut >> first) & 1) != ((cut >> second) & 1);
}

std::vector<uint16_t> triples(const std::array<uint16_t, N>& adjacency,
                              int p, int y) {
  std::vector<uint16_t> result;
  for (int i = 0; i < N; ++i) {
    if ((y >> i) & 1) result.push_back((1U << 10) | (1U << 11) | (1U << i));
    for (int j = i + 1; j < N; ++j) {
      if (edge(adjacency, i, j)) {
        result.push_back((1U << 11) | (1U << i) | (1U << j));
      }
      if (edge(adjacency, i, j) != cut_edge(p, i, j)) {
        result.push_back((1U << 10) | (1U << i) | (1U << j));
      }
    }
  }
  for (int i = 0; i < N; ++i) {
    for (int j = i + 1; j < N; ++j) {
      for (int k = j + 1; k < N; ++k) {
        const int parity = edge(adjacency, i, j) + edge(adjacency, i, k) +
                           edge(adjacency, j, k);
        if (parity & 1) result.push_back((1U << i) | (1U << j) | (1U << k));
      }
    }
  }
  std::sort(result.begin(), result.end());
  result.erase(std::unique(result.begin(), result.end()), result.end());
  return result;
}

std::vector<uint16_t> blocker_masks(const std::vector<uint16_t>& clutter) {
  std::array<bool, 1 << LOCAL_N> hits{};
  for (int value = 0; value < (1 << LOCAL_N); ++value) {
    hits[value] = std::all_of(clutter.begin(), clutter.end(),
                              [value](uint16_t row) { return value & row; });
  }
  if (hits[0]) return {0};
  std::vector<uint16_t> blocker;
  for (int value = 1; value < (1 << LOCAL_N); ++value) {
    if (!hits[value]) continue;
    bool minimal = true;
    for (int bits = value; bits; bits &= bits - 1) {
      const int bit = bits & -bits;
      if (hits[value ^ bit]) {
        minimal = false;
        break;
      }
    }
    if (minimal) blocker.push_back(value);
  }
  return blocker;
}

int allowed_z_mask(int x, int delta_p, int delta_q, int yi, int yj) {
  int result = 0;
  for (int z = 0; z < 2; ++z) {
    const int total = yi + yj + x + (x ^ delta_p) + z + (z ^ delta_q);
    if (total == 0 || total == 2) result |= 1 << z;
  }
  return result;
}

std::array<uint16_t, N> allowed_adjacency(const State& state) {
  std::array<uint16_t, N> result{};
  const int q = state.p ^ state.y;
  for (int first = 0; first < N; ++first) {
    for (int second = first + 1; second < N; ++second) {
      const int allowed = allowed_z_mask(
          edge(state.adjacency, first, second),
          cut_edge(state.p, first, second), cut_edge(q, first, second),
          (state.y >> first) & 1, (state.y >> second) & 1);
      if (allowed & 2) {
        result[first] |= uint16_t{1} << second;
        result[second] |= uint16_t{1} << first;
      }
    }
  }
  return result;
}

std::array<int, 4> profile(int p, int y) {
  std::array<int, 4> result{};
  for (int vertex = 0; vertex < N; ++vertex) {
    ++result[2 * ((p >> vertex) & 1) + ((y >> vertex) & 1)];
  }
  return result;
}

std::array<int, 4> target_q_profile(const State& state, bool complement) {
  return profile((state.p ^ state.y) ^ (complement ? (1 << N) - 1 : 0),
                 state.y);
}

uint16_t transform_local_row(uint16_t row, uint64_t permutation_code) {
  uint16_t result = row & ((1U << 10) | (1U << 11));
  for (int vertex = 0; vertex < N; ++vertex) {
    if ((row >> vertex) & 1U) {
      const int image = (permutation_code >> (4 * vertex)) & 15;
      result |= uint16_t{1} << image;
    }
  }
  return result;
}

std::vector<uint16_t> full_blocker(const State& p_state, const State& q_state,
                                   uint64_t permutation_code) {
  std::vector<uint16_t> p_blocker{uint16_t{1} << 0, uint16_t{1} << 1};
  for (uint16_t row : p_state.local_blocker) {
    uint16_t global = (row & ((1U << N) - 1)) << 2;
    if (row & (1U << 10)) global |= 1U << 12;
    if (row & (1U << 11)) global |= 1U << 13;
    p_blocker.push_back(global);
  }
  std::vector<uint16_t> q_blocker{uint16_t{1} << 12, uint16_t{1} << 13};
  for (uint16_t source_row : q_state.local_blocker) {
    const uint16_t row = transform_local_row(source_row, permutation_code);
    uint16_t global = (row & ((1U << N) - 1)) << 2;
    if (row & (1U << 10)) global |= 1U << 0;
    if (row & (1U << 11)) global |= 1U << 1;
    q_blocker.push_back(global);
  }

  std::vector<uint16_t> candidates;
  candidates.reserve(p_blocker.size() * q_blocker.size());
  for (uint16_t first : p_blocker) {
    for (uint16_t second : q_blocker) candidates.push_back(first | second);
  }
  std::sort(candidates.begin(), candidates.end(), [](uint16_t first, uint16_t second) {
    const int first_size = __builtin_popcount(static_cast<unsigned>(first));
    const int second_size = __builtin_popcount(static_cast<unsigned>(second));
    return first_size != second_size ? first_size < second_size : first < second;
  });
  candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
  std::vector<uint16_t> result;
  for (uint16_t row : candidates) {
    if (std::none_of(result.begin(), result.end(),
                     [row](uint16_t prior) { return (prior & row) == prior; })) {
      result.push_back(row);
    }
  }
  return result;
}

bool blocker_schedule_passes(const State& p_state, const State& q_state,
                             uint64_t permutation_code,
                             int facet_count, int* a_out, int* b_out) {
  const auto blocker = full_blocker(p_state, q_state, permutation_code);
  std::array<int, GLOBAL_N + 1> histogram{};
  std::vector<uint16_t> twos;
  for (uint16_t row : blocker) {
    ++histogram[__builtin_popcount(static_cast<unsigned>(row))];
    if (__builtin_popcount(static_cast<unsigned>(row)) == 2) twos.push_back(row);
  }
  std::sort(twos.begin(), twos.end());
  const std::vector<uint16_t> expected_twos{
      uint16_t((1U << 0) | (1U << 12)), uint16_t((1U << 1) | (1U << 12)),
      uint16_t((1U << 0) | (1U << 13)), uint16_t((1U << 1) | (1U << 13))};
  if (twos != expected_twos) return false;
  for (int size = 0; size <= GLOBAL_N; ++size) {
    if (histogram[size] && size != 2 && size != 4 && size != 5) return false;
  }
  const int a = histogram[4];
  const int b = histogram[5];
  // Do not identify b with the archive's ``missing_five_count_offset``: the
  // latter counts all missing independent five-sets, whereas b counts only
  // inclusion-minimal five-nonfaces.  The complete face-vector audit is a
  // downstream check.  At this stage only the exact facet/a schedule is used.
  if (a < 0 || a > 7 || facet_count != 130 - 2 * a) {
    return false;
  }
  *a_out = a;
  *b_out = b;
  return true;
}

struct EmbeddingSearch {
  const State& source;
  const std::array<uint16_t, N>& allowed;
  int source_cut;
  int target_cut;
  int target_y;
  std::array<int, N> order{};
  std::array<int, N> image{};
  int used = 0;
  std::unordered_map<uint64_t, uint64_t> images;

  EmbeddingSearch(const State& source_value,
                  const std::array<uint16_t, N>& allowed_value,
                  int source_cut_value, int target_cut_value, int target_y_value)
      : source(source_value), allowed(allowed_value), source_cut(source_cut_value),
        target_cut(target_cut_value), target_y(target_y_value) {
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [this](int first, int second) {
      const int first_degree = __builtin_popcount(source.adjacency[first]);
      const int second_degree = __builtin_popcount(source.adjacency[second]);
      return first_degree != second_degree ? first_degree > second_degree : first < second;
    });
    image.fill(-1);
  }

  void recurse(int index) {
    if (index == N) {
      uint64_t result = 0;
      for (int first = 0; first < N; ++first) {
        for (int second = first + 1; second < N; ++second) {
          if (!edge(source.adjacency, first, second)) continue;
          result |= uint64_t{1} << pair_index(image[first], image[second]);
        }
      }
      uint64_t permutation_code = 0;
      for (int vertex = 0; vertex < N; ++vertex) {
        permutation_code |= uint64_t(image[vertex]) << (4 * vertex);
      }
      images.emplace(result, permutation_code);
      return;
    }
    const int vertex = order[index];
    const int source_p = (source_cut >> vertex) & 1;
    const int source_y = (source.y >> vertex) & 1;
    for (int target = 0; target < N; ++target) {
      if ((used >> target) & 1) continue;
      if (source_p != ((target_cut >> target) & 1) ||
          source_y != ((target_y >> target) & 1)) continue;
      bool possible = true;
      for (int prior_index = 0; prior_index < index; ++prior_index) {
        const int prior = order[prior_index];
        if (edge(source.adjacency, vertex, prior) &&
            !edge(allowed, target, image[prior])) {
          possible = false;
          break;
        }
      }
      if (!possible) continue;
      image[vertex] = target;
      used |= 1 << target;
      recurse(index + 1);
      used ^= 1 << target;
      image[vertex] = -1;
    }
  }
};

}  // namespace

int main() {
  std::vector<State> states;
  std::string line;
  while (std::getline(std::cin, line)) {
    if (line.empty()) continue;
    std::istringstream input(line);
    State state;
    input >> state.word >> state.p >> state.y >> state.multiplicity;
    state.adjacency = decode_graph6(state.word);
    state.edge_mask = edge_mask(state.adjacency);
    state.triple_count = static_cast<int>(triples(state.adjacency, state.p, state.y).size());
    state.profile = profile(state.p, state.y);
    state.local_blocker = blocker_masks(triples(state.adjacency, state.p, state.y));
    states.push_back(state);
  }

  long long type_pairs = 0;
  long long embedded_graphs = 0;
  long long blocker_passing = 0;
  for (int p_index = 0; p_index < static_cast<int>(states.size()); ++p_index) {
    const State& p_state = states[p_index];
    const auto allowed = allowed_adjacency(p_state);
    const int target_q = p_state.p ^ p_state.y;
    for (int q_index = 0; q_index < static_cast<int>(states.size()); ++q_index) {
      const State& q_state = states[q_index];
      const int facet_count = p_state.triple_count + q_state.triple_count -
                              __builtin_popcount(static_cast<unsigned>(p_state.y));
      if (facet_count < 116 || facet_count > 130 || facet_count % 2) continue;
      std::unordered_map<uint64_t, uint64_t> pair_images;
      bool profile_compatible = false;
      for (int complemented = 0; complemented < 2; ++complemented) {
        if (q_state.profile != target_q_profile(p_state, complemented)) continue;
        profile_compatible = true;
        const int source_cut = q_state.p ^ (complemented ? (1 << N) - 1 : 0);
        EmbeddingSearch search(q_state, allowed, source_cut, target_q, p_state.y);
        search.recurse(0);
        for (const auto& [image_mask, permutation_code] : search.images) {
          pair_images.emplace(image_mask, permutation_code);
        }
      }
      if (!profile_compatible) continue;
      ++type_pairs;
      embedded_graphs += pair_images.size();
      for (const auto& [z_mask, permutation_code] : pair_images) {
        int a = -1;
        int b = -1;
        if (!blocker_schedule_passes(
                p_state, q_state, permutation_code, facet_count, &a, &b)) continue;
        ++blocker_passing;
        std::cout << "PAIR " << p_index << ' ' << q_index << ' '
                  << p_state.edge_mask << ' ' << z_mask << ' ' << p_state.p << ' '
                  << p_state.y << ' ' << facet_count << ' ' << a << ' ' << b << '\n';
      }
    }
    if ((p_index + 1) % 100 == 0) {
      std::cerr << "progress=" << (p_index + 1) << '/' << states.size()
                << " type_pairs=" << type_pairs << " embeddings=" << embedded_graphs
                << " passing=" << blocker_passing << '\n';
    }
  }
  std::cerr << "states=" << states.size() << " type_pairs=" << type_pairs
            << " embeddings=" << embedded_graphs << " passing=" << blocker_passing
            << '\n';
}
