// Exploratory exhaustive branch filter for the record-12 graph shadow.
//
// Feed `geng -q -t 10` on stdin.  For every unlabelled triangle-free graph X
// and every normalized cut p with X^p triangle-free, construct the twelve-
// vertex branch triple system.  The optional command-line arguments fix the
// cut size and y mask.  The special y value -2 prints every state satisfying
// the coupled ridge equations.  The special y value -3 applies only the safe
// part of the branch-blocker filter: blocker rows meeting either branch
// special must have order four or five.  Rows supported entirely on the ten
// isolates are deliberately retained because the opposite branch can make
// them non-transversals of the coupled clutter.

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr int N = 10;
constexpr int BRANCH_N = 12;

std::vector<std::pair<int, int>> decode_graph6(const std::string& word) {
  if (word.empty() || static_cast<unsigned char>(word[0]) - 63 != N) {
    throw std::runtime_error("expected an order-ten short graph6 word");
  }
  std::vector<int> bits;
  for (std::size_t k = 1; k < word.size(); ++k) {
    const int value = static_cast<unsigned char>(word[k]) - 63;
    for (int shift = 5; shift >= 0; --shift) bits.push_back((value >> shift) & 1);
  }
  std::vector<std::pair<int, int>> edges;
  int at = 0;
  // graph6 scans the strict lower triangle by columns:
  // (0,1), (0,2), (1,2), (0,3), ...
  for (int second = 1; second < N; ++second) {
    for (int first = 0; first < second; ++first) {
      if (bits.at(at++)) edges.emplace_back(first, second);
    }
  }
  return edges;
}

bool edge(const std::array<uint16_t, N>& adjacency, int i, int j) {
  return (adjacency[i] >> j) & 1U;
}

bool cut_edge(int p, int i, int j) {
  return ((p >> i) & 1) != ((p >> j) & 1);
}

bool triangle_free(const std::array<uint16_t, N>& adjacency, int p) {
  for (int i = 0; i < N; ++i) {
    for (int j = i + 1; j < N; ++j) {
      for (int k = j + 1; k < N; ++k) {
        if ((edge(adjacency, i, j) != cut_edge(p, i, j)) &&
            (edge(adjacency, i, k) != cut_edge(p, i, k)) &&
            (edge(adjacency, j, k) != cut_edge(p, j, k))) {
          return false;
        }
      }
    }
  }
  return true;
}

bool locally_pair_feasible(const std::array<uint16_t, N>& adjacency,
                           int p, int y) {
  const int q = p ^ y;
  for (int first = 0; first < N; ++first) {
    for (int second = first + 1; second < N; ++second) {
      const int x = edge(adjacency, first, second);
      const int dp = cut_edge(p, first, second);
      const int dq = cut_edge(q, first, second);
      const int yi = (y >> first) & 1;
      const int yj = (y >> second) & 1;
      bool possible = false;
      for (int z = 0; z < 2; ++z) {
        const int total = yi + yj + x + (x ^ dp) + z + (z ^ dq);
        possible |= total == 0 || total == 2;
      }
      if (!possible) return false;
    }
  }
  return true;
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

std::vector<uint16_t> blocker_masks(const std::vector<uint16_t>& triples) {
  std::array<bool, 1 << BRANCH_N> hits{};
  for (int value = 0; value < (1 << BRANCH_N); ++value) {
    hits[value] = std::all_of(triples.begin(), triples.end(),
                              [value](uint16_t triple) { return value & triple; });
  }
  std::vector<uint16_t> blocker;
  // The blocker of the empty clutter is {empty}.  Retaining that row is
  // essential in the safe-filter modes: after y-edges are added it produces
  // the two singleton specials and the y support.
  if (hits[0]) return {0};
  for (int value = 1; value < (1 << BRANCH_N); ++value) {
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

std::vector<uint16_t> add_y_edges(const std::vector<uint16_t>& base_blocker,
                                  int y) {
  std::vector<uint16_t> candidates;
  const uint16_t specials = (1U << 10) | (1U << 11);
  for (uint16_t row : base_blocker) {
    if (row & specials) {
      candidates.push_back(row);
    } else {
      candidates.push_back(static_cast<uint16_t>(row | y));
      candidates.push_back(static_cast<uint16_t>(row | (1U << 10)));
      candidates.push_back(static_cast<uint16_t>(row | (1U << 11)));
    }
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

std::array<int, BRANCH_N + 1> blocker_histogram(
    const std::vector<uint16_t>& blocker) {
  std::array<int, BRANCH_N + 1> histogram{};
  for (uint16_t row : blocker) {
    ++histogram[__builtin_popcount(static_cast<unsigned>(row))];
  }
  return histogram;
}

std::string histogram_string(const std::array<int, BRANCH_N + 1>& histogram) {
  std::ostringstream out;
  bool first = true;
  for (int size = 0; size <= BRANCH_N; ++size) {
    if (!histogram[size]) continue;
    if (!first) out << ',';
    first = false;
    out << size << ':' << histogram[size];
  }
  return out.str();
}

}  // namespace

int main(int argc, char** argv) {
  const int wanted_cut_size = argc >= 2 ? std::atoi(argv[1]) : -1;
  const int wanted_y = argc >= 3 ? std::atoi(argv[2]) : 0;
  std::string word;
  long long graph_count = 0;
  long long cut_count = 0;
  long long survivor_count = 0;
  while (std::getline(std::cin, word)) {
    if (word.empty()) continue;
    ++graph_count;
    std::array<uint16_t, N> adjacency{};
    for (const auto [i, j] : decode_graph6(word)) {
      adjacency[i] |= 1U << j;
      adjacency[j] |= 1U << i;
    }
    if (!triangle_free(adjacency, 0)) throw std::runtime_error("non-triangle-free X");
    // Normalize complementary cuts by insisting 0 is outside p.
    for (int p = 0; p < (1 << N); ++p) {
      if (p & 1) continue;
      if (wanted_cut_size >= 0 && __builtin_popcount(static_cast<unsigned>(p)) != wanted_cut_size) continue;
      if (!triangle_free(adjacency, p)) continue;
      ++cut_count;
      if (wanted_y == -2) {
        for (int y = 0; y < (1 << N); ++y) {
          if (!locally_pair_feasible(adjacency, p, y)) continue;
          ++survivor_count;
          std::cout << word << ' ' << p << ' ' << y << '\n';
        }
        continue;
      }
      const auto base_blocker = blocker_masks(triples(adjacency, p, 0));
      if (wanted_y == -3) {
        const uint16_t specials = (1U << 10) | (1U << 11);
        for (int y = 0; y < (1 << N); ++y) {
          if (!locally_pair_feasible(adjacency, p, y)) continue;
          const auto blocker = add_y_edges(base_blocker, y);
          bool allowed = true;
          for (uint16_t row : blocker) {
            const int size = __builtin_popcount(static_cast<unsigned>(row));
            if ((row & specials) && size != 4 && size != 5) {
              allowed = false;
              break;
            }
          }
          if (!allowed) continue;
          ++survivor_count;
          std::cout << word << ' ' << p << ' ' << y << '\n';
        }
        continue;
      }
      const int y_first = wanted_y < 0 ? 0 : wanted_y;
      const int y_last = wanted_y < 0 ? (1 << N) : wanted_y + 1;
      for (int y = y_first; y < y_last; ++y) {
        const auto histogram = blocker_histogram(add_y_edges(base_blocker, y));
        bool allowed = histogram[4] || histogram[5];
        for (int size = 0; size <= BRANCH_N; ++size) {
          if (histogram[size] && size != 4 && size != 5) allowed = false;
        }
        if (!allowed) continue;
        ++survivor_count;
        std::cout << word << ' ' << p << ' ' << y << ' '
                  << histogram_string(histogram) << '\n';
      }
    }
  }
  std::cerr << "graphs=" << graph_count << " valid_cuts=" << cut_count
            << " survivors=" << survivor_count << '\n';
}
