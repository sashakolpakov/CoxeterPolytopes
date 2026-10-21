// Canonicalize an unordered family of eight record-30 A5 masks under
// Aut(3K2) x S8.  The colored incidence graph has six active vertices,
// eight isolates, three matching-edge nodes, and eight A5-row nodes.

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <nauty.h>

namespace {
constexpr int N = 14;
constexpr int ROWS = 8;
using Mask = std::uint16_t;

Mask Transform(Mask value, const std::array<int, N>& permutation) {
  Mask result = 0;
  for (int vertex = 0; vertex < N; ++vertex) {
    if (value & (Mask{1} << vertex)) result |= Mask{1} << permutation[vertex];
  }
  return result;
}
}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != ROWS + 1) throw std::runtime_error("expected eight decimal A5 masks");
    std::array<Mask, ROWS> family{};
    for (int row = 0; row < ROWS; ++row) {
      const unsigned long value = std::stoul(argv[row + 1]);
      if (value >= (1u << N) || std::popcount(static_cast<unsigned>(value)) != 5) {
        throw std::runtime_error("each argument must be a 14-bit five-set");
      }
      family[row] = static_cast<Mask>(value);
    }
    std::sort(family.begin(), family.end());
    if (std::adjacent_find(family.begin(), family.end()) != family.end()) {
      throw std::runtime_error("A5 masks must be distinct");
    }

    constexpr int node_count = N + 3 + ROWS;
    const int words = SETWORDSNEEDED(node_count);
    std::vector<graph> input(static_cast<std::size_t>(words) * node_count);
    std::vector<graph> canonical(static_cast<std::size_t>(words) * node_count);
    EMPTYGRAPH(input.data(), words, node_count);
    for (int edge = 0; edge < 3; ++edge) {
      ADDONEEDGE(input.data(), N + edge, edge, words);
      ADDONEEDGE(input.data(), N + edge, 11 + edge, words);
    }
    for (int row = 0; row < ROWS; ++row) {
      for (int vertex = 0; vertex < N; ++vertex) if (family[row] & (Mask{1} << vertex)) {
        ADDONEEDGE(input.data(), N + 3 + row, vertex, words);
      }
    }
    std::array<int, node_count> lab{}, ptn{}, orbits{};
    int position = 0;
    auto add_cell = [&](const std::vector<int>& cell) {
      for (int value : cell) {
        lab[position] = value;
        ptn[position++] = 1;
      }
      ptn[position - 1] = 0;
    };
    add_cell({0, 1, 2, 11, 12, 13});
    add_cell({3, 4, 5, 6, 7, 8, 9, 10});
    add_cell({14, 15, 16});
    add_cell({17, 18, 19, 20, 21, 22, 23, 24});
    DEFAULTOPTIONS_GRAPH(options);
    statsblk stats{};
    options.defaultptn = FALSE;
    options.getcanon = TRUE;
    densenauty(input.data(), lab.data(), ptn.data(), orbits.data(), &options,
               &stats, words, node_count, canonical.data());
    if (stats.errstatus != 0) throw std::runtime_error("nauty canonicalization failed");

    // lab[canonical_position] is the old vertex at that position.  Convert
    // the two original-vertex color cells back to the record-30 numbering.
    constexpr std::array<int, 6> active_labels = {0, 1, 2, 11, 12, 13};
    constexpr std::array<int, 8> isolate_labels = {3, 4, 5, 6, 7, 8, 9, 10};
    std::array<int, N> permutation{};
    for (int canonical_position = 0; canonical_position < 6; ++canonical_position) {
      permutation[lab[canonical_position]] = active_labels[canonical_position];
    }
    for (int canonical_position = 6; canonical_position < 14; ++canonical_position) {
      permutation[lab[canonical_position]] = isolate_labels[canonical_position - 6];
    }
    std::array<Mask, ROWS> representative{};
    for (int row = 0; row < ROWS; ++row) representative[row] = Transform(family[row], permutation);
    std::sort(representative.begin(), representative.end());

    std::cout << "{\"canonical\":[";
    for (int row = 0; row < ROWS; ++row) {
      if (row) std::cout << ',';
      std::cout << representative[row];
    }
    std::cout << "],\"permutation\":[";
    for (int vertex = 0; vertex < N; ++vertex) {
      if (vertex) std::cout << ',';
      std::cout << permutation[vertex];
    }
    std::cout << "]}\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
