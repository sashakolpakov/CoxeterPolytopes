// Lexicographically least unordered eight-A5 family under Aut(3K2) x S8.
// Also emit one original-vertex permutation attaining the least family.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>

namespace {
constexpr int N = 14;
constexpr int ROWS = 8;
constexpr std::array<int, 6> ACTIVE = {0, 1, 2, 11, 12, 13};
constexpr std::array<int, 8> ISOLATE = {3, 4, 5, 6, 7, 8, 9, 10};
using Mask = std::uint16_t;
using Family = std::array<Mask, ROWS>;

bool PreservesMatching(const std::array<int, 6>& permutation) {
  std::set<std::pair<int, int>> moved;
  for (auto [left, right] : std::array<std::pair<int, int>, 3>{{{0, 3}, {1, 4}, {2, 5}}}) {
    moved.insert(std::minmax(permutation[left], permutation[right]));
  }
  return moved == std::set<std::pair<int, int>>{{0, 3}, {1, 4}, {2, 5}};
}

std::uint8_t ReverseBits(std::uint8_t value) {
  std::uint8_t result = 0;
  for (int bit = 0; bit < ROWS; ++bit) result |= ((value >> bit) & 1u) << (7 - bit);
  return result;
}

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
    Family input{};
    for (int row = 0; row < ROWS; ++row) {
      const unsigned long value = std::stoul(argv[row + 1]);
      if (value >= (1u << N) || std::popcount(static_cast<unsigned>(value)) != 5) {
        throw std::runtime_error("each argument must be a 14-bit five-set");
      }
      input[row] = static_cast<Mask>(value);
    }
    std::sort(input.begin(), input.end());
    if (std::adjacent_find(input.begin(), input.end()) != input.end()) {
      throw std::runtime_error("A5 masks must be distinct");
    }
    std::array<std::uint8_t, 8> original_columns{};
    for (int column = 0; column < 8; ++column) for (int row = 0; row < ROWS; ++row) {
      original_columns[column] |= ((input[row] >> ISOLATE[column]) & 1u) << row;
    }

    Family best{};
    std::array<int, N> best_permutation{};
    bool first = true;
    std::array<int, 6> active_permutation = {0, 1, 2, 3, 4, 5};
    do {
      if (!PreservesMatching(active_permutation)) continue;
      std::array<Mask, ROWS> active_parts{};
      for (int row = 0; row < ROWS; ++row) for (int index = 0; index < 6; ++index) {
        if (input[row] & (Mask{1} << ACTIVE[index])) {
          active_parts[row] |= Mask{1} << ACTIVE[active_permutation[index]];
        }
      }
      std::array<int, ROWS> order = {0, 1, 2, 3, 4, 5, 6, 7};
      do {
        std::array<std::pair<std::uint8_t, int>, 8> columns{};
        for (int column = 0; column < 8; ++column) {
          std::uint8_t pattern = 0;
          for (int row = 0; row < ROWS; ++row) {
            pattern |= ((original_columns[column] >> order[row]) & 1u) << row;
          }
          columns[column] = {ReverseBits(static_cast<std::uint8_t>(~pattern)), column};
        }
        std::sort(columns.begin(), columns.end());
        Family candidate{};
        for (int row = 0; row < ROWS; ++row) candidate[row] = active_parts[order[row]];
        for (int target = 0; target < 8; ++target) {
          const std::uint8_t pattern = static_cast<std::uint8_t>(
              ~ReverseBits(columns[target].first));
          for (int row = 0; row < ROWS; ++row) {
            candidate[row] |= ((pattern >> row) & 1u) << ISOLATE[target];
          }
        }
        if (!std::is_sorted(candidate.begin(), candidate.end())) continue;
        if (!first && !(candidate < best)) continue;
        first = false;
        best = candidate;
        std::iota(best_permutation.begin(), best_permutation.end(), 0);
        for (int index = 0; index < 6; ++index) {
          best_permutation[ACTIVE[index]] = ACTIVE[active_permutation[index]];
        }
        for (int target = 0; target < 8; ++target) {
          best_permutation[ISOLATE[columns[target].second]] = ISOLATE[target];
        }
      } while (std::next_permutation(order.begin(), order.end()));
    } while (std::next_permutation(active_permutation.begin(), active_permutation.end()));
    if (first) throw std::runtime_error("no canonical form found");

    Family replay{};
    for (int row = 0; row < ROWS; ++row) replay[row] = Transform(input[row], best_permutation);
    std::sort(replay.begin(), replay.end());
    if (replay != best) throw std::runtime_error("canonical permutation replay failed");
    std::cout << "{\"canonical\":[";
    for (int row = 0; row < ROWS; ++row) {
      if (row) std::cout << ',';
      std::cout << best[row];
    }
    std::cout << "],\"permutation\":[";
    for (int vertex = 0; vertex < N; ++vertex) {
      if (vertex) std::cout << ',';
      std::cout << best_permutation[vertex];
    }
    std::cout << "]}\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
