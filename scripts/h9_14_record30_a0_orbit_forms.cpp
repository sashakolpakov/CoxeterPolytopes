// Emit every row/isolate-double-lex matrix form in the full
// Aut(3K2) x S8(isolates) x S8(rows) orbit of one eight-A5 family.

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>

namespace {
constexpr std::array<int, 6> A = {0, 1, 2, 11, 12, 13};
constexpr std::array<int, 8> I = {3, 4, 5, 6, 7, 8, 9, 10};
using Family = std::array<std::uint16_t, 8>;

bool PreservesMatching(const std::array<int, 6>& p) {
  std::set<std::pair<int, int>> moved;
  for (auto [a, b] : std::array<std::pair<int, int>, 3>{{{0, 3}, {1, 4}, {2, 5}}}) {
    moved.insert(std::minmax(p[a], p[b]));
  }
  return moved == std::set<std::pair<int, int>>{{0, 3}, {1, 4}, {2, 5}};
}
}

int main(int argc, char** argv) {
  try {
    if (argc != 9) throw std::runtime_error("expected eight row masks");
    Family input{};
    for (int row = 0; row < 8; ++row) input[row] = std::stoi(argv[row + 1]);
    std::array<std::uint8_t, 8> isolate_columns{};
    for (int j = 0; j < 8; ++j) for (int row = 0; row < 8; ++row) {
      isolate_columns[j] |= ((input[row] >> I[j]) & 1u) << row;
    }
    std::set<Family> forms;
    std::array<int, 6> active_permutation = {0, 1, 2, 3, 4, 5};
    do {
      if (!PreservesMatching(active_permutation)) continue;
      std::array<std::uint16_t, 8> active_parts{};
      for (int row = 0; row < 8; ++row) for (int j = 0; j < 6; ++j) {
        if ((input[row] >> A[j]) & 1u) {
          active_parts[row] |= std::uint16_t{1} << A[active_permutation[j]];
        }
      }
      std::array<int, 8> order = {0, 1, 2, 3, 4, 5, 6, 7};
      do {
        std::array<std::uint8_t, 8> columns{};
        for (int j = 0; j < 8; ++j) for (int row = 0; row < 8; ++row) {
          columns[j] |= ((isolate_columns[j] >> order[row]) & 1u) << row;
        }
        std::sort(columns.begin(), columns.end());
        Family family{};
        for (int row = 0; row < 8; ++row) family[row] = active_parts[order[row]];
        for (int j = 0; j < 8; ++j) for (int row = 0; row < 8; ++row) {
          family[row] |= ((columns[j] >> row) & 1u) << I[j];
        }
        if (std::is_sorted(family.begin(), family.end()) &&
            std::adjacent_find(family.begin(), family.end()) == family.end()) {
          forms.insert(family);
        }
      } while (std::next_permutation(order.begin(), order.end()));
    } while (std::next_permutation(active_permutation.begin(), active_permutation.end()));
    std::cout << "{\"form_count\":" << forms.size() << ",\"forms\":[";
    bool first = true;
    for (const auto& family : forms) {
      if (!first) std::cout << ',';
      std::cout << '[';
      for (int row = 0; row < 8; ++row) {
        if (row) std::cout << ',';
        std::cout << family[row];
      }
      std::cout << ']';
      first = false;
    }
    std::cout << "]}\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
