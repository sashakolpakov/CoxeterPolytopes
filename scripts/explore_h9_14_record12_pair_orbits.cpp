// Exploratory full graph-shadow quotient for balanced record-12 branch pairs.

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

namespace {

using Matrix = uint32_t;  // 5x5 matrix, row-major.

constexpr Matrix K = (Matrix{1} << 25) - 1;

int pair_index(int first, int second) {
  int result = 0;
  for (int i = 0; i < first; ++i) result += 9 - i;
  result += second - first - 1;
  return result;
}

Matrix from_pair_mask(uint64_t value) {
  Matrix result = 0;
  for (int first = 0; first < 5; ++first) {
    for (int second = 5; second < 10; ++second) {
      if (value & (uint64_t{1} << pair_index(first, second))) {
        result |= Matrix{1} << (5 * first + second - 5);
      }
    }
  }
  return result;
}

Matrix transpose(Matrix value) {
  Matrix result = 0;
  for (int row = 0; row < 5; ++row) {
    for (int column = 0; column < 5; ++column) {
      if (value & (Matrix{1} << (5 * row + column))) {
        result |= Matrix{1} << (5 * column + row);
      }
    }
  }
  return result;
}

uint64_t layer_min(Matrix x, Matrix z) {
  uint64_t result = ~uint64_t{0};
  for (int complement_x = 0; complement_x < 2; ++complement_x) {
    const Matrix xx = complement_x ? (x ^ K) : x;
    for (int complement_z = 0; complement_z < 2; ++complement_z) {
      const Matrix zz = complement_z ? (z ^ K) : z;
      result = std::min(result, uint64_t{xx} | (uint64_t{zz} << 25));
      result = std::min(result, uint64_t{zz} | (uint64_t{xx} << 25));
    }
  }
  return result;
}

std::vector<std::array<int, 5>> permutations5() {
  std::array<int, 5> value{0, 1, 2, 3, 4};
  std::vector<std::array<int, 5>> result;
  do {
    result.push_back(value);
  } while (std::next_permutation(value.begin(), value.end()));
  return result;
}

uint64_t canonical(Matrix x, Matrix z) {
  const auto permutations = permutations5();
  std::array<std::array<int, 32>, 120> row_maps{};
  for (int index = 0; index < 120; ++index) {
    for (int row = 0; row < 32; ++row) {
      int image = 0;
      for (int bit = 0; bit < 5; ++bit) {
        if (row & (1 << bit)) image |= 1 << permutations[index][bit];
      }
      row_maps[index][row] = image;
    }
  }
  uint64_t result = ~uint64_t{0};
  for (int transposed = 0; transposed < 2; ++transposed) {
    Matrix source_x = transposed ? transpose(x) : x;
    Matrix source_z = transposed ? transpose(z) : z;
    for (const auto& row_permutation : permutations) {
      for (int column_index = 0; column_index < 120; ++column_index) {
        Matrix image_x = 0;
        Matrix image_z = 0;
        for (int source_row = 0; source_row < 5; ++source_row) {
          const int target_row = row_permutation[source_row];
          const int x_row = (source_x >> (5 * source_row)) & 31;
          const int z_row = (source_z >> (5 * source_row)) & 31;
          image_x |= Matrix(row_maps[column_index][x_row]) << (5 * target_row);
          image_z |= Matrix(row_maps[column_index][z_row]) << (5 * target_row);
        }
        result = std::min(result, layer_min(image_x, image_z));
      }
    }
  }
  return result;
}

}  // namespace

int main() {
  std::string line;
  std::map<uint64_t, int> orbit_multiplicity;
  int input_count = 0;
  while (std::getline(std::cin, line)) {
    if (!line.starts_with("BALANCED_PAIR ")) continue;
    std::istringstream input(line);
    std::string marker;
    uint64_t x_mask = 0;
    uint64_t z_mask = 0;
    input >> marker >> x_mask >> z_mask;
    const Matrix x = from_pair_mask(x_mask);
    const Matrix z = from_pair_mask(z_mask);
    if (__builtin_popcountll(x_mask) != __builtin_popcount(x) ||
        __builtin_popcountll(z_mask) != __builtin_popcount(z)) {
      throw std::runtime_error("a balanced graph has a non-cross edge");
    }
    ++orbit_multiplicity[canonical(x, z)];
    ++input_count;
  }
  std::cerr << "input=" << input_count << " full_orbits="
            << orbit_multiplicity.size() << '\n';
  for (const auto& [representative, multiplicity] : orbit_multiplicity) {
    const Matrix x = representative & K;
    const Matrix z = (representative >> 25) & K;
    std::cout << "ORBIT " << x << ' ' << z << ' ' << multiplicity << '\n';
  }
}
