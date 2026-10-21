// Independent exact quotient replay for the corrected record-12 census.
//
// This deliberately shares no record-12 production code.  It canonicalizes
// the two 5x5 binary layers under independent layer complementation, layer
// swap, S5 x S5 row/column relabelling, and bipartition transpose.  The
// column permutations are eliminated exactly by sorting the five combined
// two-layer column words for each of the 120 row permutations.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

using u64 = std::uint64_t;

constexpr u64 kMatrixAll = (u64{1} << 25) - 1;
constexpr u64 kBalancedCut = ((u64{1} << 5) - 1) << 5;

[[noreturn]] void fail(const std::string& message) {
  throw std::runtime_error(message);
}

std::vector<std::array<int, 5>> permutations5() {
  std::vector<std::array<int, 5>> result;
  std::array<int, 5> permutation{0, 1, 2, 3, 4};
  do {
    result.push_back(permutation);
  } while (std::next_permutation(permutation.begin(), permutation.end()));
  if (result.size() != 120) fail("S5 enumeration failed");
  return result;
}

const std::vector<std::array<int, 5>> kPermutations5 = permutations5();

u64 transpose_matrix(u64 value) {
  u64 result = 0;
  for (int row = 0; row < 5; ++row) {
    for (int column = 0; column < 5; ++column) {
      if (value & (u64{1} << (5 * row + column))) {
        result |= u64{1} << (5 * column + row);
      }
    }
  }
  return result;
}

u64 canonical_bipartite_pair(u64 x, u64 z) {
  std::array<unsigned, 5> x_rows{};
  std::array<unsigned, 5> z_rows{};
  for (int row = 0; row < 5; ++row) {
    x_rows[row] = static_cast<unsigned>((x >> (5 * row)) & 31);
    z_rows[row] = static_cast<unsigned>((z >> (5 * row)) & 31);
  }

  u64 best = std::numeric_limits<u64>::max();
  for (const auto& row_order : kPermutations5) {
    std::array<u64, 5> columns{};
    for (int source_column = 0; source_column < 5; ++source_column) {
      u64 value = 0;
      for (int target_row = 0; target_row < 5; ++target_row) {
        const int source_row = row_order[target_row];
        value |= u64((x_rows[source_row] >> source_column) & 1)
                 << (5 * target_row);
        value |= u64((z_rows[source_row] >> source_column) & 1)
                 << (25 + 5 * target_row);
      }
      columns[source_column] = value;
    }
    std::sort(columns.begin(), columns.end(), std::greater<u64>());
    u64 encoded = 0;
    for (int target_column = 0; target_column < 5; ++target_column) {
      encoded |= columns[target_column] << target_column;
    }
    best = std::min(best, encoded);
  }
  return best;
}

u64 canonical_pair(u64 x, u64 z) {
  u64 best = std::numeric_limits<u64>::max();
  for (int swapped = 0; swapped < 2; ++swapped) {
    const u64 first = swapped ? z : x;
    const u64 second = swapped ? x : z;
    for (int complement_first = 0; complement_first < 2;
         ++complement_first) {
      for (int complement_second = 0; complement_second < 2;
           ++complement_second) {
        const u64 xx = first ^ (complement_first ? kMatrixAll : 0);
        const u64 zz = second ^ (complement_second ? kMatrixAll : 0);
        best = std::min(best, canonical_bipartite_pair(xx, zz));
        best = std::min(
            best,
            canonical_bipartite_pair(transpose_matrix(xx),
                                       transpose_matrix(zz)));
      }
    }
  }
  return best;
}

int pair_index(int first, int second) {
  if (first > second) std::swap(first, second);
  int index = 0;
  for (int left = 0; left < 10; ++left) {
    for (int right = left + 1; right < 10; ++right, ++index) {
      if (left == first && right == second) return index;
    }
  }
  fail("invalid pair");
}

u64 matrix_from_pair_mask(u64 value) {
  u64 result = 0;
  for (int row = 0; row < 5; ++row) {
    for (int column = 0; column < 5; ++column) {
      if (value & (u64{1} << pair_index(row, column + 5))) {
        result |= u64{1} << (5 * row + column);
      }
    }
  }
  if (std::popcount(result) != std::popcount(value)) {
    fail("balanced graph has a non-cross edge");
  }
  return result;
}

std::vector<std::string> read_nonempty_lines(const std::string& path) {
  std::ifstream input(path);
  if (!input) fail("cannot read " + path);
  std::vector<std::string> result;
  std::string line;
  while (std::getline(input, line)) {
    if (!line.empty()) result.push_back(line);
  }
  return result;
}

struct PairRow {
  u64 x;
  u64 z;
  int a;
};

PairRow parse_pair(const std::string& line) {
  std::istringstream stream(line);
  std::string marker;
  std::uint64_t p_index, q_index, raw_x, raw_z, cut, y, facets, a, b;
  if (!(stream >> marker >> p_index >> q_index >> raw_x >> raw_z >> cut >> y >>
        facets >> a >> b) || marker != "PAIR") {
    fail("malformed PAIR row");
  }
  std::string excess;
  if (stream >> excess) fail("excess field in PAIR row");
  if (cut != kBalancedCut || y != 0 || facets != 130 - 2 * a ||
      b != 12 - 2 * a || a > 4) {
    fail("passing pair left balanced Dehn--Sommerville schedule");
  }
  return {matrix_from_pair_mask(raw_x), matrix_from_pair_mask(raw_z),
          static_cast<int>(a)};
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 4) {
      std::cerr << "usage: " << argv[0]
                << " PRODUCTION_PAIRS HOSTILE_PAIRS ORBIT_STREAM\n";
      return 2;
    }
    const auto production_lines = read_nonempty_lines(argv[1]);
    const auto hostile_lines = read_nonempty_lines(argv[2]);
    if (production_lines.size() != 42000 || hostile_lines.size() != 42000) {
      fail("pair stream does not contain 42,000 rows");
    }
    const std::unordered_set<std::string> production_set(
        production_lines.begin(), production_lines.end());
    const std::unordered_set<std::string> hostile_set(hostile_lines.begin(),
                                                       hostile_lines.end());
    if (production_set.size() != 42000 || hostile_set.size() != 42000 ||
        production_set != hostile_set) {
      fail("production and hostile PAIR sets differ");
    }

    std::unordered_map<u64, std::uint64_t> observed;
    std::unordered_map<u64, int> branch;
    for (const auto& line : hostile_lines) {
      const PairRow row = parse_pair(line);
      const u64 key = canonical_pair(row.x, row.z);
      ++observed[key];
      const auto [position, inserted] = branch.emplace(key, row.a);
      if (!inserted && position->second != row.a) {
        fail("quotient mixed Dehn--Sommerville branches");
      }
    }

    std::unordered_map<u64, std::uint64_t> archived;
    for (const auto& line : read_nonempty_lines(argv[3])) {
      std::istringstream stream(line);
      std::string marker, excess;
      u64 x, z, multiplicity;
      if (!(stream >> marker >> x >> z >> multiplicity) || marker != "ORBIT" ||
          (stream >> excess)) {
        fail("malformed ORBIT row");
      }
      const u64 key = x | (z << 25);
      if (canonical_pair(x, z) != key) {
        fail("archived ORBIT representative is not canonical");
      }
      if (!archived.emplace(key, multiplicity).second) {
        fail("archived ORBIT representative repeats");
      }
    }
    if (archived != observed || observed.size() != 371) {
      fail("replayed and archived full quotients differ");
    }

    std::map<int, int> histogram;
    std::uint64_t multiplicity_sum = 0;
    for (const auto& [key, multiplicity] : observed) {
      ++histogram[branch.at(key)];
      multiplicity_sum += multiplicity;
    }
    const std::map<int, int> expected{{0, 48}, {1, 95}, {2, 144},
                                      {3, 63}, {4, 21}};
    if (histogram != expected || multiplicity_sum != 42000) {
      fail("371-orbit histogram or multiplicity sum changed");
    }

    std::cout << "{\n"
              << "  \"passing_pair_rows\": 42000,\n"
              << "  \"full_symmetry_orbits\": " << observed.size() << ",\n"
              << "  \"orbit_multiplicity_sum\": " << multiplicity_sum << ",\n"
              << "  \"schedule_orbit_histogram\": {\"0\": 48, \"1\": 95, "
                 "\"2\": 144, \"3\": 63, \"4\": 21},\n"
              << "  \"production_hostile_pair_sets_equal\": true,\n"
              << "  \"archived_orbit_multiset_equal\": true,\n"
              << "  \"all_archived_representatives_canonical\": true\n"
              << "}\n";
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
