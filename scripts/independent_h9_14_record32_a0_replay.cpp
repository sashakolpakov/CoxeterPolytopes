// Independent lossless a=0 replay for H9/14 record 32.
//
// This source shares no project code with the producer.  It deliberately uses
// numeric active-mask order (rather than order-by-cardinality), an explicit
// column-multiplicity recursion, std::bitset facet words, and a literal
// minimal-cover/common-intersection transversal postcheck.  Its canonical
// target serialization is compared semantically with the production archive.

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <bitset>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

constexpr int kN = 14;
constexpr int kRows = 6;
constexpr int kIsolates = 8;
constexpr int kFacetWordSize = 262;
constexpr int kTargetFacets = 114;
constexpr std::array<int, 6> kActiveVertices = {0, 1, 2, 11, 12, 13};
constexpr std::array<int, 8> kIsolateVertices = {3, 4, 5, 6, 7, 8, 9, 10};
constexpr std::array<std::uint16_t, 4> kMissingEdges = {
    (1u << 0) | (1u << 11),
    (1u << 1) | (1u << 12),
    (1u << 0) | (1u << 13),
    (1u << 2) | (1u << 13),
};
constexpr std::array<int, kN> kCaps = {
    1, 2, 3, 5, 5, 5, 5, 5, 5, 5, 5, 3, 2, 1,
};

using Family = std::array<std::uint16_t, kRows>;
using FacetWord = std::bitset<kFacetWordSize>;

std::vector<std::uint16_t> g_active_masks;
std::vector<std::uint16_t> g_minimal_covers;
std::vector<std::uint16_t> g_five_covers;
std::array<FacetWord, 1 << kN> g_hit_words;
std::vector<Family> g_configurations;

struct Result {
  std::uint64_t paths = 0;
  std::uint64_t repeated = 0;
  std::uint64_t literal_small = 0;
  std::uint64_t blocker_valid = 0;
  std::uint64_t targets = 0;
  std::set<Family> target_families;
};

bool HitsEveryMissingEdge(std::uint16_t value) {
  for (const auto edge : kMissingEdges) {
    if ((value & edge) == 0) return false;
  }
  return true;
}

bool IsIndependent(std::uint16_t value) {
  for (const auto edge : kMissingEdges) {
    if ((value & edge) == edge) return false;
  }
  return true;
}

bool IsMinimalCover(std::uint16_t value) {
  if (!HitsEveryMissingEdge(value)) return false;
  for (int vertex = 0; vertex < kN; ++vertex) {
    const std::uint16_t bit = std::uint16_t{1} << vertex;
    if ((value & bit) && HitsEveryMissingEdge(value ^ bit)) return false;
  }
  return true;
}

void Initialize() {
  const std::uint16_t active_mask =
      (1u << 0) | (1u << 1) | (1u << 2) |
      (1u << 11) | (1u << 12) | (1u << 13);
  for (std::uint16_t value = 0; value < (1u << kN); ++value) {
    if ((value & ~active_mask) == 0 && IsIndependent(value)) {
      // Numeric order is intentional and differs from production.
      g_active_masks.push_back(value);
    }
    if (value && IsMinimalCover(value)) g_minimal_covers.push_back(value);
    if (std::popcount(value) == 5 && HitsEveryMissingEdge(value)) {
      g_five_covers.push_back(value);
    }
  }
  if (g_active_masks.size() != 24 || g_minimal_covers.size() != 6 ||
      g_five_covers.size() != kFacetWordSize) {
    throw std::runtime_error("record-32 independent state space changed");
  }
  for (std::uint16_t face = 0; face < (1u << kN); ++face) {
    for (std::size_t index = 0; index < g_five_covers.size(); ++index) {
      if (face & g_five_covers[index]) g_hit_words[face].set(index);
    }
  }
}

bool ActiveCapsHold(const Family& family) {
  for (const int vertex : kActiveVertices) {
    int incidence = 0;
    for (const auto row : family) incidence += (row >> vertex) & 1;
    if (incidence > kCaps[vertex]) return false;
  }
  return true;
}

std::array<std::uint8_t, 6> MissedRows(const Family& active_family) {
  std::array<std::uint8_t, 6> answer{};
  for (std::size_t cover = 0; cover < g_minimal_covers.size(); ++cover) {
    for (int row = 0; row < kRows; ++row) {
      if ((active_family[row] & g_minimal_covers[cover]) == 0) {
        answer[cover] |= std::uint8_t{1} << row;
      }
    }
  }
  return answer;
}

bool HasActiveSmallExtension(
    const Family& active_family,
    const std::array<std::uint8_t, 6>& missed) {
  for (const std::uint8_t rows : missed) {
    if (rows == 0) return true;
    for (const int vertex : kActiveVertices) {
      bool contained = true;
      for (int row = 0; row < kRows; ++row) {
        if ((rows & (1u << row)) &&
            !(active_family[row] & (std::uint16_t{1} << vertex))) {
          contained = false;
          break;
        }
      }
      if (contained) return true;
    }
  }
  return false;
}

bool HasLiteralSmallTransversal(const Family& family) {
  // Every minimal graph cover has order three.  For each such cover, compute
  // the literal intersection of all A5 rows it misses.  The cover itself is a
  // transversal when no row is missed; otherwise a common vertex extends it
  // to a transversal of order at most four.
  for (const auto cover : g_minimal_covers) {
    bool missed_any = false;
    std::uint16_t common = (1u << kN) - 1;
    for (const auto row : family) {
      if ((row & cover) == 0) {
        missed_any = true;
        common &= row;
      }
    }
    if (!missed_any || common != 0) return true;
  }
  return false;
}

int FacetCount(const Family& family) {
  FacetWord word;
  word.set();
  for (const auto row : family) word &= g_hit_words[row];
  return static_cast<int>(word.count());
}

void ConsumeFamily(const Family& active_family,
                   const std::array<std::uint8_t, kIsolates>& columns,
                   Result& result) {
  ++result.paths;
  Family family = active_family;
  for (int isolate = 0; isolate < kIsolates; ++isolate) {
    const int vertex = kIsolateVertices[isolate];
    for (int row = 0; row < kRows; ++row) {
      if (columns[isolate] & (1u << row)) {
        family[row] |= std::uint16_t{1} << vertex;
      }
    }
  }
  std::sort(family.begin(), family.end());
  if (std::adjacent_find(family.begin(), family.end()) != family.end()) {
    ++result.repeated;
    return;
  }
  if (HasLiteralSmallTransversal(family)) {
    ++result.literal_small;
    return;
  }
  ++result.blocker_valid;
  if (FacetCount(family) == kTargetFacets) {
    ++result.targets;
    result.target_families.insert(family);
  }
}

void EnumerateColumnMultiplicities(
    const Family& active_family,
    const std::vector<std::uint8_t>& allowed,
    int start,
    int columns_used,
    std::array<int, kRows>& row_sums,
    std::array<std::uint8_t, kIsolates>& columns,
    Result& result) {
  if (columns_used == kIsolates) {
    if (std::all_of(row_sums.begin(), row_sums.end(),
                    [](int value) { return value == 0; })) {
      ConsumeFamily(active_family, columns, result);
    }
    return;
  }
  const int columns_left = kIsolates - columns_used;
  for (int position = start; position < static_cast<int>(allowed.size());
       ++position) {
    const std::uint8_t column = allowed[position];
    int maximum = columns_left;
    for (int row = 0; row < kRows; ++row) {
      if (column & (1u << row)) maximum = std::min(maximum, row_sums[row]);
    }
    for (int multiplicity = 1; multiplicity <= maximum; ++multiplicity) {
      const int next_used = columns_used + multiplicity;
      bool feasible = true;
      for (int row = 0; row < kRows; ++row) {
        const int next = row_sums[row] -
                         multiplicity * ((column >> row) & 1);
        if (next < 0 || next > kIsolates - next_used) {
          feasible = false;
          break;
        }
      }
      if (!feasible) continue;
      for (int row = 0; row < kRows; ++row) {
        row_sums[row] -= multiplicity * ((column >> row) & 1);
      }
      for (int offset = 0; offset < multiplicity; ++offset) {
        columns[columns_used + offset] = column;
      }
      EnumerateColumnMultiplicities(
          active_family, allowed, position + 1, next_used,
          row_sums, columns, result);
      for (int row = 0; row < kRows; ++row) {
        row_sums[row] += multiplicity * ((column >> row) & 1);
      }
    }
  }
}

void ProcessConfiguration(const Family& active_family, Result& result) {
  const auto missed = MissedRows(active_family);
  std::vector<std::uint8_t> allowed;
  for (int column = 0; column < (1 << kRows); ++column) {
    if (std::popcount(static_cast<unsigned>(column)) > 5) continue;
    bool valid = true;
    for (const auto rows : missed) {
      if ((column & rows) == rows) {
        valid = false;
        break;
      }
    }
    if (valid) allowed.push_back(static_cast<std::uint8_t>(column));
  }
  std::array<int, kRows> row_sums{};
  for (int row = 0; row < kRows; ++row) {
    row_sums[row] = 5 - std::popcount(active_family[row]);
  }
  std::array<std::uint8_t, kIsolates> columns{};
  EnumerateColumnMultiplicities(
      active_family, allowed, 0, 0, row_sums, columns, result);
}

void EnumerateActiveFamilies(int depth, int start, Family& family,
                             std::uint64_t& cap_count,
                             std::uint64_t& cover_count) {
  if (depth == kRows) {
    if (!ActiveCapsHold(family)) return;
    ++cap_count;
    const auto missed = MissedRows(family);
    if (HasActiveSmallExtension(family, missed)) return;
    ++cover_count;
    g_configurations.push_back(family);
    return;
  }
  for (int index = start; index < static_cast<int>(g_active_masks.size()); ++index) {
    family[depth] = g_active_masks[index];
    EnumerateActiveFamilies(depth + 1, index, family, cap_count, cover_count);
  }
}

void WriteTargets(const std::string& path, const std::set<Family>& families) {
  std::ofstream output(path);
  if (!output) throw std::runtime_error("could not open independent target output");
  for (const auto& family : families) {
    for (int row = 0; row < kRows; ++row) {
      if (row) output << ',';
      output << family[row];
    }
    output << '\n';
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    std::string output_path;
    if (argc == 3 && std::string(argv[1]) == "--target-output") {
      output_path = argv[2];
    } else if (argc != 1) {
      std::cerr << "usage: " << argv[0] << " [--target-output PATH]\n";
      return 2;
    }
    Initialize();
    Family family{};
    std::uint64_t cap_count = 0;
    std::uint64_t cover_count = 0;
    EnumerateActiveFamilies(0, 0, family, cap_count, cover_count);
    if (cap_count != 22513 || cover_count != 14731) {
      throw std::runtime_error("independent active configuration count changed");
    }

    const unsigned available = std::thread::hardware_concurrency();
    const unsigned thread_count = std::max(1u, std::min(12u, available ? available : 1u));
    std::atomic<std::size_t> next{0};
    std::vector<Result> partial(thread_count);
    std::vector<std::thread> workers;
    for (unsigned worker = 0; worker < thread_count; ++worker) {
      workers.emplace_back([worker, &next, &partial]() {
        while (true) {
          const auto index = next.fetch_add(1, std::memory_order_relaxed);
          if (index >= g_configurations.size()) break;
          ProcessConfiguration(g_configurations[index], partial[worker]);
        }
      });
    }
    for (auto& worker : workers) worker.join();

    Result total;
    for (const auto& row : partial) {
      total.paths += row.paths;
      total.repeated += row.repeated;
      total.literal_small += row.literal_small;
      total.blocker_valid += row.blocker_valid;
      total.targets += row.targets;
      total.target_families.insert(
          row.target_families.begin(), row.target_families.end());
    }
    if (!output_path.empty()) WriteTargets(output_path, total.target_families);
    std::cout << "{\n"
              << "  \"kind\": \"independent_h9_14_record32_a0_replay\",\n"
              << "  \"active_mask_count\": " << g_active_masks.size() << ",\n"
              << "  \"active_after_caps\": " << cap_count << ",\n"
              << "  \"active_after_cover\": " << cover_count << ",\n"
              << "  \"column_paths\": " << total.paths << ",\n"
              << "  \"repeated_row_paths\": " << total.repeated << ",\n"
              << "  \"literal_small_transversal_paths\": "
              << total.literal_small << ",\n"
              << "  \"blocker_valid_paths\": " << total.blocker_valid << ",\n"
              << "  \"target_paths\": " << total.targets << ",\n"
              << "  \"target_family_count\": "
              << total.target_families.size() << "\n"
              << "}\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
