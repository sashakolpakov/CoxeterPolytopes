// Lossless standalone a=0 completion scout for H9/14 record 32.
//
// This executable intentionally imports no project code.  It enumerates all
// six unordered A5 rows using every independent active subset of P4 + K2 and
// all multiplicity vectors of the eight S8-equivalent isolate columns.  The
// output is a search-space census, not a topology or Coxeter theorem.

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

namespace {

constexpr int kVertexCount = 14;
constexpr int kRowCount = 6;
constexpr int kIsolateCount = 8;
constexpr std::array<int, 6> kActive = {0, 1, 2, 11, 12, 13};
constexpr std::array<int, 8> kIsolates = {3, 4, 5, 6, 7, 8, 9, 10};
constexpr std::array<std::uint16_t, 4> kEdges = {
    (1u << 0) | (1u << 11),
    (1u << 1) | (1u << 12),
    (1u << 0) | (1u << 13),
    (1u << 2) | (1u << 13),
};
constexpr std::array<int, kVertexCount> kA5Caps = {
    1, 2, 3, 5, 5, 5, 5, 5, 5, 5, 5, 3, 2, 1,
};
constexpr int kScheduledFacetCount = 114;
constexpr int kExpectedFiveCoverCount = 262;
constexpr int kWordCount = (kExpectedFiveCoverCount + 63) / 64;

using Family = std::array<std::uint16_t, kRowCount>;
using Word = std::array<std::uint64_t, kWordCount>;

std::vector<std::uint16_t> g_active_rows;
std::vector<std::uint16_t> g_graph_covers;
std::vector<std::uint16_t> g_five_covers;
std::array<Word, 1 << kVertexCount> g_hit_words{};

std::uint64_t g_cap_feasible_active_configurations = 0;
std::uint64_t g_active_cover_feasible_configurations = 0;
std::uint64_t g_distribution_paths = 0;
std::uint64_t g_repeated_row_paths = 0;
std::uint64_t g_small_transversal_paths = 0;
std::uint64_t g_blocker_valid_paths = 0;
std::uint64_t g_target_paths = 0;
std::array<std::uint64_t, kExpectedFiveCoverCount + 1> g_facet_histogram{};
std::set<Family> g_target_families;
std::map<std::pair<std::uint64_t, std::uint32_t>, std::uint64_t>
    g_distribution_count_cache;
std::uint64_t g_positive_distribution_problems = 0;
bool g_count_only = false;
bool g_literal_small_transversal_postcheck = false;
std::vector<Family> g_stream_configurations;

struct StreamResult {
  std::uint64_t distribution_paths = 0;
  std::uint64_t repeated_row_paths = 0;
  std::uint64_t small_transversal_paths = 0;
  std::uint64_t blocker_valid_paths = 0;
  std::uint64_t target_paths = 0;
  std::array<std::uint64_t, kExpectedFiveCoverCount + 1> facet_histogram{};
  std::set<Family> target_families;
};

thread_local StreamResult* g_thread_result = nullptr;

bool IsIndependent(std::uint16_t value) {
  for (std::uint16_t edge : kEdges) {
    if ((value & edge) == edge) return false;
  }
  return true;
}

bool CoversGraph(std::uint16_t value) {
  for (std::uint16_t edge : kEdges) {
    if ((value & edge) == 0) return false;
  }
  return true;
}

bool IsMinimalGraphCover(std::uint16_t value) {
  if (!CoversGraph(value)) return false;
  for (int vertex = 0; vertex < kVertexCount; ++vertex) {
    const std::uint16_t bit = std::uint16_t{1} << vertex;
    if ((value & bit) && CoversGraph(value ^ bit)) return false;
  }
  return true;
}

void InitializeRowsAndCovers() {
  const std::uint16_t active_mask =
      (1u << 0) | (1u << 1) | (1u << 2) |
      (1u << 11) | (1u << 12) | (1u << 13);
  for (std::uint16_t value = 0; value < (1u << kVertexCount); ++value) {
    if ((value & ~active_mask) == 0 && IsIndependent(value)) {
      g_active_rows.push_back(value);
    }
    if (value && IsMinimalGraphCover(value)) g_graph_covers.push_back(value);
    if (std::popcount(value) == 5 && CoversGraph(value)) {
      g_five_covers.push_back(value);
    }
  }
  std::sort(g_active_rows.begin(), g_active_rows.end(), [](auto left, auto right) {
    if (std::popcount(left) != std::popcount(right)) {
      return std::popcount(left) < std::popcount(right);
    }
    return left < right;
  });
  const std::array<std::uint16_t, 6> expected_covers = {
      (1u << 0) | (1u << 1) | (1u << 2),
      (1u << 0) | (1u << 2) | (1u << 12),
      (1u << 0) | (1u << 1) | (1u << 13),
      (1u << 1) | (1u << 11) | (1u << 13),
      (1u << 0) | (1u << 12) | (1u << 13),
      (1u << 11) | (1u << 12) | (1u << 13),
  };
  if (g_active_rows.size() != 24 ||
      !std::equal(g_graph_covers.begin(), g_graph_covers.end(),
                  expected_covers.begin(), expected_covers.end())) {
    throw std::runtime_error("record-32 active masks or graph covers changed");
  }
}

void InitializeHitWords() {
  if (g_five_covers.size() != kExpectedFiveCoverCount) {
    throw std::runtime_error("record-32 five-cover count changed");
  }
  for (std::uint16_t face = 0; face < (1u << kVertexCount); ++face) {
    Word word{};
    for (std::size_t index = 0; index < g_five_covers.size(); ++index) {
      if (face & g_five_covers[index]) {
        word[index / 64] |= std::uint64_t{1} << (index % 64);
      }
    }
    g_hit_words[face] = word;
  }
}

bool HasRepeatedRows(const Family& family) {
  for (int index = 1; index < kRowCount; ++index) {
    if (family[index - 1] == family[index]) return true;
  }
  return false;
}

bool HasSmallTransversal(const Family& family) {
  // Every inclusion-minimal graph cover has order three.  It extends inside
  // the order-four budget exactly when the rows it misses have nonempty
  // common intersection (or when it misses no row).
  for (std::uint16_t cover : g_graph_covers) {
    std::uint16_t common = (1u << kVertexCount) - 1;
    bool missed = false;
    for (std::uint16_t face : family) {
      if ((face & cover) == 0) {
        missed = true;
        common &= face;
      }
    }
    if (!missed || common) return true;
  }
  return false;
}

int FacetCount(const Family& family) {
  Word word = g_hit_words[family[0]];
  for (int row = 1; row < kRowCount; ++row) {
    for (int block = 0; block < kWordCount; ++block) {
      word[block] &= g_hit_words[family[row]][block];
    }
  }
  int answer = 0;
  for (std::uint64_t block : word) answer += std::popcount(block);
  return answer;
}

void ConsumeColumns(const Family& active_family,
                    const std::array<std::uint8_t, kIsolateCount>& columns) {
  if (g_thread_result == nullptr) {
    throw std::runtime_error("stream worker lacks a result sink");
  }
  StreamResult& result = *g_thread_result;
  ++result.distribution_paths;
  Family family = active_family;
  for (int isolate_index = 0; isolate_index < kIsolateCount; ++isolate_index) {
    const int vertex = kIsolates[isolate_index];
    const std::uint8_t column = columns[isolate_index];
    for (int row = 0; row < kRowCount; ++row) {
      if (column & (1u << row)) family[row] |= std::uint16_t{1} << vertex;
    }
  }
  std::sort(family.begin(), family.end());
  if (HasRepeatedRows(family)) {
    ++result.repeated_row_paths;
    return;
  }
  if (g_literal_small_transversal_postcheck && HasSmallTransversal(family)) {
    ++result.small_transversal_paths;
    return;
  }
  ++result.blocker_valid_paths;
  const int facet_count = FacetCount(family);
  ++result.facet_histogram[facet_count];
  if (facet_count == kScheduledFacetCount) {
    ++result.target_paths;
    result.target_families.insert(family);
  }
}

void EnumerateColumns(const Family& active_family,
                      const std::vector<std::uint8_t>& allowed,
                      int start,
                      int depth,
                      std::array<int, kRowCount>& remaining,
                      std::array<std::uint8_t, kIsolateCount>& columns) {
  if (depth == kIsolateCount) {
    if (std::all_of(remaining.begin(), remaining.end(),
                    [](int value) { return value == 0; })) {
      ConsumeColumns(active_family, columns);
    }
    return;
  }
  const int columns_left_after = kIsolateCount - depth - 1;
  for (int index = start; index < static_cast<int>(allowed.size()); ++index) {
    const std::uint8_t column = allowed[index];
    bool feasible = true;
    for (int row = 0; row < kRowCount; ++row) {
      const int next = remaining[row] - ((column >> row) & 1);
      if (next < 0 || next > columns_left_after) {
        feasible = false;
        break;
      }
    }
    if (!feasible) continue;
    for (int row = 0; row < kRowCount; ++row) {
      remaining[row] -= (column >> row) & 1;
    }
    columns[depth] = column;
    EnumerateColumns(active_family, allowed, index, depth + 1, remaining, columns);
    for (int row = 0; row < kRowCount; ++row) {
      remaining[row] += (column >> row) & 1;
    }
  }
}

std::uint64_t CountColumns(const std::vector<std::uint8_t>& allowed,
                           int start,
                           int depth,
                           std::array<int, kRowCount>& remaining,
                           std::unordered_map<std::uint64_t, std::uint64_t>& memo) {
  if (depth == kIsolateCount) {
    return std::all_of(remaining.begin(), remaining.end(),
                       [](int value) { return value == 0; });
  }
  std::uint64_t key = static_cast<std::uint64_t>(start);
  key |= static_cast<std::uint64_t>(depth) << 6;
  for (int row = 0; row < kRowCount; ++row) {
    key |= static_cast<std::uint64_t>(remaining[row]) << (10 + 3 * row);
  }
  if (const auto found = memo.find(key); found != memo.end()) return found->second;
  const int columns_left_after = kIsolateCount - depth - 1;
  std::uint64_t answer = 0;
  for (int index = start; index < static_cast<int>(allowed.size()); ++index) {
    const std::uint8_t column = allowed[index];
    bool feasible = true;
    for (int row = 0; row < kRowCount; ++row) {
      const int next = remaining[row] - ((column >> row) & 1);
      if (next < 0 || next > columns_left_after) {
        feasible = false;
        break;
      }
    }
    if (!feasible) continue;
    for (int row = 0; row < kRowCount; ++row) {
      remaining[row] -= (column >> row) & 1;
    }
    answer += CountColumns(allowed, index, depth + 1, remaining, memo);
    for (int row = 0; row < kRowCount; ++row) {
      remaining[row] += (column >> row) & 1;
    }
  }
  memo.emplace(key, answer);
  return answer;
}

void ProcessCoverFeasibleConfiguration(const Family& active_family) {
  std::array<std::uint8_t, 6> missed_rows{};
  for (std::size_t cover_index = 0; cover_index < g_graph_covers.size(); ++cover_index) {
    for (int row = 0; row < kRowCount; ++row) {
      if ((active_family[row] & g_graph_covers[cover_index]) == 0) {
        missed_rows[cover_index] |= std::uint8_t{1} << row;
      }
    }
  }
  std::vector<std::uint8_t> allowed;
  for (int column = 0; column < (1 << kRowCount); ++column) {
    if (std::popcount(static_cast<unsigned>(column)) == kRowCount) continue;
    bool safe = true;
    for (std::uint8_t rows : missed_rows) {
      if ((column & rows) == rows) {
        safe = false;
        break;
      }
    }
    if (safe) allowed.push_back(static_cast<std::uint8_t>(column));
  }
  std::array<int, kRowCount> remaining{};
  std::uint64_t allowed_word = 0;
  for (std::uint8_t column : allowed) allowed_word |= std::uint64_t{1} << column;
  std::uint32_t row_sum_code = 0;
  for (int row = 0; row < kRowCount; ++row) {
    remaining[row] = 5 - std::popcount(active_family[row]);
    row_sum_code |= static_cast<std::uint32_t>(remaining[row]) << (3 * row);
  }
  if (g_count_only) {
    const auto problem = std::make_pair(allowed_word, row_sum_code);
    auto found = g_distribution_count_cache.find(problem);
    if (found == g_distribution_count_cache.end()) {
      std::unordered_map<std::uint64_t, std::uint64_t> memo;
      const std::uint64_t count = CountColumns(allowed, 0, 0, remaining, memo);
      found = g_distribution_count_cache.emplace(problem, count).first;
    }
    g_distribution_paths += found->second;
    g_positive_distribution_problems += found->second > 0;
    return;
  }
  std::array<std::uint8_t, kIsolateCount> columns{};
  EnumerateColumns(active_family, allowed, 0, 0, remaining, columns);
}

void ConsumeActiveConfiguration(const Family& active_family) {
  ++g_cap_feasible_active_configurations;
  std::array<std::uint8_t, 6> missed_rows{};
  for (std::size_t cover_index = 0; cover_index < g_graph_covers.size(); ++cover_index) {
    for (int row = 0; row < kRowCount; ++row) {
      if ((active_family[row] & g_graph_covers[cover_index]) == 0) {
        missed_rows[cover_index] |= std::uint8_t{1} << row;
      }
    }
  }
  // Exact active fourth-vertex extension test.
  for (std::uint8_t rows : missed_rows) {
    for (int vertex : kActive) {
      bool hits_all = true;
      for (int row = 0; row < kRowCount; ++row) {
        if ((rows & (1u << row)) &&
            !(active_family[row] & (std::uint16_t{1} << vertex))) {
          hits_all = false;
          break;
        }
      }
      if (hits_all) return;
    }
  }
  ++g_active_cover_feasible_configurations;
  if (g_count_only) ProcessCoverFeasibleConfiguration(active_family);
  else g_stream_configurations.push_back(active_family);
}

void RunStreamWorkers() {
  const unsigned hardware = std::thread::hardware_concurrency();
  const unsigned worker_count = std::max(1u, std::min(12u, hardware ? hardware : 1u));
  std::atomic<std::size_t> next{0};
  std::vector<StreamResult> results(worker_count);
  std::vector<std::thread> workers;
  workers.reserve(worker_count);
  for (unsigned worker = 0; worker < worker_count; ++worker) {
    workers.emplace_back([worker, &next, &results]() {
      g_thread_result = &results[worker];
      while (true) {
        const std::size_t index = next.fetch_add(1, std::memory_order_relaxed);
        if (index >= g_stream_configurations.size()) break;
        ProcessCoverFeasibleConfiguration(g_stream_configurations[index]);
      }
    });
  }
  for (std::thread& worker : workers) worker.join();
  for (const StreamResult& result : results) {
    g_distribution_paths += result.distribution_paths;
    g_repeated_row_paths += result.repeated_row_paths;
    g_small_transversal_paths += result.small_transversal_paths;
    g_blocker_valid_paths += result.blocker_valid_paths;
    g_target_paths += result.target_paths;
    for (int facets = 0; facets <= kExpectedFiveCoverCount; ++facets) {
      g_facet_histogram[facets] += result.facet_histogram[facets];
    }
    g_target_families.insert(
        result.target_families.begin(), result.target_families.end());
  }
}

void EnumerateActiveConfigurations(int depth,
                                   int start,
                                   Family& family,
                                   std::array<int, kVertexCount>& incidences) {
  if (depth == kRowCount) {
    ConsumeActiveConfiguration(family);
    return;
  }
  for (int index = start; index < static_cast<int>(g_active_rows.size()); ++index) {
    const std::uint16_t row = g_active_rows[index];
    bool feasible = true;
    for (int vertex : kActive) {
      if ((row & (std::uint16_t{1} << vertex)) &&
          incidences[vertex] == kA5Caps[vertex]) {
        feasible = false;
        break;
      }
    }
    if (!feasible) continue;
    family[depth] = row;
    for (int vertex : kActive) incidences[vertex] += (row >> vertex) & 1;
    EnumerateActiveConfigurations(depth + 1, index, family, incidences);
    for (int vertex : kActive) incidences[vertex] -= (row >> vertex) & 1;
  }
}

void PrintVertices(std::uint16_t value) {
  std::cout << '[';
  bool first = true;
  for (int vertex = 0; vertex < kVertexCount; ++vertex) {
    if (value & (std::uint16_t{1} << vertex)) {
      if (!first) std::cout << ',';
      std::cout << vertex;
      first = false;
    }
  }
  std::cout << ']';
}

void PrintResult(bool print_targets) {
  std::cout << "{\n";
  std::cout << "  \"kind\": \"h9_14_record32_a0_lossless_scout\",\n";
  std::cout << "  \"theorem_certificate\": false,\n";
  std::cout << "  \"independent_active_mask_count\": " << g_active_rows.size() << ",\n";
  std::cout << "  \"minimal_graph_cover_count\": " << g_graph_covers.size() << ",\n";
  std::cout << "  \"active_row_multisets_after_caps_before_cover_filter\": "
            << g_cap_feasible_active_configurations << ",\n";
  std::cout << "  \"active_row_multisets_after_exact_active_cover_filter\": "
            << g_active_cover_feasible_configurations << ",\n";
  std::cout << "  \"column_multiset_paths\": " << g_distribution_paths << ",\n";
  std::cout << "  \"distinct_column_distribution_problem_count\": "
            << g_distribution_count_cache.size() << ",\n";
  std::cout << "  \"positive_column_distribution_configuration_count\": "
            << g_positive_distribution_problems << ",\n";
  std::cout << "  \"count_only\": " << (g_count_only ? "true" : "false") << ",\n";
  std::cout << "  \"literal_small_transversal_postcheck\": "
            << (g_literal_small_transversal_postcheck ? "true" : "false")
            << ",\n";
  std::cout << "  \"repeated_row_paths\": " << g_repeated_row_paths << ",\n";
  std::cout << "  \"small_transversal_paths\": " << g_small_transversal_paths << ",\n";
  std::cout << "  \"blocker_valid_paths\": " << g_blocker_valid_paths << ",\n";
  std::cout << "  \"facet_count_histogram\": {";
  bool first = true;
  for (int facet_count = 0; facet_count <= kExpectedFiveCoverCount; ++facet_count) {
    const std::uint64_t count = g_facet_histogram[facet_count];
    if (!count) continue;
    if (!first) std::cout << ',';
    std::cout << "\"" << facet_count << "\":" << count;
    first = false;
  }
  std::cout << "},\n";
  std::cout << "  \"scheduled_facet_count\": " << kScheduledFacetCount << ",\n";
  std::cout << "  \"scheduled_target_path_count\": " << g_target_paths << ",\n";
  std::cout << "  \"exact_unordered_target_family_count\": "
            << g_target_families.size();
  if (print_targets) {
    std::cout << ",\n  \"target_families\": [";
    bool first_family = true;
    for (const Family& family : g_target_families) {
      if (!first_family) std::cout << ',';
      std::cout << '[';
      for (int row = 0; row < kRowCount; ++row) {
        if (row) std::cout << ',';
        PrintVertices(family[row]);
      }
      std::cout << ']';
      first_family = false;
    }
    std::cout << ']';
  }
  std::cout << "\n}\n";
}

void WriteTargetFile(const std::string& path) {
  std::ofstream output(path);
  if (!output) throw std::runtime_error("could not open target output path");
  output << "{\n  \"scheduled_facet_count\": " << kScheduledFacetCount
         << ",\n  \"scheduled_target_path_count\": " << g_target_paths
         << ",\n  \"exact_unordered_target_family_count\": "
         << g_target_families.size() << ",\n  \"target_family_masks\": [\n";
  bool first_family = true;
  for (const Family& family : g_target_families) {
    if (!first_family) output << ",\n";
    output << "    [";
    for (int row = 0; row < kRowCount; ++row) {
      if (row) output << ',';
      output << family[row];
    }
    output << ']';
    first_family = false;
  }
  output << "\n  ]\n}\n";
}

}  // namespace

int main(int argc, char** argv) {
  try {
    bool print_targets = false;
    std::string target_output_path;
    for (int index = 1; index < argc; ++index) {
      const std::string option = argv[index];
      if (option == "--targets") print_targets = true;
      else if (option == "--target-output") {
        if (++index == argc) {
          std::cerr << "--target-output requires a path\n";
          return 2;
        }
        target_output_path = argv[index];
      }
      else if (option == "--count-only") g_count_only = true;
      else if (option == "--literal-small-transversal-postcheck") {
        g_literal_small_transversal_postcheck = true;
      } else {
        std::cerr << "usage: " << argv[0]
                  << " [--count-only] [--targets] [--target-output PATH] "
                     "[--literal-small-transversal-postcheck]\n";
        return 2;
      }
    }
    if (g_count_only && (print_targets || !target_output_path.empty() ||
                         g_literal_small_transversal_postcheck)) {
      std::cerr << "--count-only cannot be combined with streaming options\n";
      return 2;
    }
    InitializeRowsAndCovers();
    InitializeHitWords();
    Family family{};
    std::array<int, kVertexCount> incidences{};
    EnumerateActiveConfigurations(0, 0, family, incidences);
    if (!g_count_only) RunStreamWorkers();
    if (g_cap_feasible_active_configurations != 22513 ||
        g_active_cover_feasible_configurations != 14731) {
      throw std::runtime_error("record-32 active-row census changed");
    }
    if (!target_output_path.empty()) WriteTargetFile(target_output_path);
    PrintResult(print_targets);
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
