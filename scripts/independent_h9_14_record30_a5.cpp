// Independent exact a=5 incidence quotient for H9/14 record 30.
//
// The missing graph is 3K2 + 8K1.  Five unordered missing A4 rows are
// represented by their six distinguished active columns and the multiset of
// their eight isolate columns.  The program derives the 742-row ground set,
// enforces every vertex/edge/triple strict-link interval, quotients by
// S5 x S8, applies the forced-A5-shadow vertex Macaulay gate for the exact
// Dehn--Sommerville total b=38, and finally quotients survivors by
// Aut(3K2) x S8 x S5.  It is compact-only and makes no Coxeter-label or
// finite-volume assertion.

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cstdint>
#include <iostream>
#include <map>
#include <mutex>
#include <numeric>
#include <set>
#include <stdexcept>
#include <thread>
#include <tuple>
#include <unordered_set>
#include <vector>

namespace {

constexpr int N = 14;
constexpr int R = 5;
constexpr int I = 8;
constexpr std::array<int, 6> ACTIVE = {0, 1, 2, 11, 12, 13};
constexpr std::array<int, 8> ISOLATE = {3, 4, 5, 6, 7, 8, 9, 10};
constexpr std::array<std::uint16_t, 3> EDGES = {
    (1u << 0) | (1u << 11),
    (1u << 1) | (1u << 12),
    (1u << 2) | (1u << 13),
};

using Mask = std::uint16_t;
using Rows = std::array<Mask, R>;
using Columns = std::array<std::uint8_t, I>;
using Word = std::array<std::uint64_t, 22>;  // ceil(1372/64)

struct Gate {
  bool strict = false;
  bool aggregate = false;
  int shadow = 0;
  int additional = 0;
  auto operator<=>(const Gate&) const = default;
};

struct LinkWitness {
  Mask subset = 0;
  std::vector<int> h;
  std::vector<int> g;
};

struct Task {
  std::array<std::uint8_t, R> active_rows{};
};

struct Partial {
  std::uint64_t column_paths = 0;
  std::uint64_t equal_active_duplicates = 0;
  std::uint64_t repeated_rows = 0;
  std::uint64_t interval_s8_orbits = 0;
  std::uint64_t interval_labelled = 0;
  std::uint64_t viable_s8_orbits = 0;
  std::uint64_t viable_labelled = 0;
  std::map<Gate, std::uint64_t> gate_histogram;
  std::map<std::vector<int>, std::uint64_t> full_orbit_weights;
  std::map<std::vector<int>, Rows> representatives;
};

std::vector<std::uint8_t> g_active_masks;
std::vector<std::array<int, N>> g_graph_aut;
std::vector<Mask> g_fours;
std::array<int, 1 << N> g_four_index;
std::vector<Mask> g_fives;
std::vector<Word> g_shadow;
std::array<Word, N> g_vertex_fives;
std::vector<Task> g_tasks;

std::uint64_t Binomial(int n, int k) {
  if (k < 0 || k > n) return 0;
  k = std::min(k, n - k);
  std::uint64_t out = 1;
  for (int j = 1; j <= k; ++j) out = out * (n - k + j) / j;
  return out;
}

int Factorial(int n) {
  int out = 1;
  for (int j = 2; j <= n; ++j) out *= j;
  return out;
}

int MacaulayRaise(int value, int degree) {
  int remaining = value;
  int next_top = -1;
  int raised = 0;
  for (int lower = degree; lower >= 1; --lower) {
    int top = lower - 1;
    while ((next_top < 0 || top + 1 < next_top) &&
           Binomial(top + 1, lower) <= static_cast<std::uint64_t>(remaining)) {
      ++top;
    }
    remaining -= static_cast<int>(Binomial(top, lower));
    raised += static_cast<int>(Binomial(top + 1, lower + 1));
    next_top = top;
  }
  if (remaining != 0) throw std::runtime_error("bad Macaulay expansion");
  return raised;
}

bool Independent(Mask value) {
  for (Mask edge : EDGES) if ((value & edge) == edge) return false;
  return true;
}

Mask Transform(Mask value, const std::array<int, N>& action) {
  Mask out = 0;
  for (int v = 0; v < N; ++v) if (value & (Mask{1} << v)) {
    out |= Mask{1} << action[v];
  }
  return out;
}

int WordCount(const Word& word) {
  int out = 0;
  for (auto block : word) out += std::popcount(block);
  return out;
}

int IntersectCount(const Word& a, const Word& b) {
  int out = 0;
  for (int j = 0; j < 22; ++j) out += std::popcount(a[j] & b[j]);
  return out;
}

void Initialize() {
  g_four_index.fill(-1);
  // Active independent masks with at most two vertices.  The latter is the
  // exact cap-zero condition on every independent active triple.
  for (int value = 0; value < 64; ++value) {
    Mask full = 0;
    for (int j = 0; j < 6; ++j) if (value & (1 << j)) full |= Mask{1} << ACTIVE[j];
    if (Independent(full) && std::popcount(static_cast<unsigned>(value)) <= 2) {
      g_active_masks.push_back(static_cast<std::uint8_t>(value));
    }
  }
  if (g_active_masks.size() != 19) throw std::runtime_error("active mask count changed");

  // Literal A4 ground set: independent rows satisfying the single-row local
  // intervals.  For this graph that is precisely active cardinality <= 2.
  int independent_four_count = 0;
  for (Mask value = 0; value < (1u << N); ++value) {
    if (std::popcount(value) != 4 || !Independent(value)) continue;
    ++independent_four_count;
    int active_count = 0;
    for (int v : ACTIVE) active_count += (value >> v) & 1u;
    if (active_count <= 2) {
      g_four_index[value] = static_cast<int>(g_fours.size());
      g_fours.push_back(value);
    }
  }
  if (independent_four_count != 806 || g_fours.size() != 742) {
    throw std::runtime_error("A4 ground set changed");
  }

  for (Mask value = 0; value < (1u << N); ++value) {
    if (std::popcount(value) == 5 && Independent(value)) g_fives.push_back(value);
  }
  if (g_fives.size() != 1372) throw std::runtime_error("A5 ground set changed");
  g_shadow.resize(g_fours.size());
  for (int j = 0; j < static_cast<int>(g_fives.size()); ++j) {
    const Mask five = g_fives[j];
    for (int v = 0; v < N; ++v) if (five & (Mask{1} << v)) {
      g_vertex_fives[v][j / 64] |= std::uint64_t{1} << (j % 64);
    }
    for (int drop = 0; drop < N; ++drop) if (five & (Mask{1} << drop)) {
      const Mask four = five ^ (Mask{1} << drop);
      const int index = g_four_index[four];
      if (index >= 0) g_shadow[index][j / 64] |= std::uint64_t{1} << (j % 64);
    }
  }

  // Aut(3K2), derived literally inside S6.
  std::array<int, 6> images = {0, 1, 2, 3, 4, 5};
  const std::set<std::pair<int, int>> expected = {{0, 3}, {1, 4}, {2, 5}};
  do {
    std::set<std::pair<int, int>> moved;
    for (auto [a, b] : expected) moved.insert(std::minmax(images[a], images[b]));
    if (moved != expected) continue;
    std::array<int, N> action{};
    std::iota(action.begin(), action.end(), 0);
    for (int j = 0; j < 6; ++j) action[ACTIVE[j]] = ACTIVE[images[j]];
    g_graph_aut.push_back(action);
  } while (std::next_permutation(images.begin(), images.end()));
  if (g_graph_aut.size() != 48) throw std::runtime_error("graph automorphism order changed");

  // All active row multisets, exactly C(19+5-1,5)=33649 tasks.  Invalid
  // active column caps/codegrees are discarded before isolate enumeration.
  for (int a = 0; a < 19; ++a) for (int b = a; b < 19; ++b)
    for (int c = b; c < 19; ++c) for (int d = c; d < 19; ++d)
      for (int e = d; e < 19; ++e) {
        Task task{{g_active_masks[a], g_active_masks[b], g_active_masks[c],
                   g_active_masks[d], g_active_masks[e]}};
        std::array<std::uint8_t, 6> cols{};
        for (int v = 0; v < 6; ++v) for (int row = 0; row < R; ++row) {
          cols[v] |= ((task.active_rows[row] >> v) & 1u) << row;
        }
        if (std::ranges::any_of(cols, [](auto x) { return std::popcount(x) > 3; })) continue;
        bool good = true;
        for (int u = 0; u < 6; ++u) for (int v = u + 1; v < 6; ++v) {
          if ((u == 0 && v == 3) || (u == 1 && v == 4) || (u == 2 && v == 5)) continue;
          good = good && std::popcount(static_cast<unsigned>(cols[u] & cols[v])) <= 1;
        }
        if (good) g_tasks.push_back(task);
      }
}

std::uint8_t PermuteColumn(std::uint8_t column, const std::array<int, R>& order) {
  std::uint8_t out = 0;
  for (int row = 0; row < R; ++row) out |= ((column >> order[row]) & 1u) << row;
  return out;
}

std::vector<std::array<int, R>> RowStabilizer(const Task& task) {
  std::vector<std::array<int, R>> out;
  std::array<int, R> order = {0, 1, 2, 3, 4};
  do {
    bool good = true;
    for (int row = 0; row < R; ++row) {
      good = good && task.active_rows[order[row]] == task.active_rows[row];
    }
    if (good) out.push_back(order);
  } while (std::next_permutation(order.begin(), order.end()));
  return out;
}

std::uint64_t PackColumns(Columns columns) {
  std::sort(columns.begin(), columns.end());
  std::uint64_t out = 0;
  for (int j = 0; j < I; ++j) out |= std::uint64_t{columns[j]} << (R * j);
  return out;
}

std::uint64_t CanonicalColumns(const Columns& columns,
                               const std::vector<std::array<int, R>>& stabilizer) {
  std::uint64_t best = ~std::uint64_t{0};
  for (const auto& order : stabilizer) {
    Columns moved{};
    for (int j = 0; j < I; ++j) moved[j] = PermuteColumn(columns[j], order);
    best = std::min(best, PackColumns(moved));
  }
  return best;
}

Rows CompileRows(const Task& task, const Columns& columns) {
  Rows rows{};
  for (int row = 0; row < R; ++row) {
    for (int v = 0; v < 6; ++v) if ((task.active_rows[row] >> v) & 1u) {
      rows[row] |= Mask{1} << ACTIVE[v];
    }
    for (int j = 0; j < I; ++j) if ((columns[j] >> row) & 1u) {
      rows[row] |= Mask{1} << ISOLATE[j];
    }
  }
  return rows;
}

bool IsolateCapsHold(const Task& task, const Columns& columns) {
  std::array<std::uint8_t, 6> active_cols{};
  for (int v = 0; v < 6; ++v) for (int row = 0; row < R; ++row) {
    active_cols[v] |= ((task.active_rows[row] >> v) & 1u) << row;
  }
  // One active plus two isolates has cap 2.
  for (int v = 0; v < 6; ++v) for (int j = 0; j < I; ++j)
    for (int k = j + 1; k < I; ++k) {
      if (std::popcount(static_cast<unsigned>(active_cols[v] & columns[j] & columns[k])) > 2) {
        return false;
      }
    }
  // Three isolates have cap 3.
  for (int i = 0; i < I; ++i) for (int j = i + 1; j < I; ++j)
    for (int k = j + 1; k < I; ++k) {
      if (std::popcount(static_cast<unsigned>(columns[i] & columns[j] & columns[k])) > 3) {
        return false;
      }
    }
  return true;
}

Gate ShadowGate(const Rows& rows) {
  Word shadow{};
  std::array<int, N> incidence{};
  for (Mask row : rows) for (int v = 0; v < N; ++v) incidence[v] += (row >> v) & 1u;
  for (Mask row : rows) {
    const int index = g_four_index[row];
    if (index < 0) throw std::runtime_error("compiled row left A4 ground set");
    for (int j = 0; j < 22; ++j) shadow[j] |= g_shadow[index][j];
  }
  const int shadow_count = WordCount(shadow);
  const int additional = 38 - shadow_count;
  bool strict = additional >= 0;
  int lower_sum = 0;
  int upper_sum = 0;
  for (int v = 0; v < N; ++v) {
    const int base = std::find(ACTIVE.begin(), ACTIVE.end(), v) == ACTIVE.end() ? 8 : 4;
    const int previous = base - incidence[v];
    const int forced = IntersectCount(shadow, g_vertex_fives[v]);
    const int before = base + 6 * incidence[v] - forced;
    strict = strict && previous >= 1 && before >= 1;
    const int lower = previous >= 0 ? std::max(0, before - MacaulayRaise(previous, 3)) : 0;
    const int upper = before - 1;
    if (upper < 0) strict = false;
    lower_sum += lower;
    upper_sum += upper;
  }
  const bool aggregate = strict && lower_sum <= 5 * additional &&
                         5 * additional <= upper_sum;
  return {strict, aggregate, shadow_count, additional};
}

bool FaceBeforeAdditionalA5(Mask value, const Rows& rows) {
  if (!Independent(value)) return false;
  for (Mask four : rows) if ((value & four) == four) return false;
  return true;
}

LinkWitness FirstResidualLinkObstruction(const Rows& rows) {
  for (int subset_size = 1; subset_size <= 3; ++subset_size) {
    for (Mask subset = 1; subset < (1u << N); ++subset) {
      if (std::popcount(subset) != subset_size ||
          !FaceBeforeAdditionalA5(subset, rows)) continue;
      const int dimension = 9 - subset_size;
      const int high_order = 5 - subset_size;
      std::vector<int> f(high_order + 1);
      f[0] = 1;  // f_{-1}
      const Mask complement = static_cast<Mask>(((1u << N) - 1) ^ subset);
      for (Mask extension = complement;; extension = (extension - 1) & complement) {
        const int order = std::popcount(extension);
        if (order >= 1 && order <= high_order &&
            FaceBeforeAdditionalA5(subset | extension, rows)) {
          ++f[order];  // f[order] stores f_{order-1}
        }
        if (extension == 0) break;
      }
      std::vector<int> h(high_order + 1);
      for (int index = 0; index <= high_order; ++index) {
        for (int order = 0; order <= index; ++order) {
          const int sign = ((index - order) & 1) ? -1 : 1;
          h[index] += sign * static_cast<int>(Binomial(
              dimension - order, dimension - index)) * f[order];
        }
      }
      std::vector<int> g;
      for (int index = 1; index <= high_order; ++index) {
        g.push_back(h[index] - h[index - 1]);
      }
      const int previous = g[g.size() - 2];
      const int high = g.back();
      if (previous < 1 || high < 1 ||
          std::max(0, high - MacaulayRaise(previous, high_order - 1)) > high - 1) {
        return {subset, h, g};
      }
    }
  }
  return {};
}

std::uint64_t S8OrbitWeight(const Columns& columns,
                            const std::vector<std::array<int, R>>& stabilizer) {
  std::map<int, int> multiplicity;
  for (auto col : columns) ++multiplicity[col];
  int column_factor = 1;
  for (auto [col, count] : multiplicity) column_factor *= Factorial(count);
  int compatible_rows = 0;
  Columns sorted = columns;
  std::sort(sorted.begin(), sorted.end());
  for (const auto& order : stabilizer) {
    Columns moved{};
    for (int j = 0; j < I; ++j) moved[j] = PermuteColumn(columns[j], order);
    std::sort(moved.begin(), moved.end());
    compatible_rows += moved == sorted;
  }
  const int stabilizer_size = column_factor * compatible_rows;
  if (stabilizer_size <= 0 || Factorial(I) % stabilizer_size) {
    throw std::runtime_error("invalid S8 stabilizer");
  }
  return Factorial(I) / stabilizer_size;
}

std::vector<int> FullOrbitKey(const Rows& input) {
  std::vector<int> best;
  bool first = true;
  for (const auto& action : g_graph_aut) {
    Rows moved{};
    for (int row = 0; row < R; ++row) moved[row] = Transform(input[row], action);
    std::array<int, R> order = {0, 1, 2, 3, 4};
    do {
      std::vector<int> key;
      for (int row = 0; row < R; ++row) {
        int active = 0;
        for (int v : ACTIVE) if (moved[order[row]] & (Mask{1} << v)) active |= 1 << v;
        key.push_back(active);
      }
      std::vector<int> columns;
      for (int v : ISOLATE) {
        int col = 0;
        for (int row = 0; row < R; ++row) col |= ((moved[order[row]] >> v) & 1u) << row;
        columns.push_back(col);
      }
      std::sort(columns.begin(), columns.end());
      key.insert(key.end(), columns.begin(), columns.end());
      if (first || key < best) { best = std::move(key); first = false; }
    } while (std::next_permutation(order.begin(), order.end()));
  }
  return best;
}

void Consume(const Task& task, const Columns& columns,
             const std::vector<std::array<int, R>>& row_stabilizer,
             std::unordered_set<std::uint64_t>& seen, Partial& out) {
  ++out.column_paths;
  const std::uint64_t key = CanonicalColumns(columns, row_stabilizer);
  if (!seen.insert(key).second) { ++out.equal_active_duplicates; return; }
  if (!IsolateCapsHold(task, columns)) return;
  Rows rows = CompileRows(task, columns);
  std::sort(rows.begin(), rows.end());
  if (std::adjacent_find(rows.begin(), rows.end()) != rows.end()) {
    ++out.repeated_rows;
    return;
  }
  for (Mask row : rows) if (g_four_index[row] < 0) return;
  ++out.interval_s8_orbits;
  const std::uint64_t weight = S8OrbitWeight(columns, row_stabilizer);
  out.interval_labelled += weight;
  const Gate gate = ShadowGate(rows);
  out.gate_histogram[gate] += weight;
  if (!(gate.strict && gate.aggregate)) return;
  ++out.viable_s8_orbits;
  out.viable_labelled += weight;
  const auto full_key = FullOrbitKey(rows);
  out.full_orbit_weights[full_key] += weight;
  auto [where, inserted] = out.representatives.try_emplace(full_key, rows);
  if (!inserted && rows < where->second) where->second = rows;
}

void ProcessTask(const Task& task, Partial& out) {
  const auto stabilizer = RowStabilizer(task);
  std::unordered_set<std::uint64_t> seen;
  Columns columns{};
  std::array<int, R> remaining{};
  for (int row = 0; row < R; ++row) remaining[row] = 4 - std::popcount(task.active_rows[row]);
  auto visit = [&](auto&& self, int pattern, int left) -> void {
    if (pattern == (1 << R)) {
      if (left == 0 && std::ranges::all_of(remaining, [](int x) { return x == 0; })) {
        Consume(task, columns, stabilizer, seen, out);
      }
      return;
    }
    int maximum = left;
    for (int row = 0; row < R; ++row) if ((pattern >> row) & 1) {
      maximum = std::min(maximum, remaining[row]);
    }
    for (int count = 0; count <= maximum; ++count) {
      const int next_left = left - count;
      bool good = true;
      for (int row = 0; row < R; ++row) {
        const int next = remaining[row] - count * ((pattern >> row) & 1);
        good = good && next >= 0 && next <= next_left;
      }
      if (!good) continue;
      for (int row = 0; row < R; ++row) remaining[row] -= count * ((pattern >> row) & 1);
      const int start = I - left;
      for (int j = 0; j < count; ++j) columns[start + j] = pattern;
      self(self, pattern + 1, next_left);
      for (int row = 0; row < R; ++row) remaining[row] += count * ((pattern >> row) & 1);
    }
  };
  visit(visit, 0, I);
}

void Merge(Partial& target, const Partial& source) {
  target.column_paths += source.column_paths;
  target.equal_active_duplicates += source.equal_active_duplicates;
  target.repeated_rows += source.repeated_rows;
  target.interval_s8_orbits += source.interval_s8_orbits;
  target.interval_labelled += source.interval_labelled;
  target.viable_s8_orbits += source.viable_s8_orbits;
  target.viable_labelled += source.viable_labelled;
  for (const auto& [key, value] : source.gate_histogram) target.gate_histogram[key] += value;
  for (const auto& [key, value] : source.full_orbit_weights) target.full_orbit_weights[key] += value;
  for (const auto& [key, value] : source.representatives) {
    auto [where, inserted] = target.representatives.try_emplace(key, value);
    if (!inserted && value < where->second) where->second = value;
  }
}

void PrintRows(const Rows& rows) {
  std::cout << '[';
  for (int row = 0; row < R; ++row) {
    if (row) std::cout << ',';
    std::cout << '[';
    bool first = true;
    for (int v = 0; v < N; ++v) if (rows[row] & (Mask{1} << v)) {
      if (!first) std::cout << ',';
      std::cout << v;
      first = false;
    }
    std::cout << ']';
  }
  std::cout << ']';
}

}  // namespace

int main() {
  try {
    Initialize();
    const unsigned hw = std::thread::hardware_concurrency();
    const unsigned threads = std::max(1u, std::min(12u, hw ? hw : 1u));
    std::atomic<std::size_t> next{0};
    std::vector<Partial> partial(threads);
    std::vector<std::thread> workers;
    for (unsigned worker = 0; worker < threads; ++worker) {
      workers.emplace_back([&, worker]() {
        while (true) {
          const auto index = next.fetch_add(1, std::memory_order_relaxed);
          if (index >= g_tasks.size()) break;
          ProcessTask(g_tasks[index], partial[worker]);
        }
      });
    }
    for (auto& worker : workers) worker.join();
    Partial total;
    for (const auto& item : partial) Merge(total, item);
    std::map<std::uint64_t, int> orbit_histogram;
    for (const auto& [key, weight] : total.full_orbit_weights) ++orbit_histogram[weight];
    bool all_residual_obstructed = true;
    std::map<std::vector<int>, LinkWitness> witnesses;
    for (const auto& [key, rows] : total.representatives) {
      const auto witness = FirstResidualLinkObstruction(rows);
      all_residual_obstructed = all_residual_obstructed && witness.subset != 0;
      witnesses.emplace(key, witness);
    }
    std::cout << "{\n"
              << "  \"kind\":\"independent_h9_14_record30_a5\",\n"
              << "  \"compact_only\":true,\n"
              << "  \"schedule\":{\"a\":5,\"b\":38,\"facets\":122},\n"
              << "  \"independent_A4_count\":806,\n"
              << "  \"interval_A4_count\":742,\n"
              << "  \"active_mask_count\":" << g_active_masks.size() << ",\n"
              << "  \"active_multisets_after_caps\":" << g_tasks.size() << ",\n"
              << "  \"raw\":" << Binomial(742, 5) << ",\n"
              << "  \"column_paths\":" << total.column_paths << ",\n"
              << "  \"equal_active_duplicate_paths\":" << total.equal_active_duplicates << ",\n"
              << "  \"repeated_row_paths\":" << total.repeated_rows << ",\n"
              << "  \"interval_s8_orbits\":" << total.interval_s8_orbits << ",\n"
              << "  \"interval_feasible\":" << total.interval_labelled << ",\n"
              << "  \"viable_s8_orbits\":" << total.viable_s8_orbits << ",\n"
              << "  \"strict_aggregate_viable\":" << total.viable_labelled << ",\n"
              << "  \"full_orbit_count\":" << total.full_orbit_weights.size() << ",\n"
              << "  \"all_viable_orbits_residual_link_obstructed\":"
              << (all_residual_obstructed ? "true" : "false") << ",\n"
              << "  \"a5_compact_schedule_excluded\":"
              << (all_residual_obstructed ? "true" : "false") << ",\n"
              << "  \"full_orbit_weight_histogram\":{";
    bool first = true;
    for (auto [weight, count] : orbit_histogram) {
      if (!first) std::cout << ',';
      std::cout << '\"' << weight << "\":" << count;
      first = false;
    }
    std::cout << "},\n  \"orbits\":[";
    first = true;
    for (const auto& [key, weight] : total.full_orbit_weights) {
      if (!first) std::cout << ',';
      std::cout << "{\"weight\":" << weight << ",\"representative\":";
      PrintRows(total.representatives.at(key));
      const auto& witness = witnesses.at(key);
      std::cout << ",\"residual_link_witness\":{"
                << "\"subset\":[";
      bool first_vertex = true;
      for (int v = 0; v < N; ++v) if (witness.subset & (Mask{1} << v)) {
        if (!first_vertex) std::cout << ',';
        std::cout << v;
        first_vertex = false;
      }
      std::cout << "],\"h\":[";
      for (int j = 0; j < static_cast<int>(witness.h.size()); ++j) {
        if (j) std::cout << ',';
        std::cout << witness.h[j];
      }
      std::cout << "],\"g\":[";
      for (int j = 0; j < static_cast<int>(witness.g.size()); ++j) {
        if (j) std::cout << ',';
        std::cout << witness.g[j];
      }
      std::cout << "]}";
      std::cout << '}';
      first = false;
    }
    std::cout << "],\n  \"gate_histogram\":[";
    first = true;
    for (const auto& [gate, count] : total.gate_histogram) {
      if (!first) std::cout << ',';
      std::cout << "{\"strict\":" << (gate.strict ? "true" : "false")
                << ",\"aggregate\":" << (gate.aggregate ? "true" : "false")
                << ",\"shadow\":" << gate.shadow
                << ",\"additional\":" << gate.additional
                << ",\"count\":" << count << '}';
      first = false;
    }
    std::cout << "],\n"
              << "  \"scope\":\"Necessary compact incidence quotient only; no topology, Coxeter-label, or finite-volume conclusion.\"\n"
              << "}\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
