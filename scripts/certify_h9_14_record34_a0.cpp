// Exact solver-independent census for the a=0 branch of H9/14 record 34.
//
// Input is the graph6 stream produced by
//
//   genbg -q -d4:0 -D5:5 6 9 e:e
//
// for e=24,...,30.  The two colour classes are the six missing A5 rows
// and nine isolated vertices.  Row degree four means that the A5 has one
// active vertex; row degree five means it has none.  The program enumerates
// every locally admissible active colouring, tests blocker purity and the
// required count 114, and finally tests every six-set ridge equation.
// It imports and trusts no SAT/SMT backend.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

using Mask = uint16_t;

static constexpr std::array<int, 3> A_VERTICES = {0, 11, 13};
static constexpr std::array<int, 2> B_VERTICES = {1, 12};
static constexpr Mask ALL_VERTICES = (1u << 14) - 1;

struct Stats {
  uint64_t graphs = 0;
  uint64_t degree_graphs = 0;
  uint64_t colourings = 0;
  uint64_t distinct = 0;
  uint64_t no_small = 0;
  uint64_t target = 0;
  uint64_t pure = 0;
  uint64_t ridge = 0;
  std::map<std::tuple<int, int, int>, uint64_t> target_profiles;
  std::map<std::tuple<int, int, int>, uint64_t> ridge_profiles;
};

static bool graph_cover(Mask subset) {
  return std::popcount(unsigned(
             subset & ((1u << 0) | (1u << 11) | (1u << 13)))) >= 2 &&
         std::popcount(unsigned(subset & ((1u << 1) | (1u << 12)))) >= 1;
}

static bool hits(Mask subset, const std::array<Mask, 6>& missing_fives) {
  if (!graph_cover(subset)) return false;
  for (Mask face : missing_fives)
    if (!(subset & face)) return false;
  return true;
}

static bool no_small_transversal(const std::array<Mask, 6>& missing_fives) {
  // Every graph transversal contains one of these six inclusion-minimal
  // covers.  A transversal of size <=4 exists exactly when a base cover
  // already hits every A5 or the A5s it misses have a common vertex.
  for (int omitted_a : A_VERTICES) {
    for (int chosen_b : B_VERTICES) {
      Mask cover = Mask(1u << chosen_b);
      for (int a : A_VERTICES)
        if (a != omitted_a) cover |= Mask(1u << a);
      Mask common = ALL_VERTICES;
      bool missed = false;
      for (Mask face : missing_fives) {
        if (!(face & cover)) {
          common &= face;
          missed = true;
        }
      }
      if (!missed || common) return false;
    }
  }
  return true;
}

static int facet_count(const std::array<Mask, 6>& missing_fives) {
  int count = 0;
  for (Mask subset = 0; subset <= ALL_VERTICES; ++subset)
    if (std::popcount(unsigned(subset)) == 5 && hits(subset, missing_fives))
      ++count;
  return count;
}

static bool pure_of_order_five(const std::array<Mask, 6>& missing_fives) {
  for (Mask subset = 0; subset <= ALL_VERTICES; ++subset) {
    if (std::popcount(unsigned(subset)) <= 4 || !hits(subset, missing_fives))
      continue;
    bool minimal = true;
    Mask remaining = subset;
    while (remaining) {
      Mask bit = Mask(remaining & Mask(-int(remaining)));
      if (hits(Mask(subset ^ bit), missing_fives)) {
        minimal = false;
        break;
      }
      remaining ^= bit;
    }
    if (minimal && std::popcount(unsigned(subset)) != 5) return false;
  }
  return true;
}

static bool ridge_equations(const std::array<Mask, 6>& missing_fives) {
  // The complement of a ridge is a six-set.  Its degree is the number of
  // five-element transversals obtained by deleting one of its vertices.
  for (Mask subset = 0; subset <= ALL_VERTICES; ++subset) {
    if (std::popcount(unsigned(subset)) != 6) continue;
    int degree = 0;
    Mask remaining = subset;
    while (remaining) {
      Mask bit = Mask(remaining & Mask(-int(remaining)));
      degree += hits(Mask(subset ^ bit), missing_fives);
      remaining ^= bit;
    }
    if (degree != 0 && degree != 2) return false;
  }
  return true;
}

static std::string audit_summary(const std::array<Mask, 6>& missing_fives) {
  std::map<int, int> histogram;
  std::vector<Mask> facets;
  std::map<Mask, int> facet_indices;
  for (Mask subset = 0; subset <= ALL_VERTICES; ++subset) {
    if (std::popcount(unsigned(subset)) == 5 && hits(subset, missing_fives)) {
      facet_indices[Mask(ALL_VERTICES ^ subset)] = int(facets.size());
      facets.push_back(Mask(ALL_VERTICES ^ subset));
    }
  }

  std::vector<int> parent(facets.size());
  for (size_t index = 0; index < parent.size(); ++index)
    parent[index] = int(index);
  auto find = [&](int value) {
    while (parent[value] != value) {
      parent[value] = parent[parent[value]];
      value = parent[value];
    }
    return value;
  };
  auto join = [&](int first, int second) {
    first = find(first);
    second = find(second);
    if (first != second) parent[second] = first;
  };
  for (Mask subset = 0; subset <= ALL_VERTICES; ++subset) {
    if (std::popcount(unsigned(subset)) != 6) continue;
    std::vector<int> incident;
    Mask remaining = subset;
    while (remaining) {
      Mask bit = Mask(remaining & Mask(-int(remaining)));
      Mask transversal = Mask(subset ^ bit);
      if (hits(transversal, missing_fives))
        incident.push_back(facet_indices.at(Mask(ALL_VERTICES ^ transversal)));
      remaining ^= bit;
    }
    if (!incident.empty()) ++histogram[int(incident.size())];
    // This is the pseudomanifold dual graph used by the generic audit: only
    // degree-two ridges provide dual edges.  Singular ridges are recorded in
    // the histogram but do not spuriously join components.
    if (incident.size() == 2) join(incident[0], incident[1]);
  }
  std::set<int> components;
  for (size_t index = 0; index < facets.size(); ++index)
    components.insert(find(int(index)));

  std::ostringstream output;
  output << " ridge_hist=";
  for (auto [degree, count] : histogram)
    output << degree << ":" << count << ",";
  output << " degree_two_dual_components=" << components.size();
  return output.str();
}

static std::array<Mask, 6> make_missing_fives(
    const std::array<Mask, 6>& isolate_parts, Mask a_rows, Mask b0_rows,
    Mask b1_rows) {
  std::array<Mask, 6> missing_fives = isolate_parts;
  int a_index = 0;
  for (int row = 0; row < 6; ++row)
    if (a_rows & (1u << row))
      missing_fives[row] |= Mask(1u << A_VERTICES[a_index++]);
  for (int row = 0; row < 6; ++row) {
    if (b0_rows & (1u << row))
      missing_fives[row] |= Mask(1u << B_VERTICES[0]);
    if (b1_rows & (1u << row))
      missing_fives[row] |= Mask(1u << B_VERTICES[1]);
  }
  return missing_fives;
}

static void emit_candidate(const std::array<Mask, 6>& missing_fives, int a,
                           int b0, int b1, Mask a_rows, Mask b0_rows,
                           Mask b1_rows, const std::string& source,
                           const char* tag) {
  std::cerr << tag << " profile=" << a << "," << b0 << "," << b1
            << " source=" << source << " colouring=" << a_rows << ","
            << b0_rows << "," << b1_rows << " rows=";
  for (Mask face : missing_fives) {
    std::cerr << "[";
    bool first = true;
    for (int vertex = 0; vertex < 14; ++vertex) {
      if (face & (1u << vertex)) {
        if (!first) std::cerr << ",";
        std::cerr << vertex;
        first = false;
      }
    }
    std::cerr << "]";
  }
  std::cerr << audit_summary(missing_fives) << "\n";
}

static void audit_colouring(const std::array<Mask, 6>& isolate_parts,
                            Mask a_rows, Mask b0_rows, Mask b1_rows, int a,
                            int b0, int b1, const std::string& source,
                            Stats& stats) {
  ++stats.colourings;
  auto missing_fives =
      make_missing_fives(isolate_parts, a_rows, b0_rows, b1_rows);
  std::set<Mask> distinct(missing_fives.begin(), missing_fives.end());
  if (distinct.size() != 6) return;
  ++stats.distinct;
  if (!no_small_transversal(missing_fives)) return;
  ++stats.no_small;
  if (facet_count(missing_fives) != 114) return;
  ++stats.target;
  ++stats.target_profiles[{a, b0, b1}];
  emit_candidate(missing_fives, a, b0, b1, a_rows, b0_rows, b1_rows,
                 source, "TARGET");
  if (!pure_of_order_five(missing_fives)) return;
  ++stats.pure;
  if (!ridge_equations(missing_fives)) return;
  ++stats.ridge;
  ++stats.ridge_profiles[{a, b0, b1}];
  emit_candidate(missing_fives, a, b0, b1, a_rows, b0_rows, b1_rows,
                 source, "RIDGE");
}

static void colour_profile(const std::array<Mask, 6>& isolate_parts,
                           Mask active_rows, int a, int b0, int b1,
                           const std::string& source, Stats& stats) {
  for (Mask a_rows = active_rows;;
       a_rows = Mask((a_rows - 1) & active_rows)) {
    if (std::popcount(unsigned(a_rows)) == a) {
      Mask remaining_rows = Mask(active_rows ^ a_rows);
      for (Mask first_b_rows = remaining_rows;;
           first_b_rows = Mask((first_b_rows - 1) & remaining_rows)) {
        if (std::popcount(unsigned(first_b_rows)) == b0 &&
            std::popcount(unsigned(remaining_rows ^ first_b_rows)) == b1 &&
            (b0 != b1 || b0 == 0 ||
             first_b_rows < Mask(remaining_rows ^ first_b_rows))) {
          audit_colouring(isolate_parts, a_rows, first_b_rows,
                          Mask(remaining_rows ^ first_b_rows), a, b0, b1,
                          source, stats);
        }
        if (first_b_rows == 0) break;
      }
    }
    if (a_rows == 0) break;
  }
}

static void process_graph6(const std::string& line, Stats& stats) {
  if (line.empty() || line[0] == '>') return;
  ++stats.graphs;
  int vertex_count = int((unsigned char)line[0]) - 63;
  if (vertex_count != 15 || line.size() != 19)
    throw std::runtime_error("expected a short graph6 record on 15 vertices");
  std::vector<int> bits;
  for (size_t index = 1; index < line.size(); ++index) {
    int value = int((unsigned char)line[index]) - 63;
    for (int bit = 5; bit >= 0; --bit) bits.push_back((value >> bit) & 1);
  }

  int position = 0;
  std::array<int, 6> degrees{};
  std::array<int, 9> column_degrees{};
  std::array<Mask, 6> isolate_parts{};
  for (int second = 1; second < vertex_count; ++second) {
    for (int first = 0; first < second; ++first) {
      int edge = bits.at(position++);
      if (!edge) continue;
      int row = -1;
      int column = -1;
      if (first < 6 && second >= 6) {
        row = first;
        column = second - 6;
      } else if (second < 6 && first >= 6) {
        row = second;
        column = first - 6;
      } else {
        throw std::runtime_error("genbg stream contains a within-colour edge");
      }
      ++degrees[row];
      ++column_degrees[column];
      isolate_parts[row] |= Mask(1u << (2 + column));
    }
  }
  if (position != 105 ||
      std::any_of(bits.begin() + position, bits.end(), [](int bit) {
        return bit != 0;
      }))
    throw std::runtime_error("malformed graph6 payload or nonzero padding");
  if (std::any_of(column_degrees.begin(), column_degrees.end(),
                  [](int degree) { return degree > 5; }))
    throw std::runtime_error("isolate-column degree exceeds the local bound");

  Mask active_rows = 0;
  for (int row = 0; row < 6; ++row) {
    if (degrees[row] == 4)
      active_rows |= Mask(1u << row);
    else if (degrees[row] != 5)
      return;
  }
  ++stats.degree_graphs;
  int active_count = std::popcount(unsigned(active_rows));
  for (int a = 0; a <= 3; ++a) {
    for (int b0 = 0; b0 <= 2; ++b0) {
      for (int b1 = b0; b1 <= 2; ++b1) {
        if (a + b0 + b1 != active_count) continue;
        // Here every row is active.  Omitting the absent third A and choosing
        // either B leaves the other B's two rows with that B as a common
        // point, producing a forbidden four-transversal.
        if (a == 2 && b0 == 2 && b1 == 2) continue;
        colour_profile(isolate_parts, active_rows, a, b0, b1, line, stats);
      }
    }
  }
}

int main() {
  Stats stats;
  std::string line;
  while (std::getline(std::cin, line)) process_graph6(line, stats);
  std::cout << "graphs " << stats.graphs << " degree_graphs "
            << stats.degree_graphs << " colourings " << stats.colourings
            << " distinct " << stats.distinct << " no_small "
            << stats.no_small << " target " << stats.target << " pure "
            << stats.pure << " ridge " << stats.ridge << "\n";
  for (auto [profile, count] : stats.target_profiles) {
    auto [a, b0, b1] = profile;
    std::cout << "target_profile " << a << "," << b0 << "," << b1 << " "
              << count << "\n";
  }
  for (auto [profile, count] : stats.ridge_profiles) {
    auto [a, b0, b1] = profile;
    std::cout << "ridge_profile " << a << "," << b0 << "," << b1 << " "
              << count << "\n";
  }
}
