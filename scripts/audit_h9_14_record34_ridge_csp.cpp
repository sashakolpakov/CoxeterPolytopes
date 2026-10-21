// Independent 64-pattern ridge-CSP audit for H9/14 record 34, a=0.
//
// This is deliberately separate from certify_h9_14_record34_a0.cpp and its
// genbg census.  A missing A5 is one of six rows.  Its intersection with the
// nine isolated vertices is encoded by the nine column words p_j in {0,...,63}.
// The active graph is K3 + K2; local link intervals give at most one active
// label per row, A-incidences <=1, and B-incidences <=2.
//
// Once absence of transversals of order <=4 identifies five-transversals with
// facet complements, the complement of every candidate ridge is a six-set C
// and
//
//     d(C) = sum_{v in C} h(C - v)  is 0 or 2,
//
// where h is one exactly when its active part covers K3+K2 and it hits all six
// rows.  Only binary and ternary isolate constraints are needed below.
//
// The search is a direct finite backtrack: column words are nondecreasing to
// quotient S9, row capacities are exact, and compatibility is checked when a
// column is appended.  It uses no SAT/SMT backend and intentionally drops the
// column cap, row distinctness, target facet count, and explicit no-small-
// transversal constraints after deriving the ridge equation.  This is a
// relaxation, so its zero-survivor result is a valid independent obstruction.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <tuple>

using Profile = std::tuple<int, int, int>;

struct Search {
  Profile profile;
  std::array<int, 6> label{};  // -1 means no active label in the row.
  std::array<int, 6> target{};
  std::array<std::array<bool, 64>, 64> pair_ok{};
  std::array<int, 9> chosen{};
  std::uint64_t nodes = 0;
  bool use_ternary = false;

  static bool graph_cover(unsigned active) {
    return std::popcount(active & 7u) >= 2 &&
           std::popcount(active & 24u) >= 1;
  }

  bool h(unsigned active, unsigned isolate_union) const {
    if (!graph_cover(active)) return false;
    for (int row = 0; row < 6; ++row) {
      bool active_hit =
          label[row] >= 0 && (active & (1u << label[row]));
      bool isolate_hit = isolate_union & (1u << row);
      if (!active_hit && !isolate_hit) return false;
    }
    return true;
  }

  int binary_degree(int first, int second) const {
    // X consists of A0,A1,A2,B1, so it omits B0.  There are C(9,2)=36
    // isolate pairs.  Deleting B1 destroys graph coverage; its h-term is zero.
    constexpr unsigned active_four = 0b10111u;
    int degree = h(active_four, unsigned(first)) +
                 h(active_four, unsigned(second));
    for (int vertex : {0, 1, 2, 4}) {
      degree += h(active_four ^ (1u << vertex),
                  unsigned(first | second));
    }
    return degree;
  }

  bool ternary_ok(int first, int second, int third) const {
    // Each graph-covering active triple is minimal (two A's and one B), so
    // deleting an active vertex contributes zero.  The three remaining terms
    // correspond to deleting one of the isolate vertices.
    for (int omitted_a = 0; omitted_a < 3; ++omitted_a) {
      for (int b_vertex : {3, 4}) {
        unsigned active_three = 1u << b_vertex;
        for (int a_vertex = 0; a_vertex < 3; ++a_vertex) {
          if (a_vertex != omitted_a) active_three |= 1u << a_vertex;
        }
        int degree = h(active_three, unsigned(first | second)) +
                     h(active_three, unsigned(first | third)) +
                     h(active_three, unsigned(second | third));
        if (degree != 0 && degree != 2) return false;
      }
    }
    return true;
  }

  explicit Search(Profile value, bool ternary)
      : profile(value), use_ternary(ternary) {
    label.fill(-1);
    auto [a_count, b0_count, b1_count] = profile;
    int row = 0;
    for (int vertex = 0; vertex < a_count; ++vertex) label[row++] = vertex;
    for (int index = 0; index < b0_count; ++index) label[row++] = 3;
    for (int index = 0; index < b1_count; ++index) label[row++] = 4;
    for (int index = 0; index < 6; ++index) {
      target[index] = label[index] < 0 ? 5 : 4;
    }
    for (int first = 0; first < 64; ++first) {
      for (int second = 0; second < 64; ++second) {
        int degree = binary_degree(first, second);
        pair_ok[first][second] = degree == 0 || degree == 2;
      }
    }
  }

  bool visit(int depth, int start, std::array<int, 6>& capacities) {
    ++nodes;
    int columns_left = 9 - depth;
    for (int row = 0; row < 6; ++row) {
      if (capacities[row] < 0 || capacities[row] > columns_left) return false;
    }
    if (depth == 9) {
      return std::all_of(capacities.begin(), capacities.end(),
                         [](int value) { return value == 0; });
    }

    for (int pattern = start; pattern < 64; ++pattern) {
      bool admissible = true;
      for (int index = 0; index < depth; ++index) {
        if (!pair_ok[pattern][chosen[index]]) {
          admissible = false;
          break;
        }
      }
      if (use_ternary) {
        for (int first = 0; first < depth && admissible; ++first) {
          for (int second = first + 1; second < depth; ++second) {
            if (!ternary_ok(chosen[first], chosen[second], pattern)) {
              admissible = false;
              break;
            }
          }
        }
      }
      if (!admissible) continue;

      for (int row = 0; row < 6; ++row) {
        capacities[row] -= bool(pattern & (1u << row));
      }
      chosen[depth] = pattern;
      if (visit(depth + 1, pattern, capacities)) return true;
      for (int row = 0; row < 6; ++row) {
        capacities[row] += bool(pattern & (1u << row));
      }
    }
    return false;
  }

  bool solve() {
    auto capacities = target;
    return visit(0, 0, capacities);
  }
};

int main() {
  int profile_count = 0;
  int binary_unsat_count = 0;
  int final_unsat_count = 0;
  for (int a_count = 0; a_count <= 3; ++a_count) {
    for (int b0_count = 0; b0_count <= 2; ++b0_count) {
      for (int b1_count = b0_count; b1_count <= 2; ++b1_count) {
        if (a_count + b0_count + b1_count > 6) continue;
        ++profile_count;
        Profile profile{a_count, b0_count, b1_count};
        Search binary(profile, false);
        bool binary_sat = binary.solve();
        if (!binary_sat) ++binary_unsat_count;

        bool final_sat = binary_sat;
        std::uint64_t ternary_nodes = 0;
        if (profile == Profile{3, 0, 0}) {
          Search binary_ternary(profile, true);
          final_sat = binary_ternary.solve();
          ternary_nodes = binary_ternary.nodes;
        }
        if (!final_sat) ++final_unsat_count;

        bool expected_binary_sat = profile == Profile{3, 0, 0};
        if (binary_sat != expected_binary_sat || final_sat) {
          throw std::runtime_error("record-34 ridge profile table changed");
        }

        std::cout << a_count << "," << b0_count << "," << b1_count
                  << " binary=" << (binary_sat ? "SAT" : "UNSAT")
                  << " nodes=" << binary.nodes;
        if (profile == Profile{3, 0, 0}) {
          std::cout << " binary+ternary="
                    << (final_sat ? "SAT" : "UNSAT")
                    << " nodes=" << ternary_nodes;
        }
        std::cout << "\n";
      }
    }
  }
  if (profile_count != 23 || binary_unsat_count != 22 ||
      final_unsat_count != 23) {
    throw std::runtime_error("record-34 ridge aggregate changed");
  }
  std::cout << "profiles 23 binary_unsat 22 binary_sat 1 final_unsat 23\n";
}
