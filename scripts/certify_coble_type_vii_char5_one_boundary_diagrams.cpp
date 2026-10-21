// Exact finite-subdiagram enumerator for the 72-wall Type-VII Coble chamber.
//
// Input is a symmetric Coxeter-relation matrix.  Entries are 0 for an
// incompatible pair, 2 for an orthogonal pair, and 3/4/6 for a finite
// dihedral pair.  The first line gives the number of walls and the two orders
// to enumerate.  Output is one finite Coxeter mask per line, represented by
// two unsigned 64-bit words.  The final SUMMARY lines give the exhaustive
// compatible-clique and finite-diagram counts.

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using Mask = unsigned __int128;

int n = 0;
std::vector<std::vector<int>> relation;
std::vector<Mask> compatible;
std::vector<Mask> finite_nonorthogonal;
std::array<unsigned long long, 2> clique_counts{};
std::array<unsigned long long, 2> finite_counts{};
std::array<int, 2> requested_orders{};

struct MaskHash {
  std::size_t operator()(Mask value) const noexcept {
    const auto low = static_cast<std::uint64_t>(value);
    const auto high = static_cast<std::uint64_t>(value >> 64);
    return std::hash<std::uint64_t>{}(low ^
                                     (high + 0x9e3779b97f4a7c15ULL +
                                      (low << 6) + (low >> 2)));
  }
};

std::unordered_map<Mask, bool, MaskHash> finite_connected_cache;

int Popcount(Mask value) {
  return __builtin_popcountll(static_cast<std::uint64_t>(value)) +
         __builtin_popcountll(static_cast<std::uint64_t>(value >> 64));
}

int TrailingZeroes(Mask value) {
  const auto low = static_cast<std::uint64_t>(value);
  if (low != 0) return __builtin_ctzll(low);
  return 64 + __builtin_ctzll(static_cast<std::uint64_t>(value >> 64));
}

std::vector<int> Vertices(Mask mask) {
  std::vector<int> result;
  while (mask != 0) {
    const Mask bit = mask & (~mask + 1);
    result.push_back(TrailingZeroes(bit));
    mask ^= bit;
  }
  return result;
}

// Classify one connected component by the finite crystallographic Coxeter
// list.  We only need yes/no here.  This is equivalent to exact positive
// definiteness for a pair-compatible crystallographic diagram.
bool FiniteConnected(const std::vector<int>& nodes) {
  const int order = static_cast<int>(nodes.size());
  if (order == 1) return true;  // A1.

  std::vector<std::vector<int>> adjacency(order);
  int edge_count = 0;
  int four_count = 0;
  int six_count = 0;
  std::pair<int, int> four_edge{-1, -1};
  for (int left = 0; left < order; ++left) {
    for (int right = left + 1; right < order; ++right) {
      const int label = relation[nodes[left]][nodes[right]];
      if (label == 3 || label == 4 || label == 6) {
        adjacency[left].push_back(right);
        adjacency[right].push_back(left);
        ++edge_count;
        if (label == 4) {
          ++four_count;
          four_edge = {left, right};
        }
        if (label == 6) ++six_count;
      }
    }
  }
  if (edge_count != order - 1) return false;

  int maximum_degree = 0;
  for (const auto& neighbours : adjacency) {
    maximum_degree = std::max(maximum_degree,
                              static_cast<int>(neighbours.size()));
  }
  if (six_count != 0) return order == 2 && six_count == 1;

  if (maximum_degree <= 2) {
    if (four_count == 0) return true;  // A_order.
    if (four_count != 1) return false;
    const auto [left, right] = four_edge;
    // B_order has its unique four-edge at an end.  F4 has it in the middle.
    if (adjacency[left].size() == 1 || adjacency[right].size() == 1)
      return true;
    return order == 4;
  }

  if (four_count != 0 || maximum_degree > 3) return false;
  int branch = -1;
  for (int vertex = 0; vertex < order; ++vertex) {
    if (adjacency[vertex].size() == 3) {
      if (branch != -1) return false;
      branch = vertex;
    }
  }
  if (branch == -1) return false;
  std::vector<int> arms;
  for (const int neighbour : adjacency[branch]) {
    int previous = branch;
    int current = neighbour;
    int length = 1;
    while (adjacency[current].size() == 2) {
      const int following = adjacency[current][0] == previous
                                ? adjacency[current][1]
                                : adjacency[current][0];
      previous = current;
      current = following;
      ++length;
    }
    arms.push_back(length);
  }
  std::sort(arms.begin(), arms.end());
  if (arms.size() != 3) return false;
  if (arms[0] == 1 && arms[1] == 1) return true;  // D_order.
  return arms == std::vector<int>({1, 2, 2}) ||
         arms == std::vector<int>({1, 2, 3}) ||
         arms == std::vector<int>({1, 2, 4});  // E6, E7, E8.
}

bool FiniteDiagram(Mask mask) {
  Mask remaining = mask;
  while (remaining != 0) {
    const Mask start = remaining & (~remaining + 1);
    remaining ^= start;
    Mask component = 0;
    Mask frontier = start;
    while (frontier != 0) {
      const Mask bit = frontier & (~frontier + 1);
      frontier ^= bit;
      const int vertex = TrailingZeroes(bit);
      component |= bit;
      const Mask fresh = finite_nonorthogonal[vertex] & remaining;
      remaining ^= fresh;
      frontier |= fresh;
    }
    const auto cached = finite_connected_cache.find(component);
    if (cached != finite_connected_cache.end()) {
      if (!cached->second) return false;
      continue;
    }
    const bool finite = FiniteConnected(Vertices(component));
    finite_connected_cache.emplace(component, finite);
    if (!finite) return false;
  }
  return true;
}

void EmitMask(int order, Mask mask) {
  const auto low = static_cast<std::uint64_t>(mask);
  const auto high = static_cast<std::uint64_t>(mask >> 64);
  std::cout << order << ' ' << low << ' ' << high << '\n';
}

void Enumerate(Mask chosen, Mask candidates, int needed, int slot) {
  if (needed == 0) {
    ++clique_counts[slot];
    if (FiniteDiagram(chosen)) {
      ++finite_counts[slot];
      EmitMask(requested_orders[slot], chosen);
    }
    return;
  }
  while (Popcount(candidates) >= needed) {
    const Mask bit = candidates & (~candidates + 1);
    candidates ^= bit;
    const int vertex = TrailingZeroes(bit);
    Enumerate(chosen | bit, candidates & compatible[vertex], needed - 1,
              slot);
  }
}

}  // namespace

int main() {
  if (!(std::cin >> n >> requested_orders[0] >> requested_orders[1])) return 2;
  if (n <= 0 || n > 128 || requested_orders[0] <= 0 ||
      requested_orders[1] <= 0)
    return 3;
  relation.assign(n, std::vector<int>(n));
  compatible.assign(n, 0);
  finite_nonorthogonal.assign(n, 0);
  for (int left = 0; left < n; ++left) {
    for (int right = 0; right < n; ++right) {
      if (!(std::cin >> relation[left][right])) return 4;
      const int label = relation[left][right];
      if (left != right && (label == 2 || label == 3 || label == 4 || label == 6))
        compatible[left] |= (Mask{1} << right);
      if (left != right && (label == 3 || label == 4 || label == 6))
        finite_nonorthogonal[left] |= (Mask{1} << right);
    }
  }
  const Mask universe = n == 128 ? ~Mask{0} : ((Mask{1} << n) - 1);
  for (int slot = 0; slot < 2; ++slot)
    Enumerate(0, universe, requested_orders[slot], slot);
  for (int slot = 0; slot < 2; ++slot) {
    std::cout << "SUMMARY " << requested_orders[slot] << ' '
              << clique_counts[slot] << ' ' << finite_counts[slot] << '\n';
  }
  return 0;
}
