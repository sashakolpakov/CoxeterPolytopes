// Exact Vinberg edge/cusp verifier for the H13/50 and H14/98 sections.
//
// Input is: hyperbolic dimension, facet count, and the positive-normal Gram
// matrix.  The search enumerates every elliptic edge and finite vertex from
// the crystallographic finite Coxeter classification.  Independently it
// reaches every irreducible affine diagram as the first non-elliptic
// extension of an elliptic proper subset, checks determinant zero with
// overflow-guarded signed 128-bit Bareiss arithmetic, and forms its full
// zero-neighborhood.  A finite-volume verdict is warranted exactly when all
// neighborhoods are full rank parabolic and every edge has two finite/ideal
// endpoints.  The companion Python wrapper pins inputs and expected output.

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <queue>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

struct Bits {
  std::uint64_t lo = 0, hi = 0;
  bool operator==(const Bits& other) const {
    return lo == other.lo && hi == other.hi;
  }
};

struct BitsHash {
  std::size_t operator()(const Bits& value) const {
    return std::hash<std::uint64_t>{}(
        value.lo ^ (value.hi * 0x9e3779b97f4a7c15ULL));
  }
};

static int Popcount(Bits value) {
  return __builtin_popcountll(value.lo) + __builtin_popcountll(value.hi);
}
static bool Empty(Bits value) { return value.lo == 0 && value.hi == 0; }
static int Lowest(Bits value) {
  return value.lo ? __builtin_ctzll(value.lo)
                  : 64 + __builtin_ctzll(value.hi);
}
static void Clear(Bits& value, int bit) {
  if (bit < 64) value.lo &= ~(1ULL << bit);
  else value.hi &= ~(1ULL << (bit - 64));
}
static void Set(Bits& value, int bit) {
  if (bit < 64) value.lo |= 1ULL << bit;
  else value.hi |= 1ULL << (bit - 64);
}
static bool Contains(Bits value, int bit) {
  return bit < 64 ? value.lo & (1ULL << bit)
                  : value.hi & (1ULL << (bit - 64));
}
static Bits Intersect(Bits left, Bits right) {
  return Bits{left.lo & right.lo, left.hi & right.hi};
}
class Enumerator {
 public:
  int dimension = 0, count = 0, edge_order = 0, vertex_order = 0;
  std::vector<std::vector<long long>> gram;
  // 0 ultra, 1 infinity, and 2,3,4,6 finite Coxeter labels.
  std::vector<std::vector<int>> label;
  std::vector<Bits> compatible;
  std::vector<unsigned long long> face_counts, finite_special_counts;
  unsigned long long finite_search_nodes = 0, affine_search_nodes = 0;
  std::unordered_map<Bits, unsigned char, BitsHash> finite_endpoints;
  std::unordered_map<Bits, unsigned char, BitsHash> ideal_endpoints;
  std::unordered_set<Bits, BitsHash> affine_components;

  std::vector<int> Nodes(Bits mask) const {
    std::vector<int> result;
    while (!Empty(mask)) {
      int bit = Lowest(mask);
      Clear(mask, bit);
      result.push_back(bit);
    }
    return result;
  }

  bool Finite(const std::vector<int>& vertices) const {
    const int size = static_cast<int>(vertices.size());
    if (size <= 1) return true;
    std::vector<std::vector<int>> adjacency(size);
    for (int i = 0; i < size; ++i) {
      for (int j = 0; j < i; ++j) {
        int code = label[vertices[i]][vertices[j]];
        if (code == 3 || code == 4 || code == 6) {
          adjacency[i].push_back(j);
          adjacency[j].push_back(i);
        }
      }
    }
    std::vector<char> seen(size, false);
    for (int start = 0; start < size; ++start) {
      if (seen[start]) continue;
      std::vector<int> component;
      std::queue<int> queue;
      seen[start] = true;
      queue.push(start);
      while (!queue.empty()) {
        int current = queue.front(); queue.pop();
        component.push_back(current);
        for (int next : adjacency[current]) if (!seen[next]) {
          seen[next] = true; queue.push(next);
        }
      }
      if (component.size() == 1) continue;
      int degree_sum = 0, max_degree = 0, four_count = 0, six_count = 0;
      std::vector<int> degrees(size, 0);
      std::pair<int,int> four{-1,-1};
      for (int local : component) {
        degrees[local] = static_cast<int>(adjacency[local].size());
        degree_sum += degrees[local];
        max_degree = std::max(max_degree, degrees[local]);
        for (int next : adjacency[local]) if (local < next) {
          int code = label[vertices[local]][vertices[next]];
          if (code == 4) { ++four_count; four = {local,next}; }
          if (code == 6) ++six_count;
        }
      }
      if (degree_sum != 2 * (static_cast<int>(component.size()) - 1) ||
          max_degree > 3) return false;
      if (six_count) {
        if (component.size() != 2 || six_count != 1 || four_count) return false;
        continue;
      }
      if (four_count) {
        if (four_count != 1 || max_degree > 2) return false;
        if (component.size() == 2) continue;
        if (degrees[four.first] == 1 || degrees[four.second] == 1) continue;
        if (component.size() == 4 && degrees[four.first] == 2 &&
            degrees[four.second] == 2) continue;
        return false;
      }
      std::vector<int> branches;
      for (int local : component)
        if (degrees[local] == 3) branches.push_back(local);
      if (branches.empty()) continue;
      if (branches.size() != 1) return false;
      int center = branches[0];
      std::vector<int> arms;
      for (int neighbor : adjacency[center]) {
        int previous = center, current = neighbor, length = 1;
        while (degrees[current] == 2) {
          int following = adjacency[current][0] == previous
                              ? adjacency[current][1] : adjacency[current][0];
          previous = current; current = following; ++length;
        }
        arms.push_back(length);
      }
      std::sort(arms.begin(), arms.end());
      if ((arms[0] == 1 && arms[1] == 1) ||
          (arms == std::vector<int>{1,2,2}) ||
          (arms == std::vector<int>{1,2,3}) ||
          (arms == std::vector<int>{1,2,4})) continue;
      return false;
    }
    return true;
  }

  bool Connected(const std::vector<int>& vertices) const {
    if (vertices.empty()) return false;
    std::vector<char> seen(vertices.size(), false);
    std::queue<int> queue;
    seen[0] = true;
    queue.push(0);
    int total = 0;
    while (!queue.empty()) {
      int local = queue.front(); queue.pop(); ++total;
      for (int other = 0; other < static_cast<int>(vertices.size()); ++other) {
        if (!seen[other] && gram[vertices[local]][vertices[other]] != 0) {
          seen[other] = true; queue.push(other);
        }
      }
    }
    return total == static_cast<int>(vertices.size());
  }

  __int128 Determinant(const std::vector<int>& vertices) const {
    int size = static_cast<int>(vertices.size());
    std::vector<std::vector<__int128>> matrix(
        size, std::vector<__int128>(size));
    for (int i = 0; i < size; ++i)
      for (int j = 0; j < size; ++j)
        matrix[i][j] = gram[vertices[i]][vertices[j]];
    if (size == 0) return 1;
    __int128 previous = 1;
    int sign = 1;
    for (int pivot_index = 0; pivot_index + 1 < size; ++pivot_index) {
      int pivot_row = pivot_index;
      while (pivot_row < size && matrix[pivot_row][pivot_index] == 0)
        ++pivot_row;
      if (pivot_row == size) return 0;
      if (pivot_row != pivot_index) {
        std::swap(matrix[pivot_row], matrix[pivot_index]);
        sign = -sign;
      }
      __int128 pivot = matrix[pivot_index][pivot_index];
      for (int row = pivot_index + 1; row < size; ++row) {
        for (int column = pivot_index + 1; column < size; ++column) {
          __int128 first = 0, second = 0, numerator = 0;
          if (__builtin_mul_overflow(matrix[row][column], pivot, &first) ||
              __builtin_mul_overflow(matrix[row][pivot_index],
                                     matrix[pivot_index][column], &second) ||
              __builtin_sub_overflow(first, second, &numerator)) {
            std::cerr << "exact Bareiss arithmetic exceeded signed 128 bits\n";
            std::exit(4);
          }
          matrix[row][column] = numerator / previous;
        }
      }
      previous = pivot;
      if (previous == 0) return 0;
    }
    return sign * matrix[size - 1][size - 1];
  }

  bool Affine(const std::vector<int>& vertices) const {
    const int size = static_cast<int>(vertices.size());
    if (size < 3 || !Connected(vertices)) return false;
    std::vector<int> deleted;
    deleted.reserve(size - 1);
    for (int omit = 0; omit < size; ++omit) {
      deleted.clear();
      for (int index = 0; index < size; ++index)
        if (index != omit) deleted.push_back(vertices[index]);
      if (!Finite(deleted)) return false;
    }
    return Determinant(vertices) == 0;
  }

  void VisitFinite(std::vector<int>& vertices, Bits chosen, Bits candidates) {
    ++finite_search_nodes;
    int order = static_cast<int>(vertices.size());
    if (order > 0) ++face_counts[order];
    if (order == edge_order) finite_endpoints.try_emplace(chosen, 0);
    if (order == vertex_order) {
      Bits remaining = chosen;
      while (!Empty(remaining)) {
        int bit = Lowest(remaining); Clear(remaining, bit);
        Bits edge = chosen; Clear(edge, bit);
        ++finite_endpoints[edge];
      }
      return;
    }
    if (order < edge_order && Popcount(candidates) < edge_order - order) return;
    Bits remaining = candidates;
    while (!Empty(remaining)) {
      int vertex = Lowest(remaining); Clear(remaining, vertex);
      vertices.push_back(vertex); Set(chosen, vertex);
      if (Finite(vertices))
        VisitFinite(vertices, chosen, Intersect(remaining, compatible[vertex]));
      Clear(chosen, vertex); vertices.pop_back();
    }
  }

  void VisitAffine(std::vector<int>& vertices, Bits chosen, Bits candidates) {
    ++affine_search_nodes;
    int order = static_cast<int>(vertices.size());
    if (order > 0) ++finite_special_counts[order];
    Bits remaining = candidates;
    while (!Empty(remaining)) {
      int vertex = Lowest(remaining); Clear(remaining, vertex);
      vertices.push_back(vertex); Set(chosen, vertex);
      if (Finite(vertices)) {
        if (order + 1 <= dimension)
          VisitAffine(vertices, chosen, Intersect(remaining, compatible[vertex]));
      } else if (Affine(vertices)) {
        affine_components.insert(chosen);
      }
      Clear(chosen, vertex); vertices.pop_back();
    }
  }

  std::vector<Bits> Components(Bits mask) const {
    std::vector<Bits> result;
    while (!Empty(mask)) {
      int start = Lowest(mask);
      Bits component;
      Bits frontier;
      Set(frontier, start);
      Clear(mask, start);
      while (!Empty(frontier)) {
        int node = Lowest(frontier);
        Clear(frontier, node);
        Set(component, node);
        Bits scan = mask;
        while (!Empty(scan)) {
          int other = Lowest(scan);
          Clear(scan, other);
          if (gram[node][other] != 0) {
            Clear(mask, other);
            Set(frontier, other);
          }
        }
      }
      result.push_back(component);
    }
    return result;
  }

  Bits ZeroNeighborhood(Bits component) const {
    Bits result = component;
    for (int facet = 0; facet < count; ++facet) {
      if (Contains(component, facet)) continue;
      bool zero = true;
      Bits scan = component;
      while (!Empty(scan)) {
        int node = Lowest(scan); Clear(scan, node);
        if (gram[facet][node] != 0) { zero = false; break; }
      }
      if (zero) Set(result, facet);
    }
    return result;
  }

  void GenerateCuspEdges(const std::vector<Bits>& components, int position,
                         Bits edge) {
    if (position == static_cast<int>(components.size())) {
      auto found = finite_endpoints.find(edge);
      if (found == finite_endpoints.end()) {
        std::cerr << "parabolic deletion is not an elliptic edge\n";
        std::exit(3);
      }
      ++ideal_endpoints[edge];
      return;
    }
    Bits omitted = components[position];
    while (!Empty(omitted)) {
      int bit = Lowest(omitted); Clear(omitted, bit);
      Bits next = edge; Clear(next, bit);
      GenerateCuspEdges(components, position + 1, next);
    }
  }

  void Run() {
    std::cin >> dimension >> count;
    edge_order = dimension - 1; vertex_order = dimension;
    gram.assign(count, std::vector<long long>(count));
    for (auto& row : gram) for (auto& value : row) std::cin >> value;
    label.assign(count, std::vector<int>(count, 2));
    compatible.assign(count, Bits{});
    for (int left = 0; left < count; ++left) {
      for (int right = 0; right < count; ++right) if (left != right) {
        long long p = gram[left][right];
        long long numerator = p * p;
        long long denominator = gram[left][left] * gram[right][right];
        int code = 0;
        if (numerator == 0) code = 2;
        else if (4 * numerator == denominator) code = 3;
        else if (2 * numerator == denominator) code = 4;
        else if (4 * numerator == 3 * denominator) code = 6;
        else if (numerator == denominator) code = 1;
        else if (numerator > denominator) code = 0;
        else { std::cerr << "forbidden Coxeter ratio\n"; std::exit(2); }
        label[left][right] = code;
        if (code == 2 || code == 3 || code == 4 || code == 6)
          Set(compatible[left], right);
      }
    }

    face_counts.assign(dimension + 1, 0);
    finite_special_counts.assign(dimension + 1, 0);
    Bits all;
    for (int vertex = 0; vertex < count; ++vertex) Set(all, vertex);
    std::vector<int> vertices;
    VisitFinite(vertices, Bits{}, all);
    VisitAffine(vertices, Bits{}, all);
    for (int right = 0; right < count; ++right) {
      for (int left = 0; left < right; ++left) {
        if (label[left][right] == 1) {
          Bits pair; Set(pair, left); Set(pair, right);
          affine_components.insert(pair);
        }
      }
    }

    std::unordered_map<Bits, std::vector<Bits>, BitsHash> cusps;
    unsigned long long bad_neighborhoods = 0;
    for (Bits component : affine_components) {
      Bits neighborhood = ZeroNeighborhood(component);
      std::vector<Bits> components = Components(neighborhood);
      int rank = 0;
      bool valid = true;
      for (Bits part : components) {
        if (!affine_components.count(part)) { valid = false; break; }
        rank += Popcount(part) - 1;
      }
      if (!valid || rank != dimension - 1) {
        ++bad_neighborhoods;
      } else {
        cusps.try_emplace(neighborhood, components);
      }
    }

    for (const auto& item : cusps)
      GenerateCuspEdges(item.second, 0, item.first);

    std::array<std::array<unsigned long long, 4>, 4> endpoint_hist{};
    unsigned long long bad_edges = 0;
    for (const auto& item : finite_endpoints) {
      unsigned finite = item.second;
      unsigned ideal = ideal_endpoints[item.first];
      ++endpoint_hist[std::min<unsigned>(finite, 3)]
                     [std::min<unsigned>(ideal, 3)];
      if (finite + ideal != 2) ++bad_edges;
    }
    std::cout << "dimension " << dimension << " facets " << count
              << " finite_nodes " << finite_search_nodes
              << " affine_nodes " << affine_search_nodes
              << "\nedge_search_prefix_counts";
    for (int rank = 1; rank <= dimension; ++rank)
      std::cout << ' ' << face_counts[rank];
    std::cout << "\nfinite_special_counts";
    for (int rank = 1; rank <= dimension; ++rank)
      std::cout << ' ' << finite_special_counts[rank];
    std::cout << "\naffine_components " << affine_components.size()
              << " cusps " << cusps.size()
              << " bad_neighborhoods " << bad_neighborhoods
              << " bad_edges " << bad_edges << "\nendpoint_hist";
    for (unsigned finite = 0; finite < 4; ++finite)
      for (unsigned ideal = 0; ideal < 4; ++ideal)
        if (endpoint_hist[finite][ideal])
          std::cout << " f" << finite << "i" << ideal << ' '
                    << endpoint_hist[finite][ideal];
    std::cout << "\n";
  }
};

int main() { Enumerator enumerator; enumerator.Run(); }
