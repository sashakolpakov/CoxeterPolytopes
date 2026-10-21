// Solver-independent positive-a census for H9/14 record 23.
//
// This executable handles the Dehn--Sommerville branches with one, two, or
// three missing four-faces.  It uses only finite integer/bit-mask enumeration:
// no SAT/SMT backend is imported or trusted.
//
// Vertices 0,1,2 form the size-three side A of the missing K_{3,2}; vertices
// 12,13 form its size-two side B; vertices 3,...,11 are the nine isolates I.
// Strict local g4 forces every extra minimal nonface into W=A union I, and
// every missing four-face must meet A.  Facet complements are then counted as
// tau_3(G)+tau_2(Z), as explained in results/h9_14_record23_positive_a.md.

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <vector>

using namespace std;

namespace {

constexpr int ACTIVE_MASK = 0b111;
constexpr int W_SIZE = 12;
constexpr int FULL_W_MASK = (1 << W_SIZE) - 1;

int popcount(int value) {
    return __builtin_popcount(static_cast<unsigned>(value));
}

vector<int> masks_of_size(int ground_size, int size, int offset = 0) {
    vector<int> result;
    for (int mask = 0; mask < (1 << ground_size); ++mask) {
        if (popcount(mask) == size) result.push_back(mask << offset);
    }
    return result;
}

using Bits = array<uint64_t, 4>;

Bits all_bits(int count) {
    Bits result{};
    for (int index = 0; index < count; ++index) {
        result[index / 64] |= uint64_t{1} << (index % 64);
    }
    return result;
}

Bits bit_and(Bits first, const Bits& second) {
    for (int index = 0; index < 4; ++index) first[index] &= second[index];
    return first;
}

int bit_count(const Bits& bits) {
    int result = 0;
    for (uint64_t word : bits) result += __builtin_popcountll(word);
    return result;
}

// Canonical key of an unordered family of isolate-only five-sets.  The
// distinguished isolate block is I cap Q, where Q is the missing four-face;
// its size is q_isolate_count.  The stabilizer is therefore
// S_z(rows) x S_q x S_{9-q}.  Isolate permutations are quotiented by the two
// lists of row-membership-column multiplicities.
vector<int> isolate_family_key(
    const vector<int>& family,
    int q_isolate_count
) {
    const int row_count = static_cast<int>(family.size());
    vector<int> permutation(row_count);
    for (int index = 0; index < row_count; ++index) permutation[index] = index;

    vector<int> best;
    bool first = true;
    do {
        vector<int> key(2 * (1 << row_count));
        for (int vertex = 3; vertex < W_SIZE; ++vertex) {
            int old_column = 0;
            for (int row = 0; row < row_count; ++row) {
                if (family[row] & (1 << vertex)) old_column |= 1 << row;
            }
            int new_column = 0;
            for (int row = 0; row < row_count; ++row) {
                if (old_column & (1 << permutation[row])) {
                    new_column |= 1 << row;
                }
            }
            const int block = vertex < 3 + q_isolate_count ? 0 : 1;
            ++key[block * (1 << row_count) + new_column];
        }
        if (first || key < best) {
            best = key;
            first = false;
        }
    } while (next_permutation(permutation.begin(), permutation.end()));
    return best;
}

struct A1Case {
    int active_size = 0;
    int q = 0;
    int q_isolate_count = 0;
    array<int, W_SIZE> chosen_a5_capacity{};
    vector<int> active_fives;
    vector<int> isolate_fives;
    vector<int> triples;
    vector<int> pairs;
    vector<int> isolate_pairs;
    array<Bits, 1 << W_SIZE> hit_triples{};
    array<Bits, 1 << W_SIZE> hit_pairs{};
    Bits all_triples{};
    Bits all_pairs{};
    map<int, uint64_t> facet_histogram;
    uint64_t valid_family_count = 0;
    uint64_t target_family_count = 0;
    int maximum_facet_count = -1;
    vector<int> maximum_family;
};

void evaluate_a1_family(
    A1Case& data,
    const vector<int>& isolate_family
) {
    const int active_face_count = 5 - static_cast<int>(isolate_family.size());
    Bits surviving_triples = bit_and(data.all_triples, data.hit_triples[data.q]);
    Bits surviving_pairs = bit_and(data.all_pairs, data.hit_pairs[data.q]);
    for (int face : isolate_family) {
        surviving_triples = bit_and(surviving_triples, data.hit_triples[face]);
        surviving_pairs = bit_and(surviving_pairs, data.hit_pairs[face]);
    }

    int tau2 = 0;
    for (int pair : data.isolate_pairs) {
        bool hits = true;
        for (int face : isolate_family) {
            if (!(pair & face)) {
                hits = false;
                break;
            }
        }
        tau2 += hits;
    }

    array<int, W_SIZE> used{};
    for (int face : isolate_family) {
        for (int vertex = 0; vertex < W_SIZE; ++vertex) {
            used[vertex] += (face >> vertex) & 1;
        }
    }

    auto record = [&](const Bits& triples,
                      const Bits& pairs,
                      const vector<int>& active_family) {
        if (bit_count(pairs) != 0) return;
        const int facet_count = bit_count(triples) + tau2;
        ++data.facet_histogram[facet_count];
        ++data.valid_family_count;
        data.target_family_count += facet_count == 96;
        if (facet_count > data.maximum_facet_count) {
            data.maximum_facet_count = facet_count;
            data.maximum_family = isolate_family;
            data.maximum_family.insert(
                data.maximum_family.end(),
                active_family.begin(),
                active_family.end()
            );
        }
    };

    if (active_face_count == 0) {
        record(surviving_triples, surviving_pairs, {});
        return;
    }

    if (active_face_count == 1) {
        for (int face : data.active_fives) {
            bool capacity_ok = true;
            for (int vertex = 0; vertex < W_SIZE; ++vertex) {
                if (used[vertex] + ((face >> vertex) & 1)
                    > data.chosen_a5_capacity[vertex]) {
                    capacity_ok = false;
                    break;
                }
            }
            if (!capacity_ok) continue;
            record(
                bit_and(surviving_triples, data.hit_triples[face]),
                bit_and(surviving_pairs, data.hit_pairs[face]),
                {face}
            );
        }
        return;
    }

    for (int first_index = 0;
         first_index < static_cast<int>(data.active_fives.size());
         ++first_index) {
        const int first = data.active_fives[first_index];
        bool first_capacity_ok = true;
        for (int vertex = 0; vertex < W_SIZE; ++vertex) {
            if (used[vertex] + ((first >> vertex) & 1)
                > data.chosen_a5_capacity[vertex]) {
                first_capacity_ok = false;
                break;
            }
        }
        if (!first_capacity_ok) continue;
        const Bits first_triples = bit_and(
            surviving_triples, data.hit_triples[first]
        );
        const Bits first_pairs = bit_and(surviving_pairs, data.hit_pairs[first]);

        for (int second_index = first_index + 1;
             second_index < static_cast<int>(data.active_fives.size());
             ++second_index) {
            const int second = data.active_fives[second_index];
            bool capacity_ok = true;
            for (int vertex = 0; vertex < W_SIZE; ++vertex) {
                if (used[vertex]
                        + ((first >> vertex) & 1)
                        + ((second >> vertex) & 1)
                    > data.chosen_a5_capacity[vertex]) {
                    capacity_ok = false;
                    break;
                }
            }
            if (!capacity_ok) continue;
            record(
                bit_and(first_triples, data.hit_triples[second]),
                bit_and(first_pairs, data.hit_pairs[second]),
                {first, second}
            );
        }
    }
}

void generate_a1_isolate_orbits(
    A1Case& data,
    int family_size,
    map<vector<int>, vector<int>>& orbits
) {
    vector<int> selected;
    array<int, W_SIZE> degrees{};
    auto recurse = [&](auto&& self, int start, int left, int common) -> void {
        if (left == 0) {
            if (common != 0) return;
            vector<int> key = isolate_family_key(
                selected, data.q_isolate_count
            );
            orbits.emplace(std::move(key), selected);
            return;
        }
        for (int index = start;
             index <= static_cast<int>(data.isolate_fives.size()) - left;
             ++index) {
            const int face = data.isolate_fives[index];
            bool capacity_ok = true;
            for (int vertex = 3; vertex < W_SIZE; ++vertex) {
                if ((face & (1 << vertex))
                    && degrees[vertex] >= data.chosen_a5_capacity[vertex]) {
                    capacity_ok = false;
                    break;
                }
            }
            if (!capacity_ok) continue;
            for (int vertex = 3; vertex < W_SIZE; ++vertex) {
                degrees[vertex] += (face >> vertex) & 1;
            }
            selected.push_back(face);
            self(
                self,
                index + 1,
                left - 1,
                selected.size() == 1 ? face : (common & face)
            );
            selected.pop_back();
            for (int vertex = 3; vertex < W_SIZE; ++vertex) {
                degrees[vertex] -= (face >> vertex) & 1;
            }
        }
    };
    recurse(recurse, 0, family_size, 0);
}

void enumerate_a1_isolate_fives(A1Case& data) {
    // For three or four isolate-only faces, quotient by the exact stabilizer
    // of q before adding all remaining active five-faces.  For five faces no
    // active face remains, so a direct capacity-pruned enumeration is faster.
    for (int family_size : {3, 4}) {
        map<vector<int>, vector<int>> orbits;
        generate_a1_isolate_orbits(data, family_size, orbits);
        const array<array<int, 2>, 3> expected_orbits{{
            {{56, 1976}}, {{33, 1009}}, {{14, 334}}
        }};
        const int expected = expected_orbits[data.active_size - 1][family_size - 3];
        if (static_cast<int>(orbits.size()) != expected) {
            throw runtime_error("a=1 isolate-family orbit count changed");
        }
        cout << "a=1 active_size=" << data.active_size
             << " isolate_family_size=" << family_size
             << " isolate_orbits=" << orbits.size() << '\n';
        for (const auto& [key, representative] : orbits) {
            (void)key;
            evaluate_a1_family(data, representative);
        }
    }

    vector<int> selected;
    array<int, W_SIZE> degrees{};
    auto recurse = [&](auto&& self, int start, int left, int common) -> void {
        if (left == 0) {
            if (common == 0) evaluate_a1_family(data, selected);
            return;
        }
        for (int index = start;
             index <= static_cast<int>(data.isolate_fives.size()) - left;
             ++index) {
            const int face = data.isolate_fives[index];
            bool capacity_ok = true;
            for (int vertex = 3; vertex < W_SIZE; ++vertex) {
                if ((face & (1 << vertex))
                    && degrees[vertex] >= data.chosen_a5_capacity[vertex]) {
                    capacity_ok = false;
                    break;
                }
            }
            if (!capacity_ok) continue;
            for (int vertex = 3; vertex < W_SIZE; ++vertex) {
                degrees[vertex] += (face >> vertex) & 1;
            }
            selected.push_back(face);
            self(
                self,
                index + 1,
                left - 1,
                selected.size() == 1 ? face : (common & face)
            );
            selected.pop_back();
            for (int vertex = 3; vertex < W_SIZE; ++vertex) {
                degrees[vertex] -= (face >> vertex) & 1;
            }
        }
    };
    recurse(recurse, 0, 5, 0);
}

A1Case make_a1_case(int active_size) {
    A1Case data;
    data.active_size = active_size;
    data.q_isolate_count = 4 - active_size;
    for (int vertex = 0; vertex < active_size; ++vertex) {
        data.q |= 1 << vertex;
    }
    for (int vertex = 3; vertex < 3 + data.q_isolate_count; ++vertex) {
        data.q |= 1 << vertex;
    }

    set<int> forced_fives;
    for (int vertex = 0; vertex < W_SIZE; ++vertex) {
        if (!(data.q & (1 << vertex))) {
            forced_fives.insert(data.q | (1 << vertex));
        }
    }
    if (forced_fives.size() != 8) {
        throw runtime_error("one A4 must force eight five-nonfaces");
    }
    for (int vertex = 0; vertex < W_SIZE; ++vertex) {
        int forced_incidence = 0;
        for (int face : forced_fives) {
            forced_incidence += (face >> vertex) & 1;
        }
        const int baseline = vertex < 3 ? 5 : 6;
        data.chosen_a5_capacity[vertex] =
            baseline
            + 6 * ((data.q >> vertex) & 1)
            - forced_incidence
            - 1;
    }

    data.isolate_fives = masks_of_size(9, 5, 3);
    data.triples = masks_of_size(W_SIZE, 3);
    data.pairs = masks_of_size(W_SIZE, 2);
    data.isolate_pairs = masks_of_size(9, 2, 3);
    data.all_triples = all_bits(static_cast<int>(data.triples.size()));
    data.all_pairs = all_bits(static_cast<int>(data.pairs.size()));
    for (int mask = 0; mask <= FULL_W_MASK; ++mask) {
        for (int index = 0; index < static_cast<int>(data.triples.size()); ++index) {
            if (mask & data.triples[index]) {
                data.hit_triples[mask][index / 64]
                    |= uint64_t{1} << (index % 64);
            }
        }
        for (int index = 0; index < static_cast<int>(data.pairs.size()); ++index) {
            if (mask & data.pairs[index]) {
                data.hit_pairs[mask][index / 64]
                    |= uint64_t{1} << (index % 64);
            }
        }
        if (popcount(mask) == 5
            && (mask & ACTIVE_MASK)
            && ((mask & data.q) != data.q)) {
            data.active_fives.push_back(mask);
        }
    }
    return data;
}

void run_a1() {
    const array<uint64_t, 3> expected_valid{{
        133726, 508884, 9515916
    }};
    for (int active_size = 1; active_size <= 3; ++active_size) {
        A1Case data = make_a1_case(active_size);
        enumerate_a1_isolate_fives(data);
        cout << "a=1 active_size=" << active_size
             << " valid_families=" << data.valid_family_count
             << " target_96=" << data.target_family_count
             << " histogram=";
        for (const auto& [facet_count, multiplicity] : data.facet_histogram) {
            cout << facet_count << ':' << multiplicity << ',';
        }
        cout << '\n';
        cout << "a=1 active_size=" << active_size
             << " maximum=" << data.maximum_facet_count
             << " q=" << data.q
             << " maximum_family_masks=";
        for (int face : data.maximum_family) cout << face << ',';
        cout << '\n';
        if (data.target_family_count != 0) {
            throw runtime_error("a=1 target facet count 96 survived");
        }
        if (data.valid_family_count != expected_valid[active_size - 1]
            || data.maximum_facet_count != 92) {
            throw runtime_error("a=1 representative-state census changed");
        }
    }
}

vector<array<int, 3>> active_permutations() {
    vector<array<int, 3>> result;
    array<int, 3> permutation{0, 1, 2};
    do {
        result.push_back(permutation);
    } while (next_permutation(permutation.begin(), permutation.end()));
    return result;
}

int permute_active_mask(int mask, const array<int, 3>& permutation) {
    int result = 0;
    for (int vertex = 0; vertex < 3; ++vertex) {
        if (mask & (1 << vertex)) result |= 1 << permutation[vertex];
    }
    return result;
}

vector<int> canonical_a4_key(const vector<int>& family) {
    const int row_count = static_cast<int>(family.size());
    vector<int> active_parts(row_count);
    vector<int> column_counts(1 << row_count);
    for (int row = 0; row < row_count; ++row) {
        active_parts[row] = family[row] & ACTIVE_MASK;
    }
    for (int vertex = 3; vertex < W_SIZE; ++vertex) {
        int column = 0;
        for (int row = 0; row < row_count; ++row) {
            if (family[row] & (1 << vertex)) column |= 1 << row;
        }
        ++column_counts[column];
    }

    vector<int> row_permutation(row_count);
    for (int index = 0; index < row_count; ++index) row_permutation[index] = index;
    vector<int> best;
    bool first = true;
    const vector<array<int, 3>> permutations = active_permutations();
    do {
        for (const auto& active_permutation : permutations) {
            vector<int> key;
            for (int new_row = 0; new_row < row_count; ++new_row) {
                key.push_back(permute_active_mask(
                    active_parts[row_permutation[new_row]], active_permutation
                ));
            }
            vector<int> transformed_counts(1 << row_count);
            for (int old_column = 0;
                 old_column < (1 << row_count);
                 ++old_column) {
                int new_column = 0;
                for (int new_row = 0; new_row < row_count; ++new_row) {
                    if (old_column & (1 << row_permutation[new_row])) {
                        new_column |= 1 << new_row;
                    }
                }
                transformed_counts[new_column] += column_counts[old_column];
            }
            key.insert(key.end(), transformed_counts.begin(), transformed_counts.end());
            if (first || key < best) {
                best = std::move(key);
                first = false;
            }
        }
    } while (next_permutation(row_permutation.begin(), row_permutation.end()));
    return best;
}

bool locally_admissible_a4_family(const vector<int>& family) {
    // For a<=3, the only nonautomatic positive-g conditions are:
    //   * a pair containing an active A-vertex lies in at most two A4s;
    //   * a triple containing an active A-vertex lies in at most one A4.
    // The B-vertices have already been removed and every A4 meets A.
    for (int pair : masks_of_size(W_SIZE, 2)) {
        if (!(pair & ACTIVE_MASK)) continue;
        int incidence = 0;
        for (int face : family) incidence += (face & pair) == pair;
        if (incidence > 2) return false;
    }
    for (int triple : masks_of_size(W_SIZE, 3)) {
        if (!(triple & ACTIVE_MASK)) continue;
        int incidence = 0;
        for (int face : family) incidence += (face & triple) == triple;
        if (incidence > 1) return false;
    }
    return true;
}

map<vector<int>, vector<int>> generate_a4_orbits(int row_count) {
    map<vector<int>, vector<int>> orbits;
    vector<int> active_parts(row_count);

    auto active_recurse = [&](auto&& active_self, int row) -> void {
        if (row < row_count) {
            for (int active = 1; active <= ACTIVE_MASK; ++active) {
                active_parts[row] = active;
                active_self(active_self, row + 1);
            }
            return;
        }

        vector<int> needs(row_count);
        for (int index = 0; index < row_count; ++index) {
            needs[index] = 4 - popcount(active_parts[index]);
        }
        vector<int> column_counts(1 << row_count);
        auto column_recurse = [&](auto&& column_self,
                                  int pattern,
                                  int remaining,
                                  vector<int>& remaining_needs) -> void {
            if (pattern == (1 << row_count)) {
                if (remaining != 0) return;
                if (any_of(
                        remaining_needs.begin(), remaining_needs.end(),
                        [](int value) { return value != 0; }
                    )) {
                    return;
                }
                vector<int> family = active_parts;
                int vertex = 3;
                for (int column = 0; column < (1 << row_count); ++column) {
                    for (int copy = 0; copy < column_counts[column]; ++copy) {
                        for (int family_row = 0;
                             family_row < row_count;
                             ++family_row) {
                            if (column & (1 << family_row)) {
                                family[family_row] |= 1 << vertex;
                            }
                        }
                        ++vertex;
                    }
                }
                if (vertex != W_SIZE) {
                    throw runtime_error("A4 column compression lost isolates");
                }
                set<int> distinct(family.begin(), family.end());
                if (distinct.size() != family.size()) return;
                if (!locally_admissible_a4_family(family)) return;
                vector<int> key = canonical_a4_key(family);
                orbits.emplace(std::move(key), family);
                return;
            }

            int maximum = remaining;
            for (int family_row = 0; family_row < row_count; ++family_row) {
                if (pattern & (1 << family_row)) {
                    maximum = min(maximum, remaining_needs[family_row]);
                }
            }
            for (int count = 0; count <= maximum; ++count) {
                bool valid = true;
                for (int family_row = 0; family_row < row_count; ++family_row) {
                    if (pattern & (1 << family_row)) {
                        remaining_needs[family_row] -= count;
                        valid &= remaining_needs[family_row] >= 0;
                    }
                }
                column_counts[pattern] = count;
                if (valid) {
                    column_self(
                        column_self,
                        pattern + 1,
                        remaining - count,
                        remaining_needs
                    );
                }
                for (int family_row = 0; family_row < row_count; ++family_row) {
                    if (pattern & (1 << family_row)) {
                        remaining_needs[family_row] += count;
                    }
                }
            }
            column_counts[pattern] = 0;
        };
        column_recurse(column_recurse, 0, 9, needs);
    };
    active_recurse(active_recurse, 0);
    return orbits;
}

struct LaterCase {
    int missing_four_count = 0;
    vector<int> a4_family;
    vector<int> forced_fives;
    int chosen_a5_count = 0;
    array<int, W_SIZE> chosen_a5_capacity{};
    vector<int> active_fives;
    vector<int> isolate_fives;
    vector<int> triples;
    vector<int> pairs;
    vector<int> isolate_pairs;
    array<Bits, 1 << W_SIZE> hit_triples{};
    array<Bits, 1 << W_SIZE> hit_pairs{};
    Bits all_triples{};
    Bits all_pairs{};
    map<int, uint64_t> facet_histogram;
    uint64_t valid_family_count = 0;
    uint64_t target_family_count = 0;
};

LaterCase make_later_case(int missing_four_count, const vector<int>& family) {
    LaterCase data;
    data.missing_four_count = missing_four_count;
    data.a4_family = family;
    set<int> forced;
    for (int face : family) {
        for (int vertex = 0; vertex < W_SIZE; ++vertex) {
            if (!(face & (1 << vertex))) forced.insert(face | (1 << vertex));
        }
    }
    data.forced_fives.assign(forced.begin(), forced.end());
    data.chosen_a5_count = 7 + 6 * missing_four_count
        - static_cast<int>(data.forced_fives.size());

    for (int vertex = 0; vertex < W_SIZE; ++vertex) {
        int a4_incidence = 0;
        int forced_incidence = 0;
        for (int face : family) a4_incidence += (face >> vertex) & 1;
        for (int face : data.forced_fives) {
            forced_incidence += (face >> vertex) & 1;
        }
        const int baseline = vertex < 3 ? 5 : 6;
        data.chosen_a5_capacity[vertex] =
            baseline + 6 * a4_incidence - forced_incidence - 1;
    }

    data.isolate_fives = masks_of_size(9, 5, 3);
    data.triples = masks_of_size(W_SIZE, 3);
    data.pairs = masks_of_size(W_SIZE, 2);
    data.isolate_pairs = masks_of_size(9, 2, 3);
    data.all_triples = all_bits(static_cast<int>(data.triples.size()));
    data.all_pairs = all_bits(static_cast<int>(data.pairs.size()));
    for (int mask = 0; mask <= FULL_W_MASK; ++mask) {
        for (int index = 0; index < static_cast<int>(data.triples.size()); ++index) {
            if (mask & data.triples[index]) {
                data.hit_triples[mask][index / 64]
                    |= uint64_t{1} << (index % 64);
            }
        }
        for (int index = 0; index < static_cast<int>(data.pairs.size()); ++index) {
            if (mask & data.pairs[index]) {
                data.hit_pairs[mask][index / 64]
                    |= uint64_t{1} << (index % 64);
            }
        }
        if (popcount(mask) != 5 || !(mask & ACTIVE_MASK)) continue;
        bool contains_a4 = false;
        for (int face : family) contains_a4 |= (mask & face) == face;
        if (!contains_a4) data.active_fives.push_back(mask);
    }
    return data;
}

void evaluate_later_family(
    LaterCase& data,
    const vector<int>& isolate_family
) {
    const int remaining = data.chosen_a5_count
        - static_cast<int>(isolate_family.size());
    if (remaining < 0 || remaining > 1) {
        throw runtime_error("a=2,3 census expects zero or one active A5");
    }
    Bits triples = data.all_triples;
    Bits pairs = data.all_pairs;
    for (int face : data.a4_family) {
        triples = bit_and(triples, data.hit_triples[face]);
        pairs = bit_and(pairs, data.hit_pairs[face]);
    }
    for (int face : isolate_family) {
        triples = bit_and(triples, data.hit_triples[face]);
        pairs = bit_and(pairs, data.hit_pairs[face]);
    }

    int tau2 = 0;
    for (int pair : data.isolate_pairs) {
        bool hits = true;
        for (int face : isolate_family) hits &= bool(pair & face);
        tau2 += hits;
    }
    array<int, W_SIZE> used{};
    for (int face : isolate_family) {
        for (int vertex = 0; vertex < W_SIZE; ++vertex) {
            used[vertex] += (face >> vertex) & 1;
        }
    }
    auto record = [&](const Bits& final_triples, const Bits& final_pairs) {
        if (bit_count(final_pairs) != 0) return;
        const int facet_count = bit_count(final_triples) + tau2;
        ++data.facet_histogram[facet_count];
        ++data.valid_family_count;
        const int target = 98 - 2 * data.missing_four_count;
        data.target_family_count += facet_count == target;
    };
    if (remaining == 0) {
        record(triples, pairs);
        return;
    }
    for (int face : data.active_fives) {
        bool capacity_ok = true;
        for (int vertex = 0; vertex < W_SIZE; ++vertex) {
            if (used[vertex] + ((face >> vertex) & 1)
                > data.chosen_a5_capacity[vertex]) {
                capacity_ok = false;
                break;
            }
        }
        if (!capacity_ok) continue;
        record(
            bit_and(triples, data.hit_triples[face]),
            bit_and(pairs, data.hit_pairs[face])
        );
    }
}

void enumerate_later_a5(LaterCase& data) {
    if (data.chosen_a5_count < 3) return;
    for (int isolate_count = 3;
         isolate_count <= data.chosen_a5_count;
         ++isolate_count) {
        vector<int> selected;
        array<int, W_SIZE> degrees{};
        auto recurse = [&](auto&& self, int start, int left, int common) -> void {
            if (left == 0) {
                if (common == 0) evaluate_later_family(data, selected);
                return;
            }
            for (int index = start;
                 index <= static_cast<int>(data.isolate_fives.size()) - left;
                 ++index) {
                const int face = data.isolate_fives[index];
                bool capacity_ok = true;
                for (int vertex = 3; vertex < W_SIZE; ++vertex) {
                    if ((face & (1 << vertex))
                        && degrees[vertex] >= data.chosen_a5_capacity[vertex]) {
                        capacity_ok = false;
                        break;
                    }
                }
                if (!capacity_ok) continue;
                for (int vertex = 3; vertex < W_SIZE; ++vertex) {
                    degrees[vertex] += (face >> vertex) & 1;
                }
                selected.push_back(face);
                self(
                    self,
                    index + 1,
                    left - 1,
                    selected.size() == 1 ? face : (common & face)
                );
                selected.pop_back();
                for (int vertex = 3; vertex < W_SIZE; ++vertex) {
                    degrees[vertex] -= (face >> vertex) & 1;
                }
            }
        };
        recurse(recurse, 0, isolate_count, 0);
    }
}

void run_later_branch(int missing_four_count) {
    map<vector<int>, vector<int>> orbits = generate_a4_orbits(missing_four_count);
    map<int, int> overlap_histogram;
    map<int, uint64_t> aggregate_histogram;
    uint64_t aggregate_valid = 0;
    uint64_t aggregate_target = 0;

    for (const auto& [key, family] : orbits) {
        (void)key;
        int shared_triple_pairs = 0;
        for (int first = 0; first < missing_four_count; ++first) {
            for (int second = 0; second < first; ++second) {
                shared_triple_pairs += popcount(family[first] & family[second]) == 3;
            }
        }
        ++overlap_histogram[shared_triple_pairs];
        LaterCase data = make_later_case(missing_four_count, family);
        enumerate_later_a5(data);
        aggregate_valid += data.valid_family_count;
        aggregate_target += data.target_family_count;
        for (const auto& [facet_count, multiplicity] : data.facet_histogram) {
            aggregate_histogram[facet_count] += multiplicity;
        }
    }

    if (missing_four_count == 2) {
        const map<int, int> expected_overlaps{{0, 16}, {1, 1}};
        const map<int, uint64_t> expected_histogram{
            {63, 120}, {65, 180}, {67, 108}, {68, 540}, {69, 1260},
            {70, 660}, {71, 1739}, {72, 1610}, {73, 1950}, {74, 1660},
            {75, 3520}, {77, 2053}, {78, 450}, {79, 15}, {81, 178},
            {82, 10}, {83, 115}, {84, 555}, {86, 45}, {87, 10}
        };
        if (orbits.size() != 17 || overlap_histogram != expected_overlaps
            || aggregate_valid != 16778
            || aggregate_histogram != expected_histogram) {
            throw runtime_error("a=2 representative-state census changed");
        }
    } else {
        const map<int, int> expected_overlaps{{0, 132}, {1, 11}, {3, 1}};
        const map<int, uint64_t> expected_histogram{
            {73, 270}, {74, 360}, {75, 540},
            {78, 180}, {84, 540}, {86, 45}
        };
        if (orbits.size() != 144 || overlap_histogram != expected_overlaps
            || aggregate_valid != 1935
            || aggregate_histogram != expected_histogram) {
            throw runtime_error("a=3 representative-state census changed");
        }
    }

    cout << "a=" << missing_four_count
         << " a4_orbits=" << orbits.size()
         << " shared_isolate_triple_pair_histogram=";
    for (const auto& [pairs, multiplicity] : overlap_histogram) {
        cout << pairs << ':' << multiplicity << ',';
    }
    cout << '\n';
    cout << "a=" << missing_four_count
         << " valid_families=" << aggregate_valid
         << " target_" << 98 - 2 * missing_four_count
         << '=' << aggregate_target
         << " histogram=";
    for (const auto& [facet_count, multiplicity] : aggregate_histogram) {
        cout << facet_count << ':' << multiplicity << ',';
    }
    cout << '\n';
    if (aggregate_target != 0) {
        throw runtime_error("positive-a target facet count survived");
    }
}

}  // namespace

int main() {
    run_a1();
    run_later_branch(2);
    run_later_branch(3);
}
