// Exact complement-incidence census for the positive-a branches of H9/14
// record 23.  This is a finite integer enumerator, not a SAT/SMT model.
//
// Usage:
//   genbg <branch-specific options> | certify_h9_14_record23_positive A B T
// where A is the number of missing four-faces, B the number of additional
// minimal missing five-faces, and T the scheduled facet count.

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr int kColumnCount = 12;
constexpr int kFullColumns = (1 << kColumnCount) - 1;
constexpr int kLeftChoiceCount = 220;
constexpr int kLeftChoiceWordCount = 4;

int population(int value) {
    return __builtin_popcount(static_cast<unsigned>(value));
}

int choose(int n, int k) {
    if (k < 0 || k > n) return 0;
    int result = 1;
    for (int index = 1; index <= k; ++index) {
        result = result * (n - k + index) / index;
    }
    return result;
}

int macaulay_raise(int value, int degree) {
    if (value < 0 || degree <= 0) {
        throw std::runtime_error("invalid Macaulay input");
    }
    int remaining = value;
    int result = 0;
    int upper = 100;
    for (int order = degree; order >= 1; --order) {
        int top = order - 1;
        while (top + 1 < upper && choose(top + 1, order) <= remaining) {
            ++top;
        }
        if (top >= order) {
            remaining -= choose(top, order);
            result += choose(top + 1, order + 1);
            upper = top;
        } else {
            upper = order - 1;
        }
    }
    if (remaining != 0) {
        throw std::runtime_error("incomplete Macaulay expansion");
    }
    return result;
}

struct Census {
    std::uint64_t generated_graph_count = 0;
    std::uint64_t distinct_row_graph_count = 0;
    std::uint64_t schedule_compatible_graph_count = 0;
    std::uint64_t locally_valid_left_choice_count = 0;
    std::uint64_t blocker_valid_left_choice_count = 0;
    std::uint64_t target_facet_count_occurrences = 0;
    int maximum_facet_count = -1;
    int maximum_three_transversal_count = -1;
    int maximum_two_transversal_count = -1;
    std::map<int, std::uint64_t> facet_count_histogram;
};

struct ParsedGraph {
    std::vector<int> complement_rows;
    std::array<int, kColumnCount> complement_columns{};
};

ParsedGraph parse_graph6(const std::string& line, int row_count) {
    int graph_order = row_count + kColumnCount;
    if (static_cast<int>(static_cast<unsigned char>(line[0])) - 63
        != graph_order) {
        throw std::runtime_error("positive-a census received wrong graph order");
    }
    int edge_bit_count = graph_order * (graph_order - 1) / 2;
    int payload_length = (edge_bit_count + 5) / 6;
    if (line.size() != static_cast<std::size_t>(payload_length + 1)) {
        throw std::runtime_error("positive-a census received malformed graph6");
    }
    ParsedGraph parsed;
    parsed.complement_rows.assign(row_count, 0);
    int bit_index = 0;
    for (int upper = 1; upper < graph_order; ++upper) {
        for (int lower = 0; lower < upper; ++lower) {
            int code = static_cast<int>(
                static_cast<unsigned char>(line[1 + bit_index / 6])
            ) - 63;
            int edge = (code >> (5 - bit_index % 6)) & 1;
            ++bit_index;
            if (!edge) continue;
            if (lower >= row_count || upper < row_count) {
                throw std::runtime_error("edge lies inside a genbg colour class");
            }
            int row = lower;
            int column = upper - row_count;
            parsed.complement_rows[row] |= 1 << column;
            parsed.complement_columns[column] |= 1 << row;
        }
    }
    for (int bit = edge_bit_count; bit < payload_length * 6; ++bit) {
        int code = static_cast<int>(
            static_cast<unsigned char>(line[1 + bit / 6])
        ) - 63;
        if ((code >> (5 - bit % 6)) & 1) {
            throw std::runtime_error("nonzero graph6 padding bit");
        }
    }
    return parsed;
}

std::vector<int> five_shadow(const std::vector<int>& missing_fours) {
    std::set<int> shadow;
    for (int face : missing_fours) {
        for (int vertex = 0; vertex < kColumnCount; ++vertex) {
            if (!(face & (1 << vertex))) shadow.insert(face | (1 << vertex));
        }
    }
    return std::vector<int>(shadow.begin(), shadow.end());
}

struct LeftChoiceTables {
    std::array<int, kLeftChoiceCount> masks{};
    std::array<
        std::array<std::uint64_t, kLeftChoiceWordCount>,
        1 << kColumnCount
    > intersects{};
};

LeftChoiceTables make_left_choice_tables() {
    LeftChoiceTables tables;
    int index = 0;
    for (int first = 0; first < kColumnCount; ++first) {
        for (int second = first + 1; second < kColumnCount; ++second) {
            for (int third = second + 1; third < kColumnCount; ++third) {
                tables.masks[index++] = (1 << first) | (1 << second)
                    | (1 << third);
            }
        }
    }
    for (int subset = 0; subset <= kFullColumns; ++subset) {
        for (int choice = 0; choice < kLeftChoiceCount; ++choice) {
            if (subset & tables.masks[choice]) {
                tables.intersects[subset][choice / 64]
                    |= std::uint64_t{1} << (choice % 64);
            }
        }
    }
    return tables;
}

const LeftChoiceTables kLeftChoices = make_left_choice_tables();

std::array<std::uint64_t, kLeftChoiceWordCount> local_link_allowed_choices(
    const std::vector<int>& missing_fours,
    const std::vector<int>& shadow,
    const std::vector<int>& missing_fives
) {
    std::array<std::uint64_t, kLeftChoiceWordCount> allowed{};
    allowed.fill(~std::uint64_t{0});
    allowed.back() &= (std::uint64_t{1} << (kLeftChoiceCount - 192)) - 1;
    for (int face : missing_fours) {
        for (int word = 0; word < kLeftChoiceWordCount; ++word) {
            allowed[word] &= kLeftChoices.intersects[face][word];
        }
    }
    for (int size = 1; size <= 3; ++size) {
        for (int subset = 1; subset <= kFullColumns; ++subset) {
            if (population(subset) != size) continue;
            int macaulay_degree;
            if (size == 1) {
                macaulay_degree = 3;
            } else if (size == 2) {
                macaulay_degree = 2;
            } else {
                macaulay_degree = 1;
            }
            int four_incidence = std::count_if(
                missing_fours.begin(), missing_fours.end(),
                [subset](int face) { return (face & subset) == subset; }
            );
            int shadow_incidence = std::count_if(
                shadow.begin(), shadow.end(),
                [subset](int face) { return (face & subset) == subset; }
            );
            int five_incidence = std::count_if(
                missing_fives.begin(), missing_fives.end(),
                [subset](int face) { return (face & subset) == subset; }
            );
            auto valid_for_type = [&](bool has_left) {
                int base_previous;
                int base_high;
                if (size == 1) {
                    base_previous = has_left ? 4 : 5;
                    base_high = has_left ? 5 : 6;
                } else if (size == 2) {
                    base_previous = has_left ? 3 : 4;
                    base_high = has_left ? 4 : 5;
                } else {
                    base_previous = has_left ? 2 : 4;
                    base_high = has_left ? 3 : 4;
                }
                int previous = base_previous - four_incidence;
                int high = base_high + 6 * four_incidence
                    - shadow_incidence - five_incidence;
                return previous >= 1 && high >= 1
                    && high <= macaulay_raise(previous, macaulay_degree);
            };
            bool valid_isolate = valid_for_type(false);
            bool valid_left = valid_for_type(true);
            for (int word = 0; word < kLeftChoiceWordCount; ++word) {
                if (!valid_isolate) {
                    allowed[word] &= kLeftChoices.intersects[subset][word];
                }
                if (!valid_left) {
                    allowed[word] &= ~kLeftChoices.intersects[subset][word];
                }
            }
        }
    }
    return allowed;
}

void audit_graph(
    const ParsedGraph& parsed,
    int missing_four_count,
    int missing_five_count,
    int target_facet_count,
    Census& census
) {
    int row_count = missing_four_count + missing_five_count;
    ++census.generated_graph_count;

    std::vector<int> missing_four_rows;
    std::vector<int> missing_five_rows;
    for (int row : parsed.complement_rows) {
        if (population(row) == 8) {
            missing_four_rows.push_back(row);
        } else if (population(row) == 7) {
            missing_five_rows.push_back(row);
        } else {
            throw std::runtime_error("unexpected complement-row degree");
        }
    }
    if (missing_four_rows.size() != static_cast<std::size_t>(missing_four_count)
        || missing_five_rows.size()
            != static_cast<std::size_t>(missing_five_count)) {
        throw std::runtime_error("genbg stream has wrong row-degree distribution");
    }
    for (int first = 0; first < kColumnCount; ++first) {
        for (int second = first + 1; second < kColumnCount; ++second) {
            if (!(parsed.complement_columns[first]
                  & parsed.complement_columns[second])) {
                throw std::runtime_error("genbg stream violates the -Y1 gate");
            }
        }
    }
    if (std::set<int>(missing_four_rows.begin(), missing_four_rows.end()).size()
            != missing_four_rows.size()
        || std::set<int>(missing_five_rows.begin(), missing_five_rows.end()).size()
            != missing_five_rows.size()) {
        return;
    }
    ++census.distinct_row_graph_count;

    std::vector<int> missing_fours;
    std::vector<int> missing_fives;
    for (int row : missing_four_rows) missing_fours.push_back(kFullColumns ^ row);
    for (int row : missing_five_rows) missing_fives.push_back(kFullColumns ^ row);
    std::vector<int> shadow = five_shadow(missing_fours);
    for (int face : missing_fives) {
        if (std::binary_search(shadow.begin(), shadow.end(), face)) return;
    }
    if (shadow.size() + missing_fives.size()
        != static_cast<std::size_t>(7 + 6 * missing_four_count)) {
        return;
    }
    ++census.schedule_compatible_graph_count;

    int three_transversal_count = 0;
    for (int first = 0; first < kColumnCount; ++first) {
        for (int second = first + 1; second < kColumnCount; ++second) {
            for (int third = second + 1; third < kColumnCount; ++third) {
                three_transversal_count += !(
                    parsed.complement_columns[first]
                    & parsed.complement_columns[second]
                    & parsed.complement_columns[third]
                );
            }
        }
    }

    auto locally_allowed = local_link_allowed_choices(
        missing_fours, shadow, missing_fives
    );

    for (int choice = 0; choice < kLeftChoiceCount; ++choice) {
                if (!(locally_allowed[choice / 64]
                      & (std::uint64_t{1} << (choice % 64)))) continue;
                int left_mask = kLeftChoices.masks[choice];
                ++census.locally_valid_left_choice_count;

                int zero_active_rows = 0;
                for (int row = 0; row < row_count; ++row) {
                    if (population(parsed.complement_rows[row]) != 7) continue;
                    int original_face = kFullColumns
                        ^ parsed.complement_rows[row];
                    if (!(original_face & left_mask)) zero_active_rows |= 1 << row;
                }
                if (!zero_active_rows) continue;
                bool no_common_isolate = true;
                for (int vertex = 0; vertex < kColumnCount; ++vertex) {
                    if (left_mask & (1 << vertex)) continue;
                    no_common_isolate &= bool(
                        parsed.complement_columns[vertex] & zero_active_rows
                    );
                }
                if (!no_common_isolate) continue;
                ++census.blocker_valid_left_choice_count;

                int two_transversal_count = 0;
                for (int first = 0; first < kColumnCount; ++first) {
                    if (left_mask & (1 << first)) continue;
                    for (int second = first + 1;
                         second < kColumnCount;
                         ++second) {
                        if (left_mask & (1 << second)) continue;
                        two_transversal_count += !(
                            parsed.complement_columns[first]
                            & parsed.complement_columns[second]
                            & zero_active_rows
                        );
                    }
                }
                int facet_count = three_transversal_count
                    + two_transversal_count;
                ++census.facet_count_histogram[facet_count];
                census.target_facet_count_occurrences +=
                    facet_count == target_facet_count;
                if (facet_count > census.maximum_facet_count) {
                    census.maximum_facet_count = facet_count;
                    census.maximum_three_transversal_count =
                        three_transversal_count;
                    census.maximum_two_transversal_count =
                        two_transversal_count;
                }
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 4) {
            throw std::runtime_error(
                "usage: positive-enumerator A B TARGET_FACETS"
            );
        }
        int missing_four_count = std::stoi(argv[1]);
        int missing_five_count = std::stoi(argv[2]);
        int target_facet_count = std::stoi(argv[3]);
        if (missing_four_count < 1 || missing_four_count > 3
            || missing_five_count < 1 || missing_five_count > 5) {
            throw std::runtime_error("unsupported record-23 branch parameters");
        }
        int row_count = missing_four_count + missing_five_count;
        Census census;
        std::string line;
        while (std::getline(std::cin, line)) {
            if (line.empty() || line.rfind(">>", 0) == 0) continue;
            audit_graph(
                parse_graph6(line, row_count),
                missing_four_count,
                missing_five_count,
                target_facet_count,
                census
            );
        }
        std::cout << "generated_graph_count "
                  << census.generated_graph_count << '\n';
        std::cout << "distinct_row_graph_count "
                  << census.distinct_row_graph_count << '\n';
        std::cout << "schedule_compatible_graph_count "
                  << census.schedule_compatible_graph_count << '\n';
        std::cout << "locally_valid_left_choice_count "
                  << census.locally_valid_left_choice_count << '\n';
        std::cout << "blocker_valid_left_choice_count "
                  << census.blocker_valid_left_choice_count << '\n';
        std::cout << "target_facet_count_occurrences "
                  << census.target_facet_count_occurrences << '\n';
        std::cout << "maximum_facet_count "
                  << census.maximum_facet_count << '\n';
        std::cout << "maximum_three_transversal_count "
                  << census.maximum_three_transversal_count << '\n';
        std::cout << "maximum_two_transversal_count "
                  << census.maximum_two_transversal_count << '\n';
        for (const auto& [facet_count, multiplicity]
             : census.facet_count_histogram) {
            std::cout << "histogram " << facet_count << ' '
                      << multiplicity << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
