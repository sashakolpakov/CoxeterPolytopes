// Exact solver-independent census for the H9/14 record-23 a=0 branch.
//
// Input is the graph6 stream produced by
//
//   genbg -q -Y1 -d7:0 -D7:7 7 12 49:49
//
// The first colour class consists of seven complement rows, one for each
// missing five-face.  Every row has seven 1s because its missing face has
// five of the twelve usable vertices.  The second colour class consists of
// those twelve vertices.  -Y1 says every pair of columns has a common 1:
// equivalently, every vertex pair misses some missing five-face, so there is
// no transversal of order at most two after adjoining {12,13}.
//
// genbg generates one representative of every bicoloured isomorphism class.
// This checker independently validates every property used from its command
// line, removes repeated rows (which would repeat a minimal nonface), expands
// every choice of the three K_{3,2}-left vertices, applies all vertex/edge/
// triple strict-link bounds, and computes the exact blocker facet count.

#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <string>

namespace {

constexpr int kRowCount = 7;
constexpr int kColumnCount = 12;
constexpr int kGraphOrder = kRowCount + kColumnCount;

int population(int value) {
    return __builtin_popcount(static_cast<unsigned>(value));
}

struct Census {
    std::uint64_t generated_graph_count = 0;
    std::uint64_t distinct_row_graph_count = 0;
    std::uint64_t locally_valid_left_choice_count = 0;
    std::uint64_t target_facet_count_occurrences = 0;
    int maximum_facet_count = -1;
    int maximum_three_transversal_count = -1;
    int maximum_two_transversal_count = -1;
    std::array<int, 3> maximizing_left_vertices{};
    std::array<int, kColumnCount> maximizing_columns{};
    std::map<int, std::uint64_t> facet_count_histogram;
};

void audit_graph6_line(const std::string& line, Census& census) {
    if (line.empty()) return;
    if (line.rfind(">>", 0) == 0) return;
    if (static_cast<int>(static_cast<unsigned char>(line[0])) - 63
        != kGraphOrder) {
        throw std::runtime_error("record-23 census received wrong graph order");
    }
    constexpr int edge_bit_count = kGraphOrder * (kGraphOrder - 1) / 2;
    constexpr int payload_length = (edge_bit_count + 5) / 6;
    if (line.size() != static_cast<std::size_t>(payload_length + 1)) {
        throw std::runtime_error("record-23 census received malformed graph6");
    }

    std::array<int, kColumnCount> columns{};
    std::array<int, kRowCount> rows{};
    int bit_index = 0;
    for (int upper = 1; upper < kGraphOrder; ++upper) {
        for (int lower = 0; lower < upper; ++lower) {
            int code = static_cast<int>(
                static_cast<unsigned char>(line[1 + bit_index / 6])
            ) - 63;
            int edge = (code >> (5 - bit_index % 6)) & 1;
            ++bit_index;
            if (!edge) continue;
            if (lower >= kRowCount || upper < kRowCount) {
                throw std::runtime_error(
                    "genbg stream contains an edge inside a colour class"
                );
            }
            int row = lower;
            int column = upper - kRowCount;
            columns[column] |= 1 << row;
            rows[row] |= 1 << column;
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

    ++census.generated_graph_count;
    for (int row : rows) {
        if (population(row) != 7) {
            throw std::runtime_error("genbg stream violates the row-sum gate");
        }
    }
    for (int first = 0; first < kColumnCount; ++first) {
        for (int second = first + 1; second < kColumnCount; ++second) {
            if (!(columns[first] & columns[second])) {
                throw std::runtime_error("genbg stream violates the -Y1 gate");
            }
        }
    }

    std::set<int> distinct_rows(rows.begin(), rows.end());
    if (distinct_rows.size() != kRowCount) return;
    ++census.distinct_row_graph_count;

    // A triple hits all seven missing faces exactly when the corresponding
    // three complement columns have empty intersection.
    int three_transversal_count = 0;
    for (int first = 0; first < kColumnCount; ++first) {
        for (int second = first + 1; second < kColumnCount; ++second) {
            for (int third = second + 1; third < kColumnCount; ++third) {
                three_transversal_count += !(
                    columns[first] & columns[second] & columns[third]
                );
            }
        }
    }

    // Conditions not depending on which three columns are declared to be L.
    // In original incidence, g4/g3/g2 strict positivity gives degree caps
    // 5,4,3 on isolate vertices, isolate pairs, and isolate triples.  In
    // complement incidence these are union-size lower bounds 2,3,4.
    bool universal_local_bounds = true;
    for (int vertex = 0; vertex < kColumnCount; ++vertex) {
        universal_local_bounds &= population(columns[vertex]) >= 2;
    }
    for (int first = 0; first < kColumnCount; ++first) {
        for (int second = first + 1; second < kColumnCount; ++second) {
            universal_local_bounds &= population(
                columns[first] | columns[second]
            ) >= 3;
        }
    }
    for (int first = 0; first < kColumnCount; ++first) {
        for (int second = first + 1; second < kColumnCount; ++second) {
            for (int third = second + 1; third < kColumnCount; ++third) {
                universal_local_bounds &= population(
                    columns[first] | columns[second] | columns[third]
                ) >= 4;
            }
        }
    }
    if (!universal_local_bounds) return;

    // A left vertex has the stronger degree cap 4; a pair or triple involving
    // one has the stronger link-codegree caps 3 and 2.  Compute the columns
    // eligible to receive the left colour.
    std::array<bool, kColumnCount> eligible_left{};
    for (int left = 0; left < kColumnCount; ++left) {
        bool eligible = population(columns[left]) >= 3;
        for (int other = 0; other < kColumnCount && eligible; ++other) {
            if (other == left) continue;
            eligible &= population(columns[left] | columns[other]) >= 4;
        }
        for (int second = 0; second < kColumnCount && eligible; ++second) {
            if (second == left) continue;
            for (int third = second + 1;
                 third < kColumnCount && eligible;
                 ++third) {
                if (third == left) continue;
                eligible &= population(
                    columns[left] | columns[second] | columns[third]
                ) >= 5;
            }
        }
        eligible_left[left] = eligible;
    }

    for (int first_left = 0; first_left < kColumnCount; ++first_left) {
        if (!eligible_left[first_left]) continue;
        for (int second_left = first_left + 1;
             second_left < kColumnCount;
             ++second_left) {
            if (!eligible_left[second_left]) continue;
            for (int third_left = second_left + 1;
                 third_left < kColumnCount;
                 ++third_left) {
                if (!eligible_left[third_left]) continue;
                std::array<bool, kColumnCount> is_left{};
                is_left[first_left] = true;
                is_left[second_left] = true;
                is_left[third_left] = true;

                // Rows in this intersection are exactly the A5 faces
                // disjoint from L.  They must be nonempty and have no common
                // isolate, or L (possibly plus one point) is too small a
                // transversal.
                int zero_active_rows = columns[first_left]
                    & columns[second_left] & columns[third_left];
                if (!zero_active_rows) continue;
                bool no_common_isolate = true;
                for (int vertex = 0; vertex < kColumnCount; ++vertex) {
                    if (is_left[vertex]) continue;
                    no_common_isolate &= bool(
                        columns[vertex] & zero_active_rows
                    );
                }
                if (!no_common_isolate) continue;

                int two_transversal_count = 0;
                for (int first = 0; first < kColumnCount; ++first) {
                    if (is_left[first]) continue;
                    for (int second = first + 1;
                         second < kColumnCount;
                         ++second) {
                        if (is_left[second]) continue;
                        two_transversal_count += !(
                            columns[first] & columns[second]
                            & zero_active_rows
                        );
                    }
                }
                int facet_count = three_transversal_count
                    + two_transversal_count;
                ++census.locally_valid_left_choice_count;
                ++census.facet_count_histogram[facet_count];
                census.target_facet_count_occurrences += facet_count == 98;
                if (facet_count > census.maximum_facet_count) {
                    census.maximum_facet_count = facet_count;
                    census.maximum_three_transversal_count =
                        three_transversal_count;
                    census.maximum_two_transversal_count =
                        two_transversal_count;
                    census.maximizing_left_vertices = {
                        first_left, second_left, third_left
                    };
                    census.maximizing_columns = columns;
                }
            }
        }
    }
}

}  // namespace

int main() {
    try {
        Census census;
        std::string line;
        while (std::getline(std::cin, line)) audit_graph6_line(line, census);
        std::cout << "generated_graph_count "
                  << census.generated_graph_count << '\n';
        std::cout << "distinct_row_graph_count "
                  << census.distinct_row_graph_count << '\n';
        std::cout << "locally_valid_left_choice_count "
                  << census.locally_valid_left_choice_count << '\n';
        std::cout << "target_facet_count_occurrences "
                  << census.target_facet_count_occurrences << '\n';
        std::cout << "maximum_facet_count "
                  << census.maximum_facet_count << '\n';
        std::cout << "maximum_three_transversal_count "
                  << census.maximum_three_transversal_count << '\n';
        std::cout << "maximum_two_transversal_count "
                  << census.maximum_two_transversal_count << '\n';
        std::cout << "maximizing_left_vertices";
        for (int value : census.maximizing_left_vertices) {
            std::cout << ' ' << value;
        }
        std::cout << '\n';
        std::cout << "maximizing_complement_columns";
        for (int value : census.maximizing_columns) {
            std::cout << ' ' << value;
        }
        std::cout << '\n';
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
