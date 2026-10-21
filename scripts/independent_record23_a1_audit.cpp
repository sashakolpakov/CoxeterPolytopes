// Hostile independent audit of the H9/14 record-23 a=1 branch.
// Archived audit checker: it shares no implementation with the primary filter.

#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

static constexpr int NR = 6;
static constexpr int NC = 12;
static constexpr int N = NR + NC;
static constexpr int FULL = (1 << NC) - 1;

static int pc(int x) { return __builtin_popcount((unsigned)x); }

static int binom(int n, int k) {
    if (k < 0 || k > n) return 0;
    int value = 1;
    for (int i = 1; i <= k; ++i) value = value * (n - k + i) / i;
    return value;
}

static int macaulay_raise(int value, int degree) {
    int remainder = value;
    int answer = 0;
    int ceiling = 100;
    for (int d = degree; d >= 1; --d) {
        int top = d - 1;
        while (top + 1 < ceiling && binom(top + 1, d) <= remainder) ++top;
        remainder -= binom(top, d);
        answer += binom(top + 1, d + 1);
        ceiling = top;
    }
    if (remainder != 0) throw runtime_error("Macaulay expansion failed");
    return answer;
}

struct Totals {
    uint64_t generated = 0;
    uint64_t structural = 0;
    uint64_t left_choices = 0;
    uint64_t target = 0;
    int maximum = -1;
    int best_t3 = -1;
    int best_t2 = -1;
    int best_A = 0;
    int best_L = 0;
    array<int, 5> best_B{};
    map<int, uint64_t> histogram;
};

static void consume(const string& line, Totals& totals) {
    if (line.empty() || line.rfind(">>", 0) == 0) return;
    if ((unsigned char)line[0] - 63 != N) throw runtime_error("wrong order");
    int bits = N * (N - 1) / 2;
    int payload = (bits + 5) / 6;
    if ((int)line.size() != payload + 1) throw runtime_error("bad graph6 length");
    array<int, NR> row_complements{};
    array<int, NC> column_complements{};
    int index = 0;
    for (int upper = 1; upper < N; ++upper) {
        for (int lower = 0; lower < upper; ++lower) {
            int code = (unsigned char)line[1 + index / 6] - 63;
            int edge = (code >> (5 - index % 6)) & 1;
            ++index;
            if (!edge) continue;
            if (lower >= NR || upper < NR) throw runtime_error("not bipartite");
            int row = lower;
            int column = upper - NR;
            row_complements[row] |= 1 << column;
            column_complements[column] |= 1 << row;
        }
    }
    ++totals.generated;
    int a_row = -1;
    vector<int> b_rows;
    for (int row = 0; row < NR; ++row) {
        int degree = pc(row_complements[row]);
        if (degree == 8) {
            if (a_row != -1) throw runtime_error("two degree-eight rows");
            a_row = row;
        } else if (degree == 7) {
            b_rows.push_back(row);
        } else {
            throw runtime_error("wrong row degree");
        }
    }
    if (a_row < 0 || b_rows.size() != 5) throw runtime_error("wrong row profile");
    for (int i = 0; i < NC; ++i) for (int j = i + 1; j < NC; ++j)
        if (!(column_complements[i] & column_complements[j]))
            throw runtime_error("-Y1 failure");

    int A = FULL ^ row_complements[a_row];
    array<int, 5> B{};
    set<int> distinct_B;
    bool structural = true;
    for (int i = 0; i < 5; ++i) {
        B[i] = FULL ^ row_complements[b_rows[i]];
        structural &= (B[i] & A) != A;  // B_i must not contain A_4.
        distinct_B.insert(B[i]);
    }
    structural &= distinct_B.size() == 5;
    if (!structural) return;
    ++totals.structural;

    int t3 = 0;
    for (int i = 0; i < NC; ++i) for (int j = i + 1; j < NC; ++j)
        for (int k = j + 1; k < NC; ++k)
            t3 += !(column_complements[i] & column_complements[j]
                    & column_complements[k]);

    array<int, 1 << NC> shadow_count{};
    for (int extra = 0; extra < NC; ++extra) {
        if (A & (1 << extra)) continue;
        int shadow = A | (1 << extra);
        for (int size = 1; size <= 3; ++size) {
            for (int subset = shadow; subset; subset = (subset - 1) & shadow) {
                if (pc(subset) == size) ++shadow_count[subset];
            }
        }
    }

    for (int x = 0; x < NC; ++x) for (int y = x + 1; y < NC; ++y)
        for (int z = y + 1; z < NC; ++z) {
            int L = (1 << x) | (1 << y) | (1 << z);
            if (!(A & L)) continue;
            int I = FULL ^ L;
            vector<int> f0;
            for (int face : B) if (!(face & L)) f0.push_back(face);
            if (f0.empty()) continue;
            int common = I;
            for (int face : f0) common &= face;
            if (common) continue;

            bool local = true;
            for (int size = 1; size <= 3 && local; ++size) {
                for (int S = 1; S <= FULL; ++S) {
                    if (pc(S) != size) continue;
                    bool with_left = bool(S & L);
                    int base_prev = 0, base_high = 0;
                    if (size == 1) {
                        base_prev = with_left ? 4 : 5;
                        base_high = with_left ? 5 : 6;
                    } else if (size == 2) {
                        base_prev = with_left ? 3 : 4;
                        base_high = with_left ? 4 : 5;
                    } else {
                        base_prev = with_left ? 2 : 4;
                        base_high = with_left ? 3 : 4;
                    }
                    int r = (A & S) == S;
                    int q = shadow_count[S];
                    int d = 0;
                    for (int face : B) d += (face & S) == S;
                    int previous = base_prev - r;
                    int high = base_high + 6 * r - q - d;
                    int degree = 4 - size;
                    if (previous < 1 || high < 1
                        || high > macaulay_raise(previous, degree)) {
                        local = false;
                        break;
                    }
                }
            }
            if (!local) continue;

            int t2 = 0;
            for (int i = 0; i < NC; ++i) {
                if (L & (1 << i)) continue;
                for (int j = i + 1; j < NC; ++j) {
                    if (L & (1 << j)) continue;
                    bool hits = true;
                    for (int face : f0) hits &= bool(face & ((1 << i) | (1 << j)));
                    t2 += hits;
                }
            }
            int facets = t3 + t2;
            ++totals.left_choices;
            ++totals.histogram[facets];
            totals.target += facets == 96;
            if (facets > totals.maximum) {
                totals.maximum = facets;
                totals.best_t3 = t3;
                totals.best_t2 = t2;
                totals.best_A = A;
                totals.best_L = L;
                totals.best_B = B;
            }
        }
}

int main() {
    try {
        Totals totals;
        string line;
        while (getline(cin, line)) consume(line, totals);
        cout << "generated " << totals.generated << '\n';
        cout << "structural " << totals.structural << '\n';
        cout << "left_choices " << totals.left_choices << '\n';
        cout << "target " << totals.target << '\n';
        cout << "maximum " << totals.maximum << '\n';
        cout << "best_t3 " << totals.best_t3 << '\n';
        cout << "best_t2 " << totals.best_t2 << '\n';
        cout << "best_A " << totals.best_A << '\n';
        cout << "best_L " << totals.best_L << '\n';
        cout << "best_B";
        for (int face : totals.best_B) cout << ' ' << face;
        cout << '\n';
        for (auto [value, count] : totals.histogram)
            cout << "histogram " << value << ' ' << count << '\n';
    } catch (const exception& error) {
        cerr << error.what() << '\n';
        return 1;
    }
}
