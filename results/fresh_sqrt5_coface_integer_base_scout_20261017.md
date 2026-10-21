# Coface-summed integer prism-base packing: exact fixed-witness diagnostic

The 7,504 unlabelled integer packing constraints below have **no violations** on the saved rank-13 restricted dual. All face and vertex masses remain unchanged. This family therefore leaves its exact obstruction

\[
\frac{66908102699388132000000000000000}
{5419878674792788039126804769429}
\approx12.3449447329>12
\]

unchanged. This is a statement about a restricted weighting model, not a dimension theorem or a globally realizable polytope.

The standalone [source](../scripts/fresh_sqrt5_coface_integer_base_scout_20261017.py) and [saved exact data](fresh_sqrt5_coface_integer_base_scout_20261017.json) use only standard-library integer/Fraction arithmetic. Generation and a fresh `--verify` reconstruction both passed. No optimizer was invoked.

## Valid inequality and scope

Fix a maximal simplex face \(S\) of dimension \(d\) and an immediate \((d+1)\)-coface \(U\). The [shared-base proof](fresh_prism_simplex_shared_capacity_20261017.md) shows that bases of the selected maximal dotted prisms in this \((S,U)\) have pairwise disjoint vertex subsets \(B\subseteq V(S)\). A prism of dimension \(k\) has a simplex base with \(k\) vertices. Consequently, for each integer \(t=2,\ldots,d+1\),

\[
\sum_B\left\lfloor\frac{|B|}{t}\right\rfloor
\leq\left\lfloor\frac{d+1}{t}\right\rfloor.
\]

There are \(13-d\) immediate cofaces of \(S\) in a simple 13-polytope. Summing the last inequality gives the zero-vertex-cost signed column

\[
\sum_Q\sum_{B\subset S\text{ a base of }Q}
\left\lfloor\frac{k_Q}{t}\right\rfloor\alpha_Q
- (13-d)\left\lfloor\frac{d+1}{t}\right\rfloor\alpha_S\leq0.
\]

The positive coefficient counts each actual maximal-prism base once, rather than each of its vertices or quadrilateral subfaces. For a square, the invariant convention counts all four edge bases. For a higher-dimensional prism, it counts its two simplex bases. The existing [weighted-capacity theorem](fresh_vertex_weighted_prism_base_capacity_20261017.md) supplies the compatible incidence framework. These statements concern the selected maximal dotted-prism aggregates; no assertion is made about arbitrary overlapping all-finite product faces.

## Literal reconstruction and result

The diagnostic binds the frozen [flag maps](fresh_sqrt5_flag_refined_base_capacity_20261017.json), [actual-provenance data](fresh_sqrt5_actual_provenance_dual_20261017.json), and [full model](fresh_sqrt5_flag_local_packing_full_model_20261017.json). It uses precisely the unchanged face masses of the [five-port diagnostic](fresh_sqrt5_joint_five_port_diagnostic_20261017.md): 2,219 positive faces, including 355 actual quad provenances.

Flag occurrences are grouped by `(provenance, deleted base normal)`. Within each group the base type and coface root agree; its canonical simplex vertices and original omitted pairs are distinct; their count equals the actual prism dimension. The groups recover all 1,858 saved actual provenances, including zero-mass ones. The four-base/two-base count is checked separately for every provenance. Every base type maps consistently to one simplex row, and all 4,969 such simplex rows are distinct. No merged-provenance representative is chosen.

| Simplex dimension \(d\) | Divisor \(t\) | Columns | Positive loads |
|---:|---:|---:|---:|
| 1 | 2 | 3,460 | 380 |
| 2 | 2 | 780 | 96 |
| 2 | 3 | 780 | 81 |
| 3 | 2 | 432 | 138 |
| 3 | 3 | 432 | 117 |
| 3 | 4 | 432 | 117 |
| 4 | 2 | 297 | 189 |
| 4 | 3 | 297 | 189 |
| 4 | 4 | 297 | 181 |
| 4 | 5 | 297 | 116 |

There are 6,016 literal bases, 924 with positive quad mass, and 1,604 positive-load columns. Every exact slack is nonnegative; the maximum capacity utilization is exactly one. The output saves every base grouping, every provenance coefficient and all 7,504 loads/capacities/slacks. Its direct SHA-256 manifest binds the input sources/data and the two proof notes.

No simplex or vertex mass changes. Therefore the previously checked face constraints, signed transfers, complete local-simplex packing covers and three-/four-/five-port repairs remain valid for this same witness. No new cover or independent normalizer combination is needed.

## Remaining direction

This test sums over immediate cofaces and uses only base cardinalities. It does **not** test labelled subset constraints within a coface or the full packing cone of actual subsets of \(V(S)\). A useful next refinement must retain simultaneous `(coface root q, vertex subset T)` orbits, with every canonical base vertex included. Independent marginal orbits for \(q\) and \(T\) need not give a valid common incidence count. Those labelled constraints remain untested here; no further repair or global optimization was run.
