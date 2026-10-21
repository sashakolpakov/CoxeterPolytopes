# Direct rank-14 simplex payment test

This package records an unresolved rank-14 partial payment test with every finite vertex budget present and no discriminant partition. HiGHS reached its 300-second limit without returning a primal objective or an exact certificate. The optimization was not extended or restarted; no feasibility or obstruction result is claimed.

The [source](../scripts/fresh_sqrt5_rank14_core_spectator_lp_20261017.py) reconstructs the frozen [connected-core and finite-spectator domain](fresh_sqrt5_rank14_core_spectator_domain_20261017.md). It includes all 178,723 simplex aggregates from that domain and all 11,751 finite rank-14 vertex forests in the Q(sqrt(5)) catalogue. Every vertex budget explicitly counts all 91 unordered root pairs.

The [saved result](fresh_sqrt5_rank14_core_spectator_lp_20261017.json) preserves all exact rows, budgets, variable mappings, and the terminal solver status. Its certificate list is empty. A fresh saved reconstruction replay passed; it checks the exact domain and does not repeat the optimizer or assert that the timed-out model is feasible.

There are 405,933 whole-forest corner-orbit variables. Exactly 152,193 occur in at least one face inequality. The other 253,740 are fixed to zero during optimization and retained in the complete certificate domain. This does not weaken the optimization: their coefficients occur only in nonnegative upper-budget sums, so decreasing them to zero preserves every lower face inequality and can only improve a vertex budget.

Each maximal simplex with `m` vertices has one demand `2 * binomial(m,3)`. All its literal vertex-corner multiplicities are retained. No individual triangle requirement is imposed in place of the justified simplex aggregate.

The sparse optimizer has a 300-second limit and is invoked once. Exact primal recovery checks every aggregate and every complete vertex budget. Exact dual recovery uses the saved face multipliers and rational column-load bounds, including all columns omitted from the numerical solve. A saved replay reconstructs the entire domain and verifies the certificates without rerunning optimization.

The intended threshold is 14: in dimension 14, `f2/f0 > 21`, and `2T+Q >= 5 f2 - 91 f0 > 14 f0`. Nevertheless, complete vertex coverage does not supply complete face coverage. The connected order-15 cores and other new connected cores are absent, as are quadrilateral and general polygonal face obligations. In particular, `A14`, `B14`, `D14`, and `A1+B13` have complete budgets but no face occurrences in this inherited-core domain.

This is a Q(sqrt(5)) calculation. It neither addresses every ground field nor proves that an exclusion in dimension 14 would automatically exclude larger dimensions. No dimension bound follows from a partial feasible rule.
