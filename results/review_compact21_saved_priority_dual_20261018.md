# Root review of the exact two-row saved-dual repair

18 October 2026. **PASS for the exact diagnostic and its stated scope.**
The [report](fresh_compact21_extreme_e8_priority_saved_dual_20261018.md)
and full streaming checker were read after the separate geometric review
of the new active/inactive charge rule.

The old dual multipliers are nonpositive, and the checker reconstructs
every touched column with exact fractions. For target forest mass mu,
splitting it into active mass lambda and inactive mass mu-lambda leaves
every other column unchanged. The extreme bonus column load increases
by2(mu-lambda). Its exact feasibility condition is therefore
lambda>=mu+L/2, together with0<=lambda<=mu. The active row has right-hand
side one larger than before, so the dual lower bound decreases by exactly
lambda. Taking the least allowed mass gives the stated optimal repair
within this fixed-multiplier family, namely7158707133073/179000000000,
strictly below40. Direct changed-row summation provides a separate exact
arithmetic check in the source.

The streaming parser bounds each object and hashes the old model without
loading it. Untouched bonus columns have zero load; every changed column
is touched and checked. The geometric interpretation of the frozen old
model is inherited, and all its old face rows are retained by the new
reviewed rule.

The author records generation and fresh exact replay; this root review
does not assert an additional execution. This is a repair of one saved
dual, not an optimum of the new model. A subsequent numerical support
scout still lies above40. Full successor feasibility and a dimension20
theorem remain open. No Lean verification is claimed.
