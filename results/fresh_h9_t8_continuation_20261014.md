# T-shaped case 8: opposite-root family and remaining state

This phase excludes another continuous metric family for the explicit eighteen-facet compact-H9 candidate. It does not produce a compact H9 polytope or exclude the complete incidence.

The [previous T8 package](fresh_h9_t8_metric_20261014.md) remains unchanged. Its two continuous exclusions were integrated into a fresh bounded necessary-label search. The [first saved continuation](fresh_h9_t8_continuation_20261014_case_08.json) reached an opposite-root assignment in 933 nodes, with both high block labels representing all finite integers at least seven.

The [opposite-root proof](fresh_h9_t8_continuation_20261014_opposite.md) and [exact certificate](fresh_h9_t8_continuation_20261014_opposite.json) now exclude that assignment and a larger family: both own nonright anchor contacts may have any finite labels at least three. All other mandatory contacts remain fixed. The proof recomputes literal actual-vertex guards, verifies 400 exact continuous quadratic inequalities over `Q(sqrt(2),sqrt(5))`, and uses checked derivative identities and positive lower bounds to cover the two unrestricted own-contact magnitudes. No dotted distances are fixed. This is a compact-realization obstruction using actual positive-definite vertex Grams.

After integrating this gate, the [next bounded search](../scripts/fresh_h9_t8_continuation_20261014_next.py) stops at the first [remaining mandatory assignment](fresh_h9_t8_continuation_20261014_next_case_08.json), after 934 nodes. Its counters are 11,842 spherical-vertex failures, 511 orthogonal-negative-block failures, and 3/1/1 applications of the three frozen metric-family gates. These are traversal counters, not distinct incidence counts; the two searches share their initial traversal.

The remaining state has support blocks

```
[1,2] = I2(m), [4,5] = B2, [7,8] = B2, [12,13] = I2(n),
m,n finite integers >=7.
```

Anchor 11 has own labels `(2,3)` at roots `(12,13)`, while anchor 17 has `(3,3)`. This differs from the newly excluded opposite-root state only at `(13,17)`, whose label rises from 2 to 3. The asymmetric one-contact/two-contact Schur calculation remains open. A next metric test should derive its own projection budgets and continuous rank constraints; the current family theorem cannot simply be transferred through that change.

The saved state satisfies the mandatory Lannér-triangle, actual spherical-vertex, and necessary bridge filters. It is **not** a complete assignment of all meeting contacts and is not a Lorentzian metric candidate yet. All category-7 occurrences retain their full finite-label tail; no cutoff or representative angle is used. The searches stop at their first survivor, so they do not exhaust case 8.

Generation and saved replays pass for both searches and the metric gate:

```
python scripts/fresh_h9_t8_continuation_20261014.py --case 8 --limit 100000 --verify
python scripts/fresh_h9_t8_continuation_20261014_opposite.py --verify
python scripts/fresh_h9_t8_continuation_20261014_next.py --case 8 --limit 100000 --verify
```

The metric producer uses exact quadratic arithmetic for fresh LDL guards and SymPy algebraic-field arithmetic for its polynomial table. The [independent review](fresh_h9_t8_continuation_20261014_opposite_review.md) reconstructs all 400 determinant/Bernstein cases and literal vertex guards with separate Fraction arithmetic, checks the general identities with a separate rational-polynomial implementation, and reviews the complete own-label extension. This independent gate review does not duplicate the CSP traversal or the earlier incidence catalogue. No further search is running in this phase.
