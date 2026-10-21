# Shared prism-base capacities: exact transfer export and repaired dual

The stronger capacity shared across all prism dimensions still leaves a
valid exact budget-12 obstruction. Analytically repairing the original dual
gives

`340741162849360 / 27183053503119 > 12`,

approximately **12.5350583889**. No optimizer was run. This is a lower bound
for the enlarged partial payment family, not its claimed optimum or a
dimension exclusion.

The [source](../scripts/fresh_sqrt5_prism_shared_transfers_20261017.py)
and [certificate](fresh_sqrt5_prism_shared_transfers_20261017.json) bind the
[shared-capacity proof](fresh_prism_simplex_shared_capacity_20261017.md),
[complete base-provenance map](fresh_sqrt5_prism_simplex_base_incidence_20261017.md),
and [separate-dimension transfer export](fresh_sqrt5_prism_base_transfers_20261017.md).

For a represented maximal simplex type S of dimension d, the shared
inequality is

`sum_(k,Q) k m(Q,S,k) f_Q <= (d+1)(13−d) f_S`.

The 2,845 exported signed columns therefore have positive coefficient k m
on each prism row, negative coefficient `−(d+1)(13−d)` on S's row, and zero
vertex cost. There are 6,316 positive quad entries. Each merged quad row has
one prism dimension k, so k times its minimum base multiplicity is exactly
the minimum shared coefficient over all represented provenances. Missing
base types remain unknown and receive no fabricated simplex row or mass.

This capacity implies each earlier separate-k capacity:

`k binomial(d+1,k) = (d+1) binomial(d,k−1) >= d+1`.

The checker verifies this comparison for each source column and explicitly
checks all 3,840 separate-k dual loads after the shared repair.

Starting from the original combined dual, all quad masses are kept fixed.
Each represented simplex mass is raised to the maximum of its old value and
its required shared load divided by `(d+1)(13−d)`. This raises 423 masses,
including 410 previously zero rows. The resulting certificate has 2,807
positive face rows and 672 positive vertex masses. It satisfies every shared
and separate transfer with nonpositive signed load. Every original pair,
simplex-clique and disjoint-port column is rechecked by exact rational
arithmetic after recomputing its necessary vertex mass.

The preceding separate-k analytic repair gave approximately 12.5623648223;
the stronger shared repair gives the displayed 12.5350583889. These are
explicit surviving dual certificates, not optimized comparisons of the
two families. The new repaired simplex masses dominate the earlier repair
componentwise, which is also checked.

Generation and a fresh saved replay both passed. The separately
studied mixed-three-port columns are outside this certificate because the
new simplex masses change their loads. Missing base types, unenumerated
faces, fixed-determinant stars and simultaneous ambient realization also
remain outside. The result shows that these represented-base transfer
capacities alone do not remove the target-12 obstruction.
