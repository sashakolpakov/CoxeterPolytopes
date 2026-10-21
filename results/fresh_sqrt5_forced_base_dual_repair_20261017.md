# The forced-base extension still has an exact bound above 12

Adding the missing base stars and their endpoint charges does not remove the
current partial-model obstruction. The exact analytically repaired dual gives

\[
 \boxed{\frac{31154121494284588800}{2547155390157199561}
 \;\approx\;12.2309465746\;>12.}
\]

The [source](../scripts/fresh_sqrt5_forced_base_dual_repair_20261017.py) and
[complete sparse certificate](fresh_sqrt5_forced_base_dual_repair_20261017.json)
reconstruct the model from frozen inputs, without a numerical optimizer or
a duplicate full-domain export. Since the exact analytic bound already
exceeds 12, the conditional optimization branch was not needed.

## Domain and repair

The domain combines the frozen
[pair, clique and disjoint-port model](fresh_sqrt5_rank13_combined_disjoint_ports_20261017.md),
[represented shared prism/base transfers](fresh_sqrt5_prism_shared_transfers_20261017.md),
and the [forced missing-base patch](fresh_sqrt5_forced_base_simplex_domain_20261017.md).
It has:

| Quantity | Count |
|---|---:|
| Complete finite vertex budgets | 6,558 |
| Positive-cost variables | 211,339 |
| Face and incidence rows | 82,530 |
| Active positive-cost variables | 85,841 |
| Original represented shared transfers | 2,845 |
| Added missing-base shared transfers | 2,848 |
| Total zero-cost signed transfers | 5,693 |

Every old variable, face row, budget entry and represented transfer remains.
The 22 old clique charges have zero on the new rows, as explicitly prescribed
by the frozen patch. Existing one-port variables are reused at the new
simplex corners, rather than receiving a separate budget.

Keep all quad masses and all represented simplex masses from the previous
shared-capacity repair. For each new missing base type S, set its mass to

\[
 \alpha_S=
 \frac{\sum_Q \underline m(Q,S)\,k_Q\alpha_Q}
 {(d_S+1)(13-d_S)},
\]

where the underlined multiplicity is the minimum over **every** provenance
of the corresponding merged quad row. There are 2,537 positive new masses:

| New base dimension | Positive masses |
|---|---:|
| 1: maximal edges | 1,904 |
| 2: triangles | 402 |
| 3: tetrahedra | 150 |
| 4: four-simplexes | 81 |

The maximal-edge rows contribute zero to the objective. Their two endpoint
port charges still increase relevant vertex loads. Higher simplexes add
their actual triangular demand, `2 C(d+1,3)`. The other selected missing-base
rows retain zero dual mass; they remain present in the model.

All 5,693 signed transfer loads are checked exactly to be nonpositive. The
repair preserves every old face mass and raises only previously absent base
rows. This is a prescribed analytic repair, not an optimized dual.

## Exact verification

For every positive-cost variable j attached to vertex type v, recompute

\[
 L_j=\sum_f a_{fj}\alpha_f,\qquad
 \beta_v=\max_{j\text{ at }v}\frac{L_j}{c_j},
\]

where `c_j` is its literal pair, clique or port budget multiplicity. Every
positive-cost column is checked against `L_j <= c_j beta_v`. Columns absent
from the positive support have zero load. Independently, every zero-cost
transfer column is checked against

\[
 \sum_Q k_Q\underline m(Q,S)\alpha_Q
 -(d_S+1)(13-d_S)\alpha_S\le0.
\]

The result has 5,344 positive face masses and 752 positive vertex masses;
491 vertex masses increase. Its exact objective numerator and normalization
are

\[
 \sum_f b_f\alpha_f=\frac{1620989504989}{123750000000},\qquad
 \sum_v\beta_v=\frac{2547155390157199561}{2378376000000000000}.
\]

The additional demand is `29615942231/123750000000`; the additional vertex
mass is `107204507717238121/2378376000000000000`. Their recomputed ratio gives
the displayed bound. The certificate saves every positive face and vertex
mass, every maximizing positive-cost column, and all 5,693 signed loads.

## Meaning and limits

This rules out a budget-12 rule in this expanded weighting family on the
saved partial face domain. It does not prove the existence of a compact
polytope or an actual global incidence distribution, and it gives no new
dimension upper bound. Further global compatibility or stronger incidence
constraints may exclude this dual distribution.

The separate mixed-three-port charges are **not** included. Their features
on the newly activated base stars must be reconstructed before their earlier
diagnostics can be combined with this certificate. Remaining absent base
stars, fixed arithmetic-form classes, other simplex cores, all-finite quad
aggregates, other faces and 231 untouched finite forests remain outside.

**FINAL/FROZEN:** generation and a fresh saved reconstruction replay passed.
Direct source manifests, Python syntax and local report links were checked.
The exact verifier reuses the frozen sparse-dual engine; no independent
second arithmetic implementation or external review is claimed.
