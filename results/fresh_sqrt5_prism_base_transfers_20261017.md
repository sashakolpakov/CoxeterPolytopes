# Represented prism-base transfers and an exact repaired dual

The represented-base incidence constraints reject the previous dual, but
they do **not** remove the budget-12 obstruction. An explicit analytic repair
gives the exact lower bound

`694360866673186440 / 55273101561404603 > 12`,

approximately **12.5623648223**, for the family enlarged by every exported
represented-base transfer. No optimizer was run and this lower bound is
not claimed optimal.

The [source](../scripts/fresh_sqrt5_prism_base_transfers_20261017.py) and
[certificate/export](fresh_sqrt5_prism_base_transfers_20261017.json) bind the
[prism-base diagnostic](fresh_sqrt5_prism_simplex_base_incidence_20261017.md)
and the [combined disjoint-port model](fresh_sqrt5_rank13_combined_disjoint_ports_20261017.md).

## Signed columns with zero vertex cost

For each represented maximal simplex type S and prism dimension k, let
`C = binomial(d+1,k)(13−d)`. The incidence inequality is

`sum_Q m(Q,S,k) f_Q − C f_S <= 0`.

For each merged quad row the exported m is the minimum multiplicity across
**all** its provenances. Consequently this inequality remains valid for
every distribution of the row mass among those stars. Unknown base types
are retained in the provenance data but are not assigned a zero mass.

A nonnegative transfer parameter has coefficient +m on each quad row,
coefficient −C on the simplex row, and vertex cost zero. It transfers part
of the required payment from the simplex aggregate to prism aggregates.
The sign is essential: its dual condition is

`sum_Q m alpha_Q − C alpha_S <= 0`.

All 68,150 quad rows and all 81,183 source provenances are processed. The
export has 3,840 transfer columns, 6,316 positive quad entries and one
negative simplex entry per column. There are 11,442 quad rows whose base
coefficients vary among provenances; taking only the first would be invalid.
The data reference the full frozen provenance maps instead of duplicating
their geometry. Missing base type/dimension pairs remain explicit unknowns.

## Analytic repair

Keep every quad mass fixed. For each represented simplex row independently,
replace its mass by the maximum of its old mass and every required ratio
`(sum_Q m alpha_Q)/C`. This is the least componentwise repair with the old
simplex masses as lower bounds; it does not optimize the resulting dual
objective.

The repair raises 419 simplex masses and makes 410 previously zero simplex
rows positive. It repairs all 434 original violations. The resulting dual
has 2,807 positive face rows and 672 positive vertex masses. The checker
recomputes each vertex mass as the exact maximum column load divided by its
positive budget cost, checking every inherited pair, clique and disjoint-port
column. Every new zero-cost column is checked separately and has nonpositive
dual load. Its exact load is saved, including slack and equality cases.

This proves the displayed lower bound for the augmented partial model.
It is stronger evidence than merely observing that the old fixed dual
violates new constraints: the obstruction survives an explicit valid repair.
It does not exclude stronger geometric conditions or missing base-star
constraints.

## Scope and validation

Generation and a fresh saved replay both passed. The certificate
includes only the inherited pair weights, 22 earlier simplex clique charges,
6,528 disjoint-simplex-port charges and the represented-base transfers.
The separate mixed-three-port columns are **not** included: activating new
simplex rows changes their loads, so their earlier check cannot be inherited.
Separate capacities for different k are used; no simultaneous antichain or
product-completion claim is assumed.

There is no new dimension bound, feasible global rule, complete face
catalogue or assertion of simultaneous compact ambient realization.
