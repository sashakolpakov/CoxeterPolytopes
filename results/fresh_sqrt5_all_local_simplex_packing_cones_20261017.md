# All 740 local simplex-packing cones retain a bound above 12

Every positive-simplex forest in the saved dual now has an explicit rational
mixture of disjoint literal port blocks covering all its required simplex
corner marginals. After the corresponding vertex-mass repair, the dual bound
is

\[
 \boxed{\quad
 {69774201720926976000000000000000
  \over5448173768070467817598921440139}
 \simeq12.80689726342>12.\quad}
\]

Consequently, allowing **all nonnegative charges depending on the finite
forest and the maximal-simplex port-block orbit**, with vertex cost determined
by disjoint-block packing, cannot bring this partial model to budget 12.
This is a complete certificate for that charge family on the current dual's
positive support. It is not a dimension theorem or an ambient realization.

## Result and exact status

The [fourteen-forest pilot](fresh_sqrt5_local_simplex_packing_cones_20261017.md)
is reused unchanged. The remaining 726 forests are processed in decreasing
current vertex mass, then decreasing occupied-orbit count. No global LP is
run.

| Check | Result |
| --- | ---: |
| Positive-simplex forests covered | 740 / 740 |
| Positive block orbits covered | 6,011 |
| Distinct packing count vectors enumerated across forests | 608,079 |
| Explicit positive mixture terms | 5,255 |
| Matching exact local primal/dual optima | 623 |
| Exact lower/upper brackets | 117 |
| Open cases or resource stops | 0 |
| New exact cuts beyond subset and clique corrections | 23 |

The 117 bracket cases are **not** reported as exact local optima. Their
rational packing covers nevertheless prove the upper bound needed for the
universal charge argument below. The largest bracket width is less than
`6.8×10^-10`; their total width is

`261530894195162461841 / 19800000000000000000000000000`

or about `1.32×10^-8`. The final global denominator uses the exact covering
masses, including these conservative increases. There are 121 increased
vertex masses: 23 have genuine positive lower-bound cuts, while the other
increases come from the conservative exact covers.

Each new local solve had a twenty-second cap. Enumeration was limited to
one million cached count-vector entries and thirty seconds per forest, with
any stop required to remain explicitly open. None occurred. The largest
single count-vector list had 14,137 entries; the largest total cache had
157,453. Generation completed in about seventeen seconds, and the fresh
optimizer-free replay in about nine seconds.

## Why the certificate covers every nonnegative orbit-local charge

The [full subset package](fresh_sqrt5_all_subset_port_packing_20261017.md)
reconstructs every corner of all 3,574 positive simplex rows in the saved
optimized dual. It identifies the complete port-block orbit within the full
finite vertex forest, and checks that these blocks reproduce the inherited
pair and one-port coefficients. Let w(F,j) be the resulting total corner mass
of orbit j at forest F.

Fix F. A packing P is any set of pairwise disjoint selected literal blocks.
Let a(P,j) count its blocks in orbit j. The saved certificate supplies
nonnegative rational coefficients lambda(P) with

\[
 \sum_P\lambda(P)a(P,j)\ge w(F,j)
 \quad\hbox{for every }j,
 \qquad
 L_F=\sum_P\lambda(P).
\]

Every positive lambda comes with explicit literal block masks, whose mutual
disjointness and orbit membership are checked. Thus the covering property
does not depend on the optimizer, numerical tolerance, or the completeness of
the count-vector enumeration.

For **any** nonnegative orbit charge q(j), write

\[
 c_F(q)=\max_P\sum_j a(P,j)q(j).
\]

Multiplying the marginal cover by q(j) and summing gives

\[
 \sum_j w(F,j)q(j)
 \le\sum_P\lambda(P)\sum_j a(P,j)q(j)
 \le L_F c_F(q).
\]

This inequality applies to the whole cone of nonnegative charges, not merely
the finitely many optimal charges exported by the solver. Setting
`mu_new(F)=max(mu_old(F),L_F)` therefore satisfies every new orbit-local dual
column simultaneously. A zero-marginal block can be omitted from the test:
a nonnegative charge on it adds no load and cannot reduce the maximum packing
cost. Forests with no positive simplex marginal likewise add no constraint.

The actual [disjoint-port theorem](fresh_simplex_disjoint_ports_20261017.md)
ensures that a vertex's maximal-simplex blocks form such a packing. Complete
automorphism orbits are retained. Averaging a whole packing over the
finite-forest automorphism group produces uniform marginals within each
orbit; symmetrizing a literal charge preserves that uniform-marginal load and
cannot increase its maximum packing cost. This also explains why orbit-count
mixtures are the appropriate exact representation.

All original face multipliers alpha remain fixed. The previous
[intersecting-clique repair](fresh_sqrt5_intersecting_block_cliques_20261017.md)
already satisfies every inherited positive-cost column and all signed
transfers. Increasing mu preserves the former, and unchanged alpha preserves
the latter. The [pointed-base check](fresh_sqrt5_pointed_simplex_base_capacity_20261017.md)
requires no further repair. Hence the unchanged weighted face demand D and
new total vertex mass M form a valid dual against the entire enlarged local
charge cone. A primal budget 12 would imply `D<=12M`, contradicting the exact
ratio displayed above.

This last step is why complete exact **covering mixtures** suffice even in the
117 cases where the best local value is only bracketed. Matching lower bounds
are useful diagnostics, but are not a premise of the universal obstruction.

## New cuts and a small explicit example

The largest new cut occurs at forest `3A1+A3+A5+I2(5)` (vertex index 947),
increasing its required mass by `85790137/165000000000`. Its optimal charge
uses half- and quarter-integral coefficients. All 23 new charge columns are
exported with exact rational coefficients on the old positive simplex rows.

A particularly small new cut occurs at forest `2A2+D4+H3+I2(5)` (index 4357).
In its canonical thirteen-port numbering the nonzero charges are:

- charge `1/3` on the four-port blocks `{0,1,11,12}` and `{2,3,11,12}`;
- charge `2/3` on the block `{4,5,6,7}`;
- charge `1/3` on every pair `{a,b}` with `a∈{0,1,2,3}` and `b∈{4,6,7}`.

Their maximum disjoint-packing charge is one. If the second block is present,
no charged pair is available and at most one first-family block is available.
Without it, at most three charged pairs fit; a first-family block uses two of
their four left ports, leaving at most two charged pairs. Thus every packing
has charge at most one, and the saved mixtures attain the resulting bound.
Its exact load `220013/1200000000` exceeds the prior mass
`1677213/10000000000` by `234343/15000000000`.

## Reproduction, limits, and validation

- [Exact source](../scripts/fresh_sqrt5_all_local_simplex_packing_cones_20261017.py)
- [All local charges, explicit packing mixtures, exact brackets, and repaired dual](fresh_sqrt5_all_local_simplex_packing_cones_20261017.json)

Run `python3 scripts/fresh_sqrt5_all_local_simplex_packing_cones_20261017.py --verify`
to reconstruct every count-vector list and verify the exact certificates.
Replay invokes no optimizer. The generation progress file lives in the
temporary directory and is not an additional mathematical input.

**FINAL/FROZEN.** Generation and a separate fresh saved replay passed, followed
by syntax, direct source hashes, and local report-link checks. The root's
[independent Fraction cover review](review_sqrt5_all_local_packing_cover_20261017.md)
also passed generation and fresh replay. Without importing producer code,
running an optimizer, or enumerating count vectors, it checks all 740 forests,
6,011 marginal orbits, 5,255 explicit packing terms, and 23,324 literal blocks,
and independently reproduces the same repaired bound. The independent check
proves the universal covering conclusion; it does not duplicate the 623 local
optimality claims. These are internal reviews, not external refereeing or
proof-assistant formalization.

The remaining structural possibilities include geometric incompatibility
among otherwise disjoint blocks, charges depending on nonincident walls of
the full simplex star, joint simplex/quadrilateral configurations, and signed
relations between faces or coface flags. These are outside the local charge
cone checked here. Missing global face coverage and common ambient
realizability also remain separate. No new ambient dimension bound follows
from this method obstruction.
