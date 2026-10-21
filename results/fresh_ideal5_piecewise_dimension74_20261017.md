# Unrestricted finite-volume dimension at most 74

The completed local three-pair replacement and the exact global comparison
give the following computer-assisted conclusion:

**Every finite-sided, finite-volume real hyperbolic Coxeter polytope has
dimension at most 74.**

There is no arithmetic, compactness, simple-cusp, cusp-factor or facet-count
restriction. This strengthens the completed
[dimension-75 proof](fresh_ideal5_piecewise_dimension75_20261017.md).
The [source](../scripts/fresh_ideal5_piecewise_dimension74_20261017.py) and
[certificate](fresh_ideal5_piecewise_dimension74_20261017.json) check the new
local interfaces and all 999 reduced profiles at dimension 75. The
[independent review](fresh_ideal5_piecewise_dimension74_20261017_review.md)
reconstructs every profile and agrees on every exact row.

## What changes locally

The [order-15 affine B reduction](fresh_order15_affine_b_pair_reduction_20261017.md)
removes exactly three invariant paid pairs, (0,11), (11,13) and (11,14), from
the canonical affine B14 diagram. Its order is 15 and rank is 14. The cost
falls from 100 to 97. The three pairs have distances 11,12,12, so no basic
short-path or elementary-triangle payment is lost.

The complete saved induced-map census assigns these pairs to exactly seven
overlap contexts. Six have a preferred whole finite D14 that no supporting
white can touch. In the seventh, the marked A12 can only remain whole or
acquire a single residual-zero endpoint extension to whole A13. That added
root lies in the Lorentzian core span and cannot attach to core-orthogonal
whites; its exact products exclude every second core-contact white. In all
seven cases the rest of the corner is an orthogonal positive-definite white
block. Thus these are actual finite corners, paid by the existing whole
component cutoff 17. An abstract finite subdiagram at an ideal corner would
not suffice; the local proof establishes the complete finite outcome.

No finite coupon or finite cost is added. Every other selected pair and
local replacement stays as in the completed finite8/ideal5 theorem and its
dimension-75 refinement. All triangle payments remain at least two and
quadrilateral payments at least one. The finite cost remains at most 8n
through dimension 80. The prior D17 singleton replacement and its finite
A4/D13 coupons are retained without change.

At order 15 the other costs are A75, C95 and D96, so the largest cost is
now 97. The normalized ideal-factor envelope at affine rank r becomes:

|Rank r|Normalized cost|
|---|---:|
|1 through 12|r/2|
|13|87/14|
|14|97/15|
|15|111/16|
|16|126/17|
|17 and larger|163/23|

The previous arbitrary-rank estimates are unchanged. Actual nonsimple-cusp
incidence multiplicities still multiply each factor cost by its link vertex
count divided by that factor's order, as in the inherited factor-fiber
theorem. The only new envelope entry is at rank 14.

## Complete concentration and exact comparison

For a cusp of total rank d=n-1 with factor ranks r_i, put

    H(t)=product_i(1+t+...+t^r_i), V=H(1),
    G=sum_j [t^j]H(t)*min(j,d-j).

The normalized curvature term is (5G-dV)/(2V). Subtracting the sum of the
costs in the table gives the ideal margin. The inherited convex-order
lemma permits concentration (a,b) to (a-1,b+1), with 1<=a<=b, without
increasing normalized curvature. The cost must simultaneously be preserved
or increased; this is checked again for the altered envelope.

Low parts concentrate until one vanishes or the other reaches 12. Duplicate
13 parts may be replaced by 12 and 14, raising cost by exactly 4/105.
Duplicate 15 parts may be replaced by 14 and 16, raising it by exactly
7/2040. Both values are strictly positive, although small. Two parts at
least 18 concentrate until one reaches 17, preserving their constant cost.
Each step strictly increases the bounded integer sum of squared parts at
fixed total rank, so moves that recreate a middle rank cannot cycle.

Every partition therefore reduces to repeated 12,14,16,17, at most one
13 and one 15, at most one part from 1 through 11, and at most one part
at least 18. The source enumerates every such partition of 74.

|Ambient dimension|Reduced profiles|Minimum ideal margin|Finite slack|
|---:|---:|---:|---:|
|75|999|63/74698|95/2|

The minimum is attained at (10,12,12,12,12,16). Every margin is an exact
positive rational number; no floating tolerance is used. The finite
coefficient is beta_75=74*70/8, so beta_75-8*75=95/2>0. The inherited
global shape inequality also has strictly positive constant 5*(75-1)=370.
Summing the retained local face payments and comparing them with these
finite and ideal budgets gives a contradiction at dimension 75.
The completed dimension-75 theorem excludes every dimension at least 76.
Together they prove the stated ceiling 74.

At ambient dimension 74, the formal rank profile (1,12,12,12,12,12,12)
has margin -168143/742586. Thus this envelope comparison does not give
a dimension-73 bound. This profile is a scalar obstruction to the present
argument; no polytope realizing it is asserted.

The separate compact arithmetic ceiling 21 and arbitrary compact ceiling
29 are unchanged. The unrestricted result makes the effective ambient
ceiling for the compact-facet alternative 74; its older auxiliary proof
still gives 79. No new finite-volume construction is claimed.

## Verification status

Primary generation and a fresh separate-process saved replay both terminated
exit0/PASS for the local interfaces and all 999 exact profiles.
The local three-pair proof has passed generation, fresh replay and two
internal mathematical/source reviews. The complete inherited overlap and
map censuses retain their frozen validation records and are not rerun here.
The separate standard-library scalar reconstruction agrees on every one
of the 999 rows, the minimum, concentration increments and boundary. It
uses a different enumeration and both polynomial and cumulative curvature
calculations. Two complete internal source/mathematical reads of the
assembled argument found no identified gap in the local interfaces,
cost envelope, concentration, global inequality or inherited higher tail.
Those reads did not repeat the inherited geometric censuses. Direct
source bindings and Python syntax pass. Source, data and report
are final/frozen.

All reviews are internal to the agent workflow. External mathematician
refereeing and Lean or other proof-assistant formalization are not claimed.

```
python3 scripts/fresh_ideal5_piecewise_dimension74_20261017.py
python3 scripts/fresh_ideal5_piecewise_dimension74_20261017.py --verify
```
