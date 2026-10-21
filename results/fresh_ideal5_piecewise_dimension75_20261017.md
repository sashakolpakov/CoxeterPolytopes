# Dimension at most75 from the D17 singleton transfer

Complete internally checked computer-assisted proof; validation scope is below.

Every finite-sided, finite-volume real
hyperbolic Coxeter polytope has dimension at most75. No arithmetic,
compactness, cusp-shape or facet-count restriction is imposed. The completed
[dimension77 proof](fresh_ideal5_piecewise_dimension77_20261017.md) supplies
the higher-dimensional tail; the new comparison excludes76 and77.

The change is one pair in the17-node affine D factor. The
[singleton-pair transfer](fresh_d17_singleton_pair_coupon_20261017.md)
replaces it by two finite half-payments, while the
[enlarged finite budget](fresh_finite8_a4_d13_neighbor_coupon_budget_20261017.md)
retains cost8n through dimension80. The complete
[finite8/ideal5 theorem](fresh_finite8_ideal5_combined_payments_20261017.md)
provides the rest of the face coverage. This composition uses the completed
2,614-context classification, not the older no-short-nine hypothesis.

## Local change and its coverage

In the standard17-node affine D diagram, the pair(1,11) on the simple
spine has distance10. It is fixed as an unordered pair by every diagram
automorphism. Removing it changes the D cost from127 to126. The other
order17 types retain costs85,111,108, so the worst normalized factor cost
changes from127/17 to126/17. Every other rank cost is unchanged.

The complete1,153-pattern map census requires this pair in exactly one
whole marked-core context, index2518 of the new2,614-context domain. Its
positive affine kernel forces every extra supporting white to be right
to the17-node affine component. The core radical then makes its contacts
with the two remaining black roots equal. Nonzero contacts supply a paid
path of length2. If they vanish, the entire core is sealed. Its two adjacent
finite corners have whole A4 and D13 components, with the A4 mark an endpoint
and the D13 mark the unique long-arm neighbor of the fork center. Each
receives the new half-payment, giving total1.

This literal whole-core argument is necessary: an abstract affine embedding
alone would not seal the finite components. The complete overlap theorem
excludes every extra cross-core bond, and the local lemma handles all extra
supporting whites. The new pair has distance10, so all elementary triangle
payments and short-path reductions remain unchanged. The pre-existing
removed D pairs retain their own sealed or universal transfers.

## Finite counting and the ideal envelope

Each pair of whole A4 and D13 components contributes at most one additional
unit: two A4 endpoints, one specified D13 node, half a unit each. If a is
the total rank of small components and b is the number of D13 components,
the extra cost is at most ab/4. Every D13 has base slack26 rather than the
generic8 already charged in the finite8 proof. The resulting excess bound
is the old bound plus b(a/4-18). When b>0 the rank constraint gives a<=67,
making this correction negative; when b=0 it vanishes. Thus the full
finite-vertex budget remains8n, including all old coupons.

The new normalized ideal-factor envelope, for affine rank r, is

    r/2 for1<=r<=12;
    87/14 at13;20/3 at14;111/16 at15;126/17 at16;
    163/23 forr>=17.

Actual cusp-factor multiplicities are those of the inherited factor-fiber
theorem: factor pair costs multiply by the link-vertex count divided by
factor order. No simple-cusp assumption is introduced.

## Exact global comparison

Use the inherited convex-order theorem: transferring(a,b) to(a-1,b+1),
with1<=a<=b, cannot increase normalized curvature. Concentrate two ranks
at most11 until one vanishes or the other reaches12; their linear budget
is preserved. Replace duplicate13 by12 and14, increasing cost by5/21.
Replace duplicate15 by14 and16, increasing cost by83/408. Concentrate two
ranks at least18 until one reaches17; their constant cost is unchanged.
Every move increases the bounded integer sum of squared ranks, so it
terminates. Terminal profiles have arbitrary multiplicities of12,14,16,17
and at most one each of13,15,a rank1..11 and a rank at least18.

The [checker](../scripts/fresh_ideal5_piecewise_dimension75_20261017.py)
enumerates every terminal profile in the two remaining dimensions. For
`H(t)=product(1+t+...+t^r)`, put V=H(1) and
`G=sum_j [t^j]H(t) min(j,d-j)`, where d is total affine rank. The inherited
normalized curvature is `(5G-dV)/(2V)`. Subtracting the new envelope gives:

|Ambient dimension|Profiles|Minimum ideal margin|Finite margin|
|---:|---:|---:|---:|
|76|1053|11615/448188|57|
|77|1111|41558/485537|68|

Both ideal and finite margins are strictly positive, and the inherited
global shape inequality retains the positive constants375 and380. Hence
both dimensions are impossible. The77 theorem excludes every larger
dimension, giving75.

At ambient75 the formal profile(12,12,12,12,12,14) has negative margin
`-16873/171366`. This is a limitation of the comparison, not a realized
hyperbolic polytope. This package makes no dimension74 claim.

## Verification status

Primary generation and a fresh separate-process saved replay terminate
exit0/PASS, including all2164 exact scalar rows and the local certificate
interfaces. Two separate full source/mathematical reads pass the complete
composition, including the new whole-component coupon, the short8-first
reduction and the inherited higher-dimensional tail.

The separate standard-library-only
[independent reconstruction](fresh_ideal5_piecewise_dimension75_20261017_review.md)
matches every primary scalar row, both minima and the boundary profile.
It independently enumerates the reduced partitions and computes curvature
using sliding convolution, cumulative inclusion-exclusion and an
absolute-deviation identity. These checks do not duplicate the upstream
2,614-context classification,1,153-pattern embedding census or every
inherited geometric proof; those retain their frozen validation records.
All reviews are internal. External refereeing and proof-assistant
formalization are not claimed. Source, data and report are final/frozen.

Compact arithmetic21 and arbitrary compact29 are unchanged. The
unrestricted result also gives an effective compact-facet ambient ceiling75,
while the older auxiliary compact-facet proof remains79.

    python3 scripts/fresh_ideal5_piecewise_dimension75_20261017.py
    python3 scripts/fresh_ideal5_piecewise_dimension75_20261017.py --verify
