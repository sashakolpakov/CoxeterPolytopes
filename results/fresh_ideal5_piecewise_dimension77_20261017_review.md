# Independent review of the ideal-five dimension comparison

The full written composition and independent exact comparison pass: the
argument proves dimension at most77 for every finite-sided, finite-volume
real hyperbolic Coxeter polytope. No arithmetic or cusp-shape restriction
is imposed. The [checker](../scripts/fresh_ideal5_piecewise_dimension77_20261017_review.py)
imports only standard-library arithmetic. It independently enumerates3,722
reduced profiles in ambient dimensions78,79,80, computes curvature in two
different ways, and derives15 all-rank affine cost branches.
The [certificate](fresh_ideal5_piecewise_dimension77_20261017_review.json)
compares these results with the [primary proof](fresh_ideal5_piecewise_dimension77_20261017.md).

## Scalar calculations

The envelope under examination is

    b(r)=r/2 for1<=r<=12;
    b(13)=87/14; b(14)=20/3; b(15)=111/16; b(16)=127/17;
    b(r)=163/23 forr>=17.

The checker chooses the optional ranks1,...,11,13,15 and at least18, then
recursively enumerates all multiplicities of12,14,16,17 at the remaining
total rank. Sorted shapes are deduplicated; no producer list is imported.
Sliding-window convolution computes

    H(t)=product_i(1+t+...+t^r_i),
    V=product_i(r_i+1),  G=sum_j h_j min(j,d-j).

A separate inclusion-exclusion calculation expands the numerator
product_i(1-t^(r_i+1))=sum_l a_l t^l and uses

    G=floor(d/2)*V-2*sum_l a_l*binom(floor(d/2)-l+s,s+1).

Both exact values agree. Polynomial symmetry, H(1)=V and the independent
absolute-deviation identity for normalized curvature are checked for every
shape. Curvature minus the proposed cost gives:

| Ambient dimension | Profiles | Minimum margin | Finite slack beta_n-9n |
|---:|---:|---:|---:|
|78|1173|409611/2171546|0|
|79|1239|63437/167042|21/2|
|80|1310|3585185/7619196|20|

The least margin occurs at(1,12,16,16,16,16). The ambient77 boundary shape
(12,16,16,16,16) has margin-16664/1085773, and its finite slack is-9. This
is a limitation of the proposed comparison, not an actual chamber.

## Complete shape reduction

Use the inherited convex-order theorem: transferring(a,b) to(a-1,b+1) for
1<=a<=b cannot increase normalized curvature, with all other factors fixed.
Two ranks at most11 can be concentrated until one vanishes or the other
reaches12; their linear budget is preserved. The new rank13 cost is less
than13/2, so the older reduction that concentrated through13 is not valid.

Instead replace duplicate13 by12 and14. This increases the budget by5/21.
Replace duplicate15 by14 and16, increasing it by107/408. Concentrate two
ranks at least18 until one becomes17; their constant budget is unchanged.
Each move therefore cannot increase curvature minus cost. It strictly
increases the bounded integer sum of squared ranks, proving termination.

Terminal partitions have arbitrary multiplicities of12,14,16,17, at most
one rank1,...,11, at most one13, at most one15 and at most one rank at least18.
The independent enumeration above exhausts that domain. Positivity for all
terminal profiles proves positivity for every partition of the same rank.

## Independent all-rank ideal counts

For factor order m>=18, the proposed ideal rule pays distance at most5 and
distance6,...,9 when a mark is within5 of a branch or nonsimple-bond feature.
The independent spine-interval formulas are

    A: 5m;
    B: 5m-6  +sum_{s=8,9,10,11} min(m-s,13);
    C: 5m-15 +sum_{s=6,7,8,9} min(m-s,14);
    D: 5m+3  +sum_{s=10,11,12,13} min(m-s,12).

For B the short-distance contribution is5m-14, and each of the four longer
distances adds two fork-leaf pairs plus the spine window of width13. For D
the short contribution is5m-13; four leaf pairs at each longer distance give
the constant16 adjustment, and the spine window has width12. C has a path
with two end feature windows of combined width14. Opposite fork ends are
too distant to produce an extra leaf pair at m>=18.

Each min term has one breakpoint. The five branches per type have cost/m
of the form a+b/m, whose extrema occur at endpoints or the limiting value5.
Independent exact maxima are160/23 for B,75/11 for C and163/23 for D.
This establishes the all-rank envelope from the stated rule, without
extrapolating a finite graph sample. Small-order preferred-pattern coverage
must separately justify the exceptional rank13,...,16 values.

## Finite budget and composition obligations

The [enlarged finite-eight budget](fresh_finite8_a17_d11_coupon_budget_20261017.md)
is at most8n through n=80, including whole-small rank5, D11 and A17 coupons.
Changing from finite-eight/leaf-nine to full finite-nine adds at most n:
for a component of rank r>=18, the additional nonleaf distance-nine pairs
number r-11 in types A/B and r-12 in type D. Smaller components were already
fully paid. Thus the enlarged finite-nine rule costs at most9n through80.

At n=78 the finite coefficient equals9n, rather than strictly exceeding it.
The inherited shape inequality nevertheless retains its strictly positive
constant5(n-1)=385. Thus nonnegative finite slack and strictly positive
ideal margins give the contradiction, even when the polytope has no ideal
vertices. The positive global slacks at79 and80 are390 and395. The completed
dimension80 theorem supplies the higher tail.

## Complete local composition review

The full source and mathematical argument of the
[combined local theorem](fresh_ideal5_combined_payments_20261017.md), the
[small-factor census](fresh_ideal5_small_factor_payments_20261017.md) and
the primary comparison were read. The following interfaces close the
geometric obligations; the scalar enumeration alone would not suffice.

The finite-nine rule contains the prior finite-eight triangle rule.
At ideal factors of order14 through17, a separate breadth-first distance
reconstruction checks that every elementary triangle pair remains in the
saved paid set. The16 target costs agree with the stated rank envelope.
No replacement paying a quadrilateral only once is mistaken for a triangle
payment of two.

For a quadrilateral with a marked path of length at most9, the retained
finite rule or the complete ideal-five local transfers supply a unit or
sealed finite half-coupons. The order14 cycle uses its separate strict
spectator-spectrum transfer. The universal affine-D13/D14 transfers have
immediate small affine alternatives of orders12/13, both wholly paid.
All other terminal alternatives are short paths or the listed whole finite
components. Their source and proof reviews are recorded in the local
packages.

After this reduction, the original no-short-nine hypothesis holds for the
complete949-case classification. Its low-core supplement has component
bound at most12, and the571 remaining large-bound patterns are exhaustively
embedded with labels and nonedges respected. Orders16/17 D deletions have
minimum distances11/12, checked independently here; thus they never require
their context-specific transfer before the no-short-nine reduction. After
classification their affected context sets match the saved sealed-core
proofs. Every nonzero terminal payment in those18 cases is finite, and a
nonorthogonal supporting root supplies a retained length-two path. No
unpaid whole affine order14 conclusion or cyclic ideal transfer is imported.

The all-rank spine counts and actual cusp-factor multiplicity
V/(r+1) justify the ideal costs also at nonsimple cusps. The inherited
shape inequality retains its positive global-extremum contribution, and
the concentration moves preserve or increase costs while decreasing
curvature. The old reduction through rank13 would fail with the new budget;
the independently reconstructed12/14/16/17 reduction avoids that error.

No gap was found in this complete composition. It uses the finite-nine
classification; complete finite-eight overlap coverage remains a separate
open problem and is not needed for77.

## Verification status

The complete read-only run, generation and fresh separate-process saved
replay all terminated exit0/PASS, comparing every primary profile field,
all three minima and the boundary profile. The checker independently
reconstructed all15 infinite-order affine branches and checked the16
small-factor pair sets and126 finite-cost increment interfaces. A separate
full source/mathematical read of this reviewer checked the enumeration,
cumulative-tail identity, branch extrema, graph inclusions, no-short-nine
reduction and positive global slack, finding no gap. That peer read did
not independently rerun the3,722-row arithmetic or the571/949 censuses.

The independent checker verifies20 direct source/data/report bindings.
It does not rerun the571-pattern map census,949-core classification,
historical contact arithmetic or all preceding geometric proofs; those
retain their explicit frozen validation records. The new composition has
a separate full source and mathematical read. Source, data and report are
final/frozen after direct binding, syntax and local-link checks. All reviews are internal;
external mathematician refereeing and proof-assistant formalization are
not claimed.
