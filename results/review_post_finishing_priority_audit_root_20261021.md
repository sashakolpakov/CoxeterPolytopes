# First mathematical audit after the six-paper delivery

21 October 2026. This is a mathematical source read, not a new certificate
replay, compilation or external referee report. The six-paper receipt
`2c924e7e` remains the delivered edition. No manuscript was changed.

## The finite57 scalar counting spine

I independently checked the main counting spine of
`manuscripts/finite_volume57.tex`, source `c300db69`, at the explicitly
stated local-geometric and finite-output interfaces. These calculations
do not certify the full local atlas or its finite enumerations.

The payment comparison sums one nonnegative simultaneous corner rule.
Triangles require two, quadrilaterals one, and larger polygons need no
positive lower bound. Summing actual face-corner incidences once therefore
gives D5 <= sum C(v). The final finite maximum is an upper bound over a
formal superset of actual finite vertices; no realization of an attaining
forest is needed.

The residual Bellman certificate proposition needs BOTH directions.
All last-item inequalities bound every reachable multiset. A tight
predecessor for every stored positive-rank finite value realizes that
value; positive item ranks force its predecessor chain to reach the
reachable rank-zero base. An unreachable positive rank must have no
reachable predecessor. With the complete finite rank domain and alphabet
as inputs, these conditions exclude omitted and spurious finite entries.
Checking only attaining paths would not suffice. The printed statement
has all of these requirements.

For the eight preliminary excess forests, the printed costs are
78+3*100=378. Removing one total unit from each gives 377. The only
new positive-coupon families have costs 170+200+1/2=741/2 and
66+217/2+200+1=751/2, both below 377. This composition still depends
on the claimed exhaustive old forest list, the actual local transfer,
the rank-capacity classification of all newly increased forests and
the simultaneous validity of its compensating coupons.

For the ideal branch I rederived the transfer direction. At fixed sum,
U_(a-1)+U_(b+1) has the same mean as U_a+U_b. The displayed integer
stop-loss differences are nonnegative: in the second range its numerator
has derivative a+b-2t >= 0 and value (a+2)(b-a+1)/3 at t=a. Symmetry
covers the upper half, piecewise linearity the noninteger thresholds,
and zero outside the support the remaining thresholds. Conditioning
proves preservation on adding the other independent uniform factors.

Consequently the absolute-deviation expectation increases, so normalized
G, curvature and margin decrease. For a,b <= 7 the envelope cost sum is
unchanged, since b+1 <= 8 and a-1=0 is omitted with cost zero. The factor
product decreases by b-a+1>0. The integer square-rank potential increases
by 2(b-a+1), so the process terminates with all but at most one positive
rank at least eight. The stopping domain has the unique small-rank
branch and sorted tail printed in the paper. Positivity of the terminal
normalized margin is essential when multiplying V>=V' and m>=m'>0.
Given the complete terminal finite-output proposition, the resulting
actual ideal surplus is at least 36. I did not rerun its 2,552 profiles.

For the global shape calculation I checked the constants and signs:
sum_(j=0)^57 min(j,57-j)=2*(1+...+28)=812;
57*58/2=1653; (5*812-57*58)/2=377;
sum_(j=29)^57(2j-57)=1+3+...+57=841.
Two distinct extrema each contribute at least 841 extra uncut edges,
and their common coefficient 5/2 gives 4205. For a product-of-simplices
link, E=57V/2 and the h-polynomial has no missing degree, so the same
lower extrema correction is applicable. These are conditional on the
stated actual-link and generic-index incidence theorems, not an assertion
that arbitrary nonsimple links are products of simplices.

The upper-minus-lower comparison is therefore
0 >= 4205 + sum(κ-C_actual) >= 4205+36*v_ideal > 0.
It also excludes the all-finite case. Exclusion of dimension58 alone
does not imply the theorem: the complete higher-dimensional tail remains
an independent obligation, assigned to the separate composition reviewer.

Verdict: PASS for this scalar/counting spine at its stated interfaces.
This read supplies no full local atlas, native census, complete-runtime
or geometric Lean claim.

## Prepared two-foot formal kernel

I read all of `formal/CoxeterProofs/TwoFootPaidKernel.lean` (source
`cb282fd0`) and its map `docs/two_foot_paid_logic_20261021.md`, together
with the used finite-sum, weighted-pressure and completed-square lemmas
in its unchanged imports. All five mathematical statements follow from
the primitive hypotheses. This verdict precedes compilation.

Write P=sum a_i Ω_i, Q=sum b_i Ω_i, w=sum β_i Ω_i,
c_i=β_i-μ*a_i-ν*b_i, r=Z-P, s=U-Q. With μ,ν,a_i,b_i >= 0,
the JOINT capacity μ*a_i+ν*b_i<=β_i gives β_i,c_i>=0 and
c_i>0 implies β_i>0. Thus the original floors are available wherever
the remaining sum has support. No separate spending of the same β_i
is permitted.

Coordinate linearity gives w=sum c_i Ω_i+μ*P+ν*Q and
W=w+μ*r+ν*s=sum c_i Ω_i+μ*Z+ν*U. On the SAME input x,
the retained floors and the two actual-foot inequalities give

    x.W >= w.w + μ*(h-w.P) + ν*(t-w.Q).

This is derived from weighted primitive floors; the combined pressure
is not assumed. The separately proved full expansion is

    W.W = w.w+2μ*w.r+2ν*w.s+μ²*r.r+2μν*r.s+ν²*s.s.

Combining 2*x.W-W.W <= x.x with this pressure gives exactly

    x.x >= w.w+2μ*(h-w.Z)+2ν*(t-w.U)
                 -μ²*r.r-2μν*r.s-ν²*s.s.

The w.P and w.Q terms cancel algebraically. No flatness is needed;
the mutual residual term r.s remains, including either possible sign.
All sums and floors are explicitly bounded by n, all dot products by d;
empty bounds are allowed. The scalar is an ordered commutative ring,
not an integer replacement for real heights.

The module contains no desired energy or dimension inequality as a
premise, and imports no sorry or new declared axiom. Actual feasible
feet, the common farthest vector and concave reference, supported floors,
capacity estimates and concrete strict numerical comparisons are still
application inputs. The Section5/7 B20 map preserves both nonorthogonal
products. It does not instantiate the geometry or prove B20 in Lean.

Verdict: PASS for source mathematics and statement scope. Actual bounded
compilation and exact axiom stdout must be obtained separately before
any new declaration is counted as kernel checked. At the first admission
observation this turn, 5.5 GiB available memory was below the unchanged
6 GiB start gate; no proof child was started.

## Previously pending finite classification and two source profiles

I read the complete high-output classification162965f7 and the complete
port-sensitive financing proof f75f72f4, including their full finite-table
interfaces in c300. Both pass at their stated written boundaries. Their
separate outside-author reviews are d4da9652 and2e784a0f. Neither supplies
an entire new local rule or a global56 conclusion.

For the classification I checked separation of the unique component
r>=21, its special D21=134 and later D_r=5r+30 sums, and every guard
that could have crossed the separated component. The remaining Y has
rank<=19. The h12..16 per-rank allowances dominate actual small marks;
for example D14/S5 gives exactly4.5 per small rank, D15/S5 gives4,
and rank16 small marks give at most3. Rank11's largest small allowance
is4.9, attained by the A5/D11 replacement after counting its own same
and designated-A17 costs. Its possible E6 cross is paid by E6's unused
baseline allowance. In the h<=10 case, reserving actual A9 saves1.39R
and A10 saves2.1R from the .21RL bound; their respective context
quadratics leave denominators .64 and .2775 in the maximum bounds.
The rank11 completion has denominator .765 as printed. Thus these
cases are not an assumption that all residual forests are connected.

The exact short residual rows preserve every A1 presence restoration,
the D15 mark11 sixteenth and the second A17's independent small cross
charges. In particular Q18=115.25 comes from A17+A1, rather than D17+A1.
For M<=10 I separately checked the at-most-one-intermediate argument
and the all-small concentration up to rank5. The necessary exception
Q6=16 is S5+A1. These checks give the complete strict hard list and
the two maxima371.75 at D23+D17 and D24+D16, with designated A17
included. Every r>=26 row is strictly below363.75.

For the financing proof I reconstructed the full deletion-block bounds
F16<=108.5, F17<=113.5 and the nonleaf-D16 bound95.625, including
the A1/A2/A3 Boolean restorations and designated-A17 cross charges.
For an added small forest of rank R<=3, small same, designated cross
and old ordinary cross cost at most(1+1+2.5)R. The worst held order11
increment is11R/64; guarded rank14..16 cases are separately below5R.
This verifies the complete +5R bound used in the two large-parent cases.
It includes disconnected outputs and the possible additional whole A17.

I counted the non-A high ports directly from the affine paths/forks:
at23 there are at most six, at24 at most eight. The second D16 can
require credit only at its four leaf ports, since its nonleaf bound
makes the recipient finite-safe. The two canonical coefficient values
are29714 and29448; multiplying their V values7344 and7225 gives
kappa159878 and156210 and canonical margins10274 and8871. The
six simultaneous worst debits in the note leave respectively8546,
6818,7171,5471,4519,6823. The minimum4519 is correct.

Finally I checked the actual incidence is into(q, its designated whole
A17), with at most one credit per recipient even if two incidences
arrive. Low-output halves are combined before charging; the A-parent
credit is1/2 once, the other high-parent credit7.5 once. The earlier
unique-A17 high-source families cannot share a new recipient and a
different ideal source. The new source profiles contain neither the
older whole-A15 nor the whole-A16/E8 priority. Thus the canonical
budgets are compatible and not spent twice. The new parent26--40
theorem has its own complete written proof and two full reviews;
it is not inferred from these two earlier selected profiles.
