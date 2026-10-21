# Exact ambient60 payment frontier under the completed dimension60 rule

The current rule does not close the ambient60 comparison. Its complete finite
rank60 maximum is **413**, against curvature405. Among all2,113 terminal cusp
profiles of rank59,77 have negative margins and none has zero margin. The
smallest normalized margin is−1823/11560 at(9,9,9,16,16). The largest absolute
terminal deficit is192250 at(9,9,9,9,9,14).

This is an exact obstruction diagnostic for the
[completed dimension60 rule](fresh_order11_order16_dimension60_20261018.md),
not a new local payment theorem or a dimension59 claim. Formal forests and
cusp products in this calculation are not asserted to occur in an actual
hyperbolic polytope. Negative terminal absolute margins are not extended
to arbitrary nonterminal cusps.

## Current finite rule and complete rank60 count

The [finite420 certificate](fresh_rank61_complete_finite420_budget_20261018.md)
provides the exact current component and cross-component rules. The final
whole A1/B10 doubled-terminal coupon is a unit, not the predecessor half:
the source explicitly adds its doubled increment1 and rechecks simultaneous
two-end domination for all33,489 typed component pairs. All four earlier
coupon increments and the D17 restoration are retained. The current exceptional
component totals are

|Component|A15|B15|D15|A16|B16|D16|A17|B17|D17|
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|Total|103|103|103|111|104|111|115|105|119|

All other base weights, typed coupons, exceptional/I2 domination and the
depth8 rule are inherited from the bound finite sources. The D16 upper
representative now has its actual cost111. Some other representative coupon
allowances remain conservative; the saved records distinguish the upper
cost from the actual cost of each reconstructed forest.

Every finite forest is dominated component by component, on both ends of
every coupon, by the following representatives. Rank6 uses E6, rank9 uses B9,
ranks10..16 use D of the same rank, and A17 remains separate because of its
cross-component roles. All other ranks at most8 use A; all other ranks at
least17 use D. Whole A17 and D17 are thus different species throughout.
The source checks base and simultaneous cross domination on the complete
finite-type list through rank60, including the exceptional types and the
arbitrary-dihedral representative.

There are1,044 possible multiplicity vectors for the nine fixed species
B9,D10,...,D16,A17 with total rank at most60. They produce792 interaction
keys. The remaining species have ranks1..8 or17..60. There are no
cross-payments among these residual representatives, or among the fixed
ones. For each residual rank the saved key gives its base cost plus every
interaction with the fixed forest. The54,288 encoded item costs are compared
directly with the current typed rule.

Generation computes new item-first unbounded knapsack tables through rank60
for all792 keys. The separate replay does not repeat that optimization: it
checks every possible last item as an upper inequality and checks one
attaining last item at every nonzero rank. These1,142,064 inequalities and
47,520 predecessors prove both the upper bound and attainment inductively.
Every residual optimum is then combined with its entire fixed forest.
The newly computed table values also agree with all792 appropriate prefixes
of the frozen rank61 tables; they were not copied from those tables during
generation.

The maximum is413. Actual formal forests attaining it include

- D60;
- D17+D43;
- 2D17+D26.

The finite curvature is

\[
G_f=\sum_{j=0}^{59}\min(j,59-j)=870,
\qquad \beta_{60}=\frac{5G_f-59\cdot60}{2}=405.
\]

All17 fixed-vector alternatives whose upper cost exceeds405 also have an
actual reconstructed witness above405. The certificate includes every such
record. Particularly useful necessary obstructions are:

|Actual forest|Cost|Deficit over405|
|---|---:|---:|
|D60 or2D17+D26|413|8|
|4D15|412|7|
|D16+D17+D27|412|7|
|D15+D17+D28|411|6|
|2D16+D28|411|6|
|D15+D16+D29|410|5|
|A2+2A5+3D16|408|3|

The4D15 witness may use any four A15/B15/D15 components: their base costs
agree and they have no mutual coupons. Thus a change only to D15 would not
resolve that rank15 obstruction. If only the rank15 total changes, it needs
at least7/4 saving per component. Other displayed witnesses require separate
savings, and new coupon costs must be counted before drawing any sufficient
finite conclusion.

One proposed additional rule would delete distance7 pairs whose two marks
have feature depth at least11, with counts max(0,n−29), max(0,n−30),
max(0,n−31) in A_n,B_n,D_n respectively. This packet does not prove that
geometric rule. Even if it holds with all other costs unchanged, its exact
rank60 maximum is still413: monotonicity gives the old upper413, while
2D17+D26 has no such pair and still attains413. This conditional statement
requires no new optimization. In particular, removing the singleton D60
obstruction is insufficient; medium component ranks remain binding.

## Complete ideal domain and exact moments

For a rank-r affine factor, its normalized cost is its total pair payment
divided by r+1. The current costs are

\[
c_r=r/2\quad(1\le r\le9),\qquad
(c_{10},...,c_{16})=
\left(\frac{53}{11},\frac{63}{12},\frac{73}{13},\frac{83}{14},
\frac{94}{15},\frac{105}{16},\frac{116}{17}\right).
\]

The source reads the completed proof's actual order11 and16 family totals,
checks all2,365 prior normalized cost rows, and compares every tail family
formula with the frozen
[ambient61 frontier](fresh_ambient61_payment_frontier_20261018.md).
For order m at least18 the exact inherited family totals are

\[
\begin{aligned}
A_m&=5m,\\
B_m&=5m-6+\sum_{k=8}^{11}\min(m-k,13),\\
C_m&=5m-15+\sum_{k=6}^{9}\min(m-k,14),\\
D_m&=5m+3+\sum_{k=10}^{13}\min(m-k,12).
\end{aligned}
\]

The maximum of these totals divided by m is used separately at each rank;
there is no replacement by a constant tail envelope.

Only the inherited low9 concentration is used. Since all costs through
rank9 are still r/2, it preserves total normalized cost and decreases the
normalized curvature term until there is at most one factor of rank1..8;
all other factors have rank at least9. No high-rank concentration or
restriction on repeated large factors is made. Sorted recursion enumerates
every partition of59 in this domain, giving2,113 profiles.

For each profile(r_i), generation and replay reconstruct

\[
H(t)=\prod_i(1+t+\cdots+t^{r_i}),\qquad
V=\prod_i(r_i+1),\qquad
G=\sum_{j=0}^{59}[t^j]H(t)\min(j,59-j).
\]

The normalized margin and terminal absolute margin are

\[
M=\frac{5G-59V}{2V}-\sum_i c_{r_i},\qquad VM.
\]

All coefficients are exact integers. Symmetry, length, total coefficient
sum, the absolute-deviation identity and a cumulative-distribution identity
are checked independently of the convolution accumulation. For this odd
link dimension the CDF calculation centers at30; symmetry makes its mean
absolute deviation equal to that centered at59/2. Thus it gives the same
G/V=59/2−deviation. The odd-dimensional constant is handled explicitly.

There are77 negative profiles, no zero profiles. The full list is saved.
Some leading normalized deficits are:

|Affine ranks|Normalized margin|
|---|---:|
|(9,9,9,16,16)|−1823/11560|
|(9,9,12,14,15)|−2311/15600|
|(9,9,11,15,15)|−1117/7680|
|(9,9,11,14,16)|−8251/61200|
|(9,9,9,9,9,14)|−769/6000|

The singleton rank59 profile is positive. Its order60 family costs are
A300,B346,C341,D351, so its normalized margin is9/10 and its absolute
surplus is54. The largest terminal absolute deficit,−192250, occurs at
(9,9,9,9,9,14), not at the normalized minimizer.

As in the inherited frontier, negative terminal absolute margins are not
universal absolute lower bounds: the monotonic product-V argument requires
positive margins before it can be applied in that direction.

## Exact conditional local targets

These scenarios change only factor envelopes, leave ranks at most9 unchanged,
and assert no new local theorem or finite coupon allowance. A total saving
delta at affine rank r raises a profile's normalized margin by
count(r)*delta/(r+1). Every saved scenario is checked on all2,113 profiles,
including the currently positive ones.

Some individual orders need a reduction independently of all the others.
While ranks at most9 are unchanged, the following profiles contain just one
affected order and force the stated strict saving at that order regardless
of improvements at other ranks at least10:

|Affine order|Profile|Strict necessary total saving|
|---|---|---:|
|13|(8,9,9,9,12,12)|>1559/15600|
|15|(9,9,9,9,9,14)|>769/400|
|16|(8,9,9,9,9,15)|>3361/3000|
|17|(9,9,9,16,16)|>1823/1360|

Reducing every order11..17 total by1 still leaves three negative profiles.
Reducing every such total by2 makes all profiles positive. A narrower
sufficient rational scenario, leaving order11 unchanged, is:

|Affine order|Current maximum|Conditional new maximum|Saving|
|---|---:|---:|---:|
|12|63|125/2|1/2|
|13|73|145/2|1/2|
|14|83|331/4|1/4|
|15|94|92|2|
|16|105|415/4|5/4|
|17|116|229/2|3/2|

This gives normalized minimum73/104000 at(9,9,9,9,11,12) and terminal
absolute minimum54 at(59). The certificate includes the exact strict
threshold for each coordinate while all other entries in this scenario
are fixed. These are sufficient all-family envelope targets, not a claim
that any particular Coxeter family can attain the savings for free. This
scenario was selected by a small numerical linear-program scout followed
by rational rounding; only the displayed exact row comparisons are used
in the certificate, with no optimality claim about the scenario.

Changing order10 would invalidate the particular cost-preserving low9
reduction used here and would require a new terminal-domain argument or
an expanded domain. This packet makes no such change.

## Global status and verification

At ambient60 the strengthened two-extrema correction is4500, because each
link contributes at least
sum(2j−59 for2j>59)=900 and the two extremal coefficients total5. It is
counted once. Current finite deficit8 and the negative ideal profiles keep
the standard pointwise comparison open; they do not contradict the
[completed dimension60 proof](fresh_order11_order16_dimension60_20261018.md),
which excludes ambient61 and higher. No new finite-volume example is
asserted.

- [Producer](../scripts/fresh_ambient60_payment_frontier_20261018.py).
- [Exact certificate](fresh_ambient60_payment_frontier_20261018.json).

Final generation session22668, chunk5b0d55, and separate fresh replay
session19716, chunkbc9b8f, both terminated exit0/PASS. The fresh replay
checks saved finite Bellman certificates without reoptimizing, while
reconstructing every ideal convolution and scalar comparison. It does
not repeat any geometric census, local payment theorem or physical forest
realization. Direct metadata checks pass13 source bindings, six local
links and Python syntax.

The root agent completed a full source/report mathematical read: PASS.
It checked the current A1/B10 increment and simultaneous typed domination,
the complete fixed/residual partition and zero residual interactions,
item identities and last-item induction, actual413 attainers and the
conditional unchanged2D17+D26 witness, the odd-link CDF/symmetry identity,
low9-only domain and exact tail formulas, conditional targets and the
restriction on negative absolute margins, and the single constant4500.
That review did not replay DP, convolution or any geometric census.

A second agent completed an independent full producer/report mathematical
read, including the final residual-pair zero-interaction guard: PASS. It
checked simultaneous typed domination and the actual A1/B10 upgrade,
the complete1,044-vector/792-key decomposition and last-item certificates,
the actual413 witnesses against405, the complete low9-only2,113-profile
domain with exact family tails and odd-link moment identity, and the
conditional-only targets and thresholds. Saved-summary inspection confirmed
the17 offending finite alternatives,77 negative ideal profiles and unchanged
413 in the hypothetical deep7-only scenario. This review did not repeat DP,
convolution, scalar arithmetic, graph-distance counts or any geometric census.
Source, certificate and report are final/frozen.

    python3 scripts/fresh_ambient60_payment_frontier_20261018.py
    python3 scripts/fresh_ambient60_payment_frontier_20261018.py --verify
