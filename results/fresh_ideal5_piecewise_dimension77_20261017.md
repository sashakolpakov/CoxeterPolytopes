# Dimension77 from the complete finite9/ideal5 payment rule

**Completed internally checked computer-assisted proof.** Every finite-sided,
finite-volume real hyperbolic Coxeter polytope has
dimension at most77, without arithmetic, simple-cusp or facet-count
assumptions. The argument below combines complete local geometric proofs
with an exact finite scalar comparison; it does not assume the uncompleted
finite8 overlap classification.

The [source](../scripts/fresh_ideal5_piecewise_dimension77_20261017.py) and
[certificate](fresh_ideal5_piecewise_dimension77_20261017.json) bind the
[complete local theorem](fresh_ideal5_combined_payments_20261017.md), its
[571-pattern census and all-rank counts](fresh_ideal5_small_factor_payments_20261017.md),
and the completed [dimension80 proof](fresh_ideal6_piecewise_dimension80_20261017.md)
with its independent review. The local theorem gives finite cost at most9n
for n≤80 and ideal cost at most V times the sum of the factor-rank envelope

    b(r)=r/2 for1≤r≤12;
    b(13)=87/14; b(14)=20/3; b(15)=111/16; b(16)=127/17;
    b(r)=163/23 forr≥17.

The key small-factor change is the selective order14 rule. It retains87
pairs in affine B,84 in C,82 in D and70 in the cycle. These are complete
pair sets with every actual preferred-component embedding included, except
where a separate universal or sealed finite-coupon transfer replaces it.
The combined proof checks that no elementary triangle payment or immediate
small-affine payment is lost by these removals. The context-specific D
removals at orders16 and17 have distances greater than9, so they are used
only after the complete no-short9 classification; they do not create a
circular short-path reduction.

## Complete reduction of cusp-link shapes

For the actual link L=product_i Delta^(r_i), put d=sum_i r_i=n−1 and

    H(t)=product_i(1+t+...+t^r_i)=sum_j h_j*t^j,
    V=product_i(r_i+1),       G=sum_j h_j*min(j,d−j).

Its normalized curvature is (5G−dV)/(2V). The inherited
[convex-order theorem](fresh_distance14_unrestricted_dimension156_20261017.md)
shows that a transfer(a,b)→(a−1,b+1), for1≤a≤b, can only decrease this
curvature while all other ranks remain fixed.

The new envelope drops below the former linear cost at rank13, so it would
be incorrect to concentrate all low parts up to13. Instead concentrate
two parts strictly below12 until one vanishes or the other reaches12;
their total cost is unchanged. Replace two13 parts by12 and14, increasing
cost by5/21. Replace two15 parts by14 and16, increasing cost by107/408.
Finally concentrate two parts strictly above17 until one reaches17;
their constant cost is unchanged. Every operation can only decrease
curvature minus cost.

Each operation strictly increases the integer sum of squared ranks at
fixed total d, and this sum is bounded by d². Thus the process terminates
even though several kinds of moves may be interleaved. Its terminal shapes
have arbitrary copies of12,14,16,17; at most one13 and at most one15;
at most one low part1–11; and at most one high part at least18. The source
enumerates every such shape, so its finite comparison covers every integer
partition of each relevant link dimension.

## Exact comparison and the zero finite margin

The primary checker computes H,G,V by integer convolution and checks the
same curvature through a separate cumulative inclusion-exclusion moment
identity. The exact reduced-shape counts and minima are:

| Ambient dimension | Profiles | Minimum normalized ideal margin | Finite slack beta_n−9n |
|---:|---:|---:|---:|
|78|1173|409611/2171546|0|
|79|1239|63437/167042|21/2|
|80|1310|3585185/7619196|20|

Here beta_n=n(n−6)/8 for even n and (n−1)(n−5)/8 for odd n.
The overall minimum occurs at ambient78 and ranks(1,12,16,16,16,16).
All3,722 ideal margins are strictly positive.

The inherited [link-shape inequality](fresh_cusp_link_shape_curvature_20261017.md)
is

    D5 ≥ beta_n*v_f + sum_cusps kappa(L) + 5(n−1).

The face-payment theorem gives an upper bound9n*v_f plus the ideal
budgets. In dimension78 the finite coefficient is exactly equal to9n,
so it must not be described as a strictly positive finite margin. The
global term5(n−1)=385 is strictly positive even when there are no cusps;
when cusps are present their margins are positive too. Dimensions79 and80
also have positive finite margins. The inequalities therefore contradict
the existence of a polytope in every dimension78,79,80. The completed80
theorem excludes all larger dimensions. Together they give n≤77.

At ambient77, the formal profile(12,16,16,16,16) has margin
−16664/1085773. This comparison does not exclude77. The profile is not a
constructed hyperbolic Coxeter polytope or an assertion of sharpness.

The former auxiliary compact-facet argument still has its proved bound79;
the unrestricted77 theorem now gives the stronger effective bound77 for
that subclass. The separate compact dimension29 and compact-arithmetic
dimension21 results are unchanged. No new finite-volume example is claimed.

## Verification scope

The separate [independent review](fresh_ideal5_piecewise_dimension77_20261017_review.md)
records its reconstruction and mathematical scope. Primary full generation
and a fresh separate-process saved replay both terminated exit0/PASS;
syntax and all19 direct bindings also pass. The separate standard-library
checker independently reconstructed and compared all3,722 rows, the three
minima including their global slack, the boundary example, all15 affine
branches and all16 small target pair sets with their triangle inclusions.
Its final bound generation and saved replay are recorded in that separate
review. A full source/mathematical read found no gap in the complete local
composition, cost directions, terminating concentration, actual cusp
multiplicities, zero finite margin at78 or inherited80 tail. The complete
local binder also passed generation, fresh replay and separate full
mathematical/source reads. Source, data and this report are final/frozen.

The proof imports the published finite/affine classification, the complete
old minimal-core and overlap proofs, and the explicitly linked local
geometry. Their arithmetic and source-review scopes are recorded in their
own frozen packages. A scalar pass alone would not prove the theorem; the
complete local composition is an essential input. All reviews here are
internal. External mathematical refereeing, publication priority and
proof-assistant formalization remain unestablished.
