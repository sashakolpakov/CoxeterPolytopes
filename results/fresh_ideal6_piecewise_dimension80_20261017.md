# Dimension80 from the finite9/ideal6 payment rule

**Completed internally checked computer-assisted proof.** Every finite-sided,
finite-volume real hyperbolic Coxeter polytope has dimension at most80,
without an arithmetic, simple-cusp or facet-count hypothesis. The local
geometric arguments, their complete composition and the scalar comparison
have passed the internal checks specified below. External refereeing and
proof-assistant formalization remain separate, outstanding checks.

The [source](../scripts/fresh_ideal6_piecewise_dimension80_20261017.py) and
[certificate](fresh_ideal6_piecewise_dimension80_20261017.json) bind the
[combined local payment theorem](fresh_ideal6_combined_payments_20261017.md),
the [completed86 theorem](fresh_order16_piecewise_dimension86_20261017.md)
and its independent review. The new local input has finite cost at most9n
for n<=108 and ideal cost V*sum_i b(r_i), with

    b(r)=r/2 for1<=r<=13;
    b(14)=20/3; b(15)=111/16; b(16)=127/17;
    b(r)=95/13 forr>=17.

The order15 improvement rests on a universal transfer of nine long affine
D14 pairs to finite half-payments. It is not an assumption that the63 saved
small-core contexts have no further supporting whites. The
[broad finite budget](fresh_broad_finite_half_coupons_20261017.md) covers
every finite coupon in the combined proof.

## Reduction of all cusp-link shapes

For a link L=product_i Delta^(r_i), put d=sum_i r_i=n-1,

    H(t)=product_i(1+t+...+t^r_i)=sum_j h_j*t^j,
    V=product_i(r_i+1),       G=sum_j h_j*min(j,d-j).

The normalized curvature is kappa(L)/V=(5G-dV)/(2V). The inherited
[convex-order argument](fresh_distance14_unrestricted_dimension156_20261017.md)
proves that transferring(a,b) to(a-1,b+1) for1<=a<=b can only decrease
this curvature, with other factors fixed.

Concentrate two ranks below13 until one vanishes or the other reaches13.
Their linear cost is preserved. Replace two14 parts by13 and15; this
increases cost by5/48. Replace two15 parts by14 and16; this increases cost
by107/408. Concentrate two ranks above17 until one reaches17; their constant
cost is preserved. In each case curvature minus cost cannot increase.

The middle two moves can create each other's exceptional ranks, so a single
pass is not a termination proof. Instead every move strictly increases
the integer sum of squared ranks at fixed total d, which is bounded by d^2.
After all available moves stop, there are arbitrary copies of13,16,17,
at most one low part1,...,12, at most one14, at most one15, and at most
one high part at least18. Checking all these reduced shapes therefore
covers every integer partition of d.

## Exact comparison

The primary checker computes H,G,V by exact integer convolution and checks
the equivalent absolute-deviation identity. The independent checker uses
a different enumeration, sliding-window convolution and a cumulative
inclusion-exclusion identity. All5,008 reduced shapes in dimensions81--86
agree and have positive curvature minus the ideal cost.

| Ambient dimension | Profiles | Minimum normalized margin | Finite slack beta_n-9n |
|---:|---:|---:|---:|
|81|742|940675/11294304|31|
|82|776|1132393/9143008|41|
|83|813|3418145/9143008|105/2|
|84|852|15382621/27429024|63|
|85|892|221295/326536|75|
|86|933|12970945/15059072|86|

Here beta_n=n(n-6)/8 for even n and(n-1)(n-5)/8 for odd n. The overall
minimum occurs in ambient81 at ranks(2,13,13,13,13,13,13).

The inherited [shape inequality](fresh_cusp_link_shape_curvature_20261017.md)
gives

    D5 >= beta_n*v_f + sum_cusps kappa(L) + 5(n-1).

The actual face payments give the upper bound9n*v_f plus the ideal budgets.
The strictly positive finite and ideal margins contradict these inequalities
for every n=81,...,86. The completed86 theorem excludes all n>=87.
Together these two ranges give n<=80.
The compact-facet ambient ceiling79 and compact/arithmetic bounds are unchanged.

At ambient80, the formal profile(1,13,13,13,13,13,13) has margin
-1131729/7529536. This particular comparison therefore does not exclude80.
It is a formal link-shape obstruction, not a constructed polytope or a claim
that the bound is sharp. Lowering order14 costs or a different global
allocation is the next required improvement.

## Verification scope

The [independent review](fresh_ideal6_piecewise_dimension80_20261017_review.md)
states its arithmetic and composition scope. Primary generation and fresh
saved replay both terminated exit0/PASS, as did generation and replay of
the combined local-payment binder. A separate agent independently compared
every field of all5,008 scalar profiles, the six dimension minima and the
boundary profile, and independently derived the12 all-rank affine branches.
The complete source and written composition received a separate mathematical
read, including the D14 profiles and sealing, the short-path/no-short-path
split, termination of concentration, actual cusp multiplicities and the86
tail. These checks identify no remaining gap in the stated argument.

The independent scalar checker does not duplicate the inherited minimal-core
classification or local geometric censuses. Those remain explicitly linked
proof inputs with their own certificates and review scopes. Final version
bindings and saved replays are recorded in the independent review and
the dimension80 checkpoint archive. All reviews here are internal; no
external mathematical refereeing, publication priority or proof-assistant
formalization is asserted.

    python3 scripts/fresh_ideal6_piecewise_dimension80_20261017.py
    python3 scripts/fresh_ideal6_piecewise_dimension80_20261017.py --verify
