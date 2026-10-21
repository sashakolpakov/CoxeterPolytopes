# Independent review of the dimension80 comparison

**FINAL/FROZEN: independent arithmetic, written/source composition review, final hash-bound generation and fresh saved replay PASS.** The [independent checker](../scripts/fresh_ideal6_piecewise_dimension80_20261017_review.py) imports no producer or project arithmetic. It reconstructs all5,008 reduced profiles using standard-library integer/Fraction arithmetic and compares every saved scalar field with the [primary certificate](fresh_ideal6_piecewise_dimension80_20261017.json). The complete [primary proof](fresh_ideal6_piecewise_dimension80_20261017.md) and [local composition](fresh_ideal6_combined_payments_20261017.md) have been read; no mathematical gap was identified within their stated inherited inputs.

## Independent arithmetic

The proposed normalized cusp-factor budget is

    b(r)=r/2 for1<=r<=13;
    b(14)=20/3; b(15)=111/16; b(16)=127/17;
    b(r)=95/13 forr>=17.

For each link dimension d=80,...,85, the independent enumeration chooses the optional high part in{0,18,...,d}, optional low part in{0,...,12}, and the presence or absence of14 and15. It then chooses the multiplicities of17 and16. The remaining nonnegative rank must be divisible by13, forcing its multiplicity. Zero parts are omitted and sorted shapes are deduplicated. This exhausts the stated reduced domain without using a producer list.

Sliding-window convolution independently computes

    H(t)=product_i(1+t+...+t^r_i)=sum_j h_j t^j,
    V=product_i(r_i+1),  G=sum_j h_j min(j,d−j).

The code checks polynomial length, symmetry and H(1)=V. It separately expands product_i(1−t^(r_i+1))=sum_l a_l t^l, combining repeated shifts. For s factors and q=floor(d/2), the cumulative-tail identity gives

    G=qV−2 sum_l a_l binom(q−l+s,s+1).

This agrees exactly with the coefficient computation for every profile and the boundary profile. The absolute-deviation identity is also checked:

    kappa/V=(5G−dV)/(2V)
           =3d/4−(5/2)E|sum_i U_i−d/2|.

| Ambient dimension | Profiles | Minimum curvature minus cost | Finite slack beta_n−9n |
|---:|---:|---:|---:|
|81|742|940675/11294304|31|
|82|776|1132393/9143008|41|
|83|813|3418145/9143008|105/2|
|84|852|15382621/27429024|63|
|85|892|221295/326536|75|
|86|933|12970945/15059072|86|

The minimum is940675/11294304 at ranks(2,13,13,13,13,13,13). The ambient80 boundary profile(1,13,13,13,13,13,13) has margin−1131729/7529536. This boundary is a limitation of the particular scalar comparison, not a realized chamber or a sharpness claim.

## Completeness of the shape reduction

The inherited convex-order theorem says that transferring(a,b) to(a−1,b+1) for1<=a<=b can only decrease normalized curvature, with the other uniform factors fixed.

Two ranks at most12 can be concentrated until one vanishes or the other becomes13; their linear budget is preserved. Two ranks14 can be replaced by13 and15, increasing the budget by5/48. Two ranks15 can be replaced by14 and16, increasing it by107/408. Two ranks at least18 can be concentrated until one becomes17, preserving the constant budget. In every case curvature minus budget cannot increase.

The two middle-rank moves can create each other's exceptional rank, so they must be allowed to repeat rather than applied only once in a fixed order. Termination still follows: every elementary transfer strictly increases the integer sum of squared ranks at fixed total d. The process is bounded and ends with at most one14 and one15. Low transfers stop at13, high transfers stop at17, and the middle moves produce only13,14,15,16. The final partition therefore consists of arbitrary copies of13,16,17, at most one low rank1,...,12, at most one14, at most one15, and at most one high rank≥18. Every possible terminal shape is enumerated above.

The helpful cost increases are checked exactly, not treated as cost preservation. Positivity at a reduced shape implies positivity at the original shape.

## Local finite-budget review

A separate narrow source/mathematical read of the broad finite half-coupon budget passes. Let a be the total rank of whole finite components of rank at most4, and let t count whole components D12,...,D16,A13,A15. The small-component base slack is at least15a/2. Each listed large component has rank at least12, slack at least24 and at most three eligible marked nodes. The total cross-component half-payment is at most3at/2, counted once per eligible pair. Thus

    cost<=9n+a*(3t/2−15/2)−24t,    a+12t<=n.

For n<=108, t<=9. The coefficient of a is nonpositive for t<=5. For t=6,7,8,9 the maximal excesses are−90,−96,−138,−216. This proves the finite budget independently of the small components' particular finite types. It is a counting theorem only: the separate local geometry must ensure that a proposed coupon refers to whole actual components.

## Independent all-rank ideal cost calculation

For factor order m>=18, a direct spine-interval count gives the following costs for the ideal distance6/feature6 rule:

    A: 6m;
    B: 6m−14 + sum_{s=9,10,11} min(m−s,15);
    C: 6m−21 + sum_{s=7,8,9} min(m−s,16);
    D: 6m−7 + sum_{s=11,12,13} min(m−s,14).

For D, the distance-at-most6 contribution is6m−19; each distance7,8,9 adds four leaf-to-main pairs and a main-spine window of length min(m−4−distance,14). For B, the corresponding short contribution is6m−20, with two leaf-to-main pairs per added distance and main windows of combined width15. For C the short contribution is6m−21 and the two end windows have combined width16. Opposite-end leaf pairs do not enter these distances at m>=18.

The independent checker constructs the twelve linear branches from their window breakpoints and compares all branch coefficients with the local certificate. On each branch cost/m is a+b/m, so its maximum occurs at an endpoint or at the limiting value6. The exact normalized maxima are36/5 for B at m=25,85/12 for C at m=24, and95/13 for D at m=26. This proves the all-rank95/13 envelope without extrapolating a finite graph sample. The exceptional orders15,16,17 use the explicitly bound small-factor census, with normalized maxima20/3,111/16,127/17. Whole factors of order at most14 cost binom(m,2)/m=(m−1)/2.

## Local geometry and composition review

The complete source and written proof of the [universal D14 long-pair transfer](fresh_affine_d14_long_pair_transfer_20261017.md) were read separately. The large finite-D extension profiles are exhaustive in the stated endpoint domains. At a fork center, the arbitrary outer contacts obey the strict edge-PD inequality u²+v²<4. Nonextreme end profiles have strictly positive determinant margin. In the remaining cases both heights are1/2 before the spectator-gap argument is applied. Its equality case permits exactly one unit contact to an isolated spectator; all other supporting whites are orthogonal to the resulting core. Thus the finite coupons refer to whole actual components. The72 literal equality cases and their corner/coupon reconstruction remain producer arithmetic inputs, not a second census by this reviewer.

The combined payment theorem retains the complete triangle payment: removed D pairs are too distant to remove the short or feature witnesses. For a quadrilateral with some corner path of length at most9, finite corners and small ideal factors pay directly. The mixed6 tree and cycle replacements cover the remaining distance7 cases and tree distances7--9. For a large cycle of distance8 or9, the preceding mixed7 transfer gives finite9 or ideal7 payment, and the latter is covered by mixed6 if necessary. The order16 antipodal cycle has its separate finite A15-center coupons. These implications are applied under the assumption that the entire face is unpaid, so no iterative termination assumption is needed.

If no corner has a path of length at most9, the hypothesis of the inherited complete preferred-component theorem is exactly satisfied. Its low and high cases give components of order at most12 and17 respectively. Finite components pay directly; ideal components of order at most14 pay in full. Every induced marked image at orders15--17 is included in the small-factor census, except the explicitly removed nine D pairs at each order. The universal D14 transfer and the frozen order16/order17 transfers pay those exceptions. The combined binder matches the nine removed order15 pairs to the D14 marked-pair domain. The broad finite budget includes all finite coupons used in these transfers and is valid throughout the needed dimensions81--86.

At an actual nonsimple cusp, a pair in a rank-r simplex factor has multiplicity V/(r+1). The inherited actual-factor theorem therefore gives the ideal budget V*sum b(r), with no simplicity, arithmetic or minimum-factor-rank restriction. The inherited shape inequality has lower bound beta_n*v_f+sum kappa(L)+5(n−1). Every finite and ideal coefficient exceeds its payment cost in dimensions81--86 by the exact positive margins above. The completed [dimension86 theorem](fresh_order16_piecewise_dimension86_20261017.md) supplies the tail n>=87. The resulting composition proves the stated unrestricted bound80; the compact-facet ambient bound79 is unchanged.

## Verification scope

All5,008 primary scalar rows, six minima and the boundary row agree exactly. Polynomial coefficients, cumulative-tail curvature and absolute-deviation curvature agree for all5,009 tested shapes. All twelve independently reconstructed affine branches agree with the local binder. The complete shape reduction and primary/local mathematical composition have received the source reads described above. The parent reviewer separately checked the D14 local proof and this independent scalar source.

The local geometric classifications, the complete inherited preferred-component census, and the actual cusp-factor/shape theorems remain explicit dependencies. This review does not duplicate their full geometric or embedding censuses. Final generation and a separate fresh saved replay both terminated with exit0/PASS after the primary report was frozen. All31 direct source hashes, source syntax and local report links were checked. No external refereeing, publication-priority claim or proof-assistant formalization is asserted.

The [saved review](fresh_ideal6_piecewise_dimension80_20261017_review.json) binds primary source SHA256 `7c160770b5c8503eda48ae723d83ef0b7c6c540f37acf9f554d7249152d96567`, primary JSON `0779c84a01fd5d9a4179db3ffcbad0dcd3b9744b6f8c2f5c13487d122bf5c133`, and primary report `918c37a7f8bcf1eae66871dcd0602698fd8095b2958b3931fa2a4037ee98e95f`.

    python3 scripts/fresh_ideal6_piecewise_dimension80_20261017_review.py
    python3 scripts/fresh_ideal6_piecewise_dimension80_20261017_review.py --verify
