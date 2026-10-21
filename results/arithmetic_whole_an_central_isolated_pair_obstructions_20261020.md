# Central isolated caps at a whole A_n vertex

20 October 2026. Own paper proof, prepared for independent review.
No numerical optimization, enumeration, import, Lean or TeX execution.
These are necessary pair restrictions, not a dimension bound.

Let P be a compact arithmetic Coxeter polytope with an actual whole
A_n vertex p. Normalize its old roots α_i to norm two, put N=n+1,
and write

    d_i=i(N−i)/N,    C_ij=i(N−j)/N for i<j.

Consider two distinct actual opposite walls, each having just one
nonzero old contact, a finite native label at its own contact index
i or j. Throughout, **central** means precisely

    i≥5,   N−j≥5,   i<j.                              (1)

It does not mean the larger interval defined only by d_i>2. Put
Δ=j−i, r=N−j, C=ir/N, and d_Δ=Δ(N−Δ)/N. Then

    d_i+d_j=2C+d_Δ,       d_i d_j=C²+CΔ,
    C=ir/(i+r+Δ) ≥25/(10+Δ).                          (2)

The last bound follows by increasing either i or r from its minimum 5.

## 1. Exact normalization and common field

The isolated finite normalization and cyclotomic restriction are derived
in
[arithmetic_full_an_isolated_finite_continuation_packing_20261020.md](arithmetic_full_an_isolated_finite_continuation_packing_20261020.md),
SHA256 `9ab34f073ea95bed2040bafe6b8573331721f9b1e75d65172f814637da046810`.
For clarity, the pieces needed here work in every dimension admitting
(1), not merely the n≥13 convention of that predecessor: (1) implies
d_i,d_j>2 and at least three old nodes on both sides of each contact.

For x=4cos²(π/m), the root is scaled OVER the ground field to

    c_i²=2x,   c_i.α_i=−x,   c_i.α_k=0 (k≠i).

The full-rank old-star Gram is Lorentzian over Q(x). A proper ground
field extension would have a nonphysical embedding fixing this Gram,
so the actual ground field equals Q(x). Conjugate positive definiteness
gives 0<x'<2/d_i<1 for every nonphysical conjugate.

The predecessor's elementary coprime-exponent argument therefore leaves
only labels3,4,5,6,8,12. Labels3,4,6 contain positive-null affine
subdiagrams: respectively arms (1,3,3), the 3/3/4 star, and the 6−3
path. These are present under (1) and impossible in a compact chamber.
Thus the common actual field and label are exactly one of

| label | x | x²=sx−ν | nonphysical x' |
|---|---|---|---|
| 5 | φ²=φ+1 | 3x−1 | φ^(-2) |
| 8 | 2+√2 | 4x−2 | 2−√2 |
| 12 | 2+√3 | 4x−1 | 2−√3 |

Different rows cannot coexist, since their full-rank stars fix distinct
ground fields. In all three rows {1,x} is an integral basis.

Writing the root as c_i=a_i p−xω_i, with ω_i the old C-dual weight,
gives a_i²=x(xd_i−2) and a_i>0. The common timelike residual therefore
fixes the normalized norm-two pair product

    R=(c_i.c_j)/x
     =xC−sqrt[(xd_i−2)(xd_j−2)] ∈ Q(x).               (3)

The positive square root is essential: both walls face the same p.

## 2. Every such central pair is ultraparallel

If two actual facets meet, acuteness and their Coxeter ridge angle give
R=−2cos(π/ell). Its membership in the indicated quadratic field leaves
exactly these finite possibilities:

    label5:  0,−1,−φ;
    label8:  0,−1,−√2;
    label12: 0,−1,−√3.                               (4)

This is an ACTUAL-facet native-angle statement. Arbitrary intersecting
arrangement mirrors can have other angles and must not be substituted.
Positive definiteness of their pair Gram implies an actual intersection
by positive acute projection, so there is no omitted finite pair whose
intersection is merely outside P.

Squaring (3), and using (2), gives the necessary exact identity

    C[x²Δ+2x(R−2)]−2x d_Δ+4−R²=0.                   (5)

Comparing rational and x coefficients gives the complete table:

| label | finite R | forced C | forced d_Δ when needed |
|---|---|---|---|
| 5 | 0 | 4/Δ | 6−8/Δ |
| 5 | −1 | 3/Δ | 9/2−9/Δ |
| 5 | −φ=1−x | 4/(Δ−2) | (11Δ−30)/(2(Δ−2)) |
| 8 | 0 | 2/Δ | 4−4/Δ |
| 8 | −1 | 3/(2Δ) | 3−9/(2Δ) |
| 8 | −√2=2−x | 1/(Δ−2) | 2 |
| 12 | 0 | 4/Δ | 8−8/Δ |
| 12 | −1 | 3/Δ | 6−9/Δ |
| 12 | −√3=2−x | 1/(Δ−2) | 2 |

No division at Δ=2 is intended: in those rows the original constant
equation is C(Δ−2)=1 or 4, already impossible for Δ≤2.

For C=a/Δ, the bound (2) leaves only Δ=1 when a is 3 or 4, and
no positive integer Δ when a is 2 or 3/2. At Δ=1 the four surviving
d_Δ entries are negative or zero, whereas d_1>0. For C=1/(Δ−2),
(2) would require Δ≤5/2, contrary to Δ>2.

Only label5, R=−φ remains. Its bound gives 2<Δ≤30/7, so Δ=3 or 4.
For Δ=3 the table gives d_3=3/2 and hence N=6, inconsistent with
i+r≥10. For Δ=4 it gives N=32 and C=2; now i+r=28, i,r≥5 implies
ir≥5·23=115, whereas C=2 would require ir=64. This is impossible.

Thus every finite pair case is excluded. A parallel pair R=−2 would
give a positive future-null combination of the two outward cap roots,
feasible against every other wall by acuteness; compactness forbids it.
Consequently

    R<−2: every pair in (1) is ultraparallel.           (6)

This proof uses no three-spacing assumption and no rank table.

## 3. An arithmetic denominator-spacing obstruction

The normalized product R is an algebraic integer. Indeed its pair
cyclic product is R²∈O_K; since R∈K, the monic equation T²−R² makes
R integral. Write R=A+Bx with A,B∈Z, and set T=Δ+2B.
The rational coefficient of (5) is exactly

    A²−ν B²+ν C T=4.                                (7)

Hence ν C T is an integer, and

    N divides ν i j T.                               (8)

This retains the factor ν=2 in the label8 case. It cannot be silently
canceled at dyadic divisors of N.

At the other embedding, positive definiteness gives |R'|<2. Since
R<−2, B<0. The arithmetic-geometric mean in (3) gives

    R≥2−x d_Δ/2.

As x/(x−x')<2 in every row, this implies

    −Δ < T=Δ+2B < Δ.                                 (9)

Put M=N/gcd(νij,N). If M≥Δ, equations (8) and (9) force T=0.
Then Δ=2k, B=−k for an integer k≥1, and (7) becomes

    A²−ν k²=4.                                      (10)

For labels5/12, ν=1; the factors A−k,A+k have the same parity and
product 4, which forces k=0, a contradiction.

For label8, ν=2. This last step is root's elementary elimination,
independently checked here. Set x'=2−√2. If A<0 then
R'=A−kx'<−(√2+x')k=−2k≤−2, impossible. If A>0 and k=1, A²=6 is
not a square integer. If k≥2, then

    A²−(2+kx')² = 4k(k−√2)/(1+√2)>0,

so R'>2, also impossible. A=0 cannot satisfy (10).

Therefore the following necessary inequality holds uniformly:

    N/gcd(νij,N) < j−i,
    ν=1 for labels5/12, ν=2 for label8.                (11)

In particular the originally requested N/gcd(ij,N)<j−i holds for
types5/12. Prime N allows at most one central cap of either of those
types. No stronger denominator cancellation is claimed for label8.

## 4. Remaining global obligation

These conclusions concern isolated finite CONTACT indices in precisely
(1). They do not assert that every central opposite wall is of this
kind, or count a shared nonisolated wall as two isolated supports.
The separate full-continuation envelope is needed to turn restrictions
on these indices into restrictions on short actual edge directions.

Moreover (6) means a positive-definite mixed replacement cannot contain
two of these caps: its two-root principal Gram is already indefinite.
The new mixed-replacement inequalities must therefore exploit long or
nonisolated actual opposite walls, or later actual vertices on the
faces between these disjoint caps. No whole-A_n farthest exclusion,
component energy floor or uniform dimension improvement follows here.
