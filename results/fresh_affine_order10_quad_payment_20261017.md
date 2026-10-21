# Ideal quadrilateral payment for every affine factor of order at least10

**Final internally reviewed proof.** Suppose a finite-sided,
finite-volume hyperbolic Coxeter polytope has no connected affine
subdiagram of order at most9. Then every quadrilateral having an ideal
vertex receives at least one under the finite Vinberg rule and the affine
rule paying distance at most3 or a fork-invariant end segment of order
at most6. Arbitrarily many other cusp factors and supporting spherical
spectators are allowed. This note alone makes no ambient dimension claim.

The frozen [large-factor proof](fresh_large_affine_factor_quad_payment_20261017.md)
established the assertion when the varying affine factor has order
m>=18. Here the complete classification and scalar arguments extend
to m>=10. There is exactly one new unpaid pair of exceptional profiles;
it is computed explicitly below. All predecessor files remain unchanged.

## Geometry and the generic arguments

Write the corners as W+a+b, W+b+c, W+a+d, W+c+d. At the chosen ideal
corner its varying affine factor F contains a,b; the other roots form
a positive definite spectator forest Z right to F. Every actual edge
diagram W+one black is positive definite. Set h=−c.nu and k=−d.nu for
the positive null combination of F. Since c,d are the opposite sides,
the actual ideal vertex belongs to neither, and h,k>0.

The order-nine prohibition makes F classical A/B/C/D of order m>=10.
If b.c, a.d or c.d is nonzero, the corresponding adjacent marked pair
is already paid. Assume all three vanish. Eliminating Z changes the
normalized four-root determinant by the nonnegative quantity

` alpha*(p/h²−2q/(hk)+r/k²), `

where [[p,q],[q,r]] is positive semidefinite and alpha>0. Thus strict
positive determinant contradictions before eliminating spectators
remain contradictions afterward. Finite-component payments retain the
separate role checks described below.

For an affine A cycle, its deleted A_(m−1) path has order at least9,
so the complete finite/affine one-root extension list in the
[cycle proof](fresh_simple_cusp_affine_a_quadrilateral_20261017.md) is
unchanged. The simple-leaf E8 witness uses two cycle arcs of lengths2,5;
m>=10 leaves them nonadjacent at their ends. All other low-contact
exclusions are smaller. At an unpaid corner both old arc lengths are
at least4; because their sum m>=10, at least one exceeds4. The original
strict reciprocal-arc monotonicity argument therefore applies. The
PSD spectator correction preserves its strict sign.

For trees, use the same cusp coefficients and resistance Green identity
as the [complete tree proof](fresh_simple_cusp_affine_tree_completion_20261017.md).
At main-path marks with separation P>=4, write the outer types B_r or
D_r, with B1 allowed and D2 denoting the two fork leaves. The order is
m=r+s+P+1. Single simple leaves are now excluded using the
[order10 E7/E8 witness lemma](fresh_affine_carriers_order10_20261017.md):
available opposite simple lengths total m−3>=7; use arms(1,3,3) or
(1,2,5), with the same explicit near-decoration exceptions. Doubled
leaves still have the small star/path witnesses. For an outer fork-leaf
attachment the long arm has at least m−4>=6 simple steps, so the
previous E8 exclusion remains available as well.

Consequently all ordinary height bounds and arbitrary outer-contact
Cauchy–Schwarz bounds remain unchanged. The largest ordinary main-path
contribution is

` g_* = 4(sqrt2−1), `

and the separated-support determinant is at least
`1−2g_*/P >= 3−2sqrt2 >0`.
At P=4 the two second-node contacts may coincide; the same Green entry
is valid, as in the inherited proof.

## Complete small-inner accounting

The following rank statements are enough; no new label grid is used.

* An extension joining both old end components has full order m>=10.
  Its finite/affine connected type is therefore classical. The complete
  B/D port list is unchanged. If it joins only one of the two D2 leaves,
  its component has order m−1>=9: a finite exceptional type is still
  impossible, and an exceptional affine type of order9 is forbidden.
* The inner B block has rank at least P+1>=5 and has no exceptional
  finite extension. The inner D block has rank at least P+2>=6.
  Its only extra finite extensions are D6→E7 and D7→E8.
  Exceptional affine extensions have at most9 nodes and are forbidden.
* D8 and larger D blocks admit only the ordinary classical profiles.
  Thus an ordinary opposite block need not have rank at least9:
  rank8 suffices for D, and rank5 for B.

Every D6→E7 case is paid at that adjacent corner. The new root is at
the end of E7's length-two arm, with inverse diagonal2. A spectator
can contact that component only through the new root; a nonzero acute
Coxeter contact would give the forbidden affine E7. Thus no spectator
can turn this paid E7 into a larger finite component.

For D7→E8 there are two main-path positions: D2/P5 and D3/P4. The
first has a marked branching vertex, paid by Vinberg's published
special-vertex lemma. In either position the new root has E8 inverse
diagonal4, so a spectator contact costs more than its square2 and is
impossible. The old finite role therefore persists.

The D3/P4 scalar proof already covers every ordinary opposite profile,
including all outer-only contacts and the added classical affine
profiles. If the opposite inner block is D6, its additional E7 profile
is paid as above. If both sides have D7 exceptions, the rank equations
give just two possibilities:

| Outer ends | P | m | Disposition |
|---|---:|---:|---|
| D2,D2 |5|10|Each E8 extension has the paid branching role |
| D3,D3 |4|11|New dual-E8 calculation below |

Indeed both exceptional inner ranks equal7, so r=s=7−P; with D ends
r,s>=2 and P>=4, only P=4 or5 is possible. There are no additional
small-small cases hidden at orders12–17. The checker records872
boundary-order main-path contexts as corroboration of this elementary
rank accounting; it is not a substitute for the finite classification.

## The affine-D11 dual-E8 case

Number the old main path0,...,6, with fork leaves7,8 at0 and9,10 at6.
The marked roots are a=1,b=5. The new root c is simple to10 and has
opposite product a.c=−x; d is simple to8 with b.d=−y. All other old
products of c,d vanish, as does c.d. Eliminating the nine old white
roots gives the exact matrix, in order a,b,c,d,

```
 1/4   -1/4   -x    -1/2
-1/4    1/4   -1/2  -y
 -x     -1/2   5/4   0
-1/2    -y      0    5/4
```

The heights are h=x+1/2 and k=y+1/2. If x=0, the new simple leaf and
the old fork give an induced affine E8 of order9. If x=1, the path
a=1,2,3,4,5,6,10,c is an induced all-simple eight-cycle, affine A7.
Both are forbidden. The same holds for y by symmetry. Since the first
two nonzero Coxeter magnitudes are1 and sqrt2, this proves
x,y>=sqrt2, including all parallel and dotted possibilities.

Put u=1/h,v=1/k, and z=1/(sqrt2+1/2). The exact normalized determinant is

` F(u,v)=1−u−v−(u²+v²)/16+uv/2. `

For 0<=u,v<=z both partial derivatives are at most−1+z/2<0. Therefore

` F(u,v)>=F(z,z)=181/98−62sqrt2/49>0. `

The checker derives this matrix and polynomial directly from the literal
13-root Gram, verifies all four affine kernels used to exclude x,y=0,1,
and confirms the two inverse-diagonal4 E8 roles. The positive determinant
contradicts the actual face's corank one, and remains positive after any
PSD spectator correction. In fact the role calculation already forces
those spectators right to c,d.

## Terminal and fork-leaf marks

Deleting a terminal doubled node or a fork leaf leaves a connected
finite B/D diagram of order m−1>=9. The complete four-profile lists are
unchanged: isolated, long-end3, long-end4 and next-node3, with the
appropriate fork-mate interpretation. All low-height exclusions use
the order10 leaf lemma or the old small stars.

A doubled terminal with a main mark has P>=6, contribution at most5/2,
and an opposite B_P inner block, which has no exceptional extension.
A fork leaf with a main mark has resistance R>=6, contribution at most4,
and opposite D_(P+1). Its E7 extension is paid; its P6/E8 exception is
exactly the inherited fork-leaf exception. The ordinary opposite bound
g_* keeps both determinant margins positive. The finite D8 and E8
payments in the exceptional branch are preserved by their original
marked-role arguments, which forbid the only possible affine spectator
extensions of order at most9.

Two unpaid terminal marks belong to opposite ends and have R=m−1>=9.
Their total contribution is at most8; their local contact supports stay
separated. Hence the normalized determinant is at least1−8/9=1/9>0.
This replaces the unnecessarily large old resistance17 bound.

The three mark-position classes exhaust all affine classical trees.
Together with the cycle case and the PSD spectator argument, every
ideal-corner quadrilateral is paid under the order-nine prohibition.

## Verification scope

Generation passes the new rank accounting, physical dual-E8 elimination,
four positive affine kernels, exact positive margin and terminal9 bound.
The proof explicitly rechecks every use of the former m>=18 assumption;
the inherited all-rank Green identity, finite/affine extension lists and
other scalar cases remain dependencies. No old large census was rerun,
and no finite label or rank cutoff is used to infer the theorem. Generation
and fresh saved replay pass. A separate root mathematical/source read
passes the grounded rank thresholds, complete exceptional accounting,
literal dual-E8 kernels and determinant, and terminal resistance9
argument. That read does not duplicate the872 contextual checks or the
inherited classifications' arithmetic. Syntax, direct source hashes and
local links pass. Source, data and report are final and frozen. This is
internal review, not external refereeing or formalization. The finite-only
quadrilateral transfer, triangle payment and global curvature composition
are separate inputs before an ambient bound can be stated.
