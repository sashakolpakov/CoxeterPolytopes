# Exact small-face completions obstruct an angle-by-angle dimension-29 strategy

Every marked pair in E7 and E8 admits an exact compact triangular or quadrilateral face of a local Coxeter chamber, with finite corners and no parabolic subdiagram. The same holds for all marked pairs in A4, B4, D4, F4 and H4. These are **separate local realizations**, not a simultaneous completion or a compact chamber in dimension29. The projected face angles need not themselves be Coxeter angles.

This gives a specific barrier to the preceding [dimension-29 target](fresh_compact_dimension_bound_20261015.md). The extremal forest $E_7+E_7+E_7+E_8$ has rank29 and Vinberg weight87. Each of its positive-weight marked pairs can occur on a compact small face in a realizable local star. Therefore one cannot eliminate this forest by proving that three units of its marked-angle weight have no triangle or quadrilateral completion using only the marked vertex, its finite neighboring corners, Lorentz inertia, and absence of parabolic subdiagrams. Simultaneous compatibility of several stars, further mandatory facets, or a global incidence/weight argument remains necessary.

The [exact source](../scripts/fresh_compact29_small_faces_20261015.py), [literal recipes](fresh_compact29_small_faces_20261015_recipes.json), and [certificate](fresh_compact29_small_faces_20261015.json) are separate from the frozen dimension-bound package.

| Vertex component | Marked pairs | Triangle witnesses | Quadrilateral witnesses |
|---|---:|---:|---:|
|A4|6|6|0|
|B4|6|6|0|
|D4|6|6|0|
|F4|6|6|0|
|H4|6|6|0|
|E7|21|13|8|
|E8|28|13|15|
|Total|79|56|23|

The numerical scouts suggested the recipes. Their unsuccessful trials are not used as exclusions. In particular, the table does not claim that the quadrilateral cases cannot also have triangle completions. Each displayed success is rebuilt exactly without numeric eigenvalues.

## Local constraints and exact realization

For an irreducible spherical vertex diagram $T$ of rank $r$, let $a,b$ be its marked nodes and let $W=T\setminus\{a,b\}$. All root squares are2. A triangle recipe adds one root $c$; a quadrilateral adds $c,d$, with cyclic black order $a,b,c,d$. Every ambient pair is a finite Coxeter pair or a strictly dotted pair.

The checker verifies positive definiteness of every $W$ plus two consecutive black roots. Thus all three or four corners are finite. Let $R$ be the Schur complement of the positive-definite white Gram matrix.

- For a triangle, $R$ has size3 and negative determinant, while each adjacent two-node block is positive definite. Hence its signature is $(2,1)$.
- For a quadrilateral, $R$ has size4, determinant zero, a nonsingular Lorentz three-node principal block, and positive-definite adjacent pairs. Its signature is $(2,1,0)$. Both opposite pairs are indefinite.

All triangle entries lie in $K=\mathbb Q(\sqrt2,\sqrt5)$. Each quadrilateral uses one positive quadratic radical over $K$. Its sole dotted entry has magnitude $y=(-B+\sqrt D)/(2A)$, with $D=B^2-4AC$, where the exact determinant equation is $Ay^2+By+C=0$. The checks $A>0,C<0,D>0$ and $4A+2B+C<0$ prove that this is the positive root and that $y>2$. No rounded scout value is treated as an exact distance.

The certificate also gives the three or four vertex rays in a Lorentz three-dimensional residual basis. Each ray is timelike, annihilates exactly its two incident black roots, and pairs strictly negatively with the other black roots. All pairwise scalar products of vertex rays are negative, so they lie in the same time component. Their geodesic convex hull is the required compact triangle or quadrilateral, with precisely the indicated sides.

Their sum is timelike and strictly behind all black walls, while orthogonal to $W$. A sufficiently small perturbation in the negative white-dual directions gives a strict ambient chamber point. Keeping one white equation zero gives a relative interior point of that white facet. Likewise, starting inside a black side and perturbing in the orthogonal complement of its own normal preserves that side while making the white inequalities strict. The positive-definite support matrices guarantee these perturbations exist. Thus these are actual local face stars, not merely formal Gram signatures.

Every principal subdiagram is checked for a parabolic obstruction. The implementation uses exact Schur pivots: a negative pivot, or a zero diagonal with a nonzero off-diagonal entry, certifies failure of positive semidefiniteness. A singular positive-semidefinite principal block would fail the certificate. There are **25,457 exact principal tests**; subsets containing the quadrilateral's strictly dotted pair are rejected as positive-semidefinite candidates by their negative two-node determinant. The certificate records those additional subset counts separately.

## Why the obstruction persists in the rank29 target

Choose one of the $E_7$ or $E_8$ components in $E_7^3+E_8$, and use its appropriate local witness. Add the remaining three spherical components orthogonally to the complete star. The full white support now has rank27. The ambient signature becomes $(29,1)$, with one additional Gram zero in the quadrilateral presentation. Every neighboring corner has the required positive-definite rank29 support. Orthogonal spherical summands cannot introduce a parabolic principal subdiagram.

Consequently all positive-weight marked pairs of this rank29 forest are individually locally feasible. We verify all28 E8 pairs, so no identification of the four zero-weight E8 positions is needed. Vinberg assigns weight1 to all21 E7 pairs and24 of the28 E8 pairs, giving $3\cdot21+24=87$. [Vinberg, §6.1–2](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt).

This result is consistent with Vinberg's separate E7 warning about uniform nonnegative local weights that vanish between components. The new certificate instead gives explicit geometric witnesses for the marked-angle exclusion proposal, including E8 and all rank-four types. It does not prove that arbitrary rank29-specific reweightings are impossible, or that the local stars coexist in a single compact polytope.

## Remaining constructive and global questions

No new finite-volume example or improved unconditional dimension bound follows. The present angle-by-angle target is blocked at the explicit $E_7^3+E_8$ forest. A next global test would couple the stars of several marked angles at the same finite vertex: their new roots must share one Lorentz ambient space, satisfy mutual Coxeter pair conditions, and preserve all claimed faces. The independent one-star constructions do not establish that joint feasibility.

Generation and saved replay use:

```sh
python scripts/fresh_compact29_small_faces_20261015.py
python scripts/fresh_compact29_small_faces_20261015.py --verify
```

The [independent review](fresh_compact29_small_faces_20261015_review.md) reconstructs all79 stars with separate SymPy field arithmetic, checks all306 candidate pair intersections (260 feasible timelike rays and46 rejected opposite pairs), and repeats the complete no-affine subset checks. Generation and saved replay pass. One E8 quadrilateral remains in degree4; the other22 quadratic extensions have degree8.

The compact H9 cursor and every earlier frozen source remain unchanged by this package.
