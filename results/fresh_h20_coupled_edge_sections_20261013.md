# Coupled two- and three-contact sections through source edges

This search constructs new local Coxeter corners using **two or three active wall pairings**, then tests the full H20 sections. All **264 bounded proposals** fail at an actual neighboring finite vertex, with exact certificates against all 210 source walls. A separate finite reduction also closes **all real transverse positions in three specified coupled F4 families**: every one of their 84 possible closure candidates has an actual obtuse ridge. No new Coxeter polytope is established.

The files are the [bounded scout](../scripts/fresh_h20_coupled_edge_sections_20261013_scout.py), its [recipes](fresh_h20_coupled_edge_sections_20261013_scout.json), the [exact checker](../scripts/fresh_h20_coupled_edge_sections_20261013.py), and the [264-case exact certificate](fresh_h20_coupled_edge_sections_20261013.json). The continuous F4 reduction has its own [candidate census](fresh_h20_coupled_edge_sections_20261013_f4_scout.json) and [84-case exact certificate](fresh_h20_coupled_edge_sections_20261013_f4.json).

## Genuinely coupled local templates

Let E be the twenty incident walls of an actual source edge, with positive-definite Gram G. Normalize a section normal e to squared norm one, and let a be its twenty pairings. Put

\[
 v=G^{-1}a,\qquad w=\sum_{i\in E}v_i r_i,\qquad
 \sigma=a^TG^{-1}a,\qquad \gamma=1-\sigma.
\]

If \(\gamma>0\), choose an interior timelike edge point P and a spacelike vector z in the Lorentzian edge plane, with \(z\perp P\), \(z^2=\gamma\). Then \(e=w+z\) is a unit spacelike normal transverse to the edge at P. The initial section vertex has projected Gram \(G-aa^T\), which is positive definite because \(\sigma<1\). Every incident wall is an actual facet there.

The following templates pass the exact local norm and angle gates. Unlisted pairings are zero. Chain entries are in chain order; the A1 factor in the two joins has old squared norm d.

| Original component(s) | Active pairings | New finite component(s) | \(\sigma\) |
|---|---|---|---|
| A_n, n=2,…,20 | first two entries \((1,-1)\) | A1 + B_(n−1) | \(n/(n+1)\) |
| A1 + A2 | \((\sqrt{d/5},1)\), second at an A2 endpoint | B3 | \(13/15\) |
| A1 + A3 | \((\sqrt{d/5},1)\), second at an A3 endpoint | F4 | \(19/20\) |
| A3 | \(\sqrt6(1/3,-1/2,1/3)\) | A2 + A1 | \(5/6\) |

The opposite signs in the first template cancel the first old edge; the remaining chain acquires an end label 4. In the last template both old edges cancel, while the two old leaves acquire a label-3 edge. Thus the three-contact search is not a restatement of the earlier one-contact A4-to-H4 construction.

After examining the saved 6,492 literal edge supports, the bounded scout selects one source edge for each template and available unordered endpoint-type pair (finite/finite, finite/ideal, ideal/ideal). This gives 44 selections. Each uses the three projective positions 1/3, 1/2, 2/3 between the two normalized source endpoints, and both signs of z: 264 recipes in total. These counts do not assert that the literal choices are inequivalent under the full source symmetry group.

## Exact global failures of the bounded proposals

Every sampled section has a valid initial finite Coxeter corner. Each nevertheless has a forbidden actual ridge at a neighboring finite vertex: **261 are obtuse**, and **three have squared cosine strictly between 1/4 and 1/2**, the consecutive Coxeter values for labels 3 and 4.

These are all exact failures, not just numerical proposals. The exact checker imports no numerical scout arithmetic and uses no saved floating-point coordinates. It reconstructs the roots and points using rational matrices and the separate sparse square-root arithmetic implementation. Only combinatorial witness indices and rational recipes are taken from the scout.

For an old edge wall j, let \(w_j\) be its old dual vector. The section-edge dual is

\[
 T_j=w_j-\frac{(G^{-1}a)_j}{\gamma}z.
\]

The checker proves \((T_j,e)=(T_j,P)=0\) and \((T_j,r_i)=\delta_{ij}\) for all twenty old incident walls. For a proposed next wall l, put \(d=(T_j,r_l)<0\) and \(p=(P,r_l)<0\). It avoids algebraic division by using the positive scaled contact

\[
 X=(-d)P+pT_j.
\]

All **210 inequalities** are checked at X. It is future timelike, lies in the section, and has exactly the nineteen retained edge walls and l as its zero set; every other inequality is strict. Convexity along the line then proves that this is the first existing-wall contact in that direction. The retained nineteen projected normals are positive definite, and \((T_j,r_l)\ne0\) proves independence of the new one. Since all twenty are orthogonal to timelike X, they form a positive-definite incident system. The bad pair therefore gives a genuine local ridge of a simple finite vertex. No conclusion relies on a Gram pair belonging to redundant or disjoint inequalities.

## Complete continuous test of three coupled F4 families

The sampled F4 families have source endpoint-orbit pairs **(0,28), (0,9), (16,173)**. Their old components, in the order A1 followed by A3, are

\[
 (9;4,48,14),\quad(4;9,29,26),\quad(14;9,29,26).
\]

For these three families both active source roots have norm two, and the active pairings are \(\sqrt{2/5}\) and 1. The residual norm is \(\gamma=1/20\). The following argument reduces **every real transverse position**, with both choices of z, to finitely many rational centers.

Suppose such a section were Coxeter and finite-volume. Its actual F4 face has dimension 16. This face has finite volume inherited from the section; **it need not be compact**, since F4 does occur in affine diagrams. At a finite point of any boundary facet of this face, the extra incident wall must be orthogonal to all four F4 normals: F4 has no proper connected finite spherical extension. A boundary facet has finite points, even when some of its vertices are ideal. Consequently all the face's bounding normals arise in this orthogonal way, and its angles are inherited Coxeter angles.

A finite-volume H16 polytope needs at least seventeen bounding hyperplanes. Only sixteen come from the complementary walls of E, so some original wall l outside E must bound the F4 face. Orthogonality to the two active projected normals says

\[
 (r_l,r_{A1})=\sqrt{2/5}\,(e,r_l),\qquad
 (r_l,r_{A3\text{ end}})=(e,r_l).
\]

The left sides are integers, while their required ratio is irrational. Both vanish, and \((e,r_l)=0\). The other two F4 orthogonalities likewise force their original pairings to vanish. Thus l is orthogonal to the entire original A1+A3 component. Since w lies in its span, \((w,r_l)=0\), and therefore \((z,r_l)=0\).

Let \(p_l\) be the rational orthogonal projection of \(r_l\) to the two-dimensional Lorentzian edge plane. It is nonzero: otherwise the outside wall would contain the whole actual source edge. The relation \(p_l\perp z\) makes \(p_l\) timelike and parallel to the crossing point P. Hence P must be the future normalization of one of the 210 rational wall projections. This is an exhaustive necessary condition, not a sampled parameter set.

The exact census checks every eligible old wall, rejecting nontimelike projections and points outside the open edge with rational signs. The three families leave **17, 11, 14 distinct centers**, respectively. Each has exactly two possible normalized z vectors, giving **84 section normals**. Their separate exact replay checks the complete projection census, both branches at every center, and an actual obtuse ridge against all 210 walls in every case. Thus these three continuous F4 families are closed. This is not a classification of every coupled F4 template in the source census, nor of other two- or three-contact local diagrams.

## Replay

```sh
python scripts/fresh_h20_coupled_edge_sections_20261013.py --verify
python scripts/fresh_h20_coupled_edge_sections_20261013.py --f4 --verify
```

The producer and saved-result replays pass. The initial local templates remain useful construction formulas, but none of the exact sections certified here is a Coxeter polytope. The other continuous coupled families remain open.

The [root mathematical/source review](fresh_h20_coupled_edge_sections_20261013_review.md)
includes fresh saved replays of both the 264-case and 84-case certificates.
