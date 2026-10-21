# Direct rank-three attachments of two complete H17 sources

**Exact family obstruction.** The [checker](../scripts/fresh_h20_h17_rank_three_coupling_20261015.py) and [certificate](fresh_h20_h17_rank_three_coupling_20261015.json) classify direct H20 attachments of either of two complete H17 Coxeter sources to an arbitrary spherical rank-three support. Every such chamber has infinite volume. There is a nonordinary51-facet construction, but its exact escape and integral Cartan matrix prevent treating it as a new finite-volume candidate.

“Direct” means that the walls consist of the three support walls and walls inducing all the fixed source facets in their orthogonal H17 face. All finite Coxeter labels, reducible supports, and the unbounded `I2(m)+A1` family are allowed. Additional walls whose projection does not induce one of those source facets, changed source facets, or omitted source facets are outside this theorem. In particular, this does not exclude a later nonintegral completion using extra projected wall directions.

The19-facet source is the classical chamber for `U+E8+E8`; its graph is printed in [Roulleau, §15.4, p.90](https://epiga.episciences.org/10377/pdf), with attribution to Vinberg. The48-facet source is the reflective `(18,4,0)` diagram. Its24 subdivided `K4,4` roots and24 perfect-matching roots, including their parity and Klein-four cosets, are described in [Alexeev, Theorem1.4 and §§3–4](https://arxiv.org/pdf/2209.09110). The checker binds the repository’s exact source data and reconstructs the48 physical roots independently.

## A general propagation criterion

Normalize support roots and source roots to squared norm2. Let `J` be the inverse positive support Gram. It is entrywise nonnegative. An inducing ambient root, also normalized to norm2, has the form

\[
r_i=\lambda_i n_i-w_i,\qquad
w_i=\sum_j c_{ij}\omega_j,\quad c_{ij}\ge0,\quad
2\lambda_i^2+w_i^2=2,\quad\lambda_i>0.
\]

Call it ordinary when `w_i=0`. On a sourceπ3 edge, an ordinary root and a nonordinary root would have ambient product `−lambda`, strictly between−1 and0. That is not a Coxeter product for norm-two roots. Consequently ordinary/nonordinary status is constant on each sourceπ3 component.

If a component contains `k+1` pairwise orthogonal source roots, all of them cannot be nonordinary in a rank-`k` support. Their ambient products would be `w_i·w_j`; these are nonnegative by `J≥0`, while acuteness makes them nonpositive. This would give `k+1` nonzero mutually orthogonal vectors in a positive `k`-space. Thus the entire component must be ordinary. The argument also covers multiple inducing walls for a source facet: choosing any representatives gives the same constraint.

For the19-facet source theπ3 graph is connected and the literal source facets1,3,6,8 are mutually orthogonal. Therefore all19 inducing walls are ordinary for every rank-three spherical support. A strict timelike source point plus a sufficiently large vector in the strict negative support chamber is a strict spacelike point. This proves infinite volume without any arithmetic assumption.

## The48-facet source: the complete remaining classification

Itsπ3 graph has components of sizes24,12,12. The first component contains the four mutually orthogonal boundary roots `B1,…,B4`; hence all24 boundary and exceptional roots are ordinary. Each remaining `C_sigma` has a sourceπ4 contact with four ordinary exceptional roots. If it is nonordinary, that contact forces `lambda=1/sqrt2`: its scaled product lies strictly between−sqrt2 and0, whose only permitted nonzero Coxeter value is−1. Thus its transverse vector has squared norm1.

All spherical rank-three supports now reduce to a small exact weight calculation. A nonzero normalized contact magnitude is1, or at least`sqrt2`; parallel and dotted magnitudes are at least2. Positive inverse entries rule out multi-contact profiles except in `A1^3`. The inverse diagonals are:

| Support | Inverse diagonals | Norm-one profiles |
|---|---|---|
| A3 | 3/4,1,3/4 | middleπ3 |
| B3 | 1,2,3/2 | firstπ3 |
| H3 | all greater than1 | none |
| A2+A1 | 2/3,2/3,1/2 | A1π4 |
| I2(4)+A1 | 1,1,1/2 | either I2π3, or A1π4 |
| A1³ | 1/2,1/2,1/2 | oneπ4, or twoπ3 |

The missing singleton costs4/3 and3/2 fall strictly between the permitted squared magnitudes1 and2. For `I2(m)+A1`, `m≥5`, the I2 inverse diagonal is at least`1+1/sqrt5>1`, leaving only the A1π4 profile. This handles the full unbounded label tail.

Within a nonordinary12-component, every sourceπ3 edge requires transverse product1/2. Only `A1^3` has any such pair: its three distinct doubleπ3 weights. The source component is complete tripartite on three sets of four roots, with parallel pairs inside each set. The latter require transverse product0 or1, so each four-set has one common weight, and the three sets have different weights. Both parity components cannot be nonordinary: a cross-parity right pair would have positive transverse product. Therefore the only possibilities are all ordinary, or `A1^3` with exactly one nonordinary parity component. There are precisely twelve assignments, all related by a column transposition of the source and permutations of the three support coordinates.

## Literal51-wall representative and strict escape

Use source coordinates `(F1,F2,E11,…,E44)` with hyperbolic form `−2F1F2+sum Eij²`, and add three positive unit coordinates. Keep the24 boundary/exceptional roots and one parity class of twelve matching roots. For the other parity class replace

\[
C_\sigma\longmapsto C_\sigma-e_i-e_j,
\]

assigning the three pairs `(0,1),(0,2),(1,2)` to its three four-element cosets. Add the three support roots `e_i`. The resulting51 roots have **985 right,152π3,120π4 and18 parallel pairs**. A common strict center is

\[
H=(5,5,-1^{19}),\qquad H^2=-31.
\]

The exact rational vector

\[
X=(5,3,0,0,-3,0,0,0,-2,-1,-1,-2,0,0,0,-3,0,0,0,0,-3)
\]

has `X²=7` and every wall product nonpositive. Consequently

\[
E=X+H/1000,\qquad E^2=6949969/1000000>0
\]

has every wall product strictly negative and is future-oriented relative to`H`. This is a strict spacelike escape. Each facet is actual: the perpendicular projection of`H` to that wall is timelike and strict at every other wall. All twelve coordinate-isometric variants and their transported escapes are checked.

Scaling each boundary root and each changed matching root by1/2 gives51 rational generators spanning the full21-dimensional space. Every Cartan entry is integral. Their integer span is therefore a full rational lattice preserved by all51 reflections. This is a literal lattice statement, not an inference from rational similarity and not a reflectivity claim. The strict escape already settles the finite-volume verdict for these chambers.

Generation and saved replay pass with exact rational and symbolic radical arithmetic:

```sh
python scripts/fresh_h20_h17_rank_three_coupling_20261015.py --verify
```

The useful next construction step would have to leave this direct-inducing family, for example by adding noninducing walls or changing the retained H17 source. No such completion is constructed here, and no extra integral cap bank was started.
