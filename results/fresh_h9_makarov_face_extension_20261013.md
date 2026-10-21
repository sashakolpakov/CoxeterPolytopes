# A compact five-dimensional face as a new H9 construction seed

The new exact construction has **12 actual Coxeter facets in H9**, contains a compact Makarov five-dimensional prism as a face, and retains ten certified finite vertices. It still has infinite volume, certified by a strict spacelike point. No new compact polytope or new dimension is claimed.

The useful change of starting point is that the five-dimensional face is compact from the outset. This differs from the current geometric H9 branch, whose compact four-dimensional prism has been extended into a still noncompact five-dimensional face, and from lifting the eleven fixed wall directions of Bugaenko's compact H8 example. It also makes no restriction to an integral reflective lattice.

## Literature and dimension scope

The July 2026 preprint by Guo, Ma, Zang and Zheng still states the known example records as compact dimension 8 and noncompact finite-volume dimension 21. This is the standard frontier in the current primary sources checked, rather than an inference from an old survey. [Introduction, *Several families of incommensurable noncompact hyperbolic Coxeter polytopes*](https://arxiv.org/pdf/2607.14715).

Dimension 20 is itself another meaningful target. Kellerhals explicitly records the absence of a known cofinite example in dimension 20 in her 2017 introduction; Felikson's author-maintained overview lists the known dimensions as at most 19 and 21. The 2026 record statement does not separately settle the missing dimension, and no later H20 example was found in this check. Thus the careful status is “H20 remains unknown in the sources checked,” not a nonexistence assertion. [Kellerhals, p.1](https://doi.org/10.3390/math5030043), [author overview](https://www.maths.dur.ac.uk/users/anna.felikson/Polytopes/polytopes.html).

The exact source face is the compact prism \(P_0=[5,3,3,3,3]\) in Emery–Kellerhals, §1 and §4. They attribute \(P_0\) and \(P_2=[5,3,3,3,4]\) to Makarov and state that both are arithmetic. These individual blocks must not be confused with nonarithmetic examples obtained by particular gluings. [Emery–Kellerhals, *The three smallest compact arithmetic hyperbolic 5-orbifolds*](https://homeweb.unifr.ch/kellerha/pub/IML-1112s-21.pdf).

Borcherds' face criterion, presented with explicit angle rules by Allcock, makes a spherical H4 face a natural anchor: H4 has no connected spherical extension, so ambient facets meeting that entire face are orthogonal to H4. This criterion runs from an ambient polytope to its face; it does **not** guarantee a reverse extension. Here the reverse extension is an explicit construction problem. Doubling or gluing a known chamber alone preserves its dimension. The dimension increase below instead comes from adjoining a four-dimensional positive normal space and then seeking its bounding walls. [Allcock, Theorems 2.1–2.2](https://web.ma.utexas.edu/users/allcock/research/doubling.pdf).

## Exact seed

Set \(\phi=(1+\sqrt5)/2\). Let \(G\) be the six-by-six Gram matrix with diagonal 2 and consecutive off-diagonal entries

\[
-\phi,-1,-1,-1,-1.
\]

All other off-diagonal entries vanish. Write \(I=G^{-1}\). The first six face normals are the coordinate vectors; the seventh is

\[
w_6=-I e_5,\qquad w_6^2=(7-\sqrt5)/11.
\]

It is orthogonal to \(w_0,\ldots,w_4\), and \((w_6,w_5)=-1\). The normalized dotted weight squared is \((7+\sqrt5)/8>1\). The exact checker enumerates all five-wall supports, obtains ten feasible future timelike rays, and verifies their positive-definite incident Grams. A strict timelike center \(h\), expressed as a strictly positive combination of all seven normals, bounds a polyhedral affine slice. Thus the vertex census proves compactness directly; it does not merely trust the literature identification.

Let \(H\) be the norm-two H4 matrix with chain labels \([5,3,3]\), and \(B=H^{-1}\). The new ambient form is

\[
Q=G\oplus H,\qquad \operatorname{signature}(Q)=(9,1).
\]

The seven old normals and four H4 simple roots give eleven acute Coxeter facets. The ten old face vertices lift to actual finite H9 vertices, with a positive-definite five-wall old block and the orthogonal H4 block. A small explicit perturbation of \(h\) supplies a common strict timelike point.

There is also an elementary all-label version of the H4 face gate. Every entry of \(B\) is positive and

\[
\min_i B_{ii}=3+\sqrt5>2.
\]

For a norm-two Coxeter normal nonorthogonal to H4, at least one negative pairing has magnitude at least one. Its H4 projection therefore has square at least \(3+\sqrt5\), and its old component is timelike. It cannot meet the entire H4 face. This includes arbitrary finite labels and parallel/dotted relations. See the independently implemented [H4 Schur gate](fresh_h9_h4_transverse_gate_20261013.md).

## A whole first family is closed, then escaped by a changed center

For each of the ten compact-face vertex vectors \(v\) and each H4 fundamental weight \(\omega_j=B e_j\), the first family is

\[
r_{v,j}(t)=\left(\sqrt{\frac{t^2 B_{jj}-2}{-v^2}}\,v,-t\omega_j\right),
\qquad t=2\cos(\pi/m),\quad m\ge3.
\]

These are norm-two walls with angle \(\pi/m\) to H4 wall \(j\) and right angles to the other three. The bounded screen tests the five labels \(3,4,5,6,10\), giving 200 proposals; 152 pass every old-wall angle and 48 do not. This count is not a mirror-distinctness or exhaustive-cap claim. For those checks, squared cosines lie in \(\mathbb Q(\sqrt5)\). The cyclotomic degree bound leaves labels \(2,3,4,5,6,8,10,12\); labels 8 and 12 generate different quadratic fields, leaving exactly the six values tested by the gate.

More strongly, one archived integer spacelike point is strict for **all forty families at every real \(t\ge1\)**. In particular, no finite selection from those families can yield finite volume. The proof checks all forty \(t=1\) pairings exactly by positive squared inequalities. After dividing a cap by \(t\), the coefficient of its negative old-point pairing is \(\sqrt{B_{jj}-2/t^2}\), which increases with \(t\). Thus no angle cutoff is involved. See also the [monotonicity proof](fresh_h9_h4_transverse_monotonicity_20261013.md).

Putting all four cap centers at one old point is another poor template: for each adjacent pair of H4 weight rays, the checker gives the exact positive number

\[
B_{ij}^2-(B_{ii}-1/2)(B_{jj}-1/2)>0.
\]

Consequently such two norm-two caps have positive pairing for every pair of finite angle labels, and cannot both be retained acute facets. Distinct centers are necessary for this four-ray template.

The subsequent small trial releases the vertex-center condition. Its accepted exact cap has old component

\[
u=-I e_0+q I e_4,
\quad
q=\frac{-b+\sqrt D}{2a},
\]

where

\[
a=(6-4\sqrt5)/11,\quad b=(34+14\sqrt5)/11,
\quad D=(336+128\sqrt5)/11.
\]

Its full normal is \(r=(u,-B e_0)\). It has norm two. Its eleven old relations are seven right angles, two \(\pi/3\) angles (old wall 0 and H4 wall 0), and two strictly dotted pairs (old walls 4 and 6). Its old center lies on the compact prism's edge with active walls \(1,2,3,5\), rather than at one of its vertices. The cap cuts the common escape of every vertex-centered family.

All of this is checked over \(\mathbb Q(\sqrt5,\sqrt{11})\) using rational coefficient arithmetic and exact signs: the checker verifies \(\sqrt D=(16\sqrt{11}+4\sqrt{55})/11\). Its internal encoding is \(A+B\sqrt D\), with \(A,B\in\mathbb Q(\sqrt5)\); the form and original seed use only the base field. Projecting the strict center onto the new wall produces a timelike point strict on every other wall, proving the new facet is actual. The compact H5 face and all ten previous finite vertices survive strictly. However, a second explicit integer spacelike point is strict on all twelve walls. The changed center solves a real obstruction, but the resulting chamber still has infinite volume.

## Concrete next experiment and comparison with existing work

The next bounded experiment is to solve several such old edge/face-center systems **jointly**, with at least one cap on each required transverse region, and impose all mutual angles before any volume test. For a fixed old codimension-three face and two finite-angle anchors, the old normal is determined by a quadratic norm equation. Every proposed coherent set then has an explicit Gram matrix; compactness requires all extreme rays of its pointed affine cone slice to be future timelike. Covering the ten old vertices or cutting a few chosen escapes is insufficient.

This is a reverse-face construction using a completely specified compact H5 source. It is separate from the [fixed H8 direction lift](fresh_h8_lift_h9_20261013.md), the [single-top H4-to-H5 prism obstruction](fresh_h9_prism_lift_20261013.md), and the existing arithmetic lattice searches. If a compact completion with fourteen facets is sought, its full nerve must additionally survive the eight open shadows in the [global H9/14 report](h9_14_facet_completion.md); ordinary combinatorial incidence is not a metric realization. No claim of worldwide novelty is made for this infinite-volume seed or cap.

## Reproduction

- [Exact standalone checker](../scripts/fresh_h9_makarov_face_extension_20261013.py)
- [Frozen source, normals, vertices, family obstruction and changed-center cap](fresh_h9_makarov_face_extension_20261013.json)
- [Independent H4 matrix gate](fresh_h9_h4_transverse_gate_20261013.md)
- [Source/proof review and fresh replays](fresh_h9_makarov_face_extension_20261013_review.md)
- [Separate sparse-radical arithmetic replay](fresh_h9_makarov_face_extension_20261013_independent.md)

```
python3 scripts/fresh_h9_makarov_face_extension_20261013.py --verify
python3 scripts/fresh_h9_h4_transverse_gate_20261013.py --verify
```

Both fresh replays pass. Numerical linear programs and the tiny anchor scout supplied proposals only; neither a tolerance nor a numerical eigenvalue is used by the certificate. The saved cap is one exactly reconstructed trial, not an exhaustive census of the face-center systems.

The separate sparse-radical replay also passes: all twelve norms, all 66 pairs, the strict center, ten finite vertices, the actual new facet point, both escape comparisons and all forty family inequalities. It imports no producer arithmetic. A second reviewer checked the compactness and all-real-parameter arguments and freshly replayed the primary certificate.
