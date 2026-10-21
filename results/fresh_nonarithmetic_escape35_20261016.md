# A compact Lambert quadrilateral and an exact 36-wall chamber

An exact new wall cuts the strict escape of the retained 35-wall H9 chamber. It closes two old finite-edge ends into a compact Lambert quadrilateral, with three right angles and one angle π/6. The resulting chamber has **36 actual facets, 100 finite vertices and 266 finite edges**.

It still has **infinite volume**. A different strict spacelike escape is certified, and the complete cone calculation finds **1,068 hyperideal rays and no ideal rays**. The finite-corner germ count increases from 356 to 368. The compact quadrilateral is a local geometric improvement; neither the wall count nor that face closure is a finite-volume claim.

The [producer](../scripts/fresh_nonarithmetic_escape35_20261016.py) and [exact certificate](fresh_nonarithmetic_escape35_20261016.json) preserve the frozen [35-wall predecessor](fresh_nonarithmetic_escape34_20261016.md) and its [independent review](fresh_nonarithmetic_escape34_20261016_review.md). The separate [complete-cone calculation](fresh_nonarithmetic_escape35_20261016_cone.md) starts from the independently checked complete 32-wall cone.

## Exact recipe and Coxeter contacts

Use the original physical form Q, retained normals n_i and common strict future timelike center H. The new wall has index 35. Its right support is

\[
S=[0,1,2,4,5,6,7,9].
\]

The nine covectors indexed by S together with wall 8 have rank nine. Their kernel is literally the last coordinate vector k. Define p by

\[
(p,n_i)=0\quad(i\in S),\qquad (p,n_8)=-1,\qquad (p,k)=0.
\]

The exact identities are

\[
p^2=-5-\frac{5\sqrt5}{2},\qquad
k^2=\frac{\sqrt5}{4},\qquad (p,k)=0.
\]

Set

\[
D=10+\frac{68\sqrt5}{15},\qquad r=p-\sqrt D\,k.
\]

Then r²=2/3. Both signs of the square root are checked; the other branch is obtuse to a retained wall. The chosen r retains H and has positive pairing with the preceding strict escape.

Its 35 new wall relations are exactly eight right contacts, one π/6 contact with wall 8, and 26 strictly dotted contacts. In particular,

\[
(r,n_8)=-1,\quad n_8^2=2,\quad
\frac{4(r,n_8)^2}{r^2n_8^2}=3=4\cos^2(\pi/6).
\]

Every unspecified contact is proved strictly dotted; no assumed finite-angle alphabet or numerical fitting accepts a wall. The new wall is dotted to old walls 3 and 12, so it is outside the preceding obstruction for walls perpendicular to the compact tetrahedral face's six-root support. Direct projection of H to every wall supplies 36 future timelike relative-interior facet witnesses.

A bounded numerical screen selected this recipe. That discovery screen makes no completeness or exclusion claim. Acceptance uses only the exact linear recipe, both quadratic branches and all old-wall signs and Coxeter relations.

The arithmetic adds a formal positive square root of D to the preceding formal tower. Each scalar is encoded as [a,b], meaning a+b√D; a and b retain their predecessor encoding. Opposing-sign terms are compared by a²−Db² recursively. Division is only by nonzero scalars of the frozen base field. No irreducibility or minimal field-degree assertion is needed.

## Complete finite-corner update

The new facet's complete finite-neighbor list is S together with wall 8. All nine eight-element subsets are checked. Exactly two give positive-definite nine-root vertex supports:

| Vertex | Incident walls |
|---:|---|
| 98 | 0,1,2,4,5,6,7,9,35 |
| 99 | 0,1,2,4,5,6,8,9,35 |

Their actual future timelike points are strict at every other wall. All 98 old finite vertices are strictly retained. This proves the complete 100-vertex census, also recovered independently of this local enumeration by the separate full-cone clipping calculation.

The three new finite edges are 6–98, 8–99 and 98–99. Together with the old edge 6–8 they form a closed four-cycle. Thus two old finite-edge ends close, while

\[
9\cdot100-2\cdot266=368
\]

finite-corner germs remain open. These germs are distinct from the full cone's hyperideal rays.

## Entire compact two-face

The four vertices have common positive-definite support

\[
F=[0,1,2,4,5,6,9].
\]

The boundary walls in F-perpendicular Lorentz three-space are 7, 8, 11 and 35. Their projected norms are respectively 5/4, 2, √5/4 and 2/3. Adjacent pairs are 7–11, 11–8, 8–35 and 35–7. Their angles are right, right, π/6 and right. The two remaining pairs, 7–8 and 11–35, are strictly dotted.

The checker constructs a nondegenerate Lorentz three-frame and enumerates all six boundary-pair intersection lines, checking both orientations. Exactly the four adjacent pairs have feasible orientations; all four are future timelike and satisfy all 36 inequalities with precisely their nine incident zeros. The two opposite-pair lines are spacelike and neither orientation is feasible. The sum of the four finite rays is a strict relative-interior timelike point.

The four boundary covectors span the three-space, so their full-dimensional feasible cone is pointed. Every extreme ray therefore occurs among the six checked intersections. Its four extreme rays are precisely the verified finite rays. The projective section is a compact quadrilateral, and all other chamber inequalities hold throughout its positive hull. This proves closure of the entire F-face, rather than merely four isolated finite corners.

## Renewed escape and the support obstruction

Take the old 32-wall exterior ray at index 9, active on

\[
[0,1,2,4,5,6,7,10,27].
\]

Reconstruct V and normalize (V,H)=−1. Its norm is positive and every pairing with the 36 walls is nonpositive. The exact point

\[
E_{36}=V+H/8192
\]

has positive norm, is future oriented relative to H, and is strict against every wall. Consequently an open patch of the future ideal boundary remains feasible, proving infinite volume.

The orthogonal projection of E36 to F-perpendicular is future timelike and strict against all four quadrilateral sides. The complete face-cone calculation therefore places it in the positive hull of the four finite vertex rays. This yields an all-real stopping rule: **no coherent acute additional wall perpendicular to F can cut E36**, regardless of its coefficients or contact labels.

To justify vertex preservation, let h be any strict future timelike interior point of the proposed extension, and v an old finite vertex. Its active Gram matrix G is positive definite with nonpositive off-diagonal entries, hence G inverse is nonnegative. Decompose h=a v+Σ b_i n_i, with a>0 and b=G inverse times ((h,n_i))≤0. For an added acute normal u, (u,h)<0 and (u,n_i)≤0 imply (u,v)<0. Thus u is negative on the positive hull of the four face vertices. If u is perpendicular to F, its pairing with E36 equals its pairing with that projection and is negative. Further escape-cutting construction must leave this right-contact support.

This obstruction concerns the single specified support. It does not rule out arbitrary compatible cap clusters or changed retained walls.

## Verification and global scope

```sh
python3 scripts/fresh_nonarithmetic_escape35_20261016.py
python3 scripts/fresh_nonarithmetic_escape35_20261016.py --verify
```

Generation and fresh saved replay pass. The [independent review](fresh_nonarithmetic_escape35_20261016_review.md) also passes generation and fresh saved replay, independently reconstructing the new normal, all 35 contacts, 36 facet witnesses, all 100 finite vertices and 266 edges, all six quadrilateral intersection lines and the renewed escape. The full 36-wall cone has its own generation and replay records, with arithmetic shared with this producer; its independent complete-cone review is separate. No frozen predecessor was edited. This source, data and report are final and frozen.
