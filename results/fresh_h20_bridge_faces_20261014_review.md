# Two exact complete H2 bridge faces

The independent certificate [JSON](fresh_h20_bridge_faces_20261014_review.json) and [checker](../scripts/fresh_h20_bridge_faces_20261014_review.py) establish that each of the two specified codimension-18 faces of the frozen 95-wall chamber is the **entire hyperbolic plane**. This is a local geometric statement, not a finite-volume result or a new-cap certificate.

The supports are

```
S61 = [1,2,3,4,5,6,7,8,9,10,12,13,14,15,17,31,61,87]
S75 = [1,2,3,4,5,6,7,8,9,10,12,13,14,15,17,31,75,87].
```

The input is the hash-bound [95-wall certificate](fresh_h20_e7_edge_cusp_repair_20261014_exact.json). The review imports no project arithmetic: it uses sparse radical products over the 16-element multiquadratic basis, a further positive square root, rational matrix elimination, and outward rational sign enclosures. It reconstructs all 95 norms, the ambient inertia `(20,1)`, both support Grams and their inverses, and every old-wall projection. The old full pair census is an explicit inherited dependency.

## Actuality and absence of old boundaries

For either support, its rational Gram `B` is positive definite and `B^-1` is entrywise nonnegative. If `H` is the old strict future center, its orthogonal projection is

`H' = H - sum_i (B^-1 (H·s))_i s_i`.

The coefficients subtracted are negative. Hence `H'` remains strictly inside every nonsupport inequality; its square is negative and its time orientation is future. It is an actual relative interior point of the face.

More strongly, **all 77 nonsupport projected normals have strictly negative square and are future oriented**. Therefore their pairings with every future timelike vector in the residual Lorentz three-space are negative. None defines a geodesic boundary in that plane. The support projections vanish. Thus the whole `H2` belongs to the face, rather than only a polygonal region in it.

In particular, adjoining any one old wall to either support gives a negative Schur remainder, so cannot produce a positive-definite rank-19 extension. This conclusion is exact for these supports; the checker does not enumerate other PD18 supports.

## Orientation of possible separating caps

Let `v_i` be the projected old normals. The checker verifies that their sum `T'` is future timelike, and that `v_0,v_11,v_28` span the residual three-space. It also verifies that the old strict escape projects to a future timelike vector `E'` on each face.

Consequently, a nonzero vector `r` satisfying `r·v_i <= 0` for every old wall obeys `r·T' < 0`: equality would force every individual pairing to vanish and hence `r=0`. If it also satisfies `r·E' > 0`, it can be neither future nor past causal. It is necessarily spacelike and can be normalized as a cap normal. This removes a potential normalization or orientation gap in an exact separating-cone construction. Numerical infeasibility of a cone LP is not certified by this review.

The same sum supplies a common strict future face center for any finite set of such acute caps. A sufficiently small inward perturbation by the old strict center makes all support inequalities strict as well.

## What closure would require

Caps perpendicular to the support induce geodesic halfplanes in this complete `H2`. Any finite-volume completion of the ambient chamber must in particular give this retained face finite area. A finite-angle mate creates a corner, but one cap or a pair of caps does not close the whole plane into a finite-area polygon. The bounded numerical search for mates is separate from this certificate; no all-family angle obstruction follows here. Caps having nonzero support components are also outside this perpendicular-cap statement.

The elementary [acute-extension preservation lemma](acute_extension_preservation_20261014.md) was source/proof reviewed with no issue: positive normal squares, exact acuteness and a common strict future center make all listed facets actual, strictly retain old finite vertices, and strictly retain proper bounded full cusp links. Those assertions do not provide a global vertex census or finite volume.

## Reproduction and scope

```
python3 scripts/fresh_h20_bridge_faces_20261014_review.py --verify
```

Generation and saved replay pass. The saved JSON binds its checker and the old input source/data. No producer files were modified. This review certifies two specified entire-plane faces, their orientations and the separator implication; it does not certify a new cap, any numerical LP outcome, or completeness of the current PD18-face list.
