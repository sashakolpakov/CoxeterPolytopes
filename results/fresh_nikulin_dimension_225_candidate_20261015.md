# A research proof candidate for the finite-volume bound 225

**Candidate theorem.** Every finite-sided finite-volume real hyperbolic Coxeter polytope in H^n satisfies **n≤225**. No arithmeticity, compactness, simple-ideal-vertex, facet-count, or prescribed-diagram hypothesis is imposed.

This is a research proof candidate assembled from two explicit arguments and their internal independent checks. It is not a claim of publication priority, an externally peer-reviewed result, or the strongest currently known bound. Nikulin's 1986 ICM report already announces the unrestricted inequality n<300; the present work reconstructs explicit constants rather than treating 995 as the only historical assertion.

The [dependency checker](../scripts/fresh_nikulin_dimension_225_candidate_20261015.py) and [certificate](fresh_nikulin_dimension_225_candidate_20261015.json) bind the exact source, data, and proof reports. They verify the case partition and scalar comparisons, and recursively validate the saved source manifests. They do not purport to turn the geometric arguments or the cited classifications into a formal proof-assistant development.

## The two cases

Choose a narrow facet F supplied by Nikulin's narrow-facet theorem. Hyperplanes of any two framing facets whose hyperplanes meet H(F) have distance strictly less than arccosh 7. The proof below uses framing facets sharing actual ridges with F, a subset of that published framing set. With outward normals normalized to square two, such framing normals x,y cannot satisfy −(x,y)≥14. [Nikulin, printed p.655, Theorem A](https://www.mathunion.org/fileadmin/ICM/Proceedings/ICM1986.1/ICM1986.1.ocr.pdf).

**If F has an ideal vertex, use the product-cusp argument.** Its horospherical link is a product of irreducible affine Coxeter simplices with total rank N=n−1. Select two outgoing edges by deleting one node in each affine component, avoiding the node of F. The terminating walls can be chosen to meet F in actual ridges, including when the second endpoints are ideal and nonsimple.

At either endpoint the positive-semidefinite Schur bound is at most two. Inverse diagonals are at least one half, inverse entries are nonnegative, and every nonright norm-two Coxeter contact has magnitude at least one. Thus at most four retained contacts occur per endpoint. At most eight affine factors are affected; allowing one additional factor for the avoided narrow-facet node gives at most nine marked factors.

The [full product proof](fresh_nikulin_narrow_20261015_products.md) proves a uniform affine energy bound 2r/9 per factor of rank r, with loss at most two per marked factor. The two null-height ratios agree across all affine factors, without any arithmetic or equal-scale assumption. If M is the larger of this ratio and its reciprocal, then M≥1 and the exact Lorentz identity yields

    −(x,y) ≥ M ((n−1)/9 − 11).

For n≥226 the right side is at least 14, contradicting narrowness. Hence this case has n≤225. The [classical formulas](fresh_nikulin_narrow_20261015.md), [product arithmetic](fresh_nikulin_narrow_20261015_products.json), and [separate independent arithmetic review](fresh_nikulin_product_independent_20261015.md) make the all-rank and exceptional constants explicit.

**If F has no ideal vertex, it is compact.** This uses finite-sidedness and finite volume: the projective closure of the polytope is the convex hull of finitely many finite and ideal vertices, and the vertices of a facet are among these. A facet with only finite vertices is therefore compact. It is simple because its vertices are finite vertices of an acute Coxeter chamber. It need not be a Coxeter polytope itself.

Apply the separate [compact-facet lemma](fresh_nikulin_narrow_20261015_compact_facet.md), which gives n≤143. Its weights are defined using the ambient spherical vertex diagrams. Minimal indefinite subdiagrams are allowed to be quasi-Lannér, of order at most ten. Compact triangles supply two positive-weight corners; compact quadrilaterals supply one by a white path of length at most 17. At each ambient finite vertex, at most 17n pairs have that graph distance. The simple-polytope h-vector inequality on the facet contradicts n≥144.

Thus the two exhaustive cases give n≤max(225,143)=225.

## What was repaired during the audit

The first proposed treatment of the compact-facet case attempted to transfer Vinberg's finer compact-chamber weights directly. That would have used a false intermediate assertion: a compact marked face need not have an ambient star free of parabolic subdiagrams. A Coxeter tetrahedron with path labels 3–3–6 supplies an explicit counterexample.

The final compact-facet proof does not use that assertion. For a compact quadrilateral, the full block consisting of its white support and an opposite pair of black nodes is indefinite: a positive-definite block would give an actual intersection of opposite sides, while a singular positive-semidefinite block would give an ideal point on the compact quadrilateral. Minimal-indefinite extraction then permits quasi-Lannér diagrams. The original primary sources state the needed order-ten classification and extraction, and the proof separately explains why disconnected affine degeneracies cause no exception. [Prokhorov, §2, printed pp.407–408](https://www.mathnet.ru/php/getFT.phtml?jrnid=im&paperid=1484&what=fullteng).

The older frozen product report describes its result as conditional because its all-finite alternative was unresolved at that checkpoint. That report remains unchanged. The separate compact-facet proof supplies the missing branch here; it does not retroactively import an unproved compact-facet Coxeter property.

## Evidence and literature status

The product review reconstructs the exceptional affine matrices and their energy constants with separate arithmetic, checks independent classical coordinate/localization examples, and verifies the all-rank scalar inequalities. The separate [hostile mathematical review](fresh_dimension225_hostile_20261015.md) checks common null-height scaling, the contact budget, predetermined omission choices, and the actual-ridge endpoint bridge. It also audits the compact-facet source and geometry, with independent scalar replay. The combined certificate records the specific review files and their hashes.

The primary literature inputs are the narrow-facet theorem, spherical and affine Coxeter classification, quasi-Lannér order at most ten, the spherical graph-distance pair bound, and the simple-polytope h-vector relations. The new geometric reductions and estimates are written out in the bound reports; the checker does not silently substitute a finite numerical search for any all-rank statement.

Nikulin's printed pp.657–659 outline the narrow-facet method, work out an affine-A example, and announce n<300. The [earlier literature audit](dimension_bound_prokhorov_20261015_proof_review.md) distinguishes that announcement from the fully explicit 995 proof and its independently audited 990 refinement. This combined candidate is presented for further mathematical scrutiny with that historical distinction preserved. No current-best or novelty assertion is made.

```sh
python scripts/fresh_nikulin_dimension_225_candidate_20261015.py --verify
```
