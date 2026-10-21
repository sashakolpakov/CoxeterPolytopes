# Root mathematical review of the small-three-face dimension theorem

18 October 2026. **PASS for the written universal proof.** If every
three-face of a finite-volume right-angled hyperbolic n-polytope has
facet-plus-cusp sum at most10, then n<=8. This does not change the
unconditional mixed ceiling12.

The complete [proof](right_angled_mixed_small_three_face_dimension_research_20261018.md)
and its inherited R7/Q8 degree classification were read. The R7 facet
data imply that the central finite vertex meets the three quadrilateral
facets, and each leaf meets a unique bad triangular facet. Its outward
normal is sinh(t_i)q-cosh(t_i)n_i. Two such bad facets share an ideal
vertex and are opposite walls in that cusp rectangle; product-1 gives
sinh(t_i)sinh(t_j)=1. All three sines equal1, fixing every finite edge
at cosh length sqrt2.

For an ambient finite edge, two ideal directions give an R7 leaf when
k>=2. With k<=1 and n>=4, two other finite directions give an R7 center.
Thus every finite edge has the same length under the carrier hypothesis.
Every triple of ideal directions must be P3. Its opposite finite vertex
is exactly2q minus the three actual orthonormal normals. Different triples
give different opposite vertices at q; reversing a P3 face gives the
same move back. The sum of the defined displacements therefore cancels
over actual finite edges and P3 faces, with their correct multiplicities.

The displayed Delta formula follows directly from C(k,3)=k*C(k-1,2)/3.
For k>=9 its terms are nonnegative and its time coefficient is positive;
k=9,m=0 is future null, not the zero vector. For k<=8 the nine base
cases and increasing finite differences cover every n>=9. A finite sum
of nonzero future causal vectors cannot cancel. At least one finite
vertex exists because each allowed three-face has finite vertices.
This proves the claimed contradiction.

The author records tiny exact sqrt2 base checks; this root review is
mathematical and does not claim a separate execution of their checker.
The earlier P3 attachment/classification facts remain explicit internal
dependencies. Full Lean verification and external refereeing are pending.
The result forces at least one three-face of facet-plus-cusp sum>=11 in
dimensions9--12. A quantitative fraction of such faces, sufficient for
a new unconditional dimension bound, has not been proved.
