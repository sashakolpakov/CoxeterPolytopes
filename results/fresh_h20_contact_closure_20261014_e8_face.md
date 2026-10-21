# An exact H12 completion problem inside the H20 hybrid

The frozen 76-wall replacement has an actual A12 + E8 corner. Its E8 face is currently an infinite-volume H12 chamber bounded by exactly the twelve A12 walls. **Any finite-volume completion retaining these eight E8 facets must add some walls orthogonal to the whole E8 support.** This is a necessary contact condition, not an obstruction to every completion.

The [exact checker](../scripts/fresh_h20_contact_closure_20261014_e8_face.py) and [certificate](fresh_h20_contact_closure_20261014_e8_face.json) bind the [76-wall construction](fresh_h20_contact_closure_20261014_replacement.json). Generation and saved replay pass. They reuse its exact scalar arithmetic. The [separate review](fresh_h20_contact_closure_20261014_e8_face_review.md) independently reconstructs the projections, intrinsic form, reversed source caps and source isometry.

All indices below refer to that 76-wall certificate. The E8 support is [13,14,15,16,19,27,34,68], and the ordered orthogonal A12 chain is [64,1,2,3,4,5,6,7,8,9,10,11].

The E8 Gram has diagonal two and edges minus one. Its inverse is entrywise positive, with diagonal [2,6,12,20,8,30,14,4]. Projecting all 76 roots to the perpendicular Lorentz space gives eight zero projections, twelve A12 roots, and **56 strictly future-timelike projections**. These last roots impose no boundary on the H12 face.

The geometric face theorem applies: a spherical subdiagram without an A-type or D5 component defines a Coxeter face. See Borcherds' theorem as stated in [Allcock, Theorem 2.1, printed page 743](https://arxiv.org/pdf/0903.0138). The calculation here also checks the particular projected normals directly.

## The all-label contact gate

Let r be a new acute-compatible Coxeter normal with square q>0, and write c_i=−r·e_i≥0 against the norm-two E8 roots. If it is not orthogonal to E8, some nonzero coefficient satisfies c_i²≥q/2: every finite nonright label is at least three, while parallel and dotted relations give a stronger bound. Positivity of the inverse E8 Gram B gives

\[
c^{\mathsf T}B^{-1}c\ge2c_i^2\ge q.
\]

The projection of r to the face has nonpositive square. Its orientation is essential: use the **proposed completion's own strict interior point**, projected onto E8-perp. The nonnegative inverse and acute pairings make its pairing with the projected new normal strictly negative. Thus the new projection is future causal and cannot create a face facet.

A finite-volume finite-faceted hyperbolic polytope has finite-volume faces. Additions all nonorthogonal to E8 therefore cannot complete this infinite-volume face. Arbitrary fields and finite labels are permitted. E8 can occur inside affine cusp diagrams, so a completed face need not be compact.

## Exact intrinsic form

For the saved finite apex V, set

\[
T=\frac{3V}{4628+2303\sqrt{79}}.
\]

The checker verifies

\[
V^2=-\frac{(4628+2303\sqrt{79})^2}{117},
\qquad T^2=-\frac1{13}.
\]

In the A12-root basis followed by T, the form is

\[
\operatorname{diag}\!\left(C_{A_{12}},-\frac1{13}\right),\qquad\det=-1.
\]

Its ambient embedding lies in Q(√79); all 76 reduced projection coordinates are archived. This rational intrinsic form does not imply arithmeticity of the ambient chamber.

The saved negative sum Y of A12 fundamental weights satisfies Y²=182, Y·a_i=−1 and Y·V=0. Hence V+sqrt(−V²/182)Y is future null and strictly satisfies all face walls. An outward perturbation proves infinite volume.

The [complete fourteen-wall I12,1 chamber](fresh_h12_a12_source_20261014.md) supplies a concrete smaller completion. Its [bounded transplant test](fresh_h20_contact_closure_20261014_h12_transfer.md) initially fails ambient gates, but a [one-wall replacement](fresh_h20_contact_closure_20261014_h12_replacement.md) completes the H12 face while leaving an ambient H20 escape.

    python scripts/fresh_h20_contact_closure_20261014_e8_face.py --verify
