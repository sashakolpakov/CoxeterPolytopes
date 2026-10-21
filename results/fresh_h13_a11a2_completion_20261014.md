# A complete H13 source with an A11+A2 corner

**2026-10-14 — exact construction and saved replay PASS.** The [standalone Fraction checker](../scripts/fresh_h13_a11a2_completion_20261014.py) and [certificate](fresh_h13_a11a2_completion_20261014.json) give a complete finite-volume H13 Coxeter polytope with **15 actual facets, 44 finite vertices and two bounded ideal vertices**. It has the required A11+A2 corner. This solves the intrinsic source problem left by the [87-wall H20 E7 face](fresh_h20_e7_face_20261014.md). It does not itself establish compatibility or finite volume of an ambient H20 extension.

The [independent review checker](../scripts/fresh_h13_a11a2_completion_20261014_review.py) and [review certificate](fresh_h13_a11a2_completion_20261014_review.json) pass. A separate enumeration using at most two positive gap coordinates reconstructs all 46 projective vertices, the two bounded links and the normalized transfer Gram.

## Source and actual facets

The construction takes the E8 face at finite-corner orbit 243 in the [complete 210-wall H21 corner census](fresh_h21_source_corner_walk_20261013.md). Its supporting E8 roots are

`[8,11,22,25,26,28,33,39]`.

Their norm-two diagram has branch arms of lengths 1,2,4. At the same corner, the other components are an A11 of norm-two roots and an A2 of norm-four roots. The latter has Gram `[[4,−2],[−2,4]]`; dividing its roots by √2 gives the usual norm-two A2 Gram.

Projecting all 210 old normals to the E8 orthogonal complement yields:

| Projection type | Count |
|---|---:|
| Positive norm, already orthogonal to E8 | 15 |
| Zero vector, the supporting E8 roots | 8 |
| Strictly future null, hence redundant in hyperbolic space | 2 |
| Strictly future timelike, hence redundant | 185 |

The 15 actual face walls, in certificate order, are

`[6,34,27,23,17,24,46,13,14,48,4,36,162,3,9]`.

The first eleven form an A11 chain; the next two form A2. All 105 unordered pairs are checked directly: 91 right, 13 at π/3 and one at π/4. A strict timelike center and its orthogonal projection to each supporting wall certify all 15 actual facets. The checker verifies the ambient signature through an explicit orthogonal E8, A11+A2 and timelike basis.

The choice of E8 face is motivated by Borcherds' theorem, stated in [Allcock, Theorem 2.1, printed page 743](https://arxiv.org/pdf/0903.0138): a spherical supporting diagram with no A-type or D5 component gives a Coxeter face. Here the complete source geometry is also certified directly, so its volume conclusion does not depend on importing the H21 volume theorem or assuming that a projected inequality is an actual facet.

## Exact global closure

The saved corner ray V has square −4. Use the unnormalized corner roots, with Gram `diag(A11,2A2)`, and let their dual weights be w_i. Every future point has positive projective coefficient along V. In the slice with that coefficient one, write

\[
p=V-\sum_{i=1}^{13}z_iw_i,\qquad z_i\geq0.
\]

The two additional facets give linear upper bounds on z. Their coefficients are nonnegative, and every coordinate occurs with a strictly positive coefficient in at least one bound. Thus this projective slice is bounded and its extreme points appear among the 105 choices of thirteen active inequalities from fifteen.

The checker exhausts those 105 choices exactly. After deduplication there are 46 feasible rays: 44 timelike and two null, with no spacelike ray. Each finite vertex has a positive-definite rank-13 incident Gram and satisfies every other inequality strictly. Its physical point is also checked against all 210 original source inequalities.

At the first null ray the affine components have 9 and 5 nodes, with quotient ranks 8 and 4. At the second there is one 13-node affine component of quotient rank 12. Every component has an explicit positive kernel, a positive-definite deleted-node minor and a positive vector sum proportional to the future cusp ray. The total quotient rank is 12 in both cases, so the full horospherical links are bounded. This establishes finite volume, rather than merely a bounded search prefix or finite vertex sample.

## Two transfer formulas

Order the source A11 chain as

`[6,34,27,23,17,24,46,13,14,48,4]`,

and the A2 pair as `[36,162]`, normalized by `1/√2`. Set T=V/√2, so T²=−2. In this frame the intrinsic form is exactly

\[
\operatorname{diag}(C_{A_{11}},C_{A_2},-2).
\]

Let λ_j, j=1,…,11, be the A11 fundamental weights and μ_1,μ_2 the normalized A2 weights. The two missing corner facets are

\[
R=-\lambda_{11}-\sqrt2\,\mu_2+\frac{\sqrt2}{4}T,
\qquad
S=-\lambda_3+\frac{\sqrt2}{4}T.
\]

Both have norm two and their mutual pairing is zero. Reversing A11 replaces indices 11 and 3 by 1 and 9; swapping A2 exchanges μ_1 and μ_2. These four corner identifications are the small finite set of direct transfer tests.

The certificate stores every source normal in this normalized frame, using coefficients `a+b√2`, and reconstructs all 225 ordered Gram entries from those coordinates. The original two A2 normals retain norm four in that list; they are √2 times the corresponding basis roots and define the same walls as their norm-two normalizations.

This is a geometric field-extension construction. The extra √2 is allowed; the source is not being claimed to be a chamber of the full integral reflection group of `A11⊕A2⊕<-2>`. The [separate integral-lattice gate](fresh_h13_a11a2_lattice_gate_20261014.md) remains valid. Arithmetic H13 sources are permitted; no noncompact arithmetic H20 construction is inferred.

## Search scope and ambient obligation

The preliminary [norm-two subgroup prefix](fresh_h13_a11a2_completion_20261014_scout.json) and [full integral-root prefix](fresh_h13_a11a2_completion_20261014_full_scout.json) are bounded explorations only. The full integral prefix starts with A11+G2, explicitly retaining the norm-six apex mirror, and orders candidates by `t²/q`. Neither prefix supplies a reflectivity or nonreflectivity theorem. The [diagonal-source scout](../scripts/fresh_h13_a11a2_completion_20261014_diagonal_scout.py) likewise makes only bounded root/corner claims. None of those prefixes is a premise of this complete face construction.

To use the source in H20, each transferred intrinsic facet must still pass every retained ambient Coxeter gate. The E7 gate allows an orthogonal lift or the unique single-π/3 lift at wall 33; it does not guarantee either lift is compatible. The source has been handed to the separate ambient worker for these tests and any explicitly checked core replacement. No ambient finite-volume conclusion is included here.

```sh
python3 scripts/fresh_h13_a11a2_completion_20261014.py --verify
```
