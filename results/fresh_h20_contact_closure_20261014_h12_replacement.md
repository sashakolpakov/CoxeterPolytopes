# A 77-wall H20 chamber with a completed H12 face

**Deleting wall0 from the frozen 76-wall chamber and adding the reversed I12,1 pair gives a coherent 77-facet H20 chamber. Its E8 face is a complete finite-volume H12 chamber. Ambient H20 volume remains infinite.**

The [checker](../scripts/fresh_h20_contact_closure_20261014_h12_replacement.py) and [certificate](fresh_h20_contact_closure_20261014_h12_replacement.json) bind the [76-wall predecessor](fresh_h20_contact_closure_20261014_replacement.json), [face coordinates](fresh_h20_contact_closure_20261014_e8_face.json), [transplant test](fresh_h20_contact_closure_20261014_h12_transfer.json), and independently constructed [H12 source](fresh_h12_a12_source_20261014.md). Generation and saved replay pass. This checker reuses the predecessor's exact scalar arithmetic; the [independent face review](fresh_h20_contact_closure_20261014_e8_face_review.md) separately checks the source transfer and strict center/escape signs.

The inserted roots are

\[
R=-\lambda_{10}+2T,\quad R^2=2;\qquad
S=-\lambda_2+3T,\quad S^2=1;\qquad R\cdot S=0.
\]

Their only failed predecessor relation was R against wall0. Removing it leaves 151 new unordered pairings: **39 right, one π/3, one π/4, and 110 dotted**. All 2,926 final pairings are checked. The old rational strict center still works; its projections onto all 77 hyperplanes are strict timelike facet points. Thus every inequality is an actual facet. The ambient field remains Q(√2,√79).

## The completed face

In new numbering the E8 support is [12,13,14,15,18,26,33,67]; the inserted roots are 75 and 76. The certificate checks an exact isometry between the fourteen orthogonal face walls and the complete I12,1 source. The remaining 55 nonorthogonal projections are future timelike and impose no face boundaries.

All 37 projective source vertices are lifted explicitly. Its **35 finite vertices become actual H20 finite vertices**, with positive-definite rank-20 incident Grams and all nonincident inequalities strict. This is not a complete ambient vertex census.

The two source ideal vertices are feasible null rays in H20, but their ambient links remain unbounded. Each has the source's rank-eleven affine link plus a finite E8 factor. A saved negative E8 fundamental weight of norm two is perpendicular to the null ray and nonpositive on every incident normal, giving a nonzero link recession direction. The finite-volume H12 face alone therefore does not close those ambient cusps.

## Global status

The exact predecessor escape survives:

\[
E^2=\frac{269344645863312077}{40000000000000000}>0.
\]

It has the correct future orientation and is strict against all 77 walls. Ambient volume is infinite.

The conditional non-quasi-arithmetic obstruction is rechecked after deletion. New indices 17 and 21 have Cartan two-cycle 12+8√2; indices 23 and 24 have rational Gram [[2,−6],[−6,2]], determinant −32. Any cofinite completion retaining these facets violates the nonidentity-embedding semidefiniteness condition because the rational indefinite block survives every embedding. Rationality of the H12 face form does not change this.

The separate [79-wall completion](fresh_h20_e8_cusp_completion_20261014.md) now repairs both ambient E8 cusp factors while preserving the face and the global escape. Its [independent review](fresh_h20_e8_cusp_completion_20261014_review.md) checks the positive affine kernels and full quotient ranks. This 77-wall certificate remains frozen.

    python scripts/fresh_h20_contact_closure_20261014_h12_replacement.py --verify
