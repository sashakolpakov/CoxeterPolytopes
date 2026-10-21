# Root review: determinant-three golden G2 edge and full facet

20 October 2026. Full paper review; no mathematical worker, build or
formal verification was run.

Reviewed source:
`results/golden_determinant3_g2_first_edge_and_full_facet_20261020.md`,
SHA256 `24870d50f72e768ea3dc4734680c49689b462da56af0850d6bb58e05e9dfa51e`.

**PASS for the complete first edge and full-facet identity.** This
review does not establish nonreflectivity or a new example.

For L_n=<-phi>+I_n+A2, n>=7, I checked the primitive norm pool
{1,2,3,6}, including the inert-prime divisibility and positive-unit
normalization. The time-zero coefficient equations exclude additional
mixed roots and yield the full B_n+G2 germ. The finite continuation
classification gives all ten profiles, including the two possible
B_n extension rows.

The edge equation is 3Z^2-2phi*x^2=R. Intersection before tau=1/2 gives
x^2 <= (3phi^2/2)R, and conjugate positivity gives the other embedding
bound. I rechecked the F9 square list, every residue exclusion, and the
remaining coefficient strips. They give exactly the source's positive
candidates. The q=6 remnants fail respectively the squareclass of 2
and the rational norm 76. The sole q=3 candidate reconstructs c.

The covector of c and c^2=3 make projection integral. In its orthogonal
complement, z_2=x and the binary Gram is [[2-phi,-1],[-1,2]]. The two
displayed vectors have norms 1 and -phi, are orthogonal and have unit
determinant, proving an integral diagonal identification.

I also checked fullness, rather than only the lattice identity. For a
crossing root v+kc, q<=3 excludes nonzero k by the integral field norm.
For q=6 the two strict conjugate bounds leave only k=+-1; reflection
integrality then forces v into 3M, incompatible with v^2=3. Orthogonal
primitive q=3 or 6 roots are excluded in the unimodular complement,
and every intrinsic norm1/2 reflection extends to the ambient lattice.
Thus the actual c facet carries the full <-phi>+I_(n+1) arrangement.

The determinant distinction uses the K-valued squareclass of 3,
certified by its odd valuation, not its square rational norm. The note
correctly leaves the exact published cutoff for the literal golden
diagonal family as an unresolved external-input check. Neither an even
E8 example nor a commensurability claim replaces that missing statement.
No uniform dimension bound or whole-class exclusion follows here.
