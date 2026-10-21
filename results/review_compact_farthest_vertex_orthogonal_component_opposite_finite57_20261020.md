# Independent review: opposite walls at a farthest vertex

20 October 2026. **PASS**, written mathematical/source review only.
No mathematical program, enumeration, Lean or TeX job was run.

Reviewed in full:
`results/compact_farthest_vertex_orthogonal_component_opposite_obstruction_20261020.md`,
SHA256 `4216d06c2e03fd336898b8dac289b9969fa096c8a1d08bd67ddeb8ab5ec2f23d`.
This reviewer did not participate in deriving that source.

The signs and actual-wall interface check independently. The edge points
`p_i=cosh(L_i)p-sinh(L_i)u_i` give
`-z.p_i=E(cosh(L_i)-sinh(L_i)x_i)`, so maximal vertex distance gives
`x_i>=tanh(L_i/2)`, with the stated direction of inequality. A shared
opposite wall is nonincident to p, hence `A=-p.c>0`.

For the retained actual walls, all-wall acuteness gives `beta_W>=0` and
the positive-definite Gram with nonpositive off-diagonal entries has
nonnegative inverse. Thus `c_F=c+sum(alpha_s n_s)` is the orthogonal
projection to the retained-face space, with **positive**, not negative,
correction coefficients. Strict interiority of z gives `z.c_F<0` even
though this projected vector is not itself an actual wall normal.

When `u_i.u_j=0`, the face space has orthonormal Lorentz basis
`p,u_i,u_j`. The actual endpoint incidences force the coefficients
`A,-A/tanh(L_i),-A/tanh(L_j)`. Therefore `z.c_F<0` gives the strict upper
bound `x_i/r_i+x_j/r_j<1`. The two farthest-vertex inequalities instead
give at least `1+(tau_i^2+tau_j^2)/2>1`. Both strict signs have the needed
actual hypotheses: z is interior and both edge lengths are positive.
No hypothetical farther point or projected point is assumed to belong
to P.

The component consequence follows because the inverse finite Gram is
block diagonal and is strictly positive within each irreducible block.
The optional two-face interpretation also checks: the convex hull of
the three actual vertices lies in that face, while its three actual
supporting half-planes contain the face in the same triangle. The
algebraic proof does not need this recognition step. The empty retained
set is harmless.

The standard inputs remain explicit: finite vertices of a compact
Coxeter polytope are simple, the all-wall acuteness property applies to
its actual facet normals, and positive-definite matrices of this sign
pattern have nonnegative inverse. The result is conditional on choosing
an actual vertex farthest from an actual interior point; it does not
assert that every vertex admits such a point. No arithmetic restriction,
number of directions per component, or dimension ceiling follows.
