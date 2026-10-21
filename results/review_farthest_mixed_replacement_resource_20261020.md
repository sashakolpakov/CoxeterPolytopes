# Outside-author review: actual mixed replacements

**PASS, source-only.** I read the complete 6,102-byte source
`results/compact_farthest_vertex_mixed_replacement_inequalities_20261020.md`,
SHA256 `3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84`,
and independently rederived its projection, replacement and correlation
inequalities. No mathematical correction is requested. This is not a
scientific replay, cap classification or dimension theorem.

I also read its farthest-energy interface
`compact_farthest_vertex_energy_obstruction_20261020.md`, SHA256
`50c64e5a521d3e4e0f277a8ece02d3a29aa4ae8b71485cf3088211e587f45bff`.

## Actual geometry and the all-rank condition

For a positive-definite actual-wall Gram G, the disjoint-positive-parts
argument correctly proves `G^{-1}h>=0` for `h>=0`. With
`q0=p+M G^{-1}h`, direct expansion gives

`M^T q0=0`, `q0^2=-(1+D)`, `p.q0=-(1+D)`.

The last identity fixes the future time orientation. Every other actual
wall has nonpositive product with q0 because both the added coefficients
and its old wall contacts have the required signs. This proves actual
feasibility, not just an abstract spherical subdiagram. If the number
of selected independent walls is n, their hyperbolic intersection has
dimension zero and the feasible point is an actual finite vertex. No
intrinsic Coxeter property of a facet is being assumed.

Compactness correctly bounds every point by the farthest vertex:
if `y=sum alpha_v p_v`, with nonnegative coefficients summing to one,
then `-y^2>=1` by the reverse Lorentz Cauchy inequality, and
`-z.(y/sqrt(-y^2))<=E`. This proves the comparison used for projected
points of arbitrary selected rank in the compact theorem.

Writing a replacement as `(p-v)/sqrt(1-R)` gives exactly
`x.v>=1-sqrt(1-R)`. Since `|x|<1` and `R<1`, the quantity `s=x.v`
satisfies `s<1`, which is essential when squaring to obtain
`R<=2s-s^2`. The zero replacement is separately covered.

## Isolated caps and the weighted certificate

For distinct actual isolated caps, their spatial directions are the
unit duals u_i, not the old simple normals. Their Gram is therefore
`diag(b)(K_J-r r^T)diag(b)`. Congruence by `K_J^{-1/2}` proves the
positive-definiteness criterion `r^T K_J^{-1}r<1` in every rank.

The retained old walls are orthogonal to all these caps. Together they
form n independent actual walls with positive-definite Gram. Section 1
therefore supplies a finite vertex for every admissible J. Solving its
orthogonality equations gives precisely
`v_J=sum (K_J^{-1}r)_i u_i`; negative coefficients cause no problem.
Taking a singleton recovers the ordinary incident-edge inequality.

For nonnegative weights, a positive weighted sum of the right sides
forces `x.V>0`, hence V is nonzero. Cauchy and `|x|<1` then give the
strict inequality (8). Native Coxeter products are still an independent
necessary condition, as the source explicitly states.

## Exact check of the B20 diagnostic

The fractions `700/821` and `12600/13969` follow from the stated
half-tangents. In the B_n ordering whose dual directions are
`(e1+...+ei)/sqrt(i)`, `K_12=1/sqrt(2)`.

Using only the rational intervals printed in the source, I obtain

`0.904883152838 < R < 0.906837105704`,

`0.6170722089 < s < 0.6177879316 < 0.618`.

Since `2s-s^2` is increasing for `s<1`, it is below `0.854076`.
Thus the claimed margin greater than `0.05` is valid, and R is also
strictly below one, so the two-cap replacement test applies. This
diagnostic uses the stated pair of saturated edge constraints; its
failure does not require the old witness's separate norm estimate.
The source correctly says that the formal caps also fail the native
Coxeter alphabet. It does not call this a realized polytope.

## Useful finite-volume scope, separate from the compact statement

The projection proof itself does not use compactness. Suppose an
arbitrary finite-sided, finite-volume acute P has at least one finite
vertex. Choose an interior future unit z and let p maximize `-z.p`
over **only the finite vertices**. This is a finite nonempty set and
still gives `|x|<1`.

Every full-rank positive-definite actual-wall replacement is a finite
vertex. Its farthest comparison, the isolated-cap inequalities (7),
and the weighted inequality (8) consequently remain valid for this p,
regardless of how many cusps P has. This bypasses the multicusp causal-
height barrier for these finite targets.

It does not justify a comparison with arbitrary lower-rank projected
points, nor with ideal endpoints. An entirely ideal polytope has no
finite vertex to select. Thus the outstanding global task is adequate
positive-definite replacement coverage, or a controlled loss for ideal
continuations; no such coverage follows merely from these inequalities.

**Review boundary:** complete written mathematical read and elementary
hand arithmetic. No candidate mutation, executable import, scientific
job, build, formal proof or independent census was performed.
