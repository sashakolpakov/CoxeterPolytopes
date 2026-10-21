# Root review: regular inball and the 120-cell volume

20 October 2026. Complete independent written review of
`compact_ra4_regular_inball_volume_comparison_20261020.md`, SHA256
`6b915db5b15a3cbcec6305eb76cfc5d7ba08b5c9793d4a64a6d91518760ebeaa`.
PASS as a conditional geometric theorem. No numerical execution,
manuscript build or Lean verification was used.

The spherical comparison is proved, rather than imported. For a positive
definite four-by-four unit-diagonal Gram matrix with off-diagonal entries
at least c, the all-ones Rayleigh quotient gives lambda_max>=1+3c.
AM--GM on the other three eigenvalues and the decreasing function
lambda*((4-lambda)/3)^3 then bound its determinant by the regular
determinant. On the positive coefficient simplex, its quadratic form is
pointwise at least the regular quadratic form. Both inequalities have
the required direction in the exact cone-volume integral.

I independently checked the integral normalization: the radial Gaussian
integral in dimension four is two; writing the positive coefficients as
s*a contributes Jacobian s^3 and radial integral 2/(a^t*G*a)^2. These
factors cancel, leaving precisely the displayed square-root determinant
times simplex integral. If any pair product exceeds c, the quadratic
form comparison is strict throughout the simplex interior. The equality
statement therefore also follows.

At an actual finite right-angled vertex, projecting the four mutually
orthogonal Lorentz normals to the Euclidean space perpendicular to the
ball center gives spatial Gram I+t*t^t before normalization. The ball
hypothesis bounds every t_i from below, so the resulting unit Gram has
all pair products at least phi/2. These are the actual Euclidean outward
normal cones of the bounded Klein polytope. Their interiors are disjoint
and they cover S3; no assumption about a prefix of walls is used.

The classical regular comparison model has 120 dodecahedral facets and
600 vertices. Its intrinsic dodecahedral height h^2=phi/2 implies center
height T0^2=phi^3, so its 600 congruent cones give the comparison angle
2*pi^2/600. Consequently the unknown polytope has at least 600 vertices.
The pinned compact RA4 source supplies W=v-5f>=0 and the exact volume
formula. Thus v-4f=(v+4W)/5>=120. Equality requires both v=600 and W=0;
the previously reviewed zero-excess recognition supplies the claimed
isometry. These uses of the regular model and recognition are explicit
inputs, not conclusions of the cone comparison.

The hypothesis is essential to the present argument: P must contain a
ball of radius asinh(phi^(3/2)). An inequality about its farthest wall
does not supply this ball. No universal minimum, second-volume theorem,
dimension improvement or unrestricted decomposition is proved. The
source's opening “awaiting independent review” records its creation
status; this hash-bound review records the subsequent PASS.
