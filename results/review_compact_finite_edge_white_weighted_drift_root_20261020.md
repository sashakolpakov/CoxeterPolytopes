# Root review: globally shared white-diagram drift

20 October 2026. Complete independent written review of
`compact_finite_edge_white_weighted_drift_20261020.md`, SHA256
`486ce698db8c705f27ca58c092f30b67e695680fcf2c7a3015641996d586b165`.
PASS. No computation, dimension theorem or Lean verification is claimed.

I independently derived the edge tangent, principal-cofactor factor and
hyperbolic half-angle identity in equation (1). The scalar coefficient
of an edge agrees at its two endpoints because its white diagram and
length are shared. Thus the sum cancels exactly over the whole finite
edge graph. Under the opposite of (2), every right-hand term is nonzero
future causal: A is positive, the dual-vector sum is tangent to the
future unit vertex, and its norm is at most A. Their finite sum cannot
vanish. This proves the strict global drift conclusion.

For A_n, the displayed inverse-Gram entries and normalized Markov kernel
are exact. Multiplication by K inverse gives the stated positive weights;
their positivity follows from strict concavity of sqrt(i*(N-i)), with
zero boundary values. Cauchy--Schwarz in the K inner product proves the
minimum ratio, and conversion from dual to simple normals independently
checks its norm and normalization. The product-kernel inverse quadratic
form gives equation (4). Its logarithmic asymptotic follows by summing
the logit increments; the sum of their cubic errors stays bounded.

The alternative determinant weight has exactly the stated common
vertex factor, and the reciprocal diagonal weights sum to H_n. Retaining
only the long directions in the positive sum A proves the strict short
weight inequality. No existence claim about a chosen short wall is
needed for this calculation.

The scope restrictions are necessary and correctly retained. Shared
white types impose the same f at A, B, D and disconnected corners. The
global drift vertex need not be A_n, and one closing wall may account
for several incident directions. No all-vertex profile bound has yet
been established, so this lemma changes no dimension ceiling.
