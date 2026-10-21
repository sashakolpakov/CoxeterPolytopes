# Source audit of the orthogonal-component opposite obstruction

20 October 2026. Paper/source audit only; no execution, enumeration,
import, build or formal verification.

**PASS.** I read the complete source
`results/compact_farthest_vertex_orthogonal_component_opposite_obstruction_20261020.md`,
SHA256
`4216d06c2e03fd336898b8dac289b9969fa096c8a1d08bd67ddeb8ab5ec2f23d`.
It correctly proves that one actual opposite wall at the vertex farthest
from an interior point cannot close edges in two orthogonal dual
directions, hence cannot close edges in distinct whole finite
components. It assumes neither arithmetic nor one contact per component.
It proves no dimension ceiling by itself.

I supplied the short Schur proof used in this source, so this record
is an audit of its written hypotheses, translation and scope rather
than a third independent origin of the argument. Root proposed the
separate geometric right-triangle proof, which I also checked.

The correction coefficients are exactly
alpha=B_W^(-1)*beta_W with beta_W=-c.n_W>=0; their sign is nonnegative,
so c_F=c+sum alpha*n_W is orthogonal to W and satisfies z.c_F<0.
This uses STRICT interiority of z and all-wall acuteness. It does not
claim c_F is a facet, integral vector or unit normal. No reflection
group is assigned to the projected face.

The old dual vectors u_i,u_j are orthogonal to W and to p. When their
mutual product is zero, they and p have diagonal Gram(-1,1,1) and
span W^perp. The two actual endpoint equations determine
c_F=A*p-(A/r_i)*u_i-(A/r_j)*u_j with A=-p.c>0. Thus the interior
inequality is exactly x_i/r_i+x_j/r_j<1, with no omitted spatial
component or sign change.

The SAME farthest vertex gives x_i>=tau_i and x_j>=tau_j.
The exact identity tau/r=(1+tau^2)/2 makes their sum strictly greater
than1. Both edge lengths are positive and finite, so no zero-length,
ideal-endpoint or equality case is lost. The argument also covers an
empty W without a correction term.

The component interpretation follows from block-diagonal inversion
and strict inverse positivity on each irreducible finite Coxeter
component. The optional geometric explanation is correct: the two
old outgoing sides and their common opposite line bound the exact
triangle in the actual two-face. Feasibility of all three vertices
and convexity give the contained triangle; the three actual wall
halfspaces give the reverse inclusion. The main algebraic proof
does not depend on this explanation or an additional projection lemma.

The theorem is a new restriction at a selected farthest vertex, not
at every vertex. Same-wall multiplicity within one component remains
possible and still needs native profile/length control for a global
bound. The source states this limitation explicitly and does not
substitute a formal positive-definite corner for an actual endpoint.
