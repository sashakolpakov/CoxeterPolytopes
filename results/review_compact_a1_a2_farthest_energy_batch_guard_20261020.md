# Independent review of the A1/A2 farthest-vertex exclusion

20 October 2026. Reviewer: batch_guard. Paper/source review only: no
execution, imports, enumeration, optimizer, build or formal verification.

**Conclusion: PASS within the stated pointwise scope.** No mathematical
flaw was found in the complete proof
`results/compact_a1_a2_farthest_vertex_energy_exclusion_20261020.md`,
SHA256
`3bcf38496cd0aa7ba922dd13101595a479c0c14937d5c08bb368e2df435707e2`.
The theorem excludes an A1^a+A2^r whole diagram at the vertex farthest
from an interior point in any compact Coxeter n-polytope, n>=14.
It assumes no arithmetic hypothesis. It does not exclude other finite
vertex diagrams and therefore does not prove a new dimension ceiling.

## Interfaces independently checked

The farthest-vertex signs agree with
`results/compact_farthest_vertex_energy_obstruction_20261020.md`:
tau_i<=x.u_i, ||x||<1. Orthogonal component energies add. In an A2
component both the individual squared constraints and
(tau_1+tau_2)^2/3 are valid on every active-set branch; the proof never
misuses the unconstrained quadratic when an inverse coefficient is
negative.

The shadow bound is valid only for DISTINCT outside walls. I
independently rederived the inverse-Gram Schur identity, its entrywise
nonnegative remainder, and tanh L>=t_i/(A sqrt(d_i)). The uses in the
dotted isolated and N7..10 cases correctly establish that the known
wall cannot close the mate. The same mechanism is recorded separately
in `results/compact_actual_wall_reciprocal_shadow_20261020.md`.

Actual finite endpoint recognition is justified by projection onto the
orthogonal complement of a positive definite set of retained walls.
The inverse-Gram coefficients are nonnegative, so acuteness preserves
every other halfspace. The applications to the finite isolated,
F=(A2_3,A2_5), and N5 mate directions therefore give actual endpoints,
not just abstract finite diagrams.

## Inventories and compatibility

I independently checked both release-cost windows. For an A1 short
direction they are

    1-s<Q<1-(2/3)s;

with an A2 contact the only surviving profile is M_m=(A1_m,A2_3),
finite m>=7. The stated boundary candidates are eliminated by the
displayed positive affine null diagrams, including the exact equality
cases. The original A1-only inventory is unchanged by orthogonal A2
spectators. Its compatibility argument counts short DIRECTIONS,
including the two directions of a single (3,7) or (4,5) wall.

For the mixed M profiles, a common A1 contact or a same-node A2 contact
forces positive product; disjoint A2 blocks give product in (-1/2,0).
Only opposite nodes of the same A2 block can pair, so at most two M
walls coexist. An M wall cannot coexist with an A1-only short wall.
These arguments justify the claimed two-direction capacity without
assuming the two directions have distinct opposite walls.

For A2 short directions I rederived

    max(0,1-(4/3)B0)<Q<1-(200/243)B0.

The low native sums listed for labels3 and4 are exhaustive, and every
listed candidate contains the stated affine path or star. Label5
leaves exactly the A1_3 and A2_3 recipients. Label6 contains the old
affine3--6 path. Higher labels leave only A1_3, with B0<729/800;
labels7,8,9,10 satisfy the necessary cutoff and11 or higher do not.
A double contact in an A2 block is either affine or has shortest
sinh-square at least4sqrt2/9>31/50. Thus (11) is complete, including a
short mate release whose contacted release is shorter still.

## Energy and final arithmetic

The rational comparisons at tau*, tau0 and the isolated7 directions
check. In particular the two A1-short directions have squared
half-tangent sum>2/25 in every surviving case, giving
E_A1>a/10-3/25. In the finite/dotted two-wall case one uses the
known common-white factor delta=3/4 as well as the height bounds.

For finite isolated labels>=8 the actual mate contributes>4/25.
The exact F mate has cosh L=phi^3 and squared half-tangent1/phi.
The N5 mate has cosh L=phi^2/sqrt3 and squared half-tangent>1/5.
For dotted isolated and N7..10 profiles the shadow estimates give
tanh^2 L_mate>16/25 and hence squared half-tangent>1/4.
All are strictly greater than e=319/2000.

The only deficient block is isolated7, contributing>9/100. Two such
walls are impossible both within one A2 component and in distinct
components. Its height also forbids EVERY A1-short profile: disjoint
support gives a forbidden product in (-1/2,0), while an M wall sharing
that block gives positive product even at the smaller opposite-node
overlap. Thus the two exceptional debits are not charged together.

The final lower bounds check exactly. With isolated7, at n14 the
uniform lower bound is319*14/4000-7/100=2093/2000>1.
Without it, n>=15 is immediate. At n14 and a>=2 the least displayed
bound is1037/1000>1; when a=0 the A2 energy exceeds2233/2000>1.
These contradict the SAME farthest vertex's ||x||^2<1, rather than an
average over different vertices or separately chosen extrema.

The inherited geometric inputs are compact simplicity, all-wall
acuteness, the native Coxeter gap, and compactness exclusion of positive
affine null combinations. The report states these inputs and supplies
the local computations. No whole-polytope finite-volume enumeration,
published dimension bound, or unproved facet-descent assumption is
used.
