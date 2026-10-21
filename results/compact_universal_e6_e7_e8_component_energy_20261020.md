# Universal component energy from an extended D4 star

20 October 2026. Own paper proof, no computation, build or formalization.
Other finite components are arbitrary. These are component bounds, not
a new uniform dimension ceiling.

## Statement and actual-edge inputs

At an actual vertex of a compact hyperbolic Coxeter polytope, let T be
an E6, E7 or E8 component of its WHOLE finite vertex diagram. Use unit
outward normals, inverse old Gram C, dual unit directions
u_i=omega_i/sqrt(C_ii), and actual finite outgoing edge lengths L_i.
Put tau_i=tanh(L_i/2) and define the component covering energy by

    E_T=min { ||y||^2 : y in span(T), y.u_i>=tau_i for every i in T }.

Then, respectively,

    E_E6 > 9/20,       E_E7 > 5/9,       E_E8 > 2/3.       (1)

Each bound exceeds rank(T)/14. At a vertex farthest from a strict
interior point, the component energies add to less than1. This last
interface is the previously proved
`compact_farthest_vertex_energy_obstruction_20261020.md`, SHA256
`50c64e5a521d3e4e0f277a8ece02d3a29aa4ae8b71485cf3088211e587f45bff`.

The actual reciprocal-shadow input is
`compact_actual_wall_reciprocal_shadow_20261020.md`, SHA256
`c913ec3bd8e0f26a8c679a1ede2b98b251ad3f36efeae60486213ddfa3e8badc`.
It implies that DISTINCT actual opposites of directions i,j satisfy

    tanh(L_i)*tanh(L_j) >= K_ij,
    K_ij=C_ij/sqrt(C_ii*C_jj).                         (2)

Indeed the inverse-Gram Schur identity gives
T_j(c_i)>=K_ij/tanh(L_i), and the actual opposite of j, if different
from c_i, has tanh(L_j)>=T_j(c_i). No half tangent is substituted in
(2). If E_T<=e<1, each tau_i<=sqrt(e), and consequently

    tanh(L_i)*tanh(L_j) <= 4e/(1+e)^2.                (3)

Thus K_ij>4e/(1+e)^2 forces the SAME actual opposite for i and j.

A sharper version will also be used below. Write k=K_ij and
z=e(1+k)/2<1. For a covering vector of squared norm at most e,
Cauchy--Schwarz gives `(tau_i+tau_j)^2<=2e(1+k)=4z`, hence
q=tau_i*tau_j<=z. Since tau_i^2+tau_j^2>=2q,

    tanh(L_i)*tanh(L_j)
       =4q/[1+tau_i^2+tau_j^2+q^2]
       <=4q/(1+q)^2 <=4z/(1+z)^2.                  (3a)

The last function is increasing for0<=q<=z<1. Thus this stronger
upper bound being below k likewise forces the same actual opposite.

## An extended D4 star cannot share one opposite

Let S consist of a central old node o and its three simple neighbors
a,b,c, so the induced old diagram is D4. Suppose at least one of its
three arms extends to another old node. Then no actual outside wall
d can be the actual opposite of ALL four directions in S.

Here the explicitly imported classical classification input is:
minimal indefinite Coxeter diagrams with all proper subdiagrams
positive definite (Lanner diagrams) have order at most5. Compactness
also excludes actual affine principal subdiagrams. These are the only
non-elementary diagram-classification facts used in this argument.

To prove the claim, the full old finite vertex diagram together with
d spans the Lorentz space and has a nonpositive principal subdiagram.
Choose a minimal one. It contains d because all old walls are finite.
It contains EACH node of S: deleting any such node leaves the actual
finite replacement vertex and hence a positive definite Gram matrix.
The minimal subdiagram is not affine in a compact polytope, so it is
Lanner and has at most5 nodes. It therefore consists EXACTLY of S+d.
This step also handles arbitrary contacts to orthogonal spectators.

We now classify this D4+d Lanner diagram directly. If d contacts both
o and a leaf, it has a nonpositive proper triangle. If it contacts two
leaves, those leaves, o and d have a nonpositive proper four-cycle.
Hence d contacts either just o or just one leaf within S. Contact to
o alone must be3: a4 or larger already makes a nonpositive proper
four-node star. But the four simple leaves then give affine D4,
not Lanner. Contact to a leaf alone must be3,4 or5, because the
proper four-node chain with the two subsequent simple edges must be
finite (A4, B4 or H4). Contact3 gives finite D5, and contact4 gives
an affine diagram. Thus the ONLY possibility is a5 to one leaf q,
and no other contact to S.

This possibility cannot close the four actual directions in the
extended old diagram:

* If q itself extends to q', the four nodes d,q,o,q' contain the
  star with incident labels5,3,3. Its unit Schur complement is
  `1-cos^2(pi/5)-1/4-1/4<0`. It omits another node of S.
* If q does not extend, choose a different extending leaf r and its
  next node r'. The five nodes d,q,o,r,r' contain the chain
  `5,3,3,3`. In raw norm2 convention its consecutive determinants
  are `2,3-phi,4-2phi,5-3phi,6-4phi`; the last is negative. This
  chain omits the third leaf of S.

In either case deleting the omitted node of S therefore does NOT
leave a positive definite diagram, contradicting its actual finite
replacement. Additional contacts of d to q' or r' cannot repair the
obstruction: the base nonpositive connected Gram has a positive
nonpositive quadratic witness, and extra nonpositive products only
decrease that quadratic form. The argument never treats a general
nonpositive subdiagram as forbidden by compactness; it uses the
specific positive definiteness of the replacement missing that S-node.

## The E-type central correlations

Write the E6,E7,E8 trees as three arms of lengths `(1,2,L)` from o,
where L=2,3,4 respectively. Their central four-node S is D4 and has
an extending arm, so the preceding obstruction applies.

For completeness, the correlations below follow without a matrix
census. In RAW norm2 convention the central Schur complement is

    delta=2-[1/2+2/3+L/(L+1)],
    (C_raw)_oo=1/delta.

For the node next to o on an arm of length m, put t=m/(m+1).
Its diagonal inverse entry is `t+t^2/delta`, and its covariance with
o is `t/delta`. Thus its squared central correlation is

    K^2=t/(delta+t).

The minimum is on the arm of length1. For E6,E7,E8 the respective
delta values are1/6,1/12,1/30, giving

    min_(j adjacent to o) K_oj
       =sqrt(3)/2,       sqrt(6/7),       sqrt(15)/4.   (4)

Scaling every old root to unit norm does not change K.

Suppose the corresponding bound in (1) failed. At e=9/20,5/9,2/3,
the right side of (3) is respectively

    720/841,           45/49,           24/25.

Each is strictly smaller than its value in (4). The exact squared
comparisons reduce to

    2073600<2121843,    14175<14406,    9216<9375.

Consequently the opposite of o is also the actual opposite of each
of its three neighbors. This contradicts the extended-D4 obstruction
and proves (1), independently of all spectator types and contacts.

## D5 and D6 corollary

The same extended-star obstruction applies to D_r, r>=5, whose arms
have lengths `(1,1,r-3)`. The preceding arm inverse formula gives

    min_(j adjacent to o) K_oj=sqrt((r-2)/r).

Consequently

    E_D5 > 5/14,             E_D6 > 3/7.             (5)

For D5 use e=5/14 and k>=sqrt(3/5). Clearing the positive denominator
in `k-4z/(1+z)^2`, with z=e(1+k)/2, gives

    25k^3+330k^2+529k-560.

It is strictly increasing for k>=0. At k=sqrt(3/5) it equals
`544k-362>0`, as sealed by `887808>655220` after squaring positive
sides. For D6 use e=3/7 and k>=sqrt(2/3). The corresponding polynomial
is `9k^3+102k^2+121k-168`, again strictly increasing; at the lower k
it equals `127k-100>0`, sealed by `32258>30000`. Equation(3a) forces
the center and its three neighbors to have one actual opposite,
contradicting the extended-star obstruction. This proves (5) with
arbitrary spectators as well.

## Combined restricted-forest corollary

If n>=15, a farthest vertex cannot have all its finite components in

    {all irreducible types of rank<=4, D5,D6,E6,E7,E8}.  (6)

For this aggregation additionally pin
`compact_universal_small_a_component_energy_20261020.md`, SHA256
`801581c4b855289aca7982b2d0f4bb6f5849f151a78f40b86c9d403e84f415e6`,
and `compact_farthest_rank_four_forest_exclusion_20261020.md`, SHA256
`68ad62a619603eb7e6c2acc5689c901dde3f9d0897a7c1580b3deb948f5832b6`.
Their component estimates and short-wall classifications hold with
arbitrary finite spectators. The following simpler accounting avoids
extending their separate direction-propagation rates to E components.

Every non-A1 block in (6) pays strictly more than rank/14, except
possibly two named ISOLATED rank2 defects:

* OLD A2, NEW single7 contact: squared height h_A in(0,1/12);
* OLD I2(7), NEW single3 contact: squared height h_D in(0,1/3).

Each defect nevertheless pays more than9/100. These facts, including
the universal dihedral completeness, are proved in the pinned
rank-four source Sections2--4. Section5 of the pinned small-A source
gives at most ONE short A1 direction; all other A1 directions pay
strictly more than1/10>1/14. The only short A1 wall types are the
isolated dotted I and dotted-plus-A1_3 J, both with squared height
less than1/2.

An isolated rank2 defect cannot coexist with another such defect
block, of either type, or with an I/J short A1 wall. Their spatial
supports are orthogonal, and their actual wall product is minus the
product of their positive heights, strictly between -1/2 and0. This
violates the native Coxeter gap. In particular at most one defect
block occurs, and in its case there is no short A1 direction. All
other ranks then pay strictly more than rank/14, giving

    total energy >9/100+(n-2)/14 >=713/700>1.

If neither rank2 defect occurs, discard the energy of the possible
single short A1 direction entirely. The remaining rank is at least
n-1 and every remaining component/direction pays strictly more than
rank/14. Thus total energy is greater than `(n-1)/14>=1`. If no short
A1 direction exists, the still simpler bound is greater than n/14.
Every case contradicts the strict total farthest energy less than1,
proving (6).

The restricted family (6) still omits A/B components of rank at
least5 and D components of rank at least7. No ambient dimension
ceiling for unrestricted compact Coxeter polytopes is asserted.
