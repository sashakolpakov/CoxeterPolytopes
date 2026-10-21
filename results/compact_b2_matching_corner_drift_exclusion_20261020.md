# Repeated B2 corners cannot be drift witnesses in dimension at least twelve

20 October 2026. Direct geometric proof by root; independent review
pending. No enumeration, computation, manuscript build or Lean job.

Let P be a compact hyperbolic Coxeter polytope with a finite vertex p
whose diagram is a disjoint union of r copies of I2(4)=B2. Its dimension
is n=2r. For ANY globally shared positive white-diagram weight, p cannot
satisfy the strict drift inequality A_f(p)<||rho_f(p)|| when r>=6.
This excludes one more vertex family; it is not a dimension ceiling.

Write s=sqrt2, alpha=1+s, and c=1/sqrt2. The unit-normal Gram of one
block and its inverse are

    B=[[1,-c],[-c,1]],       B^-1=2[[1,c],[c,1]].

## 1. Actual closing walls and their supports

An outward unit wall m not through p has coordinates

    m=a*p-sum_i beta_i*omega_i,   a=-<p,m>>0,
    beta_i=-<m,n_i>>=0,
    a^2=beta^t*B_total^-1*beta-1.

Every nonzero beta is cos(pi/k), k>=3, or exceeds one. Equality one
is excluded by compactness: the positive null sum of that pair of unit
normals would give an ideal point feasible against all other walls.
The native positive values therefore start at1/2,1/sqrt2,phi/2,...,
where phi=(1+sqrt5)/2.

Suppose m closes an edge in the distinguished block. Its common-white
Gram, consisting of m, the retained mate and all other B2 blocks, is
positive definite. A contact to one complete other B2 must be a single
simple3 contact. Indeed two contacts have Schur cost at least1+c>1;
a single contact has cost2*beta^2<1, forcing beta=1/2. Each such
off-block contact therefore has cost1/2. At most one exists, since two
would already give white Schur cost at least one. This conclusion uses
positive definiteness alone, without a diagram classification.

Call the two distinguished contacts beta,gamma and let E be the
off-block cost, either zero or1/2. The exact closing-edge formula is

    cosh L_beta=sqrt2*(beta+c*gamma)/sqrt(1-E-gamma^2).
                                                               (1)

The denominator is positive at the actual finite endpoint.

If beta and gamma are both nonzero, (1) gives

    cosh L_beta >= alpha/sqrt3.                         (2)

If E=1/2 and gamma=0, positive height forces beta>1/2, hence
beta>=1/sqrt2 and cosh L_beta=2beta>=sqrt2>alpha/sqrt3.
If E=1/2 and beta=0, the finite endpoint condition forces gamma=1/2;
the total spatial cost is then one, giving a=0, impossible. Thus every
edge whose closing wall has an off-block contact also satisfies (2).

The only remaining case is E=0 with a single distinguished contact.
If it is dotted, every edge it can close has cosh L>sqrt2, so again (2).
If it is finite, positive height requires beta>1/sqrt2, hence its
label k is at least five. This wall is supported ONLY in that one block.

## 2. A finite one-contact wall closes both block directions

Here is the actual-face fact needed to avoid assigning incompatible
independent lower bounds to the two edges. In an acute polytope every
positive-definite subset S of actual unit wall normals has nonempty
intersection with P. For an interior future vector q, put

    z=q-N_S*B_S^-1*(N_S^t*q).

B_S has nonpositive off-diagonal entries and positive definite inverse
with nonnegative entries. Thus z=q+sum_s lambda_s*n_s with lambda_s>0.
Every other wall has negative product with z by acuteness. Moreover
z^2=q^2-(N_S^t*q)^t*B_S^-1*(N_S^t*q)<0, and
<q,z>=z^2<0, so z is future timelike. Normalizing z gives a point of P
on all walls of S.

For a one-contact finite wall and either choice of retained mate, the
n-wall Gram is positive definite: its only new two-node factor has
angle pi/k or pi/2. The preceding fact makes the corresponding endpoint
actual. Consequently the SAME wall closes BOTH old block directions;
another wall cannot truncate either endpoint first.

Their exact lengths satisfy

    cosh L_1=sqrt2*beta,
    cosh L_2=beta/sqrt(1-beta^2),       beta=cos(pi/k).    (3)

Both increase with beta. For k>=6, (3) implies

    tanh(L_1/2)+tanh(L_2/2)
       >=sqrt3-sqrt2+sqrt(2-sqrt3)>4/5.                 (4)

For the strict rational estimate, sqrt3-sqrt2>3/10 follows from
sqrt3+sqrt2<13/4<10/3; and sqrt(2-sqrt3)>1/2 follows from sqrt3<7/4.

For k=5, (3) instead gives

    cosh L_1=phi/sqrt2,   cosh L_2=phi/sqrt(3-phi),
    tanh(L_1/2)>1/4,     tanh(L_2/2)>3/8.               (5)

Squaring the equivalent positive cosh inequalities reduces them to
phi>353/225 and phi>6481/4177, respectively; both follow from phi>8/5.
Thus this exceptional block still has half-angle sum greater than5/8.

There is at most ONE exceptional block at p. The height of a k=5
one-contact wall satisfies a_5^2=(sqrt5-1)/4, strictly between zero and
one half. Two such walls supported in different blocks would have
product -a_5^2 in(-1/2,0), forbidden by the Coxeter gap. This counts
actual walls, with their two already established edge directions.

## 3. Summing the blocks

Put t_F=tanh(arcosh(alpha/sqrt3)/2). Then t_F>2/5, since
alpha/sqrt3>29/21 follows from sqrt2>200/147. A block with no finite
one-contact closing wall has BOTH edges bounded by (2), so its
half-angle sum exceeds4/5. A block with such a wall has the sum (4),
except for the at-most-one k=5 block covered by (5). Therefore

    sum_edges_at_p tanh(L_i/2) > (r-1)*4/5+5/8.

Every old inverse-Gram diagonal is two. Every released white diagram
is the same A1+B2^(r-1), so any shared positive white function has the
same value f at all directions. It follows that

    A_f(p)>f*sqrt2*(4r/5-7/40),
    ||rho_f(p)||^2=f^2*r*(4+2sqrt2).                    (6)

At r=6, the squared lower bound in (6) divided by f^2 is1369/32,
which exceeds42>24+12sqrt2. The ratio (4r/5-7/40)^2/r increases for
r>=1: after multiplication by1600 its consecutive difference is
1024-49/(r*(r+1))>0. Thus (6) gives A_f>||rho_f|| for every r>=6.

The global identity used here is in
`results/compact_finite_edge_white_weighted_drift_20261020.md`.
Its strict witness can therefore never have diagram B2^r in dimension
at least twelve, for ANY common-white positive function f. No arithmetic
assumption is used. Mixed labels, A1/A2 factors and higher-rank finite
factors still require separate compatible estimates. No new uniform
dimension bound follows from this family exclusion alone.
