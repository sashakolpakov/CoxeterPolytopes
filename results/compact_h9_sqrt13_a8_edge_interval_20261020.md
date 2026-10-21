# The quadratic13 host: a complete A8 edge to an actual A9 vertex

20 October 2026. Own paper proof, before execution. This supplies a third
complete finite edge of the literal H9 host. It proves neither finite
volume nor reflectivity of the full chamber.

Use O=Z[omega], omega²=omega+3, eta=omega+1, eta²=3eta+1, and the literal
L=<-eta>+I9. The complete primitive norm pool {1,2}, the actual B9 germ
at p=(1;0^9), the first two finite edges and the chamber involution G
are proved in `compact_h9_sqrt13_b9_two_edge_route_20261020.md`, SHA256
`b922647902d32beccd98fbe5a2781ddcbb9cafaeb10a12902eb3ebbce809547f`.
Native outward roots are alpha_i=-e_i+e_(i+1) for i=1,...,8,
and alpha_9=-e_9. No arbitrary-height root prefix is used below.

## Complete finite inventory and a small upper

Release alpha_9 at p, retaining A8. The edge is
p+tau*(e_1+...+e_9), with light endpoint tau=sqrt(eta/9).
A new simple normal is either A8-orthogonal, or has one product -1
with a retained root. Finite classification gives precisely:

|New full diagram|Root norm q|Contact node k|R|
|---|---:|---|---:|
|A8+A1|1|none|9|
|A8+A1|2|none|18|
|B9|1|1 or8|1|
|D9|2|2 or7|4|
|A9|2|1 or8|10|

Thus there are eight literal placements. B9 needs the short new end
root, whereas A9 and D9 need a long root and simple bonds. Other
contacts give a cycle, a nonfinite fork, or a forbidden finite bond for
the two norm classes. The finite rank-nine list has no exceptional
component containing A8.

For a contact at k, the spatial entries are z+1 in the first k positions
and z in the others. Orthogonal profiles have all entries z and k=0.
In particular x,z belong to O, where x is the time coordinate. Put

    Z=9z+k,       Z²-9eta*x²=R=9q-k(9-k).

The forward orientation has x>0,Z>0; its contact is tau=eta*x/Z.
Spatial integrality is exactly Z≡k modulo9O. It is a necessary condition
before any old-wall compatibility test.

The explicit norm-two mirror

    a=(2eta-1;eta+1,eta,eta,eta,eta,eta,eta,eta,eta)

is integral and reflective. Its variables are x*=2eta-1 and Z*=9eta+1,
with Z*²-9eta*x*²=10. Its contact is strictly before the light endpoint.
It is therefore a proved finite upper, without assuming it is first.
For any profile contacting by this upper, comparison of the positive
contact ratios gives the exact bounds

    10x²<=R(2eta-1)²,       x'²<=R*eta/9.                 (1)

The second follows from eta'=-1/eta and conjugate positivity. The root
alphabet and local finiteness ensure that a first simple wall occurs by
this ordinary finite point and has one of the displayed profiles.

## All-height isolated exclusions

The constant-spatial equations are 9z²-eta*x²=q. At the split prime
(3,omega), eta=1; q=1 would force x²=2 modulo3. At the other split
prime (3,omega-1), eta=2; q=2 would again force x²=2 modulo3. Both are
impossible. These use different places above3, and exclude both
isolated profiles at every height. They do not follow by reducing the
Z equation alone and forgetting that Z=9z.

## Tiny complete domains

Write x=A+Bomega. The inequalities in(1), together with
3<eta<10/3 and 18/5<sqrt13, give these exhaustive domains.

* For R=1, x<9/5 and |x'|<5/8. Thus
  B=(x-x')/sqrt13 is an integer strictly between-1 and1. Then x is
  a positive integer with |x'|<1, impossible.
* For R=4, x<18/5 and |x'|<5/4, so B is0 or1. At B=0 only x=1
  survives. At B=1 the physical and conjugate bounds give A=1,
  hence x=eta. Both meet the exact bounds; no other integer pair does.
* For R=10, x<=2eta-1<6 and |x'|<2, hence B is0,1 or2.
  At B=0 only x=1 survives. At B=1, A=0,1,2,3. At B=2 the physical
  upper gives A<=1 and the conjugate bound gives A>=1. These are the
  six entries in the second table below.

For reference, the coarse R=4 bounds follow from
(2/5)(5+8eta)<38/3<(18/5)² and 4eta/9<40/27<(5/4)².
The R=1 bounds follow from (5+8eta)/10<19/6<(9/5)² and
eta/9<10/27<(5/8)². The indicated coefficient strips include every
possible root, not merely small-unit or small-height samples.

For E=A+Beta, its rational norm is A²+3AB-B². A field square must have
square rational norm. The D9 tests are:

|x|N(9eta*x²+4)|Outcome|
|---|---:|---|
|1|43|Nonsquare|
|eta|1231|Nonsquare|

The A9 tests are:

|x|N(9eta*x²+10)|Positive square root, if any|
|---|---:|---|
|1|289|(eta+3), which fails Z≡1 and Z≡8 modulo9|
|omega|901|None|
|eta|3259|None|
|eta+1|4861|None|
|eta+2|1819|None|
|2eta-1|2809|9eta+1, permitting k=1 only|

The first radicand equals (eta+3)². Its omega coefficient is1, so its
positive square root is congruent to neither integer1 nor8 modulo9O.
The last radicand is82+261eta=(9eta+1)². All other printed norms are
nonsquares. A nonzero square has only the two signed square roots in
a field; positive contact selects the displayed positive ones.

Consequently the entire interval has exactly one primitive candidate:
the proposed long root a, at the k=1 placement. It is the first SIMPLE
wall of the full chamber. No finite cap list or old-wall product was
used to obtain this uniqueness.

## Actual endpoint and scope

Postselection products with the native walls are
(a,alpha_1)=-1, (a,alpha_i)=0 for2<=i<=8, and
(a,alpha_9)=-eta, the latter strictly ultraparallel. With the two earlier
actual caps v,d, the products are -2eta-2 and -7eta-2 respectively,
also strictly ultraparallel. Thus the new wall is acute/native to all
eleven old walls.

Its future apex is

    Q=(9eta+1;(5eta+2)^9),       Q²=-10eta.

It vanishes on a and alpha_1,...,alpha_8; these nine normals have the
positive A9 Gram. The other old inequalities are strict:

    (Q,alpha_9)=-5eta-2,
    (Q,v)=-16eta-6,       (Q,d)=-44eta-12.

The complete first-contact argument gives the full actual finite A9
inward germ. The chamber involution G already proved for L supplies
another actual A9 corner as its image. Neither a finite-volume claim
nor an infinite-order chamber symmetry follows from these two corners:
G remains the known involution.

The primitive alphabet, finite extension classification, local
finiteness and full-chamber interpretation are written inputs to a
future exact arithmetic checker. The literal host remains open. No
additional manuscript, Lean result or compact-H9 example is claimed.
