# The quadratic13 host: the far A9 edge and an infinite chamber period

20 October 2026. **Complete paper derivation; no new worker has been run.**
This closes the literal lattice `L=<-eta>+I9` by a complete eight-profile
finite interval, followed by an integral infinite-order chamber isometry.
It does not exclude its K-similarity or commensurability class, or compact
hyperbolic9 Coxeter polytopes generally.

## Inputs and conventions

Put `omega=(1+sqrt(13))/2`, `eta=omega+1`. Thus
`eta²=3eta+1`, `eta'=-1/eta`, and `O=Z[omega]=Z[eta]`.
Coordinates are `(time; nine spatial coordinates)` and the form is
`-eta*time² + sum spatial²`. The chamber inequalities are `(x,r)<=0`.

The complete primitive root alphabet `{1,2}`, the full initial B9 germ,
two first finite edges, and the integral chamber involution G are proved
in `compact_h9_sqrt13_b9_two_edge_route_20261020.md`, SHA256
`b922647902d32beccd98fbe5a2781ddcbb9cafaeb10a12902eb3ebbce809547f`.
The subsequent full A9 corner is proved in
`compact_h9_sqrt13_a8_edge_interval_20261020.md`, SHA256
`e3e62a6ec4d834e2e0c7d4abc4d28961839f1a467ac06966ef08854757814228`.
These are written geometric/arithmetic inputs, not conclusions of a
future scalar checker.

In particular, for `alpha_i=-e_i+e_(i+1)` (1<=i<=8) and `alpha_9=-e_9`,
the chamber C has the actual simple wall

    a=(2eta-1;eta+1,eta,eta,eta,eta,eta,eta,eta,eta)

and the full actual A9 vertex

    Q=(9eta+1;(5eta+2)^9),       Q²=-10eta.

Its outward simple normals are `a,alpha_1,...,alpha_8`. Release alpha_8,
retaining the A8 chain

    beta_1=a, beta_2=alpha_1, ..., beta_8=alpha_7.

The word "first" below always means first in this full integral mirror
arrangement, starting in the actual chamber germ at Q. It does not mean
first among a finite cap list.

## 1. An integral diagonal basis for the retained orthogonal plane

A vector perpendicular to all beta_i has the form

    (x;z,z,z,z,z,z,z,z, (2eta-1)x-(eta+5)z).

In its `(x,z)` coordinates the following are an O-basis:

    P0=(5eta+1;(eta²)^8,2eta+1),       P0²=-eta,
    R0=(7eta+1;(4eta+2)^8,3eta-1),     R0²=9,
    (P0,R0)=0.

Indeed their two-column determinant is eta, a unit. Write `P=P0/eta`;
explicitly

    P=(eta+2;eta^8,omega),      P²=-1/eta.

Both endpoints used below lie on this plane:

    Q  =(eta+3)P0-eta R0,
    Q' =(eta+3)P0+eta R0.

They are future timelike, with square -10eta. Their midpoint ray is P0.
The second expression defines Q' initially; no chamber membership of P
or Q' is assumed in using this plane.

## 2. A proved finite upper, before any first-wall assumption

Let w0 be the longest Weyl element of the retained A8. It fixes its
orthogonal plane pointwise and sends `beta_j` to `-beta_(9-j)`.
It is integral, being a product of reflections in norm-two roots.
The point inversion

    J_P(x)=-x-2eta*(P,x)P

is integral because P is integral and its squared norm is the unit
`-1/eta`. It is an involutive Lorentz isometry preserving the future
cone, fixing P and negating P-perp. It commutes with w0.

Set `H=J_P w0`, another integral future-preserving involution. For the
released root,

    b0=w0 alpha_8=a+alpha_1+...+alpha_8
      =(2eta-1;eta^8,eta+1),       (P,b0)=-1,
    b=H alpha_8=-b0+2eta P
      =(8eta+3;(5eta+2)^8,eta²).

Thus b is a primitive integral norm-two root. Direct products are

    (b,a)=-1,       (b,alpha_i)=0  (1<=i<=7),
    (b,alpha_8)=-2eta-1,          (b,alpha_9)=-eta².

The earlier walls `v=(1;omega,0^8)` and
`d=(eta;eta+1,eta+1,1,0^6)` have products
`(b,v)=-15eta-5`, `(b,d)=-40eta-11`. Hence these are native finite
bonds or strict ultraparallel pairs, not a numerical compatibility test.
Moreover `H Q=Q'`, `(Q',b)=0`, while `(Q,b)<0`. The displayed Q' is
an ordinary finite point forward on the released edge. Consequently b
gives a finite upper even before it is known to be a wall of C.

## 3. Complete profiles and their exact lattice glue

A first new wall at an ordinary point on this edge has a positive
definite one-root extension of the retained A8. The finite crystallographic
classification and the complete norm alphabet give exactly:

|Full corner diagram|New norm q|Contact node k|R|
|---|---:|---|---:|
|A8+A1|1|none|9|
|A8+A1|2|none|18|
|B9|1|1 or8|1|
|D9|2|2 or7|4|
|A9|2|1 or8|10|

The six attached placements have product -1 with beta_k and zero with
the other retained roots. There is no exceptional finite rank-nine
diagram containing A8. Multiple contacts or other placements are
nonfinite. No spectator-contact restriction is being imposed.

For an attached profile its integral coordinates are necessarily

    r=(x; (z+1)^(k-1), z^(9-k),
              (2eta-1)x-(eta+5)z-(k+eta-4)),

where the first two spatial blocks together have length eight. Put

    A=-(P0,r)/eta=-(P,r)
      =2omega*x-(eta+4)z-(k+eta-5),
    Z=(R0,r)=-9eta*x+9(eta+2)z+(eta+3)k-9.

The orthogonal projection onto the retained plane is
`A P0+(Z/9)R0`; the white projection is the negative fundamental weight
at k and has square `k(9-k)/9`. Thus exactly

    Z²-9eta*A²=9q-k(9-k)=R,
    A in O,              Z congruent (eta+3)k modulo9O.       (1)

The congruence is also sufficient for lattice reconstruction: the
orthogonal vector with the same `(x,z)` has coefficients
`A+(k+eta-5)` on P0 and `(Z-(eta+3)k)/9+1` on R0. Both are integral
precisely as asserted. This proves the glue condition, rather than
assuming an abstract A8 reversal extends to L.

For isolated profiles the root itself is in the orthogonal plane,
so it is `A P0+B R0` with A,B in O. Its norm equation is

    9B²-eta*A²=q.

For q=1 reduction at `(3,omega)` gives `A²=2` over F3. For q=2,
reduction at `(3,omega-1)` gives the same contradiction. Thus both
isolated profiles are excluded at all heights. Reducing only (1) and
forgetting Z=9B would not prove this conclusion.

## 4. The entire finite interval has one new candidate

Orient any prospective outgoing wall by `(Q,r)<0`. Contact by Q' means
`(Q',r)>=0`. Since both points are perpendicular to the white roots,

    (Q,r) =-eta*((eta+3)A+Z),
    (Q',r)=-eta*((eta+3)A-Z).

These imply `Z>0` and `Z>=(eta+3)|A|`, where equality on the starting
side is excluded by the strict first inequality. As
`(eta+3)²-9eta=10`, equation (1) gives the exact interval bounds

    10A²<=R,              A'²<=R*eta/9.                    (2)

The second follows from conjugate positivity. For every attached
profile R<=10, so `|A|<=1`, `|A'|<2`. If A=c+d*omega, then
`|d|=|A-A'|/sqrt(13)<3/sqrt(13)<1`. Hence A is an integer.
Consequently A=0 for R=1,4, and A is one of -1,0,1 for R=10.
No height prefix or finite numerical search remains.

* R=1: Z=1, but neither k=1 nor k=8 satisfies the congruence in (1).
* R=4: Z=2, but neither k=2 nor k=7 satisfies that congruence.
* R=10,A=0: this requires sqrt(10) in Q(sqrt(13)), impossible.
* R=10,A=+/-1: the positive square root is Z=eta+3. The congruence
  permits k=1 only. For k=8 its eta coefficient would require
  `1=8 modulo9`, which is false.

For clarity, the preceding R1/R4 failures already follow from their
eta coefficients: k is nonzero modulo9, whereas Z=1 or2 has eta
coefficient zero. At R10,A=-1 the reconstruction is exactly `-b0`;
it has `(Q,-b0)=0` and is excluded by the strict outgoing condition.
It is a known non-simple root through the starting vertex, not an
earlier edge wall. At A=+1 the reconstructed root is exactly b.

Thus b is the unique new primitive candidate in the whole interval.
Local finiteness of the arithmetic mirror arrangement and the finite
upper imply that it is the actual first simple wall. The new vertex
Q' has the full actual A9 inward germ `beta_8,...,beta_1,b`.
This is a full-arrangement statement, with no old-cap filter and no
incomplete simultaneous-candidate selection.

## 5. An integral chamber period of infinite order

The integral isometry H permutes the entire root arrangement. It sends
Q to Q', sends beta_j to beta_(9-j), and sends the released alpha_8 to
the now-proved first wall b. Hence it carries the full inward germ of C
at Q to the full inward germ of C at Q'. Unique continuation of an
arrangement cell proves `H(C)=C`. The earlier two-edge proof similarly
gives `G(C)=C` and `G²=I`. No assertion that arbitrary partial corner
maps preserve a chamber is used here.

We compute the trace of HG without a numerical matrix program. The four
mutually orthogonal norm-two roots

    gamma_1=a+alpha_1+...+alpha_7,
    gamma_2=alpha_1+...+alpha_6,
    gamma_3=alpha_2+...+alpha_5,
    gamma_4=alpha_3+alpha_4

give `w0=I-sum gamma_j*(gamma_j,.)`. The explicit G columns in the
first input give:

|j| (gamma_j,G gamma_j) |
|---|---|
|1| -63eta-19 |
|2| 1-2eta |
|3| -eta |
|4| 2 |

Since `trace(G)=6`, this gives `trace(w0 G)=66eta+22`. Independently,

    G P=(8eta+1;11eta+3,8eta+3,eta^6,omega),
    (P,G P)=-51eta-14.

Because w0 fixes P and H=J_P w0,

    trace(HG)=-trace(w0G)-2eta*(P,G P)
             =268eta+80 > 10.                              (3)

The ambient vector space has dimension ten. A finite-order real matrix
has complex eigenvalues of modulus one and hence trace of absolute
value at most ten. Thus HG has infinite order. This conclusion does
not follow merely from composing two different involutions; (3) is the
needed certificate.

## 6. Literal-host nonreflectivity and exact scope

Let W be the group generated by every integral root reflection of L.
Its arrangement is locally finite and C is a fundamental chamber;
W acts simply transitively on its chambers. Every power of HG preserves
C, and no nonidentity element of W preserves its open chamber. Therefore
these powers give infinitely many distinct cosets of W in O^+(L).
The latter is an arithmetic lattice, so infinite index implies W is not
a lattice: its chamber C has infinite hyperbolic volume. Equivalently,
the literal integral lattice L is nonreflective.

In particular there is no compact hyperbolic9 Coxeter polytope all of
whose reflecting walls preserve this literal L: its reflection group
would lie in W, and a discrete overgroup of a lattice has finite index
and finite covolume. This is not a classification of K-forms, of
commensurable lattices, or of all compact9 reflection groups.

All new reductions in this note are exact paper algebra. A future tiny
checker may verify (1)--(3), primitive norms, the column Gram identities
and the finite list of congruences. Such a checker does not replace the
written primitive-root theorem, actual first-contact argument, chamber
germ continuation or arithmetic-lattice/index reasoning above.
