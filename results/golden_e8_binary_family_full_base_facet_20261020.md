# A full Bugaenko base for the entire unit-coupled golden binary family

20 October 2026. Paper proof only. No scientific execution, root search,
Lean, TeX or frozen-source change. This supplies a uniform construction
reduction, not uniform nonreflectivity or a new compact example.

## Statement

Let O=Z[phi], phi^2=phi+1. For every integer k>=-1 put

    u=phi^k,   a=phi*u,
    b=4+phi*u^2=4+phi^(2k+1),
    L_k=E8(O) orthogonal-sum [[-2phi,a],[a,2]].          (1)

The full integral reflection arrangement has a mirror

    r=(u,2;0),      r^2=2b,      Gr=(0,b;0).             (2)

Every other mirror meeting r is orthogonal to it and has a primitive
normal of square2 after unit rescaling. Its full intrinsic arrangement
is exactly the literal compact Bugaenko arrangement on
`E8(O)+<-2phi>`. The inherited Bugaenko reversal extends integrally,
fixes r with its orientation, and stabilizes the full adjacent chamber.

**No squarefreeness condition on b is required.** Primitive norm2b
mirrors are mutually ultraparallel. Proper divisor norm classes can
exist and are not discarded from the ambient chamber problem.

## 1. Admissibility and coprimality hold uniformly

Write n=2k+1>=-1. Then b=4+phi^n>4 and
`b'=4-phi^(-n)` lies in `[4-phi,4)`. The binary determinant is
`-phi*b`, negative physically and positive at the conjugate embedding;
the conjugate first pivot `-2phi'` is positive. Thus (1) is admissible
of signature(9,1), with conjugate positive definite. It has no nonzero
K-rational isotropic vector. Both b and b' increase with n, so
`N(b)>=(3+phi)(4-phi)=11`; in particular b is not a unit.

The ideal(b) is coprime to(30), for ALL such k:

* Modulo2 the ring is F4 and b=phi^n is a unit.
* Modulo3 it is F9 and phi^4=-1, so phi has order8. The equation
  b=0 would require the odd exponent n to be4 modulo8.
* At the unique prime over5 the residue is F5 with phi=3. Its odd
  powers are2 or3, and b=4+phi^n is respectively1 or2, not0.

Negative exponent n=-1 causes no exception since phi is a unit.
This proof neither assumes nor tries to decide squarefreeness of b.

## 2. Primitive norms and all finite intersection angles

The lattice is even. If m is a primitive integral normal of a lattice
reflection, write m^2=2d, with d in O totally positive. Primitivity
and reflection integrality give `Gm/d` integral. Applying the adjugate
and a Bezout combination of its primitive coordinates gives

    d | det(G),       hence d | b.                       (3)

The standard golden PID/unit facts are explicit inputs: every primitive
mirror has such an O-normal, and totally positive units are squares.
Consequently when(b/d) is the unit ideal one may rescale the primitive
normal by an O-unit and arrange d=b exactly. Equation(3) covers all
divisor exponents when b has repeated factors; it is not a squarefree
norm alphabet assertion.

For arbitrary INTERSECTING mirrors of this discrete arithmetic group,
the angle is a rational multiple of pi. Since its double-angle cosine
lies in Q(phi), the complete possible values x=4cos^2(theta) are

    0, 1, 2, 3, phi^2, 2-phi, 2+phi, 3-phi.             (4)

This includes the non-native arrangement angles2pi/5 and3pi/10;
one cannot use only the chamber's pi/m list. In the order displayed,
the rational norms of x are

    0, 1, 4, 9, 1, 1, 5, 5,

and those of x-4 are

    16, 9, 4, 1, 5, 5, 1, 1.                            (5)

Thus every nonzero value x and every value x-4 is a unit at each prime
dividing b. The classical finite-order/cyclotomic angle reduction in(4)
is the same explicit input used in the pinned prime-isolation proof.
Parallel distinct mirrors would give a nonzero K-null vector and are
excluded by conjugate positive definiteness.

## 3. Partial-divisor roots cannot meet a full-b mirror obliquely

Let v have primitive norm2b and m primitive norm2d with b/d a nonunit.
Reflection integrality gives `(v,L_k) subset bO`, so
`(v,m)=b*gamma`, gamma in O. At a nonorthogonal finite intersection,

    x=(v,m)^2/(b*d)=(b/d)*gamma^2.                      (6)

Choose any prime dividing b/d. The right side has positive valuation,
whereas every nonzero value in(4) is a unit there. This is impossible.
Thus such a mirror either is orthogonal to v or is ultraparallel.
Repeated factors in b/d do not change this argument.

## 4. Whole-modulus congruence removes the squarefactor boundary

Put c=(0,-1;0), c^2=2. The binary basis(r,c) has determinant -u,
an O-unit, and its Gram is

    [[2b,-b],[-b,2]].                                   (7)

For any primitive norm2b normal m, unimodularity of E8 and reflection
integrality give m_E8=b*z, z in E8(O). In the unit basis write

    m=lambda*r+mu*c+b*z,    lambda,mu in O.

Integrality of `(m,c)/b=-lambda+2mu/b` and gcd(b,2)=1 imply
mu=b*nu with nu in O. Since E8 is even, its norm equation becomes

    1=lambda^2-b*lambda*nu+b*nu^2+b*z^2/2,
    lambda^2=1 modulo b.                               (8)

This is an exact congruence modulo the ENTIRE principal ideal(b).
It does not assume that O/b is reduced. For two such primitive roots,

    (m1,m2)/b = 2lambda1*lambda2 modulo b.

Squaring and using(8) gives x=4cos^2(theta)=4 modulo b. Equation(5)
excludes every finite intersection, including orthogonality, since b
is a nonunit supported away from2,3,5. Coincident mirrors are omitted;
parallel mirrors were excluded above. Therefore distinct full-b mirrors
are mutually ultraparallel.

A rank-one-radical argument modulo a prime would require a valuation-one
factor of b. Formula(8), using the unit binary basis and the even E8
summand, is why that squarefactor boundary is absent for this family.
It is not a claim about arbitrary E8-plus-binary lattices with nonunit
coupling or bad-prime discriminant.

## 5. The complete base arrangement and inherited automorphism

Equations(3),(6),(8) show that every mirror crossing r is orthogonal.
For a primitive reflector m=(x,y;z) orthogonal to r, equation(2) gives
y=0. Reflection integrality against Y, whose product with m is a*x,
and the fact that a is an O-unit imply d|x. E8 unimodularity likewise
gives d|z. A nonunit d would divide every coordinate of m, contrary
to primitivity. Hence d is a unit and m can be normalized to square2.

The full orthogonal lattice is exactly `E8(O)+<-2phi>` on y=0.
Every one of its norm2 reflectors extends to L_k because its integral
covector has 2/q=1. The calibrated Bugaenko norm2 chamber is the full
intrinsic arrangement. Nonorthogonal ambient mirrors are disjoint
from r, so cannot insert extra cuts into it. Choose a generic interior
point of this intrinsic chamber and either adjacent ambient chamber C;
the entire compact Bugaenko chamber is the actual r facet of C.

Let S0 be its calibrated future diagram reversal. Its native time
vector X satisfies `S0X-X in 2(E8+OX)`. The extension

    S|_(E8+OX)=S0,
    SY=Y-(u/2)(S0X-X)                                  (9)

is integral, isometric, involutive and fixes r positively. It preserves
the intrinsic full facet chamber and its chosen inward side. The generic
facet-germ theorem therefore gives S(C)=C. This supplies a genuine
full-chamber involution uniformly, not merely an integral lattice map.

## Scope and next construction obligation

All k>=-1 now have the FULL Bugaenko base and the inherited full-chamber
involution without separate squarefree tests or base-root searches.
This does not supply a second independent chamber automorphism. Nor
does it prove that c is the first transverse wall when proper divisor
norm classes occur: such classes remain available away from r and may
change edge continuations. The known unit-bisector mechanism is special,
and its completed k1 golden31 exclusion does not extend to all k here.

The infinite family is therefore not excluded. Any further route must
control transverse actual mirrors or derive a new global obstruction,
rather than repeat the already settled full-base calibration. No other
literal lattice, similarity class, commensurability class or dimension
ceiling is excluded by this note.

## Bound mathematical inputs

* `results/bugaenko_facet_binary_extension_and_collar_audit_20261020.md`,
  SHA256 `774e87f7abb33680a06879a17a5be9f06e2218bf12a95023bcff5c7046ba1227`:
  golden ring facts, primitive normalization, family and S0 glue.
* `results/compact_h9_golden_discriminant59_edge_route_20261020.md`,
  SHA256 `725d8775cebec51d8c18b4b0a0d3be0255e8b1472621adb4d3cefa656adb429e`:
  complete finite intersection-angle list and intrinsic calibration.
  Its prime59 conclusion is not used as a family exclusion.
* `data/known_examples/bugaenko_compact_h8.json`, SHA256
  `c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c`,
  and its complete transcript, SHA256
  `c4178fc4678bacb48c430b2cec18e16577078824fdd6dea39375e63ed6ae015a`.

The new whole-modulus proof and its scope await independent review.
No executable certificate or worldwide novelty claim is made.
