# The golden discriminant31 member: a unit-bisector obstruction

20 October 2026. Complete paper proof and fixed-check preparation, awaiting
outside mathematical review. No scientific execution, prefix search, import,
Lean or TeX job was performed. This is the GOLDEN norm31 host, not the older
closed norm31 lattice over Q(sqrt2). No global status is changed here.

## 1. Literal host, norm pool, and full base facet

Put O=Z[phi], phi^2=phi+1, and use the native E8 order of the calibrated
Bugaenko fixture. Consider

    L31=E8(O) orthogonal-sum B,
    B=[[-2phi,phi^2],[phi^2,2]],
    b=5+2phi,    delta=phi*b=2+7phi,    det B=-delta.      (1)

Coordinates are (x,y;z). At the physical embedding B is Lorentzian;
at the conjugate embedding -2phi'>0 and -delta'>0. Thus L31 is
admissible of signature(9,1), with no nonzero K-rational null vector.
Since N(b)=31, (b) is prime and the determinant has valuation one.

The golden PID/unit and primitive-reflector argument in the pinned
binary audit applies: a primitive reflector of square2d has d dividing
det(L31), and totally positive units are squares. After unit rescaling,
the complete necessary norm pool is

    q=2 or 2b.                                          (2)

Both occur: c=(0,-1;0) has square2, and r=(phi,2;0) has square2b and
covector Gr=(0,b;0). A long reflector has all pairings divisible by b.

The whole-arrangement prime isolation proof of the pinned L59 route
applies verbatim to prime31: every long mirror meets every short mirror
orthogonally, and distinct long mirrors are ultraparallel. Its finite
angle norms involve only primes2,3,5, including the non-native5/10
arrangement angles; neither31 nor a parallel K-null pair is possible.
Consequently the FULL r facet is y=0, with exactly the calibrated
E8(O)+<-2phi> reflection arrangement. Conversely every intrinsic
norm2 reflection extends, since 2/q=1. In particular all eleven
Bugaenko facet walls are actual, rather than a finite seed.

The full zero chamber at p0=(1,0;0) is E8+A1(r), by the same prime
isolation/norm argument as for L59. Orient the chamber by r<=0, so
y<=0. The inherited Bugaenko reversal S extends integrally fixing r:
on the time basis vector X and second binary vector Y,

    SX=(9+16phi)X+(6+10phi)omega1,
    SY=Y-(phi/2)(SX-X).                                (3)

The half is integral since SX-X is twice an intrinsic lattice vector.
This isometry preserves the inward generic germ of the FULL r facet,
so S stabilizes the full ambient chamber C. It is a future involution
with trace2. The already closed families19/59 are not being reopened.

## 2. The first transverse edge is complete without a search

The binary basis(r,c) has determinant -phi, a unit, and Gram
`[[2b,-b],[-b,2]]`. Its two mirrors bound the finite segment from p0
to p1=(2,-phi^2;0), where p1^2=-2delta and p1.c=0.

A first intervening mirror along this retained E8 edge must be binary:
a finite spherical diagram cannot properly extend E8. Write its
primitive normal lambda*r+mu*c, lambda,mu in O. Intervening contact
forces lambda*mu<0 at the physical embedding.

For q2, the physical norm equation gives lambda^2<=1/b. Conjugate
positive definiteness gives lambda'^2<=4/[b'(4-b')]. Hence

    N(lambda)^2 <= 4/[31(4-b')]<1,
    4-b'=2phi-3>1/5.                                  (4)

This excludes nonzero lambda. For q2b, reflection integrality against c
gives b|2mu, hence mu=b*gamma. The identical two-embedding argument
applies to gamma and excludes it. Thus no mirror intervenes, c is
actual, and p1 has full E8+A1(c) germ. This proof does not filter out
incompatible roots by an arbitrary finite seed.

## 3. An actual E7 germ at a unit timelike bisector

In the native fixture let

    omega7=(2,3,4,6,5,4,3,2),    omega7^2=2,
    g=(1,0;phi*omega7),         g^2=2.

This g is actual Bugaenko wall8. It has c.g=-phi^2. Put

    t=c+g=(1,-1;phi*omega7),     t^2=-2/phi.              (5)

The vector is future, since p0.t=-(1+3phi)<0. It is feasible for ALL
actual walls: for any wall other than c,g this follows by summing
its two nonpositive products; for c,g the product is -1/phi<0.
Every actual wall through t must therefore be orthogonal to BOTH c,g.
The seven old walls -e0,...,-e6 do so and form E7.

Here is a complete exclusion of any further actual zero wall. For a
primitive reflector m orthogonal to c,g, integrality and c-orthogonality
give

    m=(2u,-phi^2*u;z),   u in O,
    z.omega7=b*u,
    m^2=z^2-2phi*b*u^2.                                (6)

For q2, positive Cauchy in E8 at BOTH embeddings gives

    b*(b-4phi)*u^2<=4,
    N(u)^2<=16/[N(b)*N(b-4phi)]=16/341<1.

Thus u=0. Then z.omega7=0 and z^2=2, so z is an ordinary E7 root.
For clarity, pure E8(O) norm2 vectors are ordinary E8 roots: writing
z=A+phi B gives A^2+B^2=2 and 2(A,B)+B^2=0; evenness and minimum
norm2 force B=0. This also identifies the zero-root system, not just
its rank.

For q2b, E8 unimodularity and reflection integrality give z=b*z0.
Equations(6) become

    b*z0^2-2phi*u^2=2,      z0.omega7=u.

Cauchy at both embeddings now yields

    u^2<=4/(5-2phi)<phi^2,
    u'^2<=4/(5-2phi')<1,
    N(u)^2<=16/11<2.                                  (7)

A nonzero u must be a unit, hence +/-phi^k. Its conjugate bound forces
k>0, while its physical bound excludes k>=1. If u=0, the equation
b*z0^2=2 is impossible since z0^2 is even and b is a nonunit. Hence
there is no long root orthogonal to both c,g.

We have proved that the actual germ at the normalized t is EXACTLY
E7. The argument does not assume that a merely integral involution
preserves a chamber, or that a short metric facet is itself Coxeter.

## 4. A second full-chamber involution and its infinite product

Let P8 be the orthogonal projector to the E8 summand. The integral
involution

    W=I-2P8+omega7*(omega7,. )

fixes the binary block and omega7, and acts as -I on its orthogonal
E7. It equals the usual longest E7 action; the displayed formula
already proves its integral extension and isometry. Define

    H=(-I-phi*t*(t,. )) W
     =-I+2P8-omega7*(omega7,. )-phi*t*(t,. ).             (8)

The first factor is the future point inversion at t, since t^2=-2/phi.
Thus H is an integral future isometry, H^2=I, Ht=t, fixing each of the
seven actual E7 wall normals individually. It therefore preserves
the chamber's inward local germ at t. Local finiteness and uniqueness
of the chamber with this germ imply H(C)=C. No firstness of a new
unknown wall is required here.

The trace calculation uses the same native fixture, avoiding a basis
permutation. Write k10 for its last Bugaenko wall. Then

    Se0=e7, Se1=-k10, Se2=e6, Se3=e5, Se4=e4,
    Se5=e3, Se6=e2, Se7=e0,
    omega1=(5,8,10,15,12,9,6,3).

Its relevant exact contractions are

    tr(S)=2,                  tr(S P8)=-8-16phi,
    (omega7,S omega7)=-10-18phi,
    (t,St)=-20-36phi.                                  (9)

For the last identity, equations(3) and Sg=wall0 give

    St=(8+13phi,-1;(5+8phi)omega1+phi*omega0),
    omega0=(4,5,7,10,8,6,4,2),
    (omega7,omega1)=3, (omega7,omega0)=2.

Its binary contribution is -46-77phi and its E8 contribution26+41phi.
Also tr(S P8)=1-(9+16phi) from the two nonzero native diagonal
contributions, and the omega7 contraction is2-3(4+6phi). Formula(8)
therefore gives

    tr(SH)=-tr(S)+2tr(S P8)-(omega7,S omega7)-phi*(t,St)
          =28+42phi>10.                                (10)

A finite-order real10-dimensional isometry has trace of absolute value
at most10. Thus SH has infinite order and stabilizes C. The full
reflection group acts freely on chambers, so its intersection with
Stab(C) is trivial. Powers of SH give infinitely many distinct cosets
modulo that full reflection group. Hence L31 is nonreflective.
The explicitly imported Borel--Harish-Chandra arithmetic lattice
theorem, together with conjugate-positive anisotropic compactness,
then gives infinite chamber volume. This excludes this literal host,
not its whole similarity or commensurability class.

## 5. Infinite-family scope and fixed-check plan

The infinite audited family is B(u)=[[-2phi,phi*u],[phi*u,2]], u=phi^k,
k>=-1. The present argument closes its k1 member. It does NOT give
an obstruction for every k. For an inherited Bugaenko wall of time
phi^j, j=0,1,4, the analogous bisector square is
2(2-phi^(k+j+1)). Its integral point inversion requires
phi^(k+j+1)-2 to be a unit: the vector has a unit y-coordinate and
its E8 covector has a unit entry. In the timelike range the Lucas-norm
identity in the pinned audit makes this possible only at exponent2.
This recovers exactly k0,j1 (closed H19) and k1,j0 (here), not all
larger determinants. Other mirrors and other architectures remain open.

If independently approved, ONE tiny fixed checker may reconstruct
the native E8 Gram and fixture, B/r/c/g/t/S/W/H, all integral adjoint
inverses, first-edge scalar identities, both germ Cauchy norm seals,
the E7 Gram/orthogonal complement, future signs, and(9)--(10). There
is no height loop, Pell search, or cone enumeration. The geometric
actualness, prime isolation and group-index implications remain written
premises. Root alone may schedule bounded generation and fresh replay.

## Fixed input bindings

* `results/bugaenko_facet_binary_extension_and_collar_audit_20261020.md`,
  SHA256 `774e87f7abb33680a06879a17a5be9f06e2218bf12a95023bcff5c7046ba1227`:
  PID/unit facts, primitive pool, infinite family, inherited extension.
* `results/compact_h9_golden_discriminant59_edge_route_20261020.md`,
  SHA256 `725d8775cebec51d8c18b4b0a0d3be0255e8b1472621adb4d3cefa656adb429e`:
  prime isolation proof including arbitrary arrangement angles and full
  base-facet calibration; its prime59 outcome is not reused as a31 proof.
* `data/known_examples/bugaenko_compact_h8.json`, SHA256
  `c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c`;
  its complete transcript `data/known_examples/bugaenko_vinberg_transcript.json`,
  SHA256 `c4178fc4678bacb48c430b2cec18e16577078824fdd6dea39375e63ed6ae015a`.

These are inherited mathematical inputs. No claim of worldwide novelty,
new compact example, new dimension ceiling, or completed executable
verification accompanies this preparation edition.
