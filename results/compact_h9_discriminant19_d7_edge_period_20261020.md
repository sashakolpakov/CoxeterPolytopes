# The discriminant19 H9 host: a complete D7 edge and chamber period

20 October 2026. **Own complete paper derivation; no new checker or search
has run.** The old 31-wall seed was only a finite, infinite-volume prefix.
Here a different, tiny first-edge argument gives a full-chamber symmetry.
Its product with the known compact-facet reversal has trace `10+14phi>10`.
This proves nonreflectivity of the literal H19 lattice, subject to the
explicit inherited complete Bugaenko-chamber input below. It is not a
K-similarity, commensurability-class or general compact-H9 exclusion.

## 1. Full arrangement, actual base facet and first transverse vertex

Use the native form

    G=[[-2phi,phi],[phi,2]]+C,     phi²=phi+1,
    b=4+phi, D=phi*b=1+5phi,      N(b)=19.

C is the native E8 Cartan matrix with edges `01,12,23,34,45,56,27`.
Write `omega_j=C^-1 e_j`. In particular

    omega0=(4,7,10,8,6,4,2,5),       omega0²=4,
    omega7=(5,10,15,12,9,6,3,8),     C omega7=e7.

The complete primitive norm pool is `{2,2b}`, up to unit squares. This
is the determinant-prime/primitive-covector theorem of the pinned
mixed-norm report: a norm `2d` primitive reflector has `Gr in dO^10`,
hence `d|det G`; O=Z[phi] is a PID and totally positive units are squares.
Both classes occur. In particular this is not a norm-two-only argument.

We first apply the whole-arrangement isolation lemma used in the reviewed
prime59 route, now at the prime `(b)` over19. A long root has all pairings
in bO. Between short and long intersecting mirrors the squared-angle
quantity would have the form `b*a²`. All nonzero possibilities over K are

    1, 2, phi², 3, 2+phi, 2-phi, 3-phi.

Their rational norms are `1,4,1,9,5,1,5`, none divisible by19. This list
includes arbitrary arrangement angles, not just fundamental angles.
For two distinct long roots their reductions are nonzero on the same
one-dimensional radical modulo b. Writing `v2=a*v1+b*u`, their canonical
norms imply `a²=1 mod b` and `(v1,v2)/b=2a mod b`. Thus an intersection
would require its squared-angle quantity to be4 modulo b. The rational
norms of the listed quantities minus4 are `9,4,5,1,1,5,1`; for the
orthogonal value0 the norm is16. Again none is divisible by19. Parallel
distinct mirrors are excluded by anisotropy, since the conjugate form is
positive definite. Consequently distinct long mirrors are ultraparallel,
and a short mirror meeting a long one is perpendicular to it.

Take the long root and its inward side

    r=(1,2;0), r²=2b, Gr=(0,b;0),       (x,r)<=0 iff y<=0.

Its perpendicular lattice is exactly `M=<-2phi>+E8` on y=0. Every ambient
mirror meeting r is a perpendicular short mirror; conversely every
norm-two root of M extends to L. Thus the entire restricted arrangement
is the complete norm-two M arrangement. The inherited, audited Vinberg
calibration identifies its chamber F with Bugaenko's compact eleven-wall
H8 chamber. Specifically, `bugaenko_vinberg_transcript.md` and its JSON
prove the full primitive norm-two alphabet for M, exhaust every earlier
height through the last wall, and combine that order with the independent
41-finite-vertex compactness certificate and Vinberg's stopping criterion.
This is a **full-chamber input**, not the assertion that any compact
root prefix exhausts an arrangement.

Choose the adjacent ambient chamber C on y<0, at
`p_r=(1,0;0)`, whose full inward germ is E8+A1(r). There is no other long
zero mirror by isolation; the short zero germ is the known E8 germ of F.
Let

    c=(0,-1;0), c²=2,      p=(2,-phi;0), p²=-2D.

The edge retaining E8 goes from p_r towards p. The binary roots r,c form
an O-basis, with Gram `[[2b,-b],[-b,2]]`. Its first new mirror must be
E8-orthogonal, because E8 has no connected positive-definite rank-nine
extension. A binary root crossing the open segment is `lambda*r+mu*c`
with lambda*mu<0. For norm2, positivity of the two binary norms gives

    lambda²<=1/b,   lambda'²<=4/[b'(4-b')],
    N(lambda)²<=4/[19(4-b')]<1.

Here `b'=5-phi`, `4-b'=phi-1>3/5`. Thus lambda cannot be nonzero.
For norm2b the pairing with c forces `mu=b*gamma`. Dividing the norm
equation by b gives the identical two bounds on gamma, so gamma cannot
be nonzero. Neither class contains an intervening mirror. Therefore c
is actually the first wall and p has the full actual E8+A1(c) germ.
This proves the full-arrangement starting vertex independently of any
later appended seed walls.

## 2. A complete eight-profile interval at p

Release native E8 node0, retaining `D7(-e1,...,-e7)+A1(c)`. Its edge is

    V(tau)=p+tau*omega0,      tau_light=sqrt(D/2).

The D7 arms have lengths1,1,4. A new short root may be D7-orthogonal and
make angle2,3,5 with c, extend the long arm to D8, extend either fork to
E8, or join the long arm to c to form D9. A new long root is orthogonal
to all whites by isolation. These give exactly eight literal profiles:

|New component|k=-(new,c)|D7 contact j|R|
|---|---:|---|---:|
|isolated short|0|none|8|
|A2 on c|1|none|6|
|I2(5) on c|phi|none|6-2phi|
|isolated long|0|none|8b|
|D8|0|6|4|
|E8|0|1|1|
|E8|0|7|1|
|D9 bridge|1|6|2|

No finite connected component containing D7 has a 5-bond; the finite
classification gives D8,E8,D9 as the only possibilities here. This is a
complete finite first-contact list, not an assumption of sparse contacts.

Write the binary part as `(X/2)p-(k/2)c`, so x=X and
`2y=k-phi*X`. For a contact at j write the spatial part as
`((Z-m_j)/4)*omega0+omega_j`, with m_j=(omega0,omega_j); otherwise it is
`(Z/4)*omega0`. The exact E8 data are

    m_6=2, m_1=7, m_7=5,
    omega_6²-m_6²/4=1,
    omega_1²-m_1²/4=omega_7²-m_7²/4=7/4.

All rows have the equation and collision parameter

    Z²-2D*X²=R,      tau=D*X/Z,       X,Z>0.             (1)

At `(b)`, phi=15 in F19. The R-values8,6-2phi,2 reduce to8,14,2,
all nonsquares. Thus the isolated short, golden c-bond and D9 bridge
rows are impossible at every height. The four short rows left are
A2, D8, and the two E8 placements; the long row is separate.

The known polar root from the pinned 25-wall report supplies the finite
upper

    z=(2phi^4,-phi^5; b*phi²*omega0),       z²=2b,
    Gz=b*(-phi^5,0;phi²*e0).

It is primitive and reflection-integral. Its variables are
`Xz=2phi^4`, `Zz=4b*phi²`, hence its contact is

    tau_z=phi³/2 < sqrt(D/2).

For the strict inequality note `V(tau_z)²=3-2phi<0`. No compatibility
with a cap list is used to treat z as an upper. First-contact comparison
and conjugate positivity in (1) give

    X²<=R*K,     X'²<=R'/a,
    K=phi^8/(2b)=(44+71phi)/38 <17/4,
    a=-2D'=10phi-12 >4.                              (2)

These coarse rational bounds use only `8/5<phi<13/8`.

## 3. All other profiles fail in tiny exact domains

For the D8/E8 rows k=0, binary integrality gives X=2Y, Y in O.

For E8, (2) yields `0<Y<21/20`, `|Y'|<1/4`. In Y=A+Bphi these imply
B=0, since `sqrt5>11/5`. No positive integer has conjugate absolute
value below1/4. Both placements are excluded, without needing their
different spatial congruences.

For D8 the bounds are `0<Y<21/10`, `|Y'|<1/2`. Thus B=0 or1. B=0
fails the conjugate bound. At B=1 that bound forces A=1, giving
`Y=phi²>21/10`, again impossible. This proves the entire D8 exclusion;
there is no unbounded Pell remainder.

For A2, the bounds are `0<X<51/10`, `|X'|<5/4`, so X=A+Bphi has
B=0,1,2. The exact binary congruence `phi*X=1 modulo2O` says both A,B
are odd, because `phi^-1=phi+1 modulo2`. Hence B=1. The conjugate strip
permits A=0 or1, so A=1 and X=phi². Equation(1) now requires

    Z²=2D*phi^4+6=40+56phi,
    N(40+56phi)=704,

which is not an integer square (26²<704<27²). Thus this last short
profile also fails even before its spatial divisibility test.

For the long row reflection integrality against E8 forces
`Z/(4b)=t in O`; binary integrality gives X=2A, A in O. Its norm equation
is exactly

    2b*t²-phi*A²=1.                                   (3)

The interval bound is A<=phi^4, equivalently `0<t<=phi²`. Conjugate
positivity gives `2b'*t'²<=1`, hence `|t'|<2/5` since b'>27/8.
If t=u+vphi, the two strips force v=0 or1. The integer case fails;
v=1 forces u=1. Thus t=phi². Equation(3) gives A²=phi^8, and A>0
selects A=phi^4. Reconstruction is exactly the upper z.

Every profile in the finite interval has now been addressed. Local
finiteness of the arrangement and the finite upper imply that z is the
unique actual first simple wall. The new future endpoint

    V=p+(phi³/2)*omega0,       V²=3-2phi<0

has full inward germ `D7+c+z`. Its positive definite rank-nine diagram
has no unresolved simultaneous new wall. This establishes an actual
full-chamber edge, not just an admissible appended polar cap.

## 4. A new chamber involution from two full actual germs

Define

    P=(phi²,-phi;phi²*omega0),       P²=-2phi,
    J_P(x)=-x-(P,x)P/phi.

Since phi is a unit and P is integral, J_P is an integral Lorentz
involution fixing P and preserving the future cone. Direct products give

    (P,r)=-D,    (P,c)=-1,
    J_P r=z,    J_P c=g,
    g=old11=(phi,0;phi*omega0).

The root g is one of the eleven actual Bugaenko side walls of r. The
rank-nine positive system `D7+r+g` therefore defines an actual ambient
finite vertex, not merely a seed vertex. More explicitly its future ray is

    B=(1,0;(phi/2)*omega0),       B²=1-phi<0.

It vanishes on r,g and D7. These are actual walls, and projection of an
interior point onto their positive definite Gram span stays feasible
(inverse-Gram entries are nonnegative); its remaining wall inequalities
are strict. Thus the complete inward germ at B is `D7+r+g`.

Let w0 be the longest element of the retained D7 Weyl group. It fixes
its perpendicular space and satisfies

    w0 e1=-e7, w0 e7=-e1, w0 ei=-ei (2<=i<=6),
    w0 e0=e0+3e1+5e2+4e3+3e4+2e5+e6+3e7.

It fixes the two binary basis vectors. The last line follows by fixing
omega0; all columns are integral. Put H=J_P*w0. Since w0 fixes P, the
two involutions commute. On the outward D7 whites H induces the fork
swap1<->7, and it maps z to r and c to g. Moreover

    (P,V)=-phi,      H V=P-V=(phi-1)B.

Hence H maps the full inward actual germ at V to the full inward actual
germ at B. It preserves the entire primitive mirror arrangement, so
unique continuation of the same chamber proves `H(C)=C`. This is not
an unsupported symmetry of an incomplete 31-wall seed. H has genuine
E8/binary mixing and is not the earlier binary Pell search seed A19.

## 5. The facet reversal and an infinite-order product

The inherited Bugaenko reversal of F extends integrally to L fixing r.
Here are explicit native columns, eliminating any index-two glue ambiguity:

    A e_x=(9+16phi,0;(6+10phi)*omega7),
    A e_y=(-4-8phi,1;-(3+5phi)*omega7),
    A e0=e6, A e1=e5, A e2=e4, A e3=e3,
    A e4=e2, A e5=e1, A e6=e0, A e7=-k.

The actual Bugaenko wall k is

    k=(2+3phi,0;
       6+10phi,12+20phi,18+30phi,15+24phi,
       12+18phi,8+12phi,4+6phi,9+16phi).

These columns are the known facet isometry on y=0 and fix r=e_x+2e_y;
thus they preserve the full Gram form. They are integral and A²=I.
They permute all eleven outward walls of F: the three simple pairs
0<->6,1<->5,2<->4; node3 fixed; `-e7<->k`; and `g<->old9`, where
`old9=(1,0;phi*(2,4,6,5,4,3,2,3))`. The y coordinate is unchanged,
so the occupied side y<0 is preserved. Whole-arrangement facet adjacency
then proves `A(C)=C`. The explicit columns verify the extension rather
than assuming that every intrinsic facet symmetry extends to the host.

We calculate the trace from these columns, without a numerical program:

    trace(A w0)=4+6phi,
    (AP)_time=8+13phi,  (AP)_y=-phi,
    (AP)_spatial0=27+42phi,
    (P,AP)=-20-34phi.

For the first trace, the only nonzero diagonal contributions are
`9+16phi`, `1`, `-17-30phi`, `12+20phi`, `-1`. For the last pairing,
`GP=(-phi^5,1;phi²*e0)`, so it is
`-phi^12-phi+phi²*(27+42phi)=-20-34phi`.
Since w0 fixes P,

    trace(AH)=-trace(Aw0)-(P,AP)/phi
             =10+14phi > 10.                            (4)

An order-finite real ten-dimensional matrix has every eigenvalue on the
unit circle and trace of absolute value at most10. Thus AH has infinite
order, and both of its factors preserve the full chamber C.

The full integral reflection group W acts simply transitively on its
chambers. Its intersection with the chamber stabilizer is trivial, so
the powers of AH represent infinitely many cosets of W in O^+(L).
The latter is a uniform arithmetic lattice (admissibility and anisotropy
are the standard arithmetic-lattice input). Infinite index forces C to
have infinite hyperbolic volume. In particular no compact Coxeter
polytope with all wall reflections preserving this literal L exists.
This does not exclude other integral lattices on similar or commensurable
forms, and does not change any general dimension bound.

## 6. Provenance and exact scope

The former polar z is reused as a proved finite upper, but its status as
the first full-arrangement wall is established here. The later 31-wall
seed, its 89-vertex lower bound, its last escape and its bounded 62-profile
search are not premises of this proof. Thus there is no height-loop or
append-only completeness assumption.

The inherited sources and current byte pins are:

|Source|SHA256|
|---|---|
|`compact_h9_discriminant19_mixed_norm_seed_20261019.md`|`70fa010732d837d5c044701b4c18ed5e68cb36c18cc0954dc7d7f369250a5c87`|
|`compact_h9_discriminant19_polar_caps_research_20261019.md`|`792634e935866ed78d95d0a93c04fdd7a8484f3605b5c50a66c1d376c9599f1b`|
|`compact_h9_golden_discriminant59_edge_route_20261020.md`|`725d8775cebec51d8c18b4b0a0d3be0255e8b1472621adb4d3cefa656adb429e`|
|`og10_golden_h9_infinite_stabilizer.md`|`1d9a71240d499e498395f5c1c7eade7c2f307f2fed5e44c2578e02b1357392a1`|
|`data/known_examples/bugaenko_compact_h8.json`|`c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c`|
|`bugaenko_vinberg_transcript.md`|`113e54694c6e2c05a2de4b8bdd7a5a072aabd2572aba760ecb8d4e039ebda719`|
|`data/known_examples/bugaenko_vinberg_transcript.json`|`c4178fc4678bacb48c430b2cec18e16577078824fdd6dea39375e63ed6ae015a`|
|`data/candidates/h19_vinberg_prefix.json`|`91c16b581e5763d6a283b6818f0e087e941fbb987b6b05a7e283febf1e32d0f1`|

The first-contact classification, full facet/germ continuation, and
arithmetic index-to-volume implications are written mathematical layers.
A future tiny exact checker may verify the fixed field identities,
profiles, finite strips, columns, Gram and trace, but does not formalize
these universal geometric/group-theoretic layers. No checker, census,
Lean, TeX or numerical job was run in producing this note.
