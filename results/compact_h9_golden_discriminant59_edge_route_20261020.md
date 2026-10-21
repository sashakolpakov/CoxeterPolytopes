# A discriminant-59 extension of the compact Bugaenko facet

20 October 2026. **New paper construction route, not a reflective lattice
or a compact H9 example.** No search, import, checker, optimizer, Lean job
or build was run. The first two transverse edges below have complete
paper proofs; they do not require an increasing-height prefix.

Put O=Z[phi], phi^2=phi+1, and use the literal E8 Cartan matrix and native
node order of `data/known_examples/bugaenko_compact_h8.json`. Consider

    L59 = E8(O) + B,
    B = [[-2phi,-1],[-1,4+2phi]],
    delta=5+12phi,       b=7+5phi,       delta=phi*b.

The coordinates are (x,y;z), with z in the eight-dimensional E8 module.
The physical binary block has negative determinant -delta; its conjugate
is positive definite. Hence L59 is admissible of signature (9,1), and
has no nonzero K-rational isotropic vector. Since N(b)=59 and
phi=34 mod(b), (b) is prime. The determinant has valuation one there, so
the lattice has no proper integral overlattice.

This is separated by its odd prime-59 determinant valuation from the
inspected excluded golden prime-11, -55 and -131 hosts, the excluded
unimodular/old local-Witt hosts, and the live discriminant-19 host. It is
still an E8-plus-binary construction, not a claim to escape that entire
architecture. A bounded text screen found no earlier literal
discriminant-59 host in the inspected repository notes. No literature-wide
absence or novelty claim is made.

## 1. Complete norm pool and isolated long mirrors

The lattice is even. If a primitive integral reflector has square 2d,
reflection integrality gives Gr/d integral. The adjugate and a Bezout
relation among its primitive coordinates give d|det(G). Since O is a PID
and all totally positive units are squares, the complete necessary
primitive norm pool, modulo unit squares, is

    2, 2b.

Both classes occur below. These squarefree norm ideals also imply that
every integral vector of either displayed norm is primitive. The ring
facts here are the same elementary golden-ring facts used by the pinned
earlier norm-pool proofs; no incomplete root shell is used.

There is a useful uniform observation. In an even golden lattice with
determinant ideal (b) of prime norm p not in {2,3,5}, and with this norm
pool, all long mirrors (norm 2b) are mutually ultraparallel and meet short
mirrors (norm 2) only orthogonally. Indeed a long root v satisfies
(v,L) subset bO. With a short root, 4cos^2(theta)=b a^2. For adjacent
walls of a fundamental chamber, the complete finite-angle possibilities
over Q(sqrt(5)) are

    4cos^2(pi/m) in {0,1,2,phi^2,3,2+phi}.

For arbitrary intersecting arrangement mirrors one must also include
2-phi and 3-phi, arising from angles 2pi/5 and 3pi/10. Their rational
norms are respectively 1 and 5. Thus all seven nonzero possibilities
have norms indivisible by p. For two
long primitive roots v,w, the reduction of each lies in the same
one-dimensional radical modulo b. Write w=a v+b u. Reducing w^2/b=2
modulo b gives a^2=1, while (v,w)/b=2a modulo b. Thus their finite angle
would have 4cos^2(theta)=4 modulo b. The norms of each displayed value
minus 4 involve only 2,3,5, a contradiction. Parallel mirrors would
give a nonzero K-isotropic vector and are excluded by admissibility.
This proves a whole-arrangement statement for L59; it does not assert
that every long orthogonal lattice in every such host is isometric.

## 2. An actual compact Bugaenko facet

Take

    r=(phi-1,-2;0),       r^2=2b,       Gr=(0,-b;0).

Its reflection is integral: (x,y;z) maps to
(x+(phi-1)y,-y;z). Its mirror is y=0, whose entire orthogonal lattice is
the literal calibrated E8+<-2phi>. Every ambient mirror crossing r is
orthogonal by Section 1. Conversely every norm-two reflector of this
orthogonal lattice extends to L59. Thus its full intrinsic arrangement
is exactly the known compact eleven-facet Bugaenko arrangement, rather
than a finite wall prefix.

Let p=(1,0;0). Its full zero chamber is E8+A1(r). For a short zero root,
the scalar part has energy 2b a^2 at both embeddings and cannot be
nonzero since N(b)>1. For a long zero root, reflection integrality puts
its E8 part in bE8; dividing the norm equation by b gives an E8 energy
with rational norm at most 4/59 unless it is zero. The two scalar roots
are then +/-r. The remaining zero roots are exactly the calibrated
ordinary E8 norm-two system. This fixes the full-chamber sector C at p.

All eleven Bugaenko side walls are therefore actual walls of C, and all
41 vertices of the compact r facet are actual finite vertices of C.

## 3. The complete first transverse edge: r to a short wall

Put

    c=(1,1;0),       c^2=2,       Gc=(-1-2phi,3+2phi;0).

The determinant of the binary basis (r,c) is phi^2, a unit. Its Gram is

    [[2b,-b],[-b,2]].

It bounds a nonempty finite segment in H1, from p to

    h=(3+2phi,1+2phi;0),       h^2=-2delta,       (h,c)=0.

A first mirror before h must be orthogonal to the retained whole E8,
because E8 has no finite rank-nine connected extension. It is therefore
binary. Write its primitive normal as lambda*r+mu*c. A mirror cutting
the interior of this segment must have lambda*mu<0.

For norm 2, the physical norm equation gives lambda^2<=1/b. The positive
conjugate Gram gives

    lambda'^2 <= 4/[b'(4-b')].

Consequently N(lambda)^2<=4/[59(4-b')]<1, impossible for nonzero lambda.
Here b'=12-5phi, and 59(4-b')>4 follows already from phi>21/13.

For norm 2b, integrality of its pairing with c implies mu=b gamma.
The same argument gives gamma^2<=1/b and
gamma'^2<=4/[b'(4-b')], again impossible. Thus no intervening integral
mirror exists. The first endpoint is exactly h, of type E8+A1(c).
This is an all-height first-edge result, not a finite search cutoff.

The positive-height Bugaenko walls have x=1, phi, 2+3phi. Their products
with c are -(1+2phi)x, all strictly ultraparallel; the eight zero walls
are orthogonal to c. Thus the eleven old walls plus r,c form a valid
thirteen-wall seed. The compact r facet is retained, and h is an
additional finite vertex. A strict interior follows by perturbing an
interior point of y=0 slightly toward y>0; acute projection verifies all
thirteen actual facets independently of the full-arrangement argument.

## 4. One further COMPLETE finite edge at h

Release native E8 end node 7, retaining E7 on nodes 0,...,6 and c. Put

    u=(phi^2,phi;phi^2*omega7),
    omega7=(2,3,4,6,5,4,3,2),       omega7^2=2.

Then u^2=2, (u,c)=-1, and u is orthogonal to E7. The edge is
h+tau*omega7. Its proposed first parameter is

    tau_star=delta/(2phi^2)=(7phi-2)/2 < sqrt(delta).

The complete finite-extension inventory is only five profiles. A short
root can be isolated, join c by a 3- or 5-bond, or extend E7 to E8 by the
unique 3-bond at node 6 while staying orthogonal to c. A long root must
be orthogonal to all eight short whites by Section 1. No connected
rank-nine spherical component can contain E7 and the c component.

Write a new root's binary part as (X/2)h-(k/2)c, with k=0,1,phi. For
the E7-orthogonal profiles write its E8 part as (Z/2)omega7. For the E8
extension write it as ((Z-3)/2)omega7+omega6, where omega6^2=6 and
(omega6,omega7)=3. These are integral variables, rather than merely
field coordinates: since h_x-h_y=2 and c_x-c_y=0, X=x-y belongs to O.
The unimodular E8 dual coordinates put Z/2 in O in the orthogonal
profiles and (Z-3)/2 in O in the extension profile, so Z belongs to O
in each case. Every profile reduces to

    Z^2-delta*X^2=R,       X>0, Z>0,

with respective R=4,3,3-phi,1,4b. Contact by tau_star implies

    X^2 <= R/3,       X'^2 <= R'/(-delta'),

because 4phi^4-delta=3. Here 2<-delta'<5/2 and delta<25.

For R=1 or 3-phi these bounds give N(X)^2<1, excluding the profile.
For R=3 or 4 they give |N(X)|<2; X is therefore a unit. Both |X| and
|X'| are strictly less than phi except that X<=1 for R=3. The golden
unit group consequently forces X=1 in either case. At R=4, Z^2=9+12phi
has rational norm 45, so is not a field square. At R=3, Z=2phi^2 is the
unique positive square root, giving exactly the displayed u.

For R=4b, reflection integrality against E8 gives t=Z/(2b) in O. The
two embedding bounds imply

    N(t)^2 <= (3+delta)/(3*59) < 1,

contradicting t nonzero. This exhausts all five profiles, without a Pell
orbit computation or coefficient enumeration. Hence u is the unique
first simple wall, and its endpoint has whole type E7+A2.

Post-selection compatibility is explicit: u pairs with r by -delta,
with c by -1, with the eight E8 walls by zero except node 7 (-phi^2),
and with the three positive-height Bugaenko walls x=1,phi,2+3phi by

    -phi,       -3-3phi,       -9-15phi.

These are right, pi/3, pi/5 or strictly ultraparallel relations. No old
cap compatibility was used to discard a first-contact candidate.

## 5. A compact actual pentagonal two-face

The old Bugaenko wall v=(1,0;phi*omega7) is orthogonal to E7 and r.
The calculation above gives (u,v)=-phi, while c,v and r,u are
ultraparallel. The E7 face is therefore the pentagon bounded in order
by the projection of native node 7, c,u,v,r. Its angles are

    pi/2, pi/3, pi/5, pi/2, pi/2.

Four vertices are already actual by the Bugaenko facet and the two
complete edges. The fifth is actual because the retained E7 plus u,v
has positive-definite Gram E7+I2(5). Projection of an actual interior
point onto these nine actual walls stays feasible and strict on every
other wall: the inverse Gram is nonnegative and all wall products are
nonpositive. Thus the fifth vertex is not a prefix-only proposal.
Convexity then identifies the whole E7 face with this compact pentagon.
In particular the present construction supplies at least 44 actual
finite vertices (41 on r, plus h and the two off-r new pentagon vertices).
It does not close the other outgoing edges of the full H9 chamber.

## 6. The inherited reversal extends, but supplies no period

The Bugaenko diagram reversal swaps indices i<->8-i for i=0,...,8 and
9<->10. Its intrinsic isometry S is integral. In the literal E8 basis,

    Sx=(9+16phi)x+(6+10phi)omega1,
    omega1=(5,8,10,15,12,9,6,3),

and Se0=e7, Se2=e6, Se3=e5, Se4=e4, Se5=e3, Se6=e2, Se7=e0,
while Se1 is minus the archived eleventh Bugaenko wall. These identities
also follow directly from the two E8 corner frames and their equal
norm -2phi. Since Sx-x is in twice the intrinsic lattice, S extends to
L59 fixing r through the actual glue

    y=((phi-1)x-r)/2.

It preserves C by its inward germ along the complete r facet. However
S^2=1. This is one chamber involution, not an infinite-order stabilizer.
There is no K-isometry exchanging r and c up to mirror scaling because
b is not a K-square (N(b)=59). The new E7+A2 endpoint has no integral
germ swap: exchanging c,u would reflect in c-u, whose norm is 6 and
whose E8 covector includes the unit -phi^2 at node 7, not divisible by 3.
Thus the obvious A2 swap fails reflection integrality. These checks
prevent importing an old two-involution translation without its glue.

## 7. Exact scope of the remaining obstruction

The short c facet cannot be decided merely by replacing it with its
orthogonal integral lattice E8+<-2delta>. There are ambient integral
norm-two mirrors meeting c nonorthogonally: for every E8 root a of norm
two,

    w=(-2-phi,-phi^2;phi^2*a),       w^2=2,       (w,c)=-1.

Their restrictions have norm 3/2 and are not the full orthogonal-lattice
root pool. An intrinsic-facet nonreflectivity calculation that omitted
them would not exclude L59. The current proof neither proves nor
assumes that this short metric facet is itself Coxeter.

This is the reason to proceed by complete actual finite edges, rather
than declaring a scalar-lattice face nonreflective. The next useful
check is a small exact algebra/first-edge certificate for the identities
and five-profile elimination above, followed only by a genuinely new
edge direction if the parent authorizes it. No height-prefix loop is
proposed. Reflectivity and compactness of the full host remain open.

## Source boundary

The calibrated Bugaenko data and complete transcript are inherited:

* `data/known_examples/bugaenko_compact_h8.json`, SHA256
  `c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c`;
* `data/known_examples/bugaenko_vinberg_transcript.json`, SHA256
  `c4178fc4678bacb48c430b2cec18e16577078824fdd6dea39375e63ed6ae015a`.

The inspected earlier host exclusions are the prime-11, -55 and -131
corner-symmetry reports and the local-Witt symmetry review. Their scope
does not extend to this determinant valuation. The focused primary
screen of [McLeod, Section 3.2](https://etheses.durham.ac.uk/id/eprint/7743/1/thesis31072013.pdf)
covers the displayed diagonal phi and +/-3+2sqrt(5) families, not this
integral prime-59 form. The [Bugaenko 1984 bibliographic record](https://www.mathnet.ru/eng/vmumm3408)
states a unimodular scope, which likewise does not classify L59. These
observations screen two specific potential dead-host transfers; they
are not a complete published classification search.
