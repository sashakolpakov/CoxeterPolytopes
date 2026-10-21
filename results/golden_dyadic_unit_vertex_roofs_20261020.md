# Forced roofs over unit vertices of the full M6 face

20 October 2026. Paper proof and a finite next construction interface.
No new worker, enumeration, compactness claim or frozen-source edit.

Import the frozen full M6 face from
`results/golden_dyadic_d4_research_20261020_freeze.json`, the fixed
sixteen-wall coordinates from
`golden_dyadic_sixteen_wall_caps_20261020.md`
(SHA256 `c288b5de619471772c8381d98afa3996ac8851b4d2d146ad520d6474ee8ba92b`),
and the exact P3 basis/finite-angle argument of
`golden_dyadic_triple_h3_compact_face_20261020.md`
(SHA256 `665a9861ae70634b20a80fa93007276fb607e254e84a5c3d637d337adfee46a8`).
Notation: `R5=M6 orthogonal-sum P3`, the actual norm4 H3 roots are
(a3,s,r), and `wP=-2omega1=a2-V0` has square phi+2. All coordinates
below use this fixed decomposition.

## 1. A uniform actual-wall lemma

Let T be ANY integral future actual M6 vertex vector with `T^2=-phi`.
Then

    b_T=T+wP                                              (1)

is an ACTUAL ambient facet normal: it is the unique first simple wall
on the edge releasing a3 from the full vertex `(M6 germ at T)+H3`.
This conclusion does not need a classification of the spectator germ.

Indeed T has unit norm, giving an integral orthogonal splitting
`M6=O*T+N`, with N positive definite at both real embeddings. The
norm2 integral vector(1) is an upper on the ray `T-tau*omega1` at

    tau_b=2phi/(phi+2)=2/sqrt(5),

strictly before light. Consider any first actual wall no later than
this upper. Write its canonical primitive norm and integral splitting
as

    n=a*T+z+p,  q=n^2 in{1,2,4},
    a in O, a>0, z in N, p in P3, epsilon=z^2.

The strict positivity of a follows from the full actual starting germ.
Put L=p.wP. Collision before b gives `L/a>=phi+2`. Positive Cauchy
in P3 and the norm equation give the COMPLETE interval bounds

    2a^2+epsilon<=q,
    a'^2<=phi*(q-epsilon').                           (2)

No root-height cutoff is imposed. For a=A+Bphi,
`Tr(a^2)=2A^2+2AB+3B^2`.

For q1, (2) implies trace at most `1/2+phi<3`; only a=0,+1,-1
have such a trace, and none has `0<a^2<=1/2`. For q4, primitive
reflection integrality gives `n.T in2O`, hence a in2O. Its a/2
satisfies the same impossible q1 bounds. These arguments use the
entire imported primitive norm pool, not selected finite profiles.

For q2, the first inequality(2), with a>0, gives `0<=epsilon<2`.
At the other embedding the nonzero a' has a strictly positive
contribution, giving `0<=epsilon'<2`. An element of O in both these
intervals is0 or1: its phi-coefficient B satisfies
`|B|*sqrt(5)=|epsilon-epsilon'|<2`, hence B0. If epsilon1, (2)
again gives the impossible q1 bounds on a. Therefore epsilon0 and
z0. This is the crucial uniform spectator exclusion.

A new norm2 root has a finite nonorthogonal product with either
retained norm4 I2(5) root only at angle pi/4, product-2. A4-bond
adjoining that5-bond is already nonfinite. Thus p is orthogonal to
s,r. The exact P3 basis(a3,s,wP) gives `p=m*wP`, m in O.
Now (2) leaves only `a=1,phi^-1` in the q2 strip, since its trace
is below5. The norm equation reads

    (phi+2)m^2=2+phi*a^2.

The second a-value gives `m^2=phi^2/(phi+2)`, of nonintegral field
norm1/5. The first gives m^2=1, and positive L forces m1. Hence
n=b_T is the unique first simple wall. All simultaneous candidates
and arbitrary finite spectator contacts were included in(2).

The exact same statement holds for any unimodular golden M_d face
with this integral P3 splitting and primitive norm pool. For this
note only the already identified full M6 chamber is used.

## 2. The seventeen-wall prefix still escapes

Let b be the actual roof of the triple-H3 time T0 from the imported
face proof. Add b to the sixteen fixed walls. At the different
actual M6 vertex

    B=(1;(phi^-1)^4,phi^-4,0), B^2=-2phi^-7,

take the H3 component `-phi^-3*omega3`. The resulting vector Y has

    Y^2=-2phi^-7+(3/4)phi^-4=(183-113phi)/4>0.

For example phi<89/55<183/113 gives the strict sign exactly. Its
M6 active diagram is H4(u,beta1,beta2,beta3)+A1(w)+A1(beta6).
All sixteen old wall products and the new b product are:

    u=w=beta1=beta2=beta3=beta6=0;
    beta4=-2phi^-3, beta5=-phi^-4, h=-2;
    a3=s=0, r=-phi^-3;
    a2=1/2-phi, A(a2)=1/2-phi^4;
    W=-phi, Z=-phi^2+phi/2,
    b=1/2-2phi^-1.

They are all nonpositive. The strict interior point of the fixed
sixteen-wall model is also strict for b, so the17-wall cone is
noncompact; a full DD enumeration is unnecessary. This statement
concerns that finite prefix only. No full-lattice exclusion follows.

## 3. A new actual wall cuts that escape

The same actual M6 chamber has the different vertex time

    p0=(phi^3;(phi^2)^4,0,0), p0^2=-phi.

Its active full M6 diagram is
`H4(u,beta1,beta2,beta3)+B2(beta5,beta6)`; beta4, w and h have
strict values `-phi^2,-1,-phi^5`. Thus the uniform lemma applies.
The forced new actual roof is

    b0=p0+wP,  b0^2=2.                              (3)

Its only nonzero products with the17 known normals are

    beta4:-phi^2; w:-1; h:-phi^5; a3:-2;
    a2:-2phi; A(a2):-8phi^3;
    W:-2phi^3; Z:-3phi^2; b:-2phi.

These are all native or strictly dotted. Most importantly,

    b0.Y=1/2-phi^-2>0.

Thus(3) removes the exact remaining escape by a proved ACTUAL edge,
rather than a guessed cap. The integral involution A supplies the
partner A(b0) automatically; A fixes wP and exchanges M6 vertices.

The H4 at B is fully sealed because no finite irreducible extension
of H4 exists in rank5. Its root lattice is unimodular, so its exact
perpendicular lattice splits integrally. In the M6 coordinates it is
`M2=O*p0+O*e5+O*e6`, with full native triangle
`w=(phi-1;phi,0), beta5=-e5+e6, beta6=-e6`. Adding P3 identifies
the full actual5-face lattice with literal R1. This lattice identity
does not by itself identify all walls of that5-face or import a
metric chamber reversal as an ambient symmetry.

## 4. The complete unit-vertex list is finite and explicit

There are exactly FIVE integral norm-minus-phi times among the23
vertices of the full M6 chamber:

    U={V0,A(V0),p0,A(p0),T0},
    A(p0)=(10phi+7;(6phi+4)^2,(5phi+3)^4).

Here is a complete hand inventory, rather than a vertex enumeration
job. A finite6-subset of the nine native walls cannot contain both
beta6,h (they are dotted). If it contains beta6 but not h, it must
omit exactly two further walls from(u,beta1,...,beta5,w). If beta5
is present then w must be omitted, since the(4,5)3-chain is already
nonfinite. The complete finite cases are therefore the following10:

|two omitted besides h|determinant of the six active roots|unit time|
|---|---|---|
|w,u|1|V0|
|w,beta1|2|none|
|w,beta2|3-phi|none|
|w,beta3|2phi^-2|none|
|w,beta4|phi^-4|p0|
|beta5,u|10|none|
|beta5,beta1|16|none, as below|
|beta5,beta2|6(3-phi)|none|
|beta5,beta3|8phi^-2|none|
|beta5,beta4|2phi^-4|none|

Their finite types follow by cutting the native path, respectively:
B6; A1+B5; I2(5)+B4; H3+B3; H4+B2; A4+A1+A1_short;
A3+2A1+A1_short; I2(5)+A2+A1+A1_short;
H3+2A1+A1_short; H4+A1+A1_short.

All determinant rows except1,phi^-4,16 have an odd prime-ideal
valuation. If an integral norm-minus-phi time existed, its unimodular
positive complement would contain these six roots with determinant
equal to an index SQUARE times a unit. The odd valuations forbid it.
For determinant16 the ray is explicitly

    (1;phi^-1,1/2,1/2,1/2,1/2,0), norm=-phi^-3.

Its only positive norm-minus-phi normalization scales by phi^2 and
still has half-integral spatial entries, so is not an integral time.

The case h present and beta6 absent gives the10 reflected vertices
under the known integral involution A. If both beta6,h are omitted,
the remaining diagram is the seven-node path

    u --5-- beta1 --3-- beta2 --3-- beta3
        --3-- beta4 --3-- beta5 --5-- w.

A six-subset is finite only when the omitted path node is beta2,
beta3 or beta4. The central case is H3+H3 and gives T0. The two
others are H4+I2(5), determinant `phi^-4(3-phi)`, with odd valuation
at the prime above5. They admit no integral unit time. This proves
all23 cases and the complete list U without any unbounded search.

Thus there are exactly five forced roofs. Their ten products, in
the order(a2,A(a2),b,b0,A(b0)), are

    a2: -18phi-10,-4phi-1,-2phi,-8phi^3;
    A(a2): -4phi-1,-8phi^3,-2phi;
    b: -2phi,-2phi;
    b0: -6phi^3.

They are pairwise strictly ultraparallel, as required for consistency
with their proved actualness. This supplies a complete finite roof
set, not an invitation to append further unit-height guesses.

## 5. Even the complete five-roof finite cone has an exact escape

Add b,b0,A(b0) to the original sixteen walls, for19 in total. No DD
job is required to decide this gate either. In the original R5
coordinates the integral spacelike vector is especially small:

    X=(2phi-1;(phi^2)^3,0,0;-1,0,0,-1)
     =(V0+p0)-4omega3,
    X^2=4(phi+2)>0.

The products with all19 walls are

    beta1=beta2=beta3=beta5=beta6=a3=s=0,
    beta4=-phi^2, u=-1, w=-2phi^2, h=-11phi-7,
    r=-4, a2=b0=0,
    A(a2)=-4phi^6, b=-4phi-3, A(b0)=-2phi^7,
    W=-4phi^2, Z=-2phi^4.

Thus even ALL forced unit-vertex roofs leave an exact spacelike
escape. This is a negative answer for the finite19-wall candidate,
not for the full R5 arrangement.

The surviving support has an actual geometric interpretation. V0
and p0 are adjacent actual M6 vertices, with common rank5 white
system `A3(beta1,beta2,beta3)+B2(beta5,beta6)`. Together with the
actual norm4 A2(a3,s), this defines an actual2-face. Its old side
r joins those two vertices; the known side walls beta4,u meet r
at its ends, and their farther corners with a2,b0 are finite by
the actual positive acute projection. The two farther roofs a2,b0
are ultraparallel, leaving a genuine missing closure arc in this
finite model.

The common white A3+B2+A2 is NOT claimed sealed. New actual walls
may have finite contacts to its components, and the diagonal lattice
of its perpendicular K-space alone does not describe its induced
arrangement. Its new closure must be studied as a complete actual
edge or by another full geometric argument. This is a concrete
next interface, rather than an unbounded cap-prefix proposal.
No new example, global period or full-host obstruction is claimed.
