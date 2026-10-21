# An integral chamber period for the golden dyadic dimension-nine host

20 October 2026. Complete paper candidate, pending root/outside review
and fixed exact verification. No worker, build, prefix or new example.

This concerns ONLY the literal lattice

    R5=<-phi>+I5+D4(O), O=Z[phi], phi^2=phi+1.

It is not a statement about every lattice in its K-similarity or
commensurability class. Import the frozen actual M6/H3 packet and
the following two source-bound complete actual-wall proofs:

* `golden_dyadic_unit_vertex_roofs_20261020.md`, SHA256
  `e0c6b78e7c9c9ece0b77ef75cdfdcdc54d25c5177b1b474a3bcbead9515653e6`;
* `golden_dyadic_a3_b2_b3_actual_edge_20261020.md`, SHA256
  `400ba9e7fb1369816a9cdfa10d51a22a75a5babec0a3b5aa7f128e2fbb9327dc`.

All paths above are in `results/`. Their actualness and full-germ
claims, not just their coordinate identities, are used here. Use
the integral orthogonal splitting R5=M6+P3 and notation of those
sources. In particular wP=-2omega1, b0=p0+wP is actual, and

    N=(phi^3;phi,phi,phi,phi,0,0)-2phi*omega3

is the new actual norm2 wall from the complete17-profile edge.

## 1. Two actual H4 faces

The old positive norm2 H4 frame and a new frame are

    E=(u,beta1,beta2,beta3),
    E'=(N,beta4,beta3,beta2),

both in the chain order5,3,3. Their common Gram C has diagonal2
and successive off-diagonal entries(-phi,-1,-1), with
`det C=phi^-4`. Thus both split integrally from R5. All four roots
in each frame are ACTUAL and their Gram is positive definite; the
actual positive projection guarantees their nonempty actual face.

An H4 cannot have a nonorthogonal finite extension: the finite
classification has no irreducible component of rank at least5
containing H4. This includes nonparabolic overgroups. Hence each
H4 face is sealed for the FULL arrangement. Its exact perpendicular
lattice gives its full intrinsic arrangement, not just a projected
native prefix. The known roots

    B3(a2,a3,s)+A1(beta6)

are orthogonal to the new E' and are actual. Their norms are(2,4,4)
and1. Together with E' they form a positive definite rank8 system W.
The corresponding actual1-face will be used below.

## 2. A primitive binary plane in the new face

For any original integral vector x, write its residue perpendicular
to E' as `x-proj_E'(x)`. The projection is integral since C is
unimodular. Put

    zeta=e1-proj_E'(e1), eta=e5-proj_E'(e5),
    P1=eta-zeta,
    g=2phi^-1*eta-zeta.

These vectors are orthogonal to E', a2,a3,s and e6. Exact contractions
give

    zeta^2=-32phi-19, eta^2=-18phi-11,
    zeta.eta=-24phi-15,
    P1^2=-2phi, g^2=1, P1.g=0.

The signs are fixed by `(V0,P1)=-phi^9<0`, so P1 is future.
No primitivity or time normalization is inferred from a determinant.
The following explicit formulas give integral coordinates:

    L=3beta4+2beta3+beta2,
    P1=-e1+e5+phi^5*N+phi^3*L,
    g=(2phi^4;(phi^3)^5,0)-4phi^2*omega3.

In original(time;I5;D4) coordinates,

    g=(5phi+2;(phi^3)^4,0;-phi^4,0,0,-phi^4).

Its unit norm proves primitivity and reflection integrality. The
norm2 primitive reflector

    Q=g+e6                                              (1)

has product-1 with beta6=-e6 and is orthogonal to E' and a2,a3,s.
Its unit last I-coordinate proves primitivity.

For an independent projection calculation, C^-1 has upper triangle

    4phi^4,3phi^5,2phi^5,phi^5;
    6phi^4,4phi^4,2phi^4;
    2(3-phi)phi^4,(3-phi)phi^4;
    2phi^2.

Using the E' pairing vectors `(phi,0,0,0)` and `(0,1,0,0)` for e1,e5
gives every displayed contraction and the two residue formulas.

## 3. A complete finite interval from an actual feasible point

Project the known actual vertex V0 orthogonally to W. The actual
positive acute projection is future timelike and feasible for EVERY
actual wall. In the displayed binary basis it is exactly

    S=(phi^8/2)*P1-2phi^5*g,
    S proportional to P1-t_s*g, t_s=4phi^-3.             (2)

The segment from S to P1 is timelike. It is a concrete finite segment,
not an unbounded height interval. Any actual wall cutting it before
P1 makes a positive definite completion of the whole W.

By the full H4 seal its normal is orthogonal to E'. There are exactly
SIX finite weighted profiles at B3(a2,a3,s)+A1(beta6):

|role|canonical norm q|c=(r,e6)|B3 contact|positive white square kappa|
|---|---:|---:|---|---:|
|isolated|1|0|none|0|
|isolated|2|0|none|0|
|isolated|4|0|none|0|
|A1 to B2|2|1|none|0|
|B3 to F4|2|0|a2:-1|3/2|
|B3 to B4|4|0|s:-2|2|

These use all primitive norm classes1,2,4. The B3 end extensions
are the same complete weighted inventory as in the preceding edge.
The isolated norm1 A1 can only have a norm2 new4-bond. A bridge
between it and B3 has two4-bonds in rank5, so is nonfinite. Multiple
contacts or a middle contact cannot be finite. Thus no other finite
native assignment or spectator contact is omitted.

Write a hypothetical first normal as

    r=a*P1+b*g+c*e6+v,
    v in the B3 K-space, v^2=kappa,
    q=-2phi*a^2+b^2+c^2+kappa.                         (3)

Here A=2a is in O since P1 is integral and has norm-2phi; b,c
are in O since g,e6 have unit norm. For q4 reflection integrality
further gives a in O and b in2O. The values of c and kappa are
those in the table; no assertion that v is integral separately is
needed.

If a>0 then the wall has a negative value at P1; feasibility at S
implies feasibility on the whole segment. Only a<0 can cut earlier.
Such a wall has b>0 and collision

    0<t_*=-2phi*a/b<=t_s.

Consequently b>=(phi^4/2)|a|. Applying this and the other positive
embedding to(3) gives the complete necessary bounds

    a^2<=4q/(13phi^2), a'^2<=q*phi/2.                 (4)

The first constant follows exactly from
`phi^8/4-2phi=(13/4)phi^2`; dropping c^2,kappa only enlarges the
candidate domain. Both positive projections have nonnegative square
at the other embedding.

## 4. Exhausting every possible earlier mirror

For A=2a<0, the complete lists allowed by(4) are

|q|complete negative A list|
|---|---|
|1|-phi^-1|
|2|-phi^-1|
|4|-1,-phi^-1,-phi^-2,-2phi^-1|

These are direct two-embedding coefficient strips. For explicit
verification write A=m+nphi. At q4, the widest strip gives
`|A|<11/8`, `|A'|<18/5`; hence n=-2,-1,0 or1. Substitution gives
the four listed elements (the nearest rejected m+nphi are phi-3
and3-2phi). The q1 and q2 exact bounds delete all but-phi^-1.
For example q2 excludes-1 physically and-phi^-2 at the conjugate
embedding, while q1 is narrower still. These are necessary strips,
so including an extra element would be harmless, but none is omitted.

For q1, a=-phi^-1/2 makes `b^2=phi^2/2`, not integral, contradicting
b in O. At q2 the same a gives, in the three respective profiles,

    b^2=(phi+3)/2, phi^2/2, phi/2.

Each is nonintegral. Thus no norm1 or2 earlier wall exists.

For q4, a in O forces A in2O and leaves only `a=-phi^-1`. Its
isolated and B4 profiles respectively require

    b^2=2phi^2 or2phi,

which cannot be a square of an element of2O. The prime2 is inert,
and both right sides have dyadic valuation1, instead of at least2.
Thus every possible earlier mirror has been excluded.

At a0, the same six-profile table gives only g(q1) and Q(q2) as
primitive positive candidates. Isolated q2 would need b^2=2;
the F4 row b^2=1/2; the B4 row b^2=2. The isolated q4 root2g is
imprimitive. None is missing by normalization. The two surviving
vectors meet simultaneously at P1, but

    g=Q+beta6.

Both g and Q are integral reflectors with the halfspaces oriented
by their strictly negative products at the old initial vertex.
The equality g0 on the full chamber requires BOTH independent
equalities Q0 and beta6_0. Hence g is a nonsimple positive root,
whereas Q is the unique new SIMPLE wall. The positive definite full
corner at P1 is therefore ACTUAL with complete germ

    H4(E')+B3(a2,a3,s)+B2(Q,beta6).                  (5)

## 5. A full integral lattice isometry, with both actual germs

There is also an old actual full corner

    H4(E)+B3(b0,a3,s)+B2(beta5,beta6),                (6)

since all these walls are already actual and their whole Gram is
positive definite. Its future normal is

    P_old=phi^2*p0-2omega3, P_old^2=-2phi.

Define D by the following TEN independent vector images:

    (u,beta1,beta2,beta3) -> (N,beta4,beta3,beta2),
    b0 -> a2, a3 -> a3, s -> s,
    e5 -> -g, e6 -> e6, P_old -> P1.                 (7)

The old H4, old B3, the two unit vectors e5,e6 and P_old are
orthogonal blocks. Their image blocks have EXACTLY the same full
Gram. Thus(7) defines a K-isometry, not a similitude or a normalized
diagram match. In particular beta5=-e5+e6 maps to Q, and beta6
is fixed. It maps all positive roots in the FULL actual germ(6)
to those in the FULL actual germ(5), with the future normal mapped
to the future normal. The B2 simple-root signs are fixed explicitly.

Here is its entire10x10 matrix in an explicit INTEGRAL basis of the
original lattice R5, not in a merely similar lattice. Use

    B=(V0,e1,e2,e3,e4,e5,e6,a3,s,wP).

The first seven vectors are the unimodular M6 basis; the last three
form the exact P3 basis. In original(time;I5;D4) coordinates its
nontrivial columns are

    V0=(phi^3;0^5;-phi^2,-phi,0,-phi^3),
    e1=(-phi;0^5;phi,0,0,phi), e2,...,e6=I5 units,
    a3=(0;0^5;0,0,2,0), s=(0;0^5;-1,-1,-1,1),
    wP=(-phi^3;0^5;phi^2,phi^2,-1,phi^3).

Thus this displayed basis also specifies an explicit original-coordinate
matrix by `D_original=B*[D]_B*B^-1`. The actual integral-basis matrix
is the following; p abbreviates phi only inside the array:

```text
 60p+37 -22p-13 -22p-13 -22p-13 -22p-13 -6p-4 0 0 0 -26p-16
 22p+14  -8p-5   -8p-5   -8p-5   -8p-5  -2p-1 0 0 0 -10p-6
 19p+12  -7p-4   -7p-4   -7p-4   -7p-5  -2p-1 0 0 0  -8p-5
 19p+12  -7p-4   -7p-4   -7p-5   -7p-4  -2p-1 0 0 0  -8p-5
 19p+12  -7p-4   -7p-5   -7p-4   -7p-4  -2p-1 0 0 0  -8p-5
 19p+12  -7p-5   -7p-4   -7p-4   -7p-4  -2p-1 0 0 0  -8p-5
   0       0       0       0       0       0   1 0 0    0
 28p+16 -10p-6  -10p-6  -10p-6  -10p-6   -4p  0 1 0 -12p-8
 14p+8   -5p-3   -5p-3   -5p-3   -5p-3   -2p  0 0 1  -6p-4
 42p+24 -15p-9  -15p-9  -15p-9  -15p-9   -6p  0 0 0 -18p-11
```

For independent derivation rather than trusting the array, the only
dyadic glue identity in the old complement is

    p0=(phi^2/2)P_old-(phi/2)(3b0+2a3+s).

Its image under(7) is integral: modulo2,
`P1=phi^2(a2+s)` and phi^4=phi, so the numerator cancels. More
explicitly

    Dp0=(26phi+17;2phi^5,(phi^6)^4,0;
           4phi^4,2phi^4,6phi^4),
    DV0=(60phi+37;22phi+14,(19phi+12)^4,0;
           28phi+16,14phi+8,42phi+24).

With `F=phi^-1*N-phi^-2*DV0`, the other images are

    De1=F, De2=F+beta4, De3=F+beta4+beta3,
    De4=F+beta4+beta3+beta2, De5=-g, De6=e6,
    Da3=a3, Ds=s, DwP=a2-Dp0.

These verify the array and EVERY original integral basis image.
Gram preservation gives `(det D)^2=1`, hence det D is+1 or-1.
Its integral adjugate proves that D^-1 is integral as well. Thus
D is a full lattice automorphism, not a one-sided embedding.

## 6. Chamber matching and the obstruction, conditional on review

The integral isometry D preserves the FULL root arrangement. By the
actual future vertex and positive full-germ correspondence(5)--(7),
it takes an inward local chamber germ to an inward germ of the SAME
global chamber C. Connectedness of an arrangement chamber therefore
gives D(C)=C. Matching H4 root spans alone would not justify this step;
the complete preceding actual edge and B2 sign choice are essential.

The full displayed matrix has

    trace D=10+11phi>10.

A finite-order real10x10 isometry has all eigenvalues on the unit
circle and absolute trace at most10. Thus D has infinite order.
If C were compact, its locally finite arrangement would give finitely
many facets and a spanning finite vertex set. Its discrete stabilizer
would then be finite (equivalently it acts faithfully on that finite
spanning set). This contradicts the infinite-order D stabilizing C.

For the finite-volume conclusion use the standard arithmetic input
explicitly: the integral orthogonal group of this admissible form is
an arithmetic uniform lattice, because the other real embedding is
positive definite. If its reflection subgroup were cofinite it would
have finite index, and hence compact quotient. The noncompact full
reflection chamber therefore has infinite volume. This excludes the
literal R5 integral arrangement as a compact H9 construction host.

This is a COMPLETE PAPER CANDIDATE pending independent review and
fixed exact checks of the finite strips, full matrix/Gram/inverse,
future signs and trace. It is not yet a central-status promotion,
uniform compact-Coxeter bound, commensurability-class exclusion, or
new example. No successor root search is required if it passes.
