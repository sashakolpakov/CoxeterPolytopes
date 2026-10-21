# Golden k4: a bounded transverse edge and a candidate chamber period

20 October 2026. Complete PAPER candidate, pending independent full
review. No scientific execution, factorization search, Pell search,
worker, Lean or TeX job. All accepted sources and freezes are unchanged.
Firstness and the resulting literal-host exclusion are NOT promoted
by this awaiting-review edition.

The proper-divisor residue obstruction was supplied by root and is
rederived below with its lattice and normalization prerequisites. The
explicit upper, finite hand tables and period were derived here.

## 1. Literal lattice, actual starting edge, and a reusable upper

Use O=Z[phi], phi^2=phi+1, with native E8 basis e0,...,e7 and edges
02,13,23,34,45,56,67. Its Gram has diagonal2 and edge products-1.
Let omega_j be its raw duals. In particular

    omega7=(2,3,4,6,5,4,3,2), omega7^2=2,
    omega6=(4,6,8,12,10,8,6,3), omega6^2=6,
    omega6.omega7=3.

The literal k4 lattice is E8 orthogonal-sum B, where

    u=phi^4=2+3phi, a=phi^5=3+5phi,
    P=phi^9=21+34phi, b=P+4=25+34phi,
    B=[[-2phi,a],[a,2]], D=phi*b=34+59phi.

Its full chamber C has the Bugaenko base r=(u,2;0), norm2b, and
the ACTUAL E8+A1(c) endpoint

    c=(0,-1;0), p1=(2,-a;0), p1^2=-2D, p1.c=0.

These are the accepted full-base and uniform first-binary-edge inputs,
pinned at the end. The old E8 outward normals are -e_j. We release
-e7 and retain E7 on -e0,...,-e6 together with c. The outgoing axis
is p1+tau*omega7, 0<tau<sqrt(D); no assertion that a metric face
is Coxeter intrinsically is needed.

Here is a useful ALL-k upper before specializing. Put P=phi^(2k+1),
u=phi^k, b=P+4 and p1=(2,-phi*u;0), for any k>=-1. The equation

    b*V^2-phi*X^2=1                                  (1)

has an explicit positive integral solution:

    P=1 mod2: V=(P+1)/2, X=u*(P+3)/2;
    P!=1 mod2: V=(P^2+3P+1)/2,
               X=u*(P^2+5P+5)/2.                    (2)

Indeed O/2=F4 and P is a unit. In the second case P^2+P+1=0
there, so both numerators are even. The two polynomial identities are

    (P+4)(P+1)^2-P(P+3)^2=4,
    (P+4)(P^2+3P+1)^2-P(P^2+5P+5)^2=4.

Thus m_b=X*p1+b*V*omega7 has norm2b, integral reflection covector
Gm_b/b, and is primitive: a common prime of X and bV contradicts(1),
while the first two coordinates have content X since a is a unit.
It is orthogonal to all eight retained walls. It cuts the outgoing
axis at tau=phi*X/V<sqrt(D). This is a uniform finite UPPER mirror,
not a uniform firstness or chamber-symmetry theorem.

For k4 the first line gives the particularly small upper

    Xb=75+121phi, Vb=11+17phi,
    b*Vb^2-phi*Xb^2=1,
    m_b=Xb*p1+b*Vb*omega7,
    tau_b=phi*Xb/Vb.                                 (3)

Local finiteness of the full arrangement gives a first actual simple
opposite at or before this finite axis point. We classify that actual
wall, not an assumed fixed seed.

## 2. Primitive classes and COMPLETE finite first-simple profiles

The even golden lattice and primitive reflection integrality give
norm2d with d totally positive and d|b. Here

    b=(3+phi)(6+7phi)=d11*d29, N(d11)=11,N(d29)=29.

Both factors are totally positive and have distinct prime norm. The
golden PID and totally-positive-unit-square facts therefore normalize
the primitive square to exactly2d for

    d=1,d11,d29,b.

Rescaling by a unit changes no mirror and preserves every integral
coordinate condition below. No unidentified unit norm class is lost.

If d is nonunit, every retained norm2 wall must be orthogonal to the
first mirror. Indeed its pairing with that mirror is d*gamma by
reflection integrality, hence the squared angular invariant is
d*gamma^2. Every nonzero finite intersection value in the accepted
golden angle list has rational norm supported on2,3,5. A prime over
11 or29 dividing d forbids such a value. This argument also covers
nonnative rational arrangement angles, although the actual first
simple wall already has native angles.

For an integral vector orthogonal to c, its binary coordinates obey
2y=-a*x. Since a is a unit, x=2X,y=-aX with X in O. Orthogonality
to the retained E7 makes its E8 part Y*omega7, Y in O. Integrality
of the reflecting covector against e7 gives d|Y. Consequently every
nonunit first wall has precisely

    m=X*p1+d*V*omega7,
    d*V^2-phi*(b/d)*X^2=1,        X,V in O.            (4)

For d11, reduction modulo d11 has phi=8 in F11 and phi*(b/d)=1.
Equation(4) forces X^2=-1 mod11, impossible. For d29 the residues
are phi=24, b/d=27, phi*(b/d)=10 in F29; thus X^2=26 mod29.
The nonzero squares modulo29 are
1,4,9,16,25,7,20,6,23,13,5,28,24,22, so this is impossible too.
These exclusions hold at all heights, not only before(3).

For d=1, finite classification of the ACTUAL rank9 endpoint is short.
The eight retained facets have diagram E7+A1(c). A new actual simple
wall can be orthogonal to both; join only c by a finite bond; or
extend E7 to E8 by simple3 at its unique length-three-arm tip6,
remaining orthogonal to c. It cannot connect both old components:
no connected rank9 finite diagram contains E7. The E7 extension is
unique by finite E8 classification. A norm2 wall joining the norm2 c
has raw product -eta with eta in O. The finite golden angle list
leaves only native eta=1 or phi, besides zero. Values2 and3 are not
squares in Q(phi); the two label10 values have norm5; phi^-1 gives
the nonnative angle2pi/5 rather than a simple chamber angle.

Thus ALL actual first-simple possibilities, with outward orientation
and all-wall acuteness against the released -e7, are

    U0: m=X*p1+Y*omega7, Y^2-DX^2=1;
    Ueta: m=(U/2)*p1-(eta/2)*c+Y*omega7,
          4Y^2-DU^2=Reta=4-eta^2,
          U in O, eta-aU in2O, eta=1 or phi;
    E8: m=X*p1+omega6+Y*omega7,
        Z=2Y+3, Z^2-4DX^2=1;
    B: m=X*p1+bV*omega7, bV^2-phiX^2=1.              (5)

In every row the axis time coefficient is positive, and Y or V is
positive when needed for the outgoing finite collision. In the E8
row Y>=0 by acuteness with -e7. The unit rows are automatically
integral subject to the displayed congruence. The Ueta congruence
is exactly integrality of its original binary second coordinate;
it is not a change to a merely similar lattice.

Only actual first SIMPLE profiles are being classified. A nonnative
or nonsimple mirror at a simultaneous endpoint cannot create an
earlier endpoint without an actual simple opposite there, covered by
this list. The final selected wall's simplicity is checked below.

## 3. Complete finite scalar strips before the upper

For U0 and E8, comparison of their collision with tau_b gives

    U0: X^2<=Xb^2/b, |X'|<=1/sqrt(-D');
    E8: X^2<=Xb^2/(4b), |X'|<=1/(2sqrt(-D')).

For Ueta it gives

    U^2<=Reta*Xb^2/b, U'^2<=Reta'/(-D').

For B it gives 0<V<=Vb and |V'|<=1/sqrt(b'). These follow by
writing the squared collision ratio tau^2/D as1 minus, respectively,
1/Y^2,1/Z^2,Reta/(4Y^2),1/(bV^2). No limiting-height assumption
or root-prefix filter occurs.

The elementary exact interval1.618<phi<1.6181 proves

    Xb/sqrt(b)<30.3, -D'=59phi-93>2.462,
    b'=59-34phi>3.9846, Vb<39.

For the first inequality,100[(30.3)^2*b-Xb^2]
=268625-157594phi>0 even using phi<1.7. Convenient conservative
strips, with all variables positive physically, are therefore

    U0: X<30.3, |X'|<.638;
    E8: X<15.15, |X'|<.319;
    U1: U<54, |U'|<9/8;
    Uphi: U<20, |U'|<5/4;
    B: V<39, |V'|<.502.                              (6)

Write the searched variable as A+Bphi. Since
B=(variable-variable')/sqrt5, these strips give respectively
0<=B<=13,6,24,9,17. The value B=0 gives no positive variable in
these strips except when barred by parity anyway. Using
1.618<phi<1.6181 in |A+(1-phi)B| gives exactly the candidate pairs
in the following tables. For U1, A is even and B odd; for Uphi,
both A and B are odd, since a=phi^5=phi^2 modulo2.
The tables use the CONSERVATIVE strips, so a few extra values past
the exact upper are harmlessly retained.

Each entry is the rational norm of the required field square. A
non-square integer norm excludes it. These are finite HAND tables;
no program or scientific calculation was run.

### E8: N(1+4D X^2)

| (A,B) | Norm |
| --- | ---: |
| (1,2) | 4189 |
| (2,3) | 19225 |
| (3,5) | 58589 |
| (4,6) | 15649 |

All are nonsquares, so no E8 extension occurs in the interval.

### U0: N(1+D X^2)

| (A,B) | Norm |
| --- | ---: |
| (0,1) | 20 |
| (1,1) | 569 |
| (1,2) | 2005 |
| (2,3) | 5764 |
| (2,4) | 4189 |
| (3,4) | 3641 |
| (3,5) | 15605 |
| (4,6) | 19225 |
| (4,7) | 22436 |
| (5,8) | 41369 |
| (5,9) | 10945 |
| (6,9) | 28900 |
| (6,10) | 58589 |
| (7,11) | 71641 |
| (7,12) | 51745 |
| (8,12) | 15649 |
| (8,13) | 108820 |

Only28900=170^2 needs more than its norm. Here X=3phi^4 and the
radicand is15130+24480phi. Its trace is54740. An integral square
root would have integer trace whose square is54740+/-340, namely
55080 or54400. Neither is an integer square. Thus U0 is absent.

### U1: N((D U^2+3)/4)

| (A,B) | Norm |
| --- | ---: |
| (0,1) | 44 |
| (2,3) | 1121 |
| (2,5) | 59 |
| (4,5) | 1136 |
| (4,7) | 5204 |
| (6,9) | 8649 |
| (6,11) | 6551 |
| (8,13) | 20444 |
| (10,15) | 16049 |
| (10,17) | 26951 |
| (12,19) | 41936 |
| (12,21) | 10944 |
| (14,21) | 8009 |
| (14,23) | 61259 |

Only8649=93^2 needs a trace test. Its radicand is3783+6120phi,
trace13686. The two possible squared traces of its integral square
root are13872 and13500; neither is a square. Thus U1 is absent.

### Uphi: N((D U^2+3-phi)/4)

| (A,B) | Norm |
| --- | ---: |
| (1,1) | 181 |
| (1,3) | 505 |
| (3,3) | 191 |
| (3,5) | 3581 |
| (5,7) | 5391 |
| (5,9) | 8791 |

All are nonsquares, so Uphi is absent.

### B: N((b V^2-1)/phi)

| (A,B) | Norm |
| --- | ---: |
| (1,1) | 229 |
| (1,2) | 1116 |
| (2,3) | 3439 |
| (2,4) | 639 |
| (3,5) | 9521 |
| (4,6) | 9931 |
| (4,7) | 10819 |
| (5,8) | 25444 |
| (6,9) | 7991 |
| (6,10) | 34259 |
| (7,11) | 41229 |
| (7,12) | 17236 |
| (8,13) | 67131 |
| (9,14) | 41564 |
| (9,15) | 62729 |
| (10,16) | 97951 |
| (11,17) | 3481 |

The sole square is3481=59^2. Its field square is exactly
(75+121phi)^2, and positive orientation gives X=Xb,V=Vb.

For reproducibility of these hand tables, if w=A+Bphi then

    Nw=A^2+AB-B^2,
    Tr(Dw^2)=127A^2+422AB+338B^2,
    N(1+Dw^2)=1+Tr(Dw^2)-319(Nw)^2,
    N(1+4Dw^2)=1+4Tr(Dw^2)-5104(Nw)^2,
    N((Dw^2+3)/4)
      =[381A^2+1266AB+1014B^2+9-319(Nw)^2]/16,
    N((Dw^2+3-phi)/4)
      =[465A^2+1520AB+1225B^2+5-319(Nw)^2]/16,
    N((bw^2-1)/phi)
      =84A^2+254AB+211B^2-1-319(Nw)^2.

If a field radicand in O had a square root in K, that root would
be in O by integral closedness, justifying both norm and integer-trace
tests. There are58 retained scalar candidates in total.

## 4. Candidate firstness and the FULL actual germs

The preceding complete finite first-simple inventory has no earlier
winner and exactly the upper(3) as its surviving simple possibility.
Hence, SUBJECT TO INDEPENDENT REVIEW OF THIS PAPER CENSUS, m_b is
the actual first wall and the new vertex is

    Pnew=(Vb/phi)*p1+Xb*omega7, Pnew^2=-2/phi.        (7)

It has whole actual germ E7+A1(c)+A1(m_b). A distinct simultaneous
non-simple mirror cannot change this conclusion without another
actual simple winner from the exhausted inventory. In particular
full-b mirrors are mutually ultraparallel, and a unit-norm mirror
cannot insert an extra finite root length into the full-b A1 factor.
This is not firstness inferred merely from the displayed upper.

The old Bugaenko norm2 wall g=(1,0;phi*omega7) is an ACTUAL ambient
wall by the full base calibration. The PD set of ACTUAL walls

    -e0,...,-e6, r, g

has a finite actual vertex by positive acute projection, with vector

    Pold=(1,0;omega7), Pold^2=-2/phi.

Its whole germ is E7+A1(r)+A1(g). This proves the old vertex from
actual walls rather than treating a metric fixture as the whole
ambient chamber by analogy.

## 5. Integral period conditional on the firstness gate

Define T to fix E7 pointwise and map the orthogonal frame

    (r,g,Pold) -> (m_b,c,Pnew).

Both frames have Gram diag(2b,2,-2/phi). Their determinants in the
three-space(time,Y,omega7) are2/phi, so T is an isometry of
determinant1. It is future because Pold and Pnew are future.

Here X0,Y0 denote the original two binary BASIS vectors, not the
scalars Xb,Vb. Put Z0=omega7. The formulas sufficient to recover
the full native matrix are

    T X0=phi*Vb*p1-phi*c+phi^2*Xb*omega7,
    T Y0=(m_b-u*T X0)/2,
    T omega7=phi*c-Vb*p1-phi*Xb*omega7,
    T e_j=e_j (j<=6),
    T e7=e7+(T omega7-omega7)/2.                    (8)

For k4, uVb=1 and Xb=phi^-1 modulo2. These two congruences make
both halves in(8) integral. Explicitly Tomega7-omega7 has binary
coordinates(-2Vb,aVb-phi) and E8 coefficient-(phiXb+1), all even.
For TY0, its first binary coordinate is integral and its second
numerator is -a(Xb-aVb+1); modulo2 this vanishes because
aVb=phi and Xb=phi+1. Its E8 numerator is bVb-u phi^2 Xb, also
zero modulo2. Thus T is integral on the ORIGINAL full L_k4 basis,
not just a sublattice or a K-similar orthogonal form. Determinant1
gives an integral inverse.

The trace on E7 is7. On the orthogonal three-space, direct expansion
of(8) gives

    tr(T|3)=phi*b*Vb/2-phi*Xb-a*(Xb+1)/2
            =151+245phi,
    tr(T)=158+245phi>10.                            (9)

Under the actual-firstness conclusion in Section4, T maps the full
oriented old vertex germ to the full oriented new vertex germ. An
integral isometry preserves the FULL reflection arrangement, so this
local chamber agreement proves T(C)=C. Equation(9) proves infinite
order, since a finite-order10-dimensional real matrix has trace
at most10.

The full reflection group W acts freely on chambers. Its intersection
with their setwise stabilizer is trivial. Distinct powers of T then
give infinitely many cosets in O^+(L_k4)/W, proving nonreflectivity.
With the explicit external Borel--Harish-Chandra arithmetic-lattice
input, and anisotropy from the positive nonphysical embedding, this
also gives infinite full-chamber volume. This conclusion remains a
CANDIDATE pending independent review, not a current status change.
It concerns only this literal k4 module, not its whole similarity or
commensurability class, and supplies no dimension ceiling.

## Bound inputs and remaining review gate

* Full-base theorem:
  `golden_e8_binary_family_full_base_facet_20261020.md`, SHA256
  `ea5d96e331ac91cb9efc457db2ddaa24bac3341e9e6781038150765895675afd`.
* Actual first binary edge:
  `golden_e8_binary_family_uniform_first_edge_20261020.md`, SHA256
  `71cfa11bcadf5fffae017d712d168eeec44eecf7d5b74f595bcebca7ad92af2c`.
* Native E8 and full BG calibration:
  `data/known_examples/bugaenko_compact_h8.json`, SHA256
  `c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c`.
* Golden norm/angle/ring and full-germ conventions are those of the
  pinned full-base theorem and its inherited audit774e87f7. Standard
  finite Coxeter classification and the actual positive-projection
  lemma are explicit geometric inputs, not checked by a program here.

The complete58-value hand census, both field-square exceptions,
finite-profile classification and integral period require root and
outside-author review before any firstness or host-exclusion promotion.
No checker is proposed or run by this paper-only draft.
