# Norm239: an original-integral upper on the second H2 edge

21 October2026. **New manual supporting candidate, pending root and
outside-author review.** The starting transverse germ below is itself a
separately frozen candidate. This note proves a finite mirror upper and
an exhaustive finite first-simple-wall contract conditional on that germ.
It does not identify the next actual wall, complete the chamber, or prove
a chamber period. No worker, import, enumeration process, test, Lean or
TeX job ran; no frozen paper, guide, worker or earlier source was edited.

## 1. Literal module and the actual starting boundary

Work over O=Z[phi], phi²=phi+1, in the ORIGINAL module

    L=E8(O) perpendicular-sum [[-2phi,1],[1,8]],
    D=1+16phi, b=phi D=16+17phi, N(D)=-239.

The native E8 Gram has diagonal2 and products-1 on
02,13,23,34,45,56,67. Its actual old walls are -e_j. Use

    q=(21+33phi)X-(17+25phi)Y, q²=-2D,
    m=(3+7phi)X-(3+5phi)Y, m²=2, q.m=0,
    H=Omega7=(2,3,4,6,5,4,3,2), H²=2,
    J=Omega6=(4,6,8,12,10,8,6,3), J²=6, H.J=3.

Here H.e_j is1 at7 and0 elsewhere. Put A=phi²,
E=7-4phi>0 and F=3-phi=A E. The preceding transverse source gives

    R=phi(q-m)/2+phi³H, R²=2, R.m=-phi,
    Q=2Aq+DH, Q²=-2DE,
    Q perpendicular to E7,m,R.

Its claimed actual full germ is E7+I2(5), with E7 on0..6. The literal
original binary coordinates of R are (13+22phi,-10-17phi).

| Direct input | SHA-256 | Interface |
|---|---|---|
| results/compact_h9_prime239_first_e7_golden_transverse_edge_candidate_20261021.md | 2a207fb82414203e34625d323238abc4c297541b6e34f23c4d28b16b7c3efa15 | Actual Q, R, full E7+I2(5) germ, original parity and all finite one-node profiles; pending separate review |
| results/compact_h9_prime239_full_base_and_first_binary_edge_candidate_20261021.md | 879a3feba630072a4b127e68a3729c4da3fca4f3c44630d9aae623ddd899de02 | Actual original base, q,m, reflection covectors and the original dyadic module |
| results/golden_odd_glue_full_base_and_discriminant_reduction_20261021.md | b107f0de4f76da942100e0c82867cd6d4eeee40d833304fee67cebd84d3e41ce | Full integral arrangement, primitive norm alphabet2/2b, golden ring and local finiteness |

Release OLD m and retain R and E7. With

    W=m+(phi/2)R, W²=F/2, Q.W=0,
    P(z)=Q-zW, 0<z<zmax=2sqrt(D)/phi,                 (1)

we have P(z).m=-zF/2<0 and P(z)²=-2DE+z²F/2.
Thus (1) is the actual outgoing edge germ. No claim that its full
positive interval remains in the chamber is used to construct an upper.

## 2. A diagram swap is not an original-lattice symmetry

The bisector m-R has square2(2+phi), whose halfnorm has norm5.
Its product with old e7 is -phi³, a unit. The reflection coefficient
on e7 is therefore phi³/(2+phi), which is not integral. More explicitly,
the E8 component of its correction is -phi⁶H/(2+phi); H is primitive
and has a coordinate3 coprime to the prime2+phi. This reflection fails
original-module integrality. The I2(5) diagram involution is not thereby
available in O(L).

Positive finite H2 Weyl elements also do not establish a chamber
stabilizer. For example the longest element sends m to -R and R to -m;
it sends q to the branch with positive product with m. It does not send
the actual incoming edge to (1). These sign and module issues must not
be bypassed by a similarity or by treating q,m as a unimodular basis.

## 3. A third-power conjugated boost is original-integral

Let r0=X+2phi Y be the old norm2b base wall. Its reflection and that
of m are original-integral, future-preserving isometries. Their product
B=s_r0 s_m fixes E8 and has original binary matrix

    B=[[336+545phi,436+704phi],
       [176+285phi,227+369phi]], detB=1,
    t=trB=563+914phi>2.

Writing ell=109+176phi, its q,m matrix is

    Bq=(tq-D ell m)/2, Bm=(-ell q+tm)/2,
    t²-D ell²=4.                                      (2)

These identities follow directly from the two displayed reflection
covectors in the binary-edge input. No new arithmetic-group criterion
or computational output is assumed.

Set N=H+R. Then N²=8A and the spacelike reflection

    S(v)=v-N(v.N)/(4A)

is future-preserving, fixes E7, and satisfies S(H)=-R, S(R)=-H.
Since A is a unit and N is original-integral, S has entries in(1/4)O.
It need not preserve L. Define V=SBS; V fixes E7 and R.

The following ORIGINAL-coordinate identities improve the generic
denominator16 estimate to V in(1/4)End_O(L). Put r=Vm and
K=Ve7-e7. Their E8 components in this table multiply H:

| Vector | X coefficient | Y coefficient | H coefficient |
|---|---|---|---|
| r | -2023-13099phi/4 | (6241+10097phi)/4 | -218-352phi |
| K | -44349/2-35879phi | 17095+55321phi/2 | -2389-3864phi |

For direct verification, in q,m coordinates r has coefficients
(-369-595phi)/8 and(777+1267phi)/8. Also, if
s=21+33phi, u=17+25phi, a=3+7phi, v=3+5phi, then

    Ve_j=e_j for j0..6, Ve7=e7+K,
    VX=X-(7+10phi)(r-m)-2v A K,
    VY=Y-(9+13phi)(r-m)-2a A K.                        (3)

Equations(2), S(H)=-R and the exact inverse
X=(vq-um)/2, Y=(aq-sm)/2 prove the table and(3) by expansion.

Let Pi be orthogonal projection onto the positive fixed subspace
E7+R. It has half-integral ORIGINAL entries: PiX=(X.R)R/2 and
PiY=(Y.R)R/2, while

    Pi e_j=e_j (j0..6),
    Pi e7=e7-H/2+(phi³/2)R.

On its complementary Lorentz plane V has trace t and determinant1.
On the fixed subspace it is the identity. Cayley--Hamilton therefore
gives the FULL ten-coordinate identity

    V³=(t²-1)V-tI+(2+t-t²)Pi.                         (4)

The Pi term is essential; dropping it would be a false full-space
identity. Since t=3+2phi modulo4, t²=1 modulo4, whereas
2+t-t²=0 modulo2. Equations(3)--(4) show V³ is original-integral.
Also V^-1=tI-V+(2-t)Pi has entries in(1/4)O, so the same argument
shows V^-3 is integral. Thus

    T=V³=S B³ S is in O^+(L), fixes R and E7.           (5)

This is an original-module isometry, not a similarity. The smaller
orbit vectors Vm and V²m are not integral: Vm has the fractional
binary coordinates above, and the H coefficient of V²m is phi/2
moduloO. This only concerns this particular orbit, not all next mirrors.

## 4. Its norm2 mirror is a genuine finite upper, not a simple-wall claim

In the active plane spanned by Q,W, T is a nontrivial future boost.
Its trace there is t³-3t>2; in particular it has no fixed future point.
The signs can be verified without choosing a numerical branch:

    Q.r=D(92+151phi)/2>0,
    W.r=(775+1265phi)/4>0,
    Q.(Tm)=(t²-1)Q.r>0,
    Q.(T^-1 m)=-(t²-1)Q.r<0.                         (6)

Thus TQ is future, perpendicular to R and E7, and has negative
product with m. It lies on the positive ray(1). It is also perpendicular
to the original-integral norm2 reflector Tm. Consequently that reflector
is a finite mirror upper on the actual outgoing edge, independently of
whether it is simple or whether any earlier facet intervenes.

For an entirely explicit upper parameter, put C3=t²-1. Equations(6)
and W.Pi m=0 give

    zup = [2D C3(92+151phi)] /
          [C3(775+1265phi)-2tF],
    0<zup<2sqrt(D)/phi.                               (7)

Equivalently the hyperbolic distance of this upper from Q obeys
coshL=(t³-3t)/2. Finiteness follows already from the future vector
TQ, not from a decimal comparison at the ideal endpoint.

The upper is deliberately NOT called a facet or a full actual germ.
The orientation of Tm that is negative at Q has positive product with
R, so it cannot by itself be the next simple facet. Formula(5) also
does not prove T(C)=C: no actual full germ at TQ has been supplied.
The size of coshL is about4 billion. Even this improved third-power
upper gives a very large finite coefficient box, not an economical
enumeration proposal or authorization to launch one.

## 5. Exhaustive finite contract for the next first simple wall

The compact segment from Q to the upper meets finitely many arrangement
mirrors. Starting from the actual germ, a first new simple facet exists
on it. Finite spherical classification retains E7+R and gives exactly
the same FIVE profile types, freshly expressed in the old original frame:

    r_new=(T0 q-K0 m)/2+hH+epsilon J,
    norm2:  (epsilon,kappa)=(0,0),(0,1),(0,phi),(1,0),
    norm2b: (epsilon,kappa)=(0,0),
    r_new.R=-kappa.                                   (8)

Indeed an E7 completion is only E8 at native6; that row is orthogonal
to R. A connected finite rank9 diagram cannot contain E7. For norm2
the exact strip1<=kappa<2, |kappa'|<2 in O gives only1 orphi.
For norm2b every nonzero finite product with a retained norm2 root
would give a totally positive angular invariant b eta² of norm>=239,
although both embeddings are<4; hence all such products vanish.

The parameters T0,K0,h belong to O. Original binary integrality is
EXACTLY T0=K0 modulo2O, not an invented unimodular q,m condition.
The contact and norm equations are

    K0=D T0-2A(2h+3epsilon)-2kappa/phi,
    K0²-D T0²+4h²+12epsilon h+12epsilon²=4d,           (9)

where d=1 or b. The first equation automatically implies the parity.
For norm2b require h=bj, j in O. Then epsilon=kappa=0 and
K0=D(T0-4phi³j); its ORIGINAL binary reflection covectors divided
by b are

    [-v T0+u(K0/D), -a T0+s(K0/D)]/(2phi),

which are integral because u-v and s-a are even. Its E8 covectors
are integral after dividing by b as well. Norm2 reflections are
automatically integral on the original even module. These norm classes
are primitive: neither the prime above2 nor that above b can divide
the vector twice in its norm. No norm79 lattice identity is used.

For the orientation supporting the chamber at Q, require

    Q.r_new=D(2h+3epsilon-2A T0)<0,
    K0+phi kappa/2>0,
    z(r_new)=D(2A T0-2h-3epsilon)/(K0+phi kappa/2)
              in(0,zup].                             (10)

Acuteness to the actual old -e7 and m additionally gives h>=0 and
K0>=0. Every known old facet may be checked explicitly; omitting such
extra necessary filters only enlarges the finite domain below.

Here is a complete finite bound, including all coefficient tails. Let
U=Q/sqrt(2DE), C=(t³-3t)/2 and S0=sqrt(C²-1). If a norm2d mirror
meets the segment of lengthL=arcoshC, Lorentz decomposition and ordinary
Cauchy--Schwarz give, for any fixed vector g,

    |r_new.U|<=sqrt(2d)S0,
    |r_new.g|<=sqrt(2d)[S0|g.U|+C sqrt(g²+(g.U)²)].   (11)

For the first inequality write the collision unit vector as
cosh(l)U+sinh(l)V, 0<=l<=L, and use
r_new²>=(r_new.U)²/sinh(l)². The second follows by decomposing g
into its U and U-perpendicular parts. At the other embedding the
form is positive definite, so |r_new'.g'|<=sqrt(2d'g'²).

Taking g=-q/D and g=e7 gives the explicit physical/conjugate bounds

    |T0|<=sqrt(2d)[S0*4A/sqrt(2DE)+C sqrt(2/E)],
    |T0'|<=2sqrt(d'/(-D')),
    |h|<=sqrt(2d)[S0*D/sqrt(2DE)+C sqrt(2+D²/(2DE))],
    |h'|<=2sqrt(d').                                  (12)

For z=a0+a1phi with |z|<=B, |z'|<=B', these imply

    |a1|<=(B+B')/sqrt5,
    |a0|<=((phi-1)B+phi B')/sqrt5.                    (13)

Outward integer ceilings of(12)--(13) give a finite FOUR-integer box
for each profile. It can be streamed in constant storage. For each
pair T0,h compute K0 by(9), then check the exact norm, original
coordinates/covectors, primitivity, the signs and interval(10), and
any known-wall acuteness filters. Exact Q(phi) signs compare collisions
by cross multiplication, retaining ALL equal minima and both labels
where applicable. There is no row-order winner rule.

Completeness is conditional only on the stated actual starting germ,
norm alphabet and finite-germ classification: every first simple wall
is in(8)--(13). Every admitted vector is itself an integral reflector,
so none can cut strictly before the actual first collision while the
segment stays inside the full chamber. Thus the least admitted positive
collision equals that first collision. All minimum ties and the actual
positive root cone must then be checked to identify the full simple
germ; a bare minimum label does not establish it. This is a mathematical
contract, not an implemented or executed checker and not a report of
its result. A substantially sharper upper/domain would be needed for
a sensible guarded census. Actual next-wall identification, a chamber
period and whole-host reflectivity remain open.
