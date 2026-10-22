# Norm239 after E7+2A1: an original integral upper and a594-row domain

22 October 2026. **Complete prospective manual candidate, UNREVIEWED.**
Root supplied the skew-operator/half-integral boost idea. This source
independently derives its full original-basis integrality, orientation,
finite-simple profiles and finite domain. The author contributed the
previous domain/worker and fourth-edge review. No worker is prepared,
imported or executed here; no enumeration, symbolic process, test,
Lean, TeX or Git operation occurs. No frozen source or paper is changed.

**Conditional proposed theorem.** Starting from the actual E7+2A1 germ
below, release S and retain E7 and the norm2b wall T. There is a proved
finite original-integral mirror upper on that outgoing edge. Every
possible first SIMPLE facet before it folds into three phases and five
profiles, with at most594 scalar rows and two field-square signs per
row. This supplies no new first wall, completed chamber, period,
finite-volume example or nonreflectivity claim.

## 1. Literal module and exact starting boundary

Put O=Z[phi], phi²=phi+1, D=1+16phi, b=phi D=16+17phi,
Norm(b)=239. The ORIGINAL module is

    L=E8(O) perpendicular-sum [[-2phi,1],[1,8]],
    H=Omega7=(2,3,4,6,5,4,3,2), H²=2.

The native E8 bonds are02,13,23,34,45,56,67. E7 is spanned by
e0,...,e6 and is orthogonal to X,Y,H. In particular the active
projection of the ORIGINAL e7 is H/2, not H. Let

    S=(55+90phi)X-(43+69phi)Y+(6+10phi)H,
    T=aX+yY+hH,
    a=12529+20274phi, y=-9655-15623phi, h=1446+2343phi,
    P=[(21786+873phi)X-(9285+5310phi)Y
                                    +(-5475+5040phi)H]/1441.

The starting input states that -e0,...,-e6,S,T are the full ACTUAL
simple germ, S²=2, T²=2b, S.T=0, and P is its future timelike
corner ray. Release S. The outgoing line is P-zS for z>0.
The initial segment lies in the SAME full integral-arrangement chamber.

| Exact input | SHA-256 | Interface |
| --- | --- | --- |
| `results/compact_h9_prime239_fourth_edge_and_e7_two_a1_germ_candidate_20261022.md` | `258a6e99703e4d65ed7d71020fc642edd7c36a644e59acf63bde94ab5e9d49ef` | Actual P,S,T and full starting germ, conditional on its reviewed execution boundary |
| `results/prime239_released_r_runtime_manifest_20261022.json` | `e3e6a54e34ddab2d4213bd355c1f21c6534b275e6496b54d8b0326d9ad3bd85f` | All twenty actual preceding units plus separate combiner; NOT an execution of this new domain |
| `results/review_prime239_fourth_edge_e7_two_a1_composition_and_runtime_replay_reviewer_20261022.md` | `4aad6a513df98f142385584cc7091eae8e5a7c8bf04f616009df3ac4a7fc6dbb` | Full preceding composition/runtime review, including original reflector T |
| `results/compact_h9_prime239_full_base_and_first_binary_edge_candidate_20261021.md` | `879a3feba630072a4b127e68a3729c4da3fca4f3c44630d9aae623ddd899de02` | Original admissible module and full arrangement |
| `results/golden_odd_glue_full_base_and_discriminant_reduction_20261021.md` | `b107f0de4f76da942100e0c82867cd6d4eeee40d833304fee67cebd84d3e41ce` | Complete primitive norm2/2b alphabet and original-module normalization |

Discreteness/local finiteness, finite spherical classification and
finite local Coxeter-chamber recognition remain the inherited geometric
inputs. O has its two real embeddings; the conjugate form is positive
definite. The ideal(b) is prime since its rational norm is239. Chamber
compactness or finite covolume of its reflection group is NOT assumed.

## 2. A skew operator on the active3-space

Use G=diag(binary,2) in basis X,Y,H, so det G=-2D. The usual cross
matrix C(T) has rows(0,-h,y),(h,0,-a),(-y,a,0). Set A=G^-1 C(T).
Then AT=0, A is G-skew, and

    A=[ h/D,       8h/D,       (-8y-a)/D;
        2phi h/D, -h/D,        (y-2phi a)/D;
        -y/2,      a/2,         0 ].

Let P_T=T(GT)^t/(2b). Direct multiplication, or the three-dimensional
cross identity, gives

    A²=-(T²/det G)*(I-P_T)=phi*(I-P_T).            (1)

The relevant ORIGINAL integral coefficients are

    h/D=54+87phi,
    (-8y-a)/D=2407+3894phi,
    (y-2phi a)/D=-1867-3021phi,
    GT/b=(-1154-1867phi,-1487-2407phi,0,...,0,33+54phi).

Each follows by multiplying by D or b and using phi²=phi+1.
The last vector is the full ten-coordinate covector, not only its
compressed active entries.

Extend A by zero on E7. Let P_active be the orthogonal projection onto
span(X,Y,H), and write E=P_active-P_T. Then E is the orthogonal
projection onto the Lorentz plane perpendicular to E7 and T. On the
FULL space, A²=phi E, AE=EA=A and A P_T=P_T A=0.

The half-integrality statements require the old e7 correction:

* 2P_active sends X,Y to2X,2Y, sends e0..e6 to0 and sends e7 to H.
* 2P_T=T tensor(GT/b), hence sends EVERY original basis vector into L.
* 2AX and2AY are integral, since the only half entries of the displayed
  A are its last-row entries-y/2,a/2. Also A e_i=0 for i<=6 and
  2A e7=A H=(2407+3894phi)X-(1867+3021phi)Y is integral.

Thus2E L is contained in L and2A L is contained in L. This proves
integrality on the ORIGINAL E8 basis; active coefficient integrality
alone would not have sufficed.

## 3. A short integral finite mirror upper, with the correct direction

Define on the full space

    V=I+(phi-1)E+B, B=A or-A.                     (2)

Choose the sign so that BP=-alpha S with alpha>0. This is possible
because A preserves the Lorentz plane, AP is perpendicular to P and
nonzero, and S spans that perpendicular line. Explicitly
alpha=|AP.S|/2 belongs to K=Q(phi), and skew-adjointness with(1) gives

    alpha²=-phi P²/2>0,  BS=-(phi/alpha)P.

Formula(2) is the IDENTITY on E7 and T. On the Lorentz plane it is
phi I+B; phi²-phi=1 proves that its inverse changes the sign of B.
Since E is self-adjoint and B is skew, V is an isometry. Also
P.VP=phi P²<0, so it preserves the future cone. Its translation
length ell on this line satisfies cosh ell=phi.

Although V need not be original-integral, both signs are half-integral.
Its square and inverse square are original-integral by the FULL maps
just established:

    V^2=I+2phi E+2phi B,
    V^-2=I+2phi E-2phi B.                          (3)

They are mutually inverse, hence preserve L onto itself. More generally
V^j=I+(lambda_j-1)E+mu_j B, lambda_j,mu_j in O: multiplication gives
lambda_(j+1)=phi(lambda_j+mu_j), mu_(j+1)=lambda_j+phi mu_j;
negative j reverses mu. Thus all required signed powers are half-integral.

In the outgoing direction the exact future finite point is

    V²P=(1+2phi)P-2phi alpha S,
    z_upper=2phi alpha/(1+2phi)>0.

It lies strictly before the ideal end: writing z_ideal²=-P²/2,

    z_upper²/z_ideal²=4phi³/(1+2phi)²<1,
    (1+2phi)²-4phi³=1.

The ORIGINAL integral norm2 mirror r_upper=-V²S has incidence at V²P.
Indeed r_upper=(2phi²/alpha)P-(1+2phi)S, so
P.r_upper<0 and -S.r_upper=2(1+2phi)>0. It supplies a genuine
finite arrangement-mirror upper. Neither its simplicity nor chamber
stabilization by V² is assumed. Local finiteness and the actual starting
germ therefore give a first new simple facet at or before this upper.

## 4. COMPLETE finite-simple profiles with a norm2b retained wall

Normalize a possible first simple wall to primitive square2d, d1 or b,
using the pinned alphabet. Its finite endpoint retains E7 and T.
Write epsilon for its possible E7 contact and kappa=-r.T>=0.

First suppose d1. Since GT/b is integral, r.T=b k with k in O.
A nonzero finite contact would obey0<b k²<4 in BOTH embeddings,
by the positive definite two-root principals. Its norm would then be
239 Norm(k)²<16, impossible. Thus kappa=0. Finite classification of
E7 plus one node leaves exactly no E7 contact or the E8 extension at
native node6, with unit product-1. There is no other rank8 finite
completion containing E7. Thus the norm2 profiles are

    (epsilon,kappa)=(0,0),(1,0).

For d=b, original reflection integrality and primitivity imply that
every covector of r is divisible by b. At the prime(b), a coordinate
of r is nonzero, since r²=2b has valuation1; its reflection matrix
I-r(Gr)^t/b then forces this divisibility. Any contact with a norm2
E7 root would again have239 Norm(k)²<16, so epsilon=0.

Now kappa=b k for k in O, with0<=k<2 and |k'|<2. Write k=a+c phi.
The difference of embeddings gives-2<c sqrt5<4, hence c0 or1.
For c0, a0 or1. For c1, a=-1 or0. The resulting candidates are
0,1,phi-1,phi; the actual native gap k=2cos(pi/m)>=1 for a
nonzero finite bond removes phi-1. Consequently the norm2b profiles are

    (epsilon,kappa)=(0,0),(0,b),(0,b phi).

This proves the five-profile list without transporting the previous
norm2-retained list or silently replacing T by a norm2 root.

## 5. A594-row scalar domain after nearest-phase folding

A first simple collision has distance in[0,2ell] along the outgoing
edge. Choose a nearest phase j in{0,1,2} and seed=V^-j r. Its signed
collision has distance at most ell/2 from P; it need not be positive
at P or acute to released walls. V fixes E7 and T, so its profile is
unchanged and its WHOLE original coordinates lie in(1/2)L.

Let g=seed.S and Z0=Omega6-(3/2)H. Here Z0²=3/2 and it is
orthogonal to P,S,T, with the required native node6 unit contact.
The unique decomposition is

    seed=tau P+(g/2)S-(kappa/(2b))T+epsilon Z0,
    tau²=[g²/2+kappa²/(2b)+(3/2)epsilon²-2d]/(-P²). (4)

Both exact field-square signs and zero must be retained. Coefficients
tau or g/2 are NOT individually presumed half-integral; test the
WHOLE reconstructed original seed. Since

    cosh(ell/2)=phi/sqrt2<3/2,

the usual Lorentz-space Cauchy bound at the signed collision gives
|g|<3sqrt(d). Conjugate positive definiteness gives |g'|<=2sqrt(d').

For norm2, g belongs to(1/2)O. Write2g=a+c phi. From |g|<3 and
|g'|<=2, we obtain |c|sqrt5<10, hence c=-4,...,4. The exact
conjugate interval |a-c(phi-1)|<=4 contains at most nine integers a.
There are at most81 scalar rows for each norm2 profile/phase.

For norm2b we need the STRONGER divisibility g/b in(1/2)O. Put
f=2seed in L, f²=8b. A coordinate of f is nonzero modulo(b).
Conjugating its original root reflection by the signed half-integral
phase powers gives a matrix in(1/4)O. But that matrix equals

    I-f(Gf)^t/(4b).

Multiplication by4 and the nonzero coordinate force b to divide every
covector of f. In particular b divides f.S=2g. Thus jg=g/b is half
integral, with

    |jg|<3/sqrt(b)<1/2, |jg'|<=2/sqrt(b')<7/8.

Here b>36; phi<13/8 gives b'=33-17phi>43/8>256/49, proving the
second strict bound. Writing2jg=a+c phi now gives
|c|sqrt5<11/4, hence c=-1,0,1. The conjugate interval
|a-c(phi-1)|<7/4 contains at most four integers a. There are at
most12 scalar rows per norm2b profile/phase. The whole finite domain is

    3*(2*81+3*12)=594 scalar rows.                  (5)

This is a proved coefficient overlist, not a measured runtime or an
observed first wall. There are at most two exact roots of(4) per row.

## 6. Exact continuation/firstness obligations; nothing executed

An eventual implementation must enforce the physical/conjugate scalar
bounds, solve(4) exactly, reconstruct the WHOLE seed and test its original
half membership. For nonzero g the centered collision is

    z_seed=(P.seed)/g,
    z_seed²<=(-P²/2)*(phi-1)/(phi+1).               (6)

Both signs and equality are retained. When g0 the mirror gives no
isolated collision with this line and is skipped. Reconstruct
r=V^j seed; only now impose original INTEGRAL coordinates and every
reflection-matrix entry, exact square2d, original outgoing orientation,
positive collision before the finite upper and acuteness to every
known actual wall -e0,...,-e7,r0,m,R,S,T. E7 and T stay retained.

Keep exact minima and ALL ties, rejecting overflow rather than
truncating. Complete all three phases/five profiles under identical
current pins; absent actual units earn no coverage. Only after their
actual terminal ledgers and a complete exact combination exist may
firstness be inferred: every actual first simple wall occurs in the
domain, while every admitted root is an ORIGINAL arrangement mirror
and cannot cut the chamber interior first. A positive definite nine-root
Gram, full-coordinate tie cone and finite local chamber interpretation
are then needed for a new full-germ claim.

No worker or output exists for this domain in this source. A future
implementation and its complete mathematical contract require root and
other-author review BEFORE any guarded launch. Root retains execution
ownership, unchanged40GiB/6GiB admission and8GiB/3GiB abort floors,
59s/1GiB/one-thread bounds and the shared serial lock. Formula(5) alone
does not certify those resource bounds will be met.

In particular neither the half-integral V nor the original-integral V²
has been proved to preserve a full actual chamber germ. No period,
host exclusion, finite volume, compact H9 example, form-class exclusion,
global dimension change or new execution/formalization is asserted.
