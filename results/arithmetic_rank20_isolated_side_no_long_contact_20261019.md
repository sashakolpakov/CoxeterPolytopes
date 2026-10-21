# Isolated opposite side: complete exclusion with no long-arm contact

19 October 2026. **Written local theorem candidate, awaiting independent
review.** The uniform ideal-divisibility arguments and their native-label
completion below were supplied by the coordinating root agent and checked
independently here. The finite and zero opposite-product reductions are
written out at the smaller ranks. No scientific process, enumeration,
optimizer, Lean or build ran. Earlier pinned sources are unchanged.

This is a restricted actual-quadrilateral exclusion. It does not yet
prove the two half refunds at any remaining D+D/E whole forest, or a
new dimension ceiling.

## 1. Actual geometry, complete contacts, and the arithmetic interface

Let P be a compact arithmetic simplest-type20-polytope over one of

    Q(sqrt3), Q(sqrt6), Q(sqrt2,sqrt3),
    K9,K11,K15,K16,K20,       Kj=Q(cos(2pi/j)).

At an actual vertex suppose the ENTIRE diagram is D_r+T_s, r+s=20,
r>=10, with T_s=D_s (4<=s<=10) or E6,E7,E8. Thus10<=r<=16.
Write a for the D_r fork, b,e for its two short leaves, and
a-l1-...-lL for its long arm. Put

    N=r-2 in {8,...,14},       L=N-1,       s0=r/(4N)<=5/16.

Consider an actual quad with cyclic sides a,b,c,d and complete white
support W=A1(e)+A_L+T_s. Assume:

* d is perpendicular to every old root except b;
* c is perpendicular to b,d, the entire spectator T_s, and every li;
* c.e=-gamma, where gamma=0 or cos(pi/ell), ell>=3;
* a.c=-x is initially unrestricted by opposition.

**Theorem.** This actual quadrilateral cannot occur.

The c restrictions are the precise no-long-contact, spectator-orthogonal
part of the isolated-d row. No contacts to an unspecified old factor
have been discarded: D_r+T_s is the entire old vertex. The roots
c,d are actual facets of the ambient compact polytope, not formal
Gram entries or a proposed completion.

Project off W. Set t=1-gamma^2 and X=x+gamma/2. The four-root Gram
is

    [ s0   -1/2   -X    0   ]
    [-1/2   1      0   -beta]
    [-X     0      t    0   ]
    [ 0    -beta   0    1   ],      beta=-b.d.

Here t>0 by the actual finite corner. The unchanged b,d normals
are opposite sides of the acute compact polygon, so beta>1:
positive-definite two-normal projection would give an actual common
corner, and the positive-null case gives a feasible ideal ray.
The22 actual normals in Lorentz dimension21 give determinant zero:

    (X^2/t-s0)*(beta^2-1)=1/4.

Put

    R=X^2/t,       z=4*beta^2,       delta=z-4=1/(R-s0)>0.          (1)

The old full vertex plus c and plus d have costs

    qc=2N*x^2+2N*x*gamma+(r/2)*gamma^2,
    qd=r*z/8.                                                    (2)

Each is an independent physically Lorentzian block: an additional
actual wall is strict at the simple old vertex and has a nonzero
timelike perpendicular component. Nonphysical admissibility therefore
requires qc'<1 and qd'<1. In particular

    0<z'<8/r<=4/5,       -4<delta'<-4+8/r<0.                      (3)

Cyclic integrality makes z integral. Since qd depends rationally on
z, an embedding fixing z would preserve its negative physical
determinant. Consequently K=Q(z)=Q(delta). All other conjugates
in(3) are genuine embeddings of this same K.

Squared contacts and actual cycle products belong to K; unit-normal
x or gamma is not assumed to do so. With gamma>0 put

    w=4*x^2,       q=4*x*gamma,       v=4*gamma^2.

These are integral: w,v are length-two cycles and -q is the actual
a-e-c triangle cycle in the norm-two Gram. In particular q^2=wv.
At every nonphysical embedding0<v'<4, since the actual e,c
finite block is positive definite after arithmetic normalization.
All formulas are determinant-ratio and cyclic identities, invariant
under positive diagonal changes of arithmetic root normalization.

## 2. No e contact: the root's uniform denominator lemma

First let gamma=0. Then w=4*x^2 is integral, qc=N*w/2 and
K=Q(w). Nonphysical positivity gives0<w'<2/N. The case x=0
is physically impossible, since c would be orthogonal to the whole
positive20-dimensional old vertex.

If x is a finite native cosine,1<=w<4. Degree d>=2 gives

    0<Norm(w)<4*(2/N)^(d-1)<=1,

a contradiction. Thus K=Q, outside the eight fields.

For x>=1, use (1) to obtain

    (N*w-r)*delta=4N.

Set g=gcd(N,r)=gcd(r,2), M=N/g, k=r/g and A=M*w-k.
Then A is integral and

    A*delta=4M,       gcd(M,k)=1.

The ideals(A) and(M) are coprime: modulo M, A=-k, and integer
Bezout for M,k gives a ring-of-integers Bezout identity for A,M.
Thus(A) divides the ideal(4). Equivalently4/A is integral;
this also follows by multiplying that Bezout identity by4 and
using A*delta=4M. Hence

    |Norm(A)|<=4^d.

Physically A>4, since w>=4 gives A>=(3r-8)/g>4. At all other
embeddings0<w'<2/N gives

    -r/g<A'<-(r-2)/g,       |A'|>N/g>=4.

Therefore |Norm(A)|>4^d, the required contradiction. This dotted
subcase uses no ground-field list and no assumption that x itself
belongs to K. The argument is uniform over all r>=10 for which the
stated actual pattern occurs.

## 3. gamma>0: zero and finite x, with the small-rank checks

For x=0, qc=r*gamma^2/2>1 because r>=10 and gamma>=1/2.
Thus K=Q(gamma^2), with every other gamma-square in(0,2/r),
contained in(0,1/5). The exact cosine-square trace

    Tr(cos^2(pi/ell))=d/2+mu(ell)/4,
    d=phi(ell)/2,

forces d<=3: for d>=4 the lower trace d/2-1/4 exceeds
1+(d-1)/5. The complete phi(ell)<=6 labels are
3,4,5,6,7,8,9,10,12,14,18. The degree-three labels and10
have a nonphysical cosine-square greater than1/5. Rational labels
3,4,6 force Q; labels5,8 force the excluded exact fields
Q(sqrt5),Q(sqrt2). Only ell=12, K=Q(sqrt3), survives at this
stage.

Now let x=cos(pi/j) be finite nonzero. Completing the square in
(2), using q^2=wv at each embedding, gives

    qc'>=(4N/r)*(x^2)'.

For j>=4 any embedding fixing x^2 would give qc'>=2N/r>1.
Hence K=Q(x^2), and all its other conjugates satisfy

    0<(x^2)'<r/(4N)<=5/16.

The same exact trace gives degree at most4: at degree5 the lower
bound is9/4 and the strict upper bound is also9/4, already
impossible. The complete phi(j)<=8 list and the explicit conjugate
test leave only j=4,5,6,8,12. To spell out that finite test, the
other labels7,9,10,14,15,16,18,20,24,30 have a primitive exponent
a with a/j<=3/10, respectively a=2,2,3,3,2,3,5,3,5,7. Their
conjugate square is at least(5-sqrt5)/8>5/16.

Among j=4,5,6,8,12, the eight-field equality K=Q(x^2) leaves
only j=12 and K=Q(sqrt3). Finite gamma-square in K permits
ell=3,4,6,12. Since x=cos(pi/12) is in the nontrivial square
class sqrt2*K and x*gamma is in K, only ell=4 or12 remain.
Equal label12 has nonphysical cost

    qc'=(9r-16)*(2-sqrt3)/8>1,

already at r=10. For example(2-sqrt3)/4>1/16 and the remaining
factor is at least37. Thus only (j,ell)=(12,4) remains here.

If j=3, x=1/2 and the cycle puts gamma in K. Formula(2) factors

    qc-1=(gamma+1)*(r*gamma+r-4)/2.

All nonphysical gamma conjugates therefore lie in
(-1,-1+4/r), a subset of(-1,-3/5). The cosine trace is at
least-1/2, whereas for degree d>=4 its strict upper bound is
1-(3/5)*(d-1)<=-4/5. Hence d<=3, giving ell=3,4,5,6,7,9.
Label3 is the affine a-e-c triangle; label4 forces Q(sqrt2).
Labels5,7,9 have a conjugate greater than-3/5. Only ell=6,
K=Q(sqrt3), remains.

The three residual profiles

    (x,e-label)=(0,12),(1/2,6),(cos(pi/12),4)

have exactly the same R in(1):

    R0=(7+4sqrt3)/4.

Rationalizing z=4+1/(R0-r/(4N)) at its two embeddings gives

    Tr(z)=4+4*(7r-15)/(3r^2-6r-1).

Its added fraction is strictly between0 and1 for every r>=10,
since the denominator minus numerator is3r^2-34r+59>0.
Thus z is not integral. This closes ALL zero and finite nonzero
x cases with an e contact, without importing an old n=19 threshold.

## 4. Dotted x: a general prime-ideal norm inequality

Assume x>=1 and gamma>0. Let u=v-4, whose conjugates are all
strictly negative, and let U=|Norm(u)|. Define the integral scalar

    B=N*(4w+4q+2v-4)+2u.

Expanding (1) exactly gives

    B*delta=-4N*u,       B=N*C+2u,
    C=4w+4q+2v-4 in O_K.                            (4)

**Divisibility lemma.** The ideal(B) divides the ideal(8u^2).
Here no principal-ideal-domain assumption is used. At any prime
ideal write a=valuation(2), n=valuation(N), t0=valuation(2u),
b=valuation(B). If n>t0, (4) gives b=t0. If n<=t0, integrality
of delta and (4) give

    b<=valuation(4Nu)=a+n+t0<=a+2t0=valuation(8u^2).

The first case also satisfies the last inequality. This proves the
claim at every prime ideal.

Taking norms in(4) now gives

    |Norm(delta)|>=(N/2)^d/U.

The other d-1 conjugates of delta have absolute value below4
by(3). Consequently

    delta>4*(N/8)^d/U.                              (5)

On the physical side x>=1 and gamma>=1/2 give
R>25/12 unless at the boundary, and in all cases

    R-s0>=25/12-5/16=85/48>1,
    0<delta<=48/85<1.                               (6)

We next compose(5) with the COMPLETE native e-label possibilities
in the eight fields, not a guessed finite label cutoff.

## 5. Native labels and exact norms of u

The real cyclotomic field F_ell=Q(cos(2pi/ell)) is contained in K
because gamma^2 is a cyclic scalar. Since [K:Q]<=5, phi(ell)<=10.
The complete labels after excluding7,14 (their K7 field is not a
subfield of any of the eight) are

    3,4,5,6,8,9,10,11,12,15,16,18,20,22,24,30.

This is just the elementary phi bound plus the quadratic subfields
and conductors of the eight named fields, not a new arithmetic
classification. The rational labels3,4 are handled in Section6.

For u=-4*sin^2(pi/ell),

    |Norm_Fell/Q(u)|=Phi_ell(1),
    U=Phi_ell(1)^[K:Fell].                          (7)

Pairing conjugate factors1-zeta,1-zeta^(-1) proves(7) directly.
The cyclotomic value is1 unless ell is a prime power p^j, when
it is p. The only nontrivial remaining entries are therefore:

|ell|possible K|d|U|
|---:|---|---:|---:|
|5|K15 or K20|4|25|
|8|K16 or Q(sqrt2,sqrt3)|4|4|
|9|K9|3|3|
|11|K11|5|11|
|16|K16|4|2|

For every remaining non-prime-power label, including6, U=1.
Equations(5)--(6) exclude U=1, and also the ell8,9,16 rows:
their lower bounds are respectively at least1,4/3,2.

For ell11, gamma>2/3 gives R>16/5. Hence

    delta<80/231<4/11,

contradicting(5). No approximate angle evaluation is used.

For ell5, gamma=(1+sqrt5)/4>21/26. Indeed sqrt5>29/13 by
845>841. Monotonicity of(1+gamma/2)^2/(1-gamma^2) gives

    R>5329/940,       R-s0>16/3,       delta<3/16.                 (8)

If N>=9, (5) gives

    delta>4*(9/8)^4/25=6561/25600>1/4,

contrary to(8). If N=8, then r=10,d=4,U=25. The norm bounds are

    |Norm(delta)|>=256/25>10,
    |Norm(delta)|<64*(3/16)=12.

All three nonphysical delta conjugates are negative, so
Norm(delta)=-11. This is impossible in either K15 or K20:
rational11 is unramified and its residue degree in each real
cyclotomic field is2, since11^2=1 modulo15 or20 but11 is not
congruent to either1 or-1. Every prime ideal above11 thus has norm
11^2. An integral element with absolute norm11 would give a
principal ideal of norm11, which cannot exist. This is the usual
cyclotomic Frobenius/decomposition law; the needed finite orders
have been given explicitly. It closes the last irrational label.

## 6. The rational e labels3 and4

These short arguments are also due to the root. They avoid a
potentially weak extension-field norm bound for rational u.

For ell4, let A=2*sqrt2*x=q and Z=A+1. The ACTUAL triangle cycle
makes A integral, and

    R=Z^2/4,       Z'^2<2/N.

Put g=gcd(N,r), B=(N/g)*Z^2-r/g. Equation(1) gives
B*delta=4N/g. The same integer Bezout argument as Section2 gives
(B) dividing(4). Physically B>4; every other conjugate has
|B'|>N/g>=4. Thus its norm is simultaneously greater than and
at most4^d, impossible.

For ell3, let Z=4x+1=2q+1, an integral scalar from the same
actual cycle. Then

    R=Z^2/12,       Z'^2<6/N.

Put g=gcd(N,3r)=gcd(N,6) and B=(N/g)*Z^2-3r/g.
Now B*delta=12N/g and the coprimality argument gives(B)|(12).
For r=10,12,13,15,16, every nonphysical |B'|>3N/g>=12,
and the physical B>12, contradicting |Norm(B)|<=12^d.
For r=11, B=3Z^2-11 is coprime to3; thus(B)|(4), while all
absolute conjugates exceed4. For r=14, B=2Z^2-7 is coprime
to2; thus(B)|(3), while all absolute conjugates exceed3.
These cover exactly r=10,...,16. Native ell6 was already covered
by U=1 in Section5. The cycle scalings are explicit; no claim
that a unit-normal x itself is integral has been made.

## 7. Scope and the remaining hard cases

Sections2--6 exclude every actual isolated-d quad in which c has
no D_r long-arm contact and no spectator contact, including zero,
finite and dotted a.c. This is a complete restricted consumer
branch, with the actual whole-factor hypotheses stated in Section1.

Still required for the ten proposed half refunds are the isolated-d
branches where c touches the long arm, with or without e or a
spectator bridge, and the final all-row/triangle/simultaneous-source
composition. The separate spectator-only and bridge-d reports have
their own review status. No old n=19/20 theorem has been asserted
at these smaller ranks merely by analogy.

Mathematical inputs are finite Coxeter classification, the compact
actual-face and positive-null arguments, arithmetic cyclic integrality
and strict nonphysical positive definiteness, elementary cosine
traces and cyclotomic integral/norm and residue-degree facts. The
existing native eight-field source is
`fresh_arithmetic_high_label_20261016.md`; the old finite-c
identities are also displayed in
`arithmetic_dn_fork_whisker_quad_attachment_classification_20261019.md`.
The smaller-rank thresholds and the new ideal arguments are supplied
here explicitly. This source has no execution or Lean claim and
does not change any dimension table or installed weight rule.
