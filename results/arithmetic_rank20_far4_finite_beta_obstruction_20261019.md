# Finite opposite products in the far4 first-d row

19 October 2026. Root-agent survivor and Schur arguments, independently
checked and combined with a complete native-square reduction by
batch_guard. Paper only: no enumeration, scientific job, calibration,
Lean, TeX or build was run. Existing sources remain unchanged.

This note excludes finite beta in first-d row B over the eight target
fields. It does not exclude beta>=1 in B, or the remaining F row,
and does not install a simultaneous refund or a dimension bound.

## 1. Scope, actual cyclic variable and full-star interface

The entire old compact finite vertex is D_r+T_s, r+s=20,
10<=r<=16. Here T_s is D_s,4<=s<=10, or E6,E7,E8. Let a be
the D_r fork, b,e its short leaves and l_1,...,l_L its long arm,
where N=r-2>=8,L=N-1. In the ACTUAL quadrilateral a,b,c,d with
complete white support W=A1(e)+A_L+T_s, assume

    d.l_L=-sqrt2/2,
    d is orthogonal to every other old root except b,
    c.b=c.d=0.

Write b.d=-beta, a.c=-x. The actual acute Coxeter-wall rule allows
beta=0, a finite native cosine cos(pi/m),m>=3, or beta>=1;
the present statement concerns beta<1 only. Similarly x=0 or
x>=1/2, and 0<x<1 is a native finite cosine.

The simple old path from b to l_L and the terminal4 bond to d
give the ACTUAL integral cyclic variable

    A=2sqrt2*beta in O_K.

We do NOT assert beta in K. Projecting d against the FULL old
rank20 vertex gives

    q_d=r*beta^2/2+sqrt2*beta+1
        =r*A^2/16+A/2+1.                           (1)

Because d is strict at the actual simple old vertex, this
old-plus-d block is physically independent and Lorentzian, so
q_d>1. Its nonphysical blocks are positive definite. Thus beta=0
is impossible and, at every nonphysical embedding,

    -8/r<A_i<0.                                   (2)

Also K=Q(A): if Q(A) were a proper subfield, a nonphysical K
embedding extending its physical embedding would repeat q_d>1.
This uses the ENTIRE old star; an unspecified embedded component
would not suffice. These are the inherited actual arithmetic
admissibility interfaces, not a numerical Gram-feasibility test.

## 2. The cosine-square field cannot be a proper subfield

Suppose beta=cos(pi/m),m>=3, and put

    F=Q(A^2)=Q(cos(2pi/m)).

Since A^2 belongs to F, [K:F]<=2. If A were not in F, the
quadratic automorphism of K=F(A) would send A to-A. Its physical
value is then -A_0<=-sqrt2, because beta>=1/2. That nonphysical
conjugate violates(2), since r>=10 gives8/r<=4/5<sqrt2.
Consequently

    K=F=Q(A)=Q(cos(2pi/m)).                         (3)

This step does not assume a square root lies in its square field;
it PROVES that necessary condition from the conjugate bound.
In particular it allows the possibility that a weighted cosine,
such as the m=20 scalar, does belong to a smaller field than beta.
The subsequent square-conjugate test still checks that case.

The eight target fields are
Q(sqrt3),Q(sqrt6),K9,K11,K15,K16,K20,K24. They are nonrational
and have degree at most5. The standard real-cyclotomic degree
formula and(3) give phi(m)<=10. The COMPLETE possible list is

    m=3,4,5,6,7,8,9,10,11,12,14,15,16,18,20,22,24,30.

For completeness, every prime divisor has p-1<=10, so it is
2,3,5,7 or11. Their maximal exponents are4,2,1,1,1;
the multiplicative Euler factors give precisely the displayed
integers m>=3. This is an elementary finite list, not a worker
enumeration. For m=3,4,6 the square field F is Q, contradicting
the target field list via(3).

## 3. All nonrational labels except12 fail a square conjugate

For each F embedding given by a residue k coprime to m,

    sigma_k(A^2)=4+4cos(2k*pi/m).

Every nonphysical value would have to satisfy, by(2),

    sigma_k(A^2)<64/r^2<=16/25.                     (4)

But the following complete table violates(4) for every remaining
label other than12:

|m|k|lower bound for sigma_k(A^2)|
|---|---|---|
|5|2|3-sqrt5>3/4|
|7|2|greater than2|
|8|3|4-2sqrt2>1|
|9|2|greater than2|
|10|3|greater than2|
|11|2|greater than2|
|14|3|greater than2|
|15|2|greater than2|
|16|3|greater than2|
|18|5|greater than2|
|20|3|greater than2|
|22|3|greater than2|
|24|5|greater than2|
|30|7|greater than2|

Each k is coprime to m and not congruent to +/-1, so the
embedding is nonphysical. For the entries bounded by2, k/m<1/3,
hence cos(2k*pi/m)>-1/2. The other two bounds use sqrt5<9/4
and sqrt2<3/2. All comparisons are exact. In particular m=20
is explicitly retained and rejected by its square conjugate; it
is not discarded from an incorrect beta-field inventory.

The sole remaining label is m=12. Then

    beta=cos(pi/12),
    A=1+sqrt3,       K=Q(sqrt3),
    A'=1-sqrt3.

For r>=11, (2) fails because sqrt3-1>8/11, equivalent to
363>361 after squaring. Thus the ONLY necessary finite-beta
profile left is

    N=8,r=10,T_s=D10, beta=cos(pi/12),K=Q(sqrt3).    (5)

This is an arithmetic profile, not an existence assertion.

## 4. Complete finite c inventory in the exceptional profile

In(5), W+d=A1(e)+B8+D10. The terminal4 bond of B8 is at
the d end. Finite W+c+d with c.d=0 leaves c either

* with no old long-path contact; it may contact the spectator,
  e, or neither, subject to the full finite spectator table;
* simply attached to the near long endpoint l_1, with no
  spectator contact, and e contact absent or simple.

Indeed finite extension of the large B component remains a B
path. An internal contact makes a fork with a4 bond, a new near4
bond gives two4 bonds, and the far endpoint d is prohibited.
A c simultaneously bridging B8 to D10 also leaves a4 bond
and a fork in connected rank19; no finite diagram has that form.
If c contacts D10 but not the path, it attaches to its long
endpoint simply; an additional e contact is simple too. Thus
its spectator residual t_T is rational. With no spectator
contact set t_T=1. For any no-path c,

    t_c=t_T-gamma^2,       gamma=-c.e,
    X=-a_bar.c_bar=x+gamma/2.

If the spectator is contacted, gamma is0 or1/2. If not, gamma
is zero or a finite native cosine. This keeps all actual white
contacts rather than silently assuming spectator orthogonality.

## 5. Exact physical Schur scalars

Set z=sqrt3 and epsilon=sqrt2/16. After projecting off W,

    a_bar^2=5/16, b^2=1, d_bar^2=1/8,
    a_bar.b=-1/2, a_bar.d_bar=-epsilon,
    b.d_bar=-beta,
    Delta=beta^2-1/8=(3+2z)/8>0.

Eliminate the b,d block [[1,-beta],[-beta,1/8]], whose
determinant is-Delta. The remaining a scalar is

    R0=5/16+[(1/8)/4+beta*epsilon+epsilon^2]/Delta
       =(2+z)/8=beta^2/2<1/2.                      (6)

These are PHYSICAL real unit-normal calculations; they do not
pretend that the unit-normal beta itself belongs to Q(sqrt3).
Only the displayed simplified expressions are used below.

## 6. No path contact

Here c is orthogonal to b,d even after projection. Actual corank
therefore gives

    (x+gamma/2)^2=R0*t_c.                          (7)

Since t_c<=1 and R0<1/2, the native gap permits only x=0
or x=1/2: every higher positive native value is at least
cos(pi/4), with square1/2. At x=1/2, any gamma>0 would make
the left side at least9/16>1/2. Thus gamma=0. But then

    t_c=1/(4R0)=4-2sqrt3,

which is irrational, contrary to the rational spectator residual.

At x=0 with a spectator contact, gamma=0 makes(7) impossible.
The only other allowed value gamma=1/2 would force

    t_c=1/(16R0)=1-sqrt3/2,

again irrational, while t_T-1/4 is rational.

At x=0 without a spectator contact, t_c=1-gamma^2. Equation(7)
forces

    4gamma^2=(20+8sqrt3)/13.

But 4gamma^2 is an actual integral two-cycle, whereas this
quadratic number has trace40/13, not an integer. It is not an
algebraic integer. This also covers gamma=0, which already fails(7).
Every no-path continuation is excluded.

## 7. Near path contact

Now c.l_1=-1/2, gamma is0 or1/2, and there is no spectator
contact. The projected data are

    t_c=9/16-gamma^2,
    X=x+7/16+gamma/2,
    -c_bar.d_bar=epsilon.

The same b,d elimination gives the exact cross correction and
remaining c scalar

    J=(beta*epsilon/2+epsilon^2)/Delta=1/16,
    t_c+epsilon^2/Delta=(12+sqrt3)/24-gamma^2.

Actual corank is therefore

    (x+3/8+gamma/2)^2
       =R0*[(12+sqrt3)/24-gamma^2].                (8)

If x>0, the native gap makes the left side at least49/64.
The right side is less than(1/2)*(14/24)=7/24, impossible.
At x=0,gamma=0, the right side is

    (27+14sqrt3)/192>27/192=9/64,

strictly greater than the left side. At x=0,gamma=1/2, the
left side is25/64=75/192, while the right side is

    (15+8sqrt3)/192<31/192.

No equality in(8) is possible. This closes the exceptional
profile(5) with the full finite c inventory.

## 8. Exact conclusion and remaining work

Every finite beta in the far4 first-d row B is excluded over
the eight target fields. The proof uses the ACTUAL weighted cyclic
variable A, a complete cosine-SQUARE field inventory including
label20, and a full corner check for the only remaining profile.
No explicit integral basis or enumeration result is assumed.

The B branch beta>=1 and F branch beta>=1 remain outside this
theorem. Global triangle consumers and simultaneous half-refund
composition are also not settled by it. No rank20 allowance or
compact arithmetic dimension ceiling is changed here.
