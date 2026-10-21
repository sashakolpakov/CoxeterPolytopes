# A global nine-direction positive-part inequality

20 October 2026. **Complete analytic proof; awaiting independent review.**
This replaces the rational interval step in the separate nine-direction
proof by one polynomial inequality. The earlier proof and certificate
remain unchanged and give separate support. No scientific calculation,
source import, numerical experiment, enumeration, optimizer, Lean process
or build ran while deriving this proof. The companion exact-polynomial
checker is source only until the coordinating agent schedules it.

Let p1,...,p11 be arbitrary real numbers satisfying

    sum p_i=3,       sum p_i^2=1,

and let I contain nine coordinates. For every triple T in I put
r_T=2-sum_(i in T)p_i. Then

    sum_(T subset I, |T|=3) (r_T^-2-37/100)_+ <=1323/25.       (1)

More precisely, with delta=3-sum_I p_i,

    (1/84) sum_T (r_T^-2-37/100)_+
         <=63/100-delta/18.                                  (2)

Equality in (1) holds exactly when the nine coordinates in I are1/3
and the omitted two coordinates are0. This is a scalar frame theorem,
not a new dimension theorem; other subset sizes and the geometric
transport composition are separate obligations.

## 1. Exact frame bounds

Write S=sum_I p_i, q=sum_I p_i^2, and V=q-S^2/9. The sum of the
omitted coordinates is delta=3-S, so their squared sum is at least
delta^2/2. Thus

    0<=V<=2delta/3-11delta^2/18,
    0<=delta<=12/11,
    V<=2/11.                                                  (3)

The interval for delta follows by requiring the quadratic upper bound
for V to be nonnegative. Put

    mu=1+delta/3,       y_i=p_i-S/9 (i in I).

Then sum y_i=0, sum y_i^2=V and 1<=mu<=15/11. For a uniform random
triple in I,

    r=mu+d,       d=-sum_T y_i,       E d=0.

All r are positive, with the uniform bound r>11/20. Indeed after
centering all eleven coordinates by3/11, their squared norm is2/11.
The centered indicator of any triple has squared norm24/11. Cauchy
therefore gives

    sum_T p_i <=(9+4sqrt3)/11,
    r >=(13-4sqrt3)/11 >11/20.                                 (4)

The last comparison is equivalent to139>80sqrt3, and its squared
integer comparison is19321>19200. Negative individual p_i are allowed
throughout; no clipping or unproved monotonicity in the coordinates
has been used.

## 2. The three centered moments used below

Let T3=sum y_i^3 and T4=sum y_i^4. The probabilities that a given
set of one, two or three distinct coordinates belongs to the random
triple are1/3,1/12,1/84, respectively. Expanding powers and using
sum y_i=0 gives

    E d^2=V/4,
    E d^3=-3T3/28,
    E d^4=5V^2/28-3T4/28<=V^2/6.                              (5)

For the last inequality T4>=V^2/9 by Cauchy. For example the fourth
moment follows by substituting

    sum_(i!=j) y_i^3 y_j=-T4,
    sum_(i<j) y_i^2 y_j^2=(V^2-T4)/2,
    sum_i y_i^2 sum_(j<k, j,k!=i)y_j y_k=T4-V^2/2

in the expansion. No independence assumption is made for the triples.

We also need the elementary sharp cubic estimate

    |T3|<=7V^(3/2)/sqrt72.                                    (6)

For V=0 it is immediate. Otherwise each y_i<=M=sqrt(8V/9), since
the other eight entries sum to -y_i. Put b=-M/8. Then

    0<=sum_i(M-y_i)(y_i-b)^2
      =-T3+7MV/8.

Applying the same argument to -y gives(6). Consequently

    (3/7)T3 <= V^(3/2)/(2sqrt2) <=V/(2sqrt11)<V/6              (7)

when V>0, using V<=2/11. The weak inequality with V/6 also holds
at V=0. This deliberately uses a loose bound; a joint cubic/quartic
optimization is unnecessary here.

## 3. A global fourth-order clipping bound

Put c=37/100 and b=10/sqrt37, so b^-2=c. Since

    b>18/11,       mu<=15/11,

the number a=b-mu satisfies a>3/11. The first comparison follows
from100*121>324*37. For r>=b, the concave function c-r^-2 lies
below its tangent at b. Hence for all positive r,

    (c-r^-2)_+ <=(2/b^3)(r-b)_+.

For every real d and a>0,

    (d-a)_+ <=27d^4/(256a^3).                                 (8)

For d<=a this is trivial; for d>a the ratio(d-a)/d^4 attains its
maximum at d=4a/3. Substituting d=r-mu and the bounds for a,b gives

    (c-r^-2)_+
       <=[54/(256b^3 a^3)]d^4
        <[11^6/(128*18^3)]d^4 <3d^4.                          (9)

The last integer comparison is1771561<2239488. The estimate holds
globally on the full frame domain, including triples for which the
positive part vanishes. It is not an upper estimate only on a
selected high-r subset.

## 4. A single reciprocal remainder bound

For r=mu+d>0 the exact identity is

    r^-2=mu^-2-2d*mu^-3+3d^2*mu^-4-4d^3*mu^-5
                +[(mu+4r)/(mu^5 r^2)]d^4.                  (10)

Because mu>=1 and r>11/20, its remainder coefficient is at most

    r^-2+4r^-1 <400/121+80/11=1280/121.                    (11)

Use the identity

    (r^-2-c)_+=r^-2-c+(c-r^-2)_+.

Averaging(10), and then using(5),(7),(9),(11), proves

    E(r^-2-c)_++c
      <=mu^-2+3V/(4mu^4)+V/(6mu^5)+(1643/726)V^2.           (12)

Here1643/121=1280/121+3 combines the reciprocal remainder and
clipping coefficients, and E d^4<=V^2/6. Every replacement is in
the upper-bound direction. In particular no full reciprocal estimate
has been mistaken for the desired positive-part estimate.

## 5. The remaining polynomial inequality

Set x=delta/3, so 0<=x<=4/11, mu=1+x, and

    0<=V<=v=x(4-11x)/2.

The right side of(12) increases with V. Substituting v, multiply its
difference from1-x/6 by(1+x)^5. Exact expansion gives

    (1+x)^5*[1-x/6-( (1+x)^-2
                      +3v/(4(1+x)^4)+v/(6(1+x)^5)
                      +(1643/726)v^2 )]=x^2 R(x),

where

    R(x)=233/24+275x/24+10x^2/3+x^3/6-x^4/6
                  -(1643/2904)(4-11x)^2(1+x)^5.             (13)

On this interval x<1, so the first five terms are at least233/24.
Moreover

    d/dx[(4-11x)^2(1+x)^5]
       =(4-11x)(1+x)^4*(-2-77x)<=0.

The product is therefore at most its value16 at zero. Consequently

    R(x)>=233/24-(1643/2904)*16=635/968>0.                   (14)

Equations(12)--(14) prove

    E(r^-2-c)_++c <=1-x/6=1-delta/18,

which is(2). Multiplying by84 gives(1), since84*(63/100)=1323/25.
The proof is strict for delta>0. If delta=0, the inside Cauchy
inequality forces all nine entries to be1/3 and leaves zero squared
norm for the omitted two. Conversely that frame has all84 ratios1
and attains(1). This proves the asserted equality characterization.

## Verification and scope

The companion
`scripts/check_mixed_eleven_nine_direction_global_20261020.py`
checks the exact rational polynomial identities and constants, not
the geometry of opposition frames. Its source was prepared without
execution. The earlier interval certificate is independently retained
at `results/mixed_eleven_nine_direction_subset_20261020.json`.

Only the full eleven-coordinate sum and norm hypotheses and a
nine-element index set were used. The statement is therefore available
for actual P3 opposition frames once their earlier geometric construction
is supplied. It neither asserts the ten-direction case nor replaces the
separate actual-subset/unweighted-summation proof needed for a dimension
conclusion.
