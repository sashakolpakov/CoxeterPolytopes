# Arithmetic D16 plus actual I: the quartic norm-two fork branch is impossible

22 October 2026. **Complete manual candidate, UNREVIEWED pending full
outside-author review and root acceptance.** No worker, enumeration,
symbolic or numerical process, mathematical import, test, Lean or TeX
ran. The six paper editions, frozen inputs and public bounds are unchanged.

**Theorem at the exact written inputs below.** At an actual farthest
vertex with ENTIRE forest D16 plus actual isolated short I in a compact
arithmetic Coxeter17-polytope, the ACTUAL field cannot have degree4.
Together with the accepted preceding reduction, its degree is therefore
2 or3, and both fork opposites are isolated within D16 with dotted own
contacts. The two surviving parameters need not be equal. No old-I or
I-opposite contact is removed in these remaining degrees. The entire
forest is not excluded, and no dimension ceiling or example is claimed.

## 1. Exact input boundary

All paths in this table are under `results/`.

| Source | SHA-256 | Exact use |
|---|---|---|
| `arithmetic17_d16_actual_i_all_forks_isolated_and_quartic_units_candidate_20261022.md` | `45e56bd3cb9e9b01593134206c81a98268b1f1f4d3677572b4c61db7758abdfd` | Entire-forest degrees2/3/4; in degree4 each fork gives the actual-field generator and both norms in(1) |
| `arithmetic17_d16_actual_i_quartic_isolated_fork_contact_and_norm_two_reduction_candidate_20261022.md` | `7b7d4fa273b2322ab7953666394cefd5a94874bffa50d5aebfed317c0e67ee3d` | Direct proof of the quartic strip, zero OLD-I contact and joint norm-two conclusion |
| `review_arithmetic17_d16_all_nonisolated_forks_and_quartic_units_root_acceptance_20261022.md` | `ff9d259455913f3e291cb7836a8adc6f0f95503b0bb2dea4e2a60733ac9d9bfa` | Recorded root acceptance, including the separately reviewed cubic input and the unit corollary |

The third row is an acceptance record, not an additional mathematical
premise. The preparation headers of the first two sources are historical;
the record binds their later source-bound review and acceptance. Their
actual-geometry, K-root-line and integral-cycle interfaces remain explicit.
No relative D/I root-length squareclass, arithmeticity of a face or
assumed list of quartic fields is added here. The quadratic-form model
is automatic for arithmetic reflection lattices, not an extra subclass.

I authored the independently reviewed norm-two input. Root supplied the
subsequent unit corollary and proposed using its coefficient relation.
The complete coefficient bounds and interval-root contradiction below
are my new manual argument. Root and the composition agent supplied
independent crosschecks before freezing; this is not an outside-author
review of my new proof. Independent full review remains required.

## 2. A self-contained algebraic lemma

There is no algebraic integer T generating a totally real quartic field
whose physical value and other three conjugates x1,x2,x3 satisfy

    4<T<1369/32,       0<x_i<1/2,
    N(T)=2,           N(2T-1)=-1.                 (1)

All norms in this lemma are over that quartic field. In particular, the
generator hypothesis makes its minimal polynomial have degree exactly4.
We prove the lemma without enumerating fields or running a polynomial
search.

### 2.1. The trace and two intermediate coefficients have only nine possibilities

Write the monic minimal polynomial as

    P(X)=X^4-A X^3+B X^2-C X+2,   A,B,C in Z.

Evaluating at1/2, with the signs fixed by the even degree, gives

    -1=N(2T-1)=16P(1/2)=33-2A+4B-8C.

Consequently

    A-2B+4C=17,       B=(A+4C-17)/2,            (2)

and A is odd. The product of the small conjugates is2/T. Since each is
less than1/2, T>16. Apply AMGM first to the three x_i and then to their
three positive deficits1-2x_i. It gives

    1=(2T-1)*product(1-2x_i)
      <2*[T^(1/3)-2*2^(1/3)]^3.                 (3)

The bracket is positive because T>16. Put r=2^(1/3)>5/4. Equation(3)
implies

    T^(1/3)>2r+1/r>33/10,
    T>35937/1000.                                (4)

Here the function2r+1/r is strictly increasing for r>1, and its value
at5/4 is33/10. A second AMGM application gives

    x1+x2+x3>=3*(2/T)^(1/3)>27/25.              (5)

The last strict bound follows from2/T>64/1369>(9/25)^3: the integer
comparison is1000000>998001. Combining(4),(5) with the upper strip,

    37<37017/1000<A=T+x1+x2+x3
      <1369/32+3/2=1417/32<45.

Thus the odd integer A belongs to{39,41,43}.

For the coefficient C, the elementary symmetric formula and AMGM give

    C=T*(x1*x2+x1*x3+x2*x3)+2/T
      >3*2^(2/3)*T^(1/3)
      >3*(25/16)*(33/10)=2475/160>15.           (6)

For the upper bound set Y_i=1/x_i>2. Their product is T/2. For any
three positive numbers y_i>a, one has

    y1+y2+y3<2a+y1*y2*y3/a^2.

Indeed, replacing y1,y2 by a,y1*y2/a strictly increases the sum by
(y1-a)(y2-a)/a; merging the last two factors once more proves the
claim. At a=2 this yields

    C=2*(Y1+Y2+Y3)+2/T
      <8+T/4+2/T
      <8+1369/128+1/8=2409/128<19.              (7)

Here T>16 justifies2/T<1/8. Therefore C is one of16,17,18. Together
with(2), this leaves exactly nine possible integer polynomials. It
does not assert that any of the nine is irreducible or geometrically
realized.

### 2.2. Move the small-root interval to the positive half-line

Define

    Q(y)=16(1+y)^4*P(y/[2(1+y)]).

Expanding and using(2) gives exactly

    Q(y)=-y^4+(2A-8C+60)y^3
              +(2A-16C+158)y^2+(128-8C)y+32.   (8)

The substitution y=2x/(1-2x) maps(0,1/2) bijectively to(0,infinity).
Since the quartic generator is separable, its three small conjugates
would give THREE DISTINCT POSITIVE roots of Q. The following complete
table lists all its remaining possible coefficients; its leading and
constant coefficients are always-1 and32.

| A | C | B from(2) | coefficient of y^3 | coefficient of y^2 | coefficient of y |
|---:|---:|---:|---:|---:|---:|
|39|16|43|10|-20|0|
|41|16|44|14|-16|0|
|43|16|45|18|-12|0|
|39|17|45|2|-36|-8|
|41|17|46|6|-32|-8|
|43|17|47|10|-28|-8|
|39|18|47|-6|-52|-16|
|41|18|48|-2|-48|-16|
|43|18|49|2|-44|-16|

We show by exact inequalities that each row has only ONE positive root.

### 2.3. The three rows with C=16

At fixed C and y>=0, Q increases with A. On0<=y<=3 its smallest row
is therefore

    Q_(39,16)(y)=32-y^2*(y^2-10y+20).

Since y^2<=3y on this interval,

    y^2*(y^2-10y+20)<=y^2*(20-7y)
      <=32000/1323<32.                          (9)

The middle cubic has derivative y*(40-21y), so its maximum on[0,3]
is the stated value at40/21. Thus all three Q are strictly positive
throughout[0,3].

For any of these rows write Q=-y^4+a y^3+b y^2+32. Its derivative is

    Q'(y)=-y*L(y),       L(y)=4y^2-3a y-2b.

In the three rows L(0) is positive and L(3) is respectively-14,-58,-102.
Since L is an upward quadratic, it has exactly two positive roots
s1 in(0,3) and s2>3. Hence Q decreases to a positive local minimum at
s1, increases to a positive local maximum at s2, then decreases strictly
to minus infinity. It has exactly one positive root.

### 2.4. The three rows with C=17

The derivative is largest when A=43. For that row

    Q_(43,17)'(y)=-2*f(y),
    f(y)=2y^3-15y^2+28y+4.

The following exact identities prove f>0 on the whole nonnegative axis:

    f(y)=y*(2y-7)*(y-4)+4>=4       if0<=y<=7/2;
    f(y)=2*(y-4)^2*(y+1/2)+4*(y-3)>0  ify>=3.

For the second identity its two nonnegative summands cannot both vanish.
These intervals cover every y>=0. For the two smaller A values,

    Q_(A,17)'=Q_(43,17)'-2*(43-A)*(3y^2+2y)<0.

Thus each row decreases strictly from32 to minus infinity, giving
exactly one positive root.

### 2.5. The three rows with C=18

Again the derivative is largest when A=43, where

    Q_(43,18)'(y)=-2y*(2y^2-3y+44)-16<0

for y>=0: the quadratic has positive leading coefficient and discriminant
9-352=-343<0. Decreasing A only decreases the derivative by the same
nonnegative term2*(43-A)*(3y^2+2y). Every row therefore decreases
strictly from32 to minus infinity and has exactly one positive root.

The nine cases are exhausted. None can have the three distinct positive
roots required by Section2.2, proving the algebraic lemma.

## 3. Apply the lemma to the actual fork invariant

Suppose the ENTIRE actual D16+I farthest forest has quartic field K.
The accepted reduction in Section1 supplies, separately for each actual
fork, the integral generator T_i=4gamma_i^2 satisfying every condition
of(1). The lemma excludes even one such generator. Thus the quartic
case is impossible. Applying the accepted entire-forest degree bound
and nonisolated-fork exclusions leaves actual degree2 or3 and both
forks isolated within D16 with dotted own contact.

The application adds no pressure inequalities or energies and therefore
does not combine incompatible projected points or charge I twice. All
geometric and field descent work is at the explicit accepted input
boundary; the new contradiction is elementary algebraic number theory.
The two fork parameters were never identified. Positive old-I contacts
in degrees2/3 and all possible contacts with the I opposite remain
subject to their earlier actual-wall constraints. Those degrees, the
whole forest, other arithmetic17 families and the global ceiling17
remain unresolved by this note. The new quartic exclusion itself
requires full outside-author review and root acceptance before promotion.
