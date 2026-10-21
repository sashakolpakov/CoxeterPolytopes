# A far 4-bond quarter refund at whole B20 vertices

19 October2026. **New written consumer theorem candidate, awaiting
independent review.** Source and paper work only; no computation,
enumeration, optimizer, build or Lean process was run. Earlier stable
A/B/D evidence is unchanged.

## 1. Scope and the consumer that must actually be paid

Let P be a compact arithmetic simplest-type Coxeter polytope of
dimension20 over one of

    Q(sqrt(3)), Q(sqrt(6)), Q(sqrt(2),sqrt(3)),
    K9,K11,K15,K16,K20,       Km=Q(cos(2pi/m)).

Begin with the eight-field distance-at-most-two charge plus its
entire-I2(m>=6) unit, and the reviewed whole-A20, whole-D20 and
whole-B20 simple-end reductions. At a vertex whose ENTIRE diagram
is B20, subtract1/4 at the pair forming its unique 4-bond. The
baseline at this pair becomes3/4. The simple-end changes occur at
different marked pairs. This would lower the whole B20 cost from
141/4 to35.

The claim concerns this far-4 quarter alone and its composition with
the stated rules. No global allowance35 or bound19 follows: other
two-component forests still have cost36. The spectator half-refund
at whole B19+A1, if separately reviewed, has different white graph
and does not alter any face studied here.

The white support of the present selected angle is A18. This agrees
with the D20 SHORT-LEAF-PAIR support. Thus D-whisker mixing must be
included; it cannot be dismissed by separate consumer proofs. The
other changed roles have white graphs A1+A17, A1+B17, B18 or
B16+2A1 and hence do not occur in this face. The B18 positive coupons
likewise do not occur here.

Any quad with at least two selected far-4 angles retains at least
3/2. If another corner has an unchanged positive charge, it alone
retains at least1. A newly unpaid face therefore has exactly one
selected far-4 corner and all three other corners are either original
zeros or deleted D20 short-leaf pairs. The proof below excludes
precisely those configurations. Triangles retain their entire high
I2 weight2, and all remaining weights are nonnegative.

## 2. Notation and the complete 21-row finite inventory

Put

    R=19,  t=10/19,  k=9/19,  e=1/38,  h=1/sqrt(2).

Write the old B20 corner as short root a and long root b, with
a.b=-h. Its complete white graph is the path

    W={v1,...,v18}=A18,

and b.v1=-1/2; a is white-isolated. Take cyclic sides a,b,c,d.
Thus a.c=-x and b.d=-y are the opposite native products. In a
potential unpaid face all other adjacent native products vanish:

    b.c=c.d=d.a=0.                                  (1)

This holds for original zeros and for deleted D-whisker corners,
whose two short leaves are orthogonal. A/D fork roles have different
white graphs and cannot supply additional new zeros.

The c corner next to b is either originally zero or a D-whisker
deletion. Complete finite one-root continuation of W+b=A19 gives

    I: no W contact;
    H: v1 simple (the NEW zero D-whisker corner);
    F: v18 simple;
    B: v18 with label4;
    T: v17 simple.

The original zeros omit contacts to v1, since that would have marked
distance2 from b. The new H case is essential. The (a,d) zero must
be original, since a has no white neighbor. Finite continuation of
W+d gives seven options

    I,H,J,K,F,B,T,

where J contacts v1 with label4 and K contacts v2 simply. The large
A path permits only endpoint3/4 or next-to-end3 extensions. Multiple
contacts create a cycle; exceptional finite diagrams are too small.

Finiteness of W+c+d and the zero marked charge at (c,d) reduce the
possible rows exactly to

|c profile|permitted d profiles|
|---|---|
|I|I,H,J,K,F,B,T|
|H|I,H,F,B,T|
|F|I,H,J,K,F|
|B|I,H|
|T|I,H|

There are21 rows. For H+H and F+F the opposite corner is a deleted
D20 short-leaf pair. Every other retained row has original-zero
opposite mark. To check completeness, after adjoining H or F the
large component is A19; a same-end simple attachment yields D20,
an opposite-end extension gives A20/B20/D20, and the remaining
attachments give a branched4-bond or arms1,2,and a long arm. After
B or T the component is B19 or D19, whose only further nonisolated
finite extension here is simple at the far end. These are ordinary
finite-classification statements with all18 whites present.

Eliminating W gives the following squares and cross products with
old b. The native b/root product is -y for a d wall and0 for c.

|profile|square u|white cross magnitude rho with b|
|---|---:|---:|
|I|1|0|
|H|t|k|
|J|1/R|sqrt(2)k|
|K|2/R|17/38|
|F|t|e|
|B|1/R|sqrt(2)e|
|T|2/R|1/R|

The old residual roots have a^2=1, b^2=t, a.b=-h.

## 3. Arithmetic and actual-wall conventions

Every new actual wall is independent from the positive20-dimensional
old vertex span: compact Coxeter simplicity prevents an extra wall
through that vertex. Its old-span projection cost q therefore exceeds1
physically; at every nonphysical embedding the corresponding rescaled
independent block is positive definite and q'<1. Formulas below are
determinant-ratio identities in cyclic products. They do not embed
arbitrary unit-normal square roots as though they belonged to K.

Any displayed field equality follows in the usual way: a nonidentity
extension embedding fixing all quantities in q would preserve q>1.
Every actual native product is0, a finite negative cosine, or dotted.
Native magnitude1 gives a feasible positive-null pair by acuteness
and is forbidden in a compact chamber.

We use the elementary cosine-moment lemma from Section5 of
`arithmetic_dn_fork_whisker_quad_attachment_classification_20261019.md`:
if x=cos(pi/m) is in K and all its nonphysical conjugates satisfy
|x'+kappa|<r, then its degree is at most3 when

    1/2+3kappa^2-2kappa-3r^2>0.

The resulting labels are3,4,5,6,7,9. This is an exact cyclotomic
trace argument, not an empirical label cutoff.

## 4. Five d-wall arithmetic facts

The facts below depend only on the old B20 and d; c is unrestricted.

**I.** If d=I, its old cost is q_d=2R y^2. Thus Z=4y^2 is integral,
K=Q(Z), and 0<Z'<2/19. Finite y would give Z<4 and nonzero integer
norm below4(2/19)<1 in degree at least2; degree1 gives the excluded
ground field Q. Hence y>1.

**H.** If d=H, q_d=k+2R(y+k)^2, so

    |y'+k|<sqrt(5)/19<1/8.                       (2)

The cycle b-v1-d makes 2y integral. Here K=Q(y). For finite y the
moment test at kappa=k and r^2=5/361 has coefficient133/722>0.
The six resulting labels all fail (2): label3 is rational; the
negative conjugates for4,6,9 and the positive conjugate for7 fail;
for5 its conjugate is greater than-5/16, so adding k>7/16 leaves
more than1/8. Native y=0 makes q_d=9 rational and likewise forces
K=Q. Thus y>1.

**F.** If d=F, q_d=k+2R(y+e)^2. The long simple cycle makes V=2y
integral, K=Q(V), and

    |V'|<(1+sqrt(20))/19<3/10.                  (3)

Finite y has V<2 and would give norm below2(3/10)<1 in every
nontrivial degree. The value y=0 makes q_d=1/2, contradicting the
physical independence of the new wall. Consequently y>1.

**B.** If d=B, U=2sqrt(2)y is integral from the native long cycle,
and q_d-1=(R/4)U^2+U. Thus U' lies in(-4/R,0). Finite U<3 has
impossible norm below12/19; rational U forces K=Q. Zero y makes
q_d=1 and a positive-null actual diagram. Hence y>1.

**T.** If d=T, V=2y is integral and q_d-1=2y(Ry+2). Its other
conjugates lie in(-4/R,0). Finite V<2 gives impossible norm below
8/19, and rational V forces K=Q. Zero y is positive-null. Thus y>1.

The last two norm estimates use K generated by U or V, respectively,
which follows from the physically negative determinant ratio.

## 5. Every non-I c profile has x>1

This finite-product reduction is needed before applying the uniform
Lorentz bounds. Put w=sqrt(2)x. Its integrality multiple is always
derived from an actual cycle, rather than guessed from a unit entry.

**c=H.** The cycle a-b-v1-c makes 2w integral. The old cost satisfies

    q_c-1=(w+1)(20w+16)/2.

Hence K=Q(w), and A=2(w+1) has other conjugates in(0,2/5). If x
is finite, A<2+2sqrt(2)<5, so degree at least3 would give norm below
5(2/5)^2<1. Only quadratic K=Q(sqrt(3)) or Q(sqrt(6)) can remain
among the eight. Native x^2 in such a field allows rational labels
3,4,6, and additionally12 in Q(sqrt(3)). Their w values are
sqrt(2)/2,1,sqrt(6)/2,(1+sqrt(3))/2. The first has the wrong field,
the second forces Q, and the other conjugates of the last two fail
-1<w'<-4/5. Zero x also gives rational q_c=9>1 and forces Q.
Thus x>1.

**c=F.** Again 2w is integral, and

    q_c-1=10w^2+w-1/2.

Thus K=Q(w) and |(2w)'|<(1+sqrt(21))/10<3/5. If x is finite,
2w<3, so degree at least4 gives norm below3(3/5)^3<1. Only the
two quadratic fields above and the cubic K9 remain. In the quadratic
case the same four native values just listed fail the conjugate
bound or the field equality. In K9 a native cosine-square must be
rational or have label9 or18; labels7,14 give K7, not K9. Rational
labels give quadratic w or rational w. For9 and18, respectively,

    N_K9/Q(2x^2)=1/8, 3/8.

Neither is a rational square, whereas w^2=2x^2 would make it the
square of N(w). These norm values follow from the cubic minimal
polynomials of 4x^2, whose constant terms are-1 and-3. Thus finite
x is impossible. Zero x makes q_c=1/2, violating physical positive
index. Therefore x>1.

**c=B.** The cycle has TWO 4-bonds, so its integral quantity is4x,
not2x. Here

    q_c-1=20x^2+2x,
    -2/5<(4x)'<0,       K=Q(x).

For finite x, 2<=4x<4. Degree at least3 gives norm below4(2/5)^2<1.
In degree2 the negative norm is-1 and the trace2 or3, giving fields
Q(sqrt(2)) or Q(sqrt(13)), neither allowed. Degree1 gives Q. Zero
x is positive-null. Hence x>1.

**c=T.** The cycle makes 2w integral and

    q_c-1=w(20w+4)/2,
    -2/5<(2w)'<0,       K=Q(w).

For finite x, 2w<3. Degree at least3 gives norm below3(2/5)^2<1.
In degree2 the negative norm is-1 and the trace can only be1 or2,
giving Q(sqrt(5)) or Q(sqrt(2)), outside the eight. Zero x is
positive-null. Thus x>1.

For c=I, the stronger x>1 is geometric: a and c are unchanged
white-isolated unit opposite sides. Finite or parallel products
contradict acute projection or compactness.

## 6. The seven rows c=I

If c=I, its old cost is q_c=20x^2. Thus Z=4x^2 is integral and
primitive in K, with 4<Z and 0<Z'<1/5. Let d have residual square
u and white cross rho from the table, and put Y=y+rho. The rank
condition is

    (x^2-1)(Y^2-tu)=u/2,
    Z=4+2u/(Y^2-tu).                              (4)

The d-wall facts and y>=0 give the following explicit bounds:

|d|upper bound for Z|
|---|---:|
|I|Z<74/9|
|H|Z<41/9<5|
|J|Z<=17/4<5|
|K|Z<=60/11|
|F|Z<316/59|
|B|Z<37/9<5|
|T|Z<72/17<5|

For instance H uses (1+k)^2-t^2=4k; J uses
2k^2-t/R=8/19; K uses (17/38)^2-2t/R=11/76; F uses
(1+e)^2-t^2=59/76. For B,T, replace Y^2 by its strict lower
bound1 and use t<1. No finite/dotted opposition premise was omitted:
the d-wall facts already proved y>1 where used.

An upper bound Z<5 makes a nonzero integer norm impossible in any
nontrivial degree. In row I, degree at least3 gives norm below
(74/9)(1/5)^2<1; degree2 has norm1 and trace5,6,7 or8, giving
quadratic fields21,2,5 or15, all outside the eight. In rows K,F,
degree at least3 is similarly excluded, and degree2 has norm1,
trace5, field Q(sqrt(21)). Degree1 is Q throughout. Thus all seven
rows are excluded.

## 7. A shared positive-complement bound for the other rows

Let the residual c square be u, its white cross with b be rho,
the residual d square v, and its white cross with c be sigma.
Write Y=y+H, where H is the b/d white cross. Because x>1 and u<1,
the a,c span has Lorentz signature(1,1). Its orthogonal complement
in the white perpendicular Lorentz3-space is positive. Project b,d
there. The residual Gram is positive semidefinite, so

    Y <= sqrt(B_x D_x)+sigma*(h*x+rho)/(x^2-u),
    B_x=t+(h^2*u+2h*rho*x+rho^2)/(x^2-u),
    D_x=v+sigma^2/(x^2-u).                         (5)

Every positive correction on the right decreases with x>1.
For example the derivative numerator for the first fraction is

    -2h*rho*x^2-2x(h^2*u+rho^2)-2h*rho*u<0;

the linear numerator fraction also has negative derivative. Thus
one may bound the right side by its value at x=1. These inequalities
use actual Lorentz geometry and finite corner data, not sampled
directions or an assumed polytope decomposition.

## 8. The five rows c=H

Here u=t,rho=k. The maximum B_x at x=1 is

    B_1=1+sqrt(2)+5/9<3.                           (6)

### 8.1 d=H

The actual neighboring corner W+b+c is whole D20. For d, its
old-D20 determinant ratio gives

    (y+1)(20y+16)/2<0

at every nonphysical embedding, hence -1<y'<-4/5. But its old-B20
ratio in(2) gives y'>-(9+sqrt(5))/19>-12/19>-4/5.
The intervals are disjoint. The eight fields are all nontrivial,
so a nonidentity embedding exists. This row is impossible.

### 8.2 d=F,B or T

Here sigma=H=rho_d and v is the d residual square. Equations(5)-(6)
give

    y<sqrt(3*(v+rho_d^2/k))+rho_d*h/k.              (7)

For F, the square inside the parentheses is19/36. Thus y<3/2.
But V=2y is integral primitive, all its other conjugates have
absolute value<3/10 by(3), and 0<V<3. Its nonzero integer norm
would be below9/10 in every nontrivial degree. This excludes F.

For B the parentheses equal1/18, and the final term is1/18.
Consequently y<1/2+1/18<1, contradicting the d-wall fact B.
For T the parentheses equal1/9 and the final term h/9<1/12;
therefore y<2/3+1/12<1, contradicting fact T.

### 8.3 d=I: the one surviving quadratic field fails a discriminant test

Here (5) gives y^2<B_1<3. The d-I fact supplies integral primitive
Z=4y^2 with

    4<Z<12,       0<Z'<2/19.

Degree at least3 has norm below12(2/19)^2<1. In degree2 the norm
is1 and trace5 through12. The only such quadratic field among the
eight is Q(sqrt(6)), at trace10. Hence

    K=Q(sqrt(6)),       Z=5+2sqrt(6).

Put V=y^2=Z/4 and w=sqrt(2)x, which belongs to K by the actual
a-b-v1-c-a cycle. The exact Schur determinant, before applying
the upper bound, is

    (V-t)w^2-2k*w-2t*V-t+2/R=0.                  (8)

Its discriminant divided by4 is

    k^2+(V-t)(2t*V+t-2/R)
       =(4419+1804sqrt(6))/1444.

The numerator norm is1065, strictly between32^2 and33^2. A square
in Q(sqrt(6)) has rational-square norm. Thus the discriminant is
not a square, equation(8) has no root w in K, and this last H row
is excluded. This is an exact obstruction, not a claim that the
formal real quadratic equation lacks real roots.

## 9. The five rows c=F

Now u=t,rho=e. Its maximum complement square has the sharp bound

    B_1=t+(t/2+2h*e+e^2)/k
       <5567/4788<7/6,                            (9)

using h<5/7. The five d rows are I,H,J,K,F.

For I, (5) gives y^2<B_1<7/6. Fact I provides primitive integral
Z=4y^2, other conjugates below2/19, and Z<14/3. Its norm in any
nontrivial degree is below(14/3)(2/19)<1. Thus I is excluded.

For H, the maximum D_x is19/36 and the final correction in(5) is
(h+e)/18. This gives Y<4/5+1/20. Since H=k>7/15, it implies
y<1/2, contrary to fact H.

For J, the maximum D_x is1/18 and the final correction is
sqrt(2)(h+e)/18. Divide the upper bound for Y by sqrt(2). It is
less than1/5+1/20=1/4<k=H/sqrt(2), forcing y<0.

For K, the maximum D_x is1/9 and the correction is(h+e)/9.
Since sqrt(7/54)<13/36 and h<5/7,

    Y<13/36+197/2394=2123/4788
                     <2142/4788=17/38=H.

Again y<0, impossible.

For F, sigma=k, v=t, and the maximum D_x is1. The correction at
x=1 is h+e, while H=e. Thus

    y<sqrt(7/6)+h<11/10+3/4=37/20.

Fact F says V=2y>2 is integral primitive and |V'|<3/10. Degree
at least3 has norm below(37/10)(3/10)^2<1. In degree2 the norm
is1 or-1 and the integer trace is2 or3. These give only rational
degeneration or fields Q(sqrt(2)),Q(sqrt(5)),Q(sqrt(13)), none in
the eight. This completes all F rows, including the opposite
D-whisker deletion F+F.

## 10. The two rows each for c=B and c=T

For both profiles u<=2/19<1/9 and rho<=1/19<1/18. At x=1 use
t<8/15,h<3/4 and 1-u>8/9 in(5). This gives

    B_1<8/15+23/144<7/10.                         (10)

For d=I, y^2<B_1<1 contradicts fact I.

For d=H, sigma=rho and H=k. The maximum D_x is less than11/20,
and the final correction is less than29/576<1/18. Consequently

    y<2/3+1/18-7/15=23/90<1,

contrary to fact H. These are all four remaining rows.

## 11. Consumer conclusion and what is still not a dimension theorem

Every one of the21 finite contact rows has been excluded under the
actual-wall arithmetic hypotheses. Both possible D-whisker-zero
placements, including the opposite F+F one, were retained in the
inventory and checked. Thus a single far-4 quarter source cannot
lose its last payment through old zeros or the approved D20
whisker deletion. Two or more far-4 sources retain at least3/2.
All other modified supports are disjoint as stated in Section1.

This proves the complete WRITTEN quarter-refund candidate at the
unique4-bond of whole B20, subject to independent source review.
At the whole B20 vertex its cost is141/4-1/4=35. The B19+A1
spectator-half candidate, if reviewed, has cost139/4 and distinct
white support. Neither result controls the remaining D19+A1,
D18+high-I2, D+D or D+E forests of cost36. No global allowance35,
new dimension ceiling, global realization or Lean verification is
claimed by this report.
