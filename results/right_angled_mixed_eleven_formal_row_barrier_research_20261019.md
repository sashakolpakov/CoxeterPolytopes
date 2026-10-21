# Eleven-dimensional moments and a barrier to the same scalar deletion method

19 October 2026. Paper-only research. No scientific calculation, optimizer,
graph census, Lean or build ran. This note does not weaken the separate
dimension12 exclusion and supplies no geometric example in dimension11.
It identifies a precise limitation of a direct all-formal-row successor.

## Fresh dimension11 identities

For the actual P3 opposition graph in dimension11 the invariant
h=3q-sum n_i has square2. Choosing a future null Z perpendicular to
h on each opposition component again gives sum p_i=3, sum p_i^2=1.
Put x_i=p_i-3/11. Then

    sum x_i=0,      sum x_i^2=2/11,
    r_S=13/11+d_S, d_S=-sum_(i in S) x_i,
    |d_S|<=4sqrt3/11,

for the165 formal triples. The one-, two- and three-index inclusion
probabilities are3/11,3/55,1/165. Repeating the symmetric expansion
before imposing sum x=0 gives

    average d=0,
    average d^2=24/605,
    average d^3=-(4/33) sum x_i^3,
    average d^4=28/6655-(2/55) sum x_i^4 <=28/6655.

The elementary square-sum bound further gives sum x_i^4>=4/1331,
so the fourth upper bound improves to60/14641. The basic third
bound is -average d^3<=8sqrt22/3993. These are pointwise formal-row
identities, not corner counts or height-maximum restrictions.

## Exact meaning of the limited method

Fix a power a>0 and deletion constant c>0. Suppose one seeks a
UNIFORM pointwise bound M on

    sum_all165 max(0,r_S^(-a)-1+c)

over every formal frame sum p=3,sum p^2=1. A nonnegative polynomial
majorant plus moments is one way to obtain such M. Undirected
reciprocal-pair cancellation would then give 2B/f<=M/c.

The following two exact formal frames force

    M/c >1455/17=85.588235...,                     (1)

for EVERY such power and deletion constant. This already exceeds
the current strict incidence lower bound117825/1379=85.44235...;
their cross-product difference is3420. Thus this exact uniform
scalar route cannot exclude dimension11 using that lower bound,
even if all its interval and moment estimates were made sharp.

This statement does not cover graph-dependent row bounds, additional
geometric exclusions, an improved incidence lower bound, or a proof
that chooses Z in a controlled componentwise way to avoid these
formal rows. It is not an obstruction to dimension at most10.

## Two admissible formal frames

Frame A has two coordinates6/11 and nine coordinates7/33. Their
sum is3 and squared sum is1. Its formal ratios are exactly

|Triple type|Ratio r|Multiplicity|
|---|---|---:|
|two large coordinates|23/33|9|
|one large coordinate|34/33|72|
|no large coordinate|15/11|84|

Frame B has one coordinate(3-2sqrt5)/11 and ten coordinates
(15+sqrt5)/55. Their deviations from3/11 are -2sqrt5/11 and
sqrt5/55 respectively, hence again sum p=3 and squared sum1.
Its ratios are

|Triple type|Ratio r|Multiplicity|
|---|---|---:|
|contains the distinguished coordinate|(65+8sqrt5)/55|45|
|does not contain it|(65-3sqrt5)/55|120|

These are valid points of the formal Lorentz frame sphere. No
actual bad-triangle graph, chamber or opposition component is
asserted to realize any prescribed subset of their165 triples.

## Uniform power/deletion obstruction

Define

    L_A(a)=9*((33/23)^a-1)+72*((33/34)^a-1),
    L_B(a)=120*((55/(65-3sqrt5))^a-1)<0.

Since a positive part is at least its argument, retaining the
first81 terms of Frame A and the last120 terms of Frame B gives

    M/c >= max(81+L_A(a)/c, 120+L_B(a)/c).          (2)

It does not matter if a retained term's argument is negative;
this remains a legitimate lower bound for its positive part.

Let their logarithmic derivatives at a=0 be

    L_A=9 log(33/23)-72 log(34/33),
    L_B=-120 log((65-3sqrt5)/55).

The inequalities log((1+t)/(1-t))>2t and log(1+t)<t imply

    L_A >9*(5/14)-72*(1/33)=159/154>1.

Also 1<(65-3sqrt5)/55 and sqrt5>11/5, so

    0<-L_B<120*((65-3sqrt5)/55-1)
             <408/55<15/2.

From exp(a u)>=1+a u for every real u we get L_A(a)>=a L_A>a.
From 1-exp(-a u)<=a u for u>0 we get -L_B(a)<=a(-L_B)<15a/2.
Thus L_A(a)>0 and

    (-L_B(a))/L_A(a)<15/2.                        (3)

For fixed positive A=L_A(a) and negative B=L_B(a), the two
functions81+A/c and120+B/c in(2) respectively decrease and
increase with c>0. Their maximum has minimum at
c=(A-B)/39, with value

    81+39*A/(A-B)=81+39/(1+(-B)/A)
       >81+39/(1+15/2)=1455/17.

This proves(1) for every positive exponent, without taking a
limit or numerically optimizing a parameter. The same argument
works for the logarithmic transport identity itself, whose
reverse-edge logarithms cancel exactly.

## Research implication

The dimension12 positive-part proof gains an ORDINARY global bound
from actual reciprocal pairs; its mechanism is valid. At dimension11,
however, the formal row optimization alone already has the explicit
barrier above. The missing leverage must enter through stronger
ordinary incidence density, actual graph restrictions among the
selected triples, a better controlled potential, or another genuinely
geometric input. Maximum-height clique bounds remain useful local
information but cannot replace the ordinary degree average.
