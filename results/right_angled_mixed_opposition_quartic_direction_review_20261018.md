# A bounded direction polynomial improves the mixed reciprocal-square cap

This is a new analytic consequence of the frozen
[power-two packet](right_angled_mixed_opposition_power_two_research_20261018.md).
It does not alter that packet. It proves the stronger universal bound

    2B/f <= (5511760+23040sqrt(33))/34375
          <1128848/6875, approximately164.1961.

The sharper irrational expression is approximately164.1924; the previous
cap was165.43744. Here B counts canonical P3 opposition edges and f finite
vertices in a finite-volume right-angled twelve-polytope. The exact target
6620/59 is still smaller. No new dimension exclusion follows.

The independent-sign balancing lemma from the ideal-six cusp packet does
not directly apply: these are uniform three-subset sums with a zero-sum
constraint on twelve real coordinates. The following quartic majorant
preserves both that mean and the squared norm.

## Reduction to one direction polynomial

Use exactly the centered variables in the frozen moment proof:

    sum_(i=1)^12 x_i=0,  sum x_i^2=1/4,
    d_S=-sum_(i in S)x_i,  r_S=mu+d_S,  mu=5/4.

Average over the220 three-subsets. The proved identities are

    E d=0, E d^2=9/176,
    E d^3=-(27/220)sum x_i^3,
    E d^4=3/440-(3/220)sum x_i^4.

Keep the frozen pointwise reciprocal-square majorant, including its
remainder multiplier13. Before separately bounding the cubic and quartic
power sums, it gives

    E r^(-2) <= A0 + [108sum x_i^3-39sum x_i^4]/(220mu^5),
    A0=1/mu^2+3(9/176)/mu^4+39/(440mu^5).                       (1)

Thus it suffices to bound the single polynomial sum
P(x)=sum phi(x_i), where phi(t)=108t^3-39t^4. Combining the two power
sums before taking an upper bound preserves their correlation.

## A global quadratic majorant

For each coordinate the zero-sum constraint gives

    x_i^2=(sum_(j!=i)x_j)^2 <=11(1/4-x_i^2).

Therefore every x_i<=a, where

    a=sqrt(33)/12,  b=-a/11,  L=108-39a-78b.

For all real t<=a,

    (t-b)^2(a-t)(L-39t) >=0.                                   (2)

Indeed a-t>=0, and L-39t>=L-39a=108-78(a+b)>0. The last strict
inequality follows already from a<1/2 and a+b=10a/11.

The left side of (2) has leading terms39t^4-108t^3. Hence it equals

    -phi(t)+A t^2+D t+C

for fixed real constants A,D,C. Summing (2), using sum x_i=0 and
sum x_i^2=1/4, gives P(x)<=A/4+12C. The feasible vector
(a,b,...,b) has zero sum, squared norm1/4, and makes every term in (2)
zero. Thus this is the exact global maximum of P over the constrained
sphere, with no stationary-point sampling or optimization:

    max P=108(a^3+11b^3)-39(a^4+11b^4)
         =45sqrt(33)/22-1443/704.                               (3)

Here a^3+11b^3=5sqrt(33)/264 and a^4+11b^4=37/704.

## Improved cap and its limitation

Insert (3) into (1) and multiply by220, as in the frozen positive-function
opposition proof. The resulting cap is

    220/mu^2+660(9/176)/mu^4
       +39/(2mu^5)+(45sqrt(33)/22-1443/704)/mu^5
      =(5511760+23040sqrt(33))/34375.

Using sqrt(33)<23/4 gives the displayed rational upper1128848/6875.
It is strictly smaller than the old516992/3125 cap. Only the directional
moment estimate changed; the actual-opposition graph, invariant, positivity
and reciprocal-edge summation are inherited exactly.

This still cannot close the target112.20339 by an unrestricted full220
reciprocal-square average alone. Jensen gives, for every admissible frame,

    E r^(-2) >=(E r)^(-2)=16/25.

Consequently every universal upper obtained solely as220 times a bound
for that full average is at least704/5=140.8, already above6620/59.
Sharper all-order control of the same full average therefore cannot
reach that particular target. Information about which triples are actual,
or another geometric inequality, is needed for this route.

This note is a written mathematical derivation. No computation, optimizer,
original-packet edit, Lean compilation or new example was used for it.
The source moment and pointwise identities are the explicitly inherited
inputs linked at the start. A separate arithmetic checker is not claimed.
