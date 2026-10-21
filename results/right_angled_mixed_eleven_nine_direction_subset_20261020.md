# The nine-direction structured opposition estimate

20 October 2026. **Analytic proof with a finite rational interval step.**
This concerns the same structured subset contract as
`right_angled_mixed_eleven_structured_small_and_full_subsets_20261020.md`.
It does not assert a new dimension bound. No Lean or polytope enumeration
is used.

Let sum p_i=3, sum p_i^2=1 over eleven real coordinates, and let I
contain nine of them. Write r_T=2-sum_T p_i. We prove

    sum_(T subset I, |T|=3) (r_T^-2-37/100)_+ <=1323/25.

The right side is756/25+(27/100)*84. The common full-frame range
11/20<r_T<20/11 and the centered moment identities are derived in
the preceding source. Here is the complete remaining argument.

## One parameter and the reciprocal envelope

Put delta=3-sum_I p_i. Cauchy on the two omitted coordinates gives

    0<=delta<=12/11,
    mu=1+delta/3,
    0<=V=sum_I(p_i-sum_I(p_i)/9)^2
       <=v(delta)=2delta/3-11delta^2/18<=2/11.

The interval for delta is exact: the possible nine-coordinate sum is
between21/11 and3. Center within I, and write T3 for its cubic
power sum. The centered ratio moments are

    E d^2=V/4, E d^3=-3T3/28,
    E d^4=5V^2/28-3T4/28.

Projection of the vector of squared deviations onto the constant and
deviation vectors gives T4>=V^2/9+T3^2/V when V>0. Also
|T3|<=7V^(3/2)/sqrt72<(5/6)V^(3/2), by the elementary sharp
zero-sum cubic bound proved in the preceding full-subset argument.

Take any lower bounds mu0<=mu and 0<a0<=min_T r_T. The reciprocal
remainder factor is bounded above by

    K0=(4a0+mu0)/(mu0^5*a0^2).

The exact Taylor identity in the preceding source, nonnegativity of
E d^4 and the projection inequality imply, for any upper bound v0>=V,

    E r^-2 <=mu0^-2+(3/4)v0*mu0^-4+(K0/6)v0^2
                 +max_[0<=z<=(5/6)v0 sqrt(v0)]
                      [(3/(7mu0^5))*z-(3K0/(28v0))*z^2].       (1)

To justify the replacements together, first bound the nonnegative
fourth-moment remainder by K0 E d^4, and only then insert its joint
upper bound. The expression v^2/6-3T3^2/(28v) increases with v.
The cubic coefficient is positive; a negative T3 gives no larger
maximum than z=0. Thus the negative quadratic term is retained with
the SAME K0. At V=0 the claim follows by continuity or directly.
The maximum in (1) is at the smaller of its right endpoint and its
positive quadratic vertex, so no numerical optimization is required.

## The equality endpoint and a full neighborhood

At delta=0, Cauchy forces V=0: the nine coordinates are1/3 and the
other two are0. All84 ratios are1, giving equality1323/25 exactly.

For 0<delta<=1/100 we have mu>=1 and

    min r>=mu-sqrt(2V)>=1-sqrt(1/75)>7/8.

The Taylor remainder is therefore at most K=288/49<6. The exact
inequality (1+x)^-2<=1-2x+3x^2 for x>=0 follows on multiplying
by (1+x)^2; its difference is x^3(4+3x). Discarding the favorable
negative cubic/quartic coupling gives

    E r^-2 <=1-2delta/3+delta^2/3+delta/2
                       +(delta^(3/2)/5)+(4/9)delta^2
              <=1-delta/6+delta/50+(7/9)delta^2
              <=1-(5/36)delta <1.                        (2)

For the cubic coefficient, the bound (5/6)V^(3/2) above gives a
coefficient whose square is50/1323<1/25. The quartic coefficient
uses E d^4<=V^2/6 and V<=2delta/3. The last line of (2) uses
delta<=1/100. In this neighborhood

    max r<=301/300+sqrt(1/75)<18/11<10/sqrt37.

There is no positive-part correction: E(r^-2-37/100)_+=E r^-2-37/100.
This proves the full neighborhood, not just a finite list near zero.

## The rest of the interval

The exact checker covers [1/100,12/11] by consecutive closed
intervals of length1/10000, truncating the last interval exactly at
12/11. On each it takes the exact maximum v0 of the concave quadratic
v(delta), uses mu0 at the left endpoint, and sets

    a0=max(11/20, mu0-sqrt_up(2v0)).

Square roots are rounded upwards at denominator10^6 by integer-square
comparisons. This gives a valid positive lower ratio bound throughout
the interval. Formula (1) is then rational, including the clipped
quadratic maximum.

To pass to positive parts, add a correction to E r^-2-37/100.
The global correction is at most27/400. If the upper mean mu1 is
below18/11, the correction is also at most

    (113/500)*(v0/4)/[2*(18/11-mu1)].

If mu1+sqrt_up(2v0)<=18/11, it is zero. These are the range and
second-moment corrections proved in the preceding source. The checker
takes their minimum and proves that (1) plus this correction is
strictly less than1 on every interval. This is exactly the required
positive-part upper average63/100.

The equality endpoint, neighborhood (2), and this complete interval
cover prove the nine-direction contract. The certificate must be read
with the written moment and correction arguments, not as a classification
of polytopes. The separate ten-direction case and the final actual-graph
transport composition remain obligations until independently checked.
