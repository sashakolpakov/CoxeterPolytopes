# Rational envelopes for the ten-direction positive-part obligation

20 October 2026. **Source-only preparation; the three final polynomial
positivity obligations have not run or been proved here.** No scientific
job, optimizer, census, Lean or build ran. The earlier moment-reduction
report is unchanged at SHA256
`8c2051b3be0f692e5097eedd81cf494f15313037d36b0f6625ca78cd01aabd34`.
This note proves rational upper envelopes and identifies exactly what a
small coefficient checker must establish. It does not assert a new
dimension bound.

## 1. Domain and moment envelope

Use the exact ten-direction notation from
`right_angled_mixed_k10_subset_moment_reduction_20261020.md`:

    v in [(3-2sqrt5)/11,(3+2sqrt5)/11],
    mu=(11+3v)/10,
    V=(1+6v-11v^2)/10,
    w=sqrt(21V/10), a=mu-w,
    gamma=8/sqrt90, |T|<=gamma*V^(3/2).

The actual `v` interval lies strictly inside `[-7/50,7/10]`.
All denominators used below are positive there. The target positive-part
average is `261/500`, equivalently the corrected reciprocal average
must be at most `223/250`.

The exact remainder coefficient is

    K=(4a+mu)/(mu^5*a^2).

The joint cubic/quartic moment argument gives the upper envelope

    average r^(-2)
      <=mu^(-2)+7V/(10mu^4)
        +max_(|T|<=gamma V^(3/2))
          [7T/(15mu^5)+K*(43V^2/300-T^2/(15V))].       (1)

At `V=0` the bracket vanishes. Its expression multiplying `K` is
nonnegative throughout the actual moment domain, because

    43/300-gamma^2/15=259/2700>0.                      (2)

Thus increasing `K` increases the right side before taking the maximum.
It would be incorrect to increase its positive and negative occurrences
independently without (2).

## 2. The universal rational bound K<7

The full-frame ratio bound gives

    a>11/20,

since `139^2>3*80^2`. If `v>=0`, `mu>=11/10`, and hence

    K<=12000000/1771561<7.                            (3)

If `v<=0`, then `V<=1/10` and `V'>=3/5`. For `V>0`,
`a'=3/10-(21/20)V'/sqrt(21V/10)<0`; the endpoint with
`V=0` follows by continuity. Therefore

    a>=a(0)=(11-sqrt21)/10>16/25.

Also `v>-1/6` implies `mu>21/20`. Consequently

    K<7220000000/1045529856<7.                        (4)

All these comparisons are rational squarings or cross multiplications.
We may replace `K` by seven in (1), with its whole bracket retained.

## 3. Negative outside coordinate

For `v<=0`, `V<=1/10` and `mu<=11/10`. The vertex of the concave
quadratic in (1), after replacing `K` by seven, is `V/(2mu^5)`.
It is beyond `L=gamma V^(3/2)`, because

    2gamma mu^5 sqrtV
        <=(8/15)*(11/10)^5<1.

The maximum is therefore attained at `T=L`. Use
`gamma<17/20` and the tangent inequality

    sqrtV<=1/10+(5/2)V

to obtain the rational upper envelope

    Jneg(v)=mu^(-2)+7V/(10mu^4)
         +(119V/3000+119V^2/120)/mu^5
         +(1813/2700)V^2.                            (5)

The positive-part correction is zero here:
`mu+w<39/25<41/25<10/sqrt37`. Multiplying
`223/250-Jneg` by the positive `mu^5` gives the degree-nine
polynomial

    Pneg=(223/250)mu^5-mu^3-(7/10)Vmu
          -(119/3000)V-(119/120)V^2
          -(1813/2700)V^2 mu^5.                     (6)

It is sufficient to prove `Pneg>=0` on the larger rational interval
`[-7/50,0]`. The slightly larger interval includes points with `V<0`;
only polynomial positivity is being checked there, not validity of
the original moment interpretation at those extraneous points.

## 4. Nonnegative outside coordinate

For `v>=0`, maximize the same concave quadratic over all real `T`.
Completing the square gives the rational bound

    Jpos(v)=mu^(-2)+7V/(10mu^4)
                 +(301/300)V^2+7V/(60mu^10).        (7)

For `0<=v<=1/20`, `mu<=223/200`, `V<=509/4000` and
`w<13/25`. Thus `mu+w<327/200<41/25`, so again the correction
is zero. Clearing the positive `mu^10` gives

    Pzero=(223/250)mu^10-mu^8-(7/10)Vmu^6
              -(301/300)V^2 mu^10-(7/60)V.          (8)

It is sufficient to check this degree-fourteen polynomial on
`[0,1/20]`.

For the rest of the actual domain, use the correction already proved
in the moment-reduction note. Since

    (37/100)^(3/2)<113/500,
    10/sqrt37>41/25,

it is at most

    Cpos<=791V/[30000*(41/25-mu)].                    (9)

The denominator is positive even on the larger interval
`[1/20,7/10]`. After clearing it and `mu^10`, the remaining
degree-fifteen polynomial is

    Pcor=(41/25-mu)*Pzero-(791/30000)Vmu^10.         (10)

Thus positivity of `Pcor` on `[1/20,7/10]` suffices. Again,
extraneous points beyond the actual domain are only polynomial test
points; their moment interpretation is not asserted.

## 5. Exact finite verification contract

The prepared independent standard-library checker constructs (6), (8)
and (10) by exact `Fraction` polynomial arithmetic, then computes
every Bernstein coefficient on fixed subintervals of length `1/1000`.
Each full subinterval must have strictly positive coefficients. This
would prove strict positivity of all three polynomials on their full
closed domains, without a numerical optimization, a sampled inequality
or a root finder. It also checks all rational constants used in
(2)--(4) and the sign/endpoint comparisons above.

The worker is source-only until the parent has reviewed it and runs
generation and a separate replay through the shared resource wrapper.
It owns no inner global lock; it has a 45-second child deadline and
prints only a compact summary. Generation uses exclusive creation;
replay compares canonical JSON and writes nothing. The certificate
records the polynomial coefficients, exact global minimum Bernstein
coefficient and location, a digest over every coefficient, and a
small positive rational floor for every subinterval. The geometry and
the reduction to these polynomials are written premises, not claims
established by the finite coefficient check.
