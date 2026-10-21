# Ten-direction positive-part rows: exact reduction and a symmetry caution

20 October 2026. **Paper-only analytic reduction; the final interval
inequality is open.** No calculation program, optimizer, census, script
replay, Lean or build ran. This note supplies a precise one-variable
obligation for the proposed eleven-dimensional subset estimate. It
does not assert that obligation or a dimension-ten theorem.

## 1. The proposed contract

Let eleven real coordinates satisfy `sum p=3`, `sum p^2=1`, and put
`r_S=2-sum_S p` for a triple. For a ten-element set `I`, the proposed
analytic estimate is

    sum_(S subset I, |S|=3) max(0,r_S^(-2)-37/100)
          <=756/25+(27/100)*120=1566/25.       (1)

There are 120 triples, so the target average is `261/500`.
The corresponding candidate for general `k=|I|` replaces 120 by
`C(k,3)` on the right. If it holds for all ideal-direction subsets at
every finite vertex, then actual-subset deletion with `c=9/25`
and `H(q)=C(k,3)-t(q)` gives

    sum_actual [r^(-2)-1+c]<=84c+(27/100)*H(q).

The reviewed `sum H<=2E` and ordinary reciprocal-edge pairing would
then yield `2B<=84f+(3/2)E`. This is incompatible with the separately
proved candidate balanced lower bound
`2B>(159195/1861)f+(2794/1861)E`, since its intercept exceeds84
and `2794/1861-3/2=5/3722>0`. None of these implications proves
the missing analytic contract (1).

## 2. The exact one-variable domain

Write `v` for the one coordinate outside `I`. For the ten inside
coordinates let `S=sum_I p`, `N=sum_I p^2`, and set

    S=3-v,         N=1-v^2,
    V=N-S^2/10=(1+6v-11v^2)/10,
    mu=2-3S/10=(11+3v)/10.

The domain is exactly

    (3-2sqrt5)/11 <=v<=(3+2sqrt5)/11.          (2)

Indeed the centered squared norm `V` must be nonnegative, and any
ten-vector with this mean and centered norm gives an admissible
eleven-frame. Let `y_i=p_i-S/10` on `I`; then `sum y=0` and
`sum y^2=V`. For a triple put `d=-sum_S y`, so `r=mu+d`.
Zero-sum Cauchy gives

    |d|<=w=sqrt(21V/10),       a=mu-w>0.       (3)

Positivity also follows from the full eleven-frame bound
`a>=(13-4sqrt3)/11`; no formal triple is assumed to be an actual
chamber vertex.

The ten-index inclusion probabilities are `3/10,1/15,1/120`.
Expansion before imposing `sum y=0` gives

    average d=0,
    average d^2=7V/30,
    average d^3=-(7/60)*T,
    average d^4=(3/20)V^2-U/15,                (4)

where `T=sum y_i^3`, `U=sum y_i^4`. These formulas are pointwise
at one formal frame, not averages over geometric vertices.

## 3. Joint cubic and quartic information

For any zero-sum `m`-vector of squared norm `V`,

    |T|<=(m-2)*V^(3/2)/sqrt(m(m-1)),
    U>=V^2/m+T^2/V             when V>0.       (5)

Both have short direct proofs. Put
`M=sqrt((m-1)V/m)` and `b=-M/(m-1)`. Zero-sum Cauchy gives
`y_i<=M`. Summing the nonnegative polynomials
`(M-y_i)(y_i-b)^2` and using `sum y=0` bounds `T` by
`(m-2)V^(3/2)/sqrt(m(m-1))`; replacing `y` by `-y` bounds
the opposite sign. For the second inequality apply Cauchy to

    T=sum_i y_i*(y_i^2-V/m).

For `m=10`, these specialize to

    |T|<=L=8V^(3/2)/sqrt90,
    U>=V^2/10+T^2/V.                          (6)

Using separate worst-case bounds on `T` and `U` loses useful
information: a nearly maximal positive third moment forces a large
fourth moment, whose coefficient in (4) is negative.

## 4. A rigorous reciprocal envelope

For `V>0`, set

    K=(4a+mu)/(mu^5*a^2),
    A=7/(15mu^5),        C=K/(15V),
    Tstar=min(L,A/(2C)).

The exact reciprocal-square remainder identity and the fact that
`(4r+mu)/r^2` decreases give, using (4)--(6),

    average r^(-2) <= J(v),
    J(v)=mu^(-2)+7V/(10mu^4)+(43K/300)V^2
                   +A*Tstar-C*Tstar^2.       (7)

This is an actual proved upper envelope. The last term is the maximum
of a concave quadratic on `[-L,L]`; `A>0` puts its maximum at the
displayed nonnegative `Tstar`. At `V=0` define `J=mu^(-2)`.
There is no assumption of a uniform inside vector or a classification
of all critical frames.

To verify (7), the Taylor/moment part before applying (6) is

    mu^(-2)+7V/(10mu^4)+7T/(15mu^5)
                   +K*((3/20)V^2-U/15).

Substitution of `U>=V^2/10+T^2/V` produces the coefficient
`3/20-1/150=43/300`, exactly as in (7).

## 5. Correcting the positive part

Put `lambda=37/100` and `rho=1/sqrt(lambda)=10/sqrt37`.
Throughout (2), `mu<rho`. Write

    Cpos(v)=min(
      max(0,lambda-(mu+w)^(-2)),
      lambda^(3/2)*(7V/30)/(2*(rho-mu)) ).     (8)

Then

    average max(0,r^(-2)-lambda)
        <=J(v)-lambda+Cpos(v).                (9)

The correction is the average of `(lambda-r^(-2))_+`. Its pointwise
maximum on (3) gives the first term of (8). For the second, concavity
on `r>=rho` gives

    (lambda-r^(-2))_+<=2*rho^(-3)*(r-rho)_+.

With `b=rho-mu>0`, the elementary inequality
`(d-b)_+<=d^2/(4b)` and (4) give the second term of (8). Thus
if `mu+w<=rho` the correction is exactly zero, rather than being
unnecessarily bounded by a positive moment allowance.

A sufficient remaining one-variable statement for (1) is therefore

    J(v)+Cpos(v)<=223/250      on (2).         (10)

The domains, square roots, minimum branches and endpoint definitions
are explicit. Equation (10) has **not** been proved or checked here.
A rigorous interval certificate for it would prove this ten-direction
part of the proposed contract; it would not by itself settle the
other subset sizes or universal geometry.

## 6. Uniform inside coordinates are not a global maximum

The endpoint `v=(3-2sqrt5)/11` has `V=0`, giving the equal inside
coordinates of Frame B. Its inside reciprocal sum is

    120*[55/(65-3sqrt5)]^2.

It is a local maximum in the first nontrivial directions but not
a global maximum. Consider instead the exact admissible frame

    v=-1/10,
    one inside coordinate =31/100+3sqrt29/100,
    nine inside coordinates=31/100-sqrt29/300. (11)

The inside sum is `31/10`; its centered squared norm is `29/1000`,
and its total squared norm is `99/100`. Adding `v` verifies the
two full-frame equations exactly. The two ratio types are

    (321-7sqrt29)/300       with multiplicity36,
    (107+sqrt29)/100        with multiplicity84.

All corresponding reciprocal weights exceed `lambda`. Its reciprocal
sum is strictly greater than `106.8696`, whereas the equal-inside
endpoint sum is less than `106.848`. These comparisons follow by
direct rational squaring from

    5.3851<sqrt29<5.3852,       sqrt5<2.237:

the two reciprocal squares in (11) are respectively greater than
`1.1213` and `.7917`, while the endpoint reciprocal square is less
than `.8904`. Decimal numbers in this paragraph are terminating
rationals, not sampled computer estimates. The proposed contract
allows the larger reciprocal sum `107.04`, so (11) does not refute
it. It does refute a proof that silently substitutes the uniform
local maximum for the global one.

## 7. The small subset sizes require no moment optimization

The proposed contract for `0<=k<=6` follows by elementary bounds.
For every full eleven-frame and every triple, zero-sum Cauchy gives

    r>6/11,

because the exact lower endpoint is `(13-4sqrt3)/11` and
`sqrt3<7/4`. Therefore each positive-part summand is strictly less
than `121/36-37/100=673/225`. For `k<=5` there are at most ten
triples, so their sum is less than `6730/225<756/25`; this is
already less than the required right side, including `k=0,1,2`
where there are no triples. No sign of an individual `p_i` is used.

For `k=6`, pair each of its twenty triples with its complementary
triple. The ratios in such a pair sum to

    C=4-sum_I p>18/11.

Indeed the maximum possible six-coordinate sum in a full frame is
`(18+2sqrt15)/11<26/11`, by `sqrt15<4`. Put
`f(r)=max(0,r^(-2)-37/100)`. It is convex and decreasing for `r>0`.
Both paired ratios exceed `a=6/11`, so convexity makes their sum
of weights at most `f(a)+f(C-a)`, and monotonicity bounds this by

    f(6/11)+f(12/11)=12461/3600.

There are ten complementary pairs, giving

    sum_triples f(r)<12461/360<891/25
                   =756/25+(27/100)*20.

Strictness follows from the strict rational lower bounds on the
ratio and paired sum. This proves the complete `k<=6` part of the
proposed contract. It leaves sizes seven through eleven; size ten
is reduced to (10) above, not proved here.

## 8. Common moment and remainder formulas for other subset sizes

For any `k>=3`, center its selected coordinates at `S/k`, put
`V=sum y^2`, `mu=2-3S/k`, and define `T,U` as above. The one-,
two- and three-index inclusion probabilities in a uniform triple are

    P1=3/k,
    P2=6/[k(k-1)],
    P3=6/[k(k-1)(k-2)].

Expanding the subset powers, using `sum y=0`, gives

    average d^2=alpha2*V,
    average d^3=-alpha3*T,
    average d^4=alpha4*U+b4*V^2,

where

    alpha2=P1-P2=3(k-3)/[k(k-1)],
    alpha3=P1-3P2+2P3
          =3(k-3)(k-6)/[k(k-1)(k-2)],
    alpha4=P1-7P2+12P3
          =3(k^2-17k+54)/[k(k-1)(k-2)],
    b4=3P2-6P3=18(k-4)/[k(k-1)(k-2)].

For `7<=k<=11`, `alpha3>0` and `alpha4<0`. Let

    a=mu-sqrt(3(k-3)V/k)>0,
    K=(4a+mu)/(mu^5*a^2),
    L=(k-2)V^(3/2)/sqrt(k(k-1)).

The exact identity used in Section 4 is, with `r=mu+d`,

    r^(-2)=mu^(-2)-2d/mu^3+3d^2/mu^4-4d^3/mu^5
                    +d^4*(4r+mu)/(mu^5*r^2).

It follows by multiplying through by `mu^5*r^2`. The varying factor
`(4r+mu)/r^2=4/r+mu/r^2` has derivative
`-4/r^2-2mu/r^3<0`. Thus its value at the rigorous lower endpoint
`a` bounds the entire fourth-order remainder. The general counterpart
of (7) is consequently

    average r^(-2)
      <=mu^(-2)+3alpha2*V/mu^4
           +K*(b4+alpha4/k)*V^2
           +max_(|T|<=L) [4alpha3*T/mu^5+K*alpha4*T^2/V].

For `V>0` this is again an explicit concave-quadratic maximum;
for `V=0` it is `mu^(-2)`. The positive-part correction follows
Section 5 with `average d^2=alpha2*V` whenever `mu<rho`.
No assertion that a formal endpoint or a moment-extremizing vector
is a chamber object is needed or made.
