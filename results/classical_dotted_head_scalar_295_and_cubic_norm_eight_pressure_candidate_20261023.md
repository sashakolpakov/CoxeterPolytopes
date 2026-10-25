# Two exact estimates for the remaining classical terminal budget

23 October2026. **Complete SOURCE-ONLY mathematical CANDIDATE, pending
full different-author review and root acceptance.** Root authored both
proofs below. No scientific worker/import, enumeration, numerical or
symbolic process, optimizer, test, compiler, Lean or TeX job supplied
them. Every terminating decimal is an exact rational.

These are scalar and algebraic-number estimates, not a new geometric
dimension theorem. In particular the first lemma concerns ONLY the
displayed dotted/cleared-head formula on0<u<1. It does not strengthen
the different finite-own formula at u>1. The second lemma assumes an
element of actual degree3 with its three stated conjugates; it makes
no assumption that dividing that element by2 remains integral.

## 1. The exact dotted-head scalar is greater than .295

For0<u<1 define

    H(u)=[sqrt(2-u)-sqrt(u(1-u))]/sqrt2,
    R(u)=1-u+u^2/2,
    F(u)=(4/5)H(u)-(4/25)R(u).

**Lemma.** F(u)>.295 on the ENTIRE interval0<u<1.

This is the literal dotted row of
`results/compact_whole_b20_two_end_actual_foot_exclusion_20261021.md`,
SHA256 `5b6040b924e3895f92024600616bafaae0c760abd93afb76e6437b3dbe84eff3`.
Its earlier .2893 seal is not a premise of this proof. Any application
still has to prove that its ACTUAL feasible-point pressure has exactly
this H and R, and must pay the original reference coefficients.

### A global radical bound, proved by a positive polynomial

Put c=.681 and t=2-c^2=1.536239. We first prove H(u)>c/sqrt2.
The difference D=sqrt(2-u)-sqrt(u(1-u)) is positive, since
(2-u)-u(1-u)=2-2u+u^2>0. Also t-u^2>.536239>0. Thus it suffices
to show

    P(u)=(t-u^2)^2-4u(1-u)(2-u)>0.

Indeed this gives t-u^2>2sqrt(u(1-u)(2-u)), hence D^2>c^2.
Set r=u-.7. Exact expansion gives

    P(u)=r^2(r-.6)^2+3.107522r^2-.0094692r+.002616045121.

The last three terms are strictly greater than
3r^2-.01|r|+.0026. The elementary square
(|r|-.005)^2>=0 gives .01|r|<=r^2+.000025. Consequently

    P(u)>r^2(r-.6)^2+2r^2+.002575>0.

This proves the required radical bound on the whole interval, without
sampling or solving for its minimum. Since1/sqrt2>.7071, it also gives
H(u)>.4815351 everywhere.

### Monotonicity and complete rational interval cover

R decreases on(0,1). H decreases on(0,.65]: for u<=.5 its derivative
is immediately negative. On[.5,.65], the possibly positive derivative
part satisfies

    (2u-1)/sqrt(u(1-u))<=.3/sqrt(.2275)
                            <1/sqrt(1.5)<=1/sqrt(2-u).

The strict middle inequality is exactly .09*1.5<.2275. This proves
the derivative claim without assuming where the true minimum lies.

The following endpoint bounds are justified entirely by positive
rational squares. The L and U columns have L^2<2-v and U^2>=v(1-v).
Moreover sqrt2<1.4143, and L-U exceeds1.4143 times the last column.

|v|L|U|Strict lower bound for H(v)|
|---:|---:|---:|---:|
|.25|1.3228|.4331|.62|
|.4|1.2649|.4899|.54|
|.5|1.2247|.5|.51|
|.55|1.2041|.4975|.499|
|.6|1.1832|.4899|.4902|
|.625|1.1726|.48413|.4867|
|.65|1.16189|.47697|.4842|

These give the COMPLETE interval cover below. In the first seven
rows use H at the right endpoint and R at the left endpoint. In
the final row use the global radical bound and R(.65)=.56125.

|u interval|Strict lower bound for F(u)|
|---|---:|
|(0,.25]|.336|
|[.25,.4]|.307|
|[.4,.5]|.2992|
|[.5,.55]|.2992|
|[.55,.6]|.29596|
|[.6,.625]|.29656|
|[.625,.65]|.29611|
|[.65,1)|.29542808|

Every endpoint is covered and every displayed lower bound is>.295.
This proves the lemma on the unbounded dotted-own range through its
complete reciprocal parameter0<u<1. No finite-own branch is included.

## 2. A cubic of norm8 with two small positive conjugates is greater than42

**Lemma.** Let A be a degree3 algebraic integer, with ALL three
conjugates A,a,b positive,0<a,b<1/2, and Aab=8. Then A>42.

The ordinary algebraic-integer trace and elementary symmetric sum are
integers T=A+a+b and B=A(a+b)+ab. The monic minimal polynomial is

    p(X)=X^3-T X^2+B X-8.

No maximal-order formula, PID, field catalogue, unramifiedness or
integrality of A/2 is assumed.

Suppose A<=42. The norm gives A>32 and T in{33,...,42}. Since
a+b<1 and ab=8/A, AMGM and the strict trace interval give

    B>=2sqrt(8A)+8/A
      >2sqrt(8(T-1))+8/T.                       (1)

Because A>1/2 and0<a,b<1/2, p(1/2)<0. Direct evaluation gives

    B<T/2+63/4.                                (2)

We now exhaust ALL ten integer traces, with exact integer seals.
For the FIVE even traces, (2) gives the following Bmax; the positive
last column says already Bmax<2sqrt(8(T-1)), contradicting(1).

|T|Bmax|32(T-1)-Bmax^2|
|---:|---:|---:|
|34|32|32|
|36|33|31|
|38|34|28|
|40|35|23|
|42|36|16|

At T33, Bmax32, while (1) is strictly greater than32. For the next
THREE odd traces the complete rational-square seals are:

|T|Bmax|32(T-1)T^2-(T*Bmax-8)^2|
|---:|---:|---:|
|35|33|17191|
|37|34|14588|
|39|35|8087|

Here T*Bmax-8>0. Each positive integer therefore gives
Bmax-8/T<2sqrt(8(T-1)), again contradicting(1).

Only T41 remains. Equation(2) gives B<=36, while (1) gives
B>2sqrt320>35, so B36. But the resulting polynomial

    X^3-41X^2+36X-8

has discriminant

    41^2*36^2-4*36^3-32*41^3-27*8^2+144*41*36=-2704.

A negative cubic discriminant is incompatible with its three real
conjugates A,a,b. This excludes the final trace. Thus A>42, as claimed.

## 3. Intended use and exact verification boundary

These estimates may strengthen an original one-vector reference AFTER
actual geometry establishes its pressure and coefficient prices. For
example the scalar lemma improves ONLY the dotted .4 head graft by
more than .0057 relative to the old .2893 bound. The norm8 lemma can
raise the actual fork-length pressure in a cubic branch without an
enumeration or an unjustified A/2 integrality step.

They supply no actual strong-knot window, head/fork isolation, ground-
field bound, finite catalogue, dimension ceiling or polytope example.
Prospective applications that still require such geometric statements
remain conditional. Full different-author source review and root
acceptance are required before either new estimate becomes accepted.
No frozen paper, source, gate or successful execution has been changed.
