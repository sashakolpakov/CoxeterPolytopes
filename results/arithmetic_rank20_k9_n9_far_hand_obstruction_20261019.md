# Hand exclusion of the K9, N=9, far-endpoint3 row

19 October 2026. Root-agent proof, independently checked and written
out by batch_guard. No enumeration, scientific job, calibration,
Lean or build was run. The prepared worker and pure helper remain
unchanged as optional reproduction support, not executed evidence.

This settles ONLY the indicated actual-quad arithmetic row. It does
not settle other fields, ranks, contacts, spectator interactions or
the complete simultaneous-refund problem.

## 1. Scope and exact inputs

Use the isolated-d, spectator-orthogonal actual-quadrilateral setup
of `arithmetic_rank20_isolated_long_dotted_height_bounds_20261019.md`,
SHA256 `9fc0583124cba17afed393bb569f83835e8f2c4be8ce5cfa64839f0369254516`.
Fix N=9, r=11, the far simple long-endpoint contact, no e contact,
and the EXPLICIT opposite-wall hypothesis x=-a.c>=1. The actual
cycle gives A=2x in O_K, and the other opposite product gives an
integral delta. The row's exact necessary conditions are

    A_0>=2,
    Z=9A+1,
    C=(2Z^2-110)/9=18A^2+4A-12,
    C*delta=40,
    Z_i^2<10 at both nonphysical embeddings.

Here K=K9=Q(theta), theta=2cos(2pi/9), with minimal polynomial

    f(X)=X^3-3X+1.

The physical embedding is the largest real root. We explicitly use
the external integer-ring input O_K=Z[theta], the standard
maximal-real-cyclotomic integer-ring theorem. Thus every integral
A has a UNIQUE form

    A=a+b*theta+c*theta^2,       a,b,c in Z.

Irreducibility follows separately from the monic cubic rational-root
test: +/-1 are its only possible rational roots, and neither is a
root. Irreducibility alone would not establish the integral basis.

## 2. A short physical bound and trace restriction

Physically C>0. At either other embedding C_i<-10. Taking norms
in C*delta=40, and using the nonzero integral norm of delta, gives

    C_0<40^3/10^2=640,
    Z_0^2=(9C_0+110)/2<2935.

If A_0>=6, then Z_0>=55 and Z_0^2>=3025. Therefore

    2<=A_0<6.

Also sqrt10<19/6, since 360<361. The two conjugate Z bounds give

    -25/54<A_i<13/54           (i nonphysical).       (1)

Newton traces from f are

    Tr(1)=3, Tr(theta)=0, Tr(theta^2)=6,
    Tr(theta^3)=-3, Tr(theta^4)=18.

Consequently T=Tr(A)=3a+6c is a multiple of3. The bounds above imply

    29/27<T<175/27,

so T is3 or6.

## 3. Certified root and trace-dual intervals

Write the three roots v<u<t, with physical root t. The intervals

    -19/10<v<-15/8,
       1/3<u<7/20,
       3/2<t<31/20

are certified directly: the f-values at the respective endpoints are

    (-159/1000,17/512),
    (1/27,-57/8000),
    (-1/8,591/8000).

The derivative 3X^2-3 has fixed nonzero sign on each interval, so
each contains one root. The three disjoint intervals exhaust the
degree3. No numerical root approximation is used.

The trace-dual elements extracting c and b are

    lambda_c=(2theta^2+theta-4)/9,
    lambda_b=(theta^2+2theta-2)/9.

Indeed the displayed Newton traces verify

    Tr(A*lambda_c)=c,       Tr(A*lambda_b)=b,
    Tr(lambda_c)=Tr(lambda_b)=0.

Polynomial monotonicity on the root intervals gives these strict
rational bounds:

|element|physical t|small root u|negative root v|
|---|---|---|---|
|lambda_c|2/9<lambda_c(t)<4/15|-2/5<lambda_c(u)<-3/8|1/8<lambda_c(v)<3/20|
|lambda_b|1/3<lambda_b(t)<2/5|-7/50<lambda_b(u)<-13/100|-1/4<lambda_b(v)<-6/25|

These bounds follow by substituting the rational endpoints in each
quadratic, using its derivative sign. In particular each of the two
nonphysical values is less than the physical value, for BOTH duals.

## 4. The two coefficient bounds

For either lambda above set L=lambda(t)>0 and q=Tr(A*lambda),
so q is respectively c or b. Eliminating A_0 with the trace T gives

    q=T*L+sum_(i nonphysical) A_i*(lambda_i-L).

Both differences are negative. Since Tr(lambda)=0, their sum is
-3L. Applying(1), with its correct sign reversal, yields

    L*(T-13/18)<q<L*(T+25/18).                       (2)

For c and T=6, the physical lambda_c interval gives

    95/81<c<266/135.

Both endpoints lie strictly between1 and2; no integer c exists.
Thus T=3. Applying(2) again gives

    41/81<c<158/135,

so c=1. For b at T=3, the physical lambda_b interval gives

    41/54<b<79/45,

so b=1. Finally T=3a+6c gives a=-1. The only remaining candidate is

    A=theta^2+theta-1.

## 5. Final conjugate contradiction

The function X^2+X-1 increases on the small-root interval. Since
u<7/20,

    A(u)<(7/20)^2+7/20-1=-211/400<-25/54,

contradicting(1). Thus the exact necessary-condition set for this
fixed K9,N=9,far-endpoint3 dotted row is empty by hand proof.

The source-prepared
`scripts/fresh_rank20_k9_n9_endpoint3_far_20261019.py` and
`scripts/exact_totally_real_field_20261019.py` are preserved. There
is no producer result JSON or independent replay output from this
task. The argument above is independent of their execution and
still conditional on the named geometric/arithmetic interfaces and
the explicitly stated full integer-ring input. It does not lower a
dimension ceiling or establish the missing global charge composition.
