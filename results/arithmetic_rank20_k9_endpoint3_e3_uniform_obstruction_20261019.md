# Uniform K9 exclusion of endpoint3 and endpoint3+e3 dotted rows

19 October 2026. Root-agent hand proof, independently checked and
written out by batch_guard. No enumeration, calibration, scientific
job, Lean or build was run. The earlier fixed-N=9 proof and its
source-prepared worker are preserved unchanged.

## 1. Statement and exact scope

Let K=K9=Q(theta), theta=2cos(2pi/9), with physical embedding
theta=t the largest root of f(X)=X^3-3X+1. In the isolated-d,
spectator-orthogonal actual-quadrilateral setup, suppose the EXPLICIT
opposite-wall product x=-a.c satisfies x>=1. Neither an endpoint3-only
contact nor an endpoint3 contact together with e3 is possible, at
either endpoint, for any N>=8 for which the interfaces below hold.

In the rank20 frontier this covers N=8,...,14 (r=N+2). The algebraic
argument works for every N>=8; it does not independently extend the
geometric inventory to cases where those interfaces are unproved.
In particular c meeting the spectator is not included.

The actual-quad interfaces are those derived in
`arithmetic_rank20_isolated_long_dotted_height_bounds_20261019.md`,
SHA256 `9fc0583124cba17afed393bb569f83835e8f2c4be8ce5cfa64839f0369254516`.
Set A=2x. The actual simple-long cycle gives A in O_K. The opposite
product supplies an integral delta, and corank gives C*delta=M_0.
The exact identities are:

|row|Z|C|M_0|nonphysical condition|
|---|---|---|---|---|
|endpoint3 far|NA+1|[2Z^2-(N+2)(N+1)]/N|4(N+1)|Z_i^2<N+1|
|endpoint3 near|NA+N-1|[2Z^2-(N+2)(N+1)]/N|4(N+1)|Z_i^2<N+1|
|endpoint3+e3 far|2NA+N+2|A(NA+N+2)|N+2|Z_i^2<2(N+2)|
|endpoint3+e3 near|2NA+3N-2|(A+2)(NA+N-2)|N+2|Z_i^2<2(N+2)|

All inequalities indexed i refer to both nonphysical embeddings.
Physically A>=2. The geometric requirements that d be isolated,
c be orthogonal to the spectator, and the displayed integrality,
corank and conjugate-positivity implications are explicit inputs.
No premise that arbitrary ambient opposite walls must be dotted
is substituted for x>=1.

The EXTERNAL integer-ring input is O_K=Z[theta], the standard
maximal-real-cyclotomic integer-ring theorem. Thus every integral
U has a unique expression a+b*theta+c*theta^2 with a,b,c integers.
The cubic is irreducible by the rational-root test, but that fact
alone would not supply the full integer ring.

## 2. Reusable cubic trace-dual calculation

The three roots v<u<t satisfy

    -19/10<v<-15/8,
       1/3<u<7/20,
       3/2<t<31/20.

The f-values at the pairs of endpoints are respectively
(-159/1000,17/512), (1/27,-57/8000), and (-1/8,591/8000).
The derivative 3X^2-3 has fixed nonzero sign on each interval;
these three roots exhaust the degree. All bounds are rational.

Newton traces are

    Tr(1)=3, Tr(theta)=0, Tr(theta^2)=6,
    Tr(theta^3)=-3, Tr(theta^4)=18.

For U=a+b*theta+c*theta^2 define the trace-dual elements

    lambda_c=(2theta^2+theta-4)/9,
    lambda_b=(theta^2+2theta-2)/9.

The displayed traces directly verify

    Tr(U*lambda_c)=c,   Tr(U*lambda_b)=b,
    Tr(lambda_c)=Tr(lambda_b)=0,
    Tr(U)=3a+6c.

Monotonicity of the two quadratics on the certified root intervals
gives the following strict bounds:

|element|physical t|small root u|negative root v|
|---|---|---|---|
|lambda_c|2/9<lambda_c(t)<4/15|-2/5<lambda_c(u)<-3/8|1/8<lambda_c(v)<3/20|
|lambda_b|1/3<lambda_b(t)<2/5|-7/50<lambda_b(u)<-13/100|-1/4<lambda_b(v)<-6/25|

For either lambda its nonphysical values are strictly less than its
positive physical value L=lambda(t). If the two nonphysical values
of U lie in a common interval (ell,h), set T=Tr(U) and
q=Tr(U*lambda). Eliminating the physical U with its trace gives

    q=T*L+sum_(i nonphysical) U_i*(lambda_i-L).

The differences are negative and sum to -3L. Therefore

    L*(T-3h)<q<L*(T-3ell).                         (1)

The sign reversals here are essential. Whenever the multiplying
factors are positive, the physical rational bounds for L can be
substituted directly in(1).

## 3. Endpoint3 far

For either endpoint3-only row the nonphysical C_i are less than
-(N+1). Since C and delta are integral and C*delta=4(N+1),
taking the degree3 norm gives

    C_0<64(N+1),
    Z_0^2<(N+1)(65N+2)/2.                          (2)

In the far row, if A_0>=6, then Z_0>=6N+1. But

    (6N+1)^2-(N+1)(65N+2)/2=N(7N-43)/2>0

for N>=8. Thus 2<=A_0<6. The nonphysical bound on Z implies

    -1/2<A_i<1/4.

Indeed sqrt(N+1)<=N/2-1 and sqrt(N+1)<=N/4+1 for N>=8;
after squaring these follow from N(N-8)>=0. Consequently
1<Tr(A)<13/2. Its trace is a multiple of3, so T=3 or6.

Apply(1) with U=A, ell=-1/2, h=1/4. For lambda_c and T=6,

    7/6<c<2,

which contains no integer. For T=3 it gives

    1/2<c<6/5,

so c=1. Applying lambda_b at T=3 gives

    3/4<b<9/5,

so b=1. The trace then forces a=-1. At the small root,

    A(u)=u^2+u-1<(7/20)^2+7/20-1=-211/400<-1/2,

contrary to the required nonphysical interval.

## 4. Endpoint3 near

Set U=A+1. Then U_0>=3 and Z=NU-1. If U_0>=13/2, its
square contradicts(2), because

    ((13/2)N-1)^2-(N+1)(65N+2)/2
        =N(39N-186)/4>0.

Hence 3<=U_0<13/2. The same two square-root comparisons give

    -1/4<U_i<1/2.

Thus 5/2<Tr(U)<15/2 and T is3 or6. Apply(1) with
ell=-1/4,h=1/2 and lambda_c. At T=3 it gives

    1/3<c<1;

at T=6 it gives

    1<c<9/5.

Neither interval contains an integer. This excludes the near row.

## 5. Endpoint3 plus e3, far

Write M=N+2. In either e3 row C=(Z^2-M^2)/(4N), so its
nonphysical conjugates satisfy C_i<-M/4. The degree3 norm of
C*delta=M therefore gives

    C_0<16M.                                       (3)

In the far row C=A(NA+M). At A_0=4 this is20N+8,
strictly greater than16N+32 for N>=8; the expression increases
for A_0>=2. Hence 2<=A_0<4. Set U=A+1, so

    3<=U_0<5,
    0<U_i<3/4.

For the conjugate interval, the exact endpoints are
(N-2+/-sqrt(2N+4))/(2N). Its lower endpoint is nonnegative
because (N-2)^2-(2N+4)=N(N-6)>=0, and its upper endpoint is
less than3/4 because (N/2+2)^2-(2N+4)=N^2/4>0.

The trace is a multiple of3 and satisfies 3<T<13/2; thus T=6.
Apply(1) with ell=0,h=3/4. It yields

    5/6<c<8/5,       5/4<b<12/5.

Therefore c=1,b=2, and the trace forces a=0. But then

    U_0=t^2+2t>(3/2)^2+2*(3/2)=21/4>5,

contrary to its physical bound.

## 6. Endpoint3 plus e3, near

At A_0=4 the near expression C=(A+2)(NA+N-2) is30N-12,
strictly greater than16N+32 for N>=8. It increases for A_0>=2,
so(3) gives A_0<4. Set U=A+2. Then

    4<=U_0<6,
    1/4<U_i<1.

Here the exact conjugate endpoints are
(N+2+/-sqrt(2N+4))/(2N). The two squared comparisons in the
previous section put them above1/4 and below1 respectively.
Thus 9/2<Tr(U)<8, so again T=6.

Apply(1) with ell=1/4,h=1. The coefficient bounds are

    2/3<c<7/5,       1<b<21/10.

Hence c=1,b=2,a=0 and U=theta^2+2theta. At the negative root
v in(-2,0), this gives U(v)=v(v+2)<0, contradicting U_i>1/4.

## 7. Consequence and remaining boundary

All four endpoint3-only/e3 dotted rows are therefore excluded over
K9 whenever the stated actual-quad interfaces hold, uniformly for
N>=8. Together with the separate common-normal-form norm obstruction
for next3, endpoint4 and endpoint3+e4, this removes every classical
spectator-orthogonal dotted long-contact row over K9. The exceptional
r=10 E8 long-contact rows retain their separate explicit unit proof;
they have not been silently folded into the endpoint inventory.

No claim is made for c contacting the spectator, other ground fields
or incomplete triangle/simultaneous-refund composition. This remains
a local arithmetic exclusion, not an installed global allowance or
a new compact arithmetic dimension bound. The source-prepared fixed
K9,N=9 worker is now optional reproducibility support for a row
already excluded by hand; no worker certificate is claimed.
