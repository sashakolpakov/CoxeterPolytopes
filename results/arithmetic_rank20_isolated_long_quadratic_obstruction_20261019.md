# Quadratic exclusion of the isolated long-contact rows

19 October 2026. Own-paper argument proposed by the root agent and
independently checked by batch_guard. No enumeration, scientific job,
Lean, TeX or other proof process was run.

This note uses the exact geometric and arithmetic interfaces of
`arithmetic_rank20_isolated_long_dotted_height_bounds_20261019.md`,
SHA256 `9fc0583124cba17afed393bb569f83835e8f2c4be8ce5cfa64839f0369254516`.
It does not extend that source's spectator-orthogonality hypotheses.

## Statement and inputs

Let r=N+2 with 8<=N<=14. In the actual quadrilateral under study,
d is isolated from the old white graph, c is orthogonal to the
spectator, and c has one of the listed long-arm contacts. Suppose
the opposite product x=-a.c satisfies x>=1. Then none of those
long-contact rows occurs over K=Q(sqrt3) or Q(sqrt6).

The inputs are the following necessary arithmetic conditions already
derived in the height source. All listed A,C,delta are integral in K,
delta is positive physically, C*delta=M is a positive rational integer,
and every nonphysical conjugate has C'<-a for the stated a>0.

|row|A,Z|C and M|a|
|---|---|---|---|
|endpoint3|A=2x, Z=NA+h, h=1 or N-1|C=[2Z^2-(N+2)(N+1)]/N; M=4(N+1)|N+1|
|next3|A=2x, Z=NA+h, h=2 or N-2|C=(Z^2-2N-4)/N; M=8|2|
|endpoint4|A=2sqrt2*x, Z=NA+2h, h=1 or N-1|C=(Z^2-2N-4)/N; M=8|2|
|endpoint3+e3|A=2x; h=N+2 or3N-2|C=[(2NA+h)^2-(N+2)^2]/(4N); M=N+2|(N+2)/4|

The endpoint3+e4 row requires sqrt2 in K and is therefore absent
from both fields in this statement. The extra r=10 E8 long-contact
rows are already excluded by the separate odd-unit argument in the
height source. These omissions are named dependencies, not results
of an unperformed enumeration.

We also use O_K=Z[sqrt d] for d=3,6, and the declared physical
embedding sqrt d>0. This is the ordinary quadratic integral-ring
fact. In the endpoint4 row A is in K even though sqrt2 need not be;
one must not separately impose sqrt2 in K there.

## The norm bound

Since C and delta are integral and C*delta=M, the nonzero integer
norm of delta gives abs(Norm(C))<=M^2. The physical C is positive,
while its one other conjugate satisfies abs(C')>a. Therefore

    0<C<M^2/a.                                      (1)

Using a weak <= in place of the strict upper bound is harmless below.

## Endpoint3

Equation(1) gives C<=16(N+1), and hence

    Z^2<=(N+1)(17N+2)/2.                            (2)

Since h>=1 and

    (3N+1)^2-(N+1)(17N+2)/2=N(N-7)/2>0,

the physical bound is 2<=A<3. The nonphysical proper-block
inequality is Z'^2<N+1, so

    (-h-sqrt(N+1))/N < A' < (-h+sqrt(N+1))/N.

For h=1 or N-1 and N>=8 these imply -5/4<A'<1/4. For example,
sqrt(N+1)<=N/4+1 follows by squaring from N(N-8)>=0;
this proves the worst near-end lower bound. The far-end upper
bound follows from sqrt(N+1)<=N/4+1 as well. Thus

    3/4 < Tr(A)=A+A' < 13/4.

Write A=p+q*sqrt d with p,q integers. Its trace is the even
integer 2p, so p=1. Since 2<=A<3 and sqrt d>0, q>=1. For d=6
already 1+sqrt6>3. For d=3 the only possibility is q=1.

For the far role, (N/2-1)^2-(N+1)=N(N-8)/4>=0 gives A'>-1/2,
whereas 1-sqrt3<-1/2. For the near role, compare (2):

    ((5/2)N+N-1)^2-(N+1)(17N+2)/2
        =N(15N-66)/4>0.

Consequently A<5/2 in the near role, whereas 1+sqrt3>5/2.
Both endpoint3 roles are excluded.

## Next-to-end3 and endpoint4

For either row equation(1) gives C<=32 and

    Z^2<=34N+4.                                    (3)

In the next3 row A>=2 and h>=2, so Z>=2N+2. But

    (2N+2)^2-(34N+4)=2N(2N-13)>0.

In the endpoint4 row A>=2sqrt2 and h>=1. Therefore

    Z>=2sqrt2*N+2,
    (2sqrt2*N+2)^2-(34N+4)=N(8N+8sqrt2-34)>0

for N>=8. Both contradict(3). This endpoint4 argument uses its
actual integral cycle variable A, without adjoining a square root
that the ground field need not contain.

## Endpoint3 plus e3

Equation(1) gives C<=4(N+2). In the far role

    C=A(NA+N+2)>=2(3N+2)>4(N+2).

In the near role

    C=(A+2)(NA+N-2),

whose excess over A(NA+N+2) is (2N-4)(A+1)>=0.
The same contradiction holds.

These arguments exclude every isolated-d long-contact row of the
height source in the two quadratic fields. They do not cover c
meeting the spectator, do not establish all ten remaining consumers,
and do not by themselves justify a new simultaneous refund or a
dimension bound. Higher-degree rows remain open pending exact
arithmetic or a further analytic obstruction.
