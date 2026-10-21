# Explicit finite height bounds for the remaining isolated long-contact quads

19 October 2026. **New written necessary bounds, not an enumeration
result.** The endpoint3 denominator argument and E8 unit argument were
supplied by the coordinating root agent; the remaining scalar reductions
and uniform bounds are derived here from the same method. No scientific
program, search, optimizer, Lean or build ran. Pinned files are unchanged.

## 1. Exact scope and purpose

Use the actual D_r+T_s ambient20, eight-field setup in
`arithmetic_rank20_isolated_long_finite_product_20261019.md`.
Thus10<=r<=16,N=r-2 in{8,...,14}, d contacts only the old
leaf b, c.b=c.d=0, c is perpendicular to the ENTIRE spectator,
and c has a finite contact with the old D_r long arm. Here take
the unresolved ambient dotted product

    a.c=-x,       x>=1.

All the corner roots are actual outward facet roots. Let

    s0=(N+2)/(4N),       t=cbar^2,
    kappa=-abar.cbar-x,       R=(x+kappa)^2/t,
    delta=4*beta^2-4=1/(R-s0)>0,       beta=-b.d>1.

Cyclic integrality makes delta integral. Its other conjugates
satisfy-4<delta'<-4+8/r<0, and the old-vertex-plus-c cost

    qc=1-t+2N*(x+kappa)^2

satisfies qc'<1 at every nonphysical embedding. All cycle
variables below are integral elements of the ACTUAL K, whose
degree d is at most5.

**Result.** Each ordinary dotted long-contact row has an explicit
finite physical bound on its integral cycle variable, as listed
below. The exceptional E8 row is excluded outright. These bounds
turn a possible later exact algebraic-integer enumeration into a
finite proof obligation. No element enumeration, absence claim or
completed refund is asserted here.

## 2. The elementary divisor-to-height principle

Suppose a nonzero integral scalar C satisfies

    C*delta=M,       C>0 physically,
    C'<-a<0 at every nonphysical embedding,

where M,a are positive rational integers or rational constants as
used below. If M is a positive INTEGER, integrality of delta gives

    |Norm(C)|<=M^d,
    C<=M^d/a^(d-1).                                (1)

For d>=2 the last inequality is strict; the weak version covers
d=1 as well. This uses actual conjugate restrictions and integral
divisibility. It is not a numerical height cutoff chosen by a search.

## 3. A single endpoint3 contact

Let A=2x, integral by the actual simple old-path cycle. Put

    h=1 (far), or h=N-1 (near),
    Z=N*A+h,       t=(N+1)/(2N),
    C=[2Z^2-(N+2)*(N+1)]/N.

C is integral: Z=h modulo N, h^2=1 modulo N, and the numerator
has residue2-2=0. Directly from the quad corank relation,

    C*delta=4*(N+1).                                (2)

The proper-block condition qc'<1 is exactly Z'^2<N+1.
Hence

    -(N+2)*(N+1)/N<C'<-(N+1).

Equation(1) gives

    C<=4^d*(N+1)<=1024*(N+1),
    Z^2<=(N+1)*(1025N+2)/2.

Since h>=1 and Z=N*A+h>0, these inequalities force A<24.
Indeed

    (24N+1)^2-(N+1)*(1025N+2)/2
      =N*(127N-931)/2>0       (N>=8).

Thus every such actual row has

    2<=A=2x<24,
    (-h-sqrt(N+1))/N<A'<(-h+sqrt(N+1))/N.            (3)

No primitivity of A or arbitrary bound on an ambient unit-normal
entry is assumed.

## 4. A single next-to-end3 contact

Again A=2x is integral. Set

    h=2 (far), or h=N-2 (near),
    Z=N*A+h,       t=2/N,
    C=[Z^2-2N-4]/N.

The numerator is divisible by N because h^2=4 modulo N.
The exact identity and the nonphysical cost inequality are

    C*delta=8,       Z'^2<4,
    -2-4/N<C'<-2.

Therefore C<=2*4^d<=2048 and Z^2<=2050N+4. Since h>=2,

    (16N+2)^2-(2050N+4)=N*(256N-1986)>0

at N>=8. Thus

    2<=A=2x<16,
    (-h-2)/N<A'<(-h+2)/N.                           (4)

For N divisible by4 this row is actually excluded without a
search. In the far case C=N*A^2+4A-2, and in the near case
C=N*(A+1)^2-4*(A+1)-2. Hence C/2 is integral and odd modulo2,
while (C/2)*delta=4. It is a unit, but every nonphysical
|C'/2|>1 and its physical value exceeds1, a norm contradiction.
This excludes N=8,12, that is r=10,14, for both next-to-end rows.
The remaining bounds(4) hold whether or not this extra exclusion
is used.

## 5. A single endpoint4 contact

The integral actual-cycle variable is now A=2*sqrt2*x, not2x.
Let h=1 (far) or N-1 (near), and put

    Z=N*A+2h,       t=1/N,
    C=[Z^2-2N-4]/N.

Again C is integral, C*delta=8 and qc'<1 is exactly Z'^2<4.
Thus the SAME height estimate as Section4 gives

    2sqrt2<=A<16,
    (-2h-2)/N<A'<(-2h+2)/N.                         (5)

The formula does not assume x or sqrt2 separately belongs to K;
the cycle makes their indicated product integral. In the subcase
N divisible by4, C=N*A^2+4h*A+[4h^2-2N-4]/N is2 modulo4
for BOTH endpoints (the near constant is4N-10). Hence C/2 is
again an odd unit from (C/2)*delta=4, and its conjugate norm
contradicts unity. Thus the endpoint4 row is also excluded for
r=10,14. The other cases remain bounded by(5), not claimed empty.

## 6. Simple endpoint plus e3

Here A=2x is integral, t=(N+2)/(4N), and put

    h=N+2 (far), or h=3N-2 (near),
    Z=2N*A+h,
    C=[Z^2-(N+2)^2]/(4N).

The integral expansions are

    far:  C=A*(N*A+N+2),
    near: C=(A+2)*(N*A+N-2).

The exact corank identity is C*delta=N+2. The nonphysical
cost inequality is Z'^2<2(N+2), giving

    -(N+2)^2/(4N)<C'<-(N+2)/4.

Thus C<=4^(d-1)*(N+2)<=256*(N+2). Both physical expansions
give C>=N*A^2. Since324N>256*(N+2) for N>=8,

    2<=A=2x<18,
    (-h-sqrt(2(N+2)))/(2N)<A'
       <(-h+sqrt(2(N+2)))/(2N).                    (6)

These inequalities retain both endpoint roles. No near-end
attachment is substituted for its distinct far-end counterpart.

## 7. Simple endpoint plus e4

The old a-e-c-long-path cycle puts sqrt2 in K. The other actual
simple cycle then makes A=2x integral. Among the eight fields,
this forces K=K16 or Q(sqrt2,sqrt3), so d=4. Set

    h=1 (far), or h=N-1 (near),
    Z=2N*A+2h+N*sqrt2,
    C=[Z^2-2N-4]/N.

The quotient C is integral because Z=2h modulo N and h^2=1
modulo N. Corank and positivity again give

    C*delta=8,       Z'^2<4,
    -2-4/N<C'<-2.

Now C<=2*4^4=512 and Z^2<=514N+4. Physically sqrt2>1 and
h>=1, so A>=4 would give Z>9N+2, impossible since

    (9N+2)^2-(514N+4)=N*(81N-478)>0       (N>=8).

Thus the explicit sharper bound is

    2<=A=2x<4,
    (-2h-N*sqrt2'-2)/(2N)<A'
       <(-2h-N*sqrt2'+2)/(2N).                     (7)

The embedded sign sqrt2' is part of the actual embedding, not
an arbitrary independently chosen sign.

## 8. The exceptional E8 continuations are excluded outright

At r=10,N=8, c.l3 or c.l5 simple gives the finite E8
extensions of A7, with t=1/16 and kappa=k/16, k=5 or3.
They must not be omitted from the whole finite inventory.

The root's exact unit argument is as follows. A=2x is integral,
Z=8A+k, and

    (Z^2-5)*delta=16,
    C=(Z^2-5)/4
      =16A^2+4kA+(k^2-5)/4,
    C*delta=4.

The constant is1 or5, so C is odd modulo2. Since C divides4,
it is a unit. Nonphysical positivity gives Z'^2<1, hence
-5/4<C'<-1. Physically x>=1 gives C>1. Every absolute
conjugate therefore exceeds1, contradicting its unit norm.
This exclusion needs no ground-field list.

## 9. What the finite bounds do and do not verify

The table below records only necessary ranges, not an assertion
that any element in them exists or gives an actual quad:

|contact type|integral cycle variable|physical range|
|---|---|---|
|endpoint3 only|2x|[2,24)|
|next-to-end3 only|2x|[2,16)|
|endpoint4 only|2sqrt2*x|[2sqrt2,16)|
|endpoint3 plus e3|2x|[2,18)|
|endpoint3 plus e4|2x|[2,4)|
|r10 E8 extension|2x|excluded|

Equations(3)--(7) also bound EVERY nonphysical conjugate. For any
chosen certified integral basis beta_j of O_K with embedding matrix
M, the coordinate vector is M^(-1) times these conjugate values.
Thus |a_j| is bounded by the sum of the corresponding absolute
inverse-matrix entries times the stated finite conjugate bounds.
This gives a mathematically complete finite coefficient box; it is
not a heuristic cutoff on physical height. A later checker would
still have to certify the basis/embeddings, enumerate every box
element, and test the exact original integrality and positivity
conditions. NONE of that execution is claimed here.

All c-spectator contacts are explicitly excluded by the hypotheses.
The independent complete half-refund/triangle composition remains
open even if every finite box were later proved empty. These
necessary scalar bounds do not by themselves lower a vertex cost,
change a frontier table, or prove a dimension ceiling.
