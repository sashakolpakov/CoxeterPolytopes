# Four short universal component arguments

20 October 2026. Paper proof only; review pending. No computation,
enumeration, build or formal verification. Arithmeticity is not assumed.

At a farthest vertex of ANY compact hyperbolic Coxeter polytope,

    B3 energy > 1/4;
    B4, D4 and F4 energy > 1/3.

All other finite components are arbitrary. These bounds use the same
actual-wall method as the universal small-A source, rather than new
spectator inventories. The method was co-developed by root, batch_guard
and resource_cleanup; the four short native arguments below are by
resource_cleanup, awaiting independent review.

The supporting source is
`results/compact_universal_small_a_component_energy_20261020.md`,
SHA256 801581c4b855289aca7982b2d0f4bb6f5849f151a78f40b86c9d403e84f415e6.
Its Sections1--2 derive all tools recalled here, including the universal
no-common-opposite theorem and small-height isolation.

## Common argument

For unit old Gram B, inverse C, d_i=C_ii and unit dual directions u_i,
the farthest-point vector x satisfies ||x||^2<1 and
x.u_i >= tau_i=tanh(L_i/2). Distinct actual opposites force
tanh L_i tanh L_j >= k=u_i.u_j. If component energy is at most e,
log-concavity of H(t)=2t/(1+t^2) gives the upper bound

    tanh L_i tanh L_j <= H(sqrt(e(1+k)/2))^2.

For e=1/4 this is below k whenever k>=1/sqrt(3). For e=1/3 it is
below k whenever k>=1/sqrt(2): at the latter endpoint the cleared
inequality is 51k>34. In both cases the cleared difference is increasing
in k>0. Thus if adjacent dual correlations meet these thresholds, ONE
actual opposite c closes all directions of the component.

Write h=A^2 for c's squared height, beta for its old contacts, own cost
q=beta^T C beta, and v_i=(C beta)_i/sqrt(d_i). Then

    h >= q-1,       h <= R^2 min_i v_i^2,
    R^2=16/25 at e=1/4; R^2=3/4 at e=1/3.             (1)

If h<1/4, small-height isolation says c has NO contacts to other whole
components. Hence h=q-1. A single-contact case with q=1 is impossible
even without this height restriction: the old component plus c has a
positive-null Gram vector (C beta,1), forbidden in a compact chamber.

A multiple contact whose old connecting path omits an old node is also
impossible. Release that omitted node, which has the SAME opposite c.
The retained white diagram contains the cycle formed by c and the path;
its Gram is nonpositive on the all-ones cycle vector. This contradicts
the actual finite endpoint. In a path component only the two endpoints
may therefore be contacted. In D4 every path omits a node, so only one
contact is possible. These are actual-white arguments, not a blanket
prohibition of hyperbolic principal subdiagrams.

## B3

In path order (4,3),

    C=[[3,2sqrt(2),sqrt(2)],
       [2sqrt(2),4,2],[sqrt(2),2,2]].

Adjacent dual correlations are sqrt(2/3),1/sqrt(2). At energy<=1/4,
one c therefore closes all three. For contacts u,w>=1/2 at the two
endpoints,

    q-(16/25)v_2^2
      =(43u^2+34w^2+18sqrt(2)uw)/25
      >=(77+18sqrt(2))/100>1.

The last comparison is 18sqrt(2)>23. This contradicts (1).

For a single contact beta, the own diagonals d and the ratios
sigma=min_i(v_i^2/q) are

    contacted node:       1       2       3
    d:                    3       4       2
    sigma:               1/3     1/2     1/3
    q upper bound:      75/59   25/17   75/59.

The bounds follow from q-1 <= (16/25)sigma*q. Native contact squares
are cos^2(pi/m), m>=3, or exceed1. Node1 leaves only label3, q=3/4;
node2 only label3, q=1; node3 only labels3/4, q=1/2 or1. For the last
cutoff use 75/118<cos^2(pi/5). The q=1 cases are positive-null. The
others have h<=4/25 or h<=16/75<1/4, so isolation gives h=q-1<0.
Thus B3 energy is greater than1/4 universally.

## D4

The inverse has center diagonal4, leaf diagonal2, center/leaf entries2,
and distinct-leaf entries1. All adjacent dual correlations are1/sqrt(2),
so energy<=1/3 makes one c close all four. Multiple contacts were already
excluded by the common argument.

For a leaf contact, q=2beta^2 and min_i v_i^2=q/4. Equation(1) gives
q<=16/13 and beta^2<=8/13<cos^2(pi/5). Only labels3/4 remain, with
q=1/2 or1. Here h<=3/16<1/4, so isolation gives h=q-1<=0.
For a center contact q=4beta^2 and min_i v_i^2=q/2, giving q<=8/5
and beta^2<=2/5. Only label3 remains, with q=1, positive-null.
Thus D4 energy is greater than1/3 universally.

## B4

In path order (4,3,3),

    C=[[4,3sqrt(2),2sqrt(2),sqrt(2)],
       [3sqrt(2),6,4,2],[2sqrt(2),4,4,2],
       [sqrt(2),2,2,2]].

The adjacent dual correlations are sqrt(3)/2,sqrt(2/3),1/sqrt(2).
At energy<=1/3 one wall closes all four. For two endpoint contacts,

    q-(3/4)v_3^2
      =(5/2)u^2+(5/4)w^2+(sqrt(2)/2)uw
      >=15/16+sqrt(2)/8>1.

For single contacts, the same calculation in(1) gives

    contacted node:       1       2       3       4
    d:                    4       6       4       2
    sigma:               1/4     1/3     1/2     1/4
    q upper bound:      16/13    4/3     8/5    16/13.

At node2 even the minimum beta^2=1/4 is too large. Nodes1/3 leave
only label3 with q=1. Node4 leaves labels3/4 with q=1/2 or1 and
h<=3/16<1/4. The positive-null or isolation contradictions again apply.
Thus B4 energy is greater than1/3 universally.

## F4

In path order (3,4,3),

    C=[[4,6,4sqrt(2),2sqrt(2)],
       [6,12,8sqrt(2),4sqrt(2)],
       [4sqrt(2),8sqrt(2),12,6],[2sqrt(2),4sqrt(2),6,4]].

Adjacent dual correlations are sqrt(3)/2,2sqrt(2)/3,sqrt(3)/2.
At energy<=1/3 one wall closes all four. Endpoint contacts satisfy

    q-(3/4)v_2^2=(7/4)u^2+2w^2+sqrt(2)uw
                      >=15/16+sqrt(2)/4>1.

A single endpoint contact has d=4,sigma=1/2, so q<=8/5 and only
label3,q=1 remains. A single inner contact has d=12,sigma=2/3,
so q<=2 and beta^2<=1/6<1/4: none exists. Thus F4 energy is
greater than1/3 universally.

## What these bounds add to the mixed-forest argument

B3 has positive K^-1*1 and uniform covering coefficient

    T_B3=1^T K^-1 1=9-2sqrt(6)-2sqrt(2)>5/4.

Indeed sqrt(6)<49/20 and sqrt(2)<99/70 imply their sum<31/8.
Hence if all B3 half-tangent squares are at least1/5, its energy is
greater than1/4, or rank/12. The B4,D4,F4 uniform-cover coefficients
need not meet the earlier5/12-per-rank threshold; their INTRINSIC bounds
already give energy greater than1/3=rank/12. Thus each exceptional
propagation step must use the maximum of the intrinsic and cover bounds,
rather than incorrectly asserting a large cover coefficient for all
four types.

The isolated-A2_7 and small isolated-A1 propagation arguments in the
pinned source force every OTHER half-tangent square to be at least1/5.
They remain valid with these four additional components, since every
other block then costs more than its rank/12 by the cover or intrinsic
bound. The baseline per-rank costs also dominate1/14. At dimension14,
the case of one A1 and an isolated dotted-A1 debit still requires an
odd-rank block; among these component types that is A3 or B3, each
paying the same extra1/28 above its rank/14 rate. Thus the exact previous
bookkeeping excludes farthest forests with component types

    A1,A2,A3,A4,B3,B4,D4,F4       in n>=14.

This is NOT yet a claim about ALL finite components of rank at most4:
H3,H4 and the mixed I2(m>=4) interfaces are additional obligations.
Large connected A,B,D components also remain open. The value of the
argument is that arbitrary spectators and cross-component consumers
were handled once by the common geometric mechanism.
