# Every A component of rank at least seven has actual projection above 1/2

20 October 2026. Complete proof candidate, awaiting full independent
review. Paper only: no scientific execution, enumeration, optimizer,
arithmeticity, Lean or TeX. All predecessors remain unchanged.

**Uniform theorem.** At a vertex p farthest from a strict interior
point of a compact hyperbolic Coxeter polytope, every A_r component,
r>=7, has

    ||x_A||^2>1/2,

with arbitrary other finite components. The proof retains the SAME
actual farthest vector and actual feasible core feet. It is not an
edge-only energy theorem. A conditional dimension29 composition is
given in Section5; no accepted global bound is changed by this draft.

## 1. Reviewed base and notation

The r7 case is the corrected source
`compact_universal_a7_two_end_actual_projection_20261020.md`, SHA256
`5ae61f571761f517222bef38f35468462ec8fb9ad9d879850c419a6c316828a5`.
It states and pins the actual reciprocal-shadow, positive projection
and orthogonal-component no-common-opposite interfaces used here.
All paths in this note are under `results/`.

For the remaining r>=8 put N=r+1>=9. In the sum-zero space of R^N,
use unit normals n_i=(e_i-e_(i+1))/sqrt2, inverse Gram

    C_ij=2min(i,j)(N-max(i,j))/N,
    u_i=omega_i/sqrt(C_ii), h_i=x_A.u_i.

Write the strict interior point as E(p-x). Actual edges give
h_i>=tanh(L_i/2)>0. Distinct actual opposite walls give
F(h_i)F(h_j)>=u_i.u_j, F(t)=2t/(1+t^2). The positive-definite
replacement of an ALREADY actual wall is an actual endpoint by
positive acute projection; no new wall is inferred from Gram data.

For c=A p-v an actual own opposite, let Q be its total foreign
Schur cost. Each foreign component has at most one contact, cost
q>=1/4, with gap(1/4,1/3), and necessarily q<A^2. Otherwise the
reciprocal shadow and positive full replacement make c another
actual opposite in that orthogonal component, or give a feasible
null point in the equality case. These are the reviewed universal
cost/no-common-opposite statements; arbitrary spectator ranks are
allowed throughout.

Assume for contradiction that ||x_A||^2<=1/2.

## 2. A fixed-window forcing certificate, uniform in N

Consider the four pairs1/2,3/4,(r-3)/(r-2),(r-1)/r. Their correlations
increase with N when measured at either end. The r7 calibration in
the pinned base gives pair sums>1 for the outer pair and>19/16 for
the inner pair whenever the respective actual opposites are distinct.

Put L=2u1+2u2+u3+u4 and let R be its reflected-index counterpart,
R=2u_r+2u_(r-1)+u_(r-2)+u_(r-3). If all four pairs are distinct,

    x_A.(L+R)>4+2*(19/16)=51/8.                 (1)

For a left-side correlation i<j, u_i.u_j increases with N and is
bounded by sqrt(i/j). Hence

    L^2<12+6sqrt2+4/sqrt3+4sqrt(2/3)+sqrt3<28.

One explicit rational upper bound is2351/84<28, from
sqrt2<10/7, 4/sqrt3<7/3, sqrt(2/3)<5/6 and sqrt3<7/4.
Also, for a=(2,2,1,1),

    L.R=(sum_(i=1..4) a_i sqrt(i/(N-i)))^2,

which decreases with N. At N9 its positive square root is

    sqrt2+2sqrt(2/7)+2/sqrt5<10/7+15/14+9/10=17/5.

All three rational comparisons follow by squaring positive numbers.
Thus

    (L+R)^2<56+578/25=1978/25<80.

Cauchy contradicts(1), since (51/8)^2-(1/2)*80=41/64>0.
After reversing the numbering, a common actual pair is therefore
1/2 or3/4. This relabeling assumes no ambient polytope symmetry.

## 3. Complete rank-uniform finite inventories

Classical finite spherical diagram classification is an explicit
input. Its only irreducibles of rank>=8 are A,B,D and E8 at rank8.
In particular there are no new H/F exceptions as r grows. Affine
principal diagrams are forbidden in a compact Coxeter polytope.

### Common3/4

There is at most one own contact on each side{1,2,3}/{4,...,r}:
two on one side give a cycle retained in the other actual finite
replacement. A contact has finite native label3,4 or5, since a
larger label with an old3 gives a nonpositive retained rank-three path.

For a DOUBLE profile other than positions3,4, one deletion retains
both contacts, connecting its entire two old path components through
c into an irreducible positive rank-r diagram. Both new bonds are
internal in a possible B path; D/E have only3. Hence both bonds are3,
or the white diagram is already nonpositive. Both3 give a principal
affine cycle in the full actual diagram. This explicitly includes
the rank-eight E8 branches; none has a weighted bond.

At3,4 the contact at3 is3/4/5, and the contact at4 is3/4 (its retained
right path has length N-4>=5). The3,3 triangle is affine. Otherwise
the two minimal assignments3,4 and4,3 have respectively

    cosh L3=sqrt(3/(2N))*((N-3)+sqrt2*(N-4)),
    cosh L4=(3sqrt2+4)*sqrt((N-4)/(2N)).

Both increase with N>=9 and exceed3 already at N9. Increasing a
bond or Q only strengthens the appropriate endpoint formula. Thus
all adjacent double profiles force half-tangent square>1/2.

For SINGLE3, every position3,...,r-2 is forbidden by a principal
affine tree: at3 orr-2 there are arms1,2,5 of affine E8; at any
position4,...,r-3 there are arms1,3,3 of affine E7. Longer arms are
simply truncated to the displayed principal subdiagram. At the four
remaining positions the sufficient Q/height bounds are:

| Contact | Q bound | Height A^2 |
| --- | --- | --- |
| 1 | Q<5/8 | Q-(N+1)/(2N)<1/8 |
| 2 | Q<1/2 | Q-2/N |
| r-1 | Q<2/(N-3)<=1/3 | Q-2/N |
| r | Q<1/2+1/(2(N-3)) | Q-(N+1)/(2N)<3/(2N(N-3)) |

The end heights are below1/4. At2,r-1 there is at most one foreign
component, whose q=Q>A^2. Thus all single3 rows are excluded.

A SINGLE4 is the principal affine3,3,4 star at an internal old node.
At1 its left A3 white condition gives Q<1/4; at r the right white
condition gives Q<1/(N-3)<=1/6. Therefore Q=0 and height is negative.

A SINGLE5 can occur only at1 or3: the retained right path has at
least five nodes and cannot extend finitely through5, while the
middle of the retained left A3 cannot take5. Its H4 white cost
3phi^2/8>3/4 forces Q=0. At3, c is also the actual opposite4 with

    cosh L4=3phi^3*sqrt((N-4)/N)>3.

At1 it has actual projection>1/2 by foot F below. A zero own contact
cannot close any own edge. This exhausts common3/4 completely.

### Common1/2

There is at most one right contact on2,...,r and possibly contact1.
For double positions1,j with j>=3, deletion2 connects all retained
own nodes through c into a rank-r finite path or branch. If it is a
branch, all bonds are3 (D or the rank-eight E8 boundary). If a path,
a4 can occur only at its terminal bond1-c, and then j=3 or j=r.
All3 profiles give a principal affine cycle. The4@1+3@r profile
has the principal affine F4 path2-1-c-r-(r-1), labels3,4,3,3.
The sole4@1+3@3 profile has
Q<1/(2(r-1))<=1/14, hence Q=0, and is excluded by actual foot H below.

Adjacent contacts1,2 are retained: the bond at2 is3/4 because it
extends the old A_(r-1) at an end; the bond at1 is any finite native
label. Both3 give affine A2. Every other adjacent row is type L below.

For SINGLE3, internal positions3,...,r-2 have the same principal
affine E7/E8 obstructions. At1, Q<3/4 gives A^2<1/4. At r,
Q<1/2+1/(2r) gives A^2<1/(2rN)<1/4. At r-1 the white condition
Q<2/r<=1/4 leaves none and height is negative. At2,
Q<1/2+1/(2r)<=9/16 allows at most two foreign components.
One contradicts q<A^2=Q-2/N. With two, both first bonds are3, and
c,old2,old1,old3 plus the two contacted foreign nodes are exact
affine two-fork D5. Thus all single3 profiles disappear.

A SINGLE4 is affine internally; at1,Q<1/2 allows at most one
foreign component whose q>A^2=Q-1/N; at r,Q<1/r leaves none and
height is negative. A SINGLE label>=5 is possible only at1. Label6
gives the affine3,6 path. For native>=7, Q<1-beta1^2<1/4 leaves Q=0.
For5, Q=0 is foot F; every Q>0 has Q>=1/4 and is type L.

Thus under the assumed half bound the sole remaining common1/2
profiles are type L: adjacent contacts (not3,3), isolated native>=7,
or isolated5 with Q>0. No foreign branch is suppressed by assuming
the full star finite.

## 4. The actual feet and the uniform far-end contradiction

For clarity, the two excluded feet H,F are the positive actual
combinations explicitly proved in the predecessor sources:

    H: c+(sqrt2-1/2)n1+(sqrt2-1)n2+(sqrt2/2)n3,
       k=sqrt2-5/4,
       A_N^2=3/2+sqrt2-(11/2+3sqrt2)/N;

    F: c+sum_(j=1..4)((j+1)phi-2j)n_j,
       k=34phi-55,
       A_N^2=phi^2*(N-1)/(2N)-1.

H is Section4 of
`compact_universal_a6_actual_projection_three_sevenths_20261020.md`,
SHA256 `022c0cc77512576dfc7bed74be372a94956618498c9ed74508aff53c40a50ae2`;
F is Section3 of the corrected A7 source. Their coefficients are
positive and all actual wall products nonpositive, not merely own
products. Both have Q=0 and U in the A component. At N9,

    A_H^2-9k=(437-300sqrt2)/36>0,
    A_F^2-9k=50(89-55phi)/9>0.

For example sqrt2<10/7 and phi<16181/10000 seal the signs. Their
A_N^2 increase with N while k is fixed. Farthest comparison with
these ACTUAL points gives ||x_A||^2>1/2 in every N>=9, as claimed
where these profiles were excluded.

For the remaining type L, clear the spectators positively to
c_F=c+sum_S(B_S^-1 beta_S).n_S=A p-v_A, c_F^2=1-Q.
Put b=beta1,d=beta2,B=b+2d and

    h=4(b^2+bd+d^2)/3+Q-1, z0=h/B^2,
    lambda1=(4b+2d)/3-h/B,
    lambda2=2B/3-2h/B,
    xi=c_F+lambda1*n1+lambda2*n2=A p-U.

The explicit A7 calculation, independent of N, proves lambda1,lambda2>0,
xi.c=xi.n1=xi.S=0, xi.n2=-3h/(2B)<0, and

    8/81<z0<1/3,   k=-xi^2=B^2(z0-3z0^2)>0.

For completeness, z0>8/81 follows in the three rows as follows:
adjacent native contacts with one>=4 give z0>1/10; isolated>=7
has beta1>9/10 and z0=4/3-1/beta1^2>8/81; isolated5 and Q>=1/4
has z0>=3phi-14/3>2/15. The upper bound is exactly the positive
actual endpoint condition1-b^2-Q>0. Positivity of lambda2 is
2(1-b^2-Q)/B, and B*lambda1=2bd+1-Q>0.

Thus xi is a future timelike positive combination of actual normals,
feasible against EVERY actual wall, with pure A spatial vector U.
The full-vector identities are

    a_N=2/3-2/N>=4/9,
    A^2=B^2(z0+a_N),  U^2=B^2(a_N+3z0^2),
    U_raw=(B/sqrt2)*(2/3+z0,2/3+z0,2/3-2z0,0,...,0),
    U=U_raw-(sqrt2*B/N)*1_N.

The ACTUAL farthest pressure is x_A.U>=A-sqrtk. Its normalized
squared demand g^2 increases with a_N and satisfies g^2>9/20
already at a_N=4/9. Indeed the equivalent quadratic
2523z0^2-720z0+484/9 has negative discriminant-72976/3.

Let W=u_(r-1)+u_r be the far-end pair. Direct coordinates give

    rho=(U.W)/sqrt(U^2)
      =[2/sqrt(N(N-2))+sqrt(2/(N(N-1)))]/sqrt(a_N+3z0^2).

Its numerator decreases with N, its denominator increases, and at
N9,z0>8/81 it obeys rho<61/100 (the exact squared upper seal is
48114/129500<3721/10000). Also W^2<4.

If the far actual opposites are distinct, their pair sum is>1.
The perpendicular Cauchy bound, under ||x_A||^2<=1/2, instead gives

    x_A.W <61/(100sqrt2)+sqrt((1/2-9/20)*4)
          <11/25+9/20=89/100<1.

If the far opposites are common, the reflected finite inventory
leaves only another type-L foot U'. NORMAL reversal gives its raw
sum -sqrt2 B', disjoint from U's raw support. Hence

    U.U'=2BB'/N,
    corr(U,U')<(2/N)/a_N=3/(N-3)<=1/2.

Both normalized demands exceed sqrt(9/20), so the sum-direction
Cauchy inequality gives ||x_A||^2>(4*9/20)/3=3/5>1/2.
This excludes every remaining case and proves the uniform theorem.

## 5. Conditional global composition: exclusion of dimension29

The rates, small exceptional debits and accepted component caps are
collected, with all underlying exact pins, in the independently
written aggregation
`compact_forest_ordinary_a_parameter_and_dominant_classical_frontier_20261020.md`,
SHA256 `dd6f8b0de2d7d0bc835be3eb4e320013aa40792407fb2a747f35bcf04385d5ad`.
That source uses A5/A6 ordinary rates, B5--B8 ordinary rates, all
other ordinary finite types, actual B_r,D_r>9/14 for r>=9, and
component caps A<=20,B<=22,D<=23 in dimension29. Its proof makes
explicit the sole permitted short A1 debit OR one isolated rank-two7
defect of cost>9/100; those two exceptions cannot coexist.

Conditional on acceptance of that collected packet and the present
uniform theorem, a hypothetical compact29-polytope is impossible.
The pinned forest aggregation already forces an A_r with7<=r<=20
at every chosen farthest vertex. The present theorem gives its
projection>1/2. There cannot be a second A>=7, nor a B>=9 or D>=9,
because either would make total projection exceed1. Every other
component is therefore ordinary except the named small debits.

The spectator rank is S=29-r>=9. With no exception its cost>S/14.
With one short A1, discard that rank and retain cost>(S-1)/14>=8/14.
With one rank-two defect, its cost is>9/100+(S-2)/14>=9/100+7/14.
Adding the A component's strict1/2 exceeds1 in EVERY case. No
component rank or cost is counted twice, and no spectator is assumed
ordinary before excluding the second large classical component.

After importing the existing ambient compact ceiling29, the reviewed
composition would yield ceiling28. This statement uses NO claim that
an arbitrary metric face is Coxeter. The new global result remains
a proof candidate until all required source-bound audits pass.

## Contribution and limits

Finite57_research derived the uniform finite inventories, fixed-window
certificate and the all-r transfer. Root independently derived the
general-N core formulas and identified that the half-bound suffices
for the global aggregation, avoiding a stronger unnecessary9/14
target. Actual-wall tools and earlier low-rank work retain their
stated shared contributions. No A9+ rank/14 or actual9/14 statement
is asserted; no new manuscript, job or formal verification was made.
