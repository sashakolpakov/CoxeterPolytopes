# Universal A7: two end cores force actual projection above 1/2

20 October 2026. Complete paper candidate, awaiting independent review.
No scientific execution, optimizer, arithmeticity, Lean or TeX.
All earlier sources are unchanged.

**Theorem.** Every A7 component at a vertex farthest from a strict
interior point of a compact Coxeter polytope has

    ||x_A||^2>1/2=7/14,

with arbitrary other finite components. This is an ACTUAL farthest
projection theorem, not an edge-only energy statement. No ambient
dimension ceiling is claimed.

## 1. Accepted actual interfaces and a forced common pair

Use unit A7 normals n_i=(e_i-e_(i+1))/sqrt2 in the sum-zero space of
R8. Their inverse Gram is C_ij=min(i,j)(8-max(i,j))/4. Write the
interior point as z=E(p-x); put u_i=omega_i/sqrt(C_ii), h_i=x_A.u_i.
Actual edges give h_i>=tanh(L_i/2)>0. Distinct actual opposites obey

    F(h_i)F(h_j)>=u_i.u_j,   F(t)=2t/(1+t^2).

Positive acute projection, actual reciprocal shadows, and the theorem
forbidding a common opposite across orthogonal components are the
reviewed inputs:

* `compact_farthest_vertex_mixed_replacement_inequalities_20261020.md`,
  SHA256 `3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84`;
* `compact_actual_wall_reciprocal_shadow_20261020.md`,
  SHA256 `c913ec3bd8e0f26a8c679a1ede2b98b251ad3f36efeae60486213ddfa3e8badc`;
* `compact_farthest_vertex_orthogonal_component_opposite_obstruction_20261020.md`,
  SHA256 `4216d06c2e03fd336898b8dac289b9969fa096c8a1d08bd67ddeb8ab5ec2f23d`.

All paths are under `results/`. Classical finite spherical diagram
classification is an explicit input. No arithmetic classification is used.

Assume ||x_A||^2<=1/2. The correlations of pairs1/2,6/7 are sqrt(3/7);
those of3/4,4/5 are sqrt(3/5). Concavity of log F implies that distinct
opposites have pair sums>1 in the first case and>19/16 in the second.
Indeed F(1/2)^2=16/25<sqrt(3/7), while

    F(19/32)^2=1478656/1918225<387/500<sqrt(3/5).

The middle cross-multiplication difference is3025075 and
387^2<150000. If all four pairs are distinct, let

    w=2u1+2u2+u3+2u4+u5+2u6+2u7.

Then x_A.w>51/8. Direct Gram summation gives

    w^2=64/7+32/3+16/5+4+64/sqrt21+64/sqrt105
        +16/sqrt7+64/(3sqrt5)+16/sqrt3+8sqrt(3/5)
       <1208/15.

For the bound, the first four terms sum to less than28. The remaining
six terms are respectively less than14,32/5,32/5,10,28/3,32/5.
Each follows by squaring positive rational bounds. Cauchy contradicts
the contrary projection bound because

    (51/8)^2-(1/2)(1208/15)=359/960>0.

After reversing the numbering if necessary, a common pair is1/2
or3/4. This is relabeling, not an ambient isometry assumption.

## 2. Costs and complete central3/4 inventory

Write an actual common opposite c=A p-v, A>0, with own contacts
beta_i=-c.n_i and total spectator cost Q. A contacted spectator
component has one contact, cost q>=1/4, and no cost lies in(1/4,1/3).
Two contacts to one retained component give a finite-diagram cycle.
Necessarily every foreign q<A^2: q>A^2 gives a positive full
replacement and makes c also an actual opposite in that orthogonal
component; equality gives a positive feasible null point. Both are
forbidden by the accepted interfaces.

For common3/4 there is at most one own contact on each side
{1,2,3}/{4,5,6,7}. Each bond is3,4 or5: a larger bond together with
an old3 is nonpositive in a retained rank-three path.

For two contacts, every position pair except3,4 is either impossible
already in one white replacement or has both bonds3, yielding a
principal affine cycle in the full actual diagram. Explicitly:

| Positions | Sufficient finite-white test |
| --- | --- |
| 1,4;1,7;2,4;2,7;3,5;3,7 | deletion3 or4 is a rank-seven path; both variable bonds internal |
| 1,5;1,6;2,5;2,6 | deletion3 is D7 or E7; all bonds3 |
| 3,6 | deletion4 is D7; all bonds3 |

At3,4 both bonds may be3/4/5;3,3 is affine. Otherwise one endpoint
formula gives cosh L>3: at the two minimal assignments3,4 and4,3
these are respectively cosh L3=(15/8+3sqrt2/2)/sqrt(3/4)>4 and
cosh L4=2+3sqrt2/2>4. Thus its half-tangent square
exceeds1/2. No omitted H/F exception exists at rank seven.

For a single3, the full own costs and sufficient Q bounds are:

| j | Own cost | Q bound | A^2 |
| --- | ---: | ---: | --- |
| 1 | 7/16 | <5/8 | Q-9/16 |
| 2 | 3/4 | <1/2 | Q-1/4 |
| 3 | 15/16 | <5/8 | Q-1/16 |
| 4 | 1 | <3/5 | Q |
| 5 | 15/16 | <2/5 | Q-1/16 |
| 6 | 3/4 | <2/5 | Q-1/4 |
| 7 | 7/16 | <3/5 | Q-9/16 |

One foreign component contradicts q<A^2 (equality at j4). All rows
except j3,4 have at most one or have A^2<1/4. For j3,4 two foreign
components must both have first bond3, since Q<5/8. The six nodes
c,j,j-1,j+1 and their two contacted foreign nodes form the exact
affine two-fork D5 diagram, forbidden compactly. This also handles
arbitrary spectator ranks. Q=0 gives nonpositive height in every row.

Single4 at an internal node is the principal affine3,3,4 star. At
ends1,7 its white condition gives Q<1/4, so Q=0 and height<0.
Single5 is possible only at1 or3: the retained right A4 cannot
extend finitely through a5, and the middle of the retained left A3
cannot take a5. The left H4 cost3phi^2/8>3/4 forces Q=0. Contact3
gives the actual opposite4 with cosh L4>3. Contact1 is removed by
the following explicit actual foot, also used below.

## 3. An isolated end5 actual foot

For an isolated5 at old1 with no spectators, put phi=(1+sqrt5)/2.
The positive combination

    eta=c+sum_(j=1..4)((j+1)phi-2j)n_j

satisfies eta.c=eta.n1=eta.n2=eta.n3=0,
eta.n4=3phi-5<0 and eta^2=55-34phi=-k0<0. Every coefficient is
positive. Its products with every other actual wall are nonpositive
by all-wall acuteness; its positive time coefficient makes it future.
It is therefore an ACTUAL point after normalization. Here

    A^2=(7phi-9)/16,  k0=34phi-55,
    A^2-9k0=(7911-4889phi)/16>0.

For an explicit seal phi<16181/10000 gives numerator>1091/10000;
this upper bound follows from11181^2>5*5000^2. Also k0>0 follows
from5*17^2>38^2. Writing eta=A p-U, farthest comparison gives

    ||x_A||^2 >=(A-sqrt(k0))/(A+sqrt(k0))>1/2.

This excludes the last common3/4 row. It also excludes an isolated
end5 with Q=0 in the common1/2 inventory below. The foot is pure A;
no projection contribution from spectators is discarded.

## 4. Complete end1/2 inventory

There is at most one contact on the retained right path2..7 and
possibly a contact at1. Native labels at retained nodes are finite.
The complete double list is:

| Positions | Finite-white consequence |
| --- | --- |
| 1,2 | bond at2 is3/4, bond at1 any finite native label;3,3 affine |
| 1,3 | deletion2 path: bond at3 is3, bond at1 is3/4 |
| 1,4 or1,6 | deletion2 E7: both3 |
| 1,5 | deletion2 already contains affine E6, so impossible |
| 1,7 | deletion2 path: bond at7 is3, bond at1 is3/4 |

Both3 always give an affine cycle. The4@1+3@7 row contains the
principal affine F4 path2-1-c-7-6 with labels3,4,3,3. The sole
4@1+3@3 row has Q=0 and actual A projection>1/2 by the explicit
positive foot in Section4 of
`compact_universal_a6_actual_projection_three_sevenths_20261020.md`,
SHA256 `022c0cc77512576dfc7bed74be372a94956618498c9ed74508aff53c40a50ae2`.
That calculation states and proves its profile-specific extension
to every r>=6; no A6 theorem or additional finite inventory is used.

For a single3 the own costs are7/16,3/4,15/16,1,15/16,3/4,7/16.
Sufficient Q bounds for positions1..7 are respectively

    3/4, 4/7, 2/7, 1/7, 1/7, 2/7, 4/7.

All except j2 give A^2<1/4 and are removed by q<A^2. At j2 there
are at most two foreign components. One contradicts q<A^2; with
two, each other's cost exceeds1/4, so each is at least1/3, contrary
to Q<4/7. Single4 is affine at internal nodes. At1, Q<1/2 allows at
most one foreign component, whose q>A^2=Q-1/8; at7,Q<1/7 allows
none and height is negative. Thus neither occurs.

A single label>=5 is possible only at1: a retained A6 has no finite
extension with such a bond. Label6 gives the principal affine3,6
path and is forbidden. For label>=7, Q<1-beta1^2<1/4, so Q=0.
For label5, Q=0 is excluded by Section3; every Q>0 has Q>=1/4.

Consequently the only end-common rows still needing analysis are:

* adjacent contacts b=beta1,d=beta2, at least one bond>=4, with
  d corresponding to label3 or4;
* single native>=7 at1, with Q=0;
* single5 at1, with Q>0.

These are handled uniformly next, including all foreign costs.

## 5. A pure-A dual foot and the opposite-end contradiction

First clear all contacted spectator components positively:

    c_F=c+sum_S(B_S^-1 beta_S).n_S=A p-v_A,
    c_F^2=1-Q.

This is not asserted to be an actual native normal. Put b=beta1,
d=beta2 (d=0 in a single-contact row), B=b+2d, and

    h=4(b^2+bd+d^2)/3+Q-1,   z0=h/B^2.

All remaining rows have h>0. The positive full endpoint deleting2
gives 1-b^2-Q>0, equivalently z0<1/3. Define

    lambda1=(4b+2d)/3-h/B,
    lambda2=2B/3-2h/B,
    xi=c_F+lambda1*n1+lambda2*n2.

Both lambda are positive: lambda2=2(1-b^2-Q)/B, and
B*lambda1=2bd+1-Q>0. Direct multiplication gives

    xi.c=xi.n1=xi.S=0,  xi.n2=-3h/(2B)<0,
    xi^2=-k,  k=h-3h^2/B^2>0.

Thus xi is a positive combination of actual normals, has all actual
wall products nonpositive, and has positive time coefficient A. It
is an ACTUAL feasible future point. Write xi=A p-U with U purely
in A7. Exact full and core Schur costs give

    A^2=h+5B^2/12,
    U^2=A^2-k=B^2(5/12+3z0^2).

The farthest inequality x_A.U>=A-sqrtk yields

    g^2:=((A-sqrtk)/sqrt(U^2))^2>11/25.       (2)

Indeed A^2/k>324/49, since its cleared difference is the strictly
positive quadratic972z0^2-275z0+245/12, with discriminant-3755.

We also have z0>8/81 in EVERY remaining row. For a double row,
one native contact is at least1/sqrt2. The increasing quadratic
h(Q=0)-B^2/10 has endpoint minima
(28sqrt2-27)/120 or(14sqrt2-9)/60, both positive. Thus z0>1/10.
For single>=7, beta1>9/10, so z0=4/3-1/beta1^2>8/81.
The elementary bound cos(pi/7)>9/10 follows from
8t^3-4t^2-4t+1 at t=9/10 being-1/125, and its strict increase on
[sqrt3/2,1]. For single5 and Q>=1/4,
z0>=4/3-3/phi^2=3phi-14/3>2/15>8/81, using phi>8/5.

Put W=u6+u7 and rho=(U.W)/sqrt(U^2). The raw coordinates of U
are supported on1,2,3, followed by the subtraction of their mean.
Their raw sum is sqrt2 B. Therefore

    U.W=B(1/sqrt12+1/sqrt28),
    0<rho^2=(5/42+1/(2sqrt21))/(5/12+3z0^2)
       <603612/1170300<324/625,
    W^2=2+2sqrt(3/7)<10/3.                       (3)

Here5/42+1/(2sqrt21)<23/100, and the last fraction comparison
has positive cross-product difference1919700.

If the actual opposites6/7 are distinct, their pair law gives
x_A.W>1. On the other hand, decompose along U/sqrt(U^2). By(2),(3)
and ||x_A||^2<=1/2, Cauchy gives

    x_A.W <18/(25sqrt2)+sqrt((1/2-11/25)*10/3)
          =18/(25sqrt2)+1/sqrt5<51/100+9/20=24/25.

The rational bounds follow by squaring positive numbers. This is
the contradiction from an ACTUAL end core and the far-end walls.

If opposites6/7 are common, reverse Section4. Every discarded row
already contradicts the assumed1/2 bound; otherwise its pure-A
foot U' also satisfies(2). Reversal of the NORMALS is implemented
by e_j mapping to -e_(9-j), so its raw support is on6,7,8 with
raw sum -sqrt2 B'. The two raw supports are disjoint. Therefore

    U.U'=+BB'/4,
    (U.U')/sqrt(U^2*U'^2)
      =1/(4sqrt((5/12+3z0^2)(5/12+3z0'^2)))<3/5.

The two normalized directions have sum of squared norm less than
16/5, whereas their scalar products with x_A sum to more than
2sqrt(11/25). Cauchy now gives

    ||x_A||^2>(44/25)/(16/5)=11/20>1/2.

This exhausts both far-end alternatives and proves the theorem.

## Scope and contributions

The new mechanism is the cleared A2 actual dual foot together with
the opposite-end shadow law; it retains the full component vector.
The source was derived by finite57_research, building on the shared
actual-wall tools and root's emphasis on coupled rather than scalar
component pressures. Root's full audit caught and corrected the
reversal sign in the far-common branch; the positive-correlation
sum-direction estimate above is the corrected proof. It uses no
edge-only replacement for actual
farthest comparison. A5/A6 review status is separate; the only A6
import above is its displayed independently checkable profile foot.
A8 and large-A projection9/14 remain separate obligations. No new
ambient bound, manuscript or machine-verified theorem is asserted.
