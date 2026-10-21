# Universal A8: three actual core feet give projection above 4/7

20 October 2026. Complete paper candidate, awaiting independent review.
No scientific execution, enumeration, arithmeticity, Lean or TeX.
All accepted and submitted predecessor sources are unchanged.

**Theorem.** At a vertex farthest from a strict interior point of a
compact Coxeter polytope, every A8 component has actual projection
norm squared greater than4/7=8/14, with arbitrary finite spectators.
This is not an edge-only energy theorem or an ambient dimension bound.

## 1. Inputs and forcing one of four common pairs

Use n_i=(e_i-e_(i+1))/sqrt2 in the sum-zero space of R9,
C_ij=2min(i,j)(9-max(i,j))/9, u_i=omega_i/sqrt(C_ii), and
z=E(p-x), h_i=x_A.u_i. Assume ||x_A||^2<=4/7. The actual-wall
tools, spectator cost bounds and positive projection are those fully
stated and pinned in
`compact_universal_a7_two_end_actual_projection_20261020.md`, SHA256
`5ae61f571761f517222bef38f35468462ec8fb9ad9d879850c419a6c316828a5`.
All paths below are under `results/`. In particular distinct actual
opposites obey F(h_i)F(h_j)>=u_i.u_j, F(t)=2t/(1+t^2). A contacted
foreign component has one contact, cost q>=1/4, with gap(1/4,1/3),
and necessarily q<A^2 for an actual own opposite c=A p-v.

Pairs1/2,7/8 have correlation sqrt7/4; pairs3/4,5/6 have correlation
sqrt(5/8). Concavity of log F gives respective distinct-wall sums
greater than41/40 and1219/1000. Explicit seals are

    F(41/80)^2=43033600/65302561<33/50<sqrt7/4,
    F(1219/2000)^2=23775376000000/30095768093521
       <79/100<sqrt(5/8).

The two cross-multiplication differences are3304513 and28079388159.
If all four pairs are distinct, set

    w=2u1+2u2+u3+u4+u5+u6+2u7+2u8.

Then x_A.w>3269/500 and exact Gram summation gives

    w^2=1116/35+60/sqrt7+36/sqrt10+72/sqrt70<747/10.

The last three terms are less than227/10,57/5,87/10 respectively,
giving the stronger bound2614/35<747/10. Cauchy contradicts the
contrary projection bound because

    (3269/500)^2-(4/7)(747/10)=104527/1750000>0.

Thus one of1/2,3/4,5/6,7/8 has a common actual opposite. Reverse
the numbering if necessary: it suffices to consider1/2 or3/4.
Reversal acts on coordinates as e_j mapping to -e_(10-j); this
minus sign is retained throughout the two-foot comparisons below.

## 2. Complete common3/4 inventory

There is at most one own contact on each side{1,2,3}/{4,5,6,7,8}.
Two on one side give a cycle in a retained positive white diagram.
Each contact is finite and native. One old3 adjacent to a new bond
>=6 makes a nonpositive retained rank-three path.

For DOUBLE contacts other than3,4, at least one deletion retains both
contacts and connects its two old path components through c. This
is a positive connected rank-eight diagram. The finite possibilities
are A8,B8,D8,E8; in a B8 path the unique4 is terminal, whereas both
new bonds here are internal. Hence both new bonds must be3, or that
white diagram is already impossible. Both3 give a principal affine
cycle in the full actual diagram. This argument covers all fourteen
nonadjacent position pairs, including both E8 branch positions.

At3,4 the bonds are p=3/4/5 at3 and q=3/4 at4. The3,3 triangle is
affine. Otherwise one of the endpoint formulas is at least its
minimal3,4 or4,3 value:

    cosh L3 >=(2+5sqrt2/3)/sqrt(2/3)>11/3,
    cosh L4 >=(5sqrt2/3+20/9)/(sqrt10/3)>11/3.

They increase with either bond and foreign cost. Thus a corresponding
half-tangent square exceeds4/7 and contradicts the assumed bound.

For SINGLE3 at3 or6, the full old diagram plus c contains the exact
affine E8 tree with arms1,2,5. At4 or5 it contains the principal
affine E7 tree with arms1,3,3. These are forbidden compactly,
independently of spectators. At1,2,7,8, sufficient Q bounds are
5/8,1/2,1/3,7/12 respectively. The own costs are4/9,7/9,7/9,4/9.
The height A^2=own+Q-1 is less than1/4 at ends1,8. At2,7 there is
at most one foreign component and its q>A^2. Thus every row is
excluded by the foreign-cost/no-common-opposite interface.

A SINGLE4 is affine3,3,4 at an internal old node; at1 or8, the
white bounds Q<1/4 or Q<1/6 leave Q=0 and negative height.
A SINGLE5 can occur only at1 or3: the retained right A5 has no
finite extension through5, and the middle of the retained left A3
cannot take5. Its white H4 cost3phi^2/8>3/4 forces Q=0. At3, c is
also the actual opposite4, with cosh L4=sqrt5*phi^3>11/3, excluded.
At1 we retain it as foot type F below. Zero own contacts cannot close
an own edge, since its outgoing derivative would vanish.

Consequently every surviving common3/4 wall is an isolated5 at1,
with no spectators. Its ACTUAL foot is type F in Section4.

## 3. Complete common1/2 inventory

There is at most one right contact on2..8 and possibly a contact1.
For double contacts the exact finite-white restrictions are:

| Positions | Necessary restriction |
| --- | --- |
| 1,2 | bond at2 is3/4; bond at1 any finite native label;3,3 affine |
| 1,3 | deletion2 path B8 or A8: bond at3 is3, bond at1 is3/4 |
| 1,4 or1,7 | deletion2 E8: both3 |
| 1,5 or1,6 | deletion2 nonpositive already at both3 |
| 1,8 | deletion2 path: bond at8 is3, bond at1 is3/4 |

Both3 give an affine cycle. The4@1+3@8 row has the exact principal
affine F4 path2-1-c-8-7 with labels3,4,3,3. The sole4@1+3@3 row
has Q=0 (deletion2 own cost13/14) and is foot type H below.

For a single3, positions3..6 are forbidden by the same principal
affine E8/E7 trees as in Section2. At1 or8 the own cost is4/9 and
the Q bounds3/4 or9/16 make A^2<1/4. At7, Q<1/4 leaves no foreign
component and negative height. At2, own cost7/9 and Q<9/16 permit
at most two foreign components. One contradicts q<A^2; with two,
both first bonds must be3, and c,old2,old1,old3 plus the two foreign
contacted nodes form exact affine two-fork D5. Thus no single3 remains.

A single4 is affine internally. At1, Q<1/2 permits at most one
foreign component whose q>A^2=Q-1/9; at8,Q<1/8 permits none and
height is negative. A single label>=5 is possible only at1: the
retained A7 admits no such finite extension. Label6 is affine3,6.
For native>=7, Q<1-beta1^2<1/4, hence Q=0. For5, Q=0 is type F;
every Q>0 has Q>=1/4 and is retained in type L below.

The surviving end-common walls therefore have exactly three actual
foot types: L (adjacent contacts, isolated>=7, or isolated5 with
Q>0), H (4@1+3@3,Q0), and F (isolated5,Q0). No foreign branch has
been removed by assuming the whole star is finite.

## 4. Exact actual feet and uniform constants

All feet below are positive combinations of actual normals, have
positive time coefficient A, nonpositive products with EVERY actual
wall, and negative square. Thus their normalized future vectors are
ACTUAL points of the polytope. Write each as xi=A p-U, U purely
in the A8 component, and put

    g=(A-sqrtk)/sqrt(U^2),   k=-xi^2>0.

Farthest comparison gives x_A.(U/sqrt(U^2))>=g. For each type let
B>0 be its raw-coordinate sum divided by sqrt2. Exact constants are

| Type | Raw support | g^2 lower bound | U^2/B^2 lower bound | correlation with u7+u8 |
| --- | --- | --- | --- | --- |
| L | first3 | >9/20 | >4/9 | <61/100 |
| H | first4 | >53/100 | >1/3 | <3/4 |
| F | first5 | >11/20 | >2/9 | <9/10 |

**L.** Clear all spectators positively as in the A7 source. Put
b=beta1,d=beta2,B=b+2d,h=4(b^2+bd+d^2)/3+Q-1,z0=h/B^2.
The same explicit positive A2 dual foot has

    0<z0<1/3,  z0>8/81,
    k=B^2(z0-3z0^2),
    A^2=B^2(z0+4/9),  U^2=B^2(4/9+3z0^2).

The native minima and positivity proofs are unchanged from Section5
of the A7 source. The bound g^2>9/20 follows because
2523z0^2-720z0+484/9 is strictly positive, its discriminant being
-72976/3. The raw sum remains sqrt2 B. Consequently for W=u7+u8,

    rho^2=(U.W)^2/U^2
      =(2/sqrt63+1/6)^2/(4/9+3z0^2)
      <48114/129500<3721/10000.

Indeed (2/sqrt63+1/6)^2=23/252+2/(3sqrt63)<22/125,
using sqrt63>79/10. The final cross-product difference is729500.

**H.** Put s=sqrt2 and use the explicit positive foot

    xi=c+(s-1/2)n1+(s-1)n2+(s/2)n3,
    k=s-5/4,  A^2=(8+6s)/9,  U^2=(77-12s)/36,
    B=(3+s)/2.

Its full products and feasibility are the displayed profile calculation
in `compact_universal_a6_actual_projection_three_sevenths_20261020.md`,
SHA256 `022c0cc77512576dfc7bed74be372a94956618498c9ed74508aff53c40a50ae2`.
Here A^2>9/5, k<1/6, so A^2/k>54/5>23409/2209, giving g^2>53/100.
Also U^2/B^2=(77-12s)/(9(11+6s))>1/3 because22>15s.
The far correlation square is less than(22/125)/(1/3)=66/125<9/16.

**F.** Use the explicit positive isolated5 foot from Section3 of
the A7 source, now with N9. It has

    k=34phi-55,  A^2=(4phi-5)/9,
    U^2=(490-302phi)/9,  B=phi/2.

For g^2>11/20 it suffices that A^2/k>961/81. The elementary bounds
809/500<phi<40451/25000 give A^2>184/1125 and k<167/12500;
their ratio exceeds961/81 by the cross-product margin5752125.
Both phi bounds follow by squaring their positive sqrt5 equivalents.
Further U^2/B^2>2/9 is equivalent to979>605phi, which follows from
phi<16181/10000. The far correlation square is less than
(22/125)/(2/9)=99/125<81/100. This proves the table.

## 5. Coupling to the other end

After a numbering reversal if needed, the forced common pair in
Section1 supplies one of the left feet L,H,F. Put W=u7+u8; then
W^2=2+sqrt7/2<10/3. If the actual opposites7/8 are DISTINCT, their
pair law gives x_A.W>41/40. Each foot contradicts this:

* For L, the two-direction Cauchy expression
  rho*v+sqrt((4/7-v^2)(10/3-rho^2)), with v>=sqrt(9/20), is
  decreasing in v and increasing in rho on the indicated domain.
  Its upper endpoint is less than41/100+3/5=101/100. For the second
  term, its square is1510229/4200000<9/25; the first inequality
  follows from33489/200000<1681/10000.
* For H, the simpler separate bounds give
  x_A.W<(3/4)sqrt(4/7)+sqrt((4/7-53/100)*10/3)
  <57/100+19/50=19/20.
* For F, they give
  x_A.W<(9/10)sqrt(4/7)+sqrt((4/7-11/20)*10/3)
  <7/10+27/100=97/100.

Every displayed bound is strictly below41/40. These comparisons
use only positive rational squaring, not numerical optimization.

If opposites7/8 are COMMON, Section3 supplies a right foot of type
L,H or F. Normal reversal gives its raw sum -sqrt2 B'. For every
pair except F/F the two raw supports are disjoint; hence

    U.U'=2BB'/9,
    corr(U,U')=(2/9)/sqrt((U^2/B^2)(U'^2/B'^2)).

The table constants give these complete sum-direction Cauchy seals:

| Foot pair | Correlation upper bound | Lower bound for x_A.(unitU+unitU') | Resulting projection lower bound |
| --- | --- | --- | --- |
| L/L | 1/2 | 2sqrt(9/20) | >3/5 |
| L/H | 3/5 | 139/100 | >19321/32000>4/7 |
| L/F | 5/7 | 141/100 | >139167/240000>4/7 |
| H/H | 2/3 | 36/25 | >3888/6250>4/7 |
| H/F | 5/6 | 73/50 | >15987/27500>4/7 |

For example L/F uses corr<1/sqrt2<5/7, and H/F uses
corr<sqrt(2/3)<5/6. The scalar lower bounds follow from
sqrt(9/20)>67/100, sqrt(53/100)>18/25, sqrt(11/20)>37/50.

Finally F/F means the two actual walls are isolated end5 normals
with no spectators. Their exact product is

    c.c'=-A^2+(phi^2/4)C_18=(11-7phi)/18 in(-1/2,0).

This violates the native Coxeter gap for two distinct acute actual
walls. The signs follow from8/5<phi<2. Thus F/F is impossible too.
Every far-end alternative is now excluded, proving the theorem.

## Scope and contribution

This source extends the A7 coupled-core method, retaining the correct
positive mean correlation after reversing the normal basis. The new
three-foot inventory and its finite mixed table were derived by
finite57_research; root independently supplied the general-N core
formulas and proposed the A8 extension. Geometry/feasibility remain
written inputs, not claims of a machine certificate. Universal A9+
projection9/14 and any new ambient dimension ceiling remain open.
