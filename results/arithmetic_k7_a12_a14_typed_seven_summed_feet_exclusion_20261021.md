# K7: A12/A14 couplings with either typed7 defect

21 October 2026. Complete MANUAL supporting proof, pending full root and
outside-author review. No scientific, numerical, symbolic, enumeration,
Lean or TeX job; no manuscript, accepted source, index or ceiling change.
This is one companion argument in the existing arithmetic paper family.

**Theorem.** At a farthest vertex of a compact arithmetic first-type
Coxeter polytope containing either isolated rank-two7 defect, an A12
component has energy greater than.703101 and an A14 component greater
than.7897732, with arbitrary further finite spectators. Consequently
the following rank18 forests are impossible:

* A14 plus either typed7 defect and an ordinary irreducible T2;
* A12 plus either typed7 defect and an ordinary irreducible T4;
* A12 plus either typed7 defect and two ordinary irreducible T2 components.

Here typed7 means OLD A2/NEW isolated7 or OLD I2(7)/NEW isolated3.
Their actual isolated opposite is part of the hypothesis. A bare old
I2(7) component is not asserted to be a defect. No new global ceiling
or complete rank18 frontier follows.

## 1. Inputs and the field forced by the actual defect

All paths below are under `results/`. Only their stated local geometric,
arithmetic-normalization and component-accounting interfaces are inputs.

| Input | SHA-256 | Used interface |
|---|---|---|
| `compact_universal_small_a_component_energy_20261020.md` | `801581c4b855289aca7982b2d0f4bb6f5849f151a78f40b86c9d403e84f415e6` | Actual OLD A2/NEW isolated7 profile |
| `compact_farthest_rank_four_forest_exclusion_20261020.md` | `68ad62a619603eb7e6c2acc5689c901dde3f9d0897a7c1580b3deb948f5832b6` | Both isolated typed7 profiles, their energy>.09 and compatibility |
| `arithmetic_uniform_classical_weak_helper_obstruction_20261021.md` | `535ac0d0fa92d28ccab2976fa147a4ccff0bd3bfcd1e04ca4e134c6ae6857905` | Field-valued normalization, cyclic integrality and lifting every minimal-field embedding |
| `compact_a19_a20_central_window_actual_six_sevenths_20261020.md` | `8075f92f4766d681cdef42040925586bcba76ebdd99c059ce5d75a706dc3ee0c` | Finite-white list and actual H4/A2/singleton formulas, not its larger-rank conclusion |
| `compact_an_large_components_arbitrary_spectator_transfer_20261020.md` | `8d3b999761be7dd46e7990973b3ec4d4a746d4d682b044d53695f4aae405e3e5` | Zero-own repair and positive clearing of entire spectator blocks |
| `compact_forest_ordinary_a_parameter_and_dominant_classical_frontier_20261020.md` | `dd6f8b0de2d7d0bc835be3eb4e320013aa40792407fb2a747f35bcf04385d5ad` | Ordinary T2>1/7, T4>2/7 and orthogonal component accounting |

The standard external hypotheses are compact acute simplicity, finite/
affine classification, the actual native-wall gap, and the first-type
cyclic-integrality/strict conjugate-positive-definiteness criterion.
No field census, parked refinement, degree-specific enumerator or
geometric Lean result is assumed.

Put B=4cos²(pi/7). The actual isolated defect supplies a three-wall
principal path[7,3]. Its contact square puts B in the actual totally
real ground field K. Normalize the central root to square2 and its
3-neighbor to square2, product−1. Rescale the7-neighbor over K to
product−1; its square is2/B. Thus this actual principal Gram, up to
permutation, is

    [ 2  -1  -1 ]
    [-1   2   0 ]
    [-1   0  2/B].

Its determinant is6/B−2<0 since B>3; its two-node principal blocks
are positive definite, so it is physically Lorentzian. This construction
uses field-valued rescaling, not adjoining sqrt B. The scalar B has
irreducible polynomial X³−5X²+6X−1, with no rational root. Consequently
K7=Q(B) is the real cubic7-field. If K were a proper extension of K7,
an additional real embedding of K fixing K7 would
be nonphysical but preserve the displayed Lorentzian Gram. This
contradicts strict conjugate positivity. Therefore the ACTUAL field
is exactly K=K7, for either typed7 defect and with all spectators.

In particular K contains no quadratic subfield. Native5,8,12 contacts
would put respectively sqrt5,sqrt2,sqrt3 in K through their squared
products, so are impossible. A simple3-helper/own4 cap would put
sqrt2 itself in K by the closed-path normalization of535; it too
is impossible. These are restrictions on actual walls, not on a
normalized projected foot.

## 2. Complete small-rank geometric inventory

Set N=13 for A12 or N=15 for A14. Use UNIT old normals and duals

    D_i=2i(N-i)/N, C_ij=2min(i,j)(N-max(i,j))/N,
    h_i=x_A.omega_i, h0=hN=0, ||x||²<1,
    e_A=(1/2)sum_(i=0..N-1)(h_(i+1)-h_i)².

All pressures use the ONE original farthest vector. Call5<=i<=N-5
weak when h_i<=q_i-.61, q_i=sqrt(D_i). The minimum D is80/13 or
20/3; throughout both windows `12/5<q_i<3`.

At an actual endpoint the retained A paths have ranks at least4.
A finite extension by the opposite c has exactly the following white
roles: no contact; one simple3 at an admissible end/D/E position;
one end4; or two end3 contacts, one to each path. A path plus c has
rank at least5, excluding H/F exceptions. The joined rank is12 or14,
still>8; the tree/two-long-arm/bond argument forces the same A-path
completion for the two-white role. Every other ordinary role has the
same explicit local finite/affine principal test as in8075.

The zero-own-contact repair remains complete at these ranks. Global
internal3 positions contain affine arms1,2,5 or1,3,3, internal4 has
the squared-weight affine star1/2,1/4,1/4, and two3 contacts make an
affine cycle. Only the far-end3, node2 and far-end4 cases remain.
Their retained side parameter is m>=5; the respective height bound
<1/(2m), foreign-cost bound Q<2/m and bound Q<1/m are unchanged.
For the first two, a contacted whole spectator gives a reciprocal
shadow above one in its orthogonal component, the actual-opposite
contradiction of8d3; absence of foreign contact gives negative height.
The last has Q<1/5 and hence Q0 and negative height. Thus own contact
gamma is positive, without a whole-star-finiteness assumption.

The one scalar change below q=2.5 is explicit. The strong-edge target

    R(q)=100/61-(q-.61)/[q(2q-.61)]

increases already for q>=12/5: its derivative has numerator
`2q²-4(.61)q+.61²>0` and positive denominator. Therefore
`R(q)<R(3)<3/2`, since `3/2-R(3)=1669/197274>0`.
All ordinary nonisolated role bounds use m>=5 and exceed38/25, hence
are still strong. In the two-white role the former q>=2.5 seal is
replaced directly by `q/2+1/q>97/60>38/25`, valid for q>12/5.
The D/E, near3-own>=5 and end4 role estimates have no additional
larger-rank requirement. Own3/4/6 affine principals are retained.

For either exceptional role, a positive spectator projection Q costs
at least1/4. The exact entire-spectator clearing gives far-helper
ratio `>227/(25sqrt35)` or near4 ratio square `>121/35`. Both exceed
the present threshold. Thus every weak nonisolated exception has Q0.
Its complete list is a simple3-helper/own4 near cap, or

    v_A=gamma omega_i+omega_e/2, e=1 orN-1,
    own native5/native>=7/dotted.                    (1)

Every near cap is excluded by sqrt2 not in K7; no mate or foreign
contact has been discarded. Native5 is absent too, so a remaining
helper has t=2gamma>9/5. Its exact edge ratio satisfies
`ratio²>=4gamma-2`; R<3/2 implies t<17/8<54/25. Both helper ends
remain allowed and no adjacency statement is needed.

## 3. No remaining weak end-helper over K7

Normalize(1) by its simple3 helper as in535. Then t is an algebraic
integer in K7. Every nonidentity embedding of Q(t), lifted to the
actual field, has strict positive Schur complement

    D_i t_sigma²+2C_ie t_sigma+D_e<4.

Reverse to i<=N/2. For negative t the cheapest helper cost is
`Q_-(t)=2it²+4t+2-2(it+1)²/N`; for positive t it is
`Q_+(t)=2it²+2-2(it-1)²/N`. On the relevant sign interval each
increases with i and N. At the smallest N13,i5,

    Q_-(-83/100)-4=69/1625>0,
    Q_+(1/2)-4=2/13>0.

The quadratics at zero are below4 and increase away from their
admissible intervals. Thus EVERY other conjugate lies in(-.83,.5).

Since K7 has degree3, t has degree1 or3. In degree3 its nonzero
integer norm has magnitude `<(54/25)(83/100)²<2`. If either small
conjugate were positive, the magnitude would be
`<(54/25)(1/2)(83/100)=.8964<1`, impossible. Both are therefore
negative, with positive absolute product1/t. The integral trace is
strictly greater than `9/5-2(.83)=.14>0`. By AM–GM it is at most
`t-2/sqrt(t)<1`, since `t(t-1)²<(54/25)(29/25)²<4`. An integer
cannot lie in(0,1). In degree1 the only integer in(1.8,2.16) is2;
its actual own parallel pair has a future feasible null sum and
contradicts compactness. All weak end-helpers are excluded.

## 4. Remaining isolated feet and continuous loss absorption

For a cap isolated WITHIN A, write
`c=A p-gamma omega_i-v_S`, Q=||v_S||², u=(1-Q)/gamma².
The actual endpoint retains the entire spectator forest, so0<=Q<1.
Isolated own3/4/6 and parallel have their affine/null principals;
native5/8/12 are absent in K7. Every other own contact has gamma>.9.

For Q>0, Q>=1/4 gives u<(.75)/(.81)<1; no foreign u>1 branch
has been silently omitted. The case u=1 is an actual future feasible
null ray and is impossible. For Q0 finite native7, the actual A2
negative core gives

    loss < (100/81)/(q+sqrt(D_i-100/81))
                            +sqrt(4/3-gamma^-2).

Here .9<gamma<.902. Even at N13, q>2.48 and
`sqrt(D_i-100/81)>=sqrt(5180/1053)>2.217`; hence the first term is
`<100000/380457<.264`. The second is `<sqrt(17/150)<.337` because
gamma²<.902²<50/61. The total is<.601<.61. For every Q0 finite
m>=9, the actual singleton gives loss `<a_i+tan(pi/m)`, where
`a_i=1/[q+sqrt(D_i-1)]`. The elementary estimate
`tan(pi/9)<1386/3727<3/8` follows from pi<22/7,
sin a<a and cos a>1-a²/2. Since a_i<.212, the loss is<.587<.61.
All finite native labels have now been accounted for.

Retain EVERY isolated case with0<u<1, including positive foreign Q.
The vector

    eta=c+v_S+gamma u n_i
        =gamma[sqrt(D_i-u)p-(omega_i-u n_i)]

uses positive combinations of actual old spectator normals. It has
eta.S=0, eta.c=0, eta.n_i=-gamma(1-u)<0,
`eta²=-gamma²u(1-u)` and positive time coefficient. All other wall
products are nonpositive by acuteness. Its normalization is an
actual full-polytope point, so farthest comparison gives

    x_A.(omega_i-u n_i)>=sqrt(D_i-u)-sqrt(u(1-u)).  (2)

The loss from q is at most `a_i u+sqrt(u(1-u))`. Use the following
rational bounds for the references in Section5:

| N | a_i<a | every positive beta_i>=b | radical bound r>sqrt(1+a²) |
|---|---:|---:|---:|
|13|.212|.0305|1.023|
|15|.202|.0245|1.021|

At N13 the denominator for a_i exceeds2.48+2.27=4.75; at N15 it
exceeds2.581+2.38=4.961. Positive squares at the displayed D minima
verify these estimates and a_i<a. For ALL0<u<1,

    a u+sqrt(u(1-u))-(b/2)(2u-u²)<.61.             (3)

Indeed on u<=1/2 its unadjusted bound is a/2+1/2, respectively
.606 or.601. On u>=1/2 use2u-u²>=3/4 and the continuous maximum
`(a+sqrt(1+a²))/2`. The adjusted bounds are respectively
`.6175-.0114375=.6060625` and `.6115-.0091875=.6023125`.
The radical comparisons are strict positive-square tests. The maximum
itself follows by expanding the circle identity, or by Cauchy on
`(2u-1,2sqrt(u(1-u)))`; it is not an input about actual cap data.

## 5. One reference and all simultaneous feet

Choose symmetric rational heights with linear end tails:

| N | central H_i to the midpoint | E_H=||w||² |
|---|---|---:|
|13|H5=1.870,H6=H7=1.931|.703101|
|15|H5=1.971,H6=2.073,H7=H8=2.122|.7897732|

The values are strictly below q_i-.61: positive squares compare
2.480²,2.541² with80/13,84/13, and
2.581²,2.683²,2.732² with20/3,36/5,112/15. The reference slopes
decrease, so beta_i=(2H_i-H_(i-1)-H_(i+1))/2 are nonnegative and
`w=sum beta_i omega_i`. They are positive only in the stated central
window. At13 their positive values are.1565,.0305 and reflections;
at15 they are.1461,.0265,.0245 and reflections, verifying b.
Direct summation gives

    E_13=1.870²/5+.061²=.703101,
    E_15=1.971²/5+.102²+.049²=.7897732.

Let S contain the retained isolated indices with positive beta, and set
`Z=w-sum_S beta_i u_i n_i`. At every other positive coefficient the
actual height is strong. Summing those constraints and(2) gives

    x_A.Z>=E_H-sum_S beta_i(loss_i-.61).

Since w.n_i=beta_i and all distinct old products are nonpositive,

    ||Z||²<=E_H-sum_S beta_i²(2u_i-u_i²).

This inequality handles ADJACENT feet as well: their cross terms have
the favorable nonpositive sign. By(3) and beta_i>=b, completing the
square now gives

    e_A>=2x_A.Z-||Z||²
       >=E_H+sum_S beta_i[beta_i(2u_i-u_i²)-2(loss_i-.61)]
        >E_H.

If S is empty strictness comes from an available strong floor;
otherwise each displayed summand is strictly positive. There are no
golden omissions or independent endpoint gains in this proof. The
spectators were positively cleared before the pressures were summed.

## 6. Consumers and exact scope

The typed7 block pays>.09. Ordinary T2 pays>1/7, and T4 or two
ordinary T2 components pay>2/7. These are orthogonal whole-component
projections of the SAME original vector; count the defect once.
For A14 the total is greater than
`.7897732+.09+1/7>1`. For either A12 alternative it is greater than
`.703101+.09+2/7>1`. This proves all stated exclusions.

The theorem is conditional on an ACTUAL typed7 defect, which forces
K7; it is not an arbitrary-field A12/A14 energy theorem or an intrinsic
stronger energy claim for an ordinary I2(7). It does not address I/J,
other rank18 families, whole A/B/D components or a new global ceiling.
All six papers, accepted bounds and prior source/evidence packets
remain unchanged. This new finite-white/field/summed-foot composition
requires full written root and outside-author review before acceptance.
