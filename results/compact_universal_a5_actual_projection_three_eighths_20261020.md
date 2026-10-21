# Universal A5 actual projection above 3/8

20 October 2026. Complete paper candidate, awaiting independent review.
No scientific execution, enumeration job, optimization, arithmeticity,
Lean or TeX. Accepted sources are unchanged.

**Theorem.** At a vertex p farthest from a strict interior point of a
compact Coxeter polytope, every WHOLE A5 component satisfies

    ||x_A||^2>3/8>5/14,

with arbitrary other finite components. This is an actual farthest
projection estimate, with the no-common-opposite theorem used below;
it is not asserted for arbitrary edge-cover vectors.

Thus A5 joins the components with ordinary rank/14 cost. In dimension
n>=15, the existing structural alternative strengthens to: an A
component of rank>=6, a B component of rank>=5, or one dominant D_r
with r>=max(10,n-5) and spectator total rank<=5. This is a forest
restriction, not an unconditional dimension ceiling.

## 1. Actual tools, costs and the central pair reduction

Use the accepted actual reciprocal-shadow and farthest
orthogonal-component opposite inputs, respectively
`compact_actual_wall_reciprocal_shadow_20261020.md`, SHA256
`c913ec3bd8e0f26a8c679a1ede2b98b251ad3f36efeae60486213ddfa3e8badc`,
and `compact_farthest_vertex_orthogonal_component_opposite_obstruction_20261020.md`,
SHA256 `4216d06c2e03fd336898b8dac289b9969fa096c8a1d08bd67ddeb8ab5ec2f23d`.
Positive acute projection uses
`compact_farthest_vertex_mixed_replacement_inequalities_20261020.md`, SHA256
`3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84`.
All paths are under `results/`.

In unit path order the A5 inverse Gram is

    C=(1/3)*[[5,4,3,2,1],[4,8,6,4,2],[3,6,9,6,3],
             [2,4,6,8,4],[1,2,3,4,5]].

Put u_i=omega_i/sqrt(C_ii), h_i=x_A.u_i, and F(h)=2h/(1+h^2).
At the actual farthest point 0<h_i<1. Distinct actual opposites obey
`F(h_i)F(h_j)>=u_i.u_j`. On(0,1), F is increasing and log(F) is
strictly concave, so two directions of correlation1/sqrt2 have

    h_i+h_j>109/100.                               (1)

For an explicit seal, F((109/100)/2)^2<707/1000<1/sqrt2.
The first comparison has positive numerator2028179827 after clearing
denominators, and the second follows from2*707^2<1000^2.

Here u2.u3=u3.u4=1/sqrt2, u2.u4=1/2. If BOTH central pairs2/3
and3/4 have distinct opposites, (1) gives

    x_A.(u2+2u3+u4)>109/50,
    ||u2+2u3+u4||^2=7+4sqrt2<443/35.

The last bound uses sqrt2<99/70. Cauchy contradicts squared norm
at most3/8, because

    (109/50)^2-(3/8)*(443/35)=209/35000>0.          (2)

Consequently under the contrary projection bound a wall c closes
both2/3 or both3/4. Reverse the A5 numbering if necessary and assume
it closes2/3. This relabeling uses no ambient lattice symmetry.

Write c=A p-v, A>0, and let Q be its total spectator cost. Every
spectator component is retained in each actual finite replacement.
It has at most one c contact (two give a nonpositive cycle), with
cost q>=1/4. There is no cost in(1/4,1/3): a nontrivial component
has inverse diagonal at least4/3, while the next A1 native square
after1/4 is1/2. If Q has a contacted component of cost q>=A^2,
its contacted reciprocal shadow is at least one. Equality is
impossible compactly; strict inequality makes c also its actual
opposite, contrary to the orthogonal-component theorem. Therefore

    EVERY contacted spectator component must have q<A^2.    (3)

## 2. Complete multiple-contact inventory

Deleting2 retains A1 on1 and A3 on3,4,5. Deleting3 retains A2 on1,2
and A2 on4,5. Two c contacts within either side{1,2} or{3,4,5} give
a cycle retained by one of these positive-definite replacements.
Hence there are at most two A5 contacts; if two, they are at L<=2
and R>=3. Every individual bond is finite and native, because each
old node is retained by at least one actual replacement.

The SIX position pairs have the following complete possibilities,
from the indicated actual finite white diagram. Spherical finite
classification is an explicit classical input.

| L,R | Necessary finite white profile | Full-star consequence |
| --- | --- | --- |
| 1,3 | deletion2 is path p,q,3,3: q=3,p=3 or4 | p3 is affine cycle; p4 treated below |
| 1,4 | deletion2 is D5: p=q=3 | affine cycle |
| 1,5 | deletion3 is path3,p,q,3: p=q=3 | affine cycle |
| 2,3 | p,q each3,4 or5 | p=q3 affine triangle; others below |
| 2,4 | deletion3 is path3,p,q,3: p=q=3 | affine cycle |
| 2,5 | deletion3 is path3,p,q,3: p=q=3 | affine cycle |

Rank-five paths are A5 or B5, with a4 only at a path end; the
branched rank-five case is D5, all3. This is why no H/F exception,
internal4, extra white contact or spectator bridge enlarges the list.
All exact affine principal diagrams are forbidden by compactness.

For 4@1+3@3, the deletion2 white cost is1/2+3/8=7/8, so Q=0.
The actual direction3 then has

    cosh L3=3/2+1/sqrt2=(3+sqrt2)/2>11/5.

For contacts p at2,q at3, put b=cos(pi/p), d=cos(pi/q). The actual
deletion3 formula is

    cosh L3=(2b+3d)/sqrt(3*(1-4b^2/3-Q)).

It increases with either contact and with Q. Excluding the affine
3,3 case, its minimum among the native pairs is at3,4 (the reverse
4,3 is larger), and is again(3+sqrt2)/2. The strict threshold follows
from sqrt2>7/5. Thus every nonaffine double profile has
`tanh^2(L3/2)>3/8`, already contradicting the assumed projection.

## 3. Complete single-contact inventory and spectators

One A5 contact has label3,4 or5: a larger finite label retains a
nonpositive3,6 path in one endpoint. A contact at4 is at the middle
of the right A3 after deleting2, where only label3 is finite. There
is no zero-contact case because its A5 outgoing derivative is zero.

### Single3

The following table gives the own cost and a sufficient strict upper
bound on Q from an actual endpoint:

| Contact j | Own cost | Q upper bound | Height A^2 |
| --- | ---: | ---: | --- |
| 1 | 5/12 | 2/3 | Q-7/12<1/12 |
| 2 | 2/3 | 2/3 | Q-1/3 |
| 3 | 3/4 | 5/8 | Q-1/4 |
| 4 | 2/3 | 1/2 | Q-1/3 |
| 5 | 5/12 | 5/8 | Q-7/12<1/24 |

Rows1,5 require foreign contacts for positive height, but each cost
at least1/4 exceeds that height, contrary to(3). In row4 there is
at most one contacted component; its cost Q>A^2 again contradicts(3).
In row2 there are at most two. One is impossible; with two, (3)
forces both costs greater than1/3, hence Q>2/3, impossible. In row3,
two contacted components would both have cost strictly greater than
1/4; the cost gap then makes each at least1/3, giving Q>=2/3>5/8.
One component is again impossible. Thus all single3 profiles vanish.

### Single4

An internal contact2,3 or4 gives the exact affine star3,3,4. For an
end contact1 or5 the own cost is5/6 and the actual white condition
gives Q<1/3 (at5 the sharper bound is1/4). Thus A^2=Q-1/6<1/6.
Its positivity requires a foreign contact, whose cost at least1/4
contradicts(3). All single4 profiles are excluded.

### Single5

The middle-white contact at4 is not finite. At1 or2 the retained
white cost is phi^2/3>3/4; at3 or5 it is3phi^2/8>3/4. Thus Q=0
in all remaining cases, where phi=(1+sqrt5)/2. Their actual strong
direction and length are:

| Contact j | Direction k | cosh L_k |
| --- | ---: | --- |
| 1 | 4 | phi^3/sqrt3 |
| 2 | 3 | phi^2 |
| 3 | 2 | sqrt3*phi^3 |
| 5 | 2 | phi^3/sqrt3 |

For row1, deletion4 plus c is H4+A1 and is positive definite;
positive acute projection proves that the ALREADY actual c is the
actual opposite4. This is not a new wall inferred from a finite Gram.
All other listed directions are among the given actual2/3 endpoints.
The formulas follow directly from C and the white Schur residuals;
`8-3phi^2=5-3phi=phi^-4` gives the H4 denominator.

Each entry exceeds11/5. Indeed phi^2>5/2 and
phi^3/sqrt3>4/(7/4)=16/7>11/5; the latter uses phi>3/2 and
sqrt3<7/4. Thus the actual squared half tangent is greater than3/8.

This exhausts every common2/3 profile with arbitrary spectator costs.
Together with the central-pair reduction, it proves the theorem.

## 4. Forest consequence and boundary

The reviewed ordinary-component/exception accounting in
`compact_dn_uniform_energy_and_chain_relaxation_20261020.md`, SHA256
`1bd7c4eb372083973ca059b34121cdf0e3f575945a1552d771d7b09c84617804`,
and the strengthened actual-D theorem
`compact_dn_actual_projection_nine_fourteenths_20261020.md`, SHA256
`99bde3780f9d720d9809f5d8bea6795e4a77fcaacf8ede5871f0f831fb7823e9`,
now apply with A5 added to the ordinary classes. The other inputs,
including universal small-A bounds, remain unchanged.

If n>=15, no A>=6 and no B>=5, absence of D>=9 leaves only ordinary
rank/14 costs except at most one short A1 or one isolated rank-two7
defect of cost>9/100; those two kinds cannot coexist. With no defect,
rank n gives more than n/14>1. With one short A1, ordinary rank n-1
is at least14; with one rank-two defect, ordinary rank n-2 is at
least13 and total exceeds13/14+9/100>1. Thus a D>=9 exists.

Uniqueness and the actual9/14 bound, with the same three spectator
cases, force spectator total rank<=5 and r>=max(10,n-5). This proves
the stated alternative while retaining all old exceptional profiles.
A6--A8 rates and a reusable large-A actual-energy bound are still
separate obligations; no ambient dimension ceiling follows here.
