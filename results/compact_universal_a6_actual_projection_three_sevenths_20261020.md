# Universal A6 actual projection above 3/7

20 October 2026. Complete paper candidate, awaiting independent review.
No scientific execution, arithmeticity, enumeration, Lean or TeX.
All accepted and previously submitted sources are unchanged.

**Theorem.** At a vertex p farthest from a strict interior point of a
compact Coxeter polytope, every A6 component satisfies

    ||x_A||^2>3/7=6/14,

with arbitrary other finite components. The estimate uses an actual
feasible point of the polytope and is not asserted for arbitrary
edge-cover vectors. There is no ambient dimension-bound claim.

## 1. Normalization and a forced common pair

Use outward unit normals n_i=(e_i-e_(i+1))/sqrt2, 1<=i<=6, in the
sum-zero subspace of R7. Their inverse Gram is

    C_ij=2 min(i,j)(7-max(i,j))/7.

Write the strict interior point as z=E(p-x), E>0, and put
u_i=omega_i/sqrt(C_ii), h_i=x_A.u_i. The actual opposite of edge i
gives h_i>=tanh(L_i/2)>0. Actual reciprocal shadows and distinct-wall
acuteness imply

    F(h_i)F(h_j)>=u_i.u_j,   F(h)=2h/(1+h^2),

whenever the two actual opposites are distinct. These are the accepted
interfaces in `compact_actual_wall_reciprocal_shadow_20261020.md`, SHA256
`c913ec3bd8e0f26a8c679a1ede2b98b251ad3f36efeae60486213ddfa3e8badc`.
Positive acute projection uses
`compact_farthest_vertex_mixed_replacement_inequalities_20261020.md`, SHA256
`3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84`.
All cited paths are under `results/`.

On (0,1), log F is strictly concave. Thus for correlation kappa a
lower pair-sum seal s follows from F(s/2)^2<kappa. Here

    u2.u3=u4.u5=sqrt(8/15),   u3.u4=3/4.

The outer pairs have sum>2249/2000. Indeed
F(2249/4000)^2<7301/10000<sqrt(8/15): the numerator of the first
comparison after cross multiplication is430464052923301, and
15*7301^2=799569015<800000000. The middle pair has sum>11547/10000,
since its sharp lower sum is2/sqrt3 and
3*11547^2=399999627<400000000.

If pairs2/3,3/4,4/5 all have distinct opposites, put
w=u2+2u3+2u4+u5. Then

    x_A.w>34037/10000,
    w^2=84/5+8sqrt(8/15)+8sqrt(3/10)<33781/1250.

The last seal uses sqrt(8/15)<7303/10000 and
sqrt(3/10)<5478/10000. Cauchy contradicts ||x_A||^2<=3/7, since

    (34037/10000)^2-(3/7)(33781/1250)
      =2181583/700000000>0.

Consequently a common pair is2/3,3/4 or4/5. Reversing the A6
numbering, without assuming any ambient symmetry, leaves only
common2/3 or common3/4 to classify.

## 2. Spectator bookkeeping and complete double-contact inventory

Let c=A p-v be such an actual common opposite, A>0. Put
beta_i=-c.n_i>=0 and let Q be its total spectator Schur cost. Every
spectator component is retained in both actual positive-definite
replacements. It has at most one c contact, since otherwise the
retained finite diagram contains a cycle. Each contacted component
has cost q>=1/4, and no cost lies strictly between1/4 and1/3.

Moreover every such cost obeys q<A^2. If q>A^2, deleting its contacted
node instead gives another positive-definite full replacement, so
positive acute projection makes c an actual opposite in that
orthogonal component as well. This contradicts
`compact_farthest_vertex_orthogonal_component_opposite_obstruction_20261020.md`,
SHA256 `4216d06c2e03fd336898b8dac289b9969fa096c8a1d08bd67ddeb8ab5ec2f23d`.
Equality gives the same positive feasible null projection and is
impossible compactly. Thus

    every contacted spectator cost q is strictly less than A^2.  (1)

There is no zero A-contact case: its derivative on either A edge
vanishes, so it cannot be that edge's finite opposite. At most one
contact lies on each side of the common consecutive pair; two on
one side give a cycle retained in the other replacement. Each bond
is finite and native. A bond>=6 with an adjacent old3 is nonpositive
already in a retained rank-three path, so only3,4,5 occur.

Classical finite spherical classification now gives the following
complete DOUBLE-contact lists. A listed affine cycle is a principal
all3 cycle in the full actual diagram and is forbidden compactly.

| Common pair | Contact positions | Necessary finite white restriction |
| --- | --- | --- |
| 2/3 | 1,3 | deletion2 path p,q,3,3,3: q3, p3 or4 |
| 2/3 | 1,4 or1,5 | deletion2 E6: both3 |
| 2/3 | 1,6 | deletion3 rank-six path: both3 |
| 2/3 | 2,3 | p at2 is3/4/5; q at3 is3/4 |
| 2/3 | 2,4 or2,6 | deletion3 rank-six path: both3 |
| 2/3 | 2,5 | deletion3 D6: both3 |
| 3/4 | 3,4 | both bonds3/4/5 |
| 3/4 | all other eight pairs L<=3,R>=4 | deletion3 or4 is a rank-six path or D6: both3 |

Every all3 row gives an affine cycle. Rank-six paths have a4 only
at an end, while all displayed variable bonds are internal; E6 and
D6 are all3. This explicitly retains the E6 exception rather than
assuming every finite white component is A/B.

For common2/3, positions2,3 with at least one bond>=4 are strong:
if q>=4, monotonicity of the deletion2 formula gives

    cosh L2>=(5+4sqrt2)/sqrt7>5/2;

if p>=4, the deletion3 formula gives

    cosh L3>=(4sqrt2+6)/sqrt14>5/2.

For common3/4, the adjacent double profile is affine at3,3; otherwise
one of the reversed endpoint formulas gives

    max(cosh L3,cosh L4)>=(12+9sqrt2)/sqrt42>5/2.

Each formula increases with either bond and Q. All these bounds force
the corresponding squared half tangent>3/7, since
(cosh L-1)/(cosh L+1)>3/7 exactly when cosh L>5/2.

Only4@1+3@3 remains. Its deletion2 white cost is9/10, hence Q=0.
The NEW actual-point argument for this row follows in Section4.

## 3. Complete single-contact inventory

For a single3 contact the own cost is C_jj/4. The table records a
sufficient strict Q upper bound from the two actual endpoints.

| Common pair | j | Own cost | Q bound | A^2 |
| --- | --- | ---: | ---: | --- |
| 2/3 | 1 | 3/7 | <2/3 | Q-4/7 |
| 2/3 | 2 | 5/7 | <2/3 | Q-2/7 |
| 2/3 | 3 | 6/7 | <3/5 | Q-1/7 |
| 2/3 | 4 | 6/7 | <2/5 | Q-1/7 |
| 2/3 | 5 | 5/7 | <2/5 | Q-2/7 |
| 2/3 | 6 | 3/7 | <3/5 | Q-4/7 |
| 3/4 | 1 or6 | 3/7 | <5/8 | Q-4/7 |
| 3/4 | 2 or5 | 5/7 | <1/2 | Q-2/7 |
| 3/4 | 3 or4 | 6/7 | <5/8 | Q-1/7 |

Rows with A^2<1/4 cannot have any foreign component by(1), hence
cannot have positive height. For common2/3 at j2 there are at most
two foreign components. One contradicts(1); with two, each other's
cost exceeds2/7, so the cost gap forces Q>=2/3, also impossible.

The only remaining single3 rows are j3 in common2/3, or j3/4 in
common3/4. Their Q<5/8 allows at most two foreign components. One
again contradicts(1). With two, both first bonds are3, since a4
already costs1/2 and the other component costs at least1/4. The
six nodes consisting of c, the contacted internal old node j, its
two old neighbors, and the two foreign contacted nodes are the
exact affine two-fork D5 diagram. Compactness forbids it. Thus all
single3 cases are excluded, including arbitrary spectator ranks.

A single4 at an internal old node gives the principal affine3,3,4
star and is forbidden. At ends1 or6, common3/4 gives Q<1/4; common2/3
gives Q<1/3 at1 and Q<1/5 at6. The height is A^2=Q-1/7, always
less than1/4. Therefore(1) excludes all such profiles.

A single5 with common2/3 can occur only at1 or2: attachment to the
retained A4 on3..6 would require a finite rank-five group with a5,
which does not exist. With common3/4 it can occur only at1,3,4,6:
attachment to the middle of a retained A3 is nonpositive. The white
H3/H4 cost is at least phi^2/3>3/4, so Q=0 in every case. Already
actual c then has the following actual finite endpoint, by positive
acute projection; these are not new walls inferred from Gram data.

| Common pair | Contact j | Actual direction k | cosh L_k |
| --- | --- | --- | --- |
| 2/3 | 1 | 4 | 3phi^3/sqrt21 |
| 2/3 | 2 | 3 | 4phi^2/sqrt14 |
| 3/4 | 1 or6 | 4 or3 | 3phi^3/sqrt21 |
| 3/4 | 3 or4 | 4 or3 | 9phi^3/sqrt21 |

Here phi=(1+sqrt5)/2. Every entry exceeds5/2: using phi>8/5,
sqrt21<23/5 and sqrt14<15/4, the two smallest entries exceed63/23
and208/75 respectively. This completes the single-contact list.

## 4. The remaining double profile: an actual Lanner dual foot

For4@1+3@3 one has beta1=1/sqrt2, beta3=1/2, all other contacts
zero, including every spectator. With s=sqrt2, define

    xi=c+(s-1/2)n1+(s-1)n2+(s/2)n3.

All coefficients are positive. Direct unit-Gram multiplication gives

    xi.c=xi.n1=xi.n3=0,
    xi.n2=(s-3)/4<0,
    xi^2=5/4-s=-k,   k=s-5/4>0.

Every other actual wall has nonpositive product with this positive
combination, by all-wall acuteness. Its time coefficient A>0 makes
xi future. Thus xi/sqrtk is an ACTUAL point of the full compact
polytope, even with arbitrary spectator components.

Write xi=A p-U. The old inverse Gram gives

    A^2=(5+4s)/7,
    U=(3s/4)(e1+e2)+(1/2)(e3+e4)-((1+3/s)/7)*1_7,
    U^2=A^2-k.

The farthest-point comparison with xi/sqrtk yields
x_A.U>=A-sqrtk. Therefore

    ||x_A||^2 >= (A-sqrtk)^2/(A^2-k)
                =(A-sqrtk)/(A+sqrtk)>1/2.

For the strict last inequality A^2>9k is enough, and

    A^2-9k=(335-236sqrt2)/28>0,
    335^2=112225>2*236^2=111392.

This contradicts the contrary3/7 projection bound. Together with
Sections1--3 it proves the theorem. The same positive-foot row has
A^2=3/2+sqrt2-(11/2+3sqrt2)/(r+1), increasing with r, so its >1/2
actual projection conclusion also holds whenever that exact profile
occurs in an A_r component with r>=6. This is a profile-specific
extension, not a complete A_r theorem.

## Scope

A6 now has its ordinary rank/14 cost, conditional on independent
review of this paper candidate. A5's separate submitted >3/8 source
is unchanged and remains a separate review item. A7--A8 ordinary
rates and a reusable >=9 large-A actual9/14 bound are still open.
No new ambient dimension ceiling or edge-only energy claim follows.
