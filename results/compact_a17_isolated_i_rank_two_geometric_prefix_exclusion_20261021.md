# A17 energy and the isolated-I/rank-two remainder, without arithmetic

21 October 2026. Complete manual supporting proof, pending full root and
outside-author review. No scientific, numerical, symbolic, enumeration,
Lean or TeX job; no manuscript, frozen input, canonical index or accepted
ceiling is changed. Arithmeticity is NOT assumed.

**Theorem.** At a vertex farthest from a strict interior point of a
compact hyperbolic Coxeter polytope, a whole A17 component has actual
orthogonal projection energy greater than83/100, with arbitrary finite
spectators. Consequently it cannot coexist with an isolated short A1
of type I and an ordinary irreducible rank-two component T2.

The new issue is a complete N18 version of the geometric .61 inventory
and its adjacent-helper exclusion. In its tight reversed-helper case,
the larger-rank coarse native5/native7 floors do not suffice. The native
angle bounds below repair that step without an arithmetic restriction.
The argument permits arbitrarily many nonadjacent weak sites; there is
no field packing or at-most-one-golden-cap premise.

## 1. Precise inputs and normalization

All paths below are under `results/`. Their transitive actual-point and
classification hypotheses remain explicit. Only the indicated geometric
proofs are used, not their older rank-specific conclusions.

| Input | SHA-256 | Used interface |
|---|---|---|
| `compact_large_a_two_end_small_debit_coupling_20261021.md` | `43951a969029090579637fc98f645dac674271caac0ef3337347fe6a0002ee04` | General-N coordinates, helper-product identities, rational two-end reference and deletion proof |
| `compact_a19_a20_central_window_actual_six_sevenths_20261020.md` | `8075f92f4766d681cdef42040925586bcba76ebdd99c059ce5d75a706dc3ee0c` | Finite-white role classification and actual isolated-cap formulas before its6/7 contradiction |
| `compact_an_large_components_arbitrary_spectator_transfer_20261020.md` | `8d3b999761be7dd46e7990973b3ec4d4a746d4d682b044d53695f4aae405e3e5` | Zero-own-contact repair, entire-spectator clearing and positive actual points |
| `compact_a21_adjacent_weak_far_m_simple3_packet_20261020.md` | `1c4437c98d916d558081572dd2ecc12bc75eda1328a53633fe111120ec63dcb7` | General-N actual mixed replacement and Schur inverse/pressure formulas |
| `compact_an_projection_pressure_and_weak_cap_packing_review_20261020.md` | `8b0a59493f1fa5dc7081e5e70beff2839ac5918bfcd50e089e05ae9b5546e43b` | Far5 distance-one prohibition for N>=15, both helper ends; actual near4 mate |
| `compact_dimension20_general_forest_frontier_and_proper_classical_continuation_20261021.md` | `9f13de08c8bffd4483359ad754b34da8dc1411c67ebcd0dcc8f8c980cba60b61` | Universal isolated-I plus ordinary T2 energy>13/70 and the necessary forest scope |

Finite spherical classification, compact acute simplicity and the native
gap for distinct actual walls are the standard external geometric facts.
The listed new geometric lemmas are internal proved inputs. No arithmetic
criterion, minimal field, cyclic integrality or conjugate embedding is used.

Use unit old normals `n_i=(e_i-e_(i+1))/sqrt2`, N=18, and duals omega_i.
For the ONE original farthest vector x, `||x||²<1`, let x_A be its
orthogonal A17 projection, e=||x_A||², and put

    D_i=2i(18-i)/18, C_ij=2i(18-j)/18 (i<=j),
    h_i=x_A.omega_i, h0=h18=0,
    e=(1/2) sum_(i=0..17)(h_(i+1)-h_i)².

The actual h sequence is strictly concave. On5<=i<=13 call i weak if
`h_i<=q_i-.61`, q_i=sqrt(D_i), and otherwise strong. Throughout this
window

    D_i>=65/9, q_i>67/25, q_i<=3,
    C_(i,i+1)>=20/3 for adjacent central sites.          (1)

## 2. Complete fresh N18 weak inventory, with every spectator retained

At an actual endpoint of direction i the retained old paths are
A_(i-1) and A_(17-i), both of rank at least4. A finite component
containing its actual opposite c is a tree. Its complete white-contact
list is: no contact; one simple3 at an admissible end, D or E position;
one end4; or two end3 contacts, one on each path. A path plus c has
rank at least5, excluding H/F completions. If both paths are joined,
the connected rank17 exceeds8; the two long-arm/tree/bond argument
forces the same A-path completion as in the input. Thus shortening
the ambient A component has introduced no additional white role.

Write `c=A p-v_A-v_S`, Q=||v_S||², where S is the ENTIRE spectator
forest. Its actual finite endpoint retains S, so Q<1. Any contacted
spectator component has cost at least1/4. The zero-own-contact repair
uses only retained side parameter m>=5 and is still exhaustive here:
an internal simple3 has an affine star with arms1,2,5 or1,3,3; an
internal4 has the affine squared-weight star1/2,1/4,1/4; two simple3
contacts make an affine cycle. A far-end3 with zero own contact has
height square<1/(2m) and needs foreign cost>=1/4; the node2 role has
Q<2/m and foreign cost strictly larger than its height square; an
end4 has Q<1/m<=1/5 and hence Q0 and negative height. The first two
foreign cases force that same opposite in an orthogonal component,
by the actual orthogonal-opposite obstruction. Therefore own contact
gamma is positive, without assuming that the entire added star is finite.

Here are the fresh isolated-cap scalar seals. At D>=65/9,

    q+sqrt(D-1)>259/50, a=1/[q+sqrt(D-1)]<50/259<.194.

For example q>2687/1000 and sqrt(D-1)>2494/1000 follow by positive
squares at65/9 and56/9; their sum is>5.181>259/50.
Positive clearing of all S gives, with u=(1-Q)/gamma²,

    h_i>sqrt(D-1)-sqrt(u-1),                if u>1;
    h_i>sqrt(D-u)-sqrt(u(1-u)),             if0<u<1.

These are pure-A actual-point inequalities, not native-angle claims
about a projected facet. For Q>0,u>1, the native gap gives
`u-1<=5-3phi<73/500<(5/12)²`, so the loss is
`<50/259+5/12=1895/3108<61/100` (cross margin88).
For0<u<1 the loss is

    <a u+sqrt(u(1-u))<(.194+1.019)/2=.6065<.61,

using sqrt(1+.194²)<1.019. This includes isolated dotted caps.
The case u=1 is a nonzero future feasible null ray, impossible in a
compact polytope. For Q0 and finite native labels at least8, singleton
projection gives loss `<a+tan(pi/8)<50/259+5/12<.61`.
Native5 uses the actual H4 correction and
`q+sqrt(D-25/16)>19/4`, giving loss `<25/76+5/24<.61`.
For native7 the actual A2 negative core gives

    loss < (100/81)/(q+sqrt(D-100/81))+9/25
         <2500/10368+9/25<.61,                    (2)

since q>67/25 and sqrt(D-100/81)>=sqrt485/9>61/25.
Isolated own3/4/6 and parallel contacts have the same affine/null
principal obstructions. Thus EVERY isolated cap is strong, including
all foreign and continuous dotted cases.

For the remaining white roles the necessary strong-edge ratio is

    R(q)=100/61-(q-.61)/[q(2q-.61)].

It increases on q>=5/2 and at q=3 satisfies
`3/2-R(3)=1669/197274>0`. Hence R<3/2<761/500 on our entire window.
The complete ordinary-role inequalities of the inputs use only m>=5:
near3/own>=5 has ratio at least `(101/125+2/5)/sqrt(3/5)`;
D positions have at least `(101/125+1/5)/sqrt(2/5)`; E/end4/two-white
roles have their larger bounds. All exceed761/500, with the same
explicit affine principals disposing of the smaller own labels.
This accounts for every nonexceptional role in the finite list above.

Adding Q strengthens each edge bound. The only exceptional roles are
far-end3 and near4/simple3. If Q>0 their cleared ratios are respectively
`>227/(25sqrt35)>761/500` and square `>121/35`, by the positive
entire-spectator formulas valid for m>=5. Thus EVERY weak cap has Q0.
Near4 with helper-side parameter m>=8 is strong by the old761/500
estimate, and m7 has ratio>3/2>R. The complete remaining weak list is

    near4: v=omega_i/sqrt2+omega_(i-1)/2, i5 or6,
           or its reflected role;
    far helper: v=gamma omega_i+omega_e/2,
           e1 or17, Q0, own native5/native>=7/dotted.   (3)

The far helper has sigma=gamma+1/(2m), delta=(m+1)/(2m), m>=5.
The identity `sigma²/delta-(4gamma-2)=(gamma-1-1/(2m))²/delta`
and R<3/2 imply gamma<17/16<27/25. We will use the safe uniform
range `f=.809<gamma<g=27/25`; the lower bound is phi/2>.809.
No algebraic restriction is placed on dotted gamma.

## 3. Literal helpers are nonadjacent at N18

Let adjacent central indices be i,j=i+1, with distinct actual caps
as in the far-helper row. Put A=gamma_i, B=gamma_j and R=c_i.c_j.
The general-N AM–GM identities give

    same helper: R>1-g²>-1/2;
    outward opposite helpers: R>=1/2-AB+(A-1)(B-1)/18;
    reversed opposite helpers: R>=1/2-AB+(A+1)(B+1)/18.

Both latter bounds exceed `1/2-g²-(1-f)(g-1)/18>-.7>-1/sqrt2`.
Thus the actual native gap leaves R0 for a common helper and R0/-1/2
for opposite helpers. This includes dotted/mixed caps and BOTH opposite
orientations, rather than only a symmetric assignment.

For the outward orientation the projected parameters eta_i,eta_j<4/5
give the exact positive determinant lower bound

    (D_i-4/5)(D_j-4/5)-C_ij²
      =(2/5)C_ij-24/25+8/90>=404/225>0.

As in the actual two-cap proof, this rules out R0. For the reversed
orientation put l=i,s=18-j, so l+s=17 and l,s>=5. The exact numerator is

    18F=4ls A²B²-(l+1)(s+1)(A²+B²)-2AB
                         -(l+1)A-(s+1)B+5.           (4)

It increases in A and B on[.8,g]. Indeed the bracket in its A derivative
is at least `(39/25)ls-18>=75.6`, while the remaining subtraction is
at most `2g+13=15.16`. The corresponding B derivative is the reversal.
At least one own contact is native>=7 or dotted unless both are5.
For native7 use the sharper bound gamma>.9005: the polynomial
`8z³-4z²-4z+1` is increasing on[sqrt3/2,1] and takes value
`-.003875599` at.9005. Therefore in this case it suffices to put
(A,B)=(.809,.9005), or its reverse. As a function of l, (4) is
concave; its minimum is at l5 or12. The smaller endpoint is exactly

    .05526706486=2763353243/50000000000>0,

and the other exceeds it by `7(.9005-.809)`. Hence R0 is impossible.
If BOTH own labels are5, the literal distance-one prohibition for
every N>=15, with both far ends retained, already excludes the pair.
Thus every opposite-helper survivor has R=-1/2.

The mixed replacement is now ACTUAL. After removing old i,j, its
two-cap Schur matrix is `I-t ones`, t=(ell-1)/(2ell), ell>=5, for a
common helper. For opposite helpers its diagonals are
delta_l,delta_s in(1/2,3/5] and off-diagonal-1/2. Both matrices are
positive definite with entrywise positive inverse. Append the entire
orthogonal old spectator block; Q0 makes this an unchanged finite
positive block. Positive acute projection is future and feasible
against every unselected actual wall.

Its spatial pressure has `Z=W_i omega_i+W_j omega_j`, W_i,W_j>0;
write W=W_i+W_j. The inherited explicit inverse calculation requires
only ell>=5, gamma>.8 and both cap time coefficients>2. The latter
holds here because, using C_ie>=5/9,

    A_cap²>(16/25)(65/9)+(4/5)(5/9)+17/36-1
          =817/180>4.

Consequently `sqrt(1+Dcal)/W<7/16` for a common helper, or<3/10
for opposite helpers. These constants retain the reversed mixed
weights, not just the outward assignment. Moreover rationalization
gives

    W_i sqrt(D_i)+W_j sqrt(D_j)-||Z||
      <W/[2(sqrt(65/9)+sqrt(20/3))]<W/10.

The actual farthest-point inequality therefore bounds the weighted
average of the two height losses by<43/80 or<2/5 respectively.
Both are<.61, contradicting two weak indices. This proves literal
helper nonadjacency BEFORE near4 is removed.

## 4. Near4 removal under the contrary energy bound

Suppose contrary e<=83/100. A left near4 in(3) has an actual mate
at j<=5, of normalized demand>sqrt3/2>.86. At the right there is
either another near4 mate or a strong k in{12,13}, by the just-proved
helper nonadjacency. A strong such demand is>.77 because q>67/25.
The normalized dual correlation is bounded by

    rho<=sqrt(5*6/(13*12))=sqrt(5/26)<.44.

Therefore the two demands sum to>1.63, whereas the squared sum-direction
norm is<2.88. Cauchy contradicts e<=.83 since
`1.63²-.83*2.88=.2665>0`. Opposite near4 mates have larger demands
and smaller correlation. Reflection treats a right near4. Every
remaining weak site is thus a nonadjacent pure helper.

## 5. One rational reference, all holes and the final contradiction

Choose a symmetric rational height sequence H_(18-i)=H_i with

    H5=2.077, H6=2.218, H7=2.314, H8=2.371, H9=2.389,

linear from0 to5 and reflected at the right. These are strictly below
the strong floors: positive squaring compares `(H_i+.61)²` with
`65/9,8,77/9,80/9,9`, respectively. The successive positive slopes
are `.4154,.141,.096,.057,.018`, then reflect negatively, so H is
concave. Its reference has coordinates `(H_i-H_(i-1))/sqrt2`.

If5 is weak, replace the entire left end by the line0 to6;6 is strong.
Make the reflected replacement if13 is weak; then12 is strong.
At every remaining internal weak j replace H_j by the average of
its neighbours. This removes its dual coefficient, preserves
concavity, and costs exactly

    (H_(j-1)-2H_j+H_(j+1))²/4.

Nonadjacency makes the affected edge pairs disjoint, so these costs
add. Every positive surviving dual is paid by an actual strong floor.
Summation by parts and completing the square therefore give
`e>=2x_A.w-||w||²>=||w||²`. No separate optimized energy is added.

The following rational bounds cover EVERY weak set; the decimals
are terminating values or abbreviations of the indicated rational sums.

| Weak end5 sites | Initial reference energy | Internal holes | Total loss at most |
|---|---:|---:|---:|
|0|.8954558|at most4|4(.00050625)=.002025|
|1|.8640827333...|at most3|3(.00050625)=.00151875|
|2|.8327096666...|at most3|3(.00038025)=.00114075|

For direct verification set a=2.077,b=2.218 and
`T=.096²+.057²+.018²=.012789`. The three initial energies are
`a²/5+(b-a)²+T`, `(b²/6+a²/5+(b-a)²+2T)/2`, and `b²/6+T`.
With zero weak ends the internal interval is6..12. With one weak end
it is7..12 or6..11; with two it is7..11. Their maximum nonadjacent
cardinalities are4,3,3. The largest absolute second differences are
.045 for the first two rows and.039 for the last, explaining every
loss bound without enumerating weak subsets.

Even the last and smallest lower bound is

    .8327096666...-.00114075=.8315689166...>83/100.

This contradicts the contrary bound and proves the A17 component
theorem with arbitrary spectators.

Finally the universal isolated-I/ordinary-T2 joint bound is>13/70:
if I has height square h<=1/5, its propagated floor gives T2 energy
>1/5; if h>1/5, I costs>3/70 and T2 costs>1/7. I is counted once.
Orthogonal addition to the SAME original A17 projection yields

    ||x||²>83/100+13/70=711/700>1.

Thus A17+I+T2 is impossible, including arbitrary further finite
spectators. This removes that one necessary general dimension20 family
after review. It does not settle A16/typed7, A18/typed7, A19/I, whole
A20/B20/D20, or give a general compact ceiling19. The existing papers,
accepted bounds and all arithmetic supporting packets remain unchanged.
