# A uniform D-component energy floor and its current limitation

20 October 2026. Paper proof only; no scientific execution, optimizer,
worker import, Lean, TeX or accepted-source edits. The floor below is
universal, with arbitrary finite spectator components. The forest
consequences concern a farthest vertex. No uniform dimension ceiling
or rank/14 rate for all D components is asserted.

## 1. Every D_r, r>=8, has covering energy greater than 4/7

Use unit old roots and the covering-energy convention of
`compact_universal_e6_e7_e8_component_energy_20261020.md`, SHA256
`acf261d97f6eb1a11a64502e17e079ef7f418b9746fbf85dbdf6915750fe9a18`.
Suppose E_Dr<=4/7. The complete internal-pair classification
`compact_dn_internal_consecutive_common_opposite_20261020.md`, SHA256
`ff1aaa687171eec9ebf2669746bdc37cd753dcbcef62890ebf49be20a791bd11`,
shows that a common opposite for directions 4/5 forces
`cosh L5=4+5/sqrt2>7`, hence component energy greater than 3/4.
Thus these two opposites are distinct.

For r>=9, the same theorem applied at i=5 makes the 5/6 opposites
distinct as well. For r=8, 5/6 is the center/long-neighbor pair, and
the complete arbitrary-spectator theorem
`compact_dn_center_neighbor_common_opposite_energy_20261020.md`, SHA256
`e371611c4903ad362654ee0d4656457140cdc5e4374b4e19a2697718c2d679b1`,
likewise rules out a common opposite under this energy assumption.

For every r>=8, the unit inverse-D entries for directions 4,5,6 are
`C_ij=2 min(i,j)`. Consequently their direction Gram is independent
of r, with correlations q=sqrt(4/5), k=sqrt(5/6), and q*k. The
distinct-opposite full-tangent inequalities and the entire analytic
argument in Sections 2--4 of
`compact_d8_three_direction_component_energy_20261020.md`, SHA256
`7555245ab9ee941efbd3a6f7657a7a103d4c9c22f6ef0a63d1330c66eb9ab860`,
therefore apply verbatim. No further D8-specific contact inventory is
used in those sections: orthogonal residual coordinates only add
nonnegative squared norm. They give the incompatible bounds
`h4*h5>=p(q)>.509` and
`h4*h5<101781/200000<.509`.

This contradiction proves **E_Dr>4/7 for every r>=8**. In particular,
the statement does not transport the numerical rate r/14 from rank 8
to larger ranks.

## 2. Farthest-vertex consequences with explicit spectator conditions

The energy sum is strictly below one at a farthest vertex, by
`compact_farthest_vertex_energy_obstruction_20261020.md`, SHA256
`50c64e5a521d3e4e0f277a8ece02d3a29aa4ae8b71485cf3088211e587f45bff`.
Thus there can be **at most one D component of rank at least eight**:
two would already contribute more than 8/7.

Call a spectator block/direction ordinary here only when an accepted
estimate gives energy strictly greater than its rank divided by 14.
A large D component plus ordinary spectators of total rank at least
six is impossible, since their sum exceeds `4/7+6/14=1`.

This condition must not silently include a short A1 or isolated
rank-two 7-defect. If a single short A1 is discarded, ordinary rank
at least six is still needed beyond that direction. If an isolated
rank-two defect is present, use its accepted energy greater than
9/100 and require ordinary rank at least five; then the total exceeds
`4/7+9/100+5/14=713/700>1`. The inherited defect accounting excludes
another defect and a short A1 in that case. These are conditional
forest exclusions, not an estimate for unclassified large A/B
spectators or an unrestricted bound on ambient dimension.

Here is the resulting complete conditional frontier. **For n>=15,
if no A/B component has rank at least five, there is one and only one
D_r of rank at least eight; it satisfies
`r>=max(9,n-6)`, and its spectator total rank is at most six.**

To check the exceptional accounting explicitly, pin the universal
small-A source SHA256
`801581c4b855289aca7982b2d0f4bb6f5849f151a78f40b86c9d403e84f415e6`
and rank-four forest source SHA256
`68ad62a619603eb7e6c2acc5689c901dde3f9d0897a7c1580b3deb948f5832b6`
(the paths are `compact_universal_small_a_component_energy_20261020.md`
and `compact_farthest_rank_four_forest_exclusion_20261020.md`). Their
statements hold with arbitrary finite spectators. The two possible
rank-two defects are OLD A2/NEW7 and OLD I2(7)/NEW3, both isolated
to their own old blocks, each with energy greater than 9/100 and
squared height below 1/2. The only short-A1 walls are I or J, supported
on A1 blocks, with squared height below 1/2; at most one such direction
exists. For defects in distinct blocks, or a defect and I/J, the
spatial supports are orthogonal. Their actual wall product is
`-sqrt(h_1*h_2)` in (-1/2,0), contradicting the native wall-product gap.
Thus at most one defect occurs and it cannot coexist with a short A1.
These exclusions require no rate for the large D block.

If there were no D of rank at least eight, every other block pays
strictly more than rank/14, apart from these named exceptions: this
uses the accepted D7 bound as well as D5,D6,E6,E7,E8 and the complete
rank-at-most-four list. Without a rank-two defect discard the possible
one short A1 and obtain energy greater than `(n-1)/14>=1`. With one
defect obtain energy greater than `9/100+(n-2)/14>1`. Hence a large D
must exist. Its uniqueness follows above. All its spectators are
ordinary apart from the stated exceptions. If their total rank were
at least seven, the three preceding spectator estimates exclude,
respectively, no exception, one short A1, or one isolated rank-two
defect. Consequently their rank is at most six and `r>=n-6>=9`.

## 3. An exact D9 relaxation witness below 9/14

The higher consecutive laws alone do not prove E_D9>9/14. In the
standard Euclidean D9 coordinates set

    x=(289/750,289/750,289/750,61/250,23/100,1/5,1/5,0,0).

For i<=7 put `u_i=(e1+...+ei)/sqrt(i)`; put
`u8=(e1+...+e8-e9)/3`, `u9=(e1+...+e8+e9)/3`.
Writing S_i for the first i coordinate sum gives

    (S3,S4,S5,S6,S7)=(289/250,7/5,163/100,183/100,203/100),
    ||x||^2=478411/750000<9/14.

The cleared difference for the last comparison is
`9*750000-14*478411=52246>0`. All h_i=x.u_i are positive and below
one. Define the mock full tangents `R_i=2h_i/(1+h_i^2)`; these are
numbers, not edge lengths of an asserted polytope.

For i=3,4,5,6, the necessary distinct-opposite condition
`R_i R_(i+1)>=sqrt(i/(i+1))` is equivalent to

    4(i+1)S_i S_(i+1) >= (i+S_i^2)(i+1+S_(i+1)^2).

Its exact positive left-minus-right differences are respectively

| i | Difference |
| --- | --- |
| 3 | 77871/1562500 |
| 4 | 1219/250000 |
| 5 | 600759/100000000 |
| 6 | 4901799/100000000 |

Even adding the center/fork and fork/fork distinct-wall laws does not
remove this witness. With z=S7^2=41209/10000, their cleared
inequalities are `36z>=(7+z)(9+z)` and
`324z>=7(9+z)^2`. Their positive differences are
`243618319/100000000` and `13006548233/100000000`.

The direction-2/3 distinct-wall law fails for this vector. In its
cleared form the left side is below 107/10, whereas the right side
exceeds `(259/100)(433/100)>107/10`. Thus the witness identifies a
missing low-index constraint; it is not a counterexample to a
complete D9 theorem. No native opposite normals, completion,
compactness or geometric realization is supplied by this relaxation.

## 4. One low-index spectator branch removed without a case table

At a farthest vertex, suppose an actual c closes direction 3, has
only a simple-3 D_r contact at node 1, and otherwise contacts
spectators. Let S be their total inverse-Gram projection cost.
Since C11=2, the full timelike height is `a^2=S-1/2>0`.
The actual replacement deleting 3 retains left A2, whose white
cost is `(1/4)(4/3)=1/3`, so `S<2/3` and `a^2<1/6`.
Some spectator node j must be contacted. Its native contact
beta_j>=1/2 and C_jj>=1 give

    (C beta)_j^2/C_jj >= beta_j^2 C_jj >=1/4>a^2.

The actual reciprocal-shadow alternative therefore says that c also
closes spectator direction j. This contradicts the farthest-vertex
orthogonal-component opposite theorem. This argument is restricted
to a farthest vertex; it does not assert that such a cross-component
common opposite is impossible at every vertex.

For this last observation the additional pinned inputs are
`compact_actual_wall_reciprocal_shadow_20261020.md`, SHA256
`c913ec3bd8e0f26a8c679a1ede2b98b251ad3f36efeae60486213ddfa3e8badc`,
and `compact_farthest_vertex_orthogonal_component_opposite_obstruction_20261020.md`,
SHA256 `4216d06c2e03fd336898b8dac289b9969fa096c8a1d08bd67ddeb8ab5ec2f23d`.
This removes one branch only. Smaller-left-chain 4/5 profiles and
other simultaneous directions require further complete accounting.

In particular the following universal native profiles must not be
discarded as arithmetic long-end exclusions: isolated 5 at old node
1; isolated 5 at old node 2; and contact 4 at node 2 plus 3 at node 3
with no spectator contact. Their required replacements for a common
2/3 opposite are positive definite. The last profile's exact
two-direction covering energy is `(6sqrt2-4)/7<9/14`; its two squared
half tangents are `(1+2sqrt2)/7` and `(2sqrt2-1)/3`, whose product
is 1/3. Both inverse-direction-Gram weights are positive, so the
two-direction interior formula gives that exact energy; its comparison
with 9/14 reduces to `12sqrt2<17`, sealed by `288<289`. Thus it needs
constraints from further directions. These are
surviving profiles, not claims of compact completions or a complete
new 2/3 inventory. The separate arithmetic first-type whole-D
long-end-5 theorem does not exclude them in the present universal,
arbitrary-spectator component setting.
