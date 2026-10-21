# Complete common-2/3 D inventory and the remaining direction-1 boundary

20 October 2026. Paper proof candidate, awaiting independent review.
No scientific execution, optimization, arithmeticity, Lean or TeX.
Accepted sources are unchanged. Other finite vertex components are
arbitrary throughout.

## 1. Statement and fixed inputs

Let p be a farthest vertex of a compact Coxeter polytope and let x_D
be its common spatial vector projected to a WHOLE D_r component,
r>=9. Number the long chain 1,...,r-2, with fork leaves r-1,r.
Use unit normals `(e_i-e_(i+1))/sqrt2` and the other fork normal
`(e_(r-1)+e_r)/sqrt2`. For long indices,

    C_ij=2 min(i,j),  u_i=(e1+...+ei)/sqrt(i).

**Theorem.** If one actual wall closes BOTH outgoing directions 2
and 3, then `||x_D||^2>9/14`. Thus these two actual opposites are
distinct whenever the actual D projection has squared norm at most
9/14. This is an actual farthest-projection statement, not a claim
about all edge-covering vectors or a uniform rank/14 component rate.

The new classification below uses the complete three-profile theorem
`compact_dn_low_index_three_profile_farthest_obstruction_20261020.md`,
SHA256 `b03a84986519cf9a9c16b1c84461dd0ed52806d9ada505a1066f8f750f345c01`.
Its first two cases use edge constraints; its 4+3 case uses a further
ACTUAL feasible negative-core point.

Other fixed inputs under `results/` are:

* `compact_actual_wall_reciprocal_shadow_20261020.md`, SHA256
  `c913ec3bd8e0f26a8c679a1ede2b98b251ad3f36efeae60486213ddfa3e8badc`;
* `compact_farthest_vertex_orthogonal_component_opposite_obstruction_20261020.md`,
  SHA256 `4216d06c2e03fd336898b8dac289b9969fa096c8a1d08bd67ddeb8ab5ec2f23d`;
* `compact_dn_internal_consecutive_common_opposite_20261020.md`, SHA256
  `ff1aaa687171eec9ebf2669746bdc37cd753dcbcef62890ebf49be20a791bd11`.

Finite spherical Coxeter classification and the exclusion of affine
principal diagrams in a compact Coxeter polytope are explicit inputs.

## 2. The complete finite replacement inventory

Write the actual common outside normal as c. Deleting old node 2
leaves right D_(r-2), rank at least seven, and isolated old node 1.
In this actual positive-definite replacement, the right D component
can have no c contact, or one simple 3 at its long leaf 3 or at one
fork leaf. To verify completeness, two contacts create a nonpositive
cycle; an internal 3 creates a second fork or a degree-four star;
and any bond at least 4 has white inverse-Gram cost at least one.
The latter uses inverse unit-D leaf entries 2 and (r-2)/2.

Deleting old node 3 leaves left A2 on nodes 1,2. Its actual finite
replacement permits at most one c contact to this A2, with label
3,4 or5. Two contacts contain the nonpositive triangle; label at
least6 contains the nonpositive rank-three path 3,6. These statements
hold with every spectator retained. They also ensure every D contact
under consideration is finite and native.

### Right contact 3 at old node 3

If c has no contact to node 2, the principal tree consisting of c,
old nodes 2,...,r and no node 1 has two simple forks with leaf pairs
(2,c) and (r-1,r). It is exact affine D, so c must contact 2.
There can then be no contact to 1. The complete possibilities are:

* 3 at2: the exact affine triangle c,2,3;
* 4 at2: the third case of the fixed three-profile theorem;
* 5 at2: white A2 cost `phi^2/3>3/4` in the replacement deleting3,
  so every spectator contact is excluded. Its residual is
  `(2-phi)/3` and its actual direction-3 length satisfies
  `cosh L3=(2+5phi)/sqrt2>7`.

Here phi=(1+sqrt5)/2. The last case has
`tanh^2(L3/2)>3/4>9/14`. This exhausts the long-leaf case.

### Right contact 3 at a fork leaf

For r>=10 the retained right D_(r-2) plus c already contains the
exact affine E8 tree with arm lengths 1,2,5; it is not positive
definite. For r=9 the retained right D7 plus c is exactly finite E8.
Consequently there is no c contact to node1 or any spectator: either
would enlarge this connected E8 inside an actual finite replacement
to a connected rank-at-least-nine finite diagram, impossible by the
finite classification.

It remains to check the possible contact to node2, which was removed
from that replacement. No contact leaves a full-star affine E8 tree
with arms1,2,5. A3 gives an exact all-3 cycle. A4 gives the exact
affine star centered at2 with branches to1,3,c, labels3,3,4. A5 gives
the retained path old1,old2,c,fork with labels3,5,3 in the replacement
deleting3. Its unit-Gram determinant is
`9/16-phi^2/4<0`, contradicting positive definiteness. Thus the fork
case is impossible, including the exceptional r=9 finite-E8 case.

### No right contact

If c has no left contact either, its D derivative is zero, so it
cannot close direction3. Otherwise the left inventory is exhausted
by one contact at1 or2:

| Contact | Complete reason |
| --- | --- |
| 3 at1, with spectators | Reciprocal-shadow contradiction below |
| 4 at1 | Exact terminal-4 affine D extension |
| 5 at1 | First case of the fixed three-profile theorem |
| 3 at2 | Exact all-3 two-fork affine D diagram |
| 4 at2 | Exact affine star with branches3,3,4 at node2 |
| 5 at2 | Second case of the fixed three-profile theorem |

For the first row let Q be the total spectator inverse-Gram cost.
The full height is `A^2=Q-1/2>0`; positive definiteness after deleting3
retains white A2 cost1/3, so `Q<2/3` and `A^2<1/6`. A spectator j
must be contacted. With beta_j>=1/2, nonnegative inverse entries,
and C_jj>=1, its reciprocal shadow satisfies

    T_j^2 >= beta_j^2 C_jj/A^2 > (1/4)/(1/6)>1.

Thus c is also the actual opposite of spectator direction j. That
direction and D direction3 belong to orthogonal old components,
contradicting the fixed farthest-vertex orthogonal-opposite theorem.
This does not assume the spectators form one simple chain.

All cases are now covered. Under `||x_D||^2<=9/14` the affine and
nonpositive retained cases are impossible, the 5+3 case exceeds3/4,
and the remaining three profiles are excluded by the fixed source.
This proves the theorem.

## 3. Exact boundary of the distinct-pair approach

This classification does not itself force projection norm above9/14.
Here is a purely rational D9 relaxation satisfying every distinct-wall
pair law on directions2,...,9:

    x=(423,423,321,240,240,240,83,25,0)/1000,
    ||x||^2=641213/1000000<9/14.

It lies in the closed D9 chamber. Its long partial sums S2,...,S7
are `(846,1167,1407,1647,1887,1970)/1000`, while S8=1995/1000
and the two fork projections both equal S8/3. Put
`F(h)=2h/(1+h^2)`.

For adjacent long directions i,i+1, the law is exactly

    4(i+1) S_i S_(i+1)
       >= (i+S_i^2)(i+1+S_(i+1)^2).

The positive differences, with denominator10^12, are:

| Pair | Cleared difference |
| --- | ---: |
| 2/3 | 1732252476 |
| 3/4 | 188938803039 |
| 4/5 | 227885305759 |
| 5/6 | 850862963679 |
| 6/7 | 57148587900 |
| 7/either fork | 251045977500 |

For the last row the equivalent expression is
`36 S7 S8-(7+S7^2)(9+S8^2)`. Multiplying successive laws along
the D9 direction tree, then cancelling intermediate factors F(h)^2<1,
gives ALL remaining distinct-pair inequalities within directions2,...,9;
the old direction correlations multiply along these paths. In
particular the fork/fork correlation is7/9.

The direction1/2 inequality instead has cleared difference
`-338772348164/10^12`, so it fails. The strict positive margins and
the strict norm margin permit an arbitrarily small perturbation into
the strict old chamber if that is desired. No actual normals,
native opposite configuration or compact completion are asserted.

The next necessary input is therefore a common-1/2 native inventory
and its ACTUAL pressure, or a genuine coupling beyond these pair
laws. This witness is not a counterexample to such stronger inputs.
In particular the actual negative-core pressure of an isolated
long-end7 can reject it. No ambient dimension ceiling is claimed.
