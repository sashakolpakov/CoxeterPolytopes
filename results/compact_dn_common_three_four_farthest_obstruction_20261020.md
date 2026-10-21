# A common opposite for D directions 3 and 4 forces projection above 9/14

20 October 2026. Complete paper candidate, awaiting independent review.
No scientific execution, optimization, arithmeticity, Lean or TeX.
Previously accepted sources are unchanged.

Let p be a farthest vertex of a compact Coxeter polytope, with a WHOLE
D_r component, r>=9, and arbitrary other finite components. Number
the long chain1,...,r-2 and fork leaves r-1,r. With unit old normals,
write `C_ij=2 min(i,j)` for long indices and let x_D be the actual
farthest spatial projection to this D component. Put phi=(1+sqrt5)/2.

**Theorem.** If one actual wall c closes both directions3 and4,
then `||x_D||^2>9/14`. Consequently these actual opposites are
distinct under the contrary projection bound. This is an actual
farthest-point result, not an edge-only energy theorem or a new
ambient dimension ceiling.

## Fixed inputs

Use the actual reciprocal-shadow theorem
`compact_actual_wall_reciprocal_shadow_20261020.md`, SHA256
`c913ec3bd8e0f26a8c679a1ede2b98b251ad3f36efeae60486213ddfa3e8badc`,
and the farthest orthogonal-component opposite obstruction
`compact_farthest_vertex_orthogonal_component_opposite_obstruction_20261020.md`,
SHA256 `4216d06c2e03fd336898b8dac289b9969fa096c8a1d08bd67ddeb8ab5ec2f23d`.
The low-index three-profile theorem is
`compact_dn_low_index_three_profile_farthest_obstruction_20261020.md`,
SHA256 `b03a84986519cf9a9c16b1c84461dd0ed52806d9ada505a1066f8f750f345c01`.
All paths are under `results/`. Finite spherical classification and
the exclusion of affine principal diagrams in a compact Coxeter
polytope are explicit classical inputs.

## Complete finite replacement inventory

Deleting old3 leaves right D_(r-3), rank at least six. A finite
one-node extension has no c contact, or one simple3 at its long leaf4
or at a fork leaf. Multiple contacts create a nonpositive cycle;
internal3 gives a second fork or a degree-four star; a bond at least4
has white inverse-Gram cost at least one. The right inverse leaf
entries are2 and (r-3)/2.

The small fork exceptions are retained: D6+c can be E7 when r9,
and D7+c can be E8 when r10. For r>=11 that fork extension already
contains the affine E8 tree. We do not discard possible spectator
contacts merely because an E7 replacement occurs.

Deleting old4 instead retains left A3 on1,2,3. Its complete finite
one-node inventory allows a single3 at any node, or a single4 or5
at an endpoint1 or3. An internal4 has the nonpositive star3,3,4;
a label at least6 has a nonpositive3,6 subpath; multiple contacts
give a nonpositive cycle. These are actual finite replacements with
all spectators retained, so their positivity is mandatory.

## Right simple3 at old4

If c does not contact3, the principal tree on c and old3,...,r is
exact affine D, with leaf pairs(3,c) and(r-1,r). Thus c contacts3,
and it cannot contact another left A3 node. Its label at3 is3,4 or5.

A3 gives the exact affine triangle c,3,4. For label4 or5 let Q be
the spectator inverse-Gram cost. Every nonzero component cost is at
least1/4. The left A3 endpoint inverse entry is3/2, giving white
cost3/4 for label4, or3phi^2/8 for label5. Positive definiteness
therefore forces Q=0 in BOTH cases.

For label4 the actual direction-4 length is

    cosh L4=3+2sqrt2,
    tanh^2(L4/2)=1/sqrt2>9/14.

The last comparison squares to196>162. For label5 the residual after
deleting4 is `(8-3phi^2)/8`, and

    cosh L4=(3phi+4)/sqrt(8-3phi^2)>7.

Here `0<8-3phi^2<1`; for example `8/5<phi<5/3` proves these two
strict inequalities. The same conclusion also follows from the
smaller isolated-5 formula checked below. Thus the long-leaf case
contradicts `||x_D||^2<=9/14`.

## Right simple3 at a fork

For r>=11 this was already excluded inside the actual finite right
replacement. The following argument also handles r9,10 and every
spectator possibility in those exceptional E7/E8 replacements.

A left3 contact gives an exact all-3 cycle through c, the contacted
fork, the old fork center, and the old path back to that left node.
A left endpoint4 gives the retained path fork,c,j,j',j'' with
labels3,4,3,3, the exact affine F4 path; here j,j',j'' is1,2,3 or
3,2,1. It lies in the replacement deleting4. A left endpoint5 gives
the retained3,5,3 path, whose determinant is
`9/16-phi^2/4<0`. No internal left4 or5 belongs to the finite left
inventory. With no left contact, the full old D_r+c contains the
exact affine E8 tree of arm lengths1,2,5. The five long-arm nodes
exist for every r>=9. Spectator contacts do not change any of these
principal diagrams. The entire fork case is therefore impossible.

## No right contact

With no left contact either, c has zero D derivative and cannot close
direction4. Otherwise the complete left cases are as follows.

* Simple3 at2 or3 gives an exact two-fork affine D tree, using the
  preceding old node and c as the first pair of leaves.
* A4 at1 is the terminal-4 affine D extension; a4 at3 gives the
  exact affine star at3 with branch labels3,3,4.
* A3 at1 requires spectator contacts. If their full inverse-Gram
  cost is Q, the height square is `A^2=Q-1/2`. In the finite
  replacement deleting4, the left white cost is3/8, so Q<5/8 and
  A^2<1/8. For any contacted spectator j, nonnegative inverse entries
  give `T_j^2>=beta_j^2 C_jj/A^2>(1/4)/(1/8)>1`. Reciprocal shadow
  makes c also the actual opposite of j, contrary to the farthest
  orthogonal-component theorem. This handles arbitrary spectator
  forests, not just a spectator chain.
* A5 at1 has white cost3phi^2/8>3/4, so Q=0. Its direction-2
  shadow has `T2^2=phi^3/2>1`; hence c also closes2. The fixed
  isolated-5-at1 theorem for common2/3 applies and forces projection
  above9/14.
* A5 at3 again has Q=0. Its actual direction-4 length satisfies
  `cosh L4=3phi/sqrt(8-3phi^2)>7`.

For the last strict comparison square and simplify:
`9phi^2>49(8-3phi^2)` is equivalent to `39phi^2>98`, and
`phi^2=phi+1>13/5` proves it. Positivity of the denominator follows
from phi<5/3. Thus the length half-tangent square exceeds3/4 and
cannot lie under9/14. The larger 5+3 formula in the long-leaf case
is consequently greater than7 as well.

All finite native cases, the small right E7/E8 exceptions and every
spectator branch have been addressed. This proves the theorem. No
claim about common directions1/2 or an unconditional D_r rate follows
without their remaining native7 analysis.
