# Ambient58 D13 halves with three whole rank15 witnesses

For an actual quadrilateral in ambient dimension58, the simultaneous current
payment rule admits the following change. In a whole finite D13 component,
its long-endpoint/short-leaf pairs(0,11),(0,12) have weight1/2 instead of1
**only when the same actual finite corner has three distinct whole unmarked
A15 or D15 components**. Add a half to whole D28 pair(0,12), and a half
between an endpoint of whole A12 and the long endpoint of whole D16, each
only when its corner has two distinct whole unmarked A15/D15 components.
All witnesses are whole components, and their presence is not counted with
multiplicity. All other weights, including ideal weights, are unchanged.

This is a dimension-specific theorem. Its proof uses the finite-corner rank
capacity58. It does not assert the rule in larger ambient dimensions.
Combined with the frozen B13 endpoint theorem and the complete preceding
count, the exact formal finite charge maximum becomes **377**, equal to the
finite curvature. The local theorem and count are provided here for review;
the final dimension/tail composition is a separate certificate.

The [producer](../scripts/fresh_d13_three_rank15_context_transfer_20261018.py)
and [certificate](fresh_d13_three_rank15_context_transfer_20261018.json)
bind the [coupled local rule](fresh_d17_leaf_a16_context_coupled_transfer_20261018.md),
its [complete eight-obstruction count](fresh_d17_a16_coupled_ambient58_frontier_20261018.md),
and the [B13 endpoint repair](fresh_b13_endpoint_retention_20261018.md).
No global forest optimization or cusp convolution is repeated.

## The rank-capacity guard

A finite vertex in ambient58 has a positive-definite Coxeter forest of rank58.
D13 and three whole rank15 components already use all58 ranks. Thus a corner
satisfying the new negative guard has exactly those four components. In
particular it cannot also have a whole component of rank1 through4. This
excludes additional small components only after they have been proved whole;
a small induced subdiagram alone gives no contradiction.

Similarly, the positive D28 guard forces exactly D28 and two rank15
components. The positive A12/D16 guard forces exactly A12,D16 and two rank15
components. In the latter corner the two endpoint coupons contribute one
in total to its vertex charge, although the local face uses only its one
marked endpoint. The two other rank15 factors are presence witnesses, not
a factor multiplying the coupon.

## Complete selected finite support and arbitrary extensions

The fresh selected-domain reconstruction has801 eligible supporting patterns,
862 maps to finite D13,33 affected contexts and76 maps of the changed pair.
Every physical core uses norm2 roots. The producer rechecks its Lorentz
signature, strictly positive white and all four white-plus-one-black edge
subdiagrams, and the actual selected-corner classification. The earlier
supporting-census completeness and native-label classification remain
inherited inputs; the following arguments address arbitrary extra whites.

Use the independent current short8-to-short7-to-short6 route whenever its
adjacent-black path occurs. Its current finite6/7/8 and protected short
outputs are retained, as detailed below. Otherwise the complete selected
contexts split into four groups.

Five proper contexts have a displayed A12 in the selected whole D13. Each
map requires a genuinely new missing fork leaf. It has a forced unit core
contact whose inverse projection square is one of6,12,30,102, hence exceeds
its norm2. Existing old corner roots outside the A12 are orthogonal to it,
so none can supply this missing leaf. Same-end inverse positivity permits
arbitrary additional contacts without lowering the bound; contacts at both
ends give the independent short escape. These five completions are impossible.

Four full-core radical contexts are1985,2378,2575,2581. Once the selected
D13 is whole, every remaining possible contact position has a nonzero
radical coefficient. An additional white must pair to zero with the
radical. A nonempty native contact set therefore meets both signs, which
gives a checked adjacent-black shortcut of length at most8. Without that
escape the entire displayed core is sealed. Its unmarked A4 is then whole,
contradicting the rank-capacity guard.

In twelve further contexts at least one displayed unmarked A1/A2 block has
no possible new contact of projection square at most2. Any extension of
that block would have a first new root contacting it, so the block is
whole, regardless of possible chains farther outside the core. Its small
rank again contradicts the guard.

The remaining twelve contexts require a bound on an arbitrary attached
forest, rather than just on its first neighbor. The displayed small factor
K has rank s=1,2,3. Its unique possible contact position p has inverse
square r=4/3,3/2,2 respectively. All remaining inverse diagonals are positive;
a contact magnitude larger than1 is at least sqrt2 and is excluded. The
saved positive-end bounds also show that two same-end unit contacts have
projection square greater than2. Thus every new root contacting K has
exactly one unit core contact, at p, unless the short escape occurs.

Four of these contexts permit another contact position q at the other end.
Here J_pp=J_qq=4/3 and J_pq=−4/3. If distinct new roots contact p and q,
their positive-definite residual squares are2/3. Cauchy gives their full
inner product at most−4/3+2/3=−2/3<0. They therefore have a native bond,
and the checked core distances produce a short path of length at most8.
Consequently, if K grows at all and there is no short escape, no new root
contacts the opposite position. Roots orthogonal to the core may still
occur, and are included in the next argument.

Let Y be every new root in K's whole finite component, after removing K,
and let M be its total rank. Y is a positive-definite forest and every one
of its connected components is attached to K. All contacts of Y with the
old Lorentz core occur at p. Write u for their nonnegative unit-or-zero
contact vector and G_Y for their Gram matrix. The one-negative-direction
condition on the enlarged core gives the positive-semidefinite Schur
complement

    G_Y − r u u^T >= 0,

and therefore u^T G_Y^(-1)u <=1/r. This includes indirect chains and all
branches in Y, not merely the roots with nonzero u.

For a connected finite Coxeter Gram matrix of rank m, diagonal2 and
nonzero native magnitudes at least1, every inverse diagonal is at least
m/(m+1). Here is an all-rank proof of the bound used. Root the connected
graph at the chosen node and delete it. Each remaining connected component
of rank m_j has a nonzero contact with the root. Inverse matrices are
entrywise nonnegative. Induction and Schur complementation give a root
projection at least sum_j m_j/(m_j+1), which is at least(m−1)/m; hence
the inverse diagonal is at least1/(2−(m−1)/m)=m/(m+1). The base m=1 is
1/2. Multiple root contacts or larger labels only increase the projection.
For several components of Y, every component has a unit contact, so

    u^T G_Y^(-1)u >= sum_j m_j/(m_j+1) >= M/(M+1).

Thus r M/(M+1)<=1, or M<=1/(r−1). In the three cases M<=3,2,1 and
s+M<=4. The whole component containing K is therefore small. It cannot
be one of the three required rank15 components. This proves the selected
transfer for arbitrary white support, not just for the33 saved cores.

## The genuine affine-D12 output requiring repair

Selected support is not a universal output inventory. In particular, the
[opposite-leaf affine-D12 theorem](fresh_affine_d12_leaf_pair_transfer_20261017.md)
really uses the affected whole finite D13 unit in its quiet-end branch.
We repair that theorem with its complete all-label profiles and exact
physical determinant identity, rather than deducing validity from the
absence of a D13 hit in later saved finite corners.

Write its affine factor F of order13 with spine0 through8 and short leaves
9,10 at0 and11,12 at8. The chosen ideal marks are a=9,b=11, with reflected
choices equivalent. The opposite black roots are c,d. The spectator
forest Z is orthogonal to F; let x,y be the nonnegative contact vectors
of c,d with Z. An adjacent-black bond or a common spectator neighbor is
an independent short path, already paid. In its absence x,y have disjoint
supports. This disjointness is essential whenever a difference energy is
used below. The old theorem's complete profile reduction gives each end
cost at most8, and every nonextreme end cost at most4.

If an end is nonextreme, the determinant forces costs8+4 and x=y=0.
If the quiet extreme root attaches to the selected old leaf, the old
whole A1/D12-short-leaf unit is unchanged. If it attaches to the mate,
the old quiet corner is a whole D13 of the affected kind. If its guard
fails its unit is unchanged. If its guard holds, Z consists exactly of
three whole A15/D15 factors, all orthogonal to the complete displayed star.
There are three nonextreme equality profiles, each repaired at a distinct
actual corner:

* For `isolated4`, the other corner has whole A1 and whole D12 at its short
  leaf, and receives the retained unit.
* For `mate4`, the other corner is a whole affine B13 and is handled first
  by the independent current B13 rule. Its endpoint repair pays a whole
  A1/D12-short-leaf unit and does not use either new D13 half or coupon.
* For `center3`, the other corner is another whole D13 at the affected
  pair. The same three orthogonal rank15 factors give its guard. Two
  distinct halves, one at each actual corner, sum to1.

The producer reconstructs all24 affected quiet configurations: eight of
each kind, including old mark choices and orientation. It identifies the
physical whole components and gauge-aware affine B13 images.

## Arbitrary one-sided spectators in the energy-one boundary

When both ends are extreme their heights are exactly1/2, and the old
Schur identity says

    (x−y)^T G_Z^(-1)(x−y)=1.

For two nonempty supports the inherited complete energy-one theorem
already gives either a short path or the sealed small-factor alternatives;
those current payments are unchanged. The D13 unit at issue is used in
the arbitrary one-sided case, say y=0. The quiet d corner excludes c,
so it retains all of Z orthogonally. If d attaches to the selected leaf,
its A1/D12 unit remains. If it attaches to the mate and the new D13 guard
holds, then Z is exactly three whole A15/D15 components. This is the only
new spectator classification required.

The actual edge diagrams imply Z+c and Z+d are strictly positive definite.
For A15, the only inverse diagonals below2 occur at nodes0,1,13,14 and
are15/16,7/4,7/4,15/16. Two nonzero native contacts already have energy
at least2, by inverse positivity, and are forbidden by this strict edge
condition. The complete allowed one-root profiles are simple contacts at
those four nodes and doubled contacts at0 or14. The next native squared
magnitude(3+sqrt5)/2 is already too large. For D15, the only inverse
diagonal below2 is1 at its long endpoint0, and only a simple contact fits.
These bounds cover arbitrary labels, not a bounded label search.

The producer computes the exact inverse matrices and all allowed profiles.
Every nonzero component energy is at least15/16; no A15 profile has energy1.
For completeness it also checks all disjoint two-sided profiles in a single
A15/D15 factor: every nonzero difference energy is at least15/16 and is
not1. Equal-node contacts cannot occur because x,y have disjoint supports.
Two active components already give energy at least15/8>1. Thus energy1
on these three factors is exactly a unit contact at the long endpoint of
one whole D15 by one black root, with the other two rank15 factors inactive.
Those two factors are orthogonal to F,c,d and the active factor, and hence
remain whole in both paying corners; there is no indirect white bridge.

Assume c is active and d is the quiet mate extension. Two physical30-root
stars cover the remaining choices, with all reflected choices inherited:

* If c attaches to selected a, that old leaf is absent at the opposite
  corner. There c extends the active D15 to whole D16 at its long endpoint,
  while d extends the old path to whole A12 at an endpoint. The new guarded
  half at this distinct corner supplements the quiet D13 half.
* If c attaches to a's mate, the opposite corner is a whole D28 with marks
  (0,12). Its new guarded half again supplements the quiet D13 half.

Both guards use the two inactive whole rank15 components. The exact models
have signature(28,1,1), positive white rank26, positive edge rank27, and
positive-semidefinite corners. Their canonical marks and whole components
are reconstructed directly. Adding the two orthogonal rank15 factors gives
the actual ambient58 configuration. The proof is independent of spectator
coordinates or any unrecorded native label.

## Retention of the complete current chosen outputs

Apply the repaired affine-D12 route whenever a previously chosen output
uses its quiet D13 unit. In other routes retain the same old positive
actual contributions. This is a chosen-output argument, including later
replacement outputs; merely placing the changed D13 corner first would
not prove it. The following universal interfaces are inherited from their
complete routing proofs, with the stated marked-role checks.

The [low/protected/adaptive interface](fresh_d12_low_protected_adaptive_interfaces_20261018.md)
has66 low contexts of order bound at most12,74 protected paths and9 adaptive
paths. Every protected/adaptive marked path has length9; the changed D13
pair has length11. The universal finite8 cycle/tree and cycle12–14 outputs
are short paths or marked distance9 at D13; their second-node D_m outputs
of distance11 have m at least18. Current depth6 finite6/7/8 transfers escape
to independently paid short paths, A-type or protected whole-leaf outputs,
not this D13 long-end/short-leaf pair. Their triangle protections also remain.

The old affine-D11 one-sided output is D12. Its relevant two-sided sealed
A1+A3 boundary gives D13 marks(1,11)/(1,12), which remain unit; these are
the eight replacements already used by the current A16 theorem. The
current ideal13 center/leaf theorem gives finite D13 marks(0,10), at its
fork center; this role is reconstructed directly here. In the coupled
order11 unmixed boundary the finite D_(11+r), r=1,2,3, has marks
(r,n−2)/(r,n−1), of distance n−2−r=9. At r=2 this is D13(2,11)/(2,12),
again disjoint. Other current B/C and order14–17 outputs have the previously
identified B/A types, other D orders, or cross-component roles.

The source freshly checks physical marked roles in all216 inherited literal
configurations,564 later finite corners, the newer coupled/A13/A14/A16/D17
sealed families, all24 current D14 variants, and all332 positive chosen D12
repair outputs. None uses the changed D13 pair. The old affine-D12 quiet
and arbitrary one-sided branches above are explicitly outside that
finite-domain absence argument. The
[current D12 composition](fresh_d12_half_face_composition_20261018.md)
and the bound local packets provide the complete current route inventory;
the present source does not turn that written inherited universal proof
into a machine-certified census by checking these finite records.

Every negative use is therefore either impossible under the rank58 guard,
replaced by retained independent contributions, or paid by the two new
halves at distinct corners. New positive contributions introduce no further
negative obligations. The two new coupons are used only in the two
energy-one stars proved above. Their rank15 witness components are not
paid marks and are not altered by this rule. The independent B13 priority
uses retained D12 cross units, not a D13-dependent or newly financed route,
so the repair is noncircular. Ideal weights are unchanged. Triangle
retention uses the actual-face
[finite7/ideal5 triangle theorem](fresh_finite7_ideal5_triangle_payments_20261017.md),
not a sum of three pairs in one finite diagram. Its complete minimal-core
argument gives two distinct actual-corner payments: low cores give paths
of length at most5; all1,025 high-core triples give elementary paths of
length at most7 or one of four protected E8 alternatives. Those alternatives
have finite marked length6 or7, or a whole exceptional affine factor of
order at most9. The
[current finite6 refinement](fresh_depth6_finite6_zero_transfer_20261018.md)
retains every elementary length6 choice and all three exceptional length6
choices by their persistent nearby feature or whole leaf; its fourth
exceptional alternate has length7. The later finite7/8 changes retain
these same protected outputs. The present negative pair has distance11,
so none of those chosen finite contributions changes. Their ideal
contributions are unchanged as well, and new positive coupons cannot
decrease the triangle total2.

## Complete finite count without another optimization

The frozen coupled certificate proves that every actual rank58 forest
except its eight listed B13/D13-plus-three-rank15 forests has charge at
most377. The B13 endpoint theorem lowers its four forests by1. The present
negative D13 guard lowers the other four by1. All eight now attain377.
This uses the complete actual-type certificate, not just saved witnesses.

A new positive coupon can occur only in one of the exact rank58 families
identified at the start. Their complete A15/D15 type choices give six
forests. No A1 or other contextual witness can also fit. Their costs are

    D28 + 2 rank15:       170 + 200 + 1/2 = 370.5,
    A12 + D16 + 2 rank15: 66 + 108.5 + 200 + 1 = 375.5.

The old A12/D16 cross cost is zero, checked directly. Both values are below
377. Outside these two families the new rule is pointwise no larger than
the frozen coupled rule, including its B13 deletion. This proves the exact
formal maximum377 for all finite forests; no representative upper is
mistaken for an actual attainer and no possible family is discarded because
it was absent from an earlier list of maxima.

The ideal rule, its complete2552 low8 terminal profiles and positive-only
extension to actual cusp surplus36 are unchanged. The strict global
correction4205 is inherited once. Consequently the new local/count result
supplies the formal contradiction

    0 >= 4205 + 36 v_i > 0

when inserted into the ambient58 curvature argument. A separate final
binder must record the common hypotheses, completed58 higher-dimensional
tail, and verification scope before an unrestricted57 theorem is reported.
There is no claim here of external refereeing or end-to-end Lean verification.

## Verification status

Authoritative generation27672 passed all selected, spectator, interface and
finite-cost checks with exact maximum377. Its terminal log was read after
completion. Separate-process saved replay55090/165a81 terminated exit0/PASS.
An earlier added marked-role assertion had an off-by-one path-distance
typo, corrected from n−2−r+1 to n−2−r; that run stopped before writing a
certificate and changed no mathematical role or scope. The final source
and data include the corrected role check.

Root completed a full source/report mathematical read, with the original
affine-D12 theorem reread: the selected arbitrary-support proof, quiet and
one-sided energy classification, actual-corner coupons and rank guards,
complete finite closure and inherited current-priority inventory passed.
The review did not execute the producer or reprove the historical census.
Its triangle-prose correction replaces a same-diagram shortcut by the
explicit inherited actual-face triangle proof described above.

Construction completed an independent full mathematical source/report read,
also rereading the complete affine-D12 source/report, affine-B13 singleton
and fork reports, immediate radical/end-inverse constraints, and current
center/leaf and coupled-order11 output roles. The33 selected contexts,
arbitrary attached-forest rank bound, rank58 whole-component guard,
noncircular quiet repairs, rank15 energy-one classification and the two
complete increasing families passed. No producer, native-map, inertia,
literal, DP or arithmetic replay was performed by this reviewer; the
historical census and stated universal routing remain inherited.

Both reviewers subsequently read and passed the corrected actual-face
triangle paragraph; construction additionally reread the complete
finite7/ideal5 triangle report and current finite6 protection section.
No triple enumeration was replayed. Final integrity verification passes
484 source bindings,10 local links and one Python syntax tree. Source,
certificate and report are final/frozen. No external refereeing or
end-to-end proof-assistant verification is claimed.

    python3 scripts/fresh_d13_three_rank15_context_transfer_20261018.py
    python3 scripts/fresh_d13_three_rank15_context_transfer_20261018.py --verify
