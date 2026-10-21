# Finite rank16 and17 savings with one D14 half-coupon

The same-component finite costs can simultaneously be reduced to
**A16=116, B16=104, D16=118 and D17=125**, preserving every triangle's
payment of at least2 and every quadrilateral's payment of at least1.
The only added payment is a half-coupon between any marked node of a
whole A1 or A2 component and the unique long-arm neighbor of the fork
center of a separate whole D14 component. Its added cost is1/2 per
A1/D14 component pair and1 per A2/D14 pair. Every ideal weight and all
other finite payments remain unchanged.

The [source](../scripts/fresh_finite_d16_d17_long_pair_half_transfer_20261017.py)
and [certificate](fresh_finite_d16_d17_long_pair_half_transfer_20261017.json)
prove this local transfer. A finite-forest maximum and any dimension
consequence require a separate budget. In particular no unrestricted
ceiling64 is asserted here.

## Canonical simultaneous rule

Number finite D_n by the spine0--1--...--(n-3), with its two short
leaves n-2,n-1 attached to n-3. Change these unit weights to halves:

|Whole component|Changed marked pairs|Old total|New total|
|---|---|---:|---:|
|D16|(1,12),(1,13),(1,14),(1,15)|120|118|
|D17|(2,13),(2,14),(2,15),(2,16)|127|125|

Each set is invariant under the full labelled-diagram automorphism
group. The distances are11,12,13,13; hence no short8 or protected
whole-leaf9 payment changes. All nine old finite D17 deletions remain
deleted.

Number A16 as path0--...--15, and B16 by that path with the bond14--15
doubled. Delete the following unit payments:

* A16: (0,14),(0,15),(1,14),(1,15).
* B16: (0,11),(0,12),(0,13),(0,14),(0,15),
  (1,12),(1,13),(1,14),(1,15), (2,13),(2,14),(2,15),
  (3,14),(3,15),(4,15),(5,15).

These sets are also invariant and disjoint from the short8 and
whole-leaf9 sets. Finite C16 has the same labelled graph as B16.
The old rank17 A/B costs114/105 remain unchanged.

In D14 the new coupon's distinguished node is10 in the same numbering:
the fork center is11 and its short leaves are12,13. Node10 is uniquely
specified intrinsically. The old coupon already assigns halves at
11,12,13 to every mark in a separate whole component of rank at most5.
The present addition concerns node10 only and small components A1/A2
only; it neither duplicates the old coupon nor enlarges its rank cutoff.

## Complete selected support

Use the completed finite8/ideal5 face reduction and its current
[compatible successors](fresh_ambient65_ideal_margin_closure_20261017.md).
In the complete2,614-context high-overlap domain, a whole finite target
of rank16 needs only394 patterns:2,045 contexts have whole marked
component bound below16,92 already have a selected base larger than16,
and74 protected plus9 adaptive contexts have their separate interfaces.
A rank17 target needs140 patterns, with2,391 smaller-bound contexts
and the same83 protected/adaptive contexts handled separately.

Every remaining pattern's connected selected component is mapped
injectively and inducedly into each target, retaining every native
label and nonedge. The native label is determined from
4*g_ij^2/(g_ii*g_jj), so a positive diagonal gauge never changes an angle.
The backtracking embedding routine tries every unused target node
consistent with the assigned labels and nonedges; its degree filter is
only necessary. It thus includes every actual whole-component map,
even if some emitted maps need not be geometrically realizable.

|Target|Eligible patterns|Patterns with maps|Complete maps|
|---|---:|---:|---:|
|A16|394|70|238|
|B16|394|28|100|
|D16|394|240|578|
|D17|140|75|168|

None of the selected A16/B16 images is deleted. In contrast, every
proposed D half has actual selected support. Exactly eight contexts
and16 maps per D target are affected:

    D16: 1146,1161,1174,1185,1268,1909,2302,2491
    D17: 1017,1034,1048,1059,1264,1904,2297,2488

All16 selected bases are already full finite D16 or D17 components.
There is no proper-completion case to discard. The proof therefore
uses alternate actual corners, rather than absence from a map census.

## All-label radical sealing and the alternate halves

Every affected case has the same19-root norm-two core, constructed from
two copies of catalogue diagram25 identified at their node9. Its exact
Gram has signature(17,1,1), with radical coefficients

    (-6,-3,-4,-2,-5,-4,-3,-2,-1,0,
       6, 3, 4, 2, 5, 4, 3, 2, 1).

The checker reconstructs the white15-root positive-definite principal,
all four strict edge principals, and all four finite corner principals.
It recomputes the selected corner's complete necessary contact domain
from its actual finite blocks and compares every field with the frozen
high-overlap certificate.

In an actual Lorentz realization a radical combination of this core
represents the zero vector: it is null and orthogonal to a span containing
a timelike vector, whose orthogonal complement is positive definite.
Consequently every additional supporting white's contact vector must
annihilate the displayed radical. Every allowed contact position has a
nonzero coefficient. Contacts of only one sign cannot annihilate it,
because acute contact magnitudes are nonnegative. For every opposite-sign
pair of allowed positions, the checker exhibits a marked path of length
at most8 through that white. This uses no upper-label cutoff and allows
arbitrarily many simultaneous contacts. Shortest paths can be truncated
at consecutive black vertices to obtain the actual adjacent-corner path.

Unless the already-retained short8 branch pays, every additional white
is therefore orthogonal to the entire19-root core. This excludes all
direct attachments, indirect tails and bridges to its displayed
components. Every displayed finite component is a whole actual vertex
component. Arbitrary extra supporting whites form separate positive-
definite blocks, so the alternate corners remain finite; no simple-cusp
assumption has entered.

Each D16 context has a separate whole A2/D14 corner, and each D17
context a separate whole A1/D14 corner. The D14 mark is precisely the
counterpart of the changed large mark:

|Large mark in D16/D17|Alternate D14 mark|Payment|
|---|---:|---|
|12 /13|10, long-arm neighbor of fork|New half|
|13 /14|11, fork center|Existing half|
|14 or15 /15 or16|12 or13, short fork leaf|Existing half|

The new role is needed in only two contexts of each rank, namely
1185,2491 and1059,2488. Each configuration receives its retained
selected D half plus the alternate half, totaling1. No other finite
same-component unit is spent in this argument. All component types,
canonical marked images and coupon roles are checked from the full
corner Gram.

## Simultaneous nonselected and universal interfaces

Selected support alone would not justify these changes. In particular,
the A16 endpoint pair is used by an older literal ideal replacement.
The checker reconstructs216 complete literal configurations and debits
all four new finite rules simultaneously, after first applying the
completed A17/B17/D17 deletion:

|Source|Configurations|
|---|---:|
|Old affine-order16 and order17 replacements|9+9|
|Universal affine D13 and D14 equalities|72+72|
|Order17 partial-weight replacements|12|
|Order16 partial-weight replacements|10|
|Universal affine D11 A1+A3 spectator equalities|32|

Each finite shared component is reconstructed and mapped to its whole
finite target. Every resulting literal total is at least1. Exactly25
configurations lose payment:

* Old order16 case663 loses its A16 endpoint unit but keeps two D16
  units, at pairs(0,14)/(0,15), which are unchanged. Its total stays2.
* Twenty-four D14 equality configurations lose one D16 half. Eight of
  these also lose an A16 unit. The remaining totals are at least1.

The old A17 endpoint case662 still keeps its two D17 units. No new D17
half affects any literal unit in these216 configurations. The D11
equalities' A16 images are(1,12)/(3,14), both retained. The optional D16
units displayed in the later partial-weight proof are unchanged and
are not needed by that proof, which uses an ideal half and a finite
cross-component half.

The other branches have explicit retained interfaces:

* Finite triangles use distance at most8. The short-path quadrilateral
  branch is unchanged. Its universal ideal-cycle/tree replacements
  use finite paths of length at most8, actual whole-leaf9 pairs, or
  separate-component coupons. The small P6 equalities are whole A11
  or A7; the cycle18 equality uses separate A1/A17 coupons. The
  cycle16 one-root finite alternatives have distances7 or8.
* The complete low-core closure has bound12. The all-label high-disjoint
  domain is excluded geometrically. The74 protected cases retain their
  actual leaf9 or basic ideal payment.
* The nine adaptive cases retain their leaf9 payments. In the two
  bound16 branches the selected base is a tree with marked distance9.
  Its complete D16 induced maps avoid every changed pair; A16/B16
  completions also retain that distance9 and no deleted pair has that
  distance. A larger actual target is excluded by the bound16 proof.
* Universal D12 nonextreme alternatives are whole D13; its sealed
  equality alternatives are at most15. Universal D11 one-sided
  alternatives are whole D12 or a separate A1/D11 coupon. Its32
  remaining equalities were included above. Their finite D alternatives
  are D13/D15, with the only rank16 path being A16.
* The B15 older finite alternatives have whole order at most14. The
  later B13 and B12 alternatives have order at most13 or separate
  coupons. The newer B/C missing-terminal proofs introduce no finite
  alternative. Cycle12/13/14 outputs have short finite paths or
  separate coupons. Immediate D13/D14 ideal outputs are whole affine
  factors of orders12/13, whose own completed refinements remain.
* The new affine D15 and D17 transfers explicitly consume no finite
  D16/D17 same-component units. Their existing A4/A5--D11 halves and
  E6--D11 quarter remain unchanged. The large-finite-D half transfer
  concerns ranks at least20 and uses its separate small/D10 coupon.

The frozen all-label geometric census, universal supporting-root
classification and residual sealing arguments remain inherited inputs.
This package re-enumerates the affected finite induced maps and checks
the new radical argument and all simultaneous literal totals. It does
not equate216 literal configurations with all quadrilaterals.

Thus the new rules preserve complete face coverage. They leave all
current ideal envelopes intact, including the proved positive
dimension65 ideal margins. Counting their finite cost is separate:
large finite components can still dominate the unrestricted maximum.
The existing short8 interface is retained explicitly; this report does
not assume the still-separate finite7 quadrilateral route.

## Verification status

Final generation79090 and separate fresh saved replay45872 both
terminated exit0/PASS on the complete finite-target maps,16 radical
cores,32 changed D maps,216 simultaneous literal totals and the stated
later interfaces. Python AST, all125 direct bindings and all3 local
report links pass. The final source explicitly checks that the two
newer literal domains have no whole rank17 finite component, before
using their direct old weights without an earlier rank17 debit.

Root completed a separate full source/report mathematical read, including
the imported native-label and complete induced-embedding helpers. The
review passed canonical role invariance, all-label radical sealing,
the alternate whole A1/A2--D14 coupon, A16/B16 unsupported-pair deletion,
and all216 simultaneous literal debits and universal interfaces. A
separate saved query confirmed the absence of a whole rank17 unit in
the order16 partial-weight package. This review did not independently
re-enumerate the selected maps or repeat the upstream geometric census.
No second independent mathematical read or independent arithmetic
reconstruction is claimed. Source, data and report are final/frozen.
No external refereeing or proof-assistant verification is claimed.

    python3 scripts/fresh_finite_d16_d17_long_pair_half_transfer_20261017.py
    python3 scripts/fresh_finite_d16_d17_long_pair_half_transfer_20261017.py --verify
