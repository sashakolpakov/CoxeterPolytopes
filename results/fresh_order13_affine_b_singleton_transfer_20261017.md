# Removing the doubled-terminal / fork-center pair from affine B12

**Completed local argument; verification is recorded below.** An affine B12 factor has13 nodes. In the canonical numbering used by the [order13 support diagnostic](fresh_order13_payment_support_20261017.md), its fork center is0, its short leaves are11,12, and the spine is0–1–…–10 with doubled bond9–10. This package removes only the singleton marked pair `(0,10)`, of graph distance10, from the otherwise all-pairs rule. The factor cost changes from78 to77.

The replacement uses one additional finite half-payment: the two marks lie in a whole A2 component and the unique doubled terminal of a whole B10 component. Either A2 endpoint is eligible. The existing finite half-payment for whole A1 versus the fork center of whole D12 is retained. A separate finite-vertex budget must account for the new rule. This report asserts no ambient dimension bound.

There are two distinct reasons the removed pair was previously required. Both are discharged here: the two selected L10/L25 face-star contexts, and the universal `inner4` alternative inside the affine-D14 transfer. Discharging only the selected contexts would be insufficient. The [source](../scripts/fresh_order13_affine_b_singleton_transfer_20261017.py) and [certificate](fresh_order13_affine_b_singleton_transfer_20261017.json) record both interfaces.

## Complete support and retained short pairs

The frozen order13 diagnostic enumerated all801 selected patterns that can require a factor of order13, including patterns with larger component bounds, and their induced labelled embeddings. The removed singleton occurs in exactly contexts160 and1206 of the complete finite8 overlap certificate. These are the two orientations of L10 with marks `{3,4}` and L25 with marks `{0,3}`, identifying L10 node5 with L25 node9. The source checks the saved marked-image support, rather than inferring completeness from the two examples.

Every affine-B12 diagram automorphism fixes the fork center and the doubled terminal, so the singleton is invariant. Its distance10 exceeds every short8 payment and every elementary triangle pair in the retained ideal5 triangle rule. The new source checks the complete basic and triangle pair inclusions. Thus the short8-first reduction and triangle theorem are unaffected. The original complete overlap, no-extra-cross-edge and all-support classification remains an explicit frozen input; it is not rerun here.

The immediate D14 interface is also checked literally. Of its six saved profile/mark records, only `inner4` with the opposite fork-center mark maps to the singleton. In fact the argument below excludes all three `inner4` far-mark alternatives under the unpaid-face hypothesis, including the two still-retained distance11 B pairs. The `inner_next3` alternatives have affine D type and are not changed by this package. Whole factors of order at most12 and every other order13 pair remain paid.

## The actual radical core seals all supporting whites

Use context160 numbering. Its15-root core consists of the path

`4–3–2–1–0–5–14–13–12–11–10–6`,

with simple bonds except3–2 of label4, and branches6–7 and6–8–9. Roots3 and4 have squared gauge2 in the saved rational congruence; all other squared gauges are1. The black pairs are `{3,4}` and `{6,9}`. The white set is the remaining11 roots.

The complete core has signature `(13,1,1)`. Its rational-gauge radical is

`(-1,-2,-3,-2,-1,0,6,3,4,2,5,4,3,2,1)`.

The corner retaining3,6 omits4,9. Its13-node component F is the entire affine B12 factor, with a positive kernel. The checker verifies the full signature, the white positive definiteness, all four strict edge diagrams, the positive affine kernel, and each of the four corner signatures using exact rational arithmetic.

Let e be any additional supporting white root of the actual quadrilateral. The corner containing F is positive semidefinite. The positive kernel, together with nonpositive native off-diagonal products, therefore forces e to be orthogonal to every root of F. This is a statement about the actual corner; it does not assume simplicity of the ideal vertex.

The displayed indefinite core contains a timelike vector. Its Gram radical represents an actual zero linear combination of its roots: the represented vector is orthogonal to the core and has square zero, while the orthogonal complement of a timelike vector is positive definite. Pairing the radical with e consequently gives

`-<e,scaled_root4> + 2<e,root9> = 0`.

The only remaining possible core contacts are4 and9. They either both vanish or are both nonzero. In the latter case the adjacent marks4,9 have the actual path4–e–9 of length2, which is already paid. Hence, under the unpaid-face hypothesis, every additional supporting white is orthogonal to the entire15-root core. The reversed context1206 has the same proof, with its own saved radical and scaled omitted-root coefficients.

The actual corner diagrams therefore split as the displayed core corners plus a positive-definite supporting-white diagram. In particular each displayed finite component remains a whole component, with no indirect supporting tail or bridge. Each of the three other actual corners is finite, since its full diagram is a direct sum of positive-definite blocks.

## The two finite halves

For context160, the finite corners have these marked roles:

|Marks|Whole marked components|Payment|
|---|---|---:|
|3,9|B10 at doubled terminal3; A2 at endpoint9|new1/2|
|4,6|A1 at4; D12 at fork center6|existing1/2|
|4,9|A1 at4; A2 at endpoint9|0|

The unmarked extra A1 in the first corner and unmarked path component in the last corner do not alter these roles. The reversed context gives the same two halves. The source reconstructs each connected component, its finite type, all bond labels and the marked endpoint/fork roles from the rational-gauge Gram. Thus the two halves replace the removed ideal unit in both selected contexts.

The new rule contributes exactly1 per whole A2/B10 component pair at a finite vertex: there are two A2 endpoints and one doubled terminal. If a is the total rank of the eligible small components and v the number of B10 components, this contribution is at most `a*v/2`. This elementary count does not itself establish the joint finite-vertex budget.

## The universal D14 immediate alternative has a strict gap

The [affine-D14 long-pair theorem](fresh_affine_d14_long_pair_transfer_20261017.md) originally paid a connected affine B12 alternative when a marked fork center a is deleted and its replacement c attaches by label4 to the long endpoint of the remaining D12 block. Its positive affine kernel forces all other roots of that adjacent corner to be orthogonal to the block, so no additional outer contact is missing from this profile.

Here that immediate unit is unnecessary. Keep the original D14 factor F of order15, with marks a=0 and b in `{10,13,14}`, and use the same cusp coefficient normalization and Schur coordinates as the frozen proof. Put `x=-a.c>=0`. Direct elimination gives

`c0=0, l=sqrt(2), h=x+sqrt(2), g_c=2sqrt(2)/(x+sqrt(2))<=2`.

This bound covers every nonnegative opposite product, including all finite labels, parallel and dotted cases. The source independently reconstructs the original15-node Gram and these three exact eliminations; it does not merely copy the endpoint maximum.

At the opposite end, the complete inherited profile lists give `g_d<=6` for a fork center and `g_d<=8` for a short leaf. These bounds include the former immediate center profiles themselves: `inner4` has cost at most2 and `inner_next3` at most4. The arbitrary outer-label estimates from the original proof remain in force. Thus the argument does not depend on paying another immediate order13 alternative first.

If any adjacent native product or short spectator path already pays, there is nothing to prove. Otherwise the weighted-tree identity gives

`det(S)/(h^2*k^2) = 1-(g_c+g_d)/R`,

where R is10 for the opposite fork center and12 for a short leaf. Consequently this determinant is at least1/5 or1/6, respectively. Arbitrary positive-definite spectators add

`(1/R)*(x_Z/h-y_Z/k)^T G_Z^(-1)*(x_Z/h-y_Z/k) >= 0`.

The actual quadrilateral requires the determinant to vanish, a contradiction. Reflection gives the same conclusion when the deleted center is at the other end. Therefore the D14 theorem remains valid after removing the B12 singleton; this branch supplies an actual alternative payment rather than an unpaid new assumption.

## The distance11 orbit is recorded but retained

The same sealed core also underlies contexts166,170,1744,2137, which support the two canonical pairs `(10,11),(10,12)`. Contexts166 and1744 have a whole finite B13 corner and already receive a unit. Contexts170 and2137 instead have a whole A1 versus a short fork leaf of whole D12; the currently retained rule gives only a half there. A separately justified upgrade of that particular coupon to1 would pay these cases. The checker saves all four literal outcomes, but this package does not supply that upgrade or remove either distance11 pair.

## Verification scope

Generation43252 and a separate fresh saved replay89076 both terminated exit0/PASS for both singleton cases, all four conditional distance11 cases, their24 exact corner decompositions and strict edge tests, the three D14 immediate eliminations, and the complete saved-map and pair-set interfaces. Final checks passed for syntax, all16 direct source bindings and all4 report links. Input statuses and every recorded direct source binding are verified before use.

Two separate agents completed full source and written-mathematics reviews, finding no gap in the affine-kernel/radical sealing, whole-component coupon roles, independent immediate-D14 elimination, or the distinction between the removed singleton and conditional distance11 orbit. Neither review repeated the24-corner arithmetic or the frozen801-pattern census. The complete801-pattern embedding search,82,197-overlap construction and inherited all-label profile classifications are explicit frozen dependencies, not repeated finite censuses here. Source, data and report are final and frozen. No global budget, dimension improvement, external refereeing, proof-assistant formalization or geometric realization of a new chamber is asserted.
