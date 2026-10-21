# Removing the final two nonbasic order13 affine B pairs

Starting with the completed ceiling65 payment rule, the two remaining
nonbasic affine B12 pairs may be removed. Its order13 cost falls from75
to73. No finite weight or coupon changes. The full whole-A2/whole-D11
short-leaf coupon already proved in the order12 B/C theorem is essential
to two of the supporting cases.

The order13 A/B/C/D costs become65/73/73/74, so the normalized rank12
envelope falls from75/13 to74/13. This is a local payment refinement;
it gives no new ambient bound. Other scalar and finite obstructions remain.

The [source](../scripts/fresh_order13_affine_b_fork_projection_transfer_20261017.py)
and [certificate](fresh_order13_affine_b_fork_projection_transfer_20261017.json)
re-enumerate the complete selected-map domain and check every affected
actual-component interface.

## Intrinsic pair rule and complete domain

Number the affine B12 spine0--1--...--10, with doubled bond9--10 and
short leaves11,12 at0. Delete(9,11) and(9,12), both at distance10. They
form one orbit under the diagram's leaf swap. The earlier singleton
(0,10) and terminal/short-leaf pairs(10,11),(10,12) stay absent under
their already completed proofs. The retained73 pairs are exactly the
basic ideal5/feature5 set, all with weight1.

Every short8 and elementary-triangle pair is retained, and no distance9
payment changes. The source checks the intrinsic basic set, all
automorphisms, shortest paths and the full triangle pair set inherited
from the order13 singleton theorem.

The frozen complete order13 diagnostic contains801 eligible connected
selected patterns, including those whose actual-component upper bound
is larger than13. The new source reruns exhaustive injective induced
embeddings of every pattern into the13-node B target, checking every
edge and nonedge and every failed pattern. Its output agrees with the
complete saved maps. Exactly23 contexts, with52 affected maps, contain
one of the changed pairs:

    149,151,625,630,638,715,721,1356,1364,1556,1642,1648,1668,
    1742,1849,1857,1922,1935,2135,2242,2250,2315,2328.

Four contexts already contain the whole affine factor. Nineteen describe
proper selected subdiagrams and require additional roots to complete
the target. These are different geometric situations and have separate
proofs below.

## Nineteen proper-completion contexts: impossible doubled terminal

Every affected proper core H is nonsingular Lorentzian, all its saved
roots have physical norm2, and all its nonzero off-diagonal entries
are−1. The checker recomputes its exact inverse J and actual preferred
corner. The complete necessary contact domain of any supporting white
is partitioned into two ends. Their inverse submatrices are entrywise
positive. A root meeting both ends gives a retained marked path of
length at most8.

In each of the44 maps, target node9 is an old selected marked root i,
while the doubled terminal10 is absent. An actual whole affine B12
factor would require a norm2 root e with label4 contact to i. That root
is outside H: old omitted black roots do not belong to this corner,
and every other old corner component is orthogonal to the selected
component. It therefore cannot supply the required nonzero first bond.
Other target nodes may be supplied later by arbitrary old or new
components; the proof uses only this forced first terminal contact.

The native squared magnitude of the label4 contact is2. Since the
orthogonal complement of H is positive definite,

    projection_H(e)^2 <= 2.

If e meets both ends the short-path branch already pays. Otherwise its
contacts lie on one end, where positivity gives the lower bound2*J_ii.
Every extra contact on that end only increases the bound. The exact
table over all44 maps is:

|Forced norm2 projection lower bound|Maps|
|---:|---:|
|5/2|12|
|3|12|
|7/2|8|
|7|4|
|15|8|

Every entry exceeds2, excluding all such actual completions. No
unspecified contact is set to zero. This argument covers arbitrary
native labels, multiplicities and indirect spectator chains. It does
not require the more complicated relevant-component tail bound merely
because some saved contexts were originally classified by that method.

## Four complete affine factors: radical sealing and an existing unit

Contexts149,151,1742,2135 have15-root cores of signature(13,1,1).
Their selected13-node component F is the whole affine B12 factor. The
checker reconstructs its positive kernel, the full-core radical, the
strictly positive-definite white and edge diagrams, all four corners
and their connected finite components with native labels.

For any extra supporting white e, positive semidefiniteness of the
actual corner containing F and the positive affine kernel force every
contact with F to vanish. The full-core radical is an actual zero
linear combination: the core contains a timelike vector, and a null
vector orthogonal to it must vanish. Its two coefficients at the
omitted black roots are nonzero with opposite signs. Pairing e with
this relation says its two remaining core contacts either both vanish
or both are nonzero.

Two nonzero contacts give an unchanged adjacent marked path of length2.
Otherwise every supporting white is orthogonal to the whole old core.
Each finite corner is then the displayed positive-definite corner plus
an orthogonal positive-definite white diagram. In particular its marked
components are whole actual components: neither a direct attachment
nor an indirect bridge can enlarge them.

The reconstructed finite payments are:

|Contexts|Whole finite alternatives|Total after deleting the ideal unit|
|---|---|---:|
|149,1742|B13 unit and A2/D11 short-leaf unit|2|
|151,2135|A2/D11 short-leaf unit|1|

In every A2/D11 corner the small component is exactly A2 and the large
mark is a short fork leaf of exactly D11. The source verifies these
roles from the complete connected components and the positively scaled
native Gram. It uses the already completed
[A2/D11 unit upgrade](fresh_order12_bc_unit_coupon_transfer_20261017.md),
whose cost is already included in the ceiling65 finite budget. No new
coupon is added here. The B13 units are also unchanged.

## Nonselected and later provenance

The complete short8-first reduction keeps all its original units and
universal transfers. Their ideal alternatives have distance at most9
with the retained nearby features, or are direct/shared-neighbour paths.
Their finite outputs remain unchanged. The low-core closure has whole
marked-component bound at most12 and cannot yield an order13 factor.
The disjoint high-core branch is excluded by its inherited all-label
argument. Protected high outcomes have retained finite payments or
basic ideal payments at distance at most9. The selected cases whose
weights can change are therefore exactly the23 cases above.

Universal affine-D14 immediate B alternatives have marked images
(0,10),(10,11),(10,12), none of the newly removed pair orbit. Their
existing strict-gap replacement remains unchanged. Its other immediate
alternatives are affine D, not B. Universal D13 immediate factors have
order12. The complete saved immediate map interface is checked directly.

The later universal D12 and D11 transfers end in finite alternatives
or ideal paths of length at most4 and5. Cycle12/13/14 and longer-cycle
alternatives end in their stated finite units/coupons or retained short
ideal paths. The older same-order and higher-order partial-weight
replacements retain their stated affine D weights; their other outputs
are finite. The B14/B15 projection exclusions introduce no replacement
payment, while their finite-completion predecessors use the stated
small whole finite components.

The newest order12 B/C theorem ends in finite alternatives and consumes
no order13 B unit. The order14 B refinement is a projection exclusion.
The large-D/D10 theorem uses retained short8/basic ideal paths and
protected distance9 outcomes, so its proof does not lose a payment here.
Its finite changes affect whole D_n of order at least20 and D10 coupons;
they leave the whole B13 and A2/D11 alternatives above unchanged.
The finite rank17 deletions likewise do not affect these components.

Thus all triangle payments remain at least2 and all quadrilateral
payments at least1. Actual cusp-fiber multiplicities stay unchanged;
only the order13 B factor weight sum changes, from75 to73.

## Verification and scope

Generationbfb9cf/8dee57 and separate fresh replaybd3c11/a7cb3c both
terminated exit0/PASS on all801 pattern searches,
all52 affected maps,19 complete necessary contact domains,44 exact
projection bounds, four radical/sealed corner configurations, and
the retained basic, short8, triangle and immediate interfaces. The
complete all-label geometric core census, no-extra-cross-edge proof,
finite/affine classification and arbitrary-support contact-domain
theorems remain explicitly inherited frozen inputs. They are not
replaced by an abstract finite embedding list.

Two separate full source/report mathematical reads passed. Both reviewers
reread the complete contact-domain/projection helper and the radical,
kernel and finite-component helper. They checked the physical norm2
interpretation, positive inverse end blocks with arbitrary extra contacts,
the genuinely new forced terminal, whole-component sealing and existing
A2/D11 unit, and all retained short/basic/low/immediate/later interfaces.
Root also checked all four finite totals and the exact existing coupon
budget dependency. Neither reviewer repeated the801-pattern census,
44 projection calculations or inherited all-label geometric census.
All34 direct bindings, Python AST and three local report links pass.
Source, certificate and report are final/frozen.

No finite-forest recomputation is required by this unchanged finite rule,
but any new ambient claim still needs its complete ideal comparison and
remaining finite obligations. No new bound, construction, external
refereeing or formal verification is claimed here.

    python3 scripts/fresh_order13_affine_b_fork_projection_transfer_20261017.py
    python3 scripts/fresh_order13_affine_b_fork_projection_transfer_20261017.py --verify
