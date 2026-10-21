# One fewer paid pair in the order-17 affine D factor

The complete finite8/ideal5 payment rule can remove one additional pair from the 17-node affine D diagram if it adds a finite half-payment between a whole A4 component, marked at an endpoint, and a whole D13 component, marked at the unique long-arm neighbor of its fork center. The ideal factor cost decreases from 127 to 126. This is a local face-payment transfer; its finite-vertex budget and any resulting dimension bound are separate obligations.

Here the standard affine type is **affine D16**, with 17 nodes and affine rank 16. Labels such as “D17 target” in the inherited enumeration count diagram nodes. The new finite D13 coupon uses the ordinary finite type and therefore has 13 nodes.

The [exact checker](../scripts/fresh_d17_singleton_pair_coupon_20261017.py) and [saved certificate](fresh_d17_singleton_pair_coupon_20261017.json) bind the complete [2,614-context overlap certificate](fresh_distance8_high_overlap_20261017.md), its [1,153-pattern induced-map census](fresh_distance8_ideal5_small_factor_20261017.md), and the [complete finite8/ideal5 payment theorem](fresh_finite8_ideal5_combined_payments_20261017.md). They do not replace those completeness proofs with the single core below.

## Complete support diagnostic and the only affected context

Use the canonical affine D16 graph with spine 0–1–…–12, leaves 13,14 at 0 and leaves 15,16 at 12. The inherited paid set has 127 pairs, including 110 basic pairs. Its 17 remaining pairs split into six automorphism orbits:

|Representative|Distance|Orbit size|Selected contexts requiring the orbit|
|---|---:|---:|---:|
|0,10|10|2|2|
|0,11|11|2|2|
|1,11|10|1|1|
|1,15|12|4|4|
|2,15|11|4|4|
|3,15|10|4|16|

The certificate saves the actual context keys for every row. The singleton pair 1,11 is fixed by every one of the eight graph automorphisms. Removing it leaves an invariant 126-pair set containing the whole basic set and every elementary triangle pair. In particular, the removed pair has distance 10: it is not a short8-first payment, a distance-at-most9 feature payment, or an immediate whole-factor alternative of order at most13. The earlier literal D16/D17 transfers end in finite coupons, finite units or length-two paths and do not use it.

Across every saved marked image in the complete current census, only context 2518 requires this pair. Its full key is

`left L25 with marks{3,4}; right L25 with marks{3,4}; identify left9=right9`.

Its preferred corner is the physical pair 4,14. This exact whole-core match matters: the inherited overlap certificate excludes every additional cross-core edge, including simultaneous or larger native contacts. The local argument below therefore has the literal 19-root Gram as a proved input in this context. It is not asserted for an arbitrary abstract embedding of an affine factor.

## Literal core and all four corners

Number the left L25 by 0,…,9 and the right copy by 10,…,18, sharing 9. Every edge has label3. The left edges are

`0–1, 0–2, 2–3, 0–4, 4–5, 5–6, 6–7, 7–8, 8–9`,

with the corresponding right edges obtained by adding 10 to nodes0,…,8 and retaining the shared node 9. The opposite black pairs are{3,4} and{13,14}; all other displayed nodes are white. Arbitrary additional supporting whites E are allowed.

The checker verifies the full signature(17,1,1), positive definiteness of the 15-node white diagram and each of the four 16-node actual-edge diagrams. The four 17-node displayed corners are:

|Marks|Displayed corner|New half-payment|
|---|---|---:|
|3,13|finite A4 + A9 + A4, marks in different A4 components|0|
|3,14|finite A4 + D13; A4 endpoint 3, D13 long-arm neighbor 14|1/2|
|4,13|finite D13 + A4; D13 long-arm neighbor 4, A4 endpoint 13|1/2|
|4,14|affine D16 of order 17; marked distance 10|0|

The finite D13 fork center is 10 in the second row, with short leaves 11,12 and long-arm neighbor 14. It is 0 in the third row, with short leaves 1,2 and long-arm neighbor 4. Thus these roles are distinct from the earlier fork-center/short-leaf coupons. The two new half-payments sum to 1.

For an independent realization check, the certificate matches the actual noncompact quadrilateral with sides 3,4,13,14 and four vertices in the frozen [finite-volume H17 pyramid ledger](fresh_nonincident_carrier_incidence_20261017.md). That example is corroboration; the following proof applies to any actual quadrilateral having this literal marked core.

## Additional supporting whites either pay immediately or preserve both coupons

Let F be the 17-node affine diagram obtained by omitting 3 and 13. The actual corner 4,14 contains F and is ideal. Its positive affine kernel has coefficients 1 at its four leaves and 2 on the spine. Every extra supporting white e belongs to this corner. Its nonpositive contacts with F, together with the positive kernel, force e to be orthogonal to each root of F. Equivalently, the represented positive null vector is proportional to the ideal corner ray, to which every incident wall is orthogonal. This uses the actual corner and does not assume a simple ideal vertex.

The full 19-root Gram has radical

`(-6,-3,-4,-2,-5,-4,-3,-2,-1,0,6,3,4,2,5,4,3,2,1)`.

It is an actual linear relation among the normals. Its represented vector is null and orthogonal to the indefinite core, hence orthogonal to a timelike vector; that orthogonal complement is positive definite, so the represented vector vanishes. Pairing the relation with e and using e perpendicular to F gives

`-2(e·n3) + 2(e·n13) = 0`.

If the equal contacts are nonzero, the actual corner 3,13 contains the white-interior path 3–e–13. It has a unit payment: finite distance 2 is paid, and at an ideal corner the actual varying-factor lemma places the connected marked pair in its affine factor, where distance 2 is paid. This already supplies the quadrilateral's required unit.

Otherwise every extra supporting white is orthogonal to all 19 core roots. The extra-white diagram is positive definite as a principal subdiagram of the actual face support. At each of the three displayed finite corners the full corner diagram is therefore the orthogonal sum of its positive-definite displayed diagram and this extra-white diagram. These are actual finite corners, and the displayed A4 and D13 components cannot be enlarged or joined by extra whites. Both coupon components are whole finite components, with precisely the endpoint and long-arm-neighbor marks listed above. The two half-payments survive and give the required unit.

Consequently the local transfer is valid with arbitrary extra supporting whites and nonsimple ideal vertices. Additional ambient ideal walls do not enlarge the whole finite components used by the coupons.

## Composition interface and verification scope

Replace only the order 17 affine D paid set by the saved 126-pair set and add only the stated A4-endpoint/D13-long-arm-neighbor half-coupon. The complete predecessor's triangle coverage and short8-first replacements remain valid. Every other selected image is retained or has its existing transfer; the sole new omission is the proved context 2518 above. Thus triangles still receive at least 2 and quadrilaterals at least 1, conditional on the already frozen complete predecessor theorem. The ideal pair cost for this target is 126. Establishing the total finite cost of the enlarged coupon rule is deliberately left to a separate budget proof.

Generation and separate fresh saved replay both terminated exit0/PASS. Syntax, all19 direct input hashes and all6 report links passed. Root completed a full source/mathematical read, including the imported literal Gram, inertia and component helpers, and found no gap in the affine-corner orthogonality, actual radical relation, length-two alternative, whole-component coupon roles or complete saved-map support diagnostic. That review did not duplicate the orbit/core arithmetic or the inherited overlap and embedding searches.

This checker reconstructs the six support orbits from the saved complete maps and directly verifies the one core's rational signature, positive kernel, four strict edge diagrams, four corners and coupon roles. It does not rerun the 82,197-identification or 1,153-pattern searches. Published diagram classifications and the predecessor's complete geometric arguments remain explicit inputs. Source, data and report are final/frozen. No new ambient dimension bound is asserted here, and no frozen predecessor was changed.
