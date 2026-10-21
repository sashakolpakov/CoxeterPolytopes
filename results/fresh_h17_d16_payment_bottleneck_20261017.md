# The D16 cusp payment bottleneck in the known H17 pyramid

The known finite-volume H17 pyramid realizes the difficult17-node affine D16 factor. With the retained finite rule (unit payment for distance at most9 **or connected component order at most17**) and unit payments at its other E8 ideal corners fixed, its actual faces force all136 D16 pair weights to be at least1 under any rule invariant under the intrinsic D16 diagram automorphisms. The resulting factor cost is136=8·17. This is a local obstruction to reducing that particular diagram-based budget, not a new dimension theorem.

The [exact checker](../scripts/fresh_h17_d16_payment_bottleneck_20261017.py) and [saved certificate](fresh_h17_d16_payment_bottleneck_20261017.json) reconstruct every two-face through this cusp from the already frozen [complete H17 cone and face ledger](fresh_nonincident_carrier_incidence_20261017.md). The chamber is the two-L25 diagram sharing node9. Its82 cone rays comprise80 finite vertices and two proper ideal vertices: the apex supported on9 has factors affine E8×E8; the vertex supported on{3,13} has the17-node affine D16 factor. The finite-volume, actual-facet and complete-cone assertions are inherited from that frozen chamber package, not inferred from a partial list of proposed corners.

## Exact face deficits

The D16 cusp is simple. Every pair among its17 incident facets determines one actual two-face, giving136 faces:72 triangles and64 quadrilaterals. For each pair the checker independently reconstructs the face vertices and actual side facets from the saved cone rays, then checks agreement with the frozen exhaustive face ledger. It computes each marked distance in the actual corner star, including the correct varying E8 factor at the nonsimple apex.

The fixed comparison rule is:

* at finite corners, weight1 if the marked distance is at most9 or the connected marked component has order at most17, otherwise0;
* at the other ideal cusp, weight1 on the actual nine-node E8 varying factor;
* at the D16 cusp, leave the weight as an arbitrary nonnegative real variable for its marked pair.

A triangle requires total payment2 and a quadrilateral requires1. The exact required D16 payment is the positive part of that demand minus the fixed payments at all other corners. Every resulting deficit is0 or1.

| Face type | Faces | Pairs with deficit1 | Pairs with deficit0 |
|---|---:|---:|---:|
|Triangle|72|60|12|
|Quadrilateral|64|49|15|
|Total|136|109|27|

The49 mandatory quadrilateral pairs have distance histogram

`2:1, 3:2, 4:3, 5:4, 6:5, 7:6, 8:7, 9:6, 10:5, 11:4, 12:3, 13:2, 14:1`.

In particular **15 pairs at distance greater than9 are required by actual quadrilaterals under the retained finite rule**. All other corners of these quadrilaterals are finite, so this conclusion is independent of any payment at the E8 cusp. It is not an artifact of a catalogue containing unrealizable stars.

For comparison, omitting the retained small-component finite payments and using finite distance9 alone gives62 mandatory quadrilateral pairs, including all26 pairs at distance greater than9. There are still60 mandatory triangle pairs, hence122 total. The small-component rule removes13 of these quad deficits (11 long pairs and two at distance9). The certificate saves both rules' corner payments and face deficits separately. The26-long-pair statement must not be transferred to the retained finite rule.

Fourteen of the60 mandatory triangles also contain the E8×E8 apex. Their deficits do depend on fixing the other ideal payment at1. Allowing a different weight there may change the D16 demands; this report makes no lower-bound claim for arbitrary simultaneous changes of both cusp rules. Every face, corner payment and deficit is saved explicitly, not only as an aggregate count.

## Diagram invariance forces all136 pairs

The D16 diagram has two fork-leaf interchanges and the reversal of its unique double-fork spine. These generate its full group of8 automorphisms. The checker verifies every generator against the literal graph, constructs the group, and saves for each of the136 pairs an automorphism taking it to one of the109 primary mandatory pairs. Since those109 are also among the122 comparison demands, the same witnesses establish closure for both baselines.

Consequently any intrinsic-diagram-invariant nonnegative real weighting satisfying the stated fixed-baseline face demands must give every pair weight at least1. Its cost is therefore at least136. Conversely unit weight on every pair satisfies these demands, so this local minimum is exactly136.

Without diagram invariance, the109 primary mandatory pairs may receive weight1 and the27 optional pairs weight0 in this particular chamber. Under the distance9-only comparison the corresponding counts are122 and14. These are statements about this complete fixed chamber and the specified fixed baselines, not a universal context-dependent weighting theorem. Additional ambient incidence information is therefore a possible way around the diagram-only obstruction.

## Relation to the nonincident E8 carriers

Every one of the64 quadrilaterals through the D16 cusp has all19 roots in its original face star. Both E8 factors of the other cusp occur in that star, but that cusp is not on the quadrilateral. Thus the required long-pair faces have explicit nonincident small-factor carriers at the same other cusp. This supplies concrete data for a future cross-cusp charging argument. It does not prove a bounded multiplicity or a global curvature inequality.

The [mixed finite9/ideal7 replacement lemma](fresh_ideal_distance7_mixed_replacement_20261017.md) and its improved triangle paths remain valid. This audit explains why those local improvements alone do not remove the order17 whole-factor payment cost in a diagram-invariant scheme.

## Verification scope

Generation and fresh saved replay both terminated successfully:136 reconstructed faces,109 primary and122 comparison unit deficits,8 automorphisms and136 saved orbit demands. Syntax, all7 direct source bindings and all4 local links pass. Independent internal source/proof review found no gap, and a separate read-only inspection recomputed every saved primary/comparison corner payment and face deficit and all136 automorphism witnesses. Those reviews did not repeat the original H17 cone or volume proof. This audit uses elementary exact graph and rational cone-incidence checks; it makes no new construction or dimension claim. No frozen predecessor was edited. Source, data and report are final/frozen.
