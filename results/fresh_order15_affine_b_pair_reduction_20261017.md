# Three removable payments in the order-15 affine B factor

The order-15 affine B factor can use **97 paid pairs instead of 100**, with no new finite coupon. The three omitted pairs were required only by seven selected overlap contexts. In six contexts, the marked component is a whole finite D14 that no supporting white can touch. In the seventh it remains a whole finite A12 or grows to a whole finite A13. Every case therefore already has a finite unit payment.

This is a local improvement to the complete finite8/ideal5 face-payment theorem. It is not a new ambient dimension bound. Standard type names are used below: the affine target is **affine B14**, of rank14 and order15; the finite D14 has14 roots.

The [checker](../scripts/fresh_order15_affine_b_pair_reduction_20261017.py) and [certificate](fresh_order15_affine_b_pair_reduction_20261017.json) bind the [complete high-overlap domain](fresh_distance8_high_overlap_20261017.md), its [complete selected-pattern map census](fresh_distance8_ideal5_small_factor_20261017.md), and the [finite8/ideal5 composition](fresh_finite8_ideal5_combined_payments_20261017.md). The all-label supporting-white classification is an explicit frozen input, not inferred from a bounded search in this new checker.

## Exact support diagnostic

The cost-maximizing family at order15 is affine B14: its cost100 exceeds A75, C95 and D96. In its canonical numbering the fork center is0 with short leaves13,14; the spine runs0–1–…–12, with the last bond11–12 doubled. There are two automorphisms, exchanging13 and14.

The existing100-pair set consists of91 basic pairs and these six extra orbits. A context is counted once in a row if at least one of its saved marked images meets that orbit.

|Representative|Distance|Orbit size|Supporting contexts|
|---|---:|---:|---:|
|0,10|10|1|13|
|0,11|11|1|3|
|1,11|10|1|11|
|9,13|10|2|53|
|10,13|11|2|17|
|11,13|12|2|4|

The certificate saves every supporting context's whole marked-core key. It also records the analogous orbit diagnostics for B/C/D at orders14 and16; those are diagnostics only, with no pair deletion claimed there.

Remove the invariant union

`{(0,11), (11,13), (11,14)}`.

The exact complete saved-map support is:

|Removed orbit|Affected current context indices|
|---|---|
|0,11|1284, 1623, 1641|
|11,13 and11,14|1646, 1666, 1933, 2326|

All three pairs have distance greater than9. The full basic set, elementary triangle payments, and the short8-first replacement interface remain unchanged. The old small-factor transfers use their own existing finite outcomes or retained short paths; they do not require a new basic pair here.

The previous map search deliberately allowed every induced embedding of the bounded marked component into a target. Such embeddings need not all occur as actual affine completions. The following literal supporting-white checks prove that the seven particular preferred components are always finite. No conclusion is drawn from an abstract target embedding alone.

## Six reversed-three-node contexts: the marked D14 is untouched

Contexts1284,1623,1646,1666,1933,2326 all have two L25 cores identified along the reversed path

`left7=right9, left8=right8, left9=right7`.

Use the saved numbering with right nodes0,…,6 becoming10,…,16. The17-root Gram has signature(16,1,0). The complete frozen no-short8 supporting-white classification says that any extra supporting white with a nonzero core contact has exactly one unit contact, at node3 or node13. Both projection squares are2; every other nonzero support or larger native contact is already excluded by the exact corner and Lorentz tests. Additional core-orthogonal whites are unrestricted subject to the actual positive-definite face support.

For all six contexts, the preferred marked component is a finite D14 and contains neither3 nor13. The new checker reconstructs its literal graph and positive definiteness. Consequently every extra supporting white is orthogonal to every root of this D14. No chain of extra whites can connect to it, since any such chain would need a first direct contact. The component is therefore whole and finite, with14 roots.

More explicitly, the complement of this component in the preferred corner consists entirely of supporting whites: both black marks are in D14. That complement is positive definite as a principal subdiagram of the actual face support. The complete preferred-corner Gram is thus an orthogonal sum of positive-definite blocks, and this is an actual finite corner. The whole-component-at-most17 rule gives a unit. This argument does not presume a simple ideal vertex or treat a finite subdiagram at an ideal corner as sufficient by itself.

The absence of extra cross-core edges is supplied by the frozen complete overlap certificate for the literal whole context. It is essential to the stated contact domain and is not assumed from the preferred component alone.

## Context1641: a residual-zero endpoint can only give a whole A13

The remaining context has two L25 cores identified along the reversed path

`left6=right9, left7=right8, left8=right7, left9=right6`,

with opposite black pairs{1,2} and{11,12}. Its16-root Gram is nondegenerate Lorentzian, of signature(15,1,0). The preferred marks1,11 are endpoints of the displayed A12 path. The complete supporting-white domain consists of no core contact, or one unit contact at one of1,3,11,13. In that order, the corresponding inverse block is

```
  2       1       -3       -2
  1       4/3     -2       -4/3
 -3      -2        2        1
 -2      -4/3      1        4/3
```

Only1 and11 belong to the marked A12. A white root attached at either of these nodes has projection square2, hence residual square0. The core spans a Lorentzian subspace, whose perpendicular is positive definite. The residual vector therefore vanishes: this root lies in the core span and is orthogonal to every core-orthogonal white.

It also cannot coexist with a second core-contact white. Its actual product with such a root equals the displayed projection product. For the same endpoint or the other same-end contact, this product is positive, contrary to acuteness. At the opposite end the product is−3 or−2. Two actual supporting whites have a positive-definite principal Gram, so their norm2 product must lie strictly between−2 and2. Those opposite-end possibilities are excluded as well. This argument includes distinct roots with the same contact; no multiplicity assumption is made.

If neither endpoint-contact root exists, every extra white is orthogonal to A12, so it remains whole. If one exists, it has a single simple bond at an endpoint and gives the literal whole A13. It cannot connect to any other extra white by the preceding two arguments. The checker verifies both A13 extension matrices and every mixed product used above. Both marks remain in this component; its orthogonal complement in the corner is contained in the actual positive-definite white support. The corner is finite and the existing whole-component rule pays a unit in either case.

## Local consequence and quantitative scope

All seven newly omitted selected images have a retained finite payment under arbitrary supporting-white extensions. Every other selected image is unchanged. The complete predecessor theorem may therefore use the saved97-pair affine B14 set while preserving every triangle and quadrilateral payment. No extra finite weight is introduced and the finite8 vertex budget remains unchanged. The other order15 costs remain A75, C95 and D96, so the worst order15 cost becomes97.

The previous dimension75 boundary profile `(12,12,12,12,12,14)` has normalized margin `-16873/171366`. Reducing the order15 cost by one pair changes this margin to `-9081/285610`; two pairs change it to `29879/856830`; all three change it to `87001/856830`. These are comparisons for this one formal profile only. They neither certify all profiles nor prove that the formal profile is geometrically realizable. A new global bound requires a fresh complete scalar comparison and valid concentration rules for the changed envelope.

## Verification status

Final generation and separate fresh saved replay both terminated exit0/PASS. The source checks every input's completed status and direct recorded hashes. It reconstructs each actual corner from the four black roles and verifies that the saved marked component is its literal connected component, contains both marks, and has only supporting whites in its complement. All17 direct bindings, syntax and report links pass.

Two separate full source/mathematical reads found no gap in the six untouched D14 cases, residual-zero A12/A13 argument, strict white-PD exclusions, complete saved-map affected support or basic-pair retention. Root also read the inherited ordinary-extension classification and requested the explicit component and input-hash guards now present. Neither review duplicated the seven-case arithmetic or the inherited overlap and map searches.

The new arithmetic reconstructs graph automorphisms and exact support orbits from the saved complete map corpus, and verifies seven literal finite outcomes using rational Gram inverses and signatures. It does not repeat the82,197-overlap or1,153-pattern searches. Published diagram classifications and the complete predecessor's all-label support proof remain explicit inputs. Source, data and report are final/frozen. No frozen predecessor was edited, and no new dimension theorem is asserted here.
