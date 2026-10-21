# Independent source and normalization review of the complete coface scan

**PASS: no issue found in the stated source, normalization or scope.** The [complete scan source](../scripts/fresh_sqrt5_complete_coface_base_cones_20261017.py) applies the full small disjoint-base cones to every exported coface slice. It reports zero violations and no changes to face or vertex masses. This review did not repeat all 65,069 positive-slice load calculations; it independently reviewed their implementation and binds the saved result.

Bound inputs:

* Scan source SHA-256 `b3a80f9d4cfe111a5289e96ba8f752d8a8283c5308a839553468e3ec0553d2c9`.
* [Scan data](fresh_sqrt5_complete_coface_base_cones_20261017.json) SHA-256 `10d99e69e9ae6871f5987ce606d143fed91596e3331c7b0a59fc164c55581b5e`.

The generic cone description received a separate [independent exact reconstruction](review_simplex_base_packing_cone_20261017.md), including all active-matrix ranks and all packing/facet checks. The [labelled slice producer](fresh_sqrt5_labelled_coface_base_packing_20261017.md) already reconstructs the full actual base sets, simultaneous coface symmetries, and finite spectator-copy multiplicities; its source and normalization also received a separate narrow review.

## Mapping and arithmetic checks

The scan reconstructs the provenance-to-face map from the actual base records and uses the same common positive face masses throughout. Each nonzero block demand is recalculated from its saved rational provenance coefficients. It checks equality with the exported marginals, their sum, and the total mass of the literal source bases. All 6,016 bases and 1,858 actual provenances remain in the input, including zero-mass provenances.

The canonical-coface demand is an orbit-total measure: it sums incidence over the whole coface orbit after transport and stabilizer averaging. Consequently its capacity is the coface multiplicity times the containing simplex mass. The scan uses exactly that convention, with no extra division by the orbit size. The simplex mass and capacity are independently compared with the saved values for every slice.

For a ground set of size \(d+1\), it uses the complete cone of that size. Its common-denominator conversion is exact. Every coefficient and the cone capacity are doubled, so the maximum integer dot product divided by twice the common denominator is precisely the minimum mixture mass. Taking the maximum over all saved cone rows is valid by the complete LP-duality proof, not an assumption that a numerical optimizer has found the maximum.

The 33,650 zero-demand slices are discharged analytically: all demands and all facet coefficients are nonnegative and capacity is nonnegative, so every load is zero. The remaining 870 slices comprise 382, 101, 172 and 215 slices in simplex dimensions one through four. Their row counts give

\[
382\cdot1+101\cdot1+172\cdot8+215\cdot294=65,069
\]

positive-slice inequalities. All 34,520 slices together represent 510,656 inequalities. The saved maximum capacity utilization is exactly one, with zero violations.

## Consequence and limit

Membership in the complete marginal cone proves that **every nonnegative charge depending only on a base subset within one of these symmetrized coface classes** obeys its disjoint-packing capacity at the unchanged witness. Thus this conclusion is stronger than passing a selected finite family of subset cuts. No explicit mixture is needed for this logical implication, since the independently checked complete finite facet description and exact loads establish membership.

The common symmetry averaging is compatible with whole packings: transport all bases of a packing by the same automorphism, then mix. It is not independent transport of incompatible base instances. A mixture can also be symmetrized under the coface stabilizer without changing its invariant demands or total mass.

This exhausts the **abstract disjointness capacities for these exported, symmetrized slices**. It does not prove geometric realizability of all abstract packings, compatibility of choices across different cofaces or simplices, or existence of an ambient Coxeter polytope. Additional geometry can still forbid disjoint bases from coexisting. No source mass changes, inherited-cover reuse at changed masses, new global optimization or dimension bound is claimed.

**FINAL/FROZEN.** This is an internal mathematical/source review. The producer reports generation and fresh saved replay; this note does not misattribute them as a second independent full scan.
