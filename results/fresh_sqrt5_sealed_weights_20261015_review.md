# Independent review of the component-local additive-weight barrier

**PASS.** The [independent checker](../scripts/fresh_sqrt5_sealed_weights_20261015_review.py) and [saved certificate](fresh_sqrt5_sealed_weights_20261015_review.json) reconstruct the exact pair-orbit inequalities and all fifteen local compact quadrilaterals used by the [primary barrier](fresh_sqrt5_sealed_weights_20261015.md). Generation and saved replay pass.

The scope is essential: the baseline is fixed, corrections are nonnegative and depend only on the entire connected spherical component and its marked-pair orbit, and disconnected marked pairs receive zero. This is not an obstruction to changing the baseline, paying disconnected pairs, depending on the full vertex forest, or imposing simultaneous global face constraints. Nor does it construct a compact chamber with the rank-29 forest.

The review uses standalone Fraction arithmetic over Q(sqrt5) from the previously independent contact review. It imports no producer arithmetic, canonical tree-key routine, LP optimizer, numerical geometry, SymPy or NetworkX. Labeled graph bijections and all automorphisms reconstruct the 94 correction variables in 16 component types, their exact pair multiplicities, and every one of the 112 face inequalities. The saved rational primal is checked against all constraints.

A direct sparse lower certificate avoids reliance on an optimizer or a duality implementation. The nine E8 faces each force a different long-pair orbit weight to be at least one. The six D7 faces force six orbits whose multiplicities are 1,1,1,2,2,2. Thus the forced additional totals are exactly nine in both components. Feasible saved rational assignments attain these minima; the totals are E8 >= 23 and D7 >= 21. The E8 lower bound forces the common deficit to be at most -7, and the full checked primal attains -7. Consequently D7+D7+D7+E8 has weight at least 86, above the proposed rank-29 budget 84, within precisely this component-local rule.

For each of the fifteen faces, the review independently rescales its literal Coxeter Gram by a tree-product construction, checks physical and conjugate inertia, and projects its positive definite white support. The four residual black normals span a Lorentzian three-space. A separate rectangular-RREF calculation examines all six black-pair intersections: exactly four adjacent pairs give feasible timelike rays, and both opposite pairs have mixed-sign inequalities in either orientation. Across the package these are 60 finite corners and 30 rejected opposite intersections. The sum of the four rays is a strict timelike center, with every corner on its future sheet.

These checks establish compact quadrilateral geometry: the residual inequality cone has strict interior and no lineality, its extreme rays are exhausted by the six pair intersections, and its four extreme rays are future timelike. Perturbing a face-interior point slightly in the negative white-dual direction gives a strict center for the full local chamber. This is an actual local face statement, not a claim of compactness or finite volume for the entire chamber or its arbitrary orthogonal augmentations.

The complete union enumeration and the conditional unpaid-contact sealing theorem are bound through their frozen source manifests and earlier independent review; they are not repeated here. In particular, the 112-core sealing premise remains conditional on the absence of an H4 payment. The barrier itself needs only the fifteen explicitly reconstructed local face witnesses.

```sh
python3 scripts/fresh_sqrt5_sealed_weights_20261015_review.py --verify
```
