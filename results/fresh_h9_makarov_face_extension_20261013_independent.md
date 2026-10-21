# Independent replay of the new 12-wall H9 seed

The saved 12-wall geometry passes a separate exact arithmetic check. Its new cap cuts the previous escape, preserves all ten recorded finite vertices, and has an actual facet. The resulting chamber still has a strict spacelike escape, so this is **not** a finite-volume H9 example.

The [independent checker](../scripts/fresh_h9_makarov_face_extension_20261013_independent.py) reads the [frozen seed certificate](fresh_h9_makarov_face_extension_20261013.json) directly. Its [audit output](fresh_h9_makarov_face_extension_20261013_independent.json) binds the source and arithmetic hashes. It imports no arithmetic or geometry routines from the seed producer.

The arithmetic uses the separate sparse squarefree-radical implementation, with the exact conversion

\[
\sqrt{(336+128\sqrt5)/11}
=\frac{16\sqrt{11}+4\sqrt{55}}{11}>0.
\]

Thus the producer's two extension layers become ordinary expressions in \(\mathbb Q(\sqrt5,\sqrt{11})\). The square-root identity is checked before decoding any normal.

Verified independently:

- The ambient form has signature \((9,1)\), by exact LDL elimination.
- All 12 normal norms and all 66 pairs agree with a Coxeter chamber: 53 right pairs, eight label-3 pairs, two label-5 pairs, and three strictly dotted pairs.
- The common center is timelike and strict. The saved projection onto the new cap is timelike, orthogonal to that cap, and strict at all other walls. Acuteness and the common strict center also preserve the old facets.
- Each of the ten saved vertices is future timelike, feasible for all 12 walls, and has exactly nine incident walls with positive-definite Gram matrix. The new cap is strict at every one.
- The old escape is cut strictly. The new integer escape is future spacelike and strict at all 12 walls; its exact norm is
  \(316589058752216839496570-140761570367036752425972\sqrt5>0\).
- All 40 strict squared inequalities for the earlier vertex-center/fundamental-weight families pass. Dividing each family normal by its transverse magnitude gives a decreasing escape pairing for every magnitude at least one, proving the stated all-label obstruction within those families.

Generation and saved replay **PASS**. This check does not repeat the compact H5 extreme-ray census or the separate 200-profile screen. It makes no novelty or unrestricted completion claim.

```sh
python scripts/fresh_h9_makarov_face_extension_20261013_independent.py --verify
```
