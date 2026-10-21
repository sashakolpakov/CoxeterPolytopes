# Four polar caps give a twenty-three-facet chamber

Adding all four compatible polar walls gives an exact **23-facet H9 Coxeter chamber with 69 finite vertices and 213 finite edges**. Every old finite vertex survives, and sixteen new finite vertices appear. The old strict spacelike escape survives all four inequalities, proving infinite volume.

The [source](../scripts/fresh_h9_quartet_polar_extension_20261014.py) and [certificate](fresh_h9_quartet_polar_extension_20261014.json) verify all 253 wall pairs, including 82 new pairs; all six cap-to-cap pairs are ultraparallel. A common strict timelike center and 23 explicit relative-interior facet witnesses certify that every proposed wall is an actual facet. The complete cone census has 69 finite vertices, no ideal vertices and 413 hyperideal extreme rays. Its finite-edge graph has 195 unfinished finite-corner directions.

Applying the [complete finite-angle test](fresh_h9_quartet_polar_20261014.md) to every exterior polar hyperplane of the enlarged cone finds two further compatible candidates, indices 82 and 85. Both retain the center, both miss the strict escape, and their mutual relation is ultraparallel. The remaining 411 exterior polar candidates have a forbidden finite angle. These next two walls are pursued in a separate iteration so that this certificate stays fixed.

This is a geometric extension of an infinite-volume nonarithmetic chamber, not a new finite-volume example. No novelty claim is made for the infinite-volume chamber itself. Nonpolar walls and changed retained walls remain open construction directions.

The [independent review](fresh_h9_quartet_polar_extension_20261014_review.md) reproduces all three complete cone censuses with a separate number-field implementation and independently coded cone clipping. It verifies all 1,211 polar candidates, the final facet and finite-edge geometry, and the surviving strict escape. Generation and saved replay pass.
