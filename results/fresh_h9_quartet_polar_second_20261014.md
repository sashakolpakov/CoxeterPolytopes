# The second polar iteration: twenty-five actual facets

Two more mutually compatible polar walls give an exact H9 Coxeter chamber with **25 actual facets, 81 finite vertices and 247 finite edges**. It has twelve more finite vertices than the preceding 23-wall chamber. The strict spacelike escape remains inside every wall, so the volume is still infinite.

The [source](../scripts/fresh_h9_quartet_polar_second_20261014.py) and [certificate](fresh_h9_quartet_polar_second_20261014.json) add exterior polar candidates 82 and 85 from the [preceding iteration](fresh_h9_quartet_polar_extension_20261014.md). All 300 wall pairs satisfy exact Coxeter relations; 47 pairs involve one of the two new walls. A common strict timelike center and 25 facet-interior witnesses certify the proposed facets. All 69 preceding finite vertices survive.

The complete exact cone census has 81 finite vertices, no ideal vertices and 513 hyperideal extreme rays. The finite-edge graph leaves 235 open finite-corner directions. Testing all 513 exterior polar hyperplanes against the complete finite-label alphabet finds 509 forbidden candidates and four compatible ones, indexed 62,95,96,128. All four miss the same strict escape. Their six mutual pairs are ultraparallel, as checked in an additional diagnostic, but they have not been added in this certificate.

The phase stops at this explicit polar frontier. Repeatedly adding walls that miss the escape has increased the finite geometry without resolving the volume obstruction. The four further caps, nonpolar extensions, and changes to retained walls remain possible future directions. This result supplies no new finite-volume example and makes no novelty claim for the infinite-volume chamber itself.

The [independent review](fresh_h9_quartet_polar_extension_20261014_review.md) reproduces all three complete cone censuses with a separate number-field implementation and independently coded cone clipping. It verifies all 1,211 polar candidates, the final facet and finite-edge geometry, and the surviving strict escape. Generation and saved replay pass.
