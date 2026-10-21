# Independent review of the first-star metric obstruction

**PASS.** The [primary theorem](fresh_h9_star_cube_metric_20261014.md) excludes every metric completion of **the first saved mandatory assignment of star case 0**. It does not exclude the star incidence or its other mandatory assignments. No mathematical or scope issue was found.

The separate [review checker](../scripts/fresh_h9_star_cube_metric_20261014_review.py) and [certificate](fresh_h9_star_cube_metric_20261014_review.json) import no producer arithmetic. They use SymPy domain matrices over `QQ(sqrt(2))`, independently of the producer's scalar-pair field implementation, and reproduce:

- All 30 actual-vertex nonpositive-minor witnesses: 28 least-label-3 tests forcing right contacts and two least-label-4 tests bounding the remaining anchor contacts by labels 2 or 3.
- The positive rank-eight block and its exact inverse, plus the eight literal leaf-budget supports. The inverse is entrywise nonnegative; the own contact contributes `1/4`, and every foreign diagonal coefficient is 2. Two nonright contacts or one label at least 4 already cost at least `5/4`. Thus the seven projection patterns per leaf are exhaustive over all finite integer labels.
- All 49 negative anchor determinants, with no degenerate branch.
- All 1,372 nonzero residual three-by-three determinants: 150 negative and 1,222 positive. An additional 49 direct full eleven-by-eleven determinants check the block-determinant translation.

The monotonicity arguments are valid because all normalized off-diagonal entries are nonpositive. Replacing unspecified meeting contacts by right angles decreases the nonnegative weighted adjacency matrix. If even a least nonright contact destroys positive definiteness at an actual vertex, restoring any larger contact cannot repair that vertex. The tested anchor pairs are incident at actual vertices, so they must have finite labels; a dotted or parallel value is not an omitted option.

After eliminating the positive rank-eight block, three remaining vectors would lie in a space of dimension at most two in an H9 realization. Their nonzero Gram determinant forces rank three, hence full principal rank eleven. This is an exact rank obstruction independent of every unused dotted distance and the two unused central normals. The two-dimensional Lorentz description is consistent with the negative anchor determinants, but the final contradiction already follows from rank.

Generation and saved replay pass. The review binds the primary source and certificate and checks their input-source manifests; it does not repeat the complete 240-recipe geometry construction or every other star assignment.

```sh
python3 scripts/fresh_h9_star_cube_metric_20261014_review.py --verify
```
