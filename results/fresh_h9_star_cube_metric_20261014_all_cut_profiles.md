# Complete cut-contact metric gate for the first star incidence

The compact H9 star incidence has **no realization with all four B2 blocks labeled 4 and all four positive leaves having own contacts (3,3)**. The four cut leaves are unrestricted within this statement. This closes a metric branch, not the incidence: changing a positive-leaf contact or a B2 label remains open.

The source incidence is case 0 of the [33-case tree catalogue](fresh_h9_tree_cube_wedges_20261014.json). Its positive leaves are `[2,5,8,11]`, cut leaves `[14,15,16,17]`, and four B2 blocks `[3,4]`, `[6,7]`, `[9,10]`, `[12,13]`. Their union is an eight-dimensional positive subspace. Every Coxeter Gram for H9 has rank at most ten, so its residual Gram has rank at most two.

For each cut, its missing triangle with the B2 block has a nonright own contact. If the cut met its positive leaf nonorthogonally, an actual three-wall finite face would have a nonpositive all-ones quadratic form. Thus that pair is right. The actual three-node paths then exclude every own label at least 6. The Lannér triangle condition leaves exactly eleven ordered profiles: the nine pairs from `{3,4,5}²`, together with `(2,5)` and `(5,2)`.

Both choices of the retained own root occur in actual eight-wall spherical supports. Their Schur bounds allow a `(3,3)` leaf either no foreign contact or one of six foreign π/3 contacts. Every other profile allows none. All cross-B2 pairs are forced right. For each of the 121 profile pairs at cuts 0 and 3, the two contacts from cut 0 to the final leaf pair are bounded by label 3 using exact actual-vertex failures at label 4. In those tests every unused cut contact is relaxed to right, so no monotonicity is assumed across a newly right contact.

The resulting **8,092 cases** exhaust the relevant projections and two finite contacts. **8,088** have a nonzero residual three-by-three determinant, forcing full rank at least eleven. The remaining four contain an affine A3 four-cycle in an actual finite vertex and are also impossible. The two unused cut profiles cannot repair a principal-submatrix contradiction. Dotted distances and the two central normals are never assigned numerical values.

The [primary certificate](fresh_h9_star_cube_metric_20261014_all_cut_profiles.json) and [checker](../scripts/fresh_h9_star_cube_metric_20261014_all_cut_profiles.py) use exact arithmetic over Q(√2,√5). Generation and saved-result replay pass. A [separate Fraction arithmetic review](fresh_h9_star_cube_metric_20261014_all_cut_profiles_review.md) independently reconstructs all 8,092 determinants, the 242 upper-contact bounds, and the four actual affine closures without importing the producer arithmetic. Earlier frozen stages include the [first-assignment certificate](fresh_h9_star_cube_metric_20261014.json), [second assignment](fresh_h9_star_cube_metric_20261014_second.json), [right/label-5 assignment](fresh_h9_star_cube_metric_20261014_right_cut.json), and [nine-profile monotone gate](fresh_h9_star_cube_metric_20261014_cut_profiles.json). The latter has a [separate Fraction arithmetic review](fresh_h9_star_cube_metric_20261014_cut_profiles_review.md).

A [new continuation](fresh_h9_star_cube_metric_20261014_all_cut_search_case_00.json) applies this gate only when all its fixed-label premises hold. It reaches the next complete mandatory assignment in 13 nodes: positive leaf 11 has own contacts `(3,4)`, with all other first-assignment labels unchanged. This is an unresolved metric candidate in that saved search; it is not a Coxeter realization. A separate orientation extension is being investigated without modifying these frozen certificates.

Replay:

```sh
python scripts/fresh_h9_star_cube_metric_20261014_all_cut_profiles.py --verify
python scripts/fresh_h9_star_cube_metric_20261014_all_cut_search.py --case 0 --limit 100000 --verify
```
