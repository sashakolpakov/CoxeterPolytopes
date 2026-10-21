# Independent review of the single-facet wedge obstruction

The [local theorem](fresh_dimension_wedge_20261013.md) is sound for its stated ansatz: retain every other old wall in the old hyperplane and replace one old facet by two spacelike normals with positive proportional old projections and opposite nonzero added coordinates. Unequal slopes are included.

For a finite nonright old neighbor, the three relevant normals span a positive-definite three-space. Their Gram determinant rules out three nonright Coxeter angles, so the new pair must be orthogonal. The resulting cosine-square equation permits only replacement labels (3,3), (3,4), or (3,5), corresponding to old labels 4, 6, or 10. This is an unrestricted-label argument. An old label six forces a squared slope 1/3, which gives a forbidden angle against a retained parallel neighbor.

The exact source checker passes. A separate stdlib implementation also reconstructs the three Gram matrices and checks every saved obstructing neighbor: **11 H8 facets, 36 H19 facets, and 210 H21 facets**, 257 in total. It uses integer/Fraction operations and an independent two-coefficient golden-field implementation, without NumPy or SymPy.

The side observations have the stated scope. After deleting a facet, projections in the retained positive cone leave its spacelike normal feasible; perturbation toward an actual strict interior point handles even zero projections. An old cusp not incident with the split facet likewise admits a small added coordinate and a strict spacelike perturbation. Neither statement excludes a general dimension extension with new projection directions or several moved old facets.

- [Independent certificate](fresh_dimension_wedge_20261013_independent.json)
- [Independent checker](../scripts/fresh_dimension_wedge_20261013_independent.py)

```sh
python scripts/fresh_dimension_wedge_20261013.py --verify
python scripts/fresh_dimension_wedge_20261013_independent.py --verify
```
