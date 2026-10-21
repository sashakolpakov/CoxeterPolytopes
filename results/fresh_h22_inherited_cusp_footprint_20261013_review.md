# Review of the 860-facet inherited-cusp footprint

Source/proof review and a fresh saved replay of [the checker](../scripts/fresh_h22_inherited_cusp_footprint_20261013.py) **passed**, with no mathematical issue found. This review uses its separate sparse square-root arithmetic and the already independently certified complete old corner walk; it does not claim a second implementation of that arithmetic or a new corner enumeration.

The [exact output](fresh_h22_inherited_cusp_footprint_20261013.json) distinguishes 303 bounded old cusp rays: all 280 rays of the A5^4 orbit, 21 of the 42 D4^5 rays, and two of the D12+E8 rays. The source independently reconstructs their null projections, checks distinct primitive orientations, canonicalizes their exact active sets into the complete 18-orbit census, and verifies each opposite-height parallel pair. The old quotient rank20 plus the new orthogonal affine A1 proves bounded links of rank21; this is a positive-kernel argument, not rank alone. The other 22 positive-height types have strictly timelike old projections.

The lower-bound argument is sound for a finite acute-compatible completion retaining the 860 facets. Every added old projection lies in the retained old chamber cone, hence is future causal. A nonzero such projection can be orthogonal to an old null ray only when positively proportional to that ray. Pure vertical projections are incompatible with the retained caps of both height signs. At each of the remaining old cusps, incident caps of both height signs are necessary to bound the added-coordinate recession line. Finitely many nonincident walls disappear sufficiently deep in a horoball and cannot supply the missing constraints. Each nonzero null projection serves only one old null ray.

Therefore the **1,095,687 remaining old cusp rays require at least 2,191,374 additional facets**, if this core is retained. This is a scoped lower bound, not a nonexistence theorem or a bound on constructions with changed facets. The separate [finite-contact graph obstruction](fresh_h22_covering_cells_20261013.md) explains why merely supplying many distant cusp caps is also insufficient.

```sh
python scripts/fresh_h22_inherited_cusp_footprint_20261013.py --verify
```

Fresh replay: PASS, approximately 23 seconds.
