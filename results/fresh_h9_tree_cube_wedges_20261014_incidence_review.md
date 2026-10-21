# Non-path cube-tree incidence review

Source/proof review and independent literal incidence replay **PASS**. The [checker](../scripts/fresh_h9_tree_cube_wedges_20261014_incidence_review.py) imports no producer or graph-canonicalization code. Its [certificate](fresh_h9_tree_cube_wedges_20261014_incidence_review.json) binds the frozen geometry source and data.

All 240 recipes are reconstructed: all 24 cut orders on each of the star and T trees, with each of five negative facets left unwedged. Each edge stellar subdivision replaces precisely the simplicial facets containing that edge. The independent reconstruction then applies the four ordinary-wedge vertex formulas, recomputes the base minimal nonfaces, and checks every stored bijection on the complete resulting vertex set and minimal-nonface clutter. In total, 105,184 complete vertex incidences are mapped. All recipe multiplicities and canonical-order relabelings agree.

The rational realization argument is sound. Each cut depth is half the minimum positive vertex slack outside the specified ridge. It therefore cuts exactly the intended vertices and crosses only their outgoing edges; the asserted nonzero slacks rule out a cut through an old vertex. The ordinary wedge inequalities preserve boundedness and give the stated nine-dimensional vertices. Exact active-rank, inequality, and facet-rank checks in the producer establish simplicity and all eighteen actual facets. The connected missing-face hypergraph establishes nonproduct incidence.

This review does not repeat the rational-coordinate arithmetic for all 33 representatives. It verifies every proposed identification without `pynauty`, but the assertion that no further isomorphisms exist among the 33 representatives still uses the producer's canonical-labeling dependency. It does not audit a Coxeter-label or metric search, or exclude all other cube-tree constructions.

```sh
python scripts/fresh_h9_tree_cube_wedges_20261014_incidence_review.py --verify
```
