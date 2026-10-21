# Pointed simplex vertices: no additional restriction on the selected base domain

Refining simplex vertices from their full ambient finite forest to their
automorphism orbit in the complete simplex star does **not** strengthen the
current selected-provenance model. Every new signed capacity column is
identical to, or implied coefficient by coefficient by, an existing
forest-labelled column. This is a diagnosis of this domain, not a claim that
pointed labels never help.

The exact bounds remain

| Saved dual | Exact bound | Approximation |
| --- | --- | --- |
| Full optimized model | 2790968068837079040 / 216216020508130303 | 12.9082390023 |
| Exported-only restriction | 129790338384790720 / 10289212336025873 | 12.6142151747 |

Both exceed the rank-13 target 12. No simplex or vertex mass changes, and no
optimizer is invoked. These are obstructions to the specified partial weight
model; they are not ambient realizations or dimension bounds.

## Invariant labels and literal base maps

The [weighted capacity theorem](fresh_vertex_weighted_prism_base_capacity_20261017.md)
allows any invariant label of a vertex in its containing maximal simplex S.
For an orbit O of such vertices it gives the signed column

\[
 \sum_\rho m(\rho,S,O)\alpha_\rho
 -(13-d)|O|\alpha_S\le0.
\]

Here every prism provenance uses its actual expanded row in the
[frozen full model](fresh_sqrt5_labelled_packing_full_model_20261017.md).
The positive coefficients of these new columns are deliberately zero on
unselected merged quad rows, as in the preceding forest-labelled model. No
actual split or omitted base data is assumed for those unselected rows.

The complete normalized simplex star is an orthogonal union of one connected
Lorentz core and finite positive-definite components. A weighted graph
automorphism preserves the unique Lorentz component. Conversely, any core
automorphism extends by the identity on the spectator components. All simplex
vertices are spherical one-root deletions of the core. Consequently their
complete-star orbits are exactly the root orbits computed here. Spectator
permutations cannot merge distinct core-root orbits.

The computation retains the entire edge label `(m,a,b)`, including exact dotted
values. It performs complete pointed graph-isomorphism searches for all 501
pairs of spherical-deletion roots in the 287 relevant cores and saves the
successful bijections. The graph search checks edges and nonedges; absence of
an automorphism is determined by exhaustive backtracking, not a refinement
color alone.

The [saved forest-labelled base maps](fresh_sqrt5_vertex_refined_base_capacity_20261017.md)
are checked against the literal core isomorphisms for all 1,858 actual
provenances, 6,016 base occurrences, and 14,284 base vertices. A prism's vertices
are the saved union of its genuine quadrilateral `K2,2` vertices. Its containing
simplex uses **all** spherical one-deletions, including additional vertices
outside the base. These sets are not interchanged. Every mapped base vertex
belongs to an actual simplex orbit and has the same full ambient finite forest.

## The only split and why it is redundant

Across 4,969 containing-simplex profiles, the 11,397 forest-labelled columns
become 11,412 pointed columns. Exactly 11,382 are unchanged. The other 30 come
from 15 spectator contexts of core 281, the weighted path

\[
 0\mathrel{\stackrel{5}{-}}4
 \mathrel{\stackrel{3}{-}}1
 \mathrel{\stackrel{5}{-}}2
 \mathrel{\stackrel{3}{-}}3.
\]

Its spherical deletions are 1, 2, and 4. Deleting 1 gives `A2+I2(5)`; deleting
either 2 or 4 gives `A1+H3`. The latter roots are not related by a weighted
automorphism, so they form two distinct pointed labels even though their
finite forests agree.

For each of the 15 contexts, and for **every** selected actual provenance,
the coefficient at root 2 and the coefficient at root 4 are each at most the
coefficient at root 1. All three singleton capacities equal
`13−2=11`, with the same negative simplex row. The root-1 label was already
separate in the old model. Thus its signed inequality implies both new ones.
The JSON saves each implication and every positive coefficient difference;
this includes zero-mass provenances and proves redundancy for any nonnegative
mass distribution on the selected actual domain.

Only two of those contexts have positive incoming mass in the exported-only
witness. Their existing root-1 inequalities already supply exactly the needed
capacity (one has a rounding slack of `1/10^10`). The complete dual checks do
not rely on this observation: they recount all new columns and every inherited
positive-cost and signed column for both saved duals.

## Reproduction and status

- [Exact source](../scripts/fresh_sqrt5_pointed_simplex_base_capacity_20261017.py)
- [Saved certificate, orbit maps, column implications, and both duals](fresh_sqrt5_pointed_simplex_base_capacity_20261017.json)
- [Prior support diagnosis](fresh_sqrt5_labelled_packing_dual_support_20261017.md)

Generate with `python3 scripts/fresh_sqrt5_pointed_simplex_base_capacity_20261017.py`.
Replay every orbit search, literal map, coefficient implication, and exact
column check with the same command followed by `--verify`.

**FINAL/FROZEN.** Generation and a separate fresh saved replay both passed,
including complete model reconstruction and both exact dual recounts.
Direct source hashes, Python syntax, and local report links passed. The
frozen predecessors remain unchanged. These are internal computational
checks, not external refereeing or proof-assistant formalization.
