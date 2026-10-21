# Root review of the changed H20 assembly

The complete frozen checker passed a separate root replay after dependency
hardening. It binds 40 prior reconstruction inputs, checks the changed core
against all 102 retained old normals, checks every sequential new polar pair,
and verifies the strict center, target cuts and final strict escape.

The source review found the even-wall doubling and lift reconstruction
consistent. All finite angles at the doubled H19 wall are even; all pairs
in the resulting lifted core are checked again directly. The retained old
indices exclude the removed facet 175. Every new normal has positive norm,
every pair is nonpositive and Coxeter-compatible, and the same timelike
center is strict for all 179 normals. Orthogonal projection of that center
onto any wall remains timelike and satisfies all other inequalities strictly,
so every inequality is an actual facet.

The [separate Fraction checker](../scripts/fresh_h20_changed_assembly_20261013_root_review.py)
imports no producer arithmetic. It independently reconstructs the 49-wall
double, its connected bipartite long-wall graph, and all 71 lifted normals.
It checks their 2,485 pairs and common center, all six rational polars against
the core and previous polars, their sequential strict target cuts, and the
remaining escape against these 77 rational walls. It also contracts both
rational and `sqrt(6)` coefficients of the previous escape, confirming the
new core's exact cut `12171/4096`.

The remaining escape square is exactly
`10640902782393308787/20000000000000000000 > 0`. Together with the full primary
replay, this proves infinite volume. The 102 retained radical normals and
their cross pairs remain an explicit shared dependency, rather than being
claimed as separately reconstructed by the Fraction checker.

The [saved independent result](fresh_h20_changed_assembly_20261013_root_review.json)
records the pair counts, exact inequalities, scope and source hashes.

```sh
python3 scripts/fresh_h20_changed_assembly_20261013_exact.py --verify
python3 scripts/fresh_h20_changed_assembly_20261013_root_review.py --verify
```
