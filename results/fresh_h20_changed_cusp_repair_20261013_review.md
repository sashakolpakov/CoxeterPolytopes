# Review of the 213-wall H20 cusp repair

Source and mathematical review found no issue in the [213-wall checker](../scripts/fresh_h20_changed_cusp_repair_20261013.py) or [certificate](../data/enumerations/fresh_h20_changed_cusp_repair_20261013.json). Its saved-result replay passes. This review reuses the construction's exact arithmetic; it is not a second implementation of all field operations.

Each of the 17 specified deficient cusp links gains the two norm-three walls `5δ+z` and `5δ-z`, with `z²=3`. Their sum is `10δ`, their mutual pairing is `-3`, and they are orthogonal to the old incident walls. The checker verifies every new pair as it inserts the 34 walls, all 204 inherited cusp signs, the strict center, and the surviving strict spacelike escape.

The bounded-link conclusion uses more than rank: the old positive kernel weights cover every old incident wall, are strictly positive, and combine with the new pair to give a positive multiple of δ. The actual normal span has rank 20, so its Euclidean quotient has rank 19. Every recession vector has nonpositive pairing with each wall; the positive relation forces each pairing to vanish, and the spanning condition then forces the recession vector to be zero. A strict interior point supplies the nonempty link. The source now explicitly checks that the incident coordinates are rational before applying its rational rank calculation; I separately checked this premise on all 17 links.

The all-height statement is also sound within its stated scope. Every spacelike normal in one of these planes is a positive multiple of `tδ+z` or `tδ-z`; the pure δ direction is null. A retained old wall gives a necessary lower bound `t≥A>0`, while nonnegative pairing with the saved escape gives `t≤B`. In 32 cases the bounds are disjoint. In the two remaining cases, a retained short wall has cosine strictly between zero and one half throughout the closed interval `[A,B]`, excluding right, finite Coxeter, parallel, and dotted relations. Hence every allowable cap from those 17 real planes is strictly negative at the same escape. Any **finite** collection of such additions leaves a strict spacelike witness and cannot give finite volume. This does not exclude other directions or wall replacements.

Replay:

```sh
python scripts/fresh_h20_changed_cusp_repair_20261013.py --verify
```

The review covers the 213-wall repair and its continuous plane-family obstruction, not the subsequent gap-cutting caps or a complete cusp census.
