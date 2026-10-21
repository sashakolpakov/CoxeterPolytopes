# Joint high-baseline ADE vertex set in dimension13

There are276 simply-laced spherical forests of rank13. Exactly215 have the old baseline payment at least12, the strict h-vector target in this dimension. This package enumerates the complete finite one-wall transitions within that joint set. Different forest types at the two endpoints are allowed.

The result is116,386 accepted literal attachments among835,705 tested contact subsets. There are712,393 nonspherical choices and6,926 spherical choices outside the selected high-baseline set. The type adjacency graph connects all215 selected forests. This illustrates why separate same-type frequency bounds leave a substantial gap: the competing types admit many finite local transitions.

The selection uses the established simply-laced baseline, one unit for root pairs at graph distance one or two and zero otherwise. Every component and forest total is independently recounted in the saved output. The full forest census includes the61 lower-baseline forests, and every selected forest retains its exact margin over12.

For every selected forest and every deleted root, the checker tests all contact subsets of the twelve retained roots of sizes zero through three. These are complete because the destination is an ADE forest: the new root has degree at most three, all nonright contacts have label three, and it may be an isolated A1. The checker requires the whole destination forest to lie in the joint set. It records the complete contacts and destination index for each accepted attachment. The finite classifier checks that each connected component is a tree and uses the exact A/D/E arm-length classification.

The [source](../scripts/fresh_sqrt5_rank13_high_ade_endpoints_20261017.py) and [saved data](fresh_sqrt5_rank13_high_ade_endpoints_20261017.json) support generation and complete saved replay:

```sh
python scripts/fresh_sqrt5_rank13_high_ade_endpoints_20261017.py
python scripts/fresh_sqrt5_rank13_high_ade_endpoints_20261017.py --verify
```

These are finite endpoint possibilities, not arithmetic edge realizations. No opposite wall product has been chosen. A triangle with all three vertices in the joint set is impossible in a compact arithmetic Q(sqrt5) polytope: all entries of its normalized full-star Gram are rational, whereas its negative direction cannot survive in the positive-definite nonphysical embedding. This argument does not require the three forests to agree.

Two-special triangles and all-special quadrilaterals still need their full joint classification. Non-ADE vertices and the61 lower-baseline ADE vertices remain outside the set. There is no joint incidence inequality, dimension bound or simultaneous compact realization in this package.
