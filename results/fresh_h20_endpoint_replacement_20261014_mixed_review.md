# Independent review of the mixed endpoint construction

**PASS.** The frozen mixed construction has 1,202 actual Coxeter facets and the stated strict spacelike escape. The new caps cut the preceding square-centered escape and create finite contacts, but no finite vertex involves a mixed cap. No finite-volume or complete vertex/cusp census claim follows.

The [separate checker](../scripts/fresh_h20_endpoint_replacement_20261014_mixed_review.py) and [saved audit](fresh_h20_endpoint_replacement_20261014_mixed_review.json) bind the producer source, certificate, and inherited 810-wall certificate by SHA-256. They import no producer code. This review uses FLINT rational matrices, independently rounded rational radical intervals, and one symbolic polynomial identity. The producer's source and [report](fresh_h20_endpoint_replacement_20261014_mixed.md) were also reviewed mathematically.

The replay independently reconstructs the rational spatial Gram matrix and all product vertices by linear solves, checks its positive definiteness, and rebuilds every literal mixed center. It checks all 8,624 mixed–seed relations and all 76,636 mixed–mixed pairs. Complete cyclic Gram identities justify transferring the 392 independently checked mixed–endpoint reference pairs to all 153,664 pairs. The resulting counts agree with the certificate: 784 mixed `pi/3` pairs and 392 mixed–endpoint `pi/8` pairs, with the remaining pairs dotted.

For the endpoint radius, write `c=cos(pi/8)` and `rho/R=sqrt(c²+17/7)-c`. Substitution proves

`(rho/R)²+2c(rho/R)=17/7`.

This is an exact algebraic contact identity; overlapping intervals are not used to infer equality. Independent rational intervals prove the strict separation inequalities. The inherited 395 middle/motif balls lie strictly inside the radius-eight ball, while both new orbits lie outside it. Their mutual old pair certificate is inherited, not rerun here.

The common timelike point `H=(0,200,1)` is strict. For every wall normal `r`, the projection `H-(H,r)r/r²` remains timelike and is strict on all other walls, because all cross pairings are nonpositive. Thus every wall is an actual facet. The 392 points `P_z=(z,(196+rho²)/2,1)` have norm `-rho²`; the replay checks the nineteen independent seed incidences at each, while the orbit separation proves every other inequality strict.

The finite-neighbor argument is complete for each mixed cap. It has seventeen right seed neighbors, four mixed `pi/3` neighbors and one endpoint `pi/8` neighbor. Those last five are pairwise dotted. A finite clique containing the cap therefore has at most nineteen walls, fewer than the twenty independent normals required at a finite H20 vertex. Parallel seed partners cannot occur at a finite point. This excludes new finite vertices involving these caps without assuming a numerical vertex census.

The checker also directly tests both representative parallel-cusp incidences against both complete new orbits. Sixteen independent perpendicular seed directions together with the parallel pair give quotient rank seventeen. All 784 transported tangent points are distinct; their links remain deficient. This verifies the listed cusps, not every ideal boundary direction.

Finally, the replay checks the previous target's negative ball power `2270077/5000-176sqrt(7)`, and independently verifies the new norm-one escape against all 1,202 walls and the continuous vertex-radial relaxation. The remaining escape is strict and future-oriented, so infinite volume follows. The matched three-factor shift has positive pairing `R²/8`, correctly requiring replacement rather than an additional compatible cap.

Generation and saved-result replay both pass:

```sh
python scripts/fresh_h20_endpoint_replacement_20261014_mixed_review.py --verify
```

No source, formula, count, or scope issue was found. This audit does not independently repeat the inherited middle/motif pair computation, enumerate every vertex or cusp, or establish arithmeticity or novelty.
