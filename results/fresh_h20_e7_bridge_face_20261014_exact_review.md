# Independent review of the 96-wall bridge-face cutter

**PASS.** The separate [checker](../scripts/fresh_h20_e7_bridge_face_20261014_exact_review.py) and [saved certificate](fresh_h20_e7_bridge_face_20261014_exact_review.json) independently reconstruct the new wall and its exact geometry from the frozen [construction data](fresh_h20_e7_bridge_face_20261014_exact.json). Generation and saved replay pass. The [construction report](fresh_h20_e7_bridge_face_20261014.md) correctly states that the old escape is cut but a new strict escape proves infinite volume.

This implementation imports only the independent [whole-face review's sparse Fraction arithmetic](../scripts/fresh_h20_bridge_faces_20261014_review.py), and imports no construction arithmetic. It multiplies sparse radicals in the 16-element multiquadratic basis and two successive positive quadratic extensions. Signs are certified by outward rational enclosures, rather than the producer's recursive squaring tests. Source and data hashes are bound; reports are not, allowing review links to be added safely.

## Arithmetic and actual geometry

The independent replay checks:

- All 96 normal squares and the ambient inertia `(20,1)`.
- The rational PD18 support, its entrywise nonnegative inverse, the particular-vector anchor equalities, orthogonality of the complement, and the exact norm identity for `r=p+sqrt(K)z`.
- All 95 new pairings: 18 right, one π/3, one π/6, and 75 strictly dotted. The nonright finite contacts are exactly old walls 0 and 28.
- A common strict future timelike center, and an explicit relative interior point of the new facet with all 95 old inequalities strict.
- The actual rank-19 edge: its support is the PD18 block plus the orthogonal norm-two new root, so it is genuinely independent and positive definite. Its saved future timelike point has exactly those 19 active walls and all other inequalities strict.
- The positive new-wall pairing of the old escape, and all 96 strict negative pairings of the new rational spacelike escape.
- Strict negative new-wall pairings at all three previously bounded cusps.

The new escape square is exactly

`20181174579467302741 / 10^18`.

The checker accepts its literal rational coordinates directly. It does not depend on the numerical rounding or the cusp-edge perturbation used to discover them.

## Proof review and inherited scope

The [whole-face review](fresh_h20_bridge_faces_20261014_review.md) independently establishes that the old support face is the entire `H2`. The new perpendicular cap supplies its first geodesic boundary. The rank-19 edge is an actual edge, not a claimed finite endpoint or a finite-area polygonal completion.

The [acute-extension preservation lemma](acute_extension_preservation_20261014.md) is sound under its stated exact hypotheses. Projecting a strict center onto any listed wall gives an actual facet point. At an old finite vertex, the positive-definite incident Gram has entrywise nonnegative inverse; decomposition of the new strict center then forces every acute new wall to be strictly negative at that vertex. At a proper bounded full cusp, full-support positive affine relations and span of the null perpendicular force the same strict retention for a spacelike new wall. Thus all 96 listed facets are actual, every previously certified finite vertex remains, and the three old bounded links remain locally unchanged.

The lower bound of 116 old finite vertices and the old bounded-link certificates are inherited from the bound predecessor and cusp walk; this review does not repeat those censuses. It also does not repeat all old-old pairings. The new norm, pair, point and escape arithmetic is independently reconstructed.

No new finite vertex, complete new edge census, global compactness or finite volume is claimed. No numerical mate-search result is promoted to an all-parameter obstruction. The renewed exact strict escape settles the volume status of this particular 96-wall chamber.

## Reproduction

```
python3 scripts/fresh_h20_e7_bridge_face_20261014_exact_review.py --verify
```

No issue was found in the construction formulas, sign orientation, actual-edge argument, preservation lemma or renewed infinite-volume witness.
