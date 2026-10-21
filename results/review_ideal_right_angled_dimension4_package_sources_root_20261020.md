# Root review of the ideal dimension-four package sources

20 October 2026. Complete source review PASS. Execution and final
figure inspection are recorded separately, not inferred from this review.

## Sources reviewed in full

* Authoritative proof note: `ideal_right_angled_dimension4_proof_20261020.md`,
  SHA256 `9e6f7eac5438c3a71c8c90bc5c615dcdac51072dfb7869d9654ef3c76e9d8768`.
* Final manuscript source: `manuscripts/ideal_right_angled_dimension4.tex`,
  SHA256 `52296f3297f1e5b3696250ec729ba3524658494d8a442ee2ea2e1d0f956fb9b5`.
* Conditional logic: `formal/CoxeterProofs/IdealRightAngled4.lean`,
  SHA256 `3aad406da960989e2fe7a1e7a5e60c700243ba011d9427e333571ae32ebc4fad`.
* Exact finite model worker:
  `scripts/check_ideal_right_angled_dimension4_model_20261020.py`,
  SHA256 `bb506ec1ee98dee1429d6e668553af40316ebdeda6b59b82bbf03b8f3a6dc261`.

The complete universal geometric proof is the same argument reviewed
in `review_ideal_right_angled_maximum_height_kissing_root_20261020.md`.
The current note and paper correct Kolpakov's journal theorem number
to4 and retain all actual-wall, zero-height, metric-equality and
connectedness arguments. They use no internal small-facet or chain
classification. The four-page format omits no mathematical dependency.

## Manuscript and figure

I read the full311-line first source, then made and inspected the exact
three final figure-label changes: shorten the left explanation, remove
the overlapping path label, and lower both panel explanations slightly.
These edits change no equation, theorem, assumption or geometric drawing.
The first source, PDF and build record are byte-verified gzip archives
listed in `ideal4_dimension4_first_layout_review_20261020.json`.

The figure's radii1 and2 and center separation sqrt5 are exact:
5=1+4 gives orthogonality, and the distance formula at q=(0,1)
gives sinh-distance (5+1-4)/4=1/2. The graph panel asserts only
combinatorial propagation. The proof separately establishes that
the graph edges are actual ridges; cusp-only contacts are not edges.
Final raster layout inspection remains a separate step.

## Conditional Lean boundary

The module uses an arbitrary height type. It does not replace real
heights with integer or rational heights. Finite actual facets and
neighbor lists are supplied, with exact incidence and no duplicates.
The supplied positive maximum and zero-height cusp facet are
geometric data; their existence is not silently formalized.

`PublishedInputs` separately supplies Musin's cardinality, Kolpakov's
minimum, and his metric equality conclusion. `WrittenGeometry`
separately supplies the actual projection, normalized bounds, spherical
code reduction, ideal-face inheritance, opposite-pair center argument,
and neighboring height equality. None assumes global propagation,
nonexistence of five-polytopes or the final dimension bound.

The kernel derives the24-cardinality squeeze and regular metric at a
maximal facet, closes the maximum-height predicate under adjacency,
propagates it by induction along an actual finite path, and contradicts
the height-zero cusp facet. Five-face inheritance is explicit in the
final dimension implication. Eight actual printed-axiom commands
support a later audit; this source review is not a compilation claim.

## Exact worker boundary

The worker's displayed model is the standard D4 unit-normal model.
The24 normals,276 pair products and12 spatially opposite pairs are
checked exactly. The24 displayed null vertices have six incident walls
with three opposite parallel pairs, giving cubical links. Each wall
has a strict actual facet witness (2,r), and the scalar code-bound
identity is a polynomial identity with integer coefficients.

These are meaningful finite model checks. They do not prove complete
enumeration of vertices, classify arbitrary polytopes, prove the
published external results, or formalize the universal Lorentz
reduction. Those limits are explicit in both source and output.
