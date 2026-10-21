# Independent review of the pentagonal-radius covering obstruction

The source and mathematical argument pass review. The theorem excludes finite-volume completion of the retained weighted-box 39-wall core by a finite coherent selection from all seed-compatible radius-two caps, all seed-compatible radius-`4/phi` caps, radius-one caps at literal half-vertices, and arbitrary-radius caps at literal vertices. It does not cover arbitrary radius-one centers, other nonvertex radii, or a changed core.

Reviewed files: [producer](../scripts/fresh_h20_weighted_pentagonal_20261014.py), [certificate](fresh_h20_weighted_pentagonal_20261014.json), and [report](fresh_h20_weighted_pentagonal_20261014.md). The separate [review checker](../scripts/fresh_h20_weighted_pentagonal_20261014_review.py) and [review record](fresh_h20_weighted_pentagonal_20261014_review.json) use stdlib rational arithmetic in the basis `1,sqrt(5)`, with no producer or project arithmetic imports. Independent generation and saved replay pass. The primary saved replay also passes.

The independent checker covers the new pair classes, distance-two restriction, four-cap motif, strict motif escape, and counting constants. The previous 88-cutter radius-two classification, its 7,744-pair adjacency census, and the radius-four bounds of 512 caps and 275 probes per cap are explicit dependencies; this is not a second implementation of those earlier enumerations. File hashes bind the producer, certificate, review checker, and prior certificates. The review prose can acquire links without changing these arithmetic inputs.

## All-label completeness and new arithmetic

For radius `R=4/phi`, the smaller opposite-wall gap is at most the coordinate inradius. Since `R/2>1`, a unit coordinate has only a zero small gap. For heavy and golden coordinates, `R/sqrt(2)>phi` leaves only zero or `R/2`. The heavy `R/2` choice has complementary normalized gap strictly between `1/2` and `sqrt(2)/2`, as independently verified with rational radical enclosures. Thus it cannot be any finite Coxeter angle. The remaining coordinate catalog is complete for every label, including parallel and dotted contacts.

The central relation admits exactly two full golden coordinates at label five, or at least three with a dotted relation. The half-probe power is negative only at two or three full golden coordinates. An opposite coordinate sign makes it positive. This yields the 84 cutters per probe used in the proof. Radius-one half-center caps cut only their own half-probe strictly; a one-unit sign change gives tangency instead.

Every new normalized cap-cap cosine lies in `Q(sqrt(5))`. The finite cosine values in this field are precisely `0,1/2,phi/2`; the cyclotomic degree test excludes the other labels. The checker reduces the 7,056 adjacent `R5/R5` pairs to 42 coordinate-permutation classes and the 7,392 adjacent `R5/R2` pairs to 24 classes. Separate arithmetic rejects every class. Here `R5` abbreviates `4/phi`.

The 84 adjacent `R1/R5` pairs give 28 forbidden, 21 label-five and 35 right pairs. A compatible `R5` neighbor necessarily has three full golden coordinates. The 88 adjacent `R1/R2` pairs give 11 dotted and 77 forbidden pairs, with compatibility requiring the changed golden bit to be the `R2` cutter's full coordinate. Adjacent golden `R1/R1` pairs are forbidden.

The checker independently evaluates all 3,136 distance-two pairs from their nineteen spatial coordinates, rather than reusing the producer's compact distance formula. Exactly 240 are right and 2,896 forbidden. Compatibility requires disjoint full-coordinate triples, with exactly one changed coordinate in their union. Hence three compatible `R5` neighbors of one `R1` probe would require three disjoint triples in an eight-element set, which is impossible.

The four-cap square motif has four right, one label-five and one dotted relation. Its explicit norm-one escape also passes the separate arithmetic. The source's strict-center and acute-projection argument correctly establishes actual facets; no finite-vertex conclusion is inferred from this motif.

## Graph counting and the boundary perturbation

In each fixed unit/heavy-sign fiber, let `A` contain every golden-cube probe strictly cut by at least one selected `R1`, `R2` or `R5` cap, and choose one cutter at each probe. This definition gives a bound on total coverage, without presupposing failure of the vertex-radius family.

Every induced edge either is charged to the uniquely designated full coordinate of an `R2` endpoint, or joins `R1` to `R5` and is charged to its `R1` endpoint. The former contributes at most one charge per `R2` probe. The distance-two restriction limits the latter to at most two per `R1` probe. Therefore `e(A)<=2*n_R1+n_R2<=2|A|`.

The degree-eight cube also gives `e(A)>=8|A|-1024`. Thus `|A|<=170`; 171 already violates the two inequalities. Across 2,048 fibers, these three families cover at most 348,160 probes strictly. Arbitrary vertex radii add at most `512*275=140800`, because every non-radius-four vertex cap misses all probes strictly. The union leaves at least `524288-348160-140800=35328` probes not strictly cut. Overlaps only decrease coverage.

The source handles boundary contacts rather than equating non-strict coverage failure with a strict escape. A radius-two power cannot vanish: after its golden coefficient is set to zero, its rational coefficient is strictly positive. A radius-`4/phi` power has an odd rational coefficient and cannot vanish. The earlier vertex-radius certificate likewise excludes boundary equality.

The sole remaining boundary case is a radius-one half-center differing from the missed probe in exactly one unit sign. Moving all ten unit coordinates by `epsilon` times their current signed half-values increases every such tangent ball power by `epsilon+(5/2)*epsilon^2`. The common displacement is therefore outward for all boundary balls at once. For a finite selection, sufficiently small positive epsilon preserves every other strict inequality and stays inside the vertical box. The resulting uncovered open boundary patch, or a further small spacelike lift, proves infinite volume. Finiteness of the selected wall set is used here and is stated in the theorem.

No mathematical issue was found. The numerical scout is not needed for this obstruction.

```sh
python3 scripts/fresh_h20_weighted_pentagonal_20261014_review.py --verify
```
