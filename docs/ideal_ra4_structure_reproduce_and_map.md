# Second-volume equality: reproduction and lemma map

For current reading, see the [current second-volume manuscript](../arxiv/ideal_four_dimensional_volume_and_structure/main.tex) and the
[submission guide](arxiv_submission_editions_20261028.md). The supporting
arguments and recorded verification scopes below remain distinct from those
current submission builds. Standalone draft typesetting routes have been retired.

The supporting argument translates the frozen [universal equality proof](../results/ideal_ra4_second_volume_equality_20261018.md): every entirely ideal right-angled hyperbolic4-polytope of volume $8\pi^2/3$ is the reflected double of the ideal24-cell. This strengthens the earlier second-distinct-value theorem. It is internally checked; external refereeing and end-to-end Lean formalization are not claimed.

## Reproduction

From the repository root:

```sh
python scripts/fresh_ideal_ra4_equality_structure_20261018.py --verify
python scripts/check_ideal_ra4_equality_structure_saved_20261018.py
```

Python uses only its standard library. Default generation refuses to overwrite an existing certificate; replay uses `--verify`. The helper imports no producer and uses separate reverse-column elimination. Final generation, separate-process replay and helper passed in terminal run `bb139d`. The manuscript compiled to seven pages, session16424, final chunk `34a1f3`. Typesetting is not mathematical verification.

The finite checks cover the Figure7 incidence transcription, D9 local cusp/quad facts, all364 active triples of the cuboctahedral model and all12 width ratios, all15 labeled simple4-regular six-vertex graphs, and cube square-equation rank4. They do not reprove the historical classification or universal geometry. The [audit](../results/ideal_ra4_second_volume_equality_20261018_audit.json) is integrity/status only.

|Frozen artifact|SHA256|
|---|---|
|[Producer](../scripts/fresh_ideal_ra4_equality_structure_20261018.py)|`298b3722b6aeceecab918f3c137882bf4f672cf94fb7fe646a6ac724d353ef5f`|
|[Finite supports](../results/fresh_ideal_ra4_equality_structure_20261018.json)|`589a7901494481ddbefb8b3cf2a51854c70bd804bd78128386bd650dcf4e19fc`|
|[Independent helper](../scripts/check_ideal_ra4_equality_structure_saved_20261018.py)|`464fa7261976233865b170bb747d9fa499e498860d1b337d2fc8ae3bfd096b00`|
|[Universal proof](../results/ideal_ra4_second_volume_equality_20261018.md)|`4d227e9f7f4577f0a61b30a3ebf2790c763feeb4f1acddf7799a83f96fc1d7e2`|

## Published and inherited inputs

The [second-value map](ideal_ra4_second_volume_lemma_map.md) records the published volume formula, antiprism bound/equality, nine-vertex classification, right-angled3-dimensional rigidity and24-cell minimum. Equality additionally uses the **complete ten-vertex octahedrite entry** in Kolpakov2012 §4.3 and Figure7: its four quadrilaterals have edge-adjacency graph two disjoint edges. Three reviewers independently checked the primary text and image, without claiming a kernel-verified historical census.

The frozen gap proof supplies the internal opposite-cusp injection and antiprism/D9 shape facts. They are not relabeled as published assumptions. Acute Coxeter wall products are nonpositive. Regular24-cell flag transitivity is used; the new proof explicitly derives surjectivity onto its octahedral facet symmetry group.

## Internal equality interfaces

These are proof obligations, not proposed external axioms.

1. **AcuteOrthogonalIntersection.** Projecting a strict future interior point off mutually orthogonal walls remains future timelike and strict on every other wall. Three pairwise adjacent facets have an actual common edge; four give a forbidden finite vertex.
2. **NoAntiprismAtDefect24.** The opposite-cusp injection excludes A5 and larger. Nonsquare-box propagation forces at least three extra triangular positive neighbors of A4, excluding it too.
3. **MinimumDefectThree.** Budget24, the minimum positive defect, q≥3 outside antiprisms and the ten-vertex quad matching exclude m≥5 and m4. No4-polytope census is assumed.
4. **ExactlyEightD9.** The N7 excess partitions and absence of positive triangular neighbors at D9 eliminate defect4, leaving the cuboctahedral alternative. Its twelve nonsquare cusps contradict the opposite-cusp injection. Embedded medial recovery is proved explicitly. All eight remaining positive facets have defect3.
5. **NonsquareBoxType.** All face aspects are1 or2. D9's long-side identification excludes2:2:1; each noncubic box is2:1:1 with four D9 belt facets and two octahedral ends. There are six such cusps.
6. **SixCuspSimpleGraph.** Quadrilaterals join their two nonsquare cusps. Parallel edges would repeat an asymptotic opposite-facet pair at two cusps. The simple4-regular six-vertex graph is octahedral; eight distinct D9 triangles give its dual cube. No unproved planarity or cubic-eight-vertex classification is used.
7. **EqualWidthNormalSums.** Opposite unit normals sum to a positive multiple of their cusp null vector. Equal short widths give equal sums and all six cube square equations.
8. **AffineCubeNormals.** Square equations give $n_s=a+\sum s_i b_i$; unit and right-angle equations force mutually orthogonal mode squares $(-1/2,1/2,1/2,1/2)$. The finite rank check supports, but does not replace, the universal derivation.
9. **ActualOctahedralSection.** The Lorentz4-space cuts out a regular ideal octahedron S with six actual cusps, so convexity gives $P\cap H=S$. Its face interiors and edge interiors lie inside actual D9 facets and quadrilateral ridges respectively.
10. **NoOctahedralCrossing.** An octahedral facet meets S only at ideal points; crossing H would give a forbidden finite intersection. Cutting therefore introduces only octahedral halves and no finite vertices.
11. **TwoRegularHalves.** Both halves have W0 and are24-cells by the published theorem. Flag transitivity extends every octahedral facet symmetry, fixing the gluing placement up to reflection. This proves uniqueness, not merely numerical volume addition.

The root-owned Lean gap module currently formalizes conditional integer consequences for the earlier value theorem. It does not yet formalize these eleven equality interfaces. Future formalization must keep the actual-cusp, Lorentz-normal and convexity obligations explicit, while proving finite graph/cube algebra and defect partitions internally.

## Review scope

Root completed full written-proof, producer and helper source mathematical reviews. Order16 and construction each completed full written-proof and primary-premise reviews, without claiming full producer/helper reads. None executed the certificate or reconstructed the historical classification during review. The author's generation/replay/helper runs are separate evidence.
