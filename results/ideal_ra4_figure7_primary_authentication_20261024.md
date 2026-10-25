# Ideal RA4: primary Figure7 authentication

24 October 2026. Author: `audit_sharp_ra`. This closes only the
external-image transcription boundary left open in Section 10 of
`results/audit_ideal_ra4_fresh_20261024.md`. No mathematical worker,
enumeration, test, Lean compilation or TeX build was run. The accepted
manuscript and the historical producer/checker were preserved.

## Sources and exact scope

The publisher's original article is A. Kolpakov, *On the optimality of
the ideal right-angled 24-cell*, Algebraic & Geometric Topology 12 (2012),
1941–1960, DOI [10.2140/agt.2012.12.1941](https://doi.org/10.2140/agt.2012.12.1941).
I read Case D across printed pp. 1951–1952 in the
[publisher PDF](https://msp.org/agt/2012/12-4/agt-v12-n4-p03-s.pdf).
Its complete-list passage concerns **ideal right-angled three-dimensional
octahedrites with at most ten vertices**, having only triangular and
quadrilateral faces. It identifies Figures6–7 together with the
six-vertex octahedron as the admissible list. Figure7 is its ten-vertex
entry. This is not a classification of all ideal RA3 ten-vertex
polyhedra without the face-size hypothesis: the pentagonal antiprism is
outside that hypothesis and is treated separately in our manuscript.

For a direct image check I used the author's deposited
[arXiv v2 article](https://arxiv.org/html/1211.2944v2), dated 14 November 2012.
Its Figure 7 links directly to
[octahedrite10.png](https://arxiv.org/html/1211.2944v2/octahedrite10.png).
This distinguishes the **author-deposited source image** actually viewed
from the **published completeness text** read separately. I do not claim
that the unavailable web PDF screenshot was visually inspected.

The bounded source-image download is only at
`/tmp/coxeter_ideal_ra4_figure7_primary_20261024.png`:

- 18,390 bytes; PNG 778×780, RGBA, non-interlaced.
- SHA256 `189fb05dd2006effbd59ebef3b46a9fb7896dfb6f421080b7eff99b4c0dcf46f`.
- Retrieved from the exact author URL above and viewed directly with
  `view_image`; no raster modification or local PDF rendering was used.

The image is **not copied into the repository or release payload**.
The retrieval URL/hash and the independently written incidence data are
the reproducibility record; no redistribution licence is inferred.

## Explicit labeling and incidence comparison

The source image has no A–J labels. Assign them as follows, consistently
with `ten_vertex_incidence()` in the historical producer:

|Label|Position in the source image|Producer's schematic coordinates|
|---|---|---|
|A|outer top left|(0,0)|
|B|outer top right|(6,0)|
|C|outer bottom right|(6,6)|
|D|outer bottom left|(0,6)|
|E|upper inner left|(2,2)|
|F|upper inner right|(4,2)|
|G|middle right|(5,3)|
|H|middle center|(3,3)|
|I|lower inner left|(2,4)|
|J|middle left|(1,3)|

The four quadrilateral face boundaries, up to reversal/cyclic order, are

```
ABCD    ABFE    EHIJ    CGHI
```

The outer region is the face ABCD in this planar drawing. The eight
other regions are triangles:

```
AEJ    AJD    DJI    DIC    BFG    BGC    EFH    FGH
```

Reading the drawn edges gives

```
AB BC CD DA AE AJ BF BG CG CI DI DJ EF FG GH HF HE EJ JI IH
```

Every labeled vertex has degree four. The two quadrilateral adjacencies
are precisely `ABCD--ABFE` along AB and `EHIJ--CGHI` along HI. The
other pairs have no common edge. Sharing only C or E does not create
edge adjacency. Thus their quadrilateral edge-adjacency graph is exactly
two disjoint edges, as used in the manuscript's defect-four exclusion.

These manually read boundaries match the producer's stored vertex/edge
definitions and asserted quadrilateral cycles at lines55–86, and the
standalone checker's boundary predicate at lines42–61. I inspected that
source comparison but did not execute either program. The source image
identification is now directly checked rather than inferred from a
successful incidence checker alone.

Root independently viewed this exact local image and confirmed all four
quadrilaterals, all eight triangles and the AB/HI matching. Root also read
the original publisher CaseD complete-list passage directly. That is a
second-author visual/text review, not a rerun of the historical
classification or scripts.

## Release-premise update for the existing ideal-RA4 family

Keep two separate entries in the existing assumption/proof register:

1. **External classification:** under the ideal RA3/triangles-and-quads/
   at-most-ten-vertices hypotheses, the published CaseD list is complete;
   its ten-vertex octahedrite is the Figure7 type. The historical
   classification proof remains an external input.
2. **Authenticated incidence transcription:** the author-v2 image with
   the exact URL/hash above yields the stated A–J cycles and matching.
   This identification has now received direct visual checks by two
   authors. The existing finite checker verifies consequences of these
   data; it does not prove the classification's completeness.

This changes the fresh audit's Figure7 status from “stored incidence
checked; primary image unavailable” to “primary author-image directly
checked and independently confirmed.” It does not change any dimension,
volume, decomposition claim, or Lean verification status. Dufour's
original-source access question belongs to the sharp-RA family and is
not closed by this deliverable.

## Bound repository sources

|Path|SHA256|
|---|---|
|`manuscripts/ideal_ra4_volume_and_structure_publication.tex`|`e0c2ca24cec24ffcd736e201d3c7e06d75b28ef7be8093decd321ae3ac54bda1`|
|`scripts/fresh_ideal_ra4_equality_structure_20261018.py`|`298b3722b6aeceecab918f3c137882bf4f672cf94fb7fe646a6ac724d353ef5f`|
|`scripts/check_ideal_ra4_equality_structure_saved_20261018.py`|`464fa7261976233865b170bb747d9fa499e498860d1b337d2fc8ae3bfd096b00`|
|`results/audit_ideal_ra4_fresh_20261024.md`|`38bb7f21b0910a2dcf75d8444ce2e6bdcf223af32cf945120cba320063373640`|
