# Independent geometric review of the restricted ideal RA5 block theorem

18 October 2026. **PASS, within the stated C/D/T facet-spectrum hypothesis.**
The full [proof report](ideal_ra5_octahedral_blocks_20261018.md) and its
[producer](../scripts/fresh_ideal_ra5_octahedral_blocks_20261018.py) were read.
The local data were replayed, and a separate
[standalone checker](../scripts/check_ideal_ra5_octahedral_blocks_saved_20261018.py)
reconstructed all 104 polygon faces without importing the producer.

The following geometric steps were checked independently.

* Each 24-vertex ridge has exactly two hexagonal two-faces; those are exactly
  its intersections with the other 24-vertex ridges meeting it. An isometry
  gluing the shared ridge preserves polygon vertex counts. At a finite point
  inside such a two-face, exactly three ambient right-angled facets contain
  it. Thus the correspondence identifies the actual third facet and forces
  the neighboring 24-ridge adjacencies. It is stronger than an unlabelled
  polygon-count comparison.
* Native normals inside a right-angled facet equal the ambient neighboring
  facet normals. Their Gram matrix therefore applies unchanged. The opposite
  products -1 exclude full ridge adjacency, not only 24-ridge adjacency.
  The induced neighborhoods are exactly four-cycles. The six-vertex closure
  proof for a four-regular locally four-cycle graph is valid.
* Each opposite-pair sum in an octahedral block is nonzero: a strict interior
  timelike vector pairs negatively with it. The sums are null and mutually
  orthogonal, hence lie on one future null ray. The full six-normal span has
  signature (3,0,1). Its orthogonal complement has one positive direction
  and the same null radical. Two unit normals with the same sign along that
  positive direction have product 1, so acuteness permits at most two common
  neighboring facets.
* Five pairwise orthogonal facet normals would give an actual finite vertex
  by the written projection of a strict interior vector. The correction has
  the correct sign at every other wall. This uses the entirely ideal
  hypothesis and proves that the full T-facet adjacency graph is K5-free.
* Minimum degree 22, the six-block partition and the common-neighbor bound
  rule out 30 vertices after applying the explicit published four-color
  theorem. The three-color and four-color block cases were checked directly;
  no abstract coloring is assumed geometrically realizable.

The published input was checked in the original
[Andrasfai--Erdos--Sos paper, Theorem 1.1, p.207](https://real.mtak.hu/110445/1/1974-18.pdf).
Its r=5 specialization is precisely that K5-free graphs with minimum degree
strictly greater than 8n/11 are four-colorable. The inequality at n=30 is
strict: 22 > 240/11. The published graph theorem is an external premise;
its proof was not reproduced in Lean or independently re-proved here.

The standalone checker independently examines all 4*197 pairs consisting of
one large ridge wall and another T wall. Having at least three common ideal
vertices is equivalent to orthogonality for these pairs. Every retained
vertex set has rational span rank three and hence gives the expected
two-dimensional polygon. The counts are twelve triangles, twelve
quadrilaterals and two hexagons per ridge. All records agree with the saved
certificate. Completeness of the actual T polytope and its 198 walls / 258
vertices remains inherited from the frozen doubling proof.

The producer replay and initial independent reconstruction completed in one
sequential batch, terminal chunk `650b8c`, exit 0. After adding verification of
both inherited source pins and the certificate's producer source hash, the
final independent checker passed in terminal chunk `87bd0c`, exit 0:

```sh
python3 -B scripts/fresh_ideal_ra5_octahedral_blocks_20261018.py --verify
python3 -B scripts/check_ideal_ra5_octahedral_blocks_saved_20261018.py
```

The producer replay also verifies the saved finite color-domain arithmetic.
This review did not independently enumerate every color tuple, re-prove the
original C/D/T construction or old cusp count equations, or review the new
TeX manuscript. No Lean or external human review is claimed. The conclusion
is **at least 36 T facets and 84 D facets under the C/D/T restriction**;
there is no new unconditional dimension ceiling, spectrum nonexistence
theorem, or geometric example.

Reviewed producer SHA-256:
`3de8dce876893c0988dd33b86cc7e7c3c9ec94b3f339ba6c421a68b78edf7edc`.
Reviewed certificate SHA-256:
`3d22bb8e1f8ef6951796b3bf32412d928dfb68de79f704a1bc9bc0f192f8d91d`.
Standalone checker SHA-256:
`56f244e17d319d25cb8492c37fb91416d38a89a22aef5d6af8443dcfa8912e2b`.
