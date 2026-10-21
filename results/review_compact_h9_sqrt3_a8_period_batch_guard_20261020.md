# Independent read of the conditional A8 period

20 October 2026. **PASS within the stated conditional scope.** This is a
complete paper/source read, not execution, matrix import, a new census,
Lean verification or an independent arithmetic implementation.

Reviewed sources:

* `compact_h9_sqrt3_a8_period_isometry_20261020.md`, SHA256
  `f4d80265de3c167c82e741af40d45cb62957d63ff020211e24b7616f4229137b`.
* `../scripts/check_compact_h9_sqrt3_a8_period_isometry_20261020.py`, SHA256
  `f8f3d06aa69848779ed2739ccc0561e66de63e3e5e500944970ee15f85895d30`.

I independently rederived the two A8 frames, their timelike orthogonal
vectors, every coordinate of the index-nine glue vector, the two-equation
solve for the images of v and beta, and the trace 5+3sqrt(3). No algebraic
flaw was found. The following checks explain the source coverage.

1. The original free basis is (t,z1,...,z8), with
   z9=(4+sqrt(3))t-sum(z1,...,z8). Expanding its form gives exactly the
   displayed block Gram matrix, not the ambient ten-coordinate diagonal
   matrix. The embedding columns are seven original differences, the
   displayed beta root and v. Their dimensions and first-eight coordinate
   truncations are correct. The determinant +/-1 check, together with the
   exact transformed Gram identity, certifies the full integral basis.
2. The worker independently inverts the E8 Cartan matrix and builds
   fundamental weights from its columns. Its frames are row vectors for
   the Gram tests and are transposed into columns for constructing the
   isometry. Thus `right * inverse(left)` maps each old simple root to the
   reversed new chain and maps V0 to V1; there is no frame-transpose error.
3. The exact field operations, rational Gaussian elimination and physical
   sign rule are correct. The worker separately checks the integral
   inverse, determinant, both isometry identities, frame images, future
   cone sign, and the explicit glue vector. These are stronger checks than
   merely comparing an authored matrix with itself. Its original-basis
   conjugation uses the embedding in the correct direction.
4. In the paper calculation the seven known diagonal entries sum to
   -20-12sqrt(3). Solving the remaining two image equations gives diagonal
   entries 69+40sqrt(3) and -44-25sqrt(3). Their total is
   5+3sqrt(3)>9, so the real rank-nine isometry has infinite order.
5. The already generated A7 certificate is hash-pinned, checked for its
   E8 simple-wall result and exact endpoint, and its nested source map is
   checked recursively. Streaming pre/post hashes cover the worker and
   all its dependencies. The certificate explicitly leaves the new alpha2
   interface pending and makes no nonreflectivity claim.

There is one essential wording qualification to the paper's opening
conditional statement and Section4: “first endpoint” must mean that the
displayed r5 is the **first simple wall after resolving all simultaneous
mirrors**, with the displayed A8 germ. Equality of the scalar collision
parameter alone does not suffice: a simultaneous E8 overgroup could cut
that A8 sector. This qualification is explicit in the worker's scope and
in the separately prepared full positive-cone alpha2 census. The parent
and author have recorded it without changing the pinned historical note.

With that stronger interface, both A8 germs are in the same actual full
reflection chamber. The integral isometry permutes the full mirror
arrangement and sends the inward germ to the inward germ, so it preserves
the chamber. Its infinite order contradicts compactness. The subsequent
nonreflectivity conclusion still uses the inherited anisotropy and the
full scalar-facet transfer; neither is established by this finite matrix
checker. Until the first-simple-wall interface is certified, this remains
an exact infinite-order lattice isometry and a conditional stopping test.
