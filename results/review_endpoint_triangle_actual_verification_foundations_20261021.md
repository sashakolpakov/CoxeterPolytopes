# Endpoint triangle: actual verification and scope review

21 October2026. **PASS for the recorded final compilation, separate actual
axiom audit, source-review chain and documentation status.** This review
performed only file reads, streamed hashes, byte sizes and comparisons.
No compiler, test, worker import, mathematical replay or source/map/receipt
edit was performed. It does not repeat the written geometric reviews.

## Immutable receipt and map

* Manifest `results/lean_endpoint_triangle_verification_20261021/manifest.json`:
  SHA `f37513b4701ca1583a0c890ba224a13dce00f6845b487ca37f618b6128948c4b`,
  5813 bytes; left unchanged.
* Current map `docs/endpoint_triangle_norm_kernel_logic_20261021.md`:
  SHA `d77fb0817e837d8a90a7fa5445c99b89e68893d9aaf537a0a94404374f47fdb6`.
* Exact pending-map predecessor at
  `results/source_archive/endpoint_triangle_norm_logic_map_before_execution_20261021.md.txt`:
  SHA `c9beffd99c3d4b963573dad65fcfb28637b4d85de215c187a2e6b6886fef4a75`.

I authored the independent static source review and two precision-delta
reviews, not the kernel, actual execution or receipt. Their role and
earlier shared geometric foundations remain disclosed in those reviews.
Here I independently checked the resulting execution bindings and scope.

## Actual source, object and run evidence

The current6715-byte source is unchanged at
`ce2ad70b636880d68e0e99f345f82b8db3c7bd7197636ea9c93c0d23d98df1a4`.
Both actual ledgers bind that source, the SAME compiler/environment
fingerprints, empty LOCAL dependency-object map,2048MiB and one thread.
The three toolchain/project configuration hashes and all four named
compiler binary/library hashes were streamed and match. This is not a
separately hashed transitive closure of every Std object.

| Evidence | SHA256 prefix | Bytes | Recorded outcome |
|---|---|---:|---|
| Archived compile ledger |d9e320fa1578|1524|KERNEL_COMPILED, .88s|
| Archived separate audit ledger |9b21501438c1|1443|KERNEL_AXIOM_AUDIT, .77s|
| Compile outer stdout |1819339edd4e|1720|terminal PASS_SEQUENTIAL_MODULE_CHECK|
| Audit outer stdout |07d23171f801|3981|terminal PASS_SEQUENTIAL_MODULE_CHECK|
| Actual child axiom output |3f6974bd724a|2384|ten declaration outputs|

All manifest-bound file sizes/hashes match. The current operational
ledgers in `formal/.proof-batches/` match their archived counterparts
byte-for-byte. The preserved351088-byte object at
`formal/.lake/build/lib/lean/CoxeterProofs/EndpointTriangleNormKernel.olean`
has SHA `1336f5cf7dd3c16e88e30a52972ed5371304bb35c67d218aa480f4024b695df4`,
matching the compile ledger. No compiler was invoked to inspect it.

The actual2384-byte child output hash is
`3f6974bd724ad8fb47c6e5ed1ced5191d55b4b2f7705cd70130cf855185f49ae`.
It matches BOTH successful child-output fingerprints and the first2384
bytes of the separate audit outer stdout. The compile runner prints its
ledger rather than duplicating that child text; its stdout is therefore
correctly a different size. The audit text includes the child output
followed by its ledger and terminal PASS.

Recorded initial resources are134.63GiB free disk for both runs and
6.28/6.57GiB available RAM, respectively. They exceed40/6 admission;
they are NOT measured peak RSS. The two actual child durations are not
an aggregate whole-paper runtime. The owned temporary work directory
and operational last-failure file were absent at this review.

## Ten actual declarations, not merely requested commands

I read the actual child text and compared its declaration names with
every current `#print axioms` command and the manifest list. Exactly
these ten occur, once each:

    cycle_expansion
    bezout_triangle_identity
    physical_I_cap
    conjugate_cap
    physical_extended_cap
    product_cap
    norm_small_factor_impossible
    expanding_pair_contracts
    norm_expanding_factor_impossible
    rational_positive_bounded_product_not_multiple

The first nine list exactly propext, Classical.choice and Quot.sound;
the last lists propext and Quot.sound. No sorryAx or additional axiom
occurs. The two retained warnings concern unused ordered section
variables in cycle_expansion and bezout_triangle_identity, matching the
receipt. Compilation and audit are separate evidence, but repeated
prints of these declarations are not counted as twenty distinct results.

## Preserved failures and complete review chain

Both failure-preservation records and their bound predecessor/failure
identities match the earlier reviewed records. The e187 source failed
its original broad grind; the92c source failed because `scaled` was a
redundant explicit local grind argument. Error-recovery sorryAx outputs
from either failed batch retain ZERO verification credit. They are not
substituted for the final actual output above.

The manifest binds the full static reviewd8318de0 and explicit-cancellation
reviewd92a00ae. It does not directly include the final one-line review;
per root's instruction the sealed manifest remains unchanged, and THIS
review closes that provenance link explicitly:

`results/review_endpoint_triangle_local_parameter_precision_delta_foundations_20261021.md`,
SHA `d72c0abe320e0c96cfff23053ad85e9380defe99ba44efabf6cdfdf791bf76d2`,
2540 bytes. It was streamed and matches. The full chain binds the exact
final ce2 source with unchanged statements; the delta reviews are static
logic evidence, while only the final successful runs supply actual
verification credit.

## Map delta and exact scope

The complete map diff changes the pending status to actual compile/audit
PASS, records two warnings/two failed predecessors, links the separate
written ordinary-A1 acceptance and adds the preserved-map provenance.
The algebraic description and all external assumptions are unchanged.
The linked written acceptance explicitly binds ordinary extension17d9
and its separate full review, distinct from the short-I source917c.

The final map correctly limits Lean to conditional polynomial cancellation,
Bezout witness, ordered inequalities and finite-product norm logic.
Actual walls, cyclic integrality, conjugate Schur bounds, number-field
embeddings/norm identification, same-vector reference payments, complete
finite-white geometry and the final ambient argument remain application
inputs. Standard axiom output does not enumerate those parameters or
turn this into end-to-end geometric Lean. No global ceiling, new example,
whole-paper runtime PASS or old-module rerun is asserted by this receipt.

No remaining metadata or scope correction requested.
