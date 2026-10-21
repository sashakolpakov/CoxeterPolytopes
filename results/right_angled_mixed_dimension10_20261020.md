# Current theorem: finite-volume right-angled dimension at most ten

20 October 2026. **New unconditional internal theorem, with reviewed
written proof, exact scalar generation/replay and conditional Lean count
logic.** There are no arithmetic restrictions. External refereeing and
end-to-end formal hyperbolic geometry are not claimed. The consolidated
`manuscripts/right_angled_dimension10.tex` and its12-page PDF are complete
in the same mixed-paper family: three vector figures, zero build warnings,
two full manuscript reviews and final figure-page visual PASS.

Every finite-volume right-angled hyperbolic polytope has dimension at
most10. In particular this applies to polytopes with both finite and
ideal vertices. Known mixed examples reach dimension8, leaving a gap2.
The published unrestricted right-angled ceiling is12. This supersedes
the repository's earlier ceiling11; its source, PDF and certificates
remain preserved as predecessor evidence.

## The proof in ordinary counts

For a hypothetical mixed eleven-polytope let f>0 count finite vertices,
B count actual canonical P3 three-faces and
E=sum_non-P3(facets+cusps-10)>=0. The complete balanced incidence
argument gives

    3722B>159195f+2794E.

Actual P3 oppositions preserve h=3q-sum n_i, with h^2=2. A fixed
future null direction perpendicular to h on each component gives eleven
real frame coordinates with sum3 and squared sum1. The new structured
positive-part estimate is proved for every size0,...,11 of an ideal-
ending direction subset. It yields, by actual-subset deletion and
undirected reciprocal pairing,

    4B<=168f+3E.

The comparison would force5742f+5E<0, impossible. The noncanonical
all-ideal-edge corner allowance H<=2E used here has a direct elementary
three-face proof; metric D3 recognition is unnecessary.

Any higher-dimensional finite-volume right-angled polytope has an
eleven-dimensional face of the same kind. The published compact and
entirely ideal ceilings4 and6 make that face mixed, which was just
excluded. Thus the new proof needs no import of the earlier internal
11 theorem or Dufour's old12 tail.

## Written proof and review

The full assembly source is
`right_angled_mixed_dimension10_assembly_candidate_20261020.md`, SHA256
`70ef7de40e9fbd62dbe763d64d258c40ab777a42d7550e1608074e4c29aa0326`.
Its historical pending-check header is superseded by this current
record. The independent full assembly review is
`review_right_angled_mixed_dimension10_assembly_finite57_20261020.md`;
it includes complete proofs of the direct H budget and hereditary tail.

Independent analytic review of all subset sizes:
`review_mixed_eleven_all_subset_contracts_resource_20261020.md`, SHA256
`3b0bc185cf728c22307ff0e9cd8a3ee2cb1844a9cad0e83eb205bcf2164f95e0`.
The new lower count and one-cusp geometry were independently reviewed in
`review_mixed_eleven_one_cusp_and_corner_balance_root_20261020.md`.
Root separately reviewed the full ten-direction argument and checker in
`review_mixed_eleven_k10_envelopes_root_20261020.md`, and the global
nine-direction proof in `review_mixed_eleven_global_nine_root_20261020.md`.

## Exact and formal verification

Three current scalar workers have each passed generation and separate
read-only replay under the shared resource lock:

* `check_mixed_eleven_balanced_and_subset_rows_20261020.py`:
  balanced count arithmetic and all seven/eight/eleven-direction rows;
* `check_mixed_eleven_nine_direction_global_20261020.py`:
  global nine-direction polynomial identity and positive margin635/968;
* `check_mixed_k10_rational_envelopes_20261020.py`:
 840 closed Bernstein intervals and12550 strictly positive coefficients.

The earlier complete nine-direction interval proof has its own successful
generation/replay and remains separate supporting evidence. No optimizer,
polytope enumeration or sampled-positivity inference is used.

`formal/CoxeterProofs/RightAngledMixed10.lean`, SHA256
`41186cd79e50167a227b55bc145a3c86e06a420aeacc78fc74a166e3331a4794`,
compiled in a one-thread,2048MiB-capped batch after the written and exact
checks. Its separate printed-axiom audit passed for eight theorems.
Only propext and Quot.sound occur; the exact comparison is axiom-free,
and there is no sorryAx. It imports only Std, with no legacy dependency
gap. Compilation took1.47 seconds and audit1.21 seconds.

The kernel derives the balanced lower, cusp/defect eliminations, upper
count consequence and contradiction from explicit primitive count
identities, geometric incidence bounds and an aggregate analytic
transport premise. Actual-wall geometry, the real scalar contracts and
the published compact/ideal/82-facet inputs remain written premises;
they are not hidden geometric axioms reported as formalized results.

The reproducibility/assumption map is
`docs/right_angled_dimension10_proof_map_20261020.md`.
The final source/PDF hashes are respectively
`2fd439454de989f34590681f11e0475c137effb9faa71a72011107bf5b0922c1` and
`9b1362cb77405084b83fc6798bb9a4b951d453e32c1eda1d1f632018195dc205`.
All nine disposable figure previews were removed after review. The
preserved ceiling11 source/PDF are predecessors, not another current
primary paper. The catalogue ID for the current single-command package
is `mixed10`.
The actual single-command core run
`results/mixed10_core_verification_20261020.json` passes all six jobs:
four exact replays, the conditional kernel/axiom audit and a fresh
temporary PDF build. The declared final PDF is preserved unchanged.
