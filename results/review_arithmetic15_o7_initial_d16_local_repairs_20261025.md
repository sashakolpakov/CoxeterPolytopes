# Independent review of the arithmetic O7 / initial-D16 local repairs

25 October 2026. **Accepted at exactly the five-repair delta boundary.**
This does not accept the full arithmetic manuscript, newly discharge any
external or finite-output premise, or supply a build or runtime result.

I am the author of the preceding fresh independent audit,
`results/audit_arithmetic15_o7_initial_d16_translation_fresh_20261025.md`
(SHA-256 `3dc8e8080c93e0c0e65532c74cc47aa15c7592d3b54356c63066728a5df986ec`).
Root authored these manuscript repairs. I did not edit the reviewed TeX.

## Exact source boundary

The preserved preimage is
`results/source_archive/arithmetic15_before_o7_initial_d16_audit_repairs_20261025/compact_arithmetic_dimension_bound_consolidated.tex.gz`.
Its decompressed source is the previously reviewed 131541-byte edition,
SHA-256 `122f3b98627290a95ac9323b85cb9ae95b22376af768f6f6dfe64539d776edef`.

The successor is
`manuscripts/revisions/compact_arithmetic_dimension_bound_consolidated.tex`,
132670 bytes, SHA-256
`648e0425afb72f1c9f4a21c814a7253d6189b55762a9778430d45c847c3ced46`.
The complete 3664-byte patch is
`results/arithmetic15_o7_initial_d16_audit_repairs_20261025.patch`, SHA-256
`a58ead96a5b78930b487c46ee3376b02a8b4efc37830a2e30461438df10a493d`.
Its same-stem JSON receipt correctly keeps outside review pending.

I read the entire patch, independently compared the preserved preimage
with the live successor, reread the affected statements and checked the
signed duals against the original R proof's explicit root numbering.
The actual complete diff contains exactly the five listed local repairs.

## Five checks

1. The future-unit-point displacement now explicitly belongs to
   `p^perp`. This supplies the missing spatial hypothesis of the pressure
   inequality and agrees with the applications examined in the full audit.

2. The local dual notation retains BOTH genuine signed forks, with
   `Omega_15=(P15-e16)/sqrt(2)` and
   `Omega_16=(P15+e16)/sqrt(2)`; these signs are correctly numbered.
   The original R proof uses `n15=(e15-e16)/sqrt(2)` and
   `n16=(e15+e16)/sqrt(2)`. Direct pairing gives
   `Omega_15.n15=Omega_16.n16=1` and the crossed pairings zero.
   For `i<=14`, `Omega_i=sqrt(2)P_i` is the corresponding unit-root
   dual. Defining every `h_i=x.Omega_i` leaves the original vector x
   fixed and does not identify the two fork directions.

3. The fork table is now expressly a table of UNIT time squares.
   White-three uses raw norm2 and raw time square `2A^2`; isolated
   fork uses raw norm `2U_f` and time square `2U_f(2U_f-1)`;
   white-four uses raw norm4 and time square `4A^2`.
   The added K-rescaling justification is correct: at a native
   pi/4 contact with a norm2 old root alpha, `(r,alpha)^2=(r,r)`;
   multiplying r by `2/|(r,alpha)|` gives norm4, negative white
   product -2 and own product `-a=-2sqrt(2)gamma`.
   This does not assume sqrt(2) belongs to K.

4. Replacing signed norm1 by ABSOLUTE norm1 is the required quartic
   correction. The unit conclusion is unchanged; the other real
   conjugates make the signed norm -1 in the case under discussion.

5. The far-next-three own-three branch is now the actual FIFTEEN-node
   simple-three cycle `c,n2,n3,...,n14,n_f,c`, rather than a triangle.
   There are c, thirteen spine nodes and one genuine fork. Its unit
   Gram rows have diagonal1 and two off-diagonal entries -1/2, so
   the positive all-ones vector is a kernel and their sum is null.
   That sum has c's strictly positive time. It is orthogonal to the
   selected walls; its pairing with EVERY other original wall is
   nonpositive by all-wall acuteness. It therefore gives a nonzero
   feasible future ideal ray of the original compact chamber.
   The corrected compactness contradiction is valid and uses no
   dimension conclusion or face-arithmeticity assumption.

## Retained scope

The new O7 propagation and initial-D16 inventory retain the conditional
acceptance stated in the full audit, now without its five local wording
and hypothesis defects. The remaining A-component, row, companion and
field/census transcription obligations in that report are unchanged.
The historical source ledger is not silently rebound to these successor
bytes. No canonical edition was promoted by this review.

This review used source reads, byte hashes and a source diff only.
No scientific worker was imported or executed; no enumeration, test,
Lean, TeX, Sphinx or other proof job was launched.
