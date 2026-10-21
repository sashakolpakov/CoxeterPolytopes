# Record 31, `a=1`, source orbit 0: independent incidence replay

The missing independent-incidence row for source orbit 0 is now complete.
The production nauty quotient searched 231,773 nodes; the fresh direct
binary-matrix quotient searched 242,124 nodes.  Both give 491 pure 120-facet
targets, the same complete 42-signature ridge histogram, and exactly one
`{2:540}` target orbit.

The two traversals retain different literal representatives, as expected.
Their complete `C2 x S9` canonical family key is identical:

```text
[0,0,12,12,16,16,5,5,14,17,17,42,42,50,50].
```

Direct minimal-transversal replay on each representative independently gives
120 facets and ridge histogram `{2:540}`.  The previously frozen topology and
unrestricted compact-label audits then kill this unique incidence orbit, so
source 0 still has zero Gram candidates.

Artifacts:

- production transcript:
  `data/enumerations/h9_14_record31_a1_orbit0_corrected.json`, SHA-256
  `6f771cbaf9b1b2249217c0e25aa3c4a5217b63a6f7abd1f1328155e700bfaa39`;
- independent transcript:
  `data/enumerations/h9_14_record31_a1_orbit0_independent_matrix_quotient.json`,
  SHA-256
  `f39712bdadeee4ca1a7852f52efa10ccaa75ac92b9409949ec2285fd84384325`;
- verifier:
  `scripts/verify_h9_14_record31_a1_orbit0_independent_transcript.py`;
- regression test:
  `tests/test_h9_14_record31_a1_orbit0_independent_transcript.py`.

This closes one independent-incidence gap only.  It is not a full `a=1`
schedule theorem and says nothing about finite-volume noncompact chambers.
