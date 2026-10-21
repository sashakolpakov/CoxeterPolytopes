# Record 31, `a=1`, source orbit 5: exact incidence no-go

## Result

Source orbit 5, whose positive missing four-face is

\[
\{2,11,12,13\}\qquad(\text{mask }14340),
\]

has **no** completion satisfying the corrected compact `a=1` ridge
condition.  Two terminal exact quotient traversals reach the same
theorem-facing census:

| traversal | search nodes | pure 120-facet targets | pure ridge signatures | corrected ridge targets |
|---|---:|---:|---:|---:|
| production nauty quotient, fixed-width inline key | 41,488,368 | 735 | 231 | **0** |
| independent direct binary-matrix quotient | 6,042,612 | 735 | 231 | **0** |

The complete scheduled-blocker-size histogram and the complete
`pure_target_ridge_histogram` agree exactly between the two runs.  The
required Dehn--Sommerville ridge signature is `{2:540}`; it occurs zero times
in both transcripts.  Therefore this source dies at the incidence gate.  No
topology, Coxeter-label, Gram, or Vinberg endpoint computation is needed for
this source.

This is a source-orbit theorem, not the full Record-31 `a=1` theorem.
Production sources 2 and 4 and the remaining independent-incidence rows must
still be completed.  It makes no assertion about finite-volume noncompact
chambers and constructs no new hyperbolic Coxeter polytope.

## Exact terminal data

The production traversal has depth-node partition

```text
[1, 196, 518, 8456, 165396, 3923619, 37353399, 36783]
```

and quotient-seen partition

```text
[1, 3, 51, 1099, 41169, 1713317, 25378215, 8231].
```

It records 14,346,282 exact isomorphic prunes, 5,081,805 uncoverable deaths,
21,487,754 slot-bound deaths, 3,282 no-small-transversal leaves, and 735 pure
120-facet targets.

The matrix quotient has depth-node partition

```text
[1, 196, 644, 5068, 60477, 1032822, 4905830, 37574]
```

and quotient-seen partition

```text
[1, 4, 38, 477, 11759, 412438, 2934227, 7995].
```

It records 2,675,673 exact isomorphic prunes, 878,415 uncoverable deaths,
2,381,155 slot-bound deaths, 3,283 no-small-transversal leaves, and the same
735 pure targets.  Counts before the pure-target gate need not agree because
the two state equivalences choose different intermediate representatives.
The identical terminal scheduled-blocker and ridge histograms are the
relevant cross-check.

## Frozen certificates

- Production transcript:
  `data/enumerations/h9_14_record31_a1_orbit5_corrected_quotient.json`,
  SHA-256
  `90fb08a5889d6c3fa0dd30aae35cb33feee334695482a31a29513a64ad92d294`.
- Independent transcript:
  `data/enumerations/h9_14_record31_a1_orbit5_independent_matrix_quotient.json`,
  SHA-256
  `7dee139c3b5740635df2f7814e8e34f58983a9a6469f40fa7b183d586d9c395b`.
- Verifier:
  `scripts/verify_h9_14_record31_a1_orbit5_transcripts.py`.
- Regression test:
  `tests/test_h9_14_record31_a1_orbit5_transcripts.py`.
- The production inline key is covered by the exact dense/compact/inline
  bounded equivalence certificate
  `data/enumerations/h9_14_record31_a1_inline_key_equivalence.json`.

Run the frozen verifier with

```bash
PYTHONPATH=src:. python3 scripts/verify_h9_14_record31_a1_orbit5_transcripts.py
```

The verifier pins both transcript byte hashes and the generated production
source, checks every node/prune identity and terminal histogram total, and
requires exact equality of the two terminal pure-target censuses.
