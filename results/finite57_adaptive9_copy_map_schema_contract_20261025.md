# Adaptive9 copy-map schema repair

Source-only successor; no imports, syntax checks, tests, numerical enumeration,
builds or Git operations were performed in preparing this change. The real
failed control execution is retained at
`results/finite57_adaptive20_structural_execution_20261025.stdout`; it is not
a PASS receipt.

## Exact sources

- Unchanged loader: `scripts/check_finite57_current_adaptive9_path_binding_20261025.py`,
  SHA256 `76022d95fcc5ec6e61bf1ccaea65635adda1ab8593d0bf921d75754f14bd139c`.
  The frozen universal29dd and quiet-case9e4e preparation chains retain it.
- Successor: `scripts/check_finite57_current_adaptive9_copy_map_safe_20261025.py`,
  368 lines, 19,507 bytes, SHA256
  `a62d18107a54c512f7abcb1ceb9b02dd187eed38805928b5f0d0bf096659a0e2`.
- Existing control filename:
  `scripts/test_finite57_current_adaptive9_twenty_controls_20261025.py`,
  221 lines, 11,712 bytes, SHA256
  `b9870e7fba167ce615108e3aff33a006ef942e2efd3c4815f0d8c4dfa25fdb64`.
  Its 205-line/10,779-byte preimage, SHA256
  `1d128ef3c6bfccb0921a934b02f29126df84905e07de1399b3decfb0b2064dc9`,
  was copied and its bytes hashed before editing under
  `results/source_archive/finite57_before_adaptive_copy_map_schema_repair_20261025/`.

## Representation equivalence

The saved producer input has SHA256
`c738c6af3f910e74b24c17a10c3f38cab2ed3835f4a6661c745046abb98ec8ab`.
Its producer uses `R.combined`, then accesses `mp[i]` for the old right node
`i`. In `scripts/fresh_selective_distance10_dimension116_20261017_review.py`
lines135--146, `unite` returns `[mapping[j] for j in range(len(b))]`.
That source has SHA256
`66edc19b0b431994d5ab7f5350de26a7336684ec1d2a2f69fb84c074fed19cde`.
Thus list position `i` is precisely old node `i`, and `dict(enumerate(mp))`
preserves every lookup used in the adaptive distance calculation.

The complete native orders are 22 -> 9 and 25 -> 10, from the census records
at lines3134 and3686 of
`results/fresh_high_rank_quasilanner_census_20261017.json`, SHA256
`6fe2aa3a6b5cd9ef70e25c2efea1d484f2f2c4a17219a9e1ea1966c8bc0e2729`.
These are the only native identifiers in the already fixed nine-case domain.
They supply an independent expected old-node domain, rather than inferring
completeness from an arbitrary map's own length.

`copy_map(value,old_order,new_order)` accepts an exact indexed list or a
dictionary with integer keys or canonical decimal-string keys. It requires
the full domain `range(old_order)`, bounds both orders by64, rejects aliased
integer/string keys, requires exact integer images in the physical core,
and requires injectivity. Boolean keys or images are rejected. For any
accepted dictionary, normalization preserves each old-node lookup; for the
producer list, normalization does the same by indexed enumeration. The
subsequent `du`, `dv`, inverse, correlated profile, shortcut, seal, tail and
whole-bound checks are byte-identical to760 apart from the single replaced
copy-map expression. No numerical assertion, geometry premise, source input,
batch name, status, provider key or output assertion was weakened or changed.

## Controls and remaining actual work

The original20 control calls remain. Two additions check the complete saved
nine-case list/dictionary equivalence through all three branches and refuse
a dictionary whose size is correct but whose integer domain is incomplete.
The control count is22. Fixtures remain representation/association checks;
they do not validate inverses, the authoritative2614-domain theorem or the
whole-support geometry.

Only the root may review and run these sources through the existing bounded
wrapper. Both the22 controls and the successor's actual `adaptive9` batch
remain UNRUN at this source freeze. The final checker owner must select the
successor worker identity explicitly; no old or failed receipt may be
relabelled. Independent S8 and global finite57 remain unverified by this
adapter. The final checker and authoritative inputs were not edited here.
