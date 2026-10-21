# Root source review of the fixed R5 checker

20 October 2026. PASS by complete source inspection, before execution.

Reviewed all400 lines of
`scripts/check_golden_dyadic_r5_chamber_period_20261020.py`, SHA256
`af564d36a57dd29e1ef0b3c17e30e827aa66d5279381c566d76688fbaebe4658`,
and preparation note SHA256
`27db896fd103ff7e0397596beacfb1a45d33994223768336d2bac627ce124470`.
The separate outside-author source audit is
`results/review_golden_dyadic_r5_period_fixed_checker_batch_guard_20261020.md`.

The original lattice basis is checked in both directions against a full
D4 basis. The inverse is obtained by its Gram adjoint, avoiding the
helper's nonpivoting inverse on an arbitrary indefinite matrix. All Gram
matrices passed to that inverse have the stated nonzero block pivots.
The new residues and positive projection are reconstructed independently.
All six finite profiles remain, with the separate B2 contact square and
the simultaneous primitive normals g and Q retained. The complete written
embedding bounds imply the finite coefficient rectangle; the program
checks its exact strips and the necessary integrality exclusions.

The period is reconstructed from the old apex glue and the V0 relation,
then compared against every entry of the independently written matrix.
The restrictive parser does not evaluate code. Gram preservation, both
integral inverse identities, all nine positively oriented simple-root
images, the future apices, and trace10+11phi are checked. Recursive frozen
input checks, pre/post source hashes, an exclusive certificate write,
fresh replay, a128KiB certificate cap, and the required bounded parent
are consistent with the contract. No unbounded search occurs.

The geometric classification, actualness and completeness of the germs,
and chamber-germ implication remain explicit paper proofs. They are not
established by arithmetic assertions or by Lean. For the final obstruction,
an infinite-order full-chamber stabilizer gives infinitely many distinct
cosets of the reflection group in the arithmetic lattice: a nonidentity
reflection-group element cannot stabilize its chamber. The standard
arithmetic-lattice finite-covolume theorem is an explicit external input.
This proves a literal-lattice obstruction, not a new dimension ceiling or
an exclusion of every lattice in the commensurability class.

No scientific execution, certificate generation, replay, Lean or TeX is
claimed by this source review. Root may next schedule one bounded
generation and a separate fresh replay.
