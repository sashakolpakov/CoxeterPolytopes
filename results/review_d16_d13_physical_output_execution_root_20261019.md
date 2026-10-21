# D16/D13 physical-output audit: executed generation and separate replay

19 October2026. **Both jobs actually passed.** The reviewed fixed source
`scripts/check_d16_d13_physical_output_retention_20261019.py`, SHA256
`36c3b674e296faf6cea5da48d2c1caaf586ec503763377d27089247bf4f15d70`,
was run in mode `both` under the shared parent lock. No old producer,
forest optimizer, historical census, Lean process or TeX build ran.

## Exact execution

|Record|Seconds|Peak process-tree bytes|Start RAM GiB|Start disk GiB|
|---|---:|---:|---:|---:|
|Generation|0.506|90767360|6.48|187.31|
|Separate replay|0.472|90734592|6.27|187.32|

The59-second/1GiB wrapper,45-second child alarm, single-thread environment,
nonblocking shared lock, owner-death protection and normal resource gates
were retained. Several earlier attempts stopped below6GiB RAM before
launch; none created a scientific child or execution ledger. Both actual
jobs exited0, their locks were released, and no process remains from
these jobs. They created no build intermediates.

|Output|Bytes|SHA256|
|---|---:|---|
|`d16_d13_physical_output_retention_20261019.json`|86534|`065dfac17636aa946d495317c757262df46f2e7fc2cd058936c70be1f4323293`|
|`d16_d13_physical_output_retention_generation_20261019.json`|947|`7610e6b65a15389e044fc6103d651c1f413c78ad27d3b8f77b09ecf289e14e86`|
|`d16_d13_physical_output_retention_replay_20261019.json`|948|`2ef300de9f7f8ec0ccc3bc8787d7ec1b7a54d492c2556ac6832e9706e5738dcd`|

The54 distinct inputs total106205768 bytes. Large sections were streamed;
the certificate is under87KB. The complete coverage contract, local
models, normalization and missing-amount source semantics are in
`d16_d13_physical_output_retention_preparation_20261019.md`.

## Actual checked interface

The scan covers22 historical packets/564 corners,216 literal configurations,
19 newer physical sections, five extensions and332 positive chosen-D12
outputs. Every declared section yields physical rows. Target recognition
uses entire connected components and exact native full-model matching;
unknown finite targets fail before writing. The result has zero unknowns.

|Target classification|Occurrences|
|---|---:|
|D16 retained full guarded-unit models|20|
|D16 old full order16 cores, independent short2/rank route|4|
|A4-end/D13-short unchanged half-only requests|24|
|Complete distal A4-end/D13-short marginal geometries|8|

Repeated occurrences are saved interfaces, not distinct newly enumerated
polytopes. The four D16 guarded models and all four distal native geometries
were encountered. The fixed rational inverse/model/capacity identities,
paired ideal/distal roles and declared scalar amounts check. A hash match
alone is not substituted for a geometric proof.

## Mathematical composition now supported

Together with the separately reviewed full-seal rank lemma and actual
two-angle capacity, the audited interface supports the long-end D16
absence saving `-d16*[a1=0]` and replaces the distal finite extra by
at most `2*alpha*[a4*d13>0]`, for0<=alpha<=1/2. The old half stays
available outside the explicitly modeled marginal faces. Actual-face
routing, complete contact domains, whole-factor seals, facet projection
and historical half-only source semantics remain the written premises
already reviewed; the script is not a formal proof of those premises.

The independent sequential-cover review proves a complete SAFE finite
upper2965/8 for this refined research rule, throughout that alpha range.
It uses the saved complete original-class complement and all later
positive families, rather than only the former D13 attaining forest.
The audit has now discharged its stated saved-output gate. This is
not an exact attained maximum; the new exact maximum is undetermined.

For alpha>=9/32 the completed order11 rule's ideal envelopes only
decrease, so all2270 terminal ideal margins remain positive under
the inherited concentration framework. The unique43/128 absolute
minimum remains unchanged because its terminal profile lacks rank16.
In particular alpha=1/2 can be used with finite upper2965/8 and the
same absolute ideal surplus. No new profile or forest optimizer ran.

The finite upper is still53/8 above the curvature target364. The
general dimension ceiling therefore remains57. The subsequently
proposed D16 THIRD-node and A2/D13 rollback refinements are different
roles and are NOT certified by this audit.

Reproduce with separate guarded `--write`/`--verify` jobs as documented
in the preparation note. Existing certificates must be replayed with
`--verify`; generation refuses to overwrite them. No new Lean coverage
or external refereeing is claimed.
