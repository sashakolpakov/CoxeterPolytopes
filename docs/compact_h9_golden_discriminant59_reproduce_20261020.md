# Reproduce the discriminant-59 host exclusion

The full integral reflection arrangement of E8(Z[phi]) plus the binary
Gram matrix `[[-2phi,-1],[-1,4+2phi]]` is NONREFLECTIVE. Here phi²=phi+1;
the determinant ideal has prime norm59. A complete third edge supplies
an integral full-chamber isometry of trace40+55phi>10 and infinite order.
Its powers give infinitely many distinct cosets modulo the full
reflection group. Stop further searches of this literal host.

The [closure proof](../results/compact_h9_golden59_nonreflectivity_conclusion_20261020.md),
[complete third-edge proof](../results/compact_h9_golden59_d7_edge_interval_20261020.md)
and [period proof](../results/compact_h9_golden59_long_facet_period_20261020.md)
supply the global mathematical argument. The closure's pending-execution
sentence describes its earlier written edition; the four completed
executions below supersede that historical status. No new example,
dimension ceiling or whole commensurability-class exclusion follows.
The [final independent execution/closure review](../results/review_compact_h9_golden59_final_execution_batch_guard_20261020.md)
binds the actual certificates and ledgers to that unchanged proof.

The [paper proof](../results/compact_h9_golden_discriminant59_edge_route_20261020.md)
has SHA256
`725d8775cebec51d8c18b4b0a0d3be0255e8b1472621adb4d3cefa656adb429e`.
Root independently read the complete proof and worker. The
[independent paper/source review](../results/review_compact_h9_golden_discriminant59_batch_guard_20261020.md)
also passes, SHA256
`308646e58c366e0971d2889fd886f605c1443465a8d06325cba9a7781ffa356b`.

## Reproduce both exact interfaces

From the repository root:

```sh
(
  set -e
  for coxeter_worker in check_compact_h9_golden_discriminant59_edges_20261020 check_compact_h9_golden59_d7_first_edge_20261020
  do
    python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock -- \
      "scripts/${coxeter_worker}.py" --verify
  done
)
```

The worker SHA256 is
`01a569b66b77c38056255a94dce765e4784d809dfe86735ac5bef869755edf3d`.
It reconstructs the sixteen known normals, all120 pairs, the integral
chamber involution, a strict interior point and44 distinct finite vertices.
For all five second-edge profiles it verifies the full quadratic norm
identity, contact identity and the exact inequalities in the written
all-height exclusion. Its six evaluation points determine the complete
quadratic polynomials; this is not numerical sampling. No height loop,
optimizer or complete cone census is used.

The [certificate](../results/compact_h9_golden_discriminant59_edges_20261020.json)
is44,121 bytes, SHA256
`23c650dd70cb297af4923bbc5dbd30e19f672d6329c7c014ea90567a85831a99`.
The source pins the proof, quadratic-field helper and calibrated Bugaenko
data, and hashes all of them before and after execution. Replay compares
its newly reconstructed result with the saved certificate; it does not
overwrite the certificate. Generation and replay use the same algorithms,
so they are separate executions, not independent implementations.

The third-edge/period worker SHA256 is
`5022037cc17357e1cbe745a0fc9be40b445372ad266d00bdf7bc2236dc0d70a9`.
It exhausts two proved integer domains of52 and68 pairs, reproduces both
ten-entry norm tables, selects the unique first wall, and verifies the
full-lattice period matrix, integral inverse, actual corner frames and
trace. Its [certificate](../results/compact_h9_golden59_d7_first_edge_20261020.json)
is18,175 bytes, SHA256
`b549fcaa03c907572a38aef550d6ba11b63c65d781dc2a1300f4fb66b9dd4214`.
The [independent combined proof/source review](../results/review_compact_h9_golden59_d7_period_batch_guard_20261020.md)
and [independent interval review](../results/review_compact_h9_golden59_d7_interval_resource_20261020.md)
both pass. The latter reviewer authored the period proof and explicitly
does not count as an additional independent review of that isometry.

## Written inputs and root review

The necessary primitive norm pool {2,2b}, where b=7+5phi, uses evenness,
reflection integrality, the adjugate, and the golden ring's PID/unit facts.
It is not inferred from a bounded root list. Root checked this argument
and the whole-arrangement separation of long mirrors. The latter includes
the two nonfundamental finite-angle values 2-phi and3-phi; leaving them
out would make the original complete-list assertion incorrect, although
their norms preserve the prime59 exclusion.

The actual base facet uses the inherited COMPLETE Bugaenko arrangement,
not just eleven proposed walls. The E8 transverse edge is exhausted by
the written two-embedding norm bounds. At its E8+c endpoint, finite
classification leaves exactly five profiles on the E7+c edge. Here
X=x-y is integral and E8 unimodularity supplies integral Z. The bounds
force X to be a unit or a nonzero algebraic integer of norm below one;
the unit case gives the unique short wall u. No first-wall candidate is
discarded by checking old-cap compatibility prematurely.

The actual pentagon follows from these two complete edges, the base
facet and positive Gram projection of the remaining finite corner.
This geometric implication is a written argument. The worker checks its
algebraic interfaces and44 finite points; it does not replace the written
norm inventory, spherical classification, golden unit theorem or acute
projection argument. No geometric Lean verification is claimed.

## Recorded resource-safe executions

|Execution|Result|Time|Peak process-tree RSS|
|---|---|---:|---:|
|[Generation](../results/compact_h9_golden_discriminant59_edges_generation_20261020.json)|PASS|4.237s|56,868,864 bytes|
|[Separate replay](../results/compact_h9_golden_discriminant59_edges_replay_20261020.json)|PASS|4.428s|56,623,104 bytes|
|[Third edge/period generation](../results/compact_h9_golden59_d7_generation_20261020.json)|PASS|0.506s|57,294,848 bytes|
|[Third edge/period replay](../results/compact_h9_golden59_d7_replay_20261020.json)|PASS|0.497s|57,212,928 bytes|

Each job separately acquired the shared proof lock and enforced the
unchanged40GiB disk /6GiB available-memory launch gates,59-second timeout,
1GiB process-tree limit and one numerical thread. Launch resources were
over179GiB disk and8GiB available RAM. Python used `-B`; no temporary
build files were produced. There was no Lean, TeX or numerical-library
build. The source-preparation note's “not executed” describes its earlier
historical stage; the four ledgers above record the later executions.

All four runs took9.668 seconds in total, below58MB peak process-tree
RSS. No temporary build files were produced. The finite domains and
matrix checks now discharge the computational interfaces in the written
closure. The full-chamber/coset implication remains a reviewed geometric
proof, not a claim of geometric Lean verification. No additional TeX
manuscript was built for this single-host exclusion.
The later shared [chamber-period logic module](chamber_period_logic_verification_20261020.md)
now kernel-checks the final group-theoretic implication conditionally.
The actual matrix, geometry and coset interpretation remain explicit
inputs; the full host exclusion is not geometrically Lean verified.

The [immutable 25-input packet](../results/compact_h9_golden59_nonreflectivity_20261020_freeze.json)
pins214,487 bytes of proofs, reviews, sources, data and execution evidence.
Every listed size and SHA256 was rechecked after writing the manifest.
Mutable status documents remain outside it. Historical pending-status
records are preserved unchanged; the final review records their closure.

The historical sixteen-wall seed has uncapped c-edges. Its initial44
finite vertices and compact pentagon were valid local progress, but did
not imply full-host compactness. The third-edge infinite-period proof
now excludes the entire host without a full cone census. The earlier
square-root-of-three exclusion concerns a different literal lattice.
