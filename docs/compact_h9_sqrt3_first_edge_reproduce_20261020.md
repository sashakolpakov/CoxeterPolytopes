# Reproduce the square-root-of-three host exclusion

The full integral reflection arrangement of
`L=<-sqrt3>+<1>^9` is nonreflective. This completed search exclusion uses
a scalar facet, two exact edges and an infinite-order chamber isometry.
It supplies no new compact polytope or uniform dimension bound, and does
not exclude every lattice in the same quadratic-form class.

The [closure proof](../results/compact_h9_sqrt3_nonreflectivity_conclusion_20261020.md)
binds the computational interfaces to the geometric argument. The norm
classification, complete spherical-extension inventories and reflection
chamber facts remain written proofs. No additional TeX manuscript is
proposed. The later shared conditional Lean coverage is recorded below.

## One command for all four computational interfaces

From the repository root, run this sequential block:

```sh
(
  set -e
  for coxeter_check in diagonal_first_edge a7_opposite_edge e8_alpha2_edge a8_period_isometry
  do
    python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock -- \
      "scripts/check_compact_h9_sqrt3_${coxeter_check}_20261020.py" --verify
  done
)
```

Each separate batch acquires the shared proof lock, requires40GiB disk
and6GiB available RAM, and enforces59 seconds,1GiB process-tree RSS and
one numerical thread. Each worker has an additional45-second alarm and
256KiB certificate cap. A failed resource gate stops the sequence; do
not bypass it. Replays reconstruct certificates without replacing them.
They are separate executions of the same algorithms, not second
independent implementations or independent geometric proofs.

## What the checks establish

|Interface|Exact conclusion|Recorded generation / replay|
|---|---|---|
|Full L first A8 edge|Using the written complete norm pool, the first mirror has an actual A8+A1 corner and supplies the scalar facet.|0.104s /0.229s|
|Facet's opposite A7 edge|The original A8 vertex and the E8 vertex belong to the same full facet chamber. Three simultaneous candidates reduce to one simple E8 wall.|0.248s /0.243s|
|E8 alpha2 edge|All eighteen finite profiles are exhausted. Exactly one candidate occurs: r5, the first simple wall, giving the second actual A8 vertex.|0.217s /0.208s|
|A8 period|The reversed vertex-frame map and its inverse are integral; Gram, glue, apex images and trace5+3sqrt3>9 all check exactly.|0.314s /0.276s|

The [facet reduction](../results/compact_h9_sqrt3_scalar_facet_next_edge_20261020.md)
and [integral splitting](../results/compact_h9_sqrt3_e8_integral_splitting_20261020.md)
identify the ENTIRE lattice as E8(O) plus <1-sqrt3>, with no discarded
glue. The [A7 interval](../results/compact_h9_sqrt3_a7_opposite_interval_20261020.md)
and [alpha2 interval](../results/compact_h9_sqrt3_e8_alpha2_interval_20261020.md)
prove complete finite boxes. No arbitrary height cutoff is used.
The [period proof](../results/compact_h9_sqrt3_a8_period_isometry_20261020.md)
was written conditionally before execution. Its first-simple-wall
premise is now discharged by the alpha2 certificate. The isometry
worker deliberately reports only its algebraic scope; the closure
proof supplies the geometric composition.

The later shared [chamber-period logic module](chamber_period_logic_verification_20261020.md)
kernel-checks the last group-theoretic implication from trace excess and
full chamber preservation to infinitely many distinct cosets. Its
geometric, matrix and actual-group interpretation premises remain
explicit. This does not upgrade the whole exclusion to geometric Lean.

Because the isometry preserves the full mirror arrangement and matches
the two inward corner germs, it preserves their full chamber. Its powers
give infinitely many distinct cosets modulo the full reflection group.
A compact original host would have the proved scalar facet compact as
well, contradicting its infinite chamber stabilizer.

## Evidence and resource use

All source, interval and global-interface arguments have root and
independent written reviews. The final
[independent closure review](../results/review_compact_h9_sqrt3_nonreflectivity_conclusion_finite57_20261020.md)
also binds the actual certificates and six new generation/replay ledgers.
Original generation and all seven later
executions above pass. The latter seven took1.735 seconds in total,
with maximum process-tree RSS56,229,888 bytes. Each ran separately
under unchanged resource gates. No Lean, TeX or numerical-library build
was run for this packet. Python used `-B`; no temporary build files
were produced.

The [original first-edge packet](../results/compact_h9_sqrt3_first_edge_20261020_freeze.json)
is immutable historical evidence: its replay-deferred status describes
that earlier moment. The completed
[49-input successor packet](../results/compact_h9_sqrt3_nonreflectivity_20261020_freeze.json)
records later replay and the global conclusion. Its pinned files total
281,972 bytes; every size and hash was rechecked after writing the manifest.
Mutable status documents stay outside both packets. No worldwide novelty
claim is made for this exclusion.
