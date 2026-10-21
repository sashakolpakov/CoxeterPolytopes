# Jacobian--Kummer `H12`: complete streaming `BABB` five-vertex `P5` census

Status: exact positive `SURFACE -> POLYTOPE` certificate, frozen 2026-08-09.
Worldwide priority remains open.

## Result

The `BABB` binary edge-class shard contains exactly **36,648** seed-action
and path-reversal orbits, obtained both from an explicit orbit cover of all
14,008,320 raw unrooted proper colorings and from Burnside division by the
384-element folded seed action.  Every representative reconstructs a
779-facet ordinary Coxeter matrix.

All 36,648 matrices are pairwise nonisomorphic.  Sorted exact ordinary-label
degree profiles give 34,583 classes.  The 1,950 collision classes contain
4,015 rows; exact canonical color refinement separates every class.  Of the
4,015 refined rows, 3,990 reach a discrete canonical form.  Hashes are never
used as positive equality evidence: every same-digest degree bucket is
replayed and compared using literal bytes, and an equal refinement digest
would likewise trigger literal canonical-word replay.

The new shard is also exactly separated from all previously completed `P5`
rows and the complete five-vertex star layer.  Facet count separates it from
the `B<=2` `P5` tranche.  Against the overlapping 779-facet header, all 36,247
`ABBB` rows and all 6,117 779-facet star rows were screened from their frozen
exact per-row degree digests.  There were zero common digests.  This is a
rigorous inequality proof because equal byte strings necessarily have equal
digests; only digest equality, had it occurred, would have required literal
replay.

Consequently the completed `P5` stock now contains

`20,693 + 36,247 + 36,648 = 93,588`

pairwise nonisomorphic types.  The new project increment is 36,648 because no
original five-path pilot lies in `BABB`.  The remaining `P5` queue is exactly
the **108,159** `BBBB` orbits.  The high-`B` exact cross-shape union currently
proved is star + `ABBB` + `BABB`, namely

`16,872 + 36,247 + 36,648 = 89,767`

types.  This is not a star-plus-all-`P5` theorem: star versus `B<=2` remains
open.

## Arithmetic and geometric gates

Every row passed the frozen exact boundary engine's gates:

- five chambers with exactly four internal path mirrors;
- 779 primitive integral reflective boundary roots;
- a strict timelike interior, with global minimum facet pairing 5;
- the fixed `W1--W60` Smith diagonal `1^13` and full rank 13;
- only ordinary relations `m2,m3,m4,m6,infinity`, with common inventory
  `m2=13,442`, `m3=656`, `m4=3,909`, `m6=0`, `infinity=285,024`;
- Allcock reflection index and volume multiplier five; and
- finite volume and noncompactness.

The descendants retain the same weak full-Picard K3 carrier
`U(2)+D5(-1)^2+<-4>` on a seven-dimensional family.  No full-nef,
effectivity, or surface-action claim is made.

## Streaming and memory certificate

The producer handles one boundary at a time and retains only compact row
digests and collision indices.  It never materializes all Gram matrices and
never archives a database.  Large exact degree strings and canonical words
live only for the current replay row/group.  The successful production run
used both an internal and external 1.0 GiB RSS ceiling; sampled peak RSS was
**0.325 GiB**.  The compact output is a 7.2 KiB summary plus a 4.6 MiB
deterministic gzip JSONL transcript, rather than another 100 MiB-scale JSON.

Production command:

```text
python3 scripts/run_with_rss_limit.py --memory-gib 1.0 --timeout 3600 -- \
  python3 scripts/certify_jacobian_kummer_h12_five_vertex_p5_streaming_shard.py \
  --pattern BABB
```

The exact rebuild entry point is the same command with `--check-archive`.
The successful run completed in 1,482.259 seconds.  Its outer watchdog
reported `disposition=child_exit`, `scope=direct`, and `peak_gib=0.325`.

Frozen digests:

- summary SHA-256:
  `20f2e589fa6a216d3f52d7a92befa7df04b7d7140efa852ec6b554057bd78828`;
- gzip-row archive SHA-256:
  `04c9ff621b6ff47e690902e872204790982040cc725548ed3abd808670219e22`;
- uncompressed canonical JSONL transcript SHA-256:
  `37f3a1191c6f9f9d18cb9f929761c20372ee8d5e528644f39cf8382f0fd710d1`.

## Frozen artifacts

- [`certify_jacobian_kummer_h12_five_vertex_p5_streaming_shard.py`](../scripts/certify_jacobian_kummer_h12_five_vertex_p5_streaming_shard.py)
- [`jacobian_kummer_h12_five_vertex_p5_babb_streaming_census.json`](../data/enumerations/jacobian_kummer_h12_five_vertex_p5_babb_streaming_census.json)
- [`jacobian_kummer_h12_five_vertex_p5_babb_streaming_census_rows.jsonl.gz`](../data/enumerations/jacobian_kummer_h12_five_vertex_p5_babb_streaming_census_rows.jsonl.gz)
- [`test_jacobian_kummer_h12_five_vertex_p5_babb_streaming_census.py`](../tests/test_jacobian_kummer_h12_five_vertex_p5_babb_streaming_census.py)
