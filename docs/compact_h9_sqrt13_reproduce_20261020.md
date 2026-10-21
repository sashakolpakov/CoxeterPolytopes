# Reproduce the quadratic13 H9 host exclusion

20 October 2026. The literal lattice is
`L=<-eta>+I9` over `O=Z[omega]`, where `omega²=omega+3` and
`eta=omega+1=(3+sqrt13)/2`. Its full integral reflection arrangement is
now proved nonreflective. This excludes one literal construction host;
it supplies neither a new example nor a uniform dimension bound.
The [completed closure](../results/compact_h9_sqrt13_nonreflectivity_conclusion_20261020.md)
records the complete proof chain and later execution status.

The [first route](../results/compact_h9_sqrt13_b9_two_edge_route_20261020.md)
proves the complete primitive reflection-norm pool `{1,2}`, the full
initial B9 germ, and two complete finite edges. They lead through
`B8+A1` to another B9 vertex. The resulting integral full-chamber map
is an involution, of trace6, and does not exclude the host. Its
[independent written review](../results/review_compact_h9_sqrt13_two_edges_resource_20261020.md)
passes.

The [third-edge proof](../results/compact_h9_sqrt13_a8_edge_interval_20261020.md)
retains A8 at the initial vertex. All eight finite extension placements
are covered. Two isolated profiles are excluded at different primes
above3; exact physical/conjugate bounds leave eight norm tests. There
is one first simple wall, with norm2, leading to an actual A9 vertex.
Both [independent](../results/review_compact_h9_sqrt13_a8_edge_resource_20261020.md)
[written reviews](../results/review_compact_h9_sqrt13_a8_edge_finite57_20261020.md)
pass. Its source hash is
`e3e62a6ec4d834e2e0c7d4abc4d28961839f1a467ac06966ef08854757814228`.

The [fourth edge and period](../results/compact_h9_sqrt13_a9_far_end_and_period_20261020.md)
release the far end of the A9 chain. An integral diagonal basis and exact
mod9 reconstruction reduce all eight finite profiles to three coefficient
values. One first wall remains. A second integral chamber involution H
maps the full actual germs; the product HG has trace268eta+80>10 and
therefore infinite order. Its powers give infinitely many distinct cosets
modulo the full reflection group, proving infinite chamber volume.
Root and two independent complete mathematical reviews pass.

## Exact reproduction

From the repository root, replay the workers serially:

```sh
python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock -- \
  scripts/check_compact_h9_sqrt13_two_edges_20261020.py --verify
python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock -- \
  scripts/check_compact_h9_sqrt13_a8_edge_20261020.py --verify
python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock -- \
  scripts/check_compact_h9_sqrt13_far_end_period_20261020.py --verify
```

Generation and fresh-process replay have both passed for each worker.
Each replay recomputes exact arithmetic and compares it with the saved
certificate; it is a separate run of the same implementation, not a
second independent implementation.

|Worker|Source SHA256|Generation / replay ledgers|
|---|---|---|
|First two edges|`6839ac3b07c8bc87056712692e526deacc27bdcc13d33f7fb07a72fe8e64af31`|[generation](../results/compact_h9_sqrt13_two_edges_generation_20261020.json), [replay](../results/compact_h9_sqrt13_two_edges_replay_20261020.json)|
|A8 edge|`63dd91e189e3534024b6e2a2223adfacd5633c4a65cd87ef56177dea0a340bda`|[generation](../results/compact_h9_sqrt13_a8_edge_generation_20261020.json), [replay](../results/compact_h9_sqrt13_a8_edge_replay_20261020.json)|
|Far A9 edge and period|`ee78514c118d8054be993f86b04b7ab003ae980c8296dca4c6cf14faea6ff599`|[generation](../results/compact_h9_sqrt13_far_end_period_generation_20261020.json), [replay](../results/compact_h9_sqrt13_far_end_period_replay_20261020.json)|

All six runs used the shared exclusive proof lock, one thread, the
59-second live watchdog and a 1GiB process-tree limit. Together they took
2.918 seconds; peak process-tree RSS was 57,868,288 bytes. Launches had
over 178GiB free disk and over 8.5GiB available RAM. No build intermediates
were created. Sources, certificates and compact ledgers are retained.
The [immutable 33-input packet](../results/compact_h9_sqrt13_nonreflectivity_20261020_freeze.json)
pins 171,514 bytes of proof, source, certificate and verification inputs.
All sizes and hashes were rechecked after writing the manifest.

## Verification boundary and preserved support

The scripts check the fixed algebra and finite domains. Completeness of
the primitive alphabet, spherical extension inventory, local finiteness
and the full-chamber geometric interpretation are explicit written
inputs. These are complete written proofs, not conclusions inferred from
the scripts alone. The final infinite-volume conclusion also uses the
reviewed actual-germ and arithmetic-lattice/index arguments. The shared
[ChamberPeriod module](chamber_period_logic_verification_20261020.md)
checks the final group/coset implication conditionally; actual geometric
instantiation remains outside Lean. No external refereeing is claimed.

A [separate short-facet/C3 argument](../results/compact_sqrt13_short_facets_and_c3_compact_prefix_20261020.md)
proves exact inheritance of the coordinate-short facet within this
literal family and supplies a compact nine-wall three-dimensional
prefix. It has no executed checker yet and is not a compact-H9 example.
The [prime5 obstruction](../results/compact_h9_sqrt13_a9_reversal_prime5_obstruction_20261020.md)
shows that the automatic A9 diagram reversal fails even to permute the
integral mirror arrangement. It excludes that proposed map, not the host.
The later integral point inversion in the far-edge proof supplies a
different, valid symmetry and completes the exclusion. The prime5 note
remains correct; it is not needed for that completed proof.

Stop further edge searches of this literal H9 lattice. No additional TeX
manuscript is proposed. The uniform ceiling table is unchanged, and other
forms and lattices in its commensurability class are not excluded.
