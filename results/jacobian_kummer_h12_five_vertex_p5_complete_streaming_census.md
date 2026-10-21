# Jacobian--Kummer `H12`: complete five-vertex `P5` census

Status: exact positive `SURFACE -> POLYTOPE` certificate, frozen 2026-08-10.
Worldwide priority remains open.

## Complete result

The frozen 112-wall Jacobian--Kummer `H12` factory has exactly **201,747**
seed-action and path-reversal orbits of properly colored five-vertex paths.
Exact boundary reconstruction and ordinary-Coxeter nonisomorphism certificates
prove that all 201,747 yield pairwise nonisomorphic finite-volume noncompact
`H12` reflection polytopes of reflection index and volume multiplier five.

| exact disjoint stratum | types | facets |
|---|---:|---:|
| `B<=2` | 20,693 | 752, 761, 770 |
| `ABBB` | 36,247 | 779 |
| `BABB` | 36,648 | 779 |
| `BBBB` | 108,159 | 788 |
| **complete `P5`** | **201,747** | **752--788** |

The exact type/facet histogram is
`752^142, 761^2091, 770^18460, 779^72895, 788^108159`.
Facet count separates every pair of strata except `ABBB` versus `BABB`;
literal degree-profile comparison separates those two 779-facet strata.

The newly completed `BABB+BBBB` queue contains 144,807 census rows.  The
project increment beyond previously completed census strata is therefore
144,807.  One original five-path pilot was already a `BBBB` row and no pilot
was `BABB`, so the increment beyond all earlier individual `P5` examples is
144,806.  The complete 201,747-row census contains all five original pilots
and contributes 201,742 further types.

## `BBBB` exact census

The palindromic `BBBB` shard has 41,154,000 raw unrooted proper colorings.
Explicit orbit cover and Burnside agree on exactly **108,159** rows; the
Burnside denominator is 768 and the reversal fixed sum is 659,328.  Every row
has 788 facets.

Sorted exact ordinary-label degree profiles give 100,935 classes.  The 6,366
collision groups contain 13,590 rows.  Literal degree bytes were replayed for
every member, and exact canonical color refinement separates every group.
In total 13,591 current rows were refined: the extra row is the unique
cross-star candidate described below.  Of these, 13,232 reach a discrete
canonical form.

Every row passes the primitive/integral reflection-root, strict timelike
interior, fixed Smith-`1^13` core, full-rank, and ordinary-Coxeter relation
gates.  The common relation inventory is `m2=13,584`, `m3=672`, `m4=3,940`,
`m6=0`, `infinity=291,882`; the global minimum strict-interior facet pairing
is 5.

## Exact cross-shape accounting

The complete star layer has 8,985 rows with the same 788-facet header as
`BBBB`.  Digest filtering found exactly one candidate pair.  Both rows were
rebuilt; their literal degree-profile bytes are equal.  Exact canonical color
refinement was then run on both matrices, and the resulting invariants are
different.  Thus the pair is separated by exact refinement, not by a
hash-only inference.

Together with the frozen star--`ABBB` theorem and the new star--`BABB`
separation, this proves the exact high-`B` union

`16,872 + 36,247 + 36,648 + 108,159 = 197,926`

pairwise nonisomorphic polytopes.  It does **not** prove star plus all `P5`
pairwise nonisomorphic: star versus `B<=2` remains open.

## Streaming and memory certificate

The two large shards were produced one boundary at a time.  SHA-256 serves
only as a negative equality filter and compact linkage; every equal-digest
bucket is replayed using literal bytes, and every literal degree collision is
sent to exact canonical refinement.  No collection of Gram matrices and no
database is archived.

The successful `BBBB` command was:

```text
python3 scripts/run_with_rss_limit.py --memory-gib 1.0 --timeout 9000 -- \
  python3 scripts/certify_jacobian_kummer_h12_five_vertex_p5_streaming_shard.py \
  --pattern BBBB
```

It completed in 4,963.676 seconds.  Both the internal and external RSS ceilings
were 1.0 GiB; the sampled peak was **0.406 GiB**.  The frozen `BBBB` output is
a 7.3 KiB summary and a 14 MiB deterministic gzip JSONL transcript.

Frozen `BBBB` digests:

- summary SHA-256:
  `8e89f7450687838ab2fdff6e4fa73be0d401e928f80956ff866b84b927258000`;
- gzip-row SHA-256:
  `d367b7d208ad8c7fff294978a1857f7557d3755f33df8e4afff6689ad2aec1a7`;
- uncompressed canonical JSONL SHA-256:
  `a559e54135e539b396f4b311e04a98b495a7d7e618a1b1da895afa9ec345777a`.

The tiny combined complete-`P5` archive has SHA-256
`dd6c60a8ed02e985a889100991d2da410d8b18f0902fa314e1e6e9dc900dab16`.

## Carrier and claim boundary

All rows retain the weak full-Picard K3 carrier
`U(2)+D5(-1)^2+<-4>` on a seven-dimensional family.  No full-nef,
effectivity, or surface-action claim is made.  The degree-three fork census
and the all-five-vertex-shapes census remain open.  Allcock's prior abundance
results are respected; worldwide novelty and priority remain open.

## Frozen artifacts

- [`certify_jacobian_kummer_h12_five_vertex_p5_streaming_shard.py`](../scripts/certify_jacobian_kummer_h12_five_vertex_p5_streaming_shard.py)
- [`certify_jacobian_kummer_h12_five_vertex_p5_complete_streaming_census.py`](../scripts/certify_jacobian_kummer_h12_five_vertex_p5_complete_streaming_census.py)
- [`jacobian_kummer_h12_five_vertex_p5_bbbb_streaming_census.json`](../data/enumerations/jacobian_kummer_h12_five_vertex_p5_bbbb_streaming_census.json)
- [`jacobian_kummer_h12_five_vertex_p5_bbbb_streaming_census_rows.jsonl.gz`](../data/enumerations/jacobian_kummer_h12_five_vertex_p5_bbbb_streaming_census_rows.jsonl.gz)
- [`jacobian_kummer_h12_five_vertex_p5_complete_streaming_census.json`](../data/enumerations/jacobian_kummer_h12_five_vertex_p5_complete_streaming_census.json)
- [`test_jacobian_kummer_h12_five_vertex_p5_bbbb_complete_streaming_census.py`](../tests/test_jacobian_kummer_h12_five_vertex_p5_bbbb_complete_streaming_census.py)
