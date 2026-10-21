# C28 checker execution and boundary-excess-five composition

20 October2026. Full root source review and two separate bounded exact
executions PASS. The written proof and outside-author review are
`review_compact_ra3_c28_proof_root_20261020.md` and
`review_compact_ra3_c28_finite57_20261020.md`. They are pinned by the worker.

The stdlib worker
`scripts/check_compact_ra3_c28_graph_metric_20261020.py`, SHA256
`7356fa22ba8095e480c81cd2800f62832764925ffc6c0706be904ec8c0a6bf4e`,
was read completely before execution. Its Fraction/Q(sqrt5) operations
and physical-embedding sign test are exact. Both complexes check all
triangle incidences and connected cyclic links, degrees and Euler count.
The explicit S4 and Klein-four permutations and the regular-neighbor
closure match the written proof. The tetrahedral squared identities and
distinguished edges, general k/e normal contractions, quartic factorization,
H Gram rank/relation and positive one-dimensional unit-normal fiber are
checked by field arithmetic. These are fixed identities and two fixed
spheres, not a new exhaustive enumeration or a geometric realization test.

The source pins, worker self-hash and saved certificate are checked for
freshness. Generation refuses overwrite and replay recomputes all checks
in a separate process before comparing the complete saved JSON. The
optional provenance-only mode was not run and is not proof verification.
No old C26 worker was imported or executed. Named external and written
geometric premises remain explicit in the certificate.

Actual bounded generation and fresh replay:

|Run|Ledger|Elapsed|Peak process-tree RSS|
|---|---|---:|---:|
|Generation|`compact_ra3_c28_graph_metric_generation_20261020.json`|0.068s|56,754,176B|
|Fresh replay|`compact_ra3_c28_graph_metric_replay_20261020.json`|0.067s|56,721,408B|

Both returned PASS with the same29,156-byte certificate, SHA256
`2f9fe47902fd08af4e7627cc92dd87e4ea0dfc48a4400548c4b2c120557d9415`.
The wrapper owned the exclusive lock, capped the process tree at1GiB,
used one thread and the59-second watchdog, and checked disk/RAM gates.
No temporary build output or cache was generated. No Lean or TeX job.

## Consequences and a checked stronger charge

The two written reviews establish, using the existing geometric inputs,
that every compact RA4 polytope has f>=65 and Vol>=49pi^2/6. A positive
facet adjacent to a regular dodecahedral facet has excess d>=5. This
does not say every positive facet has excess at least5.

The unchanged global-charge source is
`compact_ra4_global_low_excess_charge_20261020.md`, SHA256
`deba2d2ed3ce5025b5c1a23de3d282bcd15762505bb15dd9c549fe068eed4549`.
Use its notation: I consists of positive facets with no D neighbor,
B the other positive facets, m=|I|, TI=sum_I d, and WII,WIB,WBB are
the corresponding nonpentagonal ridge weights. Let E5=sum_B(d-5)>=0.
Then

    |B|=(2W-TI-E5)/5,
    WBB=W-TI/2-WIB/2,
    N_irregular<=4WBB/3.

Since Nu=f-m-|B|-N_irregular, substitution gives the exact inequality

    Nu>=f-(26/15)W+(13/15)TI-m+(2/3)WIB+E5/5.

The prior universal d>=2 for positive facets gives TI>=2m, hence the
interior correction is at least11m/15>=0. With f>=65, every W<=37
therefore gives Nu>=65-(26/15)*37=13/15>0, or Nu>=1 integrally.
This substitution and consequence received a separate outside-author
algebra check from the resource agent before recording here.

At volume at most that of the120-cell, the existing uniform-D center
theorem makes its matched center interior. If the D adjacency graph
is connected, the existing propagation argument then identifies the
whole polytope with the120-cell. Thus a non120-cell of that volume
with connected D graph has W>=38 and f+W>=103. Alternatively W<=37
forces disconnected D graph and nonzero separating H2 of the full
positive-facet subcomplex, by the already reviewed Alexander-duality
argument. This is not an unproved surface assumption.

The unrestricted120-cell minimum and second compact volume remain open.
This supporting successor changes no dimension bound and adds no primary
manuscript or Lean assertion. All C26 frozen sources remain unchanged.
