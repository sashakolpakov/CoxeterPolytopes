# Golden D4 dimension9 research packet

The open candidate is `R_5=<-phi>+I_5+D4(O_phi)`.
It is **not yet a compact or finite-volume example**. The global
dimension tables are unchanged. This packet is supporting research,
with no additional TeX manuscript or Lean project.

Completed own proofs give four actual finite edges and a sealed H3
face with exact full perpendicular lattice `M6=<-phi>+I6`.
An integral projection/determinant argument excludes norm4 H4
overgroups and proves the seal for the full arrangement. The
[four-edge proof](../results/golden_dyadic_d4_fourth_vertex_and_sealed_h3_face_20261020.md)
and its linked earlier sources have complete
[root](../results/review_golden_dyadic_d4_edges_root_20261020.md) and
[outside-author](../results/review_golden_dyadic_d4_chain_resource_20261020.md)
written reviews.

One complete edge inside that face now identifies all nine actual M6
walls and an integral chamber involution. Its
[finite reduction](../results/golden_dyadic_m6_actual_short_edge_preparation_20261020.md)
has full written audits, and the
[execution review](../results/review_golden_dyadic_m6_execution_root_20261020.md)
records successful generation and fresh-process replay. The census
retains both simultaneous integral normals; one is a positive
nonsimple root. It does not silently discard that candidate.

From the repository root, reproduce the finite arithmetic certificate:

```sh
python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock -- scripts/check_golden_dyadic_m6_actual_edge_20261020.py --verify
```

This is one serial process with the shared lock, disk/RAM gates,
59-second watchdog, one thread and 1-GiB process-tree cap. The worker
has an additional45-second alarm and128-KiB certificate limit. It
refuses to overwrite a certificate. The recorded successful jobs took
0.225s total, peaking at53.8MiB; no Lean, TeX, caches or build
intermediates were produced. Full source hashes are in the certificate
and compact ledgers. The immutable packet manifest is
`results/golden_dyadic_d4_research_20261020_freeze.json`.

The checker covers all nine finite continuation profiles,924 exact
coefficient bounds, complete field-square reconstruction, lattice glue,
old-wall acuteness, the two-root identity/tie, endpoint and full reversal
matrix. Geometry, finite-type completeness and the explicitly cited
compact P6-double lineage are written inputs. The lineage identifies a
known compact M6 chamber; it is not a new lower-dimensional example or
a uniqueness classification. No geometric Lean verification is claimed.

The next task is global compactness or a genuine obstruction for R_5.
Compactness of its face and a finite-order face involution do not settle
that task. Continue complete actual edges or structural arguments;
do not reopen an unbounded root prefix or compile Lean before the
research result stabilizes.
