# Compact21 triangle constraints: proof map and reproduction

For current reading, see the [current dimension-bounds manuscript](../arxiv/hyperbolic_coxeter_dimension_bounds/main.tex) and the
[submission guide](arxiv_submission_editions_20261028.md). The supporting
arguments and recorded verification scopes below remain distinct from those
current submission builds. Standalone draft typesetting routes have been retired.

The new supporting argument
proves local triangle restrictions at the two recurring supported compact
arithmetic rank21 forests. It supplies no new dimension ceiling or charge
allowance. The compact arithmetic ceiling remains21.

|Claim|Written proof|Exact verification|
|---|---|---|
|D5+2E8 complete nonadjacent triangle domain|[Independent source/geometric review](../results/review_compact21_d5_two_e8_triangle_coupling_20261018.md)|Staged27-domain producer/output inherited; source-level read here, no independent full domain rerun|
|At most one designated D5 triangle; E8 simplex and product faces|Same review|[Independent stdlib checker](../scripts/check_compact21_e8_triangle_bundles_20261018.py):10D5pairs,16E8single-wall deletions,20product corners|
|D5 theorem beside A1+E7+E8|[Complete new proof](../results/fresh_compact21_d5_a1_triangle_coupling_20261018.md)|[Producer/replay](../scripts/fresh_compact21_d5_a1_triangle_coupling_20261018.py):24,576vectors,139eligible,6roots,15acute exclusions|
|At most seven nonadjacent triangles at A1+D5+E7+E8|[Complete new proof](../results/fresh_compact21_supported_forest_triangle_coupling_20261018.md)|[Producer/replay](../scripts/fresh_compact21_supported_forest_triangle_coupling_20261018.py):fullD5+A1spectators,11roots,55pairs,2048subsets|
|Actual ground field at compact dimension21|[Separate internal reduction](../results/fresh_compact21_single_field_20261016.md)|Pinned completed predecessor; not re-proved by the new local scripts|

Published finite and affine Coxeter classifications, compact acute simplicity,
and Vinberg's cyclic-integrality/conjugate-positivity criterion are explicit
theoretical inputs. Actual finite corners, untouched exceptional isolation,
all high-label exclusions, positive-height orientation, compact affine
obstruction, actual-face projection, and simplex/product closure remain
written geometric arguments. The manuscript and the proof reports distinguish
them from finite arithmetic checks. None of the full geometric theorems is
currently claimed Lean-verified.

Run sequentially from the repository root:

```sh
python3 -B scripts/check_compact21_e8_triangle_bundles_20261018.py
python3 -B scripts/fresh_compact21_d5_a1_triangle_coupling_20261018.py --verify
python3 -B scripts/fresh_compact21_supported_forest_triangle_coupling_20261018.py --verify
```

The new producers use the shared exclusive proof flock,59-second alarm and
periodic one-GiB RSS guard. Their complete replays take less than a second
each in the recorded runs, with peak RSS under32MiB. They import no NumPy,
optimizer or Lean process and produce no temporary proof objects; `-B`
suppresses bytecode cache. The standalone identity checker must likewise be
run sequentially with other proof work. Do not launch a proof farm.

The new producers were generated and replayed in separate processes. The
initial D5+2E8 review and independent small checker passed a root mathematical
read. The [D5+A1 review](../results/review_compact21_d5_a1_triangle_coupling_20261018.md)
and [full supported-forest and TeX review](../results/review_compact21_supported_triangle_coupling_20261018.md)
also pass. These mathematical/source reads are distinct from generation
and exact replays. The PDF
compiled in b818c0 with no final-pass warnings; unused TeX files were removed
after the batch. The first draft's long-path overfull box was fixed before
that final build. No Lean compilation was used for these new lemmas.

Any attempt to use the new triangle limits in a lower dimension theorem
must still handle actual quadrilaterals, preserve all old face-payment
consumers, introduce a valid new charge inequality, and verify it against
the complete global coverage and exact weight model. No such step has been
silently assumed here.
