# Global Coxeter LP/SDP research

The mathematical specification and actual-polytope lift are in
[global_formulation.md](../../notes/lp/global_formulation.md). The retained
[source notes](../../notes/lp/global_sources.md) identify literature inputs.
This is an initial research relaxation; no improved dimension bound or
separate author's acceptance is claimed.

The organizing objects are occupation measures for complete vertex
populations, all triangles and quadrilaterals, residual corners, and
global incidence counts. A single shared Lorentz Gram schema covers
every small face. Complete component families have symbolic rank tails.
All26 intrinsic manuscript features remain common marginals, with free
dual coefficients. This code contains no selected list of vertex forests
or marked face diagrams and performs no diagram-family pricing.

The initial sparse moment level is deliberately explicit. It can be
strengthened by raising moment/localizer degrees or adding proved
uniform identities. The LP version replaces PSD blocks by declared
necessary rational-direction tests. Its failure to find a certificate
does not resolve the SDP or the geometric problem.

From the repository root, run one monitored batch at a time:

```sh
python3 -B scripts/run_bounded_scientific_long_10gb_20261123.py --parent-lock \
  --ledger results/global_tail57_export_monitor.json -- \
  scripts/run_global_coxeter_formulation.py export \
  --threshold 57 --parity odd --upper-dimension 57 \
  --output results/global_tail57_export
```

Use fresh output names. `--upper-dimension57` is written as
`--upper-dimension 57` on the command line. It uses the completed
n>=58 theorem as an anchor. Thus one global model for the band N<=n<=57
suffices for a proposed cutoff N; separate diagram families or separate
runs for each dimension are not needed. Omitting the option keeps the
whole unbounded parity tail in the new model itself.

After the source/cleanup commit is pushed, `explore-lp` exports both the
necessary-data problem and its certificate dual and searches the latter.
`explore-sdpa` sends the complete conic certificate problem to SDPA after
exact equality elimination. Both retain numerical status separately from
exact rational identities and PSD checks. Every accepted result must be
checked in its original complete matrix; support reduction is only a
candidate-recovery device.

`explore-psd-cuts` adds rational quadratic tests from violated original
PSD blocks. It lifts any recovered LP certificate back to the original
SDP and checks it there. This adds constraints to the same global
matrices, without growing a list of diagram types. The model includes
the uniform old two-unit triangle cover, contextual rank-support rows,
and positive-coordinate support localizers. Old quadrilateral transfers
whose premises require the dimension-58 LP remain outside this first model.

The angle-square observable couples inverse-Gram sums over all finite
vertex components to the face cosines. Its A/B/D formulas and rational
unbounded-tail envelopes are proved in the specification. Its dual
coefficient is free along with the 26 old feature coefficients.

The generic compiler includes affine equalities in its Farkas alternative
and includes the factor two in off-diagonal trace pairings. It does not
fix corner coefficients, impose coefficient boxes or silently replace a
partial vertex/face list by an exhaustive one. Results store source hashes
instead of duplicating the whole source tree at every run.
