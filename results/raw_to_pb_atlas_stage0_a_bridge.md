# Raw Stage 0 to PB-Atlas Stage A

**Status (2026-08-10): proved and executable numerical interface.**

The composed compiler no longer asks the caller to supply a canonical class,
physical wall scale, boundary list, or boundary weight.  Its external input is
only

\[
     (G_{\operatorname{Pic}};r_1,\ldots,r_N),
\]

an integral unimodular Lorentzian Gram matrix and primitive oriented
full-span wall vectors.  Raw Stage 0 exhausts the smooth-rational adjunction
solutions and either routes a direct weak-anticanonical branch, rejects the
packet, or constructs one complete closed PB-Atlas input for each accepted
index-two candidate.  The four frozen Coble calibrations each have a unique
accepted candidate.

This closes the **numerical interface only**.  It does not prove that a
surface or curve configuration exists.

## Interface theorem

Let a raw Stage-0 candidate be accepted as an index-two log-boundary branch.
Write its physical boundary curves as

\[
                        B_i=d_i r_i.
\]

Stage 0 has already proved

\[
                K+\frac12\sum_iB_i=0.
\]

Consequently

\[
                    \sum_i d_i r_i=-2K.                 \tag{1}
\]

The PB-Atlas boundary indices are therefore precisely the inferred Stage-0
boundary indices, and their weights are the inferred physical scales `d_i`.
Every remaining Stage-A field is copied from the raw input or from the chosen
accepted canonical candidate:

```text
picard_gram             <- raw Gram
wall_vectors            <- raw primitive walls
canonical_class         <- derived K
boundary_indices        <- derived nonweak physical walls
boundary_weights        <- derived physical scales on those walls
```

The independent Stage-A parser then rechecks unimodularity, signature, full
wall span, primitive acute orientation, characteristic parity, `K^2=10-rho`,
and equation (1).  Thus the bridge is not a format-only copy.

## Exact calibration and its provenance

The four raw Coble packets reproduce the previously frozen PB-Atlas inputs
exactly, byte-for-structure:

| raw packet | derived Stage-A branch | boundary weights |
|---|---|---|
| `H10/30` | Type II index-two | `[2]` |
| `H11/26` | Type I index-two | `[1,1]` |
| `H14/60` | Type VI index-two | `[1,1,1,1,1]` |
| `H17/48` | Type MII index-two | eight weights `1` |

The distinction in the first column is geometric and essential: `H10` has a
twice-primitive physical boundary `B=2r`, while the other rows have primitive
square-`-4` boundary classes and require the inferred half-boundary module.

The phrase `byte-for-structure` is a regression identity, not a claim that
all six fields have independent provenance.  The raw Gram and wall arrays
were deliberately projected from these same PB-Atlas fixtures, so their
equality after the round trip is tautological.  The three fields at issue in
the interface theorem are noncircular: `canonical_class`,
`boundary_indices_zero_based`, and `boundary_weights` are reconstructed by
Raw Stage 0 before the frozen PB fixture is consulted for comparison.  The
solver and converter receive neither that fixture nor the separate
post-solve `expected` metadata in the raw calibration file.  Hostile tests
poison both sources and verify that the derived payload is unchanged while a
poisoned expected payload fails the equality check.  Thus the whole-payload
comparison is useful as a regression check, while only its three derived
fields constitute independent inverse evidence.

Two negative controls are retained.  The toric `H6/9` packet is routed to the
direct weak-anticanonical theorem without inventing a boundary.  The
reflection-maximal `Bl_p(F_4)` triangle is rejected at `K+B/2 != 0` and is
never admitted to PB-Atlas.  Hence maximality is not being smuggled in as a
log-boundary certificate.

## Exact scope

The result proves that the **combined numerical** Stage 0--A pipeline has a
genuinely raw interface and that no row-specific canonical or boundary
metadata is needed in the four calibrations.  It still does not prove:

- effectivity or irreducibility of the numerical wall classes;
- existence or nonemptiness of a rational surface incidence stratum;
- completeness of the Mori cone or absence of hidden curves;
- literal same-polytope return or essential surjectivity.

Those are the geometric obligations of PB-Atlas Stage B and its eventual
return theorem.

## Reproducibility

- bridge certifier: `scripts/certify_raw_to_pb_atlas_stage0_a_bridge.py`;
- frozen certificate: `data/enumerations/raw_to_pb_atlas_stage0_a_bridge.json`;
- tests: `tests/test_raw_to_pb_atlas_stage0_a_bridge.py`;
- upstream raw theorem: `results/raw_adjunction_log_stage0.md`;
- upstream contraction theorem: `results/pb_atlas_numerical_contraction_stage_a.md`.
