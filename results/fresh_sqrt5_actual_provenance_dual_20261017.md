# An actual-provenance dual above the rank13 target

Generation and a fresh saved replay pass. Assigning the positive quadrilateral masses to their **actual saved full-star provenances**, repairing the forced simplex masses, and recounting all current charges gives the exact ratio

\[
 \frac{278543769497020170}{22887218538457519}
 =12.170276131587661\ldots>12.
\]

No optimizer was used. This removes the coordinatewise-minimum provenance relaxation from the positive support of the previous obstruction. It is an obstruction to this partial weighting family, **not** a dimension bound or a realization of these incidence frequencies in one compact arithmetic polytope.

## The old minimum vectors are not actual provenance mixtures

The [previous union certificate](fresh_sqrt5_mixed_forced_base_transfers_20261017.md) had ratio approximately12.23095 and used coordinatewise minima in its signed base-transfer columns. For a fixed prism dimension `k`, every actual provenance has the same total number of bases: four when `k=2`, and two otherwise. If the base vector varies between provenances, its coordinatewise-minimum vector has a strictly smaller total, so it cannot belong to their convex hull. The new script proves this for all243 varying positive merged rows.

There is also a literal obstruction to lifting the **unchanged old simplex masses**. Merged row10673 has mass `2041/1000000000`. Its two provenances are source contexts488 and2764, with base sets respectively `{3904,3905,3906,3907}` and `{3908,3909,3910,3911}`. Every one of those eight base simplex masses was zero. Each provenance contributes coefficient2 to each of its four bases. Thus its nonnegative mass would have to be zero under the actual shared inequalities; the two masses cannot sum to the positive merged mass. These complete data appear in `merged_row_mass_constraints`, `provenance_masses`, and `actual_shared_base_constraints` of the [saved certificate](fresh_sqrt5_actual_provenance_dual_20261017.json).

This does not prevent repairing the dual by raising simplex masses. The certificate does precisely that.

## Exact split and repair

For every positive merged quadrilateral row `Q` having `r` saved provenances, set each provenance mass to `alpha_Q/r`. All1,858 positive provenance masses therefore preserve the original1,582 merged masses exactly. No provenance is selected or omitted.

For a base type `S` of simplex dimension `d`, the actual positive shared load is

\[
 L_S=\sum_{Q,\rho} k(Q)\,m(Q,\rho,S)\,\alpha_{Q,\rho}.
\]

The repaired mass is the maximum of its old mass and `L_S/((d+1)(13-d))`. This uses the [shared prism-base capacity](fresh_prism_simplex_shared_capacity_20261017.md), the [literal base census](fresh_sqrt5_prism_simplex_base_incidence_20261017.md), and the [forced-base domain](fresh_sqrt5_forced_base_simplex_domain_20261017.md). Maximal edges retain zero curvature demand and their two endpoint one-port charges; their shared capacity is24. Square prisms use all four edge bases.

All4,969 base types appearing in the positive support have represented rows. The repair raises1,871 simplex masses, activating1,700 previously zero rows. In particular, all1,658 previously zero missing-base rows become positive, and42 represented simplex rows become positive for the first time. Their140 corners are reconstructed directly from their core and spectator graphs, and both pair and port coefficients are compared to the inherited model before computing their mixed features.

The old pair and simplex-only clique/port coefficients are constant across each merged quad row by the defining merger. The [frozen mixed-feature certificate](fresh_sqrt5_mixed_forced_base_transfers_20261017.md) checked that the mixed coefficients also agree across all1,858 positive provenances. They therefore remain valid for the explicit split. Newly positive simplex mixed features are recomputed rather than assumed covered by that earlier exporter.

## Complete checks within this support

The saved replay checks:

- all1,858 provenance masses and1,582 exact merged-mass equalities;
- all4,969 actual shared-base constraints, with6,016 nonzero provenance incidences;
- all211,339 inherited and added positive-cost columns;
- all59,962 mixed triple orbits on the797 positive vertex forests;
- all5,693 inherited conservative signed columns, as well as the actual provenance loads;
- all7,044 positive face masses and the exact numerator and denominator.

All other provenance masses are zero, so their omitted load terms are zero. Mixed columns require **no further increase** beyond the vertex masses required by the pair, clique and one-port columns. The numerator is `21642872532791/1650000000000`, and the normalizer is `22887218538457519/21235500000000000`.

The full literal distribution is local data. Global arithmetic-form compatibility, global face incidence beyond the tested constraints, and missing face classes remain outside. The uniform split is one exact feasible allocation; no optimality is claimed.

## Reproduction and integrity

Source: [fresh_sqrt5_actual_provenance_dual_20261017.py](../scripts/fresh_sqrt5_actual_provenance_dual_20261017.py). Data: [exact certificate](fresh_sqrt5_actual_provenance_dual_20261017.json).

```sh
python3 scripts/fresh_sqrt5_actual_provenance_dual_20261017.py
python3 scripts/fresh_sqrt5_actual_provenance_dual_20261017.py --verify
```

Generation18107 and fresh saved replay95306 both terminated with `PASS_EXACT_ACTUAL_PROVENANCE_DUAL_BARRIER`. The replay reconstructs the complete result and compares it to the saved JSON. The first development run encountered a tuple/list serialization-key error before writing a result; this was corrected before the successful generation and replay.

SHA256:

- source: `7f1bf158d9b49d3d977eb3e7447f2cc8cba521b1ad3dd2860335ac71b18f9a88`;
- JSON: `b9787dd4416ca1718220d6fcf27f691e7dcc6c45a29080187ef2595501b14aa8`.

Every inherited input is hash-bound by the source manifest. No frozen predecessor was edited. The fixed-total minimum-vector argument received a separate narrow mathematical read from `construction_audit`; that read is not an independent arithmetic replay of this new certificate.
