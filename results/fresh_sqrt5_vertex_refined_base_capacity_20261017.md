# Prism-base capacity at each ambient vertex type

Generation and a fresh saved replay pass. Refining the shared prism-base capacity by the **full ambient finite vertex forest** detects1,381 violated inequalities in the preceding explicit-provenance witness. Repairing its simplex masses gives

\[
 \frac{8498530120100410}{705623875797781}
 =12.043994557995845\ldots>12.
\]

The refinement therefore strengthens the current obstruction but does not yet remove it. This is an exact partial-model dual certificate, not a dimension bound or a realization of the incidence distribution in a compact polytope. No optimizer was invoked.

## Labelled shared capacity

The [vertex-weighted theorem](fresh_vertex_weighted_prism_base_capacity_20261017.md) gives, for a maximal simplex `S` of dimension `d` and an ambient vertex-forest label `F`,

\[
 \sum_\rho m(\rho,S,F)\alpha_\rho
 \le (13-d)c(S,F)\alpha_S.
\]

Here `c(S,F)` counts vertices of `S` with full ambient forest `F`, and `m(rho,S,F)` counts vertices of all bases of the actual prism provenance `rho` that lie in `S` and have that forest. The proof applies the per-vertex, per-immediate-coface capacity from the [shared-base theorem](fresh_prism_simplex_shared_capacity_20261017.md). Summing over `F` recovers the previous unlabelled capacity; the separate inequalities retain where the base vertices lie.

Every coefficient uses the **same frozen provenance masses** from the [actual split](fresh_sqrt5_actual_provenance_dual_20261017.md). No label receives an independently selected representative or a coordinatewise minimum. If `c(S,F)=0`, an actual base vertex cannot have label `F`; the reconstruction checks this structural zero condition before constructing any inequality.

## Full vertex reconstruction

The [source](../scripts/fresh_sqrt5_vertex_refined_base_capacity_20261017.py) reconstructs all4,969 containing-simplex profiles from their connected core and full spectator forest. Every spherical one-root deletion gives a full rank13 finite forest.

For each of the1,858 actual quad provenances, its complete `K2,2` census identifies every vertex of the unique dotted prism. Each of its bases is then mapped to its containing maximal simplex by the saved literal core isomorphism. The mapped deleted root must be a spherical deletion of that simplex, and both sides must have exactly the same full ambient finite forest. All14,284 vertex comparisons pass, covering6,016 base occurrences. This includes every vertex of each base, not a chosen paid corner.

There is a useful distinction between a full wall star and its maximal prism. The wall star can have additional positive-definite two-deletion vertices outside that prism; these may enlarge one base into a larger simplex. For example, canonical quad core8 has four prism vertices and an additional PD deletion `(2,4)`. The prism uses the union of all its complete `K2,2` partitions; the containing simplex uses **all** its spherical one-deletions. The source checks both sets and their inclusion explicitly. Treating all PD deletions as prism vertices would incorrectly erase the proper-base distinction.

## A literal violation and the repair

Labelled column4939 concerns base type56033, simplex face80245, of dimension4. Its core is11 and its spectator forest is `3A1 + 3A2`. Its five vertices have model forest indices `419,199,718,718,716`. Forest419 is

`4A1 + 3A2 + A3`.

Actual provenance1660, source context103609, contributes one vertex of this forest to the relevant base. Its mass is `15699/4000000000`. Since this forest occurs once in the containing4-simplex, the labelled capacity is9. The old simplex mass was `5233/30000000000`; the refined inequality forces `5233/12000000000`, a factor of `5/2` increase. This is a concrete failure of the previous total-only capacity.

The analytic repair retains every old face mass and increases each simplex mass to the maximum of its labelled load/capacity ratios. It raises627 simplex rows and activates no new row. All quad provenance masses and merged quad totals remain unchanged.

The saved certificate rechecks:

- all11,397 labelled transfer inequalities and their1,381 old violations;
- all211,339 ordinary pair, clique and one-port positive-cost columns;
- all59,962 mixed-three-port columns;
- all5,693 inherited signed transfer columns and all4,969 actual unlabelled shared constraints;
- the7,044 positive face masses and797 positive vertex masses.

Mixed columns require no additional vertex mass beyond the other positive-cost columns. The exact numerator is `65373308616157/4950000000000`; the normalizer is `705623875797781/643500000000000`.

## Saved interface and limits

The [JSON](fresh_sqrt5_vertex_refined_base_capacity_20261017.json) stores each containing simplex's literal omitted-root/forest map, every actual base vertex and its mapped simplex vertex, every labelled signed column and its provenance coefficients, the old violations, and the repaired face and vertex masses. Complete mixed row features remain referenced through the frozen actual-provenance package; their coefficients do not change when masses change.

The capacity labels here are full vertex forests. Distinct pointed simplex-star vertex orbits with the same forest are still grouped. Other global compatibility conditions, missing face classes and realization in one arithmetic quadratic form remain outside. The positive support and explicit uniform split are unchanged; this is one analytic repair, not an optimum.

The theorem received separate narrow mathematical reads, including this agent's review of the invariant-label and zero-support scope. The new source/data have generation and saved replay checks; no independent arithmetic replay of this particular package is claimed.

## Reproduction

```sh
python3 scripts/fresh_sqrt5_vertex_refined_base_capacity_20261017.py
python3 scripts/fresh_sqrt5_vertex_refined_base_capacity_20261017.py --verify
```

Generation87371 and fresh saved replay39570 both terminated with `PASS_EXACT_VERTEX_REFINED_BASE_DUAL_BARRIER`. Development assertions first exposed the extra-PD-vertex distinction described above; the corrected source generated the saved result and was then replayed unchanged.

SHA256:

- source: `a7a06b87c0ee1d13e6ce773d03a0e681afe0ca1ddedc631b160d499446eee92d`;
- JSON: `09ca1c294750aaa39a16223b5100bc19108bf69561ff8b779d913479aa2b6e2e`.

All inherited inputs and the weighted theorem are hash-bound. Frozen predecessors were not edited.
