# Alternate released H8 lifts and coupled H9 caps

**No new finite-volume H9 example was found.** Three previously unused released-wall seeds were tested. One yields an exact twelve-wall H9 configuration with a π/10 polar cap and a quartic cap meeting it at π/10. Adaptive construction reached a 493-cap dictionary for which numerical escape sampling stopped finding a witness. Exact angle classification and five certified escape clauses then proved that **no subset of these 493 caps completes the ten fixed seed walls to a finite-volume Coxeter polyhedron**. This is a finite-bank obstruction, not an unrestricted H9 nonexistence result.

This changes the seed and solves two dotted distances together. It does not enlarge the previous 452-cap bank on the label-3, omit-4 seed.

## Three distinct seeds

Lift endpoint 0 of the compact Bugaenko H8 example as in the [earlier exact lift classification](fresh_h8_lift_h9_20261013.md), then delete one old wall. Original wall indices below are zero based. Each retained ten-wall Gram matrix has rank ten and signature `(9,1)`; the previously archived affine-subdiagram gate leaves these three unused cases.

| Lifted endpoint label | Deleted old wall | Diagram decomposition | Hyperideal dual rays | Initial exact caps | Coxeter polar caps |
|---|---:|---|---:|---:|---:|
| 3 | 5 | Lorentzian seven-node component plus H3 | 8 | 256 | 0 |
| 3 | 6 | Lorentzian eight-node component plus I2(5) | 8 | 348 | 1, angle π/10 |
| 4 | 4 | Connected ten-node Lorentzian diagram | 9 | 108 | 1, disjoint |

The other dual rays have negative norm; none has zero norm. These are necessary seed gates, not compactness claims. The initial cap proposals have square 1 or 2, at most five nonzero pairings from `−{1,φ,2,1+φ,2φ}`, and contain the saved strict center. Every accepted norm and seed angle is checked exactly. Numerical filters only propose candidates, so the cap counts are not exhaustive classifications. Each whole initial bank has an exact common spacelike escape.

For a seed Gram matrix `G`, use dual coordinates `y_i=−(x,r_i)≥0`. A cap with seed pairings `p` has coefficients `c=G⁻¹p`, square `pᵀG⁻¹p`, and inequality `c·y≥0`. The polar opposite seed wall `i` has coefficients `−G⁻¹e_i`, square `(G⁻¹)ii`, and only one nonzero seed pairing, `−1`. The three polar censuses are exact; their squared angle values lie in `Q(φ)`, whose possible finite Coxeter labels here are 2,3,4,5,6,10.

## An exact π/10 construction on the omit-6 seed

The valid polar `s` opposite original wall 5 has

\[
(s,s)=\frac{6-2\phi}{5},\qquad (s,r_5)=-1,
\qquad \cos^2\angle(s,r_5)=\frac{2+\phi}{4}.
\]

Thus its angle is π/10; it is orthogonal to the other nine seed walls. Unlike the failed natural polar of the previous released seed, this is an exact Coxeter cap.

Now prescribe a norm-two cap `v` with only two nonzero old pairings:

\[
(v,r_0)=-x,\qquad (v,r_9)=-(2+\phi),
\]

and require `(v,s)=−1`. The remaining eight old pairings are zero. The norm equation is

\[
\frac\phi2 x^2+(-7+4\phi)x-16+8\phi=0.
\]

The selected root is

\[
x=-11+7\phi+(\phi-1)\sqrt{49-24\phi}>2.
\]

The second nonzero old pairing is also strictly dotted. The cap meets `s` at π/10 and is orthogonal to eight old walls. The discriminant has field norm 649, a nonsquare, so this cap uses a genuine quartic extension of `Q(φ)`.

A separate [independent first-cap audit](fresh_h8_alternate_first_cap_independent_audit_20261013.json) verifies all 66 pair relations, Gram inertia `(9,1,2)`, a strict timelike interior, and strict relative-interior witnesses for every one of the twelve facets. It also verifies that the cap cuts the targeted earlier escape. This is a real twelve-wall H9 Coxeter configuration, **not** a compact or finite-volume example. No arithmeticity or commensurability novelty is asserted.

## Coupling two dotted distances

More generally choose two free pairing magnitudes `x,z>2`, a small set of fixed magnitudes `1,φ`, and an existing cap `a`. Prescribing `(v,a)=−M`, with `M∈{0,1,φ}`, expresses `z=L+Tx`. Substitution in the norm-two equation leaves a quadratic over `Q(φ)`. The accepted branch is checked exactly for its norm, all seed angles, the prescribed anchor relation, the strict center, and violation of the current escape.

On the connected label-4/omit-4 seed, three rounds produced 73,3,0 new caps. All 76 coupled caps together with the 108 initial caps retain a strict common spacelike escape. Its valid disjoint polar was separately checked; that polar also does not complete the initial bank.

On the omit-6 seed with the π/10 polar, four successive rounds produced 42,29,47,26 caps. Each batch cuts the preceding witness. This gives 348 initial caps, one polar, and 144 coupled caps: **493 optional normals**, in addition to the ten fixed seed walls. The last all-cap system gave no spacelike LP witness in 400 trials. This numerical outcome is recorded as `UNDECIDED_NO_NUMERIC_ESCAPE`; no vertex census, finite-volume theorem, or compactness claim follows from it.

## Compatible subsets and exact escape clauses

The initial selection graph contains 86,058 exactly accepted pairs with labels 2,3,4,5,6,10 or strictly disjoint hyperplanes. It excludes parallel pairs because compact completion is the target. Its numerical proposal filter and initially restricted finite-angle alphabet mean an absent edge is not, by itself, a proved incompatibility.

Twenty maximal greedy subsets of this graph have 20–26 total facets and exact strict escapes. A second search used only necessary positive-outward-pair conflicts, deliberately allowing unknown finite-angle pairs. It generated twenty further escaped subsets and remained satisfiable after the bounded run. These subsets need not all be Coxeter; their escapes nevertheless supply valid necessary constraints on any completion from the dictionary.

For each spacelike witness `y`, every finite-volume completion must select a cap with `c·y<0`. The verifier recomputes this complete cutter list over all 493 normals. There are 41 distinct covering clauses after four further escaped subsets in the accepted-edge graph. A standard-library Boolean search independently proves that this graph has no clique meeting every clause, in only 23 visited states. This replay uses no SAT solver.

The proof reduces to five clauses with cutter counts **29,12,12,13,14** and a union of only **54 nodes**. Every compatible set covering these clauses can be reduced to nodes in that union. The [independent complete angle audit](fresh_h8_alternate_core_angles_independent_20261013.json) classifies all 1,431 pairs among those nodes. A pair lies in a field of degree 2,4 or8, so its possible finite Coxeter labels belong to the complete safe list

```
2,3,4,5,6,8,10,12,15,16,17,20,24,30,32,34,40,48,60.
```

The audit constructs exact real-cyclotomic polynomials and isolating intervals for their largest roots. Exact comparisons in the biquadratic compositum classify 782 positive pairs, 408 forbidden finite-angle pairs, 64 valid finite pairs, 175 dotted pairs, and **two parallel pairs**. There are no missing valid compact-compatible edges. Adding both parallel pairs, `(20,114)` and `(21,114)` in the 493-node order, still gives UNSAT. The [independent cover replay](fresh_h8_alternate_core_cover_independent_20261013.json) checks all five complete cutter lists over the full dictionary and archives a 17-state Boolean proof with the parallel edges included. Thus the result excludes finite-volume completion as well as compact completion.

Redundant optional inequalities do not create an exception. The ten seed walls remain true facets in every selected subsystem: all dictionary caps pair nonpositively with each seed wall, and their common strict timelike center projects orthogonally onto that seed wall while making every other inequality strict. Discard redundant optional caps, then apply the complete angle graph and escape clauses to the actual facet normals. If an escape lies on some retained cap boundaries, adding a sufficiently small positive multiple of the strict center makes all inequalities strict while preserving its positive square. Hence failing to cut any one of the five witnesses still forces infinite volume.

All escape certificates have positive Lorentz square, negative pairing with the common strict timelike center, and satisfy their selected inequalities. The archived subset escapes are strict. Interpolating with the center crosses an open patch of the future ideal boundary, proving infinite volume for each tested subset. The normalized seed cone `Σy_i=1` lies in a bounded simplex, so a future positive compactness proof would require a complete appropriate vertex census, not merely the absence of sampled escapes.

## Replay and artifacts

The [standalone verifier](../scripts/fresh_h8_alternate_verify_20261013.py) uses only standard-library `Fraction` arithmetic. It reconstructs the three lifted seeds directly from the original H8 data, then checks accepted normals, polar angles, all coupled equations, every escaped subset, the complete cutter lists, the positive-pair conflicts, and the Boolean obstruction including both parallel edges. Forged timelike escape replacements are rejected. The separate complete-angle audit uses exact SymPy algebra only to construct cyclotomic polynomials and rational root isolations; pair signs and comparisons use exact `Fraction` arithmetic. The cover audit independently checks the decisive five clauses and Boolean proof.

```sh
python scripts/fresh_h8_alternate_verify_20261013.py --verify
PYTHONPATH=. python scripts/fresh_h8_alternate_first_cap_review_20261013.py
PYTHONPATH=. python scripts/fresh_h8_alternate_core_angles_independent_20261013.py --verify
PYTHONPATH=. python scripts/fresh_h8_alternate_core_cover_independent_20261013.py --verify
```

Artifacts:

- [Standalone replay result](../data/enumerations/fresh_h8_alternate_verify_20261013.json).
- Initial seed/cap archives `data/enumerations/fresh_h8_alternate_{3_5,3_6,4_4}_20261013.json`; [producer](../scripts/fresh_h8_alternate_caps_20261013.py).
- Their `_polar_20261013.json` polar censuses; [producer](../scripts/fresh_h8_alternate_polar_census_20261013.py).
- [76 coupled caps on the connected seed](../data/enumerations/fresh_h8_alternate_coupled_20261013.json); [producer](../scripts/fresh_h8_alternate_coupled_20261013.py).
- [144 adaptive coupled caps with the π/10 polar](../data/enumerations/fresh_h8_alternate_polar_coupled_20261013.json); [producer](../scripts/fresh_h8_alternate_polar_coupled_20261013.py).
- [Accepted-edge graph and twenty compatible subsets](../data/enumerations/fresh_h8_alternate_selection_20261013.json).
- [Positive-conflict Boolean search and 37 complete escape clauses](../data/enumerations/fresh_h8_alternate_sat_20261013.json).
- [Four further escapes and accepted-graph UNSAT](../data/enumerations/fresh_h8_alternate_graph_sat_20261013.json).
- [Five-clause, 54-node reduction](../data/enumerations/fresh_h8_alternate_graph_core_20261013.json).
- [Complete angle-audit code](../scripts/fresh_h8_alternate_core_angles_independent_20261013.py) and [independent finite-volume obstruction replay](../scripts/fresh_h8_alternate_core_cover_independent_20261013.py).
- [Independent review of the five-clause obstruction](fresh_h8_alternate_core_independent_report_20261013.md).

The search remains restricted to these changed seeds and explicit finite cap banks. More subset heuristics cannot complete the closed 493-cap bank. A further construction must introduce caps outside it—for example new anchor relations or three freed pairings—or change the retained seed. Other released seeds and unrestricted cap families remain open; no unrestricted H9 completion obstruction or new finite-volume example is claimed.
