# H9 completion with batches of exact escape constraints

**No new finite-volume H9 example was found.** A different selection method resolves the previously unfinished dictionary: no subset of its1035 entries, together with the ten fixed seed walls, gives a **compact** Coxeter polytope, even when the polar cap is optional. This is an exact finite-dictionary result; it does not exclude finite-volume completion or unrestricted new caps.

The resulting bottleneck also produces49 additional accepted cap entries outside that dictionary. Nine preserve a coherent30-facet configuration. All eight maximal compatible extensions from those nine caps have31–32 facets and exact infinite-volume certificates. Dropping or replacing facets, other caps, and other seeds remain open.

## A different completion strategy

The [preceding search](fresh_h8_threefree_20261013.md) ended with a still-satisfiable proposal graph on1035 cap entries: the old493 plus542 three-free cap entries representing540 distinct additional mirrors. Its exact root and selected-configuration audits were rerun successfully before the present work.

The saved model selects cap nodes

```
348,416,473,518,689,713,789.
```

All selected pair angles pass exact checks. Together with the ten seed walls it has17 facets, but its volume is infinite. Instead of stopping after the first random LP escape, the new method enumerates candidate vertices of the bounded affine section

\[
y_i\ge0,\quad\sum_i y_i=1,\quad c_j\cdot y\ge0.
\]

Numerical halfspace intersection supplies proposals only. A small perturbation toward the saved strict center and rational reconstruction then give exact positive-norm, future-oriented, strict witnesses. Complete cutter lists are calculated over **all1035 entries**, including caps excluded from the initial polar-retaining proposal graph. Witnesses are selected for different cutter sets, prioritizing the most restrictive clauses.

The saved model produces24 such clauses. The next selection seeks redundant coverage of difficult escape clauses and completes its clique using rarity-weighted coverage instead of graph degree alone. This gives a30-facet configuration and24 further exact clauses. Only this one additional maximal clique is needed to make the original proposal graph unsatisfiable. Allowing the polar to be omitted also leaves that numerical graph unsatisfiable. These preliminary graph outcomes are not used as exact incompatibility proofs.

## Exact compact obstruction for the1035-entry dictionary

With the polar required as an actual facet, one witness suffices. Exactly13 dictionary entries cut it:

```
0,1,4,5,8,10,17,18,43,44,58,94,345.
```

Every one is exactly incompatible with the polar. This statement concerns a retained polar **facet**, not a possibly redundant polar inequality.

The optional-polar proof removes that distinction. Thirteen complete escape clauses have cutter counts

```
19,32,6,9,10,11,13,18,20,22,24,24,23.
```

Their union contains147 nodes. All10731 pairs in that union are classified exactly using the complete safe finite-angle list for pair fields of degree at most eight. The counts are:

| Pair type | Count |
|---|---:|
| Positive outward pairing | 2364 |
| Valid finite Coxeter angle | 1006 |
| Other finite angle | 3017 |
| Parallel | 153 |
| Strictly disjoint | 4191 |

A45-state Boolean proof, replayed without a SAT solver, shows that no clique in the compact-compatible graph meets every clause. Parallel pairs are excluded here: two parallel reflection mirrors yield a parabolic isometry, which cannot occur in the reflection group of a compact hyperbolic Coxeter polytope. The argument uses all possible finite labels, not the proposal graph's omitted-edge guesses.

Redundant optional cuts do not evade the conclusion. All dictionary normals pair nonpositively with each of the ten seed normals and share a strict timelike center. Orthogonal projection of that center onto a seed wall preserves timelikeness and makes every other inequality strict, proving that every seed wall remains a genuine facet. Discard redundant optional inequalities. The actual facet set of a compact Coxeter completion must be pairwise compatible and must cut each exact spacelike witness, contradicting the Boolean certificate. If a witness lies on some optional boundaries, a sufficiently small positive center perturbation makes those inequalities strict while retaining positive Lorentz square.

The [producer-free replay](../scripts/fresh_h8_cap_completion_20261013_verify.py) checks the complete1035 cutter list for each core witness, all10731 exact pair classifications, the common center, and every branch of the45-state Boolean proof. It also rechecks the48 new strict witnesses for the17- and30-facet configurations.

This theorem is **compact-only**. Adding the153 parallel edges leaves an exact core solution with nodes

```
30,34,76,78,99,137,167,232,239,377,639,661.
```

The resulting22-wall configuration passes the finite-angle/parallel/disjointness gates but also has infinite volume, certified by eight exact strict escapes. Its failure does not prove that the full dictionary has no finite-volume completion.

## New caps that preserve a coherent partial configuration

The single polar bottleneck is used as a new target for the existing three-free, two-anchor equations. The finite parameter bounds are unchanged. There are2227 numerically accepted proposals and49 exactly checked retained cap entries. Every retained cap is polar-compatible and cuts the bottleneck witness; the thirteen old cutters are polar-incompatible. Consequently every retained cap is outside the entire1035-entry dictionary, not just the original493.

All49 are compared exactly with every wall of the30-facet configuration. Nine are compatible with all of them, at new-batch indices

```
3,4,6,9,10,12,24,32,37.
```

Their complete36-pair graph has exactly eight maximal compatible subsets:

```
{3,32}, {4,32}, {6}, {9,32}, {10}, {12}, {24,32}, {32,37}.
```

All eight resulting31–32-facet configurations have exact strict spacelike escapes. Twelve independently checked witnesses are retained for each, with complete cutter lists over the expanded1084-entry dictionary. Thus these targeted caps do close the particular old bottleneck while leaving other open ends. The compatible-extension statement preserves all30 facets of the fixed partial configuration; it does not classify constructions that discard or replace them.

The [extension replay](../scripts/fresh_h8_cap_completion_20261013_extensions_verify.py) checks the49 cap equations through the producer-free root audit, all exact compatibility decisions, all512 subsets of the nine eligible entries to establish the maximal-clique list, and96 exact escape witnesses with their complete1084 cutter lists. No compactness claim is inferred from numerical vertex counts.

## Replay and artifacts

```sh
PYTHONPATH=. python scripts/fresh_h8_cap_completion_20261013_verify.py --verify
PYTHONPATH=. python scripts/fresh_h8_cap_completion_20261013_extensions_verify.py --verify
```

- [Search with batches of escape clauses](../scripts/fresh_h8_cap_completion_20261013.py) and [two tested configurations](../data/enumerations/fresh_h8_cap_completion_20261013.json).
- [Optional-polar proposal search](../data/enumerations/fresh_h8_cap_completion_20261013_optional_polar.json).
- [Exact core extraction and pair classification](../scripts/fresh_h8_cap_completion_20261013_core.py), [polar-facet bottleneck](../data/enumerations/fresh_h8_cap_completion_20261013_core.json), and [optional-polar compact core](../data/enumerations/fresh_h8_cap_completion_20261013_optional_core.json).
- [Complete compact-obstruction replay result](../data/enumerations/fresh_h8_cap_completion_20261013_verify.json).
- [Parallel candidate and eight exact escapes](../data/enumerations/fresh_h8_cap_completion_20261013_parabolic_candidate.json).
- [New49 cap entries](../data/enumerations/fresh_h8_cap_completion_20261013_bottleneck_caps.json) and [their exact root replay](../data/enumerations/fresh_h8_cap_completion_20261013_bottleneck_caps_verify.json).
- [Exact compatibility with the fixed30-facet configuration](../data/enumerations/fresh_h8_cap_completion_20261013_target_compatibility.json).
- [All eight maximal extensions](../data/enumerations/fresh_h8_cap_completion_20261013_extensions.json), [extension producer](../scripts/fresh_h8_cap_completion_20261013_extensions.py), and [extension replay result](../data/enumerations/fresh_h8_cap_completion_20261013_extensions_verify.json).

The old1035-entry dictionary is closed only for compact completion. The expanded1084-entry dictionary has not been classified, and this bounded cap family does not exhaust possible completions of the seed.
