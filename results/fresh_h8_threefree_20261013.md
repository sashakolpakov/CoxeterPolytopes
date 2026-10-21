# Three free pairings and two anchors on the released H8 seed

**Subsequent completion:** the [next search](fresh_h8_cap_completion_20261013.md) closes this 1035-entry dictionary for compact completion, with an exact proof allowing the polar to be optional. Finite-volume completion from the dictionary remains unresolved. The bounded status and counts below describe this earlier checkpoint.

**No new finite-volume H9 example was found.** This search goes outside the previously closed 493-cap bank. It constructs **540 distinct new exact mirrors**, represented by 542 archived cap entries, then tests 153 pairwise Coxeter-compatible configurations with 26–33 facets. Every tested configuration has a strict exact spacelike escape and therefore infinite volume. The expanded compatible-subset search remains incomplete; there is no obstruction theorem for the new bank or the unrestricted construction.

## Construction

Keep the ten-wall label-3, omit-original-wall-6 seed and its valid π/10 polar from the [previous report](fresh_h8_alternate_20261013.md). The ten independent normals have norm2, Gram matrix `G` of signature `(9,1)`, and retained original indices

```
0,1,2,3,4,5,7,8,9,10.
```

Write `I=G⁻¹`, `w_i=−(v,r_i)`, and require `(v,v)=2`. The new construction frees three magnitudes `w_i>2`, so those seed hyperplanes are strictly disjoint from the new cap. At most three other magnitudes are nonzero, each chosen from `{1,φ}`. The first anchor is the polar, node348 of the old bank; the second is one of60 old `Q(φ)` caps exactly compatible with it. Prescribe the two anchor pairings as negatives of values in `{0,1,φ}` and retain only exact valid Coxeter relations or strict disjointness.

The two anchor equations leave an affine line

\[
w=w_0+xw_1,\qquad
A x^2+B x+C=0,
\]

where `A=w₁ᵀIw₁`, `B=2w₀ᵀIw₁`, and `C=w₀ᵀIw₀−2`, all in `Q(φ)`. A positive-discriminant branch gives a root over an extension of degree at most four over `Q`. Its coefficients in the seed basis are `c=−Iw`. The exact acceptance gates check the norm, all ten seed angles, both anchor angles, the three strictly dotted seed pairings, the common strict timelike center, and the targeted escape cuts. A numerical proposal filter makes no completeness claim for this family.

For each accepted root, an exact comparison of normalized squared seed pairings proves that its mirror differs from every entry of the old493 bank. A separate [mirror audit](../data/enumerations/fresh_h8_threefree_mirror_aliases_20261013.json) compares every possible equal-dotted-support pair among the new entries. It finds only two aliases, pairs `(22,23)` and `(197,198)`, both in the first batch. The archived node numbering is preserved. Thus the first280 entries represent278 new mirrors, and all262 entries from the adaptive second batch add new mirrors: **540 total**.

## A small explicit new cap

Combined cap337 has the negative seed-pairing magnitudes

\[
w=(0,0,2\phi,0,1,\phi,0,\sqrt{16-7\phi},\phi,1+2\phi)
\]

in the retained original-wall order above. Here `16−7φ>4`. Its norm-two equation is

\[
\frac{4+2\phi}{5}x^2-10+2\phi=0,
\qquad x=\sqrt{16-7\phi}.
\]

It is orthogonal to the polar and meets the norm-two old anchor194 with pairing `−φ`, hence angle π/5. Its three dotted old walls are original2,8,10; the remaining seed angles are right angles, π/3, or π/5. The field norm of `16−7φ` is95, a nonsquare, so this cap genuinely uses a quartic field. It cuts all six escape directions used to guide the adaptive batch. This is a concrete additional cap, not a finite-volume completion or a claim about arithmeticity or commensurability.

## Adaptive compatible-subset search

The first pass targets all five certified escape directions from the old finite-bank obstruction. It keeps three scored proposals per anchor and target, yielding280 accepted cap entries. A witness-driven Boolean search fixes the polar, selects compatible caps, greedily completes each selected set, and seeks a spacelike escape. Every proposed selected pair is classified exactly before its volume test.

Pairs of roots can lie in a compositum of degree eight. The exact pair gate therefore uses the full safe finite-angle list

```
2,3,4,5,6,8,10,12,15,16,17,20,24,30,32,34,40,48,60,
```

with real-cyclotomic polynomials, rational isolating intervals, and exact biquadratic sign evaluation. The preliminary proposal graph uses numerical filters and the previously accepted old-old edges; its missing edges are not an exhaustive incompatibility certificate. Parallel pairs are excluded because the target is compact completion.

The first selection round tests53 exact compatible configurations, with26–32 facets, and generates53 new escapes in addition to the five inherited ones. Its proposal graph becomes UNSAT. Six escape clauses, involving25 proposed cutter nodes, suffice to reproduce that restricted-graph failure. This is search guidance, not a complete new-bank theorem.

The second construction pass targets those six newly found escape directions with the same bounds and the same two-anchor mechanism. It produces262 additional cap entries. The combined bank has1035 entries including the old493;699 are active in the polar-retaining proposal search. One hundred further exact compatible configurations, with28–33 facets and up to11 new cap entries each, all escape. The final archive contains158 complete active-node cutter clauses. Its status is `BOUNDED_SEARCH_INCOMPLETE`; compatible combinations remain untested.

For a witness `y`, the replay checks `y_i>0`, `yᵀIy>0`, `yᵀIh<0` for the common timelike center `h`, and `c·y>0` for every chosen cap. Interpolation from `h` toward `y` crosses an open patch of the future ideal boundary, proving infinite volume. All selected normals pair nonpositively and share a strict timelike center. Orthogonally projecting that center onto any selected wall preserves timelikeness and makes every other inequality strict, so the selected walls are genuine facets. The ten seed normals retain rank ten. None of these conclusions depends on absence of numerically found escapes or an incomplete vertex census.

## Exact replay

The [main replay](../scripts/fresh_h8_threefree_verify_20261013.py) imports the previously independent exact field arithmetic, but no cap producer, floating-point library, LP code, or SAT solver. The initial run checks280 accepted entries,3576 exact pair classifications,53 compatible escaped configurations, and58 complete active-node cutter clauses. The expanded run checks542 accepted entries,6306 pair classifications,100 compatible escaped configurations, and158 cutter clauses. A [second independent audit](fresh_h8_threefree_independent_audit_20261013.json) separately replays all280 initial roots and their distinction from the old493 mirrors.

```sh
PYTHONPATH=. python scripts/fresh_h8_threefree_independent_20261013.py --verify
PYTHONPATH=. python scripts/fresh_h8_threefree_verify_20261013.py --verify
PYTHONPATH=. python scripts/fresh_h8_threefree_verify_20261013.py \
  --source data/enumerations/fresh_h8_threefree_combined_caps_20261013.json \
  --selection data/enumerations/fresh_h8_threefree_combined_selection_20261013.json \
  --output data/enumerations/fresh_h8_threefree_combined_verify_20261013.json --verify
PYTHONPATH=. python scripts/fresh_h8_threefree_mirror_audit_20261013.py --verify
```

Primary artifacts:

- [Three-free, two-anchor producer](../scripts/fresh_h8_threefree_caps_20261013.py) and [initial cap archive](../data/enumerations/fresh_h8_threefree_caps_20261013.json).
- [Initial53 selected configurations and escape clauses](../data/enumerations/fresh_h8_threefree_selection_20261013.json), [selection code](../scripts/fresh_h8_threefree_selection_20261013.py), and [exact replay](../data/enumerations/fresh_h8_threefree_verify_20261013.json).
- [Six-clause proposal-graph guidance core](../data/enumerations/fresh_h8_threefree_core_20261013.json) and [core extraction code](../scripts/fresh_h8_threefree_core_20261013.py).
- [Adaptive262 cap entries](../data/enumerations/fresh_h8_threefree_adaptive_caps_20261013.json), [combined542 entries](../data/enumerations/fresh_h8_threefree_combined_caps_20261013.json), and [merge code](../scripts/fresh_h8_threefree_merge_20261013.py).
- [Combined100 selected configurations and158 clauses](../data/enumerations/fresh_h8_threefree_combined_selection_20261013.json) and [expanded exact replay](../data/enumerations/fresh_h8_threefree_combined_verify_20261013.json).
- [Exact mirror-alias replay code](../scripts/fresh_h8_threefree_mirror_audit_20261013.py).

The next constructive opening is a new compatible combination or additional cap relations that preserve a promising partial configuration. The present bounded search neither proves that this expanded bank is impossible nor rules out an unrestricted completion of the changed H9 seed.
