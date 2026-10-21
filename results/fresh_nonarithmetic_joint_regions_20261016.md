# Exact polar screening and a twelve-cap H9 batch

The complete 36-wall cone admits exactly **16 individually compatible exterior polar walls**, with no finite-label cutoff in that classification. Twelve of them form a certified simultaneous extension with **48 actual facets**. They cut the preceding strict escape and 541 of the 1,068 old hyperideal directions.

The extension nevertheless has **infinite volume**. Its complete cone has **172 finite vertices, no ideal rays and 2,109 hyperideal rays**. Its finite graph has 395 edges and 758 open finite-corner germs. A common exact strict spacelike escape survives all 16 possible old-cone polar walls, so no subset of that one-generation polar bank can complete finite volume. This is an exact construction and a scoped obstruction, not a new finite-volume example.

The frozen [36-wall chamber](fresh_nonarithmetic_escape35_20261016.md) and its [complete cone](fresh_nonarithmetic_escape35_20261016_cone.md) are retained unchanged. The latter now also has an [independent complete-cone review](fresh_nonarithmetic_escape35_20261016_cone_review.md).

## Complete polar classification

For an exterior extreme ray V, V²>0 and its pairings with every old normal are nonpositive. Thus its polar normal is already acute to all old walls. Its arbitrary scale is harmless: finite Coxeter compatibility is tested with

\[
c_i=\frac{4(V,n_i)^2}{V^2n_i^2}=4\cos^2(\pi/m_i).
\]

The [exact polar gate](../scripts/fresh_nonarithmetic_joint_regions_20261016_polar_gate.py) and [certificate](fresh_nonarithmetic_joint_regions_20261016_polar_gate.json) cover all 1,068 hyperideal directions of the complete input cone. Exactly 16 have all contacts right, finite Coxeter, or strictly dotted. The other 1,052 each have a certified product strictly between two consecutive values 4cos²(π/m) and 4cos²(π/(m+1)). Since these values increase strictly for every integer m≥2, each such gap excludes every finite label, including arbitrarily large labels.

Thresholds are enclosed by exact rational intervals. The program uses the Machin identity for π with alternating arctangent bounds, then 50 cosine Taylor terms with an explicit 102nd-order remainder. Each stored gap has strict disjoint rational endpoints. No numerical angle tolerance, number-field degree assumption or guessed angle alphabet supplies a rejection.

The initial numerical screen found only 15 candidates. The exact audit recovered ray 8, whose contact with wall 5 is exactly π/10; its numerical value had just exceeded the discovery tolerance. This additional wall is right to wall 33, dotted to the initial fifteen candidates, and retains their common strict escape. The exact 16-wall result supersedes the numerical candidate count. It does not retroactively make the numerical screen an exact proof.

## The exact twelve-cap chamber

The [batch producer](../scripts/fresh_nonarithmetic_joint_regions_20261016_batch.py) and [certificate](fresh_nonarithmetic_joint_regions_20261016_batch.json) reconstruct each of the initial fifteen polar directions from the saved complete-cone clipping provenance. Base rays are solved in the frozen octic field. New rays use their exact positive crossing combinations. There is no division in the additional formal radical extensions.

All 15×36 old contacts and all 105 mutual pairs are checked exactly. Among the fifteen, the incompatible pairs of candidate indices are 1–12 and 1–13, which are obtuse, and 4–14, whose normalized squared pairing lies strictly between the π/3 and π/4 values. Candidate indices 0–11 are pairwise compatible and are simultaneously added as walls 36–47:

| New wall | Old cone ray | Old rays cut |
|---:|---:|---:|
| 36 | 834 | 127 |
| 37 | 452 | 110 |
| 38 | 382 | 75 |
| 39 | 973 | 75 |
| 40 | 12 | 58 |
| 41 | 161 | 40 |
| 42 | 232 | 39 |
| 43 | 14 | 37 |
| 44 | 2 | 30 |
| 45 | 233 | 25 |
| 46 | 313 | 24 |
| 47 | 1 | 18 |

The exact union of those cuts contains 541 old hyperideal rays; summing the table would double-count overlaps. Wall 40 cuts the saved 36-wall escape. Walls 43 and 47 meet recent walls 35 and 34 at right angles. Walls 36–39 have substantially richer old finite-neighbor lists than the preceding isolated caps, and their mutual finite contacts are also retained.

Every one of the 48 proposed walls is an actual facet. For a normal n with q=n²>0 and h=(H,n)<0, the stored vector qH−hn is future timelike, lies on n, and is strictly inside every other inequality. The signs follow directly from the verified common center and acuteness; the producer also checks its norm and zero pairing. All 100 old finite vertices are strictly retained.

## Complete cone and finite graph

The [complete-cone program](../scripts/fresh_nonarithmetic_joint_regions_20261016_cone.py) and [census](fresh_nonarithmetic_joint_regions_20261016_cone.json) perform twelve exhaustive halfspace clips. A negative/positive pair creates a new extreme ray precisely when its complete common-zero face consists of exactly those two rays. This tests actual two-face adjacency, rather than assuming a shared eight-wall mask suffices.

The exact new ray is a positive combination of its two parents, followed by a positive rational power-of-two rescaling. A directed acyclic expression graph avoids coefficient growth. Outward rational intervals propagate all strict pairing and norm signs; ambiguous pairings fall back to the exact formal field. There are 137 exact pairing fallbacks and no unresolved norm signs. Rays exactly on the incoming wall are retained with their enlarged active masks. Every stage preserves pointedness and full dimension, and the usual exhaustive halfspace-clipping induction supplies completeness.

The final counts are 172 finite, zero ideal, 2,109 hyperideal and zero opposite-sheet causal extreme rays. The [incidence checker](../scripts/fresh_nonarithmetic_joint_regions_20261016_incidence.py) and [edge certificate](fresh_nonarithmetic_joint_regions_20261016_incidence.json) verify all 395 finite edges against the complete cone: each eight-wall common-zero face contains exactly its two finite endpoint rays. All finite vertices have nine incident walls. Thus

\[
9\cdot172-2\cdot395=758
\]

finite-corner germs remain open. Of the new edges, 35 join an old vertex to a new vertex and 94 join two new vertices. The batch closes 35 old germs, but 333 old germs and 425 new-vertex germs remain open. The increase from 1,068 to 2,109 exterior rays and from 368 to 758 open finite germs is a concrete limitation of this polar step.

## A common strict escape

Take the 36-wall exterior ray at index 16, with active set

\[
[0,1,2,4,5,6,8,17,27],
\]

and normalize it by (V,H)=−1. The exact point E=V+H/8192 is spacelike, future oriented relative to H, and strict against all 36 old walls and all fifteen initial polar candidates. The additional exact π/10 candidate is also strictly negative at E. Hence every coherent subset of the complete 16-polar bank, including the 48-facet batch, retains an open patch of the future ideal boundary and has infinite volume.

This is a one-generation statement about polar normals of the fixed old cone. Further polar directions can arise after clipping, and nonpolar walls remain outside this obstruction.

## Discovery lanes and the next unresolved target

The source-bound [numerical discovery program](../scripts/fresh_nonarithmetic_joint_regions_20261016.py) has generation and saved replays for its [right-contact family](fresh_nonarithmetic_joint_regions_20261016_right.json), [finite-contact family](fresh_nonarithmetic_joint_regions_20261016_finite.json) and [polar screen](fresh_nonarithmetic_joint_regions_20261016_polar.json). The first family prescribes right contacts to [0,1,2,5,6] and all three recent walls 33–35, then one old anchor of label 2–12. Its finite-contact enlargement allows labels 2–5 on the three recent walls. All 23,698 isolated numerical branches of the enlarged bounded family fail acuteness or the retained-center orientation. This is only a bounded numerical result, not an all-label obstruction.

The [next-polar scout](../scripts/fresh_nonarithmetic_joint_regions_20261016_next_polar.py) and [saved output](fresh_nonarithmetic_joint_regions_20261016_next_polar.json) inspect all 2,109 exterior directions of the 48-wall cone at 80 decimal digits, with input midpoints from 160-bit outward enclosures. They find 46 proposals, including contacts to up to four new walls. All 46 miss the current E numerically; their normalized pairings range from approximately −29.37 to −0.03905. These are numerical proposals only: their exact acceptance and completeness have not been certified. The earlier recovered π/10 wall shows why a numerical candidate list must not be promoted to an exact exclusion.

No further noncutting batch is added here. The next useful target is a wall that actually cuts E, including mixed or nonpolar contacts and candidates without new finite vertices. Finite-corner counts are secondary to eliminating the remaining escapes.

## Replay and status

```sh
python3 scripts/fresh_nonarithmetic_joint_regions_20261016_batch.py --verify
python3 scripts/fresh_nonarithmetic_joint_regions_20261016_polar_gate.py --verify
python3 scripts/fresh_nonarithmetic_joint_regions_20261016_cone.py --verify
python3 scripts/fresh_nonarithmetic_joint_regions_20261016_incidence.py --verify
python3 scripts/fresh_nonarithmetic_joint_regions_20261016_next_polar.py --verify
```

All five generation runs and fresh saved replays pass, as do the three original numerical discovery modes. These new producers share explicitly bound frozen arithmetic and are not independent reviews of one another. The complete starting cone and the preceding individual caps have independent reviews. New independent batch/cone review is pending. Source and data are final; this report records the exact and numerical scopes separately. No frozen predecessor was changed.
