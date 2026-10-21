# H20 coupled polar repairs: exact98 walls, a new deficient ideal link

Two coupled local repairs extend the certified 96-wall H20 chamber to **98 exactly compatible walls** and cut two strict spacelike targets. However, the first new wall also creates a forced ideal ray whose horospherical normals have rank **17 rather than 19**, so that link is unbounded. A separate exact strict spacelike escape survives all 98 walls. This remains an **infinite-volume chamber**, not a new finite-volume example.

The [exact98 archive](../data/enumerations/fresh_h20_gap_next_20261013_exact.json) retains the entire [96-wall checkpoint](fresh_h20_edge_family_completion_20261013.md). Its dependency is the frozen [96-wall certificate](../data/enumerations/fresh_h20_edge_family_completion_20261013_cusp_repairs.json), including the form of signature (20,1), all previous pair gates, strict center, and earlier geometry. Numerical proposal files supply only active-set indices; the new polar vectors, ranks, normals, target cuts, and final escape are reconstructed exactly. A separate [source review and fresh exact replay](fresh_h20_gap_next_20261013_review.md) pass, including the 96-wall dependency.

| New wall | Exact norm | Relations to all earlier walls |
|---|---:|---|
| 96 | 2 | 17 right, 3 pi/3, 1 parallel, 75 dotted |
| 97 | 1 | 15 right, 6 pi/4, 76 dotted |

All 193 new pairs pass exactly, both walls are strictly negative at the saved timelike center, and both are strictly negative at all 107 inherited ideal rays. As with the earlier acute chamber, projection of the strict center onto any wall proves that it is a true facet. No new finite-vertex census is claimed.

## First coupled repair

Let p1 be the exact polar of the active walls surrounding the saved 96-wall escape, normalized by

```text
p1 squared = 3,
a = r3, b = r42/3,
a squared = b squared = 2,
(a,b) = 0,
(p1,a) = (p1,b) = -2.
```

The coefficient matrix has exact rank20, and p1 is a future spacelike point satisfying every old half-space inequality. Its uncorrected polar mirror has two forbidden angle values, both squared cosine 2/3. Correcting one of those pairings alone does not address the other constraints coherently.

Set

```text
t = (4+sqrt10)/3,
r96 = a+b+t p1.
```

Then

```text
r96 squared = 4-8t+3t squared = 2.
```

The complete exact gate gives pi/3 at old walls 15, 17, 21; parallelism with wall18; right angles or dotted relations elsewhere. The normal lies in Q(sqrt6,sqrt10). It cuts the **previously saved strict 96-wall escape**, rather than merely a numerical approximation to a boundary ray. This simultaneous correction is the substantive change from the earlier separately prescribed anchor profiles.

## Second coupled repair

After the first cut, another exact rank20 active set gives a future spacelike polar p2 with

```text
p2 squared = 1,
(a1,a2,a3) = (r5,r14,r36/3),
(ai,aj) = 2 if i=j, and 0 otherwise,
(p2,ai) = -1 for all three leaves.
```

Define

```text
r97 = 5 p2+a1+a2+a3.
```

Its norm is exactly `25+6-30 = 1`. All 97 preceding pairings pass; the six nonright finite angles are pi/4. Its coordinates already lie in Q(sqrt6). The certificate perturbs p2 a small exact amount toward the strict center, verifies that the resulting point is strictly feasible and spacelike for all 97 previous walls, and verifies that r97 cuts it.

Thus both additions make certified geometric cuts. They do not merely pass a formal Gram test.

## The new deficient ideal ray and the remaining escape

Since r18 and r96 have norm2 and pairing -2, their positive sum

```text
delta = r18+r96
```

is a nonzero future null ray. Every wall pairing with delta is nonpositive. Apart from the parallel pair itself, its complete incident set is

```text
[0,2,4,5,7,8,10,11,12,13,14,16,19,20,26,36].
```

These 16 normals have exact Gram rank16 and are orthogonal to the parallel pair. The pair contributes rank1, giving quotient rank17 in a 19-dimensional horosphere. There is consequently a nonzero common orthogonal horospherical direction, so the link is unbounded. Wall97 is dotted to wall96 and is not incident to this ray. The [separate diagnostic](../data/enumerations/fresh_h20_gap_next_20261013_cusp_diagnostic.json) freezes this argument without changing the exact98 input.

The complete incidence assertion uses acuteness: a sum of two nonpositive pairings vanishes exactly when both vanish. The two bounded cusp links completed at the 96-wall stage are not a certificate that every later-created ideal ray has a bounded link.

Independently of that diagnosis, the final archive reconstructs another exact rank20 kernel vector from the final numerical active set and perturbs it toward the center. The resulting strict escape satisfies

```text
(x,x) = 698631401/1887436800 > 0,
(x,h) < 0,
(x,r_i) < 0 for every i=0,...,97.
```

Thus an open patch of the ideal boundary remains accessible, proving infinite volume. No broader impossibility result for completions or alternative H20 chambers follows.

## Bounded alternatives explored

The initial changed search fixed 18 active short-wall equations and only **one** finite anchor, leaving a continuous norm-conic parameter. All other variable relations were required to be strictly dotted. The [numerical conic scout](../data/enumerations/fresh_h20_gap_next_20261013.json) tested 324 profiles with labels 2,3,4,6, found 134 numerical open-interval families, and sampled 234 cap entries. Twelve null-residual profiles were omitted. A [common numerical escape](../data/enumerations/fresh_h20_gap_next_20261013_sample_bank.json) survived the entire sample bank. This is numerical evidence, not an exact exclusion of the continuous families.

Two local wall-replacement tests were also negative numerically: dropping the two obstructing walls and inserting the polar, and projecting those two walls to the polar hyperplane. Projecting the polar separately off each obstruction produced compatible normals but did not cut the saved strict escape. These failures motivated the coupled formulas above; none is promoted to a global obstruction.

Finally, a bounded three-coefficient search considered asymmetric alternatives to the second cap while retaining the first cap and the original96 walls. From the prescribed labels 3,4,5,6,8,10,12 and a parallel option, [22 candidates](../data/enumerations/fresh_h20_gap_next_20261013_asymmetric.json) passed the numerical full pair gate and cut the next target. Each [alternative 98-wall chamber](../data/enumerations/fresh_h20_gap_next_20261013_asymmetric_volume.json) had a numerical spacelike escape within its 100-objective budget. These alternatives were not exactified and their count is not a mirror-distinctness claim.

The exact outcome of this checkpoint is therefore the two explicit coupled repair formulas, their genuine target cuts, and a verified explanation of why the retained98 chamber still fails finite volume.

## Replay

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=src:. python scripts/fresh_h20_gap_next_20261013_exact.py --verify
OPENBLAS_NUM_THREADS=1 PYTHONPATH=src:. python scripts/fresh_h20_gap_next_20261013_cusp_diagnostic.py --verify
```

Both commands recompute their exact assertions. They do not rerun numerical LP or rely on numerical signs. The loader imports numerical modules, but the proof uses rational/Q(sqrt6) arithmetic and exact quadratic extensions. The 98-wall certificate and diagnostic are frozen separately. The earlier96 source must retain its recorded hash.
