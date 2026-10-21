# H20 old-face extensions and two completed cusp links

The retained 88-wall H20 chamber has been extended to **96 exactly compatible walls**. Six successive caps cut certified strict spacelike targets; four finite vertices created along that path survive all 96 walls. Two further explicit caps complete the deficient links of two new ideal rays. Nevertheless, an exact strict future spacelike point survives every wall, so this chamber has **infinite volume**. No finite-volume example in a new dimension was found.

This changes the earlier direction family: the new caps use 18, 17, or 16 orthogonal short walls, with additional nonright or right anchor equations. The previous [all-real-height obstruction for 8,346 finite-vertex directions](fresh_h20_direction_family_audit_20261013.md) does not cover these old-face families. The present numerical searches are bounded; their rejection is not an unrestricted H20 obstruction.

| Exact checkpoint | Added walls | Certified result |
|---|---:|---|
| 93 walls | Five beyond the previous 88 | All 450 new pairs pass; five strict targets are cut; three PD20 finite vertices; exact strict escape |
| 94 walls | One from a 16-short support | All 93 new pairs pass; another strict target and PD20 finite vertex; exact strict escape |
| 96 walls | Two affine cusp completions | All 189 new pairs pass; both new cusp links bounded; all four stage vertices retained; same strict escape |

All eight added walls are strictly negative at the common timelike center and at the 107 inherited ideal rays. Each wall is an actual facet: project the strict center orthogonally to its positive-norm normal. The resulting point is timelike, lies on that wall, and remains strictly inside every other half-space because all off-diagonal pairings are nonpositive. Thus the complete pair gates apply to genuine facets, rather than merely to an unpruned list of inequalities.

## Exact old-face caps

The [five-cap archive](../data/enumerations/fresh_h20_edge_family_completion_20261013_exact.json) freezes the first five choices, their quadratic fields over Q(sqrt6), explicit normals, all previous pairings, strict target cuts, inherited-cusp signs, and the actual finite-vertex certificates where available. Wall indices are zero based.

| Wall | Norm | Orthogonal short support | All earlier pair types |
|---|---:|---:|---|
| 88 | 2 | 18 | 19 right, 1 pi/3, 68 dotted |
| 89 | 2 | 18 | 19 right, 1 pi/3, 69 dotted |
| 90 | 1 | 18 | 18 right, 1 pi/3, 1 pi/4, 70 dotted |
| 91 | 1 | 18 | 20 right, 2 pi/6, 1 parallel, 68 dotted |
| 92 | 1 | 17 | 18 right, 2 pi/4, 1 pi/6, 71 dotted |
| 93 | 1 | 16 | 19 right, 1 pi/4, 1 parallel, 72 dotted |

The [sixth cap](../data/enumerations/fresh_h20_edge_family_completion_20261013_four_exact.json), wall 93, simplifies to the rational vector

```text
(13,4,-13/2,11/2,13/4,13/4,13/2,39/4,39/4,
 9/2,11/4,3/4,-9/2,-13/4,-13/4,6,13/2,13/2,5,3/2,3).
```

It has norm 1 in the archived 21-dimensional Lorentzian form. Its incident finite vertex uses the 16 prescribed short walls, old walls 24, 30, 38, and the new wall. The vertices first created at walls 88, 89, 90, and 93 are **all retained at the final 96-wall stage**, checked against every final inequality. No complete vertex census is claimed.

The [standalone independent audit](fresh_h20_edge_family_completion_20261013_independent.md) uses separate Fraction arithmetic, replays all geometry, and rejects forged normal, norm, and cut data. Its saved [93-wall result](fresh_h20_edge_family_completion_20261013_independent.json) and [94-wall result](fresh_h20_edge_family_completion_20261013_four_independent.json) both pass.

## Completing two deficient new cusps

Parallel pairs create the exact future null rays

```text
delta1 = r91 + r41/3,
delta2 = r93 + r27/3.
```

The [96-wall archive](../data/enumerations/fresh_h20_edge_family_completion_20261013_cusp_repairs.json) checks their full incidence sets before and after the two additions. Their incident Gram rank is already 19 before repair, but that alone does not make a link bounded: the first link has a spherical rank-1 factor, and the second a spherical rank-18 factor.

For the first ray, add

```text
r94 = -r9 + 2 delta1.
```

This has norm 2, is parallel to wall 9, right to 21 earlier walls, and dotted to the other 72. It completes the spherical A1 factor. The final link has four affine components of ranks 7, 10, 1, 1.

For the second ray, order the spherical component by

```text
I = [0,8,17,5,23,14,24,10,1,9,18,3,15,38,22,7,16,30].
```

Let beta_i = r_i for short walls, and beta_i = r_i/3 for indices at least 24 in this list. These beta_i all have norm 2. Define theta with coefficients

```text
[2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,1,1]
```

in that order, and add

```text
r95 = -theta + 5 delta2.
```

The certificate checks theta squared = 2, its 18 pairings are 17 zeros and one 1, and r95 squared = 2. The new wall is right to 19 earlier walls, at pi/3 to one, and dotted to the other 75, including r94. The final link has affine components of ranks 18 and 1.

For every final component the replay gives a positive-definite principal minor, a strictly positive Gram-kernel vector, and the exact identity that the weighted sum of its normals is a positive multiple of the cusp ray. Thus every component becomes a positive relation in the horospherical quotient. Those quotient normals span all 19 dimensions. A recession vector pairing nonpositively with every normal must, by the positive relations, pair zero with all of them, hence vanish. The link is therefore bounded. The strict timelike center supplies a nonempty link interior. This establishes two completed cusps without assuming that rank alone implies boundedness.

An [independent source review and fresh saved replay](fresh_h20_edge_family_completion_20261013_cusp_review.md) also pass. This review validates the cusp-link argument and preservation checks; the separate arithmetic audits above certify the earlier stage vertices.

The same saved strict escape survives the repairs, with

```text
(x,x) = 210646311442057377 / 128000000000000000000 > 0,
(x,h) = -60151/4096 < 0,
(x,r_i) < 0 for every i = 0,...,95.
```

These are exact identities and strict inequalities. A strict spacelike escape, together with the strict timelike center, gives an open ideal boundary patch, so the chamber remains infinite volume. Completing these two cusp links does not certify all ends of the chamber.

## Bounded numerical bank search

The numerical search began with 1,593 earlier finite-vertex-direction cap entries plus 76 entries from the changed 18-short family. It repeatedly targeted an escape of the intersection of **all** bank half-spaces, using additional anchor equations. Cap counts below are archived entries, not claims of globally distinct mirrors.

| Numerical bank | Entries | Sampled compatible selections | Wall counts |
|---|---:|---:|---:|
| Initial | 1,669 | 24 | 326 |
| 1 | 1,742 | 4 | 327 |
| 2 | 1,815 | 4 | 328 |
| 3 | 1,890 | 4 | 329 |
| 4 | 1,915 | 4 | 330 |
| 5 | 2,083 | 8 | 333–334 |
| 6 | 2,142 | 4 | 333–334 |
| 7 | 2,400 | 4 | 339–340 |
| 8 | 2,840 | 8 | 344–345 |

Every whole-bank intersection and every listed greedy compatible selection had a numerical spacelike escape. These are floating-point scouts, not exact universal exclusions of the banks. The exact 93-, 94-, and 96-wall constructions above are small coherent selections certified separately.

The changed families gave 73, 73, 75, and 25 successive 18-short entries; then 168 17-short entries; then 59 18-short and 258 17-short entries. At the next target the prescribed bounded 18/17-short scouts found no cutters. A 16-short support with three active long anchors and one further anchor gave 440 entries from 42,500 profiles. These searches impose specific finite labels or parallel relations at their anchors; they do not cover continuous all-dotted anchor choices or every old face.

Reproducible numerical sources are the [initial bank](../data/enumerations/fresh_h20_edge_family_completion_20261013_bank.json), [face-walk summary](../data/enumerations/fresh_h20_edge_family_completion_20261013_face_walk.json), [four-anchor scout](../data/enumerations/fresh_h20_edge_family_completion_20261013_four_angles.json), and [final bank](../data/enumerations/fresh_h20_edge_family_completion_20261013_bank8.json). Their source hashes and source-index mappings preserve proposal provenance. The remaining strict escape requires a genuinely different wall or a changed retained chamber; growing these particular banks has not closed the volume.

## Replay

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=src:. python scripts/fresh_h20_edge_family_completion_20261013_exact.py --verify
OPENBLAS_NUM_THREADS=1 PYTHONPATH=src:. python scripts/fresh_h20_edge_family_completion_20261013_volume.py --verify
OPENBLAS_NUM_THREADS=1 PYTHONPATH=src:. python scripts/fresh_h20_edge_family_completion_20261013_four_exact.py --verify
OPENBLAS_NUM_THREADS=1 PYTHONPATH=src:. python scripts/fresh_h20_edge_family_completion_20261013_four_volume.py --verify
OPENBLAS_NUM_THREADS=1 PYTHONPATH=src:. python scripts/fresh_h20_edge_family_completion_20261013_cusp_repairs.py --verify
PYTHONPATH=src:. python scripts/fresh_h20_edge_family_completion_20261013_independent.py --verify
PYTHONPATH=src:. python scripts/fresh_h20_edge_family_completion_20261013_four_independent.py --verify
```

The exact cusp-completion replay uses no numerical search or LP. Its loader imports numerical modules but all certificate assertions use exact rational and radical arithmetic. The earlier strict-escape replays also check their saved points directly rather than rerunning LP.
