# Joint cap completion of the released H8 seed

This search did **not** produce a finite-volume H9 polytope. It produced **1,756 distinct new mirrors**, each exactly compatible with most of one actual 32-facet configuration and cutting several of its open ends. An exact surviving escape rules out adding this entire new bank to that configuration. Six deliberate replacements of obstructing old caps also retain exact escapes, even after allowing every eligible new mirror. These conclusions concern the stated banks and retained cores, not arbitrary completions of the seed.

The preceding certificates in [the cap-completion report](fresh_h8_cap_completion_20261013.md) and its independent thirteen-witness cutter audit were replayed successfully before this search. The preceding bank obstruction remains compact-only when all old caps are optional; parallel pairs were not silently treated as a finite-volume obstruction here.

## Construction directed by one configuration

The seed is the ten-wall, rank-ten Lorentzian seed obtained by changing the Bugaenko H8 endpoint label from 5 to 3 and omitting original wall 6. Its Gram matrix and strict timelike center are unchanged from [the seed archive](../data/enumerations/fresh_h8_alternate_3_6_20261013.json). All ten seed normals have square 2. Write \(K=\mathbf Q(\phi)\), \(\phi^2=\phi+1\).

The starting configuration consists of these ten seed walls and the following 22 old cap nodes:

```
31, 51, 67, 89, 112, 156, 197, 310, 348, 376, 456,
469, 496, 643, 666, 683, 717, 837, 862, 1028, 1038, 1067
```

It is case 0 in [the preceding extension archive](../data/enumerations/fresh_h8_cap_completion_20261013_extensions.json). The new independent replay verifies the seed signature \((9,1)\), its inverse, every old normal's square and orientation, and all 231 cap-pair relations in this starting configuration.

Farthest-point sampling of 3,922 numerical affine-section vertices proposed 24 separated gaps. Every saved gap was then reconstructed as an **exact strict future spacelike witness**: all ten seed coordinates and all 22 cap inequalities are strictly positive, its Lorentzian square is positive, and its pairing with the common timelike center is negative. Numerical vertex enumeration is used only to propose these witnesses; no completeness or compactness conclusion relies on it. The [target archive](../data/enumerations/fresh_h9_joint_completion_20261013_targets.json) contains all coordinates, norms, center pairings, and complete old-1,084-wall cutter lists.

For each new cap, three seed-pairing magnitudes are free and required to exceed 2. At most three other magnitudes are nonzero, with values \(1\) or \(\phi\). Two anchors are chosen from the nine \(K\)-rational normals already present in the configuration, including polar wall 348. Prescribed anchor pairing magnitudes are \(0,1,\phi\). The two anchor equations leave a line \(w=w_0+xw_1\) over \(K\); the square-2 equation is quadratic in \(x\). Each accepted normal is checked exactly in \(K(\sqrt D)\).

The numerical proposal score prioritizes cutting multiple saved gaps and retaining existing facets. It retains two proposals per anchor-pair/retained-facet pattern, so this is **not an exhaustive census of the parameter family**. Accepted proposals must cut at least two of the 24 gaps and remain compact-compatible with at least 18 of the 22 old caps. Exact pair checking covers every possible finite label in the degree-at-most-eight compositum:

```
2, 3, 4, 5, 6, 8, 10, 12, 15, 16, 17, 20, 24, 30, 32, 34, 40, 48, 60
```

There are 1,781 accepted entries: 1,756 distinct mirrors outside the old 1,084-wall dictionary, 16 further copies of new mirrors, and nine mirrors already present in the old dictionary. Equality is tested exactly using all ten normalized squared seed pairings. Nonpositive seed pairings and the independent seed basis make this equivalent to equality of oriented mirrors. The exact pattern prefilter cannot discard a possible equality.

| Number of retained old caps | Accepted entries |
|---:|---:|
| 22 | 40 |
| 21 | 259 |
| 20 | 549 |
| 19 | 587 |
| 18 | 346 |

Each entry cuts between two and eight target gaps. The [cap archive](../data/enumerations/fresh_h9_joint_completion_20261013.json) records the line, norm polynomial, discriminant, exact root, both anchors, every retained/dropped old-cap relation, and the exact cut mask. The [mirror archive](../data/enumerations/fresh_h9_joint_completion_20261013_mirrors.json) records the 1,756 representatives and all aliases.

## Deliberate replacements and exact surviving gaps

Target 8 survives **all 1,781 accepted cap entries**. Consequently no subset of the starting configuration plus this new bank can have finite volume, independently of mutual compatibility among new caps. A spacelike feasible vector with a common strict timelike center supplies an open spacelike region after a small perturbation toward the center, hence an open set of ideal directions and infinite hyperbolic volume.

To repair this particular gap using old walls, its complete old-dictionary cutter list is

```
0, 1, 5, 8, 10, 18, 44, 1035, 1036, 1037, 1039.
```

Six choices require removing at most six of the 22 old caps. For each such choice, the search retains **every old cap compatible with the new mandatory cutter**, inserts that cutter, and allows every one of the 1,756 novel mirrors compatible with the resulting fixed core. Exact coverage checks already fail before mutual compatibility among the eligible new caps needs to be imposed.

| Mandatory cutter | Removed old cap nodes | Eligible novel caps | Surviving target indices |
|---:|:---|---:|:---|
| 1039 | 1038 | 39 | 3, 5, 13, 22 |
| 1036 | 67, 1038 | 41 | 3, 5, 22 |
| 1037 | 67, 1038 | 41 | 3, 5, 22 |
| 1035 | 67, 1038, 1067 | 47 | 3, 5, 22 |
| 18 | 348, 862, 1028, 1038, 1067 | 114 | 5, 22 |
| 44 | 67, 112, 348, 1028, 1038, 1067 | 105 | 3, 5, 22 |

The nine old-mirror aliases among the accepted entries also fail to cut targets 5 and 22. Thus allowing those aliases does not repair these six cases. All coordinates and complete eligibility/coverage lists are in the [replacement archive](../data/enumerations/fresh_h9_joint_completion_20261013_replacements.json), with the eleven target-8 cutter compatibility masks in the [gap archive](../data/enumerations/fresh_h9_joint_completion_20261013_gap8.json).

Each listed fixed core is pairwise acute and contains a common strict timelike interior point; every eligible cap is acute to that core. Hence its fixed walls remain actual facets by the Lorentz-orthogonal projection argument. A relation rejected as parallel here is a **compact** incompatibility. These tables do not exclude extensions which retain parallel walls, remove additional core walls, use the other five old cutters, or introduce other mirrors. In particular, there is no claim that the combined old and new dictionary is globally insufficient.

## Direct quartic-anchor probe

The cheapest target-8 repair replaces old cap 1038 by quartic cap 1039. A further probe uses that repaired cap directly as an anchor. If its normal is \(A+\sqrt d B\), the exact decomposition satisfies

\[
A^2=\tfrac72+4\phi,
\qquad B^2=-\tfrac94-\tfrac{29}8\phi,
\qquad (A,B)=0.
\]

To keep the new norm equation quadratic over \(K\), this probe imposes the sufficient conditions \(A\cdot w=M\), \(B\cdot w=0\). With an additional rational-anchor equation, four free seed distances leave a line; without the additional anchor, three free distances leave a line. **These split equations are restrictive, not necessary for compatibility with the quartic anchor.** In particular, \(B\) is timelike, so the imposed orthogonality is a substantial extra condition.

The four-free probe considered all 210 free index sets, the same fixed-support bound, and anchor magnitudes \(M\in\{0,1,\phi,2\phi,1+2\phi,2+\phi\}\), with the second anchor magnitude in \(\{0,1,\phi\}\). It produced no numerical roots passing the seed-distance and center filters. The three-free probe considered all 120 free index sets and produced 39 numerical roots passing those filters, but none passed the requirement to cut at least two saved gaps including either target 5 or target 22. No exact accepted cap was generated in these probes. Floating-point degeneracies and omitted parameter choices are not excluded.

The separate [four-free finite-angle](../data/enumerations/fresh_h9_joint_completion_20261013_quartic_anchor.json), [four-free dotted-angle](../data/enumerations/fresh_h9_joint_completion_20261013_quartic_dotted.json), [three-free two-gap](../data/enumerations/fresh_h9_joint_completion_20261013_quartic_only.json), and [three-free either-gap](../data/enumerations/fresh_h9_joint_completion_20261013_quartic_either.json) archives preserve these bounded numerical outcomes. They are not exact absence certificates.

## Exact replay and remaining opening

The producer-free [Fraction verifier](../scripts/fresh_h9_joint_completion_20261013_verify.py) imports no search producer, numerical library, optimizer, or SAT solver. It checks the 24 exact witnesses, 26,016 old-dictionary witness signs, all 1,781 cap equations and seed angles, 39,182 retained/dropped pair classifications, 42,744 new-cap witness signs, mirror identification, the common surviving target 8, and all six replacement coverage calculations. Its [audit output](../data/enumerations/fresh_h9_joint_completion_20261013_verify.json) records these checks. The normalized-square mirror-identification argument also received an independent source-level review.

Run from the repository root:

```sh
PYTHONPATH=. python scripts/fresh_h9_joint_completion_20261013_verify.py --verify
```

The generator is [fresh_h9_joint_completion_20261013.py](../scripts/fresh_h9_joint_completion_20261013.py); the coherent replacement calculation is [fresh_h9_joint_completion_20261013_replacements.py](../scripts/fresh_h9_joint_completion_20261013_replacements.py); the restricted quartic-anchor probe is [fresh_h9_joint_completion_20261013_quartic_anchor.py](../scripts/fresh_h9_joint_completion_20261013_quartic_anchor.py).

No new configuration reached a volume-certification stage: each specified extension scheme already has an exact escape before a compatible subset can close it. A remaining constructive direction is to solve the quartic-anchor equation without coefficient splitting, or to replace a larger coherent group of old caps while explicitly cutting targets 5, 8, and 22. The present certificates do not settle either direction.
