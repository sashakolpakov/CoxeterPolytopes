# Nonright spherical endpoint search on the 860-wall H22 chamber

No new compatible cap was found. This search changes the endpoint angle profile beyond the previously excluded all-right families. The abstract spherical profile census is exact; the subsequent geometric search is a bounded numerical scout and is **not** a theorem excluding every nonright bridge.

The starting data are the four actual finite vertices of the [860-wall chamber](fresh_h22_multiquadratic_20261013.md). Releasing each of their 22 incident walls gives 88 positive-definite rank-21 edge supports. The desired new norm-two wall has one or more nonright contacts with that support and a finite contact with one of the 305 positive-height caps in the large contact component. The 860 existing walls remain fixed.

## Exact spherical profile census

The [exact producer and replay](../scripts/fresh_h22_nonright_bridge_20261013_profiles.py) reconstruct the incident Grams from the independently checked radical normals. All parent root norms are retained, including norm four. For prescribed edge pairings p and edge Gram G, the extended Gram is positive definite precisely when

`R = 2 - p^T G^-1 p > 0`.

Fraction matrix inversion and exact squarefree-radical signs give the following [certificate](fresh_h22_nonright_bridge_20261013_profiles.json):

| Nonright contacts | Prescribed labels | Positive-definite profiles |
| --- | --- | ---: |
| One | 3, 4, 5, 6 | 1,951 |
| Two | 3, 4, 5, 6 | 3,706 |
| Three | 3 on every contact | 823 |
| **Total explicit profiles** | | **6,480** |

There are also **217 parametric families**: attaching the new node only to an isolated old A1 with any label `m>=7`. Their Schur residual is `2 sin²(pi/m)>0`.

The finite Coxeter diagram classification makes this list exhaustive for abstract **nonright spherical extensions** of these supports. A node has degree at most three; a degree-three node has only label-three edges; a degree-two node has no label above five; and a label at least seven occurs only in an isolated rank-two dihedral component. The checked nonnegative inverse permits single-contact Schur tests to prune multiple-contact profiles. This completeness concerns the spherical Gram profiles, not the existence or exclusion of geometric caps satisfying all 860 inequalities. The all-right profile is covered separately by the [earlier bridge analysis](fresh_h22_contact_bridge_20261013.md).

## Numerical geometric search

Let s be the edge-span vector having pairings p. In the Lorentz plane perpendicular to the edge, use the future unit vertex v and normalized released-wall dual e. A norm-two candidate has

`r = s + x v + y e`, with `-x²+y²=R`.

The searched outgoing branch is

`x=sqrt(R)(t-t^-1)/2`, `y=-sqrt(R)(t+t^-1)/2`, `t>1`.

For each prescribed finite anchor angle, its linear pairing equation becomes a quadratic in t. The scout solves both quadratic branches, then checks positive height, the retained strict center, all 860 acuteness inequalities, and every resulting finite-angle relation. The proposed endpoint pairing formulas respect each anchor's actual norm.

All **305** positive-height anchors were tested with every integer anchor label **2 through 32**. In addition to the 6,480 explicit profiles, the 217 dihedral families were sampled at labels **7, 8, 9, 10, 12, 16, 18, 24, 30, 32, 48, 64**, giving another 2,604 profiles.

| Scout | Profiles | Outgoing quadratic proposals | Proposals passing acuteness | Accepted |
| --- | ---: | ---: | ---: | ---: |
| [One contact](fresh_h22_nonright_bridge_20261013_1link_scout.json) | 1,951 | 18,446,705 | 44,543 | 0 |
| [Two contacts](fresh_h22_nonright_bridge_20261013_2link_scout.json) | 3,706 | 35,040,230 | 49,487 | 0 |
| [Three contacts](fresh_h22_nonright_bridge_20261013_3link_scout.json) | 823 | 7,781,465 | 9,907 | 0 |
| [Selected dihedral labels](fresh_h22_nonright_bridge_20261013_i2_1link_scout.json) | 2,604 | 24,620,820 | 9,961 | 0 |
| **Total** | **9,084** | **85,889,220** | **113,898** | **0** |

Counts are parameter proposals, not distinct mirrors. Every proposal passing numerical acuteness failed another finite Coxeter-angle test. No degenerate quadratic was reported in these runs. The numerical tolerances are not rigorous root-isolation bounds; the reported failure remains numerical even within the finite label list. The script retains the closest angle failures for diagnosis rather than presenting them as valid normals.

No new cap, finite vertex, or volume improvement is claimed. Untested dihedral labels, anchor labels above 32, parallel anchor contact, other edge supports, and wall replacements remain outside the search. The exact abstract profile census does not turn these finite numerical failures into an all-label obstruction.

## Replay

The exact profile certificate passes saved-result replay:

The parent's independent source review and fresh exact replay also pass; see the [scoped review](fresh_h22_nonright_bridge_20261013_review.md). That review does not upgrade the numerical geometric scout to an exact exclusion theorem.

```sh
PYTHONPATH=src:. python scripts/fresh_h22_nonright_bridge_20261013_profiles.py --verify
```

The numerical producer is [the separate scout](../scripts/fresh_h22_nonright_bridge_20261013_scout.py). Its `--links`, `--labels`, and `--anchor-labels` arguments expose the finite alphabets. For example:

```sh
PYTHONPATH=src:. OPENBLAS_NUM_THREADS=1 python scripts/fresh_h22_nonright_bridge_20261013_scout.py --links 3 --labels 3 --anchor-labels 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32
```

For the dihedral sample, import the scout, set `S.TAG += '_i2'`, and call `S.run(1,(7,8,9,10,12,16,18,24,30,32,48,64),tuple(range(2,33)))`; this keeps its archive separate from the one-contact small-label run. All four numerical archives bind the current source hashes. No frozen 860-wall or earlier all-right certificate was modified.
