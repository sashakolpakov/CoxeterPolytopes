# Golden H3 triage and two prime-7 H9 hosts

The new bounded searches produced exact 45-wall and 56-wall Coxeter cones in H9, both with certified infinite-volume escapes. They do not give a new finite-volume polytope. The positive blocks are H3+A6, and the two quadratic spaces are distinct from the explicitly compared recent E8/scalar and E7/A2 spaces. This is a statement about the listed arithmetic hosts, not worldwide novelty of reflection groups.

Write O=Z[φ], φ²=φ+1. All field values in the JSON artifacts use exact rational coefficients in the basis (1,φ).

## Why H3 cubed was triaged first

For

```
H3 = [[2,-φ,0],[-φ,2,-1],[0,-1,2]],
det(H3) = 2φ^-2,
```

an explicit matrix A satisfies Aᵀ(2I3)A=H3. Consequently the two proposed H3³ hosts are field-similar to already considered diagonal and E8-plus-scalar spaces:

```
<-φ>  + H3³  ~  <-2φ> + I9  ~  <-2φ> + E8 + <1>,
<-2φ> + H3³  ~  <-φ>  + I9  ~  <-φ>  + E8 + <1>.
```

The matrices, similarity factor 2, and exact Gram identities are saved in [the structural certificate](fresh_h3_blocks_h9_structure_20261013.json). These are field similarities; they are not asserted to be integral lattice isometries.

Nor does a previous argument using a unique integral overclosure apply automatically. The H3 discriminant has a generator of norm 3/2. The spatial H3³ discriminant is F4³, with five maximal isotropic lines, forming block-permutation orbits of sizes 2 and 3. The even host adds a fourth discriminant coordinate; it has 21 isotropic lines and five maximal isotropic planes, again in orbits of sizes 2 and 3. Those five planes share one isotropic line, but the maximal unimodular overclosures remain nonunique.

Under the displayed similarities, explicit H3 simple reflections fail to preserve the corresponding E8 or diagonal integral lattice. Thus these displayed placements do not identify the full integral reflection groups. An unidentified conjugacy remains possible. No positive-height H3³ search was run, because this triage did not establish a new full-group search lane.

The full zero-height census was nevertheless checked: after normalization by square units, the primitive candidate norm pools are {1,2,4} for the odd host and {2,4} for the even host. There are exactly three dominant zero-height roots, all of norm 2, one in each H3 block. The finite reflection orbit in each block contains 30 roots. There are no extra mixed-block zero-height mirrors in these pools.

## The discriminant-changing pivot

The actual positive-height searches use

```
Lα = <-α> + H3 + A6,   α = φ or 2φ,
```

where A6 has diagonal 2 and neighboring entries -1. Both forms have signature (9,1) in the chosen embedding and are positive definite in the conjugate embedding. Their determinant is

```
det(Lα) = -14αφ^-2.
```

The prime 7 is inert in O, and this determinant has valuation 1 at 7. The independently checked comparison list has determinants -φ, -2φ, -√5, -3φ, -(1+4φ), 1-3φ, -4φ, -6φ and -12φ, all of valuation 0 at 7. A similarity in dimension 10 changes the determinant by a square, so it cannot change this parity. This excludes field similarity with exactly those compared hosts. It does not classify all hosts elsewhere in the repository or literature. The two new H3+A6 forms also have nonsquare determinant ratio 2, so they are not field-similar to each other.

The complete primitive candidate norm alphabets are

| Time coefficient | Candidate norms, up to square units |
|---|---|
| φ | 1, 2, 4, 7, 14, 28 |
| 2φ | 2, 4, 14, 28 |

Here “candidate” does not assert that every listed norm is realized by a primitive reflection. The ring O is a PID, so each integral mirror normal can be made primitive: in fact O is norm-Euclidean, since rounding the two coefficients leaves a,b in [-1/2,1/2], where |a²+ab-b²| ≤ 5/16 < 1. To prove the norm bound, 14B^-1 is integral. If r is primitive and its reflection is integral, then 2Br/q is integral, hence (28/q)r is integral. Primitivity gives q dividing 28. The primes 2 and 7 are inert, and every totally positive unit in O is a square. These facts give the first list. The even host has all norms in 2O, giving the second list. A possible nonunit content for these lists can only be 2, and only at norms 4 or 28; the search removes a root of either norm if every coordinate is divisible by 2. This prevents retaining the same mirror twice through imprimitive parallel normals.

For every candidate norm, the zero-height census has exactly two dominant roots, both of norm 2, one in H3 and one in A6. The independent finite root orbits have sizes 30 and 42. Since the existing finite reflection group moves any vector to its dominant cone, this census establishes the complete zero-height chamber H3+A6, including the absence of extra mixed-block zero-height mirrors.

## Bounded exact search and geometry

Both searches completed the height bound h=αx²/q ≤ 30. For each norm q, the conjugate positive-definite form also bounds x'² ≤ -q'/α'. The resulting complete time-control lists have 20 and 15 entries.

The spatial enumeration uses dominant covectors d=By. Exact Cauchy bounds give d_i² ≤ 2T at both embeddings, where T=q+αx². A Schur-complement lower bound prunes partial covectors in both embeddings; the positive entries of each block's inverse Gram matrix give an additional valid bound in the chosen embedding. Final vectors are reconstructed by B^-1d and tested for integrality. The two block energies are joined exactly. Equal-height packets are handled together, with exact pairwise angle checks before accepting the packet.

| α | Walls through height 30 | Certified spacelike escape norm | Artifact |
|---|---:|---|---|
| φ | 45 | (1908-625φ)/2527 > 0 | [45-wall cone](fresh_h3_blocks_h9_a6_phi_h30_20261013.json) |
| 2φ | 56 | 468/35-8φ > 0 | [56-wall cone](fresh_h3_blocks_h9_a6_2phi_h30_20261013.json) |

Every saved normal is primitive and has an integral lattice reflection. Every wall pair passes the exact Coxeter-angle test, including the golden-field distinction between π/5 and the inadmissible angle 2π/5. Both cones have a verified strict timelike interior. Each escape is an exact spacelike extreme ray satisfying all wall inequalities. These establish infinite volume of the finite prefixes. They do not prove nonreflectivity of either full arithmetic host.

The A6 chain reversal preserves both entire wall sets. Averaging each escape with its reversed image still gives a spacelike vector satisfying every wall inequality. The resulting symmetry-fixed norms are

```
α=φ:   (258-97φ)/361 > 0,
α=2φ:  66/5-8φ > 0.
```

The original escapes lie in the three-dimensional span of time, the first H3 fundamental weight, and the third A6 fundamental weight. The averaged escapes replace the last direction by the sum of the third and fourth A6 fundamental weights. Thus imposing the visible A6 symmetry does not remove the escape. This suggests a possible targeted arithmetic slice for future work, but no global obstruction or completion within that slice is claimed here. No larger uniform height was attempted.

## Independent checks and reproduction

The [Fraction-only audit](fresh_h3_blocks_h9_independent_audit_20261013.json) replays the field similarities, discriminant census, complete zero-height census, finite root orbits, determinant distinction, all saved primitive roots, reflection integrality, chronological heights, time-control completeness, strict interiors, every Coxeter angle, both escapes and their symmetry-fixed averages. It does not import project exact arithmetic. Positive-height spatial enumeration is reviewed from its bounds and compared with independent Cartesian enumeration on bounded samples; the entire positive-height search is not independently rerun by a second algorithm.

The 17 focused regressions compare both block enumerators with independent Cartesian searches in every norm class, test a positive-height block slice, verify finite root orbit sizes, and reject forged escape norms.

```sh
PYTHONPATH=src:. python3 scripts/fresh_h3_blocks_h9_structure.py
PYTHONPATH=src:. python3 scripts/fresh_h3_blocks_h9_a6_search.py --alpha phi --height 30 --output results/fresh_h3_blocks_h9_a6_phi_h30_20261013.json
PYTHONPATH=src:. python3 scripts/fresh_h3_blocks_h9_a6_search.py --alpha 2phi --height 30 --output results/fresh_h3_blocks_h9_a6_2phi_h30_20261013.json
PYTHONPATH=. python3 scripts/fresh_h3_blocks_h9_independent_verify.py --output results/fresh_h3_blocks_h9_independent_audit_20261013.json
PYTHONPATH=src:. python3 -m pytest -q tests/test_fresh_h3_blocks_h9.py
```

The two original search runs each took approximately 22 seconds; the 17 regression tests passed in 7.71 seconds on the current workspace machine.
