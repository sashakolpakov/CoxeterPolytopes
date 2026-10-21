# Four cusp-edge families: complete finite profiles cannot cut the escape

For four specified open edges of the frozen 95-wall chamber, **every acute cap that supplies a finite Coxeter endpoint strictly misses the saved global escape**. The exact check covers all 418 finite contact profiles, all real positions of the cap, and fifteen unrestricted `I₂(m)` families with `m ≥ 7`. It does not require a cutoff on the finite contact labels.

This is a restriction on four edge families, not a completion obstruction for the whole chamber. Other edges, ideal endpoint profiles, changes to retained walls, and useful finite endpoints that do not cut this particular escape remain outside the conclusion. No new wall was accepted or added in this experiment.

The four edges are taken from the [complete repaired-cusp census](fresh_h20_e7_repaired_cusp_walk_20261014.md), independently checked in its [review](fresh_h20_e7_repaired_cusp_walk_20261014_review.md). In the six affine components, their omitted walls are:

| Cusp-germ index | Omitted walls | Finite profiles | Unrestricted I₂ tails |
|---|---|---:|---:|
| 1675 | 89, 88, 18, 91, 79, 94 | 115 | 4 |
| 1673 | 89, 88, 18, 91, 17, 94 | 115 | 4 |
| 1677 | 89, 88, 18, 91, 83, 94 | 73 | 3 |
| 1723 | 89, 88, 90, 91, 79, 94 | 115 | 4 |

Their nineteen retained roots have positive-definite Gram matrices and square two. Each opposite null limit is strict on every wall outside that support. These four supports were selected numerically for their comparatively large pairing of the opposite null direction with the escape. The selection is heuristic; the exclusions below are exact.

## Complete finite contact profiles

Normalize a proposed cap to square two, and write `p` for its projection onto the nineteen retained roots. At a finite endpoint their twenty-root Gram matrix must be positive definite. Consequently the new root has only finite Coxeter contacts with the support, and its residual square is

`q = 2 − p² > 0`.

Finite Coxeter diagrams are forests, so nonzero contacts must lie in different connected support components. Every normalized inverse diagonal is at least `1/2`. Since each nonright contact has magnitude at least one, at most three are possible; three contacts must all have label 3. The checker therefore exhausts one-, two-, and three-contact choices and applies the exact Schur condition `p² < 2`.

Labels 3, 4, 5, and 6 are checked in exact `Q(√5)` Schur arithmetic. This includes the exceptional finite diagrams whenever their extensions are positive definite; there is no crystallographic restriction. A finite connected diagram containing a label at least 7 is only `I₂(m)`. Such a contact must therefore join the new root to one isolated old `A₁`, with every other contact right. Those fifteen families are retained symbolically and treated for every `m ≥ 7`.

The exact checker regenerates this complete census. It does not infer completeness from the bounded numerical sample.

## All real endpoint positions

Let `ν` be the proper cusp, `H` the saved strict center, and `W = H − projection_S(H)` for the selected support. Set

`h = W·ν < 0`,

`z = W − W²ν/(2h)`.

Here `z` is the opposite future null direction. The cusp has full-support positive affine relations. Thus any acute spacelike cap has negative pairing with `ν`: equality would force orthogonality to all cusp roots, leaving only a null multiple of `ν`. This proves the required orientation without assuming a favorable sheet.

Every cap with the fixed projection `p` is consequently

`r(a) = p + aW + (q − a²W²)ν/(2ah)`, with `a > 0`.

This has square two. Its intersection with the old edge is

`P(a) = z − qν/(2a²h)`,

with square `−q/a² < 0`. The inherited open-edge certificate implies that this point is future and strict on every old wall outside the nineteen-root support. Hence the parameterization concerns actual finite endpoints, rather than merely pair-compatible normal vectors.

For any tested wall or the escape `v`, the pairing can be written

`r(a)·v = p·v + a(z·v) + q(ν·v)/(2ha)`.

For every wall outside the support, `z·v < 0` and `ν·v ≤ 0`. The saved escape `E` also has `z·E < 0` and `ν·E < 0`, checked exactly. These pairings are strictly decreasing functions of `a`.

For each of the 418 finite profiles the certificate supplies a positive rational separator `a₀` and a retained wall `n` such that

`r(a₀)·n > 0`, while `r(a₀)·E < 0`.

All smaller parameters violate the old-wall inequality; all larger parameters strictly miss the escape. This proves the entire real-parameter exclusion for that profile. The checker evaluates the pairings after multiplication by the negative factor `2a₀h`, avoiding any field division. Eight certificates use bridge wall 87; its additional quadratic radical is handled exactly as well.

## The unrestricted I₂ tails

For an isolated support root, let `u` be its fundamental dual, so `u² = 1/2`. Put

`c = 2cos(π/m)`, `p = −cu`, `q = 2 − c²/2`.

All `m ≥ 7` lie in the larger closed interval `√3 ≤ c ≤ 2`. At a fixed rational separator `a₀`, multiplication by `2a₀h` gives a quadratic

`K_v(c) = A_v c² + B_v c + C_v`,

where

`A_v = −(ν·v)/2`,

`B_v = −2a₀h(u·v)`,

`C_v = 2a₀²h(W·v) + (2 − a₀²W²)(ν·v)`.

For the blocking wall, this quadratic is convex. Strict negativity at both interval endpoints makes it negative throughout, which means the original wall pairing is positive. For the escape, strict positivity is established at its exact minimum using endpoint derivative signs, with an interior-minimum discriminant test available. In the saved certificates eleven minima occur at the left endpoint and four at the right. All fifteen tails close with one rational separator each; no label sampling or interval subdivision is used in this proof. Including `c = 2`, where the finite Schur complement degenerates, is only a conservative endpoint test.

## Numerical scout and exact certificate

The [bounded scout](../scripts/fresh_h20_e7_repaired_cusp_edgewalk_20261014_scout.py) and its [JSON](fresh_h20_e7_repaired_cusp_edgewalk_20261014_scout.json) tested the 418 finite profiles and 120 sampled tail labels. All 538 failed its numerical continuous acuteness/cut interval screen, so no numerical finite-anchor trial was reached. Those numerical failures are separate from the exact theorem.

The [exact checker](../scripts/fresh_h20_e7_repaired_cusp_edgewalk_20261014_exact.py) and [certificate](fresh_h20_e7_repaired_cusp_edgewalk_20261014_exact.json) regenerate the finite-profile census, check every rational separator, and prove the fifteen full tail inequalities. Numerical data only propose separator constants and blocking walls. Acceptance uses exact Fraction arithmetic in the multiquadratic field with primes `2,79,11,3,5`, together with the frozen bridge radical.

The exact producer and saved-result replay passed. The replay command is:

```sh
PYTHONPATH=src:. python scripts/fresh_h20_e7_repaired_cusp_edgewalk_20261014_exact.py --verify
```

A finite compatible collection drawn solely from these four finite-endpoint families therefore cannot remove the current strict spacelike escape. This does not rule out such caps as intermediate steps, or walls outside these four families.
