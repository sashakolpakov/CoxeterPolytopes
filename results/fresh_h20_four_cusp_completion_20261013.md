# Four H20 cusp completions, changed gap profiles, and remaining defects

The four deficient cusp links of the frozen215-wall branch are now **exactly bounded** after nine mutually compatible additions. The resulting224-wall chamber still has infinite volume. Two changed-profile caps subsequently cut its successive escapes and produce a certified finite vertex, yielding226 coherent facets. **The 226 chamber is still infinite-volume**, with an exact strict spacelike escape and two newly diagnosed unbounded cusp links. No new finite-volume H20 polytope is claimed.

| Checkpoint | Constructive result | Remaining volume evidence |
|---|---|---|
|224|Nine minimum-shift highest-root caps close the four specified links; one cuts the 215 escape.|Strict escape of square 2505286225/33554432; it survives the entire nine-family all-real-shift class.|
|225|Rational norm-one cap cuts224 escape and preserves the four bounded links.|Strict escape of square 2104017489/33554432; a cone identity excludes every acute cap orthogonal to the fixed17-short block.|
|226|Norm-one cap overQ(sqrt7), two pi/4 anchors and19 right angles, cuts225 escape and creates an actual finite vertex.|Strict escape of square 1045865484809280029/10^20; two specific new cusp links are unbounded.|

## Exact dependencies and geometric interpretation

The starting215 configuration is the [previous rational branch](fresh_h20_changed_cusp_repair_20261013.md), with [frozen certificate](../data/enumerations/fresh_h20_changed_cusp_repair_20261013_second_exact.json). Each new producer recursively checks the predecessor certificate's declared file hashes before reconstructing its normals. Previously certified old–old relations are dependencies, while every new norm, center sign, and new–old/new–new relation is recomputed exactly. Arithmetic uses rational numbers, Q(sqrt6), and explicitly adjoined positive square roots.

All listed walls are spacelike, pair nonpositively, and have a common strict timelike point in the inherited signature(20,1) form. They therefore remain actual facets: projecting that point to one wall stays timelike, and the nonpositive pairings make every other wall inequality strict. Thus the exact Coxeter pair gates describe genuine acute Coxeter chambers; they do not supply finite volume. Positive norms of the archived strict future spacelike points prove infinite volume.

## Closing the four links at their sharp minimum

The 215 branch has new parallel cusps from old walls 30,34,46,50 and wall 214. Their finite components consist of six A9 diagrams and three D6 diagrams. The rank-six diagrams are **D6, not E6**: exact reflection closure gives 60 roots. Each A9 closure has 90 roots. The unique highest positive root is recovered directly from these finite orbits, after normalizing the original norm18 simple roots by1/3.

For each finite component let theta be its norm-two highest root and delta the cusp ray. The candidate is

`r(t) = -theta + t*delta`, with `r(t)^2 = 2`.

For each of the nine families, a saved retained-wall inequality gives the sharp necessary condition`t >= 2`. At`t=2` all nine caps are mutually compatible. The 1,971 new pair checks comprise 189 right, 81 pi/3, 34 parallel and 1,667 dotted relations. In particular, the 36 new–new checks comprise 15 right, 6 parallel and 15 dotted pairs. Wall 215 cuts the prior escape by2234047239/5000000000.

At each specified cusp the certificate enumerates all 224 signs. Its incident normals have actual rank 20 and a checked positive relation using **every** incident normal, equal to a positive multiple of delta. Their images span the 19-dimensional Euclidean quotient. Any recession direction has nonpositive pairing with these images; the positive relation forces all pairings to vanish, and full quotient span then forces that direction to be zero. This proves bounded links, rather than inferring boundedness from rank alone.

The new strict escape is negative both on every minimal cap and on its delta. Consequently

`(x,r(t)) = (x,r(2)) + (t-2)*(x,delta) < 0` for every`t >= 2`.

This excludes finite-volume completion by any finite collection from these **nine fixed negative-highest-root affine-ray families**, with arbitrary real shifts and no restriction on possible labels. It does not exclude other cusp profiles or facet replacement. The newly created parallel pairs are not silently assumed to have bounded links.

[Exact224 script](../scripts/fresh_h20_four_cusp_completion_20261013.py) · [certificate](../data/enumerations/fresh_h20_four_cusp_completion_20261013.json).

## Two changes of gap profile

At the 224 escape, fixing17 active short-wall orthogonalities leaves a4-dimensional Lorentzian residual space. A bounded three-anchor scout produced 2,670 numerical proposals. Its first all-right-anchor proposal reconstructs to a rational norm-one cap, wall 224. Exact validation gives20 right, one pi/4, one parallel and202 dotted relations; it cuts the 224 escape by 47969/16384 and is strict at all four completed cusps.

Its new feasible spacelike vector p has square 63 and an exact decomposition

`p = (9/31)r18 + (85/62)r49 + (73/62)r55 + (3/31)r181 + a vector in span(S)`,

where

`S = [0,1,2,3,4,5,8,9,13,14,15,20,21,22,24,25,26]`.

All coefficients in the displayed retained-wall combination are positive. Every additional normal orthogonal toS and acute-compatible with these four retained walls has nonpositive pairing with p. This rules out **all profiles in that fixed17-block**, including arbitrary finite labels and unrestricted dotted anchors. For a proposed finite full-dimensional completion, perturb p slightly toward that completion's own strict interior point; this gives a strict spacelike escape. The argument does not require our archived center to survive arbitrary additions, and it does not apply to replacing the retained facets.

The next construction changes to a different PD19 short block. An explicitly bounded scout released each of its 19 walls in turn, testing 417,555 finite-anchor profiles and saving 1,068 numerical entries, with repetitions between released blocks. These are proposal counts, not distinct-mirror or exact exhaustive counts. The chosen exact cap is

`r225 = u - sqrt(7)*z/3`, `u^2 = -4/3`, `z^2 = 3`, hence`r225^2 = 1`.

Its19 right neighbors are

`[0,1,2,3,4,5,8,9,10,13,15,16,18,19,20,21,22,23,24]`.

It has pairing−1 with norm-two walls 17 and25, giving pi/4 at both anchors, and is dotted to the remaining 204 walls. Its pairing with225's saved escape is`5*sqrt(7)-33171/8192 > 0`. The exact projection of the strict center onto its 19-wall block and the new wall is timelike and vanishes on exactly these20 walls, proving an actual finite vertex. The four previously completed cusp links remain unchanged because both new caps pair strictly negatively with their rays. No claim is made that earlier stage-specific finite vertices all survive.

[Exact225 script](../scripts/fresh_h20_four_cusp_completion_20261013_gap_exact.py) · [certificate and fixed17-block identity](../data/enumerations/fresh_h20_four_cusp_completion_20261013_gap_exact.json). [Exact226 script](../scripts/fresh_h20_four_cusp_completion_20261013_next_exact.py) · [certificate, vertex and escape](../data/enumerations/fresh_h20_four_cusp_completion_20261013_next_exact.json).

## Two explicit remaining cusp obligations

The final226 certificate is accompanied by a separate exact diagnosis, checking every wall sign at each ray:

* `delta0 = r224 + r48/3`. Its20 incident normals have normal-span rank 20 and quotient rank 19. However, its components consist of a positive affine rank-one factor and finite spherical factors of ranks 6 and 12. A saved positive-norm recession direction is nonpositive on every incident wall and strictly negative on at least one; the link is unbounded.
* `delta1 = r225 + (r17+r25)/2`. Its18 incident normals have quotient rank 17. Its positive affine factor has rank 2 and the spherical factors total 15. Another exact recession direction proves this link unbounded.

The positive affine dependencies show the null ray belongs to each incident span; since the Gram radical is exactly one-dimensional, normal-span rank equals Gram rank plus one. These are actual feasible future null rays, not merely singular subdiagrams. A strict global center gives full-dimensional horospherical link interiors. For either saved positive-norm recession vector w, delta+epsilon*w is spacelike and satisfies every incident inequality for epsilon>0. All nonincident inequalities remain strict for sufficiently small epsilon. A further small move toward the strict center makes every inequality strict. Thus the recession proves a geometric cusp deficit directly, as well as agreeing with the separate volume escape.

[Exact diagnostic](../scripts/fresh_h20_four_cusp_completion_20261013_new_cusps.py) · [two-cusp certificate](../data/enumerations/fresh_h20_four_cusp_completion_20261013_new_cusps.json).

This is a scoped diagnosis, not a full226 cusp census. Next constructive work must change or complete these new link profiles and remove the saved remote escape; increasing heights in the closed nine-family class or keeping the excluded17-block cannot do so.

## Reproduction

From the repository root:

```sh
python3 scripts/fresh_h20_four_cusp_completion_20261013.py --verify
python3 scripts/fresh_h20_four_cusp_completion_20261013_gap_exact.py --verify
python3 scripts/fresh_h20_four_cusp_completion_20261013_next_exact.py --verify
python3 scripts/fresh_h20_four_cusp_completion_20261013_new_cusps.py --verify
```

All four frozen exact replays passed. A separate [source/proof review and fresh four-command replay](fresh_h20_four_cusp_completion_20261013_review.md) also passed. Exact reconstruction is the claim-bearing layer; the [three-anchor scout](../scripts/fresh_h20_four_cusp_completion_20261013_gap_scout.py), [changed-block scout](../scripts/fresh_h20_four_cusp_completion_20261013_next_scout.py), and their numerical volume proposals are not volume proofs or exhaustive family classifications.
