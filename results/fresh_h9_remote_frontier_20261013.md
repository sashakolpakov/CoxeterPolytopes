# An inherited H9 face closes, followed by an H3 endpoint

Two exact additions take the chamber from 104 to **106 walls**. The first closes an inherited H2 quadrilateral; the second adds an H3+6A1 finite endpoint while preserving that quadrilateral. Each cap cuts the preceding strict spacelike escape and retains all preceding walls. A different exact strict spacelike escape survives at the final stage, so this is not a compact or finite-volume H9 example.

The construction changes the frontier component being explored. It uses an inherited two-vertex component, outside the 16 new edges covered by the [preceding all-angle obstruction](fresh_h9_component_bridge_20261013_continuous.md). It retains all 16 seed walls and the previously constructed compact H4 prism.

The [frozen cap and escape](../data/enumerations/fresh_h9_remote_frontier_20261013_caps.json) and [exact certificate](../data/enumerations/fresh_h9_remote_frontier_20261013_verify.json) are replayed by the [producer-free checker](../scripts/fresh_h9_remote_frontier_20261013_verify.py). Its dependencies are the preceding exact normal formulas and source-bound certificates; it imports no numerical search or optimizer.

## Exact cap and compact face

Let the eight orthogonal norm-two walls be

`seed1, seed2, seed5, seed6, seed7, wall22, bridge0, bridge1`.

Write their span as A, let `a = nextH4`, `e = wall22`, and set `d = (a,e) < -2`. The normal `b = a - (d/2)e` is orthogonal to A and has negative norm `q_b = 2-d²/2`. The particular solution is

`u = -(seed1 + bridge1)/2 - b/q_b`, with `u² = 1 + 1/q_b`.

Project the old strict center orthogonally off A and then off b to obtain z, whose squared norm is strictly positive. The new outward normal is

`r = u + sqrt((2-u²)/z²) z`.

The orthogonal decomposition proves `r² = 2` and all nine prescribed finite pairings exactly. The chosen positive branch has the correct strict-center orientation. The 104 new pair relations are six right angles, three angles π/3, and 95 dotted pairs. All 512 subsets of its nine finite neighbors are checked: 128 contain an already dotted pair and 384 are positive definite. None is affine. Together with the inherited audit, the 105-wall chamber has no positive affine subdiagram.

The seven common walls

`seed1, seed2, seed5, seed6, seed7, bridge0, bridge1`

have Gram matrix `2I7`. Their H2 face has sides, in cyclic order,

`seed3, wall22, remoteQuadrilateral, nextH4`.

The projected side norms are `2,2,1,2`, and its consecutive intrinsic angles are `π/2, π/2, π/4, π/2`. Both pairs of opposite sides are disjoint. The checker constructs all four vertices as timelike projections, proves the nine incident walls have positive-definite Gram matrix, and checks every other one of the 105 inequalities strictly. Thus this is an actual compact Lambert quadrilateral, not merely a compatible subdiagram. Its four vertices have four finite edges and 28 remaining open rays in H9; all possible ninth-wall endpoint tests at those edges are included in the certificate.

Every wall remains an actual facet: the chamber has a strict timelike center and all distinct outward pairings are nonpositive. Orthogonal projection of that center onto any wall is timelike and satisfies every other wall strictly.

## The second cap and final 106-wall chamber

The next exact cap uses the orthogonal edge

`seed2, seed3, seed5, seed6, seed7, wall22, bridge0, bridge1`.

It has angle π/3 with `seed5` and `nextH4`, angle π/5 with `bridge1`, six right relations, and 96 dotted relations. Use the same b and q_b as above, replace the base by `-(seed5 + φ bridge1)/2`, and project the center off this new edge and b. The base norm is `(2+φ)/2`; the same final square-root formula fixes norm two and the anchor pairing exactly.

The [second cap and escape](../data/enumerations/fresh_h9_remote_frontier_20261013_second_caps.json) and [second exact certificate](../data/enumerations/fresh_h9_remote_frontier_20261013_second_verify.json) are checked by a [separate replay entry point](../scripts/fresh_h9_remote_frontier_20261013_second_verify.py). The new finite vertex has incident diagram H3+6A1, is timelike, and satisfies all other walls strictly. The nine-edge endpoint census at this vertex has one finite edge and eight open rays. That finite edge joins a quadrilateral corner: its two endpoint walls are `seed1` and the new cap. All four quadrilateral vertices are explicitly rechecked against every one of the 106 walls. Its complete 512-subset census has 128 dotted exclusions, 352 positive-definite subsets, 32 indefinite subsets, and no affine subset. Thus the resulting whole chamber is still affine-free.

The final strict spacelike witness has norm

`(4738398474900536745175249 - 2499598252246300059056117 φ)/(16·10^24) > 0`,

and center pairing `-558175253991/250000000000 < 0`. Every one of its 106 inequalities is strict. Its supporting numerical vertex involved

`seed0, seed3, seed5, seed7, wall16, wall29, wall34, wall69, mirrorFaceGap`.

The final numerical census had 81,886 projective vertices, including 81,846 spacelike vertices. As at the preceding stage, only the exact strict escape is used to prove infinite volume.

## The intermediate 105-wall volume gate

The saved new strict future spacelike point has norm

`40108497102564153266229/10^23 - (11014835526582018390521/(5·10^22)) φ > 0`,

and center pairing `-653064609353/250000000000 < 0`. All 105 halfspace inequalities are checked strictly with rigorous rational interval bounds. It proves infinite volume independently of the numerical vertex census. Its supporting numerical vertex involved

`seed2, seed3, seed6, seed7, seed8, wall37, wall70, bridge1, facePrism`.

The numerical census proposed the witness from 80,726 projective vertices; 80,687 were spacelike. Those counts are diagnostics, not the volume proof.

## Discovery scope and replay

The numerical search examined 226 inherited open edges, with finite edge labels 2, 3, 4, 5 and a finite anchor alphabet. It found 2,513 necessary continuous intervals and 2,869 distinct rounded proposals. It does not exhaust other finite labels, the 1,075 recorded isolated-dihedral tails, or arbitrary cap families. Only the selected cap and its geometric conclusions above are certified here.

Run:

```sh
python scripts/fresh_h9_remote_frontier_20261013_verify.py --verify
python scripts/fresh_h9_remote_frontier_20261013_second_verify.py --verify
```

The [numerical scout](../scripts/fresh_h9_remote_frontier_20261013_scout.py), [proposal ranking](../scripts/fresh_h9_remote_frontier_20261013_rank.py), and [first](../scripts/fresh_h9_remote_frontier_20261013.py) and [second](../scripts/fresh_h9_remote_frontier_20261013_second.py) volume-witness producers are separate from replay. The second selected proposal is ID 2670. The neighboring proposal 2669 was rejected because it creates an affine subdiagram. That rejection is only a numerical discovery diagnostic; the accepted cap's complete affine census is exact.

The [separate source/formula review](fresh_h9_remote_frontier_20261013_review.md) passes, including a fresh 160-digit replay with all 8,255 changed predecessor/final interval enclosures nested inside their frozen 80-digit bounds. It independently reviews the orthogonal-decomposition identities and geometry but reuses the interval engine; it is not a second implementation of all root arithmetic.
