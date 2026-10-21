# Source review of the 98-wall H20 checkpoint

**PASS.** I independently reviewed the [two-cap source](../scripts/fresh_h20_gap_next_20261013_exact.py), its arithmetic dependencies and the [frozen exact artifact](../data/enumerations/fresh_h20_gap_next_20261013_exact.json). Both commands passed:

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_gap_next_20261013_exact.py --verify
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_edge_family_completion_20261013_cusp_repairs.py --verify
```

This is an independent mathematical/source review and fresh replay of the existing exact arithmetic, not a second all-root implementation. The retained 96-wall certificate is an explicitly checked dependency; its earlier certificates remain inherited dependencies. The numerical scouts supply active-set indices only. Rank-20 row reduction reconstructs each one-dimensional kernel exactly, and every reconstructed point is checked against all retained walls before it is used.

The first construction has two orthogonal norm-two leaves `a,b` and a spacelike polar vector `p` with `p²=3` and `(p,a)=(p,b)=−2`. Thus

`r = a+b+t p`, `t=(4+sqrt(10))/3`,

has square `4−8t+3t²=2`; its pairing with `p` is `sqrt(10)>0`. Its common-center pairing is `−504−116 sqrt(10)<0`. The saved previous strict escape is separately rechecked against all 96 walls and is cut strictly. The second construction has three orthogonal norm-two leaves, each pairing `−1` with a norm-one polar vector `p`. Consequently `r=5p+a+b+c` has square `25+6−30=1`, pairs positively with `p` by `2`, and has center pairing `−596`. Its strict target is reconstructed by a rational perturbation and verified against all 97 preceding walls.

All 193 new pairs are checked exactly: 32 right, three `pi/3`, six `pi/4`, 151 dotted and one parallel. The recursive compositum embeddings used for old radical-field walls preserve both coefficients and squared norms. The sign routine does not assume independent radicals, and the angle gate checks nonpositive pairings before classifying a relation. The restricted angle list is sufficient to certify these accepted pairs; it is not presented as an exhaustive search over possible caps.

Both new walls pair strictly negatively with each of the 107 inherited cusp rays, so their existing local links are unchanged. These 214 signs do not claim preservation or enumeration of every later cusp or finite vertex. A common strict timelike center and all nonpositive off-diagonal pairings imply that all 98 inequalities are actual facets: projection of the center onto any wall remains timelike and strict at every other wall.

The final saved point has exact square `698631401/1887436800>0`, future-center pairing `−768229/6144<0`, and strictly negative pairings with all 98 walls. It therefore certifies infinite volume. The two successful gap repairs neither establish finite volume nor exclude further compatible repairs.
