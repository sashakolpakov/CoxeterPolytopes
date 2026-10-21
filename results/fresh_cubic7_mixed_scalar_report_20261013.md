# Mixed dihedral blocks and scalar completion in H9

At the unchanged height bound 1000, allowing one arbitrary dihedral block, one equal-coordinate block and the scalar coordinate produced 663 mirrors compatible with the original 25-wall cone. Seventy-six cut the previous common escape. The resulting Coxeter cones have 742 or 706 walls, but a new exact spacelike vector satisfies the entire combined 855-root dictionary, the original 25 walls, and the separate Pell cap. Consequently **every compatible selection from this bounded union retains an infinite-volume escape**. No finite-volume example or full-host nonreflectivity theorem is claimed.

Write τ=2cos(2π/7), β=2−τ, C=τ²+τ−1 and D=2−τ². The host is

```
<-2τ> + B^4 + <1>,     B=[[2,-C],[-C,2]].
```

The free block has covector `(u,w)` and coordinates

```
y=((2u+Cw)/β, (Cu+2w)/β),    E=u*y1+w*y2.
```

The equal block is `(a/D,a/D)`, of norm `2a²`; D is an integral unit. The exact join solves

```
E + 2a² + z² = q + 2τx².
```

The complete normalized primitive norm pool is `{1,2,β,2β}`. The [host audit](fresh_cubic7_dihedral_independent_audit_20261013.md) proves the ring, unit normalization and norm divisibility, and excludes primitive norms 4 and 4β. Unlike the equal-block-only slice, the free block can have β-divisible covector without β-divisible coordinates, so β and 2β must be retained. The absence of norm-β solutions in this run is only a bounded result.

The enumeration covers `x>0`, `0≤u≤w`, `a,z≥0`, including equal and zero free blocks, in all twelve ordered block placements. These inequalities describe the original zero-height chamber. With H=1000, the distinguished-embedding bounds are

```
x² ≤ Hq/(2τ),     u²,w² ≤ 2q(H+1),
a² ≤ q(H+1)/2,    z² ≤ q(H+1).
```

At each other embedding they are `x'²≤−q'/(2τ')`, `u'²,w'²≤2q'`, `a'²≤q'/2`, and `z'²≤q'`. The free-block bounds follow from Cauchy's inequality for the positive block. The enumerator uses exact three-embedding alphabets, keeps all pairs sharing an energy, reconstructs integral coordinates, and checks integrality of `2Gr/q`. The [independent completeness review](fresh_cubic7_two_block_completeness_audit_20261013.md) covers these bounds, the alphabet implementation and the mixed support join. This is complete within the stated support and height, not a complete ambient Vinberg prefix.

| Norm | Original-25-compatible mirrors | Cut previous escape | Compatible with each previous cone |
|---|---:|---:|---:|
| 1 | 32 | 3 | 28 |
| 2 | 559 | 71 | 512 |
| β | 0 | 0 | 0 |
| 2β | 72 | 2 | 20 |
| Total | 663 | 76 | 560 |

All candidates are saved, including those that do not cut the previous witness. Exact ray deduplication merges this dictionary with the earlier complete two-block and equal-block/scalar dictionaries, leaving 855 roots. Retaining the previous 182-wall cone and adding its 560 compatible new walls gives 742 walls. The two previous 146-wall cones each become 706-wall cones. A fresh height-ordered selection from the original 25 walls and the full union also gives 706 walls. This last selection rejects the separate high Pell cap by the exact cross-angle gate. Thus the search also permits replacing earlier choices.

Every selected pair passes the finite Coxeter-angle test: finite `4cos²` values are restricted to `{0,1,2,3,2+τ,2+C}`. In particular, the arithmetic acute angle `2π/7` is rejected. The [Fraction-only replay](fresh_cubic7_target_independent_audit_20261013.json) independently checks all 663 candidates, the 855-root union provenance, all four selected cones, reflection integrality, angle gates, strict interiors and exact escapes. It also rejects a forged root norm.

The common escape has exact positive norm

```
(53964 - 19648τ - 18749τ²)/5041 > 0.
```

The [common-escape certificate](fresh_cubic7_mixed_scalar_common_escape_20261013.json) gives its coordinates and all 880 pairings against the original walls and union dictionary. Its pairing with the Pell cap is

```
(-762 - 2661τ - 1284τ²)/71 < 0.
```

The [separate independent replay](fresh_cubic7_common_escape_independent_audit_20261013.json) checks all 881 inequalities and a strict spacelike perturbation using a common strict timelike interior. The full dictionary is not mutually Coxeter; the conclusion concerns every compatible subset. The certificate does not exclude higher roots in these slices, other supports, or a finite-volume chamber of the full host.

The [search artifact](fresh_cubic7_mixed_scalar_h1000_20261013.json) records all alphabets, norm joins, dictionaries and four geometries. The producer took 398.4 seconds, largely in repeated exact pair checks. The common-escape replay is inexpensive:

```sh
PYTHONPATH=src:. python3 scripts/fresh_cubic7_mixed_scalar_search_20261013.py --height 1000 --output results/fresh_cubic7_mixed_scalar_h1000_20261013.json
python3 scripts/fresh_cubic7_mixed_scalar_common_verify_20261013.py
```

The `pell_cap_status` value `NOT_PROCESSED_SEPARATELY` means the routine inherited its input wall set without a separate cap test. Actual wall membership was checked: the cap occurs in mixed cases `retained182` and `retained146`, and is absent from `replacement146` and `replacement_from25`.
