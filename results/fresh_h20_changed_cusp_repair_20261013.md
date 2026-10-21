# Repairing the changed H20 cusps and moving to the remaining gap

This round **closes all seventeen specified deficient inherited cusp links** of the changed 179-wall chamber. It then leaves that cusp family and constructs an exact quadratic wall cutting the remote escape and creating a finite vertex. The resulting **214-facet chamber still has infinite volume**. A separately archived rational 215-wall branch cuts the next escape but introduces four new unbounded cusp links. No finite-volume example or new dimension has been established.

The retained input is the [changed-core 179-wall construction](fresh_h20_changed_assembly_20261013.md), not the earlier fixed core whose last facet had only eighteen finite neighbors. Its form has signature (20,1), with an added orthogonal coordinate `z²=3` and strict center `h'` of square `−7599/2`.

## Coherent repair of the seventeen links

For each of the seventeen archived primitive null rays δ, add

```
C_δ,+ = 5δ + z,       C_δ,− = 5δ − z.
```

The two normals have square 3, pairing −3, and sum `10δ`. Both are orthogonal to every old wall incident to δ. All other source cusp rays pair strictly negatively with them. The full **6,647 new pairs** pass exactly: 690 right, 17 parallel and 5,940 dotted. The same center remains strict, so all 213 inequalities are actual facets.

The previous positive affine components had quotient rank 18. The added parallel pair supplies the independent vertical direction. The verifier checks that its strictly positive relation includes **every incident normal**, sums to a positive multiple of δ, and that the actual span has rank 20. Their images positively span the 19-dimensional horospherical quotient; consequently every one of these seventeen links is bounded. This is a positive-link certificate, not a conclusion drawn from rank alone.

The [213-wall archive](../data/enumerations/fresh_h20_changed_cusp_repair_20261013.json) contains the 34 normals, all new pairings, all 204 inherited-ray signs, complete positive weights and rank checks. The [checker](../scripts/fresh_h20_changed_cusp_repair_20261013.py) recursively binds its predecessor certificates. It does not claim a complete census of other ideal vertices.

## Why changing cusp heights cannot close the remote gap

The old strict spacelike escape survives all 34 new walls. More strongly, the certificate covers **every real height** in each of these seventeen cusp planes.

Any spacelike normal in `span{δ,z}` can be positively rescaled to `r=tδ+σz`, with `σ=±1` and arbitrary real `t`. Write `a=−(δ,r_j)>0` for a retained core wall. Acuteness gives

```
t ≥ A = 3σ (r_j)_z/a.
```

Let `x` be the saved strict escape. Since `(x,δ)<0`, the necessary condition `(r,x)≥0` gives the rational upper bound

```
t ≤ T = 3σ x_z / (−(x,δ)).
```

In **32 of the 34 signed cases**, the explicitly checked acuteness bound satisfies `A>T`. In the remaining two cases—cusps 23 and 26, negative height sign—the relevant short-wall cosine is `2t/√6`. The exact inequalities force it strictly between 0 and 1/2 throughout `[A,T]`. No Coxeter angle lies there, including any arbitrarily high label; dotted and parallel relations are also impossible in that interval.

Thus every acute-compatible Coxeter wall in these seventeen planes is **strictly negative at the same escape**. Any finite collection of such additions leaves infinite volume. This theorem does not restrict general normals, other cusp directions or replacement of retained walls. It is the reason the search moved to a different residual space instead of extending a height grid.

## An exact cap at the remote gap

Seventeen active short walls at the remote gap form a positive-definite block. Their orthogonal complement has signature (3,1). Three prescribed finite-angle equations and the norm equation provided numerical proposals; one particularly simple proposal is exactified in [the 214-wall certificate](../data/enumerations/fresh_h20_changed_cusp_repair_20261013_gap_exact.json).

In the existing 21-coordinate basis set

```
u = [16,−8,−2,1,1,4,9,11,12,−7,−5,−4,−1,−1,−2,5,6,8,−3,−2,0]/16,
u² = −1/64,
s = (144+32√34)/11,
r = s u + (1/2−3s/32) z.
```

Then `11s²−288s−1280=0` and `r²=2`. The wall is orthogonal to nineteen short normals and pairs −3 with old long wall 55, giving angle π/3. Its complete 213-pair census is **19 right, one π/3 and 193 dotted**. It is strictly negative at all 204 inherited null rays, so it preserves the seventeen newly completed links.

The nineteen orthogonal short walls have a positive-definite Gram matrix. Together with `r` they form a positive-definite rank-20 system. Exact projection of the strict center gives a timelike point with precisely these twenty zero inequalities and every other inequality strict: an actual finite vertex.

This cap cuts the previous escape by

```
(−19063219806 + 16893386107√34)/110000000000 > 0.
```

Nevertheless a new strict future spacelike point of square **6026371665/33554432** survives all 214 walls. The [exact checker](../scripts/fresh_h20_changed_cusp_repair_20261013_gap_exact.py) verifies the full new pair gate, target cut, actual vertex, inherited cusp signs and new escape.

## A rational next branch and its explicit obligations

At the next gap, a rational norm-one cutter is

```
[16,−6,−4,3,2,4,9,11,12,−3,−6,−3,−2,−2,−4,6,6,8,0,−2,1].
```

It retains all 214 facets, cuts their escape by **98031/8192**, and preserves the preceding finite vertex strictly. Its new pair counts are 23 right, seven π/4, four parallel and 180 dotted. The [215-wall branch](../data/enumerations/fresh_h20_changed_cusp_repair_20261013_second_exact.json) has a new strict escape of square **23588681583123127187/50000000000000000000**.

The four parallel neighbors are old walls 30, 34, 46 and 50, each of square 9. For each neighbor `r_j`, the new null ray is `δ_j=r+r_j/3`. The [branch checker](../scripts/fresh_h20_changed_cusp_repair_20261013_second_exact.py) verifies every one of the 215 inequalities at all four rays and their complete incident diagrams.

Each link has twenty incident normals, actual span rank 20 and total quotient rank 19. **None is bounded.** Its positive affine part is just the parallel pair, of rank one. The other eighteen dimensions belong to spherical components: two rank-nine components at three rays, and three rank-six components at the fourth. For every link the certificate supplies an explicit positive-square vector perpendicular to δ and nonpositive on every incident normal. These recession vectors prove unboundedness without a rank heuristic.

This 215-wall branch is therefore recorded with concrete new cusp obligations, not treated as a successful closure. The source data provides the finite components for a later coupled completion; no such completion or exhaustive obstruction is claimed here.

## Scope and reproduction

The bounded [first residual-space scout](../data/enumerations/fresh_h20_changed_cusp_repair_20261013_gap_scout.json) and [next-gap scout](../data/enumerations/fresh_h20_changed_cusp_repair_20261013_gap_scout_next.json) supplied proposals only. A [single released-short-wall scout](../data/enumerations/fresh_h20_changed_cusp_repair_20261013_gap_scout_next_drop4.json) also remains numerical; it is not an exact bank, a completeness theorem or a volume result. Only the explicitly archived 213-, 214- and 215-wall constructions above are certified.

Fresh deterministic checks passed:

```sh
python scripts/fresh_h20_changed_cusp_repair_20261013.py --verify
python scripts/fresh_h20_changed_cusp_repair_20261013_gap_exact.py --verify
python scripts/fresh_h20_changed_cusp_repair_20261013_second_exact.py --verify
```

All three use the frozen previous pair certificates for retained-old relations and check every newly introduced pair directly. They distinguish closing specified cusp links, creating an actual finite vertex and cutting a particular escape from the separate requirement of global finite volume.

The separate [213-wall source review and fresh replay](fresh_h20_changed_cusp_repair_20261013_review.md) passed after the rational-rank and full-support positivity assertions were added. It reviews the bounded-link and all-real-family proofs using the existing field implementation; it does not claim an independent reconstruction of the later 214/215 stages.

The [root source review and fresh replay of the 214-wall step](fresh_h20_changed_cusp_repair_20261013_gap_review.md) also pass, including the actual finite vertex and separate strict escape.
