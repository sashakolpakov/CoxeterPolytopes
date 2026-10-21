# Independent review of the H22 cusp-wall replacement

The [construction](../scripts/fresh_h22_core_replacement_20261013.py) passes source and theorem review and a fresh `--verify` replay. A [separate radical-arithmetic checker](../scripts/fresh_h22_core_replacement_20261013_review.py) also passes; its [output](fresh_h22_core_replacement_20261013_review.json) binds the frozen sources. No issue was found.

This is a changed **813-facet** chamber obtained from the 812-wall rational parent branch by deleting old wall 1 and adding two caps. It does not retain the entire later 860-wall branch. The caps cut the selected inherited cusp and create an actual finite vertex, but a strict spacelike escape remains. There is no finite-volume or literature-novelty claim.

## Construction and exact inequalities

Let δ be the verified cusp ray, let E be its 20 retained incident old normals, and let p be the orthogonal projection of the old strict center h to `E^perp`. The checked values are

`p²=-2484`, `(p,δ)=-46`, `t=350`.

Then `u=-p+tδ` has positive norm `U=29716`. In the physical form `Q ⊕ <2>`, the two added normals are

`r_plus=(u/sqrt(U),1/sqrt(2))`, `r_minus=(u/sqrt(U),-1/sqrt(2))`.

Both have norm two and their mutual pairing is zero. Both pair negatively with the old strict center, and both pair positively with `(δ,0)`, so the inherited cusp is actually cut. The old projection is spacelike; this is possible because an old wall was removed, and the earlier retained-core causal-projection obstruction no longer applies.

The producer's most delicate inequality is correct. If `A=(u,c)` for a parent null projection, the larger of the signed new/parent pairings is `A/sqrt(U)+sqrt(2)`. Strict dottedness is equivalent to

`A<0`, `B=A²-6U>0`, `B²>32U²`.

These signs make both squarings legitimate. The opposite height choice has a smaller pairing. Together with the old-wall tests and the new mutual right angle, the checks cover all **1,623 new unordered pairs**: 41 right and 1,582 strictly dotted.

The separate checker reconstructs the parent normals without importing the construction producer. It represents the new normals as sparse squarefree radicals and directly checks norms and every new pairing, using exact radical signs. Thus it checks the pair inequalities through a different arithmetic route from the producer's rational squared inequalities.

## Actual vertex, facets, and remaining escape

With `d=(p,δ)`, the proposed point is

`v=p+(t-p²/d)δ=p+296δ`.

Direct calculation gives `(u,v)=0`, `v²=-U`, and `(v,h)=-16100<0`. The separate checker confirms every retained-wall inequality at `(v,0)`, with equality precisely on the 20 old normals in E. The two new cap inequalities are also equalities. Their incident Gram is the checked positive-definite rank-20 Gram of E plus `2I2`, so this is an actual finite vertex of the H22 chamber.

Every inequality is a true facet. For a positive-norm normal r and a common strict timelike center h, its orthogonal projection `h-(h,r)r/r²` is timelike, lies on that wall, and remains strictly inside every other halfspace because the distinct-wall pairings are nonpositive. This also supplies a full-dimensional chamber, rather than merely a feasible collection of normals.

The separate checker verifies a strict future spacelike point against all **813** normals. Its square is **1586** and its pairing with the center is **−814**. Hence this construction still has infinite volume despite the new finite vertex and cusp cut.

## Conditional non-quasi-arithmetic conclusion

The new cap and a retained norm-two parent have an irrational cyclic Gram product in `Q(sqrt(2U))`. Since `2U=59432` is not a rational square and its coefficient is nonzero, the cyclic field contains that quadratic field. The retained rational short-root principal Gram has a negative all-ones quadratic value. It stays indefinite under every embedding, including one nontrivial on this quadratic cyclic subfield. This contradicts the conjugate positive-semidefiniteness required for a quasi-arithmetic finite-volume chamber. Thus a cofinite completion retaining this facet core would be non-quasi-arithmetic; the currently infinite-volume chamber is not being presented as such a completion.

The separate replay checks both the cyclic identity and the rational negative witness. The D20 name comes from the previously certified classification of cusp orbit 14. The present construction itself only needs, and directly verifies, the single affine incident component and its positive kernel.

## Replay scope

```sh
python scripts/fresh_h22_core_replacement_20261013.py --verify
python scripts/fresh_h22_core_replacement_20261013_review.py --verify
```

Both pass. The first includes the entire 812-wall base pair census. The second independently reconstructs and checks the new geometry, all 1,623 new pairings, the full-chamber center/vertex/escape conditions, and the field witness. It does not repeat the complete base pair census in a second implementation.
