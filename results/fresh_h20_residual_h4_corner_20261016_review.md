# Independent 145-wall H20 corner reconstruction

The [145-wall construction](fresh_h20_residual_h4_corner_20261016.md) has
**one new finite corner**, at least 171 finite vertices in total, four retained
proper cusps, and 145 actual facets. A whole rectangle of strict spacelike
escapes remains, so the chamber has infinite volume.

The [independent reviewer](../scripts/fresh_h20_residual_h4_corner_20261016_review.py)
and [certificate](fresh_h20_residual_h4_corner_20261016_review.json)
reconstruct the physical old chamber using the earlier independent field
and interval implementations. No producer arithmetic is imported. The old
chain parameter and its exact isolation remain named reviewed dependencies.
The five affine-A17 connector normals are reconstructed from their projection
and Schur-Laplacian formulas, rather than read as producer intervals.

The four residual normals, numbered 141 through 144, have rational coordinates
in the chart (x,e,f,omega17,omega31), with Gram
diag(-1,1,1,1/2,1/2). Their exact squares are

    5/16, 309/200, 19717/22500, 6567095/1058469144.

They are pairwise orthogonal and right to the remaining sixteen support
roots. The independent review checks all 570 added pairs beyond wall 140:
71 right and 499 dotted. Exact chart identities prove the right contacts;
strict rational interval inequalities prove the dotted contacts. The common
strict center and its facet projections establish all 145 actual facets.

A fresh rational solve, reversing the four orthogonality equations and
fixing the time coefficient to one, reconstructs their unique common corner
with the sixteen support roots. Its exact norm is

    -10669433158276/103116477706827025.

The sixteen-root support Gram is positive definite, and the four positive
orthogonal squares complete a positive-definite rank-20 corner Gram. Every
nonincident wall inequality is strict. Exhausting the finite-neighbor sets
of all four new walls leaves exactly this support; 19 alternative supports
contain a certified dotted pair. Thus there is exactly one new finite vertex,
not merely one displayed vertex. Old finite vertices persist by the reviewed
acute-extension lemma; the total is a lower bound because no complete old
ambient cone census is asserted. All four old proper cusp directions remain
strict on the four new walls.

The renewed escape has frame coefficients (1,-1/100,-49/50), dual
coefficients -a and -b on omega17 and omega31, and -1/100000 on the other
sixteen duals. Every point with

    29/100 <= a <= 31/100,   3/100 <= b <= 3/50

satisfies every wall inequality strictly. Pairings are affine in a,b and
are checked at all four corners. The norm is a constant plus (a^2+b^2)/2,
whose minimum on this positive rectangle is exactly
15000061/5000000000 > 0. The specific E4 at a=3/10,b=1/20 has norm
33750061/5000000000. This certifies an entire remaining gap, not only one
numerical direction.

The independent review includes the preceding mixed wall 141, which had
previously had only primary exact checking. This is a coherent construction
and obstruction package, not a finite-volume example or a new dimension.

```sh
python3 scripts/fresh_h20_residual_h4_corner_20261016_review.py --verify
```
