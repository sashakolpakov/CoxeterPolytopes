# Exact finite-corner edge queue for the mixed D4 chamber

The [checker](../scripts/fresh_h9_mixed_d4_coupling_20261014_frontier.py) and [queue](fresh_h9_mixed_d4_coupling_20261014_frontier.json) classify all **108 directed edge releases** at the twelve simple finite vertices of the frozen twelve-wall chamber:

- **56 directions** reach another finite vertex, giving **28 distinct finite edges**.
- **52 directions** reach the light cone before any further existing wall. Each has a saved exact spacelike continuation inside its projective edge.

The finite-edge graph is connected and includes all twelve vertices. All ten old finite vertices remain; the polar addition contributes two further finite vertices, ray IDs **10 and 19** in the frozen complete cone census. The chamber still has infinite volume.

## Exact coverage and open-end witnesses

The input is the [independently checked complete fifty-ray cone](fresh_h9_mixed_d4_coupling_20261014_review.md), comprising twelve finite and thirty-eight spacelike extreme rays. This derivative reuses that independent quartic arithmetic; it is not a second scalar implementation or another full ray census.

Every finite vertex has nine incident facets with positive-definite rank-nine Gram matrix. Releasing one leaves a positive-definite rank-eight support. The checker verifies rank eight directly and finds the unique other cone ray containing that support. Their common active set is exactly those eight walls. A positive combination of the two rays has exactly these equalities and every other wall inequality strict. Thus it lies in the relative interior of the corresponding two-dimensional cone face, proving adjacency rather than merely suggesting it from shared labels.

For a finite other ray, the whole projective segment stays future timelike, giving an existing finite edge. For a spacelike other ray `q`, the queue records

`x(t)=p+t(q−p)`,

its exact norm quadratic, and a rational dyadic `0<t<1` for which `x(t)` is spacelike. All nonsupport inequalities remain strict at that point, and its time orientation agrees with the center. The unique light-cone crossing lies between the finite start and this witness. Consequently no further old facet appears at or before the crossing: these are genuinely unclosed directions, not unexamined finite or ideal cone endpoints. Every witness in this data uses a denominator at most eight.

The eighty distinct rank-eight supports comprise 28 finite edges and 52 open directions. Each finite edge appears exactly twice in the directed queue. The result covers edges incident to the twelve finite vertices; it does not enumerate vertex-free edges or claim a complete exterior face lattice.

## Concrete next supports

Most open directions release a D4 wall: wall 7 gives twelve, wall 8 gives eight, and walls 9 and 10 give twelve each. The other **eight** release source walls at the two new polar vertices. They retain the new wall 11 and provide contacts beyond the original finite-corner queue:

| Queue ID | Start ray | Released wall | Rank-eight support |
|---:|---:|---:|---|
| 18 | 10 | 0 | 1,2,3,5,7,9,10,11 |
| 19 | 10 | 1 | 0,2,3,5,7,9,10,11 |
| 20 | 10 | 2 | 0,1,3,5,7,9,10,11 |
| 21 | 10 | 3 | 0,1,2,5,7,9,10,11 |
| 45 | 19 | 0 | 1,2,3,6,7,9,10,11 |
| 46 | 19 | 1 | 0,2,3,6,7,9,10,11 |
| 47 | 19 | 2 | 0,1,3,6,7,9,10,11 |
| 48 | 19 | 3 | 0,1,2,6,7,9,10,11 |

Every row contains the exact finite start, projective direction, other cone-ray ID and spacelike witness. These are inputs for a future endpoint-cap construction. No cap search, complete arithmetic root enumeration, or first-new-mirror assertion was performed here.

Generation and saved replay pass:

```sh
python scripts/fresh_h9_mixed_d4_coupling_20261014_frontier.py --verify
```

All primary and review sources remain unchanged. This derivative binds the frozen primary JSON, independent arithmetic source/certificate and its own source; their source manifests are checked before use.
