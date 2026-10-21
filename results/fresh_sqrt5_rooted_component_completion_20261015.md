# All-rank completion by one rooted finite white component

The 458 accepted first-contact profiles on the 130 saved nonsingular unpaid cores now have a complete rank classification when exactly one connected finite white component is added and only one of its roots touches the core. There are 234 isolated finite rank candidates: 209 pass the full corner and arithmetic matrix tests, while 25 have a nonfinite corner. In addition, 37 profiles each support an unbounded B-family and an unbounded D-family. No rank cutoff is used in this classification.

The [source](../scripts/fresh_sqrt5_rooted_component_completion_20261015.py) and [certificate](fresh_sqrt5_rooted_component_completion_20261015.json) preserve the [single-chain](fresh_sqrt5_single_chain_completion_20261015.md) and [two-chain](fresh_sqrt5_two_chain_completion_20261015.md) packages. This result does **not** seal these stars against additional components or attachments at other roots. It does not cover an external component with two or more roots touching the core, multiple nonorthogonal external components, dotted cores or arbitrary global completion.

## Complete rooted finite-type reduction

Every nonzero profile contacts some node belonging to an actual corner. Thus the entire external finite component, together with a single new neighbor at its attachment root and label at least three, must be spherical. Testing the weakest such contact, label three, gives a necessary rooted-diagram classification. The complete possibilities, up to automorphisms fixing the attachment root, are:

| External type | Permitted attachment roots |
|---|---|
| A_n | an endpoint; the next node from an endpoint; additionally the third node for n=5,6,7 |
| B_n | the endpoint opposite the 4-bond; additionally the 4-bond endpoint for B3 |
| D_n | the long-arm endpoint; additionally a fork endpoint for D5,D6,D7 |
| E6,E7 | an endpoint whose extension gives E7,E8 respectively |
| H3 | its endpoint opposite the 5-bond |
| I2(5) | either endpoint, which are equivalent |

The next-to-end A root starts at rank three. B2 has equivalent endpoints, and D4 has three equivalent leaves, so these are already in the unbounded endpoint families. E8, F4 and H4 cannot be extended spherically. I2(6) and I2(10) cannot have a third nonorthogonal node.

This is an all-rank consequence of the finite Coxeter classification. For A_n an interior attachment produces a tree with arm lengths `1,i−1,n−i`; finite branching gives the D case with one further arm of length one, or the E cases with the other arm lengths `2,2`, `2,3`, `2,4`. In B_n, attachment at an interior root introduces branching in a component containing a 4-bond, which is not finite. Extending its 4-bond endpoint when n≥4 contains the affine path `3–3–4–3`. In D_n, an interior attachment produces two trivalent vertices or a vertex of degree four. Extending a fork leaf when n≥8 contains the affine tree with arm lengths `1,2,5`. These arguments cover arbitrarily long classical diagrams, rather than extrapolating a finite census.

As a consistency check, the certificate tests all 1,917 rooted finite diagrams through rank 35, plus the two isolated high-dihedral failures, and matches the list exactly. The 221 passing rooted tests include equivalent roots; they are not 221 different types. Eleven fixed rooted exceptions remain after the four classical families are separated.

## Exact rank equations

Write `H` for the nonsingular core Gram matrix, `v` for its contact column, `d` for the attachment root's squared norm, and `c=v^T H^{-1}v`. If `J` is the inverse of the external finite component's norm-two Gram matrix, the determinant lemma forces

\[
2J_{rr}\,c/d=1.
\]

For A_n, `J_ii=i(n+1−i)/(n+1)`. At an endpoint the only possible rank is `n=1/(2c/d−1)`. At the next-to-end root it is `n=(4c/d+1)/(4c/d−1)`. Each value is evaluated exactly in Q(phi); it must be a positive rational integer in the specified range. A zero denominator does not give a finite A rank. The eleven fixed rooted exceptions are checked by exact inverses, followed by full four-corner and physical/conjugate inertia checks.

For the ordinary endpoint of B_n and the long-arm endpoint of D_n, `J_rr=1` for every n. In norm-two leaf elimination, the B weighted leaf gives `2−2/2=1`, the D fork gives `2−1/2−1/2=1`, and every preceding simple bond preserves this value because `2−1/1=1`. Hence these families require exactly `c/d=1/2`.

All 37 saved profiles with this ratio are singleton label-three contacts. At each actual corner, the contacted old component is a simply-laced path with the contact at an endpoint, or there is no contact in that corner. The certificate records these four literal path witnesses for every profile. Grafting B_n or D_n therefore yields respectively a B or D component with a longer arm, so all four corners remain finite for every allowed n.

Eliminating the external component always makes the same old-core Schur correction, `2vv^T/d`, independently of n and of the B/D choice. The exact B2 and D4 representatives have identical projected four-normal Grams. A separate Fraction reconstruction gives four timelike corner rays, rejects both opposite intersections, and supplies a common strict future center. Thus every member of these 74 all-rank families has the same actual compact quadrilateral face geometry, although its ambient white support grows. No compactness or finite volume of the whole ambient chamber is claimed. The signatures are also all-rank consequences of the determinant lemma and the positive semidefinite critical rank-one update, not numerical extrapolations.

The 209 isolated finite candidates have full physical inertia `(N−2,1,1)` and conjugate inertia `(N−1,0,1)`, finite corners and no witness from the inherited affine search. Their actual quadrilateral orientation is not independently asserted here. They may be included conservatively in a sufficient weight model; infeasibility of that relaxed model alone would not prove a geometric obstruction.

## Dimension-29 input

The certificate exports 1,715 literal instances relevant to dimension 29: all 209 isolated finite stars and 1,506 members of the B/D families. Each row has at most 31 normals, the maximum for a full two-face star in dimension 29. This finite restriction belongs only to the export; the preceding rooted-family theorem has no dimension cutoff.

Each `rank29_instances` record contains its full diagram, four corner forests and marked-pair orbit keys, the remaining orthogonal spectator rank, and every possible spectator A1 count. If the remaining rank is `s`, the possible counts are the integers `j` with `0≤j≤s` and `s−j≠1`: the residual rank is either zero or can be supplied by an A-component of rank at least two. This records the A1-sensitive contexts without assuming that the component is sealed against other nonorthogonal additions.

These are new inputs for a weight stress test, not a new dimension bound. General attachments involving multiple roots of a component remain open. Generation, saved replay and syntax checks pass; this package has not yet received a separate independent reconstruction.

```sh
python3 scripts/fresh_sqrt5_rooted_component_completion_20261015.py --verify
```
