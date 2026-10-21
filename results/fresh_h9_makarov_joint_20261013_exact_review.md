# Review of the exact 27-wall Makarov-face construction

**PASS.** The frozen construction gives 27 actual Coxeter facets, retains the
10 original finite vertices, and closes all 40 original transverse edges at
40 distinct new finite vertices. Its saved strict future spacelike escape
proves that the resulting H9 chamber still has infinite volume.

This is a source and geometric proof review, followed by a fresh saved-result
replay of the constructor's exact arithmetic. It is not a second implementation
of all normal or endpoint arithmetic. The root agent's separate
[pair and escape certificate](fresh_h9_makarov_joint_20261013_exact_independent.json)
provides a distinct arithmetic check of those claims.

The checkpoint uses the 11 original seed walls and 16 new caps. It replaces
the earlier changed-center twelfth wall; it is not an extension retaining
that 12-wall chamber.

The [exact constructor](../scripts/fresh_h9_makarov_joint_20261013_exact.py)
reconstructs each cap from nine prescribed linear pairings and a norm-two
quadratic. Floating comparisons select one of the two proposed branches;
the chosen normal's norm and prescribed pairings are then checked exactly.
The exact acceptance checks cover every one of the 351 pairs: 161 right,
31 with label 3, 12 with label 5, and 147 dotted. There are no parallel
pairs in this particular checkpoint. No completeness of the bounded cap
search is inferred from these acceptance checks.

The common center is timelike and strictly satisfies every inequality.
For each outward normal r, projecting this center c to its mirror gives
`p=c-(c,r)r/(r,r)`. Its norm is strictly negative. For every other normal s,
`(p,s)=(c,s)-(c,r)(r,s)/(r,r)<0`, since both `(c,r)` and `(r,s)` are
nonpositive and the first is strict. Thus each of the 27 inequalities
supports a genuine facet. This also rules out positive proportional
duplicates.

The endpoint proof is particularly useful. At each old finite vertex v,
the inverse H4 Gram matrix supplies the vector d dual to a released H4 wall:
it pairs to one with that wall, to zero with the other incident walls,
and is orthogonal to v. The line `x(t)=v-t*d`, `t>0`, is therefore the
outgoing edge determined by the other eight incident walls. All new caps
are strict at v. The constructor proposes an owner and a positive parameter
`t=(r,v)/(r,d)`, then checks the resulting point against all 27 inequalities
using exact compositum arithmetic. Its zero set is precisely the eight
old edge walls and its new owner; every other inequality is strict.

These exact checks prove first intersection even though the proposed
minimum parameter is chosen numerically. Each wall pairing is affine in t.
For a nonincident wall, strict negativity at both ends gives strict
negativity on the whole segment. The released wall has pairing `-t`;
the owner is strict before its first zero. Hence no other retained wall
can have been crossed earlier. This is first intersection with the
specified 27-wall chamber, not a first-mirror theorem for an unspecified
reflection lattice.

Each endpoint has negative norm and the correct future orientation.
The eight old incident normals have positive-definite Gram matrix, and
the exact positive Schur complement after adjoining the owner proves
that the full nine-wall incident Gram matrix is positive definite.
Consequently these are actual simple finite vertices. The 40 incident
sets are explicitly checked to be distinct. All 10 original vertices
remain feasible, with every new cap strict there.

The saved escape is strict at all 27 walls, has negative pairing with the
future center, and has exact positive norm

```text
276679291647808184 - 123399368968947096*sqrt(5).
```

A segment from the strict timelike center toward this strict spacelike
point crosses the ideal boundary through a relatively open feasible
region. The chamber therefore has infinite hyperbolic volume. Closing
the original 40 edges does not imply a complete vertex or edge census,
compactness, or finite volume; new frontier components remain possible.

Fresh replay passed:

```sh
python3 scripts/fresh_h9_makarov_joint_20261013_exact.py --verify
```

Reviewed frozen SHA-256 bindings:

| File | SHA-256 |
| --- | --- |
| `scripts/fresh_h9_makarov_joint_20261013_exact.py` | `b1ee2821eaeccf39d3fcd567b679699d78acabf47e674989aaca0d535d76950e` |
| `results/fresh_h9_makarov_joint_20261013_exact.json` | `cf1473c23c6cfc6c6a0a9ae2ed935de0d00454909e80b2279887cbbcacce02df` |

The exact artifact also binds its seed artifact, seed arithmetic source,
selection recipes and both numerical proposal banks. No unresolved
mathematical issue was found within this stated scope.
