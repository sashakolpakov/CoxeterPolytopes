# Independent review of the split-leaf cap banks

The exact seed and both declared cap banks pass an independent replay. This
establishes an infinite-volume obstruction for these fixed banks; it does not
exclude other caps or a completion of the seed outside them.

The checker uses the repository's `Golden` arithmetic instead of the producers'
Fraction-pair arithmetic. It imports neither producer. It solves the literal
Gram system on the basis consisting of all eight `U` normals and `V_0,V_4`,
rather than using the producers' closed formula for cap norms and pairings.

The seed Gram matrix has signature `(9,1,6)`. The six exact radical relations
are `U_i+V_i=U_j+V_j` within either block of four indices. Consequently a cap's
two leaf pairing magnitudes have constant sum on each block. The six distinct
sums of values in `{0,1,phi}` give exactly 51 ordered assignments per block.
The independent enumeration of all 2,601 assignments gives exactly 192 norm-two
caps. The ambient quotient has dimension ten, so compatible leaf pairings
determine the normal uniquely; there is no omitted perpendicular parameter.

The complete exact pair graph has 4,986 edges. Separate set-based clique
enumeration recovers all 8,478 maximal cliques, with at most six caps. Every
clique has an independently checked strict spacelike escape among the 36
archived witnesses.

The 36 additional dotted caps also pass the literal Gram equations, norm,
orientation, and leaf pairing checks. In `Q(phi,sqrt(phi))`, an element
`u+v*sqrt(phi)` with `v != 0` has a nonreal conjugate under `phi -> phi' < 0`.
It therefore cannot equal a finite Coxeter pairing `-2*cos(pi/m)`, which is
totally real. The values in `Q(phi)` are classified exactly: for norm-two
normals the finite possibilities are `0,-1,-phi`; pairings at most `-2` are
also admissible. This justifies the full pair gate, rather than an acute-angle
test alone.

The combined 228-cap graph has 8,460 edges and exactly 156,990 maximal cliques,
with at most 14 caps. Their graph and clique digests agree with the producer.
All 12 saved rational witnesses pass independent exact checks of positive
norm, future orientation, strict seed inequalities, and their complete cap
cut lists. Every maximal clique admits at least one of these strict escapes.
Every subset is contained in a maximal clique, so the same obstruction covers
every compatible subset of the declared bank.

The facet interpretation also has the needed geometric hypothesis. Every
bank cap pairs nonpositively with each seed normal and the common center is
strictly timelike. Projecting that center orthogonally to a seed wall makes
all other retained inequalities strict, so each seed wall remains a facet.
After any redundant optional inequalities are discarded, a Coxeter
completion's optional facets must form one of the pair-compatible subsets
already covered. A strict spacelike point and a strict timelike point give an
open ideal boundary patch and therefore infinite hyperbolic volume.

The malformed-normal regression is checked against the Gram equation itself.
No mathematical issue was found in either producer's stated finite-bank
scope. The source's broader numerical scout is not part of this exact result.

Artifacts:

- [Independent checker](../scripts/fresh_h9_split_leaf_caps_independent_20261013.py)
- [Independent audit result](fresh_h9_split_leaf_caps_independent_20261013.json)
- [192-cap source certificate](../data/enumerations/fresh_h9_split_leaf_caps_20261013_exact.json)
- [228-cap source certificate](../data/enumerations/fresh_h9_split_leaf_caps_20261013_dotted_exact.json)

Replay from the repository root:

```sh
PYTHONPATH=src python scripts/fresh_h9_split_leaf_caps_independent_20261013.py --verify
```
