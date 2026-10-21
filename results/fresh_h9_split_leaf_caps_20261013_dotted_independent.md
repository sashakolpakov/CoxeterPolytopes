# Independent audit of the 228-cap split-leaf construction

**PASS.** An independent exact reconstruction confirms all 228 norm-two
caps, every seed and cap-pair classification, the 8,460-edge compatibility
graph, and the twelve strict spacelike escapes. A separate NetworkX clique
enumeration finds exactly **156,990 maximal cliques**, of maximum size
**fourteen**. Every maximal clique has an exact escape. Consequently every
compatible subset of this literal cap bank has infinite volume.

This is a dictionary exclusion, not an exclusion of the sixteen-wall seed
or of arbitrary additional caps. It does not exactly certify or exclude
the separate 3,040-normal numerical scout. No finite-volume polytope is
claimed.

The [independent checker](../scripts/fresh_h9_split_leaf_caps_20261013_dotted_independent.py)
imports no repository constructor, arithmetic, graph routine, or numerical
search module. It implements its own Fraction arithmetic and exact signs
in `Q(phi,sqrt(phi))`, and uses `networkx.find_cliques` in place of the
producer's custom clique recursion. Its
[machine-readable result](fresh_h9_split_leaf_caps_20261013_dotted_independent.json)
records graph and clique digests, clique-size counts, individual escape
coverage, and negative controls. It audits the
[combined certificate](../data/enumerations/fresh_h9_split_leaf_caps_20261013_dotted_exact.json)
and the literal cap normals in the
[192-cap source](../data/enumerations/fresh_h9_split_leaf_caps_20261013_exact.json).
The completeness of the source's finite-leaf assignment enumeration is
separately reviewed; this audit directly verifies all its saved normals
and pairings rather than trusting that enumeration's assertions.

## Exact geometry and angle classification

The checker reconstructs every seed Gram entry from
`a=phi+1/2`, `b=-a/(12a+1)=-5/22+phi/11`, and the prescribed two groups
of four indices. The sum/difference invariant decomposition gives inertia
`(9,1,6 zero)`. The sixteen all-one coefficients represent a strict
timelike center, with square `-128a` and pairing `-8a` with each seed wall.
All seed pairs are right or strictly dotted.

All 192 inherited cap covectors are checked against their archived leaf
pairings. The 36 additional normals are independently reconstructed using
`d=phi(1+sqrt(phi))`, with exactly two plus signs in each four-block. Every
new leaf pairing is `-1` or `-d`, and every cap has square two and contains
the same strict center. Their 228 covectors are distinct, and none equals
a seed covector. Thus the Gram radical introduces no hidden duplicate
mirrors.

For norm-two normals the finite-angle pairing is `-2cos(pi/m)`. The exact
gate covers all positive integer Coxeter labels, for the following reason.

* In `K=Q(phi)`, degree at most two restricts possible labels to
  `m=2,3,4,5,6`. The values for four and six involve `sqrt(2)` and
  `sqrt(3)`, neither in `K`. The remaining pairings are exactly
  `0,-1,-phi`.
* Write an element of `L=K(sqrt(phi))` as `u+v sqrt(phi)`. The other real
  embedding of `K` sends `phi` to a negative number. If `v` is nonzero,
  its image is nonzero and the corresponding conjugate of this element
  is nonreal. This is impossible for a finite Coxeter value, whose
  conjugates are all real. Thus no omitted label can rescue a pairing
  with a nonzero radical coefficient.
* Pairing at most `-2` is separately accepted as parallel or dotted,
  using exact signs. Every other pairing is rejected.

The independent ordered edge list and sorted maximal-clique list have
the same SHA-256 digests as the producer. The literal archived
fourteen-cap selection is among these cliques, so it gives the advertised
thirty-wall configuration together with the seed. All its walls are true
facets: for an acute collection with a strict timelike interior point,
orthogonal projection to any wall remains timelike and strict on all the
other walls. This statement concerns the configuration, not finite volume.

## Why the escapes exclude volume

Each of the twelve archived vectors has rational coefficients in the
sixteen seed normals. The checker independently verifies positive Lorentz
square, negative pairing with the strict timelike center, and strict
negative pairings with every seed wall. It checks all 228 cap pairings and
recomputes the exact list of caps cutting each vector.

Every maximal clique admits at least one of these vectors with all its
selected cap pairings strictly negative. The segment from the common
timelike center toward that spacelike vector crosses the ideal sphere
while all inequalities remain strict. An open ideal patch therefore
survives, which forces infinite hyperbolic volume. Every compatible subset
extends to a maximal clique and inherits its escape, so maximal-clique
coverage proves the stated exclusion for every subset.

The checker rejects a changed new-cap norm, a zero vector presented as an
escape, and a forged cap-cut list. It uses no floating-point geometry.

```sh
python scripts/fresh_h9_split_leaf_caps_20261013_dotted_independent.py
python scripts/fresh_h9_split_leaf_caps_20261013_dotted_independent.py --verify
```

NetworkX is the only nonstdlib dependency. The replay recomputes the exact
geometry and the entire independent clique census; it does not merely
compare saved output files.
