# An all-label gate for an H4 face extension

This is a restriction on a proposed new compact-H9 seed, not a compactness
or nonexistence certificate for that seed. A compact H5 polytope can be placed
in the Lorentzian orthogonal complement of a spherical H4 root system. Any
added wall meeting the entire H4 face must be orthogonal to H4. Walls that
bound the four transverse directions must lie away from that face.

The [standalone checker](../scripts/fresh_h9_h4_transverse_gate_20261013.py)
uses only standard-library rational pairs representing `a+b*sqrt(5)`.
Its [certificate](fresh_h9_h4_transverse_gate_20261013.json) verifies positive
LDL pivots, a literal inverse and exact inequalities; no numerical predicates
or spherical-diagram classification are needed.

Normalize the H4 simple roots and a proposed new normal r to square two.
The H4 chain has off-diagonal entries `−(1+sqrt(5))/2, −1, −1`.
Every entry of its Gram inverse B is positive, and its diagonal is

```text
14+6sqrt(5), 21+9sqrt(5), 10+4sqrt(5), 3+sqrt(5).
```

Write b for the four pairings with r. Acuteness gives `b_i<=0`.
If any pairing is nonzero and is a Coxeter relation, its magnitude is at
least one (label three); parallel and dotted relations have larger magnitude.
Consequently

```text
r_H4² = bᵀ B b >= min(diag B) = 3+sqrt(5) > 2,
r_perp² = 2-r_H4² <= -1-sqrt(5) < 0.
```

The enlarged five-wall Gram cannot be positive semidefinite. Equivalently,
the timelike perpendicular component cannot be orthogonal to a point of the
hyperbolic H4 face. Thus a nonorthogonal wall cannot meet that face, even at
infinity. This argument allows arbitrary real ambient coordinates and every
Coxeter label. It leaves caps away from the face, mixed profiles and compact
completions open.

Producer and saved replay pass:

```sh
python3 scripts/fresh_h9_h4_transverse_gate_20261013.py
python3 scripts/fresh_h9_h4_transverse_gate_20261013.py --verify
```
