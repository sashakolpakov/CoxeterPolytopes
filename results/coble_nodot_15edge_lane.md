# The 15-edge no-parallel lane isolates Type II

**Status (2026-08-09): exact finite lane theorem, complete Vinberg audit,
and exact split-K3/Coble return.  It is not the whole nodot denominator and
not the class-wide AG pushback theorem.**

The executable certificate is
[`certify_coble_nodot_15edge_lane.py`](../scripts/certify_coble_nodot_15edge_lane.py),
with [archive](../data/enumerations/coble_nodot_15edge_lane.json) and
[tests](../tests/test_coble_nodot_15edge_lane.py).
An independent implementation reproduces the complete rank histogram,
binary glue census, GMP cone counts, and both marked Type-II isometries; see
the [hostile audit](coble_nodot_15edge_lane_hostile_audit.md).
The frozen archive SHA-256 is

```text
c9e2da6bcc08736fbc9848d64ca3fd273c23d535394ef727b671f369e7b20dae
```

## 1. A P-only finite denominator

Take twelve proposed norm-minus-two Coble--Mukai roots.  Their ordinary
nonorthogonality graph is required to be connected on twelve vertices, have
exactly fifteen edges and maximum degree three, and have pairing one on every
edge.  Thus there are no pairing-two/parallel pairs.  No surface, boundary,
characteristic, cusp count, determinant, known diagram, or volume verdict is
included in the input.

Nauty's complete unlabelled enumeration followed by exact Gram arithmetic is

```text
5,741 connected supports
  -> 5,741 unweighted trials
  -> exactly two Lorentzian corank-two root presentations
  -> three even-unimodular E10 extension branches.
```

The two root presentations are nonisomorphic.  Their saturated radical
quotients have determinants `-4` and `-1`.  Full discriminant-form
enumeration gives two index-two unimodular extensions of the first and the
unique index-one extension of the second.  Every extension is rechecked to
be even, unimodular, of signature `(1,9)`, and to contain all twelve roots
integrally.

## 2. Exact volume classification

All three branches receive complete GMP double description and literal
Vinberg endpoint and full affine-neighbourhood tests.

| quotient | E10 branches | finite | ideal | hyperideal | verdict |
|---|---:|---:|---:|---:|---|
| determinant `-4` | 2 | 38 | 17 | 0 | finite volume, noncompact |
| determinant `-1` | 1 | 18 | 8 | 24 | infinite volume |

For each determinant-four branch, all twelve walls are irredundant, the 270
elliptic rank-eight edges have exactly two endpoints, all 20 connected affine
components have complete cusp neighbourhoods, and there are no failures.
Their common weighted root Gram is exactly Kondō--Shimada Type II.

The determinant-one branch is rejected by 24 explicit hyperideal extreme
rays.  Noncompactness is not being confused with finite volume: this branch
has eight ideal rays as well, but the hyperideal rays make its volume
infinite and stop the surface compiler.

## 3. Complete cusp lift and same-P return

For either finite branch, lift all seventeen primitive ideal rays by the
one-boundary formula, then add the characteristic boundary.  Both integral
markings give one norm-coloured full-wall isomorphism class:

```text
dimension / facets          H10 / 30
wall split                  1 boundary + 12 roots + 17 cusp lifts
finite / ideal / hyper      1,333 / 207 / 0
Vinberg elliptic edges      10,030
connected affine audits     436
bad edge / cusp audits      0 / 0
```

The thirty walls span `I_(1,10)` with index one, the boundary is the unique
characteristic wall, and the full Gram is literally Type II.  This activates
the already exhaustive `raw_30` compiler: its split-root anti-Gram census,
half-sum saturation, primitive equivariant K3 glue, special period,
irreducibility, and hidden-curve exhaustion construct a one-boundary strong
Coble surface and return exactly the same 30-wall polytope.  The two E10
extensions are two integral marked branches of one returned full-wall
isomorphism class.

## 4. Theoretical meaning and strict boundary

This is a genuine `P -> finite lattice packages -> exact surface -> same P`
calibration obtained without selecting Type II.  Combined with the earlier
thirteen-edge theorem, both known one-boundary Coble types now arise as the
unique finite-volume survivors of two independently stated P-only lanes.

It is still not the requested flagship theorem.  First, the lane excludes
degree-above-three graphs and every fifteen-edge pairing-two graph.  Second,
the Coxeter enumeration itself eliminates the other root presentation before
algebraic geometry is used, so Kondō's classification does not create the
polyhedral exclusion.  The new theorem is therefore a scoped exact inverse
calibration and search-boundary result, not an AG-derived bound for all of
`D_1_split`, not a compact result, and not a new abstract polytope.

The remaining strict route is to exhaust the entire predeclared nodot
denominator (or another natural P-only subdenominator whose every member is
forced to have a geometric incidence branch), then use the independent Coble
classification to push back a theorem unavailable from the input Gram alone.

## Replay

```bash
PYTHONPATH=src:. python3 scripts/certify_coble_nodot_15edge_lane.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_coble_nodot_15edge_lane.py
PYTHONPATH=src:. python3 scripts/audit_coble_nodot_15edge_lane_hostile.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_coble_nodot_15edge_lane_hostile_audit.py
```
