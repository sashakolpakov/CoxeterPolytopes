# Integral overlattices of the A11+A2 face lattice

The normalized E7 face has rational form `A11 + A2 + <-2>` and determinant
`-72`. This note closes one limited way of selecting a reflection host for
that face. It is not a nonexistence theorem for a finite-volume H13 source.

Its discriminant group is `Z/12 + Z/3 + Z/2`, with self-pairing

`11a²/12 + 2b²/3 - c²/2 (mod Z)`.

The [exact enumeration](fresh_h13_a11a2_lattice_gate_20261014.json) checks all
72 classes. Only `(0,0,0)` and `(6,0,0)` have integral self-pairing. Every
integral overlattice lies in the dual lattice, so its quotient must be a
subgroup of these two classes. Both subgroups do work. Consequently there
are precisely two integral overlattices containing this fixed lattice:

* the original even lattice of determinant `-72`;
* the odd index-two extension of determinant `-18`, obtained by adjoining
  the sixth fundamental weight of A11, of norm 3.

There is no proper even integral overlattice and no integral index-three
extension. The [checker](../scripts/fresh_h13_a11a2_lattice_gate_20261014.py)
also realizes the extension in sum-zero Euclidean coordinates and computes
its integral Gram matrix and determinant independently of the discriminant
formula.

For the two A2 simple roots `alpha,beta`, the vector `alpha-beta` has norm 6.
Its reflection swaps `alpha` and `beta` and fixes the remaining summands,
including the adjoined A11 weight. It therefore preserves both lattices.
Its mirror bisects the interior of the A2 chamber. Using the full integral
reflection group of either lattice cannot keep the original `A11+A2`
apex chamber; there is already an extra mirror through its apex.

This leaves restricted reflection subgroups and different lattice models
open. In particular, a geometric transfer of a known chamber may adjoin a
square root in the timelike direction and need not be rational over the
fixed normalized face form. The lattice gate must not be used to discard
that possibility. Neither reflectivity nor a complete hyperbolic chamber
is asserted here.

Replay: `python3 scripts/fresh_h13_a11a2_lattice_gate_20261014.py --verify`.
