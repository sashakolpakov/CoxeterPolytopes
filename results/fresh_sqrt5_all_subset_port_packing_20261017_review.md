# Independent source review of all-subset port packing

**PASS; FINAL/FROZEN.** This is an independent mathematical and source
read of the [all-subset package](fresh_sqrt5_all_subset_port_packing_20261017.md).
No inconsistency was found. It is not an independent repetition of the
17,828,080 literal arithmetic tests. The primary generation and fresh
saved replay passed separately.

The frozen [source](../scripts/fresh_sqrt5_all_subset_port_packing_20261017.py)
has SHA-256 `103f80bbfc08f08cfa96cdcbc71c4abe9c2a1fbf1e4caa7dc9c5662d418ddfc2`;
the [JSON](fresh_sqrt5_all_subset_port_packing_20261017.json) has SHA-256
`01fecabbe6b8c0621205cfcc819106b21e64b9f08a1aa5396596b6c24a877c78`.

For every automorphism g of the complete labelled vertex forest, the
images of all actual maximal-simplex port blocks remain disjoint
simultaneously. Therefore their integer subset-capacity inequality
survives averaging over this one group. Orbit-stabilizer and linearity
make that average identical to averaging each block uniformly over its
own orbit. No symmetry of the chamber or independence of its blocks is
assumed. Invariant orbit codes and literal block masks implement this
normalization.

I checked the identities, for intersection sizes k=0,...,4,

`floor(k/2)=C(k,2)-2C(k,3)+4C(k,4)`,

`floor(k/3)=C(k,3)-3C(k,4)`, and `floor(k/4)=C(k,4)`.

The subset-zeta loop sums features on subsets of each mask, which is
the direction these formulas require. Common-denominator integer loads,
the hypergeometric aggregate for each subset size, exact cross-multiplied
load/capacity comparisons, and the direct recount at each maximizing
mask are consistent. The reconstruction includes every positive simplex
row, including formerly unexported classes, and compares its complete
corner block with both inherited pair and one-port coefficients.
Actual quadrilateral-provenance copies are correctly skipped for this
new simplex-only charge.

Only vertex masses increase. The exact face and provenance masses remain
fixed, so all inherited signed loads remain unchanged and positive-cost
constraints can only improve. The resulting12.818614535... value is a
fixed-alpha lower bound within the enlarged partial payment model. It is
not an optimum, a feasible primal rule or an ambient dimension bound.

The148 exported maximizing violated charges reproduce this particular
mass repair. They are correctly distinguished from the complete
all-subset family and from the still stronger full disjoint-block packing
polytope. Assigning these charges zero payment on unexported classes is
a valid weaker model only when that restriction is stated explicitly.
The compact disjoint-port hypothesis is retained; no unrestricted
finite-volume transfer on ideal-ended ports is asserted.
