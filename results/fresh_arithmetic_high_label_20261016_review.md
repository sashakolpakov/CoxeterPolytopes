# Independent review of the eight-field compact dimension bound

The [new proof](fresh_arithmetic_high_label_20261016.md) passes an independent
mathematical reading and [standalone exact verification](../scripts/fresh_arithmetic_high_label_20261016_review.py).
Its [certificate](fresh_arithmetic_high_label_20261016_review.json) imports
no producer code, project arithmetic, symbolic algebra or earlier graph engine.
This remains an internally reviewed research proof candidate; external
validation and priority are unresolved.

The theorem is: a compact arithmetic Coxeter polytope has dimension at
most 20 whenever its ground field K is not a subfield of
F=Q(sqrt2,sqrt5) and does not contain K7=Q(cos(2pi/7)). Among the twelve
fields in Vinberg's dimension-at-least-14 list, eight satisfy this
criterion. Together with the earlier Q(sqrt2) and K7 arguments, this
leaves **Q(sqrt5) and Q(sqrt2,sqrt5)** as the only unresolved fields for
compact dimension 21. This does not exclude that dimension.

## Independent finite checks

For each of the eight fields, the verifier works in an explicit cyclotomic
compositum and supplies two automorphisms:

1. one fixing both sqrt2 and sqrt5 while moving K, disproving K contained in F;
2. one fixing K while moving K7, disproving K7 contained in K.

The quadratic characters modulo 8, 12 and 5 describe sqrt2, sqrt3 and
sqrt5; products describe sqrt6 and intersections describe the biquadratic
field. For real cyclotomic fields, fixing K_m means the residue is plus or
minus one modulo m. Degrees are checked as subgroup indices. This method
is separate from the producer's squareclass, quotient-group and degree
arguments. It checks all eight applications exactly, including the
potentially ambiguous pair of degree-three fields K7 and K9.

The local graph checker enumerates **35,904 complete binary graphs** with
selected order-two or order-three opposite minimal diagrams. This differs
from the producer's generation by connected pieces, distinguished bonds
and connecting edges. It checks all shared-white cases, every cross contact,
finite-corner forest conditions, distance-one/two marked paths, and isolation
of distinguished bonds. Exactly **nine supports** remain and agree literally
with the producer. In each, both minimal diagrams are triangles on separate
white nodes, the distinguished bonds are the opposite black pairs, and the
only connection is the white-white bond.

Independent Fraction checks establish the strict thresholds
cos(pi/9)^2>7/8 and cos(pi/6)>3/4, the two Schur-cost inequalities, and the
dimension-21 h-vector identity. Path and fork counts are reconstructed from
literal graphs; the all-rank formulas follow from counting edges and pairs
of neighbors. The two parity polynomials for all dimensions above 21 are
derived coefficientwise, not inferred from a finite dimension sample.

## Mathematical scope check

The essential field step uses **equality** of the Lanner determinant field
and the ambient arithmetic ground field. A negative physical determinant
whose other ground-field conjugates are positive cannot lie in a proper
subfield: another embedding fixing that subfield would preserve its
negative value. The arithmetic conjugate-sign criterion and determinant
field statement are the published inputs identified in the primary report.

A minimal indefinite diagram of order at least four cannot contain a
label at least six: connectedness supplies a neighboring bond, and the
resulting proper three-node principal diagram fails positive definiteness.
Its entries therefore lie in F, contradicting the field hypothesis. No
finite rank cutoff or complete Lanner catalogue is needed for this step.

For triangles, the field hypothesis excludes labels at most five and
excludes label seven. A right triangle with all labels at most eight has
determinant in F, so is excluded too. Thus every triangle has either a bond
at least nine or a bond six/eight with nonright contacts at both endpoints.
The dense alternative is necessary, in particular over K16.

In a quadrilateral with no marked path of length at most two, the complete
support argument puts these bonds on the opposite black pairs. Their
four-node Gram is positive definite. Eliminating it reduces each of the
two selected white diagonals from one by a cost strictly greater than two.
The off-diagonal is unchanged and has square strictly less than one.
Therefore the residual matrix is negative definite. This creates two
negative directions in a **principal** submatrix, contradicting Lorentz
index one. Additional white nodes cannot repair it; there is no inherited
attachment exhaustion, full-star rank hypothesis or spectator assumption.

The triangle argument also retains the white support: every minimal
indefinite subset of the triangle star must contain all three black nodes,
since each proper black subset with the whites is finite. Its order is at
most three, so it is precisely the black triangle. The distinguished bond
is isolated in its finite corner and receives weight two. Quadrilaterals
receive at least one; larger faces satisfy the required inequality by
nonnegativity. Every vertex has weight at most 2n-2. At n=21 the h-vector
forces a strict average greater than 40, contradicting the vertex bound 40.

The standard compact face correspondence, compact simplicity, absence of
affine principal diagrams, arithmetic conjugate-sign criterion, finite
Coxeter classification, and h-vector symmetry/nonnegativity remain explicit
mathematical inputs. The review found no gap in their use here. It does
not formally verify those external theorems or establish novelty.

```sh
python3 scripts/fresh_arithmetic_high_label_20261016_review.py
python3 scripts/fresh_arithmetic_high_label_20261016_review.py --verify
```
