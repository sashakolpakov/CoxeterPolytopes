# Exact principal-core filters over Q(√5)

The [checker](../scripts/fresh_sqrt5_quad_filters_20261015.py) applies the
arithmetic field and signature conditions to all 1,180 marked classes in
the [complete union census](fresh_sqrt5_quad_exhaustion_20261015.md).
The [certificate](fresh_sqrt5_quad_filters_20261015.json) records every
decision. A [separate review](fresh_sqrt5_quad_filters_20261015_review.md)
independently reconstructs the squareclass constraints and full inertias.
Generation and saved replay pass. No dimension bound follows from these
principal-subdiagram conditions alone.

| Outcome | Classes |
|---|---:|
| Two physical negative directions | 617 |
| Cyclic product outside Q(√5) | 89 |
| Conjugate Gram not positive semidefinite | 78 |
| Finite-label cores passing these necessary tests | 362 |
| Dotted-distance families, not yet tested here | 34 |

The 362 finite-label survivors consist of 124 singular and 238 nonsingular
cores. Passing a core does not prove that a full quadrilateral star, or an
ambient compact polytope, exists.

## A smaller exact field by rescaling roots

Write φ=(1+√5)/2. In norm-two normalization the nonzero finite contacts
are 1,√2,φ,√3,√(φ+2), for labels 3,4,5,6,10. The squareclasses of
2,3,φ+2 are independent over Q(φ), as proved in the independently reviewed
[Lannér catalogue](fresh_sqrt5_lanner_catalogue_20261015.md). Mark each edge
by its three-bit radical class. A cyclic product lies in Q(φ) precisely
when the xor of its edge classes is zero. A breadth-first vertex coloring
therefore tests the whole graph: edge(i,j) must have color cᵢ xor cⱼ.
A conflict is an explicit forbidden cycle, regardless of any later walls.

When this condition holds, let Rᵢ be the product of the positive radicals
whose bits occur in cᵢ, and rescale root i by Rᵢ. The new diagonal is
2Rᵢ². Every off-diagonal product −cᵢⱼRᵢRⱼ lies in Q(φ): each radical
now occurs an even number of times. The checker calculates it directly
using the radicands 2,3,φ+2, without building a degree-sixteen entry field.
All entries are integral in Z[φ]. For every edge it also checks the original
normalized two-cycle 4Gᵢⱼ²/(GᵢᵢGⱼⱼ).

Both conjugates of each radicand are positive. Thus the rescaling preserves
physical inertia, and its conjugate preserves the inertia at the other
ground-field embedding. The sign of φ itself is conjugated correctly;
it is not replaced by its absolute value. The conjugate scaled matrix is
diagonally congruent to an allowed conjugate of the original Gram matrix.
This justifies testing conjugate positive semidefiniteness in Q(φ).

The field operations use rational coefficients a+bφ, with φ²=φ+1 and
conjugation a+bφ ↦ a+b−bφ. Signs are determined exactly after rewriting
twice this expression as (2a+b)+b√5. No numerical eigenvalue or sign
tolerance enters any decision.

## Signature and rank data

For each finite-label core, eliminate the positive-definite corner
W+{0,1}, leaving the two opposite omitted black nodes. The physical
two-by-two Schur matrix has two negative diagonal entries because the
corresponding opposite Lannér subsets are present. A positive determinant
therefore gives two negative eigenvalues and excludes the core in any
Lorentzian extension. A nonpositive determinant leaves at most one.

At the other embedding, the corner is positive definite. Both residual
diagonals and the determinant must be nonnegative. The exact corner
pivots, residual matrix and signs are saved for both embeddings. A zero
physical residual determinant marks a singular core; it does not by itself
exclude extensions or certify a complete face.

Finite-contact cycle products that lie in Q(φ) are algebraic integers,
because every finite contact is an algebraic integer. Thus a separate
integrality rejection is unnecessary in this finite-label stage. It remains
essential when solving the continuous dotted-distance families.

The [following contact analysis](fresh_sqrt5_quad_contacts_20261015.md)
uses compactness, already-sufficient H4 corner weights, and all one-extra-white
profiles. Arbitrarily many external white nodes and continuous distances
require further arguments beyond the present finite principal-core filter.

Replay: `python3 scripts/fresh_sqrt5_quad_filters_20261015.py --verify`.
