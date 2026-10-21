# Triangle payments with finite distance8 and ideal distance5

Every triangular face of a finite-sided finite-volume hyperbolic Coxeter
polytope receives total payment at least2 under the following rule:

* a finite corner pays1 when its marked roots have distance at most8;
* an ideal corner pays1 when its varying affine factor has order at most9,
  or its marked distance is at most5, or its marked distance is at most8
  and a mark is within5 of a branch or a nonsimple-bond endpoint.

This includes nonsimple ideal vertices and arbitrary additional supporting
white roots. The [checker](../scripts/fresh_ideal5_finite8_triangle_payments_20261017.py)
and [certificate](fresh_ideal5_finite8_triangle_payments_20261017.json) enumerate
the finite cases and their exact E8 data. This is a triangle theorem, not
a quadrilateral theorem or a dimension bound.

## Complete core domain

The actual triangle setup and minimal-indefinite-core reduction are inherited
from the [preceding triangle proof](fresh_ideal6_triangle_payments_20261017.md).
Such a core contains the three black roots, with all other roots supporting
the face. Deleting any black root gives an actual corner subdiagram and
therefore a positive-semidefinite Gram. Graph paths avoiding that third
black root remain corner paths when all supporting whites are restored.
Likewise a branch or nonsimple bond already present in the marked component
persists. The [actual cusp-factor theorem](fresh_cusp_factor_payment_framework_20261017.md)
places those paths and features in the full varying affine factor at an
ideal corner.

For core order at most6, a Steiner tree through the three marks has at
least two marked-pair paths avoiding the third mark, each with at most5
edges. This argument allows all labels in the small-core families.

For orders7,...,10, the complete inherited high-core catalogue gives
1,025 marked triples:105 at order7,224 at order8,336 at order9 and360 at
order10. The checker computes each marked path and every feature distance
in the actual corner subdiagram. All but four triples have at least two
elementary payments. Each exception is in the simply-laced10-node L25 tree,
whose edges are

    0-1, 0-2, 2-3, 0-4, 4-5, 5-6, 6-7, 7-8, 8-9.

The exceptional triples are(0,1,9),(0,2,9),(1,2,8),(1,2,9).
Each already has one elementary payment; the next argument supplies a
different corner payment. No unexamined supporting-white extension is
discarded.

## Three positive-kernel resolutions

For the triples containing9, the corner formed by the other two marks
contains the full affine E8 diagram on nodes0,...,8. Its positive kernel is

    (6,3,4,2,5,4,3,2,1).

The source verifies this kernel and a positive-definite principal8-node
matrix, proving positive semidefiniteness. Any additional supporting white
belongs to this actual corner. Its nonpositive products with a positive
affine kernel must all vanish. Thus every extra supporting white is
orthogonal to all nine displayed nodes.

At the indicated alternate corner a simple leaf and its degree-two neighbor
are consequently fixed, including after every supporting white is restored:

| Triple | Alternate corner | Fixed simple leaf | Fixed neighbor | Marked distance |
|---|---|---:|---:|---:|
|(0,1,9)|(0,9)|3|2|6|
|(0,2,9)|(0,9)|1|0|6|
|(1,2,9)|(1,9)|1|0|7|

The fixed leaf belongs to the marked connected component; it need not
itself be a mark. None of the classical irreducible affine A/B/C/D diagrams
has a simple leaf whose neighbor has degree two. Hence this actual corner
is finite, in which case the displayed distance is at most8, or the varying
factor is exceptional affine, of order at most9. Either case pays1.

## The finite-E8 resolution

For the remaining triple(1,2,8), the corner(1,2) contains finite E8 on
nodes0,...,7, together with node9. Write its normalized Gram as E. The exact
inverse is entrywise nonnegative, with diagonal

    (30,8,14,4,20,12,6,2).

An added supporting root has square2. Its nonnegative contact magnitudes
x with E must satisfy x^T E^(-1)x<=2, because the augmented principal
diagram belongs to a positive-semidefinite actual corner. Every nonzero
native contact has magnitude at least1. Nonnegative inverse entries imply
that any contact on nodes0,...,6 alone already costs more than2. Thus all
additional supporting whites are orthogonal to those seven nodes. Contacts
at7 are allowed; the argument does not assume they vanish.

At the alternate corner(1,8), node1 remains a simple leaf and node0 remains
its degree-two neighbor. The marked path1-0-4-5-6-7-8 has length6. As above,
the entire marked component is finite or exceptional affine of order at
most9, so this corner pays. Arbitrary attachments farther along the path
cannot change its protected leaf pattern.

## Verification scope

The source verifies all1,025 triples, the exact four exceptional domains,
their distinct alternate payments, positive-definite corner pivots, the
affine E8 kernel, and every entry of the finite E8 inverse. The all-label
small-core argument, positive-kernel orthogonality, projection inequality
and affine-type exclusion are analytic parts of the proof. The catalogue
classification and actual-face/cusp reduction are explicit inherited inputs.

Generation and a separate fresh saved replay both terminated exit0/PASS.
A separate agent read the complete source and mathematical proof, including
the small-core Steiner argument, all four E8 mechanisms, native contact gap,
protection of the actual leaf pattern and actual-cusp transfer. That review
does not duplicate the1,025-case census or the inherited classification.
The parent also independently obtained the E8 inverse diagonal and affine
kernel using symbolic matrix arithmetic. This package is final/frozen.
No quadrilateral coverage, ideal-factor cost,
new ambient bound, external refereeing or formalization is asserted.

    python3 scripts/fresh_ideal5_finite8_triangle_payments_20261017.py
    python3 scripts/fresh_ideal5_finite8_triangle_payments_20261017.py --verify
