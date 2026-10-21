# Two opposite-pair families over Q(√2): complete pendant-tree classification

This is a local lemma for the compact arithmetic dimension-bound programme.
Under the topology specified below, the complete quadrilateral-star rank and
arithmetic conditions leave exactly three local stars. There is **no global
dimension bound and no finite-volume polytope claim** in this package.

The [producer](../scripts/fresh_sqrt2_dotted_tail_20261015.py) and its
[exact certificate](fresh_sqrt2_dotted_tail_20261015.json) check both symbolic
rank identities, all surviving physical and conjugate signatures, all four
corners, both opposite extensions, and 509 principal subsets for affine
obstructions. Run `python3 scripts/fresh_sqrt2_dotted_tail_20261015.py --verify`.

## Precise scope

Let the cyclic black vertices of a quadrilateral star be a,b,y,d. The white
core is u,v,w. The induced diagram on b,u,v,w,d is a five-cycle, with b–d
label 4 and its other four edges label 3. There is an edge a–v of label 3.
The contact a–y is to be determined. In the single-attachment family y–v
is right-angled; in the double-attachment family it has label 3. Every
other unlisted core pairing is right-angled.

Each additional connected white component is a spherical tree attached
once to exactly one core vertex, or is orthogonal to the whole core.
Attachments to black vertices are included. Additional core edges, shared
vertices, and external components attached to two or more core vertices are
outside this lemma. In particular, the finite core census alone does not
prove these topology assumptions for every quadrilateral.

We impose the necessary conditions for a compact arithmetic Coxeter
quadrilateral star with ground field Q(√2): its four corner diagrams are
spherical; it contains no affine principal subdiagram; its physical Gram
matrix has exactly one negative direction; its conjugate is positive
semidefinite; its normalized cyclic products are algebraic integers; and
its full star Gram matrix has corank one. Orthogonal white components may
be eliminated and do not change the corank condition on the displayed block.
The finite labels 2,3,4,6,8 and the possible minimal indefinite pieces are
recorded in the [Lannér catalogue](fresh_sqrt2_lanner_catalogue_20261015.md).

## Exhausting the permitted attachments

The corners containing a,b and a,d contain D5. A connected spherical
extension of D5 is simply laced, is of type D or E, and has at most one
branching vertex. Hence no extra component can attach at v: a four-armed
star would occur. An attachment at u or w creates a second branching vertex
in one of these two corners. Attachments at a,b,d must therefore be single
simply laced paths, with at most one path at each vertex.

A first extra vertex at b produces the affine fork with central vertex b
and arm labels 4,3,3, on b,d,u and that new vertex. The same argument at d
uses d,b,w and its new vertex. Compactness excludes both. This fork's unit
Gram matrix has inertia (3,0,1); the certificate supplies its exact positive
kernel vector.

The remaining path at a has length l∈{0,1,2,3}. Its corner is respectively
D5,E6,E7,E8. Length 4 gives affine E8, and any longer path contains it.
In the double family the same argument bounds the path at y to a length
m∈{0,1,2,3}. This proves exhaustion of all pendant-tree shapes in that family.

In the single family the finite component containing y is separate from
the C5 corner. After eliminating its white vertices, write its reduced
unit diagonal as q, where 0<q≤1. Except for I2(8), every possible such finite
tree has rational squared edge contacts, so leaf elimination gives rational
q. This includes the isolated I2(6) case. An 8-edge cannot occur in a larger
spherical connected diagram; its sole exception gives

    q = 1 − cos²(π/8) = (2−√2)/4.

For rational q, equality q=1 holds precisely when y has no attachment:
every nonzero attachment subtracts a strictly positive Schur complement.

## Rank and arithmetic in the single family

Normalize all original roots to squared length 2, and let the a–y Gram
entry be −r. Eliminating the length-l path changes the unit a diagonal to

    p = (l+2)/(2(l+1)).

The full remaining seven-by-seven determinant is

    2(1+√2) [r² − q(4p−1+√2)].

Thus corank one requires

    D := r² = q((l+3)/(l+1)+√2).

Here D is the two-cycle product in the ORIGINAL normalized Gram matrix;
Vinberg's integrality condition requires D∈Z[√2]. If q is rational,
its √2 coefficient is q, so 0<q≤1 forces q=1. The rational coefficient
(l+3)/(l+1) is integral only for l=0,1. For the exceptional I2(8) value of q,
all four values l=0,1,2,3 give nonintegral D, as checked exactly in the data.
Consequently the only single-family possibilities are:

| a-tail length | D=r² | Nature of a–y | Affected complete corner component |
|---|---|---|---|
| 0 | 3+√2 | dotted | D5 at two corners |
| 1 | 2+√2 | label 8 | E6 at two corners |

There is no y-tail in either case. The second survivor therefore belongs
to the finite-label side of the complete census, although it arises from
the same rank identity.

## Rank and arithmetic in the double family

Set p=(l+2)/(2(l+1)), q=(m+2)/(2(m+1)), and t=(1−√2)/2.
The determinant is 2(1+√2) times

    r² + 2tr − 4pq + 2t(p+q).

The two fixed edges in the triangle a–v–y have normalized contact 1.
Its cycle product therefore forces r∈Z[√2]. Write r=A+B√2 with integers
A,B. The √2 coefficient of the rank equation is

    2AB+B−A−(p+q)=0.

Thus p+q is an integer. Among p,q∈{1,3/4,2/3,5/8}, this happens only at
p=q=1. The rank polynomial then factors as

    (r−1−√2)(r+2).

The unique positive root is r=1+√2, giving a dotted a–y edge. There are
no tails. This arithmetic shortcut was independently identified during
review; the producer checks the same conclusion by exact field square-root
and integrality tests for all 16 pairs (l,m). All four affected complete
corner components in this survivor are D5.

## Exact local checks and significance for weights

For each of the three survivors the complete core-and-tail Gram matrix
has physical inertia (N−2,1,1) and conjugate inertia (N−1,0,1), where N is
its number of roots. All four indicated corner matrices are positive
definite, and both opposite-pair extensions are indefinite. None of the
509 nonempty principal subsets tested across the three matrices is positive
semidefinite and singular. For the single family the calculation rescales
y by √D, giving diagonal 2D and a–y entry −D, to work in Q(√2). Since D
and its conjugate are positive, this preserves the required signatures.
This rescaling is for calculation only; the integrality argument above
uses the original normalized cycle products.

These algebraic checks also realize the four finite corners in one
hyperbolic cone. The irreducible Gram matrix has nonpositive off-diagonal
entries; a positive eigenvector for its negative least eigenvalue gives
a strict timelike interior vector. Projection onto the orthogonal
complement of any of the positive-definite corner supports preserves
timelikeness and satisfies all other wall inequalities strictly, because
the inverse of a positive-definite matrix with nonpositive off-diagonal
entries is entrywise nonnegative. The four resulting finite rays bound
the indicated compact quadrilateral. This does not imply that the whole
ambient chamber has finite volume.

For the ongoing dimension-bound argument, the complete D5 corners have
the marked long-leaf/short-leaf pair, and the E6 corners have a marked
near node on one long arm and the far node on the other. Half-unit weights
on these two pairs per complete D5 or E6 component pay at least one unit
on every surviving quadrilateral in these families. Their total cost is
one per component. A global weight theorem still requires all other
quadrilateral cores and all external components with multiple attachments
to be covered. No claim that this has been achieved is made here.
