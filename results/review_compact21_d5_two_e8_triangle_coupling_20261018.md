# Compact arithmetic rank21: D5 triangle incompatibility and forced E8 faces

18 October 2026. **New exact local lemmas, with independent arithmetic checks.**
At a compact arithmetic H21 vertex of type D5+E8+E8, at most one of the
three potentially triangular nonadjacent D5 faces can be a triangle.
An E8 triangle of either surviving type forces respectively a simplex
three-face or a simplex four-face. If both types occur in one E8 factor,
their third walls are orthogonal and force a seven-face combinatorially
isomorphic to the product of those simplices.

These are incidence restrictions. The upper bound on triangles among the
48 nonadjacent marked pairs remains **nine**. No face-charge saving,
arithmetic ceiling20, or new geometric example is claimed.

The full 148-line
[staged producer](../scripts/fresh_compact21_d5_two_e8_triangle_coupling_20261018.py),
its [saved output](fresh_compact21_d5_two_e8_triangle_coupling_20261018.json),
the immediate finite-tree classifier, exact ordered-field/Gram arithmetic,
and the triangle-transfer, universal E7 triangle, and compact21 field
reports were read. The source-level review passes with the scope below.
No complete independent replay of its 27 triangle domains was performed
here. The new tiny
[standalone checker](../scripts/check_compact21_e8_triangle_bundles_20261018.py)
imports no producer and independently checks the additional lemmas.

## Actual geometry and completeness of the staged nonadjacent domain

Use root square two. The old vertex roots span a positive-definite space
T in signature (21,1). If marks a,b specify a triangular two-face with
third wall c, the old roots with a replaced by c and with b replaced by c
are actual positive-definite corners. A whole untouched E8 cannot have a
connected finite extension, so c is orthogonal to it.

For D5 marks this confines c to the old D5 projection and the timelike
vertex direction. For E8 marks it leaves a possible D5 spectator contact;
the staged domain includes the whole D5 alongside that E8, and a second
whole E8 is forced orthogonal. D5 contacts were not silently excluded.

Every c contact is in one of the two finite corners. A label at least six
would have to be an isolated I2(m) component in that finite corner. At an
unmarked old root this would require all old neighbors to be contained
both in {a} and in {b}, impossible in these connected old components.
At a marked root it would require its only old neighbor to be the other
mark, whereas the marks are nonadjacent. Thus the complete contact domain
uses only labels two through five.

For the first corner, a finite connected diagram is a tree with degree at
most three. The new root meets at most one root in each old component;
otherwise it creates a cycle. With degree three all contacts have label
three. These facts justify the staged enumeration of supports of size
zero through three, the distinct-component test, and its label choices.
Testing the second finite corner then gives every candidate for the two
actual corner duties. There are six nonadjacent pairs in D5 and 21 in E8:
27 quotient domains, or 6+2*21=48 actual pairs.

For each positive-height candidate, every connected principal subset
containing c is checked for a singular positive-semidefinite obstruction.
A subset not containing c is positive definite already. A disconnected
affine obstruction has a connected component of this kind. Skipping
subsets containing a label-five bond is valid by the affine Coxeter
classification, whose connected diagrams have no such bond. A positive
kernel produces a nonzero future null vector obeying all chamber
inequalities, and contradicts compactness.

The finite and affine classifications, compact acute simplicity, and the
Vinberg arithmetic criterion are explicit theoretical inputs. They are
not newly proved by a finite program in this packet.

## Field and cyclic-product gates

The completed internal
[compact21 field reduction](fresh_compact21_single_field_20261016.md)
leaves exactly k=Q(sqrt5). This is a separate theorem; the three-field
list valid only from dimension22 is not being applied at21. The older
hostile report named in the staged manifest predates this final
single-field reduction. The new checker pins the final reduction too.

Every simple cycle in the old forest plus c consists of c and the unique
old path between two contacts. The old path has only label-three bonds,
so its cyclic field membership is exactly that of the two endpoint
contact products. This justifies the staged cycle check. With that check
passed, H=c-projection-squared minus two lies in k. Under the nonidentity
real embedding of k, the unchanged rational old Gram is positive definite;
the conjugate Gram's Schur complement requires H-conjugate at most zero.
This is a necessary arithmetic condition, not sufficiency for a chamber.

Write s=sqrt5 and phi=(1+s)/2. The surviving D5 contacts and squared
heights are

|Name|Nonzero contact magnitudes|Squared height H|
|---|---|---|
|C0|phi at0|(s-1)/2|
|L3,L4|phi at3 or4|(-1+5s)/8|
|U3,U4|1 at1 and phi at3 or4|(23+13s)/8|

The only potentially triangular nonadjacent D5 pairs are (0,2), (1,3),
(1,4). The first uses C0; the second uses L3 or U3; the third uses L4 or U4.
In either E8 the only two surviving third walls are c0 and c6, with contact
phi at root0 or root6 and heights H0=4+2s, H6=1+s. The only potentially
triangular nonadjacent E8 pairs are (0,2), (3,5), (3,6), (4,6).

For roots confined to different old components the squared mutual product
is the product of their heights. The ring of integers of k is Z[phi];
a+bs belongs to it exactly when 2b and a-b are integers. For a fork root
L or U paired with either E8 type, the four values are

  (23+9s)/4,  3+s/2,  (111+49s)/4,  11+9s/2.

None is integral. They are length-two cyclic products in the normalized
Cartan matrix, so the arithmetic criterion excludes all eight literal
fork/E8 root pairings. C0 instead gives 3+s and2, which pass this necessary
test. The review found no normalization or ground-field error in this step.

## New lemma: at most one D5 triangle

Let e be the future unit timelike vector at the common old vertex. A
specified contact vector determines c=p+sqrt(H)e uniquely, since its
projection p and positive height are fixed. For any two distinct candidate
walls, their mutual product is

  c_i dot c_j = p_i dot p_j - sqrt(H_i H_j).

The standalone checker verifies the displayed rational D5 inverse by
multiplication with its Cartan matrix. It then checks all ten unordered
pairs of the five candidates above. In every case p_i dot p_j is positive
and its square exceeds H_i H_j. Every such pair therefore has strictly
positive mutual product, forbidden for distinct outward walls of an acute
Coxeter chamber.

For example the low/low and high/high pairs on different forks have
mutual product (5-s)/4>0. For a mixed pair on different forks,

  p_i dot p_j=(13+7s)/8,
  (p_i dot p_j)^2-H_i H_j=(7+5s)/4>0.

For C0 against a low or high fork candidate, the corresponding positive
differences of squares are 3(s-1)/4 and (7+5s)/4. Thus the argument excludes
the long/fork combinations too, not only simultaneous fork triangles.
Because the three marked D5 faces use disjoint candidate sets, at most
one of these faces is a triangle. No assertion is made about adjacent D5
marked pairs, or about D5 beside an A1 or other unclassified spectator.

## New lemma: E8 simplex bundles and a product face

The checker verifies two exact inverse columns for E8. In particular
J00=4, J66=2, J06=2. The projections of c0,c6 have product
2*phi^2=3+s and

  H0 H6=(3+s)^2.

Both heights and 3+s are positive, so c0 dot c6=0 whenever both actual
walls occur. This identity proves compatibility of this particular mutual
product; it does not prove a complete chamber exists.

For c0, exactly deletions 0,1,2 from the old E8 give positive-definite
eight-root corners. For c6 the corresponding deletions are 3,4,5,6.
All these corners remain positive after adding the untouched old factors.
In a compact acute chamber any positive-definite principal wall set has
an actual common face: project a strict interior timelike vector off its
positive span using the entrywise nonnegative inverse, which preserves
the strict inequalities of the remaining walls. Rank21 corners therefore
give actual vertices, not merely abstract finite diagrams.

Consequently c0 forces all three triangles on pairs from {0,1,2}, and c6
forces all six triangles on pairs from {3,4,5,6}. The respective supporting
three-face and four-face have all vertices and edges of a simplex. At
every such vertex these edges exhaust the degree of the simple face;
connectedness of a convex polytope's edge graph excludes further vertices.
The faces are therefore actual simplices.

If both walls occur, retain original root7 and all old roots outside that
E8 as supporting walls. The resulting actual face has dimension seven.
The checker verifies all20 corners obtained by omitting one member of
{0,1,2,c0} and one of {3,4,5,6,c6}. Their product-simplex adjacency gives
seven edges at every vertex, again saturating the degree of this simple
seven-face. The connected-edge argument gives the whole face, with the
incidence structure of Delta3 times Delta4.

This seven-face is **not a Coxeter polytope** in its intrinsic metric.
Projecting off root7 changes the square of root2 to3/2 and leaves root1
of square two; their product stays -1. The resulting dihedral cosine is
1/sqrt3, strictly between cos(pi/3) and cos(pi/4). Thus classifications of
few-facet Coxeter polytopes cannot be applied to this product face.

## Consequences, limitations and reproduction

Within each E8 the four designated nonadjacent triangles occur in bundles
of one and three, so their number is 0,1,3 or4. A D5 fork triangle excludes
all E8 designated triangles and, by the new D5 lemma, all other D5
designated triangles: that branch has at most one. If the long D5 triangle
occurs, both E8 factors can still contribute four in the necessary model.
The overall upper bound is therefore still 1+4+4=9. Reaching nine forces
the long D5 wall and both orthogonal E8 third walls in each E8 factor, hence
two of the product seven-faces above.

The existing charge model pays at least two on a triangle and at least one
on a quadrilateral. A nontriangle may still be a charged quadrilateral.
Neither this incompatibility nor the simplex bundles establish a new
quad payment, a saving in the D5+2E8 vertex allowance, or a cut for the
different supported forest A1+D5+E7+E8. No enlarged exact LP was run and no
claim that the saved 40.2697 dual is broken follows. A useful next target
is a compatible quadrilateral or edge-endpoint classification coupled to
the forced simplex walls, with a literal test of any new charging column.

The independent tiny checker ran twice as separate Python processes in
the sequential batch ending in chunk `87bd0c`, exit0/PASS. It checks 20
product-face corners, the two E8 inverse columns, all 16 single-wall
deletions, the D5 inverse and all ten D5 candidate pairs. The batch also
replayed the separate restricted RA5 result; all three tiny commands
together took under one second. No Lean, optimizer, attachment queue,
large file load, or additional proof process was started; no temporary
proof files or bytecode caches were produced.

```sh
python3 -B scripts/check_compact21_e8_triangle_bundles_20261018.py
```

The checker pins the staged source/output and completed compact21 field
reduction. Its source SHA-256 is
`8b2af011ec40603b69da1bbf4329605c4bec7197e7e497c7a987fc6931808aa1`.
The new lemmas are written geometric proofs plus exact independent
arithmetic; they are not Lean-verified or externally reviewed. The prior
27-domain enumeration remains an inherited internally checked component,
not a fresh independent full census performed by this review.
