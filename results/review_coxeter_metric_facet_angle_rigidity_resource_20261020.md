# Independent review: metric Coxeter facets

20 October 2026. **PASS as a written mathematical argument; no flaw found.**
I read the entire proposed proof independently, before reading any root
review. No scientific program, census, checker, Lean module or build was
run. This is not external refereeing or an end-to-end formal proof.

Reviewed source:

* `results/coxeter_metric_facet_angle_rigidity_20261020.md`
* SHA256 `f688e520f3ff7ac43225aaacd89413c03f5dd3d2d008d7e4dcdcc48326ac83fb`

## Local geometry and the exceptional label

At an actual compact codimension-three face the three parent normals are
independent and positive definite. Orthogonal projection to the facet
normal's complement gives exactly the stated cosine

    (b+c1*c2)/sqrt((1-c1^2)*(1-c2^2)).

I recalculated every marked-node case of A3, B3 and H3. The resulting
squared cosines, in the source's order, are

    1/3, 1/9; 1/2, 2/3, 1/3;
    (5+sqrt5)/10, (3+sqrt5)/6, (5+2sqrt5)/15.

The reducible rank-three cases are also covered: marking the isolated
node retains I2(m); marking a node of I2(m) leaves the two projected
black nodes orthogonal. The real cyclotomic degree phi(k)/2 for k>=3
leaves exactly k=3,4,6 among rational squared cosines, and k=5,10 among
the quadratic Q(sqrt5) possibilities. The other quadratic possibilities
k=8,12 have different fields. Only1/2 survives the table. Its B3 case
has the marked node at the4-end; one oblique parent is at angle pi/4
with the chosen facet, and the other parent is perpendicular.

Thus the assertion for EVERY integer label other than2 or4 is proved,
not merely the originally suggested labels3 and5. Dotted pairs are
correctly excluded from this finite-corner assertion.

## Realized example and H8 application

For the raw [4,3,5] Gram I independently obtained the leading minors
4,4,4,4-4phi. Its proper principal diagrams are finite; its inverse
entries have the common negative orientation at the real embedding.
At the other embedding the leading minors are positive. Each displayed
root reflection preserves the stated Z[phi] lattice, including the
norm-four root. This gives an admissible first-type arithmetic compact
simplex. Eliminating its first root gives the displayed intrinsic Gram
with diagonal1,2,2 and labels4,5,2. The label4 exception is therefore
realized with the claimed arithmetic and compact scope.

I also read the literal relation matrix in
`data/known_examples/bugaenko_compact_h8.json`. The eleven listed finite
nonright edges are exact, and every one of its eleven nodes is incident
to a3- or5-edge. The claimed conclusion is only that any parent
extension of this METRIC facet has an orthogonal collar. Its existence,
boundedness and caps are not inferred from the local theorem.

## Whole-polytope class bound

The global hypothesis is that all actual metric facets are Coxeter; it
is not assumed for arbitrary compact Coxeter polytopes. Compact Coxeter
polytopes are simple. Any connected finite vertex component of rank at
least three contains an induced connected three-node subdiagram, hence
A3, B3 or H3. The subset of vertex walls defines an actual face, so its
bad marked-node projection is genuinely a corner of one of the metric
facets. This contradicts the hypothesis. Vertex diagrams are therefore
disjoint unions of A1 and I2(m).

The inheritance step is valid without invoking a general Coxeter-facet
theorem. Eliminating arbitrary white nodes from such a matching diagram
either retains a whole black I2 pair with its native angle, leaves a
single black node whose norm alone changes, or removes a component.
Different components remain orthogonal. All faces are compact metric
Coxeter polytopes with the same matching-vertex condition.

At rank9, each vertex contributes at most four nonright two-face
corners, one for each whole I2 component. Each triangle contributes at
least two, each quadrilateral at least one, and every other k-gon at
least5-k. Counting the36 two-face corners at each simple vertex gives

    4*f0 >= B >= 5*f2-36*f0.

For the dual simplicial h-vector, f0=sum h_i and
f2=sum binom(i,2)*h_i. Pairing i with9-i in
5*f2-40*f0 gives coefficients100,60,30,10,0 for i=0,...,4.
Thus its value is100*h0+60*h1+30*h2+10*h3>=100, using h0=1 and
nonnegativity. This is the claimed contradiction. A hypothetical
higher-dimensional example has a9-face retaining the matching condition,
so it is excluded as well.

## Scope and dependencies

The bound proved is dimension<=8 for the specifically defined
all-metric-facets-Coxeter class. It does not lower the arbitrary compact
ceiling29, assert that this class has an8-dimensional example, or show
that the known H8 example belongs to it. The label4 sharpness example is
separate from any sharpness claim for the dimension bound.

External inputs are the classification of finite rank-three Coxeter
diagrams, compact Coxeter simplicity and local face geometry, elementary
cyclotomic degree, and Dehn--Sommerville symmetry/nonnegativity for the
h-vector of a convex simplicial polytope. The arithmetic example uses
the usual admissible-lattice interpretation. All new scalar identities
and their geometric applications were checked on paper in this review.
