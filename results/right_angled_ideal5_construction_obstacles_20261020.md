# Ideal-five construction log: two explicit local failures

20 October2026. Source-only derivation; no scientific job ran. This is a
short construction record, not a manuscript or a nonexistence theorem.
Neither an entirely ideal right-angled H5 example nor a mixed H9 example
has been obtained. Earlier sources remain unchanged.

## 1. A finite-volume32-cell block, and the restriction of cap-only gluing

Use the form -t^2+sum(x_i^2) in dimension6. The known16-wall P5 has
positive-orthant Klein body

    y_i>=0, y_i+y_j<=1, sum(y_i)<=2.

Reflect it across its five mutually orthogonal coordinate walls. The
32 copies have disjoint interiors and form the convex polytope Q with
exactly the following72 unit outward roots:

    (1, epsilon_i e_i+epsilon_j e_j), i<j: 40 roots;
    (2, epsilon_1,...,epsilon_5):          32 roots.

Here each epsilon is independently +/-1. Equivalently Q is

    |y_i|+|y_j|<=1, sum(|y_i|)<=2.

All distinct root products are0 or at most-1. This is a genuine
finite-volume right-angled H5 polytope: it is exactly the stated
32-cell reflected union of P5. Its vertices, rather than only its
facet-pair gate, are completely explicit:

|Vertices|Number|Euclidean square in the Klein ball|
|---|---|---|
|signed coordinate axes|10|1|
|four signed half coordinates, one zero|80|1|
|one signed2/3 coordinate, four signed1/3 coordinates|160|8/9|

To see completeness, intersect with each coordinate orthant and use
the five elementary vertex classes of the P5 inequalities. The origin
is now interior. A three-half-coordinate point has only three active
pair walls and is not a vertex of Q. The three remaining classes give
the table. At a half-four point the six pair walls have rank4 and the
two cap walls differing at its zero coordinate supply the fifth rank.
All160 finite vertices are therefore retained, not mistaken for cusps.

Call the32 roots of time coordinate2 the caps. A finite vertex belongs
to exactly ONE cap: its sign vector, with absolute-coordinate sum2.
That cap has exactly five finite vertices. Any two distinct caps have
product-1,-3,-5,-7 or-9, so their facets are disjoint at finite points.

Consequently the restricted strategy of gluing Q blocks ONLY across
these cap facets cannot produce an ideal polytope in finitely many
steps. The cap reflections generate the free product of32 order-two
groups: their halfspaces have disjoint exteriors, or equivalently the
reflection Coxeter diagram has no finite bonds. Its chamber adjacency
graph is a tree. A connected finite union of k blocks has k-1 glued
cap pairs and therefore

    32k-2(k-1)=30k+2

exposed caps. Every one retains its five finite vertices. No other
internal cap contains one of those vertices, and no pair-wall facet
is used for gluing. Thus every convex assembly in this restricted
construction has exactly150k+10 finite vertices. This is NOT a
restriction on gluing other facets or on arbitrary right-angled H5
polytopes. In dimension4 the analogous coordinate-reflection block
of P4 is the ideal24-cell; the additional cap family is the precise
new obstruction in this particular five-dimensional construction.

## 2. A74-wall24-cell lift has192 actual finite vertices

Let V be the24 unit Euclidean vertices of the regular24-cell in R4:
the eight signed coordinate axes and the sixteen signed half vectors.
Let R be its24 facet normals +/-e_i+/-e_j. In coordinates(t,x,z),
consider all74 unit roots

    (1,r,0),                         r in R;
    (sqrt2,sqrt2 v,+/-1),            v in V;
    (1,0,+/-sqrt2).

Every pair product is0 or at most-1. The intended old24-cell ideal
points and the proposed upper/lower cap ideal points therefore do
not expose the defect by a pair-product test.

For each of the96 old24-cell edges vw, v.w=1/2. Set

    X=(1,(v+w)/2,+/-1/(2sqrt2)).

Then X.X=-1/8. The three old roots common to the edge and the two
new roots associated to v,w with the chosen height sign are five
mutually orthogonal unit roots, all perpendicular to X.

All other wall products are strictly negative. For an old wall use
the fact that the midpoint lies in the relative interior of the
edge vw. For a new root, a24-cell vertex u different from v,w has
u.v,u.w<=1/2 and hence u.(v+w)/2<=1/2<3/4; only u=v,w reaches3/4.
The opposite-height roots and both caps are strict as well. These
are therefore192 distinct actual simple finite vertices, not formal
corner diagrams. This lift cannot be an entirely ideal example.
No assertion about its full vertex list or total volume is needed
for that conclusion.

## 3. A concrete filter for the next asymmetric candidate

Put an ideal cusp of a prospective ideal H5 polytope at infinity.
Its vertical facets bound a Euclidean4-box B. The other walls are
hemispheres with centers c_i and radii r_i. Right angles require
center distance squared r_i^2+r_j^2; disjoint/tangent pairs require
distance at least r_i+r_j. Against a box plane, a sphere center must
lie on it or have distance at least its radius.

At every vertex of the weighted nearest-site cells in B, check

    min_i (|x-c_i|^2-r_i^2).

A positive value gives an open gap and infinite volume. A negative
value at an actual full-rank cell/box vertex gives a finite vertex.
An entirely ideal candidate must avoid both errors. Pair angles,
selected cusp stars or one compact section alone do not certify it.
The symmetric constructions above fail the second test. An asymmetric
boundary/interior sphere pattern, or a different complete root
configuration, remains the actual construction obligation.
