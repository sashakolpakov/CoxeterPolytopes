# Root review: the eleven-vertex carrier cannot attain the +3 bound

19 October2026. **PASS after complete independent written review.**
Reviewed source `right_angled_ideal4_eleven_vertex_equality_exclusion_20261019.md`,
SHA256 `b9af9b497c69da69a25faced7067b854f0e3b8332f89da57daa12a5a968c2121`.
No scientific script, numerical classification, Lean process or build
was run. The earlier actual-wall/outside-capacity, spinor and carrier
theorems remain explicit reviewed inputs. The new proof does not require
a classification of eleven-vertex three-polyhedra or of arbitrary
ideal right-angled four-polytopes.

## Saturation is used only where justified

Assuming V(Q)=36,V(R)=11 forces q=0,N0=3,G=6 by the earlier integer
defect argument. The small-corner inequalities then force exactly
eight triangles, five quadrilaterals, and one small coefficient at
each nongood R vertex. I checked the degree-four Euler/defect identity
and the equality in each inequality; no uncounted larger polygon
remains.

The six GOOD opposite facets h_u have at least one N0 witness each,
while the total good-facet/N0 incidence count is at most2N0=6.
Thus each good h_u has exactly one vertex outside u and its four
neighbors, and hence exactly six vertices. This conclusion is not
applied to a nongood facet. The six-vertex degree-four graph is the
octahedral graph, and its eight faces are triangles.

I checked the elementary square-cusp proof: send the opposite vertices
to infinity and0. Orthogonality to a rectangle's vertical side fixes
the adjacent hemisphere's center coordinate. Passing through its two
corners and0 yields x_+^2=x_-^2=-y_-y_+ and the horizontal counterparts,
so the centered rectangle is square. Restriction preserves the actual
unit normals. The two R pair-sum coefficients therefore agree at
every good vertex, with the SAME fixed null scale.

## The quadrilateral metric and nongood orthogonality obstruction

Since q=0, all consecutive normalized null products on any R polygon
are-2. For a quadrilateral the displayed four spinors have successive
determinants1 and corner coefficients2/t,t,2/t,t. I independently
checked both the recurrence and the opposite product8/B^2 at the
corners of coefficient B.

If h_u and h_v are orthogonal, then -(t_u,t_v)=1, so u and v are
not adjacent in R. At a nongood u, write its pair coefficients
A<2<=B. If both A-pair facet walls missed v, their actual-wall
products with t_v would sum to at most-2, forcing A>=2. Thus v
belongs to one of these facets. Its polygon coefficient at u is
the OTHER pair coefficient B. It cannot be a triangle; it is a
quadrilateral with u,v opposite. The preceding metric gives B=2sqrt2
and coefficient1/sqrt2 at its two other vertices.

For clarity, the source's needed convex-face fact is that two
nonadjacent vertices of a three-polyhedron cannot share two facets:
their facet intersection would be an edge containing both. This
establishes the asserted opposite position in the quadrilateral.
The separate three-facet/opposite-pair fact below is also valid;
it need not be treated as an equivalent formulation of that stronger
two-facet assertion.

At either of the two1/sqrt2 vertices, goodness gives equal pair
coefficients, while nongoodness gives the other coefficient at least2.
Neither coefficient is1. Every triangle coefficient is1, so both
vertices have four quadrilateral facets. With only five such facets,
their incident sets share at least three. Any three of the four
facets at an ideal rectangular cusp contain an opposite pair whose
hyperplanes share only that ideal point. The two distinct vertices
cannot share those three facets. This excludes the nongood case;
the earlier good-cusp lemma excludes the good case. All eleven h_u
are consequently pairwise nonorthogonal.

## The final witness count and propagation

Only after that global pairwise exclusion is established does the
proof apply the old tangency/edge argument to ALL eleven h_u. Each
has at least six vertices and hence a vertex beyond u and its four
neighbors. If that vertex had an R-neighbor, nonorthogonality would
make the two incident h walls opposite in its box link; the hemisphere
tangency then forces it adjacent to u as well. Feasibility and the
three independent common vertical walls make this an actual edge,
contradicting the chosen nonneighbor.

Each h_u therefore supplies an N0 incidence. At an N0 cusp at most
two of the pairwise nonorthogonal h_u can occur. This gives11<=6,
the contradiction. Thus V(Q)>=37 whenever V(R)=11. No all-facets
octahedral or garland conclusion is used.

Every ideal right-angled five-polytope has a positive carrier with
a three-facet of at least11 vertices. If that three-facet has more
than11 vertices, the old +3 inequality is already stronger. Therefore
its carrier always has at least37 vertices and25 facets. The existing
injection and propagation give respectively F5>=63 and3V5>=631,
so V5>=211. Injection again gives F6>=275. With2f4>=63F6, the
unchanged transport gives272V6>=6039F6 and hence V6>=6106.
I checked6039*275=1660725 and the neighboring272 multiples.

These are independently reviewed necessary size bounds. The ideal
dimension ceiling remains6 and the mixed ceiling remains12. The
new equality geometry,63/211 and275/6106 constants are not in the
current PDF or conditional Lean modules; those still record earlier
stages. No numerical execution is claimed for these substitutions.
