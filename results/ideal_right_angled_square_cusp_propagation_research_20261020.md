# Square-cusp propagation and the cubic-cusp obstruction in dimension five

20 October2026. Paper argument only; no checker, census, optimizer,
Lean or build ran. This is a metric rigidity lemma and a restricted
obstruction. It does NOT exclude arbitrary entirely ideal5/6-polytopes.
No worldwide novelty assertion or new size-bound manuscript is intended.

## 1. Canonical square-cusp scale in dimension three

Let R be an entirely ideal finite-volume right-angled three-polyhedron.
Use outward UNIT wall normals. At a cusp with a square horospherical
rectangle, its two opposite unit-normal pairs have the same positive
sum t. This is a canonically scaled future-null vector, not a freely
rescaled ideal point. On each polygon incident to that cusp the two
intrinsic side normals therefore sum to t. Its polygon corner
coefficient in the fixed-scale spinor recurrence is exactly1.

Suppose all cusps of R are square. For adjacent vertices u,v let a,b
be the two shared wall normals, with a.b=0. Their other cusp normals
are t_u-a,t_u-b at u and t_v-a,t_v-b at v. Set A=-t_u.t_v>0.
The product of t_u-a and t_v-a is1-A. Acuteness and the actual
right-angled wall-product gap imply A=1 or A>=2. If A=1, the
three distinct actual normals

    t_u-a,  t_v-a,  b

are mutually orthogonal. Acute projection of a strict interior point
would give an actual finite vertex, contrary to entire ideality.
Thus A>=2. Every edge consequently has fixed-scale spinor gap
ell=sqrt(A/2)>=1.

The cited gap/projection fact is the full all-wall lemma, not the
unsupported inference that visible pairwise intersections have a
common face. Distinctness of the displayed normals follows from
their mutual product0 and unit norms.

## 2. A triangle forces its neighbors to be triangles

For cyclic projective spinor lifts v_i at a polygon, with
t_i=T(v_i), the exact recurrence at a square cusp is

    v_(i-1)/ell_(i-1)+v_(i+1)/ell_i=v_i.

Consider an edge whose two endpoints are square and whose gap is1.
Its two consecutive recurrences give

    v_prev/ell_prev+v_next=v_current,
    v_current+v_after/ell_after=v_next.

Hence v_prev/ell_prev=-v_after/ell_after. These two rays are the
same ideal vertex. A convex ideal polygon has distinct vertices, so
it must be a triangle: the previous and following vertex have closed
up. The sign is exactly the projective wrap sign for a triangle;
no extra vertex can be inserted on that ray.

An ideal triangular face whose three cusps are square has gap1 on
each edge. Indeed its three unit side normals have pairwise product-1,
and the canonical t at a vertex is the sum of its two adjacent side
normals. The two endpoint t vectors of an edge therefore have product-2.

Every ideal right-angled three-polyhedron has at least eight triangular
faces, by its elementary degree-four Euler identity. Start with one.
When all cusps are square, the preceding argument forces every face
across its edges to be triangular, then propagates through the connected
facet-adjacency graph of the convex polyhedron. All faces are triangles.
Degree four and Euler now give6 vertices and8 faces: the octahedron.
Its ideal right-angled realization is regular, by the usual
three-dimensional realization uniqueness (or the direct six-vertex
coordinate reconstruction in the existing manuscript).

**Lemma.** An entirely ideal right-angled three-polyhedron all of whose
cusp rectangles are square is the regular ideal octahedron.

This proof uses the existing exact spinor identity, elementary Euler
count and acute all-wall projection. It does not classify all larger
ideal three-polyhedra or presume their metrics are regular.

## 3. A boundary criterion without the all-square assumption

The propagation step only needs its two endpoint cusps to be square.
If a triangular face meets a nontriangular face along uv and both
u,v are square, the triangle's two coefficients1 give gap1 on uv
regardless of its third coefficient. To see this without a wrap-sign
convention, let a,b,c be the unit intrinsic normals of the ideal
triangle, so each distinct pair has product-1. Its fixed null vectors
at the endpoints of the shared side are (a+b)/gamma_u and
(b+c)/gamma_v. Their product is -2/(gamma_u*gamma_v); here both
coefficients are1. The neighboring polygon also has
both endpoint coefficients1, so it must close as a triangle, a
contradiction. Thus every triangle/nontriangle boundary edge contains
a nonsquare cusp.

This identifies the genuine metric barriers to octahedral propagation.
It is not a uniform bound on their density. In particular a growing
polyhedron might have many triangle patches separated by nonsquare
cusps, and no estimate here rules that out.

## 4. What this gives in dimensions four and five

If every cusp box of an entirely ideal right-angled four-polytope Q
is cubic, every three-facet inherits square cusp rectangles. The lemma
forces all three-facets to be octahedra. The already reviewed internal
octahedral-facet recognition theorem then makes Q the regular24-cell.
Its kappa=3V(Q)-4F(Q)+24 is0. This last recognition is an explicit
existing internal input, not a new classification of arbitrary Q.

If every cusp box of an entirely ideal right-angled five-polytope P
were cubic, all four-facets would have cubic cusp boxes, hence all would
be regular24-cells. The universal ideal5 identity

    sum_(four-facets Q) kappa(Q)=48

would then read0=48. Therefore any hypothetical ideal5-polytope
must have a NONCUBIC cusp. More precisely, any positive-kappa carrier
Q must contain a nonoctahedral three-facet R, and that R has a
nonsquare cusp, which is also a noncubic cusp of P.

This excludes the all-cubic hypothesis, not dimension five. The positive
kappa sixteen-cell double already constructed in the repository has
nonsquare/noncubic cusps and is not contradicted. Its repeated doubles
also rule out a universal upper bound on kappa of an individual facet.
An unconditional dimension exclusion still requires a quantitative
compatibility loss among distinct facets, beyond the single-carrier
and square-propagation lemmas.

The source interfaces used here are the fixed-scale polygon recurrence,
acute projection and ideal5 carrier identity in
`manuscripts/ideal_right_angled_six_face_constraints.tex`, and the
reviewed octahedral recognition in
`manuscripts/ideal_ra4_volume_and_structure.tex`. No prior source or
proof packet was changed.
