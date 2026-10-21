# Root review: compact right-angled cutting and farthest facets

20 October 2026. Independent written mathematical review. No numerical
job, classification search, Lean compilation or TeX build was run.

## Sources and scope

Reviewed in full:

* `compact_ra4_dodecahedral_double_spectrum_cutting_20261020.md`, SHA256
  `08129e5aa206016d758114fda97df482c78de5397f7a3631bc542f79b0cbb34b`;
* `compact_ra4_farthest_dodecahedron_wall_radius_20261020.md`, SHA256
  `d082e4a918d4504a59ac200f32bbb141f1f2cb52d701006eaea7ba5b0e5bcc22`.

Both written arguments pass, with the inherited geometric inputs and
the stated intrinsic metric hypotheses retained. Neither proves an
unrestricted minimum-volume or second-volume theorem. Neither is a new
dimension bound. They remain research notes, rather than additions to
the compiled manuscripts or verified Lean interfaces.

## Intrinsic D/D2 cutting

The distinction between a metric double and an arbitrary pentagon/hexagon
face census is essential. In the metric double D2, each of the five
hexagons has lengths (2L,L,L,2L,L,L). Its long-side midpoint segment is
therefore canonical. At a seam vertex, the point is in the interior of
an actual parent edge. Simplicity gives exactly three parent facets;
all three are D2 because each neighboring ridge at that long edge is
hexagonal. This verifies the local closure without assuming a guessed
four-dimensional face graph.

The four-normal span step supplies actual geometry. For adjacent
component facets F,G, the two other facets H,I at the seam endpoints
have ultraparallel normals: they correspond to the two nonadjacent
sides neighboring a fixed side of a compact right-angled pentagon.
The Gram matrix on F,G,H,I is an orthogonal sum of I2 and a Lorentz
two-plane. Its nondegenerate Lorentz four-space is exactly the space
attached to each of F and G. Thus the whole component lies in one
actual hyperbolic three-hyperplane.

The section argument covers all boundary strata. In a seam pentagon's
interior there is exactly one active parent wall; on a seam edge there
are exactly two; at a seam vertex there are exactly the three parent
walls of the actual parent edge. No outside wall passes through such a
point. The closed seam surface is consequently also open in the
connected boundary of the compact convex section. It is the whole
boundary, rather than an abstract surface inside a possibly larger
section. Its regular pentagons and right angles identify it with D by
metric continuation.

An outside facet cannot meet the section's interior either: a linear
support functional nonpositive on a full-dimensional section cannot
vanish at an interior point without vanishing on its whole plane.
The proposed cut therefore meets exactly the twelve component D2
facets, each in its canonical pentagon, and is orthogonal to them.
Both halves have interior, remain compact and right-angled, and retain
the D/D2 metric spectrum. Other components stay on one strict side;
their sections are disjoint from this one. Cutting removes twelve D2
facets from the total and adds one piece.

The final all-D recognition uses whole metric facets and their unique
orthogonal continuations. Facet counts alone would not suffice. It
follows that the pieces are 120-cells, N_D2=12(m-1), and volume is
m times 34*pi^2/3. For m=2 there is one common whole D facet, giving
the reflected double uniquely. No path-shaped adjacency graph is
deduced. The supplementary counts f=106m+14, v=560m+40 and W=30(m-1)
agree with v-5f=W and the stated volume formula.

## Farthest facets and wall radius

The projected point is in the actual facet: all other parent wall
inequalities remain strict by the all-wall nonpositive-product
property. Only the neighbors of the selected facet supply its
intrinsic normals. Their heights are strictly less than one, so the
radial unit directions have pair products strictly less than 1/2.

The quartic expansion was independently expanded on paper:

    (x+1)(x+2/5)^2(x-1/2)
      = x^4+(13/10)x^3+(3/50)x^2-(8/25)x-2/25.

Its printed Legendre coefficients reproduce this polynomial. The
rank-two, rank-three and rank-four trace-free tensors have the stated
contractions. In particular the rank-four trace is
u^2-(7u^2+I)/7+5I/35=0. Positivity yields N<=14. At N=14 all off-diagonal
products must be -1 or -2/5, because strictness excludes 1/2; their
vector-sum square is negative. Thus N<=13.

The inherited polygon-capacity theorem excludes a hexagon when there
are at most thirteen faces. Euler then forces twelve pentagons. The
source explicitly imports twelve-face metric recognition, as does
the current compact paper; it is not supplied by the quartic alone.
This makes every maximizing facet a metric D.

Opposite pairs of D wall normals give max intrinsic height >= a*h,
where a>=1 and h^2=phi/2. Combining this with the actual projection
bound gives T(q)^2>=phi^3. Equality forces the foot to be the center
and every neighboring parent wall to have the same maximal height.
Repeating at those neighbors propagates the all-D conclusion through
the connected facet graph. The normal coordinate then fixes q as the
120-cell center. The converse follows from symmetry and the same
projection identity.

This controls a maximum distance to supporting walls. It supplies no
bound on the minimum distance, no inscribed 120-cell, and no volume
comparison. The outgoing-edge inequality coth(s_v)coth(s_w)<=phi
is correctly kept as a product constraint, rather than a separate
lower bound on every edge.

## Exact collar reduction: separate subsequent source

Also reviewed in full after the two arguments above:
`compact_ra4_exact_dodecahedral_collar_cover_20261020.md`, SHA256
`c1a3f94c880d37a6180a7b83840eee05c27883eeb8161d77f63ddaa26b942669`.
PASS as an exact reformulation, without a volume-minimum conclusion.

Normal projection into the whole actual F follows because its neighboring
wall functionals multiply by cosh(t). Every other compact wall has the
unique displayed timelike projection, with its future orientation forced
by strict negativity on the interior of F. All-wall acuteness places its
center v in F. Its roof inequality is exactly
tanh(t)<=tanh(L)cosh d(x,v). Compactness gives a finite open-ball cover
of all of F and a uniform roof bound strictly below one. The integral
primitive of cosh^3 is sinh+(sinh^3)/3 as printed.

The pairwise wall gap is stronger than merely nonpositive products.
Product zero gives cosh(d)=cosh(rho_i)cosh(rho_j), while product<-1
gives cosh(d)>cosh(rho_i+rho_j). These are respectively orthogonal
sphere boundaries and strictly disjoint balls. The boundary-plane
condition follows from the same gap. In the converse, each cap uniquely
wins at its own center: its own profile is sech(rho_i)<1, and every
other cap's profile there is at least cosh(rho_i)>1. Continuity gives
an open region of the interior of F even when the center is on its
boundary. Hence each proposed cap is an actual facet. The base and
neighboring walls remain actual at small positive heights. This checks
the claimed construction correspondence, rather than just one direction.

One useful consequence of the exact conditions is that every vertex
of F has exactly one cap center at that vertex. A ball covering a vertex
must have its center on each incident face plane: otherwise the strict
boundary-distance condition prevents it reaching that plane. Their
intersection is the vertex. Two centers at the same point are ruled
out by pair compatibility. This is still a local fact, and gives no
minimum cost for the complete cover. The global integral inequality
needed for a universal120-cell volume minimum remains unproved.
