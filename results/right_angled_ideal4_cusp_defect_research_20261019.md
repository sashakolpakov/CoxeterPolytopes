# A positive cusp defect for ideal right-angled four-polytopes

19 October 2026. New written research proof, awaiting independent review.
No computation, Lean process, optimizer or TeX build was run. Previous
pinned evidence is unchanged. The attempted stronger +6 assertion was
withdrawn after detecting an incorrect octahedron normalization; none of
that argument is used here.

**Proposed theorem.** If Q is an entirely ideal finite-volume right-angled
four-polytope and R is any of its three-facets, with m vertices, then

    V(Q) >= 3m+3.                                    (A)

The argument forces vertices with no graph-neighbor in R, rather than
only substituting a larger facet count. With the independently reviewed
at-least-eleven-vertex carrier theorem and the unchanged ideal-six
transport interface, it would give ideal RA5 at least62 facets and189
vertices, and ideal RA6 at least252 facets and5504 vertices. These
stronger counts are pending review, and no dimension exclusion is claimed.

## 1. Actual opposite walls and exact outgoing-face count

Use outward unit normals and the Lorentz metric. Let n be R's normal.
For each ideal vertex u of R, let h_u be the unique actual ambient facet
opposite n in the rectangular cusp link at u, and put

    t_u=n+h_u.

Then t_u is future null, lies in n-perpendicular Lorentz three-space,
and represents u. Different u give different h_u. Actual right-angled
wall products are0 or at most-1. The all-wall acute projection proof
from the [capacity report](right_angled_outside_cusp_capacity_20261019.md)
is used whenever an actual finite vertex or ridge is asserted below.

Write N_i for the number of vertices outside R with exactly i graph-neighbors
in R. The reviewed capacity theorem gives i in{0,1,2}. Every vertex of R
has four outgoing edges, so

    N_1+2N_2=4m,
    V(Q)=3m+N_0+N_1/2.                              (1)

Let uv be an edge of R, supported by the two other actual facets a,b.
The four walls a,b,h_u,h_v cannot be mutually orthogonal: that would
give an actual finite vertex of Q. Since h_u,h_v are both orthogonal
to a,b, their product is at most-1. Thus

    d_uv:=-(t_u,t_v)=1-(h_u,h_v)>=2.                 (2)

The two-face G_uv cut out by a,b is an ideal polygon with side uv on R.
Its sides at u and v lie respectively on h_u and h_v. It is a triangle
exactly when (h_u,h_v)=-1, equivalently d_uv=2. In that case the third
vertex is the future null vector h_u+h_v, feasible against every actual
wall by the positive-sum argument. The three halfplanes from n,h_u,h_v
already bound the ideal triangle; all other inequalities contain its
three vertices, hence the entire triangle.

Conversely every outside vertex with two R-neighbors gives one such
triangle. In the box proof of capacity, its hemisphere center is on a
box edge and the two contact corners are that edge's endpoints. The
n-2=2 common vertical walls cut out this actual two-face, and n cuts
out the side joining the two contact points. Its three vertices and
three sides are therefore precisely this triangle. The other endpoint
pair determines it uniquely. Hence N_2 is exactly the number of edges
of R whose outgoing two-face is triangular, not just bounded by it.

The ideal right-angled three-polyhedron R has2m edges. Let q count the
edges with nontriangular outgoing two-face. Then

    q=2m-N_2=N_1/2,                                 (3)
    V(Q)=3m+N_0+q.

This equality includes all outgoing edges and all outside vertices.

## 2. Normalized cusp widths in R

The four unit facet normals of R through u have two opposite pairs.
Their sums are

    a_1+a_3=alpha_u t_u,
    a_2+a_4=beta_u t_u,      alpha_u,beta_u>0.        (4)

These coefficients use the fixed normalization t_u=n+h_u, not an
independently rescaled cusp vector. At each corner of a polygonal face
of R, the sum of its two intrinsic side normals equals gamma*t_u, where
gamma is one of alpha_u,beta_u. The same coefficient occurs on the two
opposite polygonal faces through u. Thus the four polygon-corner
coefficients at u are alpha_u,beta_u,alpha_u,beta_u.

Call u good if both alpha_u and beta_u are strictly less than2; let G
count good vertices.

### Polygon lemma

Every polygonal face of R has at least two corner coefficients gamma<2.

Here is a self-contained normalization and proof. In the Lorentz plane
of the polygon use

    T(a,b)=(a²+b²,a²-b²,2ab),
    (T(v),T(w))=-2 det(v,w)².

For its k cyclically ordered vertices choose nonzero spinors v_i with
t_i=T(v_i), choosing signs so all successive determinants are positive.
Use the endpoint conventions

    v_0=-v_k,   v_{k+1}=-v_1,
    ell_i=det(v_i,v_{i+1})>0,   ell_0=ell_k.

Equation(2) gives ell_i>=1. If l_{i-1},l_i are the outward unit side
normals at vertex i, direct two-dimensional algebra gives

    l_{i-1}+l_i=gamma_i t_i,
    gamma_i=det(v_{i-1},v_{i+1})/(ell_{i-1}ell_i),
    v_{i-1}/ell_{i-1}+v_{i+1}/ell_i=gamma_i v_i.      (5)

For example, send v_i to(1,0), and write
v_{i-1}=(a,-ell_{i-1}), v_{i+1}=(b,ell_i).
The two side normals are
(a/ell_{i-1},a/ell_{i-1},-1) and(b/ell_i,b/ell_i,1);
their sum is the asserted coefficient times T(1,0). The vector
recurrence follows from the two-dimensional determinant identity.
These formulas are invariant under determinant-one spinor changes.

The k-by-k symmetric matrix of the recurrence has diagonal gamma_i,
successive off-diagonal entries -1/ell_i, and wrap entries +1/ell_k.
Both coordinate arrays of the v_i lie in its kernel; they are independent,
because the ideal vertices are distinct. Its nullity is at least two.

If at most one gamma_i were below2, delete that vertex, or any vertex
if none is below2. The remaining principal matrix is a path matrix with
diagonal at least2 and off-diagonal magnitudes at most1. It is positive
definite: its quadratic form is at least

    x_1²+x_{k-1}²+sum_j(|x_j|-|x_{j+1}|)²>0

for every nonzero vector. Consequently the full matrix would have rank
at least k-1 and nullity at most one, a contradiction. This proves the
polygon lemma, including the triangle case and the wrap signs.

A triangular face all three of whose outgoing edge-faces are triangles
has ell_1=ell_2=ell_3=1 and gamma_1=gamma_2=gamma_3=1.
This follows either from(5), or by taking the products of the three
side-normal sums: gamma_i gamma_j=1 for each pair. Such a triangle
contributes three small corners instead of the lemma's guaranteed two.

## 3. At least six good-corner units after outgoing defects

R has m+2 polygonal faces and at least eight triangles. Indeed its degree
is4, so Euler gives F=m+2 and the polygon identity is

    number of triangles - sum_{k>=5}(k-4)*number of k-gons = 8.

Each of the q exceptional edges belongs to two polygonal faces. Thus
at least max(0,8-2q) triangular faces have no exceptional edge. Let S
count polygon corners with gamma<2. Section2 gives

    S>=2(m+2)+max(0,8-2q).                            (6)

At a good vertex all four corners are small; at any other vertex at most
two are small, by the paired coefficients(4). Hence S<=2m+2G, giving

    G>=2+max(0,4-q).                                 (7)

No bounded polygon census is used in this estimate.

## 4. Good cusps force outside vertices with no R-neighbor

First a good h_u is orthogonal to no other opposite wall h_v. If it
were, (t_u,t_v)=-1. By(2), u and v cannot be adjacent in R. Nor can
distinct nonadjacent vertices of a convex three-polyhedron lie in two
common facets: their intersection would be an edge containing both.
Therefore v is incident to at most one of the four facets of R through u.
One opposite pair in(4) consists of two facets neither incident to v.
Their products with t_v are each at most-1, using actual ambient wall
products with h_v. Taking the product of their sum with t_v then forces
alpha_u>=2 or beta_u>=2. This contradicts goodness.

The three-facet h_u is itself an ideal right-angled three-polyhedron.
At u its four edges are exactly the four outgoing edges from R, so it
contains u and their four distinct other endpoints. It has at least six
vertices: its degree4, Euler and faces of at least three sides give
2E=4V, F=V+2 and3F<=2E, hence V>=6. Thus it has an additional vertex w
not adjacent to u. Adjacency inside h_u and ambient adjacency agree for
two of its vertices, because an ambient edge with both endpoints on a
face lies entirely on that face.

We claim this w has no R-neighbor. Otherwise take such a neighbor v.
Then h_v contains the edge wv, so h_u,h_v are both incident to w.
At an ideal rectangular-box cusp their product is either0 or-1.
It cannot be0 by goodness. If it is-1, they are an opposite pair.
Send w to infinity. The hemisphere of n is tangent to these two
parallel walls, at the old points u and v respectively. Its center is
midway between them and its radius is half their separation. Since v
is a box corner (it is joined to w), the opposite tangency point u has
the same other coordinates and is also a box corner. It is then joined
to w as well, contradicting the choice of w. Therefore w lies in N_0.

Each good h_u has at least one such N_0 vertex. At any N_0 vertex no
three good h_u can meet: among any three distinct faces of its
rectangular three-box link there is an orthogonal pair. That would
contradict goodness. Consequently

    2N_0>=G.                                        (8)

This argument uses actual opposite facets, not independently chosen
normalized normals or a selected subset of outgoing edges.

## 5. The positive defect and propagated counts

Combine(3),(7),(8):

    V(Q)-3m=N_0+q
        >=q+ceil((2+max(0,4-q))/2)>=3.               (9)

For q=0,1,2,3, the displayed lower bounds are respectively3,4,4,5;
for q>=4 they are at least q+1>=5. This proves(A), conditional on
independent review of the written geometric argument.

The new large-facet carrier theorem gives a positive-kappa Q inside
every ideal RA5 polytope, with a three-facet R of m>=11 vertices.
Separately, the actual cusp-opposite/ridge-neighbor injection gives
F(Q)>=2m+3. Hence

    V(Q)>=3m+3>=36,       F(Q)>=2m+3>=25,
    F(P5)>=1+25+36=62.                              (10)

The unchanged strict ideal-five Euler argument yields

    V(P5)>=3F(P5)+3>=189.

The next-dimensional injection gives

    F(P6)>=1+62+189=252.

Retain the reviewed sharp reciprocal/actual-incidence interface
a_4-F<=136V/99. Since each ideal five-facet now has at least62 facets,

    2a_4>=62F,
    V>= (1485/68)F,
    V>=ceil(374220/68)=5504.                         (11)

Here68*5503=374204<374220 and68*5504=374272. These arithmetic
consequences do not replace the unchanged universal transport proof.
The same published NK/Dufour substitutions would improve the strict
ideal-five fraction in dimensions11/12 to2/23 and the ideal-six fraction
to1/15.

## Scope and remaining gap

The universal new assertion is(A), with the stronger structural interfaces
(3),(7),(8). The all-wall acute normal/projection fact, entirely ideal box
links, face inheritance, and the earlier capacity-two theorem are explicit
dependencies. The spinor lemma and incidence counts are proved above;
they are not sampled. The numerical applications additionally inherit the
reviewed Figure7/small-facet carrier theorem, the published ideal4 facet
minimum24 for strict ideal5 Euler, and the old136/99 transport argument.

Independent review is pending. No stronger +6 theorem is claimed here,
no new Lean coverage is claimed, and no existence or ideal-dimension6
exclusion follows. The unchanged dimension ceilings remain6 entirely
ideal and12 mixed right-angled.

