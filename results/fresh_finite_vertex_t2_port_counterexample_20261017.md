# Finite-vertex failure of two-ideal triangle-port transitivity

The proposed extension of maximal-simplex port blocks to triangles having
exactly two ideal vertices is false, even at a finite vertex of a finite-volume
hyperbolic Coxeter polyhedron. The example below has two such triangles sharing
a finite–ideal edge; their third port pair spans a pentagon. This blocks that
particular proof route to either proposed inequality
\(T_2\le(5/2)E_{fi}\) or \(T_2\le2E_{fi}\).
**It does not disprove either numerical inequality.**

The older [multi-ideal triangle example](fresh_multi_ideal_triangle_scope_20261017.md)
already shows failure at a simple ideal vertex along an ideal–ideal edge.
The present example checks the additional finite-vertex premise explicitly;
no novelty claim is made. The
[standalone exact checker](../scripts/fresh_finite_vertex_t2_port_counterexample_20261017.py)
and [saved certificate](fresh_finite_vertex_t2_port_counterexample_20261017.json)
use rational and radical SymPy arithmetic in \(\mathbb Q(\sqrt2,\sqrt3)\), with
no imported project arithmetic and no numerical sign acceptance.

## Literal chamber

In the upper half-space model \((X,Y,Z)\), \(Z>0\), put
\(c=\sqrt3/2\). Take the region

\[
 -\tfrac12\le X\le c,\qquad -\tfrac12\le Y\le c,
\]
\[
 X^2+Y^2+Z^2\ge1,\qquad
 (X-c)^2+(Y-c)^2+Z^2\ge\tfrac12.
\]

With Lorentz form \(\operatorname{diag}(-1,1,1,1)\), its outward normals of
square two are

\[
\begin{array}{c|rrrr}
 B&0&0&0&-\sqrt2\\
 L&1/\sqrt2&-\sqrt2&0&1/\sqrt2\\
 D&1/\sqrt2&0&-\sqrt2&1/\sqrt2\\
 R&\sqrt{3/2}&\sqrt2&0&\sqrt{3/2}\\
 U&\sqrt{3/2}&0&\sqrt2&\sqrt{3/2}\\
 C&2&\sqrt3&\sqrt3&0
\end{array}
\]

The inequalities are \(\langle x,N_i\rangle\le0\). The nonright finite
contacts are \(B\cdot L=B\cdot D=-1\) and
\(B\cdot R=B\cdot U=-\sqrt3\), giving labels three and six.
The opposite vertical pairs have product \(-2\). The two remaining dotted
products are \(C\cdot L=C\cdot D=-\sqrt2-\sqrt6<-2\).
All other pairs are right. Thus every pair satisfies a native Coxeter gate.
The future unit timelike point \((5/4,0,0,3/4)\) is strictly inside all six
halfspaces.

## Complete geometry and finite volume

The following are all seven vertices. The coordinates in the table are
upper-half-space coordinates; the certificate stores exact projective
Lorentz rays.

| Vertex | Coordinates | Incident facets | Type |
|---|---|---|---|
| \(A\) | \((-1/2,-1/2,1/\sqrt2)\) | \(B,L,D\) | finite |
| \(I_1\) | \((-1/2,\sqrt3/2,0)\) | \(B,L,U\) | ideal |
| \(I_2\) | \((\sqrt3/2,-1/2,0)\) | \(B,R,D\) | ideal |
| \(I_\infty\) | \(\infty\) | \(L,D,R,U\) | ideal |
| \(F\) | \((\sqrt3/2,\sqrt3/2,1/\sqrt2)\) | \(C,R,U\) | finite |
| \(G\) | \((\sqrt3/2,1/(2\sqrt3),1/\sqrt6)\) | \(B,C,R\) | finite |
| \(H\) | \((1/(2\sqrt3),\sqrt3/2,1/\sqrt6)\) | \(B,C,U\) | finite |

Completeness is checked in a bounded projective section. Let
\(D_0=\sum_iN_i\); explicitly

\[
 D_0=(2+\sqrt2+\sqrt6,\sqrt3,\sqrt3,\sqrt6),\qquad D_0^2<0.
\]

On the nonzero feasible cone, \(-\langle D_0,x\rangle>0\): equality would
force every wall product to vanish, while the six normals span the four
dimensional Lorentz space. Consequently the section
\(-\langle D_0,x\rangle=1\) is bounded. Every vertex has three independent
active wall equations. Exhausting all twenty triples, and deduplicating their
rays, gives exactly the seven vertices above, with no spacelike ray. The
nonsimple ideal vertex is explicitly retained rather than treated as simple.
Positive sums of each facet's vertex rays give exact strict relative-interior
facet witnesses.

The two mixed-coordinate cusps have positive affine relations of type
\((2,3,6)\), hence bounded Euclidean triangle links. The cusp at infinity
has the rectangular link \(\widetilde A_1\times\widetilde A_1\).
The checker verifies the full incident Gram matrices are positive
semidefinite of rank two and their positive null relations. These bounded
links, together with the complete finite/ideal vertex census, establish
finite volume. This is an actual finite-volume chamber, not merely a
Lorentzian local star or a hypothetical completion.

## Why the block premise fails

The complete faces are

\[
\begin{array}{c|l}
 B&A,I_1,H,G,I_2\quad\text{(pentagon)}\\
 L&A,I_1,I_\infty\quad\text{(triangle)}\\
 D&A,I_2,I_\infty\quad\text{(triangle)}\\
 R&I_2,G,F,I_\infty\quad\text{(quadrilateral)}\\
 U&I_1,H,F,I_\infty\quad\text{(quadrilateral)}\\
 C&F,G,H\quad\text{(triangle)}.
\end{array}
\]

The finite vertex \(A\) has the three edge neighbors
\(I_1,I_2,I_\infty\). Two of its port pairs are triangular, on \(L\) and
\(D\), and both triangles have exactly two ideal vertices. They share the
edge \(A I_\infty\). The remaining port pair lies on the pentagon \(B\).
Thus triangular adjacency is not transitive at this finite vertex.
Neither triangle extends to a larger simplex face: the only three dimensional
face is the whole chamber, which has seven vertices and six facets.
The two maximal simplex faces therefore overlap on a finite–ideal edge.

The compact/simple-polytope proof cannot be repeated here. At the shared
ideal endpoint \(I_\infty\), the degree is four, so the closed degree-three
subgraph argument fails. The frozen finite-ended theorem for triangles with
at least two finite vertices is unaffected.

For scale, this example has eleven edges, \(T_2=2\), and \(E_{fi}=6\).
Both numerical inequalities under discussion hold in this example. A proof
of either one needs an incidence argument that allows overlapping triangles
at finite–ideal edges; the six-simplex ideal-count improvement alone does
not repair the missing block decomposition.

**FINAL/FROZEN.** The exact checker exhausts twenty triple frames,
all fifteen normal pairs, all six actual facets, all eleven edges, and all
three cusp links. Generation and a fresh saved replay both passed, together
with source-hash, syntax, and report-link checks. The parent independently
read the explicit geometry and the scope distinction; no second arithmetic
implementation is claimed. No higher-dimensional bound or new example in a
record dimension is claimed.
