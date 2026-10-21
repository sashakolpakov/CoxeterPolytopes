# A uniform dodecahedral facet forces whole blocks in its neighbors

20 October 2026. Paper proof only; no jobs, search, manuscript or Lean
change. This is a conditional structural lemma for the compact volume
problem, not a universal120-cell minimum or containment assertion.

Write `L=arcosh(phi)`, phi=(1+sqrt(5))/2. Use actual outward unit
normals, with negative Lorentz time. The compact right-angled inputs
are simplicity, nonpositive products of distinct actual wall normals,
and the reviewed three-consecutive-L-edge polygon closure lemma from
`manuscripts/compact_ra4_volume_constraints.tex`, SHA256
`afe807672f3a6bbf4bb2dd95d2861b60391d7c2d7c879c6d6382e749e7ad976f`.
The last lemma identifies the ENTIRE polygon as a regular pentagon.
The three-dimensional equality input is that the compact right-angled
polyhedron with20 vertices is the regular dodecahedron D; the same
source explicitly records its combinatorial and metric recognition.

## 1. A regular pentagonal flower contains a whole dodecahedron

**Lemma.** Let Q be a compact right-angled H3 polyhedron. Suppose one
actual face F and all five actual faces A_i adjacent to it are ENTIRE
regular right-angled pentagons. Then Q contains the canonical regular
dodecahedron on this flower. Either Q is D, or

    V(Q)>=30,             (V(Q)-20)/2>=5.               (1)

Equality in(1) identifies Q with the double of D across the canonical
opposite pentagon. No assertion is made that larger Q is a garland.

Match the actual F with a face of the regular D. Its five adjacent
support planes are fixed uniquely: each is perpendicular to F's plane
along its given edge. An entire regular pentagon A_i in that plane
is uniquely fixed by this shared edge and the inward side of F.
Therefore the whole six-face flower agrees metrically with D's.

The outer edges of the five A_i force five further actual supporting
planes B_i. They are precisely the five remaining canonical D planes
other than its opposite plane U. Here actualness is essential: these
are true edges of entire A_i faces, so each has its actual neighboring
face in Q. In a simple right-angled three-polyhedron its normal is
uniquely fixed by orthogonality and the given edge. Thus eleven of
the twelve canonical D halfspaces are actual Q halfspaces. The flower
already supplies fifteen of D's twenty vertices.

Let C be the opposite regular pentagon in U, let c be U's outward
unit normal for D, and denote by b_i the five actual B_i normals.
They are orthogonal to c; within the Lorentz three-space U, their
halfspaces bound exactly the compact C.

We need the ORIENTATION of the homogeneous cone, not just an
unsupported future choice. At C's center o, write

    b_i=h*o+sqrt(1+h^2)*r_i,  h>0,

where the five Euclidean unit r_i form a regular pentagon. Their sum
is zero and they positively span o-perpendicular. If

    K={u in U: <u,b_i><=0 for all i},

then pairing with `sum b_i=5h*o` shows that every u in K has
nonnegative time `a=-<u,o>`. At a=0 the positively spanning r_i force
u=0. At a=1 the inequalities have a bounded Euclidean pentagonal
section. Its whole section lies strictly inside the Klein disk:
otherwise a segment from its interior center would approach the unit
circle within these halfspaces, contradicting compactness of C.
Consequently every nonzero u in K is FUTURE timelike. This proves
the full cone statement, including exclusion of past and spacelike
components.

Now let m be any actual Q wall other than the five B_i. All-wall
acuteness gives `<m,b_i><=0`, so its projection m_U belongs to K.
If nonzero, m_U is future timelike and hence `<m,x><0` at EVERY
future unit point x of U, in particular throughout C. If m_U=0,
the unit m is either c or-c. The known flower lies strictly on D's
side of U, so the actual orientation can only be c. Therefore all
points of C are feasible for Q. Together with the fifteen feasible
flower vertices they give all twenty vertices of the canonical D.
Their convex hull proves D is contained in Q; a few matched vertices
alone were not used to assume containment.

If c itself is actual, the eleven other canonical inequalities give
Q contained in D as well, so Q=D. Otherwise all non-B_i actual walls
are strict on C, and the center of C is an interior point of Q.
The plane U cuts Q into full-dimensional compact pieces, with the
near piece exactly D and the far piece Q_far. This cut meets only
the five B_i along C, and does so orthogonally. Thus Q_far is itself
a compact right-angled H3 polyhedron. No other actual wall can
truncate C or insert an extra seam vertex.

Euler and trivalence, with every face having at least five sides,
give `V(Q_far)>=20`. The five cut vertices occur once in each piece
and were not original Q vertices. Therefore

    V(Q)=20+V(Q_far)-10>=30.

At equality Q_far has20 vertices. Euler forces twelve pentagonal
faces, and the stated metric recognition gives D. Matching its
entire cut pentagon on the opposite side of U makes it the reflected
canonical copy. This proves the lemma and its equality statement.

## 2. Application to an actual uniform facet in H4

**Corollary.** Let P be compact and right-angled in H4, with an actual
dodecahedral facet F whose twenty outgoing edge lengths all equal L.
For every actual neighboring facet G, either G is D, or

    V(G)>=30,             d(G)=(V(G)-20)/2>=5.           (2)

If V(G)=30, G is intrinsically the double of D. This conclusion is
restricted to neighbors of the specified UNIFORM F.

The common face R=F intersection G is an entire regular pentagon.
For an edge vw of R, let H be the third parent facet through that
edge. The face A=G intersection H of the intrinsic G contains vw,
whose length is L. Its next sides at v and w are precisely the
outgoing edges from F: in a simple four-polytope those edges are
the intersections of the other three parent facets. Both have
length L by hypothesis. Thus A has three consecutive L-edges.
The actual polygon closure lemma makes the whole A a regular
pentagon. This applies to all five edge-neighbors of R in G.

Hence G satisfies the lemma's complete flower hypothesis, and(2)
follows. The argument does not replace a face by a formal local
chain, presume any other outgoing length regular, or infer a
global120-cell decomposition. For a nonuniform dodecahedral F,
its adjacent faces in G need not be regular pentagons; then this
proof supplies no excess-five bound.

The remaining global volume obligation is unchanged: these blocks
must be combined across actual facets to prove `f(P)+W(P)>=120`.
Neither that inequality nor containment of a whole120-cell follows
from this corollary alone.
