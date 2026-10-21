# Distance15 bounds for simple cusps and compact facets

Let P be a finite-sided, finite-volume hyperbolic Coxeter polytope in H^n.
The following restricted conclusions hold:

* If every ideal vertex is simple, then **n<=125**.
* If P has a compact facet, then **n<=127**, with no simplicity assumption
  on its other ideal vertices.

Finite vertices are automatically simple. Neither hypothesis is asserted
for every finite-volume chamber; the full unrestricted bound163 is
unchanged. No novelty or current-best claim is made. The first assertion
uses a simple projective polytope, not an assumption that the chamber is
compact in hyperbolic space. The compact facet in the second assertion
need not itself have Coxeter angles.

The [exact checker](../scripts/fresh_simple_cusp_distance_bound_20261017.py)
and [certificate](fresh_simple_cusp_distance_bound_20261017.json) refine
the earlier [compact-facet distance17 argument](fresh_nikulin_narrow_20261015_compact_facet.md).
They use the frozen complete high-rank quasi-Lannér catalogue to reduce
the quadrilateral distance from17 to15. The all-rank finite/affine graph
bound and the passage to ideal polygon corners are proved below.

## Local face stars, including ideal corners

Every positive-dimensional k-face belongs to exactly n−k facets. The
facet normals of its supporting walls have positive-definite Gram,
because the relative interior contains a hyperbolic point. These
incidence facts, and the product-of-simplices description of an ideal
link, are stated for acute chambers in
[Potyagailo–Vinberg, printed p.64](https://ems.press/content/serial-article-files/43032).

Let Q be an actual triangular or quadrilateral two-face. Write W for its
n−2 supporting facet normals and mark its side facets black. Each side
has exactly one original ambient facet beyond W. Its edge interior gives
W plus that black node positive definite. At a finite polygon corner,
W plus the two adjacent black nodes is positive definite. At an ideal
corner it is positive semidefinite. Under the simple-cusp hypothesis
these n normals are precisely the complete incident affine vertex
diagram. The compact-facet application has only finite corners.

The full face star is indefinite: after projecting the side normals off
the positive span W, they bound a nondegenerate finite-area polygon in a
Lorentzian three-space, which they span. Every principal block has
negative index at most one, since it is a block of the ambient Gram.

We use two acute-incidence facts from the inherited compact-facet proof.
A positive-definite set of original wall normals has an actual common
face, obtained by projecting a strict timelike interior point; the
inverse of a positive-definite matrix with nonpositive offdiagonals is
entrywise nonnegative. A singular positive-semidefinite set contains a
connected affine component whose positive null combination is a nonzero
future null vector. It is orthogonal to the whole PSD block and pairs
nonpositively with all other walls, hence gives an actual ideal point.
Nonzeroness and orientation follow from its strictly negative pairing
with the common interior point.

In a quadrilateral, the two **closed opposite sides are disjoint**.
This remains true when some or all four vertices are ideal: the four
vertices are distinct, and only consecutive sides have a common endpoint.
Consequently, for either opposite black pair U, the block W+U is
indefinite. Positive definiteness would give a finite common point of
those opposite sides; singular PSD would give an ideal common point.
Both contradict the actual closed polygon incidence. Compactness of Q
is not used in this step.

## Minimal diagrams and the radius7 certificate

Every inclusion-minimal indefinite principal diagram L is connected and
has all proper principal blocks PSD. If a proper connected affine block
C occurs, its positive kernel pairs nontrivially with an adjacent outside
node; C plus that node is already indefinite. Minimality forces that to
be all of L. Thus L has the Lannér/quasi-Lannér classification scope,
with order at most10. Strict dotted pairs give only an order-two case.
The classification input and primary references are recorded in the
[frozen high-rank census](fresh_high_rank_quasilanner_census_20261017.md)
and the inherited compact-facet proof; this package does not repeat the
large classification enumeration.

Choose minimal indefinite diagrams L in W+U and M in W+V, where U,V
are the opposite black pairs. Each contains both its black nodes. On
deleting either black node the remaining subdiagram lies in W plus one
side facet, so that deletion is **positive definite**, not merely PSD.
Therefore both marked nodes are in the catalogue's `finite_deletions`.

For every such marked pair in all26 high-rank types, the checker performs
a multi-source breadth-first search from the pair. There are384 pairs:

| Maximum distance to the marked pair | Number of pairs |
|---|---:|
|1|10|
|2|123|
|3|143|
|4|61|
|5|31|
|6|12|
|7|4|

The four radius7 cases are L23 with blacks {0,1}, farthest node8; and
L25 with blacks {1,2}, {1,3}, or {2,3}, farthest node9. Every literal
distance is saved. A connected diagram of order at most5 has radius at
most3 from a two-node set: a shortest path cannot pass the other marked
node before reaching the set. This treats all smaller types without a
label cutoff. Hence the radius7 bound is complete for every possible L
and M.

The two minimal indefinite diagrams cannot be disjoint and orthogonal,
as that would give two negative directions. They share a white vertex or
have an edge between them. In the overlap case the marked sets are joined
by a path of length at most14. In the joined case it is at most
7+1+7=15. Taking a shortest path from U to V ensures that all internal
vertices are white. Its endpoints are consecutive sides of Q, so this
path lies in their actual corner diagram.

For a triangle, every indefinite subdiagram contains all three black
nodes, since any black deletion lies in a PSD corner. A minimal one has
at most10 nodes. In a spanning tree, take the minimal subtree joining
the three marks. At least two marked leaves can be deleted separately
while the other two marks stay joined. Each resulting path has length
at most8 and belongs to one actual corner diagram. This argument uses
PSD corner blocks and is therefore valid at ideal corners as well.

Give a polygon corner weight1 when its two black nodes are at graph
distance at most15 in the full incident vertex diagram; otherwise give
weight0. Triangles receive at least2 and quadrilaterals at least1.
Larger polygons need only nonnegative payment. Thus every two-face with
k sides receives at least5−k.

## All-rank finite and affine vertex costs

For a connected graph with N vertices let p_k be the number of unordered
pairs at graph distance k. The unweighted underlying graphs of the
classical finite and affine Coxeter diagrams are paths, cycles, or the
single- and double-fork trees below. Edge multiplicities do not change
graph distance.

* A path has p_k=max(N−k,0).
* A cycle has p_k=N for 2k<N, p_k=N/2 for 2k=N, and zero otherwise.
* A path on m>=2 vertices with two extra leaves at one end has N=m+2 and
  p_k=max(m−k,0)+2[1<=k<=m]+[k=2].
* A path on m>=2 vertices with two extra leaves at each end has N=m+4 and
  p_k=max(m−k,0)+4[1<=k<=m]+2[k=2]+4[k=m+1].

Each formula follows by separating path–path, leaf–path, and leaf–leaf
pairs. In every displayed case p_k<=N. The affine D4 four-leaf star is
the m=1 exception and is treated separately. Every remaining exceptional
finite/affine diagram has N<=9, so even its total number of pairs is at
most4N. The finite dihedral diagrams and affine A1 are paths of order2;
the other small F/G/H types are paths or included exceptions. These
classification alternatives cover all connected finite and affine
diagrams.

It follows in every case that the number of pairs at distance at most15
is at most15N. The same holds for a spherical forest by adding component
orders; pairs in different components have infinite distance. The
checker independently counts245 classical graphs and seven exceptional
graphs by BFS and compares the formulas. The explicit formulas and the
finite exceptional classification, not the tested rank range, prove
the all-rank statement.

At a finite vertex the full spherical forest has order n. At a simple
ideal vertex the Euclidean link is one simplex, so its full incident
diagram is connected affine of order n. In both cases all planar-angle
weights sum to at most15n. No such order-n assertion is made for a
nonsimple ideal vertex.

## Strict h-vector thresholds

The projective closure of P is a bounded simple n-polytope when all its
ideal vertices are simple. The inherited nonnegative symmetric h-vector
argument gives the strict average corner-payment lower bound

\[
 b_d=\begin{cases}d(d-6)/8,&d\text{ even},\\
 (d-1)(d-5)/8,&d\text{ odd}.\end{cases}
\]

Its strictness follows from the positive endpoint h-coefficients. Applied
to the whole chamber, the payment bound gives b_n<15n. For n=126+2t and
n=127+2t respectively, eight times b_n−15n is

\[
 4t^2+252t,\qquad 4t^2+256t+132.
\]

Both are nonnegative for t>=0. Equality in the first at t=0 still
contradicts the strict lower bound. Thus n<=125. Without the new radius
certificate, distance17 would similarly give the weaker n<=141.

If P instead has a compact facet F, put d=n−1. F is a simple d-polytope
whose vertices are finite ambient vertices. Its two-faces are actual
two-faces of P, and their weights count only a subset of the pairs in
the ambient rank-n spherical vertex forest. Hence b_d<15(d+1). At
d=127+2t and d=128+2t, eight times the difference is

\[
 4t^2+256t+12,\qquad 4t^2+260t+136,
\]

strictly positive. Therefore d<=126 and n<=127. This replaces the older
distance17 compact-facet bound143 in this restricted branch only; the
other unrestricted finite-volume branches are not changed here.

## Comparison, verification, and scope

Felikson–Tumarkin prove the much stronger n<=8 result for
**simple all-ideal** Coxeter polytopes. Their
[primary paper, §§1–3](https://arxiv.org/pdf/math/0502413)
uses the absence of elliptic order-n diagrams, together with restrictions
on ideal triangles and quadrilaterals. That absence fails when finite
vertices are allowed, so their proof is not imported to the mixed case
here. This report makes no survey-based claim that125 is a best published
bound for the restricted class.

Generation and a fresh saved replay pass. Syntax, all nine direct source
hashes, and report links pass. Two internal mathematical/source reviews
found no gap; the [separate independent review](fresh_simple_cusp_distance_bound_20261017_review.md)
records its scope and does not claim a duplicate arithmetic replay.
The checker binds the frozen complete catalogue and
the preceding geometric proofs; its new arithmetic is the384 exact
distance calculations, independent classical graph counts, and rational
threshold identities. No coefficient field, arithmeticity, optimizer,
finite-tail cutoff, or induced-face Coxeter assumption is introduced.
