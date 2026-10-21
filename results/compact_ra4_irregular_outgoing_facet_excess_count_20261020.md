# Irregular outgoing edges consume nonpentagonal ridge excess

20 October 2026. Written structural reduction, awaiting independent review.
No jobs, searches, builds or formal verification. No universal 120-cell
minimum or containment theorem is claimed.

Let P be compact and right-angled in H4. Write f for its facet count,
W=sum_R(k(R)-5) for its pentagon excess, and N_+ for its number of positive
facets. A zero-excess facet is the regular right-angled dodecahedron D.
Put phi=(1+sqrt(5))/2 and L=arcosh(phi).

Call such a D **uniform** if all twenty edges leaving its vertices have
length L, and **irregular** otherwise. Let N_u,N_i count these two kinds
of dodecahedral facet. The word irregular refers only to these outgoing
lengths; the intrinsic facet itself remains regular.

**Proposition.**

    3 N_i <= 4 W,
    N_u >= f-N_+-floor(4W/3) >= f-7W/3.                 (1)

If W<20, then every positive facet has excess at least three, and hence

    N_u >= f-2W.                                       (2)

In particular W<=19 guarantees at least seven uniform dodecahedral
facets, since the inherited facet capacity gives f>=45.

## 1. Pentagon equations on a dodecahedral facet

Fix an actual D with outward unit normal n. For a vertex v of D, let s_v
be its outgoing edge length. Its opposite endpoint wall is

    r_v=sinh(s_v)v-cosh(s_v)n.

For an edge vw of D, let A,B be the other two parent facets through it.
Then R=A intersect B is the unique perpendicular two-face of P meeting D
along this edge. If R is a pentagon, its five boundary walls are precisely
D, the two endpoint carriers in D, and the two outgoing endpoint walls
r_v,r_w. Thus r_v and r_w meet at a right angle. Since -<v,w>=phi,

    0=<r_v,r_w>
      =cosh(s_v)cosh(s_w)-phi*sinh(s_v)sinh(s_w).

Consequently, for each pentagonal perpendicular ridge,

    coth(s_v)coth(s_w)=phi.                            (3)

There is no assumption that s_v lies above or below L. Define the real
defect at every vertex of D by

    u_v=log(coth(s_v))-(1/2)log(phi).

Then u_v=0 exactly when s_v=L, and (3) says u_v+u_w=0 across every
pentagonal perpendicular ridge. In particular a pentagonal such ridge
cannot connect a regular outgoing vertex to an irregular one. This is
also the metric propagation expressed by the two-adjacent-L-edge
right-pentagon lemma, without presuming the ridge is regular in advance.

For completeness, the five-wall assertion uses actual incidence, not
only formal normals. R is a compact right-angled polygon. Its side in D
has endpoints v,w. The next sides at those endpoints are the two outgoing
edges, whose next actual side walls are r_v and r_w. A pentagon has exactly
one remaining vertex, so these last walls are adjacent there. Compact
right-angled P is simple, and the five sides are distinct.

## 2. An irregular D needs at least three nonpentagonal incidences

Mark the vertices with u_v nonzero. If this is a nonempty proper set, every
edge crossing between marked and unmarked vertices has nonpentagonal
perpendicular ridge, by (3). There are at least three crossing edges.
Here is an elementary planar proof of the latter fact. In the dual
icosahedral graph, the subgraph formed by the crossing edges has even
degree at each vertex: going around the corresponding pentagon of D
changes between marked and unmarked an even number of times. This is a
nonempty even graph, so it contains a cycle. The icosahedral graph is
simple, and therefore that cycle has at least three edges.

If all twenty vertices are marked, the subgraph formed by pentagonal
perpendicular ridges is bipartite: along each edge a nonzero u changes
sign. Equivalently no odd cycle of such edges is possible. Each of D's
twelve pentagon boundaries must therefore contain an edge with a
nonpentagonal perpendicular ridge. Each edge belongs to two pentagons,
so at least six distinct edges have that property.

Thus every irregular D has at least three distinct nonpentagonal
perpendicular ridges. A ridge cannot account for two different edges
of this same D: its intersection with D is a face of the ridge, and
containing two edges would force the whole ridge into D, contrary to
its being nonpentagonal.

## 3. A nonpentagonal ridge has limited dodecahedral incidence

Let R be a nonpentagonal k-gonal ridge. If a dodecahedral facet meets R
along an edge, that edge has length L. Different dodecahedral facets give
different edges: in a simple four-polytope the edge of a two-face has
exactly one parent facet not containing that two-face.

The inherited three-edge closure lemma says that a compact right-angled
polygon with three consecutive L-edges is the regular pentagon. It uses
the actual orthogonal closing lines and acute projection, so another
wall cannot truncate their closing vertex. Therefore no three consecutive
edges of this R can all be carried by dodecahedral third facets.

Counting the k cyclic windows of three edges, each of which has at most
two such carriers, shows that their number c(R) satisfies

    c(R)<=floor(2k/3)<=4(k-5),         k>=6.             (4)

This counts only facet intersections that are edges. A facet meeting R
only at a vertex does not enter c(R) and contributes nothing to the count.

## 4. Global count and its current consequence

Count incidences between irregular D-facets and their nonpentagonal
perpendicular ridges. Section2 gives at least3N_i. Each is among the
incidences counted by c(R), so (4) gives

    3N_i <= sum_(R nonpentagonal)c(R)
           <=4 sum_R(k(R)-5)=4W.

All nonpositive facets are D, so N_u=f-N_+-N_i. Every positive facet has
excess at least two, and the sum of all facet excesses is2W; hence
N_+<=W. This proves (1).

If a positive facet has excess two, the small-degree lemma gives a bound
of one on each incident nonpentagonal ridge weight, and excludes a single
such ridge. It therefore has exactly two hexagonal neighbors in the
nonpentagonal adjacency graph. The reviewed degree-two lemma forces
W>=20. Thus when W<20 every positive facet has excess at least three,
N_+<=2W/3, and (1) sharpens to (2).

Combined with the new center-cut theorem, a polytope with W<=19 and volume
at most the 120-cell volume has at least seven actual uniform D-facets,
and the matched 120-cell center of EACH lies in P's interior. This does
not assert that those centers are distinct or that whole matched
120-cells lie in P. Nor does low volume currently imply W<=19. These are
the remaining global gaps; neither is replaced by a count assumption.

## Inputs, unchanged

* `manuscripts/compact_ra4_volume_constraints.tex`, SHA256
  `afe807672f3a6bbf4bb2dd95d2861b60391d7c2d7c879c6d6382e749e7ad976f`:
  actual compact wall conventions and simplicity, Lemmas `small-degree`,
  `three-edges`, and the degree-two W>=20 result, zero-excess metric
  recognition, f>=45 and the facet-excess sum2W.
* `results/compact_ra4_uniform_outgoing_center_cut_volume_cost_20261020.md`,
  SHA256 `78135c4cca5d64c583727d612d5f22cf8f99f50a6a9b0aedbae17780e01c2b2c`:
  used only for the final low-volume center consequence, not to prove (1).

The exact defect equation, planar cut argument and global incidence
inequality above are the new reduction. No extra three-dimensional
classification or global garland hypothesis is introduced.
