# Minimum-sum four-dimensional and minimum-facet five-dimensional rigidity

18 October 2026. **New written equality classifications; independent review pending.**
The following results assume no arithmetic field or restricted facet census:

* A finite-volume right-angled four-polytope with the published minimum
  facet-plus-cusp sum15 is the canonical ten-facet polytope P4.
* A finite-volume right-angled five-polytope with16 facets is the canonical
  Potyagailo–Vinberg polytope P5, with16 finite and10 ideal vertices.

No dimension bound, new example or optimality beyond the stated published
minima is claimed. The geometric inputs from the
[weighted mixed-dimension packet](right_angled_mixed_weighted_triangles_research_20261018.md)
are its elementary six-facet RA3 equality classification, the explicit
three-/four-face Lorentz models, and the distinct ridge/cusp-opposite
facet injection. Those proofs do not depend on this packet.

## A further equality rigidity lemma for dimensions four/five

We first justify the normal-product gap used below. For distinct outward
unit facet normals u,v, acuteness gives g=u dot v<=0. If -1<g<0,
their two-by-two Gram matrix G is positive definite and G inverse has
nonnegative entries. For a strict future timelike interior point x,
put b=(x dot u,x dot v), whose entries are negative, and project
p=x-[u v]G^(-1)b. Then p=x+lambda u+mu v with lambda,mu>0,
p squared=x squared-b^T G^(-1)b<0, and p has the same time orientation
as x. Every other wall r satisfies p dot r<=x dot r<0 by acuteness.
Thus p lies in the relative interior of the actual ridge cut out by u,v,
which would have a nonright angle. This contradicts the right-angled
hypothesis. Distinct facet products are consequently0 or at most-1.

The following attachment observation allows the minimum-sum equality
case in dimension4 to be classified. Suppose two triangular facets with
one finite vertex share their ideal–ideal edge inside a right-angled
three-polytope. Then that whole three-polytope
is the six-facet bipyramid.

For proof, scale the shared future null endpoints u,v so
q_1 dot u=q_1 dot v=-1. The right angle at q_1 then gives u dot v=-1.
Let n_1,n_2 be the orthogonal outward normals of the two triangles.
The finite apex of the first triangle can be written q_1=u+v-n_2;
the second is q_2=s u+s^(-1)v-n_1 for some s>0. These forms follow
from their right angles and unit timelike norm. The other four outward
normals along their sides are

\[
 u-n_2,\quad v-n_2,\quad s u-n_1,\quad s^{-1}v-n_1.
\]

Two cross products are -s and -1/s. Distinct supporting facet normals of
a right-angled polytope have product0 or at most-1, so s=1. The third
ideal vertex w=u+v-n_1-n_2 is future null and equals
(u-n_1)+(v-n_2). Every other wall has nonpositive product with this sum
by acuteness; the selected walls also do. Thus w belongs to the polytope.
The six known walls define precisely the finite-volume convex hull of
q_1,q_2,u,v,w, the earlier explicit bipyramid model. All five vertices
belong to the original polytope, so convexity identifies the whole
three-polytope without presuming any additional faces first.

We use the published bound F+c>=15 for a right-angled four-polytope,
where F counts its facets and c its cusps; see
[Potyagailo–Vinberg2005, Section4, equation(8), p.71](https://ems.press/content/serial-article-files/43032).
**Equality F+c=15 holds only for the canonical ten-facet four-polytope.**

Indeed F>=10, so c<=5. In the facet adjacency graph each three-facet
with r polygon facets and s cusps has degree r, and r+2s>=12.
Summing gives a degree sum at least12F-12c, since each ambient cusp
belongs to six facets. The graph has no five-clique, by Lorentz positive
index4, so the elementary weighted clique bound gives degree sum at
most3F^2/4. To see the graph bound, maximize the sum of t_i t_j over edges
for nonnegative weights summing to1, using a maximizer of smallest
support. Nonadjacent supported vertices can be merged onto the better
one, so the support is a clique of size at most4. Its objective is at
most(1-1/4)/2=3/8; uniform weights give the stated degree-sum bound.
For F=15-c the lower bound exceeds this upper bound for
c=0,1,2. For c=3,F=12 both equal108. Uniform vertex weights then
maximize the edge-weight quadratic form in the interior of the weight
simplex, forcing all degrees equal9 by its first derivatives. Equality
in the lower bound requires r+2s=12 in every facet, so every degree r
is even. This is impossible.

Consider c=4,F=11. Every three-facet satisfies r+s<=10 by the neighbor
and cusp-opposite injection. With r+2s>=12 this gives s>=2. The value
s=4 would require r<=6, contradicting the already proved six-facet
equality case, which has exactly three cusps. Therefore all s are2 or3.
Their sum is24, so precisely nine facets have two cusps and two have
three cusps. A two-cusp facet has r=8 and r+2s=12. Thus all its polygons
have zero excess a+b-5, where a is the number of polygon sides and b
its ideal vertices. If it had a triangle, the polygon across its
ideal–ideal edge would also have zero excess and contain both endpoints,
so it too would be a triangle with one finite vertex. The attachment
lemma would make the whole facet a six-facet bipyramid, a contradiction.
Euler gives eight finite vertices and32 polygon corners in this facet.
With triangles excluded, its eight faces are therefore all quadrilaterals.

Each of the remaining three-cusp facets has r=6 or7. The former has
six triangular faces; the latter has at least four, since its total
polygon-corner count is24. Every such triangular ridge must meet another
facet capable of containing a triangular face. There is only one such
other facet, whereas the ridge neighbors must be distinct. This excludes
c=4. Finally c=5 gives F=10, whose equality case was already identified.
This proves the stated minimum-sum rigidity without assuming that the
published attainment statement was a uniqueness theorem.

It follows that **every sixteen-facet finite-volume right-angled
five-polytope is the canonical Potyagailo–Vinberg five-polytope.**
The facet injection and F+c>=15 give the facet minimum16. At equality,
each four-facet has F+c=15 and is therefore the canonical four-polytope.
All its polygon faces have one finite vertex. Pick a finite vertex q
of the five-polytope, normalized by q squared=-1, with orthonormal outward
normals n_1,...,n_5. Its triangle graph is complete on five directions,
so its opposite triangle walls are m_ij=q-n_i-n_j. These ten normals
and the five n_i are fifteen distinct facet normals.

For each four-subset S, the forced four-face has cusp
v_S=2q-sum_(i in S)n_i. If k is omitted, exactly seven of these fifteen
walls meet v_S: n_k and the six m_ij with i,j in S. The rectangular
four-box link requires the remaining sixteenth wall r to be opposite
n_k. Thus r+n_k=L_k v_S with L_k>0. Comparing these equations for
different k gives every L_k=1: the q coefficient first makes all L_k
equal, and comparison of any n_k coefficient in two equations gives
that common value1. Therefore

\[
 r=2q-\sum_{i=1}^{5}n_i.                               \tag{1}
\]

All sixteen walls are now fixed up to isometry. In projective coordinates
their region is y_i>=0, y_i+y_j<=1 and sum_i y_i<=2. The exact certificate
checks all4368 five-wall choices, solving every independent set: its26 vertices consist
of16 finite and10 ideal vertices. This is the known five-polytope, not
a new example. Every other right-angled five-polytope consequently has
at least17 facets. If D counts canonical five-faces, this gives the
additional necessary gap

\[
 D\ge17a_5-(n-4)a_4.
\]

This higher-face gap is a further universal attachment constraint. A
complete five-clique of corner triangles alone does not force the
sixteenth wall (1); the minimum-facet hypothesis was essential there.

## Exact replay and proof boundary

The [compact certificate](right_angled_mixed_low_facet_rigidity_research_20261018.json)
and [small exact checker](../scripts/check_right_angled_mixed_low_facet_rigidity_research_20261018.py)
record the four low-cusp degree-sum cases, enumerate all4368 possible
five-wall intersections of the sixteen-wall model, and check all26
vertices exactly. The checker imports only the core packet's rational
Gaussian-elimination utility and uses Python's standard library. It
also verifies that every model facet has ten vertices, five finite and
five ideal. It does not encode the universal geometric classification.

```sh
python3 -B scripts/check_right_angled_mixed_low_facet_rigidity_research_20261018.py
```

The written Lorentz, local-to-global and graph equality proofs await
independent review. No Lean compilation or heavy numerical computation
was performed, and no frozen evidence was changed.
