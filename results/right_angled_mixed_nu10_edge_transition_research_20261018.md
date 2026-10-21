# The nu3=10 type and finite-edge ideal-degree transitions

Date: 2026-10-18. Written local rigidity and necessary incidence identities.
No optimizer, Lean verification or dimension improvement is claimed.

## The seven-facet/three-cusp type is forced

Let F be a finite-volume right-angled three-polytope with nu=facets+cusps=10.
The established three-dimensional gap and facet inequalities imply one of
two cases: eight facets/two cusps (Q8), or seven facets/three cusps (R7).
The six-facet possibility is excluded by uniqueness of P3; nine or more
facets violate facets+2cusps>=12.

Q8 has eight finite vertices, all facets are (3finite,1ideal)
quadrilaterals, and every finite vertex has one ideal incident edge, as
proved in the polygon/three-face coupling packet.

R7 has four finite vertices and three cusps. Its complete facet list is:

    three bad triangles, one all-ideal triangle,
    three (3finite,1ideal) quadrilaterals.

To prove this, the total polygon corners are 3*4+4*3=24, so seven facets
require at least four triangles. Its total polygon excess sum(u+2v-5)
is one. Thus there is exactly one positive-excess facet, of one of the
types (u,v)=(0,3),(2,2),(4,1),(6,0); every other facet has zero excess.
Every bad triangle's ideal/ideal edge must meet a nonbad facet: two bad
triangles there would force P3 by the reviewed attachment lemma. Zero
excess nonbad facets are (3,1) quads or finite pentagons and cannot contain
an ideal/ideal edge. Hence the unique positive-excess facet must carry
distinct edges from every bad triangle. At least three triangles are bad,
and only the all-ideal triangle among the four possible excess-one facets
has three such edges. It is therefore the positive-excess facet, with
exactly three bad neighbors. The remaining three facets have twelve total
corners, so they are the stated quadrilaterals.

The finite-vertex ideal-edge degrees of R7 are exactly {0,2,2,2}. Its three
ideal/ideal edges lie on the all-ideal triangle, so the total number of
finite/ideal edges is 4*3-2*3=6. A degree-k finite vertex has at most zero,
zero, one or two incident bad triangles for k=0,1,2,3 respectively; the
last bound holds because three would give a K3 and make the whole face
P3. If all four finite vertices had k>=1 and their degrees summed to six,
the only distributions would be {1,1,1,3} or {1,1,2,2}. Either supports at
most two bad triangles, contrary to three. Thus one has degree zero. Its
three edges join all other finite vertices, so they each have degree at
most two; their total of six forces degree two at all three. The finite
edge graph is the star with this central vertex and three leaves.

This proof classifies the necessary facet/edge data. It uses no assumed
isometric uniqueness of Q8 or of R7.

## A direct finite-edge transition

At a finite vertex q of an ambient n-polytope let k(q) count its incident
edges with ideal other endpoint. Let e=qq' be a finite/finite edge.
The n-1 other directions at its ends are paired by the two-faces through e.
If such a two-face is a (3finite,1ideal) quadrilateral, precisely one of
these two paired directions is ideal: following the quadrilateral from
the ends of e reaches one finite and one ideal remaining vertex.
Consequently, if all two-faces through e are of this type,

    k(q)+k(q')=n-1.

In particular this holds throughout a configuration whose three-faces
are all P3, Q8 or R7. Its two-faces are bad triangles, all-ideal triangles
or (3,1) quads; a finite/finite edge cannot lie in either triangular type.

This exposes a defect at the lower-bound equality of the preceding packet.
Equality in its local gamma bound permits only k=4 with K4 on those four
ideal directions, or k=n with a complete graph (for n=11,12). Any k=4
vertex has finite edges, whose neighbors would have k=n-5, namely six or
seven. These are outside the equality set. Thus the improved lower bound
is not a realizable equality configuration. This observation is a
transition obstruction, not a uniform numerical gap above that bound.

## A potential whose error is controlled by nu-excess

Define

    Phi_n(k)=3 binomial(n-k,3)-binomial(k,2)(n-k)
            =((n-2)/2)(n-k)(n-1-2k).

For a single three-face let A_j count its finite vertices with j ideal
incident edges. P3 has only j=3; Q8 only j=1; R7 has degrees {0,2,2,2}.
All three therefore satisfy 3A0-A2=0. For any other three-face, nu>=11.
Put e=nu-10. Then

    -6e <=3A0-A2 <=30e.

Here is a proof using only elementary graph bounds. Write p for the finite
vertices and s for cusps; p+4s=2nu-4. If L counts finite/ideal edges, then
L=4s-2E_II=A1+2A2+3A3, where E_II is the number of ideal/ideal edges.
In particular A2<=min(p,2s)<=(2nu-4)/3<=6e, proving the lower bound.

For the upper bound split by s. If s=0, facets+2s>=12 gives nu>=12 and
3A0-A2=3p=6nu-12<=30e. If s=1, each finite vertex has at most one edge
to the sole cusp, so A1=4, A2=A3=0 and the potential is
3p-12=6nu-36<=30e. If s=2, A3=0 and L>=6 since there is at most one
edge between the two cusps. Hence
3A0-A2=3p-3A1-4A2<=3p-2L<=6nu-48<=30e.
Finally for s>=3 the simple planar graph on the cusps has
E_II<=3s-6, so L>=12-2s. Since
3A1+4A2+3A3>=L,

    3A0-A2 <=3p-L<=6nu-24-10s<=6nu-54<=30e.

The coefficient thirty is attained by a compact twelve-facet
three-polytope, if one uses the known right-angled dodecahedron; attainment
is unnecessary for the inequality.

Now sum over ambient three-faces. At an ambient finite vertex with k ideal
directions, the numbers of incident triples with zero and two ideal
directions are binomial(n-k,3) and binomial(k,2)(n-k). Thus, with the
nonnegative excess

    E=sum_non-P3(nu-10)
     =B-10a3+(n-2)a2+binomial(n-1,2)2^(n-3)c,

we obtain the universal error budget

    -6E <=sum_q Phi_n(k(q)) <=30E.

The potential also has a finite-edge interpretation:

    sum_q Phi_n(k(q))
       =-(n-2) sum_finite_edges_qq'[k(q)+k(q')-(n-1)].

It follows directly from the displayed polynomial formula and finite-edge
handshaking. At E=0 the individual three-face classification gives the
stronger edge-by-edge complementation statement, not merely this summed
identity. This supplies a controlled way to compare ideal degrees across
finite edges when larger three-faces occur.

## Scope and next use

The rigid R7 edge data and the potential budget couple finite-edge
transitions to the actual three-face types. They target the surviving
multipartite/equality patterns more directly than an arbitrary graph
mixture. The next unresolved step is a sufficiently strong quantitative
cost for departures from these local transitions. No claimed dimension
exclusion or additional optimized minimum follows here.

Replay the tiny supporting arithmetic with
`python3 scripts/check_right_angled_mixed_nu10_edge_transition_research_20261018.py`.
The finite degree cases and potential coefficients are exact; the
universal graph/geometry arguments above are written proofs.
