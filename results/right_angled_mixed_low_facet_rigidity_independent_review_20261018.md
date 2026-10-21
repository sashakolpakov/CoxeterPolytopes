# Independent review of the RA4/RA5 equality classifications

18 October2026. **PASS after two requested geometric clarifications and a
primary-source citation correction were incorporated.** This is a full
mathematical/source review of the new
[rigidity report](right_angled_mixed_low_facet_rigidity_research_20261018.md)
and its
[exact checker](../scripts/check_right_angled_mixed_low_facet_rigidity_research_20261018.py).
I also read the supporting RA3/P3, RA4/P4, finite-star and facet-injection
arguments in the separate
[weighted packet](right_angled_mixed_weighted_triangles_research_20261018.md).
I did not execute either checker or review that packet's dimension11/12
Dehn–Sommerville calculations in this review.

The conclusions reviewed here are: equality in the published RA4 minimum
F+c>=15 forces the canonical P4, and an RA5 polytope with16 facets is the
canonical P5. Neither conclusion excludes a new ambient dimension or
constructs a new example. Literature priority for the equality statements
was not investigated.

## Primary input and the normal-product gap

I directly opened the publisher's full
[Potyagailo–Vinberg paper](https://ems.press/content/serial-article-files/43032).
The RA4 facet minimum10 and facet-plus-cusp minimum15 are Section4,
equations(7),(8), printedp71; the new report's initial Section5 reference
was corrected. Proposition3.1 identifies the known P4/P5 models' facet
and cusp counts. Proposition4.1(a) also supports the actual-facet
intersection property used in this subject. The published attainment
statements were not treated as uniqueness theorems.

The report now explicitly proves the necessary product gap. If two
distinct outward unit normals have product g in(-1,0), their positive
two-by-two Gram matrix has nonnegative inverse. Projecting a strict
interior future timelike point onto their common perpendicular adds a
positive combination of the normals. Acuteness preserves strict
inequalities for every other wall; the projected vector remains future
timelike. This produces an actual ridge with nonright angle, contrary to
right-angledness. Thus the use of products0 or<=-1 does not assume that
hyperplane intersection automatically gives an actual face without proof.

## The two-triangle attachment lemma

The normalizations and signs check. After scaling the shared future null
endpoints u,v, the first finite apex is q1=u+v-n2, with u.v=-1. The second
is q2=s*u+s^(-1)*v-n1, s>0. Right angles at the finite apices force the
unit coefficients on their transverse outward normals. The four further
side normals are u-n2,v-n2,s*u-n1,s^(-1)*v-n1. Their two cross products
-s and-1/s, together with the proved product gap, force s=1.

The requested actual convex-hull closure is now explicit. The third cusp
w=u+v-n1-n2 is future null and is a sum of two known outward normals.
Acuteness makes it feasible for all other walls; the known walls are
checked directly. The six known halfspaces give the convex hull of the
five feasible finite/ideal vertices q1,q2,u,v,w. Both inclusions therefore
hold, so the whole three-polytope is P3. No additional corner triangle
was assumed actual merely because its formal Gram diagram exists.

## RA4 equality

For a three-facet, Euler and right-angled vertex links give
p+2s=2r-4 and the excess identity

    sum_polygon(a+b-5)=r+2s-12>=0.

The facet-adjacency graph has no K5: five mutually adjacent facets would
give five orthogonal positive unit normals in Lorentz positive index4.
The weighted clique bound consequently gives degree sum<=3F^2/4.
Summing r+2s>=12 over facets gives degree sum>=12F-12c, because each
ambient cusp belongs to six facets. The arithmetic exclusions for
c=0,1,2 check directly.

For c=3,F=12, equality gives108 on both sides. Uniform positive weights
are an interior maximizer of the edge quadratic form. First derivatives
along the simplex force a regular graph of degree9. Equality in every
local excess inequality forces r=12-2s to be even, a contradiction. The
stationarity step uses a genuine global maximum and does not assume
regularity from an average degree alone.

For c=4,F=11, the distinct neighbor/cusp-opposite injection gives
r+s<=10. The RA3 minimum6 and its rigid equality case eliminate s=4;
the other inequalities eliminate s<2. Total cusp incidence24 therefore
forces nine facets with s=2 and two with s=3. Each two-cusp facet has
r=8, zero polygon excess, eight finite vertices and32 polygon corners.
A triangle forces a second one-finite/two-ideal triangle across its
ideal–ideal edge; the attachment lemma then contradicts r=8. With no
triangles, its32 corners force all eight polygons to be quadrilaterals.

A three-cusp facet has r=6 or7. It has six triangles in the first case,
and at least four in the second, from its24 polygon corners. Each
triangle is a ridge with a distinct neighboring ambient facet. Only
the other three-cusp facet can contain a triangle, which is impossible.
The distinct-neighbor assertion follows because the intersection of two
facets is one convex face, so two different ridges cannot have the same
neighbor. This completes the c=4 exclusion. The case c=5 has F=10 and
uses the separately proved ten-facet rigidity. Larger c is excluded by
F>=10. No cusp case remains uncovered.

## RA5 equality and the exact model

With16 ambient facets the facet injection and RA4 minimum force every
four-facet to satisfy F+c=15, hence to be P4. A finite vertex exists
because such a facet has finite vertices. All two-faces through it are
one-finite/two-ideal triangles. The five old orthonormal normals and ten
opposite triangle normals m_ij=q-n_i-n_j are distinct actual facets.

For each omitted index k, the actual P4 facet has cusp
v_k=2q-sum_(j!=k)n_j. Exactly seven known normals vanish there: n_k
and the six m_ij with i,j!=k. The cusp's four-box link has eight walls.
The six m's form three opposite pairs; n_k is orthogonal to all six,
so the only missing wall must be its opposite. Since there is only one
remaining ambient facet, a single normal r satisfies

    r+n_k=L_k*v_k, L_k>0, for all k.

The comparison used in the report is valid without a separate equal-width
assumption. In the independent basis q,n1,...,n5, the q coefficient makes
all L_k equal. Comparing an n_k coefficient in two equations then forces
that common value1. Thus r=2q-sum_i n_i. This determines all sixteen
actual supporting halfspaces, not merely a local subdiagram.

The checker considers all4368 five-wall subsets and solves every
independent set exactly. The region is bounded in its projective chart
by yi>=0 and sum yi<=2, so all its vertices occur in this enumeration.
The expected26 vertices consist of the origin, five unit-coordinate
cusps, ten half-coordinate triples, five half-coordinate quadruple
cusps, and five points with one coordinate2/3 and four coordinates1/3.
Their Lorentz norms give16 finite and10 ideal vertices. Checking every
model facet has5 finite and5 ideal vertices and every distinct normal
product is0 or-1 is consistent with the canonical P5 identification.
I read the rational elimination helper and found its sign convention and
singular-set handling appropriate. No arithmetic execution is attributed
to this review.

Finally the incidence consequence
D>=17a5-(n-4)a4 follows from the new equality classification and the
actual positive-dimensional face incidence: canonical five-faces have16
facets and every other five-face has at least17. It is a necessary
distribution inequality, not a dimension contradiction.
