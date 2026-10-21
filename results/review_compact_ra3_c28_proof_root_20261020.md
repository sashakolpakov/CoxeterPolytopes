# Root review: C28 regular-pentagon exclusion and RA4 consequence

20 October2026. Full written mathematical review PASS relative to the
named external classification and the hash-pinned existing geometric
inputs. Source `compact_ra3_twenty_eight_vertices_no_regular_pentagon_20261020.md`,
SHA256 `076864a5d836456e56ac377dffbcaa581a4738ec977dbf58a7905580d9838afc`.
This is a written review, not an executable replay or formal proof.

The statement is that no compact right-angled hyperbolic three-polyhedron
with28 vertices contains an entire regular right-angled pentagonal face.
I checked the following full argument.

Trivalence gives16 faces and42 edges. Polygon capacity rules out faces
of size8 or more. If a heptagon occurs, the complete equality-pattern
argument in the pinned primary manuscript gives the Loebell L7 lattice:
the unique complementary interior vertex joins the entire opposite
cycle, because a chord gives a forbidden ear. Actual metric uniqueness
realizes the base rotation, and every pentagon has a base edge with
cosh length1+2cos(2pi/7)>phi. Otherwise the polyhedron is a C28 fullerene.

I checked the primary classification premise directly in Brinkmann,
Goedgebeur and McKay, arXiv:1207.7010v2, Table2, printed/PDF page15:
the (28,16) row has fullerenes count2. The source's two explicit dual
spheres have12 degree-five and4 degree-six vertices, and distinct
hexagon-induced graphs. Thus identifying these two spheres suffices;
there is no unstated enumeration by the repository.

In the tetrahedral case the finite actual isometry group fixes a strict
interior point. The spatial Gram of the H normals has rank at most3;
distinct normals eliminate the rank-one alternative, giving the regular
tetrahedral directions. The ordered-pair stabilizer fixes exactly the
plane of v_i,v_j. I independently derived the two adjacency equations,
B^2=3/4, A=2B, a^2=7/4 and h^2=3/4. The alternative A=0 has negative
a^2. Each pentagon has the actual edge with endpoint third-normal
product-4/3, so its length differs from the regular value arcosh(phi).
These deductions start from an actual polyhedron, not formal-Gram
existence or a combinatorial proxy for metric regularity.

For the second sphere I checked the triangle/link data, the commuting
involutions and the three pentagon orbits. Any two regular-length edges
in a right-angled pentagon force regularity: the nonadjacent case gives
a third regular-length edge by the right-pentagon identity and reduces
to the proved adjacent case. The stated propagation sends either the
first or third orbit to all pentagons, giving five consecutive L-edges
on H0 and contradicting the existing actual three-edge closure lemma.

For the remaining orbit I independently derived A and C in the
nondegenerate basis(D,E,H0,H2), the H1 expansion and its beta image.
Their norm and mutual-product equations give
e^2=1+phi^-3*k^2 and2phi^-1*k*e=phi-k^2. Eliminating e gives the
displayed quartic with roots k^2=phi^4 and(6+7phi)/29. The actual
disjoint-wall gap forces k>=1, rejecting the second root, and gives
k=phi^2,e=-phi. The H Gram has the stated null relation and contains
a timelike vector, so its span is a nondegenerate Lorentz three-space.
A and C are distinct unit vectors with the same projection onto it;
its positive one-dimensional complement permits exactly those two
choices. The distinct actual normal G has the same projection by alpha,
an impossibility. All endpoint and disjoint-wall claims use the actual
links and the pinned acute-wall lemma.

For the RA4 consequence, every positive facet adjacent to a regular
dodecahedral facet contains their whole regular pentagonal ridge.
The existing excess-one exclusion, Loebell L6 case, frozen C26 theorem
and present C28 theorem therefore leave boundary excess d>=5, or at
least30 vertices. The existing farthest-facet theorem supplies a D;
connectedness supplies a D-positive adjacency unless all facets are D.
Facet capacity gives f>=65. In the non-all-D case the existing W>=17
and volume formula yield Vol>=pi^2*(65+17+16)/12=49pi^2/6. The recognized
all-D120-cell also satisfies this inequality.

No120-cell minimality or second compact volume follows. This is a
supporting successor in the existing compact RA4 research packet.
The primary14-page manuscript and every C26 frozen input stay unchanged.
The prepared graph/algebra checker must be reviewed and replayed
separately; it will not verify external classifications or actual geometry.

Primary classification: https://arxiv.org/pdf/1207.7010v2 .
