# Independent written review: excluding compact RA4 excess sixteen

19 October 2026. **PASS, relative to the explicitly inherited geometric
inputs.** Full source and paper-mathematics review only; no graph census,
matrix program, scientific checker, Lean or TeX job ran. No reviewed or
frozen source was edited. No new PDF or formal-proof claim is made.

Reviewed source:
`results/compact_ra4_sixteen_excess_equality_exclusion_research_20261019.md`,
SHA256
`8caf24e085b204375e3a36b108513e55d8fb0cc66559b44f3e3caf1bd17977e8`.
All five sections were read, including the exact scope of its inherited
local link exclusion. The unchanged primary manuscript remains
`manuscripts/compact_ra4_volume_constraints.tex` at SHA256
`20e5b4a3aff508520b7db3c4ae1a8904c0aa2d3299b4164a657bb6432433c49e`;
it states W=0 or W>=16. This review does not silently alter that source.

## L6 metric and actual-edge propagation

The stated reflection of the Loebell face lattice preserves every
listed adjacency: p_i goes to p_(-i), q_i to q_(-i-1), with both bases
fixed. Three-dimensional uniqueness makes it an isometry of the compact
right-angled realization. On p_0 it exchanges the two p-neighbors and
the two q-neighbors and fixes the base edge. Thus the cyclic side pattern
really is (b,a,c,c,a), with b on the hexagonal base. The same marked
pattern applies to all lateral faces by the Loebell symmetries.

The regular-base formula gives cosh(b)=2. The right-angled pentagon
identity used in the report is correct with its cyclic indices: the
two far supporting wall normals after perpendicular sides u,v have
product -sinh(u)sinh(v), and their connecting perpendicular is the
opposite side. Substitution gives sinh(c)^2=2 and
cosh(a)=sinh(c)sinh(a); hence cosh(c)=sqrt(3), coth(a)=sqrt(2), and
cosh(a)=sqrt(2). This also satisfies the remaining cyclic pentagon
relations. Therefore the hexagon-base edge b is uniquely longest.

That uniqueness is sufficient for the crucial propagation. If the
shared pentagon R=F intersect A_i belongs to two Loebell L6 facets,
the base edge selected from F and the base edge selected from A_i
must be the SAME edge of the actual pentagon. In the compact simple
four-polytope, its relative interior lies in exactly three facets.
They are F,A_i and the original base neighbor B_j. Since F and A_i
meet pentagonally, B_j must be the parent of A_i's hexagonal base;
consequently A_i intersect B_j is an actual hexagon. This uses neither
an unproved four-dimensional metric classification nor a formal graph
adjacency in place of a face.

At W=16, a degree-two G facet must have k=6. Its15 distinct selected
facets contribute at least30 to sum d=32. Every selected facet has
d<=4, and at most two of the twelve lateral neighbors can have d>2.
A d=2 facet has two unit-weight G edges, so the inherited equality
case makes it L6. Each base neighbor meets six selected lateral
neighbors; at least four have d=2. The preceding propagation gives
four distinct unit-weight G edges into that base neighbor, besides
its original edge to F. Thus d(B_j)>=5, contradicting d(B_j)<=4.
An additional positive facet outside the15 consumes the same slack
and cannot evade this argument. Every degree-two case is excluded.

## The ten-node complement lemma

I checked every case, independently of a finite graph enumeration.
Deleting the sole degree-four vertex gives the inherited subcubic
nine-node case. With two adjacent degree-four vertices, deleting one
also lowers the other's degree and gives that case.

For two nonadjacent degree-four vertices u,v, the partition X,A,B,S
has the stated sizes t,4-t,4-t,t. Missing any A--B edge gives an
induced2K2, so all such edges exist. If A is nonempty, their endpoints'
degree bound forces t>=2; hence only t=2,3,4 remain.

* At t=2, A and B are saturated, so A has no internal or S edges.
  An S edge would need to cover all six distinct vertices of
  N(u) union N(v), with only four remaining endpoint degrees.
  Thus S is independent and v together with A and S gives five.
* At t=3, the same S-edge covering argument needs five contacts but
  allows only four. Hence u,v and the three vertices of S are independent.
* At t=4, an S edge uses all four available endpoint degrees on X.
  Neither endpoint has another S neighbor, and any second S edge
  would consequently induce2K2 with it. Thus S has an independent
  triple; add u,v.

The S-edge covering step follows from the explicit four vertices
y,z,u,x or y,z,v,x, since S has no u- or v-contact. It is not a
stronger unproved common-neighbor condition. All constructed independent
sets have five distinct vertices.

## Final composition and scope

With degree-two G vertices excluded, G has minimum degree3. The exact
sum32 gives N<=10, while the inherited H-degree/small-graph theorem
excludes N<=9. Thus N=10 and at most two facets have d>=4.

The previously reviewed five-neighbor exclusion for d=3 is genuinely
local. Its proof uses only that facet's three hexagons, its five
positive neighbors, the zero-neighbor regular pentagons, actual-wall
closure and the exhaustive five-node marked graph argument. It does
not require its neighbors to have excess3, or the ambient total to be15.
Its use here is therefore valid. All other facets have H-degree>=6;
H has minimum degree5 and at most two degree-five vertices. Its
ten-node complement meets every hypothesis of the new lemma, so has
an independent five, giving the forbidden K5 in H. This closes the
remaining W=16 case.

The resulting written successor is **W=0 or W>=17**, with the same
120-cell zero-case recognition. No correction is requested. The result
retains the earlier acute-wall, compact simplicity, RA3 uniqueness and
dodecahedron-recognition inputs. It supplies no universal minimum-volume
or second-volume theorem, no garland/tessellation classification, no
new dimension ceiling and no priority claim. The manuscript20e5 source
and its existing review remain preserved separately.
