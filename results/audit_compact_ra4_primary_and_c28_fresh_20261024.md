# Fresh audit: compact right-angled four-volume and C28 successor

24 October2026. **No mathematical blocker found for the primary capacity,
excess and recognition results, the farthest-dodecahedron theorem, or the
C26/C28 exclusions and their universal bounds f>=65 and Vol>=49pi^2/6.**
This is a fresh written/source audit at the explicitly named geometric,
published-classification and recorded-computation boundaries. It is not an
external referee report, a new execution, a120-cell minimum theorem or
end-to-end Lean verification.

I read the entire1188-line primary TeX, complete C26 and C28 proofs, both
fixed checkers, the farthest-facet source and the relevant root/outside
reviews. I reconstructed the deductions below before using the older
verdicts. I did not author those sources. No scientific worker/import,
enumeration, test, Lean or TeX ran, and no canonical manuscript or frozen
input was changed.

All27 files in `results/compact_ra3_c28_successor_20261020_freeze.json`
were freshly hash/size matched. The exact inventory and four authenticated
historical execution records are in
`results/audit_compact_ra4_c28_fresh_source_inventory_20261024.json`,
SHA256 `7dcdc8e65fe6dada79d0dedf4820c814e554778ff5f241ae0c71859b5f2a137f`.
Hash authentication is not mathematical regeneration.

## 1. Actual geometry and the capacity theorem

The primary source is `manuscripts/compact_ra4_volume_constraints.tex`,
SHA256 `afe807672f3a6bbf4bb2dd95d2861b60391d7c2d7c879c6d6382e749e7ad976f`.

The normal-gap proof has the correct projection signs. For a pair with
product in(-1,0), the positive-definite Gram inverse is nonnegative.
Projecting a strict future point adds nonnegative outward-normal
multiples, keeps every other original wall inequality strict, and produces
a future timelike point on exactly the selected two walls. The nonright
actual ridge is impossible. A pair of product-1 has a nonzero feasible
future null sum, impossible in a compact polytope. Orthogonal collections
similarly give actual faces. This supplies simplicity and inherited right
angles for the faces used here, rather than assuming arbitrary Coxeter-face
heredity in a general chamber.

For a facet A, the perpendicular outgoing endpoint normal
`r_x=sinh(L_x)x-cosh(L_x)a` is correctly oriented and unit. Its products
with neighbors of A determine exactly the stated footprint, and projection
to a-perpendicular proves distinctness for different vertices. A competing
actual normal on the same perpendicular line has positive product
`cosh(L-L_x)`, contradicting acuteness. Thus a formal additional wall is not
counted as another actual facet.

The polygon-closing proof is complete. With no extra face, cyclic links
force a complementary triangulated disk on the B-cycle; an ear has total
degree four. With one extra meeting a base neighbor, it enters the disk's
boundary and every boundary vertex needs too much degree for an ear.
If the unique extra meets no base neighbor, a boundary chord cuts off an
empty disk and again gives a forbidden ear away from the chord endpoints.
The sole interior vertex consequently joins the entire boundary. Hence
equality gives the full Löbell lattice, not merely a disjoint extra face.

For a three-dimensional facet F, M=V/2+2. Its M neighbors and V outgoing
opposites are exact distinct lists. Every remaining facet has footprint at
most two among those M neighbors; a footprint of three would force it to
be a counted opposite. If a neighbor needs just one extra, polygon equality
forces that extra's footprint to be a singleton. The singleton extras
cannot be shared, and every other neighbor needs at least two incidences.
Thus there are at least M extras and

    f>=1+V+2M=2V+5.

No overlap allowance is counted twice in this argument.

## 2. Excess, volume and the complete W gap

Simplicity gives e=2v, r=v+f and six ridge corners per vertex, hence
W=6v-5r=v-5f. Within F, trivalence and Euler give
`d(F)=(V-20)/2=M-12`; summing counts each ridge twice, so sum d=2W.
The compact reflection-orbifold Euler weights give
`chi=(v-4f+16)/16`. With the stated four-dimensional Gauss--Bonnet
normalization this yields

    Vol=pi^2(f+W+16)/12.

A k-gonal shared ridge forces its neighbor's excess at least2(k-5).
Summing the original facet and its distinct neighbors gives2W>=3d(F),
including any shared ridges with exactly their multiplicity in the full
sum. This proves the neighboring-excess bound. The46/96 local cutoffs
follow with the correct parity and strictness: V=48 gives f+W>=122>120;
V=98 gives f+W>=260>256.

For W=0 every facet is the regular dodecahedron. Matching one complete
facet with the known120-cell fixes all its neighboring supporting
hyperplanes by perpendicularity along whole ridges. A whole regular
dodecahedron with that fixed pentagon is unique on the required side.
Propagation is closed under neighbors in both polytopes; connectedness
identifies their complete halfspace sets. This requires every facet to be
regular and does not propagate from one lone dodecahedral facet.

I checked the entire positive-gap dispatch, not only its final inequality:

- A nonpentagonal edge weight w forces endpoint excess at least2w.
  G-degree one is impossible. Degree two forces equal weights and a
  Löbell facet. Its regular base edges make every lateral neighbor positive.
  The L6 lateral cosines(2,sqrt2,sqrt3,sqrt3,sqrt2) have a unique longest
  edge. At an excess-two neighboring facet, that edge must have the same
  hexagonal parent. Each of the two rows contributes five extra excess
  units beyond the15-facet baseline, yielding2W>=40.
- Three consecutive regular-pentagon-length edges close a compact
  right-angled polygon to the regular pentagon: the two closing walls
  are actual orthogonal walls, and the projection lemma forbids a hidden
  truncation. Two such pentagon edges, when adjacent, also force regularity.
  Four such hexagon edges force opposite remaining positions and length2L.
- H retains **all** positive-facet adjacencies, unlike G which retains
  only nonpentagonal ridges. Actual normal geometry excludes K5 and
  induced squares in H. Eligibility in each neighbor graph uses actual
  carrier-edge/vertex incidences. The four-vertex eligibility lemma gives
  minimum H-degree five whenever G has minimum degree three.
- The subcubic nine-vertex complement lemma is exhaustive: the triangle
  branch exhausts its three outside degree slots; the trianglefree
  nonbipartite branch has an induced C5 with at least two isolated outside
  vertices. The smaller orders likewise force an independent five.
- The marked five-neighbor lemma checks universal-vertex, triangle and
  trianglefree branches. In the last branch P5 or C5 would force a
  pentagon regular while sharing a2L edge. The argument uses marked
  polygon geometry, not just an inadequate abstract graph restriction.
- The ten-vertex complement lemma handles nonadjacent exceptional
  degree-four vertices by common-neighbor count t=2,3,4. Its cross contacts
  and spare-degree bounds produce the required independent five in each
  case.

Consequently0<W<=16 forces exactly ten positive facets, with at most two
above excess three. Their H-complement satisfies precisely that last lemma
and produces a forbidden K5 in H. This proves W=0 or W>=17 over the whole
range, without a finite graph-census assumption.

## 3. Farthest dodecahedron: a universal argument

The farthest-facet source has SHA256
`d082e4a918d4504a59ac200f32bbb141f1f2cb52d701006eaea7ba5b0e5bcc22`.
For a strict interior q, projection to a maximum-height wall lies in the
relative interior of its actual facet, since all other original wall
products remain strictly negative. Its intrinsic wall heights are below1.
The associated unit directions in R3 therefore have pair products below1/2.

The displayed quartic has the correct positive Legendre expansion and
f(1)=49/25, constant coefficient7/50. The trace-free tensor contractions
give positivity of its degree1..4 kernels. The resulting bound is N<=14.
Equality at14 would force every off-diagonal product into{-1,-2/5},
contradicting the nonnegative squared norm of the direction sum. Thus N<=13.

Polygon capacity excludes a face of size at least six in such an RA3
facet. Euler then gives exactly twelve pentagons, hence the regular
dodecahedron by the named metric input. This proves that **every**
maximizing facet is D and in particular supplies a D facet in every compact
RA4 polytope. No classification of possible four-polytopes is used.

The intrinsic dodecahedral height satisfies h^2=phi/2. Antipodal face
directions imply maximum intrinsic height at least ah for time-coordinate
a>=1. Combining with the original projection gives T^2>=phi^3. Equality
forces the facet foot to be its center and every neighbor to maximize the
same original q. Connected propagation gives the all-D120-cell and its
center. The argument bounds a maximum wall height, not an inradius or an
inscribed120-cell; no such stronger conclusion follows.

## 4. External classification and metric premises

I checked the official primary [Brinkmann--Goedgebeur--McKay paper,
Table2, printed/PDF page15](https://arxiv.org/pdf/1207.7010v2).
The fullerene column has one graph at26 vertices and two at28. Its page2
defines the cubic plane graphs being counted and explains why abstract
and embedded isomorphism agree. This is an external finite classification,
not a repository-generated enumeration. The official text was retrieved;
the attempted screenshot cache-missed and is not credited as visual review.

I also checked the uniqueness clause of [Roeder--Hubbard--Dunbar,
Theorem1.4, printed p832](https://aif.centre-mersenne.org/item/10.5802/aif.2279.pdf).
It applies to marked abstract types with more than four faces and
nonobtuse angles. Every use here begins with an assumed actual compact
RA3 polyhedron with at least twelve faces. The proof uses uniqueness to
realize its combinatorial symmetries; it never infers existence from a
formal Gram matrix. The standard acute-wall, Lorentz, orbifold
Gauss--Bonnet and twelve-face dodecahedron facts remain explicit published
or foundational premises, not computational conclusions.

## 5. C26: complete classification-to-metric bridge

The C26 proof has SHA256
`1bbe6adc630ff2f44424f6ef3c573b3c71dcef4f25496def9b13fd35b2e4a03f`.
Twenty-six vertices give fifteen faces. Polygon capacity allows only
pentagons and hexagons, hence twelve and three respectively. The explicit
dual sphere and count-one input therefore exhaust the combinatorics.

On an assumed actual Q, metric uniqueness supplies its finite D3xC2
isometry group, which fixes a strict interior point by orbit averaging.
The order-three action has a nonzero horizontal orbit; otherwise the three
distinct H facets would coincide. The prime-swapping involution fixes the
horizontal plane and negates the axis. The second involution fixes its
axis sign because A0 differs from A'0. These facts justify the actual
normal ansatz and all radial/axial signs; they are not assumptions about a
molecular embedding.

I rederived U^2=2/3, B^2=1/2 and the orthogonality equation
`(2+3X)(3X-1)=2`, selecting X=(sqrt17-1)/6 with X>1/3. The A-orbit
witness has squared edge cosine `cosh(length)^2=(7+sqrt17)/4>11/4>phi^2`; the B-orbit
witness has edge cosh(1+sqrt17)/4<3/2<phi. Both actual pentagon orbits
are excluded. These quadratic radicals are forced by real metric
equations; they assume no arithmetic lattice or arithmetic ground field.

## 6. C28: all three cases and the shorter final contradiction

The C28 proof has SHA256
`076864a5d836456e56ac377dffbcaa581a4738ec977dbf58a7905580d9838afc`.
Twenty-eight vertices give sixteen faces. Polygon capacity allows sizes
at most seven. A heptagon forces the full Löbell L7 equality lattice;
its regular base edge is longer than L, so none of its lateral pentagons
is regular. Otherwise the two-fullerene input applies. The two supplied
dual spheres have different induced graphs on degree-six vertices, so
are nonisomorphic and exhaust the count-two domain.

In the tetrahedral case, the spatial H Gram at the symmetry-fixed point
has rank at most three. Its rank-one alternative would identify actual
H normals. The remaining case gives regular tetrahedral directions. The
ordered-pair stabilizer fixes their two-dimensional plane, justifying
the two-coefficient P_ij expression. I rederived B^2=3/4, A=2B,
a^2=7/4 and h^2=3/4. Every pentagon has its actual distinguished edge
with cosh4/3<phi.

In the second sphere, the two commuting involutions give the three listed
pentagon orbits. Any two L-edges in a pentagon force regularity; the
nonadjacent case follows from the right-pentagon identity before applying
the adjacent case. Starting in either outer orbit forces all pentagons
regular and gives five consecutive L-edges on H0, impossible for a hexagon.

For the middle orbit, D,E,H0,H2 are a nondegenerate Lorentz basis. Its
constraints give

    A=-phi D+phi H0+H2, C=-phi E+phi H0+H2,
    e^2=1+phi^-3*k^2, 2phi^-1*k*e=phi-k^2.

Elimination yields the stated quadratic in k^2 with roots phi^4 and
(6+7phi)/29. The compact nonadjacent-wall gap rejects the second, so
k=phi^2,e=-phi. The source's Lorentz three-space and positive
one-dimensional unit-normal fiber then genuinely allow only A,C, while
the distinct actual G has the same projection. That contradiction is valid.

A shorter ending already noted in the accepted outside review is also
valid: `C.H1=-phi*e-k=0`, but C and H1 are nonadjacent in the explicit
link table. Actual orthogonal-wall projection gives an immediate
contradiction. A consolidated edition may use this ending and preserve
the longer source/checker as supporting evidence. This shortens an
already-reviewed argument; it adds no new theorem.

## 7. Fixed algorithms: halting, coverage and existing execution

Both workers were read in full. Their rational quadratic-field operations
and physical sign tests are exact. The norm-based inverses are valid for
nonzero rational pairs in the specified irrational quadratic field.
Binary exponentiation decreases its nonnegative exponent; determinants
decrease matrix order.

These are **fixed-model checks**. Triangle/edge/link loops inspect every
element of the explicit15/16-vertex complexes. Connected cyclic links,
two incidences per edge and Euler2 establish closed triangulated spheres.
The C26 closure is a finite permutation group, capped and checked at12;
C28 explicitly checks24 tetrahedral permutations and the four second-case
permutations. The regular-pentagon closure only enlarges a subset of twelve
labels and must halt. It checks all affected labels, not sampled faces.
Exact algebra checks the normal equations, factorizations and witness
edges. The complete returned JSON is compared with the saved certificate.

None of those loops enumerates all fullerenes or establishes actual compact
realization, metric uniqueness, the all-wall lemma or the universal W gap.
Those are precisely the separately written/published bridges above.
The optional C28 `--verify-code` is provenance only and cannot replace
its full fixed-model replay.

I authenticated the four saved generation/replay ledgers against current
worker/runner/certificate hashes, exact commands, stdout lengths and the
certificate digests reported by the children. Their totals remain0.131s
for C26 and0.135s for C28; maxima are56033280 and56754176-byte sampled
tree RSS. They are historical executions, not new runs from this audit.

## 8. Universal lower bound and optional charge boundary

A positive facet adjacent to D contains their entire regular pentagonal
ridge. Excess one is impossible by polygon capacity; excess two is L6 and
has no regular pentagonal face; excess three is C26; excess four is C28
or the separately resolved L7 case. Thus its excess is at least five and
it has at least thirty vertices.

The farthest-facet theorem supplies D. Unless all facets are D, connected
facet adjacency supplies a D-positive neighbor. Capacity gives f>=65.
Outside the recognized120-cell, W>=17 then gives

    Vol>=pi^2*(65+17+16)/12=49pi^2/6.

The all-D120-cell separately satisfies both inequalities. These are
unconditional constraints in the compact RA4 class, with the stated
published inputs; they are not a120-cell minimum or second-volume theorem.

I also checked the stronger charge's partition algebra and irregular-D
incidence argument. For boundary excess at least five it gives

    Nu>=f-(26/15)W+(13/15)TI-m+(2/3)WIB+E5/5,

and TI>=2m makes the last corrections nonnegative. Hence f>=65,W<=37
forces a uniform D. The later low-volume connected-D/homology conclusions
add the separate uniform-center-cut theorem. That long auxiliary metric
proof was read but not fully reconstructed in this audit; those optional
conclusions retain that explicit premise. It is not needed for f>=65 or
the49pi^2/6 volume bound and need not expand the main manuscript.

## 9. Consolidation and reproduction requirements

The mathematical argument can be consolidated inside the existing compact
RA4 paper family: capacity/excess, universal farthest-D existence, the
small regular-pentagon exclusions, then the global lower bound. The W-gap
proof and precise external inputs must remain complete. Exploratory
homology/center-cut refinements can stay outside the main TeX.

Both fixed workers pin the old canonical TeX at its afe807... bytes.
Editing that path directly would break their existing replay contract.
Preserve that historical input path and use one clearly identified new
publication edition in the same family, or supply an explicitly reviewed
snapshot/successor mechanism. Do not silently repin the frozen workers.
Their old dated reproduction guides are themselves in the freeze; new
preferred commands should be documented without rewriting frozen evidence.

Future executions should use the reviewed assertion-safe public entry,
which retains old proof gates and locks. This audit neither reruns old
successes nor claims fresh Lean, PDF compilation, visual inspection of
every page, external refereeing,120-cell minimality or a new dimension bound.
