# Independent full-source review: compact RA4 volume constraints

19 October 2026. **Full mathematical translation and figure-source PASS.**
I read all 900 lines of the new primary TeX and the
complete preparation/proof map. This was a written source review, with
primary-reference checks for the two metric inputs. No proof/checker
program, graph census, Lean, TeX build or PDF rendering ran, and neither
the new source nor any frozen predecessor was edited in this lane.

## Bound source and publication scope

* `manuscripts/compact_ra4_volume_constraints.tex`, SHA256
  `9fe28180a9097d8ee99c1898dc64e0716abf55c9c9f8cb950a6db59576d47738`.
* `docs/compact_ra4_volume_constraints_preparation_20261019.md`, read in
  full as the source's proof map and preparation record, final SHA256
  `027f063754272281b9ef47154c7b6b52b7d6f59264b47535edc426a08cb5e44d`.

The initial full read bound TeX hash `e845bc218659c59f43ae9ab9bc13fc4c0683847c6adf09436290c064a114670a`.
The author then added exactly the requested U--V edge in the fourth
figure. I reread the changed figure and refreshed preparation status;
the final 39,878-byte source bound above contains the same mathematical
text and the corrected diagram.

The source is one successor primary paper in the existing compact RA4
capacity family. It explicitly preserves the old capacity manuscript as
a predecessor. It does not claim a new compact minimum, a second-volume
gap/equality theorem, a 120-cell decomposition, a new example, or a new
dimension ceiling. The distinction between `W>=15` and the still-open
`f+W>=120` is preserved in the abstract, introduction and final section.

## Complete argument correspondence

Sections 2--4 retain the full inherited capacity proof, rather than
replace its geometric inputs by certificate assumptions. I checked:

* The two-normal M-matrix projection keeps every other wall inequality
  strict and stays future timelike. The null-sum argument excludes
  product -1 in a compact chamber. The orthogonal projection therefore
  gives actual faces, and the finite-vertex positive space gives
  simplicity.
* For disjoint facets, common-neighbor normals lie in a positive
  orthogonal complement. The normal gap makes their footprint a clique.
* The perpendicular-edge normal is
  `sinh(L_x)*x-cosh(L_x)*a`. Its incidence equivalence is exact, its
  opposites at distinct vertices are distinct, and a second actual
  normal on the same perpendicular line would have positive product
  with it. This proves the needed uniqueness, not merely existence.
* In the dual polygon argument, the complementary discs and their ears
  have the stated outside neighbors. The prescribed-edge ear version
  correctly excludes chords in the equality pattern. The unique extra
  face is disjoint from all original neighbors, and the subsequent
  chord argument gives the full Loebell face lattice without importing
  a small-polyhedron census.
* For the four-dimensional count, the lists of known faces inside each
  neighboring facet are exact. An extra facet has footprint at most
  two; an index with just one extra consumes a distinct singleton.
  Thus the incidence count gives at least M extras and `f>=2V+5`.
  No assumed disjointness of ridge contributions is used.

Section 5 derives the new neighboring-excess result directly from these
lemmas. For each distinct ridge-neighbor, its excess is at least twice
the ridge excess. The exact identity `sum_F d(F)=2W` therefore gives
`2W>=3d(F)`, with `d(F)=(V-20)/2`. The ceiling and the combination with
the facet capacity are correct.

The orbifold count in Section 5 has the correct alternating signs and
stabilizer orders:

    chi = 1-f/2+r/4-e/8+v/16 = (v-4f+16)/16.

Together with four-dimensional Chern--Gauss--Bonnet it gives
`Vol=pi^2*(f+W+16)/12`, using the independently derived `W=v-5f`.
The regular 120-cell has `(v,f,W)=(600,120,0)`. The facet double has
`(1160,226,30)`, so the two comparisons are respectively
`f+W<=120` and `f+W<256`. The exact integer/parity endpoints are:

* V48 forces `101+21=122`, hence V at most46 at or below the first
  volume. This includes both strict and equality cases, giving the
  requested 46/46 distinction correctly.
* V98 forces `201+59=260`, hence V at most96 strictly below the double.
  V96 itself only gives254, so the source correctly leaves it open.

Section 6 gives a whole-polytope zero-excess recognition. Every facet
is regular before propagation is invoked. A fixed whole pentagon and
the prescribed inward side determine the adjacent regular dodecahedron;
the perpendicular supporting normal is forced. Connectedness closes
the matched set in both facet graphs. The source explicitly rejects
propagation from just one dodecahedral facet with arbitrary neighbors.

Sections 7--8 reproduce the reviewed stronger gap proof in full. The
degree-two argument uses Loebell metric symmetry to force all lateral
neighbors positive, giving `W>=15`. The four-consecutive-edge closure
has its oriented-chain and convex-support justification. The passage
from polygon carriers to the graph on full H-neighbors uses actual
three-/four-wall intersections. The four-vertex eligibility proof then
gives minimum H-degree five under `0<W<15`.

I checked all complement cases through order nine and both branches of
the nine-vertex lemma. The triangle branch uses its three outside
degree slots and isolated vertices correctly. The triangle-free branch
uses an induced shortest odd C5, zero-or-two contacts, and the fact that
two cycle vertices cannot cover all edges of C5. The G/H distinction is
maintained throughout: induced-square exclusion applies to H, while the
edge-weight bound uses G. Nothing is silently inferred from an omitted
pentagonal adjacency. The proof excludes every integer `0<W<15` and
does not exclude or classify the endpoint W15.

## External metric references checked

The current primary text of
[Alexandrov--Bogachev--Vesnin--Egorov, Section 1.2](https://www.mathnet.ru/php/archive.phtml?jrnid=sm&option_lang=eng&paperid=9740&wshow=paper)
states the compact V20 equality recognition as the regular right-angled
dodecahedron. Trivalence converts the source's 12-face formulation to
that statement. No ideal-polyhedron assertion is being substituted.

[Roeder--Hubbard--Dunbar, Theorem 1.4, p.832](https://aif.centre-mersenne.org/article/AIF_2007__57_3_825_0.pdf#page=9)
contains the compact nonobtuse uniqueness clause. Its more-than-four-face
scope includes every Loebell polyhedron used here. Realizing the marked
cyclic automorphism by an isometry, and hence making the convex base
regular, is a valid consequence. It is not being used as an ideal or
four-dimensional uniqueness theorem.

I also checked the cited
[Potyagailo--Vinberg introduction and Proposition 4.1](https://ems.press/content/serial-article-files/43032):
they provide the stated right-angled face/simplicity and actual incidence
background. The manuscript supplies its stronger compact normal and
projection calculations explicitly. The volume argument names its
orbifold/Gauss--Bonnet input and gives the resulting face calculation;
typesetting or a saved finite computation is not presented as its proof.

## Four source-native figures

Figure `fig:ear` faithfully distinguishes the cyclic dual-neighborhood
from the complementary triangulated-disc ear. Figure `fig:capacity` is
explicitly an incidence-count illustration, not a claimed polytope or
full facet graph. Figure `fig:pentagon-chain` correctly highlights four
consecutive regular-length edges and the two adjacent positive carriers;
its caption expressly rules out Euclidean metric interpretation of the
drawn angles.

Figure `fig:nine` highlights valid independent sets in both proof
branches. The initial right panel was missing the possible edge between
the upper outside vertex U=e and the left outside vertex V. Although
its explicitly partial-diagram caption made this no proof flaw, the
drawn `U--T4` and `V--T2` otherwise formed an apparent induced 2K2.
The author added `U--V`. The final diagram is a genuine triangle-free,
maximum-degree-three, induced-2K2-free seven-vertex component with two
isolated vertices. The same five highlighted vertices are independent.
I checked this final change in the TeX source. No source edit was made
by this reviewer.

## Status

The full mathematical translation and all four figure sources pass,
with the sole requested figure clarification now applied. The preparation map correctly assigns the
proofs to their nine sections and separates external metric facts from
the new geometric and graph arguments. Page count, compiled layout,
rendered figures and any subsequent source revision require their own
build/visual record. There is no end-to-end formal geometry claim.
