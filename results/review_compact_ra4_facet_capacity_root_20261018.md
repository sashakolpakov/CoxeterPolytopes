# Root review of the compact right-angled facet-capacity argument

18 October 2026. Written mathematical and source review. This review concerns
`compact_ra4_perpendicular_facet_capacity_research_20261018.md` and its full
translation `manuscripts/compact_ra4_facet_capacity.tex`; it does not assume
a published minimum-volume classification.

The argument proves the following universal necessary condition. If a facet
of a compact right-angled hyperbolic four-polytope has V vertices, and the
polytope has f facets, then f>=2V+5. I accept the written proof, subject to
external review and the formalization limitations below.

The acute-normal gap and orthogonal projection are used with the correct
compactness hypothesis. In particular, a product -1 would give a feasible
null sum, and is excluded. The projection establishing actual faces uses
nonnegative inverse coefficients and keeps all other inequalities strict;
the proof does not confuse pairwise orthogonality with an abstract
combinatorial face.

The perpendicular-edge construction is an actual finite edge with its
other endpoint. Its opposite normal is

    r_x=sinh(L_x)x-cosh(L_x)a.

Its projection identifies x uniquely, its products with neighboring walls
identify precisely the incidences through x, and the positive product
cosh(L-L_x) excludes a second distinct opposite wall through the same
neighbor clique. These points justify both distinctness and the later
footprint restriction.

I checked the three-dimensional closing argument in the dual triangulated
two-sphere. With no extra face, the known cyclic links force the outer
B-cycle; a triangulation of its complementary polygon has an ear whose
total degree is four. With one extra face meeting a neighboring A_i, that
face is inserted in the polygon boundary and has only one outside neighbor.
Every boundary vertex again requires more than ear degree. Hence a unique
extra face must meet none of the A_i. The optional Loebell identification
is not needed in the four-dimensional count.

For a facet with M=V/2+2 ridges, the lists of M neighbors and V opposite
facets are exact and disjoint. The footprint of any remaining facet is a
clique of size at most two. Three members would define an actual vertex
and force the previously counted opposite facet. An index requiring only
one extra face consumes a distinct singleton footprint, by the equality
part of the three-dimensional lemma. All other indices need at least two
incidences, each extra serving at most two. This proves at least M extras
and f>=1+2M+V=2V+5 without rounding or unproved incidence assumptions.

The ridge excess W=v-5f is nonnegative, and the excess contributed by a
single facet is (V-20)/2. Thus f+W>=(5/2)V-5. The compact orbifold Euler
formula and four-dimensional Gauss--Bonnet give
Vol=pi^2(f+W+16)/12. The resulting limits V<=48 strictly below the regular
120-cell volume, V<=50 at or below it, and V<=104 strictly below its double
are arithmetically correct. They are local-size reductions, not a global
minimum theorem, a volume gap, or a classification of the remaining facets.

The TeX includes both geometric reductions and two explanatory incidence
figures. Its seven-page PDF compiled in two passes with no warnings;
the build ledger records hashes and removal of temporary files. Compilation
checks typesetting, not the proof. No end-to-end Lean or external human
verification is claimed. The normal-gap/projection framework, actual-edge
construction and topological ear arguments still need formalization.

**Verdict:** accept as an internally reviewed written theorem with the
scope above; do not promote either compact volume conjecture to a theorem.
