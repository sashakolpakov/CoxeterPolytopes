# Integrating the reviewed right-angled cusp-capacity theorem

19 October 2026. Read-only integration plan. No TeX, frozen report,
reproduction input, catalogue, checker or Lean source was edited, and no
build or proof job was run. Source locations refer to the versions inspected
today.

## Recommended homes

Keep the existing manuscript dependency order:

    ideal4 structure / small-facet exclusion + carriers
                         -> five_constraints -> ideal6

The complete new geometric proof and the ideal-five size theorem fit in
[the five-dimensional constraints manuscript](../manuscripts/right_angled_five_dimensional_constraints.tex),
replacing its present Section2 rather than adding a manuscript. The existing
[ideal-six manuscript](../manuscripts/ideal_right_angled_six_face_constraints.tex)
should quote that theorem and retain its own transport proof. The catalogue
already has ideal6 depending on five_constraints and carriers, so this
placement avoids reversing the dependency or introducing a citation cycle.
The universal infinite-volume attachment remains a separate theorem in
Section1; the new finite-volume lower bounds do not complete that chamber.

Canonical new written proof:
[capacity and size report](../results/right_angled_outside_cusp_capacity_20261019.md),
SHA256 ffbae765aa2594283e347c75e2a2bf41bfbe5ded217c66c3e53746c641c7ccb0.
Its [independent full review](../results/review_right_angled_outside_cusp_capacity_finite57_20261019.md)
covers the full geometric and arithmetic composition, with no execution or
full Lean claim.

## A necessary translation gap before publication

The new carrier input is:

    kappa(Q)>0 -> some three-facet R has at least11 vertices.

Its source is the new
[Figure7 exclusion](../results/ideal_ra4_ten_vertex_facet_exclusion_research_20261019.md),
the preceding all-antiprism extension, the complete stated published
classification through ten vertices, and the Oct/D-chain theorem.

The current
[antiprism appendix](../manuscripts/ideal_ra4_square_antiprism_extension_appendix.tex)
starts at line6 with the O/D/A4 extension; it is an unbuilt input fragment,
not the full new all-antiprism/Figure7 proof. The unified
[ideal4 manuscript](../manuscripts/ideal_ra4_volume_and_structure.tex)
has the Figure7 combinatorics in its equality-only reduction around lines
388-410, and its O/D chain theorem at lines84 and432 onward. Those statements
must not be silently relabeled as the full new small-facet theorem.

For the next separately checked manuscript revision, expand the existing appendix
home to contain the general-antiprism and Figure7 exclusions, with one final
at-most-ten-vertices corollary. Include or cite that completed revision from
the primary ideal4 exposition. Until this is done, the exact dependency is
the stable written research proof, not a nonexistent TeX theorem. Do not
edit the checked unified source merely to make the citation
look complete.

## Planned five-dimensional manuscript edits

Current locations in
[right_angled_five_dimensional_constraints.tex](../manuscripts/right_angled_five_dimensional_constraints.tex):

|Location|Replacement or addition|
|---|---|
|Abstract, lines22-27|State59 facets,180 vertices and strict ideal-five fraction4/43 in dimensions11/12; replace the second-volume dependency with the new small-facet bridge.|
|Section2, lines123-155|Replace the old56 proof by the compact sequence below.|
|Section3, lines177-187|Replace the old coefficients16+40 by16+43; the nonideal fraction is39/43 in dimensions11/12 and29/43 in9/10.|
|Table, lines205-209|Use ideal fractions14/43 in9/10 and4/43 in11/12. Updated vertex ratios in11/12 are13312/20669 and6656/5637; either update9/10 with a newly checked rational row or omit those optional ratios pending that check.|
|Final status, lines219 onward|Separate written universal geometry from finite arithmetic checks and report the new review scope accurately.|

A compact Section2 needs only:

1. The acute projection lemma, with all other actual wall inequalities
   strict, and its zero-or-at-most-minus-one product consequence.
2. Outside-cusp capacity at most two, including the center-at-box-vertex
   contradiction by actual finite-vertex projection. Double counting gives
   V(P)>=(1+2^(n-3))V(F).
3. The existing actual ridge-neighbor/cusp-opposite injection. For the
   carrier's R with m>=11, ideal3 Euler gives F(R)=m+2, so the two independent
   consequences are V(Q)>=3m>=33 and F(Q)>=2m+3>=25.
4. Ambient injection yields F(P5)>=59. Write the ideal-five Euler identity
   g=3v+3f-6 and strict incidence2g>24f, hence v>=3f+3>=180.

The published ideal4 facet minimum24 is used in step4. Its uniqueness is
unneeded: the actual carrier already has at least25 facets. The old internal
second-volume/equality theorem is not needed by this new size proof.

## Planned six-dimensional manuscript edits

Current locations in
[ideal_right_angled_six_face_constraints.tex](../manuscripts/ideal_right_angled_six_face_constraints.tex):

|Location|Replacement or addition|
|---|---|
|Abstract, lines15-23|State240 facets,4980 vertices, V>=ceil(5643F/272), and ideal-six fraction below5/71.|
|Section1, lines27-51|Quote the new ideal-five59/180 theorem with the precise small-facet/capacity dependencies; retain a short statement of the actual-facet injection rather than duplicating the capacity proof.|
|Section2, lines53-77|Use F6>=1+59+180=240. The simple ten-facet incidence gives V6>=18F6>=4320, a baseline subordinate to4980. It need not remain a separately advertised main theorem.|
|Sections3-4, lines79-230|Keep the universal cusp-frame, good-edge, reciprocal103/132 and actual transport proof unchanged. Its reusable endpoint is a4-F<=136V/99.|
|Final vertex proof, lines231-235|Use2a4>=59F, then(57/2)F<=136V/99 and272V>=5643F. WithF>=240, exact integer rounding gives4980.|
|Mixed boundary section, lines259-265|Replace27N6+200I6 by27N6+213I6, giving I6/N6<5/71. Cross-reference the updated five-face fraction4/43.|
|Status, lines267-279|The new carrier bridge is an internal proof. Preserve the explicit external minimum24 and NK/Dufour boundaries, plus the pending geometric Lean scope.|

The number103/132 alone is not the interface sufficient for4980: the actual
edge reciprocity, triangle count and face identities proving136/99 must
remain in the manuscript.

## Reproduction and verification updates to schedule later

Preserve all historical56/170,227/3859 and4462 certificates as valid
predecessors. Do not change their inputs or imply they reproduce the new
bounds. One new small exact arithmetic check can cover the entire new
59/180/240/4980 and mixed-fraction composition, with source pins to the
new proof and review plus the unchanged136/99 transport packet. Its ledger
should say that actual geometry is supplied by the written proof; checking
constants is not a geometric census.

The later documentation/cross-reference changes are localized:

- [Five-dimensional reproduction map](right_angled_five_dimensional_constraints_reproduce.md):
  replace the main-size dependency chain and clearly distinguish predecessor
  replay from new-composition replay.
- [Six-dimensional reproduction map](ideal_right_angled_six_face_constraints_reproduce.md):
  update the current theorem rows and the272V>=5643F interface.
- [Right-angled proof index](right_angled_proof_index.md):
  add the new current theorem, retaining older entries as predecessors.
- Catalogue entries five_constraints at line1529 and ideal6 at line1664:
  update scopes, result sources and verification ledger only after the new
  reproducibility job passes.
- [Assumptions map](manuscript_assumptions.md), lines177-205, and
  [ManuscriptRightAngled](../formal/CoxeterProofs/ManuscriptRightAngled.lean),
  lines76-106: retain old theorems; add new conditional arithmetic interfaces
  later, with the Lorentz capacity and small-facet geometry still explicit
  obligations until actually formalized. Do not transfer old4462 Lean
coverage or build status to4980. A subsequent separate module now verifies
the new arithmetic; see the [actual execution and assumption map](right_angled_cusp_capacity_verification_20261019.md).
That batch does not by itself fill the TeX translation gap above.

Only after the sources are stable: run one small arithmetic batch, then
one TeX manuscript per guarded batch, delete intermediates, visually inspect
the changed pages, and update the per-paper hashes/status. This plan is not
a request to start those jobs now.

## The next dimension-research target

Further constant substitution is unlikely to close a dimension. The new
capacity proof exposes a sharper exact quantity worth targeting: for fixed
F, let N_k count outside vertices with exactly k neighbors in F. Capacity
gives k in{0,1,2}, hence

    2V(P)-(2+2^(n-2))V(F)=2N_0+N_1.

Thus any genuine strengthening must force outside vertices with zero or one
neighbor, or constrain how the two-neighbor vertices can be arranged across
several facets. This is a concrete geometric deficit, rather than another
large generic facet census. For the ideal4 carrier with an at-least-eleven
vertex facet, a uniform positive lower bound on that deficit would improve
the first propagation step.

In parallel, the ideal-six proof still relaxes actual good edges to all32
formal sign neighbors. An additional actual-wall incompatibility or a
uniform loss in that relaxation could improve136/99. A scalar improvement
alone would still need an opposing global obstruction to exclude dimension6.

For mixed dimensions11/12, integrate the new lower ratios and ideal-face
fractions with the already existing opposition/carrier transport constraints.
No contradiction follows from4/43 or5/71 in isolation. The published ceilings
remain6 for entirely ideal and12 for mixed right-angled polytopes.

