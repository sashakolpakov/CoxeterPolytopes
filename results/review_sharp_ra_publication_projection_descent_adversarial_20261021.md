# Sharp right-angled publication: fresh adversarial mathematical read

21 October2026. Read all382 source lines of
`manuscripts/right_angled_sharp_dimension_bounds_publication.tex`,
SHA256 `0c29a7869fa46e4e1a542d8ae6f4dfc7e616136c7a1aad3984d941a95cf89722`.
**No internal proof gap or counterexample was found. No source edit is
needed for completion at the explicit published-input boundary.**
I previously reviewed this paper family and contributed editorial
clarifications. This verdict comes from rederiving the displayed argument,
not from the prior PASS records. No scientific, Lean or TeX job ran;
no manuscript, PDF or frozen verification record was changed.

1. **Actual walls and projection, lines42–45,93–149.** The assumption of
   finitely many actual facets gives a genuine attained finite maximum
   at every finite relative-interior point. With the explicitly cited
   all-wall acuteness premise, `q+t*n` has square `-1-t^2`, lies on the
   selected plane and is future. Its product with EVERY other actual
   wall is strictly negative: the original q-term is strict and the
   added term is nonpositive. Nonadjacent/ultraparallel walls therefore
   do not produce a counterexample. The result is in the relative
   interior of the actual facet, not just its supporting hyperplane.
   Every intrinsic facet comes from a ridge neighbor; orthogonality
   preserves that neighbor's unit norm and outward sign. Thus the new
   wall list is complete, finite and correctly normalized.

2. **Faces, mixed vertices and the tail, lines99–105,289–290,333–335.**
   The restriction to face dimensions>=3 is sufficient. An intrinsic
   ridge has positive dimension and a finite relative-interior point.
   At any finite stratum, incident acute normals are independent:
   splitting a linear dependence into disjoint positive/negative parts
   gives a common vector of nonpositive Euclidean square; a positive
   dependence at zero contradicts a strict interior point. The local
   face is therefore simple, and any pair of its active walls gives
   an actual ridge and hence a right angle. This justifies the short
   inheritance paragraph even when other vertices are ideal/nonsimple.
   In the finite-volume Klein closure actual faces have original
   vertices; finite triangulation into finite/ideal simplices gives
   finite volume in the dimensions used. Entire ideality is inherited.
   Actual face flags give nine- and five-faces in the respective higher
   dimensional tails. No general Coxeter-face inheritance is assumed.

3. **Descent and the code, lines128–164,197–211,277–288.** The new
   maximum is at most `A/sqrt(1+A^2)`. Induction yields
   `A_j<=A0/sqrt(1+j*A0^2)<1/sqrt(j)` for finite A0; a very large
   starting maximum or a tied maximizing facet does not lose strictness.
   Every intrinsic normal splits with a positive time height and one
   Euclidean unit direction in the CORRECT dimension. Coincident
   directions would give positive Lorentz product, so facets are neither
   lost nor double-counted. With `t_i<1/sqrt3`, each factor
   `t_i/sqrt(1+t_i^2)` is<1/2, giving the strict1/4 code in R6.
   No centering, full-span, antipodal or combinatorial code assumption
   is used by the next lemma.

4. **Cubic bound and its endpoint, lines220–274.** Direct contractions
   give `<S(u),S(v)>=24x`, the cubic/S cross contraction3x and
   `<T(u),T(v)>=x^3-3x/8`; the rank-two kernel is `x^2-1/6`.
   The displayed constant1/16 and all three nonnegative square terms
   therefore reproduce F exactly for arbitrary unit vectors in R6.
   Its sign and diagonal give N<=27. At N27, EVERY off-diagonal term
   must be zero; strict pair products exclude the1/4 root, leaving
   only-1/2. Four such vectors have squared sum-2. Thus the necessary
   conclusion really is N<=26, not merely27. This strict endpoint
   applies to the actual radial configuration from step3.

5. **Ideal equality, lines295–336.** The displayed24-cell normals have
   Lorentz norm1. Opposite spatial pairs have heights `a-r.w` and
   `a+r.w`; bounding both by1 gives a<=1, while future unit
   timelikeness gives a>=1 and forces w=0. Thus no finite point of
   that metric cell can have every height<1. Musin is used only for
   cardinality. The equality identification is the METRIC assertion
   of Kolpakov, so no kissing-code uniqueness or merely combinatorial
   rigidity is being substituted for it.

I also reopened the primary sources for the three numerical inputs:
[Alexandrov, Section4](https://arxiv.org/html/2303.09533v1#S4) explicitly
states the27-facet lower bound for every finite-volume right-angled
six-polytope, citing Dufour; its Section6 also supplies the two cited
relations. I used that universal statement rather than silently applying
a noncompact-only definition of nu to all polytopes. The original Dufour
full text was not independently retrieved in this read.
[Musin, printed p4](https://annals.math.princeton.edu/wp-content/uploads/annals-v168-n1-p01.pdf)
gives the needed cardinality bound for angular separation at least pi/3.
[Kolpakov, Theorem4, printed p1945](https://msp.org/agt/2012/12-4/agt-v12-n4-p03-s.pdf)
gives both the ideal four-dimensional24-facet minimum and isometry to
the regular metric cell at equality. These scopes match lines63–88.
[Potyagailo--Vinberg, Section3, pp66–67](https://ems.press/content/serial-article-files/43032)
also provides actual finite-volume right-angled examples through8,
supporting the cited sharpness rather than supplying an unproved example.

The substantive external foundations remain all-wall acuteness and
the named geometric facts/theorems above. They are not proved by the
scripts or conditional Lean modules. No new finite-data enumeration
is needed for the written tensor/projection argument. The publication
proof map correctly distinguishes these assumptions, historical runs,
the later editorial source and conditional formal coverage; it makes
no fresh combined-core or end-to-end geometric Lean claim. Subject to
the stated external inputs, the text supports sharp8 and entirely
ideal4. This review adds no build, visual-audit or external-refereeing claim.
