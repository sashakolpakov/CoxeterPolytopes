# Fresh mathematical audit of the sharp right-angled dimension paper

24 October 2026. Reviewer: `/root/audit_sharp_ra`.

**Outcome: the complete written arguments for ceilings eight and four
survive this fresh source audit at their stated external-input boundary.**
I found no internal geometric, algebraic, endpoint or coverage defect in
the 382-line publication source. One reference check remains incomplete:
Dufour's original article was not obtained in full. The universal
six-dimensional 27-facet statement is therefore retained as a precisely
named external premise, rather than reported as directly reverified in
that original article.

This is a manual mathematical and source-code audit, not a new runtime
check, publication, external referee report or full geometric Lean proof.
I supplied no derivation or edit to these frozen sources before this
audit. The manuscript, PDF, checkers, certificates and Lean modules were
left unchanged. No worker was imported; no mathematical enumeration,
optimizer, test, Lean or TeX process ran. Hashing and reading files were
the only local programs used. Earlier execution records are inherited
evidence with their original source bindings.

## 1. Exact sources and scope

I read the ENTIRE publication TeX, the complete 178-line descending-height
source, both complete exact checker sources, and the three selected Lean
sources. The Lean review checks the stated interfaces and their match to
the written proof; it does not replace the previously recorded kernel
compilations or actual axiom audits. The PDF hash is bound below, but this
audit does not claim a new visual inspection of every PDF page.

|Source|Lines|SHA256|
|---|---:|---|
|`manuscripts/right_angled_sharp_dimension_bounds_publication.tex`|382|`0c29a7869fa46e4e1a542d8ae6f4dfc7e616136c7a1aad3984d941a95cf89722`|
|`manuscripts/build/right_angled_sharp_dimension_bounds_publication.pdf`|—|`609a9fc025b43ae5c899b1b37b1051e592de6a8142f3ffd275fcd9dca360d443`|
|`results/right_angled_dimension8_descending_height_proof_20261020.md`|178|`9470fa1fa3a973a12c412fbc1257537655536dba4ffc4e9ac23c320615babdbe`|
|`scripts/check_right_angled_sharp_dimension_sos_20261020.py`|154|`bd0ed1c4c1b7e18a219a2a8aca2b178dbe1d559d7e9ec90ad609150aedf8b2c3`|
|`scripts/check_ideal_right_angled_dimension4_model_20261020.py`|96|`bb506ec1ee98dee1429d6e668553af40316ebdeda6b59b82bbf03b8f3a6dc261`|
|`formal/CoxeterProofs/RightAngledSharpDimensions.lean`|261|`e1844d73a3f07de6cc512e514b948dafef4253e272b3f98de36148543554e323`|
|`formal/CoxeterProofs/RightAngledCubicScalarKernel.lean`|510|`f55568035b8852a6284f1bd404c897758f30b301b1d85f38f09e4d9de365e6ee`|
|`formal/CoxeterProofs/RightAngledProjectionHeightKernel.lean`|282|`2fac8028e432a2104a58e423c65fc38c9e214fed5f52210e408aabe3f15f0e90`|

The present proof map is
`docs/right_angled_sharp_dimensions_publication_proof_map_20261021.md`,
SHA256 `ff9087dfc4851c6f86e7891452df50849e634296f67cff004b619da09fd12cc3`.
Its earlier PASS labels were not used to infer mathematical correctness.

## 2. Geometry before the projection

The manuscript assumes a finite-sided, finite-volume convex hyperbolic
polytope, with right angles at every actual ridge. It does not assume
arithmeticity, compactness, simple ideal vertices or a restricted facet
spectrum.

The all-wall sign is essential. Orthogonality of ridge neighbors alone
would not justify the inequalities against unknown or nonadjacent walls.
The original English [Vinberg survey, printed p.42, Section 2.1](https://www.mathnet.ru/php/getFT.phtml?jrnid=rm&option_lang=eng&paperid=2140&what=fullteng)
was retrieved directly. It states the nonpositive off-diagonal Gram
condition and records Andreev's implication from nonobtuse actual ridge
angles. This supplies exactly the global sign used in the paper.

The face-inheritance passage also handles the relevant noncompact cases.
At a finite relative-interior point of an intrinsic ridge of a face, all
original supporting walls through that point are pairwise orthogonal.
One way to justify that statement is the same nonadjacent-plane theorem:
intersecting original wall planes must be ridge neighbors. Independent
orthogonal normals give a local orthant. Restricting the normals to the
supporting plane of the face consequently preserves the right angles.
The argument does not assert simplicity at an ideal vertex.

In the closed Klein model a finite-volume polytope is the convex hull of
its finitely many finite or ideal vertices. An exposed face has vertices
among those vertices. Thus an entirely ideal polytope cannot acquire a
finite vertex in an actual face. A finite triangulation of a face uses
finite or ideal simplices and has finite hyperbolic volume in the
dimensions used here. Every selected face has dimension at least four;
the infinite length of an ideal edge and the zero angles at ideal
polygon vertices are irrelevant. No arbitrary Coxeter-face inheritance
is used: right angles are the special hypothesis making this valid.

## 3. The projection and complete height chain

Let (q^2=-1), (n^2=1), (q\cdot n=-t<0). Direct expansion gives

\[
(q+tn)^2=-1-t^2,\qquad (q+tn)\cdot n=0,
\qquad q\cdot(q+tn)=-1-t^2<0.
\]

The normalized vector is therefore a finite future unit point on the
selected supporting plane. For EVERY other original wall (m),

\[
(q+tn)\cdot m=q\cdot m+t(n\cdot m)<0.
\]

The first term is strict, the second nonpositive. The point lies in the
relative interior of the ACTUAL facet, with no untested wall. There is
no Gram-realization converse or formal candidate-face assumption.

Each intrinsic facet of that facet is its ridge with one actual
neighboring wall. The unit normal of that neighboring wall is already
orthogonal to (n), so restriction preserves its norm. The new height
is exactly the old height divided by (\sqrt{1+t^2}). Choosing (t)
maximal among the finitely many actual walls gives

\[
A_{j+1}\le {A_j\over\sqrt{1+A_j^2}},\qquad
A_j\le {A_0\over\sqrt{1+jA_0^2}}<{1\over\sqrt j}.
\]

The strict final inequality uses finite (A_0), rather than genericity,
uniqueness of the maximum or a limiting argument. Repeating intrinsically
does not require the selected points to project to vertices or the same
facet to maximize at successive stages.

For dimension nine, the complete chain is (9\to8\to7\to6). It produces
an actual six-face with every intrinsic unit-normal height strictly
below (1/\sqrt3). For the entirely ideal case the chain is (5\to4),
and all new heights are strictly below one. These two chains use the
same projection theorem, with no dependence between their final
dimension conclusions.

## 4. Facet counts really are spherical-code cardinalities

At the final point write each unit wall normal as

\[
m_i=(t_i,\sqrt{1+t_i^2}\,u_i),\qquad |u_i|=1.
\]

The all-wall sign gives

\[
u_i\cdot u_j\le
 {t_it_j\over\sqrt{(1+t_i^2)(1+t_j^2)}}.
\]

The function (s/\sqrt{1+s^2}) is increasing for (s>0).
Consequently heights (<1/\sqrt3) give products (<1/4), and heights

\(<1\) give products (<1/2\). If two radial directions coincided,
their Lorentz product would be

\[
\sqrt{(1+t_i^2)(1+t_j^2)}-t_it_j>0,
\]

because the square of the first positive term exceeds the square of
the second by (1+t_i^2+t_j^2). This contradicts acuteness. There is
exactly one distinct code point for each actual facet, with no lost
multiplicity or change in ambient Euclidean dimension.

## 5. The strict six-dimensional tensor bound

I independently expanded the contractions, including their dimensions.
For (x=u\cdot v) and (C(u)=u^{\otimes3}),

\[
\begin{aligned}
\langle u\otimes u,v\otimes v\rangle&=x^2,&
\langle u\otimes u,I\rangle&=1,& \langle I,I\rangle&=6,\\
\langle C(u),S(v)\rangle&=3x,&
\langle S(u),C(v)\rangle&=3x,&
\langle S(u),S(v)\rangle&=24x.
\end{aligned}
\]

The last equality has three matching-slot contributions of (6x)
and six cross-slot contributions of (x). Thus the subtractions

\(A=u\otimes u-I/6\) and (T=C-S/8\) give

\[
\langle A(u),A(v)\rangle=x^2-1/6,
\qquad \langle T(u),T(v)\rangle=x^3-3x/8.
\]

Substitution gives precisely the paper's cubic and sum of squares:

\[
F(x)=(x+1/2)^2(x-1/4)=x^3+3x^2/4-1/16,
\]

\[
\sum_{i,j}F(u_i\cdot u_j)
 ={N^2\over16}+{3\over8}\Big|\sum_i u_i\Big|^2
  +{3\over4}\Big|\sum_i A(u_i)\Big|^2
  +\Big|\sum_i T(u_i)\Big|^2.
\]

These are ordinary finite Euclidean tensor norms, so all additional
terms are nonnegative. The diagonal value is (27/16), and every
off-diagonal value is nonpositive when its product is (<1/4). Hence

\(N^2\le27N\). The empty code causes no division issue; for a nonempty
code (N\le27).

At (N=27), equality forces every nonpositive off-diagonal summand to
vanish. The strict upper product excludes the root (1/4); every
distinct product must equal (-1/2). Four selected vectors would have
squared sum (4+12(-1/2)=-2). This proves (N\le26). The strictness
cannot be relaxed to (\le1/4), and the manuscript never does so.

The actual six-face is therefore bounded above by 26 facets. The
external 27-facet premise supplies the contradiction. A higher
dimensional polytope has an actual face flag down to dimension nine,
and the established inheritance applies. There is no missing
higher-dimensional or mixed-vertex branch.

## 6. The entirely ideal endpoint

One projection produces an actual entirely ideal right-angled
four-dimensional facet, with a strict (1/2) code in (S^3).
[Musin's original journal PDF, printed p.4, Main theorem](https://annals.math.princeton.edu/wp-content/uploads/annals-v168-n1-p01.pdf)
was retrieved directly. Its kissing-number bound supplies cardinality
at most 24 for the weaker closed product condition, so it applies here.
No uniqueness of a kissing configuration is used.

[Kolpakov's original journal PDF, printed p.1945, Theorem 4](https://msp.org/agt/2012/12-4/agt-v12-n4-p03-p.pdf)
was also retrieved directly. It gives both the facet minimum 24 and
hyperbolic isometry to the regular ideal right-angled 24-cell at
equality. This is the needed FACET theorem, rather than an inference
from the separate minimum-volume theorem. Its metric conclusion is
stronger than mere combinatorial recognition and is exactly what the
remaining calculation needs.

The displayed normals ((1,r)), (r=\pm e_i\pm e_j), have Lorentz
square one. In Klein coordinates their inequalities are

\(|y_i|+|y_j|\le1\), giving the usual ideal 24-cell. The eight axial
vertices and sixteen half-sign vertices lie on the ideal sphere. For
each spatial pair (r,-r), a future unit point ((a,w)) has heights

\(a-r\cdot w\) and (a+r\cdot w\). If both are at most one, their
sum gives (a\le1); timelikeness gives (a\ge1). Thus (a=1,w=0),
and all heights equal one. The projection's strict heights are
impossible. The normals are opposite SPATIALLY, not negatives as
Lorentz vectors, and their unit scaling agrees with the projection.
Isometry preserves that scaling and the heights.

Actual five-face inheritance excludes higher entirely ideal dimensions.
For attainment, [Potyagailo–Vinberg, Section 3 and Proposition 3.1,
printed pp.66–67](https://ems.press/content/serial-article-files/43032)
was read directly for the dimension-eight construction. The regular
ideal 24-cell attains dimension four. No existence or realization claim
is extracted from an unclassified list of candidate Gram diagrams.

## 7. The external premise still requiring original-source access

The required statement is exactly:

> Every six-dimensional finite-volume right-angled hyperbolic polytope,
> allowing mixed finite/ideal vertices, has at least 27 facets.

The [Dufour publisher record](https://link.springer.com/article/10.1007/s10711-009-9454-2)
was retrieved and confirms the article and finite-volume subject, but
offers subscription access; the original PDF was not accessible. No
readable original copy was located in this pass or in the local file
inventory.

[Alexandrov's Section 4](https://arxiv.org/html/2303.09533v1#S4) explicitly
prints this universal statement and attributes it to Dufour Lemma 1
and Proposition 4; Section 6 prints its numerical consequence again.
This is a useful retrieval pointer and the manuscript discloses that
route, but it is not substituted for a direct rereading of Dufour's
original result in this audit. The direct-source check of that article's
internal numbering, hypotheses and proof remains OPEN. No hidden
ideal-only condition is being silently inserted into the theorem:
the external premise is recorded with the exact universal scope above.

The other named numerical inputs and the all-wall passage have the
direct original-source checks recorded in Sections 2 and 6. Published
proofs themselves remain external theorems, not newly re-refereed
inside this report.

## 8. Exact checker correctness and termination

The SOS worker is a finite polynomial computation over `Fraction`.
Its monomial keys are 12-component exponent tuples for six (u) and
six (v) coordinates. Addition, scaling and multiplication are the
usual coefficient operations; zero coefficients are removed. The
tensor loops cover all (6^2=36) and (6^3=216) coordinates, with no
symmetry reduction that could miss a term. Every loop is bounded by
these fixed sizes or by finite polynomial dictionaries. There is no
candidate search, completeness heuristic or floating-point tolerance.

The worker checks the contractions as polynomial identities before
imposing unit norms. It explicitly writes the remaining error as a
multiple of (u^2-1) plus a multiple of (v^2-1). Consequently the
unit-sphere specialization is universal, not tested on samples. I
checked the two norm factorizations, the scalar factorization and the
remainder formula against the written tensor expansions. They match.
The worker also checks coefficient positivity, the ratio 27 and the
four-vector endpoint value (-2).

`--verify` recomputes this finite result and compares it with the saved
certificate, including source bindings checked before and after the
computation. That is an appropriate replay specification. The saved
certificate alone is not evidence that the geometry or a code census
was verified. This audit did not execute the worker again.

The 24-cell worker generates all (6\cdot4=24) displayed normals,
checks all 276 unordered pairs, pairs every spatial opposite and
checks the 24 displayed null vertices and their six-wall cubic links.
It provides a strict finite witness for every wall to be an actual
facet. Its combinations and sign-product loops are finite and
complete for THESE DISPLAYED model data. It correctly declines to
claim complete vertex enumeration or arbitrary-polytope classification.
Neither such stronger claim is needed for the paper's opposite-pair
height argument. This worker too was read in full and was not rerun.

## 9. Conditional Lean boundary and inherited executions

The selected Lean sources retain honest layered interfaces:

* `RightAngledSharpDimensions` composes actual projection witnesses,
  the independent published inputs and typed code consequences. Only
  actual counts use `Nat`; real heights are abstract predicates.
* `RightAngledCubicScalarKernel` derives the cubic count and endpoint
  consequences from coordinate Gram and tensor data over an arbitrary
  ordered ring. Its preferred `CoordinateWitness` still supplies the
  actual unit coordinate vectors and TWO tensor contraction identities.
  Coordinate sum-of-squares arguments derive the positivity facts;
  the desired cardinality bound is not a field of that witness.
* `RightAngledProjectionHeightKernel` derives the projection algebra,
  all-wall strict signs, orthogonal heights, cleared recurrence and
  radial strict gap. It assumes a positive normalizer and primitive
  coordinate facts. Actual face interpretation, finite maxima,
  completeness of intrinsic ridge neighbors and inheritance remain
  geometric obligations.

These coordinate layers are not instantiated into a concrete model of
all actual hyperbolic polytopes. The ideal metric interpretation and
the three named published inputs remain external. This is a genuine
conditional logic boundary, not a complete formal geometric theorem.
Printed axioms alone would not list these ordinary theorem parameters.

The inherited six-job receipt
`results/right_angled_sharp_projection_height_core_verification_20261020.json`,
SHA256 `2c290d5b7bf2b94221ce66420f12b936172b89f9b3fd990338c43ddb1f925680`,
binds the source hashes above for the selected three modules and the
exact workers. It records earlier capped checks and audits; its PDF
job belongs to the predecessor manuscript edition. The CURRENT
five-page source/PDF build is separately bound by
`results/sharp_ra_audit_edition_readability_build_20261021.json`,
SHA256 `572b14642e9732653b75617c456765c53aa2e012b37921b395142d03321495b1`.
I checked those distinctions directly in the records. No combined
fresh runtime, recompiled module or additional formalization is
credited by this audit.

## 10. Manuscript and consolidation recommendation

No mathematical source repair is required by this audit. The current
five-page proof is already one focused paper with a shared projection
theorem and a clear split between the universal tensor argument and
the entirely ideal metric equality argument. Adding exploratory mixed
vertex, opposition or higher-dimensional case tables would obscure
the proof and supply no missing hypothesis.

For consolidation, retain this paper's theorem scope and displayed
calculations. Keep the original-source access limitation for the
27-facet premise visible in its proof map or external-input register.
List the still external geometric bridges next to the existing
conditional Lean results; do not rerun successful modules merely to
make a new audit date. Any later manuscript edit should have a
separate delta review and source-bound build receipt.

Root acceptance of this audit is a separate coordinating step. No
dimension table, global claim or proof-verification status was changed
by the reviewer.
