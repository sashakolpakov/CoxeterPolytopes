# SharpRA final core steps and input scope

21 October 2026. Read-only manual finishing check of
`manuscripts/right_angled_sharp_dimension_bounds_publication.tex`, SHA-256
`0c29a7869fa46e4e1a542d8ae6f4dfc7e616136c7a1aad3984d941a95cf89722`,
and its current publication proof map, SHA-256
`ff9087dfc4851c6f86e7891452df50849e634296f67cff004b619da09fd12cc3`.
No build, scientific worker, Lean job, new research or source edit.

**Finding: no remaining completion blocker found at the stated geometric
and published-input boundaries. No further editorial change requested.**
The current five-page text supplies the core algebra and the strict
endpoint arguments needed for sharp8 and entirely-ideal4.

|Current TeX lines|Manual check and result|
|---|---|
|93–105|The all-wall implication is named and cited rather than silently substituted for actual right ridge angles. The new face-inheritance explanation is correct in dimensions at least three: a ridge has finite relative-interior points; incident walls restrict orthogonally, and the Klein face has finite/ideal vertices and a finite simplex triangulation. Entire ideality and the higher-dimensional face reductions are covered.|
|118–150|For `v=q+tn`, `v²=q.v=-1-t²`, `v.n=0`, giving the correct future unit normalization. Every other wall has strictly negative product because its first numerator term is strict and the second nonpositive. All intrinsic facet walls come from orthogonal ridge neighbors and remain unit. This proves actual relative-interior feasibility, including nonneighbor walls.|
|152–164|The maximum exists in the finite wall family and is now explicitly denoted `A_j`. Monotonicity gives the local recurrence and induction gives `A_j<=A_0/sqrt(1+jA_0²)<1/sqrt(j)`. Ties do not remove strictness.|
|197–211,277–290|Unit Lorentz splitting gives the displayed radial bound. Identical radial directions would give positive Lorentz product, so facets are counted injectively. Three steps produce strict products below1/4, on an actual finite-volume right-angled six-face. Higher dimensions reduce to an actual nine-face.|
|219–274|The cubic factorization and contractions check: cubic/S cross products are `3x`, `S.S=24x`, hence `T.T=x³-3x/8`. Their positive square sum gives `N²/16<=sum F<=27N/16`. At27, strict pair products below1/4 force every off-diagonal product to be-1/2; four vectors then have squared sum-2. The strict cap26 is proved, not supplied.|
|295–317|The displayed24 normal vectors have unit Lorentz norm and opposite spatial pairs. Two height bounds imply `a<=1`; the future unit equation implies `a>=1`, and equality forces the center. This is the required metric calculation.|
|320–347|The projected four-dimensional facet inherits all three geometric hypotheses and has strict radial products below1/2. Musin bounds cardinality; Kolpakov supplies facet equality and hyperbolic isometry. No kissing-configuration uniqueness is used. Five-face inheritance closes the higher-dimensional ideal tail.|

The external statements are used with the following precise scopes:

* **Facet27, TeX63–72:** every finite-volume right-angled six-polytope,
  including mixed vertices; no arithmeticity, classification or equality
  assumption. The retained direct-primary
  [facet/acuteness audit](review_sharp_ra_published_inputs_dufour_20261020.md)
  checked Alexandrov Section4's explicit consequence and Section6's
  `nu5>=26`, `facets6>=1+nu5`, with `nu=facets+cusps`.
  Alexandrov's noncompact definition includes mixed vertices; the compact
  six-case is already excluded by the cited compact-four bound. Dufour's
  original full text was not obtained. The TeX accurately says **via
  Alexandrov**, and the current map discloses the access limit.
* **Musin, TeX74–80:** finite subsets of `S³` with pair products at most
  one half have at most24 members. The retained
  [ideal-input audit](review_sharp_ra_published_inputs_ideal_branch_20261020.md)
  checked the original definitions and Main theorem, printed pp2/4.
  The paper's strict code satisfies this closed bound.
* **Kolpakov, TeX82–88:** Theorem4, printed p1945, gives at least24
  facets for entirely ideal right-angled four-polytopes and **isometry**
  to the regular ideal24-cell at equality. The retained original-source
  audit distinguishes this from the separate volume theorem. It also
  checks the actual sharp examples; no existence claim is inferred from
  an upper bound.

**Local-source limitation.** A repository file search, including ignored
PDF/HTML paths, found no local original of Dufour, Alexandrov, Musin,
Kolpakov or Vinberg. The existing primary-source audit reports expressly
record that no PDF assets were saved. Therefore this finishing pass
checks the current statements against those preserved direct-primary
reports; it does **not** claim a fresh comparison with locally stored
original literature or a new verification of the published proofs.
No new browsing was needed to resolve a discrepancy.

The current map, lines80–106 and194–205, correctly separates exact
algebra/model checks and conditional Lean logic from actual geometry and
published theorems. Historical four-/six-job core records retain their
original bindings; the current audit edition has its separate bounded
PDF build. It claims neither a new combined successor core run nor
end-to-end geometric Lean verification. No misleading verification
language requiring correction was found.

Arithmetic research remains held. The unreviewed whole-A19 candidate
`d15ef2ac` was not mutated, reviewed to completion or promoted by this
finishing check.
