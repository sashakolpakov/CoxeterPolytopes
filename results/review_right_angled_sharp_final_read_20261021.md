# Final read: sharp right-angled dimension bounds

21 October 2026. Read the complete five-page source
`manuscripts/right_angled_sharp_dimension_bounds.tex`, SHA-256
`3c3ddd8540593004bcd5fb57832afc77ffc639268e66cb0ab29302b286654310`
(378 lines), against the current proof/assumption maps and existing
foundations, primary-source and execution records. This is a fresh
source read by an earlier foundations/primary-input reviewer, not an
additional independent outside-author review. No worker, Lean, TeX or
literature search ran; no manuscript or PDF changed.

**Finding:** the text supports the advertised universal finite-volume
right-angled ceiling eight and entirely ideal ceiling four, including
sharpness, at its explicitly named published/geometric boundaries.
No mathematical or essential exposition gap was found. One short
terminology correction is proposed below.

The two proofs have a clear common spine. All-wall nonpositivity makes
`q'=(q+t n)/sqrt(1+t^2)` strict against every other actual wall
(lines 115–147), so the construction reaches the relative interior of
an actual facet. Its intrinsic walls are precisely the orthogonal
ridge neighbors; the maximum-height recurrence therefore controls
every facet of the next face. Finitely many walls at a finite point
give the strict bound `A_j<1/sqrt(j)` (149–160). This explains both
the geometric legitimacy and the strictness needed at the endpoints.

For the general theorem, three descents from dimension nine produce
an actual finite-volume right-angled six-face. Its distinct radial
wall directions have products strictly below `1/4` (193–207,
273–287). The displayed tensor identity gives `N<=27`; equality would
force every off-diagonal product to be `-1/2`, and four vectors then
have squared sum `-2`. Thus `N<=26` (211–270), contradicting the
published universal six-dimensional minimum 27. The strict endpoint
is proved, rather than inferred from a weak spherical bound.

For the ideal theorem, one descent from dimension five gives an
entirely ideal four-dimensional facet with heights strictly below
one. Musin gives at most 24 facets; Kolpakov gives at least 24 and
metric identification with the regular ideal 24-cell. Its opposite
unit-normal pairs force any point with all heights at most one to
be the center with every height exactly one (291–333). This correctly
uses metric rigidity, with no uniqueness assertion about kissing
configurations. Face inheritance supplies both higher-dimensional
reductions. Published examples supply attainment in dimensions eight
and four (57–59).

|Boundary|Exact use and verification limit|
|---|---|
|Acute-polytope geometry, lines 93–102|Vinberg's all-wall implication, complete intrinsic ridge-neighbor normals, finite-volume/right-angle inheritance in dimensions at least three, and entire-ideality inheritance. The written projection checks all walls, including nonneighbors; nonsimple ideal vertices are allowed. No descent reaches dimension two.|
|27-facet input, lines 63–72|Every finite-volume right-angled six-polytope, including mixed finite/ideal vertices. The existing [facet/acuteness source audit](review_sharp_ra_published_inputs_dufour_20261020.md) directly checked Alexandrov's explicit consequence and scope. Dufour's original full text was not obtained; the manuscript accurately cites the consequence via Alexandrov. This limitation is bibliographic verification, not a remaining internal inequality.|
|Ideal inputs, lines 74–89|Musin's arbitrary `S^3` code cardinality bound and Kolpakov Theorem 4's facet minimum plus hyperbolic-isometry equality. Their originals and the examples were checked in the [ideal-input audit](review_sharp_ra_published_inputs_ideal_branch_20261020.md).|
|Finite data|No unpublished facet census, code enumeration, optimizer or forest table is needed. The SOS worker checks universal exact polynomial identities; the 24-cell worker checks the displayed finite model and witnesses. Neither supplies an arbitrary-polytope classification. All calculations needed by the mathematical proof appear in the text.|

The recorded [six-job combined core](right_angled_sharp_projection_height_core_verification_20261020.json)
passes two exact workers, three conditional Lean roots and a temporary
five-page PDF build. Its [30 actual printed declarations](right_angled_sharp_projection_height_core_axiom_validation_20261020.json)
are the original dimension interface's 12, cubic kernel's 10 and
projection/height kernel's eight; these repeated prints do not enlarge
the separate publication-kernel count. This read reran none of them.
The kernels derive the scalar/coordinate algebra and count implications;
actual real wall and face realization, radial/tensor realization and
the pointwise contraction bridge, positive normalizers, geometric
interpretation and published theorems remain explicit application
inputs. The older abstract geometry records are not thereby
instantiated. The package supports a complete written proof using
named published inputs, not an end-to-end geometric Lean theorem or
an external-referee verdict.

One concrete readability edit is advisable: at source line 319 replace
“four-facet `Q`” by “four-dimensional facet `Q`.” The intended dimension
is clear from the surrounding proof, but the current phrase can be
misread as a facet with four facets. This is terminological and changes
no argument. Leave the source unchanged until the root coordinates the
edit and any source/PDF receipt update. No other expansion is needed
for the five-page paper's declared scope.
