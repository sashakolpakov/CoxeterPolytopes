# Six-paper formal assumption-boundary review

Date: 20 October 2026. Reviewer: `/root/resource_cleanup`.

**Finding: PASS for disclosure and separation of assumption classes.** No new
internal theorem was found disguised as a published fact or as a consequence
of numerical replay in the six currently selected publication scopes. This
is a source-only audit of interfaces and verification claims, not external
refereeing, a new mathematical audit of every supporting proof, or an
end-to-end geometric Lean verification. The reviewer has contributed to
earlier project proofs and source exposition. No proof/checker, Lean or TeX
job was run for this audit. Only this review was written.

## Source binding and selected scopes

The selection and core/full descriptions are those of
`manuscripts/catalogue.json`, SHA256
`a60424402416fc42b64dca0c12d0bb87da954fdd96cd7d39332efd3b6ca732ad`.
The assumptions register is `docs/manuscript_assumptions.md`, SHA256
`58f4276eea135c909d125e97fcf880f4b9e60761a669dbf8e54b543c09b84a83`.

All manuscript paths below are relative to `manuscripts/`.

| Selected ID | Source | SHA256 |
|---|---|---|
| finite57 | finite_volume57.tex | c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0 |
| compact20 | compact_arithmetic_dimension_bound.tex | 1f09df2976085fd1ad1aa56ed35f1ab182ff95a8332dcfa0b0b80d4849b82203 |
| compact_farthest_bounds | compact_farthest_vertex_bounds.tex | 2718192d0d34e60418f0d0c7e00c3ed211ed0cdc0613cf5978771963a4a19705 |
| right_angled_sharp | right_angled_sharp_dimension_bounds.tex | 3c3ddd8540593004bcd5fb57832afc77ffc639268e66cb0ab29302b286654310 |
| ideal_ra4_volume_and_structure | ideal_ra4_volume_and_structure_publication.tex | e0c2ca24cec24ffcd736e201d3c7e06d75b28ef7be8093decd321ae3ac54bda1 |
| compact_ra4 | compact_ra4_volume_constraints.tex | afe807672f3a6bbf4bb2dd95d2861b60391d7c2d7c879c6d6382e749e7ad976f |

The actual selected Lean roots are the following `CoxeterProofs` modules
(shared roots appear only once in the binding):

* finite57: GuardedRank58, Concentration58, ContextCover59,
  ShapeCurvatureArithmetic, ManuscriptDataBounds, ManuscriptFinite57,
  ConcentrationDescent58.
* compact20: CompactHVector, CompactFieldBudget, DimensionComposition,
  CompactArithmetic20, ManuscriptDataBounds, ManuscriptCompact20.
* compact_farthest_bounds: CompactFarthestComposition23.
* right_angled_sharp: RightAngledSharpDimensions,
  RightAngledCubicScalarKernel.
* ideal_ra4_volume_and_structure: IdealRA4Gap, ManuscriptIdeal4.
* compact_ra4: no selected Lean root or mathematical replay checker.

The 17-root source-set binding is
`0e039ef0bd096acda4f050f1535cf1c6232ebec46186f0a595664255bff43934`.
Its exact definition is SHA256 of the concatenation of lines
`<file SHA256><two spaces><repository-relative path>\n`, in lexicographic
module-name order; module `CoxeterProofs.X` has path
`formal/CoxeterProofs/X.lean`. This binds source bytes, not imported objects
or a fresh compilation.

## Highest-priority boundary: finite57

`ManuscriptFinite57.lean` explicitly retains `WrittenGeometry` (lines38–48),
`Earlier58Tail` (lines71–74), and `Rank58Bridge.chargesAt58` (lines77–80).
These are substantial **internal mathematical inputs**, not published facts
or consequences of the selected finite-data checks:

* `finiteDomination` and `idealSurplusDomination` connect actual charges to
  checked offers and terminal profiles.
* `shapeLower` supplies the actual shape/index lower bound with the strict
  4205 correction; `paymentUpper` supplies the simultaneous actual payment
  inequality. `PublishedIndex.exactDeficit` alone does not identify the
  actual geometric index.
* `chargesAt58` realizes those data for an actual rank58 chamber.
* `Earlier58Tail.atMost58` imports the earlier internal higher-dimensional
  theorem. It is not the final57 conclusion and is not relabeled published.

The core `fresh_d13_context_dimension57_20261018.py` checks exact ideal
polynomial/margin/transfer arithmetic and finite interface composition.
Its line75 tests saved universal-domination flags; line134 reads the saved
tail bound. These assertions do not themselves prove the universal
geometric interpretation. Its explicit `verification_scope` says so.
The full tier adds inherited exact-data/domain and tail replay, without
turning the above geometric interfaces into Lean theorems.
`ConcentrationDescent58` additionally proves the discrete descent, terminal
recognition and lifting logic; the real normalized comparison and actual
cusp interpretation remain inputs. The assumptions register correctly
distinguishes all these boundaries.

## SharpRA: stronger coordinate route, remaining actual geometry

`RightAngledCubicScalarKernel.lean` does more than accept the old
`StrictCodeSixConsequences`: its preferred `CoordinateWitness` supplies
unit six-coordinate vectors, the strict product gap, 36/216 coordinate
arrays and two pointwise contraction identities. Finite sums of coordinate
squares derive Gram positivity; the cubic SOS, count bound and strict26
endpoint are then proved. `cubic_geometry_from_coordinate_bridge`
constructs the older `WrittenCubicGeometry` interface. Thus the new
spherical-code conclusion is not merely a premise of this preferred route.

The realization of actual radial-code vectors and tensors remains an
input. `ProjectionGeometry.farthestProjection` in
`RightAngledSharpDimensions.lean` also supplies the actual facet, projected
interior point and height descent; the three projection steps are composed,
not geometrically constructed in Lean. Face inheritance and the stated
Dufour, Musin and Kolpakov facts remain typed inputs. Those published facts
are separated from the internal projection and center arguments.
The two selected Python checks establish fixed24-cell model facts and
universal polynomial/tensor identities, not universal hyperbolic geometry.

## Other four scopes

* **compact20:** `EarlierInternalProofs.higherBranches` and
  `nongoldenAt21` are expressly inherited internal proofs. Whole-star
  coverage, D5/H4 consumer geometry, actual h-vector/charge interpretation
  and the golden rank21 realization remain internal interfaces. The
  published29/13-field inputs are separate. The three core checks concern
  exact certificates and composition; parked19 refinements are not selected
  Lean imports or mathematical inputs to this paper.
* **compact23:** the real-accounting bridge derives the scaled Nat budget
  from local scalar inequalities. Actual energy decomposition, component
  bounds, caps/debit extraction and the coupled D23/I obstruction remain
  internal inputs. Source-integrity checking does not prove these theorems.
* **ideal4:** the finite-model checks do not identify every actual equality
  configuration. Equality-spectrum reduction, actual cuts, zero-case
  recognition, `OctDGeometry.pathRecognition` and
  `canonicalTwoCellGluing` are explicitly internal geometric inputs. The
  recursive accounting is proved conditional on them; chain recognition
  and canonical geometric placement are not supplied by a fixed model.
* **compactRA4:** the selected scope is written geometric proof and PDF
  reproduction, with no corresponding geometric Lean or universal numerical
  checker. It does not claim a proved120-cell minimum from the supporting
  conditional bounds.

## Predicate-source binding and conclusion

The directly inspected script predicates, besides catalogue selection, were:

* `scripts/fresh_d13_context_dimension57_20261018.py`:
  `9a41aa1d0395087a41692fb1d3a0d2b7f93f667446f4724d5b7920449b6d1795`;
* `scripts/check_right_angled_sharp_dimension_sos_20261020.py`:
  `bd0ed1c4c1b7e18a219a2a8aca2b178dbe1d559d7e9ec90ad609150aedf8b2c3`;
* `scripts/check_ideal_right_angled_dimension4_model_20261020.py`:
  `bb506ec1ee98dee1429d6e668553af40316ebdeda6b59b82bbf03b8f3a6dc261`.

Successful standard-axiom printing does not inventory proof-parameter
assumptions. Exact script outputs certify their declared computations,
not universal geometric interpretation. The current selected scopes and
assumptions register respect that distinction. There is no identified
misclassification requiring a documentation patch or another kernel merely
to increase coverage counts. These remain verified conditional logic plus
declared exact-data checks, with the written geometric bridges explicit.
