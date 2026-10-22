# Full static review: consolidated thirteen-group arithmetic16 frontier

22 October2026. **Full source/logic/assumption-map PASS; no correction
requested. UNCOMPILED, zero actual axiom-print outputs.** No Lean,
compiler, import, worker, enumeration, test, TeX or Git job was run by
this reviewer. Only source/metadata reading and archive decoding occurred.
The frozen source, map and prior gates were not edited.

ENTIRE reviewed current files:

* `formal/CoxeterProofs/ArithmeticSixteenFrontierRefinement.lean`,
  SHA256 **11dd035a5556b7a8627889b98dd9ba6f835844503302007b05a57c03c9714e33**,
  10,920 bytes/284 lines;
* `docs/arithmetic_sixteen_frontier_refinement_logic_20261022.md`,
  SHA256 **f41bf48414e500400fee889475ed637fea186e9c305eb3c848c30c4b61c6d328**,
  13,190 bytes/196 lines.

I did not author this formal module or map. I authored several external
written premises, notably whole A16, the new ordinary-dihedral/D14 consumer
and the new typed-D14 proof, and contributed inherited projection/field
foundations. This is an independent review of the NEW formal composition,
its binders, finite case logic and premise mapping, not an independent
self-review of those geometric proofs. Their separate full outside-author
reviews/root gates remain necessary external authorities.

## 1. Complete data and exact refinement

Source16–140. `Group` has the original nineteen constructors in the
fresh a8cc frontier's order. `ForestCase.group` covers every constructor,
erasing only its appropriate spectator label. Its ordinary grouped cases
are intentionally labels; the formal file does not construct diagrams
or enforce numerical ranks/angles. The map accurately requires the
application's `ActualForest` relation to carry the ENTIRE corresponding
written meaning.

The two ACTUAL typed-seven meanings remain separate in `TypedSeven`.
A14 has five separate spectator constructors, and its refined list omits
ONLY `actualIOrdinaryA1`. D14 now has the complete three-case enum:
ordinary I2(m>=4), N7 and O7. Ordinary A2 and the other decomposed rank2
spectators were already excluded in a8cc; they need not be fabricated
as new D14 cases. `ofTypedSeven` maps N7 and O7 to their exact matching
D14 constructors, without changing their actual isolated-opposite meaning.

The six whole groups removed are precisely original3,5,12,13,16,17:
ordinary A11, A11+O7+A3, whole A16, B14+either typed7, D13, D14.
The thirteen surviving original numbers are

    1,2,4,6,7,8,9,10,11,14,15,18,19.

The A14 group10 remains, with ordinary T2, two ordinary A1, N7 and O7.
It is not counted as a seventh whole-group deletion. I compared the
map's entire nineteen-row table to a8cc Sections2/3/7 and to the current
root thirteen-group table; the ordinary T3/T4 expansions, H4 omission
under typed K7, actual-I cases and nonshort-A1 meanings agree. No typed
defect is identified with a bare old7 or ordinary dihedral component.

Source143–161, `exact_enumerations`, proves constructor completeness of
the three ORIGINAL enums, then closed list length/distinctness facts:
19/13 groups,5/4 A14 spectators and3 D14 spectators. The `cases ...;
decide` obligations reduce to explicit finite lists with derived decidable
equality. No desired length, maximum or coverage assertion is supplied
as a hypothesis. `decide` is ordinary kernel-checked finite reduction,
not `native_decide` or a supplied external answer. Actual elaboration
has nevertheless not occurred and is not certified by this static review.

## 2. Exactly nine EXTERNAL written interfaces

Source172–193. The `WrittenInputs` structure contains exactly

    freshComplete;
    d13Exclusion; a11OrdinaryExclusion; b14TypedSevenExclusion;
    wholeA16Exclusion; a11OldSevenA3Exclusion;
    a14ActualIOrdinaryA1Exclusion;
    d14OrdinaryDihedralExclusion; d14TypedSevenExclusion.

That is one original classification and EIGHT exclusion interfaces.
The last two are SEPARATE ordinary and typed D14 results. There is NO
primitive whole-D14 exclusion, no `Survives` field, no supplied thirteen-
group conclusion, no desired count and no dimension<=15 hypothesis.

All eight exclusion interfaces require exact dimension16, the SAME p,
`Selected p` and the matching `ActualForest p ...`. The fresh interface
produces a description for that supplied p. It does not select a different
farthest vertex, assert the result at a replacement point, require a
particular root-length squareclass or assume the surviving list.

I freshly matched ALL NINE literal mathematical source pins in the map:
a8cc202c,6cd78a94,2e69ab9c,237080ad,7ba7e0dd,ca4760e3,60bd9822,
cc2b3869 and7331b2bb. The complete hashes and paths are those of the
immutable reviewed map, not inferred aliases. The inspected statements
support the exact interfaces: D13 covers BOTH original remainders,
ordinary A11 is a uniform component consumer, B14 excludes BOTH actual
defects, A16 is WHOLE, A3's consumer is EXACTLY A11+O7+A3, and A14's
exclusion is ONLY actual I plus ordinary A1. The ordinary-D14 theorem
does not cover O7 by fiat; the separate7331 theorem covers N7/O7.

The current root whole-D14 gate is
`results/review_arithmetic16_whole_d14_root_acceptance_and_thirteen_group_frontier_20261022.md`,
SHA256 **5931aa492bf7561f3b276f4ec3e675cbabd473ba1e151ef7c6aaff074b2e5cef**.
I read it ENTIRELY and checked its ordinary/typed conjunction, original
row exhaustiveness and thirteen-row conclusion. It binds final7331,
full outside report411b and the preserved816d one-sentence predecessor.
The combined A16/A3 gate46990774, fresh/A14 gate7b796609, D13 gate7d9e43fc,
A11 gate1f0aa0f8, corrected B14 gateed81f16c, ordinary-dihedral gatee5e01fd6
and full typed-D14 review411b7177 all freshly match the map's hashes.
Every one of the map's ten local Markdown link targets exists.

## 3. Every conditional proof, including same-point witnesses

`a14_actual_i_ordinary_a1_exclusion` (202–206) exposes only the exact
external subcase and passes through the original p/selected/forest proof.
It does not claim to derive the underlying geometric inequality.

`d14_all_spectators_exclusion` (209–219) genuinely derives the WHOLE
row exclusion. It cases on ALL THREE D14 constructors, applies the
ordinary field to its one case and the typed field separately to each
typed meaning. `ofTypedSeven` reduces definitionally to the matching
constructor, so the original `forest` term has exactly the needed type.
No nontrivial rewrite, transport across a point or supplied whole-row
contradiction is hidden in these branches.

`classified_forest_refinement` (222–257) exhausts every one of the
nineteen forest constructors. It dispatches the five earlier whole-group
exclusions to exactly their matching branches, uses the DERIVED D14
theorem on `.d14 spectator`, and cases on ALL FIVE A14 alternatives.
Only A14's actual-I/ordinary-A1 branch is contradictory. Each retained
branch is a closed explicit-list membership after the case split; its
typed parameter is erased only where intended by `.group`. `False.elim`
is used solely after the matching external contradiction has been applied.

`same_point_complete_refinement` (260–266) obtains `family, forest`
from `freshComplete at16 p selected` and returns those SAME witnesses
with the derived `Survives family`. There is no replacement of p,
change of `Selected`, independent forest witness or existence assumption
for another polytope. No uniqueness of the original description is needed:
the preceding refinement applies to every supplied actual description.

`a14_spectators_refine_to_four` (269–274) takes the second conjunct of
that proven refinement at `.a14 spectator`. The definition of `Survives`
reduces to the intended four-label membership. It proves a necessary
membership condition, not the realizability of those four alternatives.

All conditional statements have `written` EXPLICITLY in their binders.
The point, selected-vertex proof, exact16 equality and actual-forest proof
are explicit as appropriate; the section parameters occur in those types.
The old missing-`include`/body-only-section-variable failure mode is not
present. The finite `cases`, `exact`, `obtain` and closed `decide` proofs
have no apparent mathematical, motive or binder defect. No custom tactic,
nonlinear solver, local module import or external worker is invoked.

## 4. Archive chain and honest execution boundary

The four-entry archive manifest at
`results/source_archive/arithmetic16_frontier_before_d14_thirteen_group_consolidation_20261022/manifest.json`
freshly hashes to **c5d1c36a09c0057ae2cccb77b9c65551ac066887a28d58c2866eb8234ec2280e**.
All four gzip byte lengths/hashes match. Independent decoding yields:

| Preserved object | Decoded bytes | Exact decoded SHA256 |
|---|---:|---|
| Earlier source |8877|18cf9215ff771f6290ef23bce950d80e55eecf86c089174dc74d9c8fa1fa938c|
| Earlier map |10037|815ea993d22694048bc2c8bb62c9ccb896ccf769d72cb24a233772787ea69d7f|
| Earlier different-author static gate |8756|f4d730903b16a2a27e83e08048b3809f51e4f0d15c8cea32babc0c601352dfc2|
| Earlier root precompile gate |3661|da6ef32e9ff4cb072c5f15db0467025f66c3aa89f771bb4e556f33bffab253bc|

The latter two live historical files remain at those exact hashes. The
source diff adds the explicit D14 spectator type, its separate inputs and
derived exclusion, removes only D14 from the survivor list, and updates
the corresponding enumerations/print request. The old static gates are
correctly restricted to18cf/815e and are NOT automatically credited to
11dd/f41b. The repaired B14/A3/D14 mathematical predecessors likewise
remain historical; the map binds only their current corrected versions.

Exactly SIX actual theorem declarations have matching fully qualified
`#print axioms` requests at279–284. Those commands have not been run:
there is no current object, actual stdout or verified axiom list. The
module imports Std only and introduces no axiom, `sorry` or hidden
mathematical implementation. Abstract theorem PARAMETERS remain external
application obligations even if future prints contain no extra axioms.

The map correctly leaves actual geometry, finite classification, field
and angle restrictions, energy inequalities, manual tables and all prior
runtime evidence EXTERNAL. The result refines a necessary exact16 domain;
it proves no dimension-fifteen bound or existence/nonexistence of every
remaining forest, and changes no frozen paper or global ceiling.

**Verdict:** the consolidated source and map are ready for root's static
gate and its separately admitted capped compilation/actual axiom audit.
This report awards source-review credit ONLY, zero execution credit.
