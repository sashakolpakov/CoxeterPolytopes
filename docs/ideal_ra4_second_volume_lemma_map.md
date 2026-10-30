# Ideal RA4 second value: precise lemma and assumption map

For current reading, see the [current second-volume manuscript](../arxiv/ideal_four_dimensional_volume_and_structure/main.tex) and the
[submission guide](arxiv_submission_editions_20261028.md). The supporting
arguments and recorded verification scopes below remain distinct from those
current submission builds. Standalone draft typesetting routes have been retired.

The supporting argument proves the second **distinct value** $8\pi^2/3$. The [frozen proof](../results/ideal_ra4_second_volume_20261018.md) is authoritative. Equality uniqueness is not part of this packet. All polytopes are finite-sided, finite-volume, entirely ideal and right-angled in dimension4.

## Published inputs

|ID|Exact input|Source and use|
|---|---|---|
|E1|An ideal right-angled vertex link is a Euclidean rectangular3-box; facets are ideal right-angled3-polyhedra.|Standard Coxeter cusp geometry, used in [Kolpakov2012](https://msp.org/agt/2012/12-4/agt-v12-n4-p03-s.pdf), preliminaries and Lemma1.|
|E2|$e=4v$, $r=3v+f$, $\sum_Rk(R)=12v$, and $mathrm{Vol}=(v-f+4)\pi^2/3$.|Kolpakov Lemmas1–2.|
|E3|The ideal right-angled24-cell is the unique minimum-volume ideal RA4 polytope, of volume $4\pi^2/3$.|Kolpakov Theorem3. Used when all defects vanish and for the original cell's volume.|
|E4|A k-face of an ideal RA3 polyhedron forces at least2k vertices; equality is the k-antiprism.|Kolpakov Proposition5, including its equality statement.|
|E5|The complete admissible nine-vertex octahedrite list has one member, the double octahedron.|Kolpakov §4.3, CaseD2 and Figure6. Only the nine-vertex entry is required for the gap.|
|E6|The indicated ideal right-angled3-polyhedral combinatorial types have unique geometric realizations.|Standard Andreev rigidity, as used in Kolpakov's proof. It turns E5 into the stated geometric double and identifies the octahedron/antiprism models.|

No new four-dimensional bounded census is assumed. The opposite-cusp injection below is **internal**, not an external input.

## Internal quantified statements

|ID|Statement and proof obligation|Status/interface|
|---|---|---|
|I1 `RidgeDefect`|For every P in the stated class, $W=\sum_R(k_R-3)=3(v-f)$, $W\ge0$, and its volume is $4\pi^2/3+W\pi^2/9$.|Algebra from E2; no extra divisibility assumption.|
|I2 `FacetDefect`|Every facet has $d=v_F-6=\sum_{R\subset F}(k_R-3)$; $\sum_Fd=2W$. Defect0 is octahedral; every positive defect is at least2.|Valency4/Euler, double counting, E4 and E6.|
|I3 `MinimumPositiveFacet`|If P is not24-cell, its finite facet set has a smallest positive defect m.|I2 plus E3: absence of positive defect would force W0.|
|I4 `OppositeCuspInjection`|For every positive facet F, its q nontriangle neighbors and its s nonsquare-cusp opposites are distinct positive-defect facets, disjoint from F and from each other.|**New internal geometry.** Opposite rectangular-box faces are congruent; asymptotic supporting hyperplanes meet at only one ideal point and cannot also meet along a ridge. This constructs the finite injection, rather than assuming a numerical inequality.|
|I5 `SelectedDefectBudget`|If m is minimum positive defect, I4 implies $d(F)+m(q+s)\le2W$, in particular $m(1+q+s)\le2W$ when d(F)=m.|Finite disjoint-list partition and nonnegative defects. The root-owned Lean module proves this list-to-budget algebra.|
|I6 `AntiprismShape`|For every integer k≥4 the ideal RA k-antiprism has2k nonsquare cusps of aspect $2\cos(\pi/k)>1$.|**Internal Lorentz/horosphere argument**, seven exact rational identities plus the unit-normal-sum interpretation of width. The checker verifies identities, not the universal interpretation by itself.|
|I7 `DoubleOctShape`|The nine-vertex facet has q3 and exactly3 nonsquare cusps of aspect2.|E5–E6 plus internal unfolding of two cusp squares. The example checks this shape independently but does not establish classification completeness.|
|I8 `SmallNontriangleCount`|A positive facet has q≥2; q2 forces equal face sizes and the antiprism equality case.|I2 and E4.|
|I9 `ExhaustiveMinimumCases`|The minimum facet is antiprism; or m3/q3/s3; or m4/q≥4; or m≥5/q≥3. The non-antiprism m2 branch is impossible.|I6–I8 and E4. Every positive m is covered.|
|I10 `DefectGap`|The preceding budget and cases force2W≥20 and W divisible3, hence W≥12; W3,6,9 are excluded.|Complete integer arithmetic, formalized conditionally in `IdealRA4Gap.lean`.|
|I11 `ExactDouble`|The reflected38-halfspace intersection equals the union of two24-cells and has f-vector(42,168,164,38), W12 and volume $8\pi^2/3$.|Internal halfspace-union proof plus exact rational construction/face checker. No hidden escape or vertex domain is omitted.|
|I12 `SecondDistinctValue`|E3, I10 and I11 give the second distinct volume value.|Written theorem. Equality classification is an additional problem.|

## Lean interfaces and remaining geometric formalization

[IdealRA4Gap.lean](../formal/CoxeterProofs/IdealRA4Gap.lean) defines `FacetCase` for the conclusion of I9. `defect_gap_of_partition` accepts two finite lists `selected` and `rest` with total sum `2*W`, selected length `1+q+s`, and all selected entries at least m. Its conclusions include `defect_gap`, `intermediate_defects_excluded` and `volume_numerator_at_least_eight`. Exact double numerators are also checked.

This is conditional theorem composition, not an axiom declaring I4 or I9 true for geometric polytopes. An end-to-end formal development must construct the facet/cusp types and the injections in I4, prove the horospherical width calculation I6, identify the E5 type using the explicitly declared published input, and connect the exact halfspace data to I11. Those internal obligations remain visible rather than being relabeled as external assumptions.

## Computational scope

The [producer](../scripts/fresh_ideal_ra4_double24_volume_20261018.py) and [standalone checker](../scripts/check_ideal_ra4_second_volume_saved_20261018.py) establish exact construction/arithmetic claims. They do not classify all4-polytopes, prove the new cusp injection merely by testing one example, or decide every W12 equality case. The checker uses separate elimination and does not import the producer. See the [reproduction guide](ideal_ra4_second_volume_reproduce.md) for hashes, commands and actual final run IDs.
