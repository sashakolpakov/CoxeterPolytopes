import CoxeterProofs.IdealRA4Gap

/-!
Conditional ideal RA4 interfaces. Volume numerator means volume in units
pi²/3, not in24-cell units. Divisibility by4 is DERIVED only in the
explicit Oct/D class. No universal decomposition is assumed or declared.
The published minimum, internal gap partition, equality spectrum reduction,
single-cut geometry and flag-transitive gluing are separate parameters.
-/
namespace CoxeterProofs.ManuscriptIdeal4

universe u

structure GapPartition (W : Nat) where
  m : Nat
  q : Nat
  s : Nat
  selected : List Nat
  rest : List Nat
  partition : selected.sum + rest.sum = 2 * W
  count : selected.length = 1 + q + s
  minimum : ∀ a ∈ selected, m ≤ a
  shape : IdealRA4.FacetCase m q s

theorem gap_from_partition (W : Nat) (divisibility : W % 3 = 0)
    (w : GapPartition W) : 12 ≤ W :=
  IdealRA4.defect_gap_of_partition w.selected w.rest divisibility
    w.partition w.count w.minimum w.shape

structure Data (P : Type u) where
  numerator : P → Nat
  ridgeDefect : P → Nat
  doubledFacets : P → Nat
  is24Cell : P → Prop
  isDouble : P → Prop
  octDSpectrum : P → Prop
  cut : P → P → P → Prop
  chain : P → Nat → Prop

/-- These are the published24-cell minimum and its normalization.
They must be instantiated for the ideal/right-angled finite-volume class. -/
structure PublishedMinimum {P : Type u} (d : Data P) : Prop where
  minimum : ∀ p, 4 ≤ d.numerator p
  minimumUnique : ∀ p, d.numerator p = 4 → d.is24Cell p
  cellNumerator : ∀ p, d.is24Cell p → d.numerator p = 4

/-- Written internal incidence/injection and facet-case reduction.
The selected lists are actual disjoint facet defects, not an assumed
second-volume conclusion. The volume identity uses the published formula
plus the internal ridge-defect double count. -/
structure GapGeometry {P : Type u} (d : Data P) : Prop where
  divisibility : ∀ p, d.ridgeDefect p % 3 = 0
  volumeIdentity : ∀ p, d.ridgeDefect p = 3 * (d.numerator p - 4)
  selectedFacetPartition :
    ∀ p, ¬ d.is24Cell p → Nonempty (GapPartition (d.ridgeDefect p))

theorem nonminimum_numerator_at_least8 {P : Type u} (d : Data P)
    (geometry : GapGeometry d) (p : P) (noncell : ¬ d.is24Cell p) :
    8 ≤ d.numerator p := by
  obtain ⟨w⟩ := geometry.selectedFacetPartition p noncell
  have gap := gap_from_partition _ (geometry.divisibility p) w
  exact IdealRA4.volume_numerator_at_least_eight gap (geometry.volumeIdentity p)

theorem minimum_or_second_gap {P : Type u} (d : Data P)
    (published : PublishedMinimum d) (geometry : GapGeometry d) (p : P) :
    d.numerator p = 4 ∨ 8 ≤ d.numerator p := by
  classical
  by_cases cell : d.is24Cell p
  · exact Or.inl (published.cellNumerator p cell)
  · exact Or.inr (nonminimum_numerator_at_least8 d geometry p cell)

theorem all_equal_of_minimum_sum (m : Nat) (xs : List Nat)
    (minimum : ∀ a ∈ xs, m ≤ a) (total : xs.sum = m * xs.length) :
    ∀ a ∈ xs, a = m := by
  induction xs with
  | nil => simp
  | cons a xs ih =>
    have head := minimum a (by simp)
    have tailMin : ∀ b ∈ xs, m ≤ b := by
      intro b hb
      exact minimum b (by simp [hb])
    have tailLower := CoxeterProofs.sum_lower_of_minimum m xs tailMin
    simp only [List.sum_cons, List.length_cons, Nat.mul_add, Nat.mul_one] at total
    have exactHead : a = m := by omega
    have exactTail : xs.sum = m * xs.length := by omega
    intro b hb
    rcases List.mem_cons.mp hb with same | inTail
    · omega
    · exact ih tailMin exactTail b inTail

theorem equality_positive_count_seven_or_eight (defects : List Nat)
    (total : defects.sum = 24) (minimum : ∀ a ∈ defects, 3 ≤ a)
    (injection : 7 ≤ defects.length) :
    defects.length = 7 ∨ defects.length = 8 := by
  have lower := CoxeterProofs.sum_lower_of_minimum 3 defects minimum
  omega

theorem equality_exactly_eight_all_three (defects : List Nat)
    (total : defects.sum = 24) (minimum : ∀ a ∈ defects, 3 ≤ a)
    (injection : 7 ≤ defects.length) (sevenExcluded : defects.length ≠ 7) :
    defects.length = 8 ∧ ∀ a ∈ defects, a = 3 := by
  have cases := equality_positive_count_seven_or_eight defects total minimum injection
  have eight : defects.length = 8 := by omega
  refine ⟨eight, all_equal_of_minimum_sum 3 defects minimum ?_⟩
  omega

/-- The Figure7 matching geometry supplies these four distinct neighbors
of defect at least5 and one further positive facet of defect at least4. -/
theorem minimum_four_budget_impossible (neighbors : List Nat) (extra : Nat)
    (count : neighbors.length = 4) (minimum : ∀ a ∈ neighbors, 5 ≤ a)
    (extraPositive : 4 ≤ extra) (budget : 4 + neighbors.sum + extra ≤ 24) :
    False := by
  have lower := CoxeterProofs.sum_lower_of_minimum 5 neighbors minimum
  omega

theorem minimum_at_least_five_forces_four (minimum count : Nat)
    (large : 5 ≤ minimum) (neighborInjection : 4 ≤ count)
    (budget : minimum * count ≤ 24) : count = 4 := by
  have lower := Nat.mul_le_mul_right count large
  omega

/-- Exact finite facts about the cuboctahedral N7 branch.
The twelve nonsquare cusps are model data; identifying an arbitrary
remaining branch with that model is still written geometry. -/
structure SevenBranchModel where
  nontriangles : Nat
  nonsquareCusps : Nat

structure SevenBranchScriptChecks (model : SevenBranchModel) : Prop where
  nontriangles : model.nontriangles = 6
  nonsquareCusps : model.nonsquareCusps = 12

theorem seven_branch_injection_impossible (model : SevenBranchModel)
    (script : SevenBranchScriptChecks model) (count : Nat)
    (seven : count = 7)
    (actualInjection : 1 + model.nontriangles + model.nonsquareCusps ≤ count) :
    False := by
  have q := script.nontriangles
  have s := script.nonsquareCusps
  omega

/-- A finite binary tree of actual octahedral cuts. Its leaf count is
recorded, but it need not be a path before the separate path interface. -/
inductive Decomposition {P : Type u} (d : Data P) : P → Nat → Prop where
  | cell {p : P} : d.is24Cell p → Decomposition d p 1
  | split {p l r : P} {a b : Nat} :
      d.cut p l r → Decomposition d l a → Decomposition d r b →
      Decomposition d p (a + b)

theorem decomposition_positive {P : Type u} {d : Data P}
    {p : P} {count : Nat} (h : Decomposition d p count) : 0 < count := by
  induction h with
  | cell _ => decide
  | split _ _ _ left right => omega

theorem one_cell_recognition {P : Type u} {d : Data P}
    {p : P} {count : Nat} (h : Decomposition d p count) (one : count = 1) :
    d.is24Cell p := by
  cases h with
  | cell cell => exact cell
  | split cut left right =>
    have lp := decomposition_positive left
    have rp := decomposition_positive right
    omega

/-- Internal geometric recognition, all at the explicit Oct/D spectrum.
The cut premise includes actual separating octahedral sections, unchanged
outside facets, right angles, and the decrease of exactly eight D facets.
No statement about arbitrary ideal RA4 polytopes occurs here. -/
structure OctDGeometry {P : Type u} (d : Data P) : Prop where
  zeroRecognition :
    ∀ p, d.octDSpectrum p → d.doubledFacets p = 0 → d.is24Cell p
  positiveCut :
    ∀ p, d.octDSpectrum p → 0 < d.doubledFacets p →
    ∃ l r, d.octDSpectrum l ∧ d.octDSpectrum r ∧ d.cut p l r ∧
      d.doubledFacets l + d.doubledFacets r + 8 = d.doubledFacets p
  numeratorAdditive :
    ∀ p l r, d.cut p l r → d.numerator p = d.numerator l + d.numerator r
  pathRecognition :
    ∀ p count, d.octDSpectrum p → Decomposition d p count → d.chain p count
  canonicalTwoCellGluing :
    ∀ p l r, d.cut p l r → d.is24Cell l → d.is24Cell r → d.isDouble p

theorem octD_decomposition_and_counts {P : Type u} (d : Data P)
    (published : PublishedMinimum d) (geometry : OctDGeometry d)
    (p : P) (spectrum : d.octDSpectrum p) :
    ∃ k, d.doubledFacets p = 8 * k ∧
      d.numerator p = 4 * (k + 1) ∧
      Decomposition d p (k + 1) ∧ d.chain p (k + 1) := by
  have aux : ∀ fuel, ∀ p : P, d.doubledFacets p = fuel →
      d.octDSpectrum p →
      ∃ k, d.doubledFacets p = 8 * k ∧
        d.numerator p = 4 * (k + 1) ∧ Decomposition d p (k + 1) := by
    intro fuel
    induction fuel using Nat.strongRecOn with
    | ind fuel ih =>
      intro p measure admissible
      by_cases zero : d.doubledFacets p = 0
      · have cell := geometry.zeroRecognition p admissible zero
        refine ⟨0, ?_, ?_, ?_⟩
        · omega
        · have volume := published.cellNumerator p cell
          omega
        · exact Decomposition.cell cell
      · obtain ⟨l, r, hl, hr, cut, decrease⟩ :=
          geometry.positiveCut p admissible (by omega)
        have lowerLeft : d.doubledFacets l < fuel := by omega
        have lowerRight : d.doubledFacets r < fuel := by omega
        obtain ⟨a, ha, va, da⟩ := ih _ lowerLeft l rfl hl
        obtain ⟨b, hb, vb, db⟩ := ih _ lowerRight r rfl hr
        have volume := geometry.numeratorAdditive p l r cut
        refine ⟨a + b + 1, ?_, ?_, ?_⟩
        · omega
        · omega
        · have joined := Decomposition.split cut da db
          have count : (a + 1) + (b + 1) = a + b + 1 + 1 := by omega
          rw [count] at joined
          exact joined
  obtain ⟨k, hk, hv, hd⟩ := aux (d.doubledFacets p) p rfl spectrum
  exact ⟨k, hk, hv, hd, geometry.pathRecognition p (k + 1) spectrum hd⟩

theorem two_cell_decomposition_is_double {P : Type u} (d : Data P)
    (geometry : OctDGeometry d) {p : P} {count : Nat}
    (h : Decomposition d p count) (two : count = 2) : d.isDouble p := by
  cases h with
  | cell cell => omega
  | split cut left right =>
    have lp := decomposition_positive left
    have rp := decomposition_positive right
    exact geometry.canonicalTwoCellGluing _ _ _ cut
      (one_cell_recognition left (by omega))
      (one_cell_recognition right (by omega))

/-- This is the INTERNAL first half of the equality proof: defect12
forces the Oct/D facet spectrum. It is not the final equality theorem,
and is not a published assumption or a universal decomposition claim. -/
structure EqualitySpectrum {P : Type u} (d : Data P) : Prop where
  spectrumAtNumerator8 : ∀ p, d.numerator p = 8 → d.octDSpectrum p

/-- More granular construction of the spectrum interface. Minimum3 uses
the antiprism and minimum4/large-minimum exclusions. The N7 geometry uses
the no-positive-triangle rule, the Figure7 matching, and embedded medial
recovery before the script model facts and cusp injection apply. -/
structure WrittenEqualityStages {P : Type u} (d : Data P)
    (model : SevenBranchModel) where
  positiveDefects : P → List Nat
  total24 : ∀ p, d.numerator p = 8 → (positiveDefects p).sum = 24
  minimum3 : ∀ p, d.numerator p = 8 →
    ∀ a ∈ positiveDefects p, 3 ≤ a
  sevenFacetInjection : ∀ p, d.numerator p = 8 →
    7 ≤ (positiveDefects p).length
  remainingSevenBranchInjection : ∀ p, d.numerator p = 8 →
    (positiveDefects p).length = 7 →
    1 + model.nontriangles + model.nonsquareCusps ≤ (positiveDefects p).length
  classifyDefectZeroAndThree : ∀ p, d.numerator p = 8 →
    (positiveDefects p).length = 8 →
    (∀ a ∈ positiveDefects p, a = 3) → d.octDSpectrum p

theorem equality_spectrum_from_stages {P : Type u} (d : Data P)
    (model : SevenBranchModel) (script : SevenBranchScriptChecks model)
    (geometry : WrittenEqualityStages d model) : EqualitySpectrum d := by
  refine ⟨?_⟩
  intro p volume
  have sevenExcluded : (geometry.positiveDefects p).length ≠ 7 := by
    intro seven
    exact seven_branch_injection_impossible model script _
      seven (geometry.remainingSevenBranchInjection p volume seven)
  obtain ⟨eight, allThree⟩ := equality_exactly_eight_all_three
    (geometry.positiveDefects p) (geometry.total24 p volume)
    (geometry.minimum3 p volume) (geometry.sevenFacetInjection p volume) sevenExcluded
  exact geometry.classifyDefectZeroAndThree p volume eight allThree

theorem equality8_is_two_cell_chain_and_double {P : Type u} (d : Data P)
    (published : PublishedMinimum d) (geometry : OctDGeometry d)
    (reduction : EqualitySpectrum d) (p : P) (volume : d.numerator p = 8) :
    d.chain p 2 ∧ d.isDouble p := by
  obtain ⟨k, hk, hv, hd, hc⟩ :=
    octD_decomposition_and_counts d published geometry p
      (reduction.spectrumAtNumerator8 p volume)
  have count : k + 1 = 2 := by omega
  have chain : d.chain p 2 := by simpa only [count] using hc
  exact ⟨chain, two_cell_decomposition_is_double d geometry hd count⟩

/-- Here units means an explicitly supplied conversion to24-cell volume
units; no universal integrality in those larger units is assumed. -/
theorem two_volume_units_is_double {P : Type u} (d : Data P)
    (published : PublishedMinimum d) (geometry : OctDGeometry d)
    (reduction : EqualitySpectrum d) (p : P) (units : Nat)
    (conversion : d.numerator p = 4 * units) (two : units = 2) :
    d.chain p 2 ∧ d.isDouble p :=
  equality8_is_two_cell_chain_and_double d published geometry reduction p (by omega)

theorem minimum_gap_and_equality {P : Type u} (d : Data P)
    (published : PublishedMinimum d) (gap : GapGeometry d)
    (cuts : OctDGeometry d) (reduction : EqualitySpectrum d) (p : P) :
    4 ≤ d.numerator p ∧
      (d.numerator p = 4 ↔ d.is24Cell p) ∧
      (¬ d.is24Cell p → 8 ≤ d.numerator p) ∧
      (d.numerator p = 8 → d.chain p 2 ∧ d.isDouble p) := by
  refine ⟨published.minimum p, ?_, ?_, ?_⟩
  · exact ⟨published.minimumUnique p, published.cellNumerator p⟩
  · exact nonminimum_numerator_at_least8 d gap p
  · exact equality8_is_two_cell_chain_and_double d published cuts reduction p

end CoxeterProofs.ManuscriptIdeal4

#print axioms CoxeterProofs.ManuscriptIdeal4.gap_from_partition
#print axioms CoxeterProofs.ManuscriptIdeal4.minimum_or_second_gap
#print axioms CoxeterProofs.ManuscriptIdeal4.equality_exactly_eight_all_three
#print axioms CoxeterProofs.ManuscriptIdeal4.equality_spectrum_from_stages
#print axioms CoxeterProofs.ManuscriptIdeal4.octD_decomposition_and_counts
#print axioms CoxeterProofs.ManuscriptIdeal4.equality8_is_two_cell_chain_and_double
#print axioms CoxeterProofs.ManuscriptIdeal4.two_volume_units_is_double
#print axioms CoxeterProofs.ManuscriptIdeal4.minimum_gap_and_equality
