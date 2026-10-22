import Std

/-!
SOURCE PREPARATION ONLY: uncompiled, zero actual axiom-print outputs.

Conditional logic for the written global arithmetic16 composition. The
fourteen forest groups, actual-field labels and endpoint/census predicates
are EXTERNAL mathematical interfaces, not geometric instances. All inputs
concern one fixed candidate and the SAME actual farthest vertex and pair.
Only Std is imported. No field arithmetic or census is executed here.
-/

namespace CoxeterProofs.ArithmeticDimensionSixteenComposition

inductive TypedSeven where
  | oldA2NewIsolatedSeven
  | oldI2SevenNewIsolatedThree

/-- An ordinary irreducible dihedral, not a typed-seven defect. -/
structure OrdinaryT2 where
  label : Nat
  labelAtLeastThree : 3 ≤ label

inductive OrdinaryT3 where
  | a3
  | b3
  | h3

inductive OrdinaryT4 where
  | a4
  | b4
  | d4
  | f4
  | h4

/-- H4 is absent only in this exact A11 typed-seven row. -/
inductive A11TypedSevenT4 where
  | a4
  | b4
  | d4
  | f4

inductive A11ActualIRemainder where
  | a5
  | d5

inductive A12OrdinaryRemainder where
  | a5
  | d5
  | a4OrdinaryA1
  | h4OrdinaryA1
  | t3T2 (rankThree : OrdinaryT3) (rankTwo : OrdinaryT2)
  | twoT2OrdinaryA1 (first second : OrdinaryT2)

inductive A12ActualIRemainder where
  | t4 (rankFour : OrdinaryT4)
  | t3OrdinaryA1 (rankThree : OrdinaryT3)
  | twoT2 (first second : OrdinaryT2)

inductive A14Remainder where
  | t3 (rankThree : OrdinaryT3)
  | t2OrdinaryA1 (rankTwo : OrdinaryT2)
  | actualIT2 (rankTwo : OrdinaryT2)

inductive A15Remainder where
  | ordinaryT2 (rankTwo : OrdinaryT2)
  | twoOrdinaryA1
  | actualIOrdinaryA1
  | typedSeven (meaning : TypedSeven)

inductive A16Remainder where
  | ordinaryA1
  | actualI

inductive D16Remainder where
  | ordinaryNonshortA1
  | actualI

/-- The FOURTEEN original groups of the fresh dimension17 frontier. -/
inductive Forest17 where
  | a11ActualI (remainder : A11ActualIRemainder)
  | a11TypedSeven (meaning : TypedSeven) (rankFour : A11TypedSevenT4)
  | a12Ordinary (remainder : A12OrdinaryRemainder)
  | a12ActualI (remainder : A12ActualIRemainder)
  | a14 (remainder : A14Remainder)
  | a15 (remainder : A15Remainder)
  | a16 (remainder : A16Remainder)
  | wholeA17
  | b15TypedSeven (meaning : TypedSeven)
  | b16ActualI
  | wholeB17
  | d15TypedSeven (meaning : TypedSeven)
  | d16 (remainder : D16Remainder)
  | wholeD17

/-- Labels for the actual-field conclusion, not a classification axiom. -/
inductive A17Field where
  | sqrtTwo
  | sqrtThree
  | nativeSeven
  | nativeEighteen

/-- The two finite cubic anchors of the complete conditional census. -/
inductive CubicAnchor where
  | nativeSeven
  | nativeEighteen

universe u v

/--
The old GLOBAL17 tail, the fresh FOURTEEN-group frontier and the thirteen
already excluded groups. The desired GLOBAL16 bound and whole-A17
impossibility are not fields of this record. All source-level geometric,
arithmetic and inherited exact-data premises remain external.
-/
structure WrittenForestInputs
    (Vertex : Type u) (dimension : Nat)
    (ActualFarthest : Vertex → Prop)
    (ActualForest : Vertex → Forest17 → Prop) : Prop where
  G17 : dimension ≤ 17
  F17 : dimension = 17 → ∃ p forest,
    ActualFarthest p ∧ ActualForest p forest
  A11I : dimension = 17 → ∀ p remainder, ActualFarthest p →
    ActualForest p (.a11ActualI remainder) → False
  AL11 : dimension = 17 → ∀ p meaning rankFour, ActualFarthest p →
    ActualForest p (.a11TypedSeven meaning rankFour) → False
  AL12ordinary : dimension = 17 → ∀ p remainder, ActualFarthest p →
    ActualForest p (.a12Ordinary remainder) → False
  AL12I : dimension = 17 → ∀ p remainder, ActualFarthest p →
    ActualForest p (.a12ActualI remainder) → False
  A14 : dimension = 17 → ∀ p remainder, ActualFarthest p →
    ActualForest p (.a14 remainder) → False
  A15T : dimension = 17 → ∀ p rankTwo, ActualFarthest p →
    ActualForest p (.a15 (.ordinaryT2 rankTwo)) → False
  A15twoA1 : dimension = 17 → ∀ p, ActualFarthest p →
    ActualForest p (.a15 .twoOrdinaryA1) → False
  A15I : dimension = 17 → ∀ p, ActualFarthest p →
    ActualForest p (.a15 .actualIOrdinaryA1) → False
  A15seven : dimension = 17 → ∀ p meaning, ActualFarthest p →
    ActualForest p (.a15 (.typedSeven meaning)) → False
  A16 : dimension = 17 → ∀ p remainder, ActualFarthest p →
    ActualForest p (.a16 remainder) → False
  B15 : dimension = 17 → ∀ p meaning, ActualFarthest p →
    ActualForest p (.b15TypedSeven meaning) → False
  B16I : dimension = 17 → ∀ p, ActualFarthest p →
    ActualForest p .b16ActualI → False
  WB : dimension = 17 → ∀ p, ActualFarthest p →
    ActualForest p .wholeB17 → False
  D15 : dimension = 17 → ∀ p meaning, ActualFarthest p →
    ActualForest p (.d15TypedSeven meaning) → False
  D16ordinary : dimension = 17 → ∀ p, ActualFarthest p →
    ActualForest p (.d16 .ordinaryNonshortA1) → False
  D16I : dimension = 17 → ∀ p, ActualFarthest p →
    ActualForest p (.d16 .actualI) → False
  WD : dimension = 17 → ∀ p, ActualFarthest p →
    ActualForest p .wholeD17 → False

/--
Separately named endpoint, degree, finite-anchor/actual-field, domain and
exhaustive-empty interfaces. ActualGroundField labels arise only AFTER
the finite-anchor theorem. No field list is assumed for every candidate.
The empty premises concern complete conditional algebraic domains, not
an assumed absence of whole A17 or dimension17.
-/
structure WrittenA17Inputs
    (Vertex : Type u) (EndpointPair : Type v) (dimension : Nat)
    (ActualFarthest : Vertex → Prop)
    (ActualForest : Vertex → Forest17 → Prop)
    (ActualEndpoints : Vertex → EndpointPair → Prop)
    (BothIsolated : EndpointPair → Prop)
    (actualDegree : Vertex → Nat)
    (ActualGroundField ActualFiniteAnchor : Vertex → EndpointPair → A17Field → Prop)
    (QuadraticDomain : Vertex → EndpointPair → Prop)
    (CubicDomain : Vertex → EndpointPair → CubicAnchor → Prop) : Prop where
  endpointIsolation : dimension = 17 → ∀ p, ActualFarthest p →
    ActualForest p .wholeA17 →
    ∃ pair, ActualEndpoints p pair ∧ BothIsolated pair
  endpointDegree : dimension = 17 → ∀ p pair, ActualFarthest p →
    ActualForest p .wholeA17 → ActualEndpoints p pair → BothIsolated pair →
    2 ≤ actualDegree p ∧ actualDegree p ≤ 5
  finiteAnchorAndActualField : dimension = 17 → ∀ p pair,
    ActualFarthest p → ActualForest p .wholeA17 →
    ActualEndpoints p pair → BothIsolated pair →
    2 ≤ actualDegree p → actualDegree p ≤ 5 →
    ∃ field, ActualGroundField p pair field ∧ ActualFiniteAnchor p pair field
  quadraticDomain : dimension = 17 → ∀ p pair field,
    ActualFarthest p → ActualForest p .wholeA17 →
    ActualEndpoints p pair → BothIsolated pair → ActualGroundField p pair field →
    (field = .sqrtTwo ∨ field = .sqrtThree) → QuadraticDomain p pair
  quadraticExhaustiveEmpty : ∀ p pair, QuadraticDomain p pair → False
  nativeSevenDomain : dimension = 17 → ∀ p pair,
    ActualFarthest p → ActualForest p .wholeA17 →
    ActualEndpoints p pair → BothIsolated pair →
    ActualGroundField p pair .nativeSeven → ActualFiniteAnchor p pair .nativeSeven →
    CubicDomain p pair .nativeSeven
  nativeEighteenDomain : dimension = 17 → ∀ p pair,
    ActualFarthest p → ActualForest p .wholeA17 →
    ActualEndpoints p pair → BothIsolated pair →
    ActualGroundField p pair .nativeEighteen →
    ActualFiniteAnchor p pair .nativeEighteen → CubicDomain p pair .nativeEighteen
  cubicExhaustiveEmpty : ∀ p pair anchor, CubicDomain p pair anchor → False

section Composition

variable {Vertex : Type u} {EndpointPair : Type v} {dimension : Nat}
variable {ActualFarthest : Vertex → Prop}
variable {ActualForest : Vertex → Forest17 → Prop}
variable {ActualEndpoints : Vertex → EndpointPair → Prop}
variable {BothIsolated : EndpointPair → Prop} {actualDegree : Vertex → Nat}
variable {ActualGroundField ActualFiniteAnchor : Vertex → EndpointPair → A17Field → Prop}
variable {QuadraticDomain : Vertex → EndpointPair → Prop}
variable {CubicDomain : Vertex → EndpointPair → CubicAnchor → Prop}
variable (forestInputs : WrittenForestInputs Vertex dimension ActualFarthest ActualForest)
variable (a17Inputs : WrittenA17Inputs Vertex EndpointPair dimension
  ActualFarthest ActualForest ActualEndpoints BothIsolated actualDegree
  ActualGroundField ActualFiniteAnchor QuadraticDomain CubicDomain)

include a17Inputs

/-- The same actual pair follows one of FOUR field routes to an empty census domain. -/
theorem whole_a17_exclusion
    (at17 : dimension = 17) (p : Vertex)
    (farthest : ActualFarthest p) (forest : ActualForest p .wholeA17) : False := by
  obtain ⟨pair, actualPair, isolated⟩ :=
    a17Inputs.endpointIsolation at17 p farthest forest
  obtain ⟨degreeLower, degreeUpper⟩ :=
    a17Inputs.endpointDegree at17 p pair farthest forest actualPair isolated
  obtain ⟨field, actualField, finiteAnchor⟩ :=
    a17Inputs.finiteAnchorAndActualField at17 p pair farthest forest
      actualPair isolated degreeLower degreeUpper
  cases field with
  | sqrtTwo =>
    exact a17Inputs.quadraticExhaustiveEmpty p pair
      (a17Inputs.quadraticDomain at17 p pair .sqrtTwo farthest forest
        actualPair isolated actualField (Or.inl rfl))
  | sqrtThree =>
    exact a17Inputs.quadraticExhaustiveEmpty p pair
      (a17Inputs.quadraticDomain at17 p pair .sqrtThree farthest forest
        actualPair isolated actualField (Or.inr rfl))
  | nativeSeven =>
    exact a17Inputs.cubicExhaustiveEmpty p pair .nativeSeven
      (a17Inputs.nativeSevenDomain at17 p pair farthest forest
        actualPair isolated actualField finiteAnchor)
  | nativeEighteen =>
    exact a17Inputs.cubicExhaustiveEmpty p pair .nativeEighteen
      (a17Inputs.nativeEighteenDomain at17 p pair farthest forest
        actualPair isolated actualField finiteAnchor)

include forestInputs

/-- The FRESH fourteen-group frontier closes without ordinary/I contact transfer. -/
theorem dimension_seventeen_exclusion (at17 : dimension = 17) : False := by
  obtain ⟨p, forest, farthest, actualForest⟩ := forestInputs.F17 at17
  cases forest with
  | a11ActualI remainder =>
    exact forestInputs.A11I at17 p remainder farthest actualForest
  | a11TypedSeven meaning rankFour =>
    exact forestInputs.AL11 at17 p meaning rankFour farthest actualForest
  | a12Ordinary remainder =>
    exact forestInputs.AL12ordinary at17 p remainder farthest actualForest
  | a12ActualI remainder =>
    exact forestInputs.AL12I at17 p remainder farthest actualForest
  | a14 remainder => exact forestInputs.A14 at17 p remainder farthest actualForest
  | a15 remainder =>
    cases remainder with
    | ordinaryT2 rankTwo =>
      exact forestInputs.A15T at17 p rankTwo farthest actualForest
    | twoOrdinaryA1 => exact forestInputs.A15twoA1 at17 p farthest actualForest
    | actualIOrdinaryA1 => exact forestInputs.A15I at17 p farthest actualForest
    | typedSeven meaning =>
      exact forestInputs.A15seven at17 p meaning farthest actualForest
  | a16 remainder => exact forestInputs.A16 at17 p remainder farthest actualForest
  | wholeA17 => exact whole_a17_exclusion a17Inputs at17 p farthest actualForest
  | b15TypedSeven meaning =>
    exact forestInputs.B15 at17 p meaning farthest actualForest
  | b16ActualI => exact forestInputs.B16I at17 p farthest actualForest
  | wholeB17 => exact forestInputs.WB at17 p farthest actualForest
  | d15TypedSeven meaning =>
    exact forestInputs.D15 at17 p meaning farthest actualForest
  | d16 remainder =>
    cases remainder with
    | ordinaryNonshortA1 =>
      exact forestInputs.D16ordinary at17 p farthest actualForest
    | actualI => exact forestInputs.D16I at17 p farthest actualForest
  | wholeD17 => exact forestInputs.WD at17 p farthest actualForest

/-- The prior GLOBAL17 tail plus the new exact-dimension17 contradiction gives16. -/
theorem arithmetic_dimension_at_most_sixteen : dimension ≤ 16 := by
  by_cases bound : dimension ≤ 16
  · exact bound
  · have at17 : dimension = 17 := by
      have prior : dimension ≤ 17 := forestInputs.G17
      omega
    exact False.elim (dimension_seventeen_exclusion forestInputs a17Inputs at17)

end Composition
end CoxeterProofs.ArithmeticDimensionSixteenComposition

#print axioms CoxeterProofs.ArithmeticDimensionSixteenComposition.whole_a17_exclusion
#print axioms CoxeterProofs.ArithmeticDimensionSixteenComposition.dimension_seventeen_exclusion
#print axioms CoxeterProofs.ArithmeticDimensionSixteenComposition.arithmetic_dimension_at_most_sixteen
