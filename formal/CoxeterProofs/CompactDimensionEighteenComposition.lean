import Std

/-!
SOURCE PREPARATION ONLY: uncompiled, with no actual axiom-print output.

Conditional logic for the written global compact18 composition. The typed
forest, terminal and head descriptions are labels for EXTERNAL mathematical
interfaces. This file does not construct a polytope, classify a diagram,
derive a geometric energy inequality or instantiate those interfaces.
All predicates refer to one fixed candidate and the SAME actual vertex.
-/

namespace CoxeterProofs.CompactDimensionEighteenComposition

/-- The two different ACTUAL typed-seven defects of the written frontier. -/
inductive TypedSeven where
  | oldA2NewIsolatedSeven
  | oldI2SevenNewIsolatedThree

/-- The complete rank-three spectator alternatives beside A16. -/
inductive A16Spectator where
  | ordinaryA3
  | ordinaryT2OrdinaryA1
  | actualIOrdinaryT2
  | typedSevenOrdinaryA1 (meaning : TypedSeven)

inductive A17Spectator where
  | ordinaryT2
  | typedSeven (meaning : TypedSeven)

inductive A18Spectator where
  | ordinaryNonshortA1
  | actualI

/-- Nine grouped forest families, retaining their actual spectator types. -/
inductive Forest19 where
  | a16 (spectator : A16Spectator)
  | a17 (spectator : A17Spectator)
  | a18 (spectator : A18Spectator)
  | wholeA19
  | b18ActualI
  | wholeB19
  | d17TypedSeven (meaning : TypedSeven)
  | d18ActualI
  | wholeD19

/-- The two actual alternatives supplied by the whole-B19 reduction. -/
inductive B19Terminal where
  | shortThree
  | shortFour

/-- Six residual head classes AFTER the written complete D19 reduction. -/
inductive D19Head where
  | nativeFive
  | nativeSeven
  | nativeEight
  | nativeNine
  | finiteAtLeastTen (label : Nat) (atLeastTen : 10 ≤ label)
  | dotted

/-- The flags mean the SAME nonstrict .61 weakness convention at5 and15. -/
def weakEndCount (atFive atFifteen : Bool) : Nat :=
  (if atFive then 1 else 0) + (if atFifteen then 1 else 0)

/-- Equality is not lost: two Boolean sites have exactly0,1 or2 weak ends. -/
theorem weak_end_count_cases (atFive atFifteen : Bool) :
    weakEndCount atFive atFifteen = 0 ∨
      weakEndCount atFive atFifteen = 1 ∨
        weakEndCount atFive atFifteen = 2 := by
  cases atFive <;> cases atFifteen <;> decide

universe u

/--
Exactly23 named WRITTEN theorem interfaces, in the order of the frozen
composition's source table. This record contains no dimension<=18 premise
and no supplied whole-A19/B19/D19 impossibility. Geometry, source hashes,
actual-wall incidence, continuous inequalities and prior classifications
are not proved by declaring this record.
-/
structure WrittenInputs
    (Vertex : Type u) (dimension : Nat)
    (ActualFarthest : Vertex → Prop)
    (ActualForest : Vertex → Forest19 → Prop)
    (WeakAtFive WeakAtFifteen : Vertex → Bool)
    (ActualB19Terminal : Vertex → B19Terminal → Prop)
    (ActualD19Head : Vertex → D19Head → Prop) : Prop where
  G19 : dimension ≤ 19
  F19 : dimension = 19 → ∃ p forest,
    ActualFarthest p ∧ ActualForest p forest
  A16 : dimension = 19 → ∀ p spectator, ActualFarthest p →
    ActualForest p (.a16 spectator) → False
  A17o : dimension = 19 → ∀ p, ActualFarthest p →
    ActualForest p (.a17 .ordinaryT2) → False
  A17t : dimension = 19 → ∀ p meaning, ActualFarthest p →
    ActualForest p (.a17 (.typedSeven meaning)) → False
  A18o : dimension = 19 → ∀ p, ActualFarthest p →
    ActualForest p (.a18 .ordinaryNonshortA1) → False
  A18I : dimension = 19 → ∀ p, ActualFarthest p →
    ActualForest p (.a18 .actualI) → False
  B18I : dimension = 19 → ∀ p, ActualFarthest p →
    ActualForest p .b18ActualI → False
  D17t : dimension = 19 → ∀ p meaning, ActualFarthest p →
    ActualForest p (.d17TypedSeven meaning) → False
  D18I : dimension = 19 → ∀ p, ActualFarthest p →
    ActualForest p .d18ActualI → False
  B19r : dimension = 19 → ∀ p, ActualFarthest p →
    ActualForest p .wholeB19 → ∃ terminal, ActualB19Terminal p terminal
  B19four : dimension = 19 → ∀ p, ActualFarthest p →
    ActualForest p .wholeB19 → ActualB19Terminal p .shortFour → False
  B19three : dimension = 19 → ∀ p, ActualFarthest p →
    ActualForest p .wholeB19 → ActualB19Terminal p .shortThree → False
  A19two : dimension = 19 → ∀ p, ActualFarthest p →
    ActualForest p .wholeA19 →
    weakEndCount (WeakAtFive p) (WeakAtFifteen p) = 2 → False
  A19one : dimension = 19 → ∀ p, ActualFarthest p →
    ActualForest p .wholeA19 →
    weakEndCount (WeakAtFive p) (WeakAtFifteen p) = 1 → False
  A19zero : dimension = 19 → ∀ p, ActualFarthest p →
    ActualForest p .wholeA19 →
    weakEndCount (WeakAtFive p) (WeakAtFifteen p) = 0 → False
  DH : dimension = 19 → ∀ p, ActualFarthest p →
    ActualForest p .wholeD19 → ∃ head, ActualD19Head p head
  D5 : dimension = 19 → ∀ p, ActualFarthest p →
    ActualForest p .wholeD19 → ActualD19Head p .nativeFive → False
  D7 : dimension = 19 → ∀ p, ActualFarthest p →
    ActualForest p .wholeD19 → ActualD19Head p .nativeSeven → False
  D8 : dimension = 19 → ∀ p, ActualFarthest p →
    ActualForest p .wholeD19 → ActualD19Head p .nativeEight → False
  D9 : dimension = 19 → ∀ p, ActualFarthest p →
    ActualForest p .wholeD19 → ActualD19Head p .nativeNine → False
  Dlarge : dimension = 19 → ∀ p label (lower : 10 ≤ label),
    ActualFarthest p → ActualForest p .wholeD19 →
    ActualD19Head p (.finiteAtLeastTen label lower) → False
  DD : dimension = 19 → ∀ p, ActualFarthest p →
    ActualForest p .wholeD19 → ActualD19Head p .dotted → False

section Composition

variable {Vertex : Type u} {dimension : Nat}
variable {ActualFarthest : Vertex → Prop}
variable {ActualForest : Vertex → Forest19 → Prop}
variable {WeakAtFive WeakAtFifteen : Vertex → Bool}
variable {ActualB19Terminal : Vertex → B19Terminal → Prop}
variable {ActualD19Head : Vertex → D19Head → Prop}
variable (written : WrittenInputs Vertex dimension ActualFarthest ActualForest
  WeakAtFive WeakAtFifteen ActualB19Terminal ActualD19Head)

include written

/-- The complete three-case weak-end conjunction, including both orientations. -/
theorem whole_a19_exclusion
    (at19 : dimension = 19) (p : Vertex)
    (farthest : ActualFarthest p) (forest : ActualForest p .wholeA19) : False := by
  cases weak_end_count_cases (WeakAtFive p) (WeakAtFifteen p) with
  | inl zero => exact written.A19zero at19 p farthest forest zero
  | inr oneOrTwo =>
    cases oneOrTwo with
    | inl one => exact written.A19one at19 p farthest forest one
    | inr two => exact written.A19two at19 p farthest forest two

/-- The actual terminal witness is obtained once, then its alternatives close. -/
theorem whole_b19_exclusion
    (at19 : dimension = 19) (p : Vertex)
    (farthest : ActualFarthest p) (forest : ActualForest p .wholeB19) : False := by
  obtain ⟨terminal, actualTerminal⟩ := written.B19r at19 p farthest forest
  cases terminal with
  | shortThree => exact written.B19three at19 p farthest forest actualTerminal
  | shortFour => exact written.B19four at19 p farthest forest actualTerminal

/-- Complete native5/7/8/9/finite>=10/dotted dispatch at the SAME actual vertex. -/
theorem whole_d19_exclusion
    (at19 : dimension = 19) (p : Vertex)
    (farthest : ActualFarthest p) (forest : ActualForest p .wholeD19) : False := by
  obtain ⟨head, actualHead⟩ := written.DH at19 p farthest forest
  cases head with
  | nativeFive => exact written.D5 at19 p farthest forest actualHead
  | nativeSeven => exact written.D7 at19 p farthest forest actualHead
  | nativeEight => exact written.D8 at19 p farthest forest actualHead
  | nativeNine => exact written.D9 at19 p farthest forest actualHead
  | finiteAtLeastTen label lower =>
    exact written.Dlarge at19 p label lower farthest forest actualHead
  | dotted => exact written.DD at19 p farthest forest actualHead

/--
The old GLOBAL19 tail plus the NEW exact-dimension19 forest exclusions imply
GLOBAL18. No dimension heredity, Coxeter-face assumption, new optimizer or
whole-no19 conclusion is supplied as a premise.
-/
theorem compact_dimension_at_most_eighteen : dimension ≤ 18 := by
  by_cases bound : dimension ≤ 18
  · exact bound
  · have at19 : dimension = 19 := by
      have prior : dimension ≤ 19 := written.G19
      omega
    obtain ⟨p, forest, farthest, actualForest⟩ := written.F19 at19
    have impossible : False := by
      cases forest with
      | a16 spectator => exact written.A16 at19 p spectator farthest actualForest
      | a17 spectator =>
        cases spectator with
        | ordinaryT2 => exact written.A17o at19 p farthest actualForest
        | typedSeven meaning => exact written.A17t at19 p meaning farthest actualForest
      | a18 spectator =>
        cases spectator with
        | ordinaryNonshortA1 => exact written.A18o at19 p farthest actualForest
        | actualI => exact written.A18I at19 p farthest actualForest
      | wholeA19 => exact whole_a19_exclusion written at19 p farthest actualForest
      | b18ActualI => exact written.B18I at19 p farthest actualForest
      | wholeB19 => exact whole_b19_exclusion written at19 p farthest actualForest
      | d17TypedSeven meaning =>
        exact written.D17t at19 p meaning farthest actualForest
      | d18ActualI => exact written.D18I at19 p farthest actualForest
      | wholeD19 => exact whole_d19_exclusion written at19 p farthest actualForest
    exact False.elim impossible

end Composition
end CoxeterProofs.CompactDimensionEighteenComposition

#print axioms CoxeterProofs.CompactDimensionEighteenComposition.weak_end_count_cases
#print axioms CoxeterProofs.CompactDimensionEighteenComposition.whole_a19_exclusion
#print axioms CoxeterProofs.CompactDimensionEighteenComposition.whole_b19_exclusion
#print axioms CoxeterProofs.CompactDimensionEighteenComposition.whole_d19_exclusion
#print axioms CoxeterProofs.CompactDimensionEighteenComposition.compact_dimension_at_most_eighteen
