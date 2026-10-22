import Std

/-!
CORRECTED SOURCE PREPARATION: no successful compilation or axiom audit yet.
The first failed compile is preserved; retained membership proofs now
unfold Survives explicitly instead of invoking decide on an opaque goal.

Conditional refinement of the FRESH rank-sixteen necessary forest list.
The geometry, the complete nineteen-group frontier and eight written
exclusion interfaces are EXTERNAL, not instances supplied by this file.
The D14 whole-row exclusion is DERIVED from its three spectator cases.
All conditional deductions retain the SAME selected vertex and forest.
This module proves no dimension-fifteen or geometric energy theorem.
-/

namespace CoxeterProofs.ArithmeticSixteenFrontierRefinement

/-- Both alternatives include their ACTUAL isolated new opposite wall. -/
inductive TypedSeven where
  | oldA2NewIsolatedSeven
  | oldI2SevenNewIsolatedThree
  deriving DecidableEq, Repr

/-- The original nineteen GROUPS; a constructor may group spectator types. -/
inductive Group where
  | a10ActualI
  | a10TypedSeven
  | a11Ordinary
  | a11ActualI
  | a11OldSevenA3
  | a12Ordinary
  | a12ActualI
  | a12TypedSeven
  | a13
  | a14
  | a15
  | wholeA16
  | b14TypedSeven
  | b15
  | wholeB16
  | d13
  | d14
  | d15
  | wholeD16
  deriving DecidableEq, Repr

/-- Five original A14 alternatives, with the two typed meanings separate. -/
inductive A14Spectator where
  | ordinaryT2
  | twoOrdinaryA1
  | actualIOrdinaryA1
  | oldA2NewIsolatedSeven
  | oldI2SevenNewIsolatedThree
  deriving DecidableEq, Repr

/-- The complete original D14 row; ordinary means I2(m), m >= 4. -/
inductive D14Spectator where
  | ordinaryDihedral
  | oldA2NewIsolatedSeven
  | oldI2SevenNewIsolatedThree
  deriving DecidableEq, Repr

def D14Spectator.ofTypedSeven : TypedSeven → D14Spectator
  | .oldA2NewIsolatedSeven => .oldA2NewIsolatedSeven
  | .oldI2SevenNewIsolatedThree => .oldI2SevenNewIsolatedThree

/--
Labels for the written whole-forest alternatives, not constructed diagrams.
Ordinary spectator subtypes are those of the pinned fresh frontier.
The actual-I and ordinary rows, and both typed-seven meanings, stay distinct.
-/
inductive ForestCase where
  | a10ActualI
  | a10TypedSeven (meaning : TypedSeven)
  | a11Ordinary
  | a11ActualI
  | a11OldSevenA3
  | a12Ordinary
  | a12ActualI
  | a12TypedSeven (meaning : TypedSeven)
  | a13
  | a14 (spectator : A14Spectator)
  | a15
  | wholeA16
  | b14TypedSeven (meaning : TypedSeven)
  | b15
  | wholeB16
  | d13
  | d14 (spectator : D14Spectator)
  | d15
  | wholeD16
  deriving DecidableEq, Repr

def ForestCase.group : ForestCase → Group
  | .a10ActualI => .a10ActualI
  | .a10TypedSeven _ => .a10TypedSeven
  | .a11Ordinary => .a11Ordinary
  | .a11ActualI => .a11ActualI
  | .a11OldSevenA3 => .a11OldSevenA3
  | .a12Ordinary => .a12Ordinary
  | .a12ActualI => .a12ActualI
  | .a12TypedSeven _ => .a12TypedSeven
  | .a13 => .a13
  | .a14 _ => .a14
  | .a15 => .a15
  | .wholeA16 => .wholeA16
  | .b14TypedSeven _ => .b14TypedSeven
  | .b15 => .b15
  | .wholeB16 => .wholeB16
  | .d13 => .d13
  | .d14 _ => .d14
  | .d15 => .d15
  | .wholeD16 => .wholeD16

def initialGroups : List Group :=
  [.a10ActualI, .a10TypedSeven, .a11Ordinary, .a11ActualI,
   .a11OldSevenA3, .a12Ordinary, .a12ActualI, .a12TypedSeven,
   .a13, .a14, .a15, .wholeA16, .b14TypedSeven, .b15,
   .wholeB16, .d13, .d14, .d15, .wholeD16]

/-- Remove exactly six WHOLE groups; retain A14 with four spectator cases. -/
def survivingGroups : List Group :=
  [.a10ActualI, .a10TypedSeven, .a11ActualI,
   .a12Ordinary, .a12ActualI, .a12TypedSeven, .a13, .a14,
   .a15, .b15, .wholeB16, .d15, .wholeD16]

def initialA14Spectators : List A14Spectator :=
  [.ordinaryT2, .twoOrdinaryA1, .actualIOrdinaryA1,
   .oldA2NewIsolatedSeven, .oldI2SevenNewIsolatedThree]

def survivingA14Spectators : List A14Spectator :=
  [.ordinaryT2, .twoOrdinaryA1,
   .oldA2NewIsolatedSeven, .oldI2SevenNewIsolatedThree]

def initialD14Spectators : List D14Spectator :=
  [.ordinaryDihedral, .oldA2NewIsolatedSeven, .oldI2SevenNewIsolatedThree]

def Survives (family : ForestCase) : Prop :=
  family.group ∈ survivingGroups ∧
    match family with
    | .a14 spectator => spectator ∈ survivingA14Spectators
    | _ => True

/-- Closed finite data: complete original enums and exact19/13,5/4,3 counts. -/
theorem exact_enumerations :
    (∀ family : Group, family ∈ initialGroups) ∧
    (∀ spectator : A14Spectator, spectator ∈ initialA14Spectators) ∧
    (∀ spectator : D14Spectator, spectator ∈ initialD14Spectators) ∧
    initialGroups.length = 19 ∧ initialGroups.Nodup ∧
    survivingGroups.length = 13 ∧ survivingGroups.Nodup ∧
    initialA14Spectators.length = 5 ∧ initialA14Spectators.Nodup ∧
    survivingA14Spectators.length = 4 ∧ survivingA14Spectators.Nodup ∧
    initialD14Spectators.length = 3 ∧ initialD14Spectators.Nodup := by
  constructor
  · intro family
    cases family <;> decide
  · constructor
    · intro spectator
      cases spectator <;> decide
    · constructor
      · intro spectator
        cases spectator <;> decide
      · decide

universe u

/--
Exactly NINE named external written interfaces: the fresh classification,
SIX earlier exclusions and TWO separate D14 spectator exclusions.
The whole D14 group is not an input. No surviving-group conclusion,
geometry instance, dimension bound or cardinality claim is assumed here.
Selected and ActualForest must be instantiated at the ORIGINAL vertex.
-/
structure WrittenInputs
    (Vertex : Type u) (dimension : Nat)
    (Selected : Vertex → Prop)
    (ActualForest : Vertex → ForestCase → Prop) : Prop where
  freshComplete : dimension = 16 → ∀ p, Selected p →
    ∃ family, ActualForest p family
  d13Exclusion : dimension = 16 → ∀ p, Selected p →
    ActualForest p .d13 → False
  a11OrdinaryExclusion : dimension = 16 → ∀ p, Selected p →
    ActualForest p .a11Ordinary → False
  b14TypedSevenExclusion : dimension = 16 → ∀ p meaning, Selected p →
    ActualForest p (.b14TypedSeven meaning) → False
  wholeA16Exclusion : dimension = 16 → ∀ p, Selected p →
    ActualForest p .wholeA16 → False
  a11OldSevenA3Exclusion : dimension = 16 → ∀ p, Selected p →
    ActualForest p .a11OldSevenA3 → False
  a14ActualIOrdinaryA1Exclusion : dimension = 16 → ∀ p, Selected p →
    ActualForest p (.a14 .actualIOrdinaryA1) → False
  d14OrdinaryDihedralExclusion : dimension = 16 → ∀ p, Selected p →
    ActualForest p (.d14 .ordinaryDihedral) → False
  d14TypedSevenExclusion : dimension = 16 → ∀ p meaning, Selected p →
    ActualForest p (.d14 (D14Spectator.ofTypedSeven meaning)) → False

section Refinement

variable {Vertex : Type u} {dimension : Nat}
variable {Selected : Vertex → Prop}
variable {ActualForest : Vertex → ForestCase → Prop}

/-- The removed A14 subcase concerns exactly this selected p and forest. -/
theorem a14_actual_i_ordinary_a1_exclusion
    (written : WrittenInputs Vertex dimension Selected ActualForest)
    (at16 : dimension = 16) (p : Vertex) (selected : Selected p)
    (forest : ActualForest p (.a14 .actualIOrdinaryA1)) : False :=
  written.a14ActualIOrdinaryA1Exclusion at16 p selected forest

/-- Derive the ENTIRE original D14-row removal from all three spectator cases. -/
theorem d14_all_spectators_exclusion
    (written : WrittenInputs Vertex dimension Selected ActualForest)
    (at16 : dimension = 16) (p : Vertex) (selected : Selected p)
    (spectator : D14Spectator) (forest : ActualForest p (.d14 spectator)) : False := by
  cases spectator with
  | ordinaryDihedral =>
    exact written.d14OrdinaryDihedralExclusion at16 p selected forest
  | oldA2NewIsolatedSeven =>
    exact written.d14TypedSevenExclusion at16 p .oldA2NewIsolatedSeven selected forest
  | oldI2SevenNewIsolatedThree =>
    exact written.d14TypedSevenExclusion at16 p .oldI2SevenNewIsolatedThree selected forest

/-- Every classified forest at the SAME selected p lies in the refined domain. -/
theorem classified_forest_refinement
    (written : WrittenInputs Vertex dimension Selected ActualForest)
    (at16 : dimension = 16) (p : Vertex) (selected : Selected p)
    (family : ForestCase) (forest : ActualForest p family) : Survives family := by
  cases family with
  | a10ActualI => simp [Survives, ForestCase.group, survivingGroups, survivingA14Spectators]
  | a10TypedSeven _ => simp [Survives, ForestCase.group, survivingGroups, survivingA14Spectators]
  | a11Ordinary =>
    exact False.elim (written.a11OrdinaryExclusion at16 p selected forest)
  | a11ActualI => simp [Survives, ForestCase.group, survivingGroups, survivingA14Spectators]
  | a11OldSevenA3 =>
    exact False.elim (written.a11OldSevenA3Exclusion at16 p selected forest)
  | a12Ordinary => simp [Survives, ForestCase.group, survivingGroups, survivingA14Spectators]
  | a12ActualI => simp [Survives, ForestCase.group, survivingGroups, survivingA14Spectators]
  | a12TypedSeven _ => simp [Survives, ForestCase.group, survivingGroups, survivingA14Spectators]
  | a13 => simp [Survives, ForestCase.group, survivingGroups, survivingA14Spectators]
  | a14 spectator =>
    cases spectator with
    | ordinaryT2 => simp [Survives, ForestCase.group, survivingGroups, survivingA14Spectators]
    | twoOrdinaryA1 => simp [Survives, ForestCase.group, survivingGroups, survivingA14Spectators]
    | actualIOrdinaryA1 =>
      exact False.elim (a14_actual_i_ordinary_a1_exclusion written at16 p selected forest)
    | oldA2NewIsolatedSeven => simp [Survives, ForestCase.group, survivingGroups, survivingA14Spectators]
    | oldI2SevenNewIsolatedThree => simp [Survives, ForestCase.group, survivingGroups, survivingA14Spectators]
  | a15 => simp [Survives, ForestCase.group, survivingGroups, survivingA14Spectators]
  | wholeA16 =>
    exact False.elim (written.wholeA16Exclusion at16 p selected forest)
  | b14TypedSeven meaning =>
    exact False.elim (written.b14TypedSevenExclusion at16 p meaning selected forest)
  | b15 => simp [Survives, ForestCase.group, survivingGroups, survivingA14Spectators]
  | wholeB16 => simp [Survives, ForestCase.group, survivingGroups, survivingA14Spectators]
  | d13 => exact False.elim (written.d13Exclusion at16 p selected forest)
  | d14 spectator =>
    exact False.elim (d14_all_spectators_exclusion written at16 p selected spectator forest)
  | d15 => simp [Survives, ForestCase.group, survivingGroups, survivingA14Spectators]
  | wholeD16 => simp [Survives, ForestCase.group, survivingGroups, survivingA14Spectators]

/-- The complete fresh classification refines without replacing p or its forest. -/
theorem same_point_complete_refinement
    (written : WrittenInputs Vertex dimension Selected ActualForest)
    (at16 : dimension = 16) (p : Vertex) (selected : Selected p) :
    ∃ family, ActualForest p family ∧ Survives family := by
  obtain ⟨family, forest⟩ := written.freshComplete at16 p selected
  exact ⟨family, forest,
    classified_forest_refinement written at16 p selected family forest⟩

/-- The A14 GROUP remains, with precisely the four permitted spectator labels. -/
theorem a14_spectators_refine_to_four
    (written : WrittenInputs Vertex dimension Selected ActualForest)
    (at16 : dimension = 16) (p : Vertex) (selected : Selected p)
    (spectator : A14Spectator) (forest : ActualForest p (.a14 spectator)) :
    spectator ∈ survivingA14Spectators := by
  exact (classified_forest_refinement written at16 p selected (.a14 spectator) forest).2

end Refinement
end CoxeterProofs.ArithmeticSixteenFrontierRefinement

#print axioms CoxeterProofs.ArithmeticSixteenFrontierRefinement.exact_enumerations
#print axioms CoxeterProofs.ArithmeticSixteenFrontierRefinement.a14_actual_i_ordinary_a1_exclusion
#print axioms CoxeterProofs.ArithmeticSixteenFrontierRefinement.d14_all_spectators_exclusion
#print axioms CoxeterProofs.ArithmeticSixteenFrontierRefinement.classified_forest_refinement
#print axioms CoxeterProofs.ArithmeticSixteenFrontierRefinement.same_point_complete_refinement
#print axioms CoxeterProofs.ArithmeticSixteenFrontierRefinement.a14_spectators_refine_to_four
