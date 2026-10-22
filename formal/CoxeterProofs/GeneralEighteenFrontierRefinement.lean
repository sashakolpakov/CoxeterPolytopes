import Std

/-!
SOURCE PREPARATION ONLY: uncompiled; six requested axiom prints, zero actual outputs.

Conditional refinement of the FRESH eleven-group general compact18 frontier.
Nine written interfaces are EXTERNAL. No geometry instance is supplied.
Whole A15 and A16 removal is derived from their separate spectator interfaces.
Every deduction retains the SAME selected vertex and actual forest witness.
The six surviving groups are necessary cases, not constructed examples.
This module proves no dimension-seventeen bound or geometric energy estimate.
-/

namespace CoxeterProofs.GeneralEighteenFrontierRefinement

/-- Each meaning includes its ACTUAL isolated new opposite wall. -/
inductive TypedSeven where
  | oldA2NewIsolatedSeven
  | oldI2SevenNewIsolatedThree
  deriving DecidableEq, Repr

/-- The two rank-five factors allowed beside actual I in the original A12 row. -/
inductive RankFive where
  | a5
  | d5
  deriving DecidableEq, Repr

/-- Ordinary irreducible rank-four factors in the original typed-seven A12 row. -/
inductive OrdinaryRankFour where
  | a4
  | b4
  | d4
  | f4
  | h4
  deriving DecidableEq, Repr

inductive A12Spectator where
  | actualI (factor : RankFive)
  | typedSeven (meaning : TypedSeven) (factor : OrdinaryRankFour)
  deriving DecidableEq, Repr

/-- T2 means ordinary irreducible I2(m), m >= 3; each unmarked A1 is nonshort. -/
inductive A15OrdinarySpectator where
  | ordinaryA3
  | ordinaryT2OrdinaryA1
  deriving DecidableEq, Repr

inductive A15PaidSpectator where
  | actualIOrdinaryT2
  | oldI2SevenNewIsolatedThreeOrdinaryA1
  deriving DecidableEq, Repr

inductive A15Spectator where
  | ordinary (spectator : A15OrdinarySpectator)
  | paid (spectator : A15PaidSpectator)
  deriving DecidableEq, Repr

inductive A16OrdinarySpectator where
  | ordinaryT2
  | twoOrdinaryA1
  | actualIOrdinaryA1
  deriving DecidableEq, Repr

inductive A16Spectator where
  | ordinary (spectator : A16OrdinarySpectator)
  | typedSeven (meaning : TypedSeven)
  deriving DecidableEq, Repr

/-- Used separately in the A17 and D17 rows; only A17 is narrowed here. -/
inductive RankOneSpectator where
  | ordinaryNonshortA1
  | actualI
  deriving DecidableEq, Repr

/-- The original ELEVEN groups, with no rank decrement or arithmetic premise. -/
inductive Group where
  | a12
  | a15
  | a16
  | a17
  | wholeA18
  | b16TypedSeven
  | b17ActualI
  | wholeB18
  | d16TypedSeven
  | d17
  | wholeD18
  deriving DecidableEq, Repr

/-- Labels for ENTIRE actual forests, not constructed Coxeter diagrams. -/
inductive ForestCase where
  | a12 (spectator : A12Spectator)
  | a15 (spectator : A15Spectator)
  | a16 (spectator : A16Spectator)
  | a17 (spectator : RankOneSpectator)
  | wholeA18
  | b16TypedSeven (meaning : TypedSeven)
  | b17ActualI
  | wholeB18
  | d16TypedSeven (meaning : TypedSeven)
  | d17 (spectator : RankOneSpectator)
  | wholeD18
  deriving DecidableEq, Repr

def ForestCase.group : ForestCase → Group
  | .a12 _ => .a12
  | .a15 _ => .a15
  | .a16 _ => .a16
  | .a17 _ => .a17
  | .wholeA18 => .wholeA18
  | .b16TypedSeven _ => .b16TypedSeven
  | .b17ActualI => .b17ActualI
  | .wholeB18 => .wholeB18
  | .d16TypedSeven _ => .d16TypedSeven
  | .d17 _ => .d17
  | .wholeD18 => .wholeD18

def initialGroups : List Group :=
  [.a12, .a15, .a16, .a17, .wholeA18, .b16TypedSeven,
   .b17ActualI, .wholeB18, .d16TypedSeven, .d17, .wholeD18]

/-- Remove precisely A12, A15, A16, B16 and D16; keep the A17 GROUP. -/
def survivingGroups : List Group :=
  [.a17, .wholeA18, .b17ActualI, .wholeB18, .d17, .wholeD18]

def initialA12Spectators : List A12Spectator :=
  [.actualI .a5, .actualI .d5,
   .typedSeven .oldA2NewIsolatedSeven .a4,
   .typedSeven .oldA2NewIsolatedSeven .b4,
   .typedSeven .oldA2NewIsolatedSeven .d4,
   .typedSeven .oldA2NewIsolatedSeven .f4,
   .typedSeven .oldA2NewIsolatedSeven .h4,
   .typedSeven .oldI2SevenNewIsolatedThree .a4,
   .typedSeven .oldI2SevenNewIsolatedThree .b4,
   .typedSeven .oldI2SevenNewIsolatedThree .d4,
   .typedSeven .oldI2SevenNewIsolatedThree .f4,
   .typedSeven .oldI2SevenNewIsolatedThree .h4]

def initialA15Spectators : List A15Spectator :=
  [.ordinary .ordinaryA3, .ordinary .ordinaryT2OrdinaryA1,
   .paid .actualIOrdinaryT2, .paid .oldI2SevenNewIsolatedThreeOrdinaryA1]

def initialA16Spectators : List A16Spectator :=
  [.ordinary .ordinaryT2, .ordinary .twoOrdinaryA1,
   .ordinary .actualIOrdinaryA1, .typedSeven .oldA2NewIsolatedSeven,
   .typedSeven .oldI2SevenNewIsolatedThree]

def initialRankOneSpectators : List RankOneSpectator :=
  [.ordinaryNonshortA1, .actualI]

def survivingA17Spectators : List RankOneSpectator := [.actualI]

/-- Both inclusion in the original domain and the exact surviving refinement. -/
def Survives (family : ForestCase) : Prop :=
  family.group ∈ initialGroups ∧ family.group ∈ survivingGroups ∧
    match family with
    | .a17 spectator => spectator ∈ survivingA17Spectators
    | _ => True

/-- Closed finite coverage, distinctness and exact11/6,12/4/5,2/1 counts. -/
theorem exact_enumerations :
    (∀ family : Group, family ∈ initialGroups) ∧
    (∀ spectator : A12Spectator, spectator ∈ initialA12Spectators) ∧
    (∀ spectator : A15Spectator, spectator ∈ initialA15Spectators) ∧
    (∀ spectator : A16Spectator, spectator ∈ initialA16Spectators) ∧
    (∀ spectator : RankOneSpectator, spectator ∈ initialRankOneSpectators) ∧
    (initialGroups.length = 11 ∧ initialGroups.Nodup ∧
     survivingGroups.length = 6 ∧ survivingGroups.Nodup ∧
     initialA12Spectators.length = 12 ∧ initialA12Spectators.Nodup ∧
     initialA15Spectators.length = 4 ∧ initialA15Spectators.Nodup ∧
     initialA16Spectators.length = 5 ∧ initialA16Spectators.Nodup ∧
     initialRankOneSpectators.length = 2 ∧ initialRankOneSpectators.Nodup ∧
     survivingA17Spectators.length = 1 ∧ survivingA17Spectators.Nodup) := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩
  · intro family
    cases family <;> decide
  · intro spectator
    cases spectator with
    | actualI factor => cases factor <;> decide
    | typedSeven meaning factor => cases meaning <;> cases factor <;> decide
  · intro spectator
    cases spectator with
    | ordinary ordinary => cases ordinary <;> decide
    | paid paid => cases paid <;> decide
  · intro spectator
    cases spectator with
    | ordinary ordinary => cases ordinary <;> decide
    | typedSeven meaning => cases meaning <;> decide
  · intro spectator
    cases spectator <;> decide
  · decide

universe u

/--
Exactly NINE external written interfaces. The fresh eleven-group list is
an input, but neither surviving membership nor the five whole-row removals
is assumed as a combined conclusion. A15 and A16 are split at their precise
written spectator interfaces. Selected and ActualForest are uninstantiated
predicates on the ORIGINAL vertex of one fixed hypothetical candidate.
-/
structure WrittenInputs
    (Vertex : Type u) (dimension : Nat)
    (Selected : Vertex → Prop)
    (ActualForest : Vertex → ForestCase → Prop) : Prop where
  freshComplete : dimension = 18 → ∀ p, Selected p →
    ∃ family, ActualForest p family
  a12Exclusion : dimension = 18 → ∀ p spectator, Selected p →
    ActualForest p (.a12 spectator) → False
  a15OrdinaryExclusion : dimension = 18 → ∀ p spectator, Selected p →
    ActualForest p (.a15 (.ordinary spectator)) → False
  a15PaidExclusion : dimension = 18 → ∀ p spectator, Selected p →
    ActualForest p (.a15 (.paid spectator)) → False
  a16OrdinaryExclusion : dimension = 18 → ∀ p spectator, Selected p →
    ActualForest p (.a16 (.ordinary spectator)) → False
  a16TypedSevenExclusion : dimension = 18 → ∀ p meaning, Selected p →
    ActualForest p (.a16 (.typedSeven meaning)) → False
  b16TypedSevenExclusion : dimension = 18 → ∀ p meaning, Selected p →
    ActualForest p (.b16TypedSeven meaning) → False
  d16TypedSevenExclusion : dimension = 18 → ∀ p meaning, Selected p →
    ActualForest p (.d16TypedSeven meaning) → False
  a17OrdinaryExclusion : dimension = 18 → ∀ p, Selected p →
    ActualForest p (.a17 .ordinaryNonshortA1) → False

section Refinement

variable {Vertex : Type u} {dimension : Nat}
variable {Selected : Vertex → Prop}
variable {ActualForest : Vertex → ForestCase → Prop}

/-- Derive the entire four-case A15 removal from BOTH separate written inputs. -/
theorem a15_all_spectators_exclusion
    (written : WrittenInputs Vertex dimension Selected ActualForest)
    (at18 : dimension = 18) (p : Vertex) (selected : Selected p)
    (spectator : A15Spectator) (forest : ActualForest p (.a15 spectator)) : False := by
  cases spectator with
  | ordinary ordinary =>
    exact written.a15OrdinaryExclusion at18 p ordinary selected forest
  | paid paid =>
    exact written.a15PaidExclusion at18 p paid selected forest

/-- Derive the entire five-case A16 removal, retaining BOTH actual typed meanings. -/
theorem a16_all_spectators_exclusion
    (written : WrittenInputs Vertex dimension Selected ActualForest)
    (at18 : dimension = 18) (p : Vertex) (selected : Selected p)
    (spectator : A16Spectator) (forest : ActualForest p (.a16 spectator)) : False := by
  cases spectator with
  | ordinary ordinary =>
    exact written.a16OrdinaryExclusion at18 p ordinary selected forest
  | typedSeven meaning =>
    exact written.a16TypedSevenExclusion at18 p meaning selected forest

/-- Every actual classified forest at the SAME selected p has the refined scope. -/
theorem classified_forest_refinement
    (written : WrittenInputs Vertex dimension Selected ActualForest)
    (at18 : dimension = 18) (p : Vertex) (selected : Selected p)
    (family : ForestCase) (forest : ActualForest p family) : Survives family := by
  cases family with
  | a12 spectator =>
    exact False.elim (written.a12Exclusion at18 p spectator selected forest)
  | a15 spectator =>
    exact False.elim (a15_all_spectators_exclusion written at18 p selected spectator forest)
  | a16 spectator =>
    exact False.elim (a16_all_spectators_exclusion written at18 p selected spectator forest)
  | a17 spectator =>
    cases spectator with
    | ordinaryNonshortA1 =>
      exact False.elim (written.a17OrdinaryExclusion at18 p selected forest)
    | actualI =>
      simp [Survives, ForestCase.group, initialGroups, survivingGroups,
        survivingA17Spectators]
  | wholeA18 =>
    simp [Survives, ForestCase.group, initialGroups, survivingGroups,
      survivingA17Spectators]
  | b16TypedSeven meaning =>
    exact False.elim (written.b16TypedSevenExclusion at18 p meaning selected forest)
  | b17ActualI =>
    simp [Survives, ForestCase.group, initialGroups, survivingGroups,
      survivingA17Spectators]
  | wholeB18 =>
    simp [Survives, ForestCase.group, initialGroups, survivingGroups,
      survivingA17Spectators]
  | d16TypedSeven meaning =>
    exact False.elim (written.d16TypedSevenExclusion at18 p meaning selected forest)
  | d17 _ =>
    simp [Survives, ForestCase.group, initialGroups, survivingGroups,
      survivingA17Spectators]
  | wholeD18 =>
    simp [Survives, ForestCase.group, initialGroups, survivingGroups,
      survivingA17Spectators]

/-- Obtain the original witness and refine it without changing p or its forest. -/
theorem same_point_complete_refinement
    (written : WrittenInputs Vertex dimension Selected ActualForest)
    (at18 : dimension = 18) (p : Vertex) (selected : Selected p) :
    ∃ family, ActualForest p family ∧ Survives family := by
  obtain ⟨family, forest⟩ := written.freshComplete at18 p selected
  exact ⟨family, forest,
    classified_forest_refinement written at18 p selected family forest⟩

/-- A17 remains, but ONLY with actual I; the D17 spectator enum is unchanged. -/
theorem a17_spectators_refine_to_actual_i
    (written : WrittenInputs Vertex dimension Selected ActualForest)
    (at18 : dimension = 18) (p : Vertex) (selected : Selected p)
    (spectator : RankOneSpectator) (forest : ActualForest p (.a17 spectator)) :
    spectator = .actualI := by
  have retained : spectator ∈ survivingA17Spectators :=
    (classified_forest_refinement written at18 p selected (.a17 spectator) forest).2.2
  simpa [survivingA17Spectators] using retained

end Refinement
end CoxeterProofs.GeneralEighteenFrontierRefinement

#print axioms CoxeterProofs.GeneralEighteenFrontierRefinement.exact_enumerations
#print axioms CoxeterProofs.GeneralEighteenFrontierRefinement.a15_all_spectators_exclusion
#print axioms CoxeterProofs.GeneralEighteenFrontierRefinement.a16_all_spectators_exclusion
#print axioms CoxeterProofs.GeneralEighteenFrontierRefinement.classified_forest_refinement
#print axioms CoxeterProofs.GeneralEighteenFrontierRefinement.same_point_complete_refinement
#print axioms CoxeterProofs.GeneralEighteenFrontierRefinement.a17_spectators_refine_to_actual_i
