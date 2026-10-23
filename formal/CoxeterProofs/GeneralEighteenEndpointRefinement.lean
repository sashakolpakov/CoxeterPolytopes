import CoxeterProofs.GeneralEighteenFrontierRefinement

/-!
SOURCE PREPARATION: no compile or axiom-audit result is asserted here.

Two separately written exclusions refine the earlier eleven-to-six kernel
to FIVE actual forest cases at the SAME selected vertex. Its nine written
inputs remain external; the two additional geometric exclusions are explicit.

The second section proves literal integer certificates for the new A18
paired-finite endpoint argument. Machin's identity, alternating-series
bounds, actual-wall reduction and monotonicity remain written inputs:
no real cosine model or full geometric formalization is supplied.
-/

namespace CoxeterProofs.GeneralEighteenEndpointRefinement

open CoxeterProofs.GeneralEighteenFrontierRefinement (ForestCase RankOneSpectator)

def fiveCases : List ForestCase :=
  [.wholeA18, .b17ActualI, .wholeB18, .d17 .actualI, .wholeD18]

theorem exact_five_cases : fiveCases.length = 5 ∧ fiveCases.Nodup := by decide

universe u

/-- Two NEW external written geometric statements; no model is instantiated. -/
structure AdditionalInputs
    (Vertex : Type u) (dimension : Nat)
    (Selected : Vertex → Prop)
    (ActualForest : Vertex → ForestCase → Prop) : Prop where
  a17ActualIExclusion : dimension = 18 → ∀ p, Selected p →
    ActualForest p (.a17 .actualI) → False
  d17OrdinaryExclusion : dimension = 18 → ∀ p, Selected p →
    ActualForest p (.d17 .ordinaryNonshortA1) → False

section Forests

variable {Vertex : Type u} {dimension : Nat}
variable {Selected : Vertex → Prop}
variable {ActualForest : Vertex → ForestCase → Prop}

theorem a17_all_spectators_exclusion
    (written : GeneralEighteenFrontierRefinement.WrittenInputs Vertex dimension Selected ActualForest)
    (additional : AdditionalInputs Vertex dimension Selected ActualForest)
    (at18 : dimension = 18) (p : Vertex) (selected : Selected p)
    (spectator : RankOneSpectator) (forest : ActualForest p (.a17 spectator)) : False := by
  cases spectator with
  | ordinaryNonshortA1 => exact written.a17OrdinaryExclusion at18 p selected forest
  | actualI => exact additional.a17ActualIExclusion at18 p selected forest

theorem classified_five_case_refinement
    (written : GeneralEighteenFrontierRefinement.WrittenInputs Vertex dimension Selected ActualForest)
    (additional : AdditionalInputs Vertex dimension Selected ActualForest)
    (at18 : dimension = 18) (p : Vertex) (selected : Selected p)
    (family : ForestCase) (forest : ActualForest p family) : family ∈ fiveCases := by
  have retained := GeneralEighteenFrontierRefinement.classified_forest_refinement written at18 p selected family forest
  cases family with
  | a17 spectator =>
    exact False.elim (a17_all_spectators_exclusion written additional at18 p selected spectator forest)
  | d17 spectator =>
    cases spectator with
    | ordinaryNonshortA1 =>
      exact False.elim (additional.d17OrdinaryExclusion at18 p selected forest)
    | actualI => decide
  | wholeA18 => decide
  | b17ActualI => decide
  | wholeB18 => decide
  | wholeD18 => decide
  | a12 spectator =>
    have impossible : False := by
      simpa [GeneralEighteenFrontierRefinement.Survives,
        GeneralEighteenFrontierRefinement.ForestCase.group,
        GeneralEighteenFrontierRefinement.initialGroups,
        GeneralEighteenFrontierRefinement.survivingGroups] using retained
    exact impossible.elim
  | a15 spectator =>
    have impossible : False := by
      simpa [GeneralEighteenFrontierRefinement.Survives,
        GeneralEighteenFrontierRefinement.ForestCase.group,
        GeneralEighteenFrontierRefinement.initialGroups,
        GeneralEighteenFrontierRefinement.survivingGroups] using retained
    exact impossible.elim
  | a16 spectator =>
    have impossible : False := by
      simpa [GeneralEighteenFrontierRefinement.Survives,
        GeneralEighteenFrontierRefinement.ForestCase.group,
        GeneralEighteenFrontierRefinement.initialGroups,
        GeneralEighteenFrontierRefinement.survivingGroups] using retained
    exact impossible.elim
  | b16TypedSeven meaning =>
    have impossible : False := by
      simpa [GeneralEighteenFrontierRefinement.Survives,
        GeneralEighteenFrontierRefinement.ForestCase.group,
        GeneralEighteenFrontierRefinement.initialGroups,
        GeneralEighteenFrontierRefinement.survivingGroups] using retained
    exact impossible.elim
  | d16TypedSeven meaning =>
    have impossible : False := by
      simpa [GeneralEighteenFrontierRefinement.Survives,
        GeneralEighteenFrontierRefinement.ForestCase.group,
        GeneralEighteenFrontierRefinement.initialGroups,
        GeneralEighteenFrontierRefinement.survivingGroups] using retained
    exact impossible.elim

/-- Preserve both the original vertex and the actual forest witness. -/
theorem same_vertex_five_case_refinement
    (written : GeneralEighteenFrontierRefinement.WrittenInputs Vertex dimension Selected ActualForest)
    (additional : AdditionalInputs Vertex dimension Selected ActualForest)
    (at18 : dimension = 18) (p : Vertex) (selected : Selected p) :
    ∃ family, ActualForest p family ∧ family ∈ fiveCases := by
  obtain ⟨family, forest⟩ := written.freshComplete at18 p selected
  exact ⟨family, forest,
    classified_five_case_refinement written additional at18 p selected family forest⟩

end Forests

namespace Arithmetic

def scale : Int := 1000000

inductive AngleIndex where
  | j7 | j8 | j9 | j10 | j11 | j12 | j13 | j14
  | j15 | j16 | j17 | j20 | j21 | j31 | j32
  deriving DecidableEq, Repr

/-- Denominator one million; lower and upper are STRICT cosine bounds. -/
def angleRow : AngleIndex → Nat × Int × Int
  | .j7 => (7, 900968, 900969)
  | .j8 => (8, 923879, 923880)
  | .j9 => (9, 939692, 939693)
  | .j10 => (10, 951056, 951057)
  | .j11 => (11, 959492, 959493)
  | .j12 => (12, 965925, 965926)
  | .j13 => (13, 970941, 970942)
  | .j14 => (14, 974927, 974928)
  | .j15 => (15, 978147, 978148)
  | .j16 => (16, 980785, 980786)
  | .j17 => (17, 982973, 982974)
  | .j20 => (20, 987688, 987689)
  | .j21 => (21, 988830, 988831)
  | .j31 => (31, 994869, 994870)
  | .j32 => (32, 995184, 995185)

def q10 (a p : Int) : Int :=
  3628800*a^10 - 1814400*p^2*a^8 + 151200*p^4*a^6
    - 5040*p^6*a^4 + 90*p^8*a^2 - p^10

def angleSeal (i : AngleIndex) : Prop :=
  let row := angleRow i
  let a : Int := 10000000 * row.1
  let pLower : Int := 31415926
  let pUpper : Int := 31415927
  0 < scale*q10 a pUpper - 3628800*row.2.1*a^10 ∧
  0 < 3628800*row.2.2*a^10 - scale*(q10 a pLower + pLower^10)

instance (i : AngleIndex) : Decidable (angleSeal i) := by
  unfold angleSeal
  infer_instance

/-- Thirty literal Taylor-table inequalities, proved by kernel reduction. -/
theorem all_angle_integer_seals (i : AngleIndex) : angleSeal i := by
  cases i <;> decide

/-- Unnormalized exact fractions; positivity of the final denominators is checked. -/
def addFraction (x y : Int × Int) : Int × Int :=
  (x.1*y.2+y.1*x.2, x.2*y.2)

def scaleFraction (k : Int) (x : Int × Int) : Int × Int := (k*x.1, x.2)

def atanPartial (base last : Nat) : Int × Int :=
  (List.range (last+1)).foldl (fun acc r =>
    let sign : Int := if r % 2 = 0 then 1 else -1
    addFraction acc (sign, ((2*r+1)*base^(2*r+1) : Nat))) (0, 1)

def machinLower : Int × Int :=
  addFraction (scaleFraction 16 (atanPartial 5 5))
    (scaleFraction (-4) (atanPartial 239 2))

def machinUpper : Int × Int :=
  addFraction (scaleFraction 16 (atanPartial 5 6))
    (scaleFraction (-4) (atanPartial 239 1))

theorem machin_integer_seals :
    0 < machinLower.2 ∧ 0 < machinUpper.2 ∧
    31415926*machinLower.2 < 10000000*machinLower.1 ∧
    10000000*machinUpper.1 < 31415927*machinUpper.2 := by decide

theorem radical_integer_seals :
    (1414213 : Int)^2 < 2*scale^2 ∧
    2*scale^2 < (1414214 : Int)^2 := by decide

/-- scale^4 times H(g/scale,h/scale) in the written proof. -/
def hNumerator (g h : Int) : Int :=
  272*g^2*h^2 - 144*(g^2+h^2)*scale^2 - 8*g*h*scale^2 + 57*scale^4

/-- scale^4 times F_(s/scale)(g/scale,h/scale). -/
def fNumerator (g h s : Int) : Int :=
  136*g^2*h^2 - 72*(g^2+h^2)*scale^2 + 19*scale^4 - 4*s*g*h*scale

theorem mutual_three_integer_signs :
    hNumerator 900969 923880 < 0 ∧
    0 < hNumerator 900968 939692 := by decide

theorem mutual_four_endrange_signs :
    fNumerator 951057 scale 1414213 < 0 ∧
    (19 : Int)*11025^2*(140+99) < 34*10783^2*140 := by decide

inductive GapIndex where
  | m11 | m12 | m13 | m14
  deriving DecidableEq, Repr

def gapRows : GapIndex → AngleIndex × AngleIndex × AngleIndex
  | .m11 => (.j11, .j31, .j32)
  | .m12 => (.j12, .j20, .j21)
  | .m13 => (.j13, .j16, .j17)
  | .m14 => (.j14, .j15, .j16)

def gapSeal (i : GapIndex) : Prop :=
  let row := gapRows i
  let m := angleRow row.1
  let below := angleRow row.2.1
  let above := angleRow row.2.2
  above.1 = below.1+1 ∧
  fNumerator m.2.2 below.2.2 1414213 < 0 ∧
  0 < fNumerator m.2.1 above.2.1 1414214

instance (i : GapIndex) : Decidable (gapSeal i) := by
  unfold gapSeal
  infer_instance

theorem all_mutual_four_gap_seals (i : GapIndex) : gapSeal i := by
  cases i <;> decide

/-- Integer discreteness closes a gap once both monotone half-lines are excluded. -/
theorem consecutive_gap_excludes_root (Root : Nat → Prop) (lastBelow : Nat)
    (below : ∀ n, n ≤ lastBelow → ¬ Root n)
    (above : ∀ n, lastBelow+1 ≤ n → ¬ Root n) (n : Nat) : ¬ Root n := by
  by_cases earlier : n ≤ lastBelow
  · exact below n earlier
  · exact above n (by omega)

end Arithmetic
end CoxeterProofs.GeneralEighteenEndpointRefinement

#print axioms CoxeterProofs.GeneralEighteenEndpointRefinement.exact_five_cases
#print axioms CoxeterProofs.GeneralEighteenEndpointRefinement.a17_all_spectators_exclusion
#print axioms CoxeterProofs.GeneralEighteenEndpointRefinement.classified_five_case_refinement
#print axioms CoxeterProofs.GeneralEighteenEndpointRefinement.same_vertex_five_case_refinement
#print axioms CoxeterProofs.GeneralEighteenEndpointRefinement.Arithmetic.all_angle_integer_seals
#print axioms CoxeterProofs.GeneralEighteenEndpointRefinement.Arithmetic.machin_integer_seals
#print axioms CoxeterProofs.GeneralEighteenEndpointRefinement.Arithmetic.radical_integer_seals
#print axioms CoxeterProofs.GeneralEighteenEndpointRefinement.Arithmetic.mutual_three_integer_signs
#print axioms CoxeterProofs.GeneralEighteenEndpointRefinement.Arithmetic.mutual_four_endrange_signs
#print axioms CoxeterProofs.GeneralEighteenEndpointRefinement.Arithmetic.all_mutual_four_gap_seals
#print axioms CoxeterProofs.GeneralEighteenEndpointRefinement.Arithmetic.consecutive_gap_excludes_root
