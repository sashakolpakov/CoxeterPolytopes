import CoxeterProofs.CompactHVector
import CoxeterProofs.ManuscriptDataBounds

/-!
Final conditional logic for the reviewed eight-field dimension19 theorem.
Actual Coxeter faces, the complete local consumer proofs, finite-forest
coverage and the earlier universal compact20 theorem are explicit inputs.
No geometry, field classification or script execution is declared proved
by this module. There are no new axioms or admitted proofs.
-/
namespace CoxeterProofs.ArithmeticEightField19

open CompactHVector ManuscriptDataBounds

/-- These mean equality of the actual ground field, not containment. -/
inductive GroundField where
  | sqrtThree | sqrtSix | cyclotomicNine | cyclotomicEleven
  | cyclotomicFifteen | cyclotomicSixteen | cyclotomicTwenty
  | biquadraticTwoThree
  deriving DecidableEq, Repr

def multiplicity20 (i : Nat) : Nat := if i = 10 then 1 else 2

def coefficient20 (i : Nat) : Nat :=
  if i = 10 then 35
  else 5 * (i * (i - 1) / 2 + (20 - i) * (19 - i) / 2) - 380

theorem paired_coefficients20 :
    (List.range 11).map coefficient20 =
      [570, 475, 390, 315, 250, 195, 150, 115, 90, 75, 35] := by decide

theorem checked_coefficients20 :
    (List.range 11).all (fun i =>
      decide (35 * multiplicity20 i ≤ coefficient20 i)) = true := by decide

theorem coefficient_lower (i : Nat) (hi : i < 11) :
    35 * multiplicity20 i ≤ coefficient20 i := by
  exact of_decide_eq_true
    ((List.all_eq_true.mp checked_coefficients20) i (List.mem_range.mpr hi))

def endpoint20 : Term := { multiplicity := 2, coefficient := 570, mass := 1 }

theorem endpoint_surplus20 (tail : List Term)
    (coverage : ∀ t ∈ tail, ∃ i < 11,
      t.multiplicity = multiplicity20 i ∧ t.coefficient = coefficient20 i) :
    35 * vertices (endpoint20 :: tail) + 500 ≤
      weighted (endpoint20 :: tail) := by
  have coefficients : ∀ t ∈ tail, 35 * t.multiplicity ≤ t.coefficient := by
    intro t ht
    obtain ⟨i, hi, hm, hc⟩ := coverage t ht
    rw [hm, hc]
    exact coefficient_lower i hi
  have rest := weighted_lower 35 tail coefficients
  simp only [vertices, weighted, List.map_cons, List.sum_cons, endpoint20,
    Nat.mul_one] at *
  omega

/-- Quarter units clear every earlier quarter rule and the new halves.
The complete prior forest theorem and actual ten-family selection supply
oldAllowance. The two distinct changed marks remove four quarter units. -/
structure VertexEntry where
  target : Bool
  oldCost : Nat
  oldAllowance : if target then oldCost = 144 else oldCost ≤ 140

def revisedCost (v : VertexEntry) : Nat :=
  if v.target then v.oldCost - 4 else v.oldCost

theorem revised_vertex_upper (v : VertexEntry) : revisedCost v ≤ 140 := by
  cases h : v.target
  · have hold : v.oldCost ≤ 140 := by simpa [h] using v.oldAllowance
    simpa [revisedCost, h] using hold
  · have hold : v.oldCost = 144 := by simpa [h] using v.oldAllowance
    simp [revisedCost, h, hold]

/-- Pure payment logic. The geometric triangle theorem supplies the first
premise; the actual six-row exclusion and role separation supply the
second. This does not prove those geometric premises. -/
theorem triangle_payment (outside changed : Nat) (geometry : 2 ≤ outside) :
    8 ≤ 4 * outside + 2 * changed := by omega

theorem quadrilateral_payment (outside changed : Nat)
    (geometry : 2 ≤ changed ∨ 1 ≤ outside) :
    4 ≤ 4 * outside + 2 * changed := by omega

/-- Written actual-face and h-vector interfaces. The list contains EVERY
finite vertex. IncidenceLower is supplied by actual polygon payments,
their complete double count, and the stated nonnegative symmetric h-vector
identities. It is not the desired dimension bound. -/
structure Rank20Data where
  actualField : GroundField
  tail : List Term
  entries : List VertexEntry
  coefficientCoverage : ∀ t ∈ tail, ∃ i < 11,
    t.multiplicity = multiplicity20 i ∧ t.coefficient = coefficient20 i
  completeVertexCount : entries.length = vertices (endpoint20 :: tail)
  incidenceLower :
    4 * weighted (endpoint20 :: tail) ≤ (entries.map revisedCost).sum

theorem no_rank20_data (d : Rank20Data) : False := by
  have upper := sum_le_cap_mul_length (d.entries.map revisedCost) 140 (by
    intro a ha
    obtain ⟨v, _, rfl⟩ := List.mem_map.mp ha
    exact revised_vertex_upper v)
  rw [List.length_map, d.completeVertexCount] at upper
  have strict := Nat.mul_le_mul_left 4
    (endpoint_surplus20 d.tail d.coefficientCoverage)
  have lower := d.incidenceLower
  omega

/-- The earlier universal20 result supplies the higher tail. The new
geometric bridge is required separately for EACH actual listed field. -/
theorem dimension_at_most19 (n : Nat) (actualField : GroundField)
    (earlierUniversal20 : n ≤ 20)
    (geometry : n = 20 → ∃ d : Rank20Data, d.actualField = actualField) : n ≤ 19 := by
  have excluded : n ≠ 20 := by
    intro h
    obtain ⟨d, _⟩ := geometry h
    exact no_rank20_data d
  omega

end CoxeterProofs.ArithmeticEightField19

#print axioms CoxeterProofs.ArithmeticEightField19.paired_coefficients20
#print axioms CoxeterProofs.ArithmeticEightField19.endpoint_surplus20
#print axioms CoxeterProofs.ArithmeticEightField19.revised_vertex_upper
#print axioms CoxeterProofs.ArithmeticEightField19.triangle_payment
#print axioms CoxeterProofs.ArithmeticEightField19.quadrilateral_payment
#print axioms CoxeterProofs.ArithmeticEightField19.no_rank20_data
#print axioms CoxeterProofs.ArithmeticEightField19.dimension_at_most19
