import Std

namespace CoxeterProofs.ShapeCurvature

theorem weighted_minimum (indices : List Nat) (weight h : Nat → Nat)
    (positive : ∀ j ∈ indices, 1 ≤ h j) :
    (indices.map weight).sum ≤ (indices.map fun j => weight j * h j).sum := by
  induction indices with
  | nil => simp
  | cons j js ih =>
    have first := Nat.mul_le_mul_left (weight j) (positive j (by simp))
    have rest := ih (by intro k hk; exact positive k (by simp [hk]))
    simp only [List.map_cons, List.sum_cons, Nat.mul_one] at *
    omega

def extremalIndices : List Nat := (List.range 29).map (·+30)
def extremalWeight (j : Nat) : Nat := 2*j-58
def slack (h : Nat → Nat) : Nat :=
  (extremalIndices.map fun j => extremalWeight j * h j).sum

theorem extremal_coefficient_sum :
    (extremalIndices.map extremalWeight).sum = 870 := by decide

theorem extremal_slack_at_least870 (h : Nat → Nat)
    (positive : ∀ j, j ≤ 58 → 1 ≤ h j) : 870 ≤ slack h := by
  rw [← extremal_coefficient_sum]
  apply weighted_minimum
  intro j member
  obtain ⟨k, hk, rfl⟩ := List.mem_map.mp member
  exact positive (k+30) (by have := List.mem_range.mp hk; omega)

theorem two_extrema_scaled_correction (minimum maximum : Nat)
    (low : 870 ≤ minimum) (high : 870 ≤ maximum) :
    8700 ≤ 5*(minimum+maximum) := by omega

theorem finite_index_sum :
    ((List.range 59).map fun j => min j (58-j)).sum = 841 := by decide

theorem finite_curvature_identity : 5*841 = 58*59+783 := by decide

theorem total_payment_contradiction (finite ideal lower upper surplus : Nat)
    (payment : lower ≤ upper)
    (shape : 783*finite + 8700 + surplus ≤ 2*lower)
    (budget : 2*upper + 91*ideal ≤ 783*finite + surplus) : False := by omega

end CoxeterProofs.ShapeCurvature
