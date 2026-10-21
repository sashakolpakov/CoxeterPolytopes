import Std

/-! Small finite-sum lemmas used by conditional manuscript interfaces.
No geometric or enumeration result is declared as an axiom. -/
namespace CoxeterProofs.ManuscriptDataBounds

theorem sum_le_cap_mul_length (xs : List Nat) (cap : Nat)
    (bound : ∀ x ∈ xs, x ≤ cap) : xs.sum ≤ cap * xs.length := by
  induction xs with
  | nil => simp
  | cons x xs ih =>
    have head := bound x (by simp)
    have tail := ih (by intro y hy; exact bound y (by simp [hy]))
    simp only [List.sum_cons, List.length_cons, Nat.mul_add, Nat.mul_one]
    omega

structure Charge where
  cost : Nat
  curvature : Nat

def costs (xs : List Charge) : Nat := (xs.map Charge.cost).sum
def curvatures (xs : List Charge) : Nat := (xs.map Charge.curvature).sum

theorem sum_charge_gaps (xs : List Charge) (gap : Nat)
    (bound : ∀ x ∈ xs, x.cost + gap ≤ x.curvature) :
    costs xs + gap * xs.length ≤ curvatures xs := by
  induction xs with
  | nil => simp [costs, curvatures]
  | cons x xs ih =>
    have head := bound x (by simp)
    have tail := ih (by intro y hy; exact bound y (by simp [hy]))
    simp only [costs, curvatures, List.map_cons, List.sum_cons,
      List.length_cons, Nat.mul_add, Nat.mul_one] at *
    omega

theorem sum_nonpositive (xs : List Int)
    (bound : ∀ x ∈ xs, x ≤ 0) : xs.sum ≤ 0 := by
  induction xs with
  | nil => simp
  | cons x xs ih =>
    have head := bound x (by simp)
    have tail := ih (by intro y hy; exact bound y (by simp [hy]))
    simp only [List.sum_cons]
    omega

end CoxeterProofs.ManuscriptDataBounds
