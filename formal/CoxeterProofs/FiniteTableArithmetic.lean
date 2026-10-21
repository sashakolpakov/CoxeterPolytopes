import CoxeterProofs.CertificateBounds

namespace CoxeterProofs.FiniteTable

structure Table where
  cap : Nat
  allowed : List ChargedItem
  values : List Nat
  deriving DecidableEq, Repr

def check (t : Table) : Bool :=
  decide (t.values.length = t.cap+1) &&
  decide (tableValue t.values 0 = 0) &&
  checkLastItems t.cap t.allowed t.values

theorem upper_of_check (t : Table) (checked : check t = true)
    (xs : List ChargedItem)
    (membership : ∀ a ∈ xs, a ∈ t.allowed)
    (capacity : totalRank xs ≤ t.cap) :
    totalCharge xs ≤ tableValue t.values (totalRank xs) := by
  simp only [check, Bool.and_eq_true, decide_eq_true_eq] at checked
  exact additive_upper_of_last_item t.cap t.allowed (tableValue t.values)
    checked.1.2 (last_item_steps_of_check _ _ _ checked.2) xs membership capacity

end CoxeterProofs.FiniteTable
