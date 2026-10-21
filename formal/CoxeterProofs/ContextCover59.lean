import Std

/-! Presence guards do not count witness multiplicity. Fixing all A2
components and reserving at most two A1 components preserves both coupled
context charges. This is the combinatorial cover, not a geometric lemma. -/
namespace CoxeterProofs.ContextCover59

def charge (oldWitness newWitness : Bool) (a1 a2 k5 : Nat) : Nat :=
  (if oldWitness && 2 ≤ a1 then 4*a2*(a2-1) else 0) +
  (if newWitness && 1 ≤ a1 then 16*a2*k5 else 0)

def reserve (oldWitness newWitness : Bool) (a1 : Nat) : Nat :=
  if oldWitness && 2 ≤ a1 then 2
  else if newWitness && 1 ≤ a1 then 1 else 0

theorem reserve_capacity (oldWitness newWitness : Bool) (a1 : Nat) :
    reserve oldWitness newWitness a1 ≤ a1 ∧
    reserve oldWitness newWitness a1 ≤ 2 := by
  simp only [reserve]
  split <;> rename_i h
  · simp only [Bool.and_eq_true, decide_eq_true_eq] at h
    omega
  · split <;> rename_i g
    · simp only [Bool.and_eq_true, decide_eq_true_eq] at g
      omega
    · omega

theorem charge_preserved (oldWitness newWitness : Bool) (a1 a2 k5 : Nat) :
    charge oldWitness newWitness (reserve oldWitness newWitness a1) a2 k5 =
      charge oldWitness newWitness a1 a2 k5 := by
  by_cases h2 : 2 ≤ a1
  · have h1 : 1 ≤ a1 := by omega
    cases oldWitness <;> cases newWitness <;> simp [reserve, charge, h1, h2]
  · by_cases h1 : 1 ≤ a1
    · cases oldWitness <;> cases newWitness <;> simp [reserve, charge, h1, h2]
    · cases oldWitness <;> cases newWitness <;> simp [reserve, charge, h1, h2]

theorem zero_branch_has_zero_charge (oldWitness newWitness : Bool)
    (a1 a2 k5 : Nat) (zero : reserve oldWitness newWitness a1 = 0) :
    charge oldWitness newWitness a1 a2 k5 = 0 := by
  have same := charge_preserved oldWitness newWitness a1 a2 k5
  rw [zero] at same
  simpa [charge] using same.symm

theorem complete_three_branch_cover (oldWitness newWitness : Bool)
    (a1 a2 k5 : Nat) :
    ∃ r : Nat, (r = 0 ∨ r = 1 ∨ r = 2) ∧ r ≤ a1 ∧
      charge oldWitness newWitness r a2 k5 = charge oldWitness newWitness a1 a2 k5 ∧
      (r = 0 → charge oldWitness newWitness a1 a2 k5 = 0) := by
  let r := reserve oldWitness newWitness a1
  have cap := reserve_capacity oldWitness newWitness a1
  refine ⟨r, by omega, cap.1, charge_preserved _ _ _ _ _, ?_⟩
  exact zero_branch_has_zero_charge _ _ _ _ _

end CoxeterProofs.ContextCover59
