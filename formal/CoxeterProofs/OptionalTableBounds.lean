import CoxeterProofs.CertificateBounds

/-! Bellman upper bounds with explicit unreachable ranks. No numeric
sentinel or zero-filled unreachable entry is used. -/
namespace CoxeterProofs.OptionalTable
open CoxeterProofs

def lookup (values : List (Option Nat)) (n : Nat) : Option Nat :=
  values[n]?.getD none

def stepCheck (values : List (Option Nat)) (n : Nat) (a : ChargedItem) : Bool :=
  match lookup values n with
  | none => true
  | some w => match lookup values (n+a.rank) with
    | none => false
    | some z => decide (w+a.charge ≤ z)

def check (cap : Nat) (allowed : List ChargedItem) (values : List (Option Nat)) : Bool :=
  decide (lookup values 0 = some 0) &&
  (List.range (cap+1)).all (fun n => allowed.all (fun a =>
    if n+a.rank ≤ cap then stepCheck values n a else true))

structure Table where
  cap : Nat
  allowed : List ChargedItem
  values : List (Option Nat)
  deriving DecidableEq, Repr

def tableCheck (t : Table) : Bool := check t.cap t.allowed t.values

theorem step_of_check (cap : Nat) (allowed : List ChargedItem)
    (values : List (Option Nat)) (checked : check cap allowed values = true)
    (a : ChargedItem) (member : a ∈ allowed) (n w : Nat)
    (capacity : n+a.rank ≤ cap) (previous : lookup values n = some w) :
    ∃ z, lookup values (n+a.rank) = some z ∧ w+a.charge ≤ z := by
  simp only [check, Bool.and_eq_true, decide_eq_true_eq] at checked
  have row := (List.all_eq_true.mp checked.2) n (List.mem_range.mpr (by omega))
  have cell := (List.all_eq_true.mp row) a member
  simp only [capacity, ↓reduceIte, stepCheck, previous] at cell
  cases result : lookup values (n+a.rank) with
  | none => simp [result] at cell
  | some z => exact ⟨z, rfl, of_decide_eq_true (by simpa [result] using cell)⟩

theorem forest_upper (cap : Nat) (allowed : List ChargedItem)
    (values : List (Option Nat)) (checked : check cap allowed values = true)
    (xs : List ChargedItem) (membership : ∀ a ∈ xs, a ∈ allowed)
    (capacity : totalRank xs ≤ cap) :
    ∃ z, lookup values (totalRank xs) = some z ∧ totalCharge xs ≤ z := by
  induction xs with
  | nil =>
    have zero : lookup values 0 = some 0 := by
      simp only [check, Bool.and_eq_true, decide_eq_true_eq] at checked
      exact checked.1
    exact ⟨0, by simpa [totalRank] using zero, by simp [totalCharge]⟩
  | cons a xs ih =>
    have members : ∀ b ∈ xs, b ∈ allowed := by
      intro b hb; exact membership b (by simp [hb])
    have tailCapacity : totalRank xs ≤ cap := by
      simp only [totalRank, List.map_cons, List.sum_cons] at capacity ⊢
      omega
    obtain ⟨w, hw, bound⟩ := ih members tailCapacity
    have nextCapacity : totalRank xs+a.rank ≤ cap := by
      simp only [totalRank, List.map_cons, List.sum_cons] at capacity ⊢
      omega
    obtain ⟨z, hz, upper⟩ := step_of_check cap allowed values checked a
      (membership a (by simp)) (totalRank xs) w nextCapacity hw
    refine ⟨z, ?_, ?_⟩
    · simpa [totalRank, Nat.add_comm] using hz
    · simp only [totalCharge, List.map_cons, List.sum_cons] at *
      omega

theorem unreachable_excludes_forest (cap : Nat) (allowed : List ChargedItem)
    (values : List (Option Nat)) (checked : check cap allowed values = true)
    (xs : List ChargedItem) (membership : ∀ a ∈ xs, a ∈ allowed)
    (capacity : totalRank xs ≤ cap)
    (unreachable : lookup values (totalRank xs) = none) : False := by
  obtain ⟨z, hz, _⟩ := forest_upper cap allowed values checked xs membership capacity
  simp [unreachable] at hz

end CoxeterProofs.OptionalTable
