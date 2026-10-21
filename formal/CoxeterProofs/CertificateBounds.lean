import Std

/-!
Kernel-checked arithmetic used by additive finite-forest certificates.
This module proves the induction, rather than assuming the resulting upper
bound. The finite table and its item comparisons must still be supplied.
-/
namespace CoxeterProofs

structure ChargedItem where
  rank : Nat
  charge : Nat
  deriving DecidableEq, Repr

def totalRank (xs : List ChargedItem) : Nat := (xs.map ChargedItem.rank).sum
def totalCharge (xs : List ChargedItem) : Nat := (xs.map ChargedItem.charge).sum

def tableValue (table : List Nat) (n : Nat) : Nat := table[n]?.getD 0

def checkLastItems (cap : Nat) (allowed : List ChargedItem)
    (table : List Nat) : Bool :=
  (List.range (cap + 1)).all fun n =>
    allowed.all fun a =>
      if n + a.rank ≤ cap then
        decide (tableValue table n + a.charge ≤ tableValue table (n + a.rank))
      else true

theorem last_item_steps_of_check
    (cap : Nat) (allowed : List ChargedItem) (table : List Nat)
    (checked : checkLastItems cap allowed table = true) :
    ∀ a ∈ allowed, ∀ n, n + a.rank ≤ cap →
      tableValue table n + a.charge ≤ tableValue table (n + a.rank) := by
  intro a ha n hn
  unfold checkLastItems at checked
  have row := (List.all_eq_true.mp checked) n (List.mem_range.mpr (by omega))
  have cell := (List.all_eq_true.mp row) a ha
  simp only [hn, ↓reduceIte, decide_eq_true_eq] at cell
  exact cell

theorem additive_upper_of_last_item
    (cap : Nat) (allowed : List ChargedItem) (table : Nat → Nat)
    (zero : table 0 = 0)
    (steps : ∀ a ∈ allowed, ∀ n, n + a.rank ≤ cap →
      table n + a.charge ≤ table (n + a.rank))
    (xs : List ChargedItem)
    (membership : ∀ a ∈ xs, a ∈ allowed)
    (capacity : totalRank xs ≤ cap) :
    totalCharge xs ≤ table (totalRank xs) := by
  induction xs with
  | nil => simp [totalRank, totalCharge, zero]
  | cons a xs ih =>
    have tailMembership : ∀ b ∈ xs, b ∈ allowed := by
      intro b hb
      exact membership b (List.mem_cons_of_mem a hb)
    have tailCapacity : totalRank xs ≤ cap := by
      simp only [totalRank, List.map_cons, List.sum_cons] at capacity ⊢
      omega
    have tail := ih tailMembership tailCapacity
    have step := steps a (membership a (by simp)) (totalRank xs) (by
      simp only [totalRank, List.map_cons, List.sum_cons] at capacity ⊢
      omega)
    rw [Nat.add_comm (totalRank xs) a.rank] at step
    simp only [totalRank, totalCharge, List.map_cons, List.sum_cons] at *
    omega

theorem sum_lower_of_minimum (m : Nat) (xs : List Nat)
    (minimum : ∀ a ∈ xs, m ≤ a) : m * xs.length ≤ xs.sum := by
  induction xs with
  | nil => simp
  | cons a xs ih =>
    have ha := minimum a (by simp)
    have ht := ih (by
      intro b hb
      exact minimum b (List.mem_cons_of_mem a hb))
    simp only [List.length_cons, List.sum_cons, Nat.mul_add, Nat.mul_one]
    omega

/- A selected/rest partition is what the geometric injections must construct.
   No geometric injection is being postulated as a published theorem here. -/
theorem selected_defect_budget
    (W m q s : Nat) (selected rest : List Nat)
    (partition : selected.sum + rest.sum = 2 * W)
    (count : selected.length = 1 + q + s)
    (minimum : ∀ a ∈ selected, m ≤ a) :
    m * (1 + q + s) ≤ 2 * W := by
  have lower := sum_lower_of_minimum m selected minimum
  rw [count] at lower
  omega

end CoxeterProofs
