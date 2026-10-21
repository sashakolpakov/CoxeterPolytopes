import Std

/-! All-rank finite-component deficit argument for the sqrt-two arithmetic
rule. Interpreting these explicit costs as actual corner weights and proving
all face payments are separate geometric obligations. -/
namespace CoxeterProofs.CompactFieldBudget

inductive FiniteType where
  | A (rank : Nat)
  | B (rank : Nat)
  | D (rank : Nat)
  | E6 | E7 | E8 | F4 | I2
  deriving DecidableEq, Repr

def rank : FiniteType → Nat
  | .A n | .B n | .D n => n
  | .E6 => 6 | .E7 => 7 | .E8 => 8 | .F4 => 4 | .I2 => 2

def valid : FiniteType → Prop
  | .A n => 1 ≤ n
  | .B n => 2 ≤ n
  | .D n => 4 ≤ n
  | _ => True

def charge : FiniteType → Nat
  | .A n => if n = 1 then 0 else 8*n-12
  | .B n => 8*n-12
  | .D n => if n = 5 then 36 else 8*n-8
  | .E6 => 44 | .E7 => 49 | .E8 => 57 | .F4 => 24 | .I2 => 4

theorem component_deficit (t : FiniteType) (allowed : valid t) :
    charge t+4 ≤ 8*rank t := by
  cases t <;> simp only [valid, rank, charge] at *
  all_goals first | (split <;> omega) | omega

theorem large_component_deficit (t : FiniteType) (allowed : valid t)
    (large : 9 ≤ rank t) : charge t+8 ≤ 8*rank t := by
  cases t <;> simp only [valid, rank, charge] at *
  all_goals first | (split <;> omega) | omega

def forestRank (xs : List FiniteType) : Nat := (xs.map rank).sum
def forestCharge (xs : List FiniteType) : Nat := (xs.map charge).sum

theorem accumulated_deficit (xs : List FiniteType)
    (allowed : ∀ t ∈ xs, valid t) :
    forestCharge xs+4*xs.length ≤ 8*forestRank xs := by
  induction xs with
  | nil => simp [forestCharge, forestRank]
  | cons t ts ih =>
    have first := component_deficit t (allowed t (by simp))
    have rest := ih (by intro u hu; exact allowed u (by simp [hu]))
    simp only [forestCharge, forestRank, List.map_cons, List.sum_cons,
      List.length_cons, Nat.mul_add, Nat.mul_one] at *
    omega

theorem all_rank_vertex_budget (xs : List FiniteType)
    (allowed : ∀ t ∈ xs, valid t) (large : 9 ≤ forestRank xs) :
    forestCharge xs+8 ≤ 8*forestRank xs := by
  cases xs with
  | nil => simp [forestRank] at large
  | cons t ts =>
    cases ts with
    | nil =>
      have bound := large_component_deficit t (allowed t (by simp))
        (by simpa [forestRank] using large)
      simpa [forestCharge, forestRank] using bound
    | cons u us =>
      have bound := accumulated_deficit (t::u::us) allowed
      simp only [List.length_cons] at bound
      omega

end CoxeterProofs.CompactFieldBudget
