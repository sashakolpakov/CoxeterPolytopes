import Std

/-! Complete low8 absolute-moment transfers at link rank57. The finite
checks recompute the two-uniform distributions, not saved stop-loss data.
Their nonnegative signed moment remains nonnegative after convolution
with any nonnegative remaining coefficients supported in degrees0..57. -/
namespace CoxeterProofs.Concentration58

def pairCount (a b k : Nat) : Nat :=
  ((List.range (a+1)).filter fun i => i ≤ k && k-i ≤ b).length

def moment (a b shift : Nat) : Int :=
  ((List.range (a+b+1)).map fun k =>
    (pairCount a b k : Int)*Int.natAbs (2*((k+shift : Nat) : Int)-57)).sum

def momentDelta (a b shift : Nat) : Int :=
  ((a+1)*(b+1) : Nat)*moment (a-1) (b+1) shift -
    (a*(b+2) : Nat)*moment a b shift

def checkTransfer (ab : Nat × Nat) : Bool :=
  let a := ab.1; let b := ab.2
  decide (0 < a*(b+2)) && decide (a*(b+2) < (a+1)*(b+1)) &&
  (List.range 58).all (fun shift => decide (0 ≤ momentDelta a b shift))

def transfers : List (Nat × Nat) :=
  (List.range 8).flatMap fun a =>
    ((List.range 8).filter fun b => 1 ≤ a && a ≤ b).map (a,·)

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem transfers_checked : transfers.all checkTransfer = true := by decide

theorem transfer_count : transfers.length = 28 := by decide

theorem transfer_member (a b : Nat) (positive : 1 ≤ a)
    (ordered : a ≤ b) (bound : b < 8) : (a,b) ∈ transfers := by
  apply List.mem_flatMap.mpr
  refine ⟨a, List.mem_range.mpr (by omega), ?_⟩
  apply List.mem_map.mpr
  refine ⟨b, ?_, rfl⟩
  apply List.mem_filter.mpr
  exact ⟨List.mem_range.mpr bound, by simp [positive, ordered]⟩

theorem shift_moment_nonnegative (a b shift : Nat)
    (positive : 1 ≤ a) (ordered : a ≤ b) (bound : b < 8)
    (shiftBound : shift < 58) : 0 ≤ momentDelta a b shift := by
  have checked := (List.all_eq_true.mp transfers_checked) (a,b)
    (transfer_member a b positive ordered bound)
  simp only [checkTransfer, Bool.and_eq_true, decide_eq_true_eq] at checked
  exact of_decide_eq_true ((List.all_eq_true.mp checked.2) shift
    (List.mem_range.mpr shiftBound))

theorem weighted_sum_nonnegative (xs : List Nat) (weight value : Nat → Int)
    (weights : ∀ i ∈ xs, 0 ≤ weight i)
    (values : ∀ i ∈ xs, 0 ≤ value i) :
    0 ≤ (xs.map (fun i => weight i*value i)).sum := by
  induction xs with
  | nil => simp
  | cons i xs ih =>
    have head := Int.mul_nonneg (weights i (by simp)) (values i (by simp))
    have tail := ih (by intro j hj; exact weights j (by simp [hj]))
      (by intro j hj; exact values j (by simp [hj]))
    simp only [List.map_cons, List.sum_cons]
    omega

theorem convolution_moment_nonnegative (a b : Nat) (weight : Nat → Int)
    (positive : 1 ≤ a) (ordered : a ≤ b) (bound : b < 8)
    (weights : ∀ shift < 58, 0 ≤ weight shift) :
    0 ≤ ((List.range 58).map (fun shift =>
      weight shift*momentDelta a b shift)).sum := by
  apply weighted_sum_nonnegative
  · intro shift member; exact weights shift (List.mem_range.mp member)
  · intro shift member
    exact shift_moment_nonnegative a b shift positive ordered bound
      (List.mem_range.mp member)

end CoxeterProofs.Concentration58
