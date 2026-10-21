import CoxeterProofs.Generated.Ideal59

namespace CoxeterProofs.IdealProfile

def OrderedFrom : Nat → List Nat → Prop
  | _, [] => True
  | least, r::rs => least ≤ r ∧ OrderedFrom r rs

theorem ordered_length_le_sum (least : Nat) (rs : List Nat)
    (positive : 1 ≤ least) (ordered : OrderedFrom least rs) :
    rs.length ≤ rs.sum := by
  induction rs generalizing least with
  | nil => simp
  | cons r rs ih =>
    have hr : 1 ≤ r := Nat.le_trans positive ordered.1
    have rest := ih r hr ordered.2
    simp only [List.length_cons, List.sum_cons]
    omega

theorem highParts_complete (fuel left least : Nat) (rs : List Nat)
    (positive : 1 ≤ least) (ordered : OrderedFrom least rs)
    (sum : rs.sum = left) (capacity : rs.length ≤ fuel) :
    rs ∈ highParts fuel left least := by
  induction fuel generalizing left least rs with
  | zero =>
    have empty : rs = [] := List.eq_nil_of_length_eq_zero (by omega)
    subst rs
    have hz : left = 0 := by simpa using sum.symm
    simp [highParts, hz]
  | succ fuel ih =>
    cases rs with
    | nil =>
      have hz : left = 0 := by simpa using sum.symm
      simp [highParts, hz]
    | cons r rs =>
      have hr : 1 ≤ r := Nat.le_trans positive ordered.1
      have hleft : left ≠ 0 := by simp only [List.sum_cons] at sum; omega
      simp only [highParts, hleft, ↓reduceIte, List.mem_flatMap]
      refine ⟨r, ?_, ?_⟩
      · apply List.mem_filter.mpr
        exact ⟨List.mem_range.mpr (by simp only [List.sum_cons] at sum; omega),
          by simpa using ordered.1⟩
      · apply List.mem_map.mpr
        refine ⟨rs, ?_, rfl⟩
        apply ih (left-r) r rs hr ordered.2
        · simp only [List.sum_cons] at sum
          omega
        · simp only [List.length_cons] at capacity
          omega

theorem terminalDomain_complete (low : Nat) (rs : List Nat)
    (lowBound : low < 9) (ordered : OrderedFrom 9 rs)
    (sum : low + rs.sum = 58) :
    (if low = 0 then rs else low::rs) ∈ terminalDomain := by
  have capacity : rs.length ≤ 59 := by
    have bound := ordered_length_le_sum 9 rs (by omega) ordered
    omega
  have member := highParts_complete 59 (58-low) 9 rs (by omega) ordered
    (by omega) capacity
  apply List.mem_flatMap.mpr
  exact ⟨low, List.mem_range.mpr lowBound, List.mem_map.mpr ⟨rs, member, rfl⟩⟩

theorem canonical_terminal_profile_covered (low : Nat) (rs : List Nat)
    (lowBound : low < 9) (ordered : OrderedFrom 9 rs)
    (sum : low + rs.sum = 58) :
    ∃ p ∈ CoxeterProofs.Generated.Ideal59.profiles,
      p.ranks = (if low = 0 then rs else low::rs) ∧
      check CoxeterProofs.Generated.Ideal59.envelopes p = true := by
  have member := terminalDomain_complete low rs lowBound ordered sum
  rw [← CoxeterProofs.Generated.Ideal59.domain_complete] at member
  obtain ⟨p, hp, ranks⟩ := List.mem_map.mp member
  exact ⟨p, hp, ranks,
    CoxeterProofs.Generated.Ideal59.every_profile_checked p hp⟩

end CoxeterProofs.IdealProfile
