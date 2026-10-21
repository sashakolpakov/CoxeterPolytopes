import Std

namespace CoxeterProofs.FiveFaceCarrier

theorem ideal5_aggregate_identity (v e r g f D A K : Int)
    (edges : e = 8*v) (threeFaces : 3*r = 24*v+2*g)
    (euler : v-e+r-g+f = 2)
    (deltaSum : D = 8*v-2*g) (facetExcess : A = 2*g-24*f)
    (kappaSum : K = 3*D-A) : K = 48 := by omega

theorem positive_member_of_positive_sum (xs : List Int)
    (positive : 0 < xs.sum) : ∃ a ∈ xs, 0 < a := by
  induction xs with
  | nil => simp at positive
  | cons a xs ih =>
    by_cases ha : 0 < a
    · exact ⟨a, by simp, ha⟩
    · have rest : 0 < xs.sum := by simp only [List.sum_cons] at positive; omega
      obtain ⟨b, hb, hp⟩ := ih rest
      exact ⟨b, List.mem_cons_of_mem a hb, hp⟩

theorem positive_carrier (kappas : List Int) (identity : kappas.sum = 48) :
    ∃ k ∈ kappas, 0 < k :=
  positive_member_of_positive_sum kappas (by omega)

theorem carrier_is_beyond_two_known_values (delta : Nat) (kappa : Int)
    (gap : delta = 0 ∨ delta = 4 ∨ 5 ≤ delta)
    (minimum : delta = 0 → kappa = 0)
    (second : delta = 4 → kappa = -2)
    (positive : 0 < kappa) : 5 ≤ delta := by
  rcases gap with low | next | high
  · have := minimum low; omega
  · have := second next; omega
  · exact high

theorem mixed_weighted_incidence (fiveFaces nonideal ideal incidence finiteVertices
    cofaces positiveKappa : Nat)
    (partition : fiveFaces = nonideal+ideal)
    (finiteCover : nonideal ≤ incidence*finiteVertices)
    (idealCover : 48*ideal ≤ cofaces*positiveKappa) :
    48*fiveFaces ≤ 48*(incidence*finiteVertices)+cofaces*positiveKappa := by
  have scaled := Nat.mul_le_mul_left 48 finiteCover
  rw [partition, Nat.mul_add]
  omega

theorem known_model_kappas : (3*(42:Int)-4*38+24 = -2) ∧
    (3*(24:Int)-4*24+24 = 0) := by decide

end CoxeterProofs.FiveFaceCarrier
