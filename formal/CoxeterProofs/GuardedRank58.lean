import Std

/-! Exact rank-capacity consequences used by dimension-specific contextual
rules. These do not assert that the local geometric replacements are valid. -/
namespace CoxeterProofs.GuardedRank58

theorem d13_three_large (d13 large rest : Nat)
    (rank : 13*d13+15*large+rest = 58)
    (component : 1 ≤ d13) (witnesses : 3 ≤ large) :
    d13 = 1 ∧ large = 3 ∧ rest = 0 := by omega

theorem d28_two_large (d28 large rest : Nat)
    (rank : 28*d28+15*large+rest = 58)
    (component : 1 ≤ d28) (witnesses : 2 ≤ large) :
    d28 = 1 ∧ large = 2 ∧ rest = 0 := by omega

theorem a12_d16_two_large (a12 d16 large rest : Nat)
    (rank : 12*a12+16*d16+15*large+rest = 58)
    (small : 1 ≤ a12) (big : 1 ≤ d16) (witnesses : 2 ≤ large) :
    a12 = 1 ∧ d16 = 1 ∧ large = 2 ∧ rest = 0 := by omega

theorem small_whole_factor_precludes_guard (small : Nat)
    (positive : 1 ≤ small) (rankCapacity : 13+3*15+small ≤ 58) : False := by omega

theorem repaired_forest_costs :
    (77+3*100 : Nat) = 377 ∧
    (2*170+2*200+1 : Nat) ≤ 2*377 ∧
    (2*66+217+2*200+2 : Nat) ≤ 2*377 := by decide

theorem global_contradiction (idealVertices finiteSurplus : Nat)
    (comparison : 4205+36*idealVertices+finiteSurplus ≤ 0) : False := by omega

theorem ceiling57_from_tail_and_exclusion (n : Nat)
    (tail : n ≤ 58) (exclude58 : n ≠ 58) : n ≤ 57 := by omega

end CoxeterProofs.GuardedRank58
