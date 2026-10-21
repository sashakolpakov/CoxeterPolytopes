import Std

/-! Count consequences for the explicit C/D/T facet spectrum. Cusp
incidence equations and the orthogonality graph bound remain explicit
geometric premises. No abstract count assignment is asserted realizable. -/
namespace CoxeterProofs.RestrictedCarrierCounts

theorem six_divides_carriers (b c d doubles carriers : Nat)
    (shapeA : 6*b+4*c = 6*doubles+49*carriers)
    (shapeB : 4*c+6*d = 14*carriers) : carriers % 6 = 0 := by omega

theorem thirty_carriers (carriers : Nat) (positive : 0 < carriers)
    (divisibility : carriers % 6 = 0)
    (graphBound : 55*carriers ≤ 2*(carriers*carriers)) : 30 ≤ carriers := by
  have lower : 28 ≤ carriers := by
    by_cases h : 28 ≤ carriers
    · exact h
    · have cap : carriers ≤ 27 := by omega
      have sq := Nat.mul_le_mul_left carriers cap
      have contradiction : carriers = 0 := by omega
      omega
  omega

theorem sixty_six_doubles (doubles carriers : Nat)
    (carrierMinimum : 30 ≤ carriers)
    (kappaIdentity : 6*carriers = 48+2*doubles) : 66 ≤ doubles := by omega

theorem abstract_count_witness :
    6*(30:Nat) = 48+2*66 ∧
    8*(2051:Nat)+2*241+2*0 = 24*361+36*66+195*30 ∧
    6*(241:Nat)+4*105 = 6*66+49*30 ∧
    4*(105:Nat)+6*0 = 14*30 ∧
    (361+66+30:Nat) = 457 ∧ (2051+241+105+0:Nat) = 2397 := by decide

end CoxeterProofs.RestrictedCarrierCounts
