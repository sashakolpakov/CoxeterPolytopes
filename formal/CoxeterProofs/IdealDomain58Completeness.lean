import CoxeterProofs.IdealDomainCompleteness
import CoxeterProofs.Generated.Ideal58

namespace CoxeterProofs.Ideal58
open CoxeterProofs.IdealProfile

theorem terminalDomain_complete (low : Nat) (rs : List Nat)
    (lowBound : low < 8) (ordered : OrderedFrom 8 rs)
    (sum : low+rs.sum = 57) :
    (if low = 0 then rs else low::rs) ∈ terminalDomain := by
  have capacity : rs.length ≤ 58 := by
    have bound := ordered_length_le_sum 8 rs (by omega) ordered
    omega
  have member := highParts_complete 58 (57-low) 8 rs (by omega) ordered
    (by omega) capacity
  apply List.mem_flatMap.mpr
  exact ⟨low, List.mem_range.mpr lowBound, List.mem_map.mpr ⟨rs, member, rfl⟩⟩

theorem canonical_terminal_profile_covered (low : Nat) (rs : List Nat)
    (lowBound : low < 8) (ordered : OrderedFrom 8 rs)
    (sum : low+rs.sum = 57) :
    ∃ p ∈ CoxeterProofs.Generated.Ideal58.profiles,
      p.ranks = (if low = 0 then rs else low::rs) ∧
      check CoxeterProofs.Generated.Ideal58.envelopes p = true := by
  have member := terminalDomain_complete low rs lowBound ordered sum
  rw [← CoxeterProofs.Generated.Ideal58.domain_complete] at member
  obtain ⟨p, hp, ranks⟩ := List.mem_map.mp member
  exact ⟨p, hp, ranks, CoxeterProofs.Generated.Ideal58.every_profile_checked p hp⟩

end CoxeterProofs.Ideal58
