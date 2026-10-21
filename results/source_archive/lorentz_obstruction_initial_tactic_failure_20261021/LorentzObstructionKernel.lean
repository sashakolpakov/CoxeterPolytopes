import Std

/-!
Prepared internal algebra for two Lorentz-space obstructions.
Compilation and a separate actual axiom audit remain pending.

Derived: a general norm-three cutoff from the finite-orbit Cauchy budget,
its norm-two instance, the explicit golden pair comparison, and the
large-D four-by-four Schur numerator sign obstruction.

External application inputs: actual isometry and invariant timelike point,
the Euclidean decomposition/Cauchy realization, actual Gram elimination,
complete white support/native gap, and original-lattice/chamber interfaces.
No actual polytope, reflection lattice or dimension bound is formalized.
-/

namespace CoxeterProofs.LorentzObstructionKernel

open Lean.Grind

universe u
variable {R : Type u} [CommRing R] [LE R] [LT R]
  [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [OrderedRing R]

attribute [local instance] Semiring.natCast Ring.intCast

/-- Finite-order Cauchy data force second-orbit product greater than -3c. -/
theorem second_orbit_gap (c q W second : R)
    (norm_positive : 0 < c) (time_square_nonnegative : 0 ≤ q)
    (spatial_norm : W = q + c)
    (cauchy : 2 * q ^ 2 ≤ W * (W + second + q)) :
    -3 * c < second := by
  have spatial_positive : 0 < W := by grind
  have norm_square_positive := OrderedRing.mul_pos norm_positive norm_positive
  by_cases bad : second ≤ -3 * c
  · have bracket : W + second + q ≤ 2 * q - 2 * c := by grind
    have upper := OrderedRing.mul_le_mul_of_nonneg_left bracket
      (show 0 ≤ W by grind)
    have expansion : (q + c) * (2 * q - 2 * c) =
        2 * q ^ 2 - 2 * c ^ 2 := by grind only
    rw [spatial_norm, expansion] at upper
    grind
  · grind

theorem norm_two_second_orbit_gap (q W second : R)
    (time_square_nonnegative : 0 ≤ q)
    (spatial_norm : W = q + 2)
    (cauchy : 2 * q ^ 2 ≤ W * (W + second + q)) :
    -6 < second := by
  have result := second_orbit_gap (2 : R) q W second
    (by grind) time_square_nonnegative spatial_norm cauchy
  grind

theorem golden_pair_beyond_cutoff (phi second : R)
    (golden_lower : 1 < phi) (pairing : second = -80 - 132 * phi) :
    second < -6 := by
  grind

/-- Averaging, coordinate realization and Cauchy remain explicit inputs. -/
theorem no_finite_orbit_data (finiteOrder : Prop) (phi second : R)
    (golden_lower : 1 < phi) (pairing : second = -80 - 132 * phi)
    (finite_order_bridge : finiteOrder → ∃ q W : R,
      0 ≤ q ∧ W = q + 2 ∧ 2 * q ^ 2 ≤ W * (W + second + q)) :
    ¬ finiteOrder := by
  intro finite
  obtain ⟨q, W, q_nonnegative, spatial_norm, cauchy⟩ := finite_order_bridge finite
  have lower := norm_two_second_orbit_gap q W second q_nonnegative spatial_norm cauchy
  have upper := golden_pair_beyond_cutoff phi second golden_lower pairing
  grind

/-- Cleared determinant identity for the marked-D opposite-sum block. -/
theorem schur_numerator_identity (m g h : R) :
    (2 - (m + 1) * g) * (-2 * m - (m + 1) * h) - 4 =
      (m + 1) * ((m + 1) * g * h + 2 * m * g - 2 * h - 4) := by
  grind only

theorem large_d_schur_numerator_nonzero (m g h : R)
    (rank_floor : 3 ≤ m) (h_nonnegative : 0 ≤ h)
    (native_gap : g = 0 ∨ 1 ≤ g) :
    (m + 1) * g * h + 2 * m * g - 2 * h - 4 ≠ 0 := by
  rcases native_gap with zero | positive
  · rw [zero]
    grind
  · have m_nonnegative : 0 ≤ m := by grind
    have mg_lower := OrderedRing.mul_le_mul_of_nonneg_left positive m_nonnegative
    have coefficient : 0 ≤ (m + 1) * g - 2 := by
      have factor := OrderedRing.mul_le_mul_of_nonneg_left positive
        (show 0 ≤ m + 1 by grind)
      grind
    have term := OrderedRing.mul_nonneg h_nonnegative coefficient
    have rearrange : (m + 1) * g * h + 2 * m * g - 2 * h - 4 =
        h * ((m + 1) * g - 2) + 2 * m * g - 4 := by grind only
    rw [rearrange]
    grind

theorem full_white_rank_contradiction (n rank : Nat)
    (actual_rank : rank ≤ n + 1) (schur_rank : rank = n + 2) : False := by
  omega

end CoxeterProofs.LorentzObstructionKernel

#print axioms CoxeterProofs.LorentzObstructionKernel.second_orbit_gap
#print axioms CoxeterProofs.LorentzObstructionKernel.norm_two_second_orbit_gap
#print axioms CoxeterProofs.LorentzObstructionKernel.golden_pair_beyond_cutoff
#print axioms CoxeterProofs.LorentzObstructionKernel.no_finite_orbit_data
#print axioms CoxeterProofs.LorentzObstructionKernel.schur_numerator_identity
#print axioms CoxeterProofs.LorentzObstructionKernel.large_d_schur_numerator_nonzero
#print axioms CoxeterProofs.LorentzObstructionKernel.full_white_rank_contradiction
