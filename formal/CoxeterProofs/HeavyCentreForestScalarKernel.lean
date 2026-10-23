import Std

/-!
PREPARED SOURCE ONLY: uncompiled, with zero actual axiom-audit outputs.

Ordered-ring scalar algebra supporting the written finite57 heavy-centre
forest theorem. Primitive square equations imply the long-edge cone gaps;
cleared cosine/cotangent rules then imply strict distance and cone growth.
Strong Nat induction derives no-return from two local geometric bridges.
The final Nat declarations count incidences, not real lengths or heights.

External: the actual heavy graph and its three centre edges, complete
finite-centre classification, projected-normal interpretation, hyperbolic
trigonometric rules, local angle/triangle bridges, graph cycle-to-path
interpretation, forest Euler identity and actual receiver-port bound.
No global no-return or no-cycle premise or geometric instance is supplied.
Only Std is imported; no existing proof module is imported or modified.
-/

namespace CoxeterProofs.HeavyCentreForestScalarKernel

open Lean.Grind

universe u
variable {R : Type u} [CommRing R] [LE R] [LT R]
  [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [OrderedRing R]

attribute [local instance] Semiring.natCast Ring.intCast

/-- Square equations, rather than supplied chords, give all cone bounds. -/
theorem cone_primitive_bounds (c s t : R)
    (length_lower : 16 ≤ c) (s_nonnegative : 0 ≤ s)
    (hyperbolic_square : s ^ 2 = c ^ 2 - 1)
    (t_nonnegative : 0 ≤ t) (angular_square : 2 * t ^ 2 = 1) :
    256 ≤ c ^ 2 ∧ 4 * t < 3 ∧ 15 * c < 16 * s ∧ s < c := by
  have c_square_lower : 256 ≤ c ^ 2 := by
    have first : 0 ≤ c - 16 := by grind
    have second : 0 ≤ c + 16 := by grind
    have product := OrderedRing.mul_nonneg first second
    have expansion : (c - 16) * (c + 16) = c ^ 2 - 256 := by grind only
    rw [expansion] at product
    grind
  have t_cap : 4 * t < 3 := by
    by_cases target : 4 * t < 3
    · exact target
    · have first : 0 ≤ 4 * t - 3 := by grind
      have second : 0 ≤ 4 * t + 3 := by grind
      have product := OrderedRing.mul_nonneg first second
      have expansion : (4 * t - 3) * (4 * t + 3) = 16 * t ^ 2 - 9 := by
        grind only
      rw [expansion] at product
      grind
  have s_cap : s < c := by
    by_cases target : s < c
    · exact target
    · have first : 0 ≤ s - c := by grind
      have second : 0 ≤ s + c := by grind
      have product := OrderedRing.mul_nonneg first second
      have expansion : (s - c) * (s + c) = s ^ 2 - c ^ 2 := by grind only
      rw [expansion] at product
      grind
  have s_lower : 15 * c < 16 * s := by
    by_cases target : 15 * c < 16 * s
    · exact target
    · have first : 0 ≤ 15 * c - 16 * s := by grind
      have second : 0 ≤ 15 * c + 16 * s := by grind
      have product := OrderedRing.mul_nonneg first second
      have expansion :
          (15 * c - 16 * s) * (15 * c + 16 * s) =
            225 * c ^ 2 - 256 * s ^ 2 := by grind only
      rw [expansion] at product
      grind
  exact ⟨c_square_lower, t_cap, s_lower, s_cap⟩

/-- The two strict cone gaps follow from the primitive scalar data. -/
theorem long_edge_cone_gaps (c s t : R)
    (length_lower : 16 ≤ c) (s_nonnegative : 0 ≤ s)
    (hyperbolic_square : s ^ 2 = c ^ 2 - 1)
    (t_nonnegative : 0 ≤ t) (angular_square : 2 * t ^ 2 = 1) :
    4 < c - s * t ∧ 3 < s - c * t := by
  obtain ⟨_c_square_lower, t_cap, s_lower, s_cap⟩ :=
    cone_primitive_bounds c s t length_lower s_nonnegative
      hyperbolic_square t_nonnegative angular_square
  have angular_product := OrderedRing.mul_pos
    (show 0 < c by grind) (show 0 < 3 - 4 * t by grind)
  have norm_product := OrderedRing.mul_nonneg
    (show 0 ≤ c - s by grind) t_nonnegative
  have first_gap : 0 < 4 * (c - s * t) - c := by grind
  have second_gap : 0 < 16 * (s - c * t) - 3 * c := by grind
  constructor <;> grind

/-- The cosine rule gives strict growth, including negative values of cosB. -/
theorem cosine_rule_growth (c s t cR sR cosB cD : R)
    (length_lower : 16 ≤ c) (s_nonnegative : 0 ≤ s)
    (hyperbolic_square : s ^ 2 = c ^ 2 - 1)
    (t_nonnegative : 0 ≤ t) (angular_square : 2 * t ^ 2 = 1)
    (old_distance_positive : 0 < cR)
    (old_sinh_nonnegative : 0 ≤ sR) (old_sinh_cap : sR ≤ cR)
    (angle_cap : cosB ≤ t)
    (cosine_rule : cD = cR * c - sR * s * cosB) :
    4 * cR < cD := by
  have gaps := long_edge_cone_gaps c s t length_lower s_nonnegative
    hyperbolic_square t_nonnegative angular_square
  have first := OrderedRing.mul_nonneg
    (show 0 ≤ cR - sR by grind) t_nonnegative
  have second := OrderedRing.mul_nonneg old_sinh_nonnegative
    (show 0 ≤ t - cosB by grind)
  have angle_difference : 0 ≤ cR * t - sR * cosB := by grind
  have projected_difference := OrderedRing.mul_nonneg s_nonnegative angle_difference
  have strict_growth := OrderedRing.mul_pos old_distance_positive
    (show 0 < c - s * t - 4 by grind)
  have decomposition :
      cD - 4 * cR =
        cR * (c - s * t - 4) + s * (cR * t - sR * cosB) := by
    grind only
  have target_difference : 0 < cD - 4 * cR := by
    rw [decomposition]
    grind
  grind

/-- A cleared cotangent rule preserves a strict cone without division. -/
theorem cotangent_rule_growth (c s t cothR cosB sinB cotC : R)
    (length_lower : 16 ≤ c) (s_nonnegative : 0 ≤ s)
    (hyperbolic_square : s ^ 2 = c ^ 2 - 1)
    (t_nonnegative : 0 ≤ t) (angular_square : 2 * t ^ 2 = 1)
    (coth_lower : 1 ≤ cothR) (angle_cap : cosB ≤ t)
    (sine_positive : 0 < sinB) (sine_cap : sinB ≤ 1)
    (cotangent_rule : sinB * cotC = cothR * s - c * cosB) :
    3 < cotC := by
  have gaps := long_edge_cone_gaps c s t length_lower s_nonnegative
    hyperbolic_square t_nonnegative angular_square
  have first := OrderedRing.mul_nonneg
    (show 0 ≤ cothR - 1 by grind) s_nonnegative
  have second := OrderedRing.mul_nonneg
    (show 0 ≤ c by grind) (show 0 ≤ t - cosB by grind)
  have numerator_lower : 3 < cothR * s - c * cosB := by grind
  by_cases target : 3 < cotC
  · exact target
  · have product := OrderedRing.mul_nonneg
      (show 0 ≤ 3 - cotC by grind) (show 0 ≤ sinB by grind)
    grind

/-- Both consequences use the SAME primitive long-edge and angular data. -/
theorem triangle_cone_step (c s t cR sR cosB cD cothR sinB cotC : R)
    (length_lower : 16 ≤ c) (s_nonnegative : 0 ≤ s)
    (hyperbolic_square : s ^ 2 = c ^ 2 - 1)
    (t_nonnegative : 0 ≤ t) (angular_square : 2 * t ^ 2 = 1)
    (old_distance_positive : 0 < cR)
    (old_sinh_nonnegative : 0 ≤ sR) (old_sinh_cap : sR ≤ cR)
    (angle_cap : cosB ≤ t)
    (cosine_rule : cD = cR * c - sR * s * cosB)
    (coth_lower : 1 ≤ cothR)
    (sine_positive : 0 < sinB) (sine_cap : sinB ≤ 1)
    (cotangent_rule : sinB * cotC = cothR * s - c * cosB) :
    4 * cR < cD ∧ 3 < cotC := by
  exact ⟨cosine_rule_growth c s t cR sR cosB cD length_lower
      s_nonnegative hyperbolic_square t_nonnegative angular_square
      old_distance_positive old_sinh_nonnegative old_sinh_cap angle_cap cosine_rule,
    cotangent_rule_growth c s t cothR cosB sinB cotC length_lower
      s_nonnegative hyperbolic_square t_nonnegative angular_square coth_lower
      angle_cap sine_positive sine_cap cotangent_rule⟩

/-- Local angle control excludes the next return; the nondegenerate
triangle then propagates the cone. Strong induction handles every prefix.
The Nat arguments are path indices, never metric coordinates. -/
theorem finite_chain_cone_invariant (V : Type u) (path : Nat → V)
    (n : Nat) (cone : Nat → Prop)
    (_length_positive : 1 ≤ n)
    (initial_cone : cone 1) (initial_no_return : path 1 ≠ path 0)
    (local_no_return : ∀ i, 1 ≤ i → i < n →
      cone i → path (i + 1) ≠ path 0)
    (local_cone_step : ∀ i, 1 ≤ i → i < n →
      cone i → path (i + 1) ≠ path 0 → cone (i + 1)) :
    ∀ k, 1 ≤ k → k ≤ n → cone k ∧ path k ≠ path 0 := by
  intro k
  induction k using Nat.strongRecOn with
  | ind k ih =>
    intro lower upper
    by_cases first : k = 1
    · subst k
      exact ⟨initial_cone, initial_no_return⟩
    · have predecessor_positive : 1 ≤ k - 1 := by omega
      have predecessor_smaller : k - 1 < k := by omega
      have predecessor_before_last : k - 1 < n := by omega
      obtain ⟨previous_cone, _previous_no_return⟩ :=
        ih (k - 1) predecessor_smaller predecessor_positive (by omega)
      have next_no_return : path ((k - 1) + 1) ≠ path 0 :=
        local_no_return (k - 1) predecessor_positive predecessor_before_last previous_cone
      have next_cone : cone ((k - 1) + 1) :=
        local_cone_step (k - 1) predecessor_positive predecessor_before_last
          previous_cone next_no_return
      have successor : (k - 1) + 1 = k := by omega
      simpa only [successor] using And.intro next_cone next_no_return

/-- A finite chain cannot return to its start, using only the two LOCAL bridges. -/
theorem finite_chain_no_return (V : Type u) (path : Nat → V)
    (n : Nat) (cone : Nat → Prop) (length_positive : 1 ≤ n)
    (initial_cone : cone 1) (initial_no_return : path 1 ≠ path 0)
    (local_no_return : ∀ i, 1 ≤ i → i < n →
      cone i → path (i + 1) ≠ path 0)
    (local_cone_step : ∀ i, 1 ≤ i → i < n →
      cone i → path (i + 1) ≠ path 0 → cone (i + 1)) :
    path n ≠ path 0 := by
  have invariant := finite_chain_cone_invariant V path n cone length_positive
    initial_cone initial_no_return local_no_return local_cone_step
  exact (invariant n length_positive (by omega)).2

/-- The projected-normal equations imply the length bound; none is supplied. -/
theorem projected_centre_length_lower (rank delta shadow wall coshL : R)
    (rank_lower : 14 ≤ rank) (delta_positive : 0 < delta)
    (projected_norm : (rank + 1) * delta = 1)
    (shadow_lower : 1 + delta ≤ shadow) (wall_acute : wall ≤ 0)
    (normalized_product : delta * coshL = shadow - wall) :
    rank + 2 ≤ coshL ∧ 16 ≤ coshL := by
  have lower : rank + 2 ≤ coshL := by
    by_cases target : rank + 2 ≤ coshL
    · exact target
    · have product := OrderedRing.mul_pos delta_positive
        (show 0 < rank + 2 - coshL by grind)
      have expansion :
          delta * (rank + 2 - coshL) = 1 + delta - shadow + wall := by
        grind only
      rw [expansion] at product
      grind
  exact ⟨lower, by grind⟩

/-- Pure incidence algebra; the graph and its forest Euler identity are external. -/
theorem forest_boundary_identity (H E components boundary : Nat)
    (degree_three : 3 * H = 2 * E + boundary)
    (forest_euler : E + components = H) :
    boundary = H + 2 * components := by
  omega

/-- A connected component therefore has exactly H+2 outside incidences. -/
theorem connected_component_exits (H E boundary : Nat)
    (degree_three : 3 * H = 2 * E + boundary)
    (tree_euler : E + 1 = H) :
    boundary = H + 2 := by
  have count := forest_boundary_identity H E 1 boundary degree_three tree_euler
  omega

/-- With an ACTUAL three-port bound, the count gives a receiver-count bound only. -/
theorem finite_receiver_count (H E components boundary receivers : Nat)
    (degree_three : 3 * H = 2 * E + boundary)
    (forest_euler : E + components = H)
    (actual_port_bound : boundary ≤ 3 * receivers) :
    H + 2 * components ≤ 3 * receivers := by
  have count := forest_boundary_identity H E components boundary degree_three forest_euler
  omega

end CoxeterProofs.HeavyCentreForestScalarKernel

#print axioms CoxeterProofs.HeavyCentreForestScalarKernel.cone_primitive_bounds
#print axioms CoxeterProofs.HeavyCentreForestScalarKernel.long_edge_cone_gaps
#print axioms CoxeterProofs.HeavyCentreForestScalarKernel.cosine_rule_growth
#print axioms CoxeterProofs.HeavyCentreForestScalarKernel.cotangent_rule_growth
#print axioms CoxeterProofs.HeavyCentreForestScalarKernel.triangle_cone_step
#print axioms CoxeterProofs.HeavyCentreForestScalarKernel.finite_chain_cone_invariant
#print axioms CoxeterProofs.HeavyCentreForestScalarKernel.finite_chain_no_return
#print axioms CoxeterProofs.HeavyCentreForestScalarKernel.projected_centre_length_lower
#print axioms CoxeterProofs.HeavyCentreForestScalarKernel.forest_boundary_identity
#print axioms CoxeterProofs.HeavyCentreForestScalarKernel.connected_component_exits
#print axioms CoxeterProofs.HeavyCentreForestScalarKernel.finite_receiver_count
