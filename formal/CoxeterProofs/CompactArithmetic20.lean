import CoxeterProofs.CompactHVector

/-!
The final rank21 counting logic for the internally proved compact20 bound.
The universal geometric reduction and rational angle certificate are explicit
premises here. This module does not formalize those premises or identify
them with a hyperbolic polytope, and declares no new axioms.
-/
namespace CoxeterProofs.CompactArithmetic20

open CompactHVector

def pairedCoefficients21 : List Nat :=
  [630, 530, 440, 360, 290, 230, 180, 140, 110, 90, 80]

/-- Coefficient of h_i after pairing i with 21-i in 5*f₂-210*f₀.
The formula is used only at indices 0 through 10. -/
def pairedCoefficient21 (i : Nat) : Nat :=
  5 * (i * (i - 1) / 2 + (21 - i) * (20 - i) / 2) - 420

theorem paired_coefficients21_formula :
    (List.range 11).map pairedCoefficient21 = pairedCoefficients21 := by decide

theorem checked_paired_coefficients21 :
    pairedCoefficients21.all (fun c => decide (80 ≤ c)) = true := by decide

def endpoint21 : Term :=
  { multiplicity := 2, coefficient := 630, mass := 1 }

theorem endpoint_surplus21 (tail : List Term)
    (coefficients : ∀ t ∈ tail, 40 * t.multiplicity ≤ t.coefficient) :
    40 * vertices (endpoint21 :: tail) + 550 ≤
      weighted (endpoint21 :: tail) := by
  have rest := weighted_lower 40 tail coefficients
  simp only [vertices, weighted, List.map_cons, List.sum_cons, endpoint21,
    Nat.mul_one] at *
  omega

theorem no_rank21_weight_model (tail : List Term) (denominator total : Nat)
    (positiveDenominator : 0 < denominator)
    (coefficients : ∀ t ∈ tail, 40 * t.multiplicity ≤ t.coefficient)
    (incidenceLower : denominator * weighted (endpoint21 :: tail) ≤ total)
    (certificateUpper : total ≤
      denominator * (40 * vertices (endpoint21 :: tail))) : False := by
  have lower := Nat.mul_le_mul_left denominator (endpoint_surplus21 tail coefficients)
  rw [Nat.mul_add] at lower
  have strict : 0 < denominator * 550 := Nat.mul_pos positiveDenominator (by decide)
  omega

theorem ceiling20_from_remaining_golden_branch (n : Nat) (goldenBranch : Prop)
    (previousCeiling : n ≤ 21)
    (fieldReduction : n = 21 → goldenBranch)
    (goldenExcluded : ¬ goldenBranch) : n ≤ 20 := by
  have excluded21 : n ≠ 21 := by
    intro equality
    exact goldenExcluded (fieldReduction equality)
  omega

end CoxeterProofs.CompactArithmetic20

#print axioms CoxeterProofs.CompactArithmetic20.paired_coefficients21_formula
#print axioms CoxeterProofs.CompactArithmetic20.endpoint_surplus21
#print axioms CoxeterProofs.CompactArithmetic20.no_rank21_weight_model
