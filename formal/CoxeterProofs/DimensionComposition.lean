import Std

/-!
Only the final arithmetic/field composition is formalized here. The internal
geometric and certificate premises are explicit arguments, not axioms and
not claims that the full dimension theorems have been formalized.
-/
namespace CoxeterProofs.Dimension

theorem no_ambient59_from_scaled_global_inequality
    (idealVertices finiteSurplus : Nat)
    (comparison : 8700 + 91 * idealVertices + finiteSurplus ≤ 0) : False := by
  omega

theorem ceiling58_from_tail_and_exclusion
    (n : Nat) (tail : n ≤ 59) (excluded59 : n ≠ 59) : n ≤ 58 := by
  omega

inductive HighArithmeticField where
  | sqrtTwo
  | sqrtFive
  | cubicSeven
  deriving DecidableEq, Repr

def fieldCeiling : HighArithmeticField → Nat
  | .sqrtTwo => 20
  | .sqrtFive => 21
  | .cubicSeven => 20

theorem compact_arithmetic_field_composition
    (n : Nat)
    (compactCeiling : n ≤ 29)
    (highFieldBounds : 22 ≤ n → n ≤ 29 →
      ∃ k : HighArithmeticField, n ≤ fieldCeiling k) : n ≤ 21 := by
  by_cases low : n ≤ 21
  · exact low
  · obtain ⟨k, bound⟩ := highFieldBounds (by omega) compactCeiling
    cases k <;> simp [fieldCeiling] at bound <;> omega

end CoxeterProofs.Dimension
