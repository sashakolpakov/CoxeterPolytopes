import Std

/-!
Prepared exact arithmetic source. Compilation and actual axiom-output
verification are recorded separately; preparation alone is not verification.

The integers below are combinatorial scalar numerators for M = 0,...,14.
They are not real heights, coordinates, curvature or geometry variables.
The written proof supplies coefficient extraction, canonical envelopes,
convex order, actual incidence and charge allocation. No geometry or
dimension conclusion is a theorem of this module.
-/

namespace CoxeterProofs.LargeInertParentArithmetic

/-- The zero convention includes negative upper indices. -/
def binom4 (n : Int) : Int :=
  if n < 4 then 0 else n * (n - 1) * (n - 2) * (n - 3) / 24

/-- The coefficient [t^27] H(t)/(1-t)^2 for ranks (16,40-M,M). -/
def coefficient (M : Nat) : Int :=
  let m : Int := M
  30464 - binom4 (30 - m) + binom4 (13 - m) - binom4 (m - 10)

def linkVertexCount (M : Nat) : Int :=
  let m : Int := M
  17 * (41 - m) * (m + 1)

/-- Twice the canonical single-small-factor surplus before source debt. -/
def canonicalSurplus2 (M : Nat) : Int :=
  let m : Int := M
  (84 - m) * linkVertexCount M - 227 * (41 - m) * (m + 1) -
    34 * (256 - 5 * m) * (m + 1) - 10 * coefficient M

/-- Twice the full all-port debt: p=40-M is at least36 iff M<=4. -/
def sourceDebt2 (M : Nat) : Int :=
  let m : Int := M
  if M ≤ 4 then 64 * (41 - m) * (m + 1)
  else 4 * (41 - m) * (m + 1)

def remaining2 (M : Nat) : Int :=
  canonicalSurplus2 M - sourceDebt2 M

theorem coefficient_table :
    (List.range 15).map coefficient =
      [3774, 7208, 10319, 13124, 15640, 17884, 19873, 21624,
       23154, 24480, 25619, 26588, 27404, 28084, 28643] := by
  decide

theorem canonical_surplus2_table :
    (List.range 15).map canonicalSurplus2 =
      [2797, 5572, 8257, 10784, 13085, 15092, 16737, 17952,
       18669, 18820, 18337, 17152, 15197, 12404, 8725] := by
  decide

theorem remaining2_table :
    (List.range 15).map remaining2 =
      [173, 452, 769, 1056, 1245, 14228, 15757, 16864,
       17481, 17540, 16973, 15712, 13689, 10836, 7105] := by
  decide

/-- Every allowed M retains twice the required absolute surplus 173/2. -/
theorem remaining2_lower_bound :
    ∀ M : Fin 15, (173 : Int) ≤ remaining2 M.val := by
  decide

end CoxeterProofs.LargeInertParentArithmetic

#print axioms CoxeterProofs.LargeInertParentArithmetic.coefficient_table
#print axioms CoxeterProofs.LargeInertParentArithmetic.canonical_surplus2_table
#print axioms CoxeterProofs.LargeInertParentArithmetic.remaining2_table
#print axioms CoxeterProofs.LargeInertParentArithmetic.remaining2_lower_bound
