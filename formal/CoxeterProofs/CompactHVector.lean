import Std

/-! The strict compact h-vector weighting argument. Geometric f/h and
incidence identities are supplied explicitly; the strict inequality is
proved here, rather than taking the new field bounds as assumptions. -/
namespace CoxeterProofs.CompactHVector

structure Term where
  multiplicity : Nat
  coefficient : Nat
  mass : Nat
  deriving Repr

def vertices (terms : List Term) : Nat :=
  (terms.map fun t => t.multiplicity * t.mass).sum

def weighted (terms : List Term) : Nat :=
  (terms.map fun t => t.coefficient * t.mass).sum

theorem weighted_lower (threshold : Nat) (terms : List Term)
    (coefficients : ∀ t ∈ terms, threshold * t.multiplicity ≤ t.coefficient) :
    threshold * vertices terms ≤ weighted terms := by
  induction terms with
  | nil => simp [vertices, weighted]
  | cons t ts ih =>
    have tail := ih (by intro x hx; exact coefficients x (by simp [hx]))
    have head := Nat.mul_le_mul_right t.mass (coefficients t (by simp))
    simp only [vertices, weighted, List.map_cons, List.sum_cons,
      Nat.mul_add, Nat.mul_assoc] at *
    omega

theorem weighted_strict (threshold : Nat) (head : Term) (tail : List Term)
    (positive : 0 < head.mass)
    (strict : threshold * head.multiplicity < head.coefficient)
    (coefficients : ∀ t ∈ tail, threshold * t.multiplicity ≤ t.coefficient) :
    threshold * vertices (head :: tail) < weighted (head :: tail) := by
  have rest := weighted_lower threshold tail coefficients
  have first := Nat.mul_lt_mul_of_pos_right strict positive
  simp only [vertices, weighted, List.map_cons, List.sum_cons,
    Nat.mul_add, Nat.mul_assoc] at *
  omega

theorem strict_deficit_from_paired_identity
    (n betaFour deficit : Nat) (head : Term) (tail : List Term)
    (positive : 0 < head.mass)
    (strict : (2*n*(n-1)+betaFour) * head.multiplicity < head.coefficient)
    (coefficients : ∀ t ∈ tail,
      (2*n*(n-1)+betaFour) * t.multiplicity ≤ t.coefficient)
    (pairedIdentity : weighted (head::tail) =
      4*deficit + 2*n*(n-1)*vertices (head::tail)) :
    betaFour * vertices (head::tail) < 4*deficit := by
  have bound := weighted_strict _ _ _ positive strict coefficients
  rw [pairedIdentity, Nat.add_mul] at bound
  omega

def betaFour (n : Nat) : Nat :=
  if n % 2 = 0 then n*(n-6)/2 else (n-1)*(n-5)/2

def multiplicity (n i : Nat) : Nat := if 2*i = n then 1 else 2

def coefficient (n i : Nat) : Nat :=
  5 * multiplicity n i * (i*(i-1)+(n-i)*(n-i-1))

def checkDimension (n : Nat) : Bool :=
  decide ((2*n*(n-1)+betaFour n)*multiplicity n 0 < coefficient n 0) &&
  (List.range (n/2+1)).all fun i =>
    decide ((2*n*(n-1)+betaFour n)*multiplicity n i ≤ coefficient n i)

theorem checked_compact_dimensions :
    ((List.range 9).map (·+21)).all checkDimension = true := by decide

theorem coefficients_in_compact_range (n : Nat) (low : 21 ≤ n) (high : n ≤ 29) :
    checkDimension n = true := by
  apply (List.all_eq_true.mp checked_compact_dimensions) n
  exact List.mem_map.mpr ⟨n-21, List.mem_range.mpr (by omega), by omega⟩

end CoxeterProofs.CompactHVector
