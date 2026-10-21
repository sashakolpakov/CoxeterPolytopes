import CoxeterProofs.CertificateBounds

namespace CoxeterProofs.ForestDomination

theorem getD_property {α : Type} (xs : List α) (fallback : α)
    (P : α → Prop) (all : ∀ x ∈ xs, P x) (base : P fallback) :
    ∀ i : Nat, P (xs[i]?.getD fallback) := by
  induction xs with
  | nil => intro i; simpa using base
  | cons a xs ih =>
    intro i
    cases i with
    | zero => simpa using all a (by simp)
    | succ i =>
      simpa using ih (by intro x hx; exact all x (by simp [hx])) i

def cost (base : Nat → Nat) (cross : Nat → Nat → Nat) : List Nat → Nat
  | [] => 0
  | a::xs => base a + (xs.map (cross a)).sum + cost base cross xs

theorem sum_map_le (xs : List Nat) (f g : Nat → Nat)
    (pointwise : ∀ a ∈ xs, f a ≤ g a) :
    (xs.map f).sum ≤ (xs.map g).sum := by
  induction xs with
  | nil => simp
  | cons a xs ih =>
    have head := pointwise a (by simp)
    have tail := ih (by intro b hb; exact pointwise b (by simp [hb]))
    simp only [List.map_cons, List.sum_cons]
    omega

theorem ordinary_cost_domination
    (base upperBase : Nat → Nat) (cross upperCross : Nat → Nat → Nat)
    (xs : List Nat)
    (bases : ∀ a ∈ xs, base a ≤ upperBase a)
    (pairs : ∀ a ∈ xs, ∀ b ∈ xs, cross a b ≤ upperCross a b) :
    cost base cross xs ≤ cost upperBase upperCross xs := by
  induction xs with
  | nil => simp [cost]
  | cons a xs ih =>
    have base := bases a (by simp)
    have pair := sum_map_le xs (cross a) (upperCross a) (by
      intro b hb
      exact pairs a (by simp) b (by simp [hb]))
    have rest := ih
      (by intro b hb; exact bases b (by simp [hb]))
      (by intro b hb c hc; exact pairs b (by simp [hb]) c (by simp [hc]))
    simp only [cost]
    omega

end CoxeterProofs.ForestDomination
