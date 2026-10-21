import Std

namespace CoxeterProofs.DoublingKappa

def kappa (vertices facets : Int) : Int := 3*vertices-4*facets+24

theorem doubling_recurrence (v f b : Int) :
    kappa (2*v-b) (2*f-b-4) = 2*kappa v f+b-8 := by
  simp only [kappa]
  omega

theorem exact_path :
    kappa 42 38 = -2 ∧ kappa 75 63 = -3 ∧
    kappa 136 108 = 0 ∧ kappa 258 198 = 6 := by decide

theorem lower_recurrence (k next b : Int) (size : 6 ≤ b)
    (step : next = 2*k+b-8) : 2*k-2 ≤ next := by omega

theorem exponential_lower (values : Nat → Int)
    (initial : 6 ≤ values 0)
    (steps : ∀ t, 2*values t-2 ≤ values (t+1)) :
    ∀ t, 4*(2:Int)^t+2 ≤ values t := by
  intro t
  induction t with
  | zero => simpa using initial
  | succ t ih =>
    have next := steps t
    simp only [Int.pow_succ]
    omega

end CoxeterProofs.DoublingKappa
