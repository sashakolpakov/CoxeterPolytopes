import CoxeterProofs.IdealProfileArithmetic

namespace CoxeterProofs.Ideal58
open CoxeterProofs.IdealProfile

def indexSum (h : List Nat) : Nat :=
  ((List.range h.length).map fun j => min j (57-j)*h[j]?.getD 0).sum

def check (envelopes : List FractionData) (p : Profile) : Bool :=
  decide (p.ranks.sum = 57) &&
  decide (polynomial p.ranks = p.coefficients) &&
  decide (p.coefficients.length = 58) &&
  decide (vertexCount p.ranks = p.vertices) &&
  decide (p.coefficients.sum = p.vertices) &&
  decide (indexSum p.coefficients = p.indexTotal) &&
  decide (0 < p.cost.den) &&
  decide (0 < (envelopeCost envelopes p.ranks).den) &&
  decide ((envelopeCost envelopes p.ranks).num*p.cost.den =
    p.cost.num*(envelopeCost envelopes p.ranks).den) &&
  decide (72*p.cost.den+57*p.vertices*p.cost.den+
    2*p.vertices*p.cost.num ≤ 5*p.indexTotal*p.cost.den) &&
  decide (291927*p.vertices*p.cost.den+10240*p.vertices*p.cost.num ≤
    25600*p.indexTotal*p.cost.den)

theorem absolute_gap_of_check (envelopes : List FractionData) (p : Profile)
    (checked : check envelopes p = true) :
    72*p.cost.den+57*p.vertices*p.cost.den+2*p.vertices*p.cost.num ≤
      5*p.indexTotal*p.cost.den := by
  simp only [check, Bool.and_eq_true, decide_eq_true_eq] at checked
  exact checked.1.2

def terminalDomain : List (List Nat) :=
  (List.range 8).flatMap fun low =>
    (highParts 58 (57-low) 8).map fun rs => if low = 0 then rs else low::rs

end CoxeterProofs.Ideal58
