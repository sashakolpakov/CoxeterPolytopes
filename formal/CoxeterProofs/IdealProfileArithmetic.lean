import Std

/-! Exact polynomial and rational arithmetic for the terminal ideal-link
certificate. The geometric concentration and envelope lemmas are separate
internal obligations; this module does not assume them as published axioms. -/
namespace CoxeterProofs.IdealProfile

structure FractionData where
  num : Nat
  den : Nat
  deriving DecidableEq, Repr

def add (a b : FractionData) : FractionData :=
  ⟨a.num * b.den + b.num * a.den, a.den * b.den⟩

def simplexProduct (h : List Nat) (r : Nat) : List Nat :=
  (List.range (h.length + r)).map fun k =>
    ((List.range (r + 1)).map fun j =>
      if j ≤ k then h[k-j]?.getD 0 else 0).sum

def polynomial (ranks : List Nat) : List Nat :=
  ranks.foldl simplexProduct [1]

def vertexCount (ranks : List Nat) : Nat :=
  (ranks.map (· + 1)).prod

def indexSum (h : List Nat) : Nat :=
  ((List.range h.length).map fun j => min j (58-j) * h[j]?.getD 0).sum

def envelopeCost (envelopes : List FractionData) (ranks : List Nat) : FractionData :=
  (ranks.map fun r => envelopes[r]?.getD ⟨0, 0⟩).foldl add ⟨0, 1⟩

structure Profile where
  ranks : List Nat
  coefficients : List Nat
  vertices : Nat
  indexTotal : Nat
  cost : FractionData
  deriving DecidableEq, Repr

def check (envelopes : List FractionData) (p : Profile) : Bool :=
  decide (p.ranks.sum = 58) &&
  decide (polynomial p.ranks = p.coefficients) &&
  decide (p.coefficients.length = 59) &&
  decide (vertexCount p.ranks = p.vertices) &&
  decide (p.coefficients.sum = p.vertices) &&
  decide (indexSum p.coefficients = p.indexTotal) &&
  decide (0 < p.cost.den) &&
  decide ((envelopeCost envelopes p.ranks).den > 0) &&
  decide ((envelopeCost envelopes p.ranks).num * p.cost.den =
    p.cost.num * (envelopeCost envelopes p.ranks).den) &&
  decide (91 * p.cost.den + 58 * p.vertices * p.cost.den +
    2 * p.vertices * p.cost.num ≤ 5 * p.indexTotal * p.cost.den) &&
  decide (2 * p.vertices * p.cost.den + 19720 * p.vertices * p.cost.den +
    680 * p.vertices * p.cost.num ≤ 1700 * p.indexTotal * p.cost.den)

theorem absolute_gap_of_check (envelopes : List FractionData) (p : Profile)
    (checked : check envelopes p = true) :
    91 * p.cost.den + 58 * p.vertices * p.cost.den +
      2 * p.vertices * p.cost.num ≤ 5 * p.indexTotal * p.cost.den := by
  simp only [check, Bool.and_eq_true, decide_eq_true_eq] at checked
  exact checked.1.2

theorem polynomial_of_check (envelopes : List FractionData) (p : Profile)
    (checked : check envelopes p = true) : polynomial p.ranks = p.coefficients := by
  simp only [check, Bool.and_eq_true, decide_eq_true_eq] at checked
  exact checked.1.1.1.1.1.1.1.1.1.2

def highParts : Nat → Nat → Nat → List (List Nat)
  | 0, left, _ => if left = 0 then [[]] else []
  | fuel+1, left, least =>
    if left = 0 then [[]] else
      ((List.range (left+1)).filter fun r => least ≤ r).flatMap fun r =>
        (highParts fuel (left-r) r).map (r :: ·)

def terminalDomain : List (List Nat) :=
  (List.range 9).flatMap fun low =>
    (highParts 59 (58-low) 9).map fun rs => if low = 0 then rs else low :: rs

end CoxeterProofs.IdealProfile
