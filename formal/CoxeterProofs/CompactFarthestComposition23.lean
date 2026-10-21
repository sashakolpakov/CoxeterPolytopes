import Std
import Init.Grind.Ordered

/-!
Conditional final accounting for the compact farthest-vertex paper.

Charges are rational LOWER BOUNDS multiplied by 700, not actual
energies and not a discreteness assertion. The strict budget is the
explicit bridge from the sum of real orthogonal projection energies.
Geometry, exhaustive component estimates, debit exclusivity, local caps,
the D23/short-A1 obstruction and the published ceiling remain inputs.
No dimension-23 conclusion is an input. Source-bound compilation and
axiom-audit status are recorded in formal/.proof-batches.
-/

namespace CoxeterProofs.CompactFarthestComposition23

inductive ClassicalKind where
  | A | B | D
  deriving DecidableEq, Repr

structure DominantBlock where
  kind : ClassicalKind
  rank : Nat
  large : match kind with
    | .A => 7 ≤ rank
    | .B | .D => 9 ≤ rank

/-- The A20 bound supersedes, but is compatible with, five sevenths. -/
def aCharge (r : Nat) : Nat :=
  if r = 20 then 550 else if 17 ≤ r then 500 else 350

/-- The three prefix charges are deliberately weakened rational bounds:
0.88 at rank19, 0.90 at ranks20/21, and 0.94 at ranks22/23. -/
def classicalCharge (r : Nat) : Nat :=
  if 22 ≤ r then 658 else if 20 ≤ r then 630
  else if 19 ≤ r then 616 else 450

def charge (b : DominantBlock) : Nat :=
  match b.kind with
  | .A => aCharge b.rank
  | .B | .D => classicalCharge b.rank

def rankSum : List DominantBlock → Nat
  | [] => 0
  | b :: bs => b.rank + rankSum bs

def chargeSum : List DominantBlock → Nat
  | [] => 0
  | b :: bs => charge b + chargeSum bs

def withinLocalCap (b : DominantBlock) : Prop :=
  match b.kind with
  | .A => b.rank ≤ 20
  | .B => b.rank ≤ 22
  | .D => b.rank ≤ 23

/-- Ordinary rank includes paid small components and any type-J pair
whose two ranks have already been funded together. `shortA1` counts
only the unfunded isolated type I. `sevenDefects` counts either allowed
isolated rank-two seven defect. Construction of this ledger from an
actual farthest forest is an explicit, unformalized geometric input. -/
structure ForestLedger (n : Nat) where
  ordinaryRank : Nat
  shortA1 : Nat
  sevenDefects : Nat
  dominants : List DominantBlock
  rankBalance : n = ordinaryRank + shortA1 + 2 * sevenDefects + rankSum dominants
  oneDebit : shortA1 + sevenDefects ≤ 1
  budget : 50 * ordinaryRank + 63 * sevenDefects + chargeSum dominants < 700
  localCaps : ∀ b, b ∈ dominants → b.rank ≤ 29 → withinLocalCap b

def spectatorRank {n : Nat} (l : ForestLedger n) : Nat :=
  l.ordinaryRank + l.shortA1 + 2 * l.sevenDefects

theorem aCharge_half (r : Nat) : 350 ≤ aCharge r := by
  by_cases h20 : r = 20 <;> by_cases h17 : 17 ≤ r <;>
    simp [aCharge, h20, h17]

theorem aCharge_five_sevenths (r : Nat) (large : 17 ≤ r) :
    500 ≤ aCharge r := by
  by_cases h20 : r = 20 <;> simp [aCharge, h20, large]

theorem classicalCharge_nine_fourteenths (r : Nat) :
    450 ≤ classicalCharge r := by
  by_cases h22 : 22 ≤ r <;> by_cases h20 : 20 ≤ r <;>
    by_cases h19 : 19 ≤ r <;> simp [classicalCharge, h22, h20, h19]

theorem classicalCharge_prefix19 (r : Nat) (large : 19 ≤ r) :
    616 ≤ classicalCharge r := by
  by_cases h22 : 22 ≤ r <;> by_cases h20 : 20 ≤ r <;>
    simp [classicalCharge, h22, h20, large]

theorem classicalCharge_prefix20 (r : Nat) (large : 20 ≤ r) :
    630 ≤ classicalCharge r := by
  by_cases h22 : 22 ≤ r <;> simp [classicalCharge, h22, large]

theorem charge_half (b : DominantBlock) : 350 ≤ charge b := by
  cases hk : b.kind with
  | A => simpa [charge, hk] using aCharge_half b.rank
  | B =>
    have h := classicalCharge_nine_fourteenths b.rank
    simp only [charge, hk]
    omega
  | D =>
    have h := classicalCharge_nine_fourteenths b.rank
    simp only [charge, hk]
    omega

/-- Uniqueness follows from the budget; it is not a classification input. -/
theorem at_most_one_dominant {n : Nat} (l : ForestLedger n) :
    l.dominants.length ≤ 1 := by
  cases hs : l.dominants with
  | nil => simp
  | cons b bs =>
    cases bs with
    | nil => simp
    | cons c cs =>
      have hb := charge_half b
      have hc := charge_half c
      have budget := l.budget
      rw [hs] at budget
      simp only [chargeSum] at budget
      omega

/-- Existence also follows from the ordinary/debit budget. -/
theorem unique_dominant_of_large_dimension {n : Nat} (l : ForestLedger n)
    (large : 24 ≤ n) : ∃ b, l.dominants = [b] := by
  have lengthBound := at_most_one_dominant l
  cases hs : l.dominants with
  | nil =>
    have ranks := l.rankBalance
    have budget := l.budget
    have debit := l.oneDebit
    rw [hs] at ranks budget
    simp only [rankSum, chargeSum, Nat.add_zero] at ranks budget
    omega
  | cons b bs =>
    cases bs with
    | nil => exact ⟨b, rfl⟩
    | cons c cs =>
      simp only [hs, List.length_cons] at lengthBound
      omega

/-- This conservative common debit estimate covers both alternatives:
an isolated short A1 loses 50 units; a rank-two seven defect loses only
37. Keeping 50 for either is sufficient for the final proof. -/
theorem singleton_accounting {n : Nat} (l : ForestLedger n)
    (b : DominantBlock) (only : l.dominants = [b]) :
    n = b.rank + spectatorRank l ∧ charge b + 50 * spectatorRank l < 750 := by
  have ranks := l.rankBalance
  have budget := l.budget
  have debit := l.oneDebit
  rw [only] at ranks budget
  simp only [rankSum, chargeSum, Nat.add_zero] at ranks budget
  unfold spectatorRank
  constructor <;> omega

theorem singleton_cap {n : Nat} (l : ForestLedger n) (b : DominantBlock)
    (only : l.dominants = [b]) (published29 : n ≤ 29) : withinLocalCap b := by
  have ranks := (singleton_accounting l b only).1
  apply l.localCaps b
  · simp [only]
  · omega

theorem single_A_impossible {n : Nat} (l : ForestLedger n) (b : DominantBlock)
    (only : l.dominants = [b]) (large : 24 ≤ n) (published29 : n ≤ 29)
    (kind : b.kind = .A) : False := by
  obtain ⟨ranks, budget⟩ := singleton_accounting l b only
  have cap : b.rank ≤ 20 := by
    simpa [withinLocalCap, kind] using singleton_cap l b only published29
  have half := charge_half b
  have spectator7 : spectatorRank l ≤ 7 := by omega
  have r17 : 17 ≤ b.rank := by omega
  have raised : 500 ≤ charge b := by
    simpa [charge, kind] using aCharge_five_sevenths b.rank r17
  have spectator4 : spectatorRank l ≤ 4 := by omega
  have r20 : b.rank = 20 := by omega
  have finalCharge : charge b = 550 := by simp [charge, kind, aCharge, r20]
  omega

theorem single_B_impossible {n : Nat} (l : ForestLedger n) (b : DominantBlock)
    (only : l.dominants = [b]) (large : 24 ≤ n) (published29 : n ≤ 29)
    (kind : b.kind = .B) : False := by
  obtain ⟨ranks, budget⟩ := singleton_accounting l b only
  have cap : b.rank ≤ 22 := by
    simpa [withinLocalCap, kind] using singleton_cap l b only published29
  have floor : 450 ≤ charge b := by
    simpa [charge, kind] using classicalCharge_nine_fourteenths b.rank
  have spectator5 : spectatorRank l ≤ 5 := by omega
  have r19 : 19 ≤ b.rank := by omega
  by_cases h22 : 22 ≤ b.rank
  · have paid : charge b = 658 := by simp [charge, kind, classicalCharge, h22]
    omega
  · by_cases h20 : 20 ≤ b.rank
    · have paid : 630 ≤ charge b := by
        simpa [charge, kind] using classicalCharge_prefix20 b.rank h20
      omega
    · have paid : 616 ≤ charge b := by
        simpa [charge, kind] using classicalCharge_prefix19 b.rank r19
      omega

/-- The surviving D pattern is derived before invoking its geometric
coupled obstruction. No dimension24 frontier is supplied as a premise. -/
theorem single_D_frontier {n : Nat} (l : ForestLedger n) (b : DominantBlock)
    (only : l.dominants = [b]) (large : 24 ≤ n) (published29 : n ≤ 29)
    (kind : b.kind = .D) :
    n = 24 ∧ b.rank = 23 ∧ l.shortA1 = 1 ∧
      l.ordinaryRank = 0 ∧ l.sevenDefects = 0 := by
  obtain ⟨ranks, budget⟩ := singleton_accounting l b only
  have cap : b.rank ≤ 23 := by
    simpa [withinLocalCap, kind] using singleton_cap l b only published29
  have floor : 450 ≤ charge b := by
    simpa [charge, kind] using classicalCharge_nine_fourteenths b.rank
  have spectator5 : spectatorRank l ≤ 5 := by omega
  have r19 : 19 ≤ b.rank := by omega
  by_cases h22 : 22 ≤ b.rank
  · have paid : charge b = 658 := by simp [charge, kind, classicalCharge, h22]
    have r23 : b.rank = 23 := by omega
    have originalBudget := l.budget
    rw [only] at originalBudget
    simp only [chargeSum, Nat.add_zero] at originalBudget
    have debit := l.oneDebit
    unfold spectatorRank at ranks
    omega
  · by_cases h20 : 20 ≤ b.rank
    · have paid : 630 ≤ charge b := by
        simpa [charge, kind] using classicalCharge_prefix20 b.rank h20
      omega
    · have paid : 616 ≤ charge b := by
        simpa [charge, kind] using classicalCharge_prefix19 b.rank r19
      omega

/-- Conditional final ceiling. The last hypothesis is precisely the
local D23/isolated-type-I obstruction, not an ambient exclusion. -/
theorem ceiling23_from_local_component_bounds {n : Nat} (l : ForestLedger n)
    (published29 : n ≤ 29)
    (coupledD23 : ∀ b, l.dominants = [b] → b.kind = .D → b.rank = 23 →
      n = 24 → l.shortA1 = 1 → l.ordinaryRank = 0 → l.sevenDefects = 0 →
      False) : n ≤ 23 := by
  by_cases small : n ≤ 23
  · exact small
  · have large : 24 ≤ n := by omega
    obtain ⟨b, only⟩ := unique_dominant_of_large_dimension l large
    cases hk : b.kind with
    | A => exact (single_A_impossible l b only large published29 hk).elim
    | B => exact (single_B_impossible l b only large published29 hk).elim
    | D =>
      obtain ⟨n24, r23, short, ordinary, seven⟩ :=
        single_D_frontier l b only large published29 hk
      exact (coupledD23 b only hk r23 n24 short ordinary seven).elim

/-- A typed interpretation interface. Its fields are mathematical
obligations about actual compact polytopes; this file does not construct
their geometry or assert that a paper review discharges them in Lean. -/
structure GeometricBridge (P : Type) (dimension : P → Nat) (compact : P → Prop) where
  ledger : ∀ p, compact p → 24 ≤ dimension p → ForestLedger (dimension p)
  published29 : ∀ p, compact p → dimension p ≤ 29
  coupledD23 : ∀ p hp large b, (ledger p hp large).dominants = [b] →
    b.kind = .D → b.rank = 23 → dimension p = 24 →
    (ledger p hp large).shortA1 = 1 →
    (ledger p hp large).ordinaryRank = 0 →
    (ledger p hp large).sevenDefects = 0 → False

theorem compact_dimension_le23 {P : Type} {dimension : P → Nat} {compact : P → Prop}
    (g : GeometricBridge P dimension compact) (p : P) (hp : compact p) :
    dimension p ≤ 23 := by
  by_cases small : dimension p ≤ 23
  · exact small
  · have large : 24 ≤ dimension p := by omega
    exact ceiling23_from_local_component_bounds (g.ledger p hp large)
      (g.published29 p hp) (g.coupledD23 p hp large)

/-- The numeric ledger alone permits the final frontier; its actual
geometric exclusion is genuinely required by the conditional theorem. -/
theorem d23_short_numeric_frontier : 658 < 700 ∧ 23 + 1 = 24 := by decide

namespace RealAccounting

open Lean.Grind
universe u
variable {R : Type u} [CommRing R] [LE R] [LT R]
  [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [OrderedRing R]

attribute [local instance] Semiring.natCast Ring.intCast

/-- A position is retained for every component, so repeated isomorphic
blocks need not have equal actual real energies. -/
def energySum : List (DominantBlock × R) → R
  | [] => 0
  | z :: zs => z.2 + energySum zs

theorem dominant_lower_sum (zs : List (DominantBlock × R))
    (lower : ∀ z ∈ zs, (charge z.1 : R) ≤ 700*z.2) :
    (chargeSum (zs.map Prod.fst) : R) ≤ 700*energySum zs := by
  induction zs with
  | nil => simp [chargeSum, energySum, Semiring.natCast_zero, Semiring.mul_zero]
  | cons z zs ih =>
    have head := lower z (by simp)
    have tail := ih (by intro w hw; exact lower w (by simp [hw]))
    simp only [List.map_cons, chargeSum, energySum, Semiring.natCast_add]
    grind

/-- Actual energies may be arbitrary real scalars. Only local rational
lower bounds are integer-scaled; no energy-discreteness premise occurs.
The strict real budget comes from orthogonal projections of one vector.
Its geometric construction and the local estimates remain explicit. -/
structure RealLedgerWitness (n : Nat) where
  ordinaryRank : Nat
  shortA1 : Nat
  sevenDefects : Nat
  blocks : List (DominantBlock × R)
  ordinaryEnergy : R
  sevenEnergy : R
  rankBalance : n = ordinaryRank + shortA1 + 2*sevenDefects +
    rankSum (blocks.map Prod.fst)
  oneDebit : shortA1 + sevenDefects ≤ 1
  ordinaryLower : ((50*ordinaryRank : Nat) : R) ≤ 700*ordinaryEnergy
  sevenLower : ((63*sevenDefects : Nat) : R) ≤ 700*sevenEnergy
  blockLower : ∀ z ∈ blocks, (charge z.1 : R) ≤ 700*z.2
  energyBudget : ordinaryEnergy + sevenEnergy + energySum blocks < 1
  localCaps : ∀ b, b ∈ blocks.map Prod.fst → b.rank ≤ 29 → withinLocalCap b

/-- The integer budget is derived from the real inequalities, rather
than postulated or obtained by rounding actual component energies. -/
theorem strict_scaled_budget {n : Nat} (w : RealLedgerWitness (R := R) n) :
    50*w.ordinaryRank + 63*w.sevenDefects +
      chargeSum (w.blocks.map Prod.fst) < 700 := by
  have blocks := dominant_lower_sum w.blocks w.blockLower
  have ordinary := w.ordinaryLower
  have seven := w.sevenLower
  have total := w.energyBudget
  have scalar : ((50*w.ordinaryRank + 63*w.sevenDefects +
      chargeSum (w.blocks.map Prod.fst) : Nat) : R) < 700 := by
    simp only [Semiring.natCast_add]
    grind
  apply OrderedRing.lt_of_natCast_lt_natCast (R := R)
    (50*w.ordinaryRank + 63*w.sevenDefects + chargeSum (w.blocks.map Prod.fst)) 700
  simpa only [Semiring.natCast_eq_ofNat] using scalar

def RealLedgerWitness.toLedger {n : Nat} (w : RealLedgerWitness (R := R) n) :
    ForestLedger n := {
  ordinaryRank := w.ordinaryRank
  shortA1 := w.shortA1
  sevenDefects := w.sevenDefects
  dominants := w.blocks.map Prod.fst
  rankBalance := w.rankBalance
  oneDebit := w.oneDebit
  budget := strict_scaled_budget w
  localCaps := w.localCaps
}

theorem ceiling23_from_real_component_bounds {n : Nat}
    (w : RealLedgerWitness (R := R) n) (published29 : n ≤ 29)
    (coupledD23 : ∀ b, w.blocks.map Prod.fst = [b] →
      b.kind = .D → b.rank = 23 → n = 24 → w.shortA1 = 1 →
      w.ordinaryRank = 0 → w.sevenDefects = 0 → False) : n ≤ 23 := by
  exact ceiling23_from_local_component_bounds w.toLedger published29 coupledD23

end RealAccounting

end CoxeterProofs.CompactFarthestComposition23

#print axioms CoxeterProofs.CompactFarthestComposition23.at_most_one_dominant
#print axioms CoxeterProofs.CompactFarthestComposition23.unique_dominant_of_large_dimension
#print axioms CoxeterProofs.CompactFarthestComposition23.single_D_frontier
#print axioms CoxeterProofs.CompactFarthestComposition23.ceiling23_from_local_component_bounds
#print axioms CoxeterProofs.CompactFarthestComposition23.compact_dimension_le23
#print axioms CoxeterProofs.CompactFarthestComposition23.RealAccounting.dominant_lower_sum
#print axioms CoxeterProofs.CompactFarthestComposition23.RealAccounting.strict_scaled_budget
#print axioms CoxeterProofs.CompactFarthestComposition23.RealAccounting.ceiling23_from_real_component_bounds
