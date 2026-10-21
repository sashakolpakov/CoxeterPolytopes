import CoxeterProofs.ArithmeticEightField19

/-!
Conditional logic for the independently reviewed actual-K7 extension.
The geometric consumer theorems, actual ground-field interpretation,
complete forest classification and h-vector realization are written
interfaces. This module checks all coupon/source budget arithmetic,
the ordinary sum, the exact endpoint surplus and the dimension logic.
-/
namespace CoxeterProofs.ArithmeticK7Field19

open CompactHVector ManuscriptDataBounds

inductive VertexFamily where
  | a20 | b20 | d20 | a19Isolate | a18TwoIsolates
  | b19Isolate | b18TwoIsolates | a2B17Isolate | a2B16TwoIsolates
  | d19Isolate | d18HighDihedral | twoComponentDE | other
  deriving DecidableEq, Repr

def allFamilies : List VertexFamily :=
  [.a20, .b20, .d20, .a19Isolate, .a18TwoIsolates,
   .b19Isolate, .b18TwoIsolates, .a2B17Isolate, .a2B16TwoIsolates,
   .d19Isolate, .d18HighDihedral, .twoComponentDE, .other]

/-- Quarter units; `otherCost` already includes the bounded active-A4
charge. Actual finite-forest coverage and the meaning of the labels
are supplied by the written composition, not inferred from this list. -/
def baseline : VertexFamily → Nat → Nat
  | .a20, _ => 148 | .b20, _ => 148 | .d20, _ => 152
  | .a19Isolate, _ => 140 | .a18TwoIsolates, _ => 132
  | .b19Isolate, _ => 140 | .b18TwoIsolates, _ => 132
  | .a2B17Isolate, _ => 128 | .a2B16TwoIsolates, _ => 120
  | .d19Isolate, _ => 144 | .d18HighDihedral, _ => 144
  | .twoComponentDE, _ => 144 | .other, cost => cost

def debit : VertexFamily → Nat
  | .a20 => 8 | .b20 => 8 | .d20 => 12 | .a19Isolate => 8
  | .b19Isolate => 4 | .d19Isolate => 4
  | .d18HighDihedral => 4 | .twoComponentDE => 4 | _ => 0

def receivingCap : VertexFamily → Nat
  | .a19Isolate => 4 | .a18TwoIsolates => 8
  | .b19Isolate => 3 | .b18TwoIsolates => 5
  | .a2B17Isolate => 8 | .a2B16TwoIsolates => 16 | _ => 0

def envelope (family : VertexFamily) (otherCost : Nat) : Nat :=
  baseline family otherCost - debit family + receivingCap family

theorem complete_family_labels (family : VertexFamily) : family ∈ allFamilies := by
  cases family <;> decide

theorem exact_quarter_budget_table :
    allFamilies.map (fun family => envelope family 140) =
      [140, 140, 140, 136, 140, 139, 137, 136, 136, 140, 140, 140, 140] := by decide

theorem envelope_upper (family : VertexFamily) (otherCost : Nat)
    (h : otherCost ≤ 140) : envelope family otherCost ≤ 140 := by
  cases family <;> simp_all [envelope, baseline, debit, receivingCap]

structure VertexEntry where
  family : VertexFamily
  otherCost : Nat
  otherAllowance : otherCost ≤ 140
  charge : Nat
  actualRoleBound : charge ≤ envelope family otherCost

theorem vertex_upper (entry : VertexEntry) : entry.charge ≤ 140 := by
  exact Nat.le_trans entry.actualRoleBound
    (envelope_upper entry.family entry.otherCost entry.otherAllowance)

/-- Quarter charges are actual sums of the simultaneous marked-angle
rule. The list contains every finite vertex. Geometry supplies the
actual polygon incidence lower and h-vector coefficient coverage. -/
structure Rank20Data where
  tail : List Term
  entries : List VertexEntry
  coefficientCoverage : ∀ t ∈ tail, ∃ i < 11,
    t.multiplicity = ArithmeticEightField19.multiplicity20 i ∧
    t.coefficient = ArithmeticEightField19.coefficient20 i
  completeVertexCount : entries.length =
    vertices (ArithmeticEightField19.endpoint20 :: tail)
  incidenceLower :
    4 * weighted (ArithmeticEightField19.endpoint20 :: tail) ≤
      (entries.map VertexEntry.charge).sum

theorem no_rank20_data (d : Rank20Data) : False := by
  have upper := sum_le_cap_mul_length (d.entries.map VertexEntry.charge) 140 (by
    intro charge member
    obtain ⟨entry, _, rfl⟩ := List.mem_map.mp member
    exact vertex_upper entry)
  rw [List.length_map, d.completeVertexCount] at upper
  have surplus := Nat.mul_le_mul_left 4
    (ArithmeticEightField19.endpoint_surplus20 d.tail d.coefficientCoverage)
  have lower := d.incidenceLower
  omega

/-- The geometry interface is required for actual K7, never mere
containment of K7. The previous all-field20 theorem supplies the tail. -/
theorem dimension_at_most19 (n : Nat) (earlierUniversal20 : n ≤ 20)
    (actualK7Geometry : n = 20 → Nonempty Rank20Data) : n ≤ 19 := by
  have excluded : n ≠ 20 := by
    intro equal
    obtain ⟨data⟩ := actualK7Geometry equal
    exact no_rank20_data data
  omega

end CoxeterProofs.ArithmeticK7Field19

#print axioms CoxeterProofs.ArithmeticK7Field19.complete_family_labels
#print axioms CoxeterProofs.ArithmeticK7Field19.exact_quarter_budget_table
#print axioms CoxeterProofs.ArithmeticK7Field19.envelope_upper
#print axioms CoxeterProofs.ArithmeticK7Field19.vertex_upper
#print axioms CoxeterProofs.ArithmeticK7Field19.no_rank20_data
#print axioms CoxeterProofs.ArithmeticK7Field19.dimension_at_most19
