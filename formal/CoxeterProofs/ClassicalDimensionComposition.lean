import Std

/-!
Prepared, uncompiled conditional logic for the written compact19 and
first-type arithmetic17 composition. No actual geometric instance is
constructed here. All predicates concern ONE fixed candidate polytope;
the witness passed through a proof is its SAME actual farthest vertex.

The Nat variables in the first theorem are positive integer field norms,
not real coordinates, heights or curvature. The remaining Nat variables
are a dimension and an actual field degree. Only Std is imported.
-/

namespace CoxeterProofs.ClassicalDimensionComposition

/-- The quartic norm arithmetic, after the written field-norm bridge. -/
theorem quartic_norm_contradiction
    (a b n : Nat) (a_positive : 0 < a) (b_positive : 0 < b)
    (n_positive : 0 < n) (n_cap : n ≤ 6)
    (first_divisibility : 256 ∣ a * n)
    (second_divisibility : 16 ∣ b * n)
    (head_norm_cap : 9 * (a * b) < 616) : False := by
  obtain ⟨i, first_witness⟩ := first_divisibility
  obtain ⟨j, second_witness⟩ := second_divisibility
  have first_positive : 0 < a * n := Nat.mul_pos a_positive n_positive
  have second_positive : 0 < b * n := Nat.mul_pos b_positive n_positive
  have first_lower : 256 ≤ a * n := by omega
  have second_lower : 16 ≤ b * n := by omega
  have product_lower : 4096 ≤ (a * b) * (n * n) := by
    calc
      4096 = 256 * 16 := by decide
      _ ≤ (a * n) * (b * n) := Nat.mul_le_mul first_lower second_lower
      _ = a * (n * (b * n)) := Nat.mul_assoc a n (b * n)
      _ = a * (b * (n * n)) :=
        congrArg (fun k : Nat => a * k) (Nat.mul_left_comm n b n)
      _ = (a * b) * (n * n) := (Nat.mul_assoc a b (n * n)).symm
  have square_cap : n * n ≤ 36 := Nat.mul_le_mul n_cap n_cap
  have scaled_cap : (a * b) * (n * n) ≤ 36 * (a * b) := by
    calc
      _ ≤ (a * b) * 36 := Nat.mul_le_mul_left (a * b) square_cap
      _ = 36 * (a * b) := Nat.mul_comm (a * b) 36
  have lower : 4096 ≤ 36 * (a * b) := Nat.le_trans product_lower scaled_cap
  omega

/-- The complete residual degree alternatives after Meyer and degree<=4. -/
theorem degree_two_three_four_cases
    (d : Nat) (lower : 2 ≤ d) (upper : d ≤ 4) :
    d = 2 ∨ d = 3 ∨ d = 4 := by
  omega

universe u v

/-- Conditional assembly from the prior global20 tail and written branches. -/
theorem compact_dimension_at_most_nineteen
    {Vertex : Type u} (dimension : Nat)
    (ActualFarthest WholeA20 WholeB20 WholeD20 : Vertex → Prop)
    (IsolatedNativeFiveHead IsolatedFiniteLargeHead IsolatedDottedHead :
      Vertex → Prop)
    (writtenPriorGlobal20 : dimension ≤ 20)
    (writtenActualForestAt20 : dimension = 20 →
      ∃ p, ActualFarthest p ∧ (WholeA20 p ∨ WholeB20 p ∨ WholeD20 p))
    (writtenWholeA20Exclusion : dimension = 20 →
      ∀ p, ActualFarthest p → WholeA20 p → False)
    (writtenWholeB20Exclusion : dimension = 20 →
      ∀ p, ActualFarthest p → WholeB20 p → False)
    (writtenD20ActualHeadCases : dimension = 20 →
      ∀ p, ActualFarthest p → WholeD20 p →
        IsolatedNativeFiveHead p ∨ IsolatedFiniteLargeHead p ∨ IsolatedDottedHead p)
    (writtenNativeFiveExclusion : dimension = 20 →
      ∀ p, ActualFarthest p → WholeD20 p → IsolatedNativeFiveHead p → False)
    (writtenFiniteLargeExclusion : dimension = 20 →
      ∀ p, ActualFarthest p → WholeD20 p → IsolatedFiniteLargeHead p → False)
    (writtenDottedExclusion : dimension = 20 →
      ∀ p, ActualFarthest p → WholeD20 p → IsolatedDottedHead p → False) :
    dimension ≤ 19 := by
  by_cases bound : dimension ≤ 19
  · exact bound
  · have at20 : dimension = 20 := by omega
    obtain ⟨p, farthest, forest⟩ := writtenActualForestAt20 at20
    have impossible : False := by
      cases forest with
      | inl a20 => exact writtenWholeA20Exclusion at20 p farthest a20
      | inr bd =>
        cases bd with
        | inl b20 => exact writtenWholeB20Exclusion at20 p farthest b20
        | inr d20 =>
          cases writtenD20ActualHeadCases at20 p farthest d20 with
          | inl nativeFive =>
            exact writtenNativeFiveExclusion at20 p farthest d20 nativeFive
          | inr finiteOrDotted =>
            cases finiteOrDotted with
            | inl finiteLarge =>
              exact writtenFiniteLargeExclusion at20 p farthest d20 finiteLarge
            | inr dotted =>
              exact writtenDottedExclusion at20 p farthest d20 dotted
    exact False.elim impossible

/--
Conditional first-type assembly for every actual field. TypedI predicates
must mean the actual isolated short A1 of the written frontier. They are
not arbitrary ordinary A1 spectators. The quartic branch supplies norm
data, so its integer contradiction is derived by the first theorem.
-/
theorem arithmetic_dimension_at_most_seventeen
    {Vertex : Type u} {Field : Type v} (dimension : Nat)
    (ActualFarthest TypedA17I TypedB17I TypedD17I WholeD18 : Vertex → Prop)
    (RemainingFiniteD18Head IsolatedDottedD18Head : Vertex → Prop)
    (ActualGroundField : Vertex → Field → Prop) (fieldDegree : Field → Nat)
    (writtenPriorGlobal18 : dimension ≤ 18)
    (writtenActualForestAt18 : dimension = 18 → ∃ p, ActualFarthest p ∧
      (TypedA17I p ∨ TypedB17I p ∨ TypedD17I p ∨ WholeD18 p))
    (writtenA17TypedIExclusion : dimension = 18 →
      ∀ p, ActualFarthest p → TypedA17I p → False)
    (writtenB17TypedIExclusion : dimension = 18 →
      ∀ p, ActualFarthest p → TypedB17I p → False)
    (writtenD17TypedIExclusion : dimension = 18 →
      ∀ p, ActualFarthest p → TypedD17I p → False)
    (writtenD18ActualHeadCases : dimension = 18 →
      ∀ p, ActualFarthest p → WholeD18 p →
        RemainingFiniteD18Head p ∨ IsolatedDottedD18Head p)
    (writtenFiniteHeadExclusion : dimension = 18 →
      ∀ p, ActualFarthest p → WholeD18 p → RemainingFiniteD18Head p → False)
    (writtenActualDottedFieldDegree : dimension = 18 →
      ∀ p, ActualFarthest p → WholeD18 p → IsolatedDottedD18Head p →
        ∃ k, ActualGroundField p k ∧ 2 ≤ fieldDegree k ∧ fieldDegree k ≤ 4)
    (writtenQuadraticExclusion : dimension = 18 →
      ∀ p k, ActualFarthest p → WholeD18 p → IsolatedDottedD18Head p →
        ActualGroundField p k → fieldDegree k = 2 → False)
    (writtenCubicExclusion : dimension = 18 →
      ∀ p k, ActualFarthest p → WholeD18 p → IsolatedDottedD18Head p →
        ActualGroundField p k → fieldDegree k = 3 → False)
    (writtenQuarticNormData : dimension = 18 →
      ∀ p k, ActualFarthest p → WholeD18 p → IsolatedDottedD18Head p →
        ActualGroundField p k → fieldDegree k = 4 →
          ∃ a b n : Nat, 0 < a ∧ 0 < b ∧ 0 < n ∧ n ≤ 6 ∧
            256 ∣ a * n ∧ 16 ∣ b * n ∧ 9 * (a * b) < 616) :
    dimension ≤ 17 := by
  by_cases bound : dimension ≤ 17
  · exact bound
  · have at18 : dimension = 18 := by omega
    obtain ⟨p, farthest, forest⟩ := writtenActualForestAt18 at18
    have impossible : False := by
      cases forest with
      | inl a17i => exact writtenA17TypedIExclusion at18 p farthest a17i
      | inr rest =>
        cases rest with
        | inl b17i => exact writtenB17TypedIExclusion at18 p farthest b17i
        | inr rest =>
          cases rest with
          | inl d17i => exact writtenD17TypedIExclusion at18 p farthest d17i
          | inr d18 =>
            cases writtenD18ActualHeadCases at18 p farthest d18 with
            | inl finiteHead =>
              exact writtenFiniteHeadExclusion at18 p farthest d18 finiteHead
            | inr dotted =>
              obtain ⟨k, actualField, degreeLower, degreeUpper⟩ :=
                writtenActualDottedFieldDegree at18 p farthest d18 dotted
              cases degree_two_three_four_cases (fieldDegree k) degreeLower degreeUpper with
              | inl quadratic =>
                exact writtenQuadraticExclusion at18 p k farthest d18 dotted
                  actualField quadratic
              | inr cubicOrQuartic =>
                cases cubicOrQuartic with
                | inl cubic =>
                  exact writtenCubicExclusion at18 p k farthest d18 dotted actualField cubic
                | inr quartic =>
                  obtain ⟨a, b, n, ha, hb, hn, cap, first, second, headCap⟩ :=
                    writtenQuarticNormData at18 p k farthest d18 dotted actualField quartic
                  exact quartic_norm_contradiction a b n ha hb hn cap first second headCap
    exact False.elim impossible

end CoxeterProofs.ClassicalDimensionComposition

#print axioms CoxeterProofs.ClassicalDimensionComposition.quartic_norm_contradiction
#print axioms CoxeterProofs.ClassicalDimensionComposition.degree_two_three_four_cases
#print axioms CoxeterProofs.ClassicalDimensionComposition.compact_dimension_at_most_nineteen
#print axioms CoxeterProofs.ClassicalDimensionComposition.arithmetic_dimension_at_most_seventeen
