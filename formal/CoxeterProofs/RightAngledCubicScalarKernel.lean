import Std
import Init.Grind.Ordered
import CoxeterProofs.RightAngledSharpDimensions

/-!
Source-bound compilation and axiom-audit status are recorded in
formal/.proof-batches. The coordinate draft required elaboration repairs
before its first successful kernel check.

This successor does not modify the frozen RightAngledSharpDimensions module.
It replaces its two conclusion-shaped cubic premises by finite scalar data.
The lower-level interface uses the Gram diagonal/strict gap, three positive
pair sums, explicit tensor contractions and four-point Gram positivity.
The preferred coordinate constructor at the end DERIVES all these positivity
facts by finite Gram-square identities. Its only geometric/algebraic inputs
are unit six-coordinate vectors, the strict gap, and the two pointwise
contraction identities for the 36- and 216-coordinate tensor arrays.

The scalar type is an arbitrary linearly ordered commutative ring, not Rat.
Thus rational scalars are one possible instance, not an integer/rational proxy
for arbitrary real spherical codes. Unit-vector/tensor realization and the
two pointwise contraction identities still need a real-geometry bridge.
No desired cardinality, dimension bound, N=27 endpoint implication, or cubic
sum-of-squares inequality is a field of the input structures below.
-/

namespace CoxeterProofs.RightAngledCubicScalarKernel

open Lean.Grind

universe u
variable {R : Type u} [CommRing R] [LE R] [LT R]
  [Std.IsLinearOrder R] [Std.LawfulOrderLT R] [OrderedRing R]

attribute [local instance] Semiring.natCast Ring.intCast

/-- A transparent finite sum; no Finset or real-analysis library is required. -/
def sumN : Nat → (Nat → R) → R
  | 0, _ => 0
  | n + 1, f => sumN n f + f n

def pairSum (n : Nat) (f : Nat → Nat → R) : R :=
  sumN n (fun i => sumN n (f i))

theorem sumN_congr (n : Nat) (f g : Nat → R)
    (h : ∀ i, i < n → f i = g i) : sumN n f = sumN n g := by
  induction n with
  | zero => rfl
  | succ n ih =>
    have lower := ih (fun i hi => h i (by omega))
    have last := h n (by omega)
    simp only [sumN, lower, last]

theorem sumN_add (n : Nat) (f g : Nat → R) :
    sumN n (fun i => f i + g i) = sumN n f + sumN n g := by
  induction n with
  | zero => simp [sumN, Semiring.add_zero]
  | succ n ih => simp only [sumN, ih]; grind

theorem sumN_sub (n : Nat) (f g : Nat → R) :
    sumN n (fun i => f i - g i) = sumN n f - sumN n g := by
  induction n with
  | zero => simp [sumN, Ring.sub_eq_add_neg, AddCommGroup.neg_zero, Semiring.add_zero]
  | succ n ih => simp only [sumN, ih]; grind

theorem sumN_scale (n : Nat) (a : R) (f : Nat → R) :
    sumN n (fun i => a * f i) = a * sumN n f := by
  induction n with
  | zero => simp [sumN, Semiring.mul_zero]
  | succ n ih => simp only [sumN, ih]; grind

theorem sumN_const (n : Nat) (a : R) :
    sumN n (fun _ => a) = (n : R) * a := by
  induction n with
  | zero => simp [sumN, Semiring.natCast_zero, Semiring.zero_mul]
  | succ n ih => simp only [sumN, ih, Semiring.natCast_succ]; grind

theorem pairSum_congr (n : Nat) (f g : Nat → Nat → R)
    (h : ∀ i, i < n → ∀ j, j < n → f i j = g i j) :
    pairSum n f = pairSum n g := by
  apply sumN_congr
  intro i hi
  exact sumN_congr n (f i) (g i) (h i hi)

theorem pairSum_add (n : Nat) (f g : Nat → Nat → R) :
    pairSum n (fun i j => f i j + g i j) = pairSum n f + pairSum n g := by
  unfold pairSum
  simp only [sumN_add]

theorem pairSum_sub (n : Nat) (f g : Nat → Nat → R) :
    pairSum n (fun i j => f i j - g i j) = pairSum n f - pairSum n g := by
  unfold pairSum
  simp only [sumN_sub]

theorem pairSum_scale (n : Nat) (a : R) (f : Nat → Nat → R) :
    pairSum n (fun i j => a * f i j) = a * pairSum n f := by
  unfold pairSum
  simp only [sumN_scale]

theorem pairSum_const (n : Nat) (a : R) :
    pairSum n (fun _ _ => a) = (n : R) * (n : R) * a := by
  unfold pairSum
  rw [sumN_const, sumN_const]
  grind

theorem sumN_nonpos (n : Nat) (f : Nat → R)
    (h : ∀ i, i < n → f i ≤ 0) : sumN n f ≤ 0 := by
  induction n with
  | zero => simp [sumN]
  | succ n ih =>
    have lower := ih (fun i hi => h i (by omega))
    have last := h n (by omega)
    simp only [sumN]
    grind

theorem sumN_zero_terms (n : Nat) (f : Nat → R)
    (h : ∀ i, i < n → f i ≤ 0) (zero : sumN n f = 0) :
    ∀ i, i < n → f i = 0 := by
  induction n with
  | zero => intro i hi; omega
  | succ n ih =>
    have lowerLe := sumN_nonpos n f (fun i hi => h i (by omega))
    have lastLe := h n (by omega)
    have lowerZero : sumN n f = 0 := by simp only [sumN] at zero; grind
    have lastZero : f n = 0 := by simp only [sumN] at zero; grind
    intro i hi
    by_cases equal : i = n
    · simpa [equal] using lastZero
    · exact ih (fun j hj => h j (by omega)) lowerZero i (by omega)

theorem sumN_delta (n i : Nat) (a : R) :
    sumN n (fun j => if i = j then a else 0) = if i < n then a else 0 := by
  induction n with
  | zero => simp [sumN]
  | succ n ih =>
    simp only [sumN, ih]
    by_cases equal : i = n
    · subst i; simp [AddCommMonoid.zero_add]
    · by_cases before : i < n
      · have next : i < n + 1 := by omega
        simp [equal, before, next, Semiring.add_zero]
      · have next : ¬ i < n + 1 := by omega
        simp [equal, before, next, Semiring.add_zero]

/-- Sixteen times the manuscript cubic, avoiding all division notation. -/
def cubic16 (x : R) : R := 16 * (x * x * x) + 12 * (x * x) - 1

theorem cubic16_factor (x : R) :
    cubic16 x = (2 * x + 1)^2 * (4 * x - 1) := by
  grind [cubic16]

theorem cubic16_nonpos (x : R) (strict : 4 * x < 1) : cubic16 x ≤ 0 := by
  have square : 0 ≤ (2 * x + 1)^2 := OrderedRing.sq_nonneg
  have negative : 4 * x - 1 < 0 := by grind
  have product := OrderedRing.mul_nonpos_of_nonneg_of_nonpos square
    (Preorder.le_of_lt negative)
  rw [cubic16_factor]
  exact product

theorem cubic16_zero_root (x : R) (strict : 4 * x < 1)
    (zero : cubic16 x = 0) : 2 * x + 1 = 0 := by
  by_cases root : 2 * x + 1 = 0
  · exact root
  · have square : 0 < (2 * x + 1)^2 := OrderedRing.sq_pos root
    have negative : 4 * x - 1 < 0 := by grind
    have product := OrderedRing.mul_neg_of_pos_of_neg square negative
    rw [← cubic16_factor] at product
    grind

/-- Real geometry supplies a unit Gram table and this one PSD consequence.
Four-point positivity is unconditional in the indices, not conditional on
N=27. Symmetry is true geometrically but is not needed by this kernel. -/
structure ScalarCode (n : Nat) where
  gram : Nat → Nat → R
  diagonal : ∀ i, i < n → gram i i = 1
  strict : ∀ i, i < n → ∀ j, j < n → i ≠ j → 4 * gram i j < 1
  fourPSD : ∀ i j k l, i < n → j < n → k < n → l < n →
    0 ≤ 4 + 2 * (gram i j + gram i k + gram i l +
      gram j k + gram j l + gram k l)

/-- Explicit remaining Euclidean tensor bridge. `rankTwo` and `rankThree`
are the inner-product kernels of A(u) and T(u), respectively. Positivity
means the squared norm of their sums, not the desired cubic inequality. -/
structure TensorContractions {n : Nat} (d : ScalarCode (R := R) n) where
  rankTwo : Nat → Nat → R
  rankThree : Nat → Nat → R
  linearPSD : 0 ≤ pairSum n d.gram
  rankTwoPSD : 0 ≤ pairSum n rankTwo
  rankThreePSD : 0 ≤ pairSum n rankThree
  rankTwoIdentity : ∀ i, i < n → ∀ j, j < n →
    6 * rankTwo i j = 6 * (d.gram i j * d.gram i j) - 1
  rankThreeIdentity : ∀ i, i < n → ∀ j, j < n →
    8 * rankThree i j = 8 * (d.gram i j * d.gram i j * d.gram i j) -
      3 * d.gram i j

theorem cubic_pair_sum_formula (n : Nat) (g : Nat → Nat → R) :
    pairSum n (fun i j => cubic16 (g i j)) =
      16 * pairSum n (fun i j => g i j * g i j * g i j) +
      12 * pairSum n (fun i j => g i j * g i j) - (n : R) * (n : R) := by
  simp only [cubic16, pairSum_sub, pairSum_add, pairSum_scale, pairSum_const]
  grind

/-- The full scalar SOS equality, derived from the pointwise contractions. -/
theorem cubic_sum_of_squares {n : Nat} (d : ScalarCode (R := R) n)
    (t : TensorContractions d) :
    pairSum n (fun i j => cubic16 (d.gram i j)) =
      (n : R) * (n : R) + 6 * pairSum n d.gram +
        12 * pairSum n t.rankTwo + 16 * pairSum n t.rankThree := by
  have two := pairSum_congr n
    (fun i j => 6 * t.rankTwo i j)
    (fun i j => 6 * (d.gram i j * d.gram i j) - 1)
    t.rankTwoIdentity
  have three := pairSum_congr n
    (fun i j => 8 * t.rankThree i j)
    (fun i j => 8 * (d.gram i j * d.gram i j * d.gram i j) - 3 * d.gram i j)
    t.rankThreeIdentity
  simp only [pairSum_scale, pairSum_sub, pairSum_const] at two three
  have formula := cubic_pair_sum_formula n d.gram
  grind

theorem cubic_lower {n : Nat} (d : ScalarCode (R := R) n)
    (t : TensorContractions d) :
    (n : R) * (n : R) ≤ pairSum n (fun i j => cubic16 (d.gram i j)) := by
  have identity := cubic_sum_of_squares d t
  have one := t.linearPSD
  have two := t.rankTwoPSD
  have three := t.rankThreePSD
  grind

def offCubic {n : Nat} (d : ScalarCode (R := R) n) (i j : Nat) : R :=
  if i = j then 0 else cubic16 (d.gram i j)

theorem offCubic_nonpos {n : Nat} (d : ScalarCode (R := R) n)
    (i : Nat) (hi : i < n) (j : Nat) (hj : j < n) : offCubic d i j ≤ 0 := by
  by_cases equal : i = j
  · simp [offCubic, equal]
  · simpa [offCubic, equal] using cubic16_nonpos _ (d.strict i hi j hj equal)

theorem cubic_diagonal_split {n : Nat} (d : ScalarCode (R := R) n) :
    pairSum n (fun i j => cubic16 (d.gram i j)) =
      27 * (n : R) + pairSum n (offCubic d) := by
  have entries := pairSum_congr n
    (fun i j => cubic16 (d.gram i j))
    (fun i j => (if i = j then 27 else 0) + offCubic d i j) (by
      intro i hi j hj
      by_cases equal : i = j
      · subst j
        rw [d.diagonal i hi]
        simp [offCubic, cubic16]
        grind
      · simp [offCubic, equal, AddCommMonoid.zero_add])
  have diag : pairSum n (fun i j => if i = j then (27 : R) else 0) =
      27 * (n : R) := by
    unfold pairSum
    have rows := sumN_congr n
      (fun i => sumN n (fun j => if i = j then (27 : R) else 0))
      (fun _ => (27 : R)) (by
        intro i hi
        rw [sumN_delta]
        simp [hi])
    rw [sumN_const] at rows
    grind
  rw [pairSum_add, diag] at entries
  exact entries

theorem off_total_nonpos {n : Nat} (d : ScalarCode (R := R) n) :
    pairSum n (offCubic d) ≤ 0 := by
  apply sumN_nonpos
  intro i hi
  exact sumN_nonpos n (offCubic d i) (offCubic_nonpos d i hi)

/-- The count inequality is now a theorem, not a geometric input. -/
theorem tensor_count {n : Nat} (d : ScalarCode (R := R) n)
    (t : TensorContractions d) : n * n ≤ 27 * n := by
  have lower := cubic_lower d t
  have split := cubic_diagonal_split d
  have off := off_total_nonpos d
  have scalar : (n : R) * (n : R) ≤ 27 * (n : R) := by grind
  apply OrderedRing.le_of_natCast_le_natCast (R := R) (n * n) (27 * n)
  simp only [Semiring.natCast_mul]
  simpa only [Semiring.natCast_eq_ofNat] using scalar

/-- Equality forces every distinct pair into the double root of the cubic. -/
theorem endpoint_half_products {n : Nat} (d : ScalarCode (R := R) n)
    (t : TensorContractions d) (endpoint : n = 27) :
    ∀ i, i < n → ∀ j, j < n → i ≠ j → 2 * d.gram i j + 1 = 0 := by
  have lower := cubic_lower d t
  have split := cubic_diagonal_split d
  have off := off_total_nonpos d
  have scalarEquality : (n : R) * (n : R) = 27 * (n : R) := by
    simp only [endpoint, Semiring.natCast_eq_ofNat]
  have totalZero : pairSum n (offCubic d) = 0 := by grind
  have rows := sumN_zero_terms n (fun i => sumN n (offCubic d i))
    (fun i hi => sumN_nonpos n (offCubic d i) (offCubic_nonpos d i hi)) totalZero
  intro i hi j hj different
  have entry := sumN_zero_terms n (offCubic d i)
    (offCubic_nonpos d i hi) (rows i hi) j hj
  have polynomialZero : cubic16 (d.gram i j) = 0 := by
    simpa [offCubic, different] using entry
  exact cubic16_zero_root _ (d.strict i hi j hj different) polynomialZero

theorem no_endpoint_twenty_seven {n : Nat} (d : ScalarCode (R := R) n)
    (t : TensorContractions d) (endpoint : n = 27) : False := by
  have half := endpoint_half_products d t endpoint
  have p01 := half 0 (by omega) 1 (by omega) (by omega)
  have p02 := half 0 (by omega) 2 (by omega) (by omega)
  have p03 := half 0 (by omega) 3 (by omega) (by omega)
  have p12 := half 1 (by omega) 2 (by omega) (by omega)
  have p13 := half 1 (by omega) 3 (by omega) (by omega)
  have p23 := half 2 (by omega) 3 (by omega) (by omega)
  have positive := d.fourPSD 0 1 2 3 (by omega) (by omega) (by omega) (by omega)
  have zeroLessOne : (0 : R) < 1 := OrderedRing.zero_lt_one
  grind

/-- Both old consequence fields are derived; neither is a premise. -/
theorem scalar_consequences {n : Nat} (d : ScalarCode (R := R) n)
    (t : TensorContractions d) : RightAngledSharpDimensions.StrictCodeSixConsequences n := by
  exact {
    tensorCount := tensor_count d t
    equalityFourGram := fun endpoint => False.elim (no_endpoint_twenty_seven d t endpoint)
  }

/-- Reuse the frozen natural-number composition instead of duplicating it. -/
theorem strict_code_count_le26 {n : Nat} (d : ScalarCode (R := R) n)
    (t : TensorContractions d) : n ≤ 26 :=
  RightAngledSharpDimensions.strict_code_count_le26 n (scalar_consequences d t)

/-- Interface to the frozen dimension theorem. The bridge must produce the
actual finite real scalar Gram/tensor data; it does not supply a count bound. -/
def cubic_geometry_from_scalar_bridge {Face : Type} {Point : Type}
    (d : RightAngledSharpDimensions.GeometricData Face Point)
    (bridge : ∀ F q, d.radialCodeSix F q →
      ∃ s : ScalarCode (R := R) (d.facetCount F), Nonempty (TensorContractions s)) :
    RightAngledSharpDimensions.WrittenCubicGeometry d := by
  constructor
  intro F q code
  obtain ⟨s, ⟨t⟩⟩ := bridge F q code
  exact scalar_consequences s t

/-! Finite coordinate witnesses: all Gram positivity is now derived. -/

theorem sumN_nonneg (n : Nat) (f : Nat → R)
    (h : ∀ i, i < n → 0 ≤ f i) : 0 ≤ sumN n f := by
  induction n with
  | zero => simp [sumN]
  | succ n ih =>
    have lower := ih (fun i hi => h i (by omega))
    have last := h n (by omega)
    simp only [sumN]
    grind

theorem sumN_scale_right (n : Nat) (a : R) (f : Nat → R) :
    sumN n (fun i => f i * a) = sumN n f * a := by
  have entries := sumN_congr n (fun i => f i * a) (fun i => a * f i)
    (by intro i hi; grind)
  rw [sumN_scale] at entries
  grind

/-- Finite Fubini, valid also for either empty index set. -/
theorem sumN_swap (m n : Nat) (f : Nat → Nat → R) :
    sumN m (fun i => sumN n (f i)) =
      sumN n (fun j => sumN m (fun i => f i j)) := by
  induction m with
  | zero =>
    simp only [sumN, sumN_const, Semiring.mul_zero]
  | succ m ih =>
    change sumN m (fun i => sumN n (f i)) + sumN n (f m) =
      sumN n (fun j => sumN m (fun i => f i j) + f m j)
    rw [sumN_add, ih]

def coordinateDot (k : Nat) (u v : Nat → R) : R :=
  sumN k (fun a => u a * v a)

theorem coordinateDot_self_nonneg (k : Nat) (v : Nat → R) :
    0 ≤ coordinateDot k v v := by
  apply sumN_nonneg
  intro a ha
  have square : 0 ≤ (v a)^2 := OrderedRing.sq_nonneg
  simpa only [Semiring.pow_two] using square

/-- The all-ones Gram quadratic form is the squared norm of a vector sum.
This is a finite algebraic identity, not a positivity assumption. -/
theorem coordinate_gram_square (n k : Nat) (v : Nat → Nat → R) :
    pairSum n (fun i j => coordinateDot k (v i) (v j)) =
      sumN k (fun a =>
        sumN n (fun i => v i a) * sumN n (fun i => v i a)) := by
  unfold pairSum coordinateDot
  have rows : ∀ i,
      sumN n (fun j => sumN k (fun a => v i a * v j a)) =
        sumN k (fun a => v i a * sumN n (fun j => v j a)) := by
    intro i
    rw [sumN_swap n k (fun j a => v i a * v j a)]
    apply sumN_congr
    intro a ha
    exact sumN_scale n (v i a) (fun j => v j a)
  calc
    sumN n (fun i => sumN n (fun j => sumN k (fun a => v i a * v j a))) =
        sumN n (fun i => sumN k (fun a => v i a * sumN n (fun j => v j a))) :=
      sumN_congr n _ _ (fun i _ => rows i)
    _ = sumN k (fun a => sumN n (fun i => v i a * sumN n (fun j => v j a))) :=
      sumN_swap n k (fun i a => v i a * sumN n (fun j => v j a))
    _ = sumN k (fun a => sumN n (fun i => v i a) * sumN n (fun i => v i a)) := by
      apply sumN_congr
      intro a ha
      exact sumN_scale_right n (sumN n (fun i => v i a)) (fun i => v i a)

theorem coordinate_gram_nonneg (n k : Nat) (v : Nat → Nat → R) :
    0 ≤ pairSum n (fun i j => coordinateDot k (v i) (v j)) := by
  rw [coordinate_gram_square]
  apply sumN_nonneg
  intro a ha
  have square : 0 ≤ (sumN n (fun i => v i a))^2 := OrderedRing.sq_nonneg
  simpa only [Semiring.pow_two] using square

/-- Expand four coordinate vectors directly. Symmetry is supplied by the
commutative scalar multiplication, including repeated vector indices. -/
theorem coordinateDot_four_square (k : Nat) (u v w z : Nat → R) :
    coordinateDot k (fun a => u a + v a + w a + z a)
      (fun a => u a + v a + w a + z a) =
      coordinateDot k u u + coordinateDot k v v + coordinateDot k w w +
        coordinateDot k z z + 2 *
          (coordinateDot k u v + coordinateDot k u w + coordinateDot k u z +
           coordinateDot k v w + coordinateDot k v z + coordinateDot k w z) := by
  have entries := sumN_congr k
    (fun a => (u a + v a + w a + z a) * (u a + v a + w a + z a))
    (fun a => u a * u a + v a * v a + w a * w a + z a * z a +
      2 * (u a * v a + u a * w a + u a * z a +
        v a * w a + v a * z a + w a * z a))
    (by intro a ha; grind)
  simp only [sumN_add, sumN_scale] at entries
  exact entries

/-- Preferred finite-coordinate input. There is NO PSD, SOS, endpoint,
cardinality or dimension conclusion in this structure. A and T are given
as arrays of all 36 and 216 ordinary Euclidean tensor coordinates.
Their two displayed contractions are the only remaining tensor identities. -/
structure CoordinateWitness (n : Nat) where
  vector : Nat → Nat → R
  rankTwoCoordinates : Nat → Nat → R
  rankThreeCoordinates : Nat → Nat → R
  unit : ∀ i, i < n → coordinateDot 6 (vector i) (vector i) = 1
  strict : ∀ i, i < n → ∀ j, j < n → i ≠ j →
    4 * coordinateDot 6 (vector i) (vector j) < 1
  rankTwoIdentity : ∀ i, i < n → ∀ j, j < n →
    6 * coordinateDot 36 (rankTwoCoordinates i) (rankTwoCoordinates j) =
      6 * (coordinateDot 6 (vector i) (vector j) *
        coordinateDot 6 (vector i) (vector j)) - 1
  rankThreeIdentity : ∀ i, i < n → ∀ j, j < n →
    8 * coordinateDot 216 (rankThreeCoordinates i) (rankThreeCoordinates j) =
      8 * (coordinateDot 6 (vector i) (vector j) *
        coordinateDot 6 (vector i) (vector j) *
        coordinateDot 6 (vector i) (vector j)) -
      3 * coordinateDot 6 (vector i) (vector j)

def coordinateScalarCode {n : Nat} (w : CoordinateWitness (R := R) n) :
    ScalarCode (R := R) n where
  gram i j := coordinateDot 6 (w.vector i) (w.vector j)
  diagonal := w.unit
  strict := w.strict
  fourPSD := by
    intro i j k l hi hj hk hl
    have nonnegative := coordinateDot_self_nonneg 6
      (fun a => w.vector i a + w.vector j a + w.vector k a + w.vector l a)
    have expansion := coordinateDot_four_square 6
      (w.vector i) (w.vector j) (w.vector k) (w.vector l)
    have ii := w.unit i hi
    have jj := w.unit j hj
    have kk := w.unit k hk
    have ll := w.unit l hl
    grind

def coordinateTensorContractions {n : Nat} (w : CoordinateWitness (R := R) n) :
    TensorContractions (coordinateScalarCode w) where
  rankTwo i j := coordinateDot 36 (w.rankTwoCoordinates i) (w.rankTwoCoordinates j)
  rankThree i j := coordinateDot 216 (w.rankThreeCoordinates i) (w.rankThreeCoordinates j)
  linearPSD := coordinate_gram_nonneg n 6 w.vector
  rankTwoPSD := coordinate_gram_nonneg n 36 w.rankTwoCoordinates
  rankThreePSD := coordinate_gram_nonneg n 216 w.rankThreeCoordinates
  rankTwoIdentity := w.rankTwoIdentity
  rankThreeIdentity := w.rankThreeIdentity

theorem coordinate_code_count_le26 {n : Nat} (w : CoordinateWitness (R := R) n) :
    n ≤ 26 :=
  strict_code_count_le26 (coordinateScalarCode w) (coordinateTensorContractions w)

/-- The real bridge now supplies only actual unit coordinate vectors, their
strict pair gaps and two tensor contraction identities. All finite Gram
positivity and the old cubic consequences are derived in this source. -/
def cubic_geometry_from_coordinate_bridge {Face : Type} {Point : Type}
    (d : RightAngledSharpDimensions.GeometricData Face Point)
    (bridge : ∀ F q, d.radialCodeSix F q →
      Nonempty (CoordinateWitness (R := R) (d.facetCount F))) :
    RightAngledSharpDimensions.WrittenCubicGeometry d := by
  apply cubic_geometry_from_scalar_bridge (R := R) d
  intro F q code
  obtain ⟨w⟩ := bridge F q code
  exact ⟨coordinateScalarCode w, ⟨coordinateTensorContractions w⟩⟩

end CoxeterProofs.RightAngledCubicScalarKernel

#print axioms CoxeterProofs.RightAngledCubicScalarKernel.cubic_sum_of_squares
#print axioms CoxeterProofs.RightAngledCubicScalarKernel.tensor_count
#print axioms CoxeterProofs.RightAngledCubicScalarKernel.endpoint_half_products
#print axioms CoxeterProofs.RightAngledCubicScalarKernel.no_endpoint_twenty_seven
#print axioms CoxeterProofs.RightAngledCubicScalarKernel.strict_code_count_le26
#print axioms CoxeterProofs.RightAngledCubicScalarKernel.cubic_geometry_from_scalar_bridge
#print axioms CoxeterProofs.RightAngledCubicScalarKernel.coordinate_gram_square
#print axioms CoxeterProofs.RightAngledCubicScalarKernel.coordinateTensorContractions
#print axioms CoxeterProofs.RightAngledCubicScalarKernel.coordinate_code_count_le26
#print axioms CoxeterProofs.RightAngledCubicScalarKernel.cubic_geometry_from_coordinate_bridge
