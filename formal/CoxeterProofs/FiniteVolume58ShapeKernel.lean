import Std

/-!
Compilation and actual axiom-output verification are recorded separately;
source preparation alone is not verification.

Internal shape-curvature assembler for finite_volume57.tex, ambient 58.
Primitive data are index multiplicities and the below/above section split.
The two distinct extrema have opposite one-sided splits. Their two 841
gains are derived without symmetry or a degree-balance assumption.

Actual link realization, index distributions, uncut-edge interpretation,
global extrema and incidence identities remain geometric inputs. The
optional normalization adapter requires 2E=57V explicitly; in particular
this module never infers U=E for arbitrary all-above data.
-/

namespace CoxeterProofs.FiniteVolume58ShapeKernel

/-- All sums and curvatures are signed integers. -/
def sumN : Nat → (Nat → Int) → Int
  | 0, _ => 0
  | n + 1, f => sumN n f + f n

theorem sumN_congr (n : Nat) (f g : Nat → Int)
    (equal : ∀ i, i < n → f i = g i) : sumN n f = sumN n g := by
  induction n with
  | zero => rfl
  | succ n ih =>
    have previous := ih (fun i hi => equal i (by omega))
    have last := equal n (by omega)
    simp only [sumN, previous, last]

theorem sumN_add (n : Nat) (f g : Nat → Int) :
    sumN n (fun i => f i + g i) = sumN n f + sumN n g := by
  induction n with
  | zero => simp [sumN]
  | succ n ih => simp only [sumN, ih]; omega

theorem sumN_sub (n : Nat) (f g : Nat → Int) :
    sumN n (fun i => f i - g i) = sumN n f - sumN n g := by
  induction n with
  | zero => simp [sumN]
  | succ n ih => simp only [sumN, ih]; omega

theorem sumN_scale (n : Nat) (a : Int) (f : Nat → Int) :
    sumN n (fun i => a * f i) = a * sumN n f := by
  induction n with
  | zero => simp [sumN]
  | succ n ih => simp only [sumN, ih, Int.mul_add]

theorem sumN_linear (n : Nat) (a b : Int) (f g : Nat → Int) :
    sumN n (fun i => a * f i - b * g i) =
      a * sumN n f - b * sumN n g := by
  rw [sumN_sub, sumN_scale, sumN_scale]

theorem sumN_mono (n : Nat) (f g : Nat → Int)
    (bound : ∀ i, i < n → f i ≤ g i) : sumN n f ≤ sumN n g := by
  induction n with
  | zero => simp [sumN]
  | succ n ih =>
    have previous := ih (fun i hi => bound i (by omega))
    have last := bound n (by omega)
    simp only [sumN]
    omega

theorem sumN_zero (n : Nat) : sumN n (fun _ => 0) = 0 := by
  induction n with
  | zero => rfl
  | succ n ih => simp [sumN, ih]

def delta (k i : Nat) : Int := if i = k then 1 else 0

theorem sumN_delta_outside (n k : Nat) (outside : n ≤ k) :
    sumN n (fun i => delta k i) = 0 := by
  calc
    _ = sumN n (fun _ => 0) := sumN_congr n _ _ (by
      intro i hi
      have different : i ≠ k := by omega
      simp [delta, different])
    _ = 0 := sumN_zero n

theorem sumN_delta (n k : Nat) (inside : k < n) :
    sumN n (fun i => delta k i) = 1 := by
  induction n generalizing k with
  | zero => omega
  | succ n ih =>
    by_cases earlier : k < n
    · have previous := ih k earlier
      have different : n ≠ k := by omega
      simp only [sumN, delta, if_neg different, Int.add_zero]
      exact previous
    · have last : k = n := by omega
      subst k
      have previous := sumN_delta_outside n n (by omega)
      change sumN n (fun i => delta n i) + delta n n = 1
      rw [previous]
      simp [delta]

/-- Data outside 0..57 are unused. Positivity is needed only at extrema. -/
structure Link where
  h : Nat → Nat
  L : Nat → Nat
  R : Nat → Nat
  split : ∀ j, j < 58 → L j + R j = h j

def Positive (p : Link) : Prop := ∀ j, j < 58 → 1 ≤ p.h j
def AllBelow (p : Link) : Prop :=
  ∀ j, j < 58 → p.L j = p.h j ∧ p.R j = 0
def AllAbove (p : Link) : Prop :=
  ∀ j, j < 58 → p.L j = 0 ∧ p.R j = p.h j

def weight (j : Nat) : Int := (min j (57 - j) : Nat)
def belowGap (j : Nat) : Int := (j : Int) - weight j
def aboveGap (j : Nat) : Int := ((57 - j : Nat) : Int) - weight j

def G (p : Link) : Int := sumN 58 (fun j => weight j * (p.h j : Int))
def E (p : Link) : Int := sumN 58 (fun j => (j : Int) * (p.h j : Int))
def V (p : Link) : Int := sumN 58 (fun j => (p.h j : Int))
def U (p : Link) : Int := sumN 58 (fun j =>
  (j : Int) * (p.L j : Int) + ((57 - j : Nat) : Int) * (p.R j : Int))

def curvature2 (p : Link) : Int := 5 * G p - 2 * E p
def uncutCurvature2 (p : Link) : Int := 5 * U p - 2 * E p

theorem gap_totals :
    sumN 58 belowGap = 841 ∧ sumN 58 aboveGap = 841 := by
  decide

theorem finite_sums :
    sumN 58 weight = 812 ∧
    sumN 58 (fun j => (j : Int)) = 1653 ∧
    sumN 58 (fun _ => (1 : Int)) = 58 := by
  decide

theorem minimum_coefficients (j : Nat) :
    weight j ≤ (j : Int) ∧ weight j ≤ ((57 - j : Nat) : Int) := by
  have left := Nat.min_le_left j (57 - j)
  have right := Nat.min_le_right j (57 - j)
  unfold weight
  omega

theorem minimum_uncut (p : Link) : G p ≤ U p := by
  unfold G U
  apply sumN_mono
  intro j hj
  have split : (p.L j : Int) + (p.R j : Int) = (p.h j : Int) := by
    have raw := p.split j hj
    omega
  have left := Int.mul_le_mul_of_nonneg_right
    (minimum_coefficients j).1 (Int.ofNat_nonneg (p.L j))
  have right := Int.mul_le_mul_of_nonneg_right
    (minimum_coefficients j).2 (Int.ofNat_nonneg (p.R j))
  have expanded : weight j * (p.h j : Int) =
      weight j * (p.L j : Int) + weight j * (p.R j : Int) := by
    rw [← split, Int.mul_add]
  omega

/-- A positive multiplicity pays at least one copy of the coefficient gap. -/
theorem coefficient_gap (a b h : Int) (ordered : a ≤ b) (positive : 1 ≤ h) :
    a * h + (b - a) ≤ b * h := by
  have nonnegative : 0 ≤ b - a := by omega
  have product := Int.mul_le_mul_of_nonneg_left positive nonnegative
  rw [Int.mul_one, Int.sub_mul] at product
  omega

theorem edge_gap_from_positive (p : Link) (positive : Positive p) :
    G p + 841 ≤ E p := by
  calc
    G p + 841 = sumN 58 (fun j => weight j * (p.h j : Int) + belowGap j) := by
      rw [sumN_add, gap_totals.1]
      rfl
    _ ≤ E p := by
      unfold E
      apply sumN_mono
      intro j hj
      have pos : 1 ≤ (p.h j : Int) := by
        have raw := positive j hj
        omega
      exact coefficient_gap (weight j) (j : Int) (p.h j : Int)
        (minimum_coefficients j).1 pos

theorem all_below_gap (p : Link) (positive : Positive p) (below : AllBelow p) :
    G p + 841 ≤ U p := by
  have uncut : U p = E p := by
    unfold U E
    apply sumN_congr
    intro j hj
    have row := below j hj
    simp [row.1, row.2]
  rw [uncut]
  exact edge_gap_from_positive p positive

/-- No assertion U=E is used for an all-above link. -/
theorem all_above_gap (p : Link) (positive : Positive p) (above : AllAbove p) :
    G p + 841 ≤ U p := by
  calc
    G p + 841 = sumN 58 (fun j => weight j * (p.h j : Int) + aboveGap j) := by
      rw [sumN_add, gap_totals.2]
      rfl
    _ ≤ U p := by
      unfold U
      apply sumN_mono
      intro j hj
      have pos : 1 ≤ (p.h j : Int) := by
        have raw := positive j hj
        omega
      have bound := coefficient_gap (weight j) ((57 - j : Nat) : Int)
        (p.h j : Int) (minimum_coefficients j).2 pos
      have row := above j hj
      simpa [aboveGap, row.1, row.2] using bound

/-- Reserve exactly two distinct vertex positions, whether finite or ideal. -/
theorem two_extrema_curvature (n : Nat) (p : Nat → Link) (lo hi : Nat)
    (lo_lt : lo < n) (hi_lt : hi < n) (different : lo ≠ hi)
    (lo_positive : Positive (p lo)) (hi_positive : Positive (p hi))
    (minimum : AllAbove (p lo)) (maximum : AllBelow (p hi)) :
    sumN n (fun v => curvature2 (p v)) + 8410 ≤
      sumN n (fun v => uncutCurvature2 (p v)) := by
  have row : ∀ v, v < n →
      841 * delta lo v + 841 * delta hi v ≤ U (p v) - G (p v) := by
    intro v hv
    by_cases at_lo : v = lo
    · subst v
      have gain := all_above_gap (p lo) lo_positive minimum
      simp only [delta, if_pos rfl, if_neg different, Int.mul_one,
        Int.mul_zero, Int.add_zero]
      omega
    · by_cases at_hi : v = hi
      · subst v
        have other : hi ≠ lo := by omega
        have gain := all_below_gap (p hi) hi_positive maximum
        simp only [delta, if_neg other, if_pos rfl, Int.mul_zero,
          Int.mul_one, Int.zero_add]
        omega
      · have ordinary := minimum_uncut (p v)
        simp only [delta, if_neg at_lo, if_neg at_hi, Int.mul_zero,
          Int.add_zero]
        omega
  have total := sumN_mono n _ _ row
  have reservations : sumN n (fun v => 841 * delta lo v + 841 * delta hi v) =
      1682 := by
    rw [sumN_add, sumN_scale, sumN_scale,
      sumN_delta n lo lo_lt, sumN_delta n hi hi_lt]
    decide
  rw [reservations, sumN_sub] at total
  unfold curvature2 uncutCurvature2
  rw [sumN_linear, sumN_linear]
  omega

theorem finite_link_values (p : Link)
    (simplex : ∀ j, j < 58 → p.h j = 1) :
    G p = 812 ∧ E p = 1653 ∧ V p = 58 ∧ curvature2 p = 754 := by
  have g : G p = 812 := by
    calc
      G p = sumN 58 weight := sumN_congr 58 _ _ (by
        intro j hj
        simp [simplex j hj])
      _ = 812 := finite_sums.1
  have e : E p = 1653 := by
    calc
      E p = sumN 58 (fun j => (j : Int)) := sumN_congr 58 _ _ (by
        intro j hj
        simp [simplex j hj])
      _ = 1653 := finite_sums.2.1
  have v : V p = 58 := by
    calc
      V p = sumN 58 (fun _ => (1 : Int)) := sumN_congr 58 _ _ (by
        intro j hj
        simp [simplex j hj])
      _ = 58 := finite_sums.2.2
  refine ⟨g, e, v, ?_⟩
  unfold curvature2
  omega

def finiteCount (n : Nat) (finite : Nat → Bool) : Int :=
  sumN n (fun v => if finite v then 1 else 0)

def idealCurvature2 (n : Nat) (finite : Nat → Bool) (p : Nat → Link) : Int :=
  sumN n (fun v => if finite v then 0 else curvature2 (p v))

theorem finite_ideal_split (n : Nat) (finite : Nat → Bool) (p : Nat → Link)
    (simplex : ∀ v, v < n → finite v = true →
      ∀ j, j < 58 → (p v).h j = 1) :
    sumN n (fun v => curvature2 (p v)) =
      754 * finiteCount n finite + idealCurvature2 n finite p := by
  calc
    _ = sumN n (fun v => 754 * (if finite v then 1 else 0) +
          (if finite v then 0 else curvature2 (p v))) := sumN_congr n _ _ (by
      intro v hv
      cases tag : finite v with
      | false => simp [tag]
      | true =>
        have value := (finite_link_values (p v) (simplex v hv tag)).2.2.2
        simp [tag, value])
    _ = _ := by
      rw [sumN_add, sumN_scale]
      rfl

/-- Incidence inputs are I=ΣE and 2f2=ΣU, not the desired shape bound.
The conclusion is twice the manuscript's 377 vf + Σideal κ + 4205 bound. -/
theorem shape_lower_from_incidence (n : Nat) (p : Nat → Link)
    (finite : Nat → Bool) (lo hi : Nat)
    (lo_lt : lo < n) (hi_lt : hi < n) (different : lo ≠ hi)
    (lo_positive : Positive (p lo)) (hi_positive : Positive (p hi))
    (minimum : AllAbove (p lo)) (maximum : AllBelow (p hi))
    (simplex : ∀ v, v < n → finite v = true →
      ∀ j, j < 58 → (p v).h j = 1)
    (f2 I : Int) (edge_incidence : I = sumN n (fun v => E (p v)))
    (uncut_incidence : 2 * f2 = sumN n (fun v => U (p v))) :
    754 * finiteCount n finite + idealCurvature2 n finite p + 8410 ≤
      10 * f2 - 2 * I := by
  have shape := two_extrema_curvature n p lo hi lo_lt hi_lt different
    lo_positive hi_positive minimum maximum
  have split := finite_ideal_split n finite p simplex
  have incidence : sumN n (fun v => uncutCurvature2 (p v)) =
      10 * f2 - 2 * I := by
    unfold uncutCurvature2
    rw [sumN_linear]
    omega
  omega

/-- Degree balance is required only to identify the usual 57V curvature. -/
theorem balanced_curvature (p : Link) (balance : 2 * E p = 57 * V p) :
    curvature2 p = 5 * G p - 57 * V p := by
  unfold curvature2
  omega

end CoxeterProofs.FiniteVolume58ShapeKernel

#print axioms CoxeterProofs.FiniteVolume58ShapeKernel.gap_totals
#print axioms CoxeterProofs.FiniteVolume58ShapeKernel.minimum_uncut
#print axioms CoxeterProofs.FiniteVolume58ShapeKernel.edge_gap_from_positive
#print axioms CoxeterProofs.FiniteVolume58ShapeKernel.all_below_gap
#print axioms CoxeterProofs.FiniteVolume58ShapeKernel.all_above_gap
#print axioms CoxeterProofs.FiniteVolume58ShapeKernel.two_extrema_curvature
#print axioms CoxeterProofs.FiniteVolume58ShapeKernel.finite_link_values
#print axioms CoxeterProofs.FiniteVolume58ShapeKernel.shape_lower_from_incidence
#print axioms CoxeterProofs.FiniteVolume58ShapeKernel.balanced_curvature
