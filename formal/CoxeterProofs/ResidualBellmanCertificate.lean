import Std

/-!
Abstract kernel for the residual certificate in finite_volume57.tex,
Proposition `dp-certificate`. Compilation and printed-axiom audit status
are recorded in results/publication_logic_verification_20261020.json.

Ranks are positive natural numbers; weights are arbitrary integers (the
application uses the genuinely integral scale-64 residual weights). The
certificate supplies only local transition inequalities and tight predecessor
witnesses. Exactness and unreachability are conclusions, not certificate fields.
No geometry, alphabet-completeness, context-splitting or actual-table claim is
formalized here. Nothing is asserted about memo entries beyond the cutoff.
-/
namespace CoxeterProofs.ResidualBellmanCertificate

def rankSum {n : Nat} (itemRank : Fin n → Nat) : List (Fin n) → Nat
  | [] => 0
  | i :: xs => itemRank i + rankSum itemRank xs

def weightSum {n : Nat} (itemWeight : Fin n → Int) : List (Fin n) → Int
  | [] => 0
  | i :: xs => itemWeight i + weightSum itemWeight xs

structure Certificate (n cutoff : Nat) where
  itemRank : Fin n → Nat
  itemWeight : Fin n → Int
  memo : Nat → Option Int
  rank_pos : ∀ i, 0 < itemRank i
  base : memo 0 = some 0
  transition : ∀ (s : Nat) (v : Int) (i : Fin n),
    s + itemRank i ≤ cutoff → memo s = some v →
    ∃ w, memo (s + itemRank i) = some w ∧ v + itemWeight i ≤ w
  tight : ∀ (s : Nat) (v : Int),
    0 < s → s ≤ cutoff → memo s = some v →
    ∃ (t : Nat) (u : Int) (i : Fin n),
      t + itemRank i = s ∧ memo t = some u ∧ v = u + itemWeight i

variable {n cutoff : Nat}

/-- Transition dominance bounds every list, not merely a chosen witness path. -/
theorem list_upper (c : Certificate n cutoff) (xs : List (Fin n)) :
    rankSum c.itemRank xs ≤ cutoff →
    ∃ v, c.memo (rankSum c.itemRank xs) = some v ∧
      weightSum c.itemWeight xs ≤ v := by
  induction xs with
  | nil =>
    intro _
    exact ⟨0, c.base, by simp [weightSum]⟩
  | cons i xs ih =>
    intro capacity
    have tailCapacity : rankSum c.itemRank xs ≤ cutoff := by
      simp only [rankSum] at capacity
      omega
    obtain ⟨v, hv, bound⟩ := ih tailCapacity
    have nextCapacity : rankSum c.itemRank xs + c.itemRank i ≤ cutoff := by
      simp only [rankSum] at capacity
      omega
    obtain ⟨w, hw, upper⟩ :=
      c.transition (rankSum c.itemRank xs) v i nextCapacity hv
    refine ⟨w, ?_, ?_⟩
    · simpa only [rankSum, Nat.add_comm] using hw
    · simp only [weightSum]
      omega

/-- Strictly smaller predecessors give an attaining list by strong rank induction. -/
theorem memo_attained (c : Certificate n cutoff) (r : Nat) :
    r ≤ cutoff → ∀ v : Int, c.memo r = some v →
      ∃ xs : List (Fin n),
        rankSum c.itemRank xs = r ∧ weightSum c.itemWeight xs = v := by
  induction r using Nat.strongRecOn with
  | ind r ih =>
    intro capacity v state
    by_cases zero : r = 0
    · subst r
      have value : (0 : Int) = v := Option.some.inj (c.base.symm.trans state)
      exact ⟨[], rfl, value⟩
    · obtain ⟨t, u, i, split, previous, tight⟩ :=
        c.tight r v (by omega) capacity state
      have positive := c.rank_pos i
      have smaller : t < r := by omega
      have tailCapacity : t ≤ cutoff := by omega
      obtain ⟨xs, ranks, weights⟩ := ih t smaller tailCapacity u previous
      refine ⟨i :: xs, ?_, ?_⟩
      · simp only [rankSum]
        omega
      · simp only [weightSum]
        omega

def ExactMaximum (c : Certificate n cutoff) (r : Nat) (v : Int) : Prop :=
  (∃ xs : List (Fin n),
    rankSum c.itemRank xs = r ∧ weightSum c.itemWeight xs = v) ∧
  ∀ xs : List (Fin n), rankSum c.itemRank xs = r → weightSum c.itemWeight xs ≤ v

/-- Every finite memo value is the attained maximum over all lists of that rank. -/
theorem memo_exact (c : Certificate n cutoff) (r : Nat) (capacity : r ≤ cutoff)
    (v : Int) (state : c.memo r = some v) : ExactMaximum c r v := by
  refine ⟨memo_attained c r capacity v state, ?_⟩
  intro xs ranks
  have listCapacity : rankSum c.itemRank xs ≤ cutoff := by omega
  obtain ⟨w, hw, bound⟩ := list_upper c xs listCapacity
  have atRank : c.memo r = some w := by simpa only [ranks] using hw
  have same : w = v := Option.some.inj (atRank.symm.trans state)
  omega

/-- No extra reachable or exactness premise is needed in the converse. -/
theorem some_iff_exactMaximum (c : Certificate n cutoff) (r : Nat)
    (capacity : r ≤ cutoff) (v : Int) :
    c.memo r = some v ↔ ExactMaximum c r v := by
  constructor
  · exact memo_exact c r capacity v
  · intro exactValue
    obtain ⟨xs, ranks, weights⟩ := exactValue.1
    obtain ⟨w, hw, bound⟩ := list_upper c xs (by omega)
    have atRank : c.memo r = some w := by simpa only [ranks] using hw
    obtain ⟨ys, yrank, yweight⟩ := memo_attained c r capacity w atRank
    have reverseBound := exactValue.2 ys yrank
    have same : w = v := by omega
    simpa only [same] using atRank

/-- The `none` entry means there is no item list of that rank. -/
theorem none_iff_no_list (c : Certificate n cutoff) (r : Nat)
    (capacity : r ≤ cutoff) :
    c.memo r = none ↔ ¬ ∃ xs : List (Fin n), rankSum c.itemRank xs = r := by
  constructor
  · intro absent witness
    obtain ⟨xs, ranks⟩ := witness
    obtain ⟨v, hv, _⟩ := list_upper c xs (by omega)
    have atRank : c.memo r = some v := by simpa only [ranks] using hv
    rw [absent] at atRank
    cases atRank
  · intro absent
    cases state : c.memo r with
    | none => rfl
    | some v =>
      obtain ⟨xs, ranks, _⟩ := memo_attained c r capacity v state
      exact (absent ⟨xs, ranks⟩).elim

/-- The predecessor characterization is only for positive ranks; zero is the base. -/
theorem none_iff_no_predecessor (c : Certificate n cutoff) (r : Nat)
    (positive : 0 < r) (capacity : r ≤ cutoff) :
    c.memo r = none ↔ ¬ ∃ (t : Nat) (u : Int) (i : Fin n),
      t + c.itemRank i = r ∧ c.memo t = some u := by
  constructor
  · intro absent witness
    obtain ⟨t, u, i, split, previous⟩ := witness
    obtain ⟨w, hw, _⟩ := c.transition t u i (by omega) previous
    have atRank : c.memo r = some w := by simpa only [split] using hw
    rw [absent] at atRank
    cases atRank
  · intro absent
    cases state : c.memo r with
    | none => rfl
    | some v =>
      obtain ⟨t, u, i, split, previous, _⟩ := c.tight r v positive capacity state
      exact (absent ⟨t, u, i, split, previous⟩).elim

/-- Additive ranks and weights descend to multisets, represented by list permutation. -/
theorem totals_perm (c : Certificate n cutoff) {xs ys : List (Fin n)}
    (permutation : xs.Perm ys) :
    rankSum c.itemRank xs = rankSum c.itemRank ys ∧
      weightSum c.itemWeight xs = weightSum c.itemWeight ys := by
  induction permutation with
  | nil => exact ⟨rfl, rfl⟩
  | cons _ _ ih =>
    rcases ih with ⟨ranks, weights⟩
    simp only [rankSum, weightSum]
    constructor <;> omega
  | swap =>
    simp only [rankSum, weightSum]
    constructor <;> omega
  | trans _ _ ih₁ ih₂ =>
    exact ⟨ih₁.1.trans ih₂.1, ih₁.2.trans ih₂.2⟩

theorem length_le_rank (c : Certificate n cutoff) (xs : List (Fin n)) :
    xs.length ≤ rankSum c.itemRank xs := by
  induction xs with
  | nil => simp [rankSum]
  | cons i xs ih =>
    have positive := c.rank_pos i
    simp only [List.length_cons, rankSum]
    omega

/-- In the rank-58 application every finite state has a witness of at most 58 items. -/
theorem attained_within_58 (c : Certificate n cutoff) (cutoff_bound : cutoff ≤ 58)
    (r : Nat) (capacity : r ≤ cutoff) (v : Int) (state : c.memo r = some v) :
    ∃ xs : List (Fin n), rankSum c.itemRank xs = r ∧
      weightSum c.itemWeight xs = v ∧ xs.length ≤ 58 := by
  obtain ⟨xs, ranks, weights⟩ := memo_attained c r capacity v state
  have lengthBound := length_le_rank c xs
  exact ⟨xs, ranks, weights, by omega⟩

#print axioms list_upper
#print axioms memo_attained
#print axioms memo_exact
#print axioms some_iff_exactMaximum
#print axioms none_iff_no_list
#print axioms none_iff_no_predecessor
#print axioms totals_perm
#print axioms attained_within_58

end CoxeterProofs.ResidualBellmanCertificate
