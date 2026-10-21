import Std

/-! Source draft: finite termination and absolute-surplus lifting for the
rank-57 low-eight concentration process. Zero ranks are harmless padded
factors: they contribute zero mass/potential and a factor one to link size.

This file does not prove the real convex-order comparison or identify
these lists with actual horospherical links. Its normalized comparison
is an explicit, denominator-cleared premise. No computation or compiled
verification is claimed merely by the presence of this source. -/

namespace CoxeterProofs.ConcentrationDescent58

def potential (p : List Nat) : Nat := (p.map (fun r => r*r)).sum
def linkSize (p : List Nat) : Nat := (p.map (fun r => r+1)).prod

/-- A complete low-eight move, up to permutation of factors. -/
def Step (p q : List Nat) : Prop :=
  ∃ a b rs, 1 ≤ a ∧ a ≤ b ∧ b < 8 ∧
    p.Perm (a::b::rs) ∧ q.Perm ((a-1)::(b+1)::rs)

def Terminal (p : List Nat) : Prop := ∀ q, ¬ Step p q

/-- Only positive ranks below eight can participate in a move. -/
def lowFactor (r : Nat) : Bool := decide (1 ≤ r ∧ r < 8)
def lowCount (p : List Nat) : Nat := (p.filter lowFactor).length

theorem step_has_two_low_factors {p q : List Nat} (step : Step p q) :
    2 ≤ lowCount p := by
  obtain ⟨a,b,rs,positive,ordered,small,hp,hq⟩ := step
  have ha : lowFactor a = true := by simp [lowFactor] <;> omega
  have hb : lowFactor b = true := by simp [lowFactor] <;> omega
  have length := (hp.filter lowFactor).length_eq
  change 2 ≤ (p.filter lowFactor).length
  rw [length]
  simp [ha, hb]

/-- Filtering and its complementary filter expose both factors with
their multiplicities. They need not have distinct rank values. -/
theorem exists_step_of_two_low_factors (p : List Nat) (many : 2 ≤ lowCount p) :
    ∃ q, Step p q := by
  have partition := (List.filter_append_perm lowFactor p).symm
  cases hl : p.filter lowFactor with
  | nil => simp [lowCount, hl] at many
  | cons a tail =>
    cases tail with
    | nil => simp [lowCount, hl] at many
    | cons b rest =>
      have ha : 1 ≤ a ∧ a < 8 := by
        have member : a ∈ p.filter lowFactor := by rw [hl]; simp
        have hit := (List.mem_filter.mp member).2
        simpa [lowFactor] using hit
      have hb : 1 ≤ b ∧ b < 8 := by
        have member : b ∈ p.filter lowFactor := by rw [hl]; simp
        have hit := (List.mem_filter.mp member).2
        simpa [lowFactor] using hit
      let rs := rest ++ p.filter (fun r => !lowFactor r)
      have front : p.Perm (a :: b :: rs) := by simpa [hl, rs] using partition
      by_cases ordered : a ≤ b
      · exact ⟨(a-1)::(b+1)::rs,
          a,b,rs,ha.1,ordered,hb.2,front,List.Perm.refl _⟩
      · exact ⟨(b-1)::(a+1)::rs,
          b,a,rs,hb.1,(by omega),ha.2,
          front.trans (List.Perm.swap b a rs),List.Perm.refl _⟩

/-- Exact stopping-domain recognition, including arbitrary zero padding. -/
theorem terminal_iff_lowCount (p : List Nat) : Terminal p ↔ lowCount p ≤ 1 := by
  constructor
  · intro terminal
    by_cases small : lowCount p ≤ 1
    · exact small
    · obtain ⟨q, step⟩ := exists_step_of_two_low_factors p (by omega)
      exact (terminal q step).elim
  · intro small q step
    have many := step_has_two_low_factors step
    omega

theorem positive_filter_sum (p : List Nat) :
    (p.filter (fun r => decide (0 < r))).sum = p.sum := by
  induction p with
  | nil => simp
  | cons r rs ih =>
    by_cases positive : 0 < r
    · simp [positive, ih]
    · have zero : r = 0 := by omega
      simp [zero, ih]

theorem linkSize_cons (r : Nat) (rs : List Nat) :
    linkSize (r :: rs) = (r+1) * linkSize rs := by
  simp only [linkSize, List.map_cons, List.prod_cons]

theorem linkSize_perm {p q : List Nat} (perm : p.Perm q) :
    linkSize p = linkSize q := by
  exact (List.Perm.map (fun r : Nat => r+1) perm).prod_nat

/-- Removing neutral zero factors preserves the actual product denominator. -/
theorem positive_filter_linkSize (p : List Nat) :
    linkSize (p.filter (fun r => decide (0 < r))) = linkSize p := by
  induction p with
  | nil => rfl
  | cons r rs ih =>
    by_cases positive : 0 < r
    · simp [positive, linkSize_cons, ih]
    · have zero : r = 0 := by omega
      simp [zero, linkSize_cons, ih]

/-- After neutral zero factors are removed, a terminal rank57 profile
is a sorted optional low factor followed by sorted ranks at least eight.
This is domain recognition, not completeness of any table enumeration. -/
theorem terminal_rank57_sorted_profile (p : List Nat) (mass : p.sum = 57)
    (terminal : Terminal p) :
    ∃ low high : List Nat,
      (p.filter (fun r => decide (0 < r))).Perm (low ++ high) ∧
      low.length ≤ 1 ∧ (∀ r ∈ low, 1 ≤ r ∧ r < 8) ∧
      (∀ r ∈ high, 8 ≤ r) ∧ (low ++ high).Pairwise (· ≤ ·) ∧
      (low ++ high).sum = 57 ∧ linkSize (low ++ high) = linkSize p := by
  let positive := p.filter (fun r => decide (0 < r))
  let low := positive.filter lowFactor
  let rawHigh := positive.filter (fun r => !lowFactor r)
  let high := rawHigh.mergeSort (fun a b => decide (a ≤ b))
  have lowFilter : positive.filter lowFactor = p.filter lowFactor := by
    dsimp [positive]
    rw [List.filter_filter]
    apply congrArg (fun f : Nat → Bool => p.filter f)
    funext r
    by_cases pos : 0 < r
    · simp [pos]
    · have notLow : ¬ (1 ≤ r ∧ r < 8) := by omega
      simp [lowFactor, pos, notLow]
  have short : low.length ≤ 1 := by
    have bound := (terminal_iff_lowCount p).mp terminal
    simpa [low, lowFilter, lowCount] using bound
  have lowBounds : ∀ r ∈ low, 1 ≤ r ∧ r < 8 := by
    intro r member
    have hit := (List.mem_filter.mp member).2
    simpa [lowFactor] using hit
  have highPerm : high.Perm rawHigh := List.mergeSort_perm _ _
  have highBounds : ∀ r ∈ high, 8 ≤ r := by
    intro r member
    have hit := List.mem_filter.mp (highPerm.mem_iff.mp member)
    have pos := (List.mem_filter.mp hit.1).2
    have notLow := hit.2
    simp [lowFactor] at notLow
    have positiveRank : 0 < r := by simpa using pos
    omega
  have highSorted : high.Pairwise (· ≤ ·) := by
    have sorted := List.pairwise_mergeSort
      (le := fun a b : Nat => decide (a ≤ b))
      (by intro a b c hab hbc; simp_all <;> omega)
      (by intro a b; simp <;> omega) rawHigh
    exact sorted.imp (by intro a b h; simpa using h)
  have lowSorted : low.Pairwise (· ≤ ·) := by
    cases hl : low with
    | nil => simp
    | cons a tail =>
      cases tail with
      | nil => simp
      | cons b rest => simp only [hl, List.length_cons] at short; omega
  have perm : positive.Perm (low ++ high) :=
    (List.filter_append_perm lowFactor positive).symm.trans
      ((List.Perm.refl low).append highPerm.symm)
  refine ⟨low,high,perm,short,lowBounds,highBounds,?_,?_,?_⟩
  · exact List.pairwise_append.mpr ⟨lowSorted,highSorted,
      (by intro a ha b hb; have := lowBounds a ha; have := highBounds b hb; omega)⟩
  · have sums := perm.sum_nat
    have positiveMass : positive.sum = 57 := (positive_filter_sum p).trans mass
    omega
  · exact (linkSize_perm perm).symm.trans (positive_filter_linkSize p)

inductive Chain : Nat → List Nat → List Nat → Prop where
  | refl (p : List Nat) : Chain 0 p p
  | cons {p q r : List Nat} {k : Nat} :
      Step p q → Chain k q r → Chain (k+1) p r

theorem mass_le_potential (p : List Nat) : p.sum ≤ potential p := by
  induction p with
  | nil => simp [potential]
  | cons r rs ih =>
    have square : r ≤ r*r := by
      cases r with
      | zero => simp
      | succ r =>
        have h := Nat.mul_le_mul_left (r+1)
          (show 1 ≤ r+1 by omega)
        simpa using h
    simp only [potential, List.map_cons, List.sum_cons] at *
    omega

theorem potential_le_mass_square (p : List Nat) :
    potential p ≤ p.sum*p.sum := by
  induction p with
  | nil => simp [potential]
  | cons r rs ih =>
    simp only [potential, List.map_cons, List.sum_cons] at *
    simp only [Nat.add_mul, Nat.mul_add]
    omega

theorem linkSize_positive (p : List Nat) : 0 < linkSize p := by
  induction p with
  | nil => simp [linkSize]
  | cons r rs ih =>
    simp only [linkSize, List.map_cons, List.prod_cons]
    exact Nat.mul_pos (by omega) ih

theorem pair_potential_increase (a b : Nat) (positive : 1 ≤ a)
    (ordered : a ≤ b) :
    a*a+b*b+2 ≤ (a-1)*(a-1)+(b+1)*(b+1) := by
  cases a with
  | zero => omega
  | succ a =>
    simp only [Nat.succ_eq_add_one, Nat.add_sub_cancel,
      Nat.add_mul, Nat.mul_add, Nat.mul_one, Nat.one_mul]
    omega

theorem pair_size_decrease (a b : Nat) (positive : 1 ≤ a)
    (ordered : a ≤ b) :
    ((a-1)+1)*((b+1)+1) < (a+1)*(b+1) := by
  cases a with
  | zero => omega
  | succ a =>
    simp only [Nat.succ_eq_add_one, Nat.add_sub_cancel,
      Nat.add_mul, Nat.mul_add, Nat.mul_one, Nat.one_mul]
    omega

theorem step_mass (p q : List Nat) (step : Step p q) : q.sum = p.sum := by
  obtain ⟨a,b,rs,positive,ordered,small,hp,hq⟩ := step
  have ep := hp.sum_nat
  have eq := hq.sum_nat
  simp only [List.sum_cons] at ep eq
  omega

theorem step_potential (p q : List Nat) (step : Step p q) :
    potential p+2 ≤ potential q := by
  obtain ⟨a,b,rs,positive,ordered,small,hp,hq⟩ := step
  have ep := (List.Perm.map (fun r : Nat => r*r) hp).sum_nat
  have eq := (List.Perm.map (fun r : Nat => r*r) hq).sum_nat
  have inc := pair_potential_increase a b positive ordered
  simp only [potential, List.map_cons, List.sum_cons] at *
  omega

theorem step_linkSize (p q : List Nat) (step : Step p q) :
    linkSize q < linkSize p := by
  obtain ⟨a,b,rs,positive,ordered,small,hp,hq⟩ := step
  have ep := (List.Perm.map (fun r : Nat => r+1) hp).prod_nat
  have eq := (List.Perm.map (fun r : Nat => r+1) hq).prod_nat
  have strict := Nat.mul_lt_mul_of_pos_right
    (pair_size_decrease a b positive ordered) (linkSize_positive rs)
  simpa only [linkSize, List.map_cons, List.prod_cons,
    Nat.mul_assoc, ep, eq] using strict

theorem chain_mass {p q : List Nat} {k : Nat} (chain : Chain k p q) :
    q.sum = p.sum := by
  induction chain with
  | refl => rfl
  | cons first rest ih => exact ih.trans (step_mass _ _ first)

theorem chain_potential {p q : List Nat} {k : Nat}
    (chain : Chain k p q) : potential p+2*k ≤ potential q := by
  induction chain with
  | refl => omega
  | cons first rest ih =>
    have inc := step_potential _ _ first
    omega

theorem chain_linkSize {p q : List Nat} {k : Nat}
    (chain : Chain k p q) : linkSize q ≤ linkSize p := by
  induction chain with
  | refl => exact Nat.le_refl _
  | cons first rest ih =>
    exact Nat.le_trans ih (Nat.le_of_lt (step_linkSize _ _ first))

theorem chain_length_at_rank57 {p q : List Nat} {k : Nat}
    (mass : p.sum = 57) (chain : Chain k p q) : k ≤ 1596 := by
  have mq : q.sum = 57 := (chain_mass chain).trans mass
  have lower := mass_le_potential p
  have upper := potential_le_mass_square q
  have growth := chain_potential chain
  rw [mass] at lower
  rw [mq] at upper
  omega

/-- No infinite choice of moves is needed: bounded strong induction
produces a terminal state, with the same bound for every possible path. -/
theorem terminal_exists_at_rank57 (p : List Nat) (mass : p.sum = 57) :
    ∃ q k, Chain k p q ∧ Terminal q ∧ k ≤ 1596 := by
  classical
  have descent : ∀ fuel p, p.sum = 57 → 3249-potential p = fuel →
      ∃ q k, Chain k p q ∧ Terminal q ∧ k ≤ 1596 := by
    intro fuel
    induction fuel using Nat.strongRecOn with
    | ind fuel ih =>
      intro p mass measure
      by_cases next : ∃ q, Step p q
      · obtain ⟨q, step⟩ := next
        have mq : q.sum = 57 := (step_mass p q step).trans mass
        have upper := potential_le_mass_square q
        rw [mq] at upper
        have growth := step_potential p q step
        have smaller : 3249-potential q < fuel := by omega
        obtain ⟨r,k,tail,terminal,bound⟩ :=
          ih (3249-potential q) smaller q mq rfl
        have whole := Chain.cons step tail
        exact ⟨r,k+1,whole,terminal,chain_length_at_rank57 mass whole⟩
      · refine ⟨p,0,Chain.refl p,?_,by omega⟩
        intro q step
        exact next ⟨q,step⟩
  exact descent (3249-potential p) p mass rfl

/-- Clear a common charge denominator first. The normalized comparison
is (kp-cp)/vp >= (kq-cq)/vq, written without truncated subtraction.
Positive terminal surplus and link-size decrease give absolute surplus. -/
theorem absolute_surplus_lift (vp vq kp cp kq cq gap : Nat)
    (terminalSize : 0 < vq) (volume : vq ≤ vp)
    (normalised : vp*kq+vq*cp ≤ vq*kp+vp*cq)
    (terminalGap : cq+gap ≤ kq) : cp+gap ≤ kp := by
  have scaled := Nat.mul_le_mul_left vp terminalGap
  have gapMono := Nat.mul_le_mul_right gap volume
  simp only [Nat.mul_add] at scaled
  apply (Nat.mul_le_mul_left_iff terminalSize).mp
  rw [Nat.mul_add]
  omega

theorem chain_absolute_surplus {p q : List Nat} {k : Nat}
    (chain : Chain k p q) (kp cp kq cq gap : Nat)
    (normalised : linkSize p*kq+linkSize q*cp ≤
      linkSize q*kp+linkSize p*cq)
    (terminalGap : cq+gap ≤ kq) : cp+gap ≤ kp := by
  exact absolute_surplus_lift (linkSize p) (linkSize q) kp cp kq cq gap
    (linkSize_positive q) (chain_linkSize chain) normalised terminalGap

end CoxeterProofs.ConcentrationDescent58

#print axioms CoxeterProofs.ConcentrationDescent58.pair_potential_increase
#print axioms CoxeterProofs.ConcentrationDescent58.terminal_iff_lowCount
#print axioms CoxeterProofs.ConcentrationDescent58.terminal_rank57_sorted_profile
#print axioms CoxeterProofs.ConcentrationDescent58.linkSize_perm
#print axioms CoxeterProofs.ConcentrationDescent58.positive_filter_linkSize
#print axioms CoxeterProofs.ConcentrationDescent58.step_linkSize
#print axioms CoxeterProofs.ConcentrationDescent58.chain_length_at_rank57
#print axioms CoxeterProofs.ConcentrationDescent58.terminal_exists_at_rank57
#print axioms CoxeterProofs.ConcentrationDescent58.absolute_surplus_lift
#print axioms CoxeterProofs.ConcentrationDescent58.chain_absolute_surplus
