import Std

/-!
The group-theoretic closure of the two completed H9 host exclusions.

The concrete lattice isometries, their exact matrix certificates, full
actual corner germs, and the standard chamber and finite-order trace
facts are explicit inputs. This module proves the intervening logic:
trace excess excludes a positive period; powers are distinct; their
chamber stabilizers have distinct reflection-group cosets.

No particular host, geometric theorem, matrix identity or complete root
enumeration is claimed to be formalized here. All group and action laws
are visible in GroupActionModel, since this project uses Std, not Mathlib.
-/

namespace CoxeterProofs.ChamberPeriod

universe u v w

structure GroupActionModel (G : Type u) (X : Type v) where
  one : G
  mul : G → G → G
  inv : G → G
  assoc : ∀ a b c, mul (mul a b) c = mul a (mul b c)
  one_mul : ∀ a, mul one a = a
  mul_one : ∀ a, mul a one = a
  inv_mul : ∀ a, mul (inv a) a = one
  mul_inv : ∀ a, mul a (inv a) = one
  act : G → X → X
  act_one : ∀ x, act one x = x
  act_mul : ∀ a b x, act (mul a b) x = act a (act b x)

def power {G : Type u} {X : Type v} (d : GroupActionModel G X)
    (g : G) : Nat → G
  | 0 => d.one
  | n + 1 => d.mul (power d g n) g

theorem left_cancel {G : Type u} {X : Type v}
    (d : GroupActionModel G X) (a b c : G)
    (equal : d.mul a b = d.mul a c) : b = c := by
  have h := congrArg (fun z => d.mul (d.inv a) z) equal
  simpa only [← d.assoc, d.inv_mul, d.one_mul] using h

theorem power_add {G : Type u} {X : Type v}
    (d : GroupActionModel G X) (g : G) (m n : Nat) :
    power d g (m + n) = d.mul (power d g m) (power d g n) := by
  induction n with
  | zero => simp only [Nat.add_zero, power, d.mul_one]
  | succ n ih =>
    change d.mul (power d g (m + n)) g =
      d.mul (power d g m) (d.mul (power d g n) g)
    rw [ih, d.assoc]

/-- The finite-order trace bound is a standard external matrix fact.
The concrete strict excess is supplied by the reviewed exact certificate
and its real-embedding argument. Neither is hidden as a Lean axiom. -/
structure TraceInterface {G : Type u} {X : Type v}
    (d : GroupActionModel G X) where
  traceAtMostDimension : G → Prop
  finiteOrderBound : ∀ g k, 0 < k → power d g k = d.one →
    traceAtMostDimension g

theorem no_positive_period {G : Type u} {X : Type v}
    (d : GroupActionModel G X) (t : TraceInterface d) (g : G)
    (excess : ¬ t.traceAtMostDimension g) (k : Nat) (positive : 0 < k) :
    power d g k ≠ d.one := by
  intro period
  exact excess (t.finiteOrderBound g k positive period)

theorem powers_distinct {G : Type u} {X : Type v}
    (d : GroupActionModel G X) (g : G)
    (noPeriod : ∀ k, 0 < k → power d g k ≠ d.one)
    (m n : Nat) (equal : power d g m = power d g n) : m = n := by
  have not_less : ∀ a b : Nat, power d g a = power d g b → ¬ a < b := by
    intro a b same less
    have split : b = a + (b - a) := by omega
    have product : d.mul (power d g a) d.one =
        d.mul (power d g a) (power d g (b - a)) := by
      rw [d.mul_one, ← power_add, ← split]
      exact same
    have period : power d g (b - a) = d.one :=
      (left_cancel d _ _ _ product).symm
    exact noPeriod (b - a) (by omega) period
  have first := not_less m n equal
  have second := not_less n m equal.symm
  omega

theorem powers_stabilize {G : Type u} {X : Type v}
    (d : GroupActionModel G X) (g : G) (C : X)
    (fixed : d.act g C = C) (n : Nat) : d.act (power d g n) C = C := by
  induction n with
  | zero => exact d.act_one C
  | succ n ih => simp only [power, d.act_mul, fixed, ih]

/-- `cosetOf g` denotes the coset gW of the full reflection subgroup W.
The difference criterion and free chamber action are standard geometric
inputs. The subgroup itself is a predicate; no subgroup law is used by
the separation proof. Correct interpretation in the actual isometry
group is a separate, explicit obligation. -/
structure ChamberCosets {G : Type u} {X : Type v}
    (d : GroupActionModel G X) (Q : Type w) (C : X) where
  reflectionGroup : G → Prop
  cosetOf : G → Q
  equalCosetDifference : ∀ a b, cosetOf a = cosetOf b →
    reflectionGroup (d.mul (d.inv b) a)
  freeChamberAction : ∀ a, reflectionGroup a → d.act a C = C → a = d.one

theorem stabilizers_have_distinct_cosets {G : Type u} {X : Type v}
    {Q : Type w} (d : GroupActionModel G X) (C : X)
    (q : ChamberCosets d Q C) (a b : G)
    (ha : d.act a C = C) (hb : d.act b C = C)
    (same : q.cosetOf a = q.cosetOf b) : a = b := by
  have inverseFixed : d.act (d.inv b) C = C := by
    have h : d.act (d.inv b) (d.act b C) = C := by
      rw [← d.act_mul, d.inv_mul, d.act_one]
    simpa only [hb] using h
  have fixed : d.act (d.mul (d.inv b) a) C = C := by
    rw [d.act_mul, ha, inverseFixed]
  have identity := q.freeChamberAction _ (q.equalCosetDifference a b same) fixed
  have equal := congrArg (fun z => d.mul b z) identity
  simpa only [← d.assoc, d.mul_inv, d.one_mul, d.mul_one] using equal

/-- An injection of Nat into the actual coset type witnesses infinite
index. Geometry and script results enter only through the named inputs. -/
theorem infinitely_many_cosets_from_trace {G : Type u} {X : Type v}
    {Q : Type w} (d : GroupActionModel G X) (t : TraceInterface d)
    (C : X) (q : ChamberCosets d Q C) (g : G)
    (fixed : d.act g C = C) (excess : ¬ t.traceAtMostDimension g) :
    ∃ f : Nat → Q, ∀ m n, f m = f n → m = n := by
  refine ⟨fun n => q.cosetOf (power d g n), ?_⟩
  intro m n same
  have equal := stabilizers_have_distinct_cosets d C q _ _
    (powers_stabilize d g C fixed m) (powers_stabilize d g C fixed n) same
  exact powers_distinct d g (no_positive_period d t g excess) m n equal

end CoxeterProofs.ChamberPeriod

#print axioms CoxeterProofs.ChamberPeriod.left_cancel
#print axioms CoxeterProofs.ChamberPeriod.power_add
#print axioms CoxeterProofs.ChamberPeriod.no_positive_period
#print axioms CoxeterProofs.ChamberPeriod.powers_distinct
#print axioms CoxeterProofs.ChamberPeriod.powers_stabilize
#print axioms CoxeterProofs.ChamberPeriod.stabilizers_have_distinct_cosets
#print axioms CoxeterProofs.ChamberPeriod.infinitely_many_cosets_from_trace
