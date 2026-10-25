import Std

/-!
PREPARED SOURCE ONLY: no compilation or separate actual axiom audit credited.

Discrete internal logic of the whole-D16 quadratic dyadic argument.
Number fields, actual valuations, norm bounds, raw K-time squares, finite
integer-strip tables and geometry remain external ordinary parameters.
The variables below are nonnegative integer valuations and parity witnesses;
none is an assumed instance of an actual Coxeter polytope.
-/

namespace CoxeterProofs.QuadraticDyadicTripleKernel

/-- Inert unit head: the quarter-cycle, norm bound and time-square parity conflict. -/
theorem inert_unit_head_contradiction (k parity : Nat)
    (quarter : 2 ≤ k) (norm_bound : 2 * k ≤ 4)
    (time_parity : 1 + k = 2 * parity) : False := by
  omega

/-- Split unit head: both primes would consume at least six norm valuations. -/
theorem split_unit_head_contradiction (k l first_half second_half : Nat)
    (first_quarter : 2 ≤ k) (second_quarter : 2 ≤ l)
    (norm_bound : k + l ≤ 4)
    (first_parity : 1 + k = 2 * first_half)
    (second_parity : 1 + l = 2 * second_half) : False := by
  omega

/-- At a ramified unit head, each bounded fork valuation is exactly four. -/
theorem ramified_unit_fork_four (k parity : Nat)
    (quarter : 4 ≤ k) (norm_bound : k ≤ 4)
    (time_parity : 2 + k = 2 * parity) : k = 4 := by
  omega

/-- Inert head valuation one forces fork valuation one. -/
theorem inert_first_head_fork_one (k parity : Nat)
    (quarter : 2 ≤ 1 + k) (norm_bound : 2 * k ≤ 4)
    (time_parity : 3 + k = 2 * parity) : k = 1 := by
  omega

/-- Inert head valuation two needs a fork divisible by four. -/
theorem inert_second_head_pair_forces_two (k l first_half second_half : Nat)
    (first_norm : 2 * k ≤ 4) (second_norm : 2 * l ≤ 4)
    (pair_quarter : 2 ≤ k + l)
    (first_parity : 4 + k = 2 * first_half)
    (second_parity : 4 + l = 2 * second_half) : k = 2 ∨ l = 2 := by
  by_cases first_zero : k = 0
  · right
    omega
  · left
    omega

/-- Ramified head valuation one forces each fork valuation four. -/
theorem ramified_first_head_fork_four (k parity : Nat)
    (quarter : 4 ≤ 1 + k) (norm_bound : k ≤ 4)
    (time_parity : 4 + k = 2 * parity) : k = 4 := by
  omega

/-- Ramified head valuation two forces each fork valuation three. -/
theorem ramified_second_head_fork_three (k parity : Nat)
    (quarter : 4 ≤ 2 + k) (norm_bound : k ≤ 4)
    (time_parity : 7 + k = 2 * parity) : k = 3 := by
  omega

/-- Ramified head valuation three requires one high-valuation fork. -/
theorem ramified_third_head_pair_forces_three (k l first_half second_half : Nat)
    (first_quarter : 4 ≤ 3 + k) (second_quarter : 4 ≤ 3 + l)
    (first_norm : k ≤ 4) (second_norm : l ≤ 4)
    (pair_quarter : 4 ≤ k + l)
    (first_parity : 7 + k = 2 * first_half)
    (second_parity : 7 + l = 2 * second_half) : k = 3 ∨ l = 3 := by
  by_cases first_high : k = 3
  · exact Or.inl first_high
  · right
    omega

/-- The split unit-prime pressure and the fork cycle force both norms sixteen. -/
theorem split_odd_trace_pair_valuations (kp kq lp lq : Nat)
    (first_unit_prime : 3 ≤ kq) (second_unit_prime : 3 ≤ lq)
    (first_norm : kp + kq ≤ 4) (second_norm : lp + lq ≤ 4)
    (pair_quarter : 2 ≤ kp + lp) :
    kp = 1 ∧ kq = 3 ∧ lp = 1 ∧ lq = 3 := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> omega

/-- An even fork valuation at the nonunit head prime rules out this split pair. -/
theorem split_odd_trace_even_valuation_contradiction
    (kp kq lp lq parity : Nat)
    (first_unit_prime : 3 ≤ kq) (second_unit_prime : 3 ≤ lq)
    (first_norm : kp + kq ≤ 4) (second_norm : lp + lq ≤ 4)
    (pair_quarter : 2 ≤ kp + lp) (even_valuation : kp = 2 * parity) :
    False := by
  have first_one := (split_odd_trace_pair_valuations kp kq lp lq
    first_unit_prime second_unit_prime first_norm second_norm pair_quarter).1
  omega

/-- In the even-trace split case one fork is divisible by four at BOTH primes. -/
theorem split_even_trace_pair_forces_global_four
    (kp kq lp lq ap aq bp bq : Nat)
    (first_head_quarter : 1 ≤ kp) (second_head_quarter : 1 ≤ lp)
    (first_norm : kp + kq ≤ 4) (second_norm : lp + lq ≤ 4)
    (pair_quarter : 2 ≤ kq + lq)
    (first_p_even : kp = 2 * ap) (first_q_even : kq = 2 * aq)
    (second_p_even : lp = 2 * bp) (second_q_even : lq = 2 * bq) :
    (kp = 2 ∧ kq = 2) ∨ (lp = 2 ∧ lq = 2) := by
  by_cases first_q_zero : kq = 0
  · right
    constructor <;> omega
  · left
    constructor <;> omega

/-- A positive fork norm below thirty-one, divisible by sixteen, is sixteen. -/
theorem bounded_norm_multiple_sixteen (n multiplier : Nat)
    (positive : 1 ≤ n) (upper : n ≤ 30)
    (divisibility : n = 16 * multiplier) : n = 16 := by
  omega

/-- Exact norm valuation three permits only norms eight and twenty-four. -/
theorem bounded_norm_valuation_three (n odd_part parity : Nat)
    (upper : n ≤ 30) (factorization : n = 8 * odd_part)
    (odd_factor : odd_part = 2 * parity + 1) : n = 8 ∨ n = 24 := by
  by_cases first_factor : odd_part = 1
  · left
    omega
  · right
    omega

/-- The final deduction uses the actual global sixteen theorem, not face heredity. -/
theorem global_dimension_fifteen (n : Nat)
    (prior_global_ceiling : n ≤ 16) (dimension_sixteen_excluded : n ≠ 16) :
    n ≤ 15 := by
  omega

end CoxeterProofs.QuadraticDyadicTripleKernel

#print axioms CoxeterProofs.QuadraticDyadicTripleKernel.inert_unit_head_contradiction
#print axioms CoxeterProofs.QuadraticDyadicTripleKernel.split_unit_head_contradiction
#print axioms CoxeterProofs.QuadraticDyadicTripleKernel.ramified_unit_fork_four
#print axioms CoxeterProofs.QuadraticDyadicTripleKernel.inert_first_head_fork_one
#print axioms CoxeterProofs.QuadraticDyadicTripleKernel.inert_second_head_pair_forces_two
#print axioms CoxeterProofs.QuadraticDyadicTripleKernel.ramified_first_head_fork_four
#print axioms CoxeterProofs.QuadraticDyadicTripleKernel.ramified_second_head_fork_three
#print axioms CoxeterProofs.QuadraticDyadicTripleKernel.ramified_third_head_pair_forces_three
#print axioms CoxeterProofs.QuadraticDyadicTripleKernel.split_odd_trace_pair_valuations
#print axioms CoxeterProofs.QuadraticDyadicTripleKernel.split_odd_trace_even_valuation_contradiction
#print axioms CoxeterProofs.QuadraticDyadicTripleKernel.split_even_trace_pair_forces_global_four
#print axioms CoxeterProofs.QuadraticDyadicTripleKernel.bounded_norm_multiple_sixteen
#print axioms CoxeterProofs.QuadraticDyadicTripleKernel.bounded_norm_valuation_three
#print axioms CoxeterProofs.QuadraticDyadicTripleKernel.global_dimension_fifteen
