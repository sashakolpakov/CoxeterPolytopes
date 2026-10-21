# Independent mathematical assessment of the full-budget tail repair

Reviewed insertion:
`/private/tmp/finite57_full_budget_classical_unit_insert_20261020.tex`,
SHA256 `906cce8dacc25790fe582e0d1cbb8e3df745b95ec91b6152112bd0d3b288f51a`.
The canonical manuscript was read but not edited. No producer, verifier,
Lean, TeX, numerical, or scientific job was run.

**Verdict: the reported omission is real, and the proposed written repair
is mathematically sound.** The residual support census with endpoint caps
`3/2` cannot classify a failing factor of the original product, whose
endpoint budgets are `2`. Nor does a coefficient-`4/3` charge imply unit
charge. This was an omitted justification, not a counterexample to the
intended estimate. The insertion supplies the necessary stronger theorem.

I independently derived the following before reading the insertion, then
checked its full text and all 24 table entries.

1. Write `ell_ik=R(i,j)+R(i,k)-R(k,j)` and
   `c_ik=delta_k^2 R(i,k)`. Polarization gives the displayed lower bound for
   the convex-mean energy, and both affine heights are at least one, so
   `W>=E`. Nonnegative inverse entries and the zero-or-at-least-one contact
   gap give actual retained cost at least the sum of the contacted unit
   diagonal costs. Thus `ell_ik<=c_ik` on every nontrivial contact class
   proves unit charge; a single `c_ik>T` instead pays the target immediately.
   The insertion's complete B/C/D table satisfies precisely this dichotomy.
   I recalculated each entry from its prefix potential, including modified
   B rank three, modified C rank two, reversed rows, both D fork entries,
   and empty integer ranges. The omission distances exceed the targets.
   Separate modified B2 also has the stated point distance and diagonal
   threshold. This argument is coefficient one and does not use height
   balance or the residual cap.

2. For A cycles with original budgets `C_0,C_1<=2`, every retained cyclic
   distance at least three is impossible at the large orders. For distances
   one and two, direct resistance substitution gives
   `ell/c <= (m+epsilon)/(m-2)`. Consequently
   `W+C_0+C_1-T >= R(0,floor(m/2))-T-4(2+epsilon)/(m-2)`.
   I recomputed both positive numerators and their substitutions:
   `128+84u+4u^2` at even order `18+2u`, and
   `360+1708u+252u^2+8u^3` at odd order `19+2u`.
   This covers every remaining order without a sampled-rank assumption.

3. At orders 2--17 the insertion defines the **full** budget-two support
   universe, including all inverse cross terms, total cost below the
   ordinary target, and both omitted-node states. This necessarily contains
   every possible failing actual row. The Euclidean potential Gram, KKT
   expansion and height lower bound prove the claimed implication for
   arbitrary real contact magnitudes. In particular, the height bound must
   be evaluated at `21/20`; the historical `31/30` result alone would not
   suffice. Rotation covers avoid-node A with the stronger ordinary target;
   repeated avoid-A1 has target zero. The D automorphisms preserve its
   stronger ordinary estimate.

The exact small-A inequalities remain an explicitly stated finite computed
proposition. I checked that its domain and meaning agree with
`fresh_nikulin_height_charge_20261015_review.py:cycles` and that its wider
height window is recomputed in
`fresh_nikulin_defect_frontier_20261015.py:widened_window` and independently
in `fresh_nikulin_defect_frontier_20261015_review.py:cycle_replay`.
The latter reconstructs the complete 580 domains and the rational
closest-point certificates. This review does not claim a new replay of
those 580 records or replace their recorded execution evidence.

Finally I checked the interface to the separate exceptional census. An
actual failure has `C_0+C_1<T`, so setting its nonzero retained contacts to
one gives a complete necessary cheap-support list, without imposing the
residual cap. The previously defined exceptional finite check and same-row
E7 exclusions then leave exactly the stated ordered floors. Subtracting
those floors from `(2,2)` leaves `(10/9,4/3)` or `(1,3/2)` and its reverse.
Only here do the residual `3/2` support census and packing theorem apply.
The insertion therefore removes the circular coverage risk and restores
the intended classical/exceptional split.
