# Ramified-seven diagonal family: proof inputs and reproduction

The internally reviewed result is nonreflectivity of the **literal family**
\(T_n=\langle-t(2-t)\rangle\oplus I_n\), n≥3, over
\(O=\mathbf Z[t]\), \(t=2\cos(2\pi/7)\).
This is not an exclusion of all cubic-field forms, all determinant-seven
lattices, or a rational commensurability class. No new example or general
dimension ceiling follows.

This is a paper-only packet. The three preparatory proofs have full root
written reviews. The final uniform argument has full root and separate
adversarial written PASS reviews. No new exact checker, Lean module,
execution ledger, scientific
run, or external refereeing is claimed.

## Read order and exact source identities

| source | role | SHA256 |
|---|---|---|
| [Arithmetic calibration](../results/cubic7_diagonal_first_edge_and_literature_gate_20261020.md) | Full ring, PID/unit normalization, exact trace form and conjugate conventions | `f6251aacd1259f43527c6aa6ebd37c579edfbb0f4b5a1349c517545f9d6ac1f6` |
| [Norm/facet reduction](../results/cubic7_ramified_diagonal_norm_and_facet_reduction_20261020.md) | Apparent norms beta,2beta absent modulo beta; full B_n zero germ; FULL coordinate-facet identity | `34cef661bd0333b633d4cf0912b096a273635d3bff54db8fa358b7222385e87a` |
| [Complete first edge](../results/cubic7_ramified_diagonal_complete_first_edge_20261020.md) | Finite trace list; unique first norm-two wall r; full B_(n−1)+A1 endpoint | `6e8057bee97987b2bd231b07dcf113f769e66c86345d3f99743eca5e7b5f94ed` |
| [Complete second edge](../results/cubic7_ramified_diagonal_second_edge_and_finite_symmetry_20261020.md) | Seven finite profiles, exact finite interval, unique simple s; integral actual chamber involution G | `e4724ccce3557bf719586db0f7a3ae50bee5adb84da36d6a0408ded07c369c0d` |
| [Uniform obstruction](../results/cubic7_ramified_diagonal_uniform_noncompactness_20261020.md) | Explicit n=3 inventories, dyadic full-A3 exclusion, forced B3 cap H, no common hyperbolic fixed point, higher-dimensional descent | `90d097fb04ffe0336a95a84aa7bd54aa03e34520667b7596fed986d04860feda` |

The final proof does not assume a computed first Pell unit. Compactness
itself supplies a finite endpoint of the remaining A2 edge. Every possible
endpoint is addressed: an isolated root fails modulo beta; full A3 fails
at2 or enlarges to B3; one B3 orientation has an earlier integral mirror;
the other supplies H. The common-fixed-vector obstruction applies to all
possible parameters of H. G alone has order two and does not prove
nonreflectivity.

The lower-rank interface is explicit in Section1 of the uniform proof.
The earlier edge statements conservatively use n≥5; the uniform proof
checks the smaller finite-diagram inventories needed for n=3 before using
their scalar conclusions there. The full coordinate-facet identity, not a
selected embedded subdiagram, carries the dimension-three obstruction
upward.

## Verification and external boundary

The new scalar identities, finite trace lists, Gram maps, local determinant
argument and fixed-plane calculation are written arguments. There is no
numerical acceptance, optimization or unbounded height enumeration in the
proof. The imported field arithmetic uses the named Minkowski and Dirichlet
theorems. The finite spherical Coxeter classification is an explicit
structural input, including the rank-three A3/B3/H3 alternatives; its exact
contact use is stated in the source inventories.

The final cofinite conclusion additionally imports the standard arithmetic
lattice/anisotropic compactness theorem of Borel–Harish-Chandra. The
precise readable PRIMARY statement is their 1961 announcement,
*Arithmetic subgroups of algebraic groups*, Bulletin of the American
Mathematical Society67 (1961), pp579–583:
[authors' paper scan](https://www.mathi.uni-heidelberg.de/~wienhard/retreat10/references/borel_harishchandra1.pdf).
Theorem1(iv), p581, gives finite Haar measure; §4 Theorem3, p582, gives
compactness and explicitly identifies anisotropic orthogonal groups as an
example; Remark1, p582, extends the conclusions to number-field integers.
Apply these to the special orthogonal group, then its finite component
extensions. Compact factors at the other real places do not change the
conclusion at the chosen hyperbolic place. The full 1962 exposition is
[Annals of Mathematics75 (1962), 485–535](https://doi.org/10.2307/1970210).

Positive definiteness at another embedding makes the form K-anisotropic.
A cofinite reflection subgroup would have finite index in the cocompact
arithmetic orthogonal group by covolume multiplicativity, hence a compact
chamber. This external theorem is not replaced by an integer-height proxy,
an unproved ideal-ray assertion, or a claim of formal geometric coverage.

The exact Bugaenko1992 **unimodular** cubic classification is not an input.
Its precise original iff scope was not retrieved during the bounded
literature screen. Neither the author's overview phrase “n≤4” nor a
bounded nonclosed archive is used to infer this theorem.

## Existing reviews and historical evidence

The final uniform proof has complete
[root adversarial review](../results/review_cubic7_ramified_uniform_noncompactness_root_20261020.md)
and [separate full review](../results/review_cubic7_ramified_uniform_noncompactness_batch_guard_20261020.md).
Both bind the final90d097fb source and retain the external arithmetic
theorem and inherited field inputs explicitly. Further edge, root-prefix
and Pell searches of this literal family in n>=3 should stop.

The completed root reviews are
[norm/facet reduction](../results/review_cubic7_ramified_diagonal_reduction_root_20261020.md)
(`560476c4224cd2d964e3770a68bd21ffb27f7bf7c6d73415f8e674223a8edf1d`),
[first edge](../results/review_cubic7_ramified_complete_first_edge_root_20261020.md)
(`7c39346da70351b8a66f4fba4ae5196e57e56dd80498fcf02a2e8f0f0c5b8323`),
and [second edge](../results/review_cubic7_ramified_second_edge_root_20261020.md)
(`b7ce7e81a97ddcb914a60d2cf61d0f016d91cc51a12bacb59840ebd79395458f`).

[The old rank-five bounded archive](../results/cubic7_timelike_scaled_rank5_exact_prefixes.md)
provided useful root calibration. Its height1000 prefixes do not establish
reflectivity or nonreflectivity and were not rerun. They are provenance,
not premises in the uniform proof. The closed golden/silver/quadratic13
host packets and the earlier field-specific dimension refinements are not
replayed or imported here.

Reproduction currently means reading the five sources in the displayed
order and checking the stated finite identities. Byte checks may be made
with `shasum -a 256` on those paths; this is metadata validation, not proof
execution. No scientific replay command is claimed: a bounded checker
has not yet been prepared, reviewed and run. Research takes priority;
this packet does not add a TeX manuscript or Lean module.
