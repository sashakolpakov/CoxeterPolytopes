# Cubic B16 short-cap kernel: independent source review

22 October2026. **PASS of the mathematical source logic and documented
assumption boundary; UNCOMPILED.** This is not an elaboration result,
kernel-verification receipt or actual axiom audit. No worker/import/test,
Lean or TeX process was run and no bound source was changed.

Reviewed complete files:

| File | SHA-256 |
|---|---|
| formal/CoxeterProofs/CubicB16ShortCapKernel.lean | ce8b53cafa9cd18a1520545618bb854578d0e740a1c2c23a811fa16ac0901d74 |
| docs/cubic_b16_short_cap_kernel_logic_20261022.md | 4f680190280224617a239cc9e3ab5f6e5c1f91c0121c9e55a48e03ad6a490181 |

The module is by the foundations agent, with root-authored written
mathematics. I did not author this source or its six declarations. I
contributed to inherited actual-foot/reference foundations; this review
does not present them as independently formalized geometric inputs.

## Complete six-declaration check

1. `cubic_coefficients_unique` is correct. The supplied linear bounds
   force0<=T<=10. The quadratic bound forces T>=9. If T9, the quarter
   bound gives S>=21 while4S<=81 gives S<=20. Thus T10 and the two
   remaining inequalities give S25. The proof covers every integer case;
   uniqueness is derived rather than assumed. Its nonlinear expression
   is ground arithmetic after substituting each T case, suitable in
   principle for `omega`; actual tactic elaboration is still untested.

2. `i_strip_coefficients_unique` is correct for the stated complete
   range -6<=n<=6. For nonnegative n the cleared intervals contain no
   integer m. For n=-1 they contain exactly m0; n=-2..-6 contain no
   integer m. The sign-dependent orientations match the actual u/v
   strips in the written source. Bounds on n, the integral basis and
   the later l6 identification are supplied externally, not proved or
   silently instantiated here.

3. `i_height_lower_from_unit` derives its inverse bound: positivity of
   `(B-24)eta=1-24eta` gives24eta<1. From4h=2-eta,
   96h=48-24eta>47. No inverse monotonicity or desired h lower bound
   is a hypothesis.

4. `i_energy_lower_from_height` correctly obtains
   `47r²-478r+47<0` from the positive product and the primitive
   half-tangent relation. Under1000r<=99, both displayed factors are
   nonnegative. Their exact product is

       1000000(47r²-478r+47)-138647,

   which is strictly negative. The constant is correct:
   `47*99²-478*99000+47000000=138647`. Positivity of the
   square uses r<1; no irrational root or monotonicity axiom is hidden.

5. `paid_reference_norm_margin` has the correct strict margin. Its
   primitive inequalities imply

       eB > .89026 + .26*(3.5-3.311) - .0338 = .90560,
       eI > .099,
       total > 1.00460.

   Equality is permitted in the reference floor, endpoint bound and
   paid inequality; the actual short pressure and I-energy bounds
   supply the needed strictness. The source derives both cleared margin
   and total>1 rather than taking either as an assumption.

6. `cubic_short_cap_contradiction` correctly composes3--5, using
   r<=eI and total=eB+eI for ONE supplied original vector. The final
   total<1 contradicts the derived strict bound. The two integer
   coefficient theorems are intentionally separate: their number-field
   applications are not implicitly constructed inside this theorem.

All ordered-ring identities and inequalities were rederived. The
generic R variables describe conditional scalar algebra, and the Int
variables describe integral coefficients. No real geometric coordinate
is replaced by an integer. Only Std is imported; the source has no
axiom declarations, `sorry` or `admit`. It requests exactly six prints.
No claim about their future actual axiom output is made by this review.

## Application-map and status check

The three written source pins in the map streamed to their stated hashes:
51910bc233e16870d2322e66c9f2c2766076d757749088e7743277b8ddbed248,
2fe5f3c6129bb1cf264d00de2f062e068ed94e66dca290be4490d79a2d328b22,
f89214ff279908fa2d1ebb1621102adab8a4c17785418c2519d00b2d3811bc45.
The cited written-acceptance receipt matches
5d009ab4bd1275a7b0839c5b07d3b457126f4f58f0e24abdd4cbc3938032a19d.
I read the cubic coefficient deductions and the full actual-I strip/paid
reference source needed to compare these interfaces. The pinned local
toolchain is leanprover/lean4:v4.32.2.

The map correctly leaves outside Lean: actual field equality and all
embeddings; integrality, norm-one and conjugate tests; root isolation,
discriminant/full integral basis and l6/J identification; complete cap
inventory; actual farthest geometry and feasible points; the half-tangent
interpretation; and realization of the same-vector reference and its
funded P14 coefficient, tail norm and orthogonality. In particular the
paid inequality is an explicit geometric application input, not a theorem
derived by this kernel. This limits the formal scope honestly.

The current status is **six requested prints, zero actual outputs,
uncompiled source**. It supplies no quadratic B16 exclusion, forest-row
removal, global ceiling, formal manuscript theorem or implementation
verification. No mathematical source correction is requested. Any future
elaboration failure must remain separate from this source-only PASS.
