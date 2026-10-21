# Cubic B16 short-cap logic kernel

## Status and exact source

Prepared source only, 22 October 2026. **UNCOMPILED: six requested axiom prints, zero actual outputs.** No compiler, Lean import execution, worker, enumeration, test or mathematical job was run during preparation. Source review is not elaboration or an actual axiom audit. Root owns any later resource-admitted single-module compilation and separate audit.

Source: [`formal/CoxeterProofs/CubicB16ShortCapKernel.lean`](../formal/CoxeterProofs/CubicB16ShortCapKernel.lean), 6,324 bytes/137 lines, SHA-256 `ce8b53cafa9cd18a1520545618bb854578d0e740a1c2c23a811fa16ac0901d74`.

The only import is `Std`; no existing local kernel was changed or imported. The intended toolchain is the repository's `formal/lean-toolchain`, `leanprover/lean4:v4.32.2`. No local `.olean` or verified declaration is claimed. The source follows the existing generic `Lean.Grind.OrderedRing` pattern; its geometric variables range over an arbitrary linearly ordered commutative ring R. Its separate Int variables are genuine integral polynomial/basis coefficients, not integer replacements for real geometry.

## Written proof interfaces

All three sources below are under `results/`:

| Source | SHA-256 | Used written interface |
| --- | --- | --- |
| `arithmetic17_b16_actual_i_short_opposite_isolation_and_degree_restriction_20261021.md` | `51910bc233e16870d2322e66c9f2c2766076d757749088e7743277b8ddbed248` | Actual compact arithmetic17 B16+isolated short I; isolated dotted short cap, actual K=Q(B), degree2/3, strong6..14 reference and actual I relation |
| `arithmetic17_b16_isolated_cubic_cap_unique_polynomial_candidate_20261021.md` | `2fe5f3c6129bb1cf264d00de2f062e068ed94e66dca290be4490d79a2d328b22` | Cubic conjugates, integral trace/coefficient bounds, F=X³−25X²+10X−1, physical B in(24,25), discriminant473 and full integral basis |
| `arithmetic17_b16_isolated_cubic_cap_actual_i_strip_exclusion_candidate_20261021.md` | `f89214ff279908fa2d1ebb1621102adab8a4c17785418c2519d00b2d3811bc45` | Complete actual I strip, J=6−1/B, height identity and same-vector paid-reference inequality |

Their written acceptance is [`results/d17_typed7_cubic_b16_and_prime239_transverse_acceptance_20261022.json`](../results/d17_typed7_cubic_b16_and_prime239_transverse_acceptance_20261022.json), SHA `5d009ab4bd1275a7b0839c5b07d3b457126f4f58f0e24abdd4cbc3938032a19d`. This module is a later prepared logic input, not an artifact already accepted or executed by that receipt.

## Six derived declarations

Names below have prefix `CoxeterProofs.CubicB16ShortCapKernel.`. Each has one requested `#print axioms`; there are no actual outputs yet.

1. **`cubic_coefficients_unique`** takes integral S,T with `17≤S≤28`, `T≥0`, `4T≤S+15` and `4S≤T²`. It derives **S25,T10**, covering all eleven integer values T0..10. The field norm-one argument and the two conjugate-polynomial tests producing these inequalities are external. No unique coefficient pair is supplied as a premise.

2. **`i_strip_coefficients_unique`** takes integral m,n, `−6≤n≤6`, and the two sign-dependent cleared strips. It derives **m0,n−1**, covering every allowed integer n. The actual roots obey `.184<u<.1845`, `.2205<v<.222`; the divided differences give `0<m−nv<.247` and `0<m−nu<.247`. For n≥0 these imply

       2205n < 10000m < 1845n+2470;

   for n<0 they imply

       1845n < 10000m < 2205n+2470.

   The root ordering, interval derivation, coefficient range |n|≤6, actual O_K basis and representation `J=l+mB+n/B` remain explicit written application inputs. Determining l6 from the physical and u strips, and hence identifying actual J, is also external to this kernel. The theorem derives the nontrivial complete integer m,n elimination; it does not formalize the whole algebraic-number-field strip.

3. **`i_height_lower_from_unit`** takes R-valued `B>24`, positive eta, `B*eta=1` and `4h=2−eta`. It derives **96h>47**. The strict inverse bound `24eta<1` is derived by multiplying the positive gap B−24 by eta, rather than supplied.

4. **`i_energy_lower_from_height`** takes `96h>47`, `0≤r<1` and the half-tangent identity `h(1−r)²=4r`. It derives **1000r>99**. It first obtains `47r²−478r+47<0`. If `1000r≤99`, the two nonnegative factors give

       (99−1000r)(478000−47(1000r+99))
         = 1000000(47r²−478r+47)−138647 ≥ 0,

   a contradiction. Thus no monotonicity-of-a-rational-function axiom is supplied. The interpretation of r as the hyperbolic half-tangent energy lower bound for the actual I wall is external.

5. **`paid_reference_norm_margin`** takes an R-valued reference base, last height f, short pressure S, orthogonal component energies eB/eI and original total, satisfying

       100000base ≥89026, 1000f ≤3311, 2S>7,
       5000eB ≥5000base+1300(S−f)−169,
       1000eI>99, total=eB+eI.

   It derives the actual cleared margin **100000total>100460**, hence total>1. The paid inequality is the written square completion with mu=.13, support-tail norm2 and an already funded P14 coefficient. The kernel derives its rational consequence; it does not assume that coefficient payment or vector orthogonality is geometrically valid.

6. **`cubic_short_cap_contradiction`** combines declarations3–5. Its primitive inputs include B/eta/height identities, the half-tangent relation, `r≤eI`, the reference floors/paid inequality, original orthogonal energy decomposition and original norm<1. It derives False. Neither an I-energy margin nor a total-energy contradiction is a premise. The coefficient results1/2 and their number-field application are not silently instantiated inside this final theorem.

## Explicit unformalized boundaries

This is conditional internal logic, not a formal proof of the B16 theorem or an arithmetic dimension bound. External inputs still include:

- The automatic admissible quadratic-form model for arithmetic reflection lattices, actual field equality, arithmetic integrality and all embedding bounds.
- Norm-one/trace reduction, derivation of the polynomial tests, root isolation and positivity, the discriminant-index identity and full O_K basis, and the l6/J identification.
- Compact acute/farthest geometry, the entire actual cap inventory, all-wall feasibility and endpoint recognition, actual-I status, and the relation of h/r to that wall.
- The SAME original vector, strong6..14 reference, individual anchor inequalities, positive residual coefficient after payment, tail orthogonality/norm2, and the derived geometric inequality relating eB to that reference. The supplied paid inequality must be established for that vector; unrelated energy estimates cannot instantiate it.

The kernel does not close the quadratic B16 branch, remove the B16+I forest row, produce a polytope, change any global ceiling, or formalize the six manuscripts. Successful future axiom prints would enumerate axioms, not these theorem parameters. Existing verified kernels and frozen papers remain unchanged.
