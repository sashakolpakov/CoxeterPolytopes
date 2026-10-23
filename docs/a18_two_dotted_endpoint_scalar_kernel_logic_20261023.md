# A18 two-dotted endpoints: conditional ordered-ring scalar kernel

23 October 2026. **PREPARED SOURCE ONLY: uncompiled, zero actual axiom
outputs.** Ten public declarations have requested prints. Author manual
algebra/source review is complete; different-author source/map review,
root preclear and guarded compilation plus a separate actual axiom audit
are still required. This is no execution or end-to-end geometric proof.

Source: `formal/CoxeterProofs/A18TwoDottedEndpointScalarKernel.lean`,
SHA256 `f46d389c5b694adb33494807e7875d632f4f03261488f22df6de5a37303c07c2`.
Only `Std` is imported. No existing module, runner, successful object,
paper, guide or proof receipt was modified or rerun. In particular the
successful `GeneralEighteenEndpointRefinement.lean` source f6630f0f is
unchanged and is NOT imported by this new module.

The NEW implementation is my work; root authored the settled geometric
theorem. I authored/reviewed some earlier actual-point, endpoint/reference
and paired-finite foundations. My present self-check is not an independent
review of this kernel. A different author must assess the frozen source.

## Written mathematical and acceptance bindings

| File | SHA256 | Restricted role |
|---|---|---|
| results/compact18_whole_a18_two_isolated_dotted_endpoints_paid_core_exclusion_candidate_20261023.md | 4ce056b22536aecfa9744e5b5ffccb1931c8f6db802eb2726ae3cea1bef5e4dd | Complete actual geometry, N19 inventory, one-reference payments and continuous theorem |
| results/review_compact18_whole_a18_two_isolated_dotted_endpoints_paid_core_full_scope_compact_reviewer_20261023.md | c611b369970905ceb721a98f9eabaf87fe226ac1735f370f339159b93e30b09e | Full different-author written review, not an imported Lean premise |
| results/review_compact18_a18_two_dotted_endpoints_and_d17_actual_i_reduction_root_acceptance_20261023.md | 14d26b6c3cebc60e8d156b28a045338a1f02111060a715798e2e0cae122e9c21 | Root acceptance of precisely the written restricted theorem |

All three literal bindings were freshly metadata-checked. The full written
arguments provide the interfaces; their hashes alone do not prove them.
The root acceptance also discusses a separate D17 reduction, which is not
part of this module.

## Primitive scalar interface and internally derived inequalities

All coordinate quantities range over ONE generic linearly ordered
commutative ring `R` with the existing `Lean.Grind.OrderedRing` interface.
They are not Nat/Int proxies for real coordinates. The actual ordered-ring
instantiation and existence/meaning of square roots remain external.
There is no `Real` square-root function or geometric instance in Std.

For an endpoint parameter `u`, primitive inputs are

    0 <= u <= 1, a>=0, b>=0,
    19*a^2=36-19*u, b^2=u*(1-u).

In the written application, `u=g^-2` for an actual isolated dotted own
contact `g>1`, so its domain is the smaller open interval `(0,1)`.
The kernel proves the entire CLOSED interval and retains a completely
independent second endpoint parameter `v`.

| Public theorem (namespace `CoxeterProofs.A18TwoDottedEndpointScalarKernel`) | Derived content |
|---|---|
| `first_radical_chords` | `10000*a` is strictly above `13760-3902u`, `14029-4440u`, `14419-4960u` on the three respective intervals |
| `second_radical_tangent` | `4b <= 1+4u-4u^2`, derived from `(2b-1)^2>=0` and the primitive square equation |
| `three_quadratic_seals` | The three cleared quadratics are positive on the needed intervals, with the middle one positive for ALL `u` |
| `continuous_core_gap` | `10000*(a-b)>7408-2500u+1250u^2` over the ENTIRE closed interval |
| `endpoint_gains` | The exact definitions of G5/G6 imply `G5>3/50`, `G6>4699/60000`, and `60000*(G6-G5)=1099` |
| `simultaneous_paid_completion` | One scalar pressure, one squared norm and one nonnegative completed square imply the joint budget `base+gainL+gainR<=total` |
| `zero_weak_end_budget` | The zero-weak reference and two G5 gains imply `10000000*total>10440352`, hence `total>1` |
| `one_weak_end_budget` | The one-weak reference and G5/G6 imply `60000*total>61921`, hence `total>1`, in either orientation |
| `two_weak_end_budget` | The two-weak reference and two G6 gains imply `150000*total>153116`, hence `total>1` |
| `paired_endpoint_budget_contradiction` | Primitive roots/gain identities at BOTH independent endpoints, `total<1` and the legal reference-pattern budget give `False`; no final scalar gain or norm margin is assumed |

The chord proof does NOT assume analytic concavity or a chord inequality.
For `L1=13760-3902u`, `L2=14029-4440u`, `L3=14419-4960u`, its exact
derived identities are

    38*(10^8*a^2-L1^2)
      =289286476*u*(1-2u)+4365478*(1-2u)+805722,
    19*(10^8*a^2-L2^2)
      =46819800*(2u-1)*(3-4u)+306280*(3-4u)+96581,
    19*(10^8*a^2-L3^2)
      =116857600*(4u-3)*(1-u)+310080*(1-u)+19061.

Each right side is positive on its interval. A factorization of the
difference of squares and `a>=0` then supplies the strict chord.
The three final quadratic identities are

    7704-22804u+17500u^2
      =677+(1-2u)*(7027-8750u),
    721175-2089500u+1531250u^2
      =2*(875u-597)^2+8357,
    36088-99680u+70000u^2
      =703+5*(4u-3)*(3500u-2359).

Interval factors have proved nonnegative signs; no square-sign premise
is supplied by the application. The interval partition is handled by
ordered cases, not finite sampling or a label cutoff. Only kernel-checkable
`grind`/ring-order reasoning is used, with no `native_decide`, axiom,
`sorry`, Mathlib dependency or external mathematical process.

## Actual gain definitions and the ONE-vector boundary

Put `K=a-b`, `H5=263/125`, `H6=451/200` and
`S_m=2-2/m-2u+u^2`. The written gain is

    G_m=(K-H_m/m)/2-S_m/16.

Its exact denominator-cleared definitions in the kernel are

    10000*G5=5000*(a-b)-3104+1250u-625u^2,
    60000*G6=30000*(a-b)-17525+7500u-3750u^2.

Thus the implementation does not supply G5/G6 lower bounds as inputs
to its complete scalar theorem. The standalone numerical budget lemmas
accept these previously DERIVED gain bounds for modular reuse; the final
theorem explicitly calls `endpoint_gains` twice from primitive roots.
It covers zero, one or two weak end sites, retaining both one-weak
orientations. Reference floors are respectively
`.9240352`, `.8937`, `.86414` (the latter two may safely be weakened
from the strict written floors to nonstrict inputs).

The standalone `simultaneous_paid_completion` derives the budget from

    4*base+K_L+K_R-h_L-h_R <= 4*pressure,
    16*normW=16*base+S_L+S_R,
    16*gain_i=8*(K_i-h_i)-S_i,
    0<=total-2*pressure+normW.

Here `pressure` means the ONE original `x.W`, `normW` means `||W||^2`
and `total` means `||x||^2`. These meanings are external; the algebraic
deduction is internal. It does not add separately optimized energies.
The final theorem's `legal_pattern_budget` is exactly this joint budget
and the appropriate reference floor, never the desired `total>1`.

The written proof must still supply ALL of the following:

* the actual compact whole-A18 farthest vertex, SAME vector x and its
  squared-norm bound `total<1`, actual isolated dotted opposites and
  reciprocal parameter/square-root interpretation;
* the all-original-wall feasible cores and the resulting actual endpoint
  pressure, including every unknown facet;
* the COMPLETE fresh N19 weak inventory and near4 removal, nonadjacent
  helper holes and exactly the allowable five-/six-coordinate end blocks;
* the actual reference w, its stated norm floor after EVERY internal
  omission, its concavity and available individually paid boundary prices
  `1/(4m)` at each end;
* the actual residual vectors `Y_m`, their squared norms `S_m`, their
  orthogonality to the SAME w and to each other, and the legal pressure
  for ONE `W=w+(Y_L+Y_R)/4`;
* the interpretation of the one nonnegative completed square as
  `||x-W||^2>=0` and the correct end-pattern budget selection.

No inventory, all-wall incidence, endpoint isolation, geometric square-root
existence, individual coefficient capacity or orthogonality is created by
these generic scalar declarations. The theorem remains conditional at
that explicit boundary. Printed Lean axioms do not list ordinary theorem
parameters; this map records them separately.

## Requested verification and exact scope

The source footer requests EXACTLY the TEN public names in the table,
fully qualified by the stated namespace. The private difference-of-squares
helper is an implementation lemma, not another requested public output.
Currently there are ZERO actual outputs, no successful object binding and
no runtime/axiom claim. Root alone may execute an admitted single-module
compile followed by a SEPARATE audit; no successful prior module requires
rerunning. Failed or deferred attempts would receive no formal credit.

Even after those checks, the result formalizes scalar logic supporting
ONLY the accepted BOTH-isolated-DOTTED endpoint exclusion. It does not
classify nonisolated endpoints, settle the mixed finite/dotted branch,
remove whole A18, lower compact18 to17 or provide a polytope example.
The current four-group written frontier, older five-group formal snapshot,
six frozen papers and all global ceilings are unchanged by preparation.
