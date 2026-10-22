# D15 typed-seven centre kernel: full pre-execution source review

22 October 2026. **PASS for the source mathematics and assumption map; no static correctness blocker found.** Compilation and actual axiom output remain pending. No Lean, compiler, import, test or mathematical worker was run by this reviewer.

Reviewed in full:

- `formal/CoxeterProofs/D15TypedSevenCentreKernel.lean`, SHA-256 `37df90e1536c5f9134e007f84e0b29d7b73032272c434408cda3574030b225b3`, 6,804 bytes/155 lines.
- `docs/d15_typed_seven_centre_kernel_logic_20261022.md`, SHA-256 `e78df38e62cc3b39e98372fba26d922f456e35d7db947b664f76bf016935803c`, 3,825 bytes/67 lines.
- Its written application source, SHA-256 `7cb0a001efbd678de6bf26b5fe4d3712690f50281e1e7a94d65c86b80ecd0418`, is the exact corrected D15 theorem reviewed in `results/review_arithmetic17_d15_typed_seven_centre_head_full_scope_foundations_20261022.md`, SHA-256 `f29341b05666f47bc3d1478d75184ea2c17a73f42a89af7cbac1b80f5fd63d55`.

The reviewer did not author this module or map. The reviewer authored an inherited dimension17 frontier and contributed/reviewed earlier actual-point and prefix arguments. This is an independent review of the new six conditional algebraic proofs and their stated application boundary, not a claim of independently authoring their shared geometric foundations.

## Six statements and proofs

1. **`sqrt_six_seals`, lines27–49.** The two difference-of-squares identities are exact: 17²−49*6=−5 and 2000²*6−4899²=−201. In the negated target branches, both factors have the required nonnegative signs. Their products cannot be negative in the specified ordered ring. Nonnegative s and its square equation therefore imply both strict seals; no square-root estimate is supplied as a hypothesis.

2. **`first_radical_chord`, lines52–74.** Expanding either side gives −3729+7024x−3136x² after a²=15−2x. On 1<=x<=5/4 the product (x−1)(5−4x) is nonnegative and 191−32x>=151, so the squared gap is strictly positive, including the closed endpoints. If the linear target failed, a>=0 would make both difference-of-squares factors nonnegative, contradicting that gap. Thus the module derives the strict chord 200a>777−56x. The interval is wider than the application's 1<x<100/81 and is valid throughout; this does not assert an unproved geometric parameter range.

3. **`second_radical_tangent`, lines77–90.** The expanded square (sb−1)²=2−36z²−2sb follows exactly from s²=6 and 6b²=1−36z². Its nonnegativity gives sb+18z²<=1. Multiplication by nonnegative s then gives 6b+18sz²<=s. A sign assumption on b is unnecessary: the asserted inequality is valid for either algebraic root. In the geometric application b is the nonnegative radical, but omitting that redundant restriction does not weaken correctness.

4. **`centre_radical_gap`, lines93–119.** I rederived the first polynomial inequality by multiplying the chord by420, the tangent by14000, and the norm equation by3000, using 6z=6x−7. This produces precisely

       223900−14000s−11520z+(252000s−36000)z²
         <6000*(14(a−b)−norm).

   The derived lower seal 7s>17 yields 252000s−612000>=0. Multiplying this by z²>=0 legally lowers the quadratic coefficient to576000. The upper seal 2000s<4899 gives 223900−14000s>189607. The exact final completion is

       5*(189607−11520z+576000z²)
         =288*(100z−1)²+947747.

   Hence 30000*(14(a−b)−norm)>947747>947700, which gives the requested 100*(14(a−b)−norm)>3159. Every strict sign is preserved. The primitive radical square equations, root sign for a, interval, substitution and norm identity suffice; neither a supplied chord/tangent nor the desired final scalar is a theorem parameter.

5. **`funded_centre_margin`, lines122–134.** Clearing the supplied legal paid-reference inequality gives the correct common factor49. Its weakest stated constants imply

       total > .84816+.02080125+.09
                   +(31.59−14*2.886+11)/49
             = .95896125+2.186/49
             >1.00357.

   The strict scalar and spectator premises supply strictness even though the base and head-gain bounds are conservatively weak. The conclusion 100000*total>100357 and then total>1 is true. The theorem does not claim to derive the geometric paid-reference inequality itself.

6. **`strong_centre_margin`, lines137–146.** The three input floors sum to exactly .88989+.02080125+.09=1.00069125, with strictness supplied by the actual defect bound. The exact integer target100069125 at scale100000000 is correct. This conditional final-reference implication assumes neither its conclusion nor a global dimension bound.

## Logic and execution boundary

The generic carrier is an abstract linearly ordered commutative ring under the existing `Lean.Grind.OrderedRing` pattern, not an integer surrogate for real geometry. The only import is Std; no project dependency or application instance is imported. The ring expansion, order multiplication, explicit nonnegative-square steps and linear combinations have the appropriate hypotheses. The proof syntax uses the same available ordered-ring APIs as existing kernels. There is no apparent unqualified-transitivity/API issue, axiom, `sorry`, `admit` or hidden conclusion premise. Actual elaboration and tactic performance are nevertheless matters for root's separately guarded compilation; this review is not execution evidence.

Exactly six public theorem declarations and six matching `#print axioms` commands appear. They request the six names in this review, with no duplicated printer or unsolicited output counted. At this source-review point the actual print count remains ZERO. No expected axiom list has been substituted for actual output.

The map accurately binds the source and corrected written proof, describes the wider chord interval, and leaves the following external: the actual ground field and complete finite/affine cap lists; original-wall feasibility and radical interpretation; actual parameter/substitution/norm identities; reference and defect floors; selection of the appropriate branch; simultaneous positive P5/P11 capacities and orthogonality on ONE original vector. Its endpoint qualification correctly refers to a reference ending at11. Neither funded theorem licenses adding independently optimized energies, and neither charges the spectator twice. The original total<1 assumption is needed externally to turn the derived total>1 into the geometric contradiction; the module intentionally proves the forward lower-bound implication rather than a polytope theorem.

The map's prepared/uncompiled status and zero-credit statement are truthful. This review supports one guarded module compile followed by a separate six-name actual axiom audit; it supplies neither of those outcomes. No existing successful kernel, frozen manuscript, guide, accepted theorem or global bound was changed.
