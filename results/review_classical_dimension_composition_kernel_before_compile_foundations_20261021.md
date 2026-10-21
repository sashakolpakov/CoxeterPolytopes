# Classical dimension composition kernel: full precompile source review

21 October 2026. **Source/mathematics/assumption-map PASS; no correction
requested. UNCOMPILED, zero actual axiom outputs.** This review read the
entire177-line module and116-line map, with source-only inspection of the
used Nat APIs in the pinned Lean4.32.2 installation. No Lean import,
compiler, lake, test, scientific program or mathematical worker ran.
The frozen source/map and all papers, guides and checkpoints are unchanged.

| Artifact | SHA256 |
|---|---|
| formal/CoxeterProofs/ClassicalDimensionComposition.lean | c035f62de642c123ea168a61871bc6d3156113343b3f5301a877844f0062e102 |
| docs/classical_dimension_composition_logic_20261021.md | 1778aabbe0d8b7d3a7ecac83f5e7841c9c7ba7260fd5406122219798ee6b6659 |
| results/compact_dimension19_and_arithmetic_dimension17_global_composition_20261021.md | 0ac802108e0ef8a9ca821c3091fca82f7bf9331a0037639d93888f2fabb37a75 |
| results/arithmetic18_whole_d18_dotted_quartic_exclusion_by_two_cycle_20261021.md | 6796f2319b0decfeae75eada3a2c979d884ef408331891e8a216e1a1e370446a |

All four bindings freshly match. **Disclosure:** I did not author this
Lean module or its map. I authored the written quadratic D18 input and
contributed earlier actual-foot/prefix and D18 arithmetic foundations,
and separately reviewed the written global composition. This is an
outside-author review of the new formal source/map; it is not an
independent revalidation of my own inherited geometric input. The map
correctly leaves those written theorems as external application premises.

## All four statements and proofs

1. `quartic_norm_contradiction`, lines17–43, has the correct INTEGER norm
   statement. For positive a,b,n, divisibility witnesses give
   a*n=256*i and b*n=16*j. Their positive products force i,j>=1,
   hence a*n>=256 and b*n>=16. The explicit monotonic multiplication
   and reassociation prove4096<=(a*b)*(n*n); n<=6 gives n*n<=36
   and therefore4096<=36*(a*b). The last `omega` is linear in the
   common monomial a*b: 9ab<616 would give36ab<2464, a contradiction.
   Thus the substantive multiplicative lower bound is derived, not
   assumed. The term proofs use the actual argument order of
   `Nat.mul_le_mul`, `Nat.mul_le_mul_left`, `Nat.mul_assoc` and
   `Nat.mul_left_comm`; the literal4096=256*16 uses ordinary `decide`,
   not native evaluation. All positivity and cap hypotheses needed for
   the divisions-as-witnesses argument are present. Nat is appropriate
   here because these quantities are positive rational integer norms;
   it is not substituted for a real coordinate or height.

2. `degree_two_three_four_cases`, lines46–49, is exactly the complete
   integer split2<=d<=4. Zero and one are excluded by the explicit
   lower premise, not silently dropped by pattern matching. This is
   what the written degree reduction plus the external rational
   anisotropy/Meyer input supplies.

3. `compact_dimension_at_most_nineteen`, lines54–96, first separates
   dimension<=19. In the other branch the supplied PRIOR global20
   bound forces dimension20. It obtains ONE p from the complete
   actual-farthest forest witness, dispatches whole A20 and B20, and
   dispatches all three D20 head alternatives at that SAME p.
   Both the farthest predicate and whole-D20 predicate are passed
   into every head exclusion. No dimension<=19, dimension!=20,
   whole-D20 contradiction, face inheritance or downward existence
   premise is supplied. The predicates need not be mutually exclusive:
   exhaustive alternatives with individually applicable exclusions
   suffice. Empty vertex types do not introduce an implicit existence
   assumption, because the dimension20 forest premise supplies p.

4. `arithmetic_dimension_at_most_seventeen`, lines104–170, similarly
   uses the PRIOR global18 bound to force18 in the nontrivial branch.
   Its four forest alternatives preserve the actual typed-I predicates.
   In whole D18 it first consumes the exact residual finite/dotted
   head split. For the dotted branch it obtains ONE k with the actual
   ground-field relation and degree2..4. Each degree consumer receives
   the SAME p,k and all geometric/field hypotheses. The quartic branch
   obtains the primitive positive norm data and CALLS the first theorem;
   it does not assume quartic impossibility. The conjunction witness
   and degree disjunction have the right nesting and cover all cases.
   An empty abstract field type likewise causes no hidden problem:
   the actual dotted-field premise supplies the required k.

The source declares exactly these four theorems, with four matching
fully qualified `#print axioms` requests at174–177. Its only import is
`Std`; there are no local project imports, added axioms, `sorry`, native
evaluation or unsafe proof shortcuts. This is manual source inspection,
not an elaboration or axiom-output result.

## Faithfulness of the map and remaining external work

Map lines57–63 faithfully identify written G0/F/GA/GB/DH/D5/DF/DD with
the compact premises. In particular DH supplies the complete actual
isolated native5 / finite>=7 / dotted split, and D5 denotes the corrected
adaptive proof c68, not the rejected individual15/16/17-floor argument.
The prior global bound handles all larger dimensions explicitly.

Map lines65–82 retain the exact arithmetic frontier, actual isolated
short-I type, full-field branch theorems and native5's all-compatible-
golden-lattice removal. `RemainingFiniteD18Head` is explicitly the
residual7/8/10/12/18 list, rather than a license to omit native5, own0
or parallel. The degree premise includes actual field existence,
positive field degree, the upper reduction, and the external rational
Meyer exclusion; the formal code does not purport to derive those.
No finite admissible-field list or D19-to-D18 field transport enters.

Map lines84–89 correctly translate the quartic statement: a=N(U)>0,
b=abs(N(2-U))>0, n=N(V)>0, n<=6. Actual UV/4 and (2-U)V/2 integral
give256|an and16|bn at degree4, and the all-conjugate head bound gives
9ab<616. Those are the precise assumptions of the derived Nat
contradiction, with the norm interpretation/divisibility bridge left
external. The abstraction `Field : Type` plus its degree function and
actual-field relation is honest: no field operations or number theory
are formalized by that label type.

The map explicitly retains actual farthest-vertex existence/simplicity,
finite spherical/affine classification, positive finite projections and
future negative cores, all-wall signs, actual reciprocal endpoints,
native inventories, same-vector coefficient payments, first-type raw
normalization/cyclic integrality, all conjugate restrictions, the full
quadratic/cubic written proofs, and inherited script propositions.
None is instantiated in Lean. The one vertex witness passed unchanged
prevents a hidden switch of farthest vector at the composition level;
the actual same-vector energy arguments remain written hypotheses.

The map's opening UNCOMPILED status is accurate: no current module
object was present at the standard local output path, no compile/audit
receipt is claimed, and the four printer commands have ZERO actual
output credit. Root alone may later run an admitted bounded compile
and separate audit. Even success would establish this conditional
internal arithmetic/composition kernel, not either end-to-end geometric
dimension theorem, a fresh mathematical census or changed paper text.
