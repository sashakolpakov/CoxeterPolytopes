# Classical dimension composition: verified conditional logic

**Bounded compilation and separate actual axiom audit PASS: four declarations.**
The [ClassicalDimensionComposition.lean](../formal/CoxeterProofs/ClassicalDimensionComposition.lean)
is SHA256 `c035f62de642c123ea168a61871bc6d3156113343b3f5301a877844f0062e102`
(177 lines,8970 bytes). Its only import is `Std`; it adds no local
dependencies, axioms, `sorry`, native evaluation, field library or geometry
framework. All four `#print axioms` commands now have actual outputs, using
only `propext` and `Quot.sound`, with no `sorryAx` or warnings. The
[completion receipt](../results/lean_publication_completion_20261021/manifest.json)
binds unchanged source, object/compiler ledgers and actual stdout. Compile
and separate audit each took0.34s under Lean4.32.2,2048MiB/one thread and
unchanged resource gates. No successful dependency was rerun. The source's
preparation header remains historical to preserve its compiled hash.

The mathematical targets are the frozen written
[global composition](../results/compact_dimension19_and_arithmetic_dimension17_global_composition_20261021.md),
SHA256 `0ac802108e0ef8a9ca821c3091fca82f7bf9331a0037639d93888f2fabb37a75`,
especially Sections3–4, and the
[quartic norm proof](../results/arithmetic18_whole_d18_dotted_quartic_exclusion_by_two_cycle_20261021.md),
SHA256 `6796f2319b0decfeae75eada3a2c979d884ef408331891e8a216e1a1e370446a`,
Section3. These written sources and their separate reviews supply the
application premises. Linking them does not instantiate those premises
in Lean or turn source hashes into proofs.

## Derived statements

All names below have namespace `CoxeterProofs.ClassicalDimensionComposition`.

| Declaration / actual axiom print | What is derived |
|---|---|
| `quartic_norm_contradiction` | Positive integer a,b,n with n<=6,256 dividing an,16 dividing bn and9ab<616 are impossible |
| `degree_two_three_four_cases` | The complete Nat case split2<=d<=4 implies d=2 or3 or4 |
| `compact_dimension_at_most_nineteen` | Prior global20 plus the complete actual dimension20 forest and its three written branch exclusions imply dimension<=19 |
| `arithmetic_dimension_at_most_seventeen` | Prior global18 plus the complete actual dimension18 forest, typed-I exclusions, finite/dotted D18 split and actual field-degree branches imply dimension<=17 |

The first theorem opens the two divisibility witnesses. Positivity gives
256<=an and16<=bn; `Nat.mul_le_mul` and explicit multiplication
reassociation give4096<=ab*n*n. A second monotonicity argument gives
n*n<=36 and hence4096<=36ab. Only then does `omega` contradict9ab<616.
It does not assume a supplied final product lower bound or contradiction.
All multiplication is exact Nat multiplication; there is no subtraction,
saturation, division or proxy for a real height.

The arithmetic composition calls this norm theorem in its quartic case.
Its `writtenQuarticNormData` premise supplies positive integers and their
two divisibilities and bounds; it does **not** supply quartic impossibility.
The quadratic and cubic impossibilities remain named written theorems.
The separate degree-case theorem is called before dispatching those cases.

## Precise conditional interfaces

Each dimension theorem concerns ONE fixed candidate polytope. `dimension`
is its dimension. `Vertex` is any type containing its actual vertices;
`ActualFarthest` means farthest from the SAME chosen strict interior point.
The forest witness p is obtained once and passed unchanged to every
applicable head, typed-I and field premise. The code neither chooses a
different optimizer nor combines energies from distinct points.

For the compact theorem the source IDs in the written composition map as
follows: G0 is `writtenPriorGlobal20`; F is `writtenActualForestAt20`;
GA/GB are the whole-A20/B20 exclusions; DH is the isolated native5 / finite
m>=7 / dotted actual D20 head split; D5/DF/DD are its three exclusions.
The complete prior bound handles ALL higher dimensions. No downward
heredity of existence is assumed. There is no premise dimension<=19,
dimension!=20 or an unclassified whole-D20 impossibility.

For the arithmetic theorem A0 is `writtenPriorGlobal18`; F supplies exactly
typed A17+I, B17+I, D17+I or whole D18. The I predicates MUST mean the actual
isolated short A1 of that accepted frontier, rather than an arbitrary
ordinary A1 spectator. AI/BI/DI supply their full actual-field exclusions.
H18 gives actual head isolation and removes native5 by its named theorem
for EVERY compatible golden invariant lattice. Accordingly
`RemainingFiniteD18Head` denotes only its residual finite list7/8/10/12/18;
F18 excludes that list, including its separately composed12/18 input.
No native5 or own-zero/parallel case is silently assigned to F18.

`ActualGroundField p k` identifies the actual first-type field at this
vertex; `fieldDegree k` is its Nat degree. The combined written field
premise supplies2<=degree<=4: K18 supplies the upper bound and actual field
equality; positivity of degree and the named rational anisotropy/Meyer
input exclude degree1 in the19-variable rational form. Thus neither zero
nor the rational case is lost in the Nat split. Q2/Q3 give the complete
quadratic/cubic exclusions over that SAME actual field. Q4 gives the norm
data consumed by the formal quartic arithmetic.

For that norm bridge, a=N(U)>0, b=abs(N(2-U))>0 and n=N(V)>0 in actual
degree4. The all-conjugate strip yields n<=6 and9ab<616. Actual integral
UV/4 and (2-U)V/2 yield256 dividing an and16 dividing bn. None of the
number-field interpretation, norm signs, integral-root normalization or
passage from algebraic-integer divisibility to these Nat statements is
formalized by this module; they are exactly `writtenQuarticNormData`.

## Application work remaining outside Lean

The abstract predicates have no geometric instance here. All of the
following remain written, published or inherited finite-data inputs:

* existence/simplicity and the complete spherical forest at an actual
  farthest vertex; the prior global ceilings and every frontier elimination;
* all actual positive finite projections and future negative cores,
  all-wall feasibility, opposite-wall recognition, reciprocal endpoints,
  native-angle inventories and legally paid same-vector energy estimates;
* the first-type admissible quadratic space, actual K-valued raw roots,
  field equality, integral Cartan cycles/two-cycles, square classes and
  EVERY nonphysical conjugate restriction;
* the Q2/Q3 complete arithmetic tables and local/number-field arguments,
  the Q4 norm bridge above, and the external rational Meyer/anisotropy
  theorem and published invariant-lattice reflectivity inputs;
* any published classifications and exact script outputs inherited by the
  prior global/frontier theorems, with their original declared boundaries.

No new script is required by the formal proof. This source does not
replay certificates, prove field-degree reduction, classify fields or
construct a geometric Coxeter polytope. The successful compile and separate
audit verify a conditional internal arithmetic/composition kernel, not an
end-to-end Lean proof of the global bounds. The six delivered TeX/PDF
editions, formal root, lake configuration, paper maps and frozen checkpoints
are unchanged by this verification and status update.
