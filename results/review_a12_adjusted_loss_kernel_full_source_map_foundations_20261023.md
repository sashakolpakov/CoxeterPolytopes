# Full static review: A12 adjusted-loss kernel and assumptions map

23 October2026. **FULL mathematical/source/map PASS at the prepared-source
level, including the final outer-coefficient strengthening.** I read the
entire237-line Lean source and entire126-line map, the entire349-line
accepted mathematical source, its entire237-line outside review and
entire119-line root acceptance. No compiler, Lean import execution,
mathematical worker, enumeration, symbolic process, test or build ran.
Hashing, gzip roundtrip and literal-delta checks below are metadata checks.
Compilation and actual axiom-output verification remain PENDING.

## 1. Exact reviewed objects and written acceptance

| Object | Bytes | SHA256 |
|---|---:|---|
|formal/CoxeterProofs/A12AdjustedLossKernel.lean|10075|5ab0165dd843e9599f7a2c4ccb68701f2bcc9f5d701d3715ba5a7eecf9f4858d|
|docs/a12_adjusted_loss_kernel_logic_20261023.md|6114|37466f322b2e0ddf1d3568d3c2dd687d034195d8d90c9793cd67b9297609a3fb|
|results/arithmetic16_a12_both_typed_seven_rank_two_adjusted_loss_exclusion_candidate_20261023.md|17209|a35d9f561fad54a5cc21af8448c884ce9528d9b63de50cb2a2b27def4f72394a|
|results/review_arithmetic16_a12_both_typed_seven_adjusted_loss_rank_two_full_scope_compact_reviewer_20261023.md|12827|6a44888d083f87ebe602cd90ab9a7491fc381ff1e98645a9067209f5c38eb63b|
|results/review_arithmetic16_a12_both_typed_seven_root_acceptance_and_twelve_group_frontier_20261023.md|5934|209a01e628115b2e6b6368d1a2bb5d3c64082450609a8f4eda1997d54250184a|
|results/a12_adjusted_loss_outer_point08_strengthening_preservation_20261023.json|2950|163af537e741e901f5012dc8bfd99a3db08ac73d35b036b5a23f88442697689a|

I also freshly matched all seven literal mathematical input pins in a35d:
K980c5730, H535ac0d0, N801581c4, S68ad62a6, Of8edeec0,
X3732b826 and Fa8cc202c. Their full hashes are exactly those in the
accepted source and its outside review. This is an interface review of
the new scalar formalization, not a new transitive regeneration of every
geometric or arithmetic input.

I did not author this new A12 scalar module, map or a35d adjusted-loss
proof. I contributed to inherited actual-point, normalization and reference
work and reviewed related component arguments. Those shared foundations
are disclosed; they are not presented as independent new results here.

## 2. Preserved predecessor and exact strengthening

The preservation receipt binds the original source69317ae63fca640319829bce5f66cf6620439a3459c29c29d949f8b8d6e2ea43
(10077B) and map497104e52b0e69c57b0452e8c0401a2b08f4b2092509a27c1bd01346744e255e
(6115B). I independently checked both compressed hashes, decompressed
byte lengths/hashes, current hashes and the ENTIRE declared replacement
sets against the live successors. Both roundtrips and exact deltas PASS.

The archives under
`results/source_archive/a12_adjusted_loss_before_outer_point08_strengthening_20261023/`
are `formal__CoxeterProofs__A12AdjustedLossKernel.lean.gz`
(2805B,9067b8364c5ddf4914cd8f68178240e5a3d55d4ce6c9853959a903b77819253c)
and `docs__a12_adjusted_loss_kernel_logic_20261023.md.gz`
(2969B,5dc271fbb024a022b15f8f57d499d9143f6e02782b708097ea5990a98ea997e5).

Only the outer coefficient floor is weakened from beta>=.15 to
beta>=.08, with the corresponding polynomial264u+80u² replacing
124u+150u² and the corresponding multiplication factor80 replacing150.
The role disjunction is changed consistently in the theorem AND `Foot`.
The inner role, all reference data, list logic, final theorem and ten
requested print names are unchanged. The actual accepted outer reference
coefficient .15745 satisfies the new floor. No proposed later geometric
application of the .08 strengthening is assumed or accepted by this review.

## 3. Primitive circle and radical proofs

`circle_upper` is valid over the displayed ordered commutative ring.
Writing X=a(2u-1)+2cb and Y=c(2u-1)-2ab, its identity is exactly
X²+Y²=(a²+c²)((2u-1)²+4b²)=a²+c². The primitive equation
b²=u(1-u) supplies the last equality. Nonnegative Y² and the strict
gap a²+c²<r² give X²<r². If the target failed, X>=r>=0 and
(X-r)(X+r)>=0 would contradict that gap. No assumed continuous
maximum, interval sampling, sign of a/c or square-root function is used.

`radical_half` follows from (2u-1)²>=0 and the same primitive equation:
4b²<=1. A hypothetical2b>1 makes (2b-1)(2b+1)>0, contradicting
this inequality. Its b>=0 hypothesis agrees with the actual radical
interpretation, and the closed endpoints cause no defect.

## 4. Both full-interval adjusted polynomials

For the strengthened OUTER polynomial, put q=2u-u²>=0. The claimed
expression is424u+2000b-80q. On0<=u<=2/5, the half-radical
bound and u²<=u already bound264u+80u²+2000b by1137.6<1201.
On2/5<=u<=1, `circle_upper` with212,1000,1023 gives
424u+2000b<1235 because1023²-212²-1000²=1585>0.
The interval product (5u-2)(8-5u)>=0 gives25q>=16, so the
adjusted expression is<1235-80*16/25=1183.8<1201.
This is precisely loss-(.08/2)q<.6005, on the WHOLE closed interval.

For the INNER polynomial698u+61u²+4000b<2402, all three ranges
are covered with strict slack:

* On0<=u<=1/2, (1-2u)(1+2u)>=0 gives u²<=1/4. Its upper
  bound is349+61/4+2000=2364.25<2402.
* On1/2<=u<=7/10, the interval product gives20u²<=24u-7.
  The circle data1928,10000,10185 have positive gap17041 and give
  3856u+20000b<12113. Combining these yields an upper bound
  48025/20=2401.25<2402, equivalent to the written .6003125<.6005.
* On7/10<u<=1, (2u-1)(13u-9)>0. Its exact expansion gives
  (3-u)²-25b²>0; since3-u>0 and b>=0, the code's product
  contradiction proves5b<3-u. Using u²<=u then bounds the target
  strictly below2400-41u<=2400<2402.

At1/2 and7/10 the preceding closed ranges apply;0 and1 are also
included. The later strict split is used only where its factors are
strictly positive. There is no missed boundary or division by a possibly
vanishing variable.

`individual_adjusted_loss` retains EACH foot's own coefficient. In the
outer case1000beta>=80, so multiplying the nonnegative q by
1000beta-80 transfers the strengthened polynomial. In the inner case
2000beta=61 and the scale4000 give exactly the inner polynomial.
Both derive2000loss-1000beta*q<1201; that conclusion is not supplied
as an input. The actual outer loss<.212u+b and inner loss<.205u+b
are the explicitly external local-pressure interfaces in a35d.

## 5. Finite list, one square and exact rational reference

`Foot` stores interval/root/individual-role data and beta>0, without an
assumed gain. For one foot the product of beta with the derived strict
scalar inequality is exactly

    credit-debit
      =1000beta²q-beta(2000loss-1201)>0.

`finite_gain` correctly proves nonnegativity for every finite list and
strict positivity for every nonempty list by structural induction.
The empty list sums to0. The algebra needs neither nonadjacency nor a
bound on list length. This does NOT manufacture a geometric list or
justify arbitrary repetitions: the actual list and its simultaneous
pressure/norm estimates must be supplied by the application.

The cleared `one_pressure` input is exactly
pressure>=base-sum beta(loss-.6005). The cleared `one_norm` input is
normZ<=base-sum beta²q. One nonnegative square then gives

    energy>=2pressure-normZ
          >=base+(sumCredit-sumDebit)/1000.

The nonempty case is therefore strict. In the empty case, the explicit
`empty_strict` pressure and normZ<=base give the required strictness.
All scaling factors1000/2000 and their signs match a35d's ONE-vector
argument. No independent local energy gains have been added.

`reference_energy` calculates the threshold from the reference data:
3759²+5(3881-3759)²=14130081+74420=14204501. The denominator
is5*2000²=20000000. Thus base=.71022505 is a DERIVED equality,
not a supplied final A12 floor. The norm identity/reference coordinates
themselves remain geometric/application premises.

## 6. Final contradiction and exact external boundary

`typed_seven_budget_contradiction` combines the strict A12 bound with
the ACTUAL defect>.09, the separately proved DEFECT-CONDITIONAL
rank-two floor>.2, a same-vector orthogonal sum and original norm<1.
The exact sum is20004501/20000000=1.00022505>1.
`complete_typed_seven_scalar_contradiction` calls the individual/list
completion and reference calculation before that contradiction. It does
not assume the final A12 threshold or the desired norm contradiction.

The map accurately retains the following external obligations: compact
acute all-wall geometry; the original farthest vertex/vector; fresh N13
finite/foreign/dotted inventory; actual K7 lock and helper field/norm/
trace exclusions; strong finite-cap and local radical-pressure bounds;
the actual reference coefficients; the pressure/norm bridges, including
the favorable Gram sign for adjacent feet; empty-list strictness; and
orthogonality of the three WHOLE original components. The actual N7
and O7 meanings and their DISTINCT propagation proofs are retained.
A bare old7 or universal ordinary rank-two>.2 assertion would not
instantiate these hypotheses. No disconnected root-length squareclass,
arbitrary face arithmeticity or extra first-type subclass is imposed.

## 7. Static Lean and execution status

Only Std is imported. Metric quantities and finite-list entries live in
the generic linearly ordered commutative ring R with the displayed
`OrderedRing` compatibility; Nat is not substituted for metric values.
The local numeric casts and OrderedRing square/product APIs follow the
existing verified scalar-kernel pattern. No axioms, sorry, native_decide,
new dependency import or hidden real/number-field instance occurs.
The List definitions/induction, projections of `Foot`, theorem calls and
ten footer names are consistent on static inspection. I found no
specific elaboration/API blocker. Actual tactic elaboration remains a
bounded-compiler obligation; static algebraic validity is not a runtime
success claim.

There are exactly TEN requested theorem axiom prints and ZERO actual
outputs credited at this precompile gate. The source/map truthfully say
prepared source only. Even after compilation, foundational axiom prints
will not enumerate the explicit typeclass and application hypotheses.
This kernel is a conditional scalar formalization of the accepted
component theorem, not an instantiated geometric proof, a formal proof
of all twelve frontier groups, arithmetic dimension15, or a new example.
No target source, map, accepted proof, paper or shared guide was edited
by this review.
