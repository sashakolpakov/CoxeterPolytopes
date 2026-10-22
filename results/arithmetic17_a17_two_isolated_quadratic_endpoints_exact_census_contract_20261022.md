# Whole A17: conditional exact census of TWO isolated quadratic endpoints

22 October2026. Prepared source/algorithm contract by
finite57_composition_audit. **UNREVIEWED and UNEXECUTED.** The endpoint
isolation source below is itself pending complete outside review and root
acceptance. This contract assumes its stated algebraic interface explicitly;
it does not accept that theorem, execute a census or change any bound.

The new standalone stdlib worker is
`scripts/check_a17_two_isolated_quadratic_endpoints_20261022.py`.
It has three disjoint field units (9,8,8 fields) and one separate combine.
All future execution belongs to root under the existing serial guard.

## 1. Conditional actual geometry, not a field-list assumption

Assume an actual whole-A17 farthest vertex of a compact arithmetic
Coxeter17-polytope, ACTUAL QUADRATIC field K, and BOTH actual endpoint
opposites isolated within the old A17. At the SAME original vertex assume
the exact endpoint bounds and normalizations below. Arithmetic reflection
lattices have the admissible quadratic-form model automatically; there
is no extra first-type subclass restriction here.

These are the two direct mathematical input pins (under results/):

| Source | SHA256 | Exact interface |
|---|---|---|
| arithmetic17_whole_a17_short_endpoints_isolation_and_golden_exclusion_candidate_20261022.md | b80c27b00dc7ecdbdf7201327d620e10b372e35f1047964c9ae204c2d6996558 | CONDITIONAL isolated raw endpoint formulas, actual field generation and strict strips; no acceptance of its geometric reduction is supplied by this census |
| arithmetic_uniform_classical_weak_helper_obstruction_20261021.md | 535ac0d0fa92d28ccab2976fa147a4ccff0bd3bfcd1e04ca4e134c6ae6857905 | Actual K-root normalization, integral invariant Cartan cycles and nonphysical positive definiteness |

For each endpoint i=1,2, put U_i=4gamma_i^2. The conditional interface is

    U_i in O_K, K=Q(U_i),
    36/17<U_i<3969/425, every other0<U_i'<36/17.

The raw endpoint root r_i has square2U_i, sole old product-U_i and
positive time coefficient with square U_i(17U_i-36)/18. The two raw
spatial projections have mutual product U_1U_2/18: the inverse raw
A17 Gram entry between its endpoints is1/18. Therefore the ACTUAL raw
mutual product is

    F_i=U_i(17U_i-36), R=F_1F_2,
    g=(U_1U_2-sqrt(R))/18.                       (1)

The square root is PHYSICALLY POSITIVE. It must be IN K because both
actual raw normalizations lie on the same one-dimensional K-time line.
This is a full field-square condition, not merely a square rational norm.
The complete long Cartan cycle through the old A17 path is +/-g:
the endpoint-to-old entry is-1, the old forward entries are-1, the
last old-to-endpoint entry is-U_2 and the closing entry is g/U_2.
Hence g is integral. The actual two-cycle gives

    L=g^2/(U_1U_2) integral.                     (2)

Acuteness requires g<=0. **Orthogonality g=L=0 is allowed.** If g is
nonzero, its conjugate is nonzero, and positive definiteness of the two
independent conjugate normals requires0<L'<4. Their different sole old
contacts ensure independence at every embedding. For g=0 use L=0 at
BOTH places, without imposing a strict positive lower bound.

Every physical U_i and L must be a native finite invariant or dotted:
L=0 is native2; the complete positive quadratic finite alphabet is

    1,2,3,(3+sqrt5)/2,(5+sqrt5)/2,2+sqrt2,2+sqrt3,

retaining the radical values only in their actual fields. Indeed
4cos^2(pi/m)=2+2cos(2pi/m) has degree phi(m)/2 for m>=3, so the
quadratic condition gives m3,4,5,6,8,10,12. The elementary totient
list is complete: its prime divisors can only be2,3,5 and the allowed
exponents/products give precisely these m. The rational m2 value is0.
Value4 would be an actual parallel pair, excluded by compactness;
values>4 are dotted. The same alphabet applies to each actual endpoint
own U_i, so isolated finite owns are retained, not silently declared
dotted. A nonnative continuous value in(0,4) is rejected.

No endpoint equality is assumed. The census retains ALL ORDERED pairs,
including equal parameters and g=0. The golden field is included in the
domain even though the pending geometric source separately excludes it.
Its inclusion makes no stronger input necessary.

## 2. Complete quadratic fields and FULL integer coefficient boxes

Write K=Q(sqrt(d)), d>1 squarefree. Its fundamental discriminant and
full integral basis are D=d, omega=(1+sqrt(d))/2 for d=1 modulo4,
or D=4d, omega=sqrt(d) otherwise. For U=a+b omega the two strips imply

    0<b sqrt(D)=U-U'<3969/425,
    D<(3969/425)^2<88, b<=4.

Thus the complete25 fundamental discriminants are

    5,8,12,13,17,21,24,28,29,33,37,40,
    41,44,53,56,57,60,61,65,69,73,76,77,85.

The source generates the list from squarefree2<=d<=87, forms D and
retains D<88, then compares it with this exact list. No larger d can
have fundamental discriminant below88. Units own the consecutive blocks
[0:9], [9:17], [17:25], so every field is assigned once.

For b1,...,4 retain only b^2D<(3969/425)^2. Put
a0=floor(-b omega'). Since the interval0<a+b omega'<36/17 has width
strictly between2 and3, its only possible integer coefficients are

    a=a0+1, a0+2, a0+3.                          (3)

THREE candidates are essential: the strip has width>2. All strict
physical/conjugate inequalities are checked exactly. The exact floors
are isqrt(b^2d) in the sqrt basis and (isqrt(b^2d)-b)//2 in the half
basis. Irrationality makes either formula exact even when b is even.
This uses the FULL ring O_K, including its half-integral basis.

Each field has at most12 strip candidates and144 ordered pair trials.
Across25 fields there are at most300 coefficient candidates and3600
pair trials. These are STATIC upper bounds, not observed counts.
Each candidate's physical own type is checked by the preceding complete
native/dotted alphabet. No extra pruning, field assumption, endpoint
distinction or omitted zero branch is used.

## 3. Exact arithmetic and complete square test

The worker uses Fraction pairs a+b sqrt(d). Signs use exact norm
comparison in mixed-sign cases. Integrality means integer coefficients
in the sqrt basis, or integer2b and integer a-b in the half basis.
Rational norm/trace alone are NOT integrality tests.

For b=0, a+b sqrt(d) is a K-square exactly if a is a rational square
or a/d is a rational square. This includes zero and pure radical roots.
For b nonzero, let n be the nonnegative rational square root of the
norm, if it exists. Test BOTH signed norms+n and-n:

    x^2=(a+/-n)/2 a nonzero rational square,
    y=b/(2x), (x+y sqrt(d))^2=a+b sqrt(d).

These conditions are necessary and sufficient. A zero half-sum cannot
produce a mixed root. Exact squaring verifies each returned root; its
physical sign is made positive. Both mixed signed-norm branches are
required. Local helper source was copied and adapted from the separately
pinned D17 worker73c44f94; it is NOT imported or used as a D17 theorem.
Prepared arithmetic fixtures also cover g=L=0 native2 and the three-a
coefficient case; they are unexecuted source until root runs a unit.

For EVERY ordered pair the fixed first-failure order is: endpoint own
not native/dotted; time not a field square; g>0; g nonintegral;
L nonintegral; nonzero L' outside(0,4); physical mutual product not
native/dotted. A survivor contains its exact U_1,U_2,F_1,F_2,R,
positive sqrt(R),g,L and both own/mutual classes. A zero mutual product
is recorded explicitly as orthogonal/native2. No survivor is truncated.

## 4. Certificate, complete combine, termination and boundedness

Each unit requires --expected-source-sha256, then validates that self
pin, both mathematical pins, this contract, copied-helper provenance
and all FOUR existing guard/monitor/resource-source pins BEFORE its
computation and AGAIN immediately before terminal PASS. All hashes stream
65536-byte chunks. A destination must be a NEW file under results/;
exclusive creation preserves earlier success/failure evidence.

Each field emits its complete coefficient-model list, one record for
EVERY ordered pair (including rejected pairs and their exact first
reason), then a field counter summary. Empty fields are emitted too.
Every certificate ends in an exact PASS summary. A bounded JSONL writer
refuses over65536-byte records or over16MiB output; such refusal leaves
the partial evidence without a PASS. These are failure bounds, not
permission to drop a row or survivor.

Combine requires exactly three distinct --part paths with their exact
--part-sha256 values, verifies each part's bytes before and after reading,
and checks exactly one unit of every index. It streams the complete
expected field/ordered-pair sequence, RECOMPUTES every exact pair record
using the same predicates, compares its full canonical JSON and copies
it to the combined certificate. It validates terminal summaries, every
counter, EOF and all25 fields. This is complete certificate replay using
the same algorithm, not an independent mathematical proof or permission
to credit outputs without root's authentic guarded execution ledgers.
Wrong/missing/reordered/duplicated records, changed pins or parts,
oversize/truncated JSON and partial/non-PASS outputs cannot be combined.

All loops have finite displayed bounds:87 possible d, at most10 trial
square divisors per d, four b and three a per field, at most3600 ordered
pair trials and a constant number of exact arithmetic/square branches.
Integer square roots terminate on bounded integer operands. The combine
does at most the same3600 tests. No symbolic routine, numerical root
finder, Pell search, repository import, subprocess, network or parallel
work is used. Memory holds at most12 candidates, one bounded record,
small25-field metadata and exact rational temporaries. The writer streams
to disk, and stdout is only its small terminal summary.

Guarded command TEMPLATE (SOURCE_SHA is the eventual externally frozen
worker SHA; root substitutes NEW immutable output/ledger paths):

    python3 -B scripts/run_bounded_research_check_20261018.py \
      --parent-lock --ledger results/NEW_A17_0_execution.json -- \
      scripts/check_a17_two_isolated_quadratic_endpoints_20261022.py \
      --expected-source-sha256 SOURCE_SHA batch --index 0 \
      --output results/NEW_A17_0.jsonl

Indices1/2 use the same syntax. A separate guarded combine uses
`--expected-source-sha256 SOURCE_SHA combine`, then three `--part PATH`
and corresponding `--part-sha256 SHA` options, and a new `--output`.
There is NO extra python3 after the guard's --. The unchanged guard owns
the shared lock, one thread,59seconds/1GiB process-tree RSS and40GiB/6GiB
admission with8GiB/3GiB abort floors. No runtime/admission is claimed here.

## 5. Scope of a future result

At preparation there are ZERO census units, outputs, compiler runs or
axiom prints. PASS will mean a complete NECESSARY algebraic domain was
checked, not that every candidate is an actual wall pair. A nonempty
result preserves every survivor; it is not an example. An EMPTY result
would exclude these TWO isolated endpoint configurations over all actual
quadratic fields, CONDITIONAL on the stated geometric interface and its
later acceptance. It would not exclude degree>=3, whole A17, arithmetic
dimension17 or a form/commensurability class. The pending endpoint source,
six manuscripts, all earlier checkers/receipts and current guides remain
unchanged. Root plus another author must fully review the source and
domain before any eventual guarded run or acceptance.
