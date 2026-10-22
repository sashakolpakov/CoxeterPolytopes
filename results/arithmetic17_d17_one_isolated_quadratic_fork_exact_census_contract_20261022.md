# Whole D17: complete exact domain for ONE isolated quadratic fork

22 October2026. **Prepared mathematical algorithm/worker contract, pending
full independent mathematical and code review. UNEXECUTED.** No worker,
import, enumeration, symbolic/scientific process, test, Lean, compiler or
TeX ran. There is no census output, empty-domain claim or bound change.

The companion standalone stdlib source is
`scripts/check_d17_one_isolated_quadratic_fork_20261022.py`. Its intended
units are batch0, batch1, batch2 (twelve fields each), then a SEPARATE
combine. Every eventual run belongs to root's existing guarded serial
runner; this contract neither launches nor authorizes an uncapped run.

## 1. Exact conditional geometry and literal input scope

Assume an ACTUAL whole-D17 farthest vertex of a compact arithmetic
Coxeter17-polytope, with ACTUAL QUADRATIC field K. Assume ONE chosen
actual outgoing fork opposite is ISOLATED and DOTTED within D17.
No assumption about the OTHER fork is made by this census.

These are the three fixed mathematical inputs (paths under results/):

| ID | Source | SHA256 | ONLY interface used |
|---|---|---|---|
| D | arithmetic17_whole_d17_dotted_head_degree_three_and_nonisolated_quadratic_candidate_20261022.md | 4baa131a7609f1af7254e154fae99834f0ab38f9e7a31bc2a49a186b0916ceba | Actual raw isolated-fork normalization, LONG and TWO cyclic integrality and fresh W=V/2 integrality, conditioned here on this fork's isolation |
| H | arithmetic17_whole_d17_fresh_long_head_and_all_finite_head_exclusion_candidate_20261022.md | 9e2cf7d75f12c73b69e0bc9c39044acb56c37be36f94fe707efb808ee5ae9eb6 | Actual isolated dotted long head, U strips, both C_f<13; no finite-head fork conclusion |
| K | arithmetic_uniform_classical_weak_helper_obstruction_20261021.md | 535ac0d0fa92d28ccab2976fa147a4ccff0bd3bfcd1e04ca4e134c6ae6857905 | Actual admissibility, K-root lines and invariant integral Cartan cycles |

The newly prepared quadratic fork-isolation theorem0c5b and cubic
exclusion851 are NOT mathematical premises of this ONE-isolated-fork
census. A later full-D17 composition would additionally need their
complete independent reviews/root acceptance and a successful census.

At the SAME original vertex, the head and chosen fork have

    U,W in O_K, K=Q(U)=Q(W),
    4<U<11, 0<U'<2,
    2<W<676/17, 0<W'<4/17.

The ACTUAL raw head has norm2U and sole old product-U. The ACTUAL raw
fork has norm4W and sole old product-2W. Its time squares and spatial
product with the head are

    H0²=U(U-2), Hf²=W(17W-4), spatial product UW.

Thus their actual mutual raw product is

    R=U(U-2)W(17W-4), g=UW-sqrt(R),               (1)

where sqrt(R) is the PHYSICAL POSITIVE root on the actual common K-time
line. In particular R must be a square IN K, not merely have rational-
square norm. The LONG cycle makes g integral. The TWO-cycle makes

    L=g²/(2UW) in O_K.                            (2)

At the other real place the two distinct actual normals are independent
positive vectors, giving0<L'<4. They cannot be proportional since their
different sole old contacts persist under embedding. In the physical
domain of(1), g is strictly negative: the difference
R-(UW)²=2UW[(8U-17)W-2U+4] is positive, since its bracket is
>(8U-17)*2-2U+4=14U-30>0. Hence L>0 too.

The physical product L is either4cos²(pi/m), m>=3, or>4 (dotted).
L=4 would be a parallel pair, incompatible with compactness. The
complete quadratic finite alphabet is1,2,3, (3+sqrt5)/2,
(5+sqrt5)/2,2+sqrt2,2+sqrt3, retaining a value only in its actual
field. Indeed degree of4cos²(pi/m) is phi(m)/2 for m>=3;
phi(m)<=4 gives exactly m3,4,5,6,8,10,12. The elementary totient
argument permits only primes2,3,5 and the exponents/combinations giving
1,2,3,4,5,6,8,10,12. No continuous value in(0,4) is silently allowed.

Every actual pair satisfying this conditional geometry must pass all
these tests. Passing is only NECESSARY, not wall/polytope realization.

## 2. Complete fields: derived fundamental discriminant D<121

Use the FULL quadratic integer basis1,omega of the ACTUAL field. If
d>1 is squarefree, the field discriminant is D=d and
omega=(1+sqrt(d))/2 for d1 modulo4; otherwise D=4d and omega=sqrt(d).
For U=a+b omega, its two strips imply b>0 and

    b sqrt(D)=U-U'<11, hence D<121.

Thus the field list is derived, not prescribed. A finite generator loops
over squarefree integers2<=d<=120, forms their fundamental D, retains
D<121 and sorts. The exact expected36 discriminants are

    5,8,12,13,17,21,24,28,29,33,37,40,
    41,44,53,56,57,60,61,65,69,73,76,77,
    85,88,89,92,93,97,101,104,105,109,113,120.

The worker compares the generated list with this explicit list, then
selects consecutive disjoint twelve-field blocks. No squarefree d outside
the finite loop can have fundamental D<121.

## 3. Complete original-ring coefficient domains, with exact floors

For each field let a0=floor(-b omega'). For U the interval
0<a+b omega'<2 contains EXACTLY the two possible integers
a=a0+1,a0+2, since b omega' is irrational. Retain precisely those
with4<U<11. Necessarily b²D<121; since D>=5, b<=4.

For W, the interval0<a+b omega'<4/17 of width<1 contains AT MOST
the integer a=a0+1. Check that upper inequality and2<W<676/17
exactly. Necessarily289b²D<676²; since D>=5, b<=17.

The exact floors need no numerical roots. If omega=sqrt(d),
a0=isqrt(b²d). If omega=(1+sqrt(d))/2,
a0=(isqrt(b²d)-b)//2. Both formulas follow from b sqrt(d) being
irrational between two consecutive integers. These are the FULL O_K
domains; they do not replace a half-integral basis by Z[sqrt(d)].

Consequently each field has at most8 U and17 W candidates, hence at
most136 pairs. The TOTAL36-field domain has at most4896 pair tests.
This is a static upper bound, not an executed candidate count.
An empty field is still emitted and counted as a completed field.

## 4. Exact arithmetic and COMPLETE K-square test

All arithmetic uses Fraction pairs a+b sqrt(d), with positive sqrt(d).
Products, conjugation and division use the usual rational formulas.
To determine a mixed-sign pair's sign, compare a² and d b²; a nonzero
rational pair in a genuine quadratic field has zero norm only when
both entries vanish. Same-sign/rational/pure-radical branches are direct.

An element a+b sqrt(d) is integral exactly when a,b are integers for
the sqrt basis, or2b and a-b are integers for the half basis.
Checking rational norm or trace alone is insufficient and is not used
as the definition of integrality.

For rational b0, a field square is either a rational square x²=a or a
pure-radical square d y²=a. Test BOTH with exact numerator/denominator
integer square roots, including zero and negative a.
For b nonzero, a field square requires N(a+b sqrt(d))=n² with a
nonnegative rational n. Test BOTH signed norms+n and-n: one must have

    x²=(a+/-n)/2 a NONZERO rational square,
    y=b/(2x), (x+y sqrt(d))²=a+b sqrt(d).

These conditions are also sufficient. If the rational half-sum is zero
while b is nonzero it cannot give a root, so that branch is skipped.
All returned roots are checked by exact squaring, then their sign is
normalized to PHYSICAL POSITIVE. That choice is mandatory in(1).
There is no norm-sign shortcut, floating arithmetic or rational-only
square test. Prepared arithmetic fixtures cover both signed-norm mixed
branches, a pure radical, a rational square, a rational nonsquare, a
negative element, zero, sign cancellation and the half-integral basis.
They are unexecuted until a root-owned worker unit actually runs.

## 5. Streaming certificate, conjunction and resource contract

The CLI requires --expected-source-sha256. Every unit validates that
self pin and the three fixed mathematical pins plus this contract pin
BEFORE computation and AGAIN before its terminal PASS summary. Hashes
are streamed in65536-byte chunks. A destination must be a NEW file
under results/; exclusive creation prevents overwriting existing evidence.

Commands have the following interface (the printed SOURCE_SHA is a
placeholder for the eventual externally frozen exact source hash):

    python3 scripts/check_d17_one_isolated_quadratic_fork_20261022.py \
      --expected-source-sha256 SOURCE_SHA batch --index 0 --output results/NEW0.jsonl

Indices1 and2 have the same syntax. The separate combine takes exactly
three --part paths and corresponding --part-sha256 values, binds their
CURRENT bytes and requires one complete batch of each index0,1,2.
It rejects wrong source/pins, wrong fields, duplicate/out-of-order fields
or survivor indices, truncated/non-PASS summaries, inconsistent counters
and malformed/unverified survivor witnesses.

Each field is processed alone. The output writes EVERY survivor as an
individual JSONL record containing its field/full-basis coefficients,
exact U,W,H0²,Hf²,R,positive sqrt(R),g,L and the native/dotted class.
Every field emits U/W/pair and rejection counts. All necessary failures
are counted by the FIRST failing stage in this fixed order:
R not a field square; g nonintegral; L nonintegral; nonphysical L outside
(0,4); compact-parallel/non-native physical L. No survivor is truncated.

The combine streams each input, rechecks every emitted survivor and its
exact domain indices/algebra and copies EVERY survivor to its own output.
It validates all36 fields and all counter sums; it does NOT pretend to
independently rerun the rejected pair tests. Their complete coverage is
provided by the pinned batch source, actual successful batch executions
and their root-owned guarded receipts. Three valid recorded outputs and
a combine are a conjunction, not a claim of independent regeneration or
a cryptographic proof that arbitrary supplied JSON was really executed.

Every input JSONL record is limited to65536 characters, including its
newline; malformed overlong/truncated records are rejected with bounded
reads. Only one field's at most8+17 candidate objects, one JSONL record and
small36-field/counter metadata are retained in memory. No repository
worker imports, scientific package, subprocess, network, parallelism,
unbounded data load or directory scan occurs. The exact rational operands
come from this bounded algebraic domain. Root must use the unchanged
59-second/1GiB/one-thread guard with40GiB/6GiB admission, and retain all
terminal receipts. Splitting into three units does not raise any cap.

## 6. Interpretation remains conditional and unexecuted

At preparation there are ZERO actual units/outputs. A future nonempty
census is a list of necessary actual pair configurations, not mirrors
or examples. A future EMPTY census would exclude ONE isolated dotted
fork over EVERY admissible actual quadratic field at these interfaces.
Whole D17 would additionally need the separately reviewed cubic reduction
and complete quadratic fork-isolation theorem plus a source-bound
mathematical composition accepted by root. No bound guide or paper is
modified by this source/contract preparation.
