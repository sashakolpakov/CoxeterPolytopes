# Whole A17: conditional finite cubic anchor and exact other-endpoint census

22 October2026. New prepared MANUAL domain/halting/code contract by
finite57_composition_audit. **UNREVIEWED and UNEXECUTED.** Root and I
independently derived the cubic square-recovery identity and trace bound;
root supplied the short Eisenstein maximal-order proof below. No worker,
import, enumeration, symbolic process, test, compiler, Lean, TeX or Git
operation was used. There are zero actual census outputs at preparation.

Companion stdlib source:
`scripts/check_a17_finite_cubic_endpoint_20261022.py`.
Two separate field units, native7 and native18, followed by one separate
combine. All eventual execution belongs to root under the unchanged guard.

## 1. Exact CONDITIONAL interface and scope

Assume an ACTUAL whole-A17 farthest vertex of a compact arithmetic
Coxeter17-polytope. Assume BOTH actual endpoint opposites isolated and
ONE chosen endpoint finite native7 or18. Relabel the old path so that
this endpoint is first; this is coordinate relabelling, not an assumed
symmetry of the actual chamber. Set U=4gamma_anchor^2, V=4gamma_other^2.
The actual field and physical anchor are respectively

    native7: K=Q(U), P7(U)=U^3-5U^2+6U-1=0;
    native18:K=Q(U), P18(U)=U^3-6U^2+9U-3=0,

where U is the LARGEST root. The other endpoint is allowed finite OR
dotted, with the conditional strict strips

    V in O_K, 36/17<V<3969/425,
    0<V_sigma<36/17 at BOTH other embeddings.

The direct mathematical input pins (paths under results/) are:

| Source | SHA256 | Exact interface |
|---|---|---|
| arithmetic17_whole_a17_short_endpoints_isolation_and_golden_exclusion_candidate_20261022.md | b80c27b00dc7ecdbdf7201327d620e10b372e35f1047964c9ae204c2d6996558 | Conditional actual isolated endpoint raw norms, strict V strips and full ACTUAL field generation |
| arithmetic_uniform_classical_weak_helper_obstruction_20261021.md | 535ac0d0fa92d28ccab2976fa147a4ccff0bd3bfcd1e04ca4e134c6ae6857905 | Integral invariant Cartan cycles, K-root normalizations and all-place admissibility |

The separate two-dot/finite-field reduction and degree-five theorem are
NOT premises of this conditional census. Any later whole-A17 composition
must separately accept them, its quadratic exclusion and this actual
execution. A field or form is not assumed excluded here.

The actual raw roots have norms2U,2V, time squares F(U)/18,F(V)/18,
where F(X)=X(17X-36), and spatial mutual product UV/18. Thus

    R=F(U)F(V),
    g=(UV-sqrt(R))/18, sqrt(R) PHYSICALLY POSITIVE,
    g integral, L=g^2/(UV) integral.               (1)

The long Cartan cycle through A17 is +/-g. The two-cycle is L, and
all roots are normalized on the SAME actual K-time line. Hence R must
be a square IN K; a rational-square norm alone is insufficient.
Acuteness requires g<=0. **g=L=0 remains allowed.** When g is nonzero,
both other conjugates satisfy0<L_sigma<4, by positive definiteness
and the distinct sole old endpoint contacts. No strict positive bound
is imposed on the zero mutual cycle.

Every physical own V and mutual L must be native finite or dotted.
The COMPLETE finite alphabet in these cubic fields is0,1,2,3 and

    K7: native7 U; native14 U^2-3U+3;
    K9: native18 U; native9 U^2-4U+4.

Value0 is native2. Value4 is parallel and forbidden by compactness;
value>4 is dotted. Indeed phi(m)/2 must divide3, giving rational
m3/4/6 or cubic m7/9/14/18; quadratic finite values cannot lie in a
cubic field. K7 and K9 are distinct, since their field discriminants
are49 and81 as proved next. The displayed expressions follow from
2cos(3theta)=(2cos(theta))^3-3(2cos(theta)) or the double-angle identity.
No continuous nonnative value below4 is admitted.

## 2. FULL integer basis proved without a field table

The polynomial discriminants are49 for P7 and81 for P18, by the
standard cubic discriminant formula. Also

    P7(X+4)=X^3+7X^2+14X+7,
    P18(X)=X^3-6X^2+9X-3

are Eisenstein at7 and3. In particular the fields are cubic.
The order index [O_K:Z[U]] can have only the corresponding prime divisor,
since its square divides the polynomial discriminant. For native7 use
pi=U-4; for native18 use pi=U. At an extension of the normalized p-adic
valuation, the Eisenstein equation forces v(pi)=1/3: first v(pi)>0,
then its constant and leading terms must share the least valuation1.

If p divided the index, O_K/Z[pi] would have an element of order p,
represented by(A+Bpi+Cpi^2)/p with A,B,C integers not all divisible
by p. The nonzero terms in its numerator have valuations with distinct
fractional parts0,1/3,2/3. The minimum is unique and strictly below1,
because some coefficient is a p-adic unit. Its numerator therefore has
valuation<1, so the quotient has negative valuation and is not integral.
Contradiction. Consequently

    O_K=Z[U] in BOTH fields.                      (2)

The elementary index/discriminant and local-integrality facts are stated
number-theoretic foundations; no maximal order is guessed from a monic
polynomial, CAS output or named field list.

## 3. Complete coefficient boxes:2432 scalar triples

Write the OTHER endpoint in the FULL basis(2) as V=a+bU+cU^2.
The three anchor roots have the following disjoint rational brackets,
ordered from smallest to largest:

| Field | Small | Middle | Physical |
|---|---|---|---|
|native7|(.19,.20)|(1.55,1.56)|(3.24,3.25)|
|native18|(.46,.47)|(1.65,1.66)|(3.87,3.88)|

All endpoints are exact hundredths. Polynomial endpoint values, multiplied
by10^6, are respectively

    P7: -33641,8000; 11375,-11584; -35776,15625;
    P18:-32264,8423; 7125,-19304; -70797,4672.

The three sign changes of a cubic prove exactly one root in each bracket.
In Lagrange interpolation the denominator at a root r_i is
P'(r_i)=product_(j!=i)(r_i-r_j). Its absolute lower bounds at the
PHYSICAL,MIDDLE,SMALL roots are respectively

    native7:5.1072,2.268,4.104;
    native18:7.514,2.6078,4.012.

Use physical0<V<10 and other0<V_sigma<9/4. The coefficient of X^2
in the interpolant is the sum of V_i/P'(r_i). The coefficient of X
is minus the sum of V_i*(sum of the other two roots)/P'(r_i), and
the constant uses the product of the other two roots. Triangle inequalities
and the displayed brackets give the following strict term bounds in
physical,middle,small order:

| Field/coefficient | Three absolute term bounds | Sum below |
|---|---|---|
|native7,c|2;1;.55|4|
|native7,b|3.45;3.43;2.64|10|
|native7,a|.62;.65;2.78|5|
|native18,c|1.34;.87;.57|3|
|native18,b|2.84;3.76;3.11|10|
|native18,a|1.04;1.58;3.62|7|

Each bound is a positive rational multiplication against its displayed
denominator lower bound. For example native7's last constant term is
below(9/4)*3.25*1.56/4.104=11.4075/4.104<2.78; native18's last
one is below(9/4)*3.88*1.66/4.012=14.4918/4.012<3.62.
Thus the COMPLETE integer boxes are

    native7: a=-4..4,b=-9..9,c=-3..3:1197 triples;
    native18:a=-6..6,b=-9..9,c=-2..2:1235 triples.

Their total is2432. Every triple is tested against the EXACT physical
and BOTH conjugate strips and the full native own condition. Negative
b and c are included. Neither trace/norm conditions alone nor separate
coordinate membership replaces the full original integer basis.

## 4. Exact field signs, arithmetic and complete cubic K-square test

Represent elements by three Fraction coefficients modulo the relevant
monic cubic. Convolution followed by descending reductions of powers4
and3 gives products. The trace is3a+5b+13c or3a+6b+18c. The norm
is the determinant of the3-by-3 multiplication matrix. Put
e2=(Tr(x)^2-Tr(x^2))/2; the inverse is
(x^2-Tr(x)*x+e2)/N(x), by Cayley-Hamilton. Zero is checked by all
three coefficients, and integrality by their integer denominators in(2).

Signs at EACH embedding use its separate isolating bracket. Bound
a+b r+c r^2 on that positive rational interval with exact rational
interval arithmetic, then bisect its unique anchor root until the value
interval excludes zero. A nonzero quadratic polynomial cannot vanish at
a cubic generator, proving termination. More quantitatively, if q is a
common coefficient denominator and M=|a|+4|b|+16|c|, its nonzero norm
has magnitude>=q^(-3), so every embedding value has magnitude at least
1/(q^3 M^2). Interval error tends to zero under bisection. No decimal
evaluation, numerical root finder or symbolic package is used.

For the square test R in O_K is positive at ALL three embeddings.
If R is rational, a square root in an odd cubic field must be rational:
otherwise it generates a degree2 subfield, contradicting the tower law.
Test its numerator and denominator integer squares, including zero.

For nonrational R compute A=Tr(R), B=e2(R), N(R). If s^2=R in K,
integral closedness and(2) make s integral. Its norm is one of
n=+/-sqrt(N(R)), an INTEGER. Its integer trace t and e2(s)=z satisfy

    z=(t^2-A)/2 integral, z^2-2nt=B.

Conversely every such n,t gives the only possible root

    s=(tR+n)/(R+z).                              (3)

Indeed s^3-ts^2+zs-n=0 implies s(R+z)=tR+n. Here R+z cannot
vanish because R is nonrational. EXACT squaring of(3), integrality and
the physical sign are still checked; the scalar identities alone are
not declared sufficient. BOTH signed norms and ALL possible integer
traces must be tested.

There is a small complete trace bound. The finite anchor U<4 gives
0<F(U)<128; V<10 gives0<F(V)<1340. Thus physical sqrt(R)<415,
since128*1340=171520<415^2. At either other place,0<U,V<36/17
and -20<F(U),F(V)<0, so0<R_sigma<400 and every square-root
conjugate has absolute value<20. Hence |t|<455. Test precisely
t=-454,...,454 and both norm signs. This is909*2 scalar trace checks
per retained endpoint, not an unbounded search. An integer-square norm
is only a prefilter; full K-square recovery and exact verification follow.

Prepared fixtures cover multiplication/reduction, all three separate
embedding signs, original basis integrality, rational and zero squares,
a nonsquare rational, and mixed cubic square roots of negative norm.
They remain UNEXECUTED until a root-owned unit runs.

## 5. Streaming complete units, combine, halting and future interpretation

The source requires --expected-source-sha256. Every unit streams hashes
of itself, the two math pins, this contract, copied I/O/helper provenance
and FOUR original guard/monitor/resource sources before work and again
before terminal PASS. No predecessor checker is imported or run.

Unit native7 or native18 emits its header, one record for EVERY scalar
triple in its lexicographic box (including exact first failures), and a
terminal exact counter summary. Successful NECESSARY pair records contain
all three embeddings' source parameters, exact positive K-square witness,
g,L and native/dotted/zero classes. A new exclusive-created file under
results/ preserves partial failures and earlier evidence. JSONL lines
are limited to65536bytes and each output to16MiB; overflow FAILS without
PASS or silent truncation. Stdout is only a small terminal summary.

The separate combine requires TWO distinct bound parts and one owner of
each field. It rehashes them before/after reading, recomputes EVERY exact
triple record using the SAME predicates, compares canonical JSON and
requires terminal PASS/EOF/counter/source identity. It copies every
survivor; no tie or candidate is selected away. This is complete replay,
not an independent mathematical proof or evidence that arbitrary input
files were genuinely run. Root's actual guarded ledgers are still required.

Every coefficient box and trace loop is explicitly bounded. Exact sign
bisection terminates by the preceding nonzero-field argument. Operations
use finite exact rational operands. Memory retains one triple, a constant
number of3-coefficient field elements, one bounded JSONL record and small
counters; no complete file, field-candidate array, scientific import,
subprocess, parallelism or network access is used.

Guarded command TEMPLATE (root substitutes frozen source SHA and NEW paths):

    python3 -B scripts/run_bounded_research_check_20261018.py \
      --parent-lock --ledger results/NEW_A17_C7_execution.json -- \
      scripts/check_a17_finite_cubic_endpoint_20261022.py \
      --expected-source-sha256 SOURCE_SHA batch --field native7 \
      --output results/NEW_A17_C7.jsonl

Native18 uses the same syntax. A separate guarded combine has
`--expected-source-sha256 SOURCE_SHA combine`, TWO `--part PATH`
and corresponding `--part-sha256 SHA`, and a new `--output`.
There is NO extra python3 after --. The unchanged shared-lock serial
guard keeps59s/1GiB/one thread,40GiB/6GiB admission,8GiB/3GiB abort.
No prediction of runtime fit, admission, actual execution or result is
made by this prepared source.

An eventual empty authenticated census would exclude ONLY this conditional
native7/18-anchored isolated pair. Whole-A17/arithmetic16 would need its
separately reviewed complete geometric/degree/finite-alphabet reductions,
the quadratic census, this actual execution and a full root/other-author
composition. A nonempty result is a necessary model list, not an example.
