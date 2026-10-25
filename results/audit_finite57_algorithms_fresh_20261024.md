# Finite57: fresh algorithm and reproduction audit

24 October 2026. This is a source-only audit of a coherent part of the
finite-volume57 proof: the complete ideal terminal algorithm, low-eight
concentration, the final arithmetic composition, the abstract residual
certificate, and the current reproduction boundary. It is **not** an audit
of the mathematics of all136 selected scripts or of the complete geometric
higher-dimensional tail.

**Verdict.** No mathematical defect was found in the audited terminal
algorithm and concentration argument at their stated geometric and
computed-output premises. Two actionable reproduction defects were found:
optimization can remove essential assertions, and the public driver's Lean
route lacks the additional process-tree RSS wrapper required by the current
execution policy. These findings do not establish that any historical
execution used optimization, and they do not turn historical receipts into
new failed mathematical computations. They do require a fail-closed public
entry before the revised package is advertised as push-button reproducible.

I did not author or change any audited mathematical source. I used the
current files rather than treating an earlier PASS label as a proof. No
worker, mathematical-module import, enumeration, parser fixture, test, Lean,
TeX, Sphinx, optimizer or symbolic computation was run. Metadata work was
limited to source reading, streaming file hashes, catalogue/receipt decoding
and authentication of stored stdout. The accompanying exact source ledger is
`results/audit_finite57_algorithms_metadata_20261024.json`.

## 1. Exact audit boundary

The primary source is `manuscripts/finite_volume57.tex`, SHA256
`c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0`.
The audited mathematical sections are the residual-certificate proposition,
the ideal envelope and terminal/concentration argument, the final finite
closure/comparison interface, and the common algorithm specification for the
bounded higher-dimensional tail. I read the final endpoint checker and the
low-eight producer in full. I read the current driver and resource runner
in full, the bounded B13/frontier worker's domain declarations and selected
replay/entry functions, and the residual Bellman kernel in full. The
concentration Lean source was inspected for its definitions and premise
boundary, not recompiled or completely re-audited here.

Geometry, the native local census, actual corner coverage, complete finite
forest domination, context isolation and the unbounded tail are separately
audited obligations. The signatures below do not discharge those obligations.
This audit imports no arithmeticity: all finite57 mathematical premises concern
general real hyperbolic Coxeter polytopes.

## 2. Terminal cusp algorithm: definition, halting and exact domain

For a cusp link of total affine rank57, the mathematical inputs are a list of
positive integer factor ranks and the actual-charge envelope

    C_actual <= V * sum_i c(r_i),
    V = product_i(r_i+1).

The interpretation of these quantities as actual cusp fibers is a geometric
premise. The algorithm does not prove that interpretation.

The recursion `high_parts(n,least)` returns sorted partitions of n with all
parts at least `least`. Each non-base call chooses its first part r from the
finite interval `[least,n]` and recursively partitions `n-r` with lower bound
r. The remaining nonnegative rank strictly decreases. A partition's first
part and tail are unique, so induction proves both completeness and absence
of duplicates. The case n=0 returns the unique empty tail.

The terminal domain has one optional part `low` in0..7, followed by a sorted
partition of `57-low` into parts at least8. Zero is omitted. A terminal
profile's positive part below8, if present, determines that outer branch
uniquely. Thus the optional-low construction neither loses nor duplicates
profiles. The producer's independent multiplicity recursion scans every
rank8..57 and leaves a remainder below8; it enumerates the same domain by
the unique multiplicity vector. These are two domain reconstructions, not
two independent proofs of the geometric cusp model.

For each profile, starting from `[1]`, the nested convolution loops implement
exact multiplication by `1+t+...+t^r`. The output length is the new degree
plus one. Induction proves every integer coefficient of

    H(t)=product_i(1+t+...+t^r_i).

Consequently `sum(h)=V`, the symmetry check and the degree check are valid
consistency checks; importantly, the final endpoint reconstructs the full
coefficient list instead of checking just those three invariants. It then
computes `G=sum_j min(j,57-j)h_j`, `kappa=(5G-57V)/2` and the exact rational
normalized and absolute margins. No floating-point arithmetic is used.

The saved computed proposition is2552 profiles, normalized minimum87/10240
and absolute minimum36. This audit proves correctness of the stated finite
algorithm but does not rerun those2552 cases. The preserved low-eight
execution receipt records a successful full producer replay in4.903seconds
with462585856-byte sampled tree RSS. Its output identifies the distinct
minimizers `[9,9,9,15,15]` and `[57]`; minima need not occur at the same
profile. The latter absolute value also follows directly: V=58, G=812,
kappa=377, `58*c(57)=341`, hence surplus36.

The final endpoint independently reconstructs every terminal polynomial and
margin. It checks `c(r)=r/2` through rank8 and c(9)=71/16, but takes the
remaining envelope entries from the pinned input. The printed entries agree
with the manuscript table and `(5r+56)/(r+1)` at ranks24..57. Thus the
endpoint is not an independent derivation of those geometric envelope
constants; that is an explicit input, rather than an extra theorem proved
by the endpoint command.

## 3. Concentration: the complete analytic and termination argument

For independent uniforms on integer intervals, the proposed move replaces
`U_a+U_b` by `U_(a-1)+U_(b+1)`, with `1<=a<=b<=7`. Both laws have the same
mean and are symmetric. On integer thresholds up to a, the difference of
lower stop-loss expectations is

    binom(t+2,3)*(1/[a(b+2)]-1/[(a+1)(b+1)]),

which is nonnegative because the old denominator exceeds the new by
`b-a+1`. On the remaining thresholds up to the midpoint the manuscript's
numerator is

    t(a+b-t)-(a-1)(a+2b+2)/3.

It increases up to the midpoint and at t=a equals
`(a+2)(b-a+1)/3>0`. Symmetry covers the upper half, affine interpolation
covers noninteger thresholds, and equal means dispose of affine terms.
This proves convex order. Conditioning proves preservation under addition
of any common independent summand. This last assertion is necessary when
there are other cusp factors; the source includes it.

The relation `G/V=(57-E|2S-57|)/2` has the correct sign: spreading the law
increases the absolute-value expectation and decreases normalized
curvature. Both new ranks remain at most8, so `sum c(r)` is unchanged.
Therefore normalized margin does not increase. The vertex product decreases
by `b-a+1`, which is positive. No move crosses the nonlinear rank9 entry.

The integer potential `sum r_i^2` increases by `2(b-a+1)>=2`. It starts at
least57 and never exceeds57 squared. There can therefore be at most1596
moves. A profile admitting no move has at most one positive rank below8,
which is exactly the terminal domain in Section2.

The absolute-surplus lift correctly uses positivity. If the original and
terminal pairs are `(V,m)` and `(V',m')`, then

    V>=V', m>=m'>0  =>  Vm>=V'm'>=36.

This would be invalid without positive terminal margin. The source checks
positivity before multiplication. I found no sign, stopping-domain or
normalization gap in this argument. The28 explicit stop-loss rows are
additional finite arithmetic evidence; the written proof is not merely an
unexplained assertion based on those rows.

## 4. Finite Bellman exactness and its application boundary

For a fixed complete context and finite alphabet of positive-rank items,
the residual recurrence considers every allowable last item. It terminates
in increasing rank and computes the exact multiset maximum: remove one item
from any proposed multiset for the upper bound, and append the selected
item to an attaining predecessor for the reverse bound.

A saved table needs base0, *every* transition inequality, and a tight
predecessor at every positive finite stored entry. Positive item rank makes
each predecessor strictly smaller; this prevents fabricated finite states.
An unreachable positive state has no reachable predecessor. Rank0 is the
reachable base and is explicitly exempt from that last equivalence. These
conditions are correctly stated in the manuscript and in
`ResidualBellmanCertificate.lean`. The kernel derives exactness and
unreachability from local conditions, rather than taking exactness itself
as an assumption.

This abstract proof does not establish that an application's residual
alphabet dominates every actual spherical type, that all context-sensitive
counts were fixed, or that no residual pair interaction survives. The
finite producer checks representative base/pair domination, zero residual
crosses and item identities. The manuscript's coarse/actual-fiber split is
the claimed completeness bridge. I inspected that interface and its saved
receipt contracts, but did not freshly audit every originating local table
or every context formula. Their precise mathematical content remains a
separate finite-input obligation; the kernel's compilation cannot replace it.

The final endpoint assumes the predecessor's complete3298 representative
cover and748 actual fibers. It matches the eight old forests and six new
increasing forests exactly, then checks the new totals against377. This is
an arithmetic composition of a supplied complete finite classification,
not a fresh Bellman replay or a demonstration that these formal forests are
realized by chambers. Its own module docstring and output accurately say so.

## 5. Simultaneous-rule and higher-tail interfaces

The endpoint checks that the B13 change has no new coupons or ideal change;
that the guarded D13 repair supplies the old eight exclusions; that the
only increasing rank58 families have been supplied; and that their full
costs are at most377. It adds the global4205 correction once. The cusp
arithmetic supplies surplus36, so the final numerical contradiction is

    0 >= 4205 + 36*v_ideal > 0.

These assertions require the *same actual corner rule* in the face lower
bound, finite maximum and ideal envelope. The endpoint checks interface
values and pins, not the geometric all-corner or retained-consumer
arguments. Root's separate geometry audit must check that identification
and the exhaustive universal routing, including every previous debt.

Similarly, the endpoint's higher-dimensional conclusion takes the complete
bound58 tail as a pinned theorem. For the general tail algorithm, a
concentration step needs nondecreasing envelope cost in addition to convex
order, and the manuscript explicitly states that condition. A decreasing
remaining rank proves finite partition recursion. However, proving that each
advertised terminal family is reached by the allowed moves, and that its
dimension-specific charge corresponds to one complete actual rule, is not
discharged by checking the returned scalar margin. This audit does not
claim a fresh audit of all those cases or of Nikulin's unbounded-tail input.

## 6. What the present catalogue and receipts establish

Metadata authentication derived136 distinct script selections,9 Lean
roots and1 PDF job:146 full jobs and11 core jobs. The actual current
catalogue hash is in the accompanying ledger. No plan command was run.

I reauthenticated all25 successful B13/frontier records in the current
continuation ledger against their retained raw stdout byte lengths/hashes,
batch IDs, exact child command, worker/runner hashes,59-second/1GiB limits,
one-thread setting and zero recorded outer exit. The one earlier
`b13:pins-1` resource refusal remains in the ledger alongside its later
success. This is authentication of recorded execution, not rerunning it.

The25 units are parser+16B13+8frontier. Their equivalence with the two
original producer outputs has its own source/domain contracts and reviews.
The parser unit is a parser contract, not mathematical enumeration. The
111 other selected commands retain inherited separate execution records;
I did not freshly authenticate all111 in this audit. The failed original
44-job run remains failed after its endpoint success and eager-D13 RSS
failure. Neither a completed146-job runtime PASS nor transitive regeneration
is established by joining these metadata statements.

## 7. Reproduction finding R1: optimization removes essential checks

This is a static source finding, not an executed exploit. The public
driver and direct documented commands invoke the frozen bounded research
runner. It copies `os.environ`, does not reject optimization and does not
remove `PYTHONOPTIMIZE` before launching its Python child. Its process-tree
RSS, explicit wall-time, output-size, owner, child-exit and source-binding
checks use `assert`. The endpoint's mathematical comparisons and several
bounded workers also use `assert`.

Under Python optimization those statements disappear. In particular
`assert child.wait()==0` no longer checks a failed child's exit before the
wrapper emits its PASS record. The independent signal alarm and explicit
host resource-floor checks remain; it would overstate the finding to say
that every execution restriction disappears. The RSS rejection and proof
comparisons nevertheless lose enforcement.

The B13/frontier worker's `assert __debug__` at its entry is not a repair:
that assertion itself disappears. D13's newer `if not __debug__: raise`
guard is non-removable, but lies in its own replay entry; importing its
projection definitions from another worker does not call that guard.

Existing receipts do not report `sys.flags.optimize` or the inherited
optimization environment. There is no evidence here that a historical run
used optimization. The correction should preserve their bytes and scope,
rather than silently editing a pinned runner or relabeling old runs.

## 8. Reproduction finding R2 and minimal successor contract

The public driver wraps script checks in the bounded research runner but
launches the Lean runner directly. The Lean runner has explicit host
resource floors, owner/time checks, shared lock and Lean `-M2048 -j1`,
but its `guarded_run` does not sample process-tree RSS. The current policy
requires that additional monitoring; Lean's own allocation option is not
the same measurement.

A minimal successor public entry can preserve every frozen implementation:

1. Before loading proof modules or old drivers, use an ordinary `if not
   __debug__: raise` check. For every spawned Python child remove the
   inherited `PYTHONOPTIMIZE` setting and record the entry's effective
   optimization mode. Do not add `-O` or `-OO` internally.
2. Preserve the old driver too: the B13/frontier worker's FIXED ledger
   pins its exact f6a0b39f... bytes. Provide a small successor manuscript
   driver which uses the same declared selection but routes its actual
   children through the new guard and splits the Lean jobs described
   below. The direct-check route still uses the unchanged bounded research
   runner. Documented direct commands must use the new entry too; changing
   only one catalogue route leaves the other public route vulnerable.
3. For a single Lean compilation, invoke the unchanged research runner
   around `run_lean_batches_signal_safe_20261018.py MODULE`, without
   `--parent-lock` because the inner runner owns the shared lock. This
   retains the59-second/1GiB outer RSS watchdog plus the inner2048MiB,
   one-thread setting and all original admission/abort floors.
4. Use a separate invocation and preserved stdout for the actual named
   axiom audit after a matching compilation. The audit entry should
   require matching prior compilation rather than unexpectedly compiling
   missing dependencies. Do not rerun successful unchanged modules merely
   to report the new entry's existence.

This is a proposed contract, not an implemented or executed repair. Root
owns implementation, its separate source review and any admitted tests.
After correction, the manuscript should state the Python non-optimization
premise or point to the entry which enforces it. Printed Lean axioms still
do not list ordinary theorem parameters; all geometric and finite-data
premises must remain visible in the logic maps.

## 9. Scope retained for the next audit batch

The local universal projection/cover proof, native-support generators,
all remaining finite context domination and simultaneous-rule geometry,
complete dimension59..163 case algorithms, and the unbounded geometric
tail have not received full fresh source audits from this author. Later
dimension57 financing research toward56 is outside the frozen theorem57
and outside this report. Neither a new bound nor a global mathematical
PASS follows from this partial audit.

The coherent low-eight mathematical subdomain has no identified blocker;
the public-entry defects require a reproducibility correction. These are
different verdicts and should remain separate in reader documentation.
