# Full outside source review: conditional quadratic D16 two-fork census

23 October2026. **FULL different-author STATIC CODE/DOMAIN/CONTRACT PASS,
no correction requested.** This is a pre-execution source review. It is
not a worker run, syntax/AST/import check, independent second enumeration,
survivor count, geometric acceptance or whole-D16 exclusion.

I read ALL530 worker lines and ALL229 contract lines. I independently
checked the finite conditional domain, full integer rings, exact signs,
complete field squares, original raw pair identities, all ordered triples,
canonical output/replay and the complete unchanged launcher/guard closure.
All six worker PINS and both target objects freshly match below.

Composition authored this worker and contract; I authored neither. Root
supplied its prospective raw domain. I contributed to the separate D16
geometric reduction and earlier B16 mathematics, so this report is NOT
an independent review of those geometric theorems. My source review is
independent of the new implementation's author. The actual geometric
reduction remains an EXTERNAL premise with its own full outside/root gate;
this report neither assumes nor supplies that gate.

No scientific worker, repository module, parser, AST, help/plan, test,
optimizer, enumeration, compiler, Lean, TeX or Git operation was executed.
Only source reads, ordinary metadata/hash I/O and writing this report
were used. No target or frozen input was edited.

## 1. Exact immutable scope

|Object|Bytes / lines|SHA256|
|---|---:|---|
|scripts/check_arithmetic16_d16_two_isolated_dotted_forks_quadratic_20261023.py|21704 / 530|1fe11f2a6bd6daeaeb3a4f83249a756bbaa359bd30b05ed7f008a9d021f668da|
|results/arithmetic16_d16_two_isolated_dotted_forks_exact_quadratic_census_contract_20261023.md|12137 / 229|596d1bfce60ae0ad2f82fe0c154f373dd1898eca359b2d03861dcf62ac5b4b70|
|scripts/check_arithmetic16_b16_isolated_quadratic_head_short_pairs_20261023.py|19961 / 506|452587a42fcc03219cc77dabfa1ff0f691c81617c776943276b39e822d8ff42f|
|scripts/run_bounded_research_check_20261018.py|5703 / 129|4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132|
|scripts/build_illustrated_ra_manuscript_20261018.py|6606 / 149|27f14501b52a38c92121266aceb58eb2e1d36bfc21c2615ab0b4e9be4f9dee15|
|scripts/run_lean_batches_signal_safe_20261018.py|1234 / 48|c0ddc98cd98e55f669252e4f138345665c2ce2b0e4adf12d8d4ce96aed65af4d|
|scripts/run_lean_batches_20261018.py|14192 / 350|97822ecf4d9aef2bf04aa429b582c4cc36cf6d1b6f0364579251278f66fae72c|

The six PINS are exactly the contract, read-only B16 helper provenance
and four guard/launcher sources. The old helper is NEVER imported; its
old geometry, fourteen old pins and previous executions are not semantic
inputs of this new conditional search. I compared its local exact helper
patterns with the new self-contained implementation. The new worker
contains no repository import or child process. The four guard sources
were read in full, including the functions reachable through their imports.

The contract avoids a circular worker self-hash. An external freeze binds
the worker, and every actual invocation must supply its expected SHA256.
The worker authenticates that hash and all six immutable inputs before
work and before completion; the runner also hashes the checker and itself.
No unaccepted geometric source is falsely pinned as an accepted premise.

## 2. Complete conditional field and maximal-order domain

The explicit EXTERNAL premises are ACTUAL quadratic K=Q(sqrtD), all
three parameters generating K, and strict strips

    4<H<98/9, 0<H'<2;
    4<U,V<121/2, 0<U',V'<1/2.

I rederived the finite boxes from these exact inequalities. Writing
x=(a+b sqrtD)/2 gives a=x+x' integral. Thus head a5..12 and
fork a5..60 cover every strict possibility. Physical-minus-conjugate is
b sqrtD>0, so b>0; omitting b0 follows from the generator premise,
not a numerical heuristic. The gap bounds give integer
b²D<=118 for heads and b²D<=3660 for forks.

The full maximal-order parity is correct: for squarefree D=1mod4,
a,b have the same parity; otherwise both are even. In particular the
half-integral fields are not replaced by Z[sqrtD]. The head bound gives
D<=118 in the first case and4D<=118 in the second. Thus `field_list`
visits exactly the necessary squarefree radicands in2..118 with
D=1mod4 or D<=29, sorted by fundamental discriminant and D.

The three fundamental-discriminant intervals [0,40),[40,80),[80,119)
are disjoint and cover every allowed field. Both boundary values40/80
have one owner. Finite trial-square-divisor testing through isqrt(D)
is sufficient for squarefreeness. It neither invokes nor relies on an
algebraic-number factoring program or a field catalogue.

`models` iterates the entire safe coefficient boxes, applies exactly
the full parity law and all four OPEN embedding inequalities using
exact signs, and retains the complete result in(a,b) order. The
physical raw times are positive and their conjugates negative by the
stated strips; their checked identities are correct. Empty model lists
remain explicit field records. No empirical model or survivor cap is
silently used to prune the mathematical domain.

## 3. Exact arithmetic, signs and full K-square completeness

I checked every operation in Q. Coefficients are Fractions, fields
cannot be mixed, multiplication is exact and division uses the nonzero
rational norm. The sign routine handles zero and equal coefficient
signs directly. For opposite signs its sign(a)*sign(a²-Db²) is the
correct physical sign. Norm zero is impossible for a nonzero element
of the visited nonsquare fields. All equality comparisons are exact.

`integral` implements the full O_K: p,q integral in the ordinary
basis, and2q plus p-q integral for D=1mod4. This is equivalent to
the full half-basis parity condition. Encodings preserve numerator and
positive denominator for both coefficients, with the field identifier.

`rational_sqrt` checks numerator and denominator separately by isqrt
and exact squaring, including zero and the negative no-root case.
For a rational field-square radicand, `field_sqrt` retains BOTH the
rational and pure-radical possibilities. For mixed p+q sqrtD, q!=0,
the norm of any root must be one of BOTH signs +/-sqrt(NR).
Testing a²=(p+n)/2 and b=q/(2a), then resquaring EXACTLY, is
complete. The discarded a0 case cannot solve2ab=q!=0. The chosen
root is the physically positive one, with no conjugate-sign restriction.
No rational norm-square is mistakenly accepted as a field square.

The complete quadratic native L=4lambda² alphabet is also correct:
0,1,2,3 and the native5/10 golden, native8 silver and native12
sqrt3 values. Degree phi(m)/2<=2 exhausts finite m3/4/5/6/8/10/12;
orthogonal native2 is separate. Every L>4 remains dotted. Compact
parallel L4 is rejected; negative physical L cannot be accepted.
No finite-ground-field assumption is borrowed from the old B16 worker.

## 4. Actual raw pair formulas, all three quarter-products and triangle

The worker uses EXACTLY the conditional raw times and spatial products

    T_H=H(H-2), T_U=2U(2U-1), T_V=2V(2V-1),
    S_HU=HU/2, S_HV=HV/2, S_UV=7UV/2.

These match the required original root squares2H,2U,2V. The physical
future time product is the positive square root of T_i*T_j, so
g=root-S_ij is the correct negative raw Gram product. This is a fresh
D16 normalization, not the old even-B short normalization.

Every head/fork pair and EVERY ordered fork/fork pair tests, in the
declared order, product/4 integrality, full K-square time product,
physical g>=0, original path-cycle g integrality, two-cycle
L=g²/(U_i U_j) integrality,0<=L'<4, and the complete physical
native-or-DOT alphabet. g0 is retained with L0. Nonzero g cannot
have a zero algebraic conjugate, but permitting the nonnegative
lower endpoint does not discard a legitimate zero mutual product.

I independently expanded both internal invariants:

    L_HU+g_HU=(15HU-32U-8H+16)/4,
    L_UV+7g_UV=15UV/4-8(U+V)+4.

The left/right parameter positions in the implementation are correct.
These are exact consistency assertions on already computed raw data,
not additional unproved geometric filters. The head/HV equation uses
the same actual normalization. The three quarter-products are ALL
tested through the two head/fork and one fork/fork pair results.

For every ordered triple the separate exact original triangle
g_HU*g_UV*g_HV/(HUV) is tested for integrality. It is not inferred
from pair squares or assumed individual Cartan entries. Its conjugate
sign is unrestricted. Equal numerical fork parameters are retained;
the fork pair still uses the two GENUINE opposite signed directions,
with spatial product7UV/2, not a self-root scalar product.

`field_records` enumerates ALL heads and all ordered forks, including
every equal pair. Its first-outcome dispatch H/U,H/V,U/V,triangle is
complete and matches the contract. Pair and triple totals are checked
against the full Cartesian-product cardinalities. Every survivor has
all three exact pair certificates and its triangle value. There is no
survivor-count-based stopping rule or unproved symmetry quotient.

## 5. Halting, memory, records and authenticated separate replay

All loops have the proved finite coefficient bounds. A head conjugate
strip has at most4 integer a values per b, and b<=7, giving at most28
heads. A fork strip has at most1 integer a per b and b<=42, giving
at most42 forks. Consequently 28*42+42²=2940 cache entries and
28*42²=49392 triples per field are SAFE upper bounds. They are not
claimed executed counts. Every Fraction/isqrt operation terminates.

Only one field's models, pair caches, counters and one output record
are held. Survivor records stream immediately. The2MiB record and
64MiB file limits cause FAILURE, never truncation or a false terminal
completion. No total-success-time or formal heap allocation cap is
promised; actual guards remain independently mandatory.

Canonical ASCII JSON uses exact integers, sorted keys, fixed separators
and a newline. The field_begin includes all exact models/times; the
field_end includes every-pair/every-triple digests, exact first-rejection
counters and totals. Digests include successful certificates and every
rejected outcome in the deterministic enumeration order. The output
does not falsely claim to store every rejected row literally.

Outputs are NEW flat results/ files. Final symlinks are rejected,
resolved parents must be results/, and exclusive xb creation prevents
overwriting an existing result. Partial files remain possible failures;
their existence is not completion credit. Only a successful terminal
record PLUS actual guarded EXIT0 evidence could support a later result.

The SEPARATE combine requires exactly three distinct part paths and
three hashes in index order. It regenerates each complete part with
the same implementation, reads bounded lines, compares EVERY canonical
byte record and requires strict EOF. It rehashes each part after replay,
checks the disjoint whole-field total, and streams a NEW combined
certificate. Duplicate/missing/reordered/changed/oversized/partial parts
or extra trailing bytes fail. No JSON coercion or trusted saved counter
can replace this full regeneration.

The correct statuses are ENUMERATION_COMPLETE and
SAME_IMPLEMENTATION_REPLAY_COMPLETE. The flags explicitly retain
independent_implementation=false, geometry_conclusion=false and no
realization claim. This is authenticated same-implementation replay,
not a second independently authored algorithm. I have not run either.

## 6. Full guard closure and exact scope of this PASS

I read the unchanged129-line bounded runner,149-line monitor,48-line
signal-safe wrapper and350-line base source. Their imported main/build/
Lean functions are protected by main guards; using the runner's monitoring
helpers does not initiate TeX or Lean. The new worker imports none of them
and acquires no lock itself, so the contract's --parent-lock is necessary.

The shown command correctly supplies the SCRIPT directly after the
runner's --, with no second python. The wrapper obtains the shared
formal/.proof-batches/exclusive.lock nonblockingly before spawning,
retains/inherits it, checks watchdog access before child launch, sets
numerical thread variables1 and monitors the whole owned descendant
tree. Its time59s and sampled RSS1GiB ceilings,40GiB/6GiB admission,
8GiB/3GiB live abort floors and owner-loss cleanup are unchanged.
The child bootstrap retains its owner alarm/watch and owned process group.
Signal cleanup targets that group, not unrelated processes. A monitoring
access refusal fails before a scientific child; no bypass is introduced.

Actual launch authority, current resource admission, separate unit stdout/
ledger preservation, four distinct actual jobs and later evidence review
belong to root. This PASS does not authorize a lower RAM gate, a raw
worker launch, parallel parts, an imported test or a successful-module rerun.

The external geometric theorem remains separately authored/reviewed.
The historical EXTERNAL_UNACCEPTED_AT_PREPARATION flag is honest and
need not be rewritten if a later root conjunction authenticates an actual
accepted geometric source. No survivor count, empty census, chamber,
whole-D16 deletion, arithmetic15 or new Lean instance follows from this
static report. Subject to those explicit boundaries, the ENTIRE prepared
worker and contract have no remaining defect found in this full review.
