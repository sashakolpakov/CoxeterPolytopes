# Independent mathematical review of the D17 actual-I transcription

25 October 2026. Reviewer: `audit_finite57_algorithms`.

**Assessment.** I read the complete 464-line root-authored D17 fragment,
both complete primary proofs (536 and 456 lines), the complete 198-line
ordinary-D17 fragment, and the shared local B17 proof interfaces. I found
one wording/scope defect, F1 below. The missing branch already has the
strict direct ratio proved in the primary source; it needs no new
mathematics and does not change the elimination. Subject to that repair
and the explicitly stated common geometric inputs, I found no remaining
mathematical defect in this D17 translation. This is source review, not
execution, compilation, Lean verification, or a new global bound.

I authored the new B17 local-kernel transcription and contributed earlier
component/algorithm audits. I therefore do **not** count this report as
different-author review of those local lemmas. Root has separately read
and reviewed that complete B17 body. The D17 fragment, its fresh weak-five
composition, centre inventory and D17 application of the local lemmas
are root-authored and are the different-author subject of this report.
Historical PASS statements and older acceptance headers were not used
as substitutes for the arguments inspected here.

## Exact source identity and reading scope

All paths are relative to the repository. These are complete file reads,
except that only the named common interfaces of the large main paper
were revisited; their complete earlier reviews remain separate evidence.

| File | Lines / bytes | SHA-256 |
|---|---:|---|
| `manuscripts/revisions/parts/compact_exact18_d17_actual_i_working.tex` | 464 / 21953 | `b777eebcc8c912b9829d30eaf8f52057ad4b42c39c9b1ca5c2e4f474bd963172` |
| `results/compact18_d17_actual_i_fresh_tails_strong_five_and_weak_centre_reduction_candidate_20261023.md` | 536 / 24082 | `bf0a988b7a4ecbbe5616e13b5a072e6ed7150a0ce0a0b11899a6cc0551e5ea8d` |
| `results/compact18_d17_actual_i_complete_fresh_centre_and_head_exclusion_candidate_20261023.md` | 456 / 22285 | `1ffd62104db23a75a2a4445cee0630dc520ede93971855ada8f92b85a99be8e4` |
| `manuscripts/revisions/parts/compact_exact18_rank17_classical_working.tex` | 198 / 9093 | `493304469774194b36d14f06e159f49432619354db03925f906bf104fef6899a` |
| mathematical B17 local-kernel preimage, archived below | 890 / 33559 | `9e304afc509a9ea3180e2fb7dcff77f43dcc8494413e68fd1397784fcdba0c33` |
| `results/compact_classical_generic_prefix_energy_20261020.md` | 212 / 9426 | `f637d7294d17e37581330ca508c914c23a08ce259845614b055e45815d2790ef` |

The D17 freeze is
`results/compact_exact18_d17_actual_i_transcription_freeze_20261025.json`,
SHA `23031071797501476657d68b847a6d0fa67d670812313e70eaeb15d8f94953e6`.
Its 53-line ledger has SHA
`3547a353c14b9158fbe44e210619efb811d9413dbe0cd3b596a7c1f99d2c9e04`.
The complete immutable D17 copy is under
`results/source_archive/compact_exact18_d17_actual_i_transcription_20261025/`.

The shared main interfaces inspected were `app:Ijoint` (actual isolated-I
height), `app:prefixinventory` (actual rank-valid bulk and strong ratio),
`app:exact18layers` (reference payments), and `app:D17ordinary18` (only
its actual D4/D3/head lists and point constructions, without its ordinary
rank-one debit). Their source was main SHA
`4afcf81a02733f0e86d09577a521e567c2e7acd46b315e218aba8ccd7271680d`.
During this review root integrated B17 into the main. The old 495425-byte
main remains byte-exact in
`results/source_archive/compact17_before_b17_actual_i_integration_20261025/compact_farthest_vertex_bounds_consolidated.tex.gz`:
gzip SHA `a84f3ff7c7c76ae308090bab56a24b3be5553129fe4ba5f61242674c2d0746d2`.
I verified its decompressed hash and size by byte I/O only.

The B17 mathematical preimage is byte-exact in
`results/source_archive/compact_exact18_b17_actual_i_transcription_20261025/compact_exact18_b17_actual_i_working.tex.gz`,
gzip SHA `b7208a5f3a2f84a7a00f49f32d955a5456c7bd0797b4b5bd497cdc3ffad0d201`.
The current 891-line B17 successor SHA
`cc9e659c2419757fbfe42d28e42d413dd09be3aa592ab0f19e270118c22108ed`
adds only the missing closing figure colour-group brace. I inspected the
complete one-hunk patch `compact_exact18_b17_i_figure_group_repair_20261025.patch`;
there is no mathematical delta. This review does not claim the current
large main path still has its earlier hash.

## F1: distinguish a simple-three cycle from the end-four direct bound

In the left-only centre paragraph, the sentence beginning “Positive
owns three/four/six and parallel have the actual cycle/star/path/null
dispositions” is too broad without the left-contact qualification.
The second primary proof, lines 158–163, explicitly restricts the
own-three all-three cycle to a **simple-three left contact**. The
transcription must preserve this restriction. A near end-four left
contact with own three has a 3/3/4 triangle, not that affine cycle.

This is a proof-presentation defect with an already available repair.
For an end-four contact at j=1 or 14, Delta=1/15 and eta=1/sqrt(2).
Every positive own has g>=1/2, so the smallest ratio, at j=1, is

    sqrt(15)/2 + 1/sqrt(30) > 2.11 > 1.554.

The near j=14 ratio is larger than four. Thus both end-four own-three
rows are directly strong. No affine-cycle argument is needed for them.
State the simple-three own-three cycle first; retain the actual
far-end-three own-four/six/parallel dispositions; then state the direct
near/next/end-four bounds. Root agreed with this exact finding before
any source modification. Preserve b777 and repair one local hunk.

## The acyclic argument and its actual-wall domain

The theorem concerns the **entire** actual farthest-vertex forest
D17+actual isolated I in compact dimension eighteen. Actual long roots
are n_i=(e_i-e_(i+1))/sqrt(2), i<=16, and the other fork is
n17=(e16+e17)/sqrt(2). Long duals sqrt(2)P_i stop at i=15.
The signed fork duals are (P16 +/- e17)/sqrt(2). The proof never
declares e17 a short normal, sqrt(2)P16 a long D dual, or x17 positive.

There is one original vector z=E(p-x-ye), with ||x||^2+y^2<1.
The actual I opposite supplies y>=s and x1>=t. It supplies no
ordinary y^2>.1 debit. The .923 real barrier is the shared scalar
lemma, not an arithmetic statement. Classical caps retain the original
I wall and every permitted finite foreign contact.

The dependency order is sound:

1. M=11 is valid for actual D17 and gives helper coverage through twelve.
2. The two initial-five branches independently remove the extra D4/D3
   tails, before any strong-five conclusion.
3. Actual prefix nonadjacency extends through fourteen; thirteen or
   fourteen is strong.
4. Retained weak-five/weak-fourteen feet and freshly composed head/I
   payments exclude weak five, without assuming strong centre fifteen.
5. Only then does the flat-five local kernel force fifteen weak.
6. The fresh actual centre inventory removes every nonhelper row.
7. Adjacent-helper exclusion makes fourteen strong, and the conditional
   helper15/head/I kernel excludes the final pure helper.

No B17 short-tail theorem travels across this DAG. No arbitrary face is
assumed Coxeter or arithmetic. All finite classifications are of actual
white endpoint diagrams; all auxiliary points use actual wall normals.

## Fresh tails and the weak-five bootstrap

The reused ordinary-D17 material precedes its energy-dependent steps.
It supplies the complete retained A12+D4 and A13+D3 finite-white lists,
both fork choices, all foreign costs, actual H4 projections and the two
actual A2 cores. The I argument supplies new payments rather than
reusing ordinary y^2>.1 eliminations or their weight 1/7.

For the H4 exceptions, the normalized displacement is flat through
i-3, has norm below i, and demand >3.48 or >3.61. A strong-five
reference plus a first-five/I local gain closes these tails. A weak-five
reference through P1,P6 and the selected actual terminal closes them
without a head gain. The quoted bare lower bounds match the primary
proof and its actual-reference domains.

For end-three finite .9<g<1 at fourteen, the **two** A2 core corrections
retain their cross 8g^2/9. The complete decreasing-domain estimate is
G>3.388; no dotted parameter is admitted to this finite interval.
Weights .16/.15 use actual terminal capacities at ten/eleven. The
strong-five totals .931255 and .9303225 become >1 with .07665.
The weak-five bare totals plus the unused-I graft become
1.00952 and 1.0085. That graft is orthogonal to the tail residual
and its remaining P6 capacity exceeds 1/60.

For weak five, the three rows of P1/P6 capacities come from actual
concave references ending at thirteen or fourteen. The increment sums
and omission losses give precisely q13>C+.1921025 and q14>C+.2103.
The terminal prices >.14144 and >.1361 are not assigned to internal
knots. The finite/dotted fourteen feet are newly evaluated at their
actual time 28g^2+2g-1/2. Their P1/P6 projection coefficients and
residual norms retain their initial components, not merely their tail.

Finite own-five at node five gives x1>.65. All weak-fourteen alternatives
then close at the stated capacities. For the tight finite-at-least-nine
row, .18(.1085351)-.0081=.011436318 and
.9889693845+.011436318>1.0004057. Thus weak five is dotted.
The retained U5 has norm residual <1, gain >.00953 and debits
<.052326/.067883 at P1/P6. Its second coordinate, rather than a
false orthogonality assertion, is retained in the later head cross.

A nonisolated head uses the stronger .643 reference. The smallest
total .9837590045+.00466+.0126789498 exceeds one; the P6 debit
<.076369 fits .07876. An isolated finite head has its **actual**
head/I native relation, h>1/4 and y^2>.0557; q13(.5965) already
excludes it. Hence only pure or positive-foreign dotted heads remain.

For weak-fourteen owns five/seven/eight, the exact U5/14 crosses are
2a/15, a/10, a/10. They have been charged with the displayed weights.
Own eight requires the separate cap **.08**, not .1: the full P6 price
is <.067883+.006365+.013334=.087582<.08806. If clipped,
s<.11, so (s-.2)^2+.00690666>.015. The uncapped branch follows
from the .923 barrier and weighted Cauchy. The graft is orthogonal to
both retained feet because P1 and P6 remain actual paid anchors.
The two surviving fourteen rows have tail-only residuals, no P6 debit,
P1 debit <.070004 and P13 debit <.141422<.14144.

## The fresh nonflat-reference head composition

For a pure dotted head, the joint point is rebuilt at the actual mutual
parameter h, not at an assumed finite head/I pair. Its denominator
D0=u/2-v0 sqrt(h/(1+h)) is positive for u>=.2 and all h<1/2.
For .1<=u<=.2,h<=.2, u(8-3u)>=.77>2/3 gives the same sign.
Its positive actual-cap combination retains both D forks and every
unknown wall. The old products, selected zero cap product and nonpositive
actual cI product are all present in the proof.

I reconstructed the identity relating X and s0, the .91 polynomial
barrier, its shift by .13 and its quadratic lower bound. The cubic and
quartic deductions use |w|<=.13; no sampled minimum is substituted.
The general paid-reference lemma applies with m=6 because the retained
feet have zero P1/P6 products and the first coordinate stays .5965.
The capped kappa<=.12 price is <=.02 at P6, so the total
<.087883 fits .08806. The multiplier tau is nonnegative. No I energy
has already been used. The unclipped gain >.02033425 and clipped
gain >.022281 close the two reference branches. The u<=.1 pure edge
and .1<=u<=.2,h>.2 branches independently cover the excluded domains.

For positive foreign contact, the actual clearing core removes the whole
I coordinate before the head pressure is used. The bound
W0.Z<=.5965-.14759v includes the .05Y5 cross. Its full P1
price <.30713 fits .312168; P6 <.076283 fits .08806. The complete
seven-interval table adds I once and yields a uniform gain >.0214.
Both remaining reference totals exceed one. Thus strong five precedes
the centre-fifteen reduction, rather than being inferred from it.

## Complete centre inventory and every funded branch

Deleting centre fifteen retains actual A14+2A1+I. I reconstructed
qL=2 eta^2 j(15-j)/15, Delta=1-qL-delta^2-epsilon^2-Q,
sigma=g+eta j/15+(delta+epsilon)/2. The fork residuals are independent;
there are no omitted left/fork residual crosses. The sufficient ratio
T(sqrt(30))<1.554 has the actual .61 threshold.

Finite tree/arm/bond classification gives precisely the listed left,
one-fork, two-fork and foreign possibilities. In particular positive-I
single-fork pairs are (3,3),(3,4),(3,5),(4,3),(5,3), and the only
two-fork foreign case is (3,3,I3). Long left branches cannot hide a
small exceptional E/F/H enlargement. No infinite native-label list is
silently replaced by finitely many chosen labels.

The isolated centre uses future clearing cores for 0<u<1, a forbidden
actual null at u=1 and an actual finite projection for u>1. Foreign
u<75/64 and pure finite-at-least-seven u<100/81 give the required
loss bounds; the actual H4 projection handles pure five. Left-only and
full bridges are exhausted by time/foreign-cost, actual affine or direct
ratio arguments, with the F1 qualification required above.

For own zero and a finite fork >=7, the actual core has normalized
signed displacement (1^14,3-2u,u-1,+/-(1-u)). Its full norm is
25-16u+6u^2; its time ratio gives the stated G. I reconstructed all
four interval bounds using endpoint positive squares and the radical
maximum 1/sqrt(6). The least K14 total is >1.0013591. Both fork
signs remain present, and terminal .14/.13 payments fit .14144/.1361.

For fork five the correction uses the **other actual fork**. Its
spacelike corrected root is not called a point: positive finite
projection supplies the point. The pure displacement is P14, not P15.
The genuinely allowed Q=1/4 branch retains the I coordinate in this
same point. The tail energy includes its mixed norm d^2<.765.
Dotted/foreign pure head gains >.04 suffice without an additional I
charge. A pure finite head instead uses its pure gain >.0186 and
only the improvement (s-mu d)^2, with s>mu d. This replaces the
existing I coefficient. The two least totals are >1.00754385 and
>1.0056781625. The head pure residual is orthogonal to the tail.

For positive-own single-fork cases, the own-three/fork-four-or-five
Q=0 endpoint is an actual finite D17 replacement. I reconstructed
C=15/2+(17/2)b and the identity B0^2-A^2=2/17. Its signed fork
direction has full norm17 and demand >sqrt(102/7)>3.817; it does
not require a short D normal. Selected, released-fork and unknown-wall
signs justify an actual edge endpoint, not an unproved untruncated edge.

The own-five/fork-three H4 point or foreign-cleared future core is
flat through twelve, has demand >3.52 and norm <15. The funded
eleven/twelve totals are >1.0240725 and >1.0229116. The remaining
whole .9<g<1.2,Q=0 domain uses the explicitly displayed core, with
square -k, zero cap/fourteen products and negative fifteen/fork products.
The derivative sign for k/alpha^2 is negative for
2g^3-3g-3/8, hence the stated maximum bound at 1.2 is valid.
The direction is flat through **thirteen only**. At fourteen its
projection (1-a/14)P14 and residual norm <1+2a are used instead.
The least fourteen total >1.0051616 uses a legal terminal debit, not
the internal P13 coefficient.

The two-fork rows retain every finite pair and allowed foreign contact.
Pair-three-four gives an actual B17 reciprocal endpoint with the
displayed exact C. Pair-three-five uses actual finite H3 projection;
q=2 phi^-6<1/8 gives the full continuous pressure >3.53. Their
least funded totals exceed 1.0648 and 1.0240 respectively. The
triangle, null-square and direct-ratio dispositions cover the remaining
own and foreign rows.

## Final helper and shared-kernel boundary

The only weak centre left is the pure far-three helper. The actual
two-helper Schur calculation uses inverse entries 2min(i,j) through
fifteen and retains both positive fork blocks. It makes fourteen strong.
The resulting q14>.90715 reference has first-five mean f5/5<.361,
remaining P5 coefficient >.146 and terminal P14 coefficient >.1361.
These are the complete shared-reference hypotheses; mere mean/energy
bounds on an arbitrary unpaid vector would not suffice.

Every imported helper pressure uses only actual n15,n1,n2 and the
actual long dual sqrt(2)P15. Both forks have nonpositive products
with each correction. P16 occurring in Euclidean coordinates is only
a spatial sum. The time, initial residual, adverse head cross and
individual prices are identical to the conditional local B17 packet;
no B17 short-normal or terminal conclusion is invoked.

Finite/nonisolated heads are excluded before the final dotted coupling.
For pure dotted u<=.2 the standalone >.11088 gain bypasses the tail.
For u>=.2 the retained helper gain >.0421 and joint head/I gain
>.0515 are one square, including the full initial cross. The P5
debit <.09873 and P14 debit <.1273 fit their individual capacities.
The foreign head branch has its adverse cross and single I energy
included in >.0521, with P5 debit <.09273. The least final total
is .90715+.0421+.0515=1.00075. Tail-elimination contradictions
are not added to this final reference.

## Status and exact next action

Preserve b777 and both complete primaries. Make the one-hunk F1
clarification and inspect its actual delta before integration. The
proof inputs remain compact acute simplicity, finite/affine Coxeter
classification, native products for distinct actual walls and the
common original-farthest comparison/positive-projection machinery.
They are not replaced by arithmeticity, field restrictions, script
enumeration or arbitrary-face inheritance.

This report supports the D17+actual-I **row theorem**, conditional on
those named common foundations and the independently reviewed shared
local kernels. It does not establish the remaining whole B18/D18/A18
rows, a new general compact dimension bound, a complete paper build or
Lean verification. No scientific worker/import, numerical script,
enumeration, syntax check, test, Lean/TeX/Sphinx build or job was run.
Only source reads, manuscript/report writing and metadata/archival
byte verification were performed by this reviewer.
