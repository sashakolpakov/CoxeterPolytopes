# The next finite and ideal obstructions at ambient63

Under the completed dimension63 payment rule, the exact finite maximum
at rank63 is **453**, exceeding the curvature899/2 by7/2. The complete
ideal comparison at link rank62 has **2,204 profiles,21 negative margins
and no zero margins**. Thus the current payment inequalities do not
exclude dimension63. The completed unrestricted bound remains63.

The [source](../scripts/fresh_ambient63_payment_frontier_20261018.py) and
[certificate](fresh_ambient63_payment_frontier_20261018.json) compute these
obstructions and a few explicitly conditional targets. They bind the
[completed proof](fresh_order11_order17_dimension63_20261018.md), its
[enlarged finite coupon budget](fresh_e6_b10_half_coupon_forest_budget_20261018.md),
and the [completed ideal envelopes](fresh_ambient64_ideal_margin_closure_20261018.md).
No local weight reduction, new dimension bound or example is proved here.

## Exact current finite maximum

All current finite payments are retained, including the depth11 zeros,
D16=116 and D17=127, the whole E6/D11 half coupon, and the whole A1/B10
doubled-terminal half. Three representative forests attain453:

|Forest|Cost|
|---|---:|
|D15+3D16|105+3*116=453|
|2D15+D16+D17|2*105+116+127=453|
|2D17+D29|2*127+199=453|

These are formal forest attainers, not assertions that corresponding
vertices occur in a real chamber. The list is the set of canonical
attaining witnesses returned by the complete component-vector comparison;
it is not a classification of every typed maximizing forest. For example,
an A16 has the same same-component cost116 as D16, so reducing D16 alone
would leave analogous A16 obstructions.

The calculation covers all finite types. It keeps E6 at rank6 because of
its coupon, and A17 separately because of its small-component interactions.
Every other type is dominated by A_r for r<=9 and D_r for r>=10. The source
checks same-component and every type/partner cross-cost dominance under
each tested scenario, including arbitrary dihedral labels through the
inherited uniform rank2 representative. There are192 types and36,864
type/partner comparisons per scenario.

All D10 through D16 and A17 multiplicities of total rank at most63 form
700 complete vectors, giving328 distinct coupon environments. The remaining
components have ranks1 through9 or at least17. The exact interaction key
and literal coupon formula are the same as in the frozen enlarged budget;
in particular E6/D11 has doubled cross cost4, and A1/B10 has doubled cost1.
The latter remains dominated once by D10, without double-counting B10.

For each environment, rank-first dynamic programming computes every
residual value through rank63. The current tables agree entry-for-entry
with the corresponding prefixes of the frozen item-first rank64 tables.
An independent saved-table verification then checks every possible
last-component inequality and an attaining predecessor for every rank.
Each global vector reconstructs a typed attaining forest whose total
same-component and cross-component cost is recounted directly.

Per scenario the certificate contains544,152 Bellman upper inequalities,
20,664 attaining predecessors,39,200 literal interaction comparisons and
700 global alternatives. Induction on residual rank proves both the upper
bound for every forest and formal attainment; no sampled forest maximum
is being extrapolated.

## Finite targets worth distinguishing

Every scenario below preserves all current coupons and every unlisted
same-component cost. A rank16 reduction is applied to both A16 and D16;
B16 remains104. A rank15 reduction applies to A15,B15,D15 together.

|Conditional costs|Exact rank63 maximum|Slack below899/2|
|---|---:|---:|
|Current rule|453|−7/2|
|A16/D16=114 only|453|−7/2|
|D17=125 only|453|−7/2|
|A/B/D15=104, A16/D16=115, D17=125|449|1/2|
|A16/D16=114, D17=125|449|1/2|

The last row is the narrowest tested family target: two more units saved
on both A16 and D16, and two on D17, would close the finite comparison
with the enlarged coupons unchanged. Geometric feasibility is unproved.
In particular the present D17 cost127 came from a deliberate restoration
of two units, so deleting them again cannot be presumed compatible with
the D16 saving.

The three current witnesses impose useful necessary inequalities. Write
delta15,delta16,delta17,delta29 for same-component savings from the current
costs, with delta16 also valid for the A16 substitute. To reach899/2 one
must have

    delta15+3delta16 >=7/2,
    2delta15+delta16+delta17 >=7/2,
    2delta17+delta29 >=7/2.

These witness inequalities are necessary, not sufficient. The complete
DP proves sufficiency only for the explicit scenarios in the table.
The D29 witness shows why work confined to ranks15/16 cannot settle the
finite side. Conversely, changing D17 alone leaves D15+3D16 intact.

## Complete ideal rank62 domain

The current normalized costs, with r the affine rank rather than factor
order, are

    c(r)=r/2 for1<=r<=9;
    c(10)=109/22, c(11)=16/3, c(12)=74/13, c(13)=6,
    c(14)=94/15, c(15)=107/16, c(16)=118/17;
    c(r)=163/23 forr>=17.

The inherited all-rank convex-order argument concentrates two ranks below9
until one vanishes or the other reaches9; two ranks at least18 concentrate
until the smaller reaches17. Both moves preserve this cost sum and
decrease normalized curvature. Thus all actual rank62 cusp profiles reduce
to arbitrary multiplicities of9 through17, at most one rank1 through8,
and at most one rank at least18. The2,204 profiles are the complete domain
at rank62, not the2,419 profiles appropriate to the previous rank63 link.

For each profile the source computes the entire product polynomial
H=product(1+t+...+t^r), its value V=H(1), and

    G=sum_j [t^j]H*min(j,62−j),
    margin=(5G−62V)/(2V)−sum c(r).

Direct convolution is used; an exact cumulative-distribution identity
checks G independently from the coefficient vector. All rows and all21
negative rows are saved. The largest normalized deficits begin as follows.

|Affine ranks|Margin|
|---|---:|
|(10,10,12,15,15)|−22543/201344|
|(10,11,11,15,15)|−5653/50688|
|(10,10,11,15,16)|−18281/197472|
|(10,10,10,16,16)|−27785/384659|
|(11,12,12,12,15)|−863/13182|
|(9,11,12,15,15)|−2237/39936|

Two further profiles show that improving only factor orders16 and17 can
never finish the ideal comparison:

    (10,10,10,10,10,12): −8607/190333,
    (10,10,10,10,11,11): −4021/95832.

They contain only factors of orders11,12,13. Changes in larger factors
cannot change their margins. A one-unit drop in the order16 envelope
alone still leaves10 negative profiles.

## Conditional ideal targets

The following three all-family total-cost targets make every one of the
2,204 margins positive, with minimum3505/201344 at(10,10,12,15,15):

|Factor order|Current maximum|Conditional maximum|Present bottleneck|
|---|---:|---:|---|
|11|109/2|54|C has54.5; B already54|
|12|64|63|B has64; C already63|
|13|74|147/2|D has74; B and C already73|

All other family maxima remain unchanged. These are total pair-cost
drops of1/2,1,1/2, not normalized cost drops. The normalized gains are
divided by factor orders11,12,13 respectively.

With the first two drops fixed, the exact additional order13 drop needed
for nonnegative margins is4239/15488, attained at(10,10,12,15,15).
A strict increase above this threshold makes every margin positive;
the proposed half-unit drop is sufficient. There is no unaffected negative
profile hidden by that one-parameter calculation.

An alternative two-order target lowers the order11 maximum by3/4 and
the order12 maximum by1. It gives minimum157/8788 at(11,12,12,12,15).
With the order12 drop fixed, the exact order11 threshold is22543/36608.
This alternative requires reducing B as well as C at order11, since its
target maximum is215/4=53.75. The three-order target may therefore be
more convenient geometrically, although no such transfer is supplied here.

All tested changes leave the low linear range through rank9 and high
constant range from rank17 unchanged. The same2,204-profile domain remains
complete for every conditional scenario; no old concentration restriction
is silently reused after changing a boundary cost.

If a complete compatible local theorem achieved a finite maximum at most
899/2 and positive ideal margins, the already proved global constant
5*floor(63²/4)=4960 would finish the ambient63 comparison. Both new local
obligations remain open. The diagnostic does not promote that conditional
observation to a dimension62 result.

## Verification scope

The first generationc98670 passed all current and conditional exact checks.
Final source-bound generation02d48e and separate fresh saved replay35f38f
both terminated exit0/PASS. Direct audit6c2cc3 passed all12 source bindings,
all5 report links and source syntax. The source imports no
project arithmetic modules. Its current finite tables are independently
reconstructed and compared with the frozen enlarged-coupon tables; each
conditional table has its own complete upper/attainment certificate.
Geometric classification, the existing payment theorem, and all-rank
concentration remain inherited written inputs. No new coupon, pair deletion,
half transfer, or global dimension theorem is proved by this computation.

A separate agent completed a full source/report mathematical and scope
review, also reading the frozen finite-budget formulas. It passed the
complete threshold9/17 domain, current-versus-conditional costs, rank-first
Bellman upper and attainment proofs, typed coupon domination, simultaneous
A16/D16 targets and exact threshold interpretation. A light independent
standard-library read of the saved summaries confirmed the stated scenario
results and absence of unaffected zero rows in the two successful threshold
claims. That reviewer did not rerun the profile arithmetic or DP optimization.
No concrete issue was identified. Source, data and report are final/frozen;
the completed dimension63 proof and its separate aggregate are unchanged.

    python3 scripts/fresh_ambient63_payment_frontier_20261018.py
    python3 scripts/fresh_ambient63_payment_frontier_20261018.py --verify
