# A half-sum parameter for the A15/D15 repair

Write x for the marked weight at whole D15 pairs(1,13),(1,14) and y
for the weight at whole A15 pairs(0,13),(1,14). Every choice

    x>=0, y>=0, x+y=1/2

preserves triangle lower2 and quadrilateral lower1 in the current
complete face rule. The previous rule chose x=y=1/4. The executable
evaluates x=1/2,y=0, giving actual A15 cost100 and D15 cost101. This
endpoint is an accounting diagnostic, not the adopted search rule.
Both marked
orbits are intrinsic under the diagram automorphisms.

All other same-component weights, ideal weights, ordinary cross
payments and contextual rules are unchanged. In particular the guarded
A1/D5 rule and rank13 parameter9/64 remain present, A16(1,14) stays
zero, and the unused small/D15 center contribution stays zero. This
is a local face theorem. A complete changed finite count is required
before asserting a new global maximum or dimension consequence.

The [producer](../scripts/fresh_rank15_repair_parameter_retention_20261018.py)
and [certificate](fresh_rank15_repair_parameter_retention_20261018.json)
start from the [current witness-context rule](fresh_a1_d5_witness_context_retention_20261018.md)
and bind the [quarter repair theorem](fresh_a15_d15_quarter_a16_zero_center_retention_20261018.md).

## The only new obligations depend on x+y

Both A15 and D15 orbits were zero in the complete baseline preceding
the D14 selected exchange. Consequently the older complete proof,
its later finite and ideal replacements and its independent short
paths never required a positive contribution at either orbit. Adding
arbitrary nonnegative values to them preserves every such output.
The small/D15 center-zero argument also spent neither orbit.

The D14 exchange created exactly four new branches that needed the
restored A16(1,14) half. The complete quarter theorem replaced that
half at four distinct corners in each branch, as follows:

    D14 long-leaf half, D15 repair x, A15 repair y, A16 zero.

Their total is1/2+x+y=1. These are contexts673 and1558 without an
extra contact,1646 with contact13, and1933 with contact3. The whole
24-variant selected proof, including the proper-completion exclusions,
radical sealing and saturated extra-root branches, is unchanged.
It establishes the whole factors under arbitrary supporting whites
before any of these fractional weights is spent.

The producer aligns every corner with the inherited canonical role
extractor and checks all24 variants at(x,y)=(0,1/2),(1/2,0) and
(1/4,1/4). Every corner is nonnegative and every total is at least1.
Each payment and total is affine on the segment, so its endpoints
certify the entire family. The four required sums are identically1;
the other20 variants retain the lower bound throughout the segment.
The new endpoint therefore does not infer validity merely from the
four favorable saved examples.

## Complete current output retention

Replay the complete older proof with its original zero values for
both rank15 orbits. Its short-path, triangle, low/disjoint, protected,
adaptive, universal ideal and later finite outputs remain available.
The D14 selected theorem's complete current output audit is inherited;
where its lost A16 half is reached, use the24-variant half-sum proof
above. This is the same full routing used by the quarter theorem,
with its original arbitrary-support completeness inputs.

The subsequent rank13 parameter theorem changes only specified
cross-component roles and its guarded A2/A2 charge. Its four sealed
repairs spend no A15 repair-orbit contribution. The selected D15
half there is(1,12), distinct from the restored(1,13)/(1,14) orbit.
Its positive chosen contributions therefore survive. The producer
checks the marked-role exclusion in those seals. The current A1/D5
witness theorem changes only a different cross-component contribution;
its original A11 and later A9/A2 seals spend neither rank15 orbit.
Its universal chosen-output replay consequently remains compatible.

The recent rank15 quarters were never a new triangle obligation:
their marked distances are long, and the complete pre-restoration
triangle proof already supplied lower2 with both orbits zero.
Thus all triangles and quadrilaterals retain their stated lower bounds.

## Exact component costs and the necessary new count

Each orbit contains two distinct marked pairs, so the family costs are

    A15:100+2y=101−2x,
    D15:100+2x.

B15 remains97. The representative rank15 allowance must be
max(101−2x,100+2x); at the chosen endpoint it rises100.5 to101.
Actual A15 decreases by1/2 while actual D15 increases by1/2. Thus
the preceding global maximum cannot be transported by monotonicity.

The source reconstructs both complete marked tables at the chosen
endpoint, verifies their exact sums, and checks all actual finite
base costs at both endpoints and the old midpoint. In scale64,
actual A15 loses32 and actual D15 gains32 from the current rule.
Every cross/context function is the identical current function.

The old formal witness2A16+D26 still costs380.5, but this is only one
unchanged forest. More decisively, the two actual formal forests
2A15+2D14 and2D15+2D14 have costs381−4x and379+4x. Their larger cost is

    380+4*abs(x−1/4).

Thus the equal quarters uniquely minimize this lower obstruction at380.
At the tested endpoint the D15 forest costs381; A13+3D15 also costs381.
Lowering the A15 cost alone would overlook actual D15 assignments in
the same complete type fibers. No full dynamic program is needed to
reject the endpoint as an improvement. The search retains equal quarters
unless another geometric change removes this obstruction.

No complete new maximum, new obstruction count, dimension57 proof or
constructed polytope is asserted here. The completed unrestricted bound
remains58.

## Verification scope

Generation99877/chunk07b347 and separate fresh replay13977/chunk84d5d2
both terminated exit0/PASS, with logs46cc96 and4b6760 reporting the
same local rule. Two peers completed full source/report mathematical
reads of the original-zero replay, whole24-variant endpoint argument,
current contextual compatibility and exact family costs; both passed
without correction. Neither executed the producer or repeated native,
supporting-root, census, forest or ideal calculations. The explicit
minimax diagnostic above was added after those local-proof reads and
checked directly from the exact family cost functions.

This executable binds the frozen inputs,
checks all24 endpoint variants, complete changed marked tables, exact
base formulas and the stated positive rank13 interface. It runs no
historical native-map search, supporting-root census, forest dynamic
program or ideal convolution. Those geometric completeness and
sealing theorems remain explicit inherited inputs.

Source, certificate and report are final/frozen after these scoped
checks. The endpoint diagnostic is retained to explain why it is not
adopted in the next complete count.

No external mathematician refereeing or proof-assistant formalization
is asserted.

    python3 scripts/fresh_rank15_repair_parameter_retention_20261018.py --verify
