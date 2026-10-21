# Saved alternatives to the A16 half restoration

All four branches that require the newly restored A16 half have the
same canonical alternatives. Restoring one quarter to each of two
D15 entries and two A15 entries, while returning the A16 singleton to
zero, passes all24 saved whole variants. Among443 saved actual-fiber
witnesses, this changes the largest value381.5 to6093/16=380.8125.
This is a **saved-witness diagnostic**, not a new face theorem, a
reoptimized forest maximum, or a dimension57 bound.

The [producer](../scripts/fresh_a16_restoration_alternative_cost_diagnostic_20261018.py)
and [data](fresh_a16_restoration_alternative_cost_diagnostic_20261018.json)
read the frozen [selected D14/A16 exchange](fresh_d14_long_leaf_a16_half_selected_transfer_20261018.md)
and [current complete count](fresh_d14_small_d15_ambient58_frontier_20261018.md).
They also bind the full current A15/D15 weighted tables and the
[small/D15 center theorem](fresh_small_d15_center_three_eighth_retention_20261018.md).
No matrix search, root-completion calculation or dynamic program is run.

## Canonical roles in the four branches

The affected branches are673 with no extra white,1558 with no extra
white,1646 with extra contact13, and1933 with extra contact3. Each has
four distinct corners with these same canonical roles:

|Whole marked factor|Current contribution|Proposed variable|
|---|---:|---:|
|D14(0,12)/(0,13)|1/2|1/2, unchanged|
|D15(1,13)/(1,14)|0|x|
|A15(0,13)/(1,14)|0|y|
|A16(1,14)|1/2|z|

The D15 and A15 image sets are complete intrinsic automorphism orbits;
each has two marked entries. The A16 entry is reflection-invariant and
has one entry. Thus their same-component cost changes relative to the
current rule are respectively2x,2y,z−1/2. The frozen full weighted
tables confirm that both restored two-entry orbits are currently zero,
and that actual A15 and D15 each cost100.

The source extracts every corner from all24 saved whole variants and
checks its proposed lower contribution. Every nonnegative triple with
x+y+z=1/2 passes those saved variants; the tested grid has step1/16.
Each saved total is affine in these variables, and the three vertices
of that simplex are checked, so this statement also covers its interior.
The four essential branches have total1/2+x+y+z. The other saved
variants retain sufficient existing contributions. In particular,
the saturated branches have the same A15/D15 alternatives as the two
radical branches; they do not require a different coupon family.

These are whole-factor conclusions supplied by the frozen arbitrary-
support proof. This diagnostic does not independently repeat that
proof, establish retention of every other universal output after
decreasing A16, or replace the selected theorem's proper-completion
exclusions. Such a successor composition remains separate.

## Exact changes to saved forest witnesses

All443 actual fixed-type rows from the20 formerly above-curvature
representative classes are read, including every maximizer of the16
remaining actual classes. Each stored attaining forest is re-evaluated
under the current full ordinary and contextual cost function before
any proposed change is applied. These443 witnesses are not all actual
forests in those classes.

For a forest with n_D15, n_A15 and n_A16 components, put

    M = n_D15 * sum(rank(C) for every component C of rank3,4,5).

If c is a prospective small/D15 center weight, its exact cost change is

    2x*n_D15 + 2y*n_A15 + (z−1/2)*n_A16 − (3/8−c)*M.

The two short-leaf halves and rank1/rank2 D15 center roles remain
unchanged. All ordinary and contextual costs outside these four
parameters are retained. Center probes are0,1/8,1/4,3/8; there are180
tested combinations in total. The c<3/8 probes are hypothetical:
the frozen theorem proves3/8, and this diagnostic proves no stronger
center-retention theorem.

|Proposal (x,y,z;c)|Largest saved-witness value|Saved classes with a counterexample above377|
|---|---:|---:|
|Current (0,0,1/2;3/8)|763/2|16|
|D15 only (1/2,0,0;3/8)|381|14|
|A15 only (0,1/2,0;3/8)|381|14|
|Equal quarters (1/4,1/4,0;3/8)|6093/16|14|
|Equal quarters plus center0 (1/4,1/4,0;0)|6093/16|13|
|Center0 alone (0,0,1/2;0)|763/2|15|

The equal-quarter proposal would give actual A15=D15=201/2 and
A16=441/4, with no new cross coupon. The unchanged formal forest
A1+9A2+3A13 costs6093/16 for every proposal in this family, because
it contains no A15, D15 or A16. Thus these changes cannot remove all
finite obstructions, even if their universal geometry is established.
The best tested saved-witness value6093/16 occurs in88 grid proposals.

An apparent improvement on saved witnesses is a lower-obstruction
comparison only. Unrecorded forests or different residual choices may
become maximizing after the changes; representative rank15 and rank16
bases would also have to be updated. No old table is asserted to remain
optimal, and the counts13/14 are not claimed to be the complete numbers
of future obstructing classes.

## Validation and scope

Generation82643/04f872 and separate replay64119/14591c both terminated
exit0/PASS. One independent full mathematical source/report read passed,
checking the four whole-role patterns, orbit multiplicities, affine
simplex argument, current witness re-evaluation, exact cost shifts,
invariant obstruction and saved-only scope. It executed no supporting-
root, census, witness-arithmetic or optimization replay. The parent also
completed a full mathematical source/report read and passed the canonical
orbit multiplicities, four-branch identity, simplex argument, witness
alignment, center coefficient, invariant A13 obstruction and lower-bound-
only scope, without execution or arbitrary-support reproof. Source,
data and report are final/frozen after both reviews and final hash/link/AST
validation. The checks cover
frozen hashes, canonical stored roles, current table values, all24 saved
whole variants and exact shifts of443 stored forest witnesses. They do
not enumerate new actual corners, prove the proposed whole-face rule,
optimize finite forests, alter ideal payments or change the completed
unrestricted dimension58 theorem.

    python3 scripts/fresh_a16_restoration_alternative_cost_diagnostic_20261018.py --verify
