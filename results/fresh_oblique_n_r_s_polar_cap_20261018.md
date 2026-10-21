# An orthogonal cap closes R in the numerical X–S section

The new cap gives28 native numerical configurations. Each preserves all23
previous finite supports and adds exactly the finite corners N/R/T and
S/R/T. Whole R now has four finite vertices in the computed S6 section;
whole J/L/X remain compact there. The section has25 finite rays and13
open ends. Whole T remains open, so this is a concrete local closure,
not a finite-volume H9 example.

The [producer](../scripts/fresh_oblique_n_r_s_polar_cap_20261018.py),
[certificate](fresh_oblique_n_r_s_polar_cap_20261018.json) and
[standalone saved checker](../scripts/check_oblique_n_r_s_polar_cap_saved_20261018.py)
bind the frozen [profile0 X–S pilot](fresh_oblique_x_s_diagonal_native_pilot_20261018.md)
and its [six-profile extension](fresh_oblique_x_s_nonzero_r_profiles_20261018.md).
Neither predecessor is changed. This test adds a wall instead of
replacing one of the closing normals already fixed in those packets.

## Actual domain and finite duties

The full domain consists of all12 native predecessor configurations:
eight with Rprofile0, at R–S labels2,3,4,5 in each of two parents,
and four with Rprofile1 or4 and R–S label5. The two parents retain all
four original selected origins. No projected-vector deduplication
identifies different white profiles or discards their provenance.

For every configuration reconstruct J67,L68,X69,N73,S72 and R70 from
the complete six-step instruction chain. Original63/64 remain only
the early recipe anchors; original30/44/63/64/65 and new M71 are absent
from the retained chamber. The N–X and S–52 labels are12, N–10 is3,
and all predecessor normal branches and profiles are fixed.

Add fresh wall T74 with each of the seven prescribed white profiles.
Solve in the physical frame N/R/S with all three prescribed products
zero, corresponding to labels2,2,2. A nonzero white profile changes
both the projected norm and the projected right-hand sides:

    projected T² = 2 − white_projection(T)²,
    projected T·V = −white_projection(T)·white_projection(V).

In particular the nonzero R profiles must not be treated as profile0.
All old native offsets and the candidate's full common-center offset
are retained. The norm problem has both branches checked whenever
the exact finite-duty filter allows the recipe.

The genuine required finite corners are N/R/T and S/R/T. N/S is dotted;
there is no finite N/S/T requirement. Before solving normals, the
inherited exact native finite-diagram classifier tests those two
corners, with the two contacted white roots and every declared profile.
For Rprofile0, all seven T profiles pass this necessary filter at
every old R–S label2–5. For each nonzero Rprofile1/4 with R–S label5,
only Tprofile0 passes. Thus the complete84 recipes split into24
finite-duty rejects and60 eligible recipes.

All60 eligible three-frames are nonsingular and have two distinct real
norm roots. Every one of the120 roots is saved. Each is tested against
all52 retained projected original normals and six fixed new normals.
The six fixed white roots are supplied by its prescribed profile, so
there are58 original facets and64 total fixed facets before adding T.
Native, norm-residual and strict full common-center gates precede
every cone computation. A surviving branch must also pass both genuine
corner checks and a mutual recheck of all seven new walls.

The outcome is:

|R profile|Old R–S label|Native T profiles, in each parent|
|---|---:|---|
|0|2|0,1,2,4,5|
|0|3 or4|0,1,4|
|0|5|0|
|1 or4|5|0|

There are28 native caps. Sixty opposite norm branches have positive
common-center product and numerous old-wall violations. The remaining32
rejected branches fail only the X69 native product. The saved record
contains every blocker and physical product; these are numerical
native tests, not directed-interval exclusions. No norm recipe remains
unresolved in this declared domain.

## New disk and complete numerical face ledgers

Append triangles N/R/T and S/R/T to the frozen ten-triangle disk.
They attach along the consecutive boundary edges N–R and R–S, replacing
that boundary path by N–T–S. The new disk has10 vertices,21 edges and12
triangles, with boundary

    8–10–N–T–S–52–8.

All interior links are cycles and all boundary links are paths. R has
cycle X–N–T–S–X; T has boundary link N–R–S. The edge N–S is absent,
and no N/S/T triangle is introduced. The ideal-agent peer independently
confirmed this literal topology only, without a geometric conclusion.

All28 admitted configurations have the complete twelve prescribed finite
corners. The full ray sets of J/L/X/R agree exactly with their disk stars.
Whole8,11,12 are also compact in the numerical section. Whole10,52,55,
N/S/T and the other inherited open walls remain open.

The selected cone uses13 actual old boundary sides and seven new walls.
Its clean cone and independently enumerated positive-definite supports
agree at25 finite rays,11 exterior rays and13 open ends. A separate
cone using every one of the52 retained old projected constraints plus
the seven new walls has25 finite rays,15 exterior rays and13 open ends,
with no warnings. Its full40-ray and edge ledger is saved for every cap.
As in the predecessor, the other39 projected original normals are
timelike with the inherited negative orientation; they are redundant
on the hyperbolic sheet but can subdivide the exterior projective region.

Every prior finite support survives. Exactly two are added:

    {R,S,T}, {R,N,T}.

The complete exterior changes have two different patterns:

|Configurations|Removed exterior supports|Added exterior supports|
|---|---|---|
|10 caps with old R–S right|{52,N,S}, {N,R,S}|{52,N,T}, {52,S,T}|
|Other18 caps|{N,R,S}|{N,S,T}|

Thus every cap removes the targeted N/R/S exterior support and closes R.
Only the first ten also remove52/N/S. In those ten the exterior region
moves to52/N/T and52/S/T; in the other eighteen52/N/S remains. T always
has exterior rays. The full open-edge calculation, not just a finite-ray
count, gives13 ends in every case. The inherited remote regions at
7/47/55,33/59/66 and58/59/60/66 are not repaired by this cap.

## Original support debt and scope

The exact deletion set is still30/44/63/64/65. The original199-support
certificate loses19 finite supports and touches12 compact two-faces
and one compact three-face; three lost supports involve old8. The
standalone checker independently reconstructs this entire debt from
the original support list. No further original facet is released.

The preservation of23 existing and addition of two finite supports is
a numerical statement about the computed section. It does not turn
the original19/12/1 deletion account into a completed ambient chamber
certificate. New normals are not promoted to exact H9 normals, all
remaining whole-face obligations remain, and the13 exterior ends
prevent a finite-volume conclusion even for the computed section.

## Verification and limits

Initial generation4560/d49f17 passed. The final producer additionally
records the new disk and exact whole-star comparisons; its source-bound
generation33574/f2e1d0 and separate replay86811/24111d terminated exit0.
Final standalone checker63090/807624 also terminated exit0/PASS.
The first checker attempt stopped on an overgeneralized assertion that
all caps remove52/N/S. The saved producer data already distinguished
the two patterns above; the checker was corrected to verify both without
changing the producer or any geometric result.

The helper imports NumPy and no project modules. It independently
verifies all12 complete predecessor chains, all84 recipe positions,
finite-duty admissibility via the physical white-plus-three-root Gram,
both distinct solutions of each nonsingular frame, all physical products,
profile offsets, native/center gates, finite-support enumeration, both
topologies and the exact original debt. It checks all28 saved full
40-ray ledgers against every projected inequality, rebuilds their edges,
and verifies the complete finite/exterior support differences and whole
J/L/X/R closures. It does not independently solve parent normals or
enumerate geometric cones, certify intervals, or reprove the inherited
scaffold, native finite catalogue or seven-profile theorem.

The root agent completed a full source/report/standalone-helper
mathematical-scope read, including a separate reread of the cap topology.
It checked all12×7 recipes, the genuine N/R/T and S/R/T finite duties,
physical profile offsets, both norm branches and native/center gates,
the12-triangle disk and whole J/L/X/R closures, independent finite-support
enumeration, the full40-ray/edge ledgers, both10-versus18 exterior
patterns and the original exact19/12/1 debt. This was a read-only review:
no producer/helper execution, directed-interval verification, promotion
of new normals to exact ones, or reproof of the historical scaffold,
native catalogue or seven-profile theorem was performed by the reviewer.
The final metadata check verifies direct bindings, source syntax and
local report links; it is not another geometric check.

Other T contacts, different parents, labels
outside the two frozen inputs and further caps are untested. No new
dimension bound, ideal-vertex certificate or finite-volume polytope is
asserted.

    python3 scripts/fresh_oblique_n_r_s_polar_cap_20261018.py --verify
    python3 scripts/check_oblique_n_r_s_polar_cap_saved_20261018.py
