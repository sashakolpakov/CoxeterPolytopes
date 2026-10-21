# Two central walls: a valid disk topology, but no Coxeter improvement in the bounded domain

The [source](../scripts/fresh_nonarithmetic_two_central_caps_20261017.py) and
[saved result](fresh_nonarithmetic_two_central_caps_20261017.json) test a different
patch topology on the twenty frozen four-wall rings. Exact combinatorics verifies
the proposed disk. The geometric search and its complete intrinsic boundary
scores are numerical diagnostics, not exact chamber certificates.

## Topology and deletion costs

With old boundary cycle `a,b,c,d = 7,12,52,30` and ring `L,P,R,Q`, use
the peripheral triangles `abL,bcL,cdP,daQ,cLP,dPR,dRQ,aQL`. The two central
vertices give `SLP,SPR,TLQ,TQR,SLT,SRT`; the other diagonal rotates the
ring once. The source checks all links, edge multiplicities, the boundary,
and Euler characteristic. There are ten vertices, 23 edges, and fourteen
triangles. The old–new and new–new edge counts are 8 and 11; the triangle
counts by number of new vertices are 4,4,6. Thus the conditional open-end
change is `2(8+11)-3(14)=-4`.

This identity assumes the indicated old boundary is genuine and retained,
the constructed disk has the stated incidences, and its vertices are simple
and finite. It establishes no geometric existence. A separate conceptual
read by the frontier agent confirmed these combinatorial checks.

The exact control has 64 walls, 198 finite ambient vertices, and fifteen
open ends on S6. The changed family deletes old47 and53 and replaces12,
as in the frozen parent ring construction. Its inherited deletion ledger
loses the four old finite supports S6+[7,8,12], S6+[7,12,47],
S6+[8,12,52], and S6+[12,47,52]. Old53 has no finite vertex and has a
future timelike S6 projection. The previously certified compact faces
touched by these changes are the 8-heptagon and 12-quadrilateral.
Intrinsic redundancy does not preserve other ambient faces. No new
ambient finite-vertex census is claimed here.

## Bounded domain and result

Both diagonals of all twenty saved rings are tested. Each central wall
uses all seven necessary A2 profiles, both norm branches, and the same
contact alphabet `2..12,16,24,32` used by the parent scout. The three
framing contacts and the two adjacent finite endpoint conditions are
imposed; the fourth ring contact, every retained native old contact
(including projection offsets), and the common center are checked.
All coherent pairs receive a complete numerical S6 projective-cone and
edge census. This is not an all-label or all-parameter search.

Among 12,996 pairs, 4,215 pass the numerical native-contact gates. Every
one has its two central walls dotted, and every resulting boundary has
35 finite vertices, seven exterior vertices, and fifteen open ends.
There is no finite central join, no complete intended disk, and no
improved Coxeter chamber in this domain. The other 8,781 pairs fail a
numerical acute/angle gate. No near-ideal ray or degenerate face was
reported in the scored cones at the stated numerical tolerances.

## Continuous-angle eleven-end diagnostic

Allowing just the mutual central-wall angle to be non-Coxeter isolates
316 pairs with both required central endpoint Grams positive definite.
Each has a complete numerical boundary with 37 finite vertices, five
exterior vertices, and eleven open ends. All other tested native contacts
pass the numerical Coxeter/dotted gates. The topology therefore produces
the desired boundary change as ordinary hyperbolic geometry in this
sample; its Coxeter angle is the remaining obstruction.

The closest saved finite-endpoint relaxation is parent stage25,
diagonal0, with mutual product approximately -0.000638488019376382,
strictly away numerically from the right-angle value0. Its boundary has
52 compact, seven one-open, two two-open, and two exterior-only edges.
The source also saves the other relaxation recipes and scores. These
near-angle gaps have not been exactified. There is no exact new wall,
eleven-end Coxeter chamber, or finite-volume conclusion.

## Verification and continuation scope

Final generation32554 and fresh saved replay52845 both terminated with
exit0/PASS. The replay recomputed the complete bounded search and all
saved boundary scores; it is a producer replay, not independent numerical
verification. Source manifests, syntax, and report links were checked.

A possible unsearched successor changes the central degrees from4/4 to3/5:
one wall meets all four ring walls while a second truncates one nonfinite
ring corner. The frozen single-center theorem requires all four of the
first wall's corners finite and therefore does not exclude this topology.
Preliminary unsaved numerical high-contact reductions are not a certificate
and are not used in this report. No further construction search was launched
after this package's final replay.
