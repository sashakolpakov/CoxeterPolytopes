# Simultaneous long-pair deletions in whole finite rank17 components

Starting from the completed [dimension68 face-payment rule](fresh_ideal5_weighted_dimension68_20261017.md),
the same-component costs of whole finite A17, B17 and D17 can be reduced
to114,105 and127 respectively. Every triangle still receives at least2
and every quadrilateral at least1. All ideal weights and all finite
cross-component coupons remain unchanged. This is a local refinement;
an exact finite-forest budget and any dimension consequence are separate
obligations.

The [source](../scripts/fresh_finite_d17_long_pair_support_20261017.py)
and [certificate](fresh_finite_d17_long_pair_support_20261017.json) reconstruct
all affected selected maps and literal replacement interfaces. The complete
all-label core classification, extension bounds and arbitrary-support sealing
arguments are explicitly inherited from the frozen predecessors.

## Canonical rules

Number A17 as the path0--1--...--16. Number B17 by the same path with
the bond15--16 doubled. Number D17 by the spine0--1--...--14 with
short leaves15 and16 attached to14. Every rule is invariant under the
full automorphism group of the corresponding labelled diagram.

Write B for all pairs at distance at most8, together with distance9 pairs
having an endpoint that is a leaf of the whole component. These basic sets
have sizes102,102,104 respectively. They are retained in their entirety.
The new paid sets, all at unit weight, are specified exactly as follows.

For A17, retain B, the following ten selected pairs, and(0,13),(3,16):

    (2,11),(2,12),(2,13),(2,14),(3,12),
    (3,13),(3,14),(4,13),(4,14),(5,14).

This gives114 pairs. The previously deleted endpoint pair(0,16) remains
deleted. Relative to the completed rule of cost135,21 further pairs are
removed. They are the reflection orbits represented by

    (0,10),(0,11),(0,12),(0,14),(0,15),
    (1,10),(1,11),(1,12),(1,13),(1,14),(1,15).

The first ten orbits have size2; the last is fixed by reflection.

For B17, retain B and(2,11),(3,12),(4,13), giving105 pairs. The labelled
B17 diagram has no nontrivial automorphism, so this is intrinsic after
choosing the doubled end as specified. Thirty-one old units are removed.

For D17, retain every old pair except

    (0,10),(0,11),(1,10),(1,11),(1,12),
    (1,13),(1,14),(1,15),(1,16).

Only the last two are exchanged by its nontrivial automorphism. The nine
deletions give cost127. They include the internal distance9 pair(1,10),
whose endpoints are not whole-component leaves. No protected leaf9 pair
is removed. The certificate also lists every retained and removed pair
and every automorphism, avoiding any dependence on a drawing convention.

## Complete routing before the finite rank17 question

The predecessor's finite triangle payments use marked distance at most8,
which is retained. Its ideal triangle payments do not change.

For quadrilaterals, the complete reduction first resolves a marked path
of length at most8. A finite such corner retains its unit. If that corner
is ideal, the universal finite8/ideal5 replacements are still available:

* For long cycles, single-contact finite alternatives have distance at most8
  or distance9 at an actual whole-component leaf. The inverse-diagonal
  sealing and fixed-side spectator arguments in the
  [finite8 replacement](fresh_finite8_cycle_replacement_20261017.md)
  establish these roles with arbitrary supporting spectators.
* The distance6 extension of that argument has the same finite-path
  alternatives. Its only small-component equality is whole A11, and its
  cycle18 equality uses separate A1/A17 halves at A17 positions5/11.
  The radius-five tree equality uses a whole A7 path of distance6.
  These are the actual conclusions of the
  [local replacement](fresh_finite8_ideal5_local_replacement_20261017.md),
  rather than an inference from its more permissive small17 statement.
* The cycle13/14 proofs use finite paths of length at most7/8, or unchanged
  separate-component coupons. The new cycle12 proof uses finite paths of
  length at most6, paths of length4, or separate A1/A11 and A2/A11 coupons.
  In particular no one-sided arbitrary-spectator branch depends on a
  long same-component rank17 unit.

After this branch, the complete low-core closure has66 chosen contexts,
each with whole marked-component bound at most12. The all-label disjoint
high-core branch is excluded geometrically. The complete high-overlap
census leaves2,614 contexts, routed as follows:

|Route|Contexts|Effect of the deletions|
|---|---:|---|
|Whole marked component at most16|2,391|No rank17 component can occur|
|Protected actual leaf9 or retained basic ideal payment|74|Payment is retained|
|Adaptive leaf9 or whole-component bound16|9|Payment is retained|
|Whole-component bound17|140|All induced finite A/B/D17 maps are enumerated|

The bound17 cases do not stand for the entire quadrilateral domain. The
other2,474 contexts, the low and disjoint cases, and the universal short-path
branch have their separate reasons above. The frozen
[A17 provenance audit](fresh_finite_a17_endpoint_payment_support_20261017.md)
already verifies these routing counts against the full82,197-identification
core census; this source binds that audit and its underlying certificates.

## All finite A17/B17/D17 selected maps

For every one of the140 bound17 cases, take its connected selected marked
component in the saved actual core Gram. A positive diagonal congruence
does not change its native angle. For distinct nonorthogonal roots define

    rho_ij = 4*g_ij^2/(g_ii*g_jj).

The values1 and2 mean native labels3 and4. Any other value cannot occur
in an induced subdiagram of finite A17, B17 or D17. Signs are checked
negative; orthogonality is kept as an explicit nonedge. Thus a rational
gauge is never mistaken for the norm-two root normalization.

The induced-map search tries every unused target vertex consistent with
every already assigned edge and nonedge. Its degree filter is only a
necessary condition. It therefore enumerates every injective labelled
induced embedding, including every possible embedding coming from an
actual whole-component extension. It need not prove that every emitted
embedding is geometrically realizable: the over-approximation is safe
because all its marked images are retained.

|Target|Eligible patterns with maps|Complete induced maps|Distinct marked pairs|
|---|---:|---:|---:|
|A17|15|48|10|
|B17|5|18|3|
|D17|75|168|20|

Every selected image is in the retained set. Patterns with no map cannot
be subdiagrams of that actual target. Finite exceptional irreducible types
have rank at most8, so these three types exhaust finite irreducible rank17.
The finite type called C17 has the same labelled Coxeter graph as B17.
Components of every other rank keep their entire preceding rule.

## Nonpreferred finite corners in literal replacements

Absence from selected maps does not suffice: an old ideal transfer may
pay at another finite corner. The checker therefore reconstructs174
complete literal configurations:

|Source|Configurations|
|---|---:|
|Old affine-order16 replacements|9|
|Old affine-order17 replacements|9|
|Universal affine D13 equality|72|
|Universal affine D14 equality|72|
|Later affine-order17 partial weights|12|

It obtains each actual corner by adding its two black marks to the saved
white set, recomputes its connected components, and enumerates full
rank17 isomorphisms into each canonical finite target. Every rank17 match
has a positive-definite full corner Gram. These are whole components:
the old literal and partial-weight cores use the actual affine kernel
and a full-core radical to show that any additional supporting contact
either gives a retained length-two marked path or vanishes on the whole
core. The universal D13/D14 equalities force all remaining spectator
contacts to vanish. This excludes both direct attachments and indirect
tails or bridges; it is not inferred merely from a displayed subdiagram.

The literal unit images are:

|Type|Canonical marked pairs|
|---|---|
|A17|(0,13),(3,16),(0,16)|
|B17|None|
|D17|(0,12),(0,13),(0,14),(0,15),(0,16),(3,14),(3,15),(3,16)|

Every listed unit is retained except the already-deleted A17 endpoint
pair. The checker simultaneously applies all three new finite rules and
recomputes all174 totals, retaining their inherited half-payments and
ideal half-weights. Every total is at least1. The sole debit is the old
affine-order17 case662: its A17 endpoint unit was already absent, and its
two D17 units both remain. Its retained total is2. In particular the new
D17 deletion does not consume the payments that justified the old A17
endpoint deletion.

The universal D13 nonextreme zero cases use whole D12, not rank17.
The immediate D13/D14 alternatives are whole affine factors of orders12/13,
whose current refinements have separate interfaces below. Thus neither
the literal equality census nor the selected-map census is being used
as a substitute for the universal nonextreme or immediate branches.

## Every later interface in the completed rule

The checker traverses the explicitly typed finite-unit corners of the
successor packages and records their actual whole-component types.
The only rank17 units are the four D17 corners of the partial-weight
package, already included among the174 configurations above.

|Later transfer|Relevant finite alternatives|
|---|---|
|Affine-order17 singleton|Separate A4/D13 halves; short path|
|Affine-order17 partial weights|Preserved D17 units; existing finite halves plus an ideal half|
|Affine-order16 partial weights|Whole D16 or separate-component halves|
|Order15 affine B reduction|Whole D14 or a finite component of order at most14 extending A12|
|Order13 affine B transfers|Whole B13 or separate A1/D12 and A2/B10 coupons|
|Universal affine D12 transfer|Whole D13/D14/A15 or separate A1/A2 versus D12 coupons|
|Order14/15 affine B zero transfers|Projection exclusions; no finite alternative|
|Cycle12 transfer|Finite paths of length at most6 or separate A11 coupons|
|Universal affine D11 transfer|Whole D12, or literal D13/D15/A16, or separate D11 coupons|
|Order12 affine B half-transfer|Whole B11/B12 or retained separate D11 halves|
|Order12 affine C deletion|No new finite alternative|

The B15 and B11 whole finite alternatives include their arbitrary-support
sealing/exclusion statements; they are not only types observed in sampled
corners. For D12 and D11, the one-sided arbitrary-spectator arguments
select the spectator-free new root's whole finite component. Their
two-sided alternatives are either retained short paths or sealed finite
corners of order at most15/16. Those conclusions cover unbounded spectator
rank without extrapolating from the saved finite tables.

All cross-component coupons remain unchanged, including the A17
positions5/11 coupons, both full A1/D11 and A1/D12 upgrades, the B10/D13
special coupons, and the two new A11 coupons. No new finite payment is
introduced. No ideal pair is altered, and none of these preserved routes
calls back to one of the newly deleted finite units.

## Verification status

Generation9167 and fresh separate-process replay5320 both terminated
exit0/PASS on the complete selected maps, all174 simultaneous literal
totals, the retained basic sets and automorphism invariance, and every
recorded later interface. All94 direct bindings, Python AST and six local
report links pass. The full geometric core census and arbitrary-support
classification/sealing remain frozen inherited inputs.

A separate full source/report and universal-interface review passed.
That reviewer reread the finite8 cycle/local proofs to check their actual
finite-path-at-most8 or protected whole-leaf9 conclusions, the A11/A7
equalities and separate A17 coupons, the cycle13/14 finite paths and
coupons, and the later arbitrary-spectator D12/D11/B11 whole components
of order at most16 or separate-component alternatives. A separate small
standard-library check validated the234 saved induced maps, basic-pair
inclusions and all174 rational totals, including the sole inherited
case662 debit retaining2. This check did not re-enumerate the maps,
recompute corner inertias, repeat the geometric census, or calculate a
finite budget.

Root completed a full source and mathematical report read, including the
imported induced-map routine's completeness, normalization in positive
gauges, and the universal finite8 cycle/local proofs. That review passed
the retained whole-leaf9/basic pairs and automorphism invariance, complete
2,614-context routing, all140 eligible finite-target interfaces, the
simultaneous174 literal debit logic with case662 retaining2, and every
newer finite output of order at most16 or between separate components.
Root did not repeat the high-core census or corner inertia calculations.
Source, certificate and report are final/frozen.
No exact finite-forest maximum, ambient dimension theorem, new polytope,
external refereeing or proof-assistant formalization is claimed here.

    python3 scripts/fresh_finite_d17_long_pair_support_20261017.py
    python3 scripts/fresh_finite_d17_long_pair_support_20261017.py --verify
