# Four more forced raw39 vertices give another arithmetic incidence constraint

The compact arithmetic H13 model over Q(sqrt5) now uses all six vertices
outside the selected pair of raw39 quadrilaterals. The four previously unused
vertices supply two additional event equalities. One gives new positive
columns and lowers the two fixed-face-mass witness ratios to
**11.985461837279875** and **12.043857645505682**. The restricted ratio
remains above the target12. No new dimension theorem or global optimum follows.

The [producer](../scripts/fresh_raw39_additional_target_frontier_20261018.py),
[certificate](fresh_raw39_additional_target_frontier_20261018.json), and
[separate saved checker](../scripts/check_raw39_additional_targets_saved_20261018.py)
extend the frozen [raw39 shared-edge theorem](fresh_raw39_shared_edge_target_incidence_20261018.md).
Every old actual Q–Q outcome and unknown pair is retained. Two small local
covering proposals supply exact rational certificates; no global optimization
is performed. All statements below concern the selected source pair, not
arbitrary quadrilateral pairs.

## The same eight actual walls force the new targets

Keep the preceding eight-wall configuration, its source type5587=A5+A5+H3,
and its passive A5+H3 of rank8. The two raw39 quadrilaterals have two common
edge endpoints of type5587. The exact 56-support census has12 positive
definite supports, including six vertices outside both quadrilaterals:

|Target type|Whole finite forest|Outside active supports|
|---|---|---|
|3324|A1+A5+H3+H4|01467,23456|
|1880|2A1+A5+2H3|02467,13456|
|2622|A1+A2+A5+H3+I2(5)|12456,12467|

The first row was already used. The other two rows are used here. Their
actual-vertex status follows from the inherited strict-center projection:
for a positive definite acute selected support, projecting a strict
timelike interior point to its orthogonal complement leaves a timelike
point strictly inside every other actual wall. Additional acute facets
cannot remove these vertices. The new source reconstructs the same complete
Gram and compares every geometric field with the frozen certificate.

The configuration also forces two compact four-simplices, fixed respectively
by wall7 and wall5, with sides01246 and12346. Each new target lies on one
of them. Four ports at that target therefore belong to an actual simplex;
the six corresponding two-faces are triangles. These are actual faces by
the inherited convexity/strict-center proof, not a hypothetical induced
subdiagram or an assumption that no extra facet exists.

## Complete event keys and unlabelled recovery

At type1880 the active frame is2A1+H3; at type2622 it isA1+A2+I2(5).
Each of the two supports of each type has four whole-component embeddings
in the canonical rank13 frame. Thus there are16 new physical maps in total.
Every map explicitly checks that the chosen five nodes are orthogonal to
the remaining eight, whose forest is exactly A5+H3.

The three outside-wall contact covectors and positive squared timelike
heights are reconstructed at each map. The full positive definite rank13
frame fixes each spacelike projection. Positive height and the common
outward future-side sign then determine the actual normal uniquely. The
three mutual native products are checked exactly.

At either target type, the two reflected presentations give the same
unordered triple of complete covectors. Deduplicating them gives **four
event keys per target type**. Their common simplex-port set is checked too.
The complete labelled native eight-wall diagram has exactly two
automorphisms: the identity and

    0<->3, 1<->2, 4<->6, 5<->7.

The latter exchanges the two source quadrilaterals, their two edge endpoints
and the two target supports. It preserves their unlabelled pair. Consequently
an actual occurrence of a key at a fixed target recovers at most one such
unlabelled source pair. Passive-component choices are not suppressed:
they are included in the complete embedding domain and distinguished by
the complete covectors.

Define an event as a pointed occurrence of an **actual source pair** at
one of its target vertices. A compatible contact pattern alone is not
declared to be an event. Every source pair forces exactly two distinct
vertices of each target type, and every actual event has one of the four
enumerated keys. If E counts the selected unlabelled source pairs, the
pointed occurrence totals therefore satisfy

    X_5587 = X_3324 = X_1880 = X_2622 = 2E.

Recovery gives uniqueness, not a claim that each algebraically admissible
key is physically realized. This distinction is sufficient for the
necessary counting constraints below.

## Necessary gates give local event capacity one

All old Q–Q edges are preserved literally. For each new event and each old
represented Q star, use the forced simplex triangle ports, or a pair of
distinct recovered outward normals. Their actual product is

    projection - sqrt(height_product).

A strictly positive value is impossible for distinct outward walls of an
acute Coxeter chamber. When projection is zero, the squared product is
height_product in Q(sqrt5). A nonintegral trace or norm violates the
necessary normalized Vinberg two-cycle integrality condition. Equal complete
covectors denote the same actual wall and are skipped after their heights
are compared. No unresolved gate is treated as an exclusion.

At type1880, all4*26=104 event/Q pairs are excluded:42 by positive product,
46 by nonintegral two-cycle and16 by triangle ports. At type2622,212 of
4*58=232 pairs are excluded:132 positive-product,72 arithmetic and8
triangle-port exclusions. The remaining20 are kept compatible, allowing
only two represented Q orbits alongside an event.

The same necessary wall gates exclude every one of the six unordered
event/event pairs at each target. At1880 there is one arithmetic exclusion
and five positive-product exclusions; at2622 all six are positive-product
exclusions. Thus at most one selected event can occur at either actual
target vertex. This capacity is proved from all possible keys, not imposed
as an optimistic graph restriction.

The complete conservative target graphs have32 and232 maximal cliques,
and18 and169 distinct count/event vectors. Every nonempty subset of every
maximal clique is included; event equalities are not inferred from a
maximal-clique monotonicity assumption. Empty local configurations satisfy
every nonnegative-cost charge inequality automatically.

## Exact four-vertex covers and the resulting columns

Combine the unchanged source5587 and target3324 records with the two new
target records. At each of the two fixed face-mass choices, a small joint
covering problem matches the common event mass across all four types and
minimizes the total nonnegative increase in their vertex multipliers.
The numerical proposal is reconstructed and checked with exact rationals.

For each type v, the saved nonnegative charge vector satisfies

    w_v*q_v <= a_v + z_v*e_v,  0<=a_v<=1,
    sum_v z_v <=0.

Summing over actual vertices and using the four equal event totals yields
a valid positive face column with vertex costs a_v. Matching rational
covering mixtures meet every corner demand, fit the increased vertex
budgets, and have exactly the same event mass. Equality of the charge lower
bound and total increase proves optimality for this local monotone repair.
It says nothing about the optimum of the full global model.

For both witnesses, the resulting useful coefficients in the type order
(5587,3324,1880,2622) are

    a=(1,0,0,1), z=(1,0,0,-1).

Thus the new exported columns charge v_5587+v_2622. The additional1880
constraint has a valid matching cover but needs no new multiplier at
these face masses. The two new columns have nine and seven positive face
coefficients, respectively, and both are checked against both repairs.

|Witness|Additional vertex mass|Updated ratio|
|---|---|---|
|Full|934391/2000000000|111428537994904130/9296974910746791|
|Named248-Q restriction|456567/2500000000|18670750401222348/1550230079993479|

Both repairs increase only the multiplier at2622. The common event mass is
934391/2000000000. All old columns, including signed-face columns with
nonnegative vertex costs and previous incidence columns, remain valid by
fixed face masses and componentwise nondecreasing vertex multipliers.
No old provenance or pair ownership is changed.

## Verification scope

Final producer generation e4e740 and separate fresh saved replay 55afaf
terminate exit0/PASS. Saved replay rebuilds the physical maps, gates and
conservative graphs, then checks the rational certificates without invoking
an optimizer. The separate checker's final run 3f1159 also terminates
exit0/PASS, including the complete embedding, event-key and pair-domain
guards. Source, data, this report and the separate checker are final/frozen.

The separate checker uses the earlier independent checker's Q(sqrt5)
arithmetic; it imports neither this producer nor the primary arithmetic
modules. It independently enumerates the complete four support-embedding
domains and the native automorphisms, reconstructs all16 saved physical
maps, and checks every new gate. For185 positive-product gates it uses
exact rational square-root intervals; all separate already at16 bits.
These are rigorous rational enclosures, not floating-point signs.
It also checks119 nonintegral two-cycles,24 triangle-port exclusions,
all264 maximal cliques and187 new count/event vectors,420 four-type charge
inequalities,62 mixture terms, both columns and both exact repairs.

Root and a separate peer completed full producer/report/helper mathematical
reads, including the inherited native recovery, source maps and reflection
interfaces. The peer checked actual source-pair event definitions, uniqueness
of the recovered unlabelled pair, two-to-two occurrence equalities, complete
event/star and event/event domains, preservation of old Q–Q outcomes, the
four-record covering/charge dual, both exported columns and monotone repairs.
The final completeness guards compare all eight map-contact keys with the
four event keys, check marker ranges and assert each complete pair domain.
These source reviews did not rerun native geometry, embeddings, gates,
graph enumeration, optimization or saved-certificate arithmetic.

Historical Q-star provenance and ownership, and the actual-vertex and
compact-face existence arguments, remain explicit inherited inputs.
Conservative compatible pairs and local covers are necessary relaxations;
they do not establish simultaneous physical realization, complete implicit
packing coverage, a new global optimum or a dimension bound. All review
is internal; no external refereeing or formal verification is claimed.
