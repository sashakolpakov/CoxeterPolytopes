# Exact global facet selection for H9/14

## Result and theorem-facing status

The 142 graph shadows surviving the strict global GLBT and ridge-support
filters have now been subjected to two further exact layers.

First, polytopal-link constraints reduce the graph-level frontier as

\[
142\longrightarrow72\longrightarrow39\longrightarrow36\longrightarrow30.
\]

The successive conditions are: every vertex figure is nonsimplex
(\(g_1>0\)); strict vertex-figure \(g_2\); feasibility of the strict
vertex-figure \(g_3\) incidence intervals; and the corresponding edge- and
triangle-figure constraints.  The corrected nonsimplex condition matters: a
\(D\)-simplex figure is impossible because its boundary has a missing face of
order \(D+1>5\).

Second, an exact Boolean model globally selects facet complements and defines
every face of the resulting complex.  A bounded 30-second-per-Z3-call baseline
and solver-independent refinements give the following deliberately
theorem-facing partition:

| status | count | meaning |
|---|---:|---|
| independently certified impossible | **123** | the previous union of 113 certificates, plus exact solver-independent certificates for records 92, 93, 54, 23, 74, 46, 37, 34, 16, and 6 |
| explicit ordinary-polytopal selections | **19** | records 0, 1, 3, 4, 5, 7, 9, 12, 14, 17, 18, 30, 31, 32, 35, 42, 43, 44, and 47, independently audited |
| **total** | **142** | complete coverage of the second-layer graph shadows |

Equivalently, after the 123 proved whole-shadow rejections the residual has
nineteen members, each with at least one explicit ordinary-polytopal
selection.  The nineteen
are **not** hyperbolic Coxeter examples: ordinary polytopality supplies only
the combinatorial boundary sphere.  All nineteen designated selections are
now excluded as compact Coxeter nerves.  The last seven previously untested
selections—records 0, 1, 3, 4, 5, 7, and 9—have a complete independently
replayed Lanner/actual-vertex zero-state certificate.  The earlier eleven are
records 14, 17, 18, 30, 31, 32, 35, 42, 43, 44, and 47; record 12 is covered
by its stronger whole-shadow theorem.  Record 17 is also stronger than
selected-only: its exact eight-orbit
completion census is exhausted by an independently replayed compact-label
no-go.  The later corrected record-12 census exhausts that shadow into 371
abstract compact residual candidates, and complete Lanner/actual-vertex joins
leave zero label states on every one, with two no-import replays.  Thus record
12 is also a whole-shadow compact no-go.  Record 47 is now the third: all four
Dehn--Sommerville schedules are closed by strict link sums, a complete
`a=1` forced-shadow test, and a lossless `a=0` census whose unique 90-facet
orbit is the already label-obstructed selected nerve.  Record 44 is the
fourth: the `a=1` branch closes `690 -> 378 -> 0`, and the
`a=0` census closes `2,782 -> 2,214 -> 502 -> 1`, again ending at the frozen
selected nerve.  Record 43 is the fifth: strict link sums close two schedules,
the `a=1` branch closes after its forced-shadow test, and the `a=0` census has
exactly two full-graph incidence orbits; unrestricted compact-label gates
exclude both.  Record 42 is the sixth: the corrected `a=0` census covers all
82 active-row configurations and fifteen exact targets, while `a=1` has one
product orbit; complete compact-label gates exclude every branch. Record 35
is the seventh: all five compact schedules are exhausted, `a=2` is empty, the
lone `a=1` orbit has an exact index-two obstruction, and a lossless
cap-admissible `a=0` census leaves eighteen target families, all rejected by
complete Lanner/actual-vertex joins. The necessary topology audits and label
census are independently replayed. Record 32 is the eighth: strict compact
link positivity removes `a=4,5`; complete censuses leave 292 `a=0` target
orbits, ten `a=1` target orbits, the selected nerve as the sole `a=2` target,
and no `a=3` target.  Independent full-action orbit comparison and fresh
additions/blocker/ridge/Lanner/cross replays close every branch. Record 18 is
the ninth: exact searches exhaust all seven schedules, and complete compact
label gates kill all surviving `a=0,1,2` target orbits after literal
production/independent set equality.  Record 30 is the tenth: strict
global/link bounds plus exact certificates
close all eight schedules.  Its final `a=0` branch has eight production
targets, all compact-label impossible; a dual CNF with a checked DRAT proof,
a no-import orbit/clause replay, and a second full proof replay certify that
the list is exhaustive.  Record 14 is the eleventh: its final `a=1` and
`a=0` branches have literal Z3/CP-SAT equality, exact full-source target
quotients, and zero unrestricted compact-label states on all ten plus three
target orbits.  For the eight open positive shadows
`0, 1, 3, 4, 5, 7, 9, 31`, another completion remains possible.  A
changed noncompact ideal-face structure is a separate question for all of
them.  The compact-facing frontier is therefore 134 wholly excluded graph
shadows and eight open; the ordinary-incidence
partition remains `123 + 19`.  The schedule-level frontier inside Record 31
is now `a=0,1`: exact independent certificates exclude `a=2,3,4,5,6,7`.
Record 14's schedule-level frontier is empty: strict summed
vertex-link `g4` closes `a=5,6`, and complete production/independent
enumerations close all 59 `a=3` plus all 13 `a=4` source orbits with zero
incidence targets.  At `a=2`, independent Z3 and CP-SAT censuses agree on 24
literal pure targets from all 36 source orbits; the targets form six exact
fixed-source orbits, and complete compact-label joins kill all six.  At
`a=1`, the solvers agree on 19 generator-slice targets in ten full orbits;
at `a=0`, they agree on 52 targets in three full `S9 x S2` orbits.  All
thirteen final orbits die at the unrestricted compact-label gate.  See
`results/h9_14_record14_upper_schedules_no_go.md` and
`results/h9_14_record14_a2_compact_schedule_no_go.md`,
`results/h9_14_record14_a1_compact_schedule_no_go.md`, and
`results/h9_14_record14_whole_shadow_compact_no_go.md`.
For the seven open stellar-multiwedge shadows 0, 1, 3, 4, 5, 7, and 9, a
separate strict summed vertex-link `g4` theorem closes the top two schedules
of every shadow, fourteen schedule branches total.  Their summed values are
12 and 6, contradicting the required minimum 14.  No whole shadow follows
from this cutoff; see `results/h9_14_open_high_schedule_g4_no_go.md`.
For `a=2`, strict face-link positivity removes 31 of 49 positive-`A4` orbits;
two exact incidence quotients agree on twelve ridge-valid targets from the
remaining eighteen, and independent unrestricted compact-label joins kill all
twelve.
Inside Record 0 at `a=1`, all four source orbits are now closed.  Source 0
exhausts the 108 labelled order-four Lanner placements in nine disjoint
orbits of sizes `12,12,24,12,12,12,6,12,6`; all nine complete joint
incidence/finite-label models are terminal `INFEASIBLE`.  Source 1 closes all
nineteen of its placement orbits, source 2 is independently
incidence-impossible, and source 3 is closed by its exact 41-core joint
CEGAR theorem.  Hence the whole compact Record-0 `a=1` schedule is
impossible.  Record 0 as a whole remains open in schedules `a=0` and
`a=2,...,16`; the global compact `H9/14` problem and every finite-volume
noncompact variant remain open.  See the
[`source-0 aggregate theorem`](h9_14_record0_a1_source0_joint_branch_certificate.md),
[`source-1 aggregate theorem`](h9_14_record0_a1_source1_joint_branch_certificate.md),
[`source-2 theorem`](h9_14_record0_a1_source2_compact_incidence_no_go.md),
and [`source-3 theorem`](h9_14_record0_a1_source3_terminal_cegar_audit.md).
At Record 0, `a=2`, the exact source reduction is
`386,760 -> 384,780 -> 36` full missing-graph orbits.  Independent CP-SAT
and Z3 models now exclude eight fixed sources
`{0,8,9,10,12,23,31,34}` at the compact-incidence gate.  The frozen scout
ledger has 24 source orbits with explicit incidence witnesses and four
solver-open sources `{24,29,30,32}`.  This is only an eight-source tranche:
the whole `a=2` schedule and Record 0 remain open, and no finite-volume
noncompact conclusion is made.  See
[`the a=2 incidence tranche`](h9_14_record0_a2_incidence_tranche.md).
The Record-31 `a=1` census is **not complete**.  Its corrected
Dehn--Sommerville vector is
`(14,87,319,770,1281,1477,1142,540,120)`, so the ridge gate is `{2:540}`;
the earlier `{2:531}` transcript is invalid.  Source orbit 0 has one exact
topology-valid incidence target, already killed by its unrestricted compact
label join; a fresh 242,124-node direct matrix quotient independently
reproduces its unique target orbit.  Production source orbit 3, mask 8220, is now terminal: its exact
469,990,577-node quotient has 4,445 pure 120-facet targets and exactly eleven
ridge-valid incidence orbits.  All eleven pass the connected mod-two
homology-sphere audit, but production and hostile-order unrestricted compact
label joins leave zero states on every orbit, hence no Gram candidate.  Its
dense canonical keys exposed a severe memory-scaling failure mode; a separate
lossless incidence-row encoding reduces each key to at most 22 bytes and has
an exact 231,773-node dense-versus-compact equivalence replay.  See
`results/h9_14_record31_a1_orbit3_corrected_quotient.md` and
`results/h9_14_record31_a1_compact_key_equivalence.md`.  Source orbit 5,
mask 14340, is now independently closed at the incidence gate: the terminal
production and direct matrix quotients search 41,488,368 and 6,042,612 nodes,
respectively, agree exactly on 735 pure 120-facet targets and all 231 ridge
signatures, and find zero `{2:540}` targets.  See
`results/h9_14_record31_a1_orbit5_corrected_quotients.md`.  Production
source 2 and independent-incidence gaps remain, so no `a=1` schedule theorem
is claimed.  Source 4 is now incidence-impossible by an exact partition of
its six canonical root classes.  Classes 1--3 and the terminal 105,066,986-node
classes-4--6 fragment both have empty canonical `{2:540}` target-orbit sets;
the latter contains 3,184 pure 120-facet targets and 818 ridge signatures.
Only the empty target-orbit union is combined, not the nonadditive internal
fragment counts.  See
`results/h9_14_record31_a1_orbit4_selective_restart.md`.
Record 30 is closed: a complete 21-orbit, 77-target incidence and label
replay excludes `a=2`, while a strict-link theorem plus independent exact
incidence enumerators and complete label joins exclude `a=1`; the dual
CNF/DRAT theorem excludes `a=0`.  The Record-30
`a=3,4` results are complete
minimal-transversal/facet-count no-gos, while Record-31 `a=3` leaves three
topology-valid incidence candidates and kills all three at the unrestricted
compact-label gate.  In Record 30's `a=2` schedule, sparse Lanner/vertex joins
kill 74 of the 77 target orbits.  The selected product target and the two
648-state exceptions are excluded by complete whole-nerve gates; the latter
materialize 1,332 exact nonzero rank-eleven principal determinants per nerve.
See `results/h9_14_record30_a2_compact_schedule_no_go.md`.
The independently replayed `a=1` theorem is
`results/h9_14_record30_a1_compact_schedule_no_go.md`.
The schedule aggregation and final `a=0` exhaustion are
`results/h9_14_record30_whole_shadow_compact_no_go.md`.
Records 43, 44, and 47 pass the
complete incidence, link-homology, and uniform rank-four oriented-matroid
gates, and an independent multiwedge factorization proves ordinary convex
polytopality.  Records 35 and 32 have separate exact wedge/join realizations.
Records 0, 1, 3, 4, 5, 7, 9, 12, 14, 17, 18, and 31 have exact rational or
integral inverse-multiwedge realizations.  The selected record-14, 17, 18,
and 31 nerves have separate unrestricted compact-label no-gos; the newer
all-orbit record-17 theorem subsumes its selected-nerve result.
Record 16 has an exact all-schedule blocker/ridge exclusion; records 46, 37,
and 34 have the earlier exact blocker/ridge exclusions.  Record 6 is excluded
by a solver-independent two-graph theorem: the complete special-link census
and facet-count identity derive ten branches, all rejected by relative-orbit
enumeration and a no-import hostile replay.

This is a substantial nerve-level reduction, but it is **not** yet the desired
14-facet nonexistence theorem.  In particular, the nineteen positive selections
show that the combinatorial conditions are not intended to reject all ordinary
simplicial 9-polytopes.  The three multiwedge examples can even be carried by
SNC anticanonical toric boundary divisors; record 35 has a projective rational
\(\mathbb Q\)-factorial toric carrier, without a smoothness or SNC claim, and
record 32 has a smooth projective toric/SNC incidence carrier.  Coxeter
labeling and Lorentzian Gram feasibility remain separate gates in general.
The seven earlier displayed selections have a subsequent unrestricted compact-Coxeter
obstruction, but that does not exclude another completion of one of their
graph shadows or a noncompact finite-volume chamber with a changed ideal-face
structure.

## The exact facet-complement model

Let \(H\) be a surviving missing-edge graph on \(V=[14]\), and let
\(\mathcal C(H)\) be its archived ridge-support core of five-element vertex
covers.  For every \(B\in\mathcal C(H)\), introduce a Boolean variable \(x_B\).
The selected \(B\)'s are complements of nine-vertex facets.

### Closed-pseudomanifold equations

For every six-set \(C\subset V\), introduce \(r_C\in\{0,1\}\) and impose

\[
 \sum_{\substack{B\in\mathcal C(H)\\B\subset C}}x_B=2r_C.
\]

Thus the complementary eight-set is either not a ridge or lies in exactly two
facets.  Zero is allowed.

Every independent triple \(T\) is required to be covered:

\[
 \bigvee_{B\cap T=\varnothing}x_B.
\]

### Faces are defined in both directions

For every independent set \(S\) with \(4\le |S|\le8\), a face variable
\(y_S\) is defined by the biconditional

\[
 y_S\quad\Longleftrightarrow\quad
 \bigvee_{B\cap S=\varnothing}x_B.
\]

Nine-faces are the selected facets themselves.  This is stronger than merely
matching the number of selected facets: the solver counts the actual faces
generated by the selection.

One complete archived Dehn--Sommerville schedule is imposed.  If \(a\) is the
number of missing four-faces, then

\[
 a=I_4(H)-\sum_{|S|=4}y_S,
\]

and every entry \(f_3,\ldots,f_8\) must equal the corresponding entry of the
same schedule.  The implementation verifies that the archived schedule family
is consecutive and affine before using the compact linear encoding.

### The Lanner order bound is global, not post hoc

For every independent \(S\) of order at least six, impose

\[
 y_S\ \lor\ \bigvee_{v\in S}\neg y_{S\setminus\{v\}},
\]

where \(y_S\) is false for \(|S|\ge10\).  Hence no inclusion-minimal nonface
has order greater than five.  A raw ridge pseudomanifold can fail this badly;
the condition cannot be replaced by facet-count matching.

### Polytopal links

The selected missing four-faces obey every exact vertex-, edge-, and
triangle-incidence interval from the local-link calculation.  Once all face
variables are present, the solver also imposes:

- vertex figures: strict \(g_2,g_3,g_4\);
- edge figures: strict \(g_2,g_3\);
- triangle figures: strict \(g_2\);
- nonsimplex vertex count \(f_0\ge D+2\) for each such \(D\)-figure; and
- every relevant Macaulay inequality
  \(g_{i+1}\le g_i^{\langle i\rangle}\), encoded by a finite exact table.

The independent witness auditor recomputes these inequalities from the
generated face table; they are not trusted merely because Z3 asserted them.

## Solver-independent certificates and positive fixtures

### The forced-\(a=0\) minimal-transversal rejection

Exactly 36 of the 142 shadows force \(a=|\mathcal A_4|=0\).  For 35 of
them Dehn--Sommerville forces one missing five-face; the remaining graph
forces none.  The complete choice space consists of 30,934 possible
\(\mathcal A_5\) families.  For every choice, the independent archive gives
an inclusion-minimal transversal of order three or four, instead of the
required order five.  It verifies all choices directly, without Z3.  The
witness-size distribution is 3,192 of order three and 27,742 of order four.

Thirty-five of these graphs are already rejected by the local-link layer;
one lies among its 30 survivors.  Hence the union contains
\(112+36-35=113\) independently rejected shadows and leaves exactly 29.

The unique non-overlapping case is record 134 and has an especially short
certificate:

it has one forced schedule with

\[
 |\mathcal A_4|=|\mathcal A_5|=0
\]

and four missing edges.  Its entire minimal-nonface clutter would therefore be
the four edges of \(H\).  Facet complements would be exactly the
inclusion-minimal vertex covers of \(H\).  Choosing one endpoint of each edge
gives a cover of order at most four, and shrinking it gives a minimal cover of
order at most four.  This contradicts the required facet-complement order
five.  The archived certificate records the numbers \(0,0,4,4,5\); no solver
is needed to verify the contradiction.

### Records 92 and 93: strict link \(g_4\) plus blocker duality

Two further members of that 29-shadow frontier now have independent exact
exclusions.  For either record, if \(a\in\{0,1,2\}\) is the number of missing
four-faces, incidence double counting over the fourteen vertex links gives

\[
 \sum_v g_4(\operatorname{lk}v)
 =5f_4-24f_3+63f_2-112f_1+14\cdot126=18-6a.
\]

Every vertex link has strict \(g_4\geq1\), by the inherited Lanner-order
bound and the generalized lower-bound equality theorem.  Thus the sums 12
and 6 exclude \(a=1,2\), leaving only \(a=0\).  Dehn--Sommerville then forces
exactly three missing five-faces and 78 facets.

For each graph the verifier classifies all triples of independent five-sets
by how they meet the three minimal graph covers, rejects every triple having
a transversal of order at most four, and counts the remaining minimal
five-transversals.  Purity requires every minimal transversal to have order
five, and the face vector requires exactly 78 of them.  The exhaustive exact
counts are:

| record | pattern-compatible triples | no-small-transversal triples | possible numbers of five-transversals |
|---:|---:|---:|---|
| 92 | 6,773,256 | 1,176 | 57 (840 triples), 75 (336 triples) |
| 93 | 12,605,586 | 32,256 | 54 (7,560), 57 (22,680), 63 (1,260), 75 (756) |

Neither row ever gives the required 78.  This is finite integer enumeration
with an independently replayable verifier and no SAT/SMT backend.  The full
classifications are in `results/h9_14_record92_certificate.md` and
`results/h9_14_record93_certificate.md`.

### Record 54: two schedules, both exhausted

For record 54 the same double count is

\[
 \sum_v g_4(\operatorname{lk}v)=24-6a,
\]

so strict positivity leaves only \(a=0,1\).  In the \(a=1\) branch, all 718
possible missing four-faces form six exact graph-automorphism types.  Local
\(g_4\) eliminates four types; each of the remaining two has exactly 840
no-small-transversal triples.  Their five-transversal counts lie in
\(\{60,62\}\) or \(\{64,65,66,67\}\), never the scheduled 94.

For \(a=0\), four missing five-faces remain.  Strict local \(g_4\) makes their
active parts pairwise disjoint and forbids an isolate from lying in all four.
A bijective four-bit column encoding reduces the complete search to 235
active-part multisets and 1,422 incidence states, or 159 true
\(\operatorname{Aut}(H)\times S_8\times S_4\) orbits.  These represent exactly
2,229,164 labelled unordered families with no transversal of size at most
four.  Their exact five-transversal histogram has twenty values from 60 to
98 but omits the scheduled 96.  Thus record 54 is independently impossible;
the full compression and histogram are in
`results/h9_14_record54_certificate.md`.

### Record 37: five schedules and a complete five-row blocker census

For record 37 the five Dehn--Sommerville schedules have missing-four counts
(a=0,1,2,3,4), facet counts (102,100,98,96,94), and summed vertex-link
(g_4) values (30,24,18,12,6). Strict positivity excludes (a=3,4).

There are 462 individually interval-feasible missing four-faces. For (a=1),
all 462 forced five-set shadows make a vertex-link (g_4) nonpositive. For
(a=2), 96,033 of the 106,491 pairs satisfy the family intervals: 51,429
already force more than the scheduled 17 missing five-faces, and all remaining
44,604 force a nonpositive link (g_4).

For (a=0), Dehn--Sommerville requires five missing five-faces and 102
facets. A complete five-row quotient by the nine isolated graph vertices
examines 2,865,450 row-normalized states, of which 522,427 have no transversal
of order at most four. Exactly 12 states attain 102 minimal
five-transversals; they deduplicate to five explicit families. Two have ridge
histogram (1:28,2:403,3:28), and three have
(1:24,2:411,3:24). Hence every last candidate violates the degree-two
ridge condition. The frozen census and an independent direct replay are
documented in `results/h9_14_record37_certificate.md`.

### Two explicit selections

Two baseline timeouts are resolved constructively.

1. Record 30 is the nerve of
   \(\Delta_1^3\times\Delta_3^2\).  It has 128 facets, missing-face orders
   \(2,2,2,4,4\), and
   \[
   f=(14,88,328,804,1352,1568,1216,576,128).
   \]
2. Record 42 is the nerve of
   \(C_5\times\Delta_3\times\Delta_4\).  It has 100 facets, missing-face
   orders \(2,2,2,2,2,4,5\), and
   \[
   f=(14,86,309,725,1165,1295,970,450,100).
   \]

For both, the refinement archive stores the complete selected complement
list.  A backend-independent audit verifies the support core, every ridge
degree, dual connectivity, triple coverage, the complete face vector, all
minimal nonfaces, all local strict \(g\)- and Macaulay constraints, spherical
Euler characteristic of every link, and exact mod-two homology of the complex
and every face link.  These are ordinary polytopal calibration objects, not
hyperbolic Coxeter examples.

### Records 43, 44, and 47: exact multiwedge polytopes

The generalized blocker screen first resolved record 43 constructively. Its
missing graph has edges

\[
\{0,11\},\{0,12\},\{1,12\},\{0,13\},\{1,13\},\{2,13\},
\]

and adjoining the four missing five-faces

\[
\{3,4,5,6,7\},\quad
\{1,2,8,9,10\},\quad
\{2,8,9,10,11\},\quad
\{8,9,10,11,12\}
\]

gives a clutter with exactly 90 inclusion-minimal transversals, all of order
five.  Their complements generate the face vector

\[
(14,85,300,690,1088,1190,880,405,90).
\]

The independent audit verifies ridge degree two, connected dual graph, every
local strict and Macaulay inequality, spherical Euler characteristic, and
mod-two sphere homology for every face link.  A frozen 1,001-sign uniform
rank-four Gale chirotope independently passes all 45,045
Grassmann--Pluecker relations and has exactly the 90 facet complements as its
positive circuits. Records 44 and 47 have independent certificates of the
same strength.

A later exact factorization proves more. Each of the three nerves is

\[
\partial\Delta^4*K^0(J),
\]

where \(K^0\) is either a hexagon boundary or an explicitly stacked
seven-vertex 2-sphere, and \(J\) replaces one vertex by three or four clones.
The dual operation is an iterated facet wedge, which preserves convex simple
polytopality; the join is dual to product with \(\Delta^4\). Thus all three
are nerves of ordinary convex simple 9-polytopes with fourteen facets. The
construction can be chosen Delzant, making each nerve the dual intersection
complex of the reduced SNC anticanonical boundary on a smooth projective
toric 9-fold. See `results/h9_14_multiwedge_polytopality.md`.

This independent realization does not assert that the archived oriented
matroid is the realized one, and it supplies no Coxeter labels or Lorentzian
Gram matrix.

### Record 35: rational multiwedge and compact-Coxeter no-go

Record 35 has five missing edges and five missing five-faces.  Its complete
blocker consists of 102 five-sets and generates

\[
(14,86,309,726,1171,1309,986,459,102).
\]

Two clone blocks of orders three and four collapse the nerve to a nine-vertex
base.  An explicit integral realization of that base in \(\mathbb R^4\) has
all 126 homogeneous five-point determinants nonzero and exactly the required
19 supporting tetrahedra.  Five simplicial wedges recover the record-35
simplicial 9-polytope.  A rational polar and rational wedges also give a
projective \(\mathbb Q\)-factorial toric 9-fold carrier.  No Delzant,
smoothness, or SNC assertion is made.  The exact construction is in
`results/h9_14_record35_polytopality.md`.

The same nerve is nevertheless excluded from the compact Coxeter track.  A
complete join of the order-five Lanner placement tables, checked against all
102 actual rank-nine vertices, leaves two intermediate symmetry orbits; each
contains an affine \(\widetilde D_4\sqcup A_1\) or \(\widetilde D_5\)
principal subdiagram in an actual vertex.  In normalized Gram form
\(G=I-A\), adding any other non-right finite Coxeter label increases the
nonnegative matrix \(A\) entrywise.  Perron--Frobenius monotonicity therefore
cannot restore positive definiteness.  This is an unrestricted compact-label
obstruction for record 35, not merely a sparse-ansatz failure.  It does not
rule out a noncompact finite-volume chamber having a changed ideal-face
structure; see `results/h9_14_record35_coxeter_label_gate.md`.

### Record 34: complete forced-shadow and ridge exclusion

Record 34 has missing-edge graph (K_3\sqcup K_2\sqcup9I). Exact local
link arithmetic and forced-A5 shadow unions exclude every schedule with at
least one missing four-face. For (a=0), a complete bicoloured six-row/nine-
column incidence census checks 1,288,153 isomorphism representatives and
26,598,204 active colourings. Only 18 states attain the required 114 pure
order-five transversals, and every one has a ridge of degree one, three, or
four. Thus no target state is a pseudomanifold. The solver-independent proof
and direct (2^{14})-subset replay are in
`results/h9_14_record34_certificate.md`.

### Record 32: convex selection and whole-shadow compact no-go

Record 32 now has a 110-facet selected nerve with

\[
f=(14,87,318,760,1242,1400,1060,495,110).
\]

Collapsing a common three-vertex clone block leaves the boundary of a
pentagonal bipyramid; two simplicial wedges and a join with
\(\partial\Delta^4\) recover the full nerve.  This proves ordinary convex
polytopality and gives a smooth projective toric 9-fold whose SNC boundary
records the incidence complex.  It does not supply a hyperbolic metric.
Every compact Coxeter labeling of this selected nerve is excluded by the
complete overlapping-A4/disjoint-A5 Lanner census and actual-vertex replay.
A distinct displayed `a=0` selection is audited separately and is also
compact-Coxeter impossible, although its convex polytopality is not claimed.
The later whole-shadow theorem closes all six Dehn--Sommerville schedules:
`a=4,5` fail strict compact link positivity, `a=3` is empty, and the exact
`a=0,1,2` target-orbit counts are 292, ten, and one.  Complete independent
replays agree on the full target-orbit set and every incidence and unrestricted
label obstruction.  This is compact-only and makes no finite-volume
noncompact inference. See `results/h9_14_record32_polytopality.md`,
`results/h9_14_record32_a0_selection.md`, and
`results/h9_14_record32_whole_shadow_compact_no_go.md`.

### Record 16: complete graph-shadow exclusion

Record 16 is excluded as a whole graph shadow, not just at one selected
nerve.  Exact ridge equations force one blocker graph to be a complete
bipartite cut.  An analytic degree bound closes `a=0,1,2`, an exhaustive
triangle-free small-graph replay closes all but one branch for `a=3,4`, and
strict vertex-link \(g_4\) closes `a=5,6`.  The final branch forces both a
`K4,6` ridge cut and a `K5,5` two-row blocker.  Direct dualization of its 102
forced facet complements supplies an independent contradiction.  The
timed-out Z3 result is marked untrusted.  See
`results/h9_14_record16_certificate.md`.

### Eight rational stellar-multiwedge completions

The former timeout records 0, 1, 3, 4, 5, 7, 9, and 12 each have an exact
ordinary convex realization.  Every construction starts with `C(8,5)`, makes
two rational radial-beyond stellar subdivisions, and applies four simplicial
wedges (allowing a repeated carrier).  Exact supporting-hyperplane
determinants recover the base boundary, and independent hostile replays
reconstruct the expanded minimal-nonface clutter and complete blocker.

| record | stellar carriers | wedge carriers | facets | missing-four `a` |
|---:|---|---|---:|---:|
| 0 | `01`, `6,8` | `0,1,6,8` | 168 | 6 |
| 1 | `012`, `013` | `1,3,4,8` | 152 | 4 |
| 3 | `01`, `4,8` | `2,3,6,9` | 162 | 4 |
| 4 | `01`, `07` | `3,4,5,6` | 146 | 2 |
| 5 | `02`, `13` | `0,2,4,8` | 136 | 3 |
| 7 | `01`, `035` | `3,4,5,6` | 130 | 0 |
| 9 | `23`, `057` | `5,6,6,7` | 126 | 2 |
| 12 | `024`, `0,1,8` | `1,2,3,4` | 126 | 2 |

For every row, all facet complements have order five, all ridges have degree
two, the dual graph is connected, the scheduled face vector is exact, and
every face link has the mod-two homology of the required sphere.  The detailed
coordinate and certificate reports are
`results/h9_14_record0_polytopality.md`,
`results/h9_14_record1_polytopality.md`,
`results/h9_14_record3_polytopality.md`,
`results/h9_14_record4_polytopality.md`,
`results/h9_14_record5_polytopality.md`,
`results/h9_14_record7_polytopality.md`,
`results/h9_14_record9_polytopality.md`, and
`results/h9_14_record12_polytopality.md`.  These are eight selected ordinary
polytopes, not Coxeter or hyperbolic examples.

### Seven selected nerves: compact Coxeter no-go

For each displayed ordinary-polytopal selection in records 30, 32, 35, 42,
43, 44, and 47, exact Lanner placement and every actual rank-nine vertex
leave either a locally nonelliptic block or a principal Gram block with two
negative directions.  Perron--Frobenius monotonicity and interlacing make the
obstruction unrestricted over all remaining finite Coxeter labels.  This is
a theorem about the seven frozen nerve completions only.  It does not reject
another completion of the same graph shadow and does not apply Vinberg's
finite-covolume criterion to an ideal-face structure.  See
`results/h9_14_selected_nerve_coxeter_no_go.md`.

Selected records 14, 17, 18, and 31 have separate unrestricted compact-label
certificates.  In record 31 the complete cross-edge census leaves only exact
ten-wall blocks of inertia `(8,2,0)` or residual states with no simultaneous
Lanner completion; see `results/h9_14_record31_coxeter_label_gate.md`.

### The remaining seven ordinary selections: compact Coxeter no-go

The selected nerves in records 0, 1, 3, 4, 5, 7, and 9 are now excluded by a
uniform exact gate.  It reconstructs every actual rank-nine vertex, joins the
complete 108 and 420 labelled Lanner placements in orders four and five, and
sets all not-yet-forced finite pairs right angled.  The seven searches reach
zero after `5,6,5,4,3,5,3` decisive supports, respectively.  A hostile replay
reads only the frozen polytopality archives, regenerates the classical Lanner
tables, uses a separate finite Coxeter classifier, and matches every stage and
state digest.  Perron--Frobenius monotonicity makes the sparse rejections
unrestricted over the remaining finite labels.  See
`results/h9_14_open_selected_coxeter_no_go.md`.

This completes the compact-label gate for the designated ordinary selection
in every one of the nineteen positive incidence shadows.  It does not change
the whole-shadow ledger: another completion remains possible for records
0, 1, 3, 4, 5, 7, 9, 14, and 31.

## Bounded backend evidence

The baseline archive records 112 immediate local-link rejections, three Z3
UNSAT dispositions, and 27 timeouts.  Independent blocker/transversal
certificates now upgrade all three Z3 cases: records 92, 54, and 23.  Thus the
backend archive remains corroborating evidence, but no theorem-facing status
depends only on it. Independently, the local-divisibility certificate for
record 74, the four-A5 certificate for record 46, the five-A5 certificate
for record 37, the six-row certificate for record 34, and the all-schedule
certificate for record 16 upgrade five of the 27 bounded timeouts to
impossibility.  Records 0, 1, 3, 4, 5, 7, 9, 12, 14, 17, 18, 31, 32, 35, 43,
44, and 47 are upgraded to ordinary polytopality.  Record 6, the last
theorem-undecided shadow, is now independently excluded by its exact ten-
branch two-graph census.  The earlier exploratory UNSAT result is not used:
a tracked core exposed its dropped xor operand, and the corrected bounded
full-model run timed out.  The theorem instead rests entirely on the direct
link/orbit enumeration and its hostile replay.  There is now no undecided
graph shadow.

A secondary schedule-split run fixes \(|\mathcal A_4|\) before solving.  For
record 92, all three branches \(a=0,1,2\) return Z3 UNSAT.  For record 93, the
\(a=0\) branch returns UNSAT while \(a=1,2\) time out in the frozen run.  These
results remain useful propagation diagnostics, but they are no longer the
proof: the link-\(g_4\) and blocker certificates above independently exclude
both records.

## Independent topological audit and its limit

Given a selected facet list, the verifier constructs the entire Boolean face
table on \(2^{14}\) subsets.  In addition to Euler checks, it builds every
face link and computes all boundary ranks over \(\mathbf F_2\) by integer-bit
Gaussian elimination.  Passing this audit is necessary for a sphere, but is
not a PL-sphere certificate: mod-two homology does not detect odd torsion,
simple connectivity, or shellability.

The global Boolean solver does not encode homology ranks internally. A SAT
witness is therefore audited after construction; if it failed, that would
reject only the witness, not the graph.  Every one of the nineteen displayed
ordinary-polytopal selections passes the same complete link-homology audit.
Homology alone would not promote them to PL spheres; their separate product,
wedge/join, Gale-coordinate, or rational stellar-multiwedge constructions
supply the independent convex-polytope proofs.

## Reproducibility and frozen artifacts

The optional theorem-search dependency is pinned as
`z3-solver==4.16.0.0`.  The pure witness, local-link, and homology verifiers do
not import Z3.

| artifact | file SHA-256 | internal certificate SHA-256 |
|---|---|---|
| local links | `651285ca7e8d531b501d800118c042cb64fdc804c5381cd6dba59d384264052f` | `946137ef8fe8f38b7f9f7df46261a0ad617060d286d28ab7e853e7b245a75b74` |
| forced-\(a=0\) transversals | `aabf0128bc71d24b272f52802e114eb338ffb5145dad2d415b622148e4239c5a` | `4abf18addd6301b15edf7b9ef99a257aa61c94f160c524127ab23f0c828af587` |
| bounded facet baseline | `2d329439e5baeffb574932ff7bbab2dd78daf706b2917a0674bdc9f01a2f5123` | `2789e320b0b23122d53fcf27b68d9198d575af64377a485dde9b973af77a8f0c` |
| solver-independent refinement | `68a9439dc4621c8247cb9c5ae625f5d62bf670cd368417fe4defaad7dff04496` | `c976ac262bb78dc4020210e3eab8ed662d5719de0d058f1a94089c60c7f4ae25` |
| schedule split | `f9e3d149f4da9bc05220c84c1ffdde5ba8cb8ed9d388df94992d20d5f200b6f5` | `2d4e6aa60912a9060c62a07a8d6997dda00fd4eec82d60e850c41fab523f2e1c` |
| record 23 blocker | `1b50d1c42f8e221c1b1d20320e2742341d7566895f6b92f9166eab17ff227539` | `7b13a6b215c4cf6c2c871361c984f150f044a7107bce5b5fb9f484096360da5c` |
| record 54 blocker | `43dcc59d73190bf9daa46e143497f715e69151741c1fd7ad289164fb335978b3` | `911955924b1ef679e5a12c0c22a0c11184f0b1ed35ca9b6d2459502ef0b2f295` |
| record 74 divisibility | `43c7c537a3a2f04030c80b83dadd30e8e803edf792a23bb67021b019dd463b19` | `0a8812c4dbade45cabe3c544c968d265581e4d3f0fdd6b48673ff7847148ceb4` |
| record 43 selection and Gale signs | `2c6463598ab4ae2073b54fc8a9345ab099acba0492c786ba163d0ad5d29d3c66` | `af521ed4e7476c95945df9b770f36a343dcf3fc79dabed557aaab0fd97a671a3` |
| record 43 whole-shadow compact-label obstruction | `413f30f6c7d521e7a816134b5ef28c0661e5efae6431b90a72172bdba0441985` | `5b606985e5cbf4938a99060941135a3655c9c282345227dd58df06add3dbe65b` |
| record 42 whole-shadow compact-label obstruction | `c9beee7a93965774bff273b89c95430a2b9504b9d094b4730f78033aa659e74c` | `4590a0f0365dd26abd416cf4f73b909028a2863913f7c6d179a2a3e248122e9d` |
| record 44 selection and Gale signs | `8110d99a901f9102906117943c23a34c6c2d434b17c8f92f9f6e015d25c81fff` | `e66e493d0c68f7ae28ffe03e63f96209276a8d83a5ceb27c498f35b42c4a1867` |
| record 44 whole-shadow compact-label obstruction | `eaa4a69c94e6f85bcf047c3fccd811f1179887370a46591b6c43b19da50cb206` | `51a91b44a01554c6423ace7fb3b31f47a6d334ce9d9581dbd295864c4a420d35` |
| record 47 selection and Gale signs | `163da338e18d9ce8d81d1957244b9b948839b949302ee09aaf27f858fbe2c357` | `960ef6564c6b8013caf0989c95a46cc8401c846d32d47c66c057091340301546` |
| record 47 whole-shadow compact-label obstruction | `58cd3c5afd8183bbb39756a7601a47ae298311fb0d2a7f991ef4f0ec0da694e4` | `324f9a1f4bc7ed7156c8b42c30886609df2339b761f9dc1f2bebca37438b17f1` |
| record 46 four-A5 blocker | `a07dcc91611e73a0dd12f3c8f54e740483ac1cc46976fa1e8890c61260f558aa` | `5f28acdd31ae4c043e799c4dfb62c63dd4b22f95af7032c6b8af83a5fa04b3f3` |
| record 37 five-A5 blocker/ridge exclusion | `cdb7df45b1a886cc1d689c34a2a2ac0cc786afd8edd19c2770eee1ddaacccdad` | `4da8617bd4dfc77889c3c3f987e7988efd87d6d366d7efd29c7403370536cc5c` |
| record 34 all-schedule blocker/ridge exclusion | `7a7dda1a0bcc44076e6c1f40e5e47e577288c39c8bb958814d52399485723830` | `79d4de8b4c9d463bc0aca44bfd78a3c8eb1ee6bfdc211f90b47ad661e5837db2` |
| records 43/44/47 multiwedge/toric polytopality | `99753f5ed95304f9b2d804020beb3c3e787898112f4189120471c8e744f4d3b9` | `d5efd409cf5d5a5a62e9784f8c9e71eebb5a08ee20e8ea7f1d05e894d305a183` |
| record 35 rational multiwedge polytopality | `4d9a53f34bd313d183352291a2273efe7e3618761bdd21b8635a2beb08f6421e` | `bc535f94895c6e151c5bd05db142a7735e1260bc736012daa6ccbe5787d00db7` |
| record 35 unrestricted compact-label obstruction | `728da7b1d2287a3ac6983900d104502f33918ff10bcb9864d8aee4a094464175` | `fcb5b0374c1ee68b32193491dcb1d905a7ac39ef61e7c56a7bbfe6bf6b56d23c` |
| records 0/1/3/4/5/7/9 selected compact-label obstruction | `0dd28552e3b0cea7211b98c8d629618dc7a08b47bd884c8731c7df8cf269243b` | `c0f693a5ec92abb0f146d7919d7c2a22910fa07c0576d402a9700b18ef142ed7` |
| records 0/1/3/4/5/7/9 high-schedule strict-link obstruction | `8dc8fe7399a24a130323bbb9e0c7b2ec8865acc12bbaee948c24b88b7d38180e` | `a7f9d6c4ca31c91d6d456555703a273a67a70bf16ad74ab54b26afeb735156f5` |
| record 14 upper-schedule compact obstruction | `a86faeb81a004c3829af37cd436929116678c545e00e536890c38e5a06da665c` | `9392e15c3d1c446150febb215573a7ecbe1a0843a31e6ed541b0d9896bd7fa8c` |
| record 14 `a=2` complete compact-label obstruction | `5c2894c59a789dffd0869c570854d5833c3342bc8bc5ac7bbb5a9fd79bf19a6f` | `d47d7c062ae1b40a3dca8bb1f3ffe6c3f87b5444e5c0d78a8515e2caa4aab7a1` |
| record 14 `a=1` complete compact-label obstruction | `e064898f4396a7db19499b8b4c9086116dec49f0bdf4b85be2e60f7a135bec16` | `7c579a3ac2260275a504d7a5a2fa2cf851cb2eb379c040069e130a07cc1da4f3` |
| record 14 `a=0`/whole-shadow compact obstruction | `827a4bc7aa224972c6fe7bdfc3cea3e9f6c6fcb53fc44c201c92121bbeb9d466` | `75fa29552f12864ca668e82550a0d63c9b643d5b114f85a15fd641c4b523cd0e` |
| record 14 `a=0` hostile downstream replay | `75e61af5bdfedc796cd5dcf80ab870c4d2b86ab5ebcde26769b697c0c33e3c37` | `e871f807999f07e3e89a213669d34fab8ff717abec74c16fe9a96a61337966c3` |
| record 32 convex wedge/join realization | `9a6f5211b1ad6972938b28e8a2c4740bf0c3ef046ac52df269621cb723dedfd7` | `94a65380b2d56d4787a3de28d5cca5dede0540cc2ec0ff23b2a3d65b4cfbd318` |
| record 32 selected compact-label obstruction | `728d3bd5640b69e7584a061766db102c11e74cb117fd1a852ebdba8433094632` | `15bf62a37b5f65fdf7118e80be347482d30b0e47b84eeab2dd831f089b5b4412` |
| record 32 `a=0` selected completion/no-go | `07fba1cec56899f5ad428fe91e2a213cdbb7ec9fcc4b06c7981768ef3e1c7f30` | `14c3f9f5a425eb6934d7f53a7a7f071467eba7cf2e4a64b7cf9a194d3490d0bb` |
| record 32 whole-shadow compact no-go | `600d107e9739085897349cf905d411df1af260afb75facb589c0edadf74a52b0` | `486b063921841de9aa3c4611ea171801c3368c1c1c7c5cdcef938fabedfb69d9` |
| record 16 all-schedule graph-shadow exclusion | `2c4377089da7487c1d0694d5abd9759e6ac1ca5dde75a66b209fc46f459b8888` | `ac7f386a36ffe8a0ddbeab74c0fdd60913eb891365735588ee7a111fc8b7c8d5` |
| record 6 ten-branch two-graph exclusion | `a2036104134ade994b0ea5ab3eb3da9c3eaa3e86c701e71c5a34245b3c0c994e` | `cc0c3186b16cafc3c3a9572ea15f6d24880ad7bb3b18da40cf8c635d064f78ff` |
| record 0 rational stellar-multiwedge | `fbc1ff944775f8017f716bedf103a2f03068b6b4af722eebc5897dda75c563e6` | `9b2ab118ebe84e31ba72a2c8b57cf24bf0be12a9314a29356ac457be26da51a4` |
| record 0 `a=1` source-3 terminal compact CEGAR | `afe350d601d76707026475c328479c99c4080809e4665a30b585b696fa61f1c4` | terminal model `7485e42496c5f9f8f03cedde468a9f7180b98c9ee92ab6f86a23e68a0b735ec3` |
| record 1 rational stellar-multiwedge | `3779b364b4c03597299c1c92967635a34d170894fda8ac57260fafc6c0a1c1aa` | `5684d466a7bca45828da56aae3bd1ea91c5390e08f966c98bccf71b427f2b976` |
| record 3 rational stellar-multiwedge | `6bcef790610c5c320b12f4fdc2a497ad7245754d34138b7367a6dfbd3146ff14` | `500c0874643b21cc5f6c0eadaf278695d9516db1832a895afb84042fa1da4685` |
| record 4 rational stellar-multiwedge | `41ee036d8d6bbb72ae3bf638999f070ffd8ec2ed74dfb7be6066fecbff57c77a` | `0625af6609333638803df22868791179f3a659c47819bd762a07f1bd069cb1fc` |
| record 5 rational stellar-multiwedge | `a6ae11237043ca541dd30c10b6a9c0b49c1afdbce85720bc28ba37ef8921f391` | `6cdd1b7951093bb6765a3d77376e3bbad50d9f470d1752d0abc80eb061cdb015` |
| record 7 rational stellar-multiwedge | `cbe18cf69dec6d57dc4862ee055691167f12672658655b4546e2fef6d6a9bc86` | `788c8d4a16456810e719e7c0f973149edc6e0a29bab067af28d89c47f6600759` |
| record 9 rational stellar-multiwedge | `e4bb73314e3e3e53f7bc900e52b6978ad2e1f8c69a57727cff28ba91b7bd548a` | `844b80aba9da0431fdefcdc9a2fbd93b5eb72e37f81fb0344b33a7c00848fdd2` |
| record 12 rational stellar-multiwedge | `943c76986b507b86faacdbbe20e978ec863514be2f90a1747f72cdccf0017ac1` | `cd98e07ccb15d04654c951e3271239be88d842b21c7bfa5f3c7f4c48d5875ab1` |
| seven selected nerves, re-audited unrestricted compact-label obstruction | `028bbc8e266145e59fa1ee9f9796a6e2c241f4a4f09637bd73f98bced9908629` | `c2dfcbfc03dc5dbae34938ab32f3a1c1111fd729d7844fd76d04dfc6cfa44de8` |
| record 30 complete compact `a=2` schedule no-go | `06759a184708d9c2b8a4ab15fc772fc2a4b7e66ac275d2ff537bebf92cb9ebe6` | `3ac5e4ea4ddd624aef17dec45f64d403165b736bab015d088e2b8d78f43a5adf` |
| record 30 complete compact `a=1` schedule no-go | `65152e8f4b687f487852932a5dae3896d189a56d1ed99c6764f7389aaea8d341` | `932d69e5646106785171562174e5423329ac09137e0224ac65e6ce436187067d` |
| record 31 selected unrestricted compact-label obstruction | `a69e933efcdb47c4e60073870b154cbfaa375a7e1274ce4f2c0acae359277e94` | `e2aa40c2647ecf6ba04b3ffebdcaf6df2b324abeea3369568054aafe9c7b111b` |
| record 17 all-eight-orbit compact-label obstruction | `adc5b678b170fca82f7d816ec1f1a9b96c07be93b31b964d80cfc068bce7ad5b` | `59a0ebdd63a28914412cda34d60289e06a539dcf3ae4c8d65debc9e7b2f7d820` |
| record 12 all-371-orbit compact-label obstruction | `fe6c0ee7f494a67ecac88b51aed301a4fdfe26fbddbf68034916f7b7aecc97eb` | `8e06523196d30b344a0b7cec5771d968d75ab253e5c4a9ac493bdf0dec21d3d1` |
| record 92 blocker | `299e1dce5040f9473b9db168730c0302ec5551607cb55f385c8c71698256b02a` | `c0e0903caa31469d3bd95d1e10cd42a62704a1de1859cfe370749065f05a5c25` |
| record 93 blocker | `696dec108c636ef5ebdc5c9018da87bd0bd126b484689d5d18bf097c92532c0b` | `e6f0ba1d325af94a0a31bb19eab9d0eb79041354a20427a04b6a606e00d8d55b` |

Commands:

```text
PYTHONPATH=src python3 scripts/h9_14_local_link_filter.py --verify-only
PYTHONPATH=src python3 scripts/h9_14_forced_a0_transversals.py --verify-only
PYTHONPATH=src python3 scripts/h9_14_facet_completion.py --verify-only
PYTHONPATH=src python3 scripts/h9_14_facet_refinement.py --verify-only
PYTHONPATH=src python3 scripts/h9_14_schedule_split.py --verify-only
PYTHONPATH=src:. python3 scripts/certify_h9_14_record23.py --verify-only
PYTHONPATH=src:. python3 scripts/certify_h9_14_record74.py --verify-only
PYTHONPATH=src:. python3 scripts/certify_h9_14_record43.py
PYTHONPATH=src:. python3 scripts/independent_h9_14_record43_audit.py
PYTHONPATH=src:. python3 scripts/certify_h9_14_record43_all_completions.py
PYTHONPATH=src:. python3 scripts/independent_h9_14_record43_all_completions_audit.py
PYTHONPATH=src:. python3 scripts/certify_h9_14_record42_all_completions.py --verify-only
python3 scripts/independent_h9_14_record42_all_completions_audit.py
PYTHONPATH=src:. python3 scripts/certify_h9_14_record44.py
PYTHONPATH=src:. python3 scripts/independent_h9_14_record44_audit.py
PYTHONPATH=src:. python3 scripts/certify_h9_14_record44_all_completions.py
PYTHONPATH=src:. python3 scripts/independent_h9_14_record44_all_completions_audit.py
PYTHONPATH=src:. python3 scripts/certify_h9_14_record47.py
PYTHONPATH=src:. python3 scripts/independent_h9_14_record47_audit.py
PYTHONPATH=src:. python3 scripts/certify_h9_14_record47_all_completions.py
PYTHONPATH=src:. python3 scripts/independent_h9_14_record47_all_completions_audit.py
PYTHONPATH=src:. python3 scripts/certify_h9_14_record46.py --verify-only --recompute
PYTHONPATH=src:. python3 scripts/certify_h9_14_record37.py --verify-only
PYTHONPATH=src:. python3 scripts/independent_h9_14_record37_audit.py
PYTHONPATH=src python3 scripts/certify_h9_14_record34.py --verify-only
PYTHONPATH=src:. python3 scripts/certify_h9_14_multiwedge.py
PYTHONPATH=src:. python3 scripts/independent_h9_14_multiwedge_audit.py
PYTHONPATH=src:. python3 scripts/certify_h9_14_record35.py
PYTHONPATH=src python3 scripts/certify_h9_14_record35_coxeter_labels.py --verify
PYTHONPATH=src python3 scripts/certify_h9_14_record92.py --verify-only
PYTHONPATH=src python3 scripts/certify_h9_14_record93.py --verify-only
PYTHONPATH=src python3 scripts/certify_h9_14_record54.py --verify-only
PYTHONPATH=src:. python3 scripts/certify_h9_14_record32.py
PYTHONPATH=src:. python3 scripts/independent_h9_14_record32_audit.py
PYTHONPATH=src:. python3 scripts/certify_h9_14_record32_coxeter_labels.py --verify
PYTHONPATH=src:. python3 scripts/certify_h9_14_record32_a0.py
PYTHONPATH=src:. python3 scripts/independent_h9_14_record32_a0_audit.py
PYTHONPATH=src:. python3 scripts/certify_h9_14_record16.py --verify-only
PYTHONPATH=src:. python3 scripts/independent_h9_14_record16_audit.py
PYTHONPATH=src:. python3 scripts/certify_h9_14_selected_nerve_coxeter_no_go.py --verify
PYTHONPATH=src:. python3 scripts/independent_h9_14_selected_nerve_coxeter_no_go_audit.py
PYTHONPATH=src:. python3 scripts/certify_h9_14_record30_a2_schedule_no_go.py --verify-only
PYTHONPATH=src:. python3 scripts/certify_h9_14_record30_whole_shadow_compact_no_go.py
PYTHONPATH=src:. python3 scripts/certify_h9_14_record14_a1_schedule_no_go.py
PYTHONPATH=src:. python3 scripts/certify_h9_14_record14_a0_whole_shadow_no_go.py
PYTHONPATH=src:. python3 scripts/independent_h9_14_record14_a0_downstream_audit.py
PYTHONPATH=src:. python3 scripts/certify_h9_14_record0.py --verify-only
PYTHONPATH=src:. python3 scripts/certify_h9_14_record1.py
PYTHONPATH=src:. python3 scripts/certify_h9_14_record3.py
PYTHONPATH=src:. python3 scripts/certify_h9_14_record4.py --verify-only
PYTHONPATH=src:. python3 scripts/certify_h9_14_record5.py
PYTHONPATH=src:. python3 scripts/certify_h9_14_record7.py
PYTHONPATH=src:. python3 scripts/certify_h9_14_record9.py
PYTHONPATH=src:. python3 scripts/certify_h9_14_record12.py --verify-only
PYTHONPATH=src:. python3 scripts/certify_h9_14_record31_coxeter_labels.py --verify
PYTHONPATH=src:. python3 scripts/independent_h9_14_record31_coxeter_audit.py
PYTHONPATH=src:. python3 scripts/certify_h9_14_record17_all_coxeter_labels.py --verify
PYTHONPATH=src:. python3 scripts/independent_h9_14_record17_all_coxeter_audit.py
PYTHONPATH=src:. python3 scripts/certify_h9_14_record6.py
PYTHONPATH=src:. python3 scripts/independent_h9_14_record6_audit.py
```

## Remaining work

The graph-shadow incidence frontier is closed.  Record 6, the final holdout,
has exactly ten numerical two-graph branches; all ten fail either special-link
compatibility, independent-triple coverage, or the forced isolate-only
missing-four count.  Both the production certificate and a no-import hostile
replay reproduce the complete failure histogram.  This closes incidence only,
not the Coxeter theorem.

The nineteen positive shadows have one certified ordinary-polytopal
completion each.  Record 17 is now settled on the compact track: its completed
eight-orbit census and unrestricted all-orbit label gate leave zero states.
See `results/h9_14_record17_all_coxeter_label_gate.md`.
Record 12, whose missing graph is `K2,2 + 10I`, now has an exact
automorphism-quotiented whole-shadow census of 371 abstract residual incidence
complexes, with `a=0,...,4` histogram `48,95,144,63,21`.  Every row passes the
complete blocker, ridge, face-table, local `g`/Macaulay, Euler-link, and
mod-two homology-link audits.  This does not prove PL or ordinary polytopality;
the separately certified rational witness is one `a=2` row.  Complete
Lanner-placement joins and every actual rank-nine vertex restriction now
leave zero compact label states on all 371 rows, independently reproduced by
two no-import implementations.  See
`results/h9_14_record12_corrected_371_census.md` and
`results/h9_14_record12_all_coxeter_label_gate.md`.  Record 47 now has an
independently replayed all-schedule census with exactly one residual incidence
orbit; equality with the selected nerve transfers its unrestricted label
no-go to the entire shadow.  See
`results/h9_14_record47_whole_shadow_compact_no_go.md`.  Record 44 now has
the analogous all-schedule uniqueness and transferred label
no-go; see `results/h9_14_record44_whole_shadow_compact_no_go.md`.  Record 43
has exactly two residual incidence orbits and both have unrestricted
compact-label no-gos; see
`results/h9_14_record43_whole_shadow_compact_no_go.md`.  Record 42 has one
`a=0` orbit after the corrected all-82-pattern census and
one `a=1` product orbit, both compact-label impossible; see
`results/h9_14_record42_whole_shadow_compact_no_go.md`. Record 35 has three
necessary-topology-audit `a=0` orbits and one `a=1` orbit, and every compact
label branch is excluded; see
`results/h9_14_record35_whole_shadow_compact_no_go.md`. Record 32 has 292
`a=0` target orbits, ten `a=1` target orbits, the selected nerve as its sole
`a=2` target, and no `a=3` target; every compact branch is independently
replayed and excluded. See
`results/h9_14_record32_whole_shadow_compact_no_go.md`. Record 18 exhausts all
seven schedules, with exact literal-family comparisons and downstream
compact-label no-gos; see
`results/h9_14_record18_whole_shadow_compact_no_go.md`. Record 30 exhausts all
eight schedules, with a checked dual CNF/DRAT proof at
`a=0`; see
`results/h9_14_record30_whole_shadow_compact_no_go.md`. Record 14 exhausts all
seven schedules with two-solver target equality and unrestricted label
extinction; see `results/h9_14_record14_whole_shadow_compact_no_go.md`.  The
eight open positive shadows are exactly `0, 1, 3, 4, 5, 7, 9, 31`.

All nineteen designated ordinary-polytopal nerves are compact-Coxeter
impossible by unrestricted exact label gates;
records 17, 12, 47, 44, 43, 42, 35, 32, 18, and 14 are additionally impossible
across every completion orbit compatible with the compact necessities, as is
record 30 by its all-schedule aggregate theorem.  The
seven newly constructed selected nerves (records 0, 1, 3, 4, 5, 7, and 9)
now have a common independently replayed zero-state Lanner/actual-vertex
certificate.  They, together with record 31, still need whole-shadow
closure.
Any surviving label system must then pass the full global Gram rank
and signature test.

A genuinely noncompact finite-volume search is a separate branch.  It must
include the ideal-face data and satisfy Vinberg's literal stopping test:
every elliptic rank-eight edge diagram has exactly two endpoints, each either
an elliptic rank-nine diagram or one full componentwise-affine zero-
neighborhood of Coxeter rank eight.  Every connected affine subdiagram must
also have such a full rank-eight zero-neighborhood; at a nonsimple cusp this
may have more than nine nodes and several affine components, but it is one
endpoint.  No present result is a complete H9/14 nonexistence theorem or a
hyperbolic Coxeter example.
