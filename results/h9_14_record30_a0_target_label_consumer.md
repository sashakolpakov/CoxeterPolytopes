# Record 30 `a=0`: downstream compact-label consumer

## Status

The production connected-ridge incidence census has now terminated.  It
visited `28,039,413` search nodes, found `19` accepted labelled paths, and
reduced them to exactly **eight** target orbits.  A standalone downstream
consumer kills all eight orbits, with terminal compact-label histogram

```text
0 states: 8 targets.
```

Thus the production computation gives a compact `a=0` no-go.  Promotion to a
publication-grade whole-schedule theorem is deliberately withheld until a
structurally independent replay establishes that the eight-orbit incidence
list is exhaustive.  This remaining qualification concerns census
exhaustivity, not a surviving Coxeter-label or Gram state: there are no such
survivors in the production target list.

The consumer uses only the Python standard library and makes no project-module
import:

```text
scripts/independent_h9_14_record30_a0_target_labels.py
```

No universal contradiction has been proved from the numerical support data
alone (`3` dotted pairs, `8` order-five minimal nonfaces, and `132` actual
vertices).  The actual vertex family and the overlaps among the eight
supports materially enter the label CSP.  The current theorem remains the
target-specific no-go for the explicit multiwedge candidate until the
eight-orbit production census receives its independent exhaustion replay;
the production result itself is already stronger and has zero label survivors.

A hostile support audit confirms why the stronger statement is not yet
available.  Requiring eight supports merely to carry mutually compatible
members of the complete 420-placement Lanner table has solutions; failure
appears only after ordinary-face ellipticity is imposed.  Two additional
pure-`5^132` support fixtures,

```text
(59,61,122,124,248,1984,3968,14080)
(59,121,122,248,1984,3968,6016,7936),
```

are both killed by the consumer, with deterministic stage counts

```text
420 -> 96  -> 6  -> 0
420 -> 126 -> 18 -> 0.
```

They expose a promising local classification problem: several supports can
contain a common triple, and the remaining two-element pieces form a small
graph whose type controls the Lanner join.  Every graph type tested so far is
label-impossible once the relevant actual vertices are imposed.  However,
neither all such graph types nor the assertion that every admissible target
contains a decisive common-triple cluster has been exhausted.  It is therefore
a search heuristic, not a theorem or a pruning rule.

## Why the production incidence search is exhaustive

The search works in the dual facet-complement language.  A facet of a simple
nine-polytope on fourteen walls has a five-element complement, and avoidance
of the fixed missing graph `3K2` says that this complement meets each of its
three edges.  There are exactly `326` such five-covers.  They split, according
to the number of active matching vertices, as

```text
3 active + 2 isolates: 224
4 active + 1 isolate:   96
5 active + 0 isolates:   6.
```

Any 132-element target therefore contains a complement of the first type, and
`Aut(3K2) x S8` is transitive on that type.  The canonical seed
`[0,1,2,3,4]` consequently loses no target.

For each selected complement `C` and each of the nine vertices outside it,
the six-set `C+v` indexes one dual ridge group.  In a closed simple boundary
every nonempty such group has exactly two selected members, and the dual
facet graph is connected.  Starting from the seed, recursively choosing the
second member of an open ridge therefore reaches every possible connected
target.  Unit ridge partners are logical consequences and are propagated to
a fixed point.  A saturated ridge forbids all further members; these are
exact degree-two implications, not heuristic pruning.

The only other early rejections are monotone upper bounds from the scheduled
face vector.  Once exactly 1,364 five-faces have appeared, the remaining eight
independent five-sets are forced.  Their literal blocker is the unique possible
completion; the code checks that it has 132 five-complements, contains the
current partial family, is closed and connected, and then rebuilds the complete
face table, minimal-nonface family, and frozen local-link caps.  Thus this
shortcut decides the whole branch rather than sampling one completion.

Partial states are quotiented by a collision-free colored incidence graph:
active vertices, isolates, the three matching-edge nodes, and selected
five-cover nodes have separate colors.  Its isomorphism group is exactly
`Aut(3K2) x S8`.  The final eight representatives were additionally checked
pairwise nonisomorphic by a separate colored-graph isomorphism implementation.
The independent exhaustion replay now in progress approaches the same target
set from the eight missing-five supports and their literal blocker, rather
than repeating this ridge recursion.

## Exact input and incidence audit

The consumer accepts either:

- the connected-ridge census object with
  `target_orbit_representatives`; or
- the frozen explicit-candidate object with `missing_A5s`.

For every supplied representative it requires exactly eight distinct
five-element supports, none containing one of the fixed missing edges

```text
(0,11), (1,12), (2,13).
```

It then independently scans all `2^14` subsets, reconstructs the minimal
transversal blocker, and requires the exact pure histogram

```text
5: 132.
```

Separately, it enumerates every nine-facet subset avoiding the three missing
edges and eight missing five-sets.  The resulting 132 literal actual vertices
must agree exactly with complements of the blocker.  Thus malformed or
off-schedule input is rejected before any Coxeter-label conclusion.

## Complete compact-label gate

The consumer regenerates the five classical order-five Lanner types and all
420 labelled placements over `{2,3,4,5}`.  It joins the eight supports in a
deterministic maximum-overlap order.  At every stage it tests every distinct
restriction induced by the 132 literal rank-nine vertices against a standalone
finite Coxeter classifier.  Compatible merged label states are deduplicated
exactly; stage counts and state-set digests are recorded.

Unassigned meeting labels are set to their minimal value `2`.  In normalized
form the vertex Gram matrix is `I-A` with `A` nonnegative.  Replacing an
unassigned `2` by `3`, `4`, or `5` only increases `A` entrywise, so a
non-positive-definite principal diagram cannot become positive definite.
Consequently a zero terminal state count is an unrestricted compact-label
no-go within the reduced alphabet `{2,3,4,5}`.  The dotted pairs occur in no
actual vertex and cannot alter this obstruction.

A nonzero terminal state count is deliberately reported only as a survivor.
It requires a later residual-label, cross-component, and exact Gram
rank/inertia analysis; the sparse consumer does not promote it to an example.

## Frozen-candidate calibration

On `data/enumerations/h9_14_record30_a0_candidate.json`, the new consumer
independently recovers:

```text
facet-complement SHA-256:
878d579b2ef3045c8948828bac4912450598fa44f1635b2743d8a27e58e0d8e6

literal actual vertices: 132
complete Lanner placements: 420
passing states: 420 -> 96 -> 6 -> 0
```

The four processed support masks are `211,242,496,8401`; the other four need
not be joined after the state set becomes empty.  A nontrivial relabelling by
`Aut(3K2) x S8` reproduces the same stage counts.

## Completed production decision path

The frozen production input and all-target label output are

```text
data/enumerations/h9_14_record30_a0_connected_ridge_census.json
data/enumerations/h9_14_record30_a0_all_target_labels.json
```

Their file SHA-256 digests are respectively

```text
3de6293068012c6213fcae03fb824d773cab6ead16c7820242482d5752cfbf0c
3bebba2d96461e2d46151f42c2258e1818267e9b6c7c13beeb1d30e1875e805b
```

The consumer records target-family digest
`9f8d03f5dbf2d6500a2b052a0fd45c82feb28641eeb697d4329063d4a5f691ac`
and summary digest
`3b6162151fcf6da220e063099da955ab2c0590f88682d643fbb710b6835ce9be`.
Replay with

```bash
python3 scripts/independent_h9_14_record30_a0_target_labels.py \
  /path/to/record30_a0_connected_ridge_census.json --progress \
  > /path/to/record30_a0_target_labels.json
```

The first branch has occurred: all eight representatives have zero states.
The remaining task is to compare this exact target set with the output of the
independent incidence enumeration.  Had any representative survived, it would
instead have required residual finite labels, connected cross-component shells,
and exact Gram rank/inertia work before any example claim.

## Scope

Everything here is compact-only.  No ideal-face, affine-cusp, or Vinberg
finite-volume stopping audit is performed.  No finite-volume noncompact
chamber is excluded, and no new hyperbolic Coxeter polytope is claimed.  Until
the independent incidence exhaustion agrees, this report records a major
production result rather than the final `a=0` theorem.

Focused tests are in
`tests/test_h9_14_record30_a0_target_labels.py`.
