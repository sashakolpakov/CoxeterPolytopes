# Record 35: whole-shadow compact-Coxeter no-go

## Theorem and scope

Under the standard compact `H9/14` missing-face reductions, **record 35 has
exactly four residual incidence orbits passing the full necessary topology
audit, and none is the nerve of a compact hyperbolic Coxeter polytope**.
Equivalently, the whole record-35
missing-edge graph shadow is excluded on the compact track, not merely the
single ordinary-polytopal nerve selected in the earlier search.

The result is deliberately compact-only.  Compactness is used through strict
ellipticity of every rank-nine vertex diagram and strict positivity of every
vertex-link `g4`.  The proof does not enumerate ideal vertices, maximal
parabolic diagrams, or endpoint extensions, and it does not apply Vinberg's
finite-volume criterion. In particular, it neither constructs nor excludes a
noncompact chamber with this coarse shadow. Such a chamber could have a
different ideal-face incidence structure, and without the relevant parabolic
census it could have either finite or infinite volume.

The record-35 graph6 string is

```text
M?????????O?O?w??
```

with missing edges

```text
(0,11), (1,12), (0,13), (1,13), (2,13).
```

Its active vertices are `(0,1,2,11,12,13)`, and `(3,...,10)` are eight graph
isolates.  The active graph has automorphism group of order two; together
with the full `S8` action on the isolates, this is the complete graph
symmetry used below.

## All five compact schedules

Let `a` be the number of missing four-subsets and let `b` be the total number
of graph-independent missing five-subsets.  When `a>0`, `b` includes the
entire forced upper shadow of the missing A4s and is not the number of
additional minimal A5s.  Dehn--Sommerville and the local-link calculation
give the complete table

| `a` | `b` | rank-nine facets | sum of the fourteen vertex-link `g4` values |
|---:|---:|---:|---:|
| 0 | 5 | 102 | 30 |
| 1 | 11 | 100 | 24 |
| 2 | 17 | 98 | 18 |
| 3 | 23 | 96 | 12 |
| 4 | 29 | 94 | 6 |

For a compact Coxeter nerve every vertex link is spherical and satisfies the
strict inequality `g4 >= 1`.  The sum must therefore be at least fourteen,
which excludes `a=3,4`.  The exact completion census below treats all of
`a=0,1,2`.

## Positive missing-four branches

There are 490 individually admissible independent A4s after the exact
vertex-, edge-, and triple-link interval tests.
All vertex- and edge-link lower incidence bounds are exactly zero, while the
remaining singleton conditions are upper bounds.  Consequently an A4 that
fails this singleton filter cannot be repaired by adding a second A4; using
these 490 A4s as the ground set of the `a=2` census is lossless.

### The `a=1` incidence census

For a single A4, the complete forced-shadow gate is

| strict link gate | aggregate Macaulay gate | forced A5 shadow | additional minimal A5s | A4 count |
|:---:|:---:|---:|---:|---:|
| no | no | 8 | 3 | 196 |
| no | no | 9 | 2 | 112 |
| no | no | 10 | 1 | 70 |
| yes | yes | 8 | 3 | 56 |
| yes | yes | 9 | 2 | 56 |

The 112 viable A4s form exactly two orbits of size 56 under the active graph
automorphism and `S8` isolate actions.

For the first representative

```text
Q = (2,3,4,5)
```

the forced shadow has size nine, so two additional minimal A5s are needed.
There are 252 individually cap-feasible A5 candidates and
`C(252,2) = 31,626` unordered pairs.  Exactly 5,626 pairs satisfy all vertex
caps.  Only one has no transversal of order at most four, and it has the
scheduled 100 facets:

```text
R = (3,4,5,11,12)
B = (6,7,8,9,10).
```

Thus this branch supplies one residual incidence orbit.  Direct blocker
dualization gives 100 minimal transversals, all of order five.  The generated
complex has face vector

```text
(14,86,309,725,1165,1295,970,450,100),
```

all 450 ridges have degree two, the dual graph is connected, and all tested
links have spherical Euler characteristic and exact mod-two sphere homology.
It also passes the strict local `g` constraints, local Macaulay constraints,
and the global selection audit. These checks establish a survivor of the
necessary topology audit, not ordinary or Coxeter realizability.

For the second representative

```text
(2,3,4,11)
```

the forced shadow has size eight and three additional A5s are needed.  The
complete counts are 462 individually cap-feasible candidates,
`C(462,3) = 16,328,620` raw triples, 1,794,605 cap-feasible triples, and 840
triples without a small transversal.  None has the scheduled 100 facets, so
this orbit contributes no residual completion.

### The `a=2` incidence census

All `C(490,2) = 119,805` unordered distinct A4 pairs are considered.  Exact
family link intervals leave 111,699, with forced-shadow histogram

| forced A5 shadow | pairs |
|---:|---:|
| 15 | 224 |
| 16 | 24,864 |
| 17 | 42,756 |
| 18 | 30,520 |
| 19 | 11,480 |
| 20 | 1,855 |

Only 28 pairs pass the strict and aggregate shadow gate.  They form one
orbit, represented by

```text
(2,3,4,11), (2,3,4,12).
```

Its shadow has size fifteen, leaving two additional A5s.  There are 56
individually cap-feasible candidates, 1,540 raw pairs, and 540 pairs passing
all vertex caps.  Every one of the 540 has a transversal of order at most
four, which is forbidden in the compact `H9/14` selection model.  Hence the
entire `a=2` branch is empty before Coxeter labels are considered.

## Lossless `a=0` incidence census

Here the five missing A5s are minimal.  A potentially dangerous shortcut is
to assume that each A5 has a singleton active part.  No such normalization is
used. The production census allows **every independent active subset admitted
by the exact incidence caps in every row** and only then quotients the genuine
symmetries.

The five inclusion-minimal covers of the missing-edge graph are

```text
(0,1,2), (0,1,13), (1,11,13), (0,12,13), (11,12,13).
```

The lossless encoding and its checks are as follows.

1. All 18 independent subsets of the six active vertices that satisfy the
   exact per-row incidence caps are admitted as possible active parts of an
   A5 row. (The cap at vertex 13 is zero, so a row containing 13 is
   impossible.) Their size histogram is
   `0:1, 1:5, 2:8, 3:4`; in particular the four non-singleton triples
   `(0,1,2)`, `(0,2,12)`, `(1,2,11)`, and `(2,11,12)` are explicitly present
   and frozen in the certificate.
2. The five rows are taken as an unordered multiset.  Exact incidence caps
   and the small-transversal prefilter leave 1,221 active-row
   configurations.
3. For fixed active rows, each isolate is represented by its five-bit row
   incidence column.  The eight isolate columns are enumerated by their
   exact multiplicity vector.  This is a lossless quotient by `S8`, not a
   sampling or a canonical singleton ansatz.
4. Exact integer recursion solves 235 distinct column-distribution problems
   and produces 375,582 paths.  Of these, 7,239 have repeated A5 rows.  The
   prefilter is exact here: none of the remaining paths has a missed
   transversal of order at most four.  Thus 368,343 blocker-valid paths are
   retained.
5. Exactly 20 paths have the scheduled 102 facets.  Deduplication as literal
   unordered A5 families leaves 18 exact target families.

For every one of the 18 targets, blocker dualization returns exactly 102
minimal transversals, all of order five.  The complete topology and face-table
audit leaves precisely families `0,1,2,3,4` in deterministic production
order.  Each has scheduled face vector

```text
(14,86,309,726,1171,1309,986,459,102),
```

ridge degree two, connected dual graph, all required strict local and
Macaulay inequalities, spherical Euler characteristic in every face link,
and exact mod-two sphere homology for the complex and every face link.

Finally, the census exhausts both elements of the active graph automorphism
group, every isolate permutation through the exact column-multiset code, and
all `5!` row permutations.  The 18 target families form eleven orbits under
this full action; the five necessary-topology-audit survivors form exactly three.
Together with the unique `a=1` completion, record 35 therefore has exactly
four residual incidence orbits passing the full necessary topology audit.

The topology audit is used only as a necessary-condition filter.  It does not
infer polytopality or certify a PL sphere: exact mod-two homology does not
detect odd torsion or simple connectivity. One of the `a=0` orbits contains the previously proved
rational multiwedge polytope, but that separate ordinary-polytopality theorem
does not supply Coxeter angles.

## Unrestricted label obstruction for all `a=0` targets

Every order-five minimal nonface of a compact Coxeter polytope must carry an
order-five Lanner diagram.  The classical table has 420 labelled placements
on a fixed five-element support.  For each of the 18 target families, the
verifier joins the complete placement tables support by support and, after
every join, checks every restriction induced by the complete list of 102
actual rank-nine vertices.  No topology assumption or selected-nerve shortcut
is used in this label pass.

In the following exact transcript, `i:n` means that after support `i` is
joined, `n` joint sparse label states still make every current actual-vertex
restriction elliptic:

| family | successive support/state counts |
|---:|:---|
| 0 | `0:420, 1:96, 2:0` |
| 1 | `0:420, 1:96, 2:1,152, 3:2,304, 4:0` |
| 2 | `0:420, 1:96, 3:0` |
| 3 | `0:420, 1:96, 2:1,152, 3:2,304, 4:0` |
| 4 | `0:420, 3:96, 2:1,152, 1:2,304, 4:0` |
| 5 | `0:420, 2:96, 3:0` |
| 6 | `0:420, 2:96, 3:1,152, 1:0` |
| 7 | `0:420, 3:96, 4:1,152, 1:0` |
| 8 | `0:420, 3:1,152, 1:1,152, 2:0` |
| 9 | `0:420, 3:1,152, 4:1,152, 1:144, 2:0` |
| 10 | `0:420, 2:1,152, 1:1,152, 3:2,304, 4:0` |
| 11 | `0:420, 2:1,152, 1:1,152, 3:2,304, 4:0` |
| 12 | `0:420, 1:1,152, 2:0` |
| 13 | `0:420, 1:1,152, 2:0` |
| 14 | `0:420, 3:1,152, 4:1,152, 1:144, 2:0` |
| 15 | `0:420, 3:1,152, 4:0` |
| 16 | `0:420, 3:96, 4:0` |
| 17 | `0:420, 4:1,152, 3:0` |

Thus the complete sparse Lanner census fails for all 18 incidence families,
including all five necessary-topology-audit target families.

This is already an unrestricted finite-label obstruction.  In normalized
Gram convention write

```text
G = I - A,          A_ij = cos(pi/m_ij) >= 0.
```

The sparse census gives every unforced meeting pair label two, hence
`A_ij=0`.  Assigning any additional finite label `m>=3` increases `A`
entrywise.  Perron--Frobenius monotonicity cannot lower its largest eigenvalue
below one, so a non-elliptic principal actual-vertex diagram cannot become
elliptic.  The zero terminal counts therefore exclude arbitrary additional
non-right finite labels, not just the sparse ansatz.

## Unrestricted label obstruction for the `a=1` orbit

The sole residual `a=1` clutter has minimal nonfaces consisting of the five
missing edges and

```text
Q = (2,3,4,5)
R = (3,4,5,11,12)
B = (6,7,8,9,10).
```

Here `Q` and `R` overlap in three facets, while `B` is disjoint.  The complete
classical tables contain 108 labelled order-four Lanner placements on `Q`
and 420 labelled order-five placements on each of `R` and `B`.  Joining the
`Q` and `R` tables gives 1,542 compatible states.  Exact restrictions from
all 100 actual vertices leave 84 states on their six-facet union.  Combining
these with all 420 states on `B` yields 35,280 sparse block states.  All of
them pass the actual-vertex restrictions of sizes seven, eight, and nine,
whose distinct-size histogram is

```text
7:20, 8:55, 9:15.
```

There are thirty possible meeting pairs between the six-facet overlapping
cluster and `B`.  Testing one label-three edge in every sparse block gives

```text
35,280 * 30 = 1,058,400 trials,
```

of which exactly 5,040 pass all relevant actual-vertex tests.  Each of these
5,040 sparse base states admits exactly one passing edge; no base state admits
two.  Replacing that edge by label four or label five leaves zero survivors.
By monotonicity, label four also excludes every larger finite label.

This exhausts an arbitrary cross-edge set.  If a locally elliptic full
labeling had two or more cross edges, weaken every other unforced edge to
label two and either retained cross edge to label three.  Each retained edge
would have to occur among the single-edge survivors of the same sparse base
state, contradicting the exact one-edge-per-state census.  Hence at most one
cross edge can occur.

With no cross edge, `R` and `B` form an orthogonal sum of two Lanner diagrams,
so their ten-node principal Gram block has negative index two.  Every one of
the 5,040 exceptional one-edge states actually joins `R` to `B`; none uses
the `Q`-only facet 2.  Exact rooted weighted-graph classification gives

| ten-node block class | multiplicity | exact inertia |
|:---|---:|:---:|
| `A5_A5_type_1_1` | 4,320 | `(8,2,0)` |
| `A5_A5_type_4_1` | 720 | `(8,2,0)` |

For both representatives, exact leading principal minors have sign string
`++++-----+`, hence two sign changes and negative index two.  Cauchy
interlacing contradicts the single negative direction of an ambient
hyperbolic Gram matrix.  The remaining unforced pairs `(2,11)` and `(2,12)`
involve facet 2 outside the bad `R union B` principal block, so arbitrary
labels on them cannot alter its inertia.  Labels elsewhere likewise cannot
repair a bad principal submatrix.  This proves the unrestricted compact
Coxeter no-go for the fourth and last residual incidence orbit.

## Certificate and independent replay

The theorem package uses the following files:

- `src/coxeter_polytopes/combinatorics/h9_14_record35_all_completions.py`;
- `scripts/certify_h9_14_record35_all_completions.py`;
- `scripts/independent_h9_14_record35_all_completions_audit.py`;
- `data/enumerations/h9_14_record35_all_completions.json`; and
- `tests/test_h9_14_record35_all_completions.py`.

Reproduce it with

```sh
PYTHONPATH=src:. python3 scripts/certify_h9_14_record35_all_completions.py
PYTHONPATH=src:. python3 scripts/independent_h9_14_record35_all_completions_audit.py
PYTHONPATH=src:. python3 -m pytest -q tests/test_h9_14_record35_all_completions.py
```

The independent replay imports no project module. It reconstructs the five
schedules, the complete positive-branch gates, the all-cap-admissible-active-mask `a=0`
column-multiset recursion, exact blockers, all eighteen `a=0` Lanner joins,
the `a=1` unrestricted cross-edge gate, and the two exact inertia witnesses.
For the five `a=0` ridge-pseudomanifold candidates and the `a=1` completion it
also independently regenerates the full face table and checks the scheduled
face vector, exact minimal nonfaces, dual connectivity, all-link Euler
characteristics, exact mod-two homology of the complex and every face link,
strict local `g` conditions, Macaulay bounds, and missing-four intervals.
Thus the exact four-orbit topology refinement and the all-completion compact
label no-go are both replayed independently.

Frozen artifact values:

```text
certificate digest:          4ebb69e2415ff97c37d62bd8347a2e44b19288cd2d3afa3893b9d0bd763227f8
certificate archive SHA-256: 4b903a7010bfab088b5fd5f0e3fb80e809e0450602fa329d6ba76b58325e8c47
hostile-summary digest:      b01b632d01569c847194a8fb592eca4e91b0cc9865e00d296d029147b3e25332
```
