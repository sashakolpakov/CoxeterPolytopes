# Exact `a=0` census for H9/14 record 34

## Result

The `a=0` branch of record 34 is impossible.  An exact bicoloured-incidence
census finds only 18 covering states with the required 114 facets.  They are
all pure blockers of order five, but every one has ridges of degree one and
three (and one family also has degree-four ridges).  None is the boundary of
a simplicial polytope.

This conclusion uses no SAT/SMT result.  The isolated Z3 `UNSAT` run is only
corroboration and is not a step of this proof.

## The row-column encoding

The missing-edge graph is

\[
K_3(0,11,13)\sqcup K_2(1,12)\sqcup9K_1.
\]

Dehn--Sommerville requires six missing five-faces and 114 facets.  Exact
local-link intervals say that a missing five-face contains at most one of the
five active vertices; each `K3` vertex occurs in at most one row, each `K2`
vertex in at most two, and an isolate in at most five rows.  Thus the six A5s
are a `6 x 9` isolate-incidence matrix whose row degrees are four (one active
label) or five (no active label).

Write an active profile as `(a,b0,b1)`, where `a<=3` is the number of used
`K3` labels and `0<=b0<=b1<=2` are the multiplicities of the two `K2`
labels.  There are 23 formal profiles.  `(2,2,2)` has an immediate
four-transversal: omit the unused `K3` label and choose one `K2` label in a
minimal graph cover, then add the other `K2` label.  The census therefore
checks the remaining 22 profiles.

For each of the six minimal graph covers `C` (two `K3` vertices and one `K2`
vertex), a transversal of size at most four exists exactly when `C` already
hits all six A5s or the rows missed by `C` have a common vertex.  This is the
exact `no_small` test in the enumerator.

## Exhaustive quotient

`genbg` generates one representative of every unlabelled bipartite graph
with six row vertices and nine isolate vertices.  The command

```sh
genbg -q -d4:0 -D5:5 6 9 24:30
```

imposes row degree at least four and both-part degree at most five.  The C++
checker rejects rows outside degrees four/five, validates bipartiteness,
payload padding, and isolate degree, then enumerates every compatible active
colouring.  Forgetting active colours can only identify states; restoring all
colourings may overcount but cannot omit a state.  The 18 target outputs
reduce to nine orbits under
`S3(K3) x S2(K2) x S9(isolates) x S6(rows)`.

| active rows | `genbg` reps | colourings | distinct rows | no `<=4` transversal | 114 facets | ridge-valid |
|---:|---:|---:|---:|---:|---:|---:|
| 6 | 25,707 | 1,542,420 | 1,510,900 | 256 | 0 | 0 |
| 5 | 133,404 | 8,671,260 | 8,491,233 | 211,677 | 0 | 0 |
| 4 | 308,755 | 9,571,405 | 9,227,208 | 936,139 | 9 | 0 |
| 3 | 393,034 | 5,109,442 | 4,771,932 | 1,179,505 | 9 | 0 |
| 2 | 289,910 | 1,449,550 | 1,289,534 | 612,265 | 0 | 0 |
| 1 | 116,784 | 233,568 | 195,042 | 146,100 | 0 | 0 |
| 0 | 20,559 | 20,559 | 15,852 | 15,852 | 0 | 0 |
| **total** | **1,288,153** | **26,598,204** | **25,501,701** | **3,101,794** | **18** | **0** |

The target-state profile counts are `(2,0,1):6`, `(2,1,1):6`,
`(3,0,0):3`, and `(3,0,1):3`.  Their nine canonical orbits have the
following positive ridge-degree histograms:

| orbit count | degree histogram | degree-two dual components |
|---:|---|---:|
| 4 | `{1:40, 2:433, 3:40}` | 3 |
| 1 | `{1:39, 2:438, 3:33, 4:3}` | 2 |
| 1 | `{1:39, 2:435, 3:39}` | 2 |
| 1 | `{1:38, 2:437, 3:38}` | 2 |
| 1 | `{1:36, 2:441, 3:36}` | 2 |
| 1 | `{1:32, 2:449, 3:32}` | 2 |

For a six-set `S`, its complementary eight-set is a ridge precisely when at
least one of the six five-sets `S-{v}` is a transversal; its ridge degree is
the number of those transversals.  A simplicial polytope boundary requires
that number to be exactly two.  The table gives an explicit obstruction for
every target orbit.

## Reproduction and independent replay

Compile and run the complete census with

```sh
clang++ -O3 -std=c++20 scripts/certify_h9_14_record34_a0.cpp -o /tmp/record34_a0
genbg -q -d4:0 -D5:5 6 9 24:30 | /tmp/record34_a0
```

The archived manifest is
`data/enumerations/h9_14_record34_a0_census.json`.  It contains all 18 raw
target outputs, their nine orbit identifiers, source `graph6` records,
colouring masks, and ridge audits.  The Python verifier independently
reconstructs all minimal transversals, the target face vector, ridge
histograms, degree-two dual components, and symmetry quotient.

- C++ source SHA-256:
  `35a8975b2f9687604af51432324063607eaff043016cbd3346e13add232aa683`
- manifest certificate SHA-256:
  `bf878a6940b324a1a142716aca6e03f3a50758f593cf4fbca092a803e21e69b7`

