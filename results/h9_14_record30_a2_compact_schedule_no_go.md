# Record 30: complete compact `a=2` schedule no-go

## Result

The compact Dehn--Sommerville schedule `a=2` for Record 30 is impossible.
This is a schedule-level theorem, not merely a selected-nerve computation.

The exact positive-`A4` quotient has 21 orbits.  Two independent incidence
enumerators exhaust their additional minimal-`A5` choices and obtain exactly
77 symmetry orbits with the scheduled pure blocker of 128 five-element facet
complements.  A complete order-four/order-five Lanner join on every actual
rank-nine vertex kills 74 targets.  Exactly three survive that sparse gate;
all three have already or newly certified unrestricted whole-nerve no-gos.

No hyperbolic Coxeter polytope has been found.  This result is compact-only:
it neither performs Vinberg's ideal-vertex stopping test nor excludes a
finite-volume noncompact chamber.

## Complete incidence quotient

The fixed missing graph is `3K2 + 8K1`, with dotted pairs

```text
(0,11), (1,12), (2,13).
```

For every one of the 21 upstream positive-`A4` representatives, the
production program forms a colored incidence graph whose colors distinguish
facet vertices, the three unordered dotted-edge rows, the two unordered
`A4` rows, and the unordered added-`A5` rows.  Nauty canonicalizes every
partial family.

The independent C++ replay does not use nauty.  It explicitly checks all
1,935,360 automorphisms of `3K2 + 8K1`, filters the stabilizer of each
unordered `A4` pair, and uses a complete colored-incidence canonical form:
it enumerates the `3!`, `2!`, and `d!` row orders and sorts the fourteen
facet-column signatures.  Thus its equivalence relation is exactly the one
encoded by the production coloring.

Both implementations agree in every orbit on stabilizer order, search-node
count, canonical-state counts, pruning counts, blocker histogram, and the
complete set of target-orbit canonical keys.

| positive-`A4` orbit | stabilizer | added `A5`s | search nodes | pure target orbits |
|---:|---:|---:|---:|---:|
| 0 | 3,456 | 1 | 57 | 1 |
| 1 | 55,296 | 0 | 1 | 1 |
| 2 | 1,152 | 2 | 1,058 | 3 |
| 3 | 384 | 2 | 525 | 0 |
| 4 | 144 | 2 | 782 | 0 |
| 5 | 384 | 2 | 2,393 | 1 |
| 6 | 2,880 | 3 | 2,990 | 11 |
| 7 | 11,520 | 2 | 576 | 2 |
| 8 | 768 | 2 | 791 | 0 |
| 9 | 384 | 2 | 935 | 0 |
| 10 | 1,152 | 2 | 1,173 | 1 |
| 11 | 480 | 3 | 2,915 | 4 |
| 12 | 96 | 3 | 2,833 | 1 |
| 13 | 144 | 3 | 3,551 | 0 |
| 14 | 2,880 | 5 | 369,607 | 37 |
| 15 | 2,880 | 4 | 9,069 | 3 |
| 16 | 240 | 4 | 23,953 | 3 |
| 17 | 192 | 4 | 21,822 | 0 |
| 18 | 11,520 | 4 | 120,626 | 1 |
| 19 | 960 | 4 | 561,409 | 5 |
| 20 | 768 | 4 | 777,951 | 3 |
| **total** |  |  |  | **77** |

Orbit 14 is the largest five-addition branch: its raw ground set has

```text
C(966,5) = 6,937,477,486,758
```

families.  Both quotient implementations reduce it to 369,607 visited search
nodes and 37 targets.  The independent literal stabilizer expansion gives
4,290 fixed-`A4` target families, with orbit-size histogram

```text
60:17, 90:1, 120:14, 180:1, 240:1, 360:3.
```

## Complete sparse label gate

The no-project replay regenerates all nine order-four and five order-five
classical Lanner types, hence all 108 and 420 labelled placements.  For each
of the 77 incidence targets it:

1. reconstructs the literal 128-element minimal blocker;
2. checks blocker/actual-vertex complement duality;
3. joins every Lanner placement on the two `A4`s and all added `A5`s; and
4. requires every induced actual rank-nine vertex diagram to be elliptic.

Every unassigned meeting label is set to two during this join.  If a principal
vertex diagram already fails, increasing an omitted finite label to 3, 4, or
5 cannot restore positive definiteness, by Perron--Frobenius monotonicity.
The resulting zero-state conclusions are therefore unrestricted within the
target-specific finite alphabet `{2,3,4,5}`.

The terminal-state histogram is

```text
0 states:      74 targets
648 states:     2 targets
11,664 states:  1 target
```

The three survivors are exactly:

- orbit 1, target 0: `A4` masks `120,1920`, with no added `A5`;
- orbit 7, target 1: `A4` masks `57,2104`, added `A5` masks `248,1984`;
- orbit 10, target 0: `A4` masks `57,2496`, added `A5` masks `1592,1984`.

The first is the frozen selected Record-30 product nerve.  Its complete
support-normalized cross-edge gate has no admissible single `m=3` edge, and
stronger labels cannot repair it.  The selected-nerve theorem therefore
excludes this entire target orbit.

## Exact whole-nerve metric closure of the two 648-state exceptions

Orbit 7 target 1 and orbit 10 target 0 have the same decisive census:

```text
actual rank-nine vertices:       128
sparse labelled states:          648
nerve automorphism order:        1,152
sparse-state orbits:             6
orbit sizes:                     72,144,144,72,144,72
connected shells per state:      56,128,38 at sizes 2,3,4
first empty connected shell:     size 5
```

All residual internal meeting pairs are forced right-angled: adding any one
label 3, 4, or 5 already violates an actual-vertex ellipticity condition.
The minimal-nonface graph has three indefinite components, so a Lorentzian
Gram graph with negative index one must connect them.  Every connected
assignment has a two-edge spanning-tree restriction; downward monotonicity
therefore makes the shell growth through the first empty size lossless.

For each nerve, every one of

```text
6 * (56 + 128 + 38) = 1,332
```

quotient shell assignments has an eleven-node dotted-free principal block on
facets `3,...,13`.  The replay computes its determinant exactly in the
integral basis

```text
1, sqrt(2), sqrt(5), sqrt(10)
```

by a division-free signed subset recurrence.  It records the four integer
coefficients, the exact integer field norm, and the real sign certified by
rational radical enclosures.  In both nerves all 1,332 determinants are
strictly positive and every exact field norm is nonzero.  Hence every literal
assignment has Gram rank at least eleven by symmetry transport, whereas a
hyperbolic chamber in `H^9` has Gram rank at most ten.

The independent modulus-1129 residues are all nonzero as well, but they are
only redundant replay checks.  The theorem rests on the materialized exact
integer field norms and rigorously certified determinant signs, not modular
nonsingularity.

## Exact scope

This closes Record 30's compact `a=2` schedule.  Together with the earlier
certificates, it closed `a=2,3,4,5,6,7`; the subsequent independently
replayed `a=1` theorem leaves only `a=0` open.  Therefore this is still not a
whole-shadow Record-30 theorem.

Nothing here classifies ideal vertices, affine cusp diagrams, or finite-volume
noncompact chambers.  Such a claim still requires Vinberg's literal
finite-volume stopping criterion.  Nothing here supplies a Lorentzian Gram
realization, so it is not a new hyperbolic example.

## Frozen artifacts

- `scripts/certify_h9_14_record30_a2_schedule_no_go.py`
- `scripts/explore_h9_14_record30_a2_orbit14_nauty.cpp`
- `scripts/independent_h9_14_record30_a2_all_orbits_explicit.cpp`
- `scripts/independent_h9_14_record30_a2_orbit14_labels.py`
- `scripts/independent_h9_14_record30_a2_all_target_labels.py`
- `scripts/independent_h9_14_record30_a2_exceptional_whole_nerves.py`
- `data/enumerations/h9_14_record30_a2_schedule_no_go.json`

Frozen certificate SHA-256:

```text
3ac5e4ea4ddd624aef17dec45f64d403165b736bab015d088e2b8d78f43a5adf
```

Frozen archive-file SHA-256:

```text
06759a184708d9c2b8a4ab15fc772fc2a4b7e66ac275d2ff537bebf92cb9ebe6
```
