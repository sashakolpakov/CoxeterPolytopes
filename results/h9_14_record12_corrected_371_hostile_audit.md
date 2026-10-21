# Corrected record-12 371-orbit hostile audit

## Result

The corrected record-12 census has been replayed independently from the
frozen order-ten graph stream through the final 371 quotient representatives.
The replay imports no production record-12 census or residual module.  Its
three independent stages agree exactly with the production streams:

| stage | exact result |
|---|---:|
| frozen triangle-free `X` graphs | 12,172 |
| locally feasible `(X, cut, y)` rows | 1,391,382 |
| rows surviving the sound special-blocker obstruction | 6,184 |
| coloured branch orbits | 1,308 |
| ordered compatible branch-type pairs | 22,858 |
| distinct embedded `(X,Z)` rows tested | 2,725,719 |
| rows passing the complete coupled blocker schedule | 42,000 |
| full symmetry orbits | 371 |

The 371-orbit histogram by the number `a` of rank-four minimal nonfaces is

```text
a = 0   48
a = 1   95
a = 2  144
a = 3   63
a = 4   21
```

All 371 representatives pass the complete incidence and link audit described
below.  This validates the corrected **combinatorial census**.  It does not
assert that any representative admits a hyperbolic Coxeter Gram matrix, nor
that it is compact or finite-volume as a hyperbolic Coxeter polytope.

## The blocker correction

For a clutter `C`, write `b(C)` for its inclusion-minimal transversals.  Two
edge cases and one one-sided obstruction are essential.

First,

```text
b(empty clutter) = { empty row }.
```

Every set is a transversal of the empty clutter, and the empty set is its
unique inclusion-minimal transversal.  Initializing an incremental blocker by
the empty list incorrectly deletes this branch.

Second, suppose a base triple clutter `T` is enlarged by the rows
`{s,t,i}` for `i` in `Y`.  From every row `r` of `b(T)`:

- if `r` meets `{s,t}`, retain `r`;
- otherwise form `r union Y`, `r union {s}`, and `r union {t}`;
- inclusion-minimalize all candidates globally.

This is exact, including `T` empty and `Y` empty.  If neither `s` nor `t` is
chosen, every element of `Y` must be chosen; otherwise choosing either special
vertex hits the entire added family.

Third, let the two coupled complement clutters have the form

```text
P = { A union t : t in T_P },
Q = { B union u : u in T_Q }.
```

If `r` is in `b(T_P)` and meets `B`, then `r` is already a globally minimal
transversal of `P union Q`: it hits `P` because it hits every `t`, it hits all
of `Q` through `B`, and the minimality witnesses for `r` in `T_P` show that
deleting any vertex makes it miss a member of `P`.  It is not one of the four
fixed `A`-by-`B` rank-two rows.  Consequently a valid completion requires
`|r|` to be four or five.  The symmetric statement holds with `P` and `Q`
interchanged.

This proves the special-row prefilter used here.  It does **not** justify
rejecting blocker rows that avoid both special vertices: such rows can change
after coupling and global minimalization.  The older isolate-row filter was
therefore not used.

The local zero-or-two ridge equation was also solved directly.  Its feasible
`y` support lies in one side of the cut, and a selected vertex has no
`X`-neighbour inside that side.  There is no forced-`z=1` case.  The paired
replay checked this independently on all `1,308 * 45 = 58,860` retained
state-edge instances.

## Independent replay stages

### 1. Branch generation and coloured quotient

[`independent_h9_14_record12_371_branches.py`](../scripts/independent_h9_14_record12_371_branches.py)
regenerates all normalized switched cuts and locally feasible `y` values from
the 12,172 frozen `X` graphs, computes blockers incrementally, applies only the
sound special-row lemma, directly recomputes every retained full local blocker,
and performs an exact colour-preserving graph-isomorphism quotient.

It obtains:

```text
locally feasible rows                         1,391,382
sound special-row survivors                       6,184
coloured branch orbits                            1,308
sum of coloured-orbit multiplicities              6,184
empty base clutters                                   1
blocker of that empty clutter                 { empty }
```

Every retained blocker agrees with a direct recomputation.  The semantic
branch-summary digest is
`de2ce0116ee97e2b9f958f6e3b616803744b0d96817f221db341931b2dfb3e3c`.

### 2. Both cut orientations and complete coupled blockers

[`independent_h9_14_record12_371_pairs.cpp`](../scripts/independent_h9_14_record12_371_pairs.cpp)
uses an independent incremental blocker and an MRV colour-preserving subgraph
embedding search.  It explicitly runs both the direct and complementary cut
orientations and deduplicates their common images before applying the complete
coupled blocker test.

```text
compatible type pairs                         22,858
direct-profile pairs                          13,516
complement-profile pairs                       9,598
pairs allowing both orientations                 256
direct embedded images                     1,381,203
complement embedded images                 1,376,516
images shared by both orientations             32,000
distinct images tested                     2,725,719
complete-blocker passes                        42,000
```

The hostile and production 42,000-row streams are equal as sets.  Their
order-independent semantic digest is
`a189faa3d1faeb1f4db21630d2fbb816f8d3d4f1421d04eb251112d37945f336`.

### 3. Full quotient and topology

[`independent_h9_14_record12_371_quotient.cpp`](../scripts/independent_h9_14_record12_371_quotient.cpp)
canonicalizes each pair of 5-by-5 binary layers under:

- independent layer complementation and layer swap (`D8`);
- `S5 x S5` permutations within the normalized balanced cut;
- interchange of the two cut sides by matrix transpose.

The last two operations are the stabilizer of a normalized balanced cut, so
this is the normalized form of the full residual-vertex relabelling action.
The replay proves exact equality, including multiplicities, with the archived
371 rows and separately verifies that every archived row is canonical.

[`independent_h9_14_record12_371_orbits.py`](../scripts/independent_h9_14_record12_371_orbits.py)
then reconstructs every facet-complement clutter and recomputes its minimal
nonfaces.  For each of all 371 representatives it verifies:

- minimal-nonface orders are exactly two, four, and five;
- the rank-two rows are exactly the four prescribed graph edges;
- `b = 12 - 2a`, and the number of facets is `130 - 2a`;
- the complete scheduled face vector;
- every ridge has degree two;
- the facet-dual graph is connected;
- the strict local `g` and Macaulay inequalities used by the census;
- the reduced Euler characteristic of the link of every face;
- the complete mod-two homology of the link of every face.

The face-vector schedules independently recovered are:

| `a` | `(f0,...,f8)` |
|---:|---|
| 0 | `(14,87,320,780,1320,1554,1224,585,130)` |
| 1 | `(14,87,320,779,1314,1540,1208,576,128)` |
| 2 | `(14,87,320,778,1308,1526,1192,567,126)` |
| 3 | `(14,87,320,777,1302,1512,1176,558,124)` |
| 4 | `(14,87,320,776,1296,1498,1160,549,122)` |

All 371 selected-complement digests are distinct.  The topology summary's
semantic digest is
`6a58b02d39324916d24de5d43a8b1091393f83cadc7963de8581023266c66c84`.

The mod-two link calculation is a homology-manifold audit, not by itself a
certificate that every link is a PL sphere.  No claim stronger than the
checks listed above is made.

## Pinned replay streams

| stream | SHA-256 |
|---|---|
| 1,391,382 local rows | `ed12edfbff40c2fab730ed3f1740bcf6e5272554f85a807468c7cb8e515961ca` |
| 6,184 corrected special-row survivors | `c91021fe69bc047c8f2eea5f43bf708d9025492c8e3cac3552205b857731e4ee` |
| 1,308 coloured orbit rows | `640217770888c33fe60fc2d1ffd6dbb05c499c9307c40d4ec9642a482fca1905` |
| production 42,000 pair rows | `6b86880f60359a6ae638acf55a730282f2b46e5cfceb74c45d8ae7f0035d635a` |
| hostile 42,000 pair rows | `f7c093c6b4ecb3ce5a29ca1989f680ae69ffda7a73dc54de8c1d324740fc29c1` |
| archived 371 orbit rows | `4d978594b6be052328b4056d6b05c6249e80e8b156a09d2eba8057f3b7daa23c` |
| hostile topology JSON | `3119558a2e170ef2b2419b7049409e800d4ad3177fa2538db65ba8fad95145fc` |

The raw pair-file hashes differ because the independent generator emits a
different row order; exact semantic set equality is checked before quotienting.
