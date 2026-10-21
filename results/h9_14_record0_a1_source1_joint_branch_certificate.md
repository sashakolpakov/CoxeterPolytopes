# H9/14 Record 0, `a=1`, source orbit 1: joint branch certificate

## Certified conclusion

For the fixed Record-0, `a=1` source represented by the missing `A4`
`{0,2,3,4}`, all 19 exact labelled-`A4` placement-orbit branches terminate
`INFEASIBLE` in the joint compact incidence/finite-Coxeter-label model.  No
branch contains a compact-label survivor.  Since these 19 disjoint orbits
contain all 108 literal order-four Lanner placements, this closes this fixed
source for compact Coxeter polytopes.

This is a fixed-source conclusion.  A terminal checkpoint closes only its own
labelled-`A4` orbit; the fixed-source inference is made only after verifying
all branch indices `0,...,18`, terminal status `INFEASIBLE`, no survivor, and
orbit-size sum `108`.  It does not close the other Record-0 `a=1` sources or
the whole Record-0 graph shadow.

The frozen machine-readable certificate is
`data/enumerations/h9_14_record0_a1_source1_joint_branch_certificate.json`.
Its internal canonical SHA-256 is

```text
da685ee5fcfbe0ee6ba1e8951b5c115dd73cb6ca6e3bccebf2f4487c24f05777
```

## Independent offline replay

The stdlib-only aggregate verifier
`scripts/certify_h9_14_record0_a1_source1_joint_branches.py` does not invoke
CP-SAT and does not repeat an expensive branch search.  Starting from the
classical nine order-four Lanner diagrams, it independently generates all
108 labelled placements.  The fixed source acts as

```text
S3({2,3,4}) x S8({5,...,12}),  order 6 * 40320 = 241920.
```

The `S3` action on the `A4` gives 19 disjoint placement orbits, with sizes
`[6,6,6,3,6,6,6,6,6,6,6,6,6,6,6,6,6,6,3]`.  Their exact partition digest is

```text
df9a46605bdf4899a5994ef4e8ce9a8002f9c1cc0d832341c547eefba433f938
```

For each representative, the verifier enumerates the exact internal code
stabilizer.  Seven adjacent transpositions generate the outside `S8`; the two
size-three orbits have one additional internal involution.  Hence the exact
labelled-`A4` stabilizer order is `40320` for a size-six orbit and `80640` for
a size-three orbit.  The verifier records every generator as a full
14-element action, checks that it preserves the source, the code, and the
754-element facet ground, and independently recovers the checkpoint lex-step
counts (`332` for each outside generator and `402` for the internal
involution).

The fixed source also independently reconstructs every frozen model-shape
count:

| object | exact count |
|---|---:|
| candidate rank-nine facet complements | 754 |
| independent five-sets | 1,617 |
| conditional minimal-`A5` encodings | 1,608 |
| global meeting-label variables | 89 |
| face-witness constraints | 1,324 |
| all-order Lanner-closure constraints | 6,796 |
| ridge zero-or-two constraints | 1,491 |

## Nonfinite-pattern and cut certificate

The earlier three-round seed contains 73 concrete nonfinite face cuts.  The
aggregate verifier canonicalizes them afresh under all diagram permutations,
deduplicates across rounds, and checks every result is nonfinite while every
one-vertex deletion is finite.  This yields exactly 23 abstract types, with
order histogram

```text
order 3: 5; order 4: 9; order 5: 8; order 6: 1.
```

The canonical pattern-list digest is

```text
f9767cc4cd2bd873de238b315c0bac6c5d5509bc73549967da55c7f101bb0971
```

Each type is then embedded into every eligible labelled face, with every
label placement in its full isomorphism orbit, and its monotone threshold cut
is streamed into SHA-256 in the production order.

| seed round | new abstract types | expanded cuts | cut-stream SHA-256 |
|---:|---:|---:|---|
| 0 | 16 | 1,707,843 | `6ff340050b5485917e88e150e1f67432a42aa529616cfea85848943ede442408` |
| 1 | 6 | 51,410 | `fe00fefeaa6b4b0bfd41a7aa25ae3b3b917b147c2e6abc3f10ddfe10de666169` |
| 2 | 1 | 10,548 | `6b97c1c26c882d15fb236143061924145141f1ab4f05490875349fad6dee8d90` |
| total | 23 | 1,769,801 | `a4c01774988940813a748b744e585fc22086fdbadc8a5796e90078f0d76262eb` |

The digest of the three-row expansion transcript is
`f3758234e404cb557a89d9caebf8e66e3a2c2065cd4f831a2a6e656076d45010`.
Every branch checkpoint is required to contain this exact pattern library,
transcript, and installed-cut total, with zero later cut rounds.

## Terminal branch ledger

Every row below has terminal status `INFEASIBLE`, survivor `null`, and branch
closure `true`.  The checkpoint-level fixed-source flag is deliberately
`false`; only the aggregate sets it to `true` after the complete partition is
accounted for.

| branch | representative code | orbit | stabilizer order | wall seconds | solver branches | conflicts |
|---:|---|---:|---:|---:|---:|---:|
| 0 | `[2, 2, 3, 2, 3, 5]` | 6 | 40320 | 124.910116 | 69665 | 10736 |
| 1 | `[2, 2, 3, 3, 2, 5]` | 6 | 40320 | 173.121919 | 61693 | 7086 |
| 2 | `[2, 2, 4, 5, 2, 3]` | 6 | 40320 | 153.505226 | 41474 | 1698 |
| 3 | `[2, 2, 5, 2, 3, 3]` | 3 | 80640 | 214.106131 | 78142 | 10597 |
| 4 | `[2, 2, 5, 4, 2, 3]` | 6 | 40320 | 168.394044 | 37842 | 818 |
| 5 | `[2, 2, 5, 5, 2, 3]` | 6 | 40320 | 152.450472 | 35652 | 386 |
| 6 | `[2, 3, 3, 3, 4, 2]` | 6 | 40320 | 202.255903 | 45655 | 1897 |
| 7 | `[2, 3, 3, 3, 5, 2]` | 6 | 40320 | 155.209818 | 55942 | 5061 |
| 8 | `[2, 3, 4, 3, 3, 2]` | 6 | 40320 | 144.580402 | 60284 | 5034 |
| 9 | `[2, 3, 4, 4, 3, 2]` | 6 | 40320 | 148.100956 | 36166 | 330 |
| 10 | `[2, 3, 4, 5, 2, 2]` | 6 | 40320 | 115.503665 | 39694 | 1188 |
| 11 | `[2, 3, 4, 5, 3, 2]` | 6 | 40320 | 114.182845 | 36361 | 501 |
| 12 | `[2, 3, 5, 2, 3, 2]` | 6 | 40320 | 159.605517 | 123123 | 26826 |
| 13 | `[2, 3, 5, 3, 3, 2]` | 6 | 40320 | 117.448416 | 65079 | 7118 |
| 14 | `[2, 3, 5, 4, 2, 2]` | 6 | 40320 | 99.077759 | 37039 | 593 |
| 15 | `[2, 3, 5, 4, 3, 2]` | 6 | 40320 | 97.702771 | 36866 | 444 |
| 16 | `[2, 3, 5, 5, 2, 2]` | 6 | 40320 | 207.038989 | 34770 | 195 |
| 17 | `[2, 3, 5, 5, 3, 2]` | 6 | 40320 | 180.250502 | 36218 | 510 |
| 18 | `[3, 3, 5, 2, 2, 2]` | 3 | 80640 | 184.439536 | 43025 | 653 |

The aggregate verifier also accepts a legacy branch-0 terminal artifact with
the three optional solver-statistics fields absent, reports those statistics
as unavailable, and still requires its explicit terminal status to be
`INFEASIBLE`.  Other branches must contain complete terminal statistics.

## Reproduction and scope

Run the frozen verification and focused tests with

```bash
PYTHONPATH=src:. python3 scripts/certify_h9_14_record0_a1_source1_joint_branches.py
PYTHONPATH=src:. pytest -q tests/test_h9_14_record0_a1_source1_joint_branches.py
```

The verification is an independent semantic replay of the partition,
stabilizers, model counts, pattern library, cut expansion, and every terminal
checkpoint field.  It authenticates the archived CP-SAT terminal results; it
does not rerun those searches.

No finite-volume noncompact conclusion is made.  Ideal vertices, affine
cusps, maximal parabolic subdiagrams, and noncompact Vinberg endpoints are not
part of this certificate.
