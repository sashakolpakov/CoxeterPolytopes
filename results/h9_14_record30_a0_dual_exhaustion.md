# Record 30, `a=0`: independently certified compact no-go

## Theorem certified

There is no compact hyperbolic Coxeter polytope in dimension nine with
fourteen facets, Record-30 face vector

```text
(14,88,328,806,1364,1596,1248,594,132),
```

and no minimal missing four-set (`a=0`).  The production connected-ridge
census has exactly eight candidate incidence orbits, and the complete compact
Coxeter-label gate kills all eight.  A structurally different dual SAT
exhaustion, followed by an independently checked DRAT proof, now certifies
that the production list is complete.

This is a compact-only theorem.  It does not test Vinberg's ideal-face
finite-volume criterion, does not exclude noncompact finite-volume chambers,
and does not produce a new hyperbolic Coxeter polytope.

The frozen theorem object is
[`h9_14_record30_a0_dual_exhaustion.json`](../data/enumerations/h9_14_record30_a0_dual_exhaustion.json).
Its internal certificate SHA-256 is

```text
436b9bc834529f7d2d94754e43c3f785bb937ec0f7157570ad9a06d56858ab2d
```

## Production result being replayed

The connected-ridge production search visited `28,039,413` nodes, accepted
`19` labelled paths, and reduced them to eight isomorphism orbits.  The
standalone complete label consumer gives zero terminal states on every orbit:

```text
0 states: 8 targets.
```

The label computation uses every classical order-five Lanner placement in
the reduced finite alphabet `{2,3,4,5}` and tests ellipticity on all 132 actual
rank-nine vertices.  Burcroff's global multiple-edge bound (Lemma 10.10) makes this
alphabet complete for the compact problem; Perron--Frobenius monotonicity
shows that increasing an unspecified label cannot repair a failed principal
vertex diagram.

The production artifacts are

- [`h9_14_record30_a0_connected_ridge_census.json`](../data/enumerations/h9_14_record30_a0_connected_ridge_census.json), containing the eight incidence orbits; and
- [`h9_14_record30_a0_all_target_labels.json`](../data/enumerations/h9_14_record30_a0_all_target_labels.json), containing their zero-state label replays.

## Independent dual exhaustion

The replay does not repeat the production ridge-growth algorithm or its
nauty quotient.  It instead chooses the eight independent five-element
minimal nonfaces directly.  From them it defines their literal blocker among
the 326 five-subsets meeting all three fixed missing edges of `3K2`.  A
selected blocker element is the complement of an actual rank-nine vertex.

The CNF imposes the following compact necessary incidence conditions:

1. exactly eight independent missing five-sets;
2. the frozen 430 vertex/edge/triple link caps;
3. no transversal of size three or four;
4. exact logical equivalence between each blocker variable and meeting all
   eight missing five-sets;
5. exactly 132 blocker elements; and
6. degree zero or two in every one of the 833 dual-ridge groups.

Connectivity and the remaining face-vector entries are deliberately omitted.
The CNF therefore describes a larger class than the production census.  Its
UNSAT result is stronger than the exhaustion statement actually needed.

### Lossless seed and symmetry breaking

The 326 possible blocker elements split by active matching vertices as

```text
3 active + 2 isolates: 224
4 active + 1 isolate:   96
5 active + 0 isolates:   6.
```

Since 132 elements are selected but only 102 can have the last two types,
every target has a blocker element of the first type.  The group
`Aut(3K2) x S8` is transitive on its 224 elements, so fixing
`{0,1,2,3,4}` loses no orbit.

The stabilizer of this seed has order 8,640.  Eight explicit generators are
used for lexicographic symmetry breaking.  This is lossless: the
characteristic-vector lexicographic maximum in each finite stabilizer orbit
is at least each of its generator images and hence satisfies every encoded
comparison.

For each of the eight production targets, every full-group image whose
blocker contains the fixed seed was reconstructed and blocked literally.
There are 63,960 distinct blocked seed forms, distributed as

```text
1,200; 14,760; 14,640; 3,600; 3,720; 14,880; 3,720; 7,440.
```

Thus an UNSAT result means that no further necessary-incidence orbit exists
beyond the eight supplied production orbits.

## Hostile encoder audit

The frozen verifier rebuilds the acting group without importing the CNF
generator's group or seed-form helpers.  It scans all permutations of the six
active vertices preserving `3K2`, traverses `S8` literally on the isolates,
and reconstructs the full target orbits, automorphism orders, seed forms, and
all blocking clauses.  In particular, it verifies the orbit--stabilizer and
seed double-count identities target by target.

The encoder was also attacked from the positive side.  Removing exactly the
63,960 final blocking clauses leaves 920,975 clauses.  For each target, the
audit forces a lex-admissible seed image and obtains a SAT assignment.  Each
full 226,358-variable assignment was then checked literally against all
920,983 clauses and decoded back to exactly eight missing five-sets and 132
blocker elements.  Hence none of the eight expected targets is accidentally
removed by the base encoding or the lex constraints.

Small exhaustive truth tables separately verify the unary exact/at-most
counters and Boolean lex encoder.  The finalized hostile-audit summary SHA is

```text
bcfa508f5a616754481e09bdc7a84d16b97245b67d1b751a4deb7b2199b337d0
```

## Proof certificate

The exact frozen formula has

```text
variables: 226,358
clauses:   984,935
raw CNF SHA-256:
531eedff4e8ba07f3750351d3416f9834abf6fbfd6e0d338fedbb7cc463f511c
```

Kissat 4.0.4 returned `UNSATISFIABLE` in 57.19 seconds after 38,437
conflicts and emitted a binary DRAT proof.  `drat-trim` at commit
`2e3b2dc0ecf938addbd779d42877b6ed69d9a985` independently returned
`s VERIFIED`; the check used 34,109,825 resolution steps and took 589.079
seconds.

Frozen proof artifacts:

- [`h9_14_record30_a0_dual_exhaustion.cnf.zst`](../data/enumerations/h9_14_record30_a0_dual_exhaustion.cnf.zst), compressed SHA-256 `368bcf7f5eef1d418bfd26601a4b46d2979d197e0188e7c07d6d6f8573075be3`;
- [`h9_14_record30_a0_dual_exhaustion.drat.zst`](../data/enumerations/h9_14_record30_a0_dual_exhaustion.drat.zst), compressed SHA-256 `480c043221dd5a6485b18ec6ac9c6d1aaab3b811873b1948c5daec374f291159`; and
- [`h9_14_record30_a0_dual_exhaustion_drat_trim.txt`](../data/enumerations/h9_14_record30_a0_dual_exhaustion_drat_trim.txt), normalized checker transcript.

## Reproduction and exact scope

The independent sources are

- [`independent_h9_14_record30_a0_dual_cnf.py`](../scripts/independent_h9_14_record30_a0_dual_cnf.py);
- [`verify_h9_14_record30_a0_dual_cnf.py`](../scripts/verify_h9_14_record30_a0_dual_cnf.py);
- [`certify_h9_14_record30_a0_dual_exhaustion.py`](../scripts/certify_h9_14_record30_a0_dual_exhaustion.py); and
- [`test_h9_14_record30_a0_dual_cnf.py`](../tests/test_h9_14_record30_a0_dual_cnf.py).

Run the frozen fast verification with

```bash
python3 scripts/certify_h9_14_record30_a0_dual_exhaustion.py
python3 -m pytest -q tests/test_h9_14_record30_a0_dual_cnf.py
```

The certifier also accepts `--deep --solver PATH --drat-trim PATH` to
decompress the frozen files, confirm their raw digests, rerun the hostile
audit, and invoke the proof checker.  Deep replay first requires the supplied
Kissat and `drat-trim` executables to have the exact binary SHA-256 digests
archived in the theorem object; it therefore cannot silently substitute a
different solver or checker while claiming to reproduce the frozen run.

This closes the entire compact Record-30 `a=0` schedule.  It does not by
itself assert a whole-shadow Record-30 theorem; that requires the separate
aggregation of the already certified `a=1,...,7` schedules.  No statement in
this certificate concerns ideal vertices, affine cusps, or finite-volume
noncompact polytopes.
