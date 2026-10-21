# Record 30, `a=0`: final hostile audit of the dual-CNF proof

## Verdict

**PASS**, with compact-only scope.  The dual CNF is a lossless necessary
model for every genuine compact Record-30 `a=0` target after the canonical
seed, its 63,960 terminal blockers are exactly the seed-containing images of
the eight production orbits, and the archived binary DRAT proof verifies
against the archived CNF.  The separate label replay has zero states on all
eight orbits.  No fatal or material proof-chain gap remains.

This proves a compact schedule no-go.  It says nothing about noncompact
finite-volume chambers and supplies no new hyperbolic Coxeter polytope.

The machine-readable final audit is
[`h9_14_record30_a0_dual_hostile_audit.json`](../data/enumerations/h9_14_record30_a0_dual_hostile_audit.json).

## Inclusion of every genuine compact target

The three fixed missing edges form `3K2`.  At `a=0`, Dehn--Sommerville forces
eight independent missing five-sets and 132 facets.  Facet complements are
therefore five-element transversals of the eleven minimal nonfaces.

The 326 possible five-covers of `3K2` split as

```text
3 active + 2 isolates: 224
4 active + 1 isolate:   96
5 active + 0 isolates:   6.
```

Since only 102 covers have the last two types, every 132-facet target has a
type-three cover.  `Aut(3K2) x S8` is transitive on the 224 type-three covers,
so fixing `{0,1,2,3,4}` is lossless.

The remaining CNF gates are necessary and were checked literally:

- exact unary counters select eight missing `A5`s and 132 blocker elements;
- no size-three/four transversal can occur, since its complement would be a
  face with more than nine vertices;
- each blocker variable is equivalent to meeting all eight selected `A5`s;
- each nonempty six-set ridge group has degree exactly two; and
- the local-cap arrays are the strict-link bounds
  `[7,3]`, `[7,3,1]`, and `[6,3,1,0]` in face sizes one, two, and three.

For the last item, the top link-(g) values before the containing-`A5`
incidence is subtracted are `[8,4]`, `[8,4,2]`, and `[7,4,2,1]`.
Strict positivity gives exactly the encoded caps.  The optional face-vector
encoding was also audited: all eight targets have exactly 56 absent
independent six-sets and 160 absent independent seven-sets.

Sequential-counter, blocker-biconditional, ridge, missing-face-indicator, and
lex encodings passed exhaustive small truth tables.  The lex tests include
the aliased-variable case in which the right vector is a permutation of the
left.  Symmetry breaking is lossless because the characteristic-vector
maximum in each finite seed-stabilizer orbit dominates every generator image.

## Structurally independent orbit and clause replay

The no-project-import replay
[`hostile_h9_14_record30_a0_dual_cnf_orbits.py`](../scripts/hostile_h9_14_record30_a0_dual_cnf_orbits.py)
closes each full family orbit by elementary permutation generators, rather
than importing or repeating the generator's `seed_forms` routine.  Its source
SHA-256 is

```text
f49d28b445adffa22b326582dc8b4bc7d6018e3699c093d2140dc21ca81ca3b7
```

The eight `(full orbit, stabilizer, type-three blocker, seed form)` tuples are

```text
( 2240, 864, 120,  1200)
(26880,  72, 123, 14760)
(26880,  72, 122, 14640)
( 6720, 288, 120,  3600)
( 6720, 288, 124,  3720)
(26880,  72, 124, 14880)
( 6720, 288, 124,  3720)
(13440, 144, 124,  7440).
```

For every row, direct seed filtering agrees with the independent incidence
double count

```text
seed forms = full orbit * type-three blockers / 224.
```

The orbits are pairwise disjoint.  Their 63,960 seed families have SHA-256
`591e61e6cb33cbd7fb039cb6313336fa022c0ce33ad15143a81cc94778f5a25e`;
the corresponding negative DIMACS clauses have the deliberately separate
SHA-256 `76c2b85f93cdc562af26bb9cee0841b872a0b94b24ed86d3926556c0d52d56d2`.
Exact set comparison—not hash comparison
alone—shows that these clauses are precisely the CNF tail.  The frozen replay
object is
[`h9_14_record30_a0_dual_orbit_hostile.json`](../data/enumerations/h9_14_record30_a0_dual_orbit_hostile.json),
file SHA-256

```text
b33237d52eb8206cad7c8c6c5bc95e43c7ec54dac5367030a0c55fc9803122b5.
```

## Independent proof replay and provenance

The archived compressed files pass zstd integrity and decompress to

```text
CNF raw SHA-256:
531eedff4e8ba07f3750351d3416f9834abf6fbfd6e0d338fedbb7cc463f511c

DRAT raw SHA-256:
9e374009efedec7b3d6c11e2c004cca9f95e3bea3a724c1f94e4686a13802a47.
```

An independent full replay used the exact checker binary recorded by the
certificate:

```text
drat-trim h9_record30_a0_dual_lex.cnf \
  h9_record30_a0_dual_lex_kissat.drat -i
```

It returned exit code zero and `s VERIFIED` after 488.306 seconds.  This is a
new hostile-rerun timing, not the archived run's 589.079 seconds.  All semantic
counters agree exactly:

```text
127,835 / 984,935 clauses in core
106,647 / 934,371 lemmas in core
34,109,825 resolution steps
21,068 RAT lemmas
94,524 redundant literals.
```

The checker is `drat-trim` commit
`2e3b2dc0ecf938addbd779d42877b6ed69d9a985`, binary SHA-256
`f58f63b0f76945d4c4c9ff6e87afaf870f579e67c0f7cca589492df8fc7ebd47`.
The proof-producing solver is Kissat 4.0.4, binary SHA-256
`7458a5dae42d7dcc0c8402bb153fe9faee321cb6f3d079b79b0fb385e55d8bb6`.
The deep certifier now checks both full binary hashes before
replaying, as well as the decompressed CNF and proof hashes.

## Corrected findings and verification

Three audit findings were resolved before this PASS:

1. The first report draft accidentally displayed the `a=1` face vector.  It
   now displays the certified `a=0` vector
   `(14,88,328,806,1364,1596,1248,594,132)`.
2. Family and DIMACS-clause hashes now have distinct field names.
3. Deep replay now rejects solver/checker binaries that do not match the
   frozen provenance.

Focused verification:

```bash
python3 scripts/certify_h9_14_record30_a0_dual_exhaustion.py
python3 -m pytest -q \
  tests/test_h9_14_record30_a0_dual_cnf.py \
  tests/test_h9_14_record30_a0_dual_orbit_hostile.py \
  tests/test_h9_14_record30_a0_deep_tool_provenance.py
```

Result: `16 passed`.  The hostile orbit/truth-table file contributes eight of
those tests.  The main certificate's internal SHA-256 is
`436b9bc834529f7d2d94754e43c3f785bb937ec0f7157570ad9a06d56858ab2d`.
