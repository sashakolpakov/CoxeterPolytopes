# Record 0, `a=1`: terminal closure of source orbit 0

## Exact result

The last previously open `a=1` source orbit of Record 0 represented by

```text
A4 = {2,3,4,5}
```

has no compact Coxeter completion.  The conclusion is exhaustive for this
fixed source: the 108 classical labelled order-four Lanner placements split
under the full source action into nine disjoint orbits of sizes

```text
12, 12, 24, 12, 12, 12, 6, 12, 6.
```

For every orbit representative, the joint exact model simultaneously chooses
the 178 rank-nine facet complements and all finite Coxeter labels.  All nine
models terminate `INFEASIBLE`; no compact label assignment and hence no Gram
candidate survives.

This closes source orbit 0, not by extrapolating from a capped target sample,
but by covering all 108 possible labels on the unique missing `A4`.  A separate
ten-target orbit CEGAR was run only as a scout and is not used in the theorem.

## What is independently replayed

The frozen verifier reconstructs from the classical Lanner table:

- all 108 labelled order-four placements and their exact nine-orbit `S4`
  quotient, with partition digest
  `bc69a04c15d0d2126e687e28e5134f232fd2f65c1dc721285840a93c1c5830cc`;
- every labelled-`A4` stabilizer inside
  `S4({2,3,4,5}) x S7({6,...,12}) x C2({0,1})`;
- the complete joint-model counts
  `719 / 1607 / 1324 / 89 / 1617 / 6669 / 1456` for candidate facets,
  conditional `A5` encodings, face witnesses, meeting labels, independent
  five-sets, all-order Lanner closure, and ridges;
- ten deletion-minimal nonfinite Coxeter patterns, independently expanded to
  exactly 237,543 monotone face cuts with stream digest
  `033efa726adbea0daa6d7963fd1547653a5ab1c9c887fd064755fe774b8ce1ec`;
  and
- all nine complete no-timeout terminal documents, including the absence of
  a survivor and their solver statistics.

The independent replay authenticates the expensive terminal CP-SAT results;
it does not rerun roughly 36 minutes of branch searches during routine tests.
The aggregate certificate digest is

```text
fc8117901ca897ee21bc1a7becd7fcc505b4218ceb8f6a4eae55005ab58c5226.
```

## Schedule corollary and exact scope

The other three `a=1` sources were already frozen independently:

- source 1: all nineteen labelled-`A4` branches are compact-label impossible;
- source 2: the complete compact incidence model is infeasible; and
- source 3: the 41-core joint incidence/label model is terminally infeasible.

Together with the present source-0 theorem, these four certificates exhaust
the four source orbits.  Therefore **the entire compact Record-0 `a=1`
schedule is impossible**.

The quantifier stops there.  Record 0 as a whole is not closed: compact
schedules `a=0` and `a=2,...,16` remain open, while `a=17,18` were already
excluded by the strict summed-link theorem.  Consequently the global compact
frontier remains **134 wholly excluded graph shadows plus eight open shadows**
`0,1,3,4,5,7,9,31`.  This result neither proves global compact `H9/14`
nonexistence nor changes the ordinary-incidence partition `123+19`.

Nothing here addresses finite-volume noncompact chambers.  Ideal vertices,
parabolic subdiagrams, and Vinberg's endpoint criterion were not modeled.  No
new hyperbolic Coxeter polytope is constructed.

## Frozen artifacts and replay

- `scripts/certify_h9_14_record0_a1_source0_joint_branches.py`
- `data/enumerations/h9_14_record0_a1_source0_joint_branch_certificate.json`
- `tests/test_h9_14_record0_a1_source0_joint_branches.py`

Run

```bash
PYTHONPATH=src:. python3 \
  scripts/certify_h9_14_record0_a1_source0_joint_branches.py --check-archive
PYTHONPATH=src:. pytest -q \
  tests/test_h9_14_record0_a1_source0_joint_branches.py
```

