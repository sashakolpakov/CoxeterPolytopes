# Record 0, `a=1`, source 2: compact incidence no-go

## Fixed-source theorem

For the Record-0 missing-edge graph

```text
M?????????????o??
```

the `a=1` source orbit represented by

```text
A4 = {0,1,2,3}
```

cannot be the nerve of a compact hyperbolic Coxeter 9-polytope with fourteen
facets.

This is a theorem for one fixed source orbit only.  The other `a=1` source
orbits `0,1,3` remain under search, so neither the complete `a=1` schedule nor
the whole Record-0 shadow is excluded.  Nothing here excludes finite-volume
noncompact chambers.

## Exact dual model

A rank-nine vertex of the simplicial nerve has a five-element complement.
The model therefore selects the scheduled 178 facet complements rather than
guessing the 23 additional missing five-sets directly.  It imposes:

- every intended set through order four has a selected disjoint facet;
- exactly 32 independent five-sets are missing;
- every six-set has selected-five degree zero or two; and
- every nonface of order at least six contains a fixed missing edge, the
  fixed `A4`, or a missing independent five-set.

The final condition is the complete compact Lanner-order closure: compact
minimal hyperbolic Coxeter diagrams have order at most five.  It rules out
hidden minimal blockers of orders six and above.

Both exact encodings independently reconstruct the same finite model:

| item | count |
|---|---:|
| candidate facet complements | 754 |
| independent five-sets | 1,617 |
| required-face constraints | 1,324 |
| all-order closure constraints | 6,796 |
| ridge degree constraints | 1,491 |
| source-stabilizer generators | 10 |

The fixed-source stabilizer is

```text
S2({2,3}) x S9({4,...,12}) x S2({0,1}), order 1,451,520.
```

The generator moved-pair counts are

```text
201, 166, 166, 166, 166, 166, 166, 166, 166, 165.
```

The lex slice is lossless: every finite group orbit has a global
lexicographic minimum, and that minimum is no larger than each generator
image.

## Independent terminal results

- Production OR-Tools CP-SAT, one worker: `INFEASIBLE`.
- Independently written Z3 reconstruction: `unsat`.

Thus there is no necessary compact incidence target for this fixed source;
the label and Gram stages are never reached.

## Frozen artifacts

- `scripts/certify_h9_14_record0_a1_source2_incidence_no_go.py`
- `scripts/explore_h9_14_record0_a1_dual_facets.py`
- `scripts/independent_h9_14_record0_a1_dual_z3.py`
- `data/enumerations/h9_14_record0_a1_source2_incidence_no_go.json`
- `tests/test_h9_14_record0_a1_source2_incidence_no_go.py`

The frozen internal certificate SHA-256 is
`9c30e26f9466338adac2acf7aa1fc033cdeee201520002baa350bdfb1b362ab5`.

Reproduce and verify with

```console
PYTHONPATH=src:. caffeinate -i python3 \
  scripts/certify_h9_14_record0_a1_source2_incidence_no_go.py --fresh \
  --write data/enumerations/h9_14_record0_a1_source2_incidence_no_go.json
PYTHONPATH=src:. python3 \
  scripts/certify_h9_14_record0_a1_source2_incidence_no_go.py
PYTHONPATH=src:. python3 -m pytest -q \
  tests/test_h9_14_record0_a1_source2_incidence_no_go.py
```
