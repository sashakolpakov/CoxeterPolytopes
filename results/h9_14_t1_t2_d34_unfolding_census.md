# TOOT: a two-diagram surface-free `H9/14` unfolding census

Status: **exact bounded-corpus `IC0/R0` census PASS** (2026-08-08).

## Result

Freeze the primary certified exact-Gram corpus consisting of the five rows in
the `H9/14` factory, the parked rank-ten stock, and the direct parked-product
audits.  It has six `H9/14` markings and four nonisomorphic Coxeter diagrams:
`T1`, `T2`, `T3`, and the one-step parked doubling product.  Define `D34` by
surface-free conditions only:

- dimension 9, fourteen walls, Gram rank/inertia `(10;1,9,4)`;
- ordinary Coxeter relations in `{2,3,4}` with exactly eight `m4` edges; and
- the minimal involutive lift alphabet: a fixed `(-2)` root or an orthogonal
  exchanged pair of `(-2)` roots.

Inside this pinned corpus, `D34` contains exactly three markings and two
nonisomorphic diagrams: `T1` and `T2` (`C01` and its duplicate `C14`).  This is
a complete denominator only inside the stated corpus, not an enumeration of
all abstract `H9/14` diagrams.

The exact compiler no longer assumes that ten facet roots form an integral
basis.  It constructs the quotient root lattice from all fourteen roots by
Hermite normal form.  This change is essential for `T1`.

| datum | scalings | compatible-lane candidates | survivor | rank-18 unfolding | K3 complement |
|---|---:|---:|---|---|---|
| `T2`, `75m2+8m3+8m4` | 18 | `16+5` | unique `U+E8(-2)` | `(18,4,0)`, determinant `-16` | `U(2)+U(2)` |
| `T1`, `74m2+9m3+8m4` | 18 | `2+5` | unique `U+E8(-2)` | **new `(18,2,0)`, determinant `-4`** | `U+U(2)` |

Thus all 36 minimal crystallographic scalings and all 28 reflection-stable
even fixed-overlattice candidates are exhausted.  Each diagram has exactly
one survivor, and both force the same fixed lattice `U+E8(-2)`.

For `T1`, the compiler constructs 24 acute square-minus-two roots in four
fixed and ten exchanged-pair orbits.  Their off-diagonal pairings are
`242*0+34*1`, and their orbit sums reproduce the selected labelled `T1` Gram
exactly.  The roots span index two in the rank-18 lattice; this is recorded,
not silently promoted to index one.  Exact graph glue embeds the `(18,2,0)`
lattice primitively in the K3 lattice with complement `U+U(2)` and extends the
involution.  This is a new weak K3 lattice/action carrier lineage for the
known `T1` polytope, not a new polytope and not yet a nef-cone theorem.

## Exact boundary

- `U`: open.  Corpus-relative arithmetic exhaustiveness is not universality
  for all `H9/14` realizations and does not construct the required intrinsic
  algebraic surface locus.
- `R`: exact labelled fold-back and lattice/action recovery (`R0`) is closed
  for both diagrams.  Full geometric marked recovery is already available
  only for `T2` through its independent marked bridge.
- `C`: the existing `T2` promotion supplies complete invariant-nef equality;
  `T1` cone completeness and hidden-wall exclusion remain open.
- `V`: casewise finite volume is imported only after exact labelled diagram
  identity with the existing archives.  `T1` has `112+15` and `T2` has
  `90+20` finite/ideal endpoints, with zero bad elliptic edges and zero bad
  parabolic neighbourhoods.  Finite volume is not inferred from the weak K3
  compiler.

## Replay

```bash
PYTHONPATH=.:src python3 scripts/certify_h9_14_t1_t2_d34_unfolding_census.py --check-archive
PYTHONPATH=.:src pytest -q tests/test_h9_14_t1_t2_d34_unfolding_census.py
```

Artifacts:

- [`certify_h9_14_t1_t2_d34_unfolding_census.py`](../scripts/certify_h9_14_t1_t2_d34_unfolding_census.py)
- [`h9_14_t1_t2_d34_unfolding_census.json`](../data/enumerations/h9_14_t1_t2_d34_unfolding_census.json)
- [`test_h9_14_t1_t2_d34_unfolding_census.py`](../tests/test_h9_14_t1_t2_d34_unfolding_census.py)
