# Record 18 is impossible in the compact H9/14 problem

## Theorem and exact scope

There is no compact hyperbolic Coxeter 9-polytope with fourteen facets whose
missing-edge graph is record 18,

```text
M???????????o?K?_
```

with missing edges

```text
(0,12), (1,12), (2,13), (3,13), (12,13).
```

This is a whole-shadow theorem: it disposes of every completion allowed by the
compact H9/14 missing-face and Lanner reductions, not only the previously
selected ordinary-polytopal nerve.

It is strictly a compact theorem.  Ideal vertices, parabolic diagrams, and
Vinberg's finite-volume endpoint criterion are not part of the computation.
Consequently it neither constructs a new hyperbolic Coxeter polytope nor
excludes a finite-volume noncompact record-18 chamber.

## Complete branch partition

Let `a` be the number of minimal missing four-faces.  Dehn--Sommerville and
the Lanner-order bound give the complete compact schedules

| `a` | missing five-faces | facets | sum of vertex-link `g4` |
|---:|---:|---:|---:|
| 0 | 8 | 110 | 42 |
| 1 | 14 | 108 | 36 |
| 2 | 20 | 106 | 30 |
| 3 | 26 | 104 | 24 |
| 4 | 32 | 102 | 18 |
| 5 | 38 | 100 | 12 |
| 6 | 44 | 98 | 6 |

The Murai--Nevo equality case, together with the Lanner bound, makes every
compact vertex-link `g4` strictly positive.  Thus `a=5,6` are impossible
because their summed values are below fourteen; `a>=7` already has nonpositive
global `g4`.  The remaining branches close as follows.

| branch | exhaustive input | exact incidence targets | final obstruction |
|---:|---:|---:|---|
| `a=4` | 13 positive-A4 orbits | 0 | incidence/blocker |
| `a=3` | 40 positive-A4 orbits | 0 | incidence/blocker |
| `a=2` | 14 positive-A4 orbits | 8 literal targets | all 8 fail compact Coxeter labels |
| `a=1` | 3 positive-A4 orbits | 1,080 literal targets in 8 orbits | all 8 orbits fail compact Coxeter labels |
| `a=0` | 196 literal targets in 3 orbits | 196 | all 3 orbits fail compact Coxeter labels |

This covers every nonnegative value of `a`.

## Independent enumeration and literal equality

The positive-A4 incidence-matrix replay independently reproduces every one of
the 70 production orbit keys and literal orbit sizes: 3, 14, 40, and 13 in
the `a=1,2,3,4` branches.

The difficult `a=0` search is reproduced by a structurally different exact
enumerator.  Both implementations return the same set of 196 literal
eight-row families, whose explicit incidence normal forms give three classes
of sizes 144, 36, and 16.

For `a=1`, a collision-free byte-matrix production search and a standalone
128-bit quotient independently produce exactly the same 1,080 literal target
families.  For all twelve `a=2/add4` representatives, the same comparison is
made representative by representative: ten have no target, `[53,198]` has
six, and `[53,202]` has two.  The two `a=2/add5` representatives have no
incidence target.  Equality here is equality of normalized families, not an
equality of counts or probabilistic fingerprints.

## Topology and Coxeter-label gates

Every target sent to the label stage is blocker-dualized directly.  Its face
table, ridge degrees, connected dual graph, local `g` and Macaulay conditions,
and the mod-two homology of every face link are checked independently.

The compact label search regenerates the classical labelled Lanner diagrams,
joins their placements on all minimal nonfaces, and tests every actual
rank-nine vertex restriction against the finite Coxeter classification.
Unassigned labels are set to 2.  This is lossless: for normalized Gram
`G=I-A`, increasing a Coxeter label increases the nonnegative matrix `A`, so
a restriction already failing ellipticity cannot be repaired later.

The terminal label counts are zero in all three `a=0` classes, all eight
`a=1` classes, and all eight literal `a=2` targets.  Hence none of the
topologically valid incidence targets is a compact Coxeter nerve.

## Frozen certificate

The aggregate verifier checks the complete schedule partition, invokes the
frozen positive-branch, `a=3`, and literal-family verifiers, checks every
downstream zero terminal, and pins every evidence archive and source file.

```console
PYTHONPATH=.:src python3 \
  scripts/certify_h9_14_record18_whole_shadow_compact_no_go.py
PYTHONPATH=.:src pytest -q \
  tests/test_h9_14_record18_whole_shadow_compact_no_go.py
```

Principal artifacts:

- `scripts/certify_h9_14_record18_whole_shadow_compact_no_go.py`
- `data/enumerations/h9_14_record18_whole_shadow_compact_no_go.json`
- `data/enumerations/h9_14_record18_a0_independent_audit.json`
- `data/enumerations/h9_14_record18_a1_ridge_raw_families.json`
- `data/enumerations/h9_14_record18_a2_add4_ridge_raw_families.json`
- `data/enumerations/h9_14_record18_a3_independent_audit.json`

The certified compact-facing H9/14 ledger therefore advances from `131 + 11`
to `132 + 10`.  The ordinary-incidence partition remains `123 + 19` because
record 18 still has an explicit ordinary-polytopal completion; ordinary
polytopality is not hyperbolic Coxeter realizability.
