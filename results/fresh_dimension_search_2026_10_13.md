# Fresh constructions for a new dimension — 13 October 2026

**No new compact or finite-volume dimension was found.** This follow-up
changes quadratic spaces, ground fields, diagrams, and gluing constructions.
It does not reinterpret an unclosed chamber as a positive result.

The principal target here is a first **compact H9** example. This is a new
compact dimension, not a dimension with no known noncompact examples.
For arbitrary finite volume, the missing dimension H20 or a dimension
above H21 would be a different frontier. The published examples listed by
[Felikson's reference page](https://www.maths.dur.ac.uk/users/anna.felikson/Polytopes/polytopes.html),
checked on 13 October 2026, reach compact dimension eight and finite-volume
dimensions at most nineteen and twenty-one.

## Seven new arithmetic host searches

Four golden-field hosts `L_alpha=<-alpha>+E8+<1>` have determinant square
classes different from each other and from the previously searched hosts.
All have a complete candidate norm alphabet, and the searches enforce
reflection integrality and exact chronological Vinberg selection.

| Golden coefficient alpha | Height bound | Walls | Finite vertices | Hyperideal rays |
|---|---:|---:|---:|---:|
| `sqrt5` | 30 | 16 | 31 | 117 |
| `3phi` | 30 | 13 | 3 | 77 |
| `1+4phi` | 30 | 16 | 3 | 170 |
| `-1+3phi` | 30 | 14 | 24 | 54 |
| `sqrt5` | 100 | 43 | not enumerated | exact spacelike witness |

The fully enumerated rows have no cusps. A standalone rational-arithmetic
verifier validates all five saved escape certificates. The
[golden report](fresh_golden_neighbor_report.md) also proves that two other
tempting overlattice/binary modifications recover earlier groups; those
were eliminated before being counted as fresh searches.

The cubic family is over `K=Q(theta)`, where
`theta=2cos(2pi/7)` and `theta^3+theta^2-2theta-1=0`. Put `delta=2-theta`.
Its hosts are again `L_a=<-a>+E8+<1>`.

| Cubic coefficient a | Complete candidate norms | Height bound | Control levels | Walls | Verdict |
|---|---|---:|---:|---:|---|
| `theta` | `1,2` | 100 | 21 | 24 | exact escape |
| `2theta` | `1,2` | 100 | 6 | 17 | exact escape |
| `theta*delta` | `1,2,delta,2delta` | 100 | 31 | 14 | exact escape |

Each cubic host has signature `(9,1)` at the chosen embedding and is
positive definite at the other two. Since the positive block is
unimodular, primitive reflection norms satisfy the same local valuation
argument as the existing cubic E8-cube implementation. For `a=2theta`,
divisibility of every spatial coordinate by two forces
`q=-2theta*x^2 mod 4`, excluding a norm-four class. At the ramified prime
`delta` the analogous valuation is one. The Minkowski bound `14/9<2`
proves class number one; three independent unit signatures imply that
totally positive units are squares. Thus the candidate alphabets above
are complete modulo rescaling, rather than guessed norm cutoffs.

An exact dominant-vector census at height zero finds precisely the E8
and scalar A1 root orbits. At positive heights, exact conjugate bounds
enumerate the control coordinate; Schur-complement bounds enumerate E8
dominant weights; an exact square join supplies the last coordinate.
The nine initial walls and all accepted positive-height walls satisfy
every mutual nonpositive-pairing condition.

The three certified spacelike norms are, respectively,

```
(6318 - 3004 theta - 1520 theta^2)/1681 > 0
(1290 - 1841 theta + 784 theta^2)/169 > 0
-14 - 51 theta + 50 theta^2 > 0.
```

Every witness has time coordinate one, satisfies every original wall
inequality, and accompanies a strictly feasible timelike interior point.
The numerical linear program only suggests the support; exact cubic
arithmetic accepts or rejects the witness. A separate verifier uses only
Python `Fraction`, polynomial reduction, and rational root isolation.
It validates all saved roots, reflection coefficients, pairwise acuteness,
chronological level records and root-to-level matching, height bounds,
interiors, and escape inequalities. It also rejects forged norms, timelike
points disguised as escapes, false heights, reordered levels, and misplaced
roots. This checks the saved
certificates independently; it does not independently repeat enumeration
exhaustion.

Artifacts:

- [Cubic search](../scripts/search_fresh_cubic7_h9_20261013.py)
- [Independent cubic verifier](../scripts/audit_fresh_cubic7_h9_no_import_20261013.py)
- [Three H100 audit results](fresh_cubic7_h9_h100_independent_audit_20261013.json)
- [Independent mathematical scope audit](fresh_cubic7_h9_independent_scope_audit_20261013.md)
- [theta host](fresh_cubic7_h9_seed0_h100_20261013.json),
  [2theta host](fresh_cubic7_h9_seed1_h100_20261013.json),
  [theta delta host](fresh_cubic7_h9_seed2_h100_20261013.json)

**These are infinite-volume finite prefixes. The full reflection groups
remain undecided; later mirrors could remove the recorded escapes.**

## Geometric deformation and gluing

The literal known Bugaenko H8 facet admits no tilted extension: every
facet adjacent to it in a Coxeter H9 polytope must meet it orthogonally.
An orthogonal pair in its connected label-three component forces an
untilted wall, and label-three edges propagate this condition. The two
remaining label-five endpoints then rule out every integral tilt label
by exact inequalities. This is a local rigidity statement, not a proof
that a right extension or a compact H9 polytope cannot exist.
The [exact check](../scripts/fresh_bugaenko_facet_tilt_rigidity_20261013.py)
and [certificate](../data/enumerations/fresh_bugaenko_facet_tilt_rigidity_20261013.json)
record the calculation.

A separate numerical construction changes the original side diagram:
it varies finite labels and branch attachment positions, solves for a
dotted edge to give a rank-ten Lorentzian Gram matrix, and tests polar
truncations of the hyperideal vertices. This is a bounded construction
screen, not an exact nonexistence theorem. Its final report records the
finite alphabet and the degenerate and cusp-bearing cases left open.
The screen covers 47,185,920 labelled templates with original finite labels
in `{3,4,5,6}`. Its 113,099 Lorentzian skeletons split into 2,465 parked
degeneracies, 43,290 parked cusp-bearing cases, and 67,344 polar/old-wall
angle failures; there are no survivors. The new polar-angle test uses the
nearest integer label without imposing a finite maximum, at floating
tolerance `1e-7`. See the
[geometric report](fresh_bugaenko_gram_extension_20261013.md) and
[numerical archive](../data/enumerations/fresh_gram_polar_h9_complete_3456_20261013.json).
The angle gate assumes each proposed polar cut survives as a facet; it
does not first remove redundant polar cuts. This is another reason the
negative screen is not an exhaustive obstruction to compact completion.
The [subsequent polar-reduction audit](polar_reduction_h9_20261013.md)
now resolves this particular gap: after removing the invisible cuts,
all 67,344 tested noncusped systems still fail a retained angle test.
Its numerical and family-scope limitations remain.

The H18 experiment examines all relative twists of the common H17 facet
and all single-parent even doublings preserving that seam. There is only
one direct-twist class. Eleven parent-double classes produce finite-volume
H18 constructions, but none repairs the targeted lower-dimensional
forbidden angles. These derived examples do not establish a new dimension
or worldwide novelty. See the
[gluing report](fresh_gluing_h18_twists_doubles_20261013.md) for the exact
finite scope and reproducible artifacts.

## Reproduction and next distinctions

For example:

```sh
PYTHONPATH=src python3 scripts/search_fresh_cubic7_h9_20261013.py \
  --seed 0 --height 100 \
  --output results/fresh_cubic7_h9_seed0_h100_20261013.json
python3 scripts/audit_fresh_cubic7_h9_no_import_20261013.py \
  results/fresh_cubic7_h9_seed0_h100_20261013.json \
  results/fresh_cubic7_h9_seed1_h100_20261013.json \
  results/fresh_cubic7_h9_seed2_h100_20261013.json
PYTHONPATH=src:. python3 -m pytest -q tests/test_fresh_golden_neighbor.py
```

The golden tests pass all eight checks, including exact calibration against
the earlier odd-A1 ordered wall list and the two explicit lattice
identifications. The existing cubic field's number-theory and height-zero
regressions also pass (two selected tests).

Future work should preserve the distinction between these local failures
and global obstructions. It should avoid the two proved duplicate golden
modifications and the exhausted direct H18 twists. Removing redundant polar
cuts before applying angle gates was an unresolved step when this report
was written; it has since been completed in the
[next follow-up](polar_blocks_h20_search_2026_10_13.md).
Constructions outside the literal H8-facet family and different
positive root blocks are also open directions; no untested construction
is claimed here.
