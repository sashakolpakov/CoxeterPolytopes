# The 15-edge, three-pairing-two lane has no finite-volume survivor

**Status (2026-08-09): complete exact finite-denominator theorem.  All 111
Lorentzian weighted root cones and all 99 integral `E10` branches have
hyperideal extreme rays.  This is a scoped finite-volume exclusion, not an
abstract nonexistence theorem outside the lane and not a new polytope.**

The executable certificate is
[`certify_coble_nodot_15edge_three_pair2_lane.py`](../scripts/certify_coble_nodot_15edge_three_pair2_lane.py),
with [frozen class ledger](../data/enumerations/coble_nodot_15edge_three_pair2_lane.json)
and [tests](../tests/test_coble_nodot_15edge_three_pair2_lane.py).  The archive
SHA-256 is

```text
581480510c340cb77ad21f017aa48ace1e072d0f07ba35c6f54eb895a5b6be2c
```

## Exact scoped theorem

Take twelve proposed norm-minus-two roots.  Require their nonorthogonality
support to be connected, have exactly fifteen edges and maximum degree three,
mark exactly three distinct support edges with pairing two, and give each of
the other twelve support edges pairing one.  No surface, quotient determinant,
cusp count, known diagram, or volume condition is supplied as input.

For every resulting Cartan matrix of exact inertia `(9,1,2)`, quotient its
saturated rank-two radical.  Every resulting twelve-wall Coxeter cone in
`H^9` has a hyperideal extreme ray, hence infinite volume.  This holds for all
quotient cones before imposing an `E10` surface-lattice condition and again in
every even-unimodular `E10` marking.  No chamber reaches the full cusp lift or
split-K3/Coble same-`P` return.

The complete denominator is

```text
5,741 unlabelled connected supports
  x C(15,3) unordered marked edge triples
= 2,612,155 exact marked trials
-> rank 9: 3; rank 10: 208; rank 11: 27,509; rank 12: 2,584,435
-> rank-ten inertia (8,2,2): 84; target inertia (9,1,2): 124
-> 124 target markings form 111 weighted-isomorphism classes
-> all 111 quotient cones have hyperideal rays
-> 92 E10-saturable classes, with 99 integral branches
-> all 99 E10 branches have hyperideal rays
-> zero finite-volume survivors and zero surface-return attempts.
```

All 2,612,155 ranks are computed over the integers.  Exact rational
congruence then separates the 208 rank-ten matrices into the 124 desired
Lorentzian corank-two cases and 84 matrices with two negative directions.
There is no floating spectral prefilter in the certificate.  The 124 count is
for unordered edge-triple markings on unlabelled support representatives, not
vertex-labelled graphs.  Weighted isomorphism gives 100 singleton classes,
ten classes represented twice, and one represented four times.

## Lattice and volume ledgers

The marked target determinant histogram is

| determinant | markings | determinant | markings |
|---:|---:|---:|---:|
| `-1` | 59 | `-4` | 21 |
| `-8` | 3 | `-9` | 3 |
| `-13` | 8 | `-16` | 9 |
| `-28` | 1 | `-32` | 1 |
| `-40` | 1 | `-52` | 4 |
| `-56` | 1 | `-64` | 7 |
| `-72` | 3 | `-76` | 2 |
| `-80` | 1 |  |  |

Nineteen of the 111 weighted classes have nonsquare absolute determinant and
cannot possess a finite-index unimodular overlattice because
`|det L| = [U:L]^2 |det U|`.  Every square-determinant class has at least one
extension.  Complete discriminant quadratic-form enumeration gives no
extension for 19 classes, one for 85 classes, and two for seven classes: 99
integral `E10` branches.  Every branch is replayed as even, unimodular, of
signature `(1,9)`, containing all roots and reproducing their Gram.

Exact GMP double description is run independently on all 111 quotient cones
and all 99 integral markings.  All twelve walls are irredundant.  Every
quotient cone contains a certified strict timelike point, but also between 11
and 61 hyperideal extreme rays.  Each integral marking reproduces its
quotient cone's exact finite/ideal/hyperideal counts.  Full extreme-ray records
and hyperideal sublists are frozen by SHA-256 for every class and branch.

Literal Vinberg stops at the hyperideal endpoint in every case.  This is the
necessary distinction: ideal rays alone may be finite-volume cusps, while a
hyperideal extreme ray makes the acute chamber infinite-volume.  Since no
root cone passes that gate, attempting a cusp lift would not be a legitimate
completion and is correctly omitted.

## Strict timelike witnesses without a finite vertex

The corrected interior gate is essential in this shard.  The complete
equal-weight sum of extreme rays is a strict interior vector but is not
timelike in twelve classes.  Six of those have a finite extreme ray, which is
already a one-ray timelike seed.  The other six have no finite extreme ray at
all, but each has a pair of extreme rays whose sum has positive square.

For a timelike seed `v` and the complete ray sum `s`, the certificate doubles
`N` until

```text
x = N v + s
```

has positive square.  Because every wall slack of `s` is positive and every
ray slack is nonnegative, `x` remains strictly inside the cone.  The seed,
multiplier, coordinates, square, and strict-slack replay are archived exactly.
Ninety-nine quotient classes need `N=0`; the other twelve use powers of two
between 1 and 64.  Thus neither absence of a finite projective vertex nor a
spacelike equal-weight ray sum is mistaken for an empty hyperbolic chamber.

## Combined 15-edge criterion through three pairing-two edges

Combining the four frozen shards with pairing-two counts zero, one, two and
three gives an exact P-first classification for

```text
15 support edges, maximum degree 3, and 0 <= pairing-two count <= 3.
```

The combined finite denominator contains

```text
3,306,816 marked trials
238 target Lorentzian marked presentations
211 weighted root classes
180 integral E10 branches.
```

Exactly one of the 211 weighted root cones has finite volume: the
zero-pairing-two Type-II class.  Its two integral markings return the same
known Type-II full polytope and pass the frozen `raw_30` split-K3/Coble
compiler.  The other 210 weighted root cones—including every class with one,
two or three pairing-two edges—have hyperideal rays.  Equivalently, within
this declared union, finite volume holds exactly for the Type-II root diagram.

This combines the positive known return with three exact negative boundaries
into a criterion.  It does not assert that Type II is unique after allowing
four or more pairing-two edges, degree above three, a different edge count,
or a different root alphabet.

## Lineage and strict boundary

The zero-pairing-two shard isolates the known Type-II return.  The one-, two-
and three-pairing-two shards contribute respectively 28, 70, and 111
hyperideal weighted classes.  None of the present 111 root diagrams matches
the frozen Type-I or Type-II diagram.  This is lineage relative to the project
ledger only, not a claim that these Lorentzian diagrams have never appeared
in the literature.

The theorem excludes finite volume only in the exact stated lane.  It proves
neither abstract nonexistence of the diagrams as other geometric objects,
strict L4, a general dimension bound, nor the desired class-wide
`POLYTOPE -> SURFACE -> AG pushback` theorem.  No new compact or finite-volume
hyperbolic Coxeter polytope is claimed.

## Replay

```bash
PYTHONPATH=src:. python3 scripts/certify_coble_nodot_15edge_three_pair2_lane.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_coble_nodot_15edge_three_pair2_lane.py
```
