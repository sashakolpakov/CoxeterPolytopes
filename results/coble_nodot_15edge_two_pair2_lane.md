# The 15-edge, two-pairing-two lane has no finite-volume survivor

**Status (2026-08-09): complete exact finite-denominator theorem.  All 70
Lorentzian weighted root cones have hyperideal extreme rays.  This is a
negative boundary result, not a new compact or finite-volume example and not
the full nodot theorem.**

The executable certificate is
[`certify_coble_nodot_15edge_two_pair2_lane.py`](../scripts/certify_coble_nodot_15edge_two_pair2_lane.py),
with [frozen class ledger](../data/enumerations/coble_nodot_15edge_two_pair2_lane.json)
and [tests](../tests/test_coble_nodot_15edge_two_pair2_lane.py).  The archive
SHA-256 is

```text
85a2a347b84954b20e3f796958a244755fd994a76508fbfa3bad84f077c15167
```

## Exact scoped theorem

Take twelve proposed norm-minus-two roots.  Require their nonorthogonality
support to be connected, have exactly fifteen edges and maximum degree three,
mark exactly two distinct support edges with pairing two, and give every other
support edge pairing one.  No surface, determinant, cusp count, known diagram,
or volume condition is supplied as input.

For every resulting Cartan matrix of exact inertia `(9,1,2)`, quotient its
saturated rank-two radical.  Every resulting twelve-wall Coxeter cone in
`H^9` has a hyperideal extreme ray, hence infinite volume.  This is true for
all quotient cones before an `E10` condition is imposed and for every
even-unimodular `E10` marking afterward.  Therefore this lane contains no
compact or finite-volume root chamber and cannot enter the one-boundary
split-K3/Coble cusp compiler.

The complete denominator is

```text
5,741 unlabelled connected supports
  x C(15,2) unordered marked edge pairs
= 602,805 exact marked trials
-> 595,360 rank 12 + 7,362 rank 11 + 83 rank 10
-> all 83 rank-ten trials have exact inertia (9,1,2)
-> 70 weighted-isomorphism classes
-> all 70 quotient cones have hyperideal rays
-> 54 E10-saturable classes, with 57 integral saturation branches
-> all 57 E10 branches have hyperideal rays
-> zero finite-volume survivors and zero surface-return attempts.
```

Every one of the 602,805 ranks is computed over the integers.  Exact rational
congruence supplies the inertia of the rank-ten cases, with no floating
spectral prefilter.  The 83 count concerns unordered edge-pair markings on
unlabelled support representatives, not vertex-labelled graphs.  Their
weighted-isomorphism multiplicities are 60 singleton classes, eight classes
with two markings, one with three, and one with four.

## Two independent obstruction layers

The quotient determinants of the 83 marked survivors are

| determinant | marked trials |
|---:|---:|
| `-1` | 51 |
| `-4` | 7 |
| `-5` | 3 |
| `-12` | 1 |
| `-13` | 2 |
| `-16` | 1 |
| `-25` | 2 |
| `-36` | 3 |
| `-37` | 2 |
| `-40` | 1 |
| `-44` | 2 |
| `-45` | 2 |
| `-48` | 2 |
| `-52` | 1 |
| `-56` | 1 |
| `-64` | 2 |

Sixteen of the 70 weighted classes have nonsquare absolute determinant and
therefore cannot have a finite-index unimodular overlattice, since
`|det L| = [U:L]^2 |det U|`.  The other 54 classes all possess at least one
even-unimodular extension.  Complete discriminant quadratic-form enumeration
gives one extension for 51 classes and two extensions for three classes: 57
integral `E10` branches.  Each is exactly replayed as even, unimodular, of
signature `(1,9)`, containing all twelve roots and reproducing their Gram.

The geometric obstruction is stronger and independent of that lattice gate.
Exact GMP double description is run on all 70 quotient cones, including the
16 nonsaturable ones, and again on all 57 integral markings.  Every cone has
all twelve input walls irredundant and has between 10 and 49 hyperideal
extreme rays.  Each integral marking reproduces exactly the finite/ideal/
hyperideal counts of its quotient cone.  The complete ray records and their
hyperideal sublists are frozen by SHA-256 for every class and branch.

Hyperideal is not being confused with ideal.  Ideal rays may be finite-volume
cusps; one hyperideal extreme ray already makes an acute Coxeter chamber
infinite-volume.  Literal Vinberg therefore stops at the explicit hyperideal
endpoint in all 70 classes.

## A corrected exact interior witness

There is a useful methodological by-product.  The sum of all extreme rays is
strictly inside every cone, but in two of the 70 classes that sum has negative
square because the hyperideal rays dominate it.  Rejecting those cones merely
because this convenient sum is spacelike would be wrong.

The certificate instead chooses a finite extreme ray `v`, lets `s` be the
complete ray sum, and doubles `N` until

```text
x = N v + s
```

has positive square.  All wall slacks of `x` remain strictly positive because
those of `s` are already positive, and the square is checked exactly.  Sixty-
eight classes use `N=0`; the exceptional classes use `N=1` and `N=2`.
Every `E10` marking receives the same exact witness construction.  This
removes a ray-sum heuristic from the polyhedral gate without changing the
hyperideal verdict.

## Combined 15-edge criterion through two pairing-two edges

Together with the frozen zero- and one-pairing-two censuses, this gives a
clean P-first classification for the union

```text
15 support edges, maximum degree 3, and 0 <= pairing-two count <= 2.
```

Across the three exhaustive shards there are 694,661 marked trials, 114
Lorentzian marked survivors, 100 weighted root classes, and 81 integral
`E10` branches.  Exactly one weighted root cone has finite volume: the
zero-pairing-two Type-II class.  Its two integral markings return the same
known Type-II full polytope and pass the frozen `raw_30` split-K3/Coble
compiler.  The other 99 root classes, including every class with one or two
pairing-two edges, have hyperideal rays.

Thus, within this declared union, finite volume is equivalent to being the
Type-II root diagram.  This is a genuine criterion combining a positive
return with two negative boundaries; it is not a claim outside the union.

## Lineage and strict boundary

The zero-pairing-two 15-edge lane isolates the known Type-II return.  The
exactly-one-pairing-two lane has 28 hyperideal classes, and the present
exactly-two-pairing-two lane adds 70 more.  None of the present root diagrams
matches the frozen Type-I or Type-II root diagram.  This comparison is only
against the project ledger; it is not a literature-wide novelty claim about
the 70 Lorentzian diagrams.

No cusp lift is run here because none is mathematically licensed by the
finite-volume gate.  The result excludes no graph of degree above three, no
15-edge graph with three or more pairing-two edges, and no other edge count.
It proves neither strict L4, a general dimension bound, nor the requested
class-wide `POLYTOPE -> SURFACE -> AG pushback` theorem.  It produces no new
finite-volume hyperbolic Coxeter polytope.

## Replay

```bash
PYTHONPATH=src:. python3 scripts/certify_coble_nodot_15edge_two_pair2_lane.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_coble_nodot_15edge_two_pair2_lane.py
```
