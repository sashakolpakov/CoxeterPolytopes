# The 15-edge, one-pairing-two lane has no finite-volume survivor

**Status (2026-08-09): complete exact finite-denominator theorem.  Every
Lorentzian root cone in this lane has a hyperideal extreme ray.  This is a
negative search-boundary result, not a new compact or finite-volume Coxeter
polytope and not the full nodot theorem.**

The executable certificate is
[`certify_coble_nodot_15edge_one_pair2_lane.py`](../scripts/certify_coble_nodot_15edge_one_pair2_lane.py),
with [frozen archive](../data/enumerations/coble_nodot_15edge_one_pair2_lane.json)
and [tests](../tests/test_coble_nodot_15edge_one_pair2_lane.py).  The archive
SHA-256 is

```text
1e3e20026c418c6857ceb7a834693fd9217479e05506ddcab93d16304c38b084
```

## Exact scoped theorem

Take twelve proposed norm-minus-two roots.  Require their nonorthogonality
support to be connected, have exactly fifteen edges and maximum degree three,
mark exactly one support edge with pairing two, and give every other support
edge pairing one.  No surface, determinant, cusp count, known diagram, or
volume condition is supplied as input.

For every resulting Cartan matrix of exact inertia `(9,1,2)`, quotient its
saturated rank-two radical.  The resulting twelve-wall Coxeter cone in
`H^9` has a hyperideal extreme ray.  Consequently no such cone is compact or
finite-volume.  This remains true before imposing the `E10` surface-lattice
condition, and hence also for every even-unimodular `E10` saturation.

This is an exhaustive statement for the declared lane:

```text
5,741 unlabelled connected supports
  x 15 possible marked pairing-two edges
= 86,115 exact marked trials
-> 84,624 rank 12 + 1,462 rank 11 + 29 rank 10
-> all 29 rank-ten trials have exact inertia (9,1,2)
-> 28 weighted-isomorphism classes
-> all 28 quotient cones have hyperideal rays
-> 18 E10-saturable classes, with 21 integral saturation branches
-> all 21 E10 branches have hyperideal rays
-> zero finite-volume survivors and zero surface-return attempts.
```

Rank is computed over the integers on every marked trial.  Exact rational
congruence supplies inertia only after the rank-ten gate; there is no floating
spectral prefilter.  The 29 count is a count of edge markings on one
unlabelled representative of each support, not a count of vertex-labelled
graphs.  Weighted isomorphism removes the sole duplicate and leaves 28
classes.

## Lattice and volume obstructions

Ten of the 28 radical quotients have nonsquare absolute determinant.  They
cannot possess a finite-index unimodular overlattice because

```text
|det L| = [U : L]^2 |det U|.
```

The other eighteen classes have square determinant.  Complete discriminant
quadratic-form enumeration gives one `E10` saturation for fifteen classes and
two for three classes, or 21 integral branches in total.  Every alleged
extension is replayed integrally: it is even, unimodular, has signature
`(1,9)`, contains all twelve roots, and reproduces the root Gram.

Exact GMP double description is run on all 28 quotient cones, including the
ten without an `E10` extension, and again in every one of the 21 integral
markings.  All twelve walls are irredundant, every cone has a positive strict
interior point, and all causal rays lie in one time cone.  Nevertheless every
cone has between 7 and 34 hyperideal extreme rays.  This is the decisive
Vinberg obstruction: ideal rays alone could describe finite-volume cusps,
whereas a hyperideal extreme ray makes the chamber infinite-volume.  The
literal Vinberg audit therefore stops at the hyperideal endpoint rather than
misreporting these chambers as finite-volume noncompact.

There is no cusp completion to perform.  The one-boundary lift is defined
only after the root chamber passes the finite-volume gate; lifting the ideal
rays while hyperideal rays remain would not produce the desired complete
polytope.  Thus the split-K3/Coble same-`P` return is correctly recorded as
unreached, not failed after an attempted construction.

## Frozen class ledger

`F/I/H` records the exact finite, ideal, and hyperideal extreme-ray counts of
the quotient cone.  The same triple is reproduced in every `E10` marking of
that class.

| class | graph6 | marked edge | marked members | quotient det | E10 branches | F/I/H |
|---:|---|---|---:|---:|---:|---:|
| 0 | ``K??ED?WT@gX?`` | `(4, 8)` | 1 | -1 | 1 | 8/15/26 |
| 1 | ``K??ED?WT@g?[`` | `(4, 8)` | 1 | -16 | 1 | 10/12/13 |
| 2 | ``K??E@_KwCoD_`` | `(3, 7)` | 1 | -64 | 1 | 6/6/22 |
| 3 | ``K??E@_Kw@WBO`` | `(3, 7)` | 1 | -56 | 0 | 10/0/24 |
| 4 | ``K??E@_Kw@WPG`` | `(3, 7)` | 1 | -80 | 0 | 14/0/20 |
| 5 | ``K??E@_Kw@W@W`` | `(3, 7)` | 1 | -72 | 0 | 10/2/22 |
| 6 | ``K??E@_Kw@WAS`` | `(2, 7)` | 2 | -12 | 0 | 13/2/19 |
| 7 | ``K??E@_Kw@WAS`` | `(3, 7)` | 1 | -48 | 0 | 12/2/20 |
| 8 | ``K??E@_Kw@W?[`` | `(3, 7)` | 1 | -64 | 2 | 10/8/16 |
| 9 | ``K??E@_Kod_GK`` | `(0, 10)` | 1 | -48 | 0 | 11/2/11 |
| 10 | ``K??E@_Kod_GK`` | `(1, 11)` | 1 | -12 | 0 | 17/2/27 |
| 11 | ``K?AAD@OJ?s[?`` | `(4, 9)` | 1 | -1 | 1 | 14/11/25 |
| 12 | ``K?AAD?oUAKAS`` | `(4, 9)` | 1 | -4 | 1 | 14/5/30 |
| 13 | ``K?AAD?W[?kOK`` | `(4, 8)` | 1 | -16 | 1 | 12/5/16 |
| 14 | ``K?AAD?WUBGOK`` | `(8, 11)` | 1 | -40 | 0 | 14/1/13 |
| 15 | ``K?AAD?WS`gX?`` | `(4, 8)` | 1 | -49 | 1 | 8/6/18 |
| 16 | ``K?AAD?WS`gL?`` | `(4, 8)` | 1 | -61 | 0 | 10/0/22 |
| 17 | ``K?AAD?WS`gWC`` | `(4, 8)` | 1 | -64 | 1 | 10/5/17 |
| 18 | ``K?AAD?WS`gKC`` | `(4, 8)` | 1 | -76 | 0 | 10/0/22 |
| 19 | ``K?AA@_geECHG`` | `(0, 9)` | 1 | -1 | 1 | 12/4/34 |
| 20 | ``K?AA@_ghAoQO`` | `(2, 7)` | 1 | -1 | 1 | 11/11/30 |
| 21 | ``K?AA@_ghAoWG`` | `(2, 8)` | 1 | -1 | 1 | 6/11/27 |
| 22 | ``K?AA@_ghAoQG`` | `(2, 8)` | 1 | -1 | 1 | 8/8/31 |
| 23 | ``K?AA@_goO[DC`` | `(2, 7)` | 1 | -9 | 2 | 11/16/15 |
| 24 | ``K?AA@_goO[DC`` | `(4, 10)` | 1 | -36 | 2 | 11/10/7 |
| 25 | ``K?AA@_gaRCX?`` | `(2, 8)` | 1 | -1 | 1 | 6/10/34 |
| 26 | ``K?AEB?oEAWAc`` | `(4, 9)` | 1 | -1 | 1 | 11/9/27 |
| 27 | ``K?AEB?oEAWDA`` | `(2, 7)` | 1 | -1 | 1 | 16/9/28 |

The archive also freezes SHA-256 digests of the complete exact ray records
and of the hyperideal sublist for every quotient cone and every integral
branch.

## Lineage and strict boundary

The earlier 15-edge/no-pairing-two lane returns the known Type-II chamber;
the 13-edge/two-or-three-pairing-two lane returns Type I.  The present census
fills the intervening **15-edge/exactly-one-pairing-two/max-degree-three**
shard.  None of its 28 root diagrams matches either archived Type-I or
Type-II root diagram, and none survives to a full-wall comparison.  This is
lineage relative to the frozen Coble ledger only; it is not a claim that the
28 Lorentzian diagrams have never appeared anywhere in the literature.

The result excludes one natural P-first combinatorial stratum before algebraic
geometry is invoked.  It is useful search leverage, but it is not the desired
class-wide `POLYTOPE -> SURFACE -> AG pushback` theorem: it excludes no graph
of degree above three, no fifteen-edge graph with two or more pairing-two
edges, and no other edge count.  It proves neither strict L4 nor a general
dimension bound.  It also produces no new compact or finite-volume
hyperbolic Coxeter polytope.

## Replay

```bash
PYTHONPATH=src:. python3 scripts/certify_coble_nodot_15edge_one_pair2_lane.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_coble_nodot_15edge_one_pair2_lane.py
```
