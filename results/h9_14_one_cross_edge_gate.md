# Exact one-cross-edge gate for the five frozen earlier H9/14 families

> **Scope update.** The later sixth ordinary-polytopal nerve, record 35, is
> outside these five Hamming shells and is instead excluded from the compact
> Coxeter track by `results/h9_14_record35_coxeter_label_gate.md`.

## Result

The first layer outside the five frozen factor-orthogonal/sparse Coxeter-label
families is now exhausted exactly.

For every full-nerve symmetry orbit of a frozen baseline labelling, change
exactly one currently right-angled cross-factor meeting pair to
`m = 3, 4, or 5`.  Quotienting by the labelled stabilizer of each baseline
gives **609,078** one-edge extension orbits.  Exactly **6,342** of them pass
all finite-vertex and Lanner-minimal-nonface tests.  Every one of those local
survivors nevertheless has a principal Gram submatrix with at least two
negative directions.

Thus one cross edge *does* often suffice locally, but it never suffices for
Lorentzian inertia.  Relative to each frozen family:

> any Gram candidate of negative index one must contain at least two new
> cross-factor non-right finite edges.

This is a complete theorem about the one-edge shell.  It is not an
unrestricted Coxeter-label census and is not a new hyperbolic polytope.

## Exact census

| record | frozen baseline orbits | right cross pairs | fully labelled one-edge trials | one-edge symmetry orbits | locally admissible labelled extensions | locally admissible orbits | metric survivors |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 30 | 45 | 76 | 2,659,392 | 1,218 | 575,424 | 180 | 0 |
| 42 | 6,120 | 65 | 9,057,484,800 | 603,840 | 78,520,320 | 6,120 | 0 |
| 43 | 15 | 45 | 1,701,000 | 1,440 | 19,440 | 16 | 0 |
| 44 | 25 | 45 | 8,164,800 | 2,220 | 95,040 | 20 | 0 |
| 47 | 5 | 45 | 1,020,600 | 360 | 16,200 | 6 | 0 |
| **total** | **6,210** | — | **9,071,030,592** | **609,078** | **79,226,424** | **6,342** | **0** |

“Fully labelled” counts actual assignments on the fixed fourteen facet
labels, not canonical representatives.  They provide an independent closure
check on the symmetry quotient.  For example, the record-43 count is

\[
  12,600\cdot45\cdot3=1,701,000.
\]

## Why the quotient is exact

The earlier certificates give one representative for every baseline orbit
under the full automorphism group of the complete minimal-nonface
hypergraph.  For each such representative `b`, the new computation:

1. generates its full labelled stabilizer `Stab(b)`;
2. computes the exact `Stab(b)`-orbits on the cross-factor pairs;
3. treats labels 3, 4, and 5 separately; and
4. weights every pair orbit by `|Aut(N)|/|Stab(b)|` times its pair-orbit
   size to recover the number of fully labelled extensions.

An extension contains exactly one non-right cross edge, so resetting that
edge to label two recovers a unique baseline.  Hence extensions from two
different baseline orbits cannot collide.  The stabilizer quotient is
therefore the full-nerve quotient of the stated one-edge shell.

## Complete local audit

Every baseline was already checked on every rank-nine vertex and every
higher minimal nonface.  Changing a cross-factor pair has two useful exact
consequences:

- only vertices containing both endpoints can change; every such vertex is
  retested with the complete finite Coxeter classifier
  `A,B,D,E,F,H,I2`, component by component;
- no cross-factor pair lies inside an order-four or order-five minimal
  nonface, so every Lanner diagram is literally unchanged.

For the product records, the affected vertices are audited through their
verified Cartesian factorization: all affected sections in the two joined
factors are tested, and their products with the unchanged orthogonal factor
sections add no further condition.  For the multiwedges, the archived list
of all 90 vertices is used directly.

The local survivor label histograms also show that the gate is genuinely
selective.  In the three multiwedges the surviving orbit counts by new label
are

| record | `m=3` | `m=4` | `m=5` |
|---:|---:|---:|---:|
| 43 | 11 | 5 | 0 |
| 44 | 14 | 6 | 0 |
| 47 | 4 | 2 | 0 |

Thus no one-edge `m=5` multiwedge extension even passes the finite-vertex
gate, although many `m=3` and `m=4` extensions do.

## Global metric obstruction

The metric rejection does not choose dotted values and does not use a
floating-point full Gram matrix.  It exhibits a principal submatrix with
negative index at least two, which excludes any full normalized Gram matrix
of Lorentzian negative index one by interlacing.

### Record 30

The orthogonal baseline has five negative factor blocks: two Lanner
four-blocks and three dotted two-blocks.  One cross edge can join at most two
blocks.  The coupled component still contains a negative principal block,
while three negative blocks remain orthogonal.  Hence every one-edge
extension has

\[
                         n_-\geq 4.
\]

This holds for every legitimate dotted value.

### Record 42

There are two Lanner blocks and a polygon block containing a dotted
two-by-two negative principal submatrix.  If the new edge joins the polygon
to one Lanner block, the other Lanner remains orthogonal.  If it joins the
two Lanner blocks, the negative polygon principal block remains orthogonal.
In either case the coupled component contains another negative principal
block, so

\[
                         n_-\geq 2.
\]

### Records 43, 44, and 47

The isolated order-five Lanner block is joined to the residual nine-facet
block.  Of the 42 locally admissible one-edge orbits:

- 40 leave some residual Lanner support untouched by the new endpoint; its
  orthogonal direct sum with the isolated Lanner block has two negative
  directions;
- one record-43 orbit and one record-44 orbit meet the common core of all
  residual Lanner supports and require an exact coupling calculation.

For two normalized Lanner Gram blocks `A` and `B`, coupled only at nodes
`a,b` by `t=-2 cos(pi/m)`, Schur complementation gives the determinant factor

\[
 q=1-t^2(A^{-1})_{aa}(B^{-1})_{bb}.
\]

Both inverse diagonal entries are negative because the full Lanner
determinant is negative and the corresponding one-node-deletion determinant
is positive.  A positive `q` leaves one negative direction in the Schur
complement, in addition to the negative direction of `A`.  For each of the
two nontrivial survivor orbits the exact factor is

\[
 q=\frac{50-22\sqrt5}{16(2-\sqrt5)^2}>0.
\]

The numerator is positive because `50^2 = 2500 > 2420 = (22 sqrt(5))^2`.
Thus these final two principal blocks also have negative index two.

## Algebro-geometric design consequence

For a divisor-lattice construction, the result sharpens the earlier
indecomposability requirement.  The incidence varieties and their local
Lanner decorations are not enough, and neither is a single extra
cross-intersection.  Within the five frozen decoration families, the
Neron--Severi/root pairing must create at least two nonzero cross-factor
couplings while keeping every rank-nine divisor configuration elliptic.

This is still a restricted inverse-design statement: in records 43, 44, and
47, simultaneous changes to residual meeting labels outside the four Lanner
supports were not enumerated.  Those changes belong to the next CSP layer.

## Reproducibility

The exact implementation is
`src/coxeter_polytopes/combinatorics/h9_14_one_cross_edge.py`.

The frozen archive, including every one of the 6,342 surviving orbit
representatives, is
`data/enumerations/h9_14_one_cross_edge.json`, with certificate digest

```text
b6b440c8649725bbc49a43662036dc8765423b252e7c8092d85f8f7436e3e4f7
```

Replay with

```bash
PYTHONPATH=src:. python3 scripts/certify_h9_14_one_cross_edge.py --verify
PYTHONPATH=src:. pytest -q tests/test_h9_14_one_cross_edge.py
```
