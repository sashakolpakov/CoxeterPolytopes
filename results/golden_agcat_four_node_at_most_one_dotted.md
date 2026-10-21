# Golden AGCat four-node packets with at most one dotted edge

## Exact scope

This certificate advances the \(\mathbf Q(\sqrt5)\), dimensions 22--29
queue from isolated three-cycles to simultaneous compatibility on connected
four-node graphs.  It exhausts two natural tranches:

1. every connected four-node assignment whose six pair slots are finite or
   orthogonal; and
2. every connected assignment with exactly one dotted edge, where that edge
   ranges over the complete 6,869-value dotted alphabet and every other slot
   ranges over all six finite/orthogonal values.

Thus zero edges, both tree topologies, triangles with a pendant vertex,
four-cycles, diamonds with two overlapping triangles, and complete \(K_4\)
packets all occur.  The computation does **not** cover packets with two or
more dotted edges, nor does it enumerate complete connected frames of ranks
23--30.

## Exact orbit reduction

The finite-only tranche is quotiented by all of \(S_4\).  It has 2,380
connected edge-labelled orbits before cycle compatibility.

When a dotted edge is unique, it can canonically be moved to \(01\).  The
remaining quotient is its order-four setwise stabilizer
\(S_{\{0,1\}}\times S_{\{2,3\}}\).  There are 2,080 finite-complement
orbits.  Expanding the distinguished slot by all dotted values gives

\[
  2{,}080\cdot6{,}869=14{,}287{,}520
\]

canonical assignments before the cycle gate.

## Simultaneous cycle compatibility

For every present graph cycle \(C\), symmetrizability requires

\[
  \prod_{e\in C}c_e=r_C^2,
  \qquad r_C\in\mathbf Z[\phi].
\]

The certifier checks the element-squareclass XOR on all four triangles and
all three Hamilton four-cycles simultaneously.  This is not merely a
necessary filter.  On a connected graph, choose a spanning tree and
propagate a totally positive symmetrizer along it.  Trivial squareclass on
each fundamental cycle makes the corresponding chord quotient a square in
\(\mathbf Q(\sqrt5)\); hence the field-level symmetrizable Cartan data exist.
The cycle roots are golden integers because \(\mathbf Z[\phi]\) is
integrally closed.  Choosing an integral directed-edge gauge and reflection
order remains a later stage (except in the explicit path family below, where
an integral gauge is given).

In the one-dotted tranche this yields a particularly cheap exact join:

- every finite cycle avoiding \(01\) must already have trivial squareclass;
- all cycles containing \(01\) must prescribe the same squareclass for its
  dotted value;
- if \(01\) is a bridge, it has no cycle constraint and all 6,869 values
  survive this stage.

The result is 327 compatible finite-only packets and 532,734 compatible
one-dotted packets.

## Two-place determinant test

For normalized diagonal entries 2, determinant expansion by permutation
cycle type gives

\[
\det G=
16-4\sum_e c_e
+\sum_{e\perp f}c_ec_f
-4\sum_T\sqrt{\prod_{e\in T}c_e}
-2\sum_Q\sqrt{\prod_{e\in Q}c_e}.                 \tag{1}
\]

The first inner sum is over the three disjoint edge pairs, the next over the
four triangles, and the last over the three unoriented Hamilton cycles.
Every root in (1) is the identity-positive golden root.  Acuteness makes an
odd triangle cycle negative and an even four-cycle positive; conjugating the
result preserves the correlated signs rather than independently choosing
positive square roots at the second place.

At the conjugate place all pair minors are already strictly positive.
Exact signs of the leading \(3\times3\) determinant and (1) therefore give
Sylvester's positive-definiteness test.  At the selected place finite-only
packets may be positive definite or Lorentzian.  A one-dotted packet already
has one negative direction on its dotted \(2\times2\) minor and must remain
Lorentzian.

No selected or conjugate zero determinant occurs in the compatible
one-dotted tranche.  Among finite-only packets, 21 have a zero leading
three-node determinant and three different packets have a zero full
four-node determinant, at both embeddings.  Strict conjugate
positive-definiteness removes all 24.  The split matters: the former are
degenerate principal subdiagrams, while only the latter are degenerate full
four-node Grams.  This handles the affine/parabolic boundary explicitly
rather than silently identifying it with a compact packet.

## Exact counts

| tranche | connected canonical inputs | cycle-compatible | two-place retained |
|---|---:|---:|---:|
| finite/orthogonal only | 2,380 | 327 | **119** |
| exactly one dotted | 14,287,520 | 532,734 | **180,402** |

The 119 finite-only survivors split into 5 selected-positive and 114
selected-Lorentzian packets.  All 180,402 one-dotted survivors are
selected-Lorentzian.  Indeed, every one of the 532,734 cycle-compatible
one-dotted packets is selected-Lorentzian; all rejection in that tranche
comes from conjugate positive-definiteness.

The topology split is:

| tranche | topology | cycle-compatible | retained |
|---|---|---:|---:|
| finite-only | star \(K_{1,3}\) | 35 | 13 |
| finite-only | path \(P_4\) | 75 | 36 |
| finite-only | triangle + pendant | 90 | 32 |
| finite-only | cycle \(C_4\) | 33 | 10 |
| finite-only | diamond \(K_4-e\) | 62 | 19 |
| finite-only | complete \(K_4\) | 32 | 9 |
| one dotted | star \(K_{1,3}\) | 103,035 | 28,845 |
| one dotted | path \(P_4\) | 274,760 | 105,919 |
| one dotted | triangle + pendant | 135,547 | 39,686 |
| one dotted | cycle \(C_4\) | 3,927 | 1,507 |
| one dotted | diamond \(K_4-e\) | 11,022 | 3,277 |
| one dotted | complete \(K_4\) | 4,443 | 1,168 |

The nonzero diamond and complete-graph counts certify genuine overlapping
triangle and longer-cycle joins; this is not just a tree enumeration.

The four compatible/retained streams are frozen by SHA-256 digests in the
JSON artifact.  Each record contains the six canonical pair products, the
leading \(012\) determinant, and the full determinant, all as exact
\(a+b\phi\) coefficient pairs.

## An exact global rank-23--30 packet family

The same computation exposes a full-rank family that can be solved without a
large join.  In native rank \(r\), take a path, put one dotted product \(c\)
on its terminal edge, and put Coxeter label 3 (product 1) on every remaining
edge.  A tree has no cycle equations, so symmetrizability is exact and
automatic.  An integral Cartan realization is

\[
 a_{01}=-1,\quad a_{10}=-c,\qquad
 a_{i,i+1}=a_{i+1,i}=-1\quad(i\ge1),
\]

with symmetrizer \(\operatorname{diag}(c,1,\ldots,1)\).

If \(D_j\) is the normalized leading determinant on the first \(j\) nodes,
the tridiagonal recurrence gives

\[
 D_j=2j-(j-1)c.                                      \tag{2}
\]

At the selected place, \(c>4\) makes \(D_1>0\) and every \(D_j<0\) for
\(j\ge2\), hence the full frame has exactly one negative direction.  At the
conjugate place, Sylvester positivity through rank \(r\) is equivalent to

\[
 c'<{2r\over r-1}.
\]

Exhausting the entire dotted alphabet gives:

| dimension | native rank | admissible dotted products |
|---:|---:|---:|
| 22 | 23 | 3,590 |
| 23 | 24 | 3,584 |
| 24 | 25 | 3,577 |
| 25 | 26 | 3,571 |
| 26 | 27 | 3,566 |
| 27 | 28 | 3,561 |
| 28 | 29 | 3,556 |
| 29 | 30 | 3,552 |

The single choice \(c=2+2\phi\), whose conjugate is
\(4-2\phi\in(0,1)\), works uniformly in all eight ranks.  These are genuine
connected, globally compatible full-rank golden frame packets with integral
Cartan reflections in admissible arithmetic Lorentzian lattices.

Moreover, these ambient arithmetic lattices are cocompact.  A nonzero
\(K\)-isotropic vector would remain isotropic after conjugation, contradicting
positive definiteness of the conjugate form.  Thus the quadratic form is
anisotropic over \(K\), and its full arithmetic orthogonal group is uniform.
Their minimal Vinberg field is genuinely \(K=\mathbf Q(\sqrt5)\): a dotted
alphabet value cannot be rational, since its selected value is greater than
4 while its conjugate is less than 4, and that value \(c\) is already a
length-two cyclic Cartan product.

They are deliberately **not** advertised as new compact examples.  The
displayed \(r\)-wall simplicial cone has a dotted proper subdiagram and is
not a finite-volume Coxeter simplex.  Its arithmetic host still needs a
complete Vinberg run: additional mirrors could close a compact chamber, end
with an infinite-index reflection subgroup, or otherwise fail reflectivity.
Because the host lattice itself is anisotropic, any finite-index reflection
chamber would be compact.  What (2) proves is that local/global frame
compatibility alone will not empty the golden queue.

## What this proves, and what it does not

This is an exhaustive theorem for connected four-node golden packets with
at most one dotted edge.  It is a strict advance beyond the triangle census:
all cycles on the four vertices must coexist in one symmetrizable packet and
the full determinant must pass both real places.

It is not a compact Coxeter polytope census.  The 180,521 retained local
packets and the explicit full-rank path seeds are inputs for ideal gauges,
reflection orders, conductor lattices, and complete Vinberg-host completion.
A rank-23--30 frame may contain four-node subsets with several dotted edges,
so the local census is not universal; and none of the displayed full-rank
seeds has been shown reflective (equivalently, to have a finite-index
reflection subgroup).  Their ambient arithmetic lattices are cocompact, as
proved above.  This result therefore proves neither golden nonexistence nor
a new compact Coxeter example.  The status of \(\mathbf Q(\sqrt5)\) in
dimensions 22--29 remains **OPEN**.

Reproduce with:

```bash
PYTHONPATH=src python3 scripts/certify_golden_agcat_four_node_at_most_one_dotted.py --verify
PYTHONPATH=src pytest -q tests/test_golden_agcat_four_node_at_most_one_dotted.py
```

A no-import hostile implementation independently reconstructs the alphabet,
squareclasses, orbit quotients, generic permutation determinants, inertia,
all 48 strata, all four stream hashes, and the rank-23--30 path family.  Its
verdict is `PASS`; see
[`golden_agcat_four_node_at_most_one_dotted_hostile_audit.md`](golden_agcat_four_node_at_most_one_dotted_hostile_audit.md).
