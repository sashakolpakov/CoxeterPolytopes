# Complete golden AGCat triangle packets, including dotted edges

## Exact scope

This certificate advances the \(\mathbf Q(\sqrt5)\), dimensions 22--29
AGCat queue from its length-two alphabet to the complete local three-node
tranche.  It uses all 6,875 Coxeter-filtered edge products in the proven
\(A=62\) frame box, including all 6,869 dotted values.

A complete three-node nonorthogonality graph has three nonzero edges, so the
single orthogonal value is removed.  Edge triples are unordered, exactly
quotienting by all triangle vertex permutations.

This is not a rank-23--30 graph census.  Compatibility among overlapping
triangles, longer cycles, ideal gauges, reflection orders, conductor
lattices, and complete Vinberg transcripts remain open.

## Scalable exact enumeration

There are

\[
\binom{6874+2}{3}=54{,}158{,}527{,}500
\]

raw unordered nonzero edge multisets, too many for naive materialization.
The certifier instead factors every principal edge ideal exactly, using
inert and ramified valuations and Hensel-lifted valuations at split primes.
It also retains the conjugate sign, which detects the non-square unit class
represented by \(\phi\).

The 6,874 values occupy 3,661 element-squareclass buckets; the largest has
109 elements.  XOR joining yields only 27,088 compatible bucket triples.
Streaming those buckets gives exactly 4,768,141 edge multisets for which

\[
c_{12}c_{23}c_{31}=r^2,\qquad r\in\mathbf Z[\phi].
\]

No list of 4.7 million records is stored.  Two SHA-256 stream digests freeze
the sorted edges, acute directed cycle, and determinant of every packet.

## Corrected golden square-root solver

For

\[
A+B\phi=(m+n\phi)^2,
\]

putting \(t=n^2\) gives

\[
5t^2-(4A+2B)t+B^2=0.
\]

The implementation solves this integer quadratic, reconstructs \(m\), and
multiplies the candidate root back exactly.  Regression tests include
\(\phi^2=1+\phi\), \((2+\phi)^2=5+5\phi\), and
\((3+2\phi)^2=13+16\phi\).

The identity-positive root \(r\) is chosen and the acute directed
three-cycle is \(b_{123}=-r\).  The normalized Gram determinant is

\[
D=8-2(c_{12}+c_{23}+c_{31})+2b_{123}.
\]

## Exact counts

| stage | exact count |
|---|---:|
| raw unordered nonzero edge multisets | 54,158,527,500 |
| element-squareclass buckets | 3,661 |
| compatible bucket triples | 27,088 |
| integral golden square-product packets | 4,768,141 |
| rejected because conjugate determinant is nonpositive | 2,134,256 |
| conjugate-positive packets retained | **2,633,885** |
| retained selected-positive packets | 0 |
| retained selected-Lorentzian packets | **2,633,885** |
| raw selected-degenerate packets | 1 |
| retained selected-degenerate packets | 0 |

At the conjugate place, positive definiteness and frame independence require
strictly positive determinant.  At the selected place an arbitrary local
packet may be positive, degenerate, or Lorentzian.  Degenerate affine
packets must be counted before the compactness filter; compact reflection
groups have no parabolic elements.  In this census no selected-degenerate
packet survives conjugate positivity anyway.  The sole raw zero-determinant
packet is also conjugate-degenerate.

The full stream digest is

```text
41386c1bfb64c497acc2ddb967c9a0d646f82e67dfc8907b8bc7d3b50e85bbe3
```

and the conjugate-positive retained stream digest is

```text
8be12b4371f49e1a0dd35f2c290c4e5735e64da56dc209655d1e0af5b49476b6
```

## Conclusion

This completely closes the local triangle stage, dotted edges included.
It does not prove the golden high-dimensional queue empty: 2,633,885 local
triangle packets survive, and the decisive work is their compatibility in
connected ranks 23--30 followed by the order/conductor/Vinberg stages.

In particular, the fact that no retained complete triangle is
selected-positive is not an obstruction to a compact chamber.  A connected
high-rank frame need not contain a complete nonorthogonal triangle at all:
its three-node induced subgraphs may be trees or may contain orthogonal
edges.  Nor has a universal triangular three-face been forced.  Conversely,
a selected-Lorentzian three-node restriction is locally compatible with a
global form of index one.  Only the overlapping-cycle and full-rank
compatibility computation can decide which of these local packets actually
extend.

Reproduce with:

```bash
PYTHONPATH=src python3 scripts/certify_golden_agcat_triangle_packets.py --verify
PYTHONPATH=src pytest -q tests/test_golden_agcat_triangle_packets.py
```
