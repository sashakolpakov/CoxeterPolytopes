# Silver `N(A2^12)` exact H24 screen through height 10

## Outcome

The two exact arithmetic hosts

\[
 \langle-\sqrt2\rangle\oplus N(A_2^{12}),\qquad
 \langle-2\sqrt2\rangle\oplus N(A_2^{12})
\]

have no compact or finite-volume survivor through the complete Vinberg
prefix at height 10.  Each prefix has 24 height-zero walls and one packet of
12 positive walls.  Its section is a product of twelve triangles, with

\[
3^{12}=531441
\]

exact rays: 25 finite, no cusps, and 531416 hyperideal.

This is a new host-family screen, not a new compact example.  It is exhaustive
only through the displayed cutoff; later walls could cut the current rays.

## Exact ternary glue

The spatial lattice is constructed from (A_2^{12}) with the extended
ternary Golay code.  The archived generator matrix is checked directly over
\(\mathbf F_3\): it has rank six, is self-orthogonal and hence self-dual, and
has weight enumerator

\[
1+264y^6+440y^9+24y^{12}.
\]

The root-lattice determinant is (3^{12}), while the glue index is (3^6),
so the overlattice has determinant one.  Every nonzero codeword has weight
divisible by three; its discriminant representative has norm (2w/3),
which is even.  Thus the construction is even unimodular.  Since the minimum
nonzero code weight is six, a glue coset has norm at least four.  The complete
height-zero root system is therefore exactly (A_2^{12}), with 72 roots and
24 simple chamber walls.

## Complete mirror norms

For a primitive root (r=(x,y)), reflection integrality is

\[
(r,r)\mid 2(a x,(y,N)).
\]

The positive block is even unimodular.  The local ramified-dyadic argument
therefore gives the complete canonical norm lists

| coefficient | mirror norms |
|---|---|
| (a=\sqrt2) | (2,2+\sqrt2) |
| (a=2\sqrt2) | (2,4+2\sqrt2) |

Totally positive units are squares, so these lists lose no unit-square root
class.

## First packet and glue exclusion

For (a=\sqrt2), the first packet has

\[
q=2+\sqrt2,\qquad x=1+\sqrt2,\qquad y^2=6+4\sqrt2.
\]

In one dominant (A_2) block the exact reflective possibilities have square
zero, one-third, two-thirds, or all of this target.  Any non-affine
decomposition of the full target occupies at most three blocks and has a
nonzero rational or radical discriminant residue there.  Such a residue word
cannot lie in the ternary Golay code, whose minimum weight is six.  Hence the
only roots are the twelve affine highest-root caps

\[
(1+\sqrt2;(1+\sqrt2)(\alpha_{i,1}+\alpha_{i,2})),
\qquad i=1,\ldots,12.
\]

The (a=2\sqrt2) proof is identical, with norm (4+2\sqrt2) and spatial
scale (2+\sqrt2).  The exact block census and code-minimum argument are
replayed in the tests.

In dominant Dynkin-weight coordinates each cap makes its block an (A_2)
alcove triangle.  A product vertex choosing a nonzero vertex in (k) blocks
has exact Lorentz norm

\[
-a+\frac{2k}{3}c^2,
\]

where (c=\sqrt2) or (2).  Exactly the (k=0,1) vertices are finite;
every (k\ge2) vertex is hyperideal.  This gives

\[
1+12\cdot2=25
\]

finite vertices and (3^{12}-25=531416) hyperideal rays.

## Remaining controls through height 10

The silver coefficient has five controls and the twice-silver coefficient
has three.  Controls with spatial targets (2+\sqrt2),
(6+3\sqrt2), and (10+7\sqrt2) would force the integral pairing
(u\cdot v) to be half-integral.  The target (20+14\sqrt2) contradicts
minimum two after the three possible values (v^2=2,4,6) are checked.

The only later arithmetically realizable target is (10+6\sqrt2).  Writing
(y=u+\sqrt2v) gives

\[
u^2+2v^2=10,\qquad u\cdot v=3.
\]

Minimum two forces (v^2=2), (u^2=6), and ((u-v)^2=2).  Thus both
(v) and (u-v) are roots in one (A_2) component.  Its dominant weights
are ((\sqrt2,3+\sqrt2)), up to reversal, and their sum
(3+2\sqrt2) exceeds the earlier cap (2+2\sqrt2).  Every such vector is
therefore rejected by the first packet.

## Reproduction

```sh
PYTHONPATH=src python3 scripts/audit_silver_a2_ternary_golay_h24.py
PYTHONPATH=src pytest -q tests/test_silver_a2_ternary_golay.py
```

The frozen archive is
`data/enumerations/silver_a2_ternary_golay_h24_h10.json`.

