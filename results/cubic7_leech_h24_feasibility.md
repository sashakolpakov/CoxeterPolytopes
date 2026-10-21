# Cubic-seven Leech-spatial H24 feasibility audit

> **Continuation note (2026-08-13).**  The first-shell statements below
> remain valid, but the stated continuation blocker has now been bypassed for
> the first two realised packets and for a complete necessary-wall obstruction
> through `H=1000`.  See
> [`cubic7_leech_h24_h1000_continuation.md`](cubic7_leech_h24_h1000_continuation.md).

Consider

\[
L=\langle-2\theta\rangle\oplus\Lambda_{24}
\quad\text{over}\quad O=\mathbf Z[\theta],
\qquad \theta=2\cos(2\pi/7),
\]

where the Leech lattice is positive definite, even, unimodular, and has
minimum four.  This audit gives a complete local mirror-norm proof, a complete
height-zero census, and a finite exact first-shell calculation.  It is not a
reflectivity, finite-volume, or compactness result.

Reproduction:

```console
python3 scripts/audit_cubic7_leech_h24_first_shell.py
PYTHONPATH=src pytest -q tests/test_cubic7_leech.py
```

## Complete mirror norms

For a primitive vector `r=(x,y)` of norm `q`, reflection integrality is

\[
q\mid 2I(r),\qquad I(r)=(2\theta x,(y,\Lambda_{24})).
\]

Because the Leech lattice is even, `q` is always divisible by the inert prime
`(2)`.  At an odd prime, unimodularity makes divisibility of the spatial
pairing ideal equivalent to divisibility of `y`; since `theta` is a unit,
divisibility of `2 theta x` also forces divisibility of `x`, contradicting
primitivity.  Hence no odd prime divides `q`.

At `(2)`, if `y` is not divisible by two then its unimodular pairing ideal is
a unit, so `q | 2 I(r)` gives `v_2(q)<=1`.  If `y=2z`, primitivity makes `x` a
unit and

\[
q/2=-\theta x^2+2(z,z)
\]

is a 2-adic unit.  Thus `v_2(q)=1` in both cases.  The field has class number
one and full unit-signature rank, so every totally positive unit is a square.
The complete canonical mirror-norm list, modulo root rescaling by unit
squares, is therefore exactly

\[
\boxed{q=2}.
\]

This is cleaner than `E8^3`: the coefficient `2 theta` does not introduce a
second norm class.

## Height zero: exactly no roots

Write a spatial vector in the power basis, then make the integral change

\[
y=-A+(-1+\theta+\theta^2)B-\theta C,
\qquad A,B,C\in\Lambda_{24}.
\]

The trace of its field-valued Leech square is the exact positive sum

\[
\begin{aligned}
\operatorname{Tr}(y,y)={}&\|A-B\|^2+\|A-C\|^2+2\|B-C\|^2\\
&+\|A\|^2+2\|B\|^2+2\|C\|^2.
\end{aligned}
\]

An exhaustive equality-pattern check shows that every nonzero triple has at
least three units of nonzero norm contribution.  Since the Leech minimum is
four,

\[
\operatorname{Tr}(y,y)\ge12,
\]

with equality only when `B=C=0` and `A` is a minimal Leech vector.  A
height-zero mirror would have `(y,y)=2`, whose trace is six.  Therefore the
height-zero reflection census is complete and empty.

In particular, the base point `p=(1,0)` is regular: it lies on no reflection
mirror.  There is no height-zero Weyl group or dominant Weyl cone to choose.

## Standard Leech minimal shell

The implementation constructs the extended binary Golay code from the cyclic
`[23,12,7]` generator and checks its weight distribution

| weight | 0 | 8 | 12 | 16 | 24 |
|---:|---:|---:|---:|---:|---:|
| codewords | 1 | 759 | 2576 | 759 | 1 |

It then streams the standard minimal Leech coordinates, scaled by `sqrt(8)`:

| coordinate type | count |
|:--|--:|
| `(4^2,0^22)` | 1,104 |
| `(2^8,0^16)` on octads, even sign parity | 97,152 |
| `(3,1^23)` with Golay sign changes | 98,304 |
| **total** | **196,560** |

Every vector has scaled square 32, hence Leech norm four.  The total also
checks the transitive Conway orbit `|Co0|/|Co2|=196,560`.

## First positive control shell

Exact three-embedding enumeration through Vinberg height 10 gives one and
only one positive control level:

\[
x=-1+\theta+\theta^2,
\quad h=2\theta+\theta^2,
\quad (y,y)=2+4\theta+2\theta^2.
\]

The target trace is 12.  Hence every possible spatial vector is in the
trace-minimum equality case and comes from one of the 196,560 minimal Leech
vectors with `B=C=0`.  Such a vector has field-valued square exactly 4, not
`2+4 theta+2 theta^2`.  The complete first positive packet is therefore
empty.  This is a genuine finite first-shell result, not sampling.

## Equal-height packets without a Weyl chamber

The absence of height-zero roots simplifies the chamber definition but
removes the dominant-weight quotient.  A correct continuation is:

1. Keep `p=(1,0)` fixed and orient every root by `(p,r)<0`.  Since `p` is
   regular, it belongs to a unique reflection chamber; no arbitrary initial
   Weyl chamber is required.
2. Enumerate every root at one exact Vinberg height, modulo the stabilizer in
   `Co0` of all lower walls already accepted.  Initially the stabilizer is
   all of `Co0`.
3. Discard roots having positive pairing with a lower wall.  Test the entire
   surviving equal-height packet simultaneously.  Complete Vinberg
   enumeration predicts that survivors are pairwise nonpositive.  A positive
   pair is an audit failure—evidence of a missing lower root or incomplete
   orbit expansion—not permission to choose one representative arbitrarily.
4. Avoid quadratic expansion by a stabilizer chain.  A one-vector minimal
   orbit uses `Co0 -> Co2`; pair orbits use `Co2` suborbits classified by exact
   Leech inner product; triple orbits continue with the exact pair stabilizer.
   Pairwise acuteness is then checked once per double-coset relation.

## Former continuation blocker and its resolution

The next field target is not determined by the ordinary Leech theta series.
Writing `y` as the triple `(A,B,C)` above prescribes six simultaneous Leech
norms/pairings.  What is needed is a lossless canonical enumerator for
diagonal `Co0` orbits on ordered triples `(A,B,C)`, together with stabilizers
and orbit sizes.  The repository currently contains neither a Leech lattice
basis/closest-vector enumerator nor a callable `Co0` permutation/matrix
representation.

Thus a generic full wall-orbit Vinberg implementation is not available from
the one-vector theta series alone.  Treating theta-series coefficients as
triple counts, or choosing one wall from a positive equal-height packet,
would still be mathematically invalid.

The later continuation avoids those invalid steps.  Exact coefficient-Gram
equations force the first two nonempty levels into the minimal shell and
ordered pairs of minimal vectors, where the complete orbit incidences can be
counted.  Above them, adjoining a fixed norm-six Leech vector and maximising
over every positive-semidefinite minimum-four Gram extension proves a strict
hyperideal obstruction without asserting that every necessary triple is
realised.  Full wall-orbit enumeration above the certified prefix remains a
separate missing primitive.
