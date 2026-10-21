# Exact stabilizers of the bounded OG10 golden wall prefixes

> **Subsequent global result.**  The caveat below is now realized explicitly:
> the complete chamber has the involution `h=s0*A`, which transports walls
> across the control cutoff and lies in the enhanced birational image.  It
> does not preserve any of the three finite prefixes.  See
> [`og10_golden_complete_chamber_involution.md`](og10_golden_complete_chamber_involution.md).

## Scope first

This is a theorem about three **control-height prefixes**, not about the
eventual complete movable or Kähler chamber.  An isometry of a complete
infinite chamber need not preserve the subset of its walls selected by a
non-invariant height cutoff.  Consequently the result below does not bound
the birational automorphism group of an OG10 realization.  Such an
application would first require either termination of the wall search or an
intrinsic, full-chamber-invariant characterization of one certified finite
wall subset.

Within that scope, the finite answer is exact and maximally rigid.

## Bounded stabilizer theorem

Let \(R_n\) be any of the three archived, oriented norm-two wall sets on

\[
H_{OG}=\begin{pmatrix}-2\phi&1\\1&2\end{pmatrix}
       \oplus E_8\oplus\langle2\rangle
\]

with \(n=171,321,457\).  Then

\[
\operatorname{Stab}_{O(H_{OG},\mathbb Z[\phi])}(R_n)=\{1\}.
\]

Here setwise stabilization includes the orientations of the wall
half-spaces.  Equivalently, an isometry must permute the recorded inward
norm-two roots, not merely their unoriented hyperplanes.

For a root \(r\in R_n\), form its exact fingerprint

\[
F_n(r)=\{\!\{(r,s):s\in R_n\}\!\}.
\]

Every setwise isometry preserves these multisets.  Exact arithmetic in
\(\mathbb Z[\phi]\), with no floating-point comparison, gives:

| walls | distinct fingerprints | largest multiplicity |
|---:|---:|---:|
| 171 | 171 | 1 |
| 321 | 321 | 1 |
| 457 | 457 | 1 |

Thus every root is fixed.  The first eleven archived roots are the ten
height-zero coordinate roots together with

\[
r_{10}=(\phi-1,1,0,\ldots,0).
\]

Their coordinate determinant is \(\pm(\phi-1)\), whose field norm is
\(-1\).  They are therefore a \(\mathbb Z[\phi]\)-basis of the entire
field lattice.  Fixing them forces the isometry itself to be the identity.
This checks both steps needed by the argument: the fingerprints use the
whole prefix, and the fixed roots recover the whole lattice integrally.
In particular, no infinite-order element stabilizes any one of these three
oriented finite wall configurations.

The ordered fingerprint digests for the three prefixes are, respectively,

```text
bc67833ef2f38bdfd43c5b2febbd5b18104e413d6d69eb83a71433e3be353e6b
f3d4c2c468abefd0404afd758e4886365b5f79b76f783fe5aef5fb79b0a08268
0358ead0c3e9457ef1ba68b3c8f67cdb2b542dbfe79af94609aacc7602b2ccbf
```

## Height-zero symmetries

The exact Gram automorphism group of the height-zero
\(A_1+E_8+A_1\) simple system has order two: the nontrivial element swaps
the two isolated \(A_1\) nodes.  Put

\[
p=(2,-1,0,\ldots,0).
\]

The height-zero span is \(p^\perp\), so an extension must send
\(p\mapsto\varepsilon p\), \(\varepsilon=\pm1\).  Since
\(e_0=(p-r_0)/2\), the isolated-node swap would send it to

\[
e_0-\frac12e_1+\frac12e_{10}
\quad\hbox{or}\quad
-e_0+\frac12e_1+\frac12e_{10},
\]

according to the sign of \(p\).  Neither extension is integral.  The
identity diagram map has two integral extensions: the identity and the
time-reversing involution \(p\mapsto-p\), whose characteristic polynomial
is \((t+1)(t-1)^{10}\).  Hence the integral height-zero stabilizer that also
preserves the selected ray \(\mathbb R_{>0}p\) is already trivial.

Scalar multiplication by a golden unit does not evade this conclusion:
\(uI\) preserves \(H_{OG}\) only if \(u^2=1\), so only \(u=\pm1\) occur,
and \(-I\) reverses the selected chamber orientation.

## An infinite-order negative witness

There are plenty of infinite-order lattice isometries; the obstruction is
chamber stabilization.  Let \(s_i\) be reflection in archived root \(r_i\).
The product

\[
T=s_0s_{10}
\]

is the identity on the last nine coordinates and has binary block

\[
T=\begin{pmatrix}
\phi&-\phi\\
-1-\phi&2\phi
\end{pmatrix}.
\]

It is an exact integral isometry with

\[
\chi_T(t)=(t-1)^9(t^2-3\phi t+1).
\]

At the identity embedding its spectral radius is

\[
\rho(T)=\frac{3\phi+\sqrt{9\phi^2-4}}2
       \approx4.638515763465525,
\]

so it has infinite order.  On the underlying rank-22 integer lattice its
characteristic polynomial is

\[
(t-1)^{18}(t^4-3t^3-7t^2-3t+1).
\]

This is a negative witness/search seed, not a chamber symmetry: among the
457 archived walls, only indices \(1,\ldots,9\), which it fixes, have images
remaining in the archived set.  The statement stays false even if wall
orientations are forgotten.  Thus this most natural golden hyperbolic
translation crosses the wall system rather than stabilizing its selected
prefix.

## Reproducibility

The exact generator is
[`scripts/certify_og10_golden_prefix_stabilizer.py`](../scripts/certify_og10_golden_prefix_stabilizer.py),
the frozen certificate is
[`data/candidates/og10_golden_prefix_stabilizer.json`](../data/candidates/og10_golden_prefix_stabilizer.json),
and the independent regression test is
[`tests/test_og10_golden_prefix_stabilizer.py`](../tests/test_og10_golden_prefix_stabilizer.py).

```sh
PYTHONPATH=.:src python3 scripts/certify_og10_golden_prefix_stabilizer.py
PYTHONPATH=.:src pytest -q tests/test_og10_golden_prefix_stabilizer.py
```
