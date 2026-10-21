# Golden odd-Leech section in H23: exact first packets

## Scope

Let

\[
S_{23}=(1^{24})^\perp\cap O_{24},\qquad
L_a=\langle-a\rangle\oplus(S_{23}\otimes_{\mathbb Z}\mathbb Z[\phi]),
\]

where \(\phi^2=\phi+1\).  The already audited section has determinant
\(3\), Smith form \(1^{22},3\), minimum \(3\), and 2,576 minimal vectors.
This certificate treats the two elementary admissible coefficients
\(a=\phi\) and \(a=2\phi\).  At the conjugate embedding both Lorentzian
coefficients change sign, so the conjugate forms are positive definite.

This is a complete calculation only through the first nonempty
positive-height packet for each coefficient.  It is not a terminal Vinberg
run and it is not a census of every coefficient ideal.

## Local mirror norms and height zero

The primes \(2\) and \(3\) are inert in \(\mathbb Z[\phi]\).  Reflection
integrality, primitivity, the Smith form of \(S_{23}\), and normalization by
totally-positive unit squares give the complete canonical norm alphabet

\[
q\in\{1,2,3,6\}
\]

for both coefficients.  When \(3\mid q\), both spatial coefficients lie in
the one-dimensional pairing radical modulo 3.  Its certified minimum is at
least 7.

Writing \(y=u+\phi v\), a height-zero root satisfies

\[
u^2+v^2=q,\qquad 2(u,v)+v^2=0.
\]

The ungated necessary Gram-type counts for \(q=1,2,3,6\) are
\(0,0,1,1\); after the radical gate they are \(0,0,0,0\).  Thus height zero
is exactly empty.

## Exact section shells

The complete norm-4 section shell has

\[
759\binom84+24\cdot23=53{,}682
\]

vectors.  These are the sign-balanced even octad vectors and the ordered
\((+4,-4)\) raw-coordinate vectors; the odd coset contributes none.

The complete norm-5 section shell lies in the odd coset and has

\[
759\binom{16}{2}+2576\cdot12^2+759\binom{16}{2}
=553{,}104
\]

vectors.  The three terms come from Golay weights 8, 12, and 16.  The JSON
certificate records the full parity and zero-sum derivation.

## First packet for \(a=\phi\)

There are exactly seven admissible control levels through the first
nonempty one.  The first six have no necessary spatial Gram type.  The last
has

\[
h=\frac32+\frac52\phi,\quad q=2,\quad x=1+\phi,\quad
y=\phi v,\quad v^2=5.
\]

It therefore contains exactly 553,104 primitive reflecting roots.  For two
distinct norm-5 vectors, the section minimum gives \((v,w)\le3\), hence the
corresponding root pairing is at most \(-2\phi<0\).  The whole packet is
accepted simultaneously.

Let \(e=(4,-4,0^{22})/\sqrt8\), so \(e^2=4\), and put
\(z=(1,\frac23e)\).  Then

\[
z^2=\frac{16}{9}-\phi>0,
\qquad (z,r)\le-1<0
\]

for every first-packet wall.  Thus this exact prefix has a strict
hyperideal patch and is neither finite-volume nor compact.

## First packet for \(a=2\phi\)

There are exactly two admissible control levels through the first nonempty
one.  The first is empty.  The second has

\[
h=1+2\phi,\quad q=2,\quad x=\phi,\quad
y=\phi v,\quad v^2=4.
\]

It contains exactly 53,682 primitive reflecting roots.  Distinct norm-4
vectors have inner product at most 2, again giving root-pairing upper bound
\(-2\phi<0\).

For any norm-3 section vector \(b\), put \(z=(1,\frac43b)\).  Since a
norm-3 vector pairs at most 2 with a norm-4 vector,

\[
z^2=\frac{16}{3}-2\phi>0,
\qquad (z,r)\le-2+\frac23\phi<0.
\]

This prefix also has a strict hyperideal patch and is neither finite-volume
nor compact.

## Exact conclusion

Neither natural golden coefficient produces a compact or finite-volume
chamber at its first complete packet.  This is a bounded negative result:
later mirrors may cut the displayed spacelike points, and other coefficient
ideals remain outside this calculation.

Reproduce with:

```bash
PYTHONPATH=src python3 scripts/audit_golden_odd_leech_section_h23.py --verify
PYTHONPATH=src pytest -q tests/test_golden_odd_leech_section.py
```
