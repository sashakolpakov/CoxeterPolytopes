# Direct H9 Vinberg continuation through \(3+6\phi\)

## Result

For the golden fixed lattice

\[
M=\langle-2\phi\rangle\oplus E_8\oplus\langle2\rangle,
\qquad \operatorname{sig}_{\mathrm{id}}(M)=(9,1),
\]

the direct norm-two Vinberg search is exhaustive through control coordinate
\(3+6\phi\).  The new level contains exactly **ten** simple roots.  All ten
are strict on the previously certified spacelike extreme ray.  Consequently
the resulting 71-wall outer cone still has the exact hyperideal ray

\[
v=\left(
1,2\phi,4\phi,6\phi,
\frac{64\phi}{13},\frac{48\phi}{13},\frac{32\phi}{13},
\frac{16\phi}{13},\frac{40\phi}{13},\frac{\phi}{13}
\right)
\]

with

\[
v^2=\frac{210-128\phi}{169}>0.
\]

Thus the norm-two outer cone is not of finite volume at this bounded stage.
Later simple roots are not excluded and could, by themselves, truncate this
ray.  The statement here is a bounded exact Vinberg theorem, not a second
proof of completeness or infinite covolume of the full chamber.

## Direct rank-ten recursion

Write a vector of \(M\) as \(r=(x,z,u)\), and put \(w=E_8z\).  Its norm-two
equation, doubled to give a convenient exhaustive coordinate box, is

\[
2w^TE_8^{-1}w+4u^2=4+4\phi x^2.                 \tag{1}
\]

The height-zero chamber is the dominant chamber for \(E_8\oplus A_1\), so
\(w_i,u\geq0\).  If \(R=(X,Z,U)\) is an earlier simple root, the Vinberg
condition is

\[
2w\mathbin{\cdot}Z+4uU\leq4\phi xX.             \tag{2}
\]

This gives a nine-variable positive-definite Schur recursion, rather than
the ten positive variables used by the ambient H10 enumeration.  The
successive conditional diagonal coefficients of the form in (1) are

\[
1,\ \frac43,\ \frac32,\ \frac85,\ \frac53,\
\frac{12}{7},\ \frac74,\ 16,\ 4.
\]

They are all at least one.  Hence testing every dominant ring coordinate
whose square is at most the right side of (1), at both real embeddings, is
an exhaustive coordinate alphabet.  Schur remainders are tested exactly at
both embeddings, and monotone partial sums of (2) provide the wall pruning.

At the conjugate embedding, the original norm equation gives

\[
2(\phi-1)(x')^2\leq2,
\]

and therefore

\[
(x')^2\leq\phi.                                  \tag{3}
\]

Together with \(0<x\leq13\), (3) gives exactly the following candidate
coordinates, in Vinberg-height order:

\[
\begin{split}
&1,\ \phi,\ 1+\phi,\ 2\phi,\ 1+2\phi,\ 2+2\phi,
1+3\phi,\ 2+3\phi,\\
&3+3\phi,\ 2+4\phi,\ 3+4\phi,\ 2+5\phi,
3+5\phi,\ 4+5\phi,\ 3+6\phi.
\end{split}
\]

Thus \(3+6\phi\) is the unique unprocessed candidate below identity bound
13.  The next admissible coordinate is \(4+6\phi\).

## From-scratch calibration

Before trusting the continuation, the direct enumerator was replayed from
the nine height-zero \(E_8\oplus A_1\) walls.  Through \(4+5\phi\), its
positive-height root counts were

\[
2,2,2,0,4,3,2,7,0,6,7,1,13,3.
\]

Together with the nine cone walls, these are 61 roots.  Their complete
ordered tuple—not merely their counts or set—is identical to the 61 fixed
roots extracted from the independently certified 1144-wall ambient OG10
prefix.  Both ordered lists have SHA-256 digest

```text
b83e9e3018debdccc15a9934e2ab735e9882f643f1f04ab155bb52f5a881eabc
```

This calibration simultaneously checks the H9-specific conjugate filter,
height ordering, positive-definite recursion, wall pruning, root
reconstruction, and orientation against an independent enumeration.

## The \(3+6\phi\) level

The exhaustive recursion visited 156,899 nodes and accepted ten roots.  Its
depth profile was

```text
1, 49, 784, 2793, 4802, 8134, 14504, 27832, 65170, 32830
```

The ten roots have digest

```text
cb22d2c82dd89c8b423a5a5e2cc78b31488b836c3d12aa642d9dd7011582498f
```

Their pairings with the old hyperideal ray, in enumerator order, are

\[
\frac1{13}\bigl(
-22-32\phi,\ -42-42\phi,\ -10+6\phi,\ -22-12\phi,
-32-32\phi,
\]

\[
-32-22\phi,\ -16,\ -16-6\phi,\ -22-22\phi,
-18-12\phi\bigr).
\]

Every value is strictly negative at the identity embedding.  The old nine
active facets therefore remain active, and the independent extreme-ray
verifier reproduces exactly the same ray and positive norm against all 71
walls.

As a separate cross-check, the complete-chamber Pell-fold involution \(g\)
maps four of the archived 61 walls to roots at this level.  These occur at
new-root indices \(1,4,5,8\), coming respectively from old-root indices
\(37,51,59,60\).  All four appear in the direct output, while the direct
search finds six further walls not forced by that symmetry.

## Reproducible artifacts

- Enumerator: `src/coxeter_polytopes/search/og10_golden_h9_vinberg.py`
- Certificate builder:
  `scripts/certify_og10_golden_h9_vinberg_through_3_plus_6phi.py`
- Frozen certificate:
  `data/candidates/og10_golden_h9_vinberg_through_3_plus_6phi.json`
- Regression test:
  `tests/test_og10_golden_h9_vinberg_through_3_plus_6phi.py`

The frozen JSON includes every one of the 71 roots, all level-recursion
counters, root and pairing digests, the Pell-fold cross-check, and the
independently verified hyperideal-ray certificate.
