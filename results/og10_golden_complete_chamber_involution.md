# A complete-chamber involution on the golden OG10 carrier

> **Subsequent decision.**  The anti-symplectic fixed H9 chamber has two
> complete-chamber involutions whose product has trace `8+8*phi>10` and hence
> infinite order.  This proves the fixed H9 chamber and the present H10 chamber
> infinite-volume, and the enhanced birational centralizer infinite.  See
> [`og10_golden_h9_infinite_stabilizer.md`](og10_golden_h9_infinite_stabilizer.md).

## Result

Let

\[
H=\begin{pmatrix}-2\phi&1\\1&2\end{pmatrix}
  \oplus E_8\oplus\langle2\rangle
\]

over \(\mathbb Z[\phi]\), and let \(C\) be the complete norm-two Vinberg
chamber selected by the control vector \(p=(2,-1,0,\ldots,0)\).  There is an
explicit nontrivial arithmetic involution

\[
h=\begin{pmatrix}1&1-\phi\\0&-1\end{pmatrix}
  \oplus I_9
\]

such that

\[
                         hC=C.                         \tag{1}
\]

This is a theorem about the **complete**, possibly infinitely sided chamber,
not a symmetry of only the archived 457-wall prefix.

The same matrix has two useful constructions.  If

\[
A=\begin{pmatrix}1&1-\phi\\-1&\phi\end{pmatrix}\oplus I_9,
\]

then \(A\) has infinite order, \(A^2=s_0s_{10}\), and the exact chamber
fold is

\[
                           h=s_0A.                     \tag{2}
\]

Equivalently, \(h\) is the integral reflection in

\[
\nu=(\phi-1,2,0,\ldots,0),\qquad H(\nu,\nu)=6+2\phi.  \tag{3}
\]

Thus the symmetry is not itself a norm-two wall reflection.  It acts
nontrivially on the discriminant group, so in particular \(h\notin W_2\).
The chamber stabilizer and the arithmetic reflection quotient therefore
contain a subgroup of order two.  Since \(h^2=1\), this does **not** prove
that the reflection index is infinite.

## Why the chamber statement is global

The global step uses a finite Lorentzian violator bound.  Put

\[
D=1+4\phi,
\quad
\rho=(0,\tfrac12,46,91,135,110,84,57,29,68,\tfrac12),
\quad
q=p+\frac{\rho}{100}.
\]

Then

\[
(p,\rho)=0,qquad \rho^2=621,qquad
q^2=-\frac{19379}{10000}-8\phi<0.
\]

The ten height-zero simple-root pairings at \(q\) are strictly negative.
For a positive-height norm-two simple root \(r\), write \(x=r_0>0\) and

\[
r=\frac{x}{2}p+z,qquad z\perp p.
\]

Since \(p^2=-2D\), one has

\[
z^2=2+\frac{D x^2}{2}.
\]

More generally, write a timelike test point as \(q'=ap+v\), with
\(v\perp p\).  If \((q',r)\geq0\), positive-definite Cauchy--Schwarz on
\(p^\perp\) gives

\[
aDx\leq(v,z)\leq
\sqrt{v^2\left(2+\frac{Dx^2}{2}\right)},
\]

and hence

\[
x^2\leq
\frac{2v^2}{D\left(a^2D-v^2/2\right)}.               \tag{4}
\]

The archived Vinberg computation is exhaustive for every simple root with
\(0<x\leq10\).  At the two endpoints used in (1), the exact right side of
(4) is respectively

\[
\beta_q=
\frac{2029167180}{49215477949}
-\frac{1186149744}{49215477949}\phi
\approx0.0022337806,
\]

and

\[
\beta_{hq}=
\frac{98378394011}{49215477949}
-\frac{98265338719}{98430955898}\phi
\approx0.3836204750.
\]

Both are positive and strictly below \(100\).  Direct exact evaluation gives
a strictly negative pairing with every one of the 457 archived roots at both
\(q\) and \(hq\).  Formula (4) then excludes every unarchived simple root:
any nonnegative pairing would force \(x<10\), where the enumeration is
already exhaustive.  Thus both points lie in the same complete open chamber.
Because \(h\) preserves the full norm-two arrangement, it follows that
\(hC=C\).

This argument is stronger than checking that a finite wall set is preserved.
Indeed, \(h\) does not preserve the control-height prefix.

## New walls transported by the involution

Applying \(h\) to the 457 archived simple roots gives:

| outcome | count |
|---|---:|
| image already among the 457 roots | 274 |
| image outside the archived prefix | 183 |
| involution-closed union | 640 |

The 183 transported roots occur at the following control coordinates:

| control coordinate | transported roots |
|---|---:|
| \(2+5\phi\) | 20 |
| \(4+4\phi\) | 3 |
| \(3+5\phi\) | 129 |
| \(4+5\phi\) | 31 |

All 640 roots have norm two and are pairwise acute, and the set is invariant
under \(h\).  Since \(h\) is now known to stabilize the complete chamber,
the 183 images are genuine complete-chamber simple walls.  They are **not**
claimed to exhaust those four height levels; the symmetry transports walls
but is not an exhaustive height enumerator.

The exact digests are

```text
new 183 roots:
3a0cdc70e695b452eaab5b119545420c042c389333626b5454d2b7c4b248ef22

closed 640-root set:
f4f4e5afbdc1eb8eaf1a64e3c3c0f23840668ec5c61623d0de237ae9b4755e8e
```

## Discriminant, monodromy, and birational scope

The involution is geometric in the precise enhanced sense needed here.
Let \(N\) be the rank-22 trace lattice and let \(n_1,n_2\) be the archived
gluing generators of orders \(2,22\) in \(A_N\).  Exact conversion of \(h\)
to the integral trace basis gives

\[
                       n_1\longmapsto n_1,
\qquad                 n_2\longmapsto-n_2.             \tag{5}
\]

On the orthogonal complement

\[
T=\langle2\rangle\oplus\langle66\rangle
\]

take \(-I_T\).  Its discriminant action matches (5), so
\(h\oplus(-I_T)\) preserves the order-44 gluing graph.  Conjugation into the
explicit integral overlattice basis verifies from scratch that the resulting
rank-24 matrix is integral, isometric, determinant one, and has order two.
It preserves the positive orientation: \(h\) preserves the chosen positive
line and \(-I_T\) has determinant \(+1\) on the positive complement plane.
It also preserves the OG10 period line, acting on it by \(-1\).

Consequently the extension is a positive Hodge-monodromy isometry by
[Onorati's computation of OG10 monodromy](https://arxiv.org/abs/2002.04129).
It commutes with the golden action on \(N\).  Combining (1) with the ambient
exceptional-chamber uniqueness argument and
[Markman's birational Torelli theorem](https://arxiv.org/abs/1101.4606)
shows

\[
                 h\in (\Gamma_\tau^{\mathrm{mon}})_C
                 =\operatorname{im}\bigl(\operatorname{Bir}(X,\tau)
                   \longrightarrow O(N)\bigr).          \tag{6}
\]

Thus the enhanced birational image contains a nontrivial involution.  More
precisely, (6) asserts the existence of an enhanced birational self-map whose
action on \(N\) is \(h\), and whose full \(H^2\)-action can be taken to be
the order-two extension above.  No assertion about the kernel of the
birational action is needed.

The discriminant calculation also proves \(h\notin W_2\), because every
norm-two reflection acts trivially on \(A_N\).  Hence this is a genuine
nontrivial class in the arithmetic reflection quotient, rather than a long
word in known wall reflections.

## What the seed search alone did not prove

The infinite-order Pell seed \(A\) folds to the finite-order symmetry \(h\).
Its inverse does the same.  A scan of fourteen additional distinct integral
reflections with normal supported on the binary block and final
\(\langle2\rangle\) line, using coefficients in the box
\([-1,1]+[-1,1]\phi\), folded globally to only \(1\) or \(h\).  Every final
membership test used (4) with bound below \(100\).

This bounded negative result by itself says nothing about whether a second
chamber symmetry exists.  At this stage it certified none of the following:

- no infinite-order same-chamber element has been certified;
- infinite reflection index has not been proved;
- finite reflection index has not been proved; and
- the complete chamber may still have infinitely many walls.

The sharp statement of this report is the order-two lower bound and its
geometric realization.  A subsequent fixed-slice construction does find the
needed second symmetry and proves infinite reflection index; it does not
invalidate the bounded ambient seed-search statement above.

## Reproducibility

The exact generator is
[`scripts/probe_og10_golden_chamber_folding.py`](../scripts/probe_og10_golden_chamber_folding.py),
the frozen certificate is
[`data/candidates/og10_golden_complete_chamber_involution.json`](../data/candidates/og10_golden_complete_chamber_involution.json),
and the independent regression test is
[`tests/test_og10_golden_complete_chamber_involution.py`](../tests/test_og10_golden_complete_chamber_involution.py).

```sh
PYTHONPATH=.:src python3 scripts/probe_og10_golden_chamber_folding.py
PYTHONPATH=.:src pytest -q tests/test_og10_golden_complete_chamber_involution.py
```

The frozen matrix digests are

```text
field involution h:
3682d8d751dedf54fdf429e4d35ef4286b8cec5eb68f37040cce5e6fbd052607

trace-lattice action:
74339270c5e8a77d2bd6dcadef2b124b59fc42872ec4ef043311c50e1e36e847

rank-24 OG10 extension:
e5b0029e0e54a4fa41cfea5a4ee36446af17d0d1c75b93d0e8b9c3829963c1df
```
