# The determinant-11 host: a known local lattice and a new symmetry proposal

**Status: exact comparison PASS; reflectivity remains undecided.** The focused
primary-source check did not establish a published exclusion for

\[
 L=E_8\perp\begin{pmatrix}-2\phi&-1\\-1&2\end{pmatrix},
 \qquad \phi=(1+\sqrt5)/2,\quad \beta=3+\phi.
\]

This is the repository's existing **H11 field lattice**, in hyperbolic
dimension nine: changing the sign of the second binary coordinate changes
the off-diagonal entries to +1. The name H11 refers to the discriminant norm,
not to hyperbolic dimension. The previous
[norm-two prefix](h11_vinberg_prefix.md) exhausts control coordinate at most
20, with 95 walls, and explicitly leaves reflectivity undecided. The
[finite-index audit](h11_h19_finite_index_obstruction.md) likewise supplies no
full-group nonreflectivity theorem. Including the primitive norm-\(2\beta\)
class is a substantive extension of that earlier norm-two calculation, not
the discovery of a different host.

## What the literature does and does not give

Bugaenko's 1984 paper explicitly concerns **unimodular** forms over this
ring. Its scope cannot be applied to this determinant-11 lattice.
[Primary bibliographic record](https://www.mathnet.ru/eng/vmumm3408).
The focused search did not verify a theorem in Bugaenko's 1992 chapter
covering this particular nonunimodular lattice; the
[publisher record](https://doi.org/10.1090/advsov/008/02) is not being treated
as access to the chapter's full theorem statements.

There is a close, relevant primary example. McLeod's thesis, §3.2.1,
printed pp.64–65, proves that the standard diagonal integral form
\(f_d^3=\langle-d,1,1,1\rangle\), \(d=3+2\sqrt5\), is nonreflective,
using an infinite-order chamber symmetry. Theorem 3.0.8, printed p.34,
requires a one-dimensional **unimodular integral orthogonal summand** to
pass a reflection chamber to a lower-dimensional face. Applied repeatedly
to the standard diagonal lattices, this rules out their reflectivity in
dimensions at least three. It does not assert the same for every lattice
in the rational quadratic space.
[McLeod, 2013](https://etheses.durham.ac.uk/id/eprint/7743/1/thesis31072013.pdf#page=72).

No newer applicable classification was identified by the focused search.
This records a search outcome, not a claim that no such theorem exists.

## Exact comparison of the integral lattices

The standalone [checker](../scripts/fresh_h9_beta11_literature_20261014.py)
and [certificate](fresh_h9_beta11_literature_20261014.json) use Fraction
arithmetic in \(\mathbb Q(\phi)\), with no project arithmetic imports.
They reconstruct the frozen host Gram and verify

\[
 \det L=-\phi\beta=-(3+2\sqrt5),\qquad N(\beta)=11.
\]

The ideal \((\beta)\) is prime: the residue map sends \(\phi\) to 8 in
\(\mathbb F_{11}\). An integral overlattice's index ideal squared must
divide the determinant ideal. The latter has a single prime factor of
valuation one, so **there is no proper integral overlattice of L**.
Likewise, rational similarity in rank ten changes determinant by a tenth
power and a square. Its odd valuation at \((\beta)\) excludes similarity
to any rank-ten integral form with unit determinant. This directly blocks
an attempted reduction to a unimodular host.

There is nevertheless an explicit integral sublattice \(M\subset L\)
with Gram \(2f_d^9\). In binary and Euclidean E8 coordinates, put

\[
 x=2z_0,\qquad t=z_0+z_1,\qquad
 z=\operatorname{diag}(H,H,H,H)(u_1,\ldots,u_8)^T,
 \quad H=\begin{pmatrix}1&1\\1&-1\end{pmatrix}.
\]

Then the norm is
\(-2d z_0^2+2z_1^2+2\sum u_i^2\). The certificate gives this map in the
literal E8 simple-root basis, verifies its determinant \(\pm32\), and
therefore computes absolute index \([L:M]=1024\). Thus the rational form
is similar to McLeod's form, but the integral lattices differ at 2.

This distinction affects mirrors. The E8 root \(r=e_1-e_3\) has norm two
and its reflection preserves L. It sends the M-generator \(e_1+e_2\) to
\(e_2+e_3\), whose M-coordinates contain halves. Hence that reflection
does not preserve M. The explicit embedding therefore supplies no
inclusion of the full reflection group of L in the nonreflective standard
diagonal integral group. Neither this example nor rational similarity
decides the reflection index of L.

## Executable symmetry candidate

The certificate transcribes the 4×4 matrix in McLeod's Appendix B.2.1,
printed p.108, then checks it independently.
[Primary matrix](https://etheses.durham.ac.uk/id/eprint/7743/1/thesis31072013.pdf#page=116).
Let it be \(g_4\), and let T be the preceding embedding. The exact checks are

\[
 g_4^T f_d^3 g_4=f_d^3,\quad \det g_4=1,\quad
 \operatorname{tr}g_4=80+128\phi>4,
 \qquad g=T(g_4\oplus I_6)T^{-1}.
\]

All entries of **g itself** lie in \(\mathbb Z[\phi]\); no power is
needed. The checker verifies \(g^T G_Lg=G_L\), determinant one, and future
orientation. Its trace bound proves infinite projective order. Both the
literal source matrix and the transported 10×10 matrix are serialized as
coefficient pairs \([a,b]=a+b\phi\). The transported order is
\((x,E8_0,\ldots,E8_7,t)\).

**Chamber stabilization has not been checked.** The extension by six
coordinates and the extra mirrors of L both require a fresh test. A
complete separator check must include both primitive norm classes, 2 and
\(2\beta\). If the image is folded back by reflections, the resulting
matrix must be tested again for infinite order; folding can remove the
useful symmetry. This is an explicit next experiment, not a
nonreflectivity certificate.

Reproduce with

```sh
python scripts/fresh_h9_beta11_literature_20261014.py --verify
```

Generation and saved replay pass. Frozen predecessor sources, data and
reports were not modified.
