# TOOT: the simple golden norm-six carrier embeds in OG10

## The positive theorem

Let

\[
\mathcal O=\mathbb Z[\phi],\qquad \phi^2-\phi-1=0,
\]

and put

\[
F=\langle-2\phi\rangle\oplus E_8\oplus\langle2\rangle,
\qquad k=(0^9,1).
\]

Adjoin a vector `e` with norm six and pairing `3k`:

\[
H_6(k)=
\begin{pmatrix}
F&3k\\
3k^t&6
\end{pmatrix}.
\]

Its signatures at the two real embeddings are `(10,1)` and `(11,0)`, and
its determinant is `-6 phi`.  The restriction-of-scalars trace form

\[
B(x,y)=\operatorname{Tr}_{\mathbb Q(\sqrt5)/\mathbb Q}
\left(-\frac{\phi}{\sqrt5}H_6(k)(x,y)\right)
\]

is an even lattice `N` of rank 22, signature `(1,21)`, determinant `-36`,
and discriminant group

\[
A_N\cong(\mathbb Z/6)^2.
\]

There is an explicit primitive embedding

\[
N\hookrightarrow
\Lambda_{OG10}=3U\oplus2E_8(-1)\oplus A_2(-1)
\]

whose orthogonal complement is

\[
T=\langle2\rangle\oplus\langle6\rangle.
\]

For the two trace lifts

\[
D=e,\qquad E=(\phi-1)e,
\]

one has

\[
B|_{\langle D,E\rangle}=\begin{pmatrix}-6&0\\0&-6\end{pmatrix},
\]

and the embedding has been chosen so that

\[
\operatorname{div}_{\Lambda_{OG10}}(D)=3,
\qquad
\operatorname{div}_{\Lambda_{OG10}}(E)=1.
\]

The last row of `H_6(k)` lies in `3 O`, so

\[
s_e(x)=x-\frac{H_6(k)(x,e)}3e
\]

is an integral field-lattice reflection.  On the identity real eigenspace,

\[
B(e,x)=-\frac{\phi}{\sqrt5}H_6(k)(e,x),
\]

so the restriction of the ambient hyperplane `D^perp` is exactly the field
reflection mirror `e^perp`.

Consequently the primitive class `D` is an OG10 wall divisor of type
`(-6,3)`.  This closes the primitive-embedding gate left open in the earlier
norm-six reduction note.

This is a genuine positive **algebro-geometric carrier and wall theorem**.
It is not yet a new finite-volume or compact hyperbolic Coxeter polytope.

## Splitting the trace lattice

An integral change of basis puts the trace Gram matrix in the form

\[
N\cong2E_8(-1)\oplus S,
\]

where

\[
S=
\begin{pmatrix}
 2& 4& 0& 0& 0& 0\\
 4& 6& 0& 0& 0& 0\\
 0& 0&-2&-2&-3&-3\\
 0& 0&-2&-4&-3&-6\\
 0& 0&-3&-3&-6&-6\\
 0& 0&-3&-6&-6&-12
\end{pmatrix}.
\]

The change-of-basis columns are the first eight `E8` coordinate vectors,
then their `(phi-1)` multiples, followed by the two trace coordinates of
each of the `x`, `u`, and `e` field directions.  Its determinant is one.
The Smith diagonal of `N` is

\[
(1^{20},6,6).
\]

## The selected discriminant gluing

Work in the displayed basis of

\[
S\oplus T,
\qquad T=\operatorname{diag}(2,6).
\]

Adjoin

\[
g_1=(0,\tfrac12,0,0,0,0;\tfrac12,0),
\]

\[
g_2=(\tfrac12,0,0,0,-\tfrac13,\tfrac13;0,\tfrac16).
\]

Modulo `S + T`, these vectors have orders two and six.  Their exact
quadratic data in `S + T` are

\[
g_1^2=2,\qquad g_2^2=0,\qquad(g_1,g_2)=1.
\]

Thus they span an order-12 isotropic graph, with injective projections to
both discriminant groups.  The complete 36-class enumeration of `A_N`, with
the standard generators of `A_T` fixed, gives four anti-isometry generator
pairs.  Exactly two are orthogonal to `D/3`; they differ by the sign of the
order-six complement generator.  The chosen graph is one of these two.

If

\[
R=[e_1,\ldots,e_6,g_1,g_2],
\]

then `det(R)=1/12`, `R` contains the original `S + T` and both gluing
generators, and

\[
Q=R^t(S\oplus T)R=
\begin{pmatrix}
2&4&0&0&0&0&2&1\\
4&6&0&0&0&0&3&2\\
0&0&-2&-2&-3&-3&0&0\\
0&0&-2&-4&-3&-6&0&-1\\
0&0&-3&-3&-6&-6&0&0\\
0&0&-3&-6&-6&-12&0&-2\\
2&3&0&0&0&0&2&1\\
1&2&0&-1&0&-2&1&0
\end{pmatrix}.
\]

The ambient determinant is therefore `-3`, but the identification with the
OG10 lattice does not rely on determinant and signature alone.

## Literal isometry to the OG10 Beauville lattice

The unimodular matrix

\[
P=
\begin{pmatrix}
0&1&-1&1&2&1&1&2\\
0&0&0&-1&0&0&0&0\\
0&0&0&0&-1&-1&0&-3\\
0&0&0&0&0&1&1&2\\
0&0&0&0&0&0&0&1\\
0&0&0&0&1&0&0&0\\
0&0&1&1&0&0&0&0\\
1&-1&0&0&-4&-2&-2&-4
\end{pmatrix},
\qquad\det P=-1,
\]

satisfies

\[
P^tQP=3U\oplus A_2(-1).
\]

Together with the two untouched `E8(-1)` blocks, this gives a literal
24-dimensional isometry with `Lambda_OG10`.  More explicitly, let `C_N` be
the trace splitting above and set

\[
A=\operatorname{diag}(C_N,I_2)\operatorname{diag}(I_{16},R).
\]

If `Pi` moves the three `U` blocks ahead of the two `E8(-1)` blocks, then

\[
Q_{24}=\operatorname{diag}(I_{16},P)\Pi
\]

is integral and unimodular, and the exact certificate verifies

\[
Q_{24}^t\bigl(A^t(N\oplus T)A\bigr)Q_{24}
=3U\oplus2E_8(-1)\oplus A_2(-1).
\]

The certificate stores `A`, `Q_24`, `A Q_24`, both primitive embedding
matrices, and every Gram matrix as literal rational or integral matrices.
The two inclusion matrices have Smith diagonals `(1^22)` and `(1,1)`, so
both `N` and `T` are primitive; their mutual pairing is identically zero.

## The selected `(-6,3)` direction

In the glued ambient basis, the two classes have coordinates

\[
D=(0^{20},1,0,0,0),
\qquad
E=(0^{20},-1,1,0,0).
\]

Their ambient pairing rows are

\[
(D,-)=(0^{18},-3,-3,-6,-6,0,0),
\]

\[
(E,-)=(0^{19},-3,0,-6,0,-2).
\]

Hence both vectors are primitive and have square `-6`, while the row gcds
are respectively three and one.  Equivalently, `D/3` is orthogonal to the
gluing graph and survives as the nonzero order-three ambient discriminant
class, whereas `E/3` does not.

Mongardi--Onorati's
[Proposition 3.1](https://arxiv.org/abs/2010.12511) says that a multiple of a
primitive OG10 class of square `-6` and divisibility three is stably prime
exceptional.  Their Proposition 5.4 says directly that such a primitive
class is a wall divisor.  Thus the hyperplane `D^perp` is a genuine numerical
OG10/MBM wall; this statement does not assert that a chosen sign of `D` is
effective on the particular marked model.

## Algebraic realization

The complement `T` is positive definite.  In its diagonal basis choose

\[
\sigma=t_1+i\,t_2/\sqrt3.
\]

Then

\[
(\sigma,\sigma)=0,
\qquad(\sigma,\bar\sigma)=4.
\]

Period-map surjectivity realizes this line by a marked OG10-type manifold.
Because the real and imaginary parts span all of `T_R`, and because `N` is
primitive, the full integral `(1,1)` lattice is exactly `N`.  The first
trace coordinate has square two, so the positive integral `(1,1)` criterion
makes the manifold projective.  References for these standard external
inputs are Huybrechts' period and deformation paper
([arXiv:alg-geom/9705025](https://arxiv.org/abs/alg-geom/9705025)) and his
projectivity criterion
([arXiv:math/0106014](https://arxiv.org/abs/math/0106014)).  The OG10
Beauville lattice is the one computed by Rapagnetta
([arXiv:math/0606409](https://arxiv.org/abs/math/0606409)).

Therefore there exists a projective OG10-type manifold `X` with

\[
NS(X)=N
\]

and with `D` as a primitive `(-6,3)` wall class.  This is an algebraic variety
whose divisor lattice was designed from the desired Coxeter/reflection
constraint, rather than found by an unconstrained search.

## What this does and does not answer

This is the first positive resolution of the precise norm-six embedding gate:

- the rank-11 golden field carrier exists at both real places;
- its full rank-22 trace lattice is the Picard lattice of a projective OG10;
- the gluing deliberately selects `D=e`, rather than the orthogonal golden
  partner `E`, as the ambient divisibility-three wall direction; and
- the selected field mirror is therefore represented by a genuine OG10
  wall hyperplane.

It does **not** prove that the resulting wall arrangement has a finite-sided,
finite-volume, or compact chamber.  No new hyperbolic Coxeter polytope is
claimed.  Those conclusions still require a complete enumeration of every
OG10 wall meeting the golden eigenslice, chamber selection, and Vinberg's
finite-volume or compactness tests.

## Reproduction

The frozen artifacts are:

- `scripts/certify_og10_golden_norm_six_embedding.py`;
- `data/candidates/og10_golden_norm_six_embedding.json`;
- `tests/test_og10_golden_norm_six_embedding.py`.

Run

```bash
PYTHONPATH=src:. python3 \
  scripts/certify_og10_golden_norm_six_embedding.py --verify-certificate
PYTHONPATH=src:. pytest -q tests/test_og10_golden_norm_six_embedding.py
```

The verifier recomputes the two field signatures, trace signature and Smith
form, all discriminant classes and anti-isometries, the selected gluing,
ambient Gram matrix, literal 24-by-24 OG10 isometry, primitive embedding
Smith forms, orthogonal complement, and the square/divisibility of both
selected trace lifts before comparing the frozen payload and source hashes.
