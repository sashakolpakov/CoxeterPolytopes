# A projective OG10 carrier for a golden hyperbolic 10-space

> **Subsequent completion.**  The wall-completeness and chamber-selection
> gaps left open in this carrier note are now closed in
> [the complete wall and movable-cone theorem](og10_golden_complete_wall_chambers.md):
> the carrier has no primitive divisibility-three classes, its complete
> Kähler arrangement has three field-norm orbits, and its marked
> eigen-movable cone is exactly a norm-two Coxeter chamber.  Its finite-volume
> problem is equivalent to finiteness of an enhanced birational centralizer;
> see [the covolume theorem](og10_enhanced_covolume_theorem.md).  The problem
> has subsequently been decided negatively: an infinite fixed-slice chamber
> stabilizer proves the H10 chamber infinite-volume and the enhanced
> birational centralizer infinite; see
> [the infinite-stabilizer theorem](og10_golden_h9_infinite_stabilizer.md).

## Result

There exists a projective irreducible holomorphic symplectic manifold `X` of
OG10 deformation type whose Neron--Severi lattice is

\[
N\cong
U\oplus T_{11}(-1)\oplus E_8(-1)^2\oplus A_1(-1)^2,
\qquad
T_{11}=\begin{pmatrix}2&1\\1&6\end{pmatrix},
\]

and on which `N` has an integral self-adjoint endomorphism

\[
\tau^2-\tau-1=0.
\]

The `phi`-eigenspace of `tau` has signature `(1,10)`, while the conjugate
eigenspace is negative definite of rank 11.  Its positive cone therefore
supplies a genuinely algebro-geometric carrier for hyperbolic 10-space.

This crosses the rank barrier that prevents the same full-NS construction on
a K3 surface.  It does **not** yet produce a Coxeter polytope: the next gate is
to determine the OG10 MBM or prime-exceptional walls and prove that their
intersection with the `phi`-eigenspace is a complete acute finite-volume or
compact Coxeter chamber.

## The rank-11 golden form

Put

\[
\mathcal O=\mathbb Z[\phi],\qquad \phi^2-\phi-1=0,
\]

and consider the symmetric `O`-lattice with Gram matrix

\[
H=
\begin{pmatrix}-2\phi&1\\1&2\end{pmatrix}
\oplus E_8\oplus\langle2\rangle.
\]

At the identity real embedding its signature is `(10,1)`; at the conjugate
embedding it is positive definite of rank 11.  Define the integral trace
pairing

\[
B(x,y)=\operatorname{Tr}_{\mathbb Q(\sqrt5)/\mathbb Q}
       \left(-\frac{\phi}{\sqrt5}H(x,y)\right).
\]

The scalar is in the codifferent of `O`.  Direct exact reduction gives

\[
(\mathcal O^{11},B)\cong N.
\]

In particular,

\[
\operatorname{rk}N=22,
\qquad \operatorname{sig}N=(1,21),
\qquad \det N=-44.
\]

Multiplication by `phi` is integral, satisfies its golden polynomial, and is
self-adjoint for `B`.  The scalar `-phi/sqrt(5)` is negative at both real
embeddings.  It reverses the `(10,1)` identity signature to `(1,10)` and the
positive-definite conjugate signature to `(0,11)`.  Thus precisely one real
eigenspace is hyperbolic.

## Primitive embedding in the OG10 Beauville lattice

The unimodular summand

\[
C=U\oplus E_8(-1)^2
\]

plays no role in the discriminant gluing.  Write the remaining part as

\[
K=-T_{11}\oplus\langle-2\rangle^2
=
\begin{pmatrix}
-2&-1&0&0\\
-1&-6&0&0\\
0&0&-2&0\\
0&0&0&-2
\end{pmatrix},
\]

and take the positive rank-two lattice

\[
T=\langle2\rangle\oplus\langle66\rangle.
\]

In coordinates for `K + T`, adjoin the two dual vectors

\[
g_1=\left(6,-1,\frac{11}{2},\frac12;
                 \frac12,\frac12\right),
\]

\[
g_2=\left(-\frac6{11},\frac1{11},-\frac12,0;
                 0,\frac9{22}\right).
\]

Modulo `K + T`, they generate a subgroup of order `2*22=44`.  Projection to
each discriminant group is injective.  The exact quadratic data are

\[
q(g_1)=-110,\qquad q(g_2)=10,
\qquad (g_1,g_2)=25,
\]

so their graph is isotropic.  The resulting overlattice has determinant

\[
\frac{\det(K)\det(T)}{44^2}
=\frac{44\cdot132}{44^2}=3.
\]

This determinant calculation alone is not used as an identification.  An
explicit basis of the overlattice, in the original `K + T` coordinates, is

\[
R=
\begin{pmatrix}
1&0&0&0&0&3/11\\
0&1&0&0&0&5/11\\
0&0&1&0&0&1/2\\
0&0&0&1&1/2&0\\
0&0&0&0&1/2&0\\
0&0&0&0&0&1/22
\end{pmatrix}.
\]

Its determinant is `1/44`, it contains `K + T`, `g1`, and `g2`, and its Gram
matrix is

\[
Q=R^t(K\oplus T)R=
\begin{pmatrix}
-2&-1&0&0&0&-1\\
-1&-6&0&0&0&-3\\
0&0&-2&0&0&-1\\
0&0&0&-2&-1&0\\
0&0&0&-1&0&0\\
-1&-3&-1&0&0&-2
\end{pmatrix}.
\]

Finally, the unimodular matrix

\[
P=
\begin{pmatrix}
0&-1&0&0&0&0\\
-1&-1&0&0&-1&-1\\
-1&-1&0&0&-1&0\\
0&0&0&-1&0&0\\
0&0&1&1&0&0\\
2&2&0&0&1&1
\end{pmatrix},
\qquad \det P=-1,
\]

satisfies

\[
P^tQP=U^2\oplus A_2(-1).
\]

Consequently the full glued lattice is explicitly

\[
C\oplus Q
\cong3U\oplus E_8(-1)^2\oplus A_2(-1).
\]

This is the OG10 Beauville lattice computed by Rapagnetta
([arXiv:math/0606409](https://arxiv.org/abs/math/0606409)).  The displayed
overlattice basis also shows directly that `N` is primitive and that its full
orthogonal complement is exactly `T`.

As a redundant discriminant-form check, the orthogonal residual to the image
of `A_N` in

\[
A_T\cong\mathbb Z/2\oplus\mathbb Z/66
\]

is

\[
\{(0,0),(0,22),(0,44)\}.
\]

The two nonzero elements have quadratic value `4/3 mod 2Z`, which is the
discriminant form of `A2(-1)`.

## Geometric realization

Let

\[
\Lambda_{OG10}=3U\oplus E_8(-1)^2\oplus A_2(-1).
\]

The construction above gives a primitive embedding

\[
N\hookrightarrow\Lambda_{OG10},
\qquad N^\perp=T=\langle2\rangle\oplus\langle66\rangle.
\]

Choose either orientation of the positive plane `T_R`; for example the period
line may be represented in a diagonal basis by

\[
\sigma=\sqrt{33}\,e_1+i e_2,
\qquad (\sigma,\sigma)=0,
\quad (\sigma,\bar\sigma)>0.
\]

Marked period-map surjectivity realizes this period by a manifold of OG10
deformation type; see Huybrechts' basic-results paper
([arXiv:alg-geom/9705025](https://arxiv.org/abs/alg-geom/9705025)) and the
global Torelli framework of Verbitsky
([arXiv:0908.4121](https://arxiv.org/abs/0908.4121)).  Since the real and
imaginary parts of `sigma` span all of `T_R`, its integral `(1,1)` lattice is
exactly `N`, not a proper overlattice.

The lattice `N` contains a positive integral vector, for example the sum of
the two standard isotropic generators of `U`.  Huybrechts' corrected
projectivity criterion then makes the resulting manifold projective
([arXiv:math/0106014](https://arxiv.org/abs/math/0106014)).

Thus the construction produces an algebraic variety with the prescribed full
divisor lattice, rather than merely a rational Hodge structure or an abstract
hyperbolic lattice.

## An automatic wall lift for every field root

The OG10 wall theorem supplies more than realization of the carrier.  Let

\[
r\in\mathcal O^{11},\qquad H(r,r)=2,
\]

be any norm-two root of the field lattice, and regard `r` as an integral
divisor `D` in the trace lattice `N=NS(X)`.  In particular, both `D` and
`tau D` are integral `(1,1)` classes.  Exact trace evaluation gives

\[
B(D,D)=-2,
\qquad B(\tau D,\tau D)=-4,
\qquad B(D,\tau D)=-2.
\]

The class `D` is primitive because an even lattice cannot contain a vector
`D0` with `-2=k^2 D0^2` for an integer `k>1`.  Multiplication by `phi` is an
integral automorphism, with inverse `tau-1`, so `tau D` is primitive as well.
Primitivity holds in the ambient Beauville lattice because `N` is primitively
embedded.  Its discriminant group has order three, so the ambient divisibility
of a primitive vector is one or three; divisibility three would divide its
square, which is impossible for `-2` and `-4`.  Thus both classes have
divisibility one.
Mongardi--Onorati's classification therefore says that `D` and `tau D` are
wall divisors: they are precisely the `(-2,div=1)` and `(-4,div=1)` wall
types.

These two ambient walls restrict to the same wall on the golden eigenspace.
Indeed, for

\[
x\in V_\phi=\ker(\tau_{\mathbb R}-\phi),
\]

self-adjointness gives

\[
(\tau D,x)=(D,\tau x)=\phi(D,x).
\]

Hence

\[
D^\perp\cap V_\phi=(\tau D)^\perp\cap V_\phi.
\]

Under the trace-field identification this is exactly the hyperplane
`H(r,-)=0`.  Thus:

> Every norm-two reflection mirror of the rank-11 golden lattice occurs in
> the restriction of the numerical OG10 wall-divisor arrangement to the
> hyperbolic eigenslice.

This is a necessary bridge that was missing from the abstract lattice
construction.  It turns the next computation into a precise Vinberg problem
on `H`, rather than a search through unrelated varieties.  It is not yet a
completeness theorem: OG10 also has `(-6,div=3)` and `(-24,div=3)` walls, and
some `(-2)` or `(-4)` divisor walls not arising from a norm-two field root may
cut the eigenslice further.  A wall divisor is a monodromy-stable numerical
wall class; the statement also does not yet prove that each such restricted
hyperplane is a facet of the Kähler cone of this particular `X`.  We have not
yet selected a birational/marked model whose Kähler cone meets the desired
field chamber.

## What is and is not new

The lattice and period-map argument give a new positive output for the
algebro-geometric search program:

- the complex K3 full-NS golden factory stops at `H^9`;
- an OG10-type variety realizes a full rank-22 golden divisor lattice with an
  `H^10` eigenslice; and
- the construction is explicit enough to start a constrained wall search on
  the exact field form `H`.

It is not yet a new hyperbolic Coxeter polytope.  In particular, multiplication
by `phi` is an endomorphism of the divisor lattice; no claim is made here that
it is induced by an automorphism of `X`.  Nor has the MBM-wall set been shown
to be invariant under it.  The required next theorem is a level-3 statement:
classify the OG10 wall-divisor classes in this particular `N`, intersect them
with the hyperbolic eigenspace, and certify completeness and compactness of the
resulting chamber.

The relevant geometric input for that step is the characterization of Kahler
classes for OG10-type manifolds by Mongardi--Onorati
([arXiv:2010.12511](https://arxiv.org/abs/2010.12511)).

## Reproducibility

Run the exact certificate and focused test from the repository root:

```sh
PYTHONPATH=.:src python3 scripts/certify_og10_golden_h10_carrier.py
PYTHONPATH=.:src pytest -q tests/test_og10_golden_h10_carrier.py
```

The frozen output is
`data/candidates/og10_golden_h10_carrier.json`; the test requires exact
equality between that artifact and a fresh certificate run.

The certificate checks the field signatures, trace-basis isometry, integral
self-adjoint golden action, order-44 isotropic gluing graph, primitive
embedding, orthogonal complement, residual discriminant form, and the explicit
isometry to the OG10 Beauville lattice.  The period-map and projectivity
theorems are mathematical inputs, not replaced by computation.
