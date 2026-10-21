# K3 completion theorem for the eight forced golden Lanner seeds

## Outcome

Let \(L_\Gamma\) be the even trace lattice obtained from one of the eight
norm-two golden Lanner blocks classified in
results/golden_lanner_global_constraints.md. Then:

1. every \(L_\Gamma\) has a primitive embedding in the K3 lattice
   \(\Lambda_{K3}=U^3\oplus E_8(-1)^2\), unique up to
   \(O(\Lambda_{K3})\);
2. no \(L_\Gamma\) has a proper even overlattice preserved by the golden
   endomorphism \(\tau\);
3. consequently, \(L_\Gamma\) is automatically primitive in every even
   \(\tau\)-stable Néron--Severi completion containing it; and
4. none of the eight seeds is obstructed at the rank-20 completion stage:
   all eight occur primitively, with negative-definite \(\tau\)-stable
   orthogonal complement, in the discriminant-19 K3 carrier \(H_{19}\).

The last assertion is an explicit integral gluing certificate, not merely an
application of an abstract embedding theorem. Thus \(H_{19}\) is a
**universal local carrier for the eight forced seeds**.

This is a lattice-completion theorem for the example track. It does not
produce a new compact hyperbolic polytope: effectivity, the orientation of a
common nef chamber, the absence of extra restricted K3 walls, and Vinberg
completeness remain separate gates.

## 1. Primitive embeddings of the seeds alone

The seed signatures and discriminant groups are as follows. The displayed
orthogonal-complement signature is the one forced inside the K3 lattice of
signature \((3,19)\).

| seed | \(\operatorname{sig}L_\Gamma\) | \(A_{L_\Gamma}\) | K3 complement signature | \(\ell(A_{L_\Gamma})\) |
|---|---:|---|---:|---:|
| \(A_4(3,5,3)\) | \((1,7)\) | \(\mathbf Z/11\) | \((2,12)\) | 1 |
| \(A_4(5,3,5)\) | \((1,7)\) | \(\mathbf Z/19\) | \((2,12)\) | 1 |
| \(D_4(5,3,3)\) | \((1,7)\) | \((\mathbf Z/2)^4\) | \((2,12)\) | 4 |
| \(C_4\), one 5 | \((1,7)\) | \(\mathbf Z/31\) | \((2,12)\) | 1 |
| \(C_4\), opposite 5s | \((1,7)\) | \(\mathbf Z/55\) | \((2,12)\) | 1 |
| \(A_5(5,3,3,3)\) | \((1,9)\) | \((\mathbf Z/2)^2\) | \((2,10)\) | 2 |
| \(A_5(5,3,3,5)\) | \((1,9)\) | \(\mathbf Z/2\oplus\mathbf Z/10\) | \((2,10)\) | 2 |
| \(D_5\), terminal 5 | \((1,9)\) | \((\mathbf Z/4)^2\) | \((2,10)\) | 2 |

For rank-eight seeds, the complement rank is 14 and
\(\ell(A_{L_\Gamma})\leq4\leq14-2\). For rank-ten seeds, the complement rank
is 12 and \(\ell(A_{L_\Gamma})\leq2\leq12-2\). Both positive and negative
indices increase strictly on passing from \(L_\Gamma\) to \(\Lambda_{K3}\).
Nikulin's stable-range primitive-embedding theorem therefore gives existence
and uniqueness of the primitive embedding orbit in every row.

This conclusion concerns the seed only. Its K3 orthogonal complement has
signature \((2,12)\) or \((2,10)\), so it gives a positive-dimensional family
of lattice-polarized K3 surfaces. It does **not** yet construct a rank-20
Néron--Severi lattice or a ten-dimensional golden eigenspace.

## 2. Exact even-overlattice and golden-saturation audit

Even overlattices of an even lattice \(L\) are classified by isotropic
subgroups \(H\leq(A_L,q_L)\); their index is \(|H|\). The golden action
extends to the overlattice precisely when \(H\) is \(\tau\)-invariant. The
implementation constructs the Smith discriminant form exactly, enumerates
every isotropic subgroup, builds the corresponding integral Gram matrix, and
tests \(\tau H=H\).

The nontrivial even overlattices are:

| seed | proper even overlattices | proper \(\tau\)-stable ones |
|---|---|---:|
| \(A_4(3,5,3)\) | none | 0 |
| \(A_4(5,3,5)\) | none | 0 |
| \(D_4(5,3,3)\) | three of index 2, each with \(A=(\mathbf Z/2)^2\) | 0 |
| \(C_4\), one 5 | none | 0 |
| \(C_4\), opposite 5s | none | 0 |
| \(A_5(5,3,3,3)\) | one even unimodular overlattice of index 2 | 0 |
| \(A_5(5,3,3,5)\) | one index-2 overlattice with \(A=\mathbf Z/5\) | 0 |
| \(D_5\), terminal 5 | one index-2 lattice with \(A=(\mathbf Z/2)^2\), and two even unimodular index-4 lattices | 0 |

For reproducibility, the action of \(\tau\) on the canonical Smith generators
is:

| seed | Smith moduli | images of the Smith basis |
|---|---|---|
| \(A_4(3,5,3)\) | \((11)\) | \((4)\) |
| \(A_4(5,3,5)\) | \((19)\) | \((5)\) |
| \(D_4(5,3,3)\) | \((2,2,2,2)\) | \((0,1,0,0),(1,1,0,0),(0,0,1,1),(0,0,1,0)\) |
| \(C_4\), one 5 | \((31)\) | \((19)\) |
| \(C_4\), opposite 5s | \((55)\) | \((8)\) |
| \(A_5(5,3,3,3)\) | \((2,2)\) | \((0,1),(1,1)\) |
| \(A_5(5,3,3,5)\) | \((2,10)\) | \((1,5),(1,8)\) |
| \(D_5\), terminal 5 | \((4,4)\) | \((2,1),(3,3)\) |

If \(L_\Gamma\subset N\) is any even \(\tau\)-stable completion, its primitive
closure

\[
\overline L_\Gamma=(L_\Gamma\otimes\mathbf Q)\cap N
\]

is an even \(\tau\)-stable self-overlattice. The enumeration proves
\(\overline L_\Gamma=L_\Gamma\). This is stronger than choosing one primitive
embedding: primitivity is forced in *every* compatible completion.

## 3. Explicit rank-20 completion in the \(H_{19}\) carrier

Let \(N_{19}\) be the trace lattice of

\[
H_{19}=
\begin{pmatrix}-2\phi&\phi\\ \phi&2\end{pmatrix}\oplus E_8
\]

over \(\mathbf Z[\phi]\). It has signature \((1,19)\), determinant \(-19\),
discriminant group \(\mathbf Z/19\), and a primitive K3 embedding with
rank-two transcendental complement \(T_{19}\).

The bounded exact Vinberg prefix through control-coordinate bound 5 contains
the following induced Lanner subdiagrams. Root indices refer to
run_h19_vinberg_search(5).roots.

| seed | root indices | \(\det C\) | \(A_C\) | \([N_{19}:L_\Gamma\perp C]\) |
|---|---|---:|---|---:|
| \(A_4(3,5,3)\) | \(0,9,14,13\) | 209 | \(\mathbf Z/209\) | 11 |
| \(A_4(5,3,5)\) | \(0,12,8,9\) | 361 | \((\mathbf Z/19)^2\) | 19 |
| \(D_4(5,3,3)\) | \(0,12,13,9\) | 304 | \((\mathbf Z/2)^3\oplus\mathbf Z/38\) | 16 |
| \(C_4\), one 5 | \(1,0,13,16\) | 589 | \(\mathbf Z/589\) | 31 |
| \(C_4\), opposite 5s | \(11,19,13,16\) | 1045 | \(\mathbf Z/1045\) | 55 |
| \(A_5(5,3,3,3)\) | \(1,0,13,9,2\) | 76 | \(\mathbf Z/2\oplus\mathbf Z/38\) | 4 |
| \(A_5(5,3,3,5)\) | \(2,4,3,16,17\) | 380 | \(\mathbf Z/2\oplus\mathbf Z/190\) | 20 |
| \(D_5\), terminal 5 | \(1,13,6,2,16\) | 304 | \(\mathbf Z/4\oplus\mathbf Z/76\) | 16 |

Here \(C=L_\Gamma^\perp\subset N_{19}\). In every row the computation
certifies:

- the selected field Gram matrix is exactly the named Lanner matrix;
- the rectangular Smith diagonal of the trace-coordinate inclusion consists
  entirely of ones, so the seed is primitive;
- \(C\) is negative definite and preserved by \(\tau\);
- \(\tau^2=\tau+I\) and \(\tau\) is self-adjoint on \(C\);
- the gluing subgroup in \(A_{L_\Gamma}\oplus A_C\) is isotropic,
  \(\tau\)-stable, and projects isomorphically onto all of \(A_{L_\Gamma}\);
  and
- the determinant identity

\[
|\det L_\Gamma|\,|\det C|
=19\,[N_{19}:L_\Gamma\perp C]^2
\]

holds exactly.

In fact the gluing index is \(|\det L_\Gamma|\) in every row: the full seed
discriminant is absorbed by the complement, leaving the ambient
\(\mathbf Z/19\) discriminant.

## 4. Why discriminants 11 and 19 are natural, but not forced

The determinant-11 path seed also has the following matching completion in
the discriminant-11 carrier:

\[
\det L=-11,\qquad \det C=121,\qquad
A_C=(\mathbf Z/11)^2,\qquad [N_{11}:L\perp C]=11.
\]

In canonical Smith coordinates, \(\tau\) acts by \(4\) on
\(A_L=\mathbf Z/11\), by the basis images \((8,0),(1,4)\) on \(A_C\), and
by \(8=1-4\) on the residual ambient \(\mathbf Z/11\).

For the determinant-19 path in \(N_{19}\), the corresponding data are

\[
\det L=-19,\qquad \det C=361,\qquad
A_C=(\mathbf Z/19)^2,\qquad [N_{19}:L\perp C]=19.
\]

Here \(\tau\) acts by \(5\) on \(A_L\), by \((15,0),(18,5)\) on the chosen
basis of \(A_C\), and by \(15=1-5\) on the residual ambient
\(\mathbf Z/19\). The gluing line in each complement is a \(\tau\)-eigenline
with the seed eigenvalue.

This makes the two prime-discriminant paths especially clean inverse-design
seeds. It does not make the carrier canonical: \(N_{19}\) contains all eight
seed types, so the local discriminant form alone does not determine the
rank-20 completion.

## 5. Geometric wall viability is still stronger

The lattice theorem does not promote the displayed root classes to a
complete list of irreducible effective curves. On a K3 surface every
\((-2)\)-class has one effective sign by Riemann--Roch, but it can be
reducible, and its restricted hyperplane can be redundant, touch a lower
face, or cut the proposed eigen-nef cell. A valid 14-wall construction must
therefore enumerate the relevant \((-2)\)-classes in the *completed* rank-20
lattice and certify the resulting restricted chamber.

The two affine \(C_4\) seeds illustrate the distinction sharply. If the four
\(\delta_i\) are irreducible effective curves, then

\[
F=\delta_1+\delta_2+\delta_3+\delta_4
\]

is a primitive nef isotropic fiber class. Exact two-Gram calculation gives

\[
F^2=0,\qquad F\cdot\tau F=(\tau F)^2=
\begin{cases}
2,&\text{one 5-edge},\\
4,&\text{two opposite 5-edges}.
\end{cases}
\]

Thus \(F\) is not a \(\tau\)-eigenvector and does not itself give an ideal ray
of the golden eigenslice. This branch survives the local compactness test,
but it carries a forced elliptic fibration that the remaining walls must
respect.

## 6. A rank-18 finite-vertex discriminant gate

There is a useful corollary for the eventual 14-wall search. Let \(V\) be the
saturation of the rank-18 trace lattice generated by the nine walls at a
finite vertex inside a rank-20 K3 Néron--Severi lattice \(N\), and put
\(K=V^\perp\), of rank two. For every prime \(p\),

\[
\ell_p(A_V)\leq \ell_p(A_N)+\operatorname{rk}K\leq4.
\]

Here is the precise length argument. Since \(V\) is primitive, the finite
index extension \(V\perp K\subset N\) is described by an isotropic graph
subgroup \(H\subset A_V\oplus A_K\), whose projections to both factors are
injective. Work with the \(p\)-primary parts and put

\[
a=\ell_p(A_V),\quad b=\ell_p(A_K),\quad
h=\ell_p(H),\quad n=\ell_p(A_N).
\]

The exact sequences

\[
0\longrightarrow H^\perp\longrightarrow A_V\oplus A_K
\longrightarrow H^*\longrightarrow0,
\qquad
0\longrightarrow H\longrightarrow H^\perp
\longrightarrow A_N\longrightarrow0
\]

give \(a+b\leq2h+n\). Injectivity into \(A_K\) gives
\(h\leq b\), while \(b\leq\operatorname{rk}K=2\). Therefore

\[
a\leq n+b\leq n+2.
\]

Finally the primitive rank-two transcendental complement of \(N\) in the
unimodular K3 lattice has discriminant form opposite to \(q_N\), so
\(n=\ell_p(A_N)\leq2\). This proves the bound.

This is a necessary prefilter for saturated finite-vertex trace types. It is
not yet an enumeration of their possible \(\tau\)-stable saturations.

## Reproducibility

The exact certificate is implemented in
coxeter_polytopes.geometric.golden_lanner_k3 and tested by
tests/test_golden_lanner_k3.py:

~~~bash
pytest -q tests/test_golden_lanner_k3.py
~~~

The implementation uses exact \(\mathbf Q\), \(\mathbf Z[\phi]\), Smith and
Hermite normal forms only. No floating-point signature, embedding, or
discriminant calculation enters the result.

## Source and scope

The primitive-embedding step uses V. V. Nikulin, *Integral symmetric bilinear
forms and some of their applications*, Math. USSR Izvestija 14 (1980),
[DOI 10.1070/IM1980v014n01ABEH001060](https://doi.org/10.1070/IM1980v014n01ABEH001060).
The even-overlattice step is the standard isotropic-subgroup correspondence
for discriminant forms from the same paper. The \(H_{19}\) and \(H_{11}\)
carrier data and their K3 embeddings are independently audited elsewhere in
this repository.

No claim is made here that the bounded \(H_{19}\) Vinberg prefix is complete,
that its chamber is finite-volume or compact, or that any of the eight local
seed embeddings extends to a compact 14-facet chamber in \(\mathbb H^9\).
