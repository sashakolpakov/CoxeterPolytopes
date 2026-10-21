# Cusp divisors and fibrations on the rational `P12`--`P16` surfaces

## Positive theorem

> **Theorem.** Every one of the 28 ideal vertices of the complete rational
> `P12`--`P16` nef chambers is semiample on each of the two explicit,
> non-isomorphic carriers `a=0` and `a=1`. Twenty-seven cusp classes already
> give base-point-free pencils in their primitive class. The remaining
> primitive genus-three class \(F=4H-E_1-\cdots-E_{16}\) is rigid, but has
> exact semiample index four: \(h^0(mF)=1\) for \(m=1,2,3\), while \(|4F|\)
> is a base-point-free pencil of genus-nine fibres.

The two carriers have respectively \(j(C_0)=0\) and
\(j(C_1)=6912/31\), so this is genuinely the same known hyperbolic Coxeter
polytope realized with the full fibration theorem on two non-isomorphic
algebro-geometric objects. It is not a claim that the abstract Vinberg
polytopes are new.

The later [symbolic cusp-family theorem](rational_i1n_symbolic_cusp_family.md)
upgrades this two-fibre statement to every parameter in the nonempty
carrier-open family.  Two coefficientwise `ZZ[a]` rank/section charts cover
each apparent exceptional value, so the additional cusp-specific
exceptional polynomial is `1`.

| chamber | all nef rays | cusps | semiample cusps | primitive pencils | index-four pencils |
|---|---:|---:|---:|---:|---:|
| \(P_{12}\) | 37 | 2 | 2 | 2 | 0 |
| \(P_{13}\) | 46 | 3 | 3 | 3 | 0 |
| \(P_{14}\) | 94 | 5 | 5 | 5 | 0 |
| \(P_{15}\) | 145 | 6 | 6 | 6 | 0 |
| \(P_{16}\) | 325 | 12 | 12 | 11 | 1 |
| **total** | 647 | **28** | **28** | **27** | **1** |

“Index four” is used in the Halphen-type sense of first moving multiple. It
does not assert that the surface is a classical anticanonical Halphen surface.

## Vinberg's finite-volume condition in the divisor algebra

Let \(C_i\) be the irreducible negative curves giving the complete Mori cone,
and use the stored Lorentzian form

\[
  \langle x,y\rangle_L=-x_0y_0+\sum_{j>0}x_jy_j.
\]

The nef inequalities are

\[
  F\cdot C_i=-\langle F,C_i\rangle_L\geq0.
\]

An independent exact GMP `cdd` conversion gives the positive-square/isotropic
splits

\[
  35+2,\quad43+3,\quad89+5,\quad139+6,\quad313+12,
\]

reproducing the finite and ideal vertices in the Vinberg volume certificates.
For each isotropic ray \(F\), take every wall with \(F\cdot C_i=0\) and split
the Gram graph into connected components \(A_j\). Exact nullspace arithmetic
proves

\[
  \operatorname{inertia}G(A_j)=(|A_j|-1,0,1),
  \qquad
  \sum_{i\in A_j}a_iC_i=F
\]

with primitive positive integral coefficients. If there are \(c\) components,
the complete zero-neighbourhood has \(n+c-1\) nodes and total parabolic rank
\(n-1\). Thus Vinberg's ideal-endpoint gate becomes the following exact
algebro-geometric test:

> A boundary isotropic nef class is a finite-volume ideal endpoint exactly
> when its full orthogonal negative-curve configuration is a full-rank union
> of affine components whose primitive null cycles span that ray.

The word “full” is essential: a sampled curve collection cannot certify a
cusp or finite volume.

## The first 22 primitive pencils

Sixteen cusps have at least two affine components. Two components give
distinct effective divisors \(D_1,D_2\) in the same primitive class \(F\).
On a rational surface numerical and linear equivalence agree in `Pic`; the
components have no common curve and

\[
  D_1\cdot D_2=F^2=0.
\]

They are therefore disjoint and span a base-point-free pencil. Five more
cusps are \(H-E_1\), the line pencil through the first blow-up centre. On
\(X_{13}\), the class

\[
  4H-2E_1-E_2-\cdots-E_{13}
\]

has the absolutely irreducible \(Q_a\) and its distinct affine wall cycle as
two disjoint members. These arguments apply to both `a=0` and `a=1`.

## The exact plane-jet mobility criterion

In the chart \(Y=1\), put

\[
 C_a=z-x^3-a xz^2-z^3,
 \qquad Q_a=z^4+(x+z)C_a,
 \qquad S_a=C_a^2+z^2Q_a.
\]

Fifteen exact strict-transform calculations recover the selected branch
coefficients, rather than assuming a numerical Taylor expansion:

\[
\begin{aligned}
 \phi_{S,0}(x)&=x^3+x^9-x^{11}+x^{13}+x^{15}+O(x^{16}),\\
 \phi_{S,1}(x)&=x^3+x^7+x^9+x^{11}+6x^{13}+O(x^{16}).
\end{aligned}
\tag{4.1}
\]

After the fifteenth blow-up the exceptional restrictions are \(Q_a=z\) and
\(S_a=z+1\) for both parameters.

For

\[
 D=dH-\sum_{i=1}^n m_iE_i,
 \qquad A_{k-1}=m_1+\cdots+m_{k-1},
\]

truncate \(\phi_{S,a}\) before degree \(k\). A degree-\(d\) homogeneous form
\(G\) is a section of \(D\) exactly when, at every centre,

\[
 x^{-A_{k-1}}
 G\!\left(x,\phi_{k-1}(x)+x^{k-1}u\right)
\tag{4.2}
\]

has neither a negative power of \(x\) nor a term of total degree below
\(m_k\). The forbidden coefficients form an explicit integral matrix
\(J_D\), and

\[
 H^0(X,D)=\ker J_D.
\tag{4.3}
\]

This complete-ideal condition correctly permits exceptional components.
Imposing plane multiplicity \(m_i\) independently at every infinitely-near
point would be a stronger and incorrect test.

## Five new primitive pencils

The six formerly open primitive systems have the following exact dimensions
on both carriers:

| surface and cusp | plane class | monomials | forbidden jets | rank | \(h^0\) |
|---|---|---:|---:|---:|---:|
| \(X_{14}\), \(g=5\) | \(11H-3E_1-\cdots-3E_{13}-2E_{14}\) | 78 | 627 | 76 | 2 |
| \(X_{15}\), \(g=5\) | same class, coefficient of \(E_{15}\) zero | 78 | 708 | 76 | 2 |
| \(X_{16}\), \(g=5\) | same class, last two coefficients zero | 78 | 786 | 76 | 2 |
| \(X_{16}\), \(g=3\) | \(4H-E_1-\cdots-E_{16}\) | 15 | 136 | 14 | 1 |
| \(X_{16}\), \(g=10\) | \(21H-7E_1-\cdots-7E_5-6E_6-4E_7-\cdots-4E_{16}\) | 253 | 2,655 | 251 | 2 |
| \(X_{16}\), \(g=11\) | \(23H-7E_1-\cdots-7E_8-5E_9-4E_{10}-\cdots-4E_{16}\) | 300 | 3,139 | 298 | 2 |

For the five moving rows, one kernel vector is the known affine wall cycle,
whose relevant plane products are

\[
 CQ^2,\qquad LQ^2S^2,\qquad CQ^2S^2.
\]

The certificates store a complete sparse homogeneous polynomial for a second
kernel vector. The `a=0` witnesses have respectively 31, 91 and 114 terms;
the `a=1` witnesses have 32, 91 and 113 terms. Every `a=1` second polynomial
is different from its `a=0` mate, even though the nullities agree.

For every row and both parameters, exact replay proves that the second form
has precisely the assigned multiplicity sequence, hence no exceptional
component, and is coprime to the corresponding wall product. The two members
therefore have no common curve. Their intersection is \(F^2=0\), so they are
disjoint and their complete pencil is base-point-free.

## The last cusp: exact semiample index four

For

\[
 F=4H-E_1-\cdots-E_{16},
\]

the unique primitive section is the affine wall cycle with plane factor
\(Q_a\). Applying (4.2) to its first four multiples gives, for both `a=0`
and `a=1`,

| class | monomials | forbidden jets | rank | \(h^0\) | conclusion |
|---|---:|---:|---:|---:|---|
| \(F\) | 15 | 136 | 14 | 1 | rigid |
| \(2F\) | 45 | 408 | 44 | 1 | rigid |
| \(3F\) | 91 | 816 | 90 | 1 | rigid |
| \(4F\) | 153 | 1,360 | 151 | 2 | pencil |

The unique sections for \(m=1,2,3\) are the wall members \(Q_a^m\), so no
smaller multiple moves. For \(m=4\), the certificate supplies a second
degree-16 form with actual multiplicity sequence exactly \((4^{16})\): it
has no exceptional component and is coprime to \(Q_a^4\). It has 64 terms
for `a=0` and 63 terms for `a=1`, with respective sparse SHA-256 values

```text
a=0  1393b70fc23154a593c95781a8990caab364716c3cbcc308296f91883657ec12
a=1  f8b43b63e6e0c560303c398d983884f737623cc732b61ead6c2d604b18866b0b
```

The two \(4F\)-members have no common component and intersection zero, hence
are disjoint and make \(|4F|\) base-point-free. For connectedness, take the
Stein factorization through a curve \(B\). Since \(X\) is rational,
\(g(B)=0\), so \(B\simeq\mathbb P^1\). If the final map
\(B\to\mathbb P^1\) had degree \(e\), then connectedness of the first map
would give
\(h^0(X,4F)=h^0(\mathbb P^1,\mathcal O(e))=e+1\). The exact value
\(h^0(X,4F)=2\) forces \(e=1\); hence the pencil's general fibre is connected.
Adjunction gives

\[
 p_a(4F)=1+\frac{K_X\cdot4F}{2}=9.
\]

Thus the formerly negative primitive result combines with the positive
fourth-multiple construction to give the sharp criterion: **the cusp is
semiample, and its exact semiample index is four.**

## Primal-dual geometric criterion

The result separates two questions that should not be conflated:

1. Vinberg's full affine-neighbourhood test certifies that an isotropic ray is
   a finite-volume ideal endpoint and reconstructs its effective wall cycle.
2. The complete jet-kernel test computes \(H^0(mF)\). A coprime second kernel
   vector with exact cluster multiplicities certifies a fibration; successive
   one-dimensional kernels give a sharp lower bound on the semiample index.

Here the two boundaries meet exactly: 27 cusps move primitively, while the
only primitive obstruction converts into a new positive pencil at the exact
fourth multiple.

## Replay

The compositional exact artifacts are:

- [`certify_rational_i1n_cusp_fibrations.py`](../scripts/certify_rational_i1n_cusp_fibrations.py)
- [`rational_i1n_a1_cusp_jet_workers.json`](../data/enumerations/rational_i1n_a1_cusp_jet_workers.json)
- [`rational_i1n_cusp_fibrations.json`](../data/enumerations/rational_i1n_cusp_fibrations.json)
- [`test_rational_i1n_cusp_fibrations.py`](../tests/test_rational_i1n_cusp_fibrations.py)

The heavy `a=1` nullspaces are regenerated serially in fresh processes, so
GMP/FLINT heaps never overlap. The aggregate certificate pins the worker-file
SHA, its payload SHA, and the current certifier/family source hashes.

```bash
PYTHONPATH=.:src caffeinate -i python3 -m scripts.certify_rational_i1n_cusp_fibrations --rebuild-a1-workers
PYTHONPATH=.:src python3 -m scripts.certify_rational_i1n_cusp_fibrations
PYTHONPATH=.:src python3 -m scripts.certify_rational_i1n_cusp_fibrations --check-a1-workers
PYTHONPATH=.:src python3 -m scripts.certify_rational_i1n_cusp_fibrations --check-archive
PYTHONPATH=.:src pytest -q tests/test_rational_i1n_cusp_fibrations.py
```

Frozen payload SHA-256 values:

```text
a=1 worker certificate  3fbe4f7225af5964dd626386d0a6b497ff1568c7fdead7c8be9f8f3dc61975ed
aggregate v3 certificate 98db5f42c3720872fb181f9eb8763d50177510d356745d1b10b2299400d3196e
```
