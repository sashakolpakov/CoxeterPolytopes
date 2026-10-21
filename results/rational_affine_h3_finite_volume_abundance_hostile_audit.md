# Hostile audit of finite-volume affine-H3 abundance

**Verdict (2026-08-10): PASS AFTER REPAIR, with semiampleness upgraded to
base-point-freeness.** The finite-volume cusp argument really does produce
a base-point-free multiple of the affine H3 null cycle. It uses the
full-nef finite-volume hypothesis in exactly one place and uses no
reflection-maximality or classification theorem. The geometry of the
nonreduced star then forces the first multiple to move.

The first draft nevertheless contained two substantive downstream errors.

1. Its rank-six big-anticanonical endpoint was empty. Nefness and
   \(K_S.F=4\) already exclude every rank \(\rho\le10\).
2. Semiampleness alone does not identify the fibre index, so the original
   primitive-lattice-vector argument was invalid. On this actual surface,
   however, the multidegree-zero Picard group of the nonreduced affine star
   is torsion-free. This forces \(\mathcal O_F(F)\) to be trivial and makes
   \(|F|\) itself base-point-free.

The owner memo and certifier were patched. The valid theorem is now:

> For a literal finite-volume full-nef rational carrier satisfying the
> affine H3 hypotheses, \(F=4E+B_1+\cdots+B_4\) is base-point-free.
> Necessarily \(\rho(S)\ge11\). Its connected general fibre has class
> \(F\), genus three, and relatively minimal rank at most \(26\). The
> vertical exceptional tail is not bounded here.

## 1. Local lattice and wall audit

In the physical ordered basis \(E,B_1,\ldots,B_4\), the star matrix is

\[
 Q_0=
 \begin{pmatrix}
 -1&1&1&1&1\\
 1&-4&0&0&0\\
 1&0&-4&0&0\\
 1&0&0&-4&0\\
 1&0&0&0&-4
 \end{pmatrix}.
\]

Exact elimination gives

\[
 \ker Q_0=\mathbf Q(4,1,1,1,1),\qquad
 \operatorname {rank}Q_0=4.                                \tag{1.1}
\]

The five classes themselves are independent in
\(\operatorname {NS}(S)\): pairing a relation with each \(B_i\) makes it a
multiple of \(F\), and an ample class sees that the nonzero effective cycle
\(F\) is not numerically zero. Thus their span is five-dimensional with
radical \(\mathbf RF\). The crude lattice conclusion \(\rho\ge6\) is
correct; Section 5 below gives the stronger geometric conclusion
\(\rho\ge11\).

In the doubled host, \(e=E\) and \(h_i=B_i/2\) all have square \(-2\), with
\(e.h_i=1\), and

\[
                         F/2=2e+h_1+h_2+h_3+h_4.            \tag{1.2}
\]

This is the positive null relation of the affine \(D_4\) star. No
root-system classification is being smuggled into the proof:
(1.1)--(1.2) directly show that its quotient normals span and bound a
four-dimensional Euclidean factor.

If another physical negative curve \(C\) has \(F.C=0\), then

\[
                   0=4E.C+\sum_{i=1}^4B_i.C.               \tag{1.3}
\]

All five summands are nonnegative because the curves are distinct and
irreducible on a smooth surface. Hence each is zero. This verifies the
termwise orthogonality used by the owner memo. The root scaling
\(C,C,C/2\) for squares \(-1,-2,-4\) does not change it.

## 2. Finite volume really makes the cusp link compact

The hypothesis here is not bare rational polyhedrality. The nef chamber is
a finite-sided, full-dimensional, finite-volume hyperbolic polyhedron, and
\([F]\) is an isotropic point in its closure.

A finite-volume polyhedron cannot have a positive-dimensional ideal face.
After moving \([F]\) to infinity in the upper-half-space model, a
sufficiently deep horosphere meets only the finitely many walls through
\([F]\). Their inequalities cut out a full-dimensional Euclidean
polyhedron. A nonzero recession direction would give it infinite Euclidean
volume and hence give the associated cusp sector infinite hyperbolic
volume. Therefore the link is bounded and \([F]\) is an ordinary
finite-volume ideal vertex.

This argument would fail for an infinite wall arrangement without a
local-finiteness input, or for a merely polyhedral infinite-volume cone.
Neither weakened setting is claimed in the theorem.

## 3. Orthogonal product and positive rational dependence

Put \(V_F=F^\perp/\mathbf RF=U_0\perp U_1\), where \(U_0\) is the
four-dimensional quotient span of the star. A cusp-wall affine inequality
depends on its quotient normal and a constant. Star normals have only an
\(U_0\)-component, while (1.3) puts every other cusp normal in \(U_1\).
Consequently the horospherical link is literally

\[
                               P_0\times P_1.               \tag{3.1}
\]

The affine constants do not mix the two variables, so orthogonality of the
linear normals is enough for this product statement. Since the whole link
and \(P_0\) are bounded, \(P_1\) is bounded.

For a nonempty full-dimensional rational Euclidean polyhedron, boundedness
is equivalent to its inward normals positively spanning the ambient vector
space. Thus the nonstar physical cusp roots admit a rational dependence

\[
                     \sum_j a_j\bar r_j=0,
                     \qquad a_j\in\mathbf Q_{>0}.           \tag{3.2}
\]

One exact construction is to express the negative of each normal as a
nonnegative rational combination of all normals and then sum. Rationality
follows from the rational defining system. Lifting (3.2) from \(V_F\) gives

\[
                            \sum_j a_jr_j=\lambda F.        \tag{3.3}
\]

Intersection with an ample divisor proves \(\lambda>0\): both sides are
nonzero positive combinations of effective curve classes. Clearing the
rational coefficients and the possible half-root scaling for a
square-minus-four curve gives an effective integral divisor

\[
                             G\equiv nF,\qquad n>0.          \tag{3.4}
\]

No saturation, reflection maximality, or hidden sign choice enters this
step.

## 4. Numerical-to-linear and base-point-free audit

Every component of \(G\) is disjoint from every component of \(F\) by
(1.3). A smooth rational surface has \(q(S)=0\), torsion-free Picard group,
and hence no difference between integral numerical and linear equivalence.
Therefore \(G\sim nF\).

The sections with zero divisors \(G\) and \(nF\) have no common zero. They
globally generate \(\mathcal O_S(nF)\), so the complete system \(|nF|\) is
base-point-free. This proves semiampleness. The step is stronger than
merely finding two numerically equivalent cycles and has been stated with
the needed Picard hypotheses.

## 5. The low-rank poison test

The first draft attempted to treat \(\rho=6\) using
Testa--Várilly-Alvarado--Velasco. The primary source locator was accurate,
but its hypotheses expose a contradiction instead of a new endpoint.

For every smooth rational surface,

\[
                           K_S^2=10-\rho(S).                \tag{5.1}
\]

If \(\rho\le10\), rational Riemann--Roch and vanishing of rational
plurigenera give

\[
 \chi(\mathcal O_S(-K_S))=1+K_S^2\ge1,\qquad
 h^2(-K_S)=h^0(2K_S)=0.                                   \tag{5.2}
\]

Hence \(-K_S\) is effective. This is impossible because \(F\) is nef but

\[
                              (-K_S).F=-4.                 \tag{5.3}
\]

Thus \(\rho\ge11\), and the positive-relation proof of Sections 2--4 covers
every nonempty case.

For reproducibility, the audited primary source is
[arXiv:0901.1094v2](https://arxiv.org/abs/0901.1094v2). Its source archive
has SHA-256
146fdf403f85ba8383edb46777a9df5fd6a61a1f4ed4718d8661b90f757e27d3,
and BigRationalSurfacesMarch25.tex has SHA-256
6a23a31d295bb0d2bfdb644aba1a2d778edc3aaa78e17058f5dfde774b1e6e94.
Lemma 2.4 says that big \(-K_X\) forces \((-K_X).N>0\) for every nonzero
nef \(N\); Lemma 2.5 then proves such a nef non-big \(N\)
base-point-free. Putting \(N=F\) violates Lemma 2.4 before Lemma 2.5 can
supply a nonvacuous case. Theorem 2.9 on Cox-ring finite generation is
correct but unnecessary here.

## 6. The fibration-index attack and its AG repair

The first hostile replay correctly observed that semiampleness by itself
only gives a minimal moving multiple \(qF\). A purely lattice-theoretic
argument cannot force \(q=1\). The missing input is the Picard group of the
actual nonreduced fibre.

Let

\[
                         C=F_{\mathrm {red}}
                           =E+B_1+B_2+B_3+B_4.              \tag{6.1}
\]

This is an SNC tree of five rational curves. Consequently a
multidegree-zero line bundle on \(C\) is trivial:

\[
                           \operatorname {Pic}^0(C)=0.      \tag{6.2}
\]

Let \(\mathcal I\) be the nilradical of \(C\subset F\). The exact units
sequence is

\[
 1\longrightarrow1+\mathcal I\longrightarrow\mathcal O_F^*
 \longrightarrow\mathcal O_C^*\longrightarrow1.            \tag{6.3}
\]

In characteristic zero the finite logarithm is an isomorphism of sheaves
of abelian groups \(1+\mathcal I\simeq\mathcal I\). Constants on \(C\)
lift to \(F\), so cohomology identifies

\[
 \ker\!\left(\operatorname {Pic}(F)\to\operatorname {Pic}(C)\right)
                         \simeq H^1(F,\mathcal I).           \tag{6.4}
\]

The right side is an additive complex vector group, hence has no
nontrivial torsion.

For \(N=\mathcal O_F(F)\), the degrees on the five reduced components are
the intersections \(F.E,F.B_1,\ldots,F.B_4\), all zero. Thus \(N|_C\) is
trivial and \(N\) lies in (6.4). The disjoint divisor \(G\sim nF\) from
Section 4 is a unit along the entire nonreduced scheme \(F\), so

\[
                              N^n\simeq\mathcal O_F.        \tag{6.5}
\]

Torsion-freeness forces \(N\simeq\mathcal O_F\). In the exact sequence

\[
 0\longrightarrow\mathcal O_S\longrightarrow\mathcal O_S(F)
 \longrightarrow\mathcal O_F(F)\longrightarrow0,           \tag{6.6}
\]

the obstruction to lifting \(1\in H^0(F,\mathcal O_F)\) lies in
\(H^1(S,\mathcal O_S)=0\). Its lift is nowhere zero along \(F\); together
with the canonical section vanishing on \(F\), it generates
\(\mathcal O_S(F)\). Therefore \(|F|\) is base-point-free.

### 6.1 Lattice poison: why the AG step is essential

Fix \(M=U\oplus D_4(-1)\), a lattice whose discriminant group has exponent
two, and a primitive isotropic vector \(e\in U\). If \(\theta\) is the
highest \(D_4\) root, then for every \(r\ge1\)

\[
                              \alpha_0=re-\theta            \tag{6.7}
\]

together with the four finite simple roots has exactly the same affine
\(\widetilde D_4\) Gram matrix, but the positive null relation is
\(\alpha_0+\theta=re\). Thus neither the local Gram matrix nor even this
fixed small discriminant bounds null-vector divisibility. The surface
argument succeeds because it sees the normal bundle
\(\mathcal O_F(F)\), which the Néron--Severi lattice does not.

### 6.2 Genus and core rank

The pencil from \(|F|\) has \(F\) as a fibre. Its Stein base is
\(\mathbf P^1\) because \(q(S)=0\). The support of \(F\) is connected and
the coefficient vector \((4,1,1,1,1)\) has gcd one. Hence the Stein degree
is one; otherwise all fibre coefficients would be divisible by that
degree. A connected general fibre \(H\) is numerically \(F\), and

\[
                           2g(H)-2=K_S.F=4,\qquad g(H)=3.   \tag{6.8}
\]

On the relatively minimal model, \(K_f=K_{S_{\min}}+2H\) is nef, so

\[
 0\le K_f^2=10-\rho(S_{\min})+8(g(H)-1)
             =26-\rho(S_{\min}).                           \tag{6.9}
\]

Therefore \(\rho(S_{\min})\le26\). The former genus-two/rank-18 branch came
from the invalid primitive-vector argument and is deleted.

## 7. Final scope

The hostile result is therefore deliberately split.

- **Proved:** finite-volume full-nef affine H3 makes \(F\) itself
  base-point-free, with no reflection maximality or classification input.
- **Also proved:** the numerical configuration is empty in ranks at most
  ten; the resulting fibration has genus three and relatively minimal rank
  at most \(26\).
- **Not proved:** a bound on vertical blowups or an absolute
  Picard/hyperbolic-dimension bound.

The independent certifier replays both Gram matrices, the termwise
orthogonality poison test, half-root denominator clearing, the five empty
low-rank rows, the rational-tree Picard gate, the lattice divisibility
poison family, and the genus-three core calculation.
