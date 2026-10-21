# The minimal H1 pair: a sharp rank-14 log core and the unresolved exceptional tail

**Status (2026-08-10): proved core theorem; no original-rank or
reflection-maximality theorem is claimed.** This note treats the smallest
standard-host H1 incidence

\[
 B_1^2=B_2^2=-4,\qquad B_1.B_2=2,
\]

with two transverse intersection points. It uses only adjunction, the
canonical exact sequence, contraction of square-minus-one curves, and exact
intersection arithmetic. It proves that the relative log-minimal core has
Picard rank at most \(14\), and gives an explicit complete-Mori example at
equality. It does **not** prove that the Picard rank before relative
minimalization is bounded. The missing step is precisely a global theorem
forcing reflection maximality to bound the exceptional tail.

No Nikulin--Vinberg reflective-lattice classification, ABSW theorem, Prasad
formula, or arithmetic volume is used.

## 1. The sharp log-minimal core theorem

Let \(S\) be a smooth projective rational surface over an algebraically
closed field and let \(B_1,B_2\) be smooth rational curves satisfying

\[
 B_1^2=B_2^2=-4,\qquad B_1.B_2=2.                 \tag{1.1}
\]

For the clean H1 statement assume that the intersection consists of two
distinct transverse points. Put

\[
                         Z=B_1+B_2,\qquad L=K_S+Z.             \tag{1.2}
\]

### Theorem 1.1

There is a sequence of blowdowns, all disjoint from \(Z\),

\[
                  S=S_m\longrightarrow S_{m-1}\longrightarrow
                  \cdots\longrightarrow S_0,                 \tag{1.3}
\]

such that every contracted curve is a smooth rational square-minus-one
curve and

\[
                         K_{S_0}+Z_0\quad\hbox{is nef}.         \tag{1.4}
\]

The pair is unchanged near \(Z\), and

\[
                         \rho(S_0)\leq14.                      \tag{1.5}
\]

Equivalently, the hyperbolic dimension of the log-minimal core is at most
\(13\).

### Proof

The curve \(Z\) is a reduced Gorenstein genus-one cycle. Indeed,

\[
 Z^2=-4,\qquad K_S.Z=4,\qquad
 p_a(Z)=1+\frac{Z^2+K_S.Z}{2}=1.               \tag{1.6}
\]

Moreover

\[
       \deg\bigl(\omega_Z|_{B_i}\bigr)=(K_S+Z).B_i=0.          \tag{1.7}
\]

The usual residue differential on the two-component nodal cycle is nowhere
zero, so

\[
                              \omega_Z\simeq\mathcal O_Z.      \tag{1.8}
\]

Use the canonical exact sequence

\[
 0\longrightarrow\mathcal O_S(K_S)
  \longrightarrow\mathcal O_S(K_S+Z)
  \longrightarrow\omega_Z\longrightarrow0.                   \tag{1.9}
\]

Rationality gives
\(H^1(S,K_S)=H^1(S,\mathcal O_S)^\vee=0\). Hence the unit section
of \(\omega_Z\simeq\mathcal O_Z\) lifts to a section of
\(\mathcal O_S(L)\). Its zero divisor \(A\in|L|\) is effective and is
disjoint from \(Z\); \(A=0\) is allowed.

If \(A\) is not nef, an irreducible component \(E\) has \(A.E<0\).
Because \(A\) is effective, \(E^2<0\); because its support is disjoint from
\(Z\),

\[
                              K_S.E=A.E<0.                     \tag{1.10}
\]

Adjunction now forces \(E\) to be a square-minus-one curve. Explicitly, if
\(E^2=-q\), then

\[
 K_S.E=2p_a(E)-2+q<0,
\]

so \(p_a(E)=0\) and \(q=1\). Thus \(E\cong\mathbf P^1\), and it can be
blown down. The contraction is an isomorphism near \(Z\), and pushing the
lifted section forward gives the same construction on the new smooth
rational surface. Repeating decreases the Picard rank and therefore stops,
exactly when \(K+Z\) is nef.

On the terminal surface, Noether's identity and (1.1) give

\[
\begin{aligned}
 (K_{S_0}+Z_0)^2
 &=K_{S_0}^2+2K_{S_0}.Z_0+Z_0^2\\
 &=(10-\rho(S_0))+8-4\\
 &=14-\rho(S_0).
\end{aligned}                                                  \tag{1.11}
\]

A nef divisor on a smooth surface has nonnegative square, proving
\(\rho(S_0)\leq14\). \(\square\)

The same calculation explains both the strength and the exact scope. For
an intersection number \(I=B_1.B_2\), one formally gets

\[
                       (K+B_1+B_2)^2=10-\rho+2I.               \tag{1.12}
\]

But when \(I>2\), the dualizing sheaf has positive degree \(I-2\) on each
component; a lifted canonical section is no longer disjoint from the pair.
Thus (1.12) is not an absolute core bound for the \(I\geq4\) H1 branches.

## 2. Contraction and discrepancy audit

The pair in (1.1) has negative-definite intersection matrix

\[
 Q=\begin{pmatrix}-4&2\\2&-4\end{pmatrix}.                    \tag{2.1}
\]

Over \(\mathbf C\), Grauert's theorem therefore contracts it to one normal
surface point in the complex-analytic category.  For the algebraic or
projective and \(\mathbf Q\)-Gorenstein discrepancy discussion below, assume
that a contraction \(f:S\to Y\) of precisely this pair has been supplied. If

\[
                  K_S=f^*K_Y+a_1B_1+a_2B_2,
\]

then adjunction gives the exact system

\[
 Q\binom{a_1}{a_2}=\binom 22,\qquad
                  (a_1,a_2)=(-1,-1).                          \tag{2.2}
\]

Over \(\mathbf C\), with two transverse nodes, this is the two-component
cusp configuration: it is log canonical but not klt and is not a quotient
singularity. Consequently the quotient-orbifold contribution
\(1/|G_p|\) and the klt orbifold-\(c_2\) theorem used in the disjoint H1-free
lane are unavailable here.

The intersection number alone does not certify log canonicity when the two
contacts coalesce. If the curves have one tangency of order two, blow up
the tangency. The two strict transforms and the first exceptional curve
have discrepancy \(-1\) and meet at one triple point. Blowing up that triple
point creates discrepancy

\[
                              1-1-1-1=-2,                      \tag{2.3}
\]

so the contraction is not log canonical. The transverse hypothesis in
Theorem 1.1 is therefore substantive.

There is also a classification-free check for a larger negative-definite
boundary packet. Standard-host integrality makes every positive boundary
intersection even. A physical intersection at least \(4\) already makes
the corresponding \(2\times2\) principal minor non-positive. Hence every
edge in a negative-definite connected packet has physical intersection
exactly \(2\). If \(G\) is its ordinary adjacency matrix, then

\[
                            Q=-2(2I-G).                         \tag{2.4}
\]

Writing \(x=(2I-G)^{-1}{\bf1}\), the discrepancies are \(-x\). At a
vertex where \(x\) is minimal, the equation

\[
                     2x_v-\sum_{w\sim v}x_w=1                 \tag{2.5}
\]

forces that vertex to be a leaf and \(x_v\geq1\). Equality propagates and
is possible only for the two-vertex packet. Thus the pair (2.1) is the
unique connected multi-boundary packet with discrepancies exactly \(-1\);
every connected packet with at least three vertices is worse than log
canonical. No ADE list is needed for this conclusion.

## 3. Equality is real: a complete-Mori rank-14 carrier

The bound (1.5) is attained, so H1 is not an empty formal obstruction.
Start with the Hirzebruch surface \(\mathbb F_2\), with negative section
\(s\) and fiber \(f\):

\[
                    s^2=-2,\qquad f^2=0,\qquad s.f=1.          \tag{3.1}
\]

Choose two smooth transverse curves

\[
                    \bar B_i\in|s+2f|.                        \tag{3.2}
\]

They have self-intersection \(2\), mutual intersection \(2\), and
\(\bar B_1+\bar B_2=-K_{\mathbb F_2}\). Choose private points
\(q_i\in\bar B_i\), away from the mutual intersections, lying on one common
fiber. Perform six successive interior boundary blowups over each \(q_i\),
each time at the point of the strict boundary on the newest exceptional
curve. Then

\[
 B_i^2=2-6=-4,\qquad B_1.B_2=2,\qquad
 K_S+B_1+B_2=0,\qquad \rho(S)=2+12=14.                        \tag{3.3}
\]

The complete cone calculation is equation-level, not a census. In the
\(n=2\) construction of Jennifer Li, the irreducible exceptional curves and
the two distinguished fibers form a basis whose dual basis is given by
effective combinations of the same curves and the boundary. Her elementary
dual-basis lemma then gives

\[
 \overline{NE}(S)=
 \left\langle
 B_1,B_2,\ E_{i,j}\ (i=1,2,\ 1\leq j\leq6),\ F_1,F_2
 \right\rangle_{\mathbf R_{\geq0}}.                           \tag{3.4}
\]

Here the negative alphabet is exactly

\[
              2\cdot(-4)+12\cdot(-2)+2\cdot(-1).             \tag{3.5}
\]

Indeed, in each exceptional arm the first five strict exceptional curves
have square \(-2\) and the last has square \(-1\); the common fiber and the
negative section have square \(-2\). Thus (3.4) is a complete finite
rational Mori alphabet, not merely a chosen curve configuration.

For the standard host

\[
 M=\langle\operatorname{NS}(S),B_1/2,B_2/2\rangle,\qquad q_M=2q_S,          \tag{3.6}
\]

the sixteen facet normals have split

\[
                     4\cdot(-2)+12\cdot(-4)                  \tag{3.7}
\]

in \(q_M\). Their intersections give only the crystallographic orders
\(2,3,4\) (and disjoint pairs), so the full nef cone is a literal finite
Coxeter chamber. It also has finite hyperbolic volume: the cone is
rational polyhedral, every nef class has nonnegative square, and in
signature \((1,13)\) a totally isotropic face is at most a ray, so its only
vertices at infinity are ordinary ideal vertices. The host has determinant
\(2^{10}\): starting from the
unimodular rank-14 Neron--Severi lattice, doubling the form gives determinant
\(2^{14}\), and adjoining the two independent half-boundaries has index
\(4\), hence divides the determinant by \(4^2\).

This construction is extracted from the primary \(n=2\) formulas in
[Jennifer Li, *A cone conjecture for log Calabi--Yau surfaces*](https://arxiv.org/abs/2207.12483),
Section 5. The downloaded v1 source archive and its main TeX member have
SHA-256 hashes

    source archive  20991983ad5b4c4ab8fbd3340f3b438e91068f8e7e81cc851aaad9e2bc85d107
    main TeX        da8ba69b8f96359f75fa95bfd5da44aa86f28dd98337f4f24030ed86cb0e4892

No claim in this note that the chamber is reflection-maximal is imported
from a reflective-lattice table.

## 4. Why the original rank does not follow

Let \(\pi:S\to S_0\) be the composite of the \(m\) blowdowns in (1.3),
and use total-transform exceptional classes \(e_1,\ldots,e_m\). Since all
centers are disjoint from the boundary,

\[
 \operatorname{NS}(S)=\pi^*\operatorname{NS}(S_0)
       \mathbin\perp\bigoplus_{j=1}^m\mathbf Z e_j,\qquad e_j^2=-1,         \tag{4.1}
\]

and the half-boundaries are pulled back from \(S_0\). Therefore the standard
host splits exactly as

\[
                    M(S)=M(S_0)\mathbin\perp\langle-2\rangle^m.             \tag{4.2}
\]

This is the numerical exceptional tail left invisible by the rank-14 core
bound.

Nor does a linear tail automatically create a missing *pure exceptional*
reflection. For successive free blowups on the newest exceptional curve,
the physical exceptional curves are

\[
 e_1-e_2,\ e_2-e_3,\ \ldots,\ e_{m-1}-e_m,\ e_m.              \tag{4.3}
\]

In \(q_M\), the first \(m-1\) have norm \(-4\) and their reflections swap
adjacent \(e_j\)'s; the last has norm \(-2\) and its reflection changes the
sign of \(e_m\). These reflections generate every signed permutation of the
total exceptional basis, in particular the reflection in each \(e_j\).
Thus the local exceptional root system is already saturated. Any
maximality obstruction must come from a **mixed** root involving the core,
or from a global failure of finite volume; it cannot be obtained by merely
noticing the tail.

Consequently the proved statement is

\[
 \boxed{\rho(S_0)\leq14,\qquad
        \rho(S)=\rho(S_0)+m\leq14+m,}                          \tag{4.4}
\]

with no theorem here bounding \(m\).

There is one useful conditional closure. If an integral ample marking \(H\)
comes with an absolute bound

\[
             T=\sum_{C\ {\rm physical\ facet}} H.C,           \tag{4.5}
\]

then every exceptional component is a distinct physical facet and has
positive integral \(H\)-degree. Hence \(m\leq T\), and

\[
                  \rho(S)\leq14+T,\qquad
                  \dim\mathbf P(\operatorname{Nef}S)\leq13+T.              \tag{4.6}
\]

This is exactly where an absolute **total marked-facet degree** bound would
finish H1. A bound on the degree of each individual facet, without a bound
on their total or number, does not give (4.6).

There is a second conditional consequence for the arithmetic trace packet.
Suppose a full trace lattice \(T\), of rank \(d(n+1)\), is carried
isometrically by \(\operatorname{NS}(S)\), where
\(d=[k:\mathbf Q]\), and suppose it actually descends
through the particular H1 tail contraction \(\pi:S\to S_0\). Equivalently,
its image is orthogonal to the exceptional lattice in (4.1). Then

\[
                  T\hookrightarrow\operatorname{NS}(S_0)
                  \quad\Longrightarrow\quad
                  d(n+1)\leq\rho(S_0)\leq14.                  \tag{4.7}
\]

For \(k=\mathbf Q\), this gives \(n\leq13\). This is useful when descent has
an independent geometric certificate, but it cannot be silently inferred
from reflection maximality. In a literal full-nef carrier the exceptional
components in (4.3) are themselves facet roots. Hence they lie in, or span
part of, the visible full trace root block and have nonzero pairing with the
exceptional lattice (already their self-pairings are negative). Such a full
trace block is generally **not** orthogonal to the tail, so the
finite-orbit/full-trace descent premise may fail or be impossible. The
unconditional target remains a bound on the tail length \(m\), not an
assumed descent around it.

## 5. Honest conclusion

The minimal transverse H1 incidence is possible and its log-minimal core has
the sharp absolute rank bound \(14\). The contraction is lc non-klt, so the
klt orbifold route genuinely stops. Intersection theory and discrepancy
control do not bound the number of disjoint-from-boundary blowups, and the
exceptional tail contains no automatic missing local reflection.

Accordingly this note proves neither an unbounded reflection-maximal family
nor an absolute original-rank theorem. The exact remaining theorem is:

> In a literal finite-volume complete-nef carrier, prove from geometric
> reflection maximality that the mixed-root exceptional tail in (4.2) has
> uniformly bounded rank (or a finite numerical contraction orbit).

Until that bridge is proved, replacing the core by a census or invoking a
reflective-lattice classification would not be an AG solution of H1.
