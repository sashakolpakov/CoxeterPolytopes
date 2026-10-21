# The golden H9 and OG10 H10 chambers have infinite volume

## Theorem

Let

\[
L_9=\langle-2\phi\rangle\oplus E_8\oplus\langle2\rangle
\]

over \(\mathcal O=\mathbb Z[\phi]\), and let \(P_9\) be the complete
norm-two chamber obtained as the anti-symplectic fixed slice of the golden
OG10 chamber.  Then

\[
\boxed{\operatorname{vol}(P_9)=\infty.}
\tag{1}
\]

In fact, the stabilizer of \(P_9\) in \(O^+(L_9)\) contains an infinite
dihedral group.  Consequently the complete golden norm-two chamber
\(C\subset\mathbb H^{10}\) of the OG10 carrier also has infinite volume:

\[
\boxed{\operatorname{vol}(C)=\infty.}
\tag{2}
\]

For the projective marked OG10 carrier \((X,\eta,\tau)\) used in the
enhanced-covolume theorem, this further implies

\[
\boxed{|\operatorname{Bir}(X,\tau)|=\infty,}
\tag{3}
\]

where

\[
\operatorname{Bir}(X,\tau)=
\{f\in\operatorname{Bir}(X):f^*\tau=\tau f^*
\text{ on }NS(X)\}.
\]

Equation (3) is an existential consequence of the finite-covolume
equivalence.  The explicit H9 transformations below are not claimed to lift
to birational maps of this marked H10 carrier.

## 1. The Bugaenko reversal

The complete \(u=0\) facet of \(P_9\) is Bugaenko's compact eleven-wall
polytope \(P_8\).  Its labeled Coxeter diagram has a nontrivial reversal

\[
\pi=(0\ 8)(1\ 7)(2\ 6)(3\ 5)(9\ 10),
\tag{4}
\]

with nodes \(4\) fixed.  The first nine walls in the archived Bugaenko order
form an \(\mathcal O\)-basis: their coordinate determinant is
\(2\phi-3\), of field norm \(-1\).  Therefore (4) determines an integral
isometry \(A_8\) of
\(\langle-2\phi\rangle\oplus E_8\).  In the actual fixed-slice coordinates
\((x,z_1,\ldots,z_8)\), it is

\[
A_8=\begin{pmatrix}
9+16\phi&0&0&0&0&0&0&0&-2-3\phi\\
30+50\phi&0&0&0&0&0&0&1&-6-10\phi\\
60+100\phi&0&0&0&0&0&1&0&-12-20\phi\\
90+150\phi&0&0&0&0&1&0&0&-18-30\phi\\
72+120\phi&0&0&0&1&0&0&0&-15-24\phi\\
54+90\phi&0&0&1&0&0&0&0&-12-18\phi\\
36+60\phi&0&1&0&0&0&0&0&-8-12\phi\\
18+30\phi&1&0&0&0&0&0&0&-4-6\phi\\
48+80\phi&0&0&0&0&0&0&0&-9-16\phi
\end{pmatrix}.
\tag{5}
\]

Exact multiplication verifies

\[
A_8^2=1,
\qquad A_8^tH_8A_8=H_8,
\tag{6}
\]

and verifies all eleven oriented wall identities
\(A_8r_i=r_{\pi(i)}\), including the two walls not in the chosen basis.
Extend it to

\[
A=A_8\oplus 1_{\langle2\rangle}\in O(L_9).
\tag{7}
\]

## 2. Why the reversal preserves the complete chamber

This step does not depend on the 61-wall prefix.  Because \(A\) is an
integral isometry, it permutes the full set of norm-two mirrors.  Equations
(4)--(7) show that it maps the complete Bugaenko facet

\[
F=P_9\cap\{u=0\}
\]

to itself, with all oriented facet inequalities preserved.  Moreover \(A\)
fixes the last coordinate, so it preserves the side \(u\geq0\) occupied by
\(P_9\).

Choose a point in the relative interior of \(F\).  Local finiteness of the
root arrangement gives a neighborhood in which \(u=0\) is the only mirror
through that point.  Hence exactly one arrangement chamber is adjacent to
\(F\) on either side.  Both \(P_9\) and \(A(P_9)\) are adjacent to \(F\) on
the \(u\geq0\) side, so

\[
A(P_9)=P_9.
\tag{8}
\]

This is the global facet-adjacency argument that upgrades the diagram
reversal to a complete-chamber symmetry.

## 3. An infinite-order chamber stabilizer

The independently certified Pell-fold symmetry is

\[
g=\begin{pmatrix}\phi&-1\\ \phi&-\phi\end{pmatrix}_{x,u}
\oplus I_8,
\qquad g^2=1,
\qquad g(P_9)=P_9.
\tag{9}
\]

Exact multiplication of (5), (7), and (9) gives

\[
\operatorname{tr}(Ag)=8+8\phi>10.
\tag{10}
\]

If a real \(10\times10\) matrix has finite order, all its eigenvalues are
roots of unity, and therefore the absolute value of its trace is at most
\(10\).  Equation (10) proves that \(Ag\) has infinite order.  Since \(A\)
and \(g\) are involutions, they generate an infinite dihedral subgroup of
the complete chamber stabilizer.

## 4. Infinite index and infinite volume

Put

\[
\Gamma_9=O^+(L_9),
\qquad
W_9=\langle s_r:r\in L_9,\ r^2=2\rangle.
\]

The conjugate form of \(L_9\) is positive definite, so the form is
anisotropic over \(\mathbb Q(\sqrt5)\).  Thus \(\Gamma_9\) is a uniform
arithmetic lattice.  The group \(W_9\) is normal in \(\Gamma_9\), and
standard reflection theory makes it simply transitive on the chambers of
the complete mirror arrangement.  Consequently

\[
\Gamma_9/W_9\cong\operatorname{Stab}_{\Gamma_9}(P_9).
\tag{11}
\]

The right side is infinite by (8)--(10), hence \(W_9\) has infinite index.
Covolume multiplicity in the finite-volume quotient
\(\Gamma_9\backslash\mathbb H^9\) now gives (1).  In particular \(P_9\) is
not compact.

## 5. OG10 consequences and scope

The fixed-facet theorem proved the one-way gate

\[
\operatorname{vol}(C)<\infty\Longrightarrow P_9\text{ compact}.
\tag{12}
\]

Combining (1) with (12) proves (2).  The enhanced-covolume theorem for the
same marked projective OG10 carrier states

\[
\operatorname{vol}(C)<\infty
\iff |\operatorname{Bir}(X,\tau)|<\infty,
\]

so (2) proves (3).

This settles the covolume question for this algebraic carrier negatively:
it does not produce a compact or finite-volume Coxeter polytope in dimension
9 or 10.  It does not prove a general nonexistence theorem for hyperbolic
Coxeter polytopes in either dimension.  Nor does it identify an explicit
infinite-order element of \(\operatorname{Bir}(X,\tau)\); that infinitude is
forced through the already proved arithmetic-monodromy equivalence.

The exact certificate is
[`og10_golden_h9_infinite_stabilizer.json`](../data/candidates/og10_golden_h9_infinite_stabilizer.json),
recomputed by
[`certify_og10_golden_h9_infinite_stabilizer.py`](../scripts/certify_og10_golden_h9_infinite_stabilizer.py).
