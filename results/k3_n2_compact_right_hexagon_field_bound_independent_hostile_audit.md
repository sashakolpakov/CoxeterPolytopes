# Independent hostile audit: compact right-angled K3 hexagons

**Verdict (2026-08-11): PASS in the exact compact right-angled scope.**
The owner theorem needs neither a bound on every hexagon modulus nor an
assumption that an arbitrary pair cycle generates the field.  Rank-three
right-hexagon closure forces at least one second-neighbour pair cycle into
the interval \((4,68)\) at the Lorentzian place, while all its definite-place
conjugates lie in \((0,4)\).  That unique-large-conjugate pattern makes this
single integral cycle a primitive generator and proves

\[
                       |\operatorname {Disc}k|\le 68^{30}.
\]

The independently pinned owner bytes are:

    report       5cbbe7fd2d88a2129bdc0b5345f0ae377c0fa02b5f0874606fc6ff76a17783f3
    certifier    ff39e56f30c74429fb93935e3ba3f5353eeb219d059b2b38aa78c5f3df07a85a
    archive      dc52fb561ac5a16b7ee0683af1999be26909820102fabf92382b30c1cf856e10
    owner tests  06b0c81155739626843543301732e1af75c56dd2c5f9bd87c7bafab1674c808c
    scope audit  ba74f1f05d7e68d0d2cf916991dd9f248ce96127febe928f54bdf1274c672365
    scope tests  83895d73bddd5e098f53373ece083b6da6407439b08947bb58bf3a9a20f56952

## 1. Independent right-hexagon closure

Let \(u,v,w>1\) be the distinguished coshes of three alternating side
lengths, and put

\[
 x_0=u^2-1,\qquad x_2=v^2-1,\qquad x_4=w^2-1.
\]

The rank-three cofactor calculation, performed directly on the three
alternating normals, gives

\[
\begin{aligned}
 S&=u^2+v^2+w^2-1+2uvw,\\
 x_1&={S\over x_0x_2},&
 x_3&={S\over x_2x_4},&
 x_5&={S\over x_4x_0}.
\end{aligned}
\]

For example,

\[
 x_1+1
 =\frac{(uv+w)^2}{(u^2-1)(v^2-1)}.
\]

Expanding the numerator after subtracting one recovers \(S/(x_0x_2)\)
exactly.  Rotation gives the other two identities, and direct cancellation
then yields every four-consecutive minor relation

\[
                         x_i x_{i+1}=x_{i+3}x_{i+4}.
\]

The square roots are used only as positive real coordinates at the
distinguished embedding.  The proof never asserts that \(u,v,w\) belong to
the returned field.

## 2. The short-pair dichotomy is exhaustive

If one of \(x_0,x_2,x_4\) is below \(16\), the required short pair already
exists.  Otherwise order them as

\[
                         16\le x\le y\le z.
\]

Choose the odd coordinate with denominator \(yz\).  With
\(A=x+1,C=y+1,E=z+1\), monotonicity of
\((t+1)/t\) and \(\sqrt{t+1}/t\) gives

\[
\begin{aligned}
x_{\mathrm{odd}}
&<{2C+E+2C\sqrt E\over yz}\\
&\le {51+34\sqrt{17}\over256}
 <{221\over256}<1.
\end{aligned}
\]

Thus some one of all six coordinates always satisfies \(0<x_i<16\).
The independent executable checks all permutations of the alternating
ordering and exact rational cofactor samples; it does not infer a bound on
the other five coordinates.

## 3. Integrality and primitive field generation

For the short second-neighbour pair, let

\[
                         p=4(x_i+1).
\]

Writing \(P_j=1-s_j\), the finite-projective module trace
\(\operatorname {Tr}_{\mathcal O}(P_iP_{i+2})\) equals \(p\).
Consequently \(p\in\mathcal O\), without primitive normal scaling or
integrality of individual directed Cartan entries.

Compact ultraparallelity gives \(4<\sigma_0(p)<68\).  At each
nondistinguished embedding, positive-definite Cauchy--Schwarz and
nonvanishing of the dotted pairing give

\[
                              0<\sigma(p)<4.
\]

Put \(F=\mathbf Q(p)\).  If \(F\subsetneq k\), the restriction of
\(\sigma_0\) to the totally real field \(F\) has exactly \([k:F]\ge2\)
extensions to \(k\).  Since \(k\) is totally real, all those extensions are
real and all take the same value \(>4\) on \(p\).  One is
nondistinguished, contradicting the preceding interval.  Hence

\[
                              \mathbf Q(p)=k.
\]

This also rules out a gratuitous proper totally real extension above this
one-Lorentz-place packet.  The owner's minimal-cyclic-ground-field
hypothesis is retained, but the primitive-pair step does not secretly
assume pair generation.

## 4. Degree, discriminant, and scope

The primitive full trace has rank \(3d=\rho(X)\le20\), so \(d\le6\).
Every conjugate of the primitive algebraic integer \(p\) lies in
\((0,68)\).  Therefore every conjugate difference is strictly below \(68\),
and the Vandermonde gives

\[
\begin{aligned}
 |\operatorname {Disc}k|
 &\le |\operatorname {Disc}\mathbf Z[p]|\\
 &<68^{d(d-1)}
 \le68^{30}.
\end{aligned}
\]

Hermite--Minkowski gives a literal finite field set.  This result is
conditional on the primitive exact full-trace complex-K3 carrier and the
exact compact right-angled hexagon eigen-nef chamber.  It does not:

- construct that carrier or prove essential surjectivity;
- cover non-right or noncompact hexagons;
- cover arbitrary six-facet or higher-facet chambers;
- bound every dotted pair, a canonical mixed intersection, integral hosts,
  complete Gram packets, or chambers; or
- use reflection maximality, saturation, or a Coxeter/reflection-lattice
  classification.

The independent replay is
[audit_k3_n2_compact_right_hexagon_field_bound_independent_hostile.py](../scripts/audit_k3_n2_compact_right_hexagon_field_bound_independent_hostile.py),
with tests in
[test_k3_n2_compact_right_hexagon_field_bound_independent_hostile.py](../tests/test_k3_n2_compact_right_hexagon_field_bound_independent_hostile.py).
