# The pure norm-two OG10 Whittaker divisor is obstructed

## Result

Let

\[
F=\mathbf Q(\sqrt5),\qquad \mathcal O=\mathbf Z[\phi],\qquad
m_0=\frac{\phi}{\sqrt5},
\]

and let \(\rho=\rho_{L_F}\) be the 176-dimensional Hilbert--Weil
representation attached to the totally-real OG10 extension constructed in
`og10_borcherds_obstruction.md`.  Then

\[
\boxed{P_{m_0,0}=\xi_{(-4,6)}(f_{m_0,0})\ne0.}
\]

Consequently there is a cusp form

\[
g\in S_{(6,6),\rho}
\quad\text{with}\quad
b_g(m_0,0)\ne0.
\]

Equivalently, Bruinier's single Whittaker basis element \(f_{m_0,0}\) is
**not** weakly holomorphic.  Thus there is no Whittaker lift whose principal
part consists only of a nonzero coefficient at \((m_0,0)\).  This rules out
the direct single-term lift of the pure norm-two special divisor targeted in
the earlier report; it does not rule out a compensating principal part whose
total special divisor simplifies to the same divisor.

The proof is direct identity-term domination.  Normalize the diagonal
Fourier coefficient of \(P_{m_0,0}\) so that its identity term is 1.  The
exact rational certificate proves

\[
\left|a_{P_{m_0,0}}(m_0,0)-1\right|
<0.931739<1.
\]

No Hilbert cusp-form basis or dimension calculation is needed.

## 1. Normalizations and primary-source inputs

Bruinier defines the unitary Hilbert--Weil representation in equations
(3.5)--(3.8), the Poincare map \(\xi_k\) in (4.18)--(4.19), and proves its
surjectivity in Proposition 4.3.  Equation (4.17) shows that, up to a
nonzero common scalar, \(\xi_k(f_{m,\mu})\) is the vector-valued Hilbert
Poincare series indexed by \((m,\mu)\).  Proposition 4.5 gives

\[
\{g,f_{m,\mu}\}=b_g(m,\mu).
\]

These statements are in J. H. Bruinier,
[*Regularized theta lifts for orthogonal groups over totally real
fields*](https://arxiv.org/abs/0908.3076).

For the scalar Hilbert Poincare coefficient normalization and its
Kloosterman--Bessel expansion we use the general displayed Fourier formula
at the end of Section 2, and its diagonal specialization in the opening
display of Section 3, of M. Zhang and Y. Zhang,
[*On the non-vanishing of Hilbert Poincare
series*](https://doi.org/10.1016/j.jnt.2025.11.003) (the freely available
[author version](https://arxiv.org/abs/2312.04169) has the same formulas).
Their proof of Theorem 3.2 records the two bounds

\[
|J_5(x)|\le1,
\qquad
|J_5(x)|\le \frac{(x/2)^5}{5!}.
\]

For our translation lattice \(\mathcal O\), their fractional-ideal parameter
is \(\mathfrak c=\mathfrak d^{-1}\), so
\(N(\mathfrak c\mathfrak d)=1\).  Since \(D_F=5\), weight 6 gives the
scalar prefactor

\[
C_F=\frac{(2\pi)^2}{\sqrt5}.
\]

The vector-valued formula follows from exactly the same unfolding: replace
each scalar Kloosterman phase by that phase times the corresponding matrix
coefficient of \(\rho(\gamma)^{-1}\).  In the diagonal \((0,0)\) component,
write this generalized sum as

\[
K_c^\rho(m_0,\epsilon m_0;0,0).
\]

With \(n=|N(c)|\), the normalized coefficient therefore has the form

\[
1+C_F
\sum_{\epsilon\in\mathcal O^{\times,+}}
\sum_{c\in(\mathcal O\setminus\{0\})/\mathcal O^\times}
\frac{K_c^\rho(m_0,\epsilon m_0;0,0)}{n}
\prod_{i=1}^2
J_5\!\left(
\frac{4\pi\sqrt{m_{0,i}(\epsilon m_0)_i}}{|c_i|}
\right).
\tag{1}
\]

Only absolute values of the generalized Kloosterman sums are used below,
so all harmless inverse/complex-conjugate convention choices disappear.

## 2. Weil bounds: the unit row and all other rows

The representation \(\rho\) is unitary by Bruinier's statement following
equation (3.9).  A generalized Kloosterman sum for a fixed nonzero \(c\)
has at most \(|\mathcal O/(c)|=n\) summands, and every matrix coefficient
has absolute value at most 1.  Hence

\[
\boxed{\frac{|K_c^\rho|}{n}\le1.}
\tag{2}
\]

The unit row is much better.  Modulo all units it is represented by \(c=1\),
there is only one residue, and its matrix coefficient is the \((0,0)\)
entry of \(\rho(S)\), up to diagonal-unit and translation factors fixing
\(\mathfrak e_0\).  Bruinier's equation (3.6) gives

\[
\boxed{|K_1^\rho|=|\rho(S)_{0,0}|=\frac1{\sqrt{|L_F'/L_F|}}
=\frac1{\sqrt{176}}.}
\tag{3}
\]

This is the only representation-specific improvement needed.

## 3. The positive-unit orbit

The totally positive units are \(\epsilon=\phi^{2t}\), \(t\in\mathbf Z\),
and

\[
N(m_0)=\frac15.
\]

For a lower-left row of absolute norm \(n\), multiplication of its generator
by a unit shifts \(t\).  Thus its two Bessel arguments can be written

\[
x_1=\frac{s}{\sqrt n}\phi^{t+\delta},qquad
x_2=\frac{s}{\sqrt n}\phi^{-t-\delta},qquad
s=\frac{4\pi}{\sqrt5},qquad |\delta|\le\frac12.
\tag{4}
\]

Put

\[
A_n=\frac1{5!}\left(\frac{s}{2\sqrt n}\right)^5.
\]

The elementary Bessel bounds imply

\[
|J_5(x_1)J_5(x_2)|
\le
\min(1,A_n\phi^{5(t+\delta)})
\min(1,A_n\phi^{-5(t+\delta)}).
\tag{5}
\]

For \(n\ge4\), \(A_n<1\).  The right side of (5) is at most \(A_n^2\)
and also at most \(A_n\phi^{-5|t+\delta|}\).  Taking
\(K=\operatorname{bit\_length}(n)\), using the first bound for
\(|t|\le K\), the second on the tails, and the rational inequalities

\[
\phi>\frac85,qquad
A_1<\frac32,
\]

gives the entirely rational majorant

\[
U(n)=
\frac{9(2K+1)}{4n^5}
+\frac3{n^2}
\frac{(5/8)^{5K+2}}{1-(5/8)^5}.
\tag{6}
\]

For the unit row \(n=1\), the imbalance in (4) is an integral shift, so the
balanced term has argument \(s\) in both factors.  The defining alternating
series for \(J_5\), enclosed using

\[
\frac{333}{106}<\pi<\frac{355}{113},qquad
2.236067<\sqrt5<2.236068,
\]

proves exactly

\[
0.3327637<J_5(s)<0.3331040<\frac13.
\tag{7}
\]

For every other unit power, one Bessel factor is bounded by 1 and the other
by the small-argument inequality.  Using the sharper certified
\(A_1<73/50\),

\[
\sum_{t\in\mathbf Z}
|J_5(s\phi^t)J_5(s\phi^{-t})|
<
\frac19+\frac{2(73/50)}{(8/5)^5-1}
=\frac{37256}{88929}
<0.418941.
\tag{8}
\]

Combining (3), (8), and rational upper bounds for \(C_F\) and
\(\sqrt{176}\) yields a unit-row contribution smaller than

\[
0.557556.
\tag{9}
\]

## 4. All nonunit rows and the exact infinite tail

Rows \(c\) modulo \(\mathcal O^\times\) give principal integral ideals.
For an upper bound we may sum over **all** integral ideals.  The number of
ideals of norm \(n\) in \(F\) is

\[
a_F(n)=\sum_{d\mid n}\chi_5(d),
\tag{10}
\]

the coefficient identity coming from
\(\zeta_F(s)=\zeta(s)L(s,\chi_5)\).  In particular
\(a_F(2)=a_F(3)=0\), so every nonunit row has \(n\ge4\).

The certificate sums \(a_F(n)U(n)\) exactly as a rational number through
\(B=1024=2^{10}\).  For completeness, its infinite tail is not estimated
numerically.  Since

\[
K\le\log_2n+1,qquad
\left((5/8)^5\right)^K<n^{-3},qquad
a_F(n)\le d(n)\le n,
\]

the remaining sum is at most

\[
\int_B^\infty x^{-4}
\left[
\frac94(2\log_2x+3)
+\frac{3(5/8)^2}{1-(5/8)^5}
\right]dx.
\tag{11}
\]

Using \(\log B=10\log2\) and
\(\log2>2/3\), this is bounded by the exact fraction

\[
\frac{273187}{15914464444416}<1.7166\cdot10^{-8}.
\tag{12}
\]

The finite sum plus (12) is less than \(0.021194\).  Applying (2) and the
prefactor \(C_F\), the contribution from **every nonunit row** is less than

\[
0.374184.
\tag{13}
\]

Equations (9) and (13) give

\[
|E|<0.557556+0.374184<0.931739<1,
\]

which proves the theorem.

## 5. Cusps, ideal classes, and what is not being assumed

There is no hidden cusp or ideal-class truncation in the argument.

1. Bruinier's Poincare series is formed for the full
   \(\widetilde{\mathrm{SL}}_2(\mathcal O)\) from
   \(\widetilde\Gamma_\infty\backslash\widetilde\Gamma\).  Formula (1)
   is the Fourier coefficient at the same infinity cusp and sums every
   possible lower row of those cosets.
2. The passage from elements \(c\) modulo units to ideals only enlarges the
   absolute-value sum: principal ideals are bounded by all ideals using
   (10).  It therefore remains valid without a class-number assumption.
3. In this particular field there is in fact only one narrow ideal class:
   \(\mathbf Q(\sqrt5)\) has class number one, and the unit \(\phi\) has
   norm \(-1\), so its narrow and ordinary class groups agree.  Hence no
   additional cusp component is present anyway.

## 6. Exact logical consequence and scope

By Bruinier's Proposition 4.3,

\[
P_{m_0,0}=\xi(f_{m_0,0})\ne0.
\]

By Proposition 4.5, if every
\(g\in S_{(6,6),\rho}\) had \(b_g(m_0,0)=0\), then
\(f_{m_0,0}\) would pair trivially with the whole cusp space and would lie
in \(\ker\xi\), a contradiction.  Thus the coefficient functional

\[
g\longmapsto b_g(m_0,0)
\]

is nonzero.  This is the exact obstruction asked for; it is not merely a
proof that the cusp space itself is nonzero.

The conclusion is deliberately narrow:

- it rules out the **single** pure norm-two principal part;
- it does not rule out a finite compensating principal part in the kernel of
  all cusp coefficient functionals;
- although \(g_0=P_{m_0,0}\ne0\), it does not prove
  \(\Lambda(g_0)\ne0\), so it does not by itself rule out an exceptional
  scalar automorphic section with divisor \(Z(m_0,0)\);
- it does not imply anything in either direction about
  \([\Gamma_\tau:W_2]\), compactness, or finite covolume.

## Exact verification

The proof uses only Python's exact `Fraction` arithmetic:

```bash
PYTHONPATH=.:src python3 scripts/certify_og10_whittaker_nonvanishing.py
PYTHONPATH=.:src pytest -q tests/test_og10_whittaker_nonvanishing.py
```

The first command records every rational enclosure and the final exact
fraction, not just floating-point approximations.

The scalar-to-Weil unfolding normalization and the special unit-row factor
\(|\rho(S)_{0,0}|=1/\sqrt{176}\) were also replayed independently in a
hostile audit; no additional metaplectic magnitude factor occurs.
