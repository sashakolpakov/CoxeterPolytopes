# The determinant-line Whittaker gate for Bugaenko's compact \(\mathbb H^8\) chamber

**Status (2026-08-08): exact theorem and machine-checkable certificate.**

## Result

Let

\[
 F=\mathbb Q(\sqrt5),\qquad \mathcal O=\mathbb Z[\phi],\qquad
 \phi=\frac{1+\sqrt5}{2},
\]

and let

\[
 H_8=\langle-2\phi\rangle\oplus E_8
\]

be the archived Bugaenko lattice of signature $((8,1),(9,0))$.  Its
canonical determinant-line extension in Bruinier's sign convention is

\[
 \widehat H=H_8\oplus\langle-4\phi\rangle,
 \qquad \operatorname{sig}(\widehat H)=((8,2),(10,0)).
 \tag{1}
\]

With

\[
 Q(x)=\frac{\phi}{2\sqrt5}\widehat H(x,x),
 \tag{2}
\]

the exact Hilbert discriminant module is

\[
 \boxed{
 A=\widehat L'/\widehat L
 \cong\mathcal O/(2)\oplus\mathcal O/(4),
 \qquad |A|=64,
 }
 \tag{3}
\]

of level $(8)$.  It is anisotropic as an $\mathcal O$-valued finite
quadratic module, so the determinant extension has no nontrivial even
overlattice and hence no unresolved gluing branch.

All eleven facet roots have the single Heegner type

\[
 (m_0,\mu_0)=\left(\frac{\phi}{\sqrt5},0\right).
 \tag{4}
\]

For $n=8$, Bruinier's Whittaker input and obstruction weights are exactly

\[
 k=(-3,5),\qquad \kappa=(5,5).
 \tag{5}
\]

Let $P_{m_0,0}\in S_{(5,5),\rho_A}$ be the Hilbert--Weil Poincare series
which is the $\xi_k$-image of the Whittaker basis element $f_{m_0,0}$.  The
certificate proves

\[
 \boxed{P_{m_0,0}\ne0.}
 \tag{6}
\]

Consequently

\[
 \boxed{
 q^{-\phi/\sqrt5}\mathfrak e_0
 \text{ is an obstructed Whittaker principal part.}
 }
 \tag{7}
\]

This decides the **pure full norm-two Heegner column**.  It does not decide
principal parts with additional compensating Heegner types.

This is a theorem-level limit on the inverse-transfer program: the canonical
determinant line does **not** automatically supply the pure facet/norm-two
Whittaker correction, even for this compact arithmetic chamber.  A viable
general transfer theorem must therefore allow bounded compensating Heegner
types, retain a genuinely non-Borcherds algebraic-carrier branch, or do both.

There is also an important support correction.  The eleven facets do not
give eleven independent Heegner columns: all have the type (4).  Moreover,
the same norm-two divisor has 309 further sheets touching proper lower faces
of the closed chamber.  Thus the exact slice count is

\[
 \boxed{
 11\text{ facet sheets}+309\text{ lower-face contacts}
 =320\text{ closed contacts},\qquad0\text{ interior cuts}.
 }
 \tag{8}
\]

Hence the eleven facets alone are **not** selectable by the standard
full-group/discriminant-module Heegner column $Z(m_0,0)$.  The pure
coefficient at (4) selects that full norm-two Heegner divisor, including all
309 lower-face contacts.  Those contacts do not spoil chamber completeness
because none cuts the interior.  This statement does not exclude orbit
splitting after passage to a finite-index level subgroup.

## 1. Why the determinant line is \(\langle-4\phi\rangle\)

The archived form $H_8$ has determinant $-2\phi$.  To use the determinant
construction in signature $(1,8)$, first negate it.  Since the rank is nine,

\[
 \det(-H_8)=(-1)^9\det(H_8)=2\phi.
\]

The determinant-line bilinear form is $2(-1)^8\det(-H_8)=4\phi$.
Negating the resulting type-IV form to Bruinier's $(8,2)$ convention gives
precisely (1).

The scalar $\phi/\sqrt5$ is totally positive and belongs to the
codifferent.  The trace bilinear lattice associated with (2) is exactly

\[
 E_8\oplus E_8\oplus\langle-2\rangle\oplus\langle2\rangle
 \oplus\langle-4\rangle\oplus\langle4\rangle.
 \tag{9}
\]

The certificate gives a determinant-$-1$ basis change to (9), and verifies

\[
 \operatorname{sig}(9)=(18,2),\qquad\det(9)=64.
 \tag{10}
\]

## 2. The precise Weil module and the absence of gluing

Let $u$ and $e$ generate the lines of coefficients $-2\phi$ and
$-4\phi$.  Since $E_8$ is unimodular and $\phi$ is a unit, representatives
of $A$ can be written

\[
 \mu(x,y)=\frac{x}{2}u+\frac{y}{4}e,qquad
 x\in\mathcal O/(2),\quad y\in\mathcal O/(4).
 \tag{11}
\]

Their finite quadratic form is

\[
 \boxed{
 \sqrt5\,Q(\mu(x,y))
 \equiv-\frac{\phi^2x^2}{4}-\frac{\phi^2y^2}{8}pmod{\mathcal O}.
 }
 \tag{12}
\]

The additive group is

\[
 A\cong(\mathbb Z/2)^2\oplus(\mathbb Z/4)^2.
\]

Equivalently, the square-valued discriminant form of the trace lattice,
modulo $2\mathbb Z$, is

\[
 \left\langle-\frac12\right\rangle
 \oplus\left\langle\frac12\right\rangle
 \oplus\left\langle-\frac14\right\rangle
 \oplus\left\langle\frac14\right\rangle.
 \tag{13}
\]

Exact enumeration of the $4\cdot16=64$ pairs in (11) finds 40 distinct
values of (12) and only one isotropic element, namely $(0,0)$.  Even
$\mathcal O$-overlattices correspond to isotropic $\mathcal O$-submodules of
$A$, so there is no nontrivial even overgluing.  Every norm-two reflection
acts trivially on $A$: for $z\in\widehat L'$ it changes $z$ by
$(z,r)r\in\widehat L$.  Thus no additional $W(P)$-stability condition
changes this conclusion.

## 3. One arithmetic column, not eleven

Every archived simple root $r_i$ has $H_8(r_i,r_i)=2$, and therefore

\[
 Q(r_i)=\frac{\phi}{\sqrt5}=m_0,qquad r_i\bmod\widehat L=0.
 \tag{14}
\]

The eleven local facet sheets consequently collapse to the one global
coefficient $c(m_0,0)$.  Multiplying that coefficient by eleven would only
give multiplicity eleven to **every** norm-two sheet; it would not select the
eleven facets.

The independently archived compact-face enumeration exhausts every
norm-two mirror meeting the closed chamber.  It gives the profile (8): 11
facets, 309 nonfacet contacts, and no interior cut.  This has two different
consequences which must not be conflated.

1. The full norm-two arrangement still has $P$ as a chamber, because no
   extra mirror enters $P^\circ$.
2. A divisor asserted to contain only the eleven facet sheets is not the
   standard full-group Heegner column $Z(m_0,0)$: that column contains all of
   its sheets.

The pure automorphic question is therefore exactly whether the one-term
principal part in (7) is weakly holomorphic.

## 4. Exact nonvanishing of the obstruction Poincare series

Bruinier's Proposition 4.5 identifies the obstruction pairing with Fourier
coefficients of cusp forms, while equations (4.17)--(4.19) identify
$\xi_k(f_{m_0,0})$, up to a nonzero common scalar, with $P_{m_0,0}$.  It is
therefore enough to prove (6).

Normalize the diagonal Fourier coefficient as

\[
 a_{P_{m_0,0}}(m_0,0)=1+E.
 \tag{15}
\]

The Hilbert Poincare coefficient formula has scalar prefactor

\[
 C_F=\frac{4\pi^2}{\sqrt5}
 \tag{16}
\]

and two Bessel factors of order $4$.  The vector-valued unfolding replaces
the scalar Kloosterman phase by the corresponding matrix coefficient of the
unitary Weil representation.

### 4.1. The unit row

Bruinier's finite Fourier transform gives

\[
 |\rho(S)_{0,0}|=|A|^{-1/2}=\frac18.
 \tag{17}
\]

The totally positive units are $\phi^{2t}$.  With
$s=4\pi/\sqrt5$, their absolute Bessel orbit is

\[
 U_1=\sum_{t\in\mathbb Z}
 \left|J_4(s\phi^t)J_4(s\phi^{-t})\right|.
\]

Arb interval arithmetic through $|t|\le10$, followed by the analytic tail

\[
 2\frac{(s/2)^4}{4!}
 \frac{\phi^{-44}}{1-\phi^{-4}},
\]

proves

\[
 U_1<0.268385051708902,qquad
 \frac{C_F}{8}U_1<0.592301824941505.
 \tag{18}
\]

The displayed decimals are outward-rounded summaries of certified balls,
not floating-point estimates.

### 4.2. Exact vanishing at ideal norm four

The prime $2$ is inert, so there is one ideal of norm four, represented by
$c=2$.  Its units $d$ modulo $2$ are $1,\phi,\phi+1$.  After harmless upper
translations and unit diagonals, the three $(0,0)$ Weil coefficients reduce
to lower-unipotent Gauss sums with multipliers

\[
 2,\qquad2\phi,\qquad2\phi^2.
\]

For each multiplier, exact enumeration of all 64 elements of $A$ gives the
same eighth-root phase histogram

\[
 \{0:16,\ 2:16,\ 4:16,\ 6:16\}.
\]

Hence every Gauss sum is

\[
 16(1+i-1-i)=0.
 \tag{19}
\]

Thus the complete norm-four Kloosterman row vanishes.  This exact
cancellation is the decisive improvement over a representation-free
unitarity bound.

### 4.3. Every remaining lower row and the infinite tail

For a lower-left entry $c$, put

\[
 A[c]=\{\mu\in A:c\mu=0\}.
\]

The explicit finite-Weil matrix formula gives

\[
 |\rho(\gamma)_{0,0}|
 \le\sqrt{\frac{|A[c]|}{|A|}}.
 \tag{20}
\]

For $A=\mathcal O/(2)\oplus\mathcal O/(4)$, exact enumeration gives

\[
 |A[c]|=1,16,64
\]

according as the $2$-adic valuation of $c$ is $0,1$, or at least $2$.
After division by the ideal norm, every generalized Kloosterman row is
therefore bounded by $1/8$, $1/2$, or $1$, respectively.

For an ideal of norm $N\ge4$, balance a generator so the two Bessel arguments
are

\[
 \frac{s}{\sqrt N}\phi^{t+\delta},\qquad
 \frac{s}{\sqrt N}\phi^{-t-\delta},qquad |\delta|\le\frac12.
\]

The primary inequalities

\[
 |J_4(x)|\le1,qquad
 |J_4(x)|\le\frac{(x/2)^4}{4!}
\]

and the rational bounds

\[
 \phi>\frac85,qquad
 \frac1{4!}\left(\frac{2\pi}{\sqrt5}\right)^4<\frac{13}{5}
\]

give, for every integer $K\ge0$, the orbit majorant

\[
 U(N,K)=
 (2K+1)\left(\frac{13}{5N^2}\right)^2
 +2\frac{13}{5N^2}
 \frac{(5/8)^{4K+2}}{1-(5/8)^4}.
 \tag{21}
\]

The certificate takes the best of finitely many valid $K$, weights (21) by
the exact ideal count

\[
 a_F(N)=\sum_{d\mid N}\chi_5(d),
\]

and by the three Weil bounds following (20), and sums through $N=1024$.
The exact rational result is less than $0.008407298361694$ before (16).

For $N>1024$, discard all Weil savings.  Using
$a_F(N)\le N$, $K=\operatorname{bit\_length}(N)$, and
$((5/8)^4)^K<N^{-2}$ gives the exact rational tail

\[
 \frac{2243027}{27996979200}<0.000080116750596.
 \tag{22}
\]

After (16), all nonunit rows other than the exactly zero norm-four row
contribute less than

\[
 0.149847824366840.
 \tag{23}
\]

Combining (18), (19), and (23), the final certified interval is

\[
 \boxed{
 |E|<0.742149649308345<1.
 }
 \tag{24}
\]

Equation (15) is nonzero, proving (6).  Equivalently, there exists
$g\in S_{(5,5),\rho_A}$ with $b_g(m_0,0)\ne0$, so Bruinier's exact pairing
criterion obstructs (7).  No cusp-space dimension computation is needed.

## 5. Exact scope of the verdict

The result proves all of the following.

1. The determinant-line extension and its even gluing are unique.
2. The eleven facets supply one, not eleven, Heegner principal-part columns.
3. A facet-only divisor is not selectable by the standard
   full-group/discriminant-module column because that same column includes
   the 309 lower-face contacts.
4. The pure full norm-two Whittaker principal part is obstructed.

It does **not** prove any of the following.

1. The **least** weakly holomorphic principal part containing
   $q^{-m_0}\mathfrak e_0$ together with additional compensating terms has
   been found, or lies in the first compensator shell.  In fact the direct
   determinant safe-span theorem now proves abstractly that some larger
   finite safe support exists, but does not locate it.
2. There is no arbitrary meromorphic modular form supported on the norm-two
   divisor; a compact totally-real converse theorem is unavailable.
3. The standard full-group column itself separates the eleven facets.  A
   sufficiently deep discriminant-kernel level is now proved to split the
   target cosets so that their closed-slice support is exactly the chosen
   facets; the present full-group computation does not determine the least
   such level.
4. Bugaenko's reflection group fails to be reflective or compact.  Its
   compactness is independently certified and is unaffected by this
   automorphic obstruction.

Thus the requested three-way verdict is:

\[
 \boxed{
 \begin{array}{rcl}
 \text{eleven facets via standard full-group column}&:&\text{not selectable},\\
 \text{pure full norm-two Whittaker term}&:&\text{obstructed},\\
 \text{larger safe compensation}&:&\text{exists abstractly; least support unknown}.
 \end{array}}
 \tag{25}
\]

## 6. Reproduction and hashes

Run

```sh
PYTHONPATH=src:. python3 scripts/certify_bugaenko_whittaker_gate.py
PYTHONPATH=src:. pytest -q tests/test_bugaenko_whittaker_gate.py
```

The test result at freeze time is `3 passed`.

The certifier reads the two frozen geometry archives rather than hard-coding
their assertions: it checks the eleven norm-two facet diagonals, compact
41-vertex/zero-cusp profile, and the `11+309=320`, zero-interior-cut contact
profile before accepting the analytic verdict.

- certificate:
  [`data/candidates/bugaenko_h8_whittaker_gate.json`](../data/candidates/bugaenko_h8_whittaker_gate.json);
- certifier:
  [`scripts/certify_bugaenko_whittaker_gate.py`](../scripts/certify_bugaenko_whittaker_gate.py),
  SHA-256 `d95a1b23854aa60017dbe09871bdd1db1b7734c6cecdffc9215398120da7c5d2`;
- tests:
  [`tests/test_bugaenko_whittaker_gate.py`](../tests/test_bugaenko_whittaker_gate.py),
  SHA-256 `caf7d96c200911d426590ea7c283cfd682506ff619331b23c7a685e9e31e71b7`;
- compact-polytope archive:
  [`data/known_examples/bugaenko_compact_h8.json`](../data/known_examples/bugaenko_compact_h8.json),
  SHA-256 `c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c`;
- 320-contact archive:
  [`data/trace_lattices/bugaenko_k3_eigen_nef.json`](../data/trace_lattices/bugaenko_k3_eigen_nef.json),
  SHA-256 `e77bb81bcc738f4f542676748dbcf09418cab830da579a600fd1e309f03de782`.

## Follow-up: the complete first compensator support

The pure-column theorem has now been extended one finite layer.  Exact
enumeration across all 64 components finds a unique globally minimal
three-component unit orbit of norm `1/80`.  Its nonempty Heegner divisor is
wholly disjoint from the real slice, but certified coefficient bounds make
the target-plus-shell `2x2` Poincare matrix invertible.  Hence there is no
target-preserving principal part on this complete first support.  Larger and
multi-shell supports remain open.  See the
[`first compensator theorem`](bugaenko_h8_whittaker_first_compensator.md).

## 7. Primary analytic sources

1. J. H. Bruinier,
   [*Regularized theta lifts for orthogonal groups over totally real fields*](https://arxiv.org/abs/0908.3076),
   J. Reine Angew. Math. **672** (2012), 177--222.  See equations
   (3.5)--(3.9), (4.17)--(4.19), Proposition 4.5, and Theorem 6.8.
2. M. Zhang and Y. Zhang,
   [*On the non-vanishing of Hilbert Poincare series*](https://arxiv.org/abs/2312.04169),
   J. Number Theory (2026),
   [doi:10.1016/j.jnt.2025.11.003](https://doi.org/10.1016/j.jnt.2025.11.003).
   The proof uses their Hilbert Kloosterman--Bessel coefficient formula and
   the stated elementary Bessel bounds.
3. F. Strömberg,
   [*Weil representations associated to finite quadratic modules*](https://arxiv.org/abs/1108.0202),
   Math. Z. **275** (2013), 509--527.  Its explicit matrix-coefficient
   formula is the classical finite-module form of the kernel-size estimate
   (20); in the present Hilbert module the same estimate is also obtained
   directly from Bruinier's finite Fourier transform, and the exceptional
   norm-four row is enumerated without an estimate.
