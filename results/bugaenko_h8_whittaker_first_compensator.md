# The first complete Bugaenko \(\mathbb H^8\) compensator shell is obstructed

**Status (2026-08-08): exact bounded theorem and machine-checkable
certificate.**

## Result

Start only from the marked compact Bugaenko polytope and its canonical
determinant-line extension

\[
 \widehat H=\langle-2\phi\rangle\oplus E_8\oplus
 \langle-4\phi\rangle,
 \qquad
 Q=\frac{\phi\widehat H}{2\sqrt5},
\]

over \(\mathcal O=\mathbf Z[\phi]\).  Its discriminant module is

\[
 A=\mathcal O/(2)\oplus\mathcal O/(4)_{\mathrm{aux}}.
\]

Among **all 64 discriminant components**, the complete globally smallest
totally positive Whittaker-index shell, modulo positive units and
\(\mu\mapsto-\mu\), is the single orbit represented by

\[
 \boxed{
 m_1=\frac{\phi}{4\sqrt5},\qquad
 N(m_1)=\frac1{80},\qquad
 \mu_1=\left(\phi,2\phi^2\right).
 }
 \tag{1}
\]

The orbit contains three discriminant components.  Every one is
self-opposite and has nonzero auxiliary coordinate modulo \(4\).  Therefore
its complete, nonempty Heegner divisor has exact real-slice disposition

\[
 \boxed{
 0\text{ facets}+0\text{ lower-face contacts}
 +0\text{ interior cuts};\quad
 \text{the entire divisor is disjoint from the slice}.}
 \tag{2}
\]

Let \(P_0\) be the weight-\((5,5)\) Hilbert--Weil Poincare series for the
target

\[
 (m_0,\mu_0)=\left(\frac{\phi}{\sqrt5},0\right),
\]

and let \(P_1\) be the series for (1).  A certified coefficient-matrix
calculation proves

\[
 \boxed{P_0\text{ and }P_1\text{ are linearly independent}.}
 \tag{3}
\]

Consequently no weakly holomorphic Whittaker principal part supported on the
target and this **complete globally minimal compensator shell** can have a
nonzero target coefficient.  This closes the first finite compensator
branch after the pure target obstruction.

It is a bounded nonexistence theorem.  The later determinant safe-span
theorem proves that some larger finite slice-disjoint compensator exists, but
this calculation neither locates it nor supplies an independent support
bound.

## 1. The finite shell is computed from the polytope carrier

Write

\[
 \mu(x,y)=\frac{x}{2}u+\frac{y}{4}e,
 \qquad x\in\mathcal O/(2),\quad y\in\mathcal O/(4).
\]

The finite quadratic form from the
[`pure-target certificate`](bugaenko_h8_whittaker_gate.md) is

\[
 \sqrt5Q(\mu(x,y))\equiv
 -\frac{\phi^2x^2}{4}-\frac{\phi^2y^2}{8}
 \pmod{\mathcal O}.
 \tag{4}
\]

For an admissible index put

\[
 w=8\sqrt5m\in\mathcal O.
\]

Then (4) is the finite congruence

\[
 w\equiv-\phi^2(2x^2+y^2)\pmod{8\mathcal O}.
 \tag{5}
\]

Exact enumeration of the \(4\cdot16=64\) pairs gives 40 distinct residues.
Their field-norm residues modulo \(8\) are exactly

\[
 \{0,1,3,4\}.
 \tag{6}
\]

Because \(m\) is totally positive while \(\sqrt5\) has signs \((+,-)\),
the element \(w\) has signs \((+,-)\), so \(N(w)<0\), and

\[
 N(m)=-\frac{N(w)}{8^2\cdot5}.
 \tag{7}
\]

The possibilities \(-N(w)=1,2,3\) are excluded immediately by (6).  The
next value is realized by

\[
 w_1=2\phi,\qquad N(w_1)=-4,
\]

which gives (1).  There is a unique ideal of norm four in
\(\mathbf Q(\sqrt5)\), namely \((2)\), and
\(\mathcal O^\times=\{\pm\phi^j\}\).  Hence every sign-\((+,-)\) element
of norm \(-4\) is

\[
 2\phi^{2t+1}.
\]

Modulo \(8\mathcal O\), this orbit has the three residues

\[
 (0,2),\quad(2,4),\quad(6,2).
\]

Each residue has exactly one discriminant component satisfying (5).  The
complete component orbit is

\[
 \begin{aligned}
  &(\phi,2\phi^2),\\
  &(\phi^2,2),\\
  &(1,2\phi),
 \end{aligned}
 \tag{8}
\]

where the first coordinate is reduced modulo \(2\) and the second modulo
\(4\).  Multiplication by \(\phi\) cycles (8) and multiplies the index by
\(\phi^2\).  All three components are 2-torsion, so quotienting by
\(\mu\mapsto-\mu\) does not split or enlarge the orbit.

Thus (1) is the complete global minimum, not a minimum chosen inside one
convenient component.

## 2. The shell is nonempty but invisible on the real slice

Let \(r\) be any \(E_8\) root, \((r,r)=2\).  The vector

\[
 \lambda=\frac\phi2u+\phi^2r+\frac{\phi^2}{2}e
 \in\mu_1+\widehat L
 \tag{9}
\]

satisfies

\[
 \widehat H(\lambda,\lambda)
 =-2\phi\left(\frac\phi2\right)^2
 +2\phi^4
 -4\phi\left(\frac{\phi^2}{2}\right)^2
 =\frac12.
\]

Therefore \(Q(\lambda)=m_1\), proving that the global Heegner divisor in
(1) is not empty.

On the determinant-line real slice, a special hyperplane can meet the slice
only when its defining vector has zero coefficient on the auxiliary line
\(e\).  Every vector in the coset \(\mu(x,y)+\widehat L\) has auxiliary
coefficient

\[
 \frac y4+a,\qquad a\in\mathcal O.
\]

This can vanish exactly when \(y=0\) in \(\mathcal O/(4)\).  Every component
in (8) has \(y\ne0\).  This proves the whole-divisor disjointness in (2),
not merely the absence of one chosen vector from the slice.

This is the most favorable geometric kind of compensator: it can affect the
global cusp obstruction without adding any real mirror to the Bugaenko
chamber.

## 3. The two-column Poincare matrix

For \(i,j\in\{0,1\}\), let

\[
 C_{ij}=b_{P_j}(m_i,\mu_i).
\]

Both diagonal entries have identity term one, while both off-diagonal
identity terms vanish because the two indices lie in distinct unit orbits.
The scalar weight-five coefficient formula has outer norm factor

\[
 \left(\frac{N(m_i)}{N(m_j)}\right)^2.
\]

Consequently

\[
 \widetilde C=
 \operatorname{diag}(N(m_i)^2)^{-1}
 C\operatorname{diag}(N(m_i)^2)
 \tag{10}
\]

has no norm-ratio factors.  Similarity preserves invertibility.

For a unit lower-left row, every matrix entry of the finite Fourier
transform has magnitude \(1/\sqrt{|A|}=1/8\).  Arb interval arithmetic at
220-bit precision, with analytic tails after \(|t|=10\), gives

| Bessel orbit | certified upper bound |
|---|---:|
| minimal diagonal | \(<0.00037\) |
| target--minimal cross | \(<0.027\) |

For all nonunit rows, unitarity gives the conservative normalized
Kloosterman-row bound one.  If \(N\ge4\) is the ideal norm and \(A\) is the
small-argument prefactor, the certificate minimizes the exact rational
majorant

\[
 (2K+1)\left(\frac{A}{N^2}\right)^2
 +2\frac{A}{N^2}
 \frac{(5/8)^{4K+2}}{1-(5/8)^4}
 \tag{11}
\]

over finitely many \(K\), sums with the exact ideal count through
\(N=1024\), and bounds the remaining infinite tail separately.  Since

\[
 \frac{N(m_1)}{N(m_0)}=\frac1{16},
\]

the target prefactor \(<13/5\) scales by \(1/16^2\) on the minimal diagonal
and by \(1/16\) on a cross entry.  The certificate uses the rational
majorants \(1/96\) and \(1/6\), respectively.

To avoid any dependence on whether rows are vector-valued or quotiented by
the metaplectic center, every estimate involving the nonzero component is
multiplied by two.  Even with that conservative convention, the complete
entry bounds are

\[
 \begin{array}{c|c}
 \text{entry}&\text{certified error}\ \\ \hline
 (0,0)&<3/4\\
 (1,1)&<1/500\\
 (0,1),(1,0)&<1/6.
 \end{array}
 \tag{12}
\]

The sharp \((0,0)\) value is inherited and independently replayed from the
pure-target certificate; it includes the exact vanishing of the full
ideal-norm-four \((0,0)\) Kloosterman row.  The other three entries use no
special finite-Weil cancellation.

The row dominance margins of \(\widetilde C\) are therefore

\[
 1-\frac34-\frac16=\frac1{12}>0,
 \qquad
 1-\frac1{500}-\frac16=\frac{1247}{1500}>0.
 \tag{13}
\]

The Levy--Desplanques theorem proves \(\det\widetilde C\ne0\), hence
\(\det C\ne0\).  Bruinier's Poincare map identifies these columns, up to
nonzero normalizing scalars, with the cusp obstructions of the two
Whittaker basis orbits.  This proves (3).

The analytic inputs are Bruinier's unitary Weil representation, Poincare
map, and exact obstruction pairing in
[Bruinier, *Regularized theta lifts for orthogonal groups over totally real
fields*](https://arxiv.org/abs/0908.3076), together with the scalar
Kloosterman--Bessel coefficient formula in
[Zhang--Zhang, *On the non-vanishing of Hilbert Poincare
series*](https://arxiv.org/abs/2312.04169).

## 4. Exact scope

What is proved:

1. the finite module extracted from the polytope determines a complete
   globally minimal index shell across all 64 components;
2. that shell is nonempty as a global Heegner divisor but entirely disjoint
   from the canonical real slice;
3. modulo units and component sign, it gives exactly one compensator
   column; and
4. the target plus that column has no target-preserving Whittaker kernel.

What is not proved:

1. no intrinsic bound on the number or norm of compensators is known;
2. a larger finite safe kernel exists abstractly, but its least index and
   support remain unknown;
3. a sufficiently deep level can separate the facet classes, but the least
   level and its first positive correction remain uncomputed; and
4. this carrier obstruction changes nothing about the independently
   certified compactness of Bugaenko's \(\mathbb H^8\) polytope.

Thus the exact updated gate is

\[
 \boxed{
 \begin{array}{rcl}
 \text{pure target column}&:&\text{obstructed},\\
 \text{target + complete minimal shell}&:&\text{obstructed},\\
 \text{larger finite safe support}&:&\text{exists; least support unknown}.
 \end{array}}
\]

## 5. Reproduction

Run

```sh
PYTHONPATH=src:. python3 scripts/certify_bugaenko_h8_whittaker_first_compensator.py
PYTHONPATH=src:. pytest -q tests/test_bugaenko_h8_whittaker_first_compensator.py
```

The certifier checks all 64 discriminant classes, the exact global minimum,
the complete three-component unit orbit, an explicit representing vector,
the slice-disjointness criterion, and every analytic matrix inequality.

- certifier:
  [`scripts/certify_bugaenko_h8_whittaker_first_compensator.py`](../scripts/certify_bugaenko_h8_whittaker_first_compensator.py),
  SHA-256 `8877b4e12f418d074294a9da7d57cc9913329e8035283ef9dfe9b210249016a9`;
- tests:
  [`tests/test_bugaenko_h8_whittaker_first_compensator.py`](../tests/test_bugaenko_h8_whittaker_first_compensator.py),
  SHA-256 `19bc060bfe4f2a49524e6ae47ac76fbfd17028ff02e2e2704fc5c7a27e623202`;
- certificate:
  [`data/candidates/bugaenko_h8_whittaker_first_compensator.json`](../data/candidates/bugaenko_h8_whittaker_first_compensator.json).
