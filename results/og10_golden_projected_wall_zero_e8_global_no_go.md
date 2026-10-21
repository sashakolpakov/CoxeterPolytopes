# Global zero-​\(E_8\) no-go for the OG10 one-root compiler

## Theorem

Keep

\[
F=\langle-2\phi\rangle\oplus E_8\oplus\langle2\rangle,
\qquad
H(\ell)=\begin{pmatrix}F&\ell\\ \ell^t&2\end{pmatrix},
\]

and restrict only to the zero-​\(E_8\) family
\(\ell=(A,0,C)\), with arbitrary \(A,C\in\mathbf Z[\phi]\).  If
\(p=F^{-1}\ell\), then no covector satisfying

\[
0<p^2<2,\qquad (p')^2<2,
\]

separates the archived invariant interior point from the archived ray
\(v_{61}\).  This is global in \(A,C\); there is no coefficient bound in the
statement.

## Exact height reduction

Put \(X=-A/\phi\).  Then

\[
q=2p^2=C^2-\phi X^2,
\qquad 0<q,q'<4.
\]

The two exact evaluations, with the archived denominator scalings cleared,
are

\[
400\ell(x_{\rm int})=-\phi(601+199\phi)X+2C,
\]

\[
13\ell(v_{61})=\phi(C-13X).
\]

For \(X=0\), these have the same sign.  For \(X\ne0\), opposite signs force

\[
13<\frac C X<\frac{\phi(601+199\phi)}2.
\]

In particular \(|C|>13|X|\), whence

\[
(169-\phi)X^2<q<4,
\qquad |X|<\frac16.
\]

At the other embedding,

\[
q'=C'^2+(\phi-1)X'^2<4,
\qquad |X'|<\frac83.
\]

Writing \(X=a+b\phi\), Minkowski reconstruction now gives

\[
|b|<\frac{17}{6\sqrt5}<2,
\qquad |a|<\frac16+\phi<2,
\]

so \(|a|,|b|\le1\).  The norm equation gives
\(|C|<13/6\) and \(|C'|<2\).  Thus for \(C=c+d\phi\),

\[
|d|<\frac{25}{6\sqrt5}<2,
\qquad |c|<\frac{13}{6}+\phi<4,
\]

and hence \(|c|\le3, |d|\le1\).  These are necessary bounds for a
separator, derived from the inequalities rather than imposed as a search
box.

The resulting 189 pairs contain exactly 26 two-place-signature-admissible
pairs.  Their exact sign histogram on `(interior,v61)` is

\[
13\ (-,-),\qquad 13\ (+,+),
\]

with no zero and no opposite-sign pair.  An independent implementation works
directly in the conservative covector box \(|A_a|,|A_b|\le2\),
\(|C_c|\le3,|C_d|\le1\).  Its 525-pair superset has 30 admissible pairs,
again split as 15 `(-,-)` and 15 `(+,+)` with no separator.

## Finite norm alphabet

Independently, the Minkowski rectangle \(0<q,q'<4\) contains exactly

\[
1,\ 1+\phi,\ 2-\phi,\ 2,\ 2+\phi,\ 3-\phi,\ 3.
\]

Up to multiplication by totally positive square units, these form four wall
classes, of rational norms \(1,4,5,9\).  Finiteness of this norm alphabet
does **not** make the mixed-​\(E_8\) vector problem finite; the theorem above
uses the cutting inequalities to obtain the extra height bound.

## Scope

This closes the zero-​\(E_8\) branch of the one-root ansatz globally.  It says
nothing global about mixed-​\(E_8\) covectors.  In particular, it neither
produces nor excludes a hyperbolic Coxeter polytope.

The exact certificate is
`data/candidates/og10_golden_projected_wall_zero_e8_global.json`, replayed by
`scripts/certify_og10_golden_projected_wall_zero_e8_global.py`; the
source-independent superset replay is
`scripts/independent_og10_golden_projected_wall_zero_e8_global.py`.
