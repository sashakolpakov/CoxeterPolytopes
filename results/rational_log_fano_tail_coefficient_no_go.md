# Log-pair coefficients do not absorb a free exceptional tail

**Status (2026-08-10): exact AG coefficient theorem and sharp no-go.**
This note asks whether the free exceptional \(B_m\) tail left by the H1 and
H2 adjoint-core theorems can be bounded by putting canonical coefficients on
the physical negative curves and applying a log-Fano, log-BMY, or BAB-style
theorem.  The answer is exact.

* On any fixed core, a tail of length \(m\) consumes at least \(m\) units of
  anti-log-canonical square.  Thus the *existence* of a nef and big
  \(-\!(K+\Delta)\) gives the sharp conditional inequality
  \(m<[-(K_{S_0}+\Delta_0)]^2\).
* Effective fractional coefficients never reduce this cost.  The unique
  square-optimal choice is to put coefficient zero on the whole tail.
* For the minimal hostile H1 and H2 incidences treated below, nefness of
  \(-\!(K+\Delta)\) already forces coefficient one on every component of
  the genus-one core.  The adjoint section then makes \(K+\Delta\)
  effective.  After one tail blowup it is nonzero, so its negative is not
  nef.  Even on a tail-free log-Calabi--Yau equality core it is zero, not
  big.

Consequently there is no effective boundary on the physical negative
curves which turns either minimal hostile pair, with a free tail away from
the pair, into a weak log-Fano pair; in particular there is no canonical
uniformly klt choice.  Higher, parallel, or dotted H1/H2 incidences are not
claimed here.  This uses only blowup formulas, adjunction, the canonical
exact sequence, and exact intersection theory.  No Nikulin--Vinberg
classification, ABSW theorem, arithmetic volume, or census is used.

## 1. The exact tail-energy identity

Let \(S_0\) be a smooth projective surface and let \(\Delta_0\) be an
effective boundary.  Choose a point outside \(\operatorname{Supp}\Delta_0\),
blow it up, and then make each new center a free point of the newest
exceptional curve, away from its intersection with the preceding curve.
Write

\[
                    \pi:S_m\longrightarrow S_0
\]

for the composite.  In the orthogonal total-transform basis
\(e_1,\ldots,e_m\),

\[
 e_i^2=-1,\qquad e_i.e_j=0\ (i\ne j),\qquad
 K_{S_m}=\pi^*K_{S_0}+\sum_{j=1}^m e_j.                    \tag{1.1}
\]

The physical exceptional curves are

\[
 R_i=e_i-e_{i+1}\quad(1\leq i<m),\qquad E=e_m.             \tag{1.2}
\]

Give them arbitrary boundary coefficients

\[
 \Delta_m=\pi^*\Delta_0+\sum_{i=1}^{m-1}a_iR_i+bE,
                   \qquad 0\leq a_i,b\leq1.                \tag{1.3}
\]

Put \(L_0=-\!(K_{S_0}+\Delta_0)\) and
\(L_m=-\!(K_{S_m}+\Delta_m)\).

### Theorem 1.1 (tail energy)

There are exact coefficients \(x_j\) such that

\[
 L_m=\pi^*L_0-\sum_{j=1}^m x_je_j,
 \qquad L_m^2=L_0^2-\sum_{j=1}^m x_j^2,                    \tag{1.4}
\]

where, for \(m\geq2\),

\[
 \begin{aligned}
 x_1&=1+a_1,\\
 x_j&=1-a_{j-1}+a_j &&(2\leq j<m),\\
 x_m&=1-a_{m-1}+b,
 \end{aligned}                                             \tag{1.5}
\]

and for \(m=1\), \(x_1=1+b\).  In every length,

\[
             \sum_{j=1}^m x_j=m+b.                         \tag{1.6}
\]

Consequently Cauchy--Schwarz gives

\[
 \boxed{\quad
 L_m^2\leq L_0^2-\frac{(m+b)^2}{m}\leq L_0^2-m.
 \quad}                                                     \tag{1.7}
\]

Equality in the resulting bound \(L_m^2\leq L_0^2-m\) holds exactly when
every tail coefficient is zero.

#### Proof

Substitute (1.2)--(1.3) into (1.1) and collect the orthogonal \(e_j\)
coordinates.  This gives (1.4)--(1.5).  All intermediate \(a_i\) cancel
when the \(x_j\) are summed, proving (1.6).  Cauchy--Schwarz gives
\(\sum x_j^2\geq(\sum x_j)^2/m\), and \(b\geq0\) gives (1.7).

If equality holds in the resulting outer bound, then \(b=0\) and equality
in Cauchy--Schwarz makes every \(x_j=1\).  Equation (1.5) successively forces
\(a_1=\cdots=a_{m-1}=0\).  Conversely the zero tail boundary has
\(x_1=\cdots=x_m=1\), so equality is attained. \(\square\)

The local nef inequalities are also transparent:

\[
             L_m.R_i=x_i-x_{i+1},\qquad L_m.E=x_m.          \tag{1.8}
\]

Thus the square-optimal zero boundary is locally nef on the tail:
\(L_m.R_i=0\) and \(L_m.E=1\).  Fractional coefficients cannot buy back
square hidden in a different tail direction.

### Corollary 1.2 (the exact conditional rank bound)

If \(L_m\) is nef and big, then \(L_m^2>0\), and hence

\[
                       \boxed{m<L_0^2}.                     \tag{1.9}
\]

Equivalently, \(m\leq\lceil L_0^2\rceil-1\).  This conclusion does not
need a uniform discrepancy gap: it holds for all effective coefficients in
the full interval \([0,1]\).  It is sharp as an intersection inequality.
What is missing in H1/H2 is not the numerical consequence of a weak
log-Fano structure; it is the weak log-Fano structure itself.

## 2. The minimal hostile H1 and H2 incidences force coefficient one

The obstruction occurs before any global BMY theorem can be applied.  Let
\(\Gamma\) be an effective divisor having none of the displayed hostile
curves as a component.  Intersections of \(\Gamma\) with those curves are
nonnegative.

### Proposition 2.1 (H1 coefficient rigidity)

Suppose

\[
 B_1^2=B_2^2=-4,\qquad B_1.B_2=2,
\]

and write

\[
                 \Delta=a_1B_1+a_2B_2+\Gamma,
                         \qquad 0\leq a_1,a_2\leq1.
\]

If \(-\!(K_S+\Delta)\) is nef, then

\[
                 a_1=a_2=1,\qquad \Gamma.(B_1+B_2)=0.       \tag{2.1}
\]

Indeed, adjunction gives \(K.B_i=2\).  Nefness, tested on \(B_1,B_2\),
gives

\[
 \begin{aligned}
 2-4a_1+2a_2+\Gamma.B_1&\leq0,\\
 2+2a_1-4a_2+\Gamma.B_2&\leq0.
 \end{aligned}                                             \tag{2.2}
\]

Adding them yields
\(4-2(a_1+a_2)+\Gamma.(B_1+B_2)\leq0\).  The coefficient caps and
nonnegativity force (2.1).

### Proposition 2.2 (H2 coefficient rigidity)

Suppose

\[
 B^2=-4,\qquad R^2=-2,\qquad B.R=2,
\]

and write

\[
                    \Delta=aB+bR+\Gamma,
                         \qquad0\leq a,b\leq1.
\]

If \(-\!(K_S+\Delta)\) is nef, then

\[
                      a=b=1,\qquad \Gamma.(B+R)=0.           \tag{2.3}
\]

Here adjunction gives \(K.B=2\) and \(K.R=0\).  The two nef tests are

\[
 \begin{aligned}
 2-4a+2b+\Gamma.B&\leq0,\\
 2a-2b+\Gamma.R&\leq0.
 \end{aligned}                                             \tag{2.4}
\]

Their sum first forces \(a=1\) and \(\Gamma.(B+R)=0\); the second then
forces \(b=1\).

These arguments remain valid when \(\Gamma\) contains every other
physical negative curve, including an exceptional tail.  Such curves can
only add nonnegative intersection to (2.2) or (2.4).

### Theorem 2.3 (no minimal-hostile weak log-Fano boundary)

Let \(S_0\) be a smooth projective rational surface.  Let
\(D_0=B_1+B_2\) be the transverse minimal H1 pair, or let \(D_0=B+R\) be
the minimal H2 pair, on \(S_0\).  Let \(S_m\to S_0\) be any free
exceptional tail blown up away from \(D_0\), and write \(D\) for its strict
transform.  There is no effective boundary on the physical negative
curves, with coefficients at most one, for which
\(-\!(K_{S_m}+\Delta)\) is nef and big.  If \(m\geq1\), it is not even
nef.

#### Proof

The canonical exact-sequence argument for both hostile pairs gives an
effective adjoint divisor

\[
                         A_0\in|K_{S_0}+D_0|,               \tag{2.5}
\]

with \(A_0=0\) allowed.  Propositions 2.1--2.2 show that anti-log-canonical
nefness forces coefficient one on every component of \(D\).  Hence
\(\Delta=D+\Gamma\) with \(\Gamma\geq0\).  The blowup formula gives

\[
 K_{S_m}+\Delta
   \sim_{\mathbf R}\pi^*A_0+\sum_{j=1}^m e_j+\Gamma.       \tag{2.6}
\]

Here each \(e_j\) denotes its effective total-transform divisor.  The right
side is effective, and it is nonzero when \(m\geq1\).  Its
intersection with every ample divisor is then positive, so its negative
cannot be nef.  When \(m=0\), nefness of the negative similarly forces
\(A_0=\Gamma=0\); the anti-log-canonical divisor is zero and is not big.
\(\square\)

In particular, the hostile coefficients are forced to be one, so the pair
cannot have any positive uniform klt gap.  The failure is structural, not
a poor choice of fractional coefficient.

## 3. Discrepancies do not see the length

The same calculation explains why discrepancy control alone cannot replace
bigness.  In the physical prime basis,

\[
 \sum_{j=1}^m e_j
   =\sum_{i=1}^{m-1}iR_i+mE,
\]

and therefore

\[
 K_{S_m}+\Delta_m
  =\pi^*(K_{S_0}+\Delta_0)
    +\sum_{i=1}^{m-1}(i+a_i)R_i+(m+b)E.                    \tag{3.1}
\]

All relative coefficients are positive.  A crepant pullback would require
\(a_i=-i\) and \(b=-m\), which is not an effective boundary.  With zero
tail boundary the successive exceptional valuations over the smooth point
have discrepancies \(1,2,\ldots,m\), rather than discrepancies tending
toward \(-1\).  On the source, taking every tail coefficient zero preserves
whatever positive discrepancy gap the core already has.  Thus arbitrarily
long smooth blowup towers are compatible with uniform discrepancy control;
they fail the weak log-Fano hypothesis through (1.7).

This is exactly why boundedness of uniformly \(\varepsilon\)-lc log-Fano
pairs is not a tail theorem here.  That boundedness requires the log-Fano
positivity hypothesis; it does not bound all smooth higher birational
models of a fixed pair from discrepancy control alone.  Equation (1.7) is
the elementary surface identity which detects the missing positivity
assumption.

Equally, contracting the whole tail returns the same pair
\((S_0,\Delta_0)\) for every \(m\).  Any orbifold Chern number or discrepancy
bound computed only on that contracted target is literally independent of
the carrier tail length.  Turning such a target bound into a carrier-rank
bound would already require the missing equivariant descent or
minimalization theorem.

## 4. Log-Chern audit

There is also no hidden upper bound in the usual log-BMY direction.  Assume
for this paragraph that \(m\geq1\), that \(D_0\) is a reduced SNC boundary,
let \(D_m\) be its strict transform, let
\(T_m=R_1+\cdots+R_{m-1}+E\), and put
\(c_{2,\log}(S,D)=e(S\setminus D)\).  If
\(c_0=e(S_0\setminus D_0)\), then the blowup and complement formulas give

\[
\begin{array}{c|c|c}
 \text{boundary}&(K+\text{boundary})^2&c_{2,\log}\\ \hline
 D_m&(K_{S_0}+D_0)^2-m&c_0+m\\
 D_m+T_m&(K_{S_0}+D_0)^2-(m+3)&c_0-1.
\end{array}                                             \tag{4.1}
\]

For the second row, removing the whole exceptional fibre identifies the
complement with \(S_0\setminus(D_0\cup\{p\})\); this proves the constant
log Euler number.  Its square formula is Theorem 1.1 with all tail
coefficients one.  In either row, increasing \(m\) decreases the log
\(c_1^2\) side while leaving the log \(c_2\) side constant or increasing
it.  Therefore an inequality of the form \(c_{1,\log}^2\leq3c_{2,\log}\)
becomes easier, not harder, along the tail.  Only an independently proved
nef-and-big hypothesis turns the falling square into the bound (1.9).

## 5. Exact conclusion

The coefficient route has one real theorem and one sharp terminal
obstruction:

\[
 \boxed{
 \begin{gathered}
 -(K+\Delta)\ \text{nef and big on a length-}m\text{ tail}
       \Longrightarrow m<[-(K_0+\Delta_0)]^2,\\
 \text{but the minimal hostile H1/H2 incidences admit no such boundary.}
 \end{gathered}}
\]

Thus log-pair coefficients, uniform discrepancies, log-BMY, and BAB do not
close the exceptional-tail gate even for the minimal H1/H2 incidences.
This note makes no conclusion for the higher, parallel, or dotted incidence
branches.  A successful dimension theorem still needs a global input which
actually contracts or bounds the physical tail (for example total marked
degree, a finite numerical contraction orbit, or an equivariant
minimalization theorem).  The exact arithmetic is replayed by
`scripts/certify_rational_log_fano_tail_coefficient_no_go.py` and frozen in
`data/enumerations/rational_log_fano_tail_coefficient_no_go.json`.
