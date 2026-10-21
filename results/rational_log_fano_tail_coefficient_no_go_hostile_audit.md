# Hostile audit of the log-Fano tail-coefficient no-go

**Verdict (2026-08-10): PASS after scope and presentation repairs.**  The
tail-energy identity, its unique global equality case, the conditional rank
bound, coefficient rigidity for the two minimal hostile incidences, the
adjoint-effectivity obstruction, the discrepancy ladder, and both log-Chern
rows survive independent reconstruction.  The result does not cover higher,
parallel, or dotted H1/H2 incidence branches and does not itself bound an
arbitrary exceptional tail.

The audit uses total-transform blowup coordinates, adjunction, the canonical
exact sequence, positivity against an ample divisor, and elementary Euler
characteristics.  It uses no Nikulin--Vinberg classification, ABSW theorem,
arithmetic volume, or polytope census.

## 1. Blowup coordinates and the sharp energy loss

Take \(m\) successive free blowups over a point outside
\(\operatorname {Supp}\Delta _0\).  If \(e_j\) is the total transform of the
exceptional divisor introduced at the \(j\)-th blowup, then

\[
 e_i.e_j=-\delta _{ij},\qquad
 K_{S_m}=\pi ^*K_{S_0}+\sum_{j=1}^m e_j.
\]

The physical prime chain is

\[
 R_i=e_i-e_{i+1}\quad(1\leq i<m),\qquad E=e_m.
\]

For

\[
 \Delta_m=\pi^*\Delta_0+\sum_{i=1}^{m-1}a_iR_i+bE,
 \qquad 0\leq a_i,b\leq1,
\]

direct collection in the orthogonal basis gives

\[
 -\!(K_{S_m}+\Delta_m)=\pi^*L_0-\sum_{j=1}^m x_je_j,
\]

with

\[
 x_1=1+a_1,\quad
 x_j=1-a_{j-1}+a_j,\quad
 x_m=1-a_{m-1}+b,
\]

and \(x_1=1+b\) when \(m=1\).  The intermediate coefficients telescope:

\[
 \sum_jx_j=m+b,\qquad
 L_m^2=L_0^2-\sum_jx_j^2.
\]

Cauchy--Schwarz therefore proves

\[
 L_m^2\leq L_0^2-\frac{(m+b)^2}{m}\leq L_0^2-m.
\]

The equality case is genuinely sharp.  For fixed \(b\), equality in the
Cauchy step occurs exactly at

\[
 a_i=\frac{ib}{m}\quad(1\leq i<m),
 \qquad x_1=\cdots=x_m=1+\frac bm.
\]

These coefficients lie in \([0,1]\).  Equality in the *outer* bound with
loss exactly \(m\) additionally forces \(b=0\), hence every \(a_i=0\).
Thus the zero tail boundary is the unique global square minimizer.  The
local formulas also have the asserted signs:

\[
 L_m.R_i=x_i-x_{i+1},\qquad L_m.E=x_m.
\]

If \(L_m\) is nef and big, then \(L_m^2>0\) on the smooth surface, so
\(m<L_0^2\).  This is a conditional intersection bound, not an existence
theorem for a weak log-Fano boundary and not an unconditional tail bound.

## 2. The coefficient-one gates

Write \(\Gamma\) for an effective divisor containing neither displayed
hostile component.  Distinct effective curves on a smooth surface have
nonnegative intersection, so all \(\Gamma\)-terms below are nonnegative.

For the minimal H1 incidence

\[
 B_1^2=B_2^2=-4,\qquad B_1.B_2=2,\qquad K.B_i=2,
\]

the two nef tests for \(-\!(K+\Delta)\), with
\(\Delta=a_1B_1+a_2B_2+\Gamma\), are precisely the two inequalities in the
owner report.  Their sum is

\[
 2(1-a_1)+2(1-a_2)+\Gamma.(B_1+B_2)\leq0.
\]

Every summand is nonnegative.  Hence \(a_1=a_2=1\) and the displayed
\(\Gamma\)-intersection vanishes.

For the minimal H2 incidence

\[
 B^2=-4,\qquad R^2=-2,\qquad B.R=2,\qquad K.B=2,\quad K.R=0,
\]

adding the two nef tests gives

\[
 2(1-a)+\Gamma.(B+R)\leq0.
\]

Thus \(a=1\) and \(\Gamma.(B+R)=0\); the \(R\)-test then reads
\(2(1-b)\leq0\), so \(b=1\).  This proof is continuous over the entire
coefficient cube, not merely a finite coefficient grid.

The H1 adjoint argument requires its two distinct transverse contacts.  The
H2 adjoint argument only needs the length-two intersection scheme.  Neither
argument licenses the same conclusion for the higher, parallel, or dotted
incidence branches.

## 3. Adjoint effectivity kills weak log-Fano positivity

For either minimal pair \(D_0\), the established adjoint-core argument gives
an effective divisor

\[
 A_0\in|K_{S_0}+D_0|,
\]

with \(A_0=0\) allowed.  Briefly, \(D_0\) is a connected Gorenstein
genus-one curve with trivial dualizing sheaf.  Rationality gives
\(H^1(S_0,K_{S_0})=0\), so the unit dualizing section lifts through the
canonical exact sequence.  Its divisor is \(A_0\).

After a free tail away from \(D_0\), coefficient rigidity forces any
hypothetical anti-log-canonically nef boundary to have the form
\(\Delta=D+\Gamma\).  Hence, as an \(\mathbb R\)-linear equivalence,

\[
 K_{S_m}+\Delta
 \sim_{\mathbb R}\pi^*A_0+\sum_{j=1}^m e_j+\Gamma.
\]

Every total-transform divisor \(e_j\) has an effective representative.  For
\(m\geq1\), the right side is a nonzero effective \(\mathbb R\)-divisor.  It
has positive intersection with every ample divisor \(H\).  If its negative
were nef, its intersection with \(H\) would be nonnegative, contradicting

\[
 -\bigl(\pi^*A_0+\sum e_j+\Gamma\bigr).H<0.
\]

For \(m=0\), the same ample-intersection argument says anti-log-canonical
nefness forces \(A_0+\Gamma=0\).  Then the anti-log-canonical class is zero
and is not big.  The no-go is therefore stronger than failure of a uniform
klt gap: for \(m\geq1\) nefness itself fails, and coefficient one on the core
already rules out klt.

## 4. Discrepancy and log-Chern checks

The change from total transforms to physical primes is

\[
 \sum_{j=1}^m e_j=\sum_{i=1}^{m-1}iR_i+mE.
\]

Consequently the relative coefficients of \(K+\Delta\) are
\(i+a_i\) on \(R_i\) and \(m+b\) on \(E\).  With zero tail boundary these
are the usual discrepancy ladder \(1,2,\ldots,m\) over the original smooth
point.  A formal crepant choice would require \(a_i=-i\) and \(b=-m\),
outside the effective coefficient cube.  Long free towers are therefore
compatible with a fixed positive discrepancy gap on the source; discrepancy
control without log-Fano positivity cannot bound their length.

For the log-Chern calculation, assume \(m\geq1\), let \(D_0\) be reduced
SNC, and keep the centers outside \(D_0\).  Omitting the exceptional chain
from the boundary gives

\[
 (K_{S_m}+D_m)^2=(K_{S_0}+D_0)^2-m,\qquad
 e(S_m\setminus D_m)=e(S_0\setminus D_0)+m.
\]

Including the full reduced chain \(T_m\) gives square loss \(m+3\), since
the all-one tail coefficients have exceptional coordinates
\((2,1,\ldots,1)\).  Removing the full exceptional fibre identifies the
complement with \(S_0\setminus(D_0\cup\{p\})\), so

\[
 (K_{S_m}+D_m+T_m)^2=(K_{S_0}+D_0)^2-(m+3),\qquad
 e(S_m\setminus(D_m+T_m))=e(S_0\setminus D_0)-1.
\]

Thus the residual \(3c_{2,\log}-c_{1,\log}^2\) grows with \(m\) in either
row.  This is only a direction check: the audit does not assert that every
pair under discussion satisfies the hypotheses of a particular log-BMY
theorem.

## 5. Repairs, scope, and reproducibility

The hostile pass made the following necessary repairs to the owner package:

1. all global H1/H2 language was narrowed to the two minimal incidences;
2. the equality sentence now distinguishes Cauchy equality from equality in
   the outer loss-\(m\) bound;
3. \(D_0\) and its strict transform \(D\) are distinguished, and the
   fractional-divisor equivalence is written \(\sim_{\mathbb R}\);
4. the BAB sentence now states the actual missing
   \(\varepsilon\)-lc log-Fano positivity premise;
5. the full-tail log-Chern row explicitly has \(m\geq1\); and
6. the frozen schema records that no higher, parallel, or dotted H1/H2 claim
   is made.

An independent test reconstructs the coefficient formulas without importing
the owner certifier, checks exact rational equality profiles, verifies both
rigidity decompositions, reconstructs the discrepancy and log-Chern changes,
checks the frozen scope, and pins the audited bytes.  The audited owner files
have SHA-256 hashes

    bf60458b0d1162fb7b1b17e6503102697c5f118119ab360ba7b67f321cef968c  results/rational_log_fano_tail_coefficient_no_go.md
    b3622c2122329d23239cc0b9e2bd6182dda080e1220fb52b90f32b0c60ef9e01  scripts/certify_rational_log_fano_tail_coefficient_no_go.py
    b64133616a6fcfd574b8da07823567f9e21ee2347db79dbea631d93ee8750baf  data/enumerations/rational_log_fano_tail_coefficient_no_go.json
    3a5de56425413c4a9d08510215bcdf68e998c6e8104765d011fea2cbae7e1c43  tests/test_rational_log_fano_tail_coefficient_no_go.py

Subject to those bytes and hypotheses, the verdict is **PASS**.  The theorem
is a sharp obstruction to one proposed AG route, not the missing absolute
dimension-bound theorem.
