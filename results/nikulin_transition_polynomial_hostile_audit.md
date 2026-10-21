# An exact hostile audit of Nikulin's transition polynomial

**Status: the arXiv-v4 certificate is false as printed; the 2011 journal
version silently contains the correct one-exponent repair.  This audit does
not import Nikulin's global theorem and does not claim to disprove it.**

The source under audit is V. V. Nikulin, *The transition constant for
arithmetic hyperbolic reflection groups*.  In
[arXiv:0910.5217v4](https://arxiv.org/abs/0910.5217), pages 15--16,
equation (24) prints

\[
 P_2(x)=(x+2)^3(x+1)^4x^4(x-1)^3(x^2+x-1)^2
\]

and assigns it the sup norm

\[
 \delta=0.0905969664=\frac{884736}{9765625}
 \qquad (-2\leq x\leq1).
\]

That assertion is false.  The peer-reviewed version, Izvestiya: Mathematics
75:5 (2011), 971--1005,
[DOI 10.1070/IM2011v075n05ABEH002561](https://doi.org/10.1070/IM2011v075n05ABEH002561),
prints exponent `3` in equation (4.19), page 986.  We independently verify
that this is exactly the repair required by the calculation.

The exact replay is
[`scripts/independent_nikulin_transition_polynomial_audit.py`](../scripts/independent_nikulin_transition_polynomial_audit.py),
with frozen output
[`data/enumerations/nikulin_transition_polynomial_audit.json`](../data/enumerations/nikulin_transition_polynomial_audit.json).

## 1. A rational counterexample to the arXiv formula

There is no need for floating point arithmetic.  At the rational point
\(x=3/4\),

\[
 |P_2(3/4)|=\frac{6471355275}{68719476736}
              >\frac{884736}{9765625}.
\]

Thus the printed value cannot be the maximum.

The exact maximum is also elementary.  The polynomial is invariant under
\(x\mapsto-1-x\).  Put \(y=x(x+1)\); as \(x\) runs through \([-2,1]\),
\(y\) runs through \([-1/4,2]\), and

\[
 |P_2(x)|=F_2(y)=y^4(y-1)^2(2-y)^3.
\]

Its derivative factors as

\[
 F_2'(y)=-y^3(y-2)^2(y-1)(9y^2-19y+8).
\]

Checking the endpoints and these five critical values gives the unique
maximal \(y\)-value

\[
 y_+=\frac{19+\sqrt{73}}{18}
\]

and

\[
 \|P_2\|_{[-2,1]}
 =\frac{42029518+2327386\sqrt{73}}{387420489}
 =0.1598126956594688196644863666\ldots.
\]

The exact minimum polynomial of this maximum is

\[
 387420489T^2-84059036T+3538944.
\]

## 2. What the printed polynomial actually proves

For this subcase Nikulin applies the inequality

\[
 [K:\mathbb Q]\leq
 1+\frac{\log |P(t)|}{-\log\delta},\qquad t=14.
\]

For the printed \(P_2\),

\[
 P_2(14)=764468783534776320000.
\]

Let \(\delta_2=\|P_2\|_{[-2,1]}\).  The two exact algebraic sign checks

\[
 P_2(14)\delta_2^{26}>1,
 \qquad
 P_2(14)\delta_2^{27}<1
\]

show that the real right-hand side lies strictly between `27` and `28`:

\[
 1+\frac{\log P_2(14)}{-\log\delta_2}
 =27.2225709066620200011\ldots.
\]

Consequently the arXiv-v4 equation, taken literally, gives only
\([K:\mathbb Q]\leq27\), not `23`.

## 3. The independently verified repair

Replace only the last exponent by `3`:

\[
 P_3(x)=(x+2)^3(x+1)^4x^4(x-1)^3(x^2+x-1)^3.
\]

This is not an ad hoc guess.  If

\[
 Q(z)=z^3(z-1)^4(z-2)^4(z-3)^3(z^2-3z+1)^3
\]

is the earlier equation (13)/(4.8), then the exact identity

\[
 P_3(x)=Q(x+2)
\]

explains the paper's phrase “similar”/“analogous to” that polynomial.  In the
same coordinate \(y=x(x+1)\), the derivative has only the new critical
factors \((2y-1)(5y-8)\).  Direct comparison of all critical values gives

\[
 \|P_3\|_{[-2,1]}
 =|P_3(y=8/5)|
 =\frac{884736}{9765625}
 =0.0905969664.
\]

Moreover

\[
 P_3(14)=159773975758768250880000,
\]

and the exact rational inequalities

\[
 P_3(14)\delta^{22}>1,
 \qquad
 P_3(14)\delta^{23}<1
\]

place the real bound strictly between `23` and `24`:

\[
 23.2493142903636456507\ldots.
\]

Thus the repaired certificate proves the integer bound `23`.  This is exactly
the exponent printed in journal equation (4.19), so this particular arXiv
fault does **not** create a gap in the final journal theorem.

## 4. Propagation and precise scope

Had the exponent-2 formula remained unrepaired, its `s=k=3` value would feed
the arXiv-v4 Theorem 11, Propositions 12--13, Theorem 4, and the applications
that quote `N(14)<=25`.  The journal chain instead uses corrected equation
(4.19), Theorem 4.6, Propositions 4.1--4.2, and closes Theorem 3.2 on page 998.

Therefore the accurate verdict is:

- arXiv v4 contains a literal, exactly refutable polynomial certificate;
- the journal silently repairs it;
- our calculation verifies that repair independently;
- this audit neither proves nor disproves the remaining global transition
  theorem.

## 5. Why this is not our dimension theorem

Nikulin is being used here as a hostile regression corpus.  Our AG/dessin
route must separately certify the following arrows:

1. a geometric probe of every chamber, including cusps;
2. a bounded, sheet-separated marked carrier for that probe;
3. a finite numerical Gram image, not merely a finite-type surface family;
4. recovery of the complete chamber and its cusp data;
5. a bound on dimension or field degree produced by those carrier invariants.

The useful translated insight is the **small probe principle**: try to force a
bounded facet/root cluster whose complete Gram data already controls the
field.  In AG language this should be a bounded marked boundary cluster or a
bounded correspondence packet with finite Gram image.  The dangerous leap is
to assume that such a local packet controls every compact and cusped chamber.
That arrow remains an object of independent proof and counterexample search;
no Nikulin/Vinberg classification result is allowed to supply it.

