# Hostile scope audit: the power-coupled K3 pentagon obstruction

**Verdict (2026-08-11): PASS, scoped.**  The argument is an exact
conditional obstruction on a visibly declared sublocus.  The exact
K3/maximal package is retained: rank-six primitive full trace, minimal
complex projective K3, exact compact right-pentagon eigen-nef chamber,
integral complete/saturated restricted reflections, and reflection
maximality in the declared order host.  The power relation and algebraic
integrality are extra hypotheses, displayed separately.  The report does
not convert them into consequences of maximality.

## 1. Gram-closure audit

For cyclic facets, normalize only the squared correlations

\[
 {h(e_i,e_j)^2\over h(e_i,e_i)h(e_j,e_j)}.
\]

They remain elements of \(k\); no square root of a facet norm is assumed to
lie in \(k\).  With the report's labels, the five principal four-by-four
Gram determinants are

\[
\begin{aligned}
 CE-C-D-E+1,&\qquad -A+BE-B-E+1,\\
 BD-B-C-D+1,&\qquad AD-A-D-E+1,\\
 AC-A-B-C+1.&
\end{aligned}
\]

The rank-three condition makes each zero.  Solving gives

\[
 B=(A-1)(C-1),\qquad
 D={A(C-1)\over AC-A-C},\qquad
 E={C(A-1)\over AC-A-C}.
\]

Exact rational replay verifies all five identities, not only the three
printed in the proof.  Since every nonadjacent pair in a compact right
pentagon is ultraparallel, \(B,D,E>1\), so the selected denominator branch
is positive.

On the additional locus \(C=A^q\),

\[
 u=A^q-A^{q-1}-1,\qquad AC-A-C=Au,
\]

and direct substitution gives

\[
 E=1+u^{-1},\qquad D=1+A^{q-1}u^{-1}.
\]

There is no sufficiency or existence claim hidden in these necessary Gram
identities.

## 2. Integrality audit

The pair cyclic product attached to the \(E\)-pair is \(4E\).  Although
individual Cartan entries need not belong to the returned order, the pair
product does: it differs by an integer from the module trace of the product
of the two integral reflections on the finite projective returned
order-module.  Thus

\[
                         v=4E-4={4\over u}\in\mathcal O.
\]

The separate ansatz \(A\in\mathcal O\) makes
\(u=A^q-A^{q-1}-1\) integral.  Therefore the identity

\[
                                  u\,v=4
\]

really is an identity between algebraic integers.  The proof correctly
warns that \(4A\in\mathcal O\) would not imply \(A\in\mathcal O\); this is
why the extra integrality premise cannot be deleted.

## 3. Signed-place and norm audit

At the distinguished Lorentzian place, ultraparallelity gives \(A_0>1\)
and \(E_0>1\), hence \(u_0,v_0>0\).  At the positive-definite conjugate
place, strict Cauchy--Schwarz gives \(0<A_1<1\).  For every \(q\ge2\),

\[
 u_1=-\bigl(1+A_1^{q-1}(1-A_1)\bigr)\in(-2,-1),
 \qquad v_1\in(-4,-2).
\]

Consequently \(N(v)=-n\) with \(n>0\).  Algebraic integrality and
\(uv=4\) give

\[
                     n\mid16,\qquad
                     u_0<{16\over n}\le16.
\]

No inequality direction reverses: \(|v_1|<4\) gives
\(n=v_0|v_1|<4v_0\), hence \(u_0=4/v_0<16/n\).  It follows that

\[
 A_0^{q-1}(A_0-1)<17,\qquad
 A_0<{1+\sqrt{69}\over2}.
\]

This bounds the field-generating cyclic product \(4A_0\).  It is a cyclic
height bound, not a claimed bound for a K3 correspondence intersection.

## 4. Quadratic enumeration audit

The rank-six full trace of a rank-three field module forces
\([k:\mathbf Q]=2\).  Since \(A\) is an algebraic integer generator and
\(0<A_1<1\), its trace and norm obey

\[
                         2\le T\le5,
                  \qquad 1\le N\le4.
\]

Exact enumeration of this \(4\times4\) box, with a positive nonsquare
discriminant and the two root intervals, leaves

\[
 (T,N,D)=(3,1,5),(4,1,12),(4,2,8),(5,2,17),(5,3,13).
\]

Their squarefree fields are precisely

\[
 \mathbf Q(\sqrt2),\ \mathbf Q(\sqrt3),\ \mathbf Q(\sqrt5),
 \ \mathbf Q(\sqrt{13}),\ \mathbf Q(\sqrt{17}).
\]

This is a necessary five-field list.  The theorem does not claim that all
five rows, or even one row, is realized by an exact K3 pentagon satisfying
every ansatz condition.

## 5. Counter-scope

For independent \(A,C\), the closure coordinate is

\[
                       E={C(A-1)\over AC-A-C}.
\]

The moving denominator has no certified fixed rational divisor.  The audit
therefore confirms only the power-coupled obstruction.  It does not solve
the general \(n=2\) automatic-height problem, does not obtain a mixed
correspondence-degree bound from reflection maximality, and does not rule
out an exact nonorthogonally coupled infinite-field family outside this
locus.  No reflective-lattice classification enters the proof or replay.
