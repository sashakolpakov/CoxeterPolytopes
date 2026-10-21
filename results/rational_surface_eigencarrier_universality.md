# Rational-surface eigen-carrier universality

## Result

There is an exact positive answer to the proposed “soft surface first”
principle for the arithmetic class.

> **Rational-surface eigen-carrier theorem.**  Every admissible totally-real
> arithmetic Coxeter datum of the lattice-form type specified below has a
> weak algebraic carrier on a blow-up of
> ​\(\mathbf P^2\).  Its hyperbolic space is one real eigenspace of an
> integral self-adjoint algebraic correspondence, and its Coxeter chamber is
> the semialgebraic section of that eigenspace cut out by integral divisor
> classes.

Here “admissible” has its standard quadratic-form meaning: at one real place
the form is Lorentzian, and at every other real place it is definite of the
opposite sign.  For the full restriction-of-scalars trace lattice, this is
also exactly the Hodge-index condition.  This last necessity statement is
about this trace-lattice mechanism; it does **not** say that an arbitrary,
unrelated surface realization of the same real chamber must remember all
Galois-conjugate forms.

This gives two things which were previously separate:

1. a universal existence theorem for a broad class of Coxeter polytopes; and
2. a sharp diagnosis of which extra geometric requirements are needed before
   surface classification can prove nonexistence of polytopes.

The theorem is constructive.  If the trace lattice has signature \((1,m)\),
the generic algorithm uses at most

\[
1+4m
\]

blow-ups.  If the field has degree \(d\) and the hyperbolic polytope has
dimension \(n\), then \(m=d(n+1)-1\), so the displayed bound is

\[
4d(n+1)-3.
\]

No K3 rank bound occurs.

## 1. Exact statement

Let \(F\) be a totally real field of degree \(d\), let \(M\) be a
torsion-free \(\mathcal O_F\)-lattice of rank \(n+1\), and let \(q\) be a
nondegenerate symmetric \(F\)-bilinear form, integral in the sense that
\(q(M,M)\subseteq\mathcal O_F\).  Fix a real embedding \(\sigma_0\) and
assume

\[
\operatorname{sig}(q_{\sigma_0})=(1,n),
\qquad
q_\sigma<0\quad(\sigma\ne\sigma_0).
\]

Choose a totally positive \(\beta\in F\) such that

\[
B(x,y)=\operatorname{Tr}_{F/\mathbf Q}(\beta q(x,y))
\]

is integral on the underlying \(\mathbf Z\)-module \(M_{\mathbf Z}\).  For an
\(\mathcal O_F\)-valued form, \(\beta=1\) already gives integrality; allowing
\(\beta\) also covers the useful codifferent normalizations.

Suppose roots \(r_1,\ldots,r_s\in M\) define a finite-sided Coxeter chamber
\(P\subset\mathbb H(q_{\sigma_0})\), with the field reflections integral on
\(M\).  Then there are

- an integer \(D>0\);
- an integer \(b\le 4d(n+1)-3\);
- a blow-up
  \[
  X=\operatorname{Bl}_{p_1,\ldots,p_b}(\mathbf P^2);
  \]
- an integral lattice injection
  \[
  j:(M_{\mathbf Z},D^2B)\hookrightarrow
  \operatorname{NS}(X)\cong I_{1,b};
  \]
- a transpose-invariant integral divisor-product correspondence
  \(\Gamma\in\operatorname{CH}^2(X\times X)\); and
- a rational divisor-product Chow projector \(\pi_U\), whose image is a
  rank-\(d(n+1)\) Tate-motive summand carrying a unital \(F\)-action; and
- integral divisor classes \(D_i=j(r_i)\);

such that one real eigenspace \(E_{\sigma_0}\) of \(\Gamma_*\) has signature
\((1,n)\) and

\[
P\cong
\left\{
[x]\in\mathbf P(E_{\sigma_0}):
x^2>0,\ (x,D_i)\ge0\text{ for }1\le i\le s
\right\}_{\text{chosen component}}.
\]

This is an isometry up to one positive common scale.  Consequently all
Coxeter angles, ultraparallel distances, incidence data, finite-volume
status, and compactness are preserved exactly.

The word *weak* is essential: the \(D_i\) are integral divisor classes, but
need not be effective, irreducible, or roots of integral reflections of the
full Néron--Severi lattice.  The algebraic object is therefore the tuple

\[
(X,\Gamma;D_1,\ldots,D_s),
\]

and the polytope is recovered from it by a finite system of quadratic and
linear inequalities.  It is literally a semialgebraic divisor-eigenslice.

## 2. Why admissibility is exactly the Hodge-index gate

Over \(\mathbf R\), restriction of scalars gives an orthogonal decomposition

\[
(M_{\mathbf Z}\otimes\mathbf R,B)
\cong
\bigoplus_{\sigma:F\hookrightarrow\mathbf R}
(M\otimes_{F,\sigma}\mathbf R,\sigma(\beta)q_\sigma).
\]

Because \(\beta\) is totally positive, the positive index of \(B\) is the
sum of the positive indices of all conjugate forms.  The distinguished place
already contributes one.  Hodge index says that the Néron--Severi space of a
projective surface has positive index exactly one.  Therefore a trace-form
divisor carrier can exist only if every other conjugate form is negative
definite.

Conversely, that condition makes \(B\) Lorentzian, of signature

\[
(1,d(n+1)-1),
\]

and the explicit embedding below produces the surface carrier.  Thus, for a
carrier required to contain the **full trace form isometrically up to positive
scale**, the one-indefinite-place condition is necessary and sufficient at
the Hodge-lattice level.  Its orthogonal complement in the Néron--Severi
space is then automatically negative definite.

For the opposite sign convention—common in Vinberg computations—the stored
form has signature \((n,1)\) at \(\sigma_0\), is positive definite elsewhere,
and is paired with a totally negative trace scalar.  Replacing both by their
negatives gives the displayed hypotheses and leaves the trace form unchanged.
This is only sign bookkeeping; \(\beta\) in the displayed theorem remains
totally positive.

## 3. Constructive embedding in a rational-surface divisor lattice

Let \(G\) be any integral Gram matrix of signature \((1,m)\).  Rational
Gram--Schmidt gives

\[
C^tGC=\operatorname{diag}(a_0,-a_1,\ldots,-a_m),
\qquad a_i\in\mathbf Q_{>0}.
\]

The positive coefficient is a difference of two rational squares:

\[
a_0=\left(\frac{a_0+1}{2}\right)^2-
    \left(\frac{a_0-1}{2}\right)^2.
\]

Write \(a_i=p_i/q_i\).  Lagrange's four-square theorem applied to
\(p_iq_i\) gives

\[
a_i=\sum_{k=1}^4\left(\frac{u_{ik}}{q_i}\right)^2.
\]

Using disjoint negative coordinate blocks for these representations embeds
\(G\) rationally in \(I_{1,1+4m}\).  Clearing a common denominator \(D\)
produces an integral matrix \(E\) satisfying

\[
E^tI_{1,1+4m}E=D^2G.
\]

Now choose \(1+4m\) distinct points and blow up \(\mathbf P^2\).  In the
basis consisting of the pullback \(H\) and the exceptional curves \(E_i\),

\[
\operatorname{NS}(X)=
\mathbf ZH\oplus\bigoplus_i\mathbf ZE_i,
\qquad
H^2=1,\quad E_i^2=-1,
\]

so the matrix above is an injection by genuine integral divisor classes.

The implementation handles isotropic starting bases as well as diagonal
ones.  Its frozen calibration starts from \(U\oplus\langle-2\rangle\), not
from an already diagonal form.

## 4. The eigenspace is induced by an algebraic correspondence

Choose a nonzero primitive element \(\alpha\in F\) (take \(\alpha=1\) when
\(F=\mathbf Q\)).  Multiplication by \(\alpha\) on \(M_{\mathbf Q}\) is
self-adjoint for \(B\):

\[
B(\alpha x,y)=B(x,\alpha y).
\]

On the embedded rational subspace \(U=j(M_{\mathbf Q})\), transport this
operator and extend it by zero over \(U^\perp\).  No conjugate of a nonzero
algebraic number is zero.  After multiplying by an integer, the extension
\(A\) preserves \(\operatorname{NS}(X)\).  Its nonzero real eigenspaces are
exactly the conjugate summands above; the
\(\sigma_0(\alpha)\)-eigenspace is the Lorentzian one.

This integral operator is algebraic in a completely explicit sense.  Let
\(e_0=H,e_i=E_i\), and let \(e_j^*\) be the intersection-dual basis.  Then

\[
\Gamma_A=\sum_j e_j^*\times A(e_j)
\in\operatorname{CH}^2(X\times X)
\]

acts on \(\operatorname{NS}(X)\) as \(A\).  In the product basis its
coefficient matrix is symmetric precisely because \(A\) is self-adjoint, so
\(\Gamma_A={}^t\Gamma_A\) already as the displayed integral divisor-product
cycle.  The integrality uses that \(H,E_1,\ldots,E_b\) is a unimodular basis,
so its intersection-dual basis \(H,-E_1,\ldots,-E_b\) is integral.  Thus no
standard conjecture and no claim about surface automorphisms is being used.

There is a stronger rational statement behind the chosen integral operator.
Let

\[
P_U=j(j^tIj)^{-1}j^tI
\]

be the rational orthogonal projector onto \(U\), and extend multiplication by
every \(a\in F\) by zero on \(U^\perp\).  These extensions give a unital
\(F\)-action on the image of \(P_U\), with \(P_U\) as its unit.  This is an
action by actual Chow correspondences, not merely by endomorphisms of
cohomology.  Indeed, if \(G=I_{1,b}\) is the intersection matrix and \(C_A\)
is the product-basis coefficient matrix of \(\Gamma_A\), then

\[
C_A=G^{-1}A^t,
\qquad
C_A G C_B=C_{BA}.
\]

The second identity is the literal composition law obtained by intersecting
the two middle divisor factors.  Consequently \(P_U^2=P_U\) proves
\(\Gamma_{P_U}\circ\Gamma_{P_U}=\Gamma_{P_U}\) in
\(\operatorname{CH}^2(X\times X)_{\mathbf Q}\).  Since the image lies in the
divisor part of the motive of a rational surface, it is a rank-\(d(n+1)\)
Tate summand.  The full \(F\)-action and projector are rational; after choosing
a primitive field generator, a positive integer multiple of its
correspondence is the integral \(\Gamma\) in the theorem.

For \(x\) in the distinguished real eigenspace, the divisor pairing with a
transported root is

\[
(x,j(r_i))=D^2\sigma_0(\beta)q_{\sigma_0}(x,r_i).
\]

The positive factor proves that the divisor half-spaces reproduce the
original oriented chamber.

## 5. Compact exact calibration: Bugaenko on a blow-up of \(\mathbf P^2\)

The general bound is intentionally crude.  Bugaenko's compact 11-facet
polytope in \(\mathbb H^8\) has trace lattice

\[
N=\langle2\rangle\oplus\langle-2\rangle
  \oplus E_8(-1)^2.
\]

There is an explicit \(20\times18\) integer matrix \(E\) with

\[
E^tI_{1,19}E=4G_N.
\]

The construction uses

- \((3,1)\) for the vector of square \(8\);
- \((2,2)\) in two negative coordinates for the vector of square \(-8\);
  and
- two disjoint scaled coordinate realizations \(R^tR=4E_8\).

Therefore take the entirely explicit rational surface

\[
X=\operatorname{Bl}_{p_0,\ldots,p_{18}}(\mathbf P^2_{\mathbf Q}),
\qquad p_k=[1:k:k^2].
\]

All points are distinct and rational, and
\(\operatorname{NS}(X)=I_{1,19}\).  The eleven Bugaenko roots transport to
eleven integral divisor classes.  Their full pairing matrix is four times
the pinned trace pairing.  More importantly, the certificate independently
projects all eleven divisor normals to the \(4\phi\)-eigenspace and verifies

\[
  ((D_i)_{4\phi},(D_j)_{4\phi})
  =4\,\sigma_0(\alpha)\,h(r_i,r_j),
  \qquad \alpha=-\phi/\sqrt5.
\]

It also transports twice the archived strict Bugaenko interior point and
checks positive square and all eleven strict oriented inequalities.  Thus the
basis conversion, metric, and chosen half-spaces are checked independently,
not inferred from an unprojected trace Gram matrix.

Multiplication by \(\phi=(1+\sqrt5)/2\) extends, after multiplier four, to an
integral self-adjoint \(20\times20\) operator with only 40 nonzero entries
and largest absolute entry 9.  The corresponding integral cycle has 40
divisor-product terms and is literally invariant under transposition.  Its
characteristic polynomial is

\[
t^2(t^2-4t-16)^9,
\]

and its real eigenspaces are

| eigenvalue | dimension | signature |
|---|---:|---:|
| \(4\phi\) | 9 | \((1,8)\) |
| \(4\phi'\) | 9 | \((0,9)\) |
| \(0\) | 2 | \((0,2)\) |

For this calibration the rational trace-space projector \(\pi\) and the
zero-complement golden action \(A\) are also frozen at the cycle level.  Their
coefficient matrices satisfy, under exact Chow composition,

\[
\pi^2=\pi,
\qquad A\pi=\pi A=A,
\qquad A^2=A+\pi.
\]

Thus the rank-18 divisor Tate summand carries an actual
\(\mathbf Q(\sqrt5)\)-action; the integral 40-term correspondence is \(4A\).

The eleven divisors and their golden translates span the full rank-18 trace
subspace.  The inherited exact Vinberg replay gives 41 finite vertices, no
cusps, no hyperideal rays, and a complete wall list.  Hence the
\(4\phi\)-eigenslice is exactly Bugaenko's compact \(\mathbb H^8\) chamber.

This is a new realization in the present ledger, not a new abstract Coxeter
polytope.  The previously frozen K3 eigen-nef realization is geometrically
stronger for this particular chamber; the point here is that the
rational-surface construction is universal and has no Picard-rank ceiling.

## 6. What this proves—and what must be strengthened for primal no-go results

The theorem proves that a bare divisor-class carrier is too soft to constrain
admissible arithmetic Coxeter polytopes: rational surfaces supply one
automatically, in every rank.  Therefore classification of surfaces alone
cannot lower the Coxeter nonexistence bound if the definition asks only for
divisor classes and a correspondence eigenspace.

The useful stronger carrier notions are now sharply identified:

1. every facet class is represented by an effective irreducible negative
   curve;
2. the Coxeter chamber is the complete nef or movable eigenslice;
3. the reflections or diagram symmetries arise from actual birational maps;
4. the correspondence belongs to a restricted natural algebra (for example,
   generated by a prescribed fibration, involution, or moduli construction);
5. the carrier is rank-minimal, full-Picard, K3, Enriques, or hyperkähler.

Those requirements can have classification content.  The universal weak
carrier theorem separates them cleanly from the part supplied for free by
Hodge index and rational blow-ups.

## Reproducible files

- theorem and compact certificate generator:
  `scripts/certify_rational_surface_eigencarrier_universality.py`;
- reusable exact construction:
  `src/coxeter_polytopes/geometric/rational_surface_carrier.py`;
- frozen certificate:
  `data/enumerations/rational_surface_eigencarrier_universality.json`;
- tests:
  `tests/test_rational_surface_eigencarrier_universality.py`.
