# Universal weak Chow-action carriers and the leverage no-go

**Status (2026-08-11): unconditional in the marked rational-Chow category.**

## 1. Result

There is an exact answer to the missing outer quantifier once the word
*carrier* is allowed to mean a divisor Tate motive with rational algebraic
correspondences.

> **Theorem 1 (universal weak divisor/Chow-action carrier).**  Let
> \(P\subset\mathbb H^n\) be a finite-sided arithmetic Coxeter chamber of
> lattice form.  Let \(k\) be its totally real Vinberg field, let \((M,q)\)
> be a \(W(P)\)-stable integral \(\mathcal O_k\)-lattice of rank \(n+1\),
> and let \(r_1,\ldots,r_f\in M\) be oriented facet normals.  At the
> distinguished embedding \(\sigma _0\), assume
>
> \[
> \operatorname {sig}(q_{\sigma _0})=(1,n),\qquad
> q_\sigma<0\quad(\sigma\ne\sigma _0).
> \]
>
> Put \(d=[k:\mathbf Q]\) and \(R=d(n+1)\).  Then one can construct
>
> 1. a rational surface \(X=\operatorname {Bl}_{b}(\mathbf P^2_{\mathbf Q})\)
>    with \(b\le4R-3\);
> 2. a rational Chow projector
>    \(\pi_U\in\operatorname {CH}^2(X\times X)_{\mathbf Q}\), whose image
>    is a rank-\(R\) divisor Tate-motive summand;
> 3. a faithful action
>
>    \[
>       W(P)\longrightarrow
>       \operatorname {Aut}_{\operatorname {ChowMot}_{\mathbf Q}}
>       (X,\pi_U)
>    \]
>
>    by divisor-product Chow correspondences, commuting with a unital
>    \(k\)-action on \((X,\pi_U)\); and
> 4. integral divisor classes \(D_i\in\operatorname {NS}(X)\),
>
> such that the distinguished real \(k\)-eigenspace \(E_{\sigma _0}\) is
> Lorentzian and
>
> \[
> P\cong
> \left\{[x]\in\mathbf P(E_{\sigma _0}):
> x^2>0,\ (x,D_i)\ge0\ (1\le i\le f)\right\}_{\rm chosen\ component}.
> \tag{1.1}
> \]
>
> The side correspondence attached to \(r_i\) is a self-transpose
> involution of the summand, induces the original reflection \(s_{r_i}\),
> and the complete Coxeter relations hold under literal composition in the
> Chow ring.  Thus the construction returns the Gram metric, all oriented
> half-spaces, and the full side group, not only the abstract Coxeter matrix.

Finite volume and reflection maximality are not needed to construct this
carrier.  They remain properties of the returned chamber and group.

The word **weak** is indispensable.  The \(D_i\) need not be effective or
irreducible; (1.1) need not be the nef, movable, or Mori cone of \(X\); and
the side correspondences need not be graphs of automorphisms or birational
self-maps of \(X\).

### Uniform published arithmetic corollary

Vinberg's published \(n<30\) cutoff for arithmetic reflection groups and the
published global degree bound \([k:\mathbf Q]\le25\) imply

\[
R=d(n+1)\le25\cdot30=750,
\qquad b\le4R-3\le2997.                    \tag{1.2}
\]

Consequently every finite-volume arithmetic Coxeter datum has an
unconditional weak rational-Chow carrier on a blow-up at no more than 2,997
points.  By adding unused exceptional directions, all of them can be put on
one fixed lattice type \(I_{1,2997}=\operatorname {NS}(\operatorname
{Bl}_{2997}\mathbf P^2)\).  This numerical corollary **imports** the
classical dimension and field-degree bounds.  It is not a new proof of
either bound.

The base embedding and field-action theorem is frozen in
[`rational_surface_eigencarrier_universality.md`](rational_surface_eigencarrier_universality.md).
The new content here is the simultaneous exact action of every side
reflection and the sharp leverage audit below.

## 2. Proof

### 2.1 The invariant trace lattice

Arithmeticity supplies an admissible quadratic space over a totally real
field.  If a chosen arithmetic lattice is only commensurable with an
integral orthogonal group, take a normal finite-index subgroup preserving an
\(\mathcal O_k\)-lattice and sum the finitely many translates of that
lattice.  The result is a full \(W(P)\)-stable lattice \(M\).  Scaling the
form and the finitely many root-line representatives clears the remaining
denominators without changing a reflection or a wall.

Choose a totally positive \(\beta\in k\) so that

\[
B(x,y)=\operatorname {Tr}_{k/\mathbf Q}(\beta q(x,y))              \tag{2.1}
\]

is integral on the underlying rank-\(R\) group \(M_{\mathbf Z}\).  Restriction
of scalars gives

\[
(M_{\mathbf Z}\otimes\mathbf R,B)
 =\bigoplus_{\sigma:k\hookrightarrow\mathbf R}(M_\sigma,
 \sigma(\beta)q_\sigma),
\]

so \(B\) has signature \((1,R-1)\).  Every \(k\)-linear reflection in
\(W(P)\) is a \(B\)-isometry and commutes with multiplication by \(k\).

### 2.2 Embedding in a rational surface

For any integral Gram matrix \(G\) of signature \((1,R-1)\), rational
Gram--Schmidt writes

\[
G\sim_{\mathbf Q}\langle a_0,-a_1,\ldots,-a_{R-1}\rangle,
\qquad a_i>0.
\]

Write \(a_0\) as a difference of two rational squares and each \(a_i\),
\(i>0\), as a sum of four rational squares.  Lagrange's theorem and clearing
one denominator give an integral injection \(j\) satisfying

\[
j^t I_{1,b}j=D^2G,
\qquad b\le1+4(R-1)=4R-3.                 \tag{2.2}
\]

The lattice \(I_{1,b}\) is the N\'eron--Severi lattice of
\(X=\operatorname {Bl}_b\mathbf P^2\).  Put

\[
U=j(M_{\mathbf Q}),\qquad
\Pi=j(j^tI_{1,b}j)^{-1}j^tI_{1,b}.        \tag{2.3}
\]

Then \(\Pi^2=\Pi\), \(\Pi\) is intersection-self-adjoint, and its image is
\(U\).  Equation (2.2) gives the integral wall divisors \(D_i=j(r_i)\).

### 2.3 The full reflection action is algebraic

For \(g\in W(P)\), define an endomorphism of
\(\operatorname {NS}(X)_{\mathbf Q}\) by

\[
\widehat g=jgj^{-1}\Pi.                                  \tag{2.4}
\]

It acts as \(g\) on \(U\) and as zero on \(U^\perp\).  Hence

\[
\widehat g\widehat h=\widehat{gh},\qquad
\widehat 1=\Pi.                                           \tag{2.5}
\]

In particular, \(\widehat s_i^2=\Pi\), and every finite Coxeter relation
\((s_is_j)^{m_{ij}}=1\) becomes

\[
(\widehat s_i\widehat s_j)^{m_{ij}}=\Pi.                  \tag{2.6}
\]

Let \(e_0=H,e_1=E_1,\ldots,e_b=E_b\) be the unimodular divisor basis and
let \(e_a^*\) be its intersection-dual basis.  Any rational endomorphism
\(A\) of \(\operatorname {NS}(X)_{\mathbf Q}\) is induced by the literal
divisor-product cycle

\[
Z_A=\sum_a e_a^*\times A(e_a)
       \in\operatorname {CH}^2(X\times X)_{\mathbf Q}.    \tag{2.7}
\]

If \(J=I_{1,b}\), its product-basis coefficient matrix is
\(C_A=J^{-1}A^t\).  Intersecting the two middle divisor factors gives

\[
C_AJC_B=C_{BA}.                                           \tag{2.8}
\]

Thus (2.5)--(2.6) are identities of Chow correspondences, not only
cohomological identities.  The projector \(Z_\Pi\) defines the Chow-motive
summand.  A reflection is self-adjoint and involutive, so \(Z_{\widehat
s_i}\) is self-transpose and is an automorphism of \((X,Z_\Pi)\).  The
faithfulness of the native Lorentzian representation proves faithfulness on
the summand.

Multiplication by every \(a\in k\) is treated by the same construction.
It commutes with (2.4), so the \(k\)-action and \(W(P)\)-action commute.

### 2.4 Exact chamber return

The real \(k\)-eigenspace decomposition of \(U_{\mathbf R}\) is the
restriction-of-scalars decomposition.  For \(x\in E_{\sigma _0}\),

\[
(x,D_i)_X=D^2\sigma _0(\beta)
q_{\sigma _0}(x,r_i).                                    \tag{2.9}
\]

The multiplier is positive and independent of \(i\).  Therefore every
oriented half-space, angle, ultraparallel distance, and incidence relation
is unchanged.  This proves (1.1).

## 3. Why this theorem cannot supply the ABSW bounds

The obstruction is structural, not a missing computation.

### Proposition 2 (rank/field leverage no-go for the weak category)

Existence of a carrier satisfying Theorem 1, without an independently
bounded carrier complexity, implies neither a dimension bound nor a finite
set of Vinberg fields.

#### Proof

An exact full-trace carrier only gives

\[
d(n+1)=R\le\rho(X).                                       \tag{3.1}
\]

But rational blow-ups have \(\rho(\operatorname {Bl}_b\mathbf P^2)=b+1\)
for arbitrary \(b\), and the construction deliberately chooses \(b\) as a
function of \(R\).  Hence (3.1) is not a ceiling.

More strongly, on the divisor Tate part of a rational surface, product
correspondences realize the complete matrix algebra

\[
\operatorname {End}_{\mathbf Q}(\operatorname {NS}(X)_{\mathbf Q}).
\]

Thus the existence of a field action or reflection action imposes no
motivic scarcity: all rational matrices occur algebraically.  Even on one
fixed sufficiently large blow-up, infinitely many number fields of a fixed
degree act on rational Tate summands, and wall coordinates and correspondence
denominators are unbounded.  Bounded degree alone is not a finite field set.

The bound \(b\le2997\) in (1.2) does not reverse this logic: it was obtained
by inserting the already-published \(n<30\) and \(d\le25\) bounds.  Using it
to rederive either input would be circular.  \(\square\)

This proposition is limited and precise.  It says that **weak Chow-carrier
existence itself** has no pruning power.  It does not say that finite-volume
geometry, maximality, arithmetic volume, or a stronger geometric carrier
cannot yield bounds.

## 4. Two stronger-action warnings

### 4.1 Actual actions alone are also cheap

Let \(T_M\) be the split torus with cocharacter lattice \(M_{\mathbf Z}\).
Every \(g\in W(P)\subset GL(M_{\mathbf Z})\) acts on \(T_M\) by a monomial
algebraic automorphism.  On the standard compactification
\((\mathbf P^1)^R\), it becomes a monomial Cremona transformation.  The
trace form and roots recover the chamber in the real cocharacter space.
Thus merely asking for a faithful algebraic-torus or Cremona action is even
weaker than Theorem 1 and cannot create a rank bound.

### 4.2 Effective wall reflections cannot all be regular automorphisms

Let \(C\) be a nonzero effective divisor on a smooth projective surface with
\(C^2<0\).  The reflection in its class sends \([C]\) to \(-[C]\).  A
regular automorphism preserves the effective cone, while \(-[C]\) is not
effective because it has negative intersection with every ample divisor.
Therefore no regular automorphism induces the reflection in an effective
negative wall class.

Any strong mixed category must consequently distinguish at least

- effective divisors supplying cone walls;
- internal divisor cuts which are not Mori walls; and
- automorphism, birational, or correspondence operators supplying group
  generators.

The characteristic-two Borcherds surface makes exactly this distinction:
42 curves supply actual walls, while 168 further generators are Cremona
actions.

## 5. Audit of the stronger positive-characteristic route

No universal effective/nef or mixed positive-characteristic theorem is
currently closed.

1. A complex K3 carrier has \(\rho\le20\); a K3 surface in arbitrary
   characteristic has \(\rho\le22\).  Therefore an injective native block
   gives \(n+1\le22\), while a full trace/order block gives
   \(d(n+1)\le22\).
2. At rank 22, a supersingular K3 N\'eron--Severi lattice is
   \(p\)-elementary with discriminant \(-p^{2\sigma}\).  An arbitrary
   full-rank arithmetic trace lattice cannot simply be declared to be such a
   Picard lattice.
3. Borcherds' absolutely maximal \(H^{21}\) chamber is a positive endpoint:
   its characteristic-two Artin-invariant-one K3 realization mixes 42
   curves with 168 Cremona actions.  This proves one carrier, not essential
   surjectivity for every maximal arithmetic chamber.
4. Vinberg's \(X_3\) chamber is another positive weak mixed-divisor example,
   but its norm-\((-6)\) cuts are not all nef-root walls and it is not known
   to settle the maximal-only quantifier.

Accordingly the universal statement proved here is exactly the rational-Chow
one.  A theorem capable of an AG-derived ABSW bound still needs a universal
promotion to a category with an independently bounded numerical lattice and
an exhaustive effective/nef/MBM return, or some equally restrictive
canonical/minimal geometry.  Neither positive characteristic nor the word
“Cremona” supplies that promotion automatically.

## 6. Reproduction

The exact \((2,3,\infty)\) triangle replay checks the projector, three
reflection correspondences, finite Coxeter relations, parabolic product,
cycle-composition law, and the rank formula:

```bash
PYTHONPATH=src:. python3 scripts/certify_universal_weak_chow_action_carrier.py --check
PYTHONPATH=src:. pytest -q tests/test_universal_weak_chow_action_carrier.py
```

The frozen archive is
[`universal_weak_chow_action_carrier.json`](../data/enumerations/universal_weak_chow_action_carrier.json).
