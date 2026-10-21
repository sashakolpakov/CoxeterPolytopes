# The isotropic weak-log class has a period-torsion gate

**Status (2026-08-10): proved over \(\mathbf C\), with an exact lattice and
cohomology-recursion replay.**  This note isolates a semiampleness datum that
is invisible to the Néron--Severi lattice even in the simplest smooth
genus-one case.

## 1. The theorem

Let \(S\) be a smooth projective rational surface over \(\mathbf C\), and let
\(F\subset S\) be a smooth connected genus-one curve, regarded as an effective
Cartier divisor whose class is primitive in
\(\operatorname {Pic}(S)=\operatorname {NS}(S)\).  Assume

\[
                  F^2=K_S\mathbin{\cdot}F=0,
        \qquad F\text{ is nef}.                              \tag{1.1}
\]

Set

\[
                  N:=\mathcal O_F(F)\in\operatorname {Pic}^0(F). \tag{1.2}
\]

Here \(N\) has degree \(F^2=0\), which explains its membership in
\(\operatorname {Pic}^0(F)\).

Then the following are equivalent:

1. \(F\) is semiample;
2. \(|rF|\) moves for some \(r>0\);
3. \(N\) is torsion in \(\operatorname {Pic}^0(F)\).

More precisely, if \(N\) has exact order \(m\), then for every integer
\(r\geq 0\),

\[
 h^0\!\left(S,\mathcal O_S(rF)\right)
       =\left\lfloor\frac r m\right\rfloor+1,
 \qquad
 h^1\!\left(S,\mathcal O_S(rF)\right)
       =\left\lfloor\frac r m\right\rfloor.                 \tag{1.3}
\]

In particular, \(|rF|\) is fixed for \(0<r<m\), while \(|mF|\) is a
base-point-free pencil.  Thus \(m\) is both the order of the normal character
and the first moving/semiample index.  If \(N\) is non-torsion, then

\[
 h^0\!\left(S,\mathcal O_S(rF)\right)=1,
 \qquad h^1\!\left(S,\mathcal O_S(rF)\right)=0
                    \quad\text{for every }r\geq0,            \tag{1.4}
\]

and no positive multiple of \(F\) moves.

The nef hypothesis is retained because it is the hypothesis needed by the
weak-log compiler, and primitivity is the normalization of its isotropic ray.
The cohomological proof below in fact uses neither nefness nor primitivity: it
uses the stated smoothness, connectedness, effectivity, rationality, and
intersection equalities.

## 2. Exact proof

For every \(r\geq1\), restriction to the Cartier divisor \(F\) gives

\[
 0\longrightarrow\mathcal O_S((r-1)F)
 \longrightarrow\mathcal O_S(rF)
 \longrightarrow N^r\longrightarrow0.                       \tag{2.1}
\]

Because \(S\) is rational, \(\chi(\mathcal O_S)=1\) and
\(H^0(S,K_S)=0\).  Riemann--Roch and (1.1) therefore give

\[
 \chi\!\left(\mathcal O_S(rF)\right)
   =1+\frac{rF\mathbin{\cdot}(rF-K_S)}2=1.                   \tag{2.2}
\]

Moreover,

\[
 H^2(S,\mathcal O_S(rF))
   \simeq H^0(S,\mathcal O_S(K_S-rF))^\vee=0.                \tag{2.3}
\]

Indeed, if \(K_S-rF\) were effective, adding the effective divisor \(rF\)
would make \(K_S\) effective, contrary to \(H^0(S,K_S)=0\).

A degree-zero line bundle on the smooth elliptic curve \(F\) has a nonzero
section if and only if it is trivial.  Suppose first that \(N\) has exact
order \(m\).  For \(1\leq r<m\), the right side of (2.1) has no section, so
induction from \(h^0(S,\mathcal O_S)=1\) gives \(h^0(rF)=1\).  Equations
(2.2)--(2.3) then give \(H^1(S,(m-1)F)=0\).  At \(r=m\), the restriction map
in (2.1) is consequently onto \(H^0(F,N^m)=H^0(F,\mathcal O_F)\).  Choose a
lift \(u\in H^0(S,mF)\) whose restriction to \(F\) is \(1\).  Its zero divisor
\(G\sim mF\) is disjoint from \(F\).  The canonical section with zero divisor
\(mF\) and the section \(u\) have no common zero.  Hence \(|mF|\) is a
base-point-free pencil.

The same argument determines every later term, not just the first one.  If
\(m\nmid r\), then \(H^0(F,N^r)=0\), so (2.1) does not change \(h^0\).  If
\(r=qm\), the section \(u^q\) restricts nontrivially to \(F\), so the
restriction map in (2.1) is onto and \(h^0\) increases by one.  Starting at
\(h^0(0)=1\) proves the first formula in (1.3); (2.2)--(2.3) prove the second.
This is a direct restriction-sequence calculation: it does not assume that
the pencil morphism already has connected fibres, nor does it silently use a
Stein factorization.  In particular, the explicit sections \(u^q\) themselves
certify the required restriction-map ranks at all multiples of \(m\).

If \(N\) is non-torsion, then \(H^0(F,N^r)=0\) for every \(r>0\).  Sequence
(2.1) gives \(h^0(rF)=h^0((r-1)F)=1\), and (2.2)--(2.3) give \(h^1(rF)=0\).

Finally, if \(|rF|\) moves, choose the least positive such \(r\).  The new
section cannot come from \(H^0((r-1)F)\), so (2.1) gives a nonzero section of
\(N^r\); hence \(N^r\simeq\mathcal O_F\).  Alternatively, if \(rF\) is
base-point-free, its restriction \(N^r\) is a globally generated degree-zero
line bundle and is therefore trivial.  This proves all three equivalences.

## 3. Two blowups with the same marked lattice and opposite answers

Fix a smooth plane cubic \(C\subset\mathbf P^2\) and nine distinct points
\(p_1,\ldots,p_9\in C\).  On

\[
             S=\operatorname {Bl}_{p_1,\ldots,p_9}\mathbf P^2
\]

let \(H\) be the pullback of a line and \(E_i\) the exceptional curves.  The
strict transform of \(C\) satisfies the following equality of divisor classes
in \(\operatorname {Pic}(S)=\operatorname {NS}(S)\):

\[
 [F]=3H-E_1-\cdots-E_9=[-K_S],\qquad F^2=K_S\cdot F=0.       \tag{3.1}
\]

Thus \(F\sim-K_S\) linearly, not merely numerically.  The displayed vector is
primitive.

It is nef.  Indeed, \(F\cdot E_i=1\).  If an irreducible plane curve \(Q\)
does not contain \(C\), Bézout and smoothness of \(C\) at every \(p_i\)
give

\[
 F\cdot\widetilde Q
   =3\deg Q-\sum_i\operatorname {mult}_{p_i}Q\geq0.          \tag{3.2}
\]

The remaining irreducible curve \(F\) has self-intersection zero.  Thus these
checks cover every irreducible curve on the blowup.

Under the identification \(F\simeq C\), its normal character is

\[
 N\simeq\mathcal O_C(3)(-p_1-\cdots-p_9)\in\operatorname {Pic}^0(C). \tag{3.3}
\]

Now make two choices.

* **Trivial character.**  Choose another cubic meeting \(C\) transversely in
  the nine points.  They are the base locus of the resulting cubic pencil.
  After blowing them up, \(|F|\) is a base-point-free pencil and (3.3) is
  trivial.
* **Non-torsion character.**  Choose eight distinct points first.  As the
  ninth point varies on \(C\), the character (3.3) varies by a translate of
  the Abel--Jacobi isomorphism \(C\simeq\operatorname {Pic}^0(C)\).  Only
  countably many choices make it torsion, and finitely many choices repeat an
  earlier point.  Choose outside those sets.  Then \(F\) is still nef, but
  every \(|rF|\) is the single divisor \(rF\).

Both surfaces have the identical marked Néron--Severi lattice

\[
 \operatorname {NS}(S)\simeq I_{1,9}
   =\langle1\rangle\oplus\langle-1\rangle^{\oplus9},
\]

with the same vectors \(K_S=-3H+\sum E_i\) and
\([F]=3H-\sum E_i=[-K_S]\), and \(F\) is nef in both.  Their semiampleness
answers are nevertheless opposite.  This is a same-\((\operatorname {NS},K,F)\)
obstruction, not a claim that the two complete nef cones are marked-isomorphic.

## 4. Exact consequence for the polytope-to-surface compiler

In the boundary-free weak-log specialization,

\[
                         D:=-K_S\sim F.
\]

The example proves that the lattice, canonical vector, square-zero condition,
effectivity of the smooth genus-one representative, and actual nefness of
\(D\) do not determine whether \(D\) is semiample.  A correct isotropic lane
must additionally return, for example,

* an exact base-point-free certificate for a multiple of \(D\);
* the finite order of \(\mathcal O_F(F)\) when a smooth genus-one
  representative is supplied; or
* another independently verified morphism certificate.

This period gate is additional to, not a replacement for, the other carrier
gates.  In particular, this theorem package makes **no claim** that

* every isotropic weak-log class has a smooth genus-one representative;
* the displayed blowups realize a prescribed Coxeter polytope;
* a complete nef chamber or all effective facet curves have been returned;
* any relevant chamber has finite hyperbolic volume;
* any reflection group is saturated or maximal; or
* the construction is exhaustive for raw arithmetic polytopes.

The accompanying
[`certifier`](../scripts/certify_rational_isotropic_weak_log_period_gate.py)
replays the integral \(I_{1,9}\) calculation, Riemann--Roch polynomial, and
the exact torsion-order recursion.  Its finite replay does not replace the
cohomological or moduli arguments proved above.
