# Genus-three fibrations do not see the vertical blowup depth

**Status (2026-08-10): exact birational-invariance no-go proved; no
unbounded complete or reflection-maximal carrier is claimed.**  This note
isolates what remains after affine-H3 abundance has produced the
base-point-free genus-three pencil

\[
                 F=4E+B_1+B_2+B_3+B_4.
\]

The relatively minimal fibration, its generic curve, its moduli map, its
relative canonical model, and its Hodge bundle are all unchanged by
vertical blowups.  On a nonminimal smooth model the two numerical terms in
the relative Noether formula change by opposite amounts.  Consequently
slope, Hodge, stable-reduction, KSBA, and ordinary relative-Chern arguments
bound the relatively minimal core but cannot bound the vertical blowup
depth.

This failure is realized by the correct surface roots.  If a different
fibre contains a smooth rational \((-1)\)-component, arbitrarily many free
blowups on its newest exceptional component give a chain of physical
\((-2)\)-roots ending in a physical \((-1)\)-root.  Its pure-tail group is
the full signed-permutation group, and the one mixed curve already has an
orbit spanning the whole tail.  Thus canonical coloring, local reflection
saturation, and the known finite-volume attachment condition add no bound.

The construction is deliberately local.  It does not prove that the full
Mori cone after every blowup still has only the declared curves, that its
nef cone has finite volume, or that its facet group is reflection-maximal.
Those global hypotheses remain the only possible source of a depth bound.

## 1. Exact vertical-birational theorem

Let

\[
                         f_0:X\longrightarrow\mathbf P^1
\]

be a genus-three fibration on a smooth projective rational surface.  Let
\(\pi:S\to X\) be a composition of \(m\) point blowups, with every centre
contained in a fibre of the induced fibration, and put
\(f=f_0\circ\pi\).  No assumption on the incidence of the centres is needed
in this section.

### Theorem 1.1 (vertical blowup erasure)

The following data are identical for \(f\) and \(f_0\):

1. the generic genus-three curve and hence the rational moduli map from the
   base;
2. the relative canonical algebra
   \(\bigoplus_{n\ge0} f_*\omega_{S/\mathbf P^1}^{\otimes n}\), and hence
   the relative canonical model;
3. the Hodge bundle \(f_*\omega_{S/\mathbf P^1}\), whose determinant has
   degree

   \[
                              \chi_f=3;                    \tag{1.1}
   \]

4. the relatively minimal model and, after any required finite base
   change, the stable model.

On the smooth total spaces the exact changes are

\[
\begin{array}{c|ccccc}
 &\rho&c_2&K^2&K_f^2&e_f\\ \hline
 S-X&m&m&-m&-m&m,
\end{array}                                                \tag{1.2}
\]

where

\[
 K_f=K_S+2F,\qquad
 e_f=e(S)-e(\mathbf P^1)e(F_{\rm gen}).                    \tag{1.3}
\]

For every smooth rational genus-three model one has

\[
 \boxed{K_f^2=26-\rho(S),\qquad e_f=\rho(S)+10,
        \qquad K_f^2+e_f=36=12\chi_f.}                    \tag{1.4}
\]

Thus the raw pair \((K_f^2,e_f)\) records a blowup by the translation
\((-1,+1)\), but every relative Noether expression cancels it.  The
nonnegative inequality \(K_f^2\ge0\) is available only after relative
minimalization: on the blown-up surface the newest vertical exceptional
curve \(A\) satisfies

\[
                              K_f.A=-1.                    \tag{1.5}
\]

In particular, applying a relatively-minimal slope or nefness inequality
directly on \(S\) would silently assume away the tail one is trying to
bound.

### Proof

A vertical blowup is an isomorphism over the generic point of the base, so
the generic curve and the moduli map do not change.  For one blowup with
exceptional curve \(A\),

\[
       \omega_{S/\mathbf P^1}
       =\pi^*\omega_{X/\mathbf P^1}\otimes\mathcal O_S(A). \tag{1.6}
\]

More generally the discrepancy divisor of a composition of point blowups
is effective and exceptional.  Normality of the smooth target gives

\[
             \pi_*\mathcal O_S(nD_{\rm exc})=\mathcal O_X
             \qquad(n\ge0):                               \tag{1.7}
\]

a rational function allowed poles only over finitely many codimension-two
points extends across those points.  Projection formula applied to (1.6)
and all its powers proves equality of the relative canonical algebras and
of the Hodge bundles.  The relative canonical model is their relative
Proj.  Uniqueness of relative minimal and stable models gives the last
birational assertions.

Blowing up a point adds one to \(\rho\) and \(c_2=e\), and subtracts one
from \(K^2\).  Since the fibre class pulls back and is orthogonal to every
exceptional class,

\[
                    K_f^2\longmapsto K_f^2-1.             \tag{1.8}
\]

The generic fibre is unchanged, so \(e_f\) increases by one.  Finally, on
a rational surface

\[
 K_S^2=10-\rho(S),\quad e(S)=\rho(S)+2,\quad
 F^2=0,\quad K_S.F=2g-2=4.                                \tag{1.9}
\]

Substitution into (1.3) proves (1.4).  Formula (1.1) is
\(\chi_f=\chi(\mathcal O_S)+(g-1)=1+2=3\).  Adjunction gives
\(K_S.A=-1\), while \(F.A=0\), proving (1.5). \(\square\)

## 2. What stable and Chern inequalities actually bound

Let \(S_{\min}\) be the relatively minimal smooth model.  All slope and
relative-canonical nefness inequalities apply to \(S_{\min}\), not to an
arbitrary vertical blowup.  At the weakest level,

\[
 0\le K_{f,\min}^2=26-\rho(S_{\min})
 \quad\Longrightarrow\quad \rho(S_{\min})\le26.           \tag{2.1}
\]

Any stronger slope inequality still changes only the bound on this core.
If \(S\to S_{\min}\) contains \(m\) vertical blowups, then

\[
 \rho(S)=\rho(S_{\min})+m,
 \quad K_f^2=K_{f,\min}^2-m,
 \quad e_f=e_{f,\min}+m.                                  \tag{2.2}
\]

On the stable side, after the usual semistable reduction, a vertical blowup
inserts an unstable rational component: a smooth free centre gives a tail,
whereas a nodal centre gives a bridge.  Stabilization contracts it.  Hence
the stabilized dual graph, the stable boundary degree, the Hodge class, and
the Kodaira--Spencer or moduli map remain those of the core.  On the smooth
side, the extra Euler contribution is exactly cancelled by the loss in
\(K_f^2\) in (1.4).

This proves a formal impossibility statement: no function only of the
generic curve, moduli map, relative canonical algebra/model, Hodge bundle,
stable fibre data, \(\chi_f\), or relatively minimal \(K_f^2\) can bound
\(m\).  Those data are literally constant on every vertical blowup tower.
Nor can the raw Noether pair supply a bound without a new lower bound for
the nonminimal \(K_f^2\); equation (1.5) explains why the standard lower
bound is unavailable.

## 3. The correctly colored arbitrary free tail

The erasure is compatible with the local physical alphabet.  Suppose a
fibre of \(f_0\), different from the distinguished affine-H3 fibre,
contains a smooth rational curve \(C_0\) with \(C_0^2=-1\).  Blow up a
smooth point of \(C_0\), away from its intersections.  Thereafter blow up
a free point of the newest exceptional curve, again away from all
intersections.  In total-transform coordinates write the exceptional
classes as \(e_1,\ldots,e_m\).  The displayed vertical curves are

\[
 C=C_0-e_1,\qquad
 R_i=e_i-e_{i+1}\ (1\le i<m),\qquad A=e_m.                \tag{3.1}
\]

They are actual smooth rational curves with

\[
 C^2=R_i^2=-2,qquad A^2=-1,                              \tag{3.2}
\]

and form a simple chain.  Their canonical colors are exact:

\[
 K_S.C=K_S.R_i=0,qquad K_S.A=-1.                         \tag{3.3}
\]

In the doubled adapted host the corresponding roots have norms \(-4\)
and \(-2\), and

\[
 q_M(K_S,C)=q_M(K_S,R_i)=0,qquad q_M(K_S,A)=-2.          \tag{3.4}
\]

Because the blowups occur in another fibre, every distinguished
square-minus-four component \(B_j\) is orthogonal to all \(e_i\).  Thus
all pure-tail reflections preserve

\[
 M=\left\langle\operatorname {NS}(S),B_1/2,\ldots,B_4/2,
                  \ldots\right\rangle,qquad q_M=2q_{\rm NS}.            \tag{3.5}
\]

The roots

\[
                  e_1-e_2,\ldots,e_{m-1}-e_m,e_m          \tag{3.6}
\]

generate the full signed-permutation group \(W(B_m)\) on the exceptional
summand.  The only displayed mixed curve is \(C=C_0-e_1\); its exceptional
projection is \(-e_1\), and therefore

\[
             \operatorname {span}_{\mathbf R}
             \bigl(W(B_m)(-e_1)\bigr)
             =\operatorname {span}_{\mathbf R}
             (e_1,\ldots,e_m).                            \tag{3.7}
\]

Consequently this tower has, for every \(m\), the right effectivity,
adjunction colors, adapted-host parity, pure-tail reflection saturation,
and the exact orbit-span condition known to be necessary for finite
volume.  None of these statements says that the remaining global walls
exist.

## 4. Exact surviving gate

The affine-H3 fibration route has now separated the problem cleanly.

* Abundance and the nonreduced-star normal-bundle argument produce \(|F|\)
  and genus three.
* Relative canonical and slope geometry bounds only the relatively minimal
  core.
* Stable, Hodge, and moduli invariants erase every vertical free tail.
* Canonical coloring, pure-tail saturation, and finite-volume orbit span
  allow the arbitrary chain (3.1).

Therefore a depth theorem must use information which is not invariant under
vertical blowup.  Two precise sufficient kinds of input remain:

1. an absolute **total marked degree** bound which charges every physical
   tail facet, rather than a separate bound for each facet; or
2. a finite **equivariant contraction orbit/descent** theorem forcing the
   reflection group to act through a bounded collection of relatively
   minimal models.

Bare fibration invariants cannot substitute for either input.  Conversely,
this note does not refute a depth theorem from the complete-nef and
reflection-maximal hypotheses; it proves exactly why that theorem must use
their global mixed-wall content rather than the stable genus-three model.
