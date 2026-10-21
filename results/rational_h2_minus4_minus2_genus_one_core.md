# The minimal H2 pair: a sharp rank-12 adjoint core and a saturated exceptional tail

**Status (2026-08-10): proved core theorem; no original-rank theorem or
reflection-maximal counterfamily is claimed.** This note treats the smallest
standard-host H2 incidence

\[
             B^2=-4,\qquad R^2=-2,\qquad B.R=2.              \tag{0.1}
\]

It proves, using only rational-surface cohomology, adjunction, contractions,
and intersection theory, that the adjoint-minimal core has Picard rank at
most \(12\). The bound is sharp for smooth rational pairs. It does **not**
bound the Picard rank before relative minimalization. Successive free
blowups give an arbitrary exceptional \(B_m\) tail whose visible reflections
already generate its full signed-permutation group, so a naive local
reflection-overgroup argument cannot remove that tail.

No Nikulin--Vinberg classification, ABSW or Prasad theorem, arithmetic
volume, or polytope census is used.

## 1. The genus-one adjoint core

Let \(S\) be a smooth projective rational surface over an algebraically
closed field. Let \(B,R\subset S\) be distinct smooth rational curves
satisfying (0.1), and put

\[
                         D=B+R,\qquad L=K_S+D.                \tag{1.1}
\]

The core theorem only uses the scheme-theoretic intersection number
\(B.R=2\). Thus it applies both when \(B\cap R\) consists of two distinct
transverse points and when the two contacts coalesce. The singularity
statement in Section 3 will distinguish those cases.

### Theorem 1.1 (minimal H2 core)

The linear system \(|K_S+B+R|\) contains a unique effective divisor \(A\),
with \(A=0\) allowed, and

\[
                         \operatorname {Supp}(A)\cap D=\varnothing.       \tag{1.2}
\]

There is a finite sequence

\[
 S=S_m\longrightarrow S_{m-1}\longrightarrow\cdots
      \longrightarrow S_0                                             \tag{1.3}
\]

in which every map blows down a smooth rational square-minus-one curve
disjoint from \(D\), the pair is unchanged near \(D\), and

\[
                         K_{S_0}+B_0+R_0\quad\text{is nef}.                \tag{1.4}
\]

The terminal core satisfies

\[
                         \boxed{\rho(S_0)\leq12}.                         \tag{1.5}
\]

Equivalently, its projectivized nef cone has hyperbolic dimension at most
\(11\).

### Proof

The reduced Cartier curve \(D\) is connected and Gorenstein. Adjunction
gives

\[
 K_S.B=2,\qquad K_S.R=0,\qquad D^2=-2,\qquad K_S.D=2,                       \tag{1.6}
\]

and hence

\[
 p_a(D)=1+\frac{D^2+K_S.D}{2}=1.                                          \tag{1.7}
\]

Moreover

\[
 \deg(\omega_D|_B)=(K_S+D).B=0,\qquad
 \deg(\omega_D|_R)=(K_S+D).R=0.                                         \tag{1.8}
\]

Both restrictions are therefore trivial line bundles on \(\mathbf P^1\).
Serre duality on the connected Gorenstein genus-one curve gives
\(h^0(D,\omega_D)=1\). A nonzero dualizing section cannot restrict to zero
on one component: compatibility on the nonempty intersection scheme would
make its degree-zero section on the other \(\mathbf P^1\) vanish as well.
Its restrictions to both components are consequently nonzero constants,
so it has no zero anywhere on \(D\). Thus

\[
                              \omega_D\simeq\mathcal O_D.                 \tag{1.9}
\]

This argument uses the length-two intersection scheme, not nodality.

The canonical exact sequence is

\[
 0\longrightarrow\mathcal O_S(K_S)
  \longrightarrow\mathcal O_S(K_S+D)
  \longrightarrow\omega_D\longrightarrow0.                              \tag{1.10}
\]

Rationality gives

\[
 H^0(S,K_S)=0,\qquad
 H^1(S,K_S)=H^1(S,\mathcal O_S)^\vee=0.                                  \tag{1.11}
\]

Therefore restriction induces an isomorphism

\[
 H^0(S,K_S+D)\simeq H^0(D,\omega_D)\simeq k.                             \tag{1.12}
\]

The zero divisor \(A\) of its generator is unique, effective, and disjoint
from \(D\), proving (1.2).

Suppose \(A\) is not nef. Then some irreducible curve \(E\) has
\(A.E<0\). Since \(A\) is effective, \(E\) is a component of \(A\) and
\(E^2<0\). It is disjoint from \(D\), so

\[
                              K_S.E=A.E<0.                                \tag{1.13}
\]

Write \(E^2=-q\), with \(q>0\). Adjunction gives

\[
                       K_S.E=2p_a(E)-2+q<0.                               \tag{1.14}
\]

The only possibility is \(p_a(E)=0\) and \(q=1\). Hence \(E\) is a smooth
rational square-minus-one curve, and it can be blown down. The contraction
is an isomorphism near \(D\), so (0.1) remains true on the new smooth
rational surface. Reapply the same construction there. Picard rank drops at
every step, so the process terminates exactly when its unique adjoint
divisor is nef. This proves (1.3)--(1.4) without any extremal-ray or
fiber-type classification.

For every smooth projective rational surface, Noether's formula gives
\(K^2=10-\rho\). On the terminal surface, (1.6) is unchanged, and therefore

\[
 \begin{aligned}
 (K_{S_0}+D_0)^2
   &=K_{S_0}^2+2K_{S_0}.D_0+D_0^2\\
   &=(10-\rho(S_0))+4-2\\
   &=12-\rho(S_0).
 \end{aligned}                                                            \tag{1.15}
\]

A nef divisor on a smooth surface has nonnegative square. Equation (1.15)
proves (1.5). \(□\)

## 2. Sharpness of the core bound

The number \(12\) is attained without a classification. On
\(X=\mathbf P^1\times\mathbf P^1\), choose two smooth \((1,1)\)-curves
\(\bar B,\bar R\) meeting transversely in two points. Then

\[
 \bar B^2=\bar R^2=2,\qquad \bar B.\bar R=2,\qquad
                         \bar B+\bar R=-K_X.                              \tag{2.1}
\]

Blow up six distinct smooth points of \(\bar B\) and four distinct smooth
points of \(\bar R\), all away from \(\bar B\cap\bar R\). For the strict
transforms on the resulting surface \(S\),

\[
 B^2=-4,\qquad R^2=-2,\qquad B.R=2,\qquad
 K_S+B+R=0,\qquad \rho(S)=2+6+4=12.                                      \tag{2.2}
\]

Thus \(S\) is already an adjoint-minimal core and equality holds in
Theorem 1.1. This is a sharp surface-pair control; no complete-Mori,
finite-volume, or reflection-maximal claim is attached to this example.

## 3. Contraction and singularity scope

The H2 pair has negative-definite matrix

\[
 Q=\begin{pmatrix}-4&2\\2&-2\end{pmatrix},\qquad \det Q=4.                \tag{3.1}
\]

If the pair contracts to one normal point and

\[
                  K_S=f^*K_Y+a_BB+a_RR,
\]

then exact adjunction gives

\[
 Q\binom{a_B}{a_R}=\binom20,qquad
                         (a_B,a_R)=(-1,-1).                              \tag{3.2}
\]

When \(B\cap R\) consists of two distinct transverse points, this is an SNC
two-component cusp configuration. Its contraction is log canonical but not
klt. In particular, the klt quotient-orbifold \(c_2\) inequality used for
the weak-incidence wrong-log lane cannot be imported here.

The intersection matrix alone does not prove log canonicity. If the entire
intersection number is an ordinary tangency of order two, the first blowup
has discrepancy \(-1\) and produces a triple point of the two strict
transforms and the exceptional curve. Blowing up that triple point creates
discrepancy \(1-1-1-1=-2\). The tangential contraction is not log canonical.
Theorem 1.1 nevertheless remains valid because its dualizing-sheaf argument
did not assume SNC.

Nor does an index cover close H2. At equality \(B+R=-K_S\), whereas the
standard Coble/K3-cover equation is a smooth-enough branch in
\(|-2K_S|\). Formally, a double cover branched along \(D=-K_S\) has
\(K_{\widetilde S}=\pi^*(K_S+D/2)=\pi^*(K_S/2)\), not trivial canonical
class. The local index-one cover of the transverse lc cusp is not a global
K3 carrier. Thus no \(n\leq19\) K3 conclusion follows from (0.1).

The Hodge-index contribution is exactly (1.15). It bounds the nef adjoint
core, but it permits an arbitrary negative-definite orthogonal exceptional
summand before contraction. Log-BMY or orbifold-\(c_2\) input therefore does
not supply the missing total-rank theorem in this setup.

## 4. The free exceptional tail is already reflection-saturated

Let \(\pi:S_m\to S_0\) be obtained by \(m\) successive free blowups away
from \(D\): after the first blowup, put each new center on the newest
exceptional curve and away from its intersection with the preceding one.
In total-transform exceptional coordinates \(e_1,\ldots,e_m\), the physical
exceptional curves are

\[
 R_1=e_1-e_2,\ \ldots,\ R_{m-1}=e_{m-1}-e_m,\qquad E_m=e_m.                \tag{4.1}
\]

They form the chain

\[
 R_1-\cdots-R_{m-1}-E_m,qquad
 R_i^2=-2,\quad E_m^2=-1,                                               \tag{4.2}
\]

with adjacent physical intersection \(1\). Since all centers are disjoint
from \(D\),

\[
 \operatorname {NS}(S_m)=\pi^*\operatorname {NS}(S_0)
       \mathbin\perp\bigoplus_{i=1}^m\mathbf Ze_i,qquad e_i^2=-1.        \tag{4.3}
\]

For the standard H2 host

\[
 M(S)=\langle\operatorname {NS}(S),B/2\rangle,
                         \qquad q_M=2q_{\operatorname {NS}},              \tag{4.4}
\]

the half-boundary is pulled back, so

\[
                  M(S_m)=M(S_0)\mathbin\perp\langle-2\rangle^m.          \tag{4.5}
\]

The roots in (4.1) have \(q_M\)-squares
\((-4,\ldots,-4,-2)\) and adjacent pairings \(2\). Reflection in
\(e_i-e_{i+1}\) swaps \(e_i,e_{i+1}\); reflection in \(e_m\) changes the
sign of \(e_m\). Adjacent swaps and one coordinate sign change generate all
signed permutations. Hence the visible physical reflections generate

\[
 W(B_m)=O(\langle-2\rangle^m),\qquad |W(B_m)|=2^m m!.                     \tag{4.6}
\]

There is no missing pure-exceptional coordinate reflection with which to
enlarge the group: every such reflection is already in (4.6).

This tail is genuine algebraic geometry, not an abstract colored diagram.
If \(H_0\) is ample on \(S_0\), then \(H=\pi^*H_0\) is nef with
\(H^2>0\), and \(H.C=0\) for an irreducible curve \(C\) exactly when \(C\)
is one of the exceptional components. Taking \(S_0=\mathbf P^2\) makes
\([H]\) the positive ray cut out by the whole chain. Nevertheless this local
fact asserts neither that an H2 surface with such tails has a finite-volume
complete nef Coxeter chamber nor that it is reflection-maximal.

Consequently Theorem 1.1 gives only

\[
             \boxed{\rho(S_0)\leq12,\qquad
                    \rho(S)=\rho(S_0)+m\leq12+m},                         \tag{4.7}
\]

and does not bound \(m\).

There is a useful conditional degree closure. Suppose an integral ample
marking \(H\) has an absolute bound on the **total** physical-facet degree

\[
                  T=\sum_{C\ {\rm physical\ facet}}H.C,                 \tag{4.8}
\]

and suppose every component contracted by \(\pi\) is itself one of those
physical facets. Each such component then has positive integral \(H\)-degree,
so \(m\leq T\). Consequently

\[
       \rho(S)\leq12+T,
       \qquad \dim\mathbf P(\operatorname {Nef}S)\leq11+T.              \tag{4.9}
\]

The total-degree premise and the physical-facet premise are both essential.
In particular, a bound on each individual facet degree, with no bound on
their number or total degree, does not imply (4.9).

## 5. Conditional trace descent and the exact open gate

Suppose a full trace lattice \(T\) of rank \(d(n+1)\) is carried
isometrically by \(\operatorname {NS}(S)\), and suppose it independently
descends through the particular H2 contraction \(\pi:S\to S_0\). Equivalently,
its image is orthogonal to the exceptional lattice in (4.3). Then

\[
 T\hookrightarrow\operatorname {NS}(S_0)
       \quad\Longrightarrow\quad
                         \boxed{d(n+1)\leq12}.                            \tag{5.1}
\]

For the rational field \(d=1\), this gives \(n\leq11\).

Descent is not automatic. In a literal full-nef carrier whose exceptional
curves in (4.2) are facets, their roots lie inside the visible full trace
root block. They have nonzero exceptional projection (already their
self-pairings are negative), so that block is generally not orthogonal to
the tail. Finite numerical orbit, full-trace descent, and reflection
maximality must not be conflated.

The exact remaining H2 theorem is therefore global:

> In a literal finite-volume complete-nef carrier, prove from geometric
> reflection maximality that the mixed core/tail geometry uniformly bounds
> the length \(m\), or supplies an independently certified finite numerical
> contraction orbit.

The \(B_m\) calculation proves that local tail saturation alone cannot do
this. Conversely, this note constructs no unbounded reflection-maximal
counterfamily. Until the global bridge is proved, the sharp rank-12 core is
a genuine AG reduction, not an absolute dimension-bound existence theorem.

## 6. Reproducibility

The exact pair arithmetic, discrepancy vector, sharp equality blowup count,
conditional trace inequality, and symbolic \(B_m\) reflection actions are
replayed by
`scripts/certify_rational_h2_minus4_minus2_genus_one_core.py` and frozen in
`data/enumerations/rational_h2_minus4_minus2_genus_one_core.json`.
The sampled tail lengths are regression controls for the closed formulas in
Section 4, not a census and not the logical source of the theorem.
