# The noncompact K3 cusp-divisibility host gate

**Status (2026-08-11): cusp-local no-go, with an exact global repair
criterion.**  After the noncompact full-trace field collapse to
\(k=\mathbf Q\), a K3 cusp does not by itself give a finite list of integral
hosts.  The affine root diagram controls the negative-definite fibre-root
part, but it does not control the divisibility of the primitive isotropic
fibre class.  That divisibility is the numerical multisection index and can
be arbitrarily large in primitive K3-embeddable lattices with one fixed
physical affine cusp packet.

Thus

\[
 \boxed{\text{cusp-local host finiteness fails.}}
\]

The exact missing global property is a **bounded spanning physical-root
packet**: in each Picard rank, select a full rational basis of actual facet
\((-2)\)-curves whose complete intersection Gram lies in one bounded
integral box.  Such a lemma immediately gives finitely many saturated Picard
lattices.  Global use of finite volume is the target hypothesis not tested
by the family below and is the natural source of such a packet; same-host
reflection maximality is downstream and cannot replace it.

This note does **not** construct infinitely many finite-volume
reflection-maximal K3 chambers.  The family below isolates the parameter
which every attempted cusp-only proof leaves free.  Its complete root
chambers are not asserted to have finite volume, so it is not a
counterfamily to the full target theorem.

## 1. An exact unbounded K3-embeddable cusp family

For every integer \(d\ge2\), let

\[
 M_d=U(d)\oplus\langle-2\rangle
\]

with basis \(u,v,a\) and Gram matrix

\[
 Q_d=
 \begin{pmatrix}
 0&d&0\\
 d&0&0\\
 0&0&-2
 \end{pmatrix}.
 \tag{1.1}
\]

It is even of signature \((1,2)\), and

\[
                 |\operatorname{disc}M_d|=2d^2.       \tag{1.2}
\]

Consequently the lattices are pairwise nonisometric.  The vector \(u\) is
primitive and isotropic, but

\[
                 \operatorname{div}_{M_d}(u)
       =\gcd\{(u,x):x\in M_d\}=d.                     \tag{1.3}
\]

Define

\[
                  \alpha=a,\qquad \beta=u-a.
\]

Then

\[
 \alpha^2=\beta^2=-2,\qquad
 (\alpha,\beta)=2,\qquad
 \alpha+\beta=u.                                    \tag{1.4}
\]

Thus every \(M_d\) has the same fixed affine
\(\widetilde A_1\) cusp Gram

\[
                    \begin{pmatrix}-2&2\\2&-2\end{pmatrix}. \tag{1.5}
\]

The transverse square-minus-two root

\[
                         c=v-a
\]

records the escaping parameter:

\[
 c^2=-2,\qquad (u,c)=d,\qquad
 (\alpha,c)=2,\qquad (\beta,c)=d-2.             \tag{1.6}
\]

In fact \(\alpha,\beta,c\) are an integral basis of \(M_d\), with Gram

\[
 \begin{pmatrix}
 -2&2&2\\
 2&-2&d-2\\
 2&d-2&-2
 \end{pmatrix},
 \qquad \det=2d^2.                                   \tag{1.7}
\]

All three reflections are integral because a square-minus-two root acts by
\(x\mapsto x+(x,r)r\).  Formula (1.7) makes the obstruction completely
visible: the fixed affine packet leaves one transverse root pairing, the
cusp width, and the host discriminant unbounded.

## 2. The family really lies in the K3 lattice universe

Let \(e_1,f_1,e_2,f_2\) be standard bases of two copies of \(U\), and let
\(a\) be a root in one \(E_8(-1)\) summand of the K3 lattice.  The map

\[
 u\longmapsto e_1,\qquad
 v\longmapsto d f_1+e_2,\qquad
 a\longmapsto a                                      \tag{2.1}
\]

is an isometric primitive embedding of \(M_d\) into
\(U^3\oplus E_8(-1)^2\).  Primitivity is elementary: the coordinate matrix
of (2.1) has a maximal minor equal to one.  Period surjectivity and avoidance
of the Noether--Lefschetz hyperplanes therefore give an exact-Picard complex
K3 surface with \(\operatorname{NS}(X)\simeq M_d\) for every \(d\).

After Weyl reframing, \(u\) lies on the nef boundary.  A primitive nef
isotropic K3 class defines a genus-one fibration.  The roots perpendicular
to \(u\) are

\[
                       m u\pm a,\qquad m\in\mathbf Z,
\]

and an affine chamber at this cusp has \(\alpha=a\) and
\(\beta=u-a\) as consecutive simple roots.  Under the K3 ample-cone
dictionary these simple roots are physical irreducible rational
\((-2)\)-curves.  The numerical multisection index of the fibration is
exactly (1.3).

For completeness, this is an actual multisection degree, not merely a gcd in
an abstract lattice.  For any \(k>0\), the class
\(D_k=v+ku\) has \(D_k^2=2kd>0\) and \(D_k.u=d\).  Riemann--Roch makes
\(D_k\) effective (the opposite sign has negative intersection with the nef
class \(u\)).  Every irreducible component has nonnegative \(u\)-degree,
and every such degree is divisible by \(\operatorname{div}(u)=d\).  Since
their sum is \(d\), exactly one horizontal component has degree \(d\).
No smaller positive multisection degree is possible by divisibility.

Hence neither of the following bounds \(d\):

- the Picard-rank ceiling, since every example has rank three;
- primitive K3 embeddability, by (2.1);
- the affine ADE type or its component multiplicities, since (1.5) is fixed;
- the fact that the incident roots are physical \((-2)\)-curves; or
- integrality of every displayed root reflection.

This is the precise genus-one torsor/multisection-index escape anticipated
by the cusp analysis.

## 3. Why finite volume has not disappeared

For a genuine finite-volume cusp, the incident facet roots span
\(u^\perp\), but that subspace is degenerate of codimension one in the
Lorentz space.  It therefore determines neither a transverse lattice vector
nor the glue of \(\mathbf Zu\) inside the full Picard lattice.  In (1.1),
the entire incident real-root arrangement in \(u^\perp\) is independent of
\(d\), while the pairing with every transverse lattice vector is divisible
by \(d\).

The family does not prove that finite volume permits arbitrary \(d\).  On
the contrary, any proof of finite host finiteness must use finite volume
*globally*, away from a single horospherical link, to find a transverse wall
at uniformly bounded arithmetic height.  Merely knowing that the
horospherical section is bounded proves

\[
                \operatorname{span}_{\mathbf Q}{\text{incident roots}\}
                 =u^\perp,
\]

which is exactly one dimension short of what host finiteness needs.

Same-host reflection maximality does not repair that local deficit.  Once a
host \(M\) and an exact finite-volume chamber have been supplied, maximality
saturates the facet group to all reflections of \(M\).  It does not compare
different \(M_d\), bound \(\operatorname{div}(u)\), or manufacture a bounded
transverse vector before the host is known.

## 4. The exact global repair lemma

The following statement is sufficient and is the sharp next theorem to
prove independently.

> **Bounded spanning-root lemma.**  For each \(3\le\rho\le20\), there is a
> constant \(B(\rho)\) such that every finite-volume exact full-Picard K3
> root chamber of rank \(\rho\) has physical facet roots
> \(E_1,\ldots,E_\rho\) which span
> \(\operatorname{NS}(X)\otimes\mathbf Q\) and satisfy
> \(0\le E_i.E_j\le B(\rho)\) for \(i\ne j\).

Here “span” is essential.  A bounded affine cusp packet spans only
\(u^\perp\), and any bounded collection of unrelated local diagrams can
leave the same divisibility escape.

Assume the lemma.  For fixed \(\rho\), the nonsingular Gram matrix

\[
                         G=(E_i.E_j)
\]

has diagonal \(-2\) and integral off-diagonal entries in a finite box.
There are therefore finitely many possible root lattices
\(R=\sum_i\mathbf ZE_i\).  Since \(R\) spans the Picard space and the Picard
lattice is integral,

\[
                 R\subseteq\operatorname{NS}(X)\subseteq R^*. \tag{4.1}
\]

The discriminant group \(R^*/R\) is finite, so it has only finitely many
subgroups and hence only finitely many integral, even K3-compatible
overlattices.  Taking the finite union over \(3\le\rho\le20\) proves finite
hosts.

After that, exact complete-root return already gives one Weyl-chamber orbit
per host.  In the more general native-host formulation, same-host reflection
maximality gives the same conclusion by saturation to the complete mirror
arrangement.  Thus same-host reflection maximality is downstream of (4.1),
not a substitute for it.

## 5. Dependency audit

A known finite-polyhedral-surface argument states a stronger uniform form
of the bounded spanning-root lemma: it selects \(\rho\) actual exceptional
curves spanning \(\operatorname{NS}_{\mathbf Q}\), with normalized pairings
strictly below \(62\).  For K3 roots this would give
\(E_i.E_j\le61\), and Section 4 would immediately finish the noncompact
finite-host theorem, including nonsimple cusps.

That spanning assertion is global narrow-part hyperbolic geometry; it is not
a consequence of the elliptic-fibre Euler bound or of the K3 period theorem.
Its published proof chain ultimately imports the classical hyperbolic
diagram method.  It has not been independently reproduced in this report,
so it is recorded as the exact repair gate rather than silently used as a
new AG proof of reflection-lattice finiteness.

The rigorous conclusion frozen here is therefore:

\[
\boxed{
 \begin{gathered}
 \text{cusp/ADE/NS/K3-embedding data alone leave }\operatorname{div}(u)
 \text{ unbounded};\\
 \text{a global bounded spanning facet-root packet would imply finite hosts};\\
 \text{deriving that packet independently from finite volume remains the gate.}
 \end{gathered}}
\]

No dimension classification, finite-volume claim for \(M_d\), raw-carrier
essential surjectivity, or Property 4 conclusion is made.
