# Independent hostile audit: the noncompact K3 cusp-divisibility host gate

**Verdict (2026-08-11): PASS in the exact stated scope.**  The family

\[
                    M_d=U(d)\oplus\langle-2\rangle,
                    \qquad d\ge2,
\]

really consists of pairwise nonisometric, primitively K3-embeddable rank-three
lattices with one fixed physical affine (\widetilde A_1) cusp packet and
unbounded primitive-fibre divisibility.  That divisibility is not merely an
abstract gcd: it is the attained minimum multisection degree of the resulting
genus-one K3 fibration.

The family is a cusp-local no-go only.  Its full root chamber is not asserted
to have finite volume or to be reflection-maximal.  The proposed global
bounded-spanning-root statement is a sufficient repair lemma, not a theorem
proved by this package.

## 1. Frozen owner artifacts

This audit pins:

- `results/noncompact_k3_cusp_divisibility_host_gate.md`, SHA-256
  `68ce633e37d9cf3a8b7e94cdc199c0886e707e9ed14504b544162aee500d2e45`;
- `scripts/certify_noncompact_k3_cusp_divisibility_host_gate.py`, SHA-256
  `347618dea4971fa38215dbe7825639ddfb015254f80ccdd87aa65447285914a2`;
- `data/enumerations/noncompact_k3_cusp_divisibility_host_gate.json`, SHA-256
  `93b95baf9948068c85d720bbb063347930c78e62cab95bc831c02bdf8ca67b41`;
- `tests/test_noncompact_k3_cusp_divisibility_host_gate.py`, SHA-256
  `c03bd6997003e3c3a9630e042f5639372d1110ebbd55dfd2db6f21a28bdb9f48`.

## 2. The K3 embedding and exact-Picard realization are sound

In the basis (u,v,a), the Gram matrix is

\[
 Q_d=\begin{pmatrix}0&d&0\\d&0&0\\0&0&-2\end{pmatrix}.
\]

It has signature ((1,2)) and determinant (2d^2).  In two hyperbolic
planes and a primitive (E_8(-1))-root direction, the declared embedding is

\[
 u=e_1,\qquad v=df_1+e_2,\qquad a=a.
\]

Its pullback is exactly (Q_d).  More importantly, the rows corresponding
to (e_1,e_2,a) form an identity maximal minor.  Thus the cokernel is
torsion-free and the embedding remains primitive after inclusion into the
full K3 lattice.  This is stronger than a discriminant-form existence check
and works for every (d).

The orthogonal period domain has positive dimension.  Choosing a generic
period outside the countable Noether--Lefschetz hyperplanes gives a projective
K3 surface with exact Picard lattice (M_d).  Riemann--Roch chooses an
effective sign of the primitive isotropic class (u), and the standard K3
Weyl reduction moves it to the nef boundary.  A primitive nef isotropic K3
class is base-point-free and defines a genus-one pencil.

## 3. The affine roots and physical effectivity check

For (x=pu+qv+ra), the equations (x.u=0) and (x^2=-2) force

\[
                         q=0,\qquad r=\pm1.
\]

Hence all roots perpendicular to the cusp are exactly (mu\pm a).  The
pair

\[
                   \alpha=a,\qquad\beta=u-a
\]

is a consecutive affine simple pair, with

\[
 \alpha^2=\beta^2=-2,
 \qquad \alpha.\beta=2,
 \qquad \alpha+\beta=u.
\]

After choosing the corresponding root chamber as the K3 ample chamber, its
simple roots are classes of irreducible smooth rational curves.  Thus the
fixed affine packet is genuinely physical; it is not merely a pair of
abstract lattice roots.  The transverse root (v-a) records the escaping
intersection (u.(v-a)=d), but the no-go does not require it to be a simple
facet of a finite-volume chamber.

## 4. Divisibility is the actual multisection index

The pairing vector of (u) is ((0,d,0)), so

\[
                         \operatorname {div}(u)=d.
\]

Every divisor, and therefore every horizontal irreducible curve, has fibre
degree divisible by (d).  To show that the gcd is attained, take

\[
                           D_k=v+ku,\qquad k>0.
\]

Then

\[
                         D_k^2=2kd>0,
                         \qquad D_k.u=d.
\]

K3 Riemann--Roch implies that (D_k) or (-D_k) is effective.  The second
choice is impossible because it has negative intersection with the nef class
(u).  Write the resulting effective divisor as a sum of irreducible
components.  Nefness makes all their (u)-degrees nonnegative, while lattice
divisibility makes every degree a multiple of (d).  Their multiplicity-
weighted sum is (d).  Therefore exactly one horizontal component occurs
with multiplicity one and degree (d); all remaining components are
vertical.  A degree-(d) multisection exists, and no smaller positive degree
can exist.  The numerical and geometric indices coincide exactly.

## 5. What the cusp packet fails to see

The affine pair spans the degenerate plane (u^\perp), one rational
dimension short of (M_d\otimes\mathbf Q).  Its Gram matrix and null
relation are independent of (d), whereas

\[
 |\operatorname {disc}M_d|=2d^2,
 \qquad \operatorname {div}(u)=d,
 \qquad u.(v-a)=d
\]

escape.  Picard rank, primitive K3 embeddability, physical affine ADE data,
and integrality of the root reflections consequently do not prove
cusp-local host finiteness.

This does not remove finite volume from the target problem.  A finite-volume
cusp forces its incident roots to span (u^\perp), but that is still the
same deficient hyperplane.  A global argument must find transverse data;
the displayed family deliberately does not decide whether finite volume
forces that data to be bounded.

## 6. The bounded-spanning-root repair is sufficient

Suppose that for each (3\le\rho\le20), every target chamber supplies
(\rho) actual facet curves (E_1,\ldots,E_\rho) which span
(\operatorname {NS}(X)_\mathbf Q) and satisfy

\[
                    0\le E_i.E_j\le B(\rho)\quad(i\ne j).
\]

Their diagonal entries are all (-2), so for fixed rank their nonsingular
integral Gram matrices lie in a finite box.  For each resulting root lattice

\[
                        R=\sum_i\mathbf ZE_i,
\]

full rational spanning and Picard integrality give

\[
                 R\subseteq\operatorname {NS}(X)\subseteq R^*.
\]

The discriminant group (R^*/R) is finite.  Its finitely many subgroups
contain every integral intermediate lattice; evenness and K3 compatibility
only remove candidates.  Thus the lemma immediately yields finitely many
Picard hosts.  Same-host reflection maximality, when needed, enters only
after this step to identify the complete mirror group and its single chamber
orbit.

Spanning is essential.  Any bounded packet confined to (u^\perp), even a
complete affine cusp diagram, leaves the divisibility and transverse glue
uncontrolled.  The owner report correctly labels the derivation of a global
bounded spanning packet from finite volume as the remaining gate and does
not claim to have reproduced the external narrow-part diagram theorem.

## 7. Reproducibility and boundary

The separate executable replay is
`scripts/audit_noncompact_k3_cusp_divisibility_host_gate_independent_hostile.py`,
with regressions in
`tests/test_noncompact_k3_cusp_divisibility_host_gate_independent_hostile.py`.
It imports no owner certificate code.  It recomputes the primitive unit
minor, the complete bounded search for perpendicular roots, the attained
multisection degree, and a finite rank-three spanning-Gram box.

The exact conclusion remains:

\[
 \boxed{
 \text{fixed physical cusp data do not bound the host, while a bounded
 full-spanning physical root packet would.}}
\]

No finite-volume counterfamily, noncompact finite-host theorem, raw-carrier
return, or Property 4 result is asserted.
