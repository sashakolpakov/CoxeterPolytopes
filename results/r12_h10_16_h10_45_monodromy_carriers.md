# Two explicit algebraic carriers for the H10/16 and H10/45 polytopes

## Positive theorem

There are two explicit algebraic carriers for already certified hyperbolic
Coxeter polytopes in dimension ten.

1. A four-modulus elliptic K3 family with very general
   \(NS=U\oplus D_8(-1)\oplus A_1(-1)^2\) has a transverse enhancement
   meridian whose parallel transport is the \(D_8\) spinor involution.  Its
   invariant nef section is the finite-volume noncompact H10/16 polytope,
   with 16 walls.
2. A four-modulus elliptic K3 family with very general
   \(NS=U\oplus D_4(-1)^2\oplus A_1(-1)^2\) has a transverse enhancement
   meridian whose parallel transport is a \(D_4\) outer transposition.  Its
   invariant nef section is the finite-volume noncompact H10/45 polytope,
   with 45 walls.

This is an equation-to-divisor-monodromy-to-Coxeter realization.  Parallel
transport in a family is the algebraic mechanism; it is not an automorphism
of a single K3 surface.  It therefore supplies a positive
algebro-geometric answer even though the previously tested Nikulin
single-fiber automorphism lifts fail.  There is no new abstract-polytope
claim: the new item is the explicit algebraic carrier of each known
certified chamber.

Concretely, the exceptional curves of the Weierstrass ADE resolutions form
an unlabelled divisor local system over the punctured parameter disk.  A
simultaneous resolution with individually labelled terminal curves requires
the corresponding root cover.  Going once around the puncture therefore
permutes the divisor classes even though the coefficient tuple returns to
itself.  No claim about the transcendental monodromy, or about a global
regular automorphism of the returned fiber, is needed.

## The 16-wall carrier

Put \(s=t^2\) and consider

\[
y^2=x^3+P(t^2)x+tQ(t^2),
\]

where

\[
\begin{aligned}
P(s)&=-3s^3+\alpha s^2+\beta s-3r^2+3-\alpha-\beta,\\
Q(s)&=2s^4-\alpha s^3+\gamma s^2+q_1s+q_0,\\
q_1&=9r-2\alpha r-\beta r-r^3-8+3\alpha-2\gamma,\\
q_0&=2r^3-2+\alpha-\gamma-q_1.
\end{aligned}
\]

The very general fiber has configuration \(I_4^*+2I_2+10I_1\), Mordell--Weil
rank zero, and no Mordell--Weil torsion.  The height bound
\(h(P)\geq4-3=1\) excludes a nonzero torsion section, so the trivial lattice
is primitive and

\[
NS=U\oplus D_8(-1)\oplus A_1(-1)^2
\]

exactly, not merely as a finite-index sublattice.

At the clean boundary point

\[
(r,\alpha,\beta,\gamma)=(1,2,1,-2/3)
\]

take the explicit transverse slice

\[
r=1,\qquad \alpha=2,\qquad \beta=1,\qquad \gamma=-2/3+z
\]

and the meridian \(z=\varepsilon e^{i\theta}\).  In the literal third
blow-up chart of the \(I_4^*\) fiber, the two terminal components are labelled
by the roots of

\[
3q^2+2q+z+\frac13=0,
\qquad
q_\pm=\frac{-2\pm\sqrt{-12z}}6.
\]

Thus one turn around \(z=0\) exchanges \(q_+\) and \(q_-\), hence exchanges
exactly the two spinor components of \(D_8\).  The common six-component
chain, \(U\), and both \(A_1\) components remain fixed.  The central fiber
enhances \(I_4^*/D_8\) to \(I_6^*/D_{10}\).

This is a genuinely transverse and isolated local event.  At the displayed
boundary point,

\[
\begin{aligned}
P(s)&=-3s^3+2s^2+s-3,\\
Q(s)&=2s^4-2s^3-\frac23s^2+\frac73s+\frac13,
\end{aligned}
\]

and the residual discriminant factor is

\[
-112s^4+16s^3+120s^2-105s-108.
\]

Its discriminant is \(2220142134528\), its resultant with \(P\) is
\(-665523\), and it is nonzero at both marked finite fibers.  These exact
checks rule out a simultaneous accidental collision.  A base point may be
chosen in an arbitrarily small generic deformation outside the countable
Noether--Lefschetz divisors; the local braid class persists and the nearby
Picard lattice is then the exact very-general lattice above.

On the full rank-12 lattice the transport is

\[
\operatorname{id}_U\oplus\operatorname{outer}(D_8)
\oplus\operatorname{id}_{A_1^2}.
\]

The certificate gives an explicit determinant \(-1\) change of basis from
this standard model to the archived Roulleau basis.  It carries the two
spinor components to archived walls 10 and 11 and conjugates the entire
integral action, not only its rational representation.  The fixed lattice
has rank/signature/determinant

\[
(11,(1,10),16),
\]

Smith factors \(1^8,2,2,4\), and coinvariant
\(\langle-4\rangle\) of divisibility two.  The induced discriminant-form
map is checked equivariantly on all 16 elements of \((\mathbb Z/2)^4\).

## The 45-wall carrier

Now consider

\[
y^2=x^3+t^2P(t^2)x+t^3Q(t^2),
\]

with

\[
\begin{aligned}
P(s)&=-3+\alpha(s-1)+\beta(s-1)^2,\\
Q(s)&=2-\alpha(s-1)+\gamma(s-1)^2+\delta(s-1)^3.
\end{aligned}
\]

The very general fiber has configuration \(2I_0^*+2I_2+8I_1\), Mordell--Weil
rank and torsion zero, and exact lattice

\[
NS=U\oplus D_4(-1)^2\oplus A_1(-1)^2.
\]

Here the height estimate is \(h(P)\geq4-(1+1+1/2+1/2)=1\), again proving
primitive saturation.

Use the boundary point

\[
(\alpha,\beta,\gamma,\delta)=(1,1,1,2)
\]

and transverse slice \(\delta=2+z\), with the other three parameters fixed.
The first \(I_0^*\) Tate cubic is exactly

\[
\rho^3-3\rho-z+2,
\]

whose discriminant is \(27z(4-z)\).  At \(z=0\),

\[
\rho^3-3\rho+2=(\rho-1)^2(\rho+2).
\]

A positive meridian swaps the two roots specializing to \(1\), while the
root specializing to \(-2\) is fixed.  Since these roots literally label
the three outer components in the first blow-up chart, parallel transport
exchanges precisely two outer components of the first \(D_4\).  It fixes the
third outer and central components, all of the second \(D_4\), \(U\), and
both \(A_1\)'s.  The central enhancement is \(I_0^*/D_4\) to \(IV^*/E_6\).

The residual discriminant is

\[
s(112s^3-328s^2+207s+216).
\]

The other cubic has discriminant \(-14242881024\), its resultant with
\(P\) is \(-275427\), and the Tate cubic of the second \(D_4\) has
discriminant \(-112\).  Thus the intended component collision is clean.

The full Picard transport is

\[
\operatorname{id}_U\oplus\operatorname{outer}_{(12)}(D_4)
\oplus\operatorname{id}_{D_4}\oplus\operatorname{id}_{A_1^2}.
\]

Another explicit determinant \(-1\) basis change carries the exchanged
outer pair to archived walls 2 and 3 and conjugates this full integral
action to the archived `S6` transposition.  Its fixed lattice has

\[
(\operatorname{rank},\operatorname{signature},|\det|)
=(11,(1,10),64),
\]

Smith factors \(1^6,2,2,2,2,4\), and coinvariant
\(\langle-4\rangle\) of divisibility two.  The exact discriminant action is
transported equivariantly and checked on all 64 elements of
\((\mathbb Z/2)^6\).

## Vinberg finite-volume gates

The algebraic identification inherits the already complete chamber and
Vinberg audits; finite volume is not inferred from a truncated diagram.

| carrier | facets | finite vertices | cusps | elliptic edges | affine-neighbourhood audits | bad audits |
|---|---:|---:|---:|---:|---:|---:|
| H10/16 | 16 | 70 | 6 | 426 | 14 | 0 |
| H10/45 | 45 | 1,048 | 93 | 7,332 | 302 | 0 |

For H10/16, every one of the 426 elliptic rank-nine edge diagrams has
exactly two finite or ideal endpoint extensions, and the canonical complete
audit digest is

    900df6d000ca91669cea117474cda559e4188600fb2744bafa0f69fd1541c226

For H10/45, the endpoint-type counts are 578 of type `(0,2)`, 3,028 of
type `(1,1)`, and 3,726 of type `(2,0)`.  All 302 connected affine
zero-neighbourhoods pass.  Hence both chambers are rigorously
finite-volume noncompact hyperbolic Coxeter polytopes.

## What has and has not been obtained

The positive advance is two explicit divisor-lattice carriers:

\[
\text{Weierstrass equation}
\longrightarrow
\text{resolved fiber components}
\longrightarrow
\text{component braid}
\longrightarrow
O(NS)
\longrightarrow
\text{fixed Coxeter chamber}.
\]

This realizes the sought algebraic-object mechanism without requiring the
monodromy to be a regular automorphism of one fiber.  It does not create a
new abstract Coxeter diagram, and it does not evade the separate
single-fiber Nikulin obstructions.  Those statements are compatible: the
same divisor action can exist as monodromy of an algebraic family even when
it cannot extend as the proposed involution of an individual K3 surface.

The proof also isolates a reusable positive criterion.  An algebraic family
automatically supplies a genuine hyperbolic Coxeter section once one has:

1. an exact, primitively saturated generic Néron--Severi lattice and a
   complete finite wall description of its nef chamber;
2. a resolved-component local system whose meridian gives an integral wall
   permutation;
3. a fixed lattice of signature \((1,n)\); and
4. the full Vinberg edge-endpoint and affine-neighbourhood audits for the
   fixed chamber.

The first two gates are algebro-geometric, the third is the Hodge-index
hyperbolicity gate, and the fourth distinguishes finite volume from an
infinite-volume chamber.  The two carriers above pass all four gates.

## Reproducibility

- certificate generator:
  [`certify_r12_h10_16_h10_45_monodromy_carriers.py`](../scripts/certify_r12_h10_16_h10_45_monodromy_carriers.py)
- exact archive:
  [`r12_h10_16_h10_45_monodromy_carriers.json`](../data/enumerations/r12_h10_16_h10_45_monodromy_carriers.json)
- regression tests:
  [`test_r12_h10_16_h10_45_monodromy_carriers.py`](../tests/test_r12_h10_16_h10_45_monodromy_carriers.py)

Replay with:

```text
PYTHONPATH=.:src python3 scripts/certify_r12_h10_16_h10_45_monodromy_carriers.py --check-archive
```

The pinned certificate payload SHA-256 is

    ce457d3a0b1718e73391eb68ec7f7ca0aadb8f2e1f32c4b3eb12a470c3aef9c0
