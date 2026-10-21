# TOOT: an equation-level K3 family for the `A15_11` `H10/21` carrier

## Positive result

The `A15_11_1_C02` finite-volume `H10/21` chamber now has an explicit
nine-parameter algebraic K3 realization, not only an abstract period point.
Put

\[
\begin{aligned}
A(t)&=-3+a_1t+a_2t^2+a_3t^3+a_4t^4,\\
B(t)&=2+(1-a_1)t+b_2t^2+b_3t^3+b_4t^4+b_5t^5+b_6t^6,
\end{aligned}
\]

and let (X_{\mathbf a,\mathbf b}) be the minimal resolution of

\[
 y^2=x^3+A(s^2)x+B(s^2).
\]

Let (D(t)=4A(t)^3+27B(t)^2=tQ(t)).  On the nonempty Zariski open

\[
 \operatorname{Disc}(Q)\ne0,
 \qquad
 4a_4^3+27b_6^2\ne0,
\]

the following statements hold.

1. (X_{\mathbf a,\mathbf b}) is an elliptic K3 surface with fibers
   (I_2+22I_1).
2. The deck involution
   (iota(s,x,y)=(-s,x,y)) is non-symplectic and fixes one elliptic
   curve and one rational curve.
3. Its invariant integral divisor lattice is exactly
   
   \[
   H^2(X,\mathbb Z)^\iota
   \cong U\oplus E_8(-2)\oplus\langle-2\rangle=:F.
   \]
4. For a very general parameter,
   (NS(X)=F).  Under the literal standard marking already archived for
   `A15_11`, all 21 wall reflections therefore act integrally on the
   **full** Néron--Severi lattice.

This is a positive equation-to-divisor-to-Coxeter construction.  It is not
a new abstract Coxeter polytope, and it does not claim that the 21 walls are
the full nef cone.

## Why the equation works

Before base change, consider the rational elliptic surface

\[
 R:\quad y^2=x^3+A(t)x+B(t).
\]

The normalization forces

\[
 D(0)=0,\qquad D'(0)=108,\qquad c_4(0)=144,
\]

so (t=0) is an (I_1) fiber.  The two open conditions say that the
other eleven zeros of (D) are simple and that the fiber at infinity is
smooth.  Thus (R) has twelve (I_1) fibers and no reducible fibers.

The map (t=s^2) is ramified over the marked (I_1) fiber and the smooth
fiber at infinity.  It changes the ramified (I_1) into (I_2), doubles
the other eleven (I_1)'s, and leaves the ramified smooth fiber smooth.
Consequently the Euler sum is

\[
 e(X)=e(I_2)+22e(I_1)=2+22=24.
\]

The degrees of (A(s^2)) and (B(s^2)) are at most 8 and 12, so the
minimal resolution is the corresponding elliptic K3.  Equivalently, this
is the resolved double cover of the blow-up of
(\mathbb P^2) defining (R), branched over a marked nodal cubic and a
marked smooth cubic.

There is a particularly small witness proving that the parameter open is
nonempty:

\[
 A(t)=t^4+t-3,
 \qquad
 B(t)=t^6+2,
\]

for which

\[
 Q(t)=31t^{11}+12t^8-36t^7+120t^5-72t^4+108t^3
      +4t^2-36t+108.
\]

The certificate checks exactly that this polynomial is squarefree, that it
has no common root with (A), and that the discriminant at infinity is
31.  This witness proves nonemptiness; it is not being advertised as an
individually Picard-certified special member.

## The integral divisor-lattice identification

This is more concrete than identifying the lattice only by
((r,a,\delta)).

- Let (f) be the K3 fiber and (O) the zero section.  The classes
  (f,O+f) have Gram (U).
- Since (R) has no reducible fibers, its Mordell--Weil height lattice is
  (E_8).  Pulling eight Shioda classes through the degree-two base change
  doubles every intersection pairing, producing (E_8(-2)) in
  (NS(X)).
- The nonidentity component (C) of the ramified (I_2) fiber has
  (C^2=-2) and is orthogonal to the preceding corrected classes.

These eleven invariant divisor classes therefore span

\[
 L=U\oplus E_8(-2)\oplus\langle-2\rangle,
 \qquad |\det L|=2^9.
\]

On the other hand, the fixed locus has (g=1) and (k=1).  Nikulin's
relations give

\[
 g=\frac{22-r-a}{2}=1,
 \qquad
 k=\frac{r-a}{2}=1,
\]

hence (r=11,a=9); the (I_0+I_1) base-change case has
(\delta=1).  Thus the full invariant lattice has determinant (2^9).
The inclusion (L\subset H^2(X,\mathbb Z)^\iota) has equal rank and equal
determinant, so its index is one.  This proves the claimed integral
identification without leaving a possible saturation gap.

Locally at the ramified node the model is (uv=s^2).  Its resolution is
the (I_2) fiber.  The strict transform of the old nodal fiber is the one
fixed rational curve; the exceptional component is invariant but not
fixed pointwise.  This is why the fixed-locus count is exactly one rational
curve, not two.

## Why the very general full Picard lattice is exactly (F)

The coefficient normalization leaves nine parameters.  It has used
`PGL(2)` to put the ordered branch points at (0,\infty), Weierstrass
homothety to impose ((A(0),B(0))=(-3,2)), and the remaining base scaling
to impose (D'(0)=108).  The nine parameters are recovered directly from
the derivatives of (A) and (B); the archived recovery Jacobian is the
identity matrix.

The moduli of 2-elementary K3 pairs of type ((11,9,1)) also has dimension
(20-r=9).  This dimension equality is accompanied by a geometric
recovery, rather than assumed from a parameter count: the unique fixed
elliptic curve gives the type-2 elliptic fibration; quotienting by
(iota), then contracting the distinguished exceptional curve over the
ramified node, recovers (R) and its ordered (I_1)/smooth branch fibers.
Thus the normalized family has zero-dimensional generic fibers over the
moduli of pairs and dominates its component.

Any additional algebraic class is anti-invariant and imposes a
Noether--Lefschetz hyperplane.  Outside their countable union,

\[
 NS(X)=H^2(X,\mathbb Z)^\iota=F.
\]

## Coxeter consequence and boundary

The pinned carrier certificate uses exactly the same standard Gram for
(F), entry by entry.  Transporting its archived walls gives 21 integral
reflections on the full (NS(X)) of a very general member.  The already
completed Vinberg audit gives 607 finite vertices, 50 ideal vertices,
4,159 elliptic edges, and no bad edge or parabolic neighborhood.  Hence the
chamber is finite-volume and noncompact.

What is established:

- an explicit algebraic family, with a checked nonempty parameter open;
- exact full (NS) for the very general member;
- an integral geometric basis for (F), not merely a period-map label;
- all 21 certified Coxeter reflections on that full divisor lattice.

What is not established:

- a new abstract polytope;
- that the 21-wall chamber is the full nef cone;
- that its ten square-(-4) walls are irreducible curve walls; or
- that the lattice reflections are automorphisms of the surface.

The deck involution supplies the algebraic construction and acts trivially
on the very general (NS); it is not being used as a substitute for the
Coxeter reflections.

## Lineage and reproduction

The quadratic base-change construction, fixed-fiber table, and invariant
lattice framework are in Garbagnati--Salgado,
[Linear systems on rational elliptic surfaces and elliptic fibrations on K3
surfaces](https://arxiv.org/abs/1703.02783).  Their later primary
classification identifies the (k=1) case, the (I_0+I_1) branch data,
and ((r,a,\delta)=(11,9,1)):
[Elliptic fibrations on K3 surfaces with a non-symplectic involution fixing
rational curves and a curve of positive genus](https://arxiv.org/abs/1806.03097).
Shioda supplies the Mordell--Weil height-lattice formalism.  The K3 class
and abstract Coxeter polytope are therefore known; the positive project
addition is the exact normalized equation, full-Picard proof, and literal
equation-to-`A15_11` carrier bridge.

Permanent files:

- [certificate](../data/enumerations/a15_11_ue8m2a1_basechange_family.json)
- [certifier](../scripts/certify_a15_11_ue8m2a1_basechange_family.py)
- [focused tests](../tests/test_a15_11_ue8m2a1_basechange_family.py)
- [underlying chamber carrier](a15_11_h10_21_x2_k3_carrier.md)

Replay from the repository root:

```bash
PYTHONPATH=.:src python3 scripts/certify_a15_11_ue8m2a1_basechange_family.py --check-archive
PYTHONPATH=.:src pytest -q tests/test_a15_11_ue8m2a1_basechange_family.py
```
