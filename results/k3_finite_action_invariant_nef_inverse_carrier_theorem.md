# Inverse-carrier theorem for finite K3 actions and quotient nef cones

**Status (2026-08-08): frozen after hostile audit.**  The audit checked the
period/action quantifiers, complete-root separation, quotient descent,
geometric-versus-integral Coxeter gates, the scoped Meyer obstruction, and
the labeled moduli quotient.  The Coxeter and Vinberg conclusions below are
deliberately independent of the surface-existence equivalence.

The supporting calibration replay is

```text
PYTHONPATH=src:. pytest -q \
  tests/test_k3_quotient_nef_slices.py \
  tests/test_k3_nef_promotion_gate_h10_h11.py \
  tests/test_roulleau_r14_c2d_special_k3_order4_lift.py \
  tests/test_a18_c18_special_k3_c4xc2_lift.py
```

and currently gives `21 passed`.

## Result

Finite invariant sections of complete K3 nef cones form a second exact
inverse-carrier class.  The surface-existence part is an equivariant period
problem; the section-completeness part is an exact root-separation problem;
and the Coxeter and volume conclusions remain independent finite gates.

Let $N$ be an even lattice of signature $(1,r-1)$, $3\leq r\leq20$, and
let

\[
 \rho:G\longrightarrow O^+(N)
\]

be an action of a finite group.  Put $F=N^G$ and assume that $F$ has
signature $(1,d)$, $d\geq2$.  An **equivariant complete-root datum** consists
of the following additional information.

1. A full-dimensional pointed rational polyhedral cone
   $Q\subset\overline{{\cal C}_F^+}$ with $Q^\circ$ in the positive cone.
2. For every facet of $Q$, a primitive root $\delta_a\in N$ with
   $\delta_a^2=-2$ whose restricted hyperplane is that facet.
3. The complete-root separation condition

   \[
    \delta^\perp\cap Q^\circ=\varnothing
    \quad\text{for every }\delta\in N,\ \delta^2=-2.
    \tag{1}
   \]

4. A primitive embedding $i:N\hookrightarrow\Lambda_{K3}$ and an
   orientation-compatible extension

   \[
    \widetilde\rho:G\longrightarrow O(\Lambda_{K3})
   \]

   preserving $i(N)$ and restricting there to $i\rho i^{-1}$.
5. A character $\chi:G\to\mathbf C^\times$.  Define its exact-Picard
   period locus by

   \[
   \begin{aligned}
    \Omega_{\mathcal D}^{\circ}={[\omega]\in
      \mathbf P(T_{\mathbf C,\chi}):{}&
      \omega^2=0,\ (\omega,\bar\omega)>0,\\
      &\omega^\perp\cap\Lambda_{K3}=i(N)},
      \qquad T=i(N)^\perp,
   \end{aligned}
   \tag{2}
   \]

   Here $T_{\mathbf C,\chi}$ is the simultaneous subspace on which every
   $g\in G$ acts as multiplication by $\chi(g)$.

Then

\[
 \boxed{
 \Omega_{\mathcal D}^{\circ}\ne\varnothing
 \quad\Longleftrightarrow\quad
 \begin{array}{c}
 \text{there is a projective K3 surface }X\text{ and a full marking }\\
 \eta:H^2(X,\mathbf Z)\xrightarrow{\sim}\Lambda_{K3},
 \text{ and }a:G\to\operatorname{Aut}(X),\\
 \eta a(g)^*\eta^{-1}=\widetilde\rho(g)\text{ for every }g\in G,\\
 \eta(H^{2,0}(X))\subset T_{\mathbf C,\chi},\\
 \eta(\operatorname{NS}(X))=i(N),\qquad
 \eta(\operatorname{Nef}(X))\cap i(F_{\mathbf R})=i(Q).
 \end{array}}
 \tag{3}
\]

Every marked projective K3 surface with a finite action and a
finite polyhedral invariant nef section supplies data satisfying (1)--(2).
Thus nonemptiness of (2), together with the finite lattice and complete-root
gates above, is necessary and sufficient.  If a faithful surface action is
required, require $\widetilde\rho$ to be faithful; K3 automorphisms act
faithfully on integral second cohomology.

For the normal projective quotient $\pi:X\to Y=X/G$, pullback gives

\[
 \boxed{
 \pi^*N^1(Y)_{\mathbf R}=F_{\mathbf R},\qquad
 \pi^*\operatorname{Nef}(Y)=Q,
 }
 \tag{4}
\]

and

\[
 (\pi^*D)^2=|\operatorname{im}(a)|D^2.
\]

For a faithful action this factor is $|G|$.  Hence the invariant section is
the actual nef cone of an algebraic quotient surface, up to a harmless global
metric rescaling; it is not merely an artificial real slice.

## Why the criterion is exact

Choose a regular rational point $x\in Q^\circ$.  Condition (1) says that
$x$ lies in a unique chamber $D$ of the complete $(-2)$-root arrangement of
$N$, that all of $Q$ lies in $\overline D$, and that no root hyperplane
contains all of $F$.  The prescribed facet roots prevent
$\overline D\cap F_{\mathbf R}$ from extending past any facet of $Q$.
Therefore

\[
 \overline D\cap F_{\mathbf R}=Q.
 \tag{5}
\]

Because $x$ is fixed by $G$ and $D$ is the unique root chamber containing
$x$, the chamber $D$ is $G$-stable.  A period in (2) produces a marked
projective K3 surface with Picard lattice exactly $N$.  The refined period
map resolves the Weyl ambiguity by choosing the marked Kähler chamber to be
$D$.  The extended action is then a Hodge action preserving the Kähler
class $x$; strong global Torelli realizes it by automorphisms.  Equation (5)
gives (3).

For necessity, mark an actual action on $X$.  Its cohomological action gives
$\widetilde\rho$, its action on $H^{2,0}(X)$ gives $\chi$, and its period
lies in (2).  The K3 nef-wall theorem gives the facet roots.  An invariant
ample class rules out a root orthogonal to all of $F$, while the nef chamber
rules out a root hyperplane cutting $Q^\circ$.

If $v_1,\ldots,v_k$ are the complete extreme rays of $Q$, condition (1) is
equivalently the exact sign test

\[
 \{(v_j,\delta):1\leq j\leq k\}
 \quad\text{has one weak sign and is not identically zero}
 \tag{6}
\]

for every square-minus-two root $\delta$.  A complete root enumeration,
an exact stopping theorem, or a published complete chamber is still needed
to make the universal quantifier in (6) effective.  A finite sample of roots
does not prove (1).

## Folding and the Coxeter gate

On $F_{\mathbf R}$ all inequalities in one root orbit coincide.  The
orthogonal projection of a root is proportional to its orbit sum

\[
 R_{\mathcal O}=\sum_{\delta\in\mathcal O}\delta\in F.
 \tag{7}
\]

After discarding orbit inequalities whose hyperplanes miss hyperbolic space,
removing redundancies, and normalizing the remaining sums, the
projectivization of $Q$ is a geometric hyperbolic Coxeter polyhedron if and
only if:

1. the folded Gram matrix has Lorentzian inertia $(1,d)$ on its span;
2. every normalized pair is a Coxeter angle, parallel pair, or allowed
   ultraparallel pair;
3. the oriented halfspaces have a common hyperbolic interior; and
4. every displayed inequality is an actual facet.

There is a separate arithmetic gate:

5. after primitive integral normalization, every folded reflection is
   integral on the stated fixed or quotient divisor lattice.

Items 1--4 define the geometric Coxeter property.  Item 5 upgrades it to an
integral-reflection Coxeter polyhedron in the arithmetic group of the stated
lattice.  It is not part of the definition of a geometric Coxeter
polyhedron and does not follow from quotient geometry.

None of these assertions follows merely from the existence of the quotient
surface.  In particular, a rational polyhedral nef cone need not be Coxeter.

If the Coxeter gate passes, finite volume may be certified by Vinberg's
literal edge-endpoint and complete parabolic-neighborhood criterion or by an
independent exact cone computation.  There is also a structural route: since
`Q` is the complete invariant-nef section and contains an invariant ample
class, finite polyhedrality of `Q` implies finite hyperbolic volume whenever
the relative hyperbolic dimension is at least two.  Polyhedrality of an
arbitrary divisor/halfspace subsystem does not imply this, and in relative
dimension one an ideal endpoint can still give infinite length.  Literal
Vinberg or exact-cone replay therefore remains an independent wall-level
certificate even when the structural route applies.

There is also a sharp compactness restriction for the subclass passing item
5.  If the folded reflections are integral on $F$ and their chamber has
finite volume, their group has finite index in the arithmetic lattice
$O^+(F)$.  For $d\geq4$, Meyer makes $F_{\mathbf Q}$ isotropic, so this
arithmetic lattice and every finite-index subgroup are nonuniform.  Therefore

\[
 \boxed{d\geq4\quad\Longrightarrow\quad
 \begin{array}{c}
 \text{a finite-volume integral-reflection rational}\\
 \text{invariant-nef fold is noncompact}.
 \end{array}}
 \tag{8}
\]

Consequently the integral-reflection subclass of finite K3 quotient-nef
folds cannot produce the desired compact $\mathbb H^9$ example.  This
argument does not exclude a geometric Coxeter quotient fold whose facet
reflections fail to preserve $F$.  A nonrational totally-real eigenspace is
one alternative route around the arithmetic obstruction.

## The moduli problem and its finite gates

For a fixed labeled lattice action and embedding orbit, let
$\Gamma_{i,\widetilde\rho,C}$ be the subgroup induced by K3-lattice
isometries which fix $i(N)$ pointwise, centralize
$\widetilde\rho(G)$, and preserve the chosen component of (2) and the
discrete Kähler-chamber choice.  The labeled, pointwise-$N$-marked carrier
locus is

\[
 \Gamma_{i,\widetilde\rho,C}\backslash
 \Omega_{\mathcal D}^{\circ,+}.
\]

For an unlabelled action the corresponding normalizer replaces the
centralizer.  The locus's nonemptiness is decided in the following order.

1. Enumerate primitive embeddings of $N$ by discriminant forms.
2. Enumerate extensions of $\rho(G)$ through the gluing to
   $\Lambda_{K3}$, including their action on the orthogonal complement.
3. For each character $\chi$, decide whether the simultaneous period
   eigenspace contains a positive isotropic line avoiding every extra
   Noether--Lefschetz hyperplane.
4. Prove the complete-root separation test (6) and identify the irredundant
   orbit walls.
5. Apply the folded Coxeter and Vinberg gates.

Steps 1--3 are the softer equivariant surface-existence problem.  Steps 4--5
are the rigid wall and hyperbolic-polyhedron problem.  Classification or
emptiness of the period strata transfers to the corresponding **marked
finite-action carrier class**.  It becomes a nonexistence theorem for a
larger abstract Coxeter class only after a separate universality theorem
forces every member of that class to admit one of these lattice actions.

The passage to unmarked abstract polytopes is generally many-to-one:
different periods, primitive embeddings, cohomological extensions, groups,
or integral lattices can yield the same folded Coxeter diagram.

## Calibrations and strict limits

The theorem recovers the two complete quotient-nef calibrations already in
the repository:

- the rank-18 K3 involution gives the ten-facet finite-volume noncompact
  simplex in $\mathbb H^9$;
- a transposition in the rank-19 K3 $S_3$ action gives the thirteen-facet
  finite-volume noncompact polyhedron in $\mathbb H^{10}$.

The following near misses show why every hypothesis is present.

- A graph or Picard-lattice symmetry need not extend to a K3 Hodge action.
- A same-order lift can fail even when a higher-order transcendental lift
  succeeds.
- Coinvariant $(-2)$-roots obstruct an invariant Kähler class.
- Hidden ambient roots can cut a proposed fixed section.
- Primitive square-minus-four divisor reflections are not K3 nef walls.
- A complete mixed-root reflection chamber preserved by an automorphism is
  not automatically the nef cone.
- A quotient-nef cone can fail the Coxeter-angle or finite-volume gate.

Concrete certificates for these boundaries are recorded in the supporting
results below.

## Supporting results

- [`k3_strong_inverse_carrier_transfer_theorem.md`](k3_strong_inverse_carrier_transfer_theorem.md)
  supplies the complete K3 root-chamber theorem and the period/nef dictionary.
- [`k3_nef_promotion_gate_h10_h11.md`](k3_nef_promotion_gate_h10_h11.md)
  proves the exact root-separation criterion and gives hidden-root failures.
- [`k3_quotient_nef_slices.md`](k3_quotient_nef_slices.md) proves finite
  quotient descent, orbit folding, and the two calibrated quotient cones.
- [`roulleau_r14_c2d_special_k3_order4_lift.md`](roulleau_r14_c2d_special_k3_order4_lift.md)
  and [`a18_c18_special_k3_c4xc2_lift.md`](a18_c18_special_k3_c4xc2_lift.md)
  give exact equivariant gluing, period, rootless-coinvariant, Kähler, and
  Torelli realizations.
- [`vinberg_divisor_volume_gate.md`](vinberg_divisor_volume_gate.md) records
  the independent terminal finite-volume criterion.
