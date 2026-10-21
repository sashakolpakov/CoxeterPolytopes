# Mixed K3/Fourier--Mukai carriers: the sharp capacity and lift obstruction

**Status (2026-08-12).**  This note gives an unconditional obstruction theorem
for the proposed universal mixed K3 carrier.  It does **not** prove that every
reflection-maximal arithmetic chamber has such a carrier.  In fact, essential
surjectivity of this species would already prove the new arithmetic inequality

\[
                  [k:\mathbf Q](n+1)\le 23.                 \tag{0.1}
\]

No Vinberg, Nikulin, Esselmann, or ABSW dimension/finiteness theorem is used
below.  The obstruction is Hodge index plus elementary lattice theory.

The main correction to the earlier rank-only discussion is that rank 24 is not
available to a Lorentzian full-trace block.  The numerical Mukai lattice of a
supersingular K3 has rank 24, but signature \((2,22)\).  A trace block of rank
24 has signature \((1,23)\), so it cannot embed.  The sharp mixed-characteristic
ceiling is 23, not 24.  In characteristic zero it is 21.

## 1. Exact carrier contract

Let \(P\subset\mathbf H^n\) be an admissible arithmetic Coxeter datum over a
totally real field \(k\) of degree \(d\).  Let \((M,q)\) be a rank-\(n+1\)
integral order module, with

\[
 \operatorname {sig}(q_{\sigma _0})=(1,n),\qquad
 q_\sigma<0\quad(\sigma\ne\sigma _0).
\]

After a positive codifferent normalization, its full restriction-of-scalars
trace lattice is

\[
 T_P=(M_{\mathbf Z},\operatorname {Tr}_{k/\mathbf Q}(\beta q)),
 \qquad R:=\operatorname {rk}T_P=d(n+1),
 \qquad \operatorname {sig}T_P=(1,R-1).                 \tag{1.1}
\]

A **full-trace mixed Mukai carrier** means, at minimum:

1. a K3 surface, twisted K3 surface, or genuine K3 category \(\mathcal X\),
   possibly in positive characteristic;
2. an injective nondegenerate similarity

   \[
       \iota:T_P\otimes\mathbf Q\hookrightarrow
       K_{\rm num}(\mathcal X)\otimes\mathbf Q           \tag{1.2}
   \]

   by one positive rational scale, intertwining the returned order action;
3. actual Fourier--Mukai equivalences (or graph kernels of automorphisms) whose
   cohomological actions restrict to the side reflections, not arbitrary Chow
   matrices; and
4. an exact extraction rule returning every oriented side and no hidden wall.

There are two inequivalent meanings of item 3 which must not be conflated.
The weak one asks for one kernel \(\Phi_i\) over each side with
\(\Phi_i^H|_{T_P}=s_i\).  The strong one asks for a homomorphism

\[
             W(P)\longrightarrow \operatorname {Aut}_{\rm FM}(\mathcal X)
                                                               \tag{1.3}
\]

whose value on the \(i\)-th Coxeter generator has the prescribed action.  The
strong version includes the functor-level relations
\(\Phi_i^2\simeq\mathrm {Id}\) and
\((\Phi_i\Phi_j)^{m_{ij}}\simeq\mathrm {Id}\).  Derived Torelli supplies only
the weak version.  Section 4 below records the resulting extra gate.

The phrase *full trace block* is important.  A quotient, a non-isometric
semialgebraic encoding, or arbitrary divisor-product correspondences falls back
to the universal weak Chow carrier and is outside (1.2).

## 2. The sharp signature theorem

### Theorem 2.1 (Mukai capacity)

Suppose (1.2) exists.

1. If \(\mathcal X\) is a complex projective K3, a finite-height K3, a twisted
   K3 of either kind, or a characteristic-zero K3 category, then

   \[
                       \boxed{R\le21}.                   \tag{2.1}
   \]

2. If \(\mathcal X\) is a supersingular K3 (or a supersingular K3 category
   with maximal algebraic Mukai lattice), then

   \[
                       \boxed{R\le23}.                   \tag{2.2}
   \]

Consequently a mixed category allowing all characteristics still satisfies

\[
                       \boxed{d(n+1)\le23}.              \tag{2.3}
\]

#### Proof

For a projective K3 surface,

\[
 \widetilde {\operatorname {NS}}(X)
   =H^0(X,\mathbf Z)\oplus\operatorname {NS}(X)\oplus H^4(X,\mathbf Z)
   \cong U\oplus\operatorname {NS}(X).
\]

If \(\rho=\operatorname {rk}\operatorname {NS}(X)\), its signature is
\((2,\rho)\).  In characteristic zero and at finite height, \(\rho\le20\).
A subspace of signature \((1,R-1)\) therefore needs \(R-1\le20\), proving
(2.1).  For a supersingular K3, \(\rho=22\), so the same argument gives
\(R-1\le22\), proving (2.2).

A twisted Mukai lattice is obtained by a rational \(B\)-field isometry; its
rank and real inertia are unchanged.  A genuine characteristic-zero K3
category has a rank-24 K3 Mukai Hodge lattice of signature \((4,20)\); its
algebraic part is orthogonal to the positive period plane and hence has
signature at most \((2,20)\).  The same calculation gives (2.1).  The
supersingular analogue has algebraic signature at most \((2,22)\).
\(\square\)

### Corollary 2.2 (rank 24 is impossible)

No admissible full trace lattice of rank 24 embeds in a rank-24
supersingular numerical Mukai lattice: their signatures are respectively

\[
                         (1,23)\quad\hbox{and}\quad(2,22).
\]

This remains true after:

- a common positive rescaling (Sylvester inertia is unchanged);
- a nonprimitive or finite-index embedding (the rational subspace is
  unchanged);
- replacing physical facets by projected ambient roots (the trace block is
  unchanged);
- a Brauer twist or rational \(B\)-field transform (it is an isometry); or
- replacing a geometric K3 by a genuine noncommutative K3 category (the
  K3-type Mukai signature is unchanged).

Allowing arbitrary 2-Calabi--Yau categories with unbounded numerical rank
would evade the theorem, but would also discard the bounded-rank leverage of
the proposed K3 species.

## 3. Primitive embedding and discriminant gates

The inequality in Theorem 2.1 is necessary, not sufficient.  Let \(A\) be
the algebraic numerical Mukai lattice and let \(T\subset A\) be a primitive
copy of the trace block with orthogonal complement \(C\).  An integral
embedding is described by an isotropic gluing subgroup

\[
                  H\subset A_T\oplus A_C,               \tag{3.1}
\]

where \(A_L=L^\vee/L\) denotes a discriminant module.  A side reflection
\(s_i\in O(T)\) extends to \(A\) only if there is \(c_i\in O(C)\) for which

\[
                         (\bar s_i,\bar c_i)H=H.          \tag{3.2}
\]

Thus rank and signature do not construct the lattice, and an embedding of
quadratic spaces does not construct the side actions.  Scaling can clear
denominators, but it does not make (3.2) automatic; it usually enlarges the
discriminant modules which have to be glued.

There is a particularly transparent full-divisor obstruction at the
supersingular endpoint.  If a rank-22 trace lattice is required to be the
full N\'eron--Severi lattice up to finite index and common rational
similarity, then

\[
                  \det T_P\equiv -1
                  \pmod {\mathbf Q^{*2}}.                \tag{3.3}
\]

Indeed a supersingular K3 of Artin invariant \(\sigma\) has determinant
\(-p^{2\sigma}\); finite index changes determinant by a square, and scaling
a rank-22 form also changes it by a square.  Integrally, the ambient
discriminant group is \(p\)-elementary (type I at \(p=2\)), which is a still
stronger gluing restriction.  Equation (3.3) is not asserted when the block
is allowed to mix with the Mukai hyperbolic plane; in that weaker category
one must use the full gluing test (3.1).

## 4. Fourier--Mukai/Hodge lift gate

Even an isometric embedding satisfying (3.1) need not give actual kernels.
For a complex K3, write

\[
 N=\widetilde {\operatorname {NS}}(X),\qquad
 T_X=N^\perp\subset\widetilde H(X,\mathbf Z).
\]

Unimodularity gives an anti-isometry \(A_N\simeq A_{T_X}\).  A Mukai Hodge
isometry restricts to isometries \(g_N,g_T\) satisfying

\[
                \delta\,\bar g_N=\bar g_T\,\delta.       \tag{4.1}
\]

It is induced by a Fourier--Mukai equivalence only after the signed
(positive-four-plane orientation) condition is also imposed.

The transcendental Hodge structure is simple when \(N\) is the exact
algebraic lattice.  Its finite Hodge-isometry group is contained in the roots
of unity of its Zarhin endomorphism field and is therefore cyclic.  Hence:

### Proposition 4.1 (complex discriminant-image obstruction)

If the trace block is the full algebraic Mukai lattice (or if its
discriminant action injects equivariantly into \(A_N\)), a simultaneous
Fourier--Mukai lift of the side group forces its induced discriminant image
to be cyclic.

#### Proof

By (4.1) the discriminant image is conjugate to an image of the finite Hodge
isometry group of \(T_X\), which is cyclic.  The orientation condition can
shrink this image but cannot make it noncyclic. \(\square\)

This explains two calibration facts without conflating them.

- The rank-20 \(X_3\) packet can lift its long reflections over \(\mathbf C\):
  their discriminant action is the scalar \(-1\), matched by
  \(-\mathrm{id}\) on the transcendental lattice.
- For the Borcherds lattice \(U\oplus D_{20}(-1)\), the displayed reflection
  group induces \(GL(2,\mathbf F_2)\simeq S_3\) on the discriminant module.
  Thus its full action cannot be routed through a complex exact-algebraic
  K3 lattice.  Its characteristic-two supersingular realization is a
  genuinely different crystalline branch, not a complex Fourier--Mukai
  deformation.

In the supersingular branch, (4.1) is replaced by preservation of the
strictly characteristic crystalline period subspace.  An abstract isometry
of the numerical Mukai lattice need not preserve that subspace.  Ogus-type
crystalline Torelli, and its derived analogue when used, makes this an actual
extra condition.  The special Borcherds \(H^{21}\) surface passes it because
the 168 long actions are honest automorphisms and the 42 short actions are
spherical twists; maximality alone does not supply such a lift.

### Proposition 4.2 (a spherical-twist lift is not a Coxeter lift)

Let \(S\) be a 2-spherical object in the bounded derived category of a K3
surface or in a genuine K3 category.  Although the spherical twist \(T_S\)
acts on the numerical Mukai lattice as the reflection in \(v(S)\), it cannot
be the image of a Coxeter involution in (1.3).  More precisely,

\[
                  T_S(S)\simeq S[-1],\qquad
                  T_S^m(S)\simeq S[-m],                 \tag{4.2}
\]

so \(T_S\) has infinite order and in particular
\(T_S^2\not\simeq\mathrm {Id}\).

#### Proof

The first isomorphism is the defining calculation for the twist in a
2-spherical object.  Exact equivalences commute with shifts, so iteration
gives the second.  A 2-spherical object cannot be isomorphic to a nonzero
shift of itself, since its graded endomorphism space is
\(\mathbf k\oplus\mathbf k[-2]\) and is not periodic.  Hence no positive
power in (4.2) is the identity.
\(\square\)

Thus the familiar spherical-twist constructions give Artin/braid-type lifts
whose **numerical quotient** is the Coxeter action; they do not by themselves
split that quotient back to a categorical \(W(P)\)-action.  The
surjectivity in derived Torelli says that each signed Hodge isometry has a
Fourier--Mukai preimage.  It neither gives a group-theoretic section nor
preserves the order of a chosen isometry.

This distinction is not merely terminological.  For complex K3 surfaces of
Picard number one, [Fan--Lai](https://arxiv.org/abs/2302.12663) prove that an involutive autoequivalence cannot
act on the numerical Mukai lattice as reflection in a \((-2)\)-vector; such a
lift is necessarily of infinite order.  In arbitrary Picard rank and in the
supersingular branch a comparable no-splitting theorem is not presently
available here, so Proposition 4.2 is deliberately stated for the actual
spherical-twist lifts, not for every possible lift of the same numerical
reflection.

Consequently the existing \(X_3\) and Borcherds \(H^{21}\) packets establish
the weak kernel-by-kernel condition: the curve walls are represented by
spherical twists and the remaining displayed walls by other
Fourier--Mukai/graph kernels.  They do **not** establish (1.3).  If ``exact
action return'' means only recovery of the cohomological lattice action, they
remain valid controls.  If it means a strict categorical Coxeter action, a
new splitting condition must be added to the carrier contract.

## 5. Why reflection maximality does not prove essential surjectivity

Reflection maximality has the following scope: it excludes a proper discrete
overgroup generated by real hyperplane reflections in the same native
\(\mathbf H^n\).  None of the following follows formally from it.

1. **Capacity:** it does not provide (2.3).  A proof that every maximal input
   has a mixed K3 carrier would itself be a proof of (2.3), before any AG
   boundedness theorem is invoked.
2. **Local integral embedding:** it says nothing about rational Hasse
   invariants, discriminant forms, or the gluing subgroup (3.1).
3. **Kernel lift:** it says nothing about (4.1), orientation, or preservation
   of a crystalline characteristic subspace.  If a strict categorical action
   is required, it also says nothing about splitting the cohomological
   representation over the Coxeter group; Proposition 4.2 shows that the
   standard curve-wall lifts do not split it.
4. **No hidden walls:** an ambient Mukai root need not preserve the trace
   eigenspace.  Its projected hyperplane can cut the selected slice without
   defining a native integral reflection, so native maximality does not make
   it redundant.
5. **Decoration boundedness:** it does not bound polarizations, Hilbert
   polynomials of graphs, numerical kernel classes, Tor amplitude, or
   regularity.

The last point is visible even on one fixed K3.  Tensoring by \(L^m\) gives
Fourier--Mukai kernels

\[
             \mathcal O_\Delta\otimes p_2^*L^m
\]

with unbounded numerical class and graph/kernel Hilbert polynomial.  Thus
boundedness of the underlying K3 moduli does not imply quasi-compactness of
the decorated Fourier--Mukai stack.  A universal theorem would need a
separate uniform kernel bound for the selected side actions.

## 6. Exact verdict on the proposed universal species

Let \(\mathsf {Cox}^{\max}_{\rm arith}\) be the desired class of all marked,
reflection-maximal, finite-volume admissible arithmetic Coxeter data.

### Theorem 6.1 (what universality would have to prove)

An exact functor

\[
 \mathsf {Cox}^{\max}_{\rm arith}
       \longrightarrow\mathsf {K3Mukai}^{\rm mixed}
       \longrightarrow\mathsf {Cox}^{\max}_{\rm arith}
\]

with a nondegenerate full trace block can be essentially surjective only if
all of the following are proved from the raw input:

\[
\begin{array}{ll}
\text{(C1)}&[k:\mathbf Q](n+1)\le23;\\
\text{(C2)}&\text{an integral embedding satisfying every local gluing gate;}\\
\text{(C3)}&\text{simultaneous signed Hodge/crystalline lifts of all sides;}\\
\text{(C3')} &\text{a functor-level Coxeter splitting, if (1.3) is required;}\\
\text{(C4)}&\text{exact no-hidden-wall chamber return;}\\
\text{(C5)}&\text{a uniform boundedness theorem for the decorated kernels.}
\end{array}                                                \tag{6.1}
\]

Conversely, (C1) alone is far from sufficient: (C2)--(C5) are independent
arithmetic, Hodge/crystalline, cone, and moduli conditions.

#### Proof

(C1) is Theorem 2.1; (C2) is (3.1)--(3.2); (C3) is Section 4; Proposition 4.2
gives the independent condition (C3'); exact return requires (C4); and finite
Coxeter image from a quasi-compact decorated stack requires (C5). \(\square\)

The current controls occupy three different points of this theorem:

| chamber | full trace rank | mixed K3 status |
|---|---:|---|
| compact Bugaenko \(H^8\) | 18 | complex K3 carrier exists |
| Vinberg \(X_3\subset H^{19}\) | 20 | mixed complex curve/Fourier--Mukai side lifts exist; no strict \(W\)-splitting is asserted |
| absolutely maximal Borcherds \(H^{21}\) | 22 | the characteristic-two branch supplies all individual side lifts; its spherical twists are not Coxeter involutions as functors |

No frozen reflection-maximal **arithmetic** input of full trace rank at least
24 is currently available in the repository, so Theorem 2.1 is not presented
as an explicit counterexample to the target class.  High-trace finite-volume
Coxeter examples do exist outside arithmeticity (for example dimension four
with trace-field degree eight), and Corollary 2.2 correctly rejects them, but
they are not substitutes for an admissible arithmetic counterexample.

The honest conclusion is therefore sharp:

- the proposed rank-24 K3 species has a rigorous universal capacity ceiling
  \(d(n+1)\le23\);
- twisted, projected, nonprimitive, and genuine K3-category variants do not
  evade it;
- reflection maximality does not prove the embedding, Hodge/crystalline lift,
  categorical splitting, return, or kernel bounds; and
- essential surjectivity over all maximal arithmetic chambers remains an
  unproved theorem, not an available route to finiteness.

## 7. Related frozen controls

- [`bounded_universal_mixed_ag_carrier_no_go_and_next_species.md`](bounded_universal_mixed_ag_carrier_no_go_and_next_species.md)
- [`borcherds_h21_outer_routing_counterexample.md`](borcherds_h21_outer_routing_counterexample.md)
- [`projected_k3_universality_rank_saturation_no_go.md`](projected_k3_universality_rank_saturation_no_go.md)
- [`supersingular_k3_golden_h10_obstruction.md`](supersingular_k3_golden_h10_obstruction.md)
- [`universal_weak_chow_action_carrier_and_leverage_no_go.md`](universal_weak_chow_action_carrier_and_leverage_no_go.md)
