# Unconditional NS/Chow encoding and maximal arithmetic chamber finiteness

**Status (2026-08-11): unconditional theorem, with provenance separated.**

There is an unconditional answer to both literal questions suggested by the
conditional AG-to-ABSW synthesis:

1. every finite-sided finite-volume arithmetic Coxeter chamber presented in
   the standard admissible Vinberg lattice-form model has an exact weak
   encoding in the Néron--Severi lattice and divisor Chow motive of a rational
   surface; and
2. the reflection-maximal chambers among them form only finitely many
   unmarked normalized classes.

These statements are logically independent.  The first is algebraic-
geometric but weak.  The second is classical arithmetic reflection-group
finiteness.  The first does **not** prove the second.

## 1. Unconditional theorem

Let \(P\subset\mathbf H^n\), \(n\ge2\), be a finite-sided finite-volume
arithmetic Coxeter chamber of lattice form.  Let \(k\) be its minimal Vinberg
ground field, let \(d=[k:\mathbf Q]\), and let \(W(P)\) be its side-reflection
group.

### Theorem 1.1 (uniform weak NS/Chow carrier)

There is a rational surface

\[
                   X=\operatorname {Bl}_b(\mathbf P^2)
\]

with

\[
                         b\le2997,                         \tag{1.1}
\]

an integral injection of a positive square rescaling of the full
restriction-of-scalars trace lattice into \(\operatorname {NS}(X)\), a
rational divisor Chow projector \(\pi_U\), a commuting faithful action of
\(k\) and \(W(P)\) on the rank-\(d(n+1)\) Tate summand \((X,\pi_U)\), and
integral divisor classes \(D_i\) indexed by the facets of \(P\), such that

\[
 P\cong
 \left\{[x]\in\mathbf P(E_{\sigma_0}):
 x^2>0,\ (x,D_i)\ge0\text{ for every }i\right\}_{\text{chosen component}}.
                                                               \tag{1.2}
\]

Every side reflection is induced by a rational divisor-product Chow
correspondence.  On the summand it is a self-transpose involution, and all
Coxeter relations hold under literal Chow composition.  Equation (1.2)
preserves the normalized Gram matrix of the projected wall normals
\(\operatorname {pr}_{E_{\sigma_0}}D_i\), oriented half-spaces, incidence,
compactness, and finite-volume status exactly.  The raw ambient trace
pairings \((D_i,D_j)_X\) need not equal the distinguished-place Gram entries.

The explicit bound follows from three proved inputs:

\[
 n<30,\qquad d\le25,\qquad
 b\le4d(n+1)-3,
\]

and therefore

\[
 d(n+1)\le750,\qquad b\le4\cdot750-3=2997.
\]

No maximality is needed for Theorem 1.1.

### Theorem 1.2 (unconditional maximal-chamber finiteness)

If \(W(P)\) is maximal among discrete reflection groups, or only maximal
among arithmetic reflection groups, then the collection of all such \(P\),
over every \(n\ge2\), contains only finitely many unmarked normalized chamber
classes.

More precisely:

\[
 2\le n\le29,\qquad [k:\mathbf Q]\le25,
\]

the set of native pairs \((k,n+1)\) is finite, and only finitely many native
reflective-lattice similarity hosts occur.

For finite-covolume arithmetic inputs the two maximality formulations agree:
any discrete reflection overgroup has finite index and is commensurable with
the original arithmetic lattice, hence is itself arithmetic.

## 2. Proof and two independent classical endings

### 2.1 The weak carrier

For the admissible \(k\)-quadratic module \(M\), choose a totally positive
trace scalar \(\beta\) such that

\[
 B(x,y)=\operatorname {Tr}_{k/\mathbf Q}(\beta q(x,y))
\]

is integral on \(M_{\mathbf Z}\).  Admissibility makes \(B\) Lorentzian of
rank \(R=d(n+1)\).  Rational diagonalization, one difference of two squares,
and four-square decompositions of the negative coefficients give an integral
embedding, after one common square scale,

\[
                    (M_{\mathbf Z},B)\hookrightarrow I_{1,b},
                    \qquad b\le4R-3.
\]

The target is the Néron--Severi lattice of a blow-up of \(\mathbf P^2\).
The orthogonal projector onto the embedded rational subspace is a rational
divisor-product Chow correspondence.  If \(j\) denotes the embedding and
\(g\in W(P)\), then

\[
                         \widehat g=jgj^{-1}\pi_U
\]

acts as \(g\) on the summand and as zero on its complement.  In an
intersection-unimodular divisor basis, every rational endomorphism \(A\) is
represented by

\[
                         Z_A=\sum_a e_a^*\times A(e_a).
\]

Middle-factor contraction gives \(Z_g\circ Z_h=Z_{gh}\) under the standard
correspondence convention.  This proves the group action and its Coxeter
relations.  Pairing an element of the distinguished real eigenspace with
\(j(r_i)\) multiplies the original facet pairing by one positive constant,
which proves (1.2).

The published arithmetic bounds \(n<30\) and \(d\le25\) then give (1.1).
They are imported inputs to the numerical uniformity statement; the Chow
carrier does not reprove them.

### 2.2 Direct global ending

Agol--Belolipetsky--Storm--Whyte prove directly that there are only finitely
many conjugacy classes of arithmetic maximal hyperbolic reflection groups.
Nikulin's 2007 Theorem 2.2 gives an independent published terminal theorem.
A complete mirror arrangement has one reflection-group orbit of chambers,
so group conjugacy finiteness is exactly the required unmarked normalized
chamber finiteness.

### 2.3 Finite-native-pair ending without the ABSW conclusion

There is also a factorized proof which does not invoke the ABSW finiteness
conclusion.

1. Vinberg gives \(n<30\), and the published transition/low-dimensional
   results give \(d\le25\).
2. For the complete reflection chamber of a native reflective lattice,
   Nikulin's bounded-frame lemma selects a connected rank-\(m=n+1\) root
   frame with normalized pair coordinates below \(62^2\).
3. Every cyclic Cartan product \(b_C\) is an algebraic integer.  At the
   distinguished embedding and all definite embeddings,
   \[
                       |\sigma(b_C)|\le62^m\le62^{30}.
   \]
   Full-rank field inheritance and finite-index invariance identify the
   field generated by these products with \(k(P)\).
4. Algebraic integers of degree at most \(25\) with every conjugate bounded
   by \(62^{30}\) form a finite set: their monic minimal-polynomial
   coefficients lie in finite Vieta boxes.  A rank-\(30\) graph has only
   finitely many simple cycles, so only finitely many fields and native pairs
   occur.
5. Nikulin's fixed-\((k,m)\) theorem gives finitely many reflective lattice
   similarity classes.  Reflection maximality identifies \(W(P)\) with the
   full native reflection group, whose complete mirror arrangement has one
   chamber orbit.

This factorization explains exactly where the unconditional finite-pair set
comes from.  It is arithmetic, not an output of the rational surface.

## 3. Why this is not an AG-only proof of finiteness

The weak carrier is universal because it imposes no effective-cone
condition.  Its classes \(D_i\) need not be effective, irreducible,
Mori-extremal, or walls of the nef or movable cone.  Its correspondences need
not be graphs of automorphisms or birational maps.  Blow-ups have unbounded
Picard rank, and divisor-product correspondences realize the full rational
matrix algebra on their Tate part.  Carrier existence therefore gives only

\[
                         d(n+1)\le\rho(X),
\]

with a right-hand side deliberately enlarged to fit the input.  It excludes
no admissible Coxeter datum and supplies no dimension, field, or chamber
bound by itself.

The stronger current categories recover leverage only after losing
universality:

- the same marked Néron--Severi lattice can have different effective and
  nef cones;
- projected physical K3 roots do not bound ambient complements or glue;
- Vinberg's \(X_3\) is one counterexample to universal exact-module
  projected-root complex-K3 return over all finite-volume inputs; and
- the absolutely maximal Borcherds \(H^{21}\) chamber lies outside every
  current complex direct-rational, weak-log, full-trace, and projected-root
  K3 branch, although it has a characteristic-two mixed curve/Cremona
  realization.

Thus the unconditional conclusion is the conjunction of Theorems 1.1 and
1.2, with separate provenance.  It would be circular to use the imported
classical bounds or finiteness theorem and then advertise (1.1) as a new
AG derivation of those same facts.

## 4. Exact answer to the attached conditional statement

The conditional statement can now be replaced by the following scope-safe
answer.

> **Unconditional result.**  Every finite-volume arithmetic Coxeter chamber
> in the standard admissible Vinberg lattice-form model has an exact weak
> NS/Chow carrier on a rational surface, uniformly with at most \(2997\)
> blow-ups.  Every reflection-maximal chamber belongs to a
> finite set of unmarked normalized classes.  Neither assertion assumes an
> algebraic carrier.  The uniform carrier bound and the finite-class theorem
> import classical arithmetic reflection results; weak carrier existence is
> not their proof.

The sharper branchwise algebraic-geometric inequalities from the attached
statement remain valuable but carrier-relative:

\[
 n\le11\quad\text{(direct rational)},\qquad
 n\le25\quad\text{(weak log)},\qquad
 [k:\mathbf Q](n+1)\le20\quad\text{(exact complex K3)}.
\]

No non-tautological effective/nef/mixed-action carrier is currently proved
essentially surjective over all maximal inputs.  An unconditional **AG-only**
proof remains open; the unconditional global theorem above is classical.

## 5. Frozen evidence

- [unconditional classical theorem and finite-pair factorization](unconditional_maximal_arithmetic_reflection_chamber_finiteness.md);
- [independent hostile audit of the classical theorem](unconditional_maximal_arithmetic_reflection_chamber_finiteness_independent_hostile_audit.md);
- [universal weak Chow-action theorem](universal_weak_chow_action_carrier_and_leverage_no_go.md);
- [independent hostile audit of the weak Chow-action theorem](universal_weak_chow_action_carrier_and_leverage_no_go_independent_hostile_audit.md);
- [hostile universality--leverage no-go](unconditional_ag_absw_universality_leverage_no_go.md);
- [finite-native-pair chamber theorem](projected_k3_maximal_abstract_chamber_finiteness.md);
- [rank-\(20\) projected-root no-go](projected_k3_universality_rank_saturation_no_go.md); and
- [absolutely maximal rank-\(22\) outer-routing counterexample](borcherds_h21_outer_routing_counterexample.md).

The companion certificates and independent hostile audits pin the exact
scope and replay the finite-field and Chow-composition steps.
