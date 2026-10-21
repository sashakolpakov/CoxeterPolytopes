# Hostile scope theorem for an unconditional AG route to arithmetic reflection finiteness

**Status (2026-08-11): unconditional no-go for the weak formulation;
positive global AG theorem not proved.**

## Verdict

No result currently frozen in this repository proves an unconditional
algebraic-geometric version of ABSW for every maximal finite-volume
arithmetic hyperbolic reflection chamber.  This is not merely a missing
calculation.  The present evidence gives a sharp universality--leverage
dichotomy.

1. If a carrier is weak enough to be universal—an occurrence of the full
   trace lattice in a Néron--Severi group, a self-adjoint divisor
   correspondence, and arbitrary integral divisor half-spaces—then every
   admissible lattice-form Coxeter datum has such a carrier on a blow-up of
   \(\mathbf P^2\).  The Picard rank is unbounded and carrier nonemptiness
   has no dimension, field, effectivity, or chamber-finiteness content.
2. If a carrier is strengthened enough to have boundedness leverage—most
   notably an exact complex-K3 root or projected-root return—then it is not
   universal.  Vinberg's \(X_3\) disproves the all-finite-volume
   exact-module projected-root assertion, and the absolutely maximal
   Borcherds \(H^{21}\) chamber lies outside every current complex
   rational/weak-log/full-trace/projected-root K3 branch.
3. Calling arbitrary divisor inequalities or arbitrary Chow
   correspondences “mixed divisor/action data” falls back into the first
   horn.  A non-tautological mixed category has to use effective prime or
   Mori walls and actual automorphisms, pseudo-automorphisms, or birational
   maps, together with an exact fundamental-domain theorem.  One
   characteristic-two K3 realizes the Borcherds endpoint in precisely this
   stronger way, but no class-wide compiler or boundedness theorem for that
   species has been proved.

Thus an unconditional *positive* global theorem can presently be obtained
only by importing the classical arithmetic results themselves: ABSW
directly, or the appropriate classical dimension/ground-field theorems
together with fixed-data reflective-lattice finiteness and reflection
maximality.  That is an unconditional classical route, not an AG-only
derivation.

## 1. What “unconditional” must mean

Let \(\mathcal C_{\max}\) be the surface-free class of finite-sided,
finite-volume arithmetic Coxeter chambers \(P\subset\mathbf H^n\),
\(n\ge2\), whose side groups are maximal in the declared reflection
category.  An honest unconditional AG theorem needs all of the following.

| gate | required assertion | present status |
|---|---|---|
| \(U0\) | The input class is all of \(\mathcal C_{\max}\), with no carrier-existence hypothesis. | This is the desired quantifier. |
| \(U1\) | From raw Coxeter/arithmetic data, construct a nonempty finite exhaustive carrier list for every \(P\). | **Open for every non-tautological proposed stack; false for the current complex stack.** |
| \(U2\) | Every accepted package returns the same normalized chamber, every side, finite volume, and the relevant side-reflection group; all hidden geometric walls and action walls are exhausted. | Proved only branchwise and on scoped inputs. |
| \(U3\) | The geometric conditions are intrinsic: prime/effective/Mori walls or genuine maps/actions, not arbitrary numerical divisor cuts or prescribed correspondences. | Missing from the universal weak carrier; present only in special strong carriers. |
| \(U4\) | An independent AG theorem bounds the carrier stack strongly enough to give finitely many native pairs \((k,n+1)\), or directly finitely many native similarity hosts. | Proved conditionally inside declared rational and complex-K3 branches, not after \(U1\). |
| \(U5\) | Convert the bounded AG output into finitely many maximal chambers. | Proved from finite native pairs only after importing Nikulin's fixed-data theorem and using reflection maximality; an AG-only replacement is absent. |
| \(U6\) | The proof survives the same-NS, golden projected-host, \(X_3\), and Borcherds \(H^{21}\) controls below. | No present positive proposal passes all four. |

The word “unconditional” cannot mean that the implication after a carrier is
supplied is unconditional.  It must include \(U1\), exact return, and the
global quantifier \(U0\).

## 2. Formal universality--leverage dichotomy

### Theorem 2.1 (weak carriers are universal and therefore give no
nonemptiness obstruction)

Let \(P\) be an admissible totally-real arithmetic Coxeter datum of
lattice-form type.  Let its full restriction-of-scalars trace lattice have
rank

\[
                         R=[k:\mathbf Q](n+1).
\]

There is a blow-up \(X\) of \(\mathbf P^2\), an integral injection of a
positive common rescaling of that trace lattice into
\(\operatorname {NS}(X)\), a rational Chow projector carrying the field
action, a self-adjoint divisor-product correspondence, and integral divisor
classes \(D_i\), such that the distinguished eigenspace and the inequalities

\[
                         (x,D_i)\ge0
\]

return \(P\) exactly as a semialgebraic chamber.  One may take at most
\(4R-3\) blow-ups.

This is the proved
[rational-surface eigen-carrier theorem](rational_surface_eigencarrier_universality.md).
The \(D_i\) need not be effective, irreducible, extremal, or walls of the
nef or movable cone.  Every rational self-adjoint operator on the divisor
space is represented by a rational divisor-product Chow correspondence.

Consequently, the assertion

\[
          \text{“there exists an NS/divisor/correspondence carrier”}
\]

eliminates no admissible input.  Since blow-ups have unbounded Picard rank,
it also supplies no uniform cohomological ceiling.  In particular,
nonemptiness of this weak carrier category alone cannot yield a dimension
bound, a finite set of ground fields, finite native hosts, or ABSW
finiteness.

This conclusion is deliberately limited to the weak carrier's existence
condition.  It does not say that no separately imposed, genuinely geometric
condition on a rational surface can ever have leverage.

### Corollary 2.2 (arbitrary mixed data are tautological)

Suppose “mixed divisor/action carrier” permits

- arbitrary integral divisor classes as facet inequalities; and
- arbitrary rational Chow correspondences as the wall or field actions.

Then it contains the universal category of Theorem 2.1.  Exact return merely
repeats the input half-spaces inside a divisor lattice.  This enlarged
terminology does not repair \(U3\) or produce a new bound.

To avoid this corollary, a divisor side must have independently checkable
geometric status—such as a prime effective negative curve, an extremal
Mori/MBM divisor, or a component of an intrinsic discriminant—and an action
side must come from an actual automorphism, pseudo-automorphism, birational
map, monodromy operator with a geometric realization theorem, or comparably
rigid structure.  Its fundamental-domain and no-hidden-wall assertions must
then be proved, not inserted into the carrier definition.

### Proposition 2.3 (a bounded-rank carrier moves all difficulty into
essential surjectivity)

Let a carrier category have a uniform Picard or divisor-cohomology ceiling
\(B\), and require an injective full trace block.  Every carrier in it
satisfies

\[
                         [k:\mathbf Q](n+1)\le B.       \tag{2.1}
\]

Therefore essential surjectivity of this category over
\(\mathcal C_{\max}\) would immediately prove the corresponding dimension
and degree bounds.  Rank counting proves (2.1); it does **not** prove
essential surjectivity.

For complex K3 surfaces \(B=20\).  The rank consequence is useful, but the
exact projected-root category is not universal.  For K3 surfaces in
arbitrary characteristic \(B=22\); the attractive consequences

\[
                         n\le21,\qquad
                  [k:\mathbf Q](n+1)\le22
\]

remain conditional on a mixed-carrier existence theorem which has not been
proved.  Presenting those inequalities as unconditional without proving
\(U1\) would simply assume the main theorem at the essential-surjectivity
gate.

## 3. Sharp hostile controls

### 3.1 Bare NS does not determine the chamber

The blow-ups of \(\mathbf P^2\) at three noncollinear points and at three
collinear points have the same marked lattice \(I_{1,3}\), the same
canonical class, and the same displayed Chern and Hodge invariants.  Only
the collinear surface has the effective class
\(H-E_1-E_2-E_3\) giving an extra wall.  Hence even
\((\operatorname {NS},K)\) does not determine the effective, ample, or nef
cone.  See the
[NS-versus-nef audit](universal_ns_realization_and_nef_leverage_audit.md).

**Refuted formulation:** a global geometric property of the Coxeter
polytope is encoded by the bare NS lattice in a way that automatically
recovers its complete chamber.

### 3.2 Exact projected walls do not finite-ize the ambient host

For the compact golden \((2,5,5)\) triangle, an all-parameter family of
rank-nine K3 packages fixes the field, order, normalized chamber, trace
eigenspace, and physical projected square-\((-2)\) facets, while

\[
 [S_s:T_s\oplus K_s]=s^3,\qquad
 |\operatorname {disc}S_s|=64s^3(s-1)^3
\]

diverge.  See
[the projected-host scaling no-go](k3_compact_projected_root_scaling_no_go.md).

**Refuted formulation:** physical projected K3 roots and exact eigen-nef
return imply finitely many ambient NS lattices, complements, or glue
packages.

This does not obstruct abstract maximal-chamber counting after finite native
pairs and reflection maximality.  It only forbids using projected ambient
host finiteness as the missing unconditional theorem.

### 3.3 \(X_3\) refutes exact-module projected-root universality

Vinberg's finite-volume \(H^{19}\) chamber on

\[
                  U\oplus2E_8(-1)\oplus A_2(-1)
\]

has 24 primitive norm-\((-2)\) facet lines and 12 primitive
norm-\((-6)\), divisibility-three lines.  At complex K3 Picard rank 20 the
projector complement vanishes.  A common rational similarity and ambient
square-\((-2)\) lifts would force

\[
                          (a/b)^2=3,
\]

which is impossible over \(\mathbf Q\).  See the
[\(X_3\) rank-saturation theorem](projected_k3_universality_rank_saturation_no_go.md)
and its
[independent audit](projected_k3_universality_rank_saturation_no_go_independent_hostile_audit.md).

The exact quantifier is important.  This is one counterexample to
universality over all finite-volume inputs.  Absolute reflection maximality
of \(X_3\) is not proved here, so it does not settle a maximal-only
conjecture.  Its 12 short walls remain valid integral divisor cuts, but they
are not asserted in the frozen carrier theorem to be irreducible curves.
Thus arbitrary mixed-divisor return survives only by weakening the geometric
gate that was supposed to supply leverage.

### 3.4 The absolutely maximal \(H^{21}\) chamber defeats the current
complex stack

Borcherds' 210-wall chamber on

\[
                         U\oplus D_{20}(-1)
\]

is finite-volume, noncompact, arithmetic, and absolutely
reflection-maximal.  Its native rank is 22.  Therefore every complex K3
full-trace or projected-root return is impossible.  The direct rational
branch and the exhaustive weak-log rational compiler are also empty; the
last surviving weak-log colourings would realize
\(\mathbf P^2(\mathbf F_4)\) as a complex point-line configuration, and its
Fano subplane has determinant \(-2\ne0\).

See the
[Borcherds outer-routing theorem](borcherds_h21_outer_routing_counterexample.md)
and its
[independent audit](borcherds_h21_outer_routing_counterexample_independent_hostile_audit.md).

This proves failure of the current complex disjunction on a maximal input.
It does **not** prove that every conceivable complex-surface mixed-action
carrier is impossible.

The positive escape is nevertheless highly informative.  Over
characteristic two, the supersingular Artin-invariant-one K3 has 42 smooth
rational curve walls and 168 genuine regular automorphisms induced by
quintic Cremona involutions; on its Picard lattice those automorphisms act
as the 168 norm-\((-4)\) reflections.  This is an actual mixed
curve/action realization, not an arbitrary correspondence.  It shows the
kind of fourth species a viable theorem must include, while proving nothing
about its class-wide existence.

### 3.5 A prime wall cannot itself be its regular reflection action

Let \(X\) be a projective variety, let \(0\ne[C]\in\operatorname{NS}(X)\)
be the class of an effective divisor, and suppose \([C]^2<0\) in a surface
intersection lattice.  No regular automorphism of \(X\) can act on
\(\operatorname{NS}(X)\) as the orthogonal reflection

\[
            s_C(x)=x-{2(x,C)\over(C,C)}C.
\]

Indeed, \(s_C(C)=-C\), whereas pullback and pushforward by a regular
automorphism preserve the effective cone.  The class \(-C\) is not
effective: an ample class has positive intersection with every nonzero
effective divisor and negative intersection with \(-C\).

**Refuted formulation:** every effective prime root wall is simultaneously
implemented by a regular automorphism acting as reflection in that same
root.  A viable strong carrier must distinguish divisor walls from action
walls, or allow a genuinely birational/monodromy mechanism with its own
geometric realization theorem.  The characteristic-two \(H^{21}\) model
does exactly the former: its 42 curve walls and 168 automorphism-reflection
vectors are different species.

## 4. Requirement-by-requirement audit of the attached conditional result

| attached assertion | mathematical status | unconditional-global verdict |
|---|---|---|
| Direct rational exact carriers give \(n\le11\). | Proved inside that carrier species. | Does not prove every raw maximal chamber has the carrier. |
| Weak-log rational carriers give \(n\le25\). | Proved inside the certified weak-log species. | The compiler may be empty; \(H^{21}\) is an absolutely maximal empty example. |
| Exact complex-K3 trace carriers give \([k:\mathbf Q](n+1)\le20\). | Immediate and proved after exact full-trace return. | Exact projected-root universality is false; \(H^{21}\) is outside complex K3 by rank. |
| Exact admissible finite-volume full-trace/order-module complex-K3 carriers in \(n\ge2\) have finitely many fields. | Proved carrier-relative, including the cusp collapse and compact discriminant bound. | Does not construct the carrier and therefore does not finite-ize the fields of all raw inputs. |
| Finite native pairs plus reflection maximality give finitely many normalized chambers. | Proved using Nikulin's fixed-data reflective-lattice theorem. | Unconditional only after a global finite-pair theorem is supplied; not AG-only because the downstream step imports Nikulin. |
| The current compiler is exhaustive. | Contradicted by \(\Delta(2,4,6)\), \(X_3\), and the absolutely maximal \(H^{21}\) endpoint at successively broader gates. | **False.** |
| Bare NS or arbitrary divisor inequalities encode the needed global geometry. | Same-NS examples refute the first; the universal rational weak carrier makes the second tautological. | **False or non-leveraging, depending on interpretation.** |

The gap is exactly the first universal arrow, not the branchwise
inequalities after that arrow.

## 5. Minimal non-tautological category not yet ruled out

The smallest plausible replacement suggested by the hostile controls is a
**bounded mixed cone/action fundamental-domain carrier**.  A package would
have to contain:

1. a projective K3 surface in arbitrary characteristic, or another
   deformation class with a uniform divisor-cohomology ceiling;
2. an injective full restriction-of-scalars trace/order module, with
   admissibility and the integral native reflection action retained;
3. the exact nef, movable, or distinguished eigencone;
4. for each raw facet, either an intrinsic effective prime/Mori/MBM wall or
   an actual automorphism, pseudo-automorphism, or birational wall-pairing
   inducing the required native reflection;
5. exact equality of \(P\) with a geometric fundamental domain, including
   every zero-projection, hidden divisor wall, hidden action wall, and
   chamber stabilizer;
6. a finite exhaustive compiler from raw \(P\), not a carrier-success
   condition in the definition; and
7. an independent boundedness theorem controlling characteristic,
   polarization/action complexity, native fields, and native similarity
   hosts—or an explicit declaration that the last host step imports
   Nikulin.

This category is “minimal” only in the hostile logical sense:

- positive characteristic is forced by the known \(H^{21}\) realization
  if one wants the demonstrated endpoint species;
- action walls are forced by its 168 long reflections;
- non-literal or internal walls are forced by \(X_3\);
- full trace/order data are forced by the field problem;
- exact fundamental-domain equality is forced by same-NS and hidden-wall
  examples; and
- bounded cohomology is what creates any AG dimension leverage.

No theorem currently proves that every \(P\in\mathcal C_{\max}\) has such a
package.  It is therefore a research target, not an unconditional result.
If “mixed” is weakened by dropping items 4 or 5, Corollary 2.2 applies and
the category becomes a universal numerical encoding with no ABSW leverage.

## 6. What is unconditionally closed

There are two different kinds of unconditional theorem.

1. **Unconditional negative AG-route statements.**  The weak
   NS/divisor/correspondence category is universal and non-leveraging; the
   exact complex projected-root K3 category is not universal; and the
   current complex outer compiler fails on an absolutely maximal chamber.
   These are genuine proved theorems.
2. **Unconditional positive global finiteness.**  The classical ABSW theorem
   already proves that there are finitely many conjugacy classes of maximal
   arithmetic hyperbolic reflection groups.  Classical Vinberg,
   Nikulin/ground-field, Esselmann in the integral rational endpoint, and
   related arithmetic results supply the dimension and field controls in
   their exact scopes.  The repository's finite-native-pair theorem can
   compose suitable classical finite-pair inputs with Nikulin fixed-data
   finiteness and reflection maximality.  This composition is frozen in the
   [unconditional classical chamber theorem](unconditional_maximal_arithmetic_reflection_chamber_finiteness.md):
   (n\le29), \([k:\mathbf Q]\le25\), and, in the noncompact case,
   (k=\mathbf Q\) with (n\le19) or (n=21), followed by finite native
   pairs and finite unmarked normalized maximal chambers.

The second item is the only presently closed positive global route.  It may
be restated in the AG vocabulary, but doing so does not turn its arithmetic
inputs into an AG-only proof.  Using ABSW itself to certify the carrier list
would be circular if the stated goal is a new derivation of ABSW.

## 7. Scope-safe final statement

The strongest honest conclusion is:

> The existing algebraic geometry proves strong dimension and ground-field
> implications after exact carriers are supplied, but no non-tautological
> carrier species is presently proved essentially surjective over all
> maximal finite-volume arithmetic reflection chambers.  Weakening the
> carrier to bare NS, arbitrary divisor inequalities, or arbitrary Chow
> correspondences makes it universal and removes the desired leverage.
> Strengthening it to exact complex K3 root geometry restores leverage and
> loses universality.  A bounded positive-characteristic mixed
> curve/action fundamental-domain compiler is the minimal surviving
> proposal, not a theorem.  Hence no unconditional AG-only ABSW theorem may
> currently be claimed; unconditional global finiteness remains classical.

## Evidence

- [fixed-data/outer-routing synthesis](ag_to_absw_fixed_data_closure_and_outer_routing_obstruction.md)
- [independent synthesis audit](ag_to_absw_fixed_data_outer_routing_synthesis_independent_hostile_audit.md)
- [four-property completion audit](canonical_chamber_four_property_completion_audit.md)
- [universal rational weak carrier](rational_surface_eigencarrier_universality.md)
- [bare-NS versus nef-cone audit](universal_ns_realization_and_nef_leverage_audit.md)
- [golden projected-host counterfamily](k3_compact_projected_root_scaling_no_go.md)
- [\(X_3\) projected-root no-go](projected_k3_universality_rank_saturation_no_go.md)
- [Borcherds \(H^{21}\) outer-routing counterexample](borcherds_h21_outer_routing_counterexample.md)
- [finite-native-pair maximal-chamber theorem](projected_k3_maximal_abstract_chamber_finiteness.md)
- [unconditional classical maximal-chamber theorem](unconditional_maximal_arithmetic_reflection_chamber_finiteness.md)

Primary external calibrations:

- I. Agol, M. Belolipetsky, P. Storm, and K. Whyte,
  [Finiteness of arithmetic hyperbolic reflection groups](https://doi.org/10.4171/GGD/47).
- I. Dolgachev and S. Kondō,
  [A supersingular K3 surface in characteristic 2 and the Leech lattice](https://arxiv.org/abs/math/0112283).
- F. Esselmann,
  [Über die maximale Dimension von Lorentz-Gittern mit coendlicher Spiegelungsgruppe](https://doi.org/10.1006/jnth.1996.0141).
