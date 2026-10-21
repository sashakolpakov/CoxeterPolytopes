# Canonical chamber rigidity: the corrected AG-to-ABSW transfer theorem

**Status (2026-08-10): proved synthesis of the project-native carrier
theorems, with one necessary correction to the proposed universal functor.**
This note uses no Nikulin--Vinberg dimension bound, reflective-lattice
classification, Prasad volume estimate, or ABSW finiteness theorem as an
input.

The correction is fundamental:

> A non-rational arithmetic Coxeter chamber cannot be the entire nef chamber
> of an ordinary projective surface.  Literal full-nef carriers have minimal
> Vinberg field \(\mathbf Q\).  A chamber over \(k\ne\mathbf Q\) must instead
> occur as an exact distinguished eigen-nef (or invariant-nef) section of a
> larger rational Neron--Severi space.

With that correction, the desired program separates into two geometrically
different branches.  The rational branch asks for the whole nef cone.  The
non-rational branch asks for the whole restriction-of-scalars block and an
exhaustive projected-wall theorem on its distinguished eigenspace.

## 1. Arithmetic input

Let \(P\subset\mathbb H^n\) be a finite-sided arithmetic Coxeter chamber of
quadratic-form type.  Throughout, \(k\) denotes its **minimal Vinberg
ground field**, assumed totally real, and \(d=[k:\mathbf Q]\).  Write

\[
 T_{\mathbf Q}=\operatorname {Res}_{k/\mathbf Q}V,
 \qquad
 \operatorname {rank}_{\mathbf Q}T_{\mathbf Q}=d(n+1).
 \tag{1.1}
\]

At the distinguished embedding the form has signature \((1,n)\); at every
other embedding it is negative definite.  A **full trace carrier** contains
an integral reflection-stable lattice \(T_{\mathbf Z}\subset T_{\mathbf Q}\),
not just the distinguished real representation or an external local system.

The universal weak rational-surface construction already embeds a common
square multiple of this trace lattice in the Neron--Severi lattice of a
blowup of \(\mathbf P^2\), and algebraic divisor-product correspondences
recover the distinguished eigenspace and all labelled halfspaces.  That
construction is not at issue here.  Its divisors need not be irreducible,
effective, extremal, or nef-complete, so it has no dimension leverage.

## 2. A full nef chamber has ground field \(\mathbf Q\)

### Proposition 2.1 (full-nef rationality)

Let \(S\) be an ordinary projective surface and suppose that its complete nef
cone is finite-sided.  Assume every hyperbolic facet is orthogonal to a
rational numerical curve class \(C_i\in N_1(S)_{\mathbf Q}\), as is automatic
for a literal Mori facet.  If the projectivized nef cone is a Coxeter chamber,
then its minimal Vinberg field is \(\mathbf Q\).

### Proof

The intersection form and every \(C_i\) are rational.  For a negative facet
normal define the Cartan entries

\[
             a_{ij}={2(C_i.C_j)\over C_i^2}\in\mathbf Q.
\tag{2.1}
\]

Replacing \(C_i\) by \(\lambda_iC_i\) changes
\(a_{ij}\) to \((\lambda_j/\lambda_i)a_{ij}\): it rescales row \(i\)
inversely and column \(i\) directly.  Consequently every cyclic product

\[
 a_{i_1i_2}a_{i_2i_3}\cdots a_{i_ri_1}
\tag{2.2}
\]

is invariant under arbitrary nonzero wall rescaling and lies in
\(\mathbf Q\).  A common scaling of the intersection form leaves every
Cartan entry unchanged.  The minimal Vinberg ground field is generated over
\(\mathbf Q\) by these cyclic products, hence is \(\mathbf Q\).
In this numerical realization the reflections themselves also have rational
matrices in \(N^1(S)_{\mathbf Q}\); the field conclusion, importantly, is
the one obtained from the rescaling-invariant cyclic products.  \(\square\)

The conclusion remains true for a normal \(\mathbf Q\)-factorial projective
surface with rational numerical facet classes.  It is not changed by using
irrationally normalized length-two wall vectors: the cyclic products remove
those individual square roots.

Thus the literal requirement

\[
                 P=\mathbf P(\operatorname {Nef}S)
\tag{2.3}
\]

is appropriate only on the \(k=\mathbf Q\) branch.  Requiring (2.3) for all
arithmetic \(P\) would make the proposed functor inconsistent before any
boundedness argument begins.

## 3. Corrected canonical chamber rigidity

The useful carrier property is a disjunction.

### 3.1 Rational full-nef rigidity

For \(k=\mathbf Q\), a strong carrier packet consists of

\[
 (S,\operatorname {NS}(S),K_S,\{C_i\},\pi)
\tag{3.1}
\]

such that:

1. the labelled primitive wall lattice is the declared integral reflection
   host in \(\operatorname {NS}(S)\);
2. every \(C_i\) is an actual irreducible physical facet curve;
3. \(P=\mathbf P(\operatorname {Nef}S)\), including all ideal rays and every
   hidden negative curve;
4. the reflection group is saturated in the declared host, or reflection
   maximality proves that saturation; and
5. a geometric MMP contraction frame, its exceptional lattice, and the
   canonical or log-canonical class are part of the returned data.

### 3.2 Non-rational eigen-nef rigidity

For \(k\ne\mathbf Q\), a strong carrier instead consists of a surface \(S\),
an isometric injection

\[
 T_{\mathbf Z}\hookrightarrow\operatorname {NS}(S),
\tag{3.2}
\]

and a functorial integral self-adjoint \(k\)-action, or an equivalent bounded
integral correspondence algebra.  If \(E_{\sigma_0}\) is its distinguished
Lorentzian real eigenspace, the cone requirement is

\[
 P=\mathbf P\bigl(\operatorname {Nef}(S)_{\mathbf R}
                         \cap E_{\sigma_0}\bigr).
\tag{3.3}
\]

Equality (3.3) means more than recovering the displayed inequalities.  Every
effective curve whose projection cuts the interior of the eigenslice must be
included, and the labelled Gram and orientation must return exactly.  The
trace block must also be orthogonal to a geometric exceptional lattice, or
the numerical orbit of that lattice under the Coxeter operators must be
finite, so that the full trace representation descends through a genuine
minimalization.

This is the corrected form of **canonical chamber rigidity**.  It records
the chamber, the field action, and the MMP geometry without asking an
integral surface cone to have an irrational field of definition.

## 4. The project-native transfer theorem

### Theorem 4.1 (dimension transfer)

Let \(n\ge2\), and let \(\mathcal C\) be a surface-free class of arithmetic
Coxeter chambers which are reflection-maximal (absolutely, relative to the
declared integral host, or among arithmetic reflection groups as specified
by the relevant branch).  Suppose an exhaustive compiler sends every
\(P\in\mathcal C\) to at least one hostile-certified carrier of one of the
following types.

1. **Direct rational carrier.**  A smooth projective rational surface has
   rational-polyhedral complete Mori cone and literal finite-volume full-nef
   chamber on \(L=\operatorname {NS}(S)\).  Every negative extremal curve is
   smooth rational of square \(-1,-2,-4\), every facet reflection preserves
   \(L\), and the facet group is reflection-maximal in one of the senses
   above.  Square-zero tangent Mori rays are allowed.
2. **Certified weak-log rational carrier.**  A smooth projective **complex**
   rational surface carries a possibly empty disjoint reduced boundary
   \(B=\sum B_i\) of smooth rational \((-4)\)-curves, effective irreducible
   negative curves whose oriented halfspaces reproduce the proposed
   finite-volume chamber, and the adapted exponent-two host
   \(M=\langle\operatorname {NS}(S),B_i/2\rangle\).  Every displayed facet
   reflection preserves \(M\), and the facet group is reflection-maximal
   relative to \(M\), or in either stronger sense used by the rational
   branch.  For

   \[
                         D=-(K_S+B/2),
   \]

   the packet supplies a finite affine cover, compatible global sections of
   a Cartier multiple \(mD\), and exact localized unit identities proving
   that those sections generate \(\mathcal O_S(mD)\).  Equivalently, it may
   supply another independently verified base-point-free certificate of the
   same strength.  No complete Mori list is assumed.
3. **Minimalizable full-trace carrier.**  Over \(\mathbf C\), the entire
   integral trace lattice (1.1) embeds isometrically in the canonical
   orthogonal of a minimal smooth projective surface; or it embeds in a
   smooth projective surface with integral pairing isometries extending the
   simple trace reflections and fixing the canonical class, and descends
   through a supplied geometric contraction whose exceptional lattice has
   finite numerical Coxeter orbit.

Then the respective conclusions are

\[
\begin{array}{c|c}
\text{carrier branch}&\text{AG bound}\\ \hline
\text{direct rational}&n\le11,\quad k=\mathbf Q,\\
\text{certified weak-log rational}&n\le25,\quad k=\mathbf Q,\\
\text{minimalizable full trace}&d(n+1)\le20.
\end{array}
\tag{4.1}
\]

Over an algebraically closed field of odd positive characteristic, only the
last-row **numerical** ceiling changes to \(22\), provided the same full
trace injection, canonical-preserving operators, and geometric
finite-orbit descent are certified in the appropriate crystalline/etale
carrier.  This does not promote the complex weak-incidence branch or the
finite-field conclusions below to positive characteristic.  In particular, an
exhaustive complex compiler into these three branches would give the global
dimension bound

\[
                              \boxed{n\le25}.          \tag{4.2}
\]

### Proof

The first row is the direct full-NS theorem: maximality eliminates
\(\mathbb F_0\) and \(\mathbb F_e\) for \(e>2\), forces one free exceptional
chain and the missing line or fibre wall, and the first mixed spacelike ray
gives \(m\le11\) over \(\mathbf P^2\) and \(m\le10\) over \(\mathbb F_2\).
Both have \(n\le11\).  Proposition 2.1 gives the field statement.

For the second row, the localized unit identities make \(mD\)
base-point-free and hence make \(D\) nef.  Adjunction then classifies every
hidden negative curve into the boundary/\((-2)\)/\((-1)\) adapted integral
reflection alphabet.  Reflection maximality gives the full reflection group,
and an ample orientation plus the Lorentz future cone proves that the
proposed chamber is the entire nef chamber, without a complete Mori list.
Contracting the disjoint \((-4)\) boundary gives klt
\(1/4(1,1)\) points and

\[
 D^2=10-\rho+b,
 \qquad
 \widehat c_2=\rho+2-{7b\over4}.
\tag{4.3}
\]

Orbifold \(c_2\)-semipositivity yields \(\rho\le26\), hence \(n\le25\).
Again Proposition 2.1 gives \(k=\mathbf Q\).

For the third row, canonical orthogonality is part of the direct minimal
subcase.  In the nonminimal subcase, canonical preservation and the
finite-orbit equivariant-minimalization theorem descend the full trace block
orthogonally to the minimal canonical class, without bounding the
exceptional complement.  The minimal-surface classification then forces the
carrier into the Kodaira-dimension-zero lane and gives

\[
             \operatorname {rank}T_{\mathbf Q}=d(n+1)\le20.
\tag{4.4}
\]

Taking the maximum of the three ceilings proves (4.2).  \(\square\)

Every premise in Theorem 4.1 is a carrier theorem, not a property inferred
from the existence of a weak NS embedding.  This is precisely where raw
essential surjectivity enters.

## 5. Field degree is not yet a finite field set

Equation (4.4) bounds \([k:\mathbf Q]\), but bounded degree alone does not
give finitely many number fields.  Nor does a bounded or finite-type surface
moduli space: algebraic points in one positive-dimensional family can have
unbounded fields of definition.

This failure persists under much stronger geometric return in hyperbolic
dimension one.  There is an explicit rank-four K3 family with exact
finite-length distinguished eigen-nef intervals, complete restricted-root
saturation, reflection maximality in the declared \(\mathbf Z[t]\)-trace
host, and a primitive ample class of fixed square \(14\), while the minimal
fields \(\mathbf Q(\sqrt{a^2-1})\) are infinite.  The assertion is relative
to the declared order; it does not cover arbitrary overorders or
commensurable hosts, and it does not settle the required \(n\ge2\) case.
Nevertheless it proves that fixed Picard rank, exact finite-volume return,
host-relative maximality, and a fixed polarization square do not by
themselves bound field height.

There are two exact project-native upgrades.

### 5.1 Finite integral Gram image

Suppose instead that the compiler lands in one bounded finite-type
polarized carrier family (equivalently, a finite list of Hilbert
polynomials and projective embedding types), bounds the **total** degrees of
all marked cycles and correspondence components, and uses one global finite
return grammar.  The relevant relative Hilbert/Chow spaces then have
finitely many connected strata after a finite stratification on which the
carrier families are smooth proper and the integral local systems and marked
actions are defined.  Integral cohomology classes, actions, and pairings are
locally constant on those strata.  Hence only finitely many
complete normalized Gram matrices and therefore only finitely many Vinberg
fields occur.  Signed cycles must have both positive and negative parts
bounded.  A polarization degree by itself, a per-component degree bound, or
a per-object return rule does not supply this finite-type premise.

### 5.2 Bounded integral spectral packet

Alternatively, suppose \(k\) is generated by at most \(t\) marked algebraic
integers \(\alpha_j\), each with one distinguished conjugate in \([2,14]\)
and all other conjugates in \([-2,1]\).  The exact project-native auxiliary
polynomial gives

\[
 [\mathbf Q(\alpha_j):\mathbf Q]\le23,
 \qquad
 [k:\mathbf Q]\le23^t.
\tag{5.1}
\]

All conjugates are bounded by \(14\), so the coefficients of every minimal
polynomial are bounded by

\[
       |c_{r-j}|\le {r\choose j}14^j.                 \tag{5.2}
\]

For uniformly bounded \(t\), equations (5.1)--(5.2) leave a finite list of
generators and therefore a finite field set.  The unresolved input is not
Fekete theory; it is proving that the exact carrier returns a bounded number
of integral primitive spectral coordinates with this conjugate envelope.

## 6. From fields to an ABSW-type finite chamber theorem

A dimension ceiling and a finite field set still do not by themselves list
the maximal chambers.  A complete AG replacement for the ABSW conclusion
needs the following finite-host and saturation assertions:

1. only finitely many **isometry classes** of integral reflection hosts
   occur, or else a finite list of complete normalized Gram packets
   (including their full facet lists and exact-return data) is already
   supplied;
2. for every host counted by the first alternative, all displayed facet
   reflections preserve that host and reflection maximality proves that the
   displayed group equals the group generated by **all** integral
   hyperbolic reflections of the host;
3. the locally finite complete mirror arrangement of that full reflection
   group has one chamber orbit.

For the direct rational branch, the rank bound and the elementary
rational-surface blowup calculation leave finitely many unimodular
Neron--Severi hosts; the declared maximality supplies item 2.  For the weak
rational branch, the rank-\(26\) ceiling and
\(L(2)\subseteq M\subseteq L(2)^*=\frac12L\) leave finitely many adapted
hosts, and the complete weak-log transfer proves items 2--3 for the stronger
certified branch.  On the non-rational branch a finite integral Gram image proves finiteness
directly only when it contains the complete return packets; if it merely
produces finitely many trace hosts, the same host-preservation, saturation,
and complete-arrangement hypotheses are still required.  Under exactly
these alternatives the number of maximal chamber classes is finite.
Without the finite integral image, surface boundedness alone is
insufficient.

## 7. What property 4 must actually prove

The exhaustive \(P\)-first theorem is therefore not a single soft map from a
Gram lattice to a surface.  It must prove the following disjunction for every
raw maximal arithmetic input.

1. **Rational canonical escape.**  A reconstructed rational canonical class
   must fail to be nef, and the raw dual cone must expose a geometric
   \(K\)-negative contraction direction.
2. **Incidence/proximity realization.**  The contraction program must lift to
   an actual configuration of points and irreducible curves.  Pairwise Gram
   compatibility does not prove simultaneous effectivity.
3. **Cone completeness.**  The resulting surface must have no unlisted
   negative curve or tangent Mori ray cutting the proposed nef cone.
4. **Non-rational eigen-nef return.**  For \(k\ne\mathbf Q\), the full trace
   block, field action, projected effective walls, and finite-orbit
   minimalization must all be intrinsic.
5. **Finite arithmetic return.**  A bounded integral correspondence packet,
   or a bounded primitive spectral packet, must make the Gram-field image
   finite.

The raw adjunction compiler, exhaustive contraction compiler, universal
cluster charts, and first exact geometric PB-Atlas packet implement parts of
items 1--3.  They do not prove their nonemptiness for every maximal input.
The universal weak eigencarrier implements the lattice and labelled-halfspace
part of item 4, but deliberately not effective or cone-complete return.

Thus the remaining global theorem is now sharply stated: prove that every
maximal arithmetic input enters one of the three branches of Theorem 4.1 and
passes the relevant effectivity, completeness, descent, and finite-return
gates.  Bare maximality cannot be substituted for those gates: the
[rank-three canonical-orientation counterexample](rational_pfirst_effectivity_realizability_no_go.md)
is finite-area and absolutely reflection-maximal, passes the integral
adjunction alphabet, but has its marked rational canonical class inside the
proposed nef chamber and therefore has no rational full-nef realization.

Property 4 is nevertheless no longer only aspirational.  It is proved on
the entire ordinary anticanonical sequence \(I_{1,n}\), \(3\le n\le9\).
For \(3\le n\le8\), host-reflection descent and reduced-characteristic
rigidity force the standard chamber with timelike vertex \(-K\), and the
first \(n\) blowups of the flex base cluster give an exact full-nef carrier.
For \(n=9\), the same sequence terminates at the primitive anticanonical cusp
and the full nine-point flex pencil realizes the unique saturated alcove.
These are variable surface-free exact-carrier classes.  The next rational
universality boundary is twice-primitive \((-4)\) facets, \(-K\) outside the
chamber, and nonunimodular adapted hosts.

## 8. Dependencies

The project-native inputs are:

- the [universal weak eigencarrier](rational_surface_eigencarrier_universality.md);
- the [P-first effectivity and canonical-orientation no-go](rational_pfirst_effectivity_realizability_no_go.md);
- the [timelike anticanonical property-4 theorem](rational_pfirst_timelike_anticanonical_vertex_closure.md)
  and its [hostile audit](rational_pfirst_timelike_anticanonical_vertex_closure_hostile_audit.md);
- the [anticanonical-cusp property-4 theorem](rational_pfirst_anticanonical_cusp_saturation_closure.md)
  and its [hostile audit](rational_pfirst_anticanonical_cusp_saturation_closure_hostile_audit.md);
- the [direct full-NS \(n\le11\) theorem](rational_direct_full_ns_n11_dimension_theorem.md)
  and its [hostile audit](rational_direct_full_ns_n11_dimension_theorem_hostile_audit.md);
- the [weak-incidence orbifold theorem](rational_wrong_log_orbifold_rank_gate.md)
  and its [hostile audit](rational_wrong_log_orbifold_rank_gate_hostile_audit.md);
- the [finite semiampleness certificate](pfirst_weak_log_semiampleness_certificate.md),
  the [complete weak-log transfer](pfirst_weak_log_complete_carrier_transfer.md),
  and its [independent hostile audit](pfirst_weak_log_complete_carrier_transfer_independent_hostile_audit.md);
- the [weak-incidence maximal-finiteness theorem](rational_wrong_log_reflection_maximal_finiteness.md)
  and its [hostile audit](rational_wrong_log_maximal_finiteness_hostile_audit.md);
- the [equivariant minimalization theorem](equivariant_minimalization_trace_descent_theorem.md)
  and its [hostile audit](equivariant_minimalization_trace_descent_hostile_audit.md);
- the [minimal canonical-carrier theorem](project_native_minimal_surface_automorphism_dimension_bound.md)
  and its [hostile audit](project_native_minimal_surface_automorphism_dimension_bound_hostile_audit.md);
- the [integral-packet bridge](integral_packet_global_dimension_bridge.md)
  and its [hostile audit](integral_packet_global_dimension_bridge_hostile_audit.md);
- the [exact interval polynomial theorem](fekete_interval_minus2_1_probe_theorem.md)
  and its [hostile audit](fekete_interval_minus2_1_probe_hostile_audit.md);
- the [exact maximal K3 eigen-interval counterfamily](rank4_k3_exact_maximal_eigen_interval_infinite_fields.md)
  and its [hostile scope audit](rank4_k3_exact_maximal_eigen_interval_infinite_fields_hostile_scope_audit.md); and
- the [polytope-to-surface moduli principle](polytope_to_surface_moduli_principle.md).
