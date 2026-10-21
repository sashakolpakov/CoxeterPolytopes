# Strict AG-action carriers: the 2-rank theorem and the remaining no-go

**Status (2026-08-12).**  There is a genuine conditional AG-native route, but
it is narrower than “put the Coxeter group in an automorphism group.”  A
faithful strict action on a bounded polarized object bounds the hyperbolic
dimension once the ambient automorphism groups have uniformly bounded
elementary (2)-rank.  A bounded integral full-trace cohomological realization
also gives

\[
                    [k:\mathbf Q](n+1)\le B.              \tag{0.1}
\]

To get **finitely many fields and chambers**, one still needs a
quasi-compact bounded stack of the *action decorations* and a
constructible, cohomologically discrete exact-return map.  Bounded dimension,
Jordan's property, bounded underlying varieties, or bounded graph Hilbert
polynomials alone do not provide that last step.

This report uses no theorem asserting finiteness of arithmetic reflection
groups.

## 1. Three inequivalent meanings of “action”

They have to be separated before applying an AG group theorem.

1. A **strict geometric action** is a faithful homomorphism

   \[
       W(P)\longrightarrow \operatorname {Aut}(X)
       \quad\hbox{or}\quad W(P)\longrightarrow\operatorname {Bir}(X),
                                                               \tag{1.1}
   \]

   taking each simple reflection to an actual involution.  For
   Fourier--Mukai actions it means kernels \(K_w\) with coherent convolution
   isomorphisms

   \[
        K_g*K_h\simeq K_{gh},\qquad K_{s_i}*K_{s_i}\simeq\mathcal O_\Delta.
                                                               \tag{1.2}
   \]

   Thus the Coxeter relations hold at map/functor level, not merely after
   applying cohomology.

2. A **kernel-by-kernel cohomological lift** chooses an actual map or
   Fourier--Mukai kernel for each wall, but only requires its action on a
   numerical or cohomological lattice to be the desired reflection.  The
   products can differ from the Coxeter relations by cohomologically trivial
   autoequivalences.  Derived Torelli supplies lifts of this kind; it does not
   supply a group-theoretic section.

3. A **numerical/Chow action** consists of matrices or algebraic
   correspondences inducing the reflections.  Divisor-product cycles can
   even satisfy all relations literally in a Chow-motive summand, but they
   need not be graphs of maps or invertible Fourier--Mukai kernels.  Jordan
   theorems for \(\operatorname {Aut}(X)\) or \(\operatorname {Bir}(X)\) do not
   apply to this category.

The universal rational-Chow carrier belongs to the third tier.  Its full
matrix algebra is exactly why it is universal and exactly why it supplies no
independent boundedness.

## 2. The local elementary-(2) subgroup forced by a chamber

For a group \(G\), put

\[
 r_2(G)=\sup\{r:(\mathbf Z/2)^r\hookrightarrow G\}.
\]

### Lemma 2.1 (the compulsory local (2)-rank)

If a finite-sided finite-volume Coxeter chamber
\(P\subset\mathbf H^n\) has side group \(W(P)\), then

\[
                         r_2(W(P))\ge \lfloor n/2\rfloor. \tag{2.1}
\]

If (P) has a finite vertex, the lower bound improves to
\(\lceil n/2\rceil\).

#### Proof

At a finite vertex the stabilizer is a finite Coxeter group of rank (n).
Its finite-type Coxeter diagram is a forest.  One of the two parts of its
bipartition has at least \(\lceil n/2\rceil\) vertices; the corresponding
simple reflections commute and are independent, hence generate
\((\mathbf Z/2)^{\lceil n/2\rceil}\).

At an ideal vertex, the horospherical stabilizer is a cocompact Euclidean
Coxeter group in dimension (n-1).  Its finite point group is a finite
Coxeter group of rank (n-1) and occurs as a subgroup of the affine group.
The same forest argument gives rank at least
\(\lceil(n-1)/2\rceil=\lfloor n/2\rfloor\).  A finite-volume polytope has a
finite or ideal vertex, so these are exhaustive. \(\square\)

### Theorem 2.2 (strict-action dimension transfer)

Let \(\mathcal G\) be a class of AG symmetry groups satisfying

\[
                 r_2(G_X)\le a\qquad(G_X\in\mathcal G).  \tag{2.2}
\]

If every chamber in a class has a faithful strict action
\(W(P)\hookrightarrow G_X\), then

\[
                         \boxed{n\le2a+1}.                \tag{2.3}
\]

For compact chambers, \(n\le2a\).

This is an AG-native dimension theorem: it only compares the forced finite
subgroup in Lemma 2.1 with finite subgroups of the geometric symmetry group.

### Corollary 2.3 (polarized regular actions)

Suppose (L) is very ample, (N=h^0(X,L)), and the strict action preserves
the polarization.  Then

\[
 \operatorname {Aut}(X,L)\hookrightarrow \operatorname {PGL}_N.
\]

The adjoint representation

\[
 \operatorname {PGL}_N\hookrightarrow
 \operatorname {GL}(\mathfrak {sl}_N)
\]

is faithful and has dimension (N^2-1).  A faithful characteristic-zero
linear representation of \((\mathbf Z/2)^r\) has dimension at least (r),
because its characters must span the dual of that elementary abelian group.
Consequently

\[
 r_2(\operatorname {Aut}(X,L))\le N^2-1,
 \qquad
 \boxed{n\le2N^2-1};                                    \tag{2.4}
\]

the compact bound is (n\le2N^2-2).

A bounded polarized family has a uniform very ample multiple and a uniform
value of (N), so (2.4) is uniform on such a family.  If the maps are not
required to preserve a common polarization, this argument does not apply
without a separate uniform (2)-rank theorem.

### Why “Jordan” is not the theorem needed

A Jordan bound says that each finite subgroup has an abelian subgroup of
bounded index.  It says nothing about the number of generators of that
abelian subgroup.  The local group can contain, and in right-angled cases can
be, \((\mathbf Z/2)^n\); its Jordan index is already (1).  Hence Jordan's
property alone gives no bound on (n).  One needs (2.2), a uniform bound on
elementary-(2) rank (or the stronger bounded-order property).

For birational actions the same Theorem 2.2 applies whenever a uniform bound
on (r_2(\operatorname {Bir}(X))) is available.  Results on Jordan index and
results on (p)-rank are logically different here.  In particular, the
large-prime rank theorem for rationally connected varieties does not address
the compulsory prime (p=2).

## 3. What bounded cohomological rank proves

Let (P) be admissible arithmetic over its Vinberg field (k), of degree

\[
                         d=[k:\mathbf Q].
\]

For an invariant rank-(n+1) (k)-module (M), its integral full-trace
lattice is

\[
 T_P=\operatorname {Res}_{\mathcal O_k/\mathbf Z}M,
 \qquad \operatorname {rk}_{\mathbf Z}T_P=d(n+1),
 \qquad \operatorname {sig}T_P=(1,d(n+1)-1).             \tag{3.1}
\]

### Theorem 3.1 (full-trace cohomological transfer)

Assume an exact carrier provides a nondegenerate (W(P))-equivariant
injection

\[
                  T_P\otimes\mathbf Q\hookrightarrow H_X\otimes\mathbf Q
                                                               \tag{3.2}
\]

into an integral cohomological or numerical lattice of rank at most (B).
Then

\[
              \boxed{[k:\mathbf Q](n+1)\le B}.           \tag{3.3}
\]

In particular (n\le B-1) and
\([k:\mathbf Q]\le B/(n+1)\).

Strictness of the action is not needed for this rank inequality.  It is
needed if the carrier is advertised as an action of (W(P)), rather than a
list of cohomological lifts of its generators.

For a complex/finite-height K3 Mukai carrier, inertia improves the raw
rank-24 ceiling to (d(n+1)\le21); for a supersingular K3 it gives
\(d(n+1)\le23).  These are the signature bounds proved in the companion
Mukai-capacity report.

### Rank does not imply a finite field list

Equation (3.3) bounds only the degree.  There are infinitely many number
fields of any fixed degree at least two, and integral matrices of fixed size
can have unbounded entries, discriminants, and trace fields.  To pass from
(3.3) to finitely many (k), one must bound the integral action/Gram data or
put it in a genuinely quasi-compact moduli problem on which exact extraction
is discrete.

## 4. The exact positive conditional finiteness theorem

The useful target is a **bounded decorated strict-action stack**, not merely
a bounded family of underlying varieties.

### Theorem 4.1 (bounded strict cohomological carrier)

Let \(\mathcal A\) be a quasi-compact finite-type algebraic stack whose
objects consist of:

1. a polarized AG object (X) in a bounded family;
2. a finite facet marking;
3. graph kernels of regular or birational maps, or perfect
   Fourier--Mukai kernels, in finitely many fixed Hilbert/numerical strata
   and with bounded Tor amplitude;
4. a strict faithful action of the returned Coxeter group as in (1.1) or
   (1.2); and
5. an integral local system (H) of rank at most (B), containing the
   full-trace block (3.2), on which every marked generator is the returned
   reflection.

Assume the exact Gram/chamber extraction is constant on a finite
constructible stratification of \(\mathcal A\) (local constancy is enough).
Then:

\[
 \begin{gathered}
 [k:\mathbf Q](n+1)\le B,\\
 \text{only finitely many marked Coxeter chambers, Vinberg fields, and
 side groups occur in }\mathcal A.
 \end{gathered}                                           \tag{4.1}
\]

If the strict symmetry groups also satisfy (2.2), then (n\le2a+1)
independently of (3.3).

#### Proof

The rank inequality is Theorem 3.1.  A quasi-compact space has only finitely
many members of a finite constructible stratification, and extraction is
constant on each member.  Hence its image in the discrete groupoid of marked
Coxeter data is finite.  The field and group are exact invariants of that
data. \(\square\)

The finite-type hypotheses are natural only after the decorations are
bounded.  Graphs with fixed Hilbert polynomial lie in a Hilbert scheme, and
perfect kernels with fixed numerical data and amplitude lie in bounded
pieces of the moduli of complexes.  The union over all graph degrees, kernel
classes, amplitudes, and numbers of facets is not quasi-compact.

When the facet set is represented by a finite flat marking over
\(\mathcal A\), its degree is locally constant and hence uniformly bounded on
the quasi-compact stack.  If exact return is computed from its marked
sections \(r_i\) of the integral local system, the form on (H), and the
cohomological actions, then the required constructible constancy is automatic
after finite étale labelling and a finite stratification: the pairings
\((r_i,r_j)\) and traces of words are sections of a discrete integral local
system.  It is **not** automatic when return is computed from continuously
varying coefficients of the maps themselves; Section 5 is the explicit
counterexample.

The difficult global statement is therefore essential surjectivity:

> every maximal finite-volume arithmetic chamber must admit such a carrier
> in one fixed \(\mathcal A\), with exact constructible return.

This is not presently proved.  Moreover, (3.3) shows that its proof would
already establish a strong dimension/degree theorem before applying the last
quasi-compactness sentence.

## 5. Fixed-dimensional strict actions still do not finite-ize fields

There is a direct hyperbolic warning.  For every (m\ge7), let

\[
                      P_m=\Delta(2,3,m)\subset\mathbf H^2.
\]

Its Tits/Lorentz reflection representation is faithful and projectivizes to
a strict regular action

\[
                    W(P_m)\hookrightarrow\operatorname {PGL}_3
                    =\operatorname {Aut}(\mathbf P^2).   \tag{5.1}
\]

One may mark the invariant conic and the three mirror lines.  Thus all
carriers have:

- the same smooth polarized surface \((\mathbf P^2,\mathcal O(1))\);
- three involutive linear automorphisms;
- the same graph Hilbert polynomial for every generator; and
- an exact return of the hyperbolic triangle from the conic, mirrors,
  invariant form, selected real hyperbolic component, and oriented mirror
  half-spaces.

Nevertheless their Vinberg fields contain, and with the standard
normalization equal,

\[
                         k_m=\mathbf Q(\cos^2(\pi/m)),    \tag{5.2}
\]

whose degrees are unbounded as (m) ranges through odd primes.  The relation
\((s_is_j)^m=1\) is detected by torsion points of unbounded order inside the
fixed finite-type group \(\operatorname {PGL}_3\).

Therefore even fixed dimension, fixed polarization, a fixed number of
generators, strict regular actions, and fixed graph Hilbert polynomials do
not imply field or chamber finiteness.  Equivalently, the map from these
action parameters to the discrete label (m) is not constant on a finite
constructible stratification; the finite-order locus with arbitrary order is
an infinite union of torsion loci.

These triangles are not asserted to be arithmetic or reflection-maximal.
Their role is a logical one: AG boundedness of strict actions alone cannot be
the arithmetic finiteness mechanism.  Arithmetic integrality plus bounded
cohomological/height data must do additional work.

## 6. The (X_3), (H21), and spherical-twist test

The two endpoint carriers give an exact separation of the three action
notions.

| proposed action species | Vinberg (X_3\subset\mathbf H^{19}) | absolutely maximal Borcherds (H21\subset\mathbf H^{21}) |
|---|---|---|
| individual geometric/FM lifts of every side reflection | yes: 24 curve twists and 12 further FM lifts | yes: 42 curve twists and 168 graph kernels of involutive automorphisms |
| strict action using those displayed lifts | **no** | **no** |
| full-trace capacity | rank 20 fits the complex K3 ceiling | rank 22 fits the supersingular ceiling |
| proof of a class-wide bounded/total species | no | no |

For a (2)-spherical object (S), its spherical twist satisfies

\[
                       T_S(S)\simeq S[-1],\qquad
                       T_S^r(S)\simeq S[-r].              \tag{6.1}
\]

Hence (T_S) has infinite order although its numerical Mukai action is a
reflection of order two.  In particular the displayed (X_3) and (H21)
kernel lists do not satisfy (1.2), so they cannot be used as strict Coxeter
actions.

The same short walls cannot be repaired by regular automorphisms of the same
K3 surface.  If (C\ne0) is effective and (C^2<0), its reflection sends
\([C]\) to \(-[C]\).  A regular automorphism preserves the effective cone,
whereas \(-[C]\) is not effective (it has negative intersection with an
ample class).  Birational selfmaps of a smooth minimal K3 are regular, so
they do not repair this obstruction either.

The conclusion must not be overstated:

- (X_3) and (H21) **pass** the weak kernel-by-kernel mixed-Mukai species;
- their standard spherical-twist choices **fail** the strict Coxeter
  relation;
- this does not prove that every alternative FM lift of the same numerical
  reflection has infinite order in these Picard ranks, nor that no strict
  action on a different AG object exists; and
- therefore the broad strict species is **unproved**, not refuted, while the
  narrower “use the displayed curve twists as Coxeter generators” species is
  refuted.

The Picard-rank-one Nielsen-realization theorem of Fan--Lai gives stronger
non-lifting statements for generic Picard-rank-one K3 surfaces, but it does
not apply directly to (X_3) or the characteristic-two (H21) surface.

## 7. Birational and numerical boundary

Two common shortcuts do not cross the missing gate.

1. **Birational maps.**  The plane has Cremona transformations of unbounded
   degree, so boundedness of (X) does not bound the graph decorations.
   A uniform Jordan index still does not bound the elementary-(2) subgroup
   in Lemma 2.1.  A birational version of Theorem 2.2 needs a genuine uniform
   (2)-rank theorem; a birational version of Theorem 4.1 additionally needs
   bounded graph strata and discrete cohomological extraction.

2. **Chow correspondences.**  On the divisor Tate motive of a rational
   surface, divisor-product cycles realize the full rational matrix algebra.
   They can implement the Coxeter matrices and all relations, but are not
   geometric involutions.  If their cohomological rank is allowed to grow
   with (d(n+1)), the construction is tautologically universal.  If that
   rank is independently bounded, Theorem 3.1 becomes useful—but proving
   that bounded realization for every input is precisely the missing global
   theorem.

## 8. Verdict

The proposed strict-action route has two real pieces of leverage:

\[
 \boxed{r_2(G_X)\le a\ \Longrightarrow\ n\le2a+1},
 \qquad
 \boxed{\operatorname {rk}H_X\le B\ \Longrightarrow\
        [k:\mathbf Q](n+1)\le B}.                         \tag{8.1}
\]

It becomes a finiteness theorem only after a uniform essential-surjectivity
and bounded-decoration theorem puts all carriers in one quasi-compact stack
with constructibly discrete exact return.  Neither maximality nor the known
(X_3/H21) lifts proves that theorem.  Strict actions are therefore a valid
conditional AG-native transfer mechanism, but not a completed global proof.

## 9. Primary and frozen dependencies

- Prokhorov--Shramov,
  [*Jordan property for Cremona groups*](https://arxiv.org/abs/1211.3563):
  Jordan index for finite subgroups of birational groups; the distinction
  from (2)-rank is essential here.
- J. Xu,
  [*A remark on the rank of finite (p)-groups of birational
  automorphisms*](https://arxiv.org/abs/1912.10865): the large-prime rank
  theorem and explicit warning that (p=2) behaves differently.
- Lieblich,
  [*Moduli of complexes on a proper morphism*](https://arxiv.org/abs/math/0502198):
  finite-presentation moduli after the numerical/boundedness data are fixed.
- Fan--Lai,
  [*Nielsen realization problem for derived automorphisms of generic K3
  surfaces*](https://arxiv.org/abs/2302.12663): finite derived actions in the
  generic Picard-rank-one case.
- [`mukai_crystalline_universal_carrier_capacity_and_lift_obstruction.md`](mukai_crystalline_universal_carrier_capacity_and_lift_obstruction.md)
- [`bounded_universal_mixed_ag_carrier_no_go_and_next_species.md`](bounded_universal_mixed_ag_carrier_no_go_and_next_species.md)
- [`universal_weak_chow_action_carrier_and_leverage_no_go.md`](universal_weak_chow_action_carrier_and_leverage_no_go.md)
