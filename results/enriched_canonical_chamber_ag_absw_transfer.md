# Enriched canonical-chamber transfer: the exact AG route to dimension and field finiteness

**Status (2026-08-11): proved conditional synthesis.**  This theorem isolates
the precise extra datum used by the present proof to turn the four-property
POLYTOPE-to-SURFACE programme into dimension and finite-field theorems.  The
first four gates construct and return the right chamber.  A fifth,
class-uniform **polarized spectral-height** gate is a sufficient arithmetic
input.  It is not known whether the first four gates force it in dimensions
\(n\ge2\).

This is not an assertion that every maximal arithmetic Coxeter chamber
passes the gates.  Essential surjectivity is proved on the standard marked
ordinary rational lane and remains open on the exponent-two
twice-\((-4)\) and non-rational eigen-nef lanes.

## 1. The five gates

Let \(P\subset\mathbb H^n\), \(n\ge2\), be a finite-volume maximal
arithmetic Coxeter chamber with minimal Vinberg field \(k\),
\(d=[k:\mathbf Q]\), and field module \(V\) of rank \(n+1\).  An enriched
exact carrier consists of the following data.

1. **Trace fidelity.**  The complete restriction-of-scalars module
   \[
              T=\operatorname {Res}_{k/\mathbf Q}V,
              \qquad \operatorname {rank}_{\mathbf Q}T=d(n+1),
   \tag{1.1}
   \]
   its acting order, and its cyclic Cartan field are returned integrally.
   The order-module is finite projective and every facet reflection preserves
   it; otherwise the compiler must separately certify cyclic-product
   integrality and every multiplication operator used below.
   On the complex K3 branch the intersection form is the one-positive trace
   form \(q=\operatorname {Tr}_{k/\mathbf Q}(\alpha h)\), with the sign of
   \(\alpha\) chosen so that the distinguished place is Lorentzian and every
   other place contributes a negative-definite summand.
2. **Chamber fidelity.**  If \(k=\mathbf Q\), the carrier returns
   \(P=\mathbf P(\operatorname {Nef}S)\).  If \(k\ne\mathbf Q\), it
   returns \(P\) as the exact distinguished eigen-nef section of the full
   trace lattice.  A literal full nef cone cannot carry a non-rational
   minimal Vinberg field.
3. **Canonical/MMP fidelity.**  The carrier returns either the canonical
   class, a nef/semiample log class with its morphism certificate, or a
   finite-orbit equivariant minimalization through which \(T\) descends.
4. **Essential surjectivity.**  Every raw chamber in the declared
   surface-free class receives an actual carrier satisfying 1--3; proposed
   wall classes are simultaneously irreducible effective, and the return is
   exhaustive rather than a list of familiar examples.
5. **Uniform arithmetic height.**  Every non-rational input used in the
   field conclusion is returned on an exact full-Picard complex K3 trace
   carrier.  The compiler returns a big-and-nef class \(H\in T\) with
   \(H^2\ge s_0>0\).  For a field-generating selection of the large affine
   signed-cycle coordinates \(\theta_C\), it returns algebraic
   correspondences acting by multiplication and a class-uniform bound
   \[
                   H\cdot\tau_{\theta_C}(H)\le C.       \tag{1.2}
   \]
   For finiteness of chambers rather than only fields, it additionally
   returns finitely many integral host isometry classes or finitely many
   complete normalized Gram packets.

Gates 1--4 are the corrected form of canonical chamber rigidity.  Gate 5
is a class-level boundedness property, not part of the existence of one
carrier.

## 2. Dimension transfer

### Theorem 2.1 (AG dimension ceiling)

Suppose an exhaustive compiler puts every input into one of the following
hostile-audited branches:

1. a direct rational full-Neron--Severi carrier whose complete negative
   hyperbolic facets are smooth rational curves of square \(-1,-2\), or
   \(-4\), whose facet reflections are integral, and whose reflection group
   is host-maximal;
2. a smooth projective complex rational weak-log carrier with disjoint
   \((-4)\) boundary, effective irreducible proposed facets, an exact
   basepoint-free certificate for \(D=-(K+B/2)\), and adapted-host
   reflection maximality; or
3. a minimal, or finite-orbit-minimalizable, complex full-trace carrier with
   the canonical-preserving trace action required by the descent theorem.

Then respectively

\[
 \begin{array}{rcl}
 \text{direct rational}&:&n\le11,\quad k=\mathbf Q,\\
 \text{weak-log rational}&:&n\le25,\quad k=\mathbf Q,\\
 \text{complex full trace}&:&d(n+1)\le20.
 \end{array}                                             \tag{2.1}
\]

Consequently every input has \(n\le25\).  In the non-rational branch,
\(d\ge2\), so in fact \(n\le9\).  On the standard marked ordinary
rational lane the stronger P-first theorem gives \(n\le9\), and every
maximal survivor in dimensions \(3\) through \(9\) has an explicit exact
flex-pencil carrier.

### Proof

The first row is the direct rational MMP/free-chain theorem.  The second is
the semiampleness, exact-return, boundary-contraction and orbifold
second-Chern theorem.  In the third row, canonical orthogonality or
finite-orbit equivariant descent puts the full Lorentzian trace block on a
minimal Kodaira-dimension-zero surface.  Over \(\mathbf C\), the K3 Picard
ceiling gives (2.1).  Taking the largest of the three row bounds gives
\(n\le25\).  The non-rational and ordinary refinements are immediate from
\(d\ge2\) and the ordinary wall-cone theorem. \(\square\)

No Nikulin--Vinberg dimension classification, ABSW theorem, arithmetic
volume formula, or reflective-lattice census enters this implication.

## 3. Why the fifth gate gives finitely many fields

The cyclic Cartan field need not be generated by pair products.  On a finite
projective returned order-module, traces of reflection words show that every
directed cyclic Cartan product is integral.  Closed walks reduce to pair
cycles and chordless simple cycles.  At any non-distinguished embedding the
normalized arithmetic form \(h^\sigma\), equivalently \(-q_\sigma\), is
positive definite; the Gram matrix of a chordless \(\ell\)-cycle is cyclic
tridiagonal, and its kernel injects into two coordinates.  Therefore

\[
                             \ell\le n+3.              \tag{3.1}
\]

For such a cycle put

\[
             \theta_C=(-1)^\ell c(C)+2^\ell.          \tag{3.2}
\]

Positive-definite Cauchy--Schwarz bounds every non-distinguished conjugate
of the signed cycle by \(2^\ell\).  Hence \(\theta_C\) is totally positive
and is either at most \(2^{\ell+1}\) at every place, or its distinguished
conjugate dominates all the others.  In the dominant case, the
one-positive-place Hodge identity gives

\[
 H\cdot\tau_{\theta_C}(H)
       \ge \sigma_0(\theta_C)H^2.                     \tag{3.3}
\]

Thus Gate 5 bounds every conjugate of every selected generator.

### Theorem 3.1 (uniform finite field set)

Assume Gates 1--5 with uniform \(s_0,C\), including the exact full-Picard K3
refinement in Gate 5 for every non-rational input.  Then the Vinberg fields
occurring in the declared class form an explicitly finite set.

For a fixed \(n\), first put

\[
                 D_n=\left\lfloor{20\over n+1}\right\rfloor.
                                                        \tag{3.4}
\]

If \(D_n<2\), then \(k=\mathbf Q\).  If \(D_n\ge2\), put

\[
 \begin{aligned}
 L_n&=n+3,\\
 r_n&=\lfloor\log_2D_n\rfloor,&
 M_n&={D_n(D_n-1)\over2}+1,\\
 R_{0,n}&=\max\{2^{L_n+1},C/s_0\},&
 R_n&=r_nM_nR_{0,n}.
 \end{aligned}                                         \tag{3.5}
\]

Then

\[
                   |\operatorname {Disc}k|
                   \le R_n^{D_n(D_n-1)}.              \tag{3.6}
\]

Since every non-rational input has \(2\le n\le9\), a dimension-uniform
coarse bound is

\[
 |\operatorname {Disc}k|
 \le \left(32\max\{8192,C/s_0\}\right)^{30}.          \tag{3.7}
\]

For \(n=2\), the sharper bound is

\[
 |\operatorname {Disc}k|
 \le \left(32\max\{64,C/s_0\}\right)^{30}.           \tag{3.8}
\]

### Proof

If \(k=\mathbf Q\), there is nothing to prove.  If \(k\ne\mathbf Q\),
Gate 5 places the carrier on the complex full-trace K3 branch, where (2.1)
gives \(d\le D_n\); in particular \(D_n<2\) excludes this case.  For
\(D_n\ge2\), equations (3.1)--(3.3) bound all conjugates of every selected
affine cycle coordinate by \(R_{0,n}\).  A field of degree at most \(D_n\)
needs at most \(r_n\) members of any generating set: each strict
intermediate extension has degree at least two.  A positive grid with
coefficients at most \(M_n\) contains a primitive linear combination,
whose conjugates are bounded by \(R_n\).  The squared Vandermonde gives
(3.6), and the elementary symmetric functions give finite integer boxes
for every monic minimal polynomial.  Equations (3.7)--(3.8) are the maxima
of these constants over the indicated finite dimension ranges. \(\square\)

The conclusion is a literal finite field set, not only bounded field
degree.

### Corollary 3.2 (bounded Hilbert/Chow grammar supplies Gate 5)

Suppose instead that the exhaustive compiler lands in a finite list of
polarized projective embedding types and uses one finite return grammar.  It
must return every class of a **complete signed marked Neron--Severi basis**
and every signed correspondence as a difference of positive and negative
effective parts, bound the total Hilbert/Chow degree of **both** parts, and
bound every rational denominator.
Assume the marked basis spans the exact Picard lattice on every accepted
stratum.  Then Gate 5 holds, including its finite-host clause.

Indeed, the mixed number in (1.2) is one coefficient of the bounded
correspondence multidegree.  After a finite smooth-proper stratification of
the relevant Hilbert/Chow spaces, the marked integral cohomology classes,
their actions, and their pairings are locally constant.  There are finitely
many strata, so only finitely many complete marked intersection/action
matrices and host lattices occur.  A bounded polarization square alone,
bounded degree per component, or a positive-dimensional finite-type family
without bounded complete markings does not imply this conclusion.

## 4. From fields to maximal-chamber finiteness

Assume in addition the finite-host clause of Gate 5, preservation of every
host by the displayed reflections, and reflection maximality relative to
that host.  Adjoining any host-integral reflection remains discrete, so
maximality gives the full reflection subgroup \(W(P)=R(M)\).  The locally
finite complete mirror arrangement has a connected chamber graph and one
\(R(M)\)-orbit of chambers.  Hence finitely many hosts, or finitely many
complete normalized packets, give finitely many maximal chambers.

This last step is logically separate from the field theorem: bounded field
discriminant does not by itself bound the discriminant or isometry class of
the integral reflection host.

## 5. Sharpness of every additional gate

The hypotheses cannot be replaced by a weak lattice embedding.

- Without Gate 3, exceptional blowup tails have arbitrary rank and the
  trace block need not descend.
- Without Gate 4, an absolutely reflection-maximal rank-three packet can
  put the proposed rational canonical class inside the chamber and admit no
  rational full-nef carrier.
- The marked lattice and a nef isotropic class do not determine its
  semiampleness: the normal character is an independent period datum.
- The first four weak axioms do not imply Gate 5 once hyperbolic dimension
  one is allowed: rank-four minimal K3 carriers have exact finite-volume
  eigen-nef return, complete restricted-mirror saturation, relative
  reflection maximality, and a primitive polarization of fixed square
  \(14\), yet realize infinitely many fields
  \(\mathbf Q(\sqrt{a^2-1})\).  Their escaping scalar is
  \[
                         H\cdot\tau_t(H)=32a-10.       \tag{5.1}
  \]
  This counterfamily does not prove that Gate 5 is logically necessary in
  the present \(n\ge2\) scope; automatic higher-dimensional height remains
  open.

Thus the unresolved implication is now exactly

\[
 \begin{split}
 &\text{maximal exact non-rational eigen-nef return}\\
 &\qquad\Longrightarrow
 \text{canonical }H\in T\text{ and bounded large-cycle bidegrees.}
 \end{split}                                           \tag{5.2}
\]

On the rational side, the unresolved essential-surjectivity boundary is the
twice-primitive \((-4)\) lane, where \((-K).B=-2\) and the ordinary
anticanonical wall-cone argument stops.

## 6. Frozen inputs

The proof composes:

- the [canonical-chamber transfer theorem](canonical_chamber_rigidity_absw_transfer_theorem.md);
- the [direct rational theorem](rational_direct_full_ns_n11_dimension_theorem.md);
- the [complete weak-log transfer](pfirst_weak_log_complete_carrier_transfer.md);
- the [equivariant minimalization theorem](equivariant_minimalization_trace_descent_theorem.md);
- the [combined ordinary property-4 theorem](rational_pfirst_ordinary_dimension_bound_and_exact_carrier.md);
- the [chordless-cycle height theorem](k3_chordless_cycle_height_field_bound.md); and
- the [rank-four exact-maximal height counterfamily](rank4_k3_exact_maximal_eigen_interval_infinite_fields.md).

The accompanying certifier pins these dependencies and replays every
dimension and discriminant constant.  It is a dependency-and-implication
certificate, not an arbitrary raw-polytope carrier compiler.
