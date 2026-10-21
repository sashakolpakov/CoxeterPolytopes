# Dessin-native field effectivity without arithmetic-volume finiteness

**Status (2026-08-09): two unconditional finite-image/Galois implications,
plus conditional geometric gates and an explicit surface criterion.**  The
mechanism in this note is not Prasad's arithmetic covolume machinery and does
not use the ABSW finiteness theorem. The separate Mostow--Prasad strong
rigidity theorem is used only to prove uniqueness of a complete realization
once its Coxeter matrix is fixed. Its open hypothesis is instead geometric:
the outward unit-normal geometric Gram image (equivalently, the Cartan image
modulo diagonal root rescaling) of the exhaustive Coxeter--Borcherds carrier locus must
be proved zero-dimensional.  Zero-dimensionality of the full carrier stack
is needed only to enumerate carrier fields and obtain descent.

## 1. The moduli object that is actually needed

Fix \(n\) and a finite bound package \(\mathcal D\).  It contains:

- an ambient projective dimension and Hilbert polynomials for a surface, its
  ordered marked curves, the zero-cycles used for ordinary faces, and the
  exceptional clusters/dual complexes used for ideal vertices;
- a bounded **universal coefficient-order grammar**: bounds on its degree,
  module presentation, denominators, and multiplication/eigenvalue
  relations, together with a bound on the number of ordered primitive roots,
  their norm/divisibility alphabet, and scalar outward unit-normal geometric
  Gram variables on a chosen finite splitting/eigencharacter cover;
- a principal level, a trivialization at that level, the ordered root lines
  and reflection matrices, finite Coxeter orders, and exact tags for parallel
  and ultraparallel pairs;
- the determinant-line/Heegner data and bounded Borcherds data: weight,
  principal-part support and coefficients, safe correction support, and the
  divisor identity; and
- either a canonical projective gauge extracted uniquely from the marked
  incidence data, or the quotient by the ambient projective group, together
  with the oriented chamber and distinguished real embedding/signature
  label.  A freely chosen frame is not allowed: its choices would introduce
  a positive-dimensional torsor rather than rigidify the moduli problem.

Let

\[
  \mathscr H^{\rm alg}_{n,\mathcal D}/\mathbb Q
\]

be the algebraic stack of those data before imposing real inequalities.  It
is constructed from Hilbert and flag-Hilbert schemes, locally closed
incidence/contact conditions, finite level/root covers, and the finite
dimensional spaces of the bounded automorphic data.  The package
\(\mathcal D\) does **not** contain a list of coefficient or Vinberg fields:
those are residue/subfields output by the moduli calculation.  Let

\[
  \mathscr H^{\dagger}_{n,\mathcal D}
\]

denote its **real enhancement**: a geometric point together with an embedding
of its residue field into \(\mathbb R\), the distinguished Lorentzian
signature, an orientation, and the semialgebraic Coxeter/Vinberg return
certificate.  This distinction is mandatory.  A distinguished archimedean
embedding and positivity are not themselves algebraic markings on a stack
over \(\mathbb Q\).

At every geometric point, the splitting/eigencharacter marking evaluates the
unit-normal geometric Gram coordinates as scalars in
\(\overline{\mathbb Q}\).  A rational
point which merely remembers an abstract algebra such as
\(\mathbb Q[t]/(t^2-5)\), without choosing either scalar eigencharacter, is
not sufficient for field recovery and is deliberately excluded from the
rigid splitting cover.

One finite congruence quotient cannot certify that a Coxeter product has
infinite order.  Thus the exact Gram entries and the parallel/ultraparallel
tags (or a terminating exact lift certificate) are part of the object; the
level marking alone is not an exact-return theorem.

Write \(\mathscr H^{\rm rig}_{n,\mathcal D}\) for the locus on which the
canonical gauge (or projective quotient), level, roots, and order data kill
every automorphism.  If a surface automorphism can act trivially on those
data, the package must be strengthened or that point is not in the rigid
locus.

There are three different output fields.  On the fully scalar-decorated
splitting cover one has

\[
 k_V\subseteq K_{\rm geom}\subseteq M^\sharp.
\]

Here \(k_V\) is Vinberg's cyclic-product ground field,
\(K_{\rm geom}\) is the field of the outward unit-normal symmetric geometric
Gram coordinates, and \(M^\sharp\) is a field of definition of the complete
carrier plus its scalar splitting decoration.  The field of an undecorated
carrier need not contain \(K_{\rm geom}\).  None of these fields is the
unquotiented field of arbitrary primitive-root scalings or asymmetric Cartan
entries: independent root rescaling can vary those while fixing all
reflections and every cyclic product.

Let
\[
 \gamma:\mathscr H^{\rm rig}_{n,\mathcal D}\longrightarrow\mathscr G_{n,\mathcal D}
\]
forget the carrier but retain the outward unit-normal geometric Gram and its
cyclic products, and let \(Z_\gamma\) be its scheme-theoretic image.  One may
equivalently target Cartan data modulo the diagonal root-rescaling torus.

## 2. Two finite-dessin field theorems

> **Weak finite-Gram-image theorem.**  If \(Z_\gamma\) is zero-dimensional
> and finite type over \(\mathbb Q\), then the returned unit-normal geometric
> Gram fields \(K_{\rm geom}\), and hence their Vinberg subfields \(k_V\), form finite effectively
> enumerable lists.  Their degrees are bounded by the geometric point count
> of \(Z_\gamma\), or by the length of an explicit finite coordinate algebra.
> No zero-dimensionality of the carrier stack, trivial-stabilizer
> hypothesis, or descent theorem is needed.

Indeed, a zero-dimensional finite-type image is finite.  On the scalar
splitting cover, its finitely many geometric points literally evaluate the
Gram coordinates.  Exact number-field arithmetic then computes the
cyclic-product subfield at each point.  Positive-dimensional fibres of
\(\gamma\) are harmless for Coxeter-field effectivity: they may parametrize
many surfaces returning the same marked Gram.

> **Strong full-dessin theorem.**  Suppose the exact-return locus
> \(\mathscr H^{\rm rig}_{n,\mathcal D}\) is a zero-dimensional separated
> Deligne--Mumford stack of finite type over \(\mathbb Q\), with trivial
> geometric stabilizers.  Let
> \[
>   N=\#\mathscr H^{\rm rig}_{n,\mathcal D}(\overline{\mathbb Q})
> \]
> count its geometric points (or use the length of an explicit finite atlas
> as a computable upper bound).  Then:
>
> 1. its geometric points form a finite continuous \(G_{\mathbb Q}\)-set;
> 2. the field of moduli \(M_x\) of every point is a number field and
>    \([M_x:\mathbb Q]=\#(G_{\mathbb Q}x)\le N\);
> 3. the field of moduli is a field of definition; and
> 4. the fully decorated carrier fields \(M_x^\sharp\), unit-normal geometric
>    Gram fields \(K_{\rm geom}\),
>    and Vinberg fields \(k_V\) form finite effectively enumerable lists,
>    each of degree at most \(N\).

**Proof.**  A zero-dimensional finite-type scheme over a field is finite.
The same conclusion for geometric isomorphism classes follows from a finite
type atlas of the stated DM stack.  Hence \(G_{\mathbb Q}\) acts on a finite
set.  The stabilizer of \(x\) has finite index equal to the orbit size, and
its fixed field is \(M_x\).  Because the object has no automorphisms, the
isomorphisms between its stabilizer-conjugates are unique; uniqueness makes
the cocycle condition automatic.  Effective descent for the moduli stack
then descends the object to \(M_x\).

For ordered roots put

\[
 a_{ij}=\frac{2(r_i,r_j)}{(r_i,r_i)}.
\]

The Vinberg ground field is generated by the rescaling-invariant cyclic
products

\[
 k_V(x)=\mathbb Q\!\left(
   a_{i_1i_2}a_{i_2i_3}\cdots a_{i_si_1}
                         \right).
\]

The scalar splitting/root/Gram marking makes these elements evaluable in the
residue field of \(x\), so \(k_V(x)\subseteq M_x\).  Applying this to the finitely
many Galois orbits proves finiteness and the degree bound.  Given explicit
zero-dimensional equations, primary decomposition and number-field
factorization compute the residue fields, and exact arithmetic in them
computes the cyclic-product subfields.  Real-root isolation then selects the
embedding with signature \((1,n)\) and checks the distinguished Coxeter
inequalities.  This is a dessin-native enumeration of fields; no volume
formula or maximal-reflection-group catalogue enters the proof.  \(\square\)

Without trivial stabilizers, finiteness of fields of **moduli** still follows,
but effective descent can be obstructed by the residual gerbe.  This is why
full rigidification is a theorem hypothesis rather than decoration.

## 3. A surface criterion for the strong hypothesis

Bounded Hilbert polynomials put the embedded surfaces, marked curves,
zero-cycles, and cluster flags in a finite-type Hilbert/flag-Hilbert ambient
space.  They do **not** make that space finite: plane curves of fixed degree
and polarized K3 surfaces are elementary positive-dimensional
counterexamples.

A sufficient surface gate is the following.

> **Surface rigidity criterion.**  Assume the complete marked
> Hilbert/incidence/Borcherds locus attached to \(\mathcal D\) is finite type,
> its coefficient-order/lattice/local-system grammar is represented by the
> stated finite-type parameter scheme, and a canonical marked-data gauge or
> projective quotient kills its automorphisms.
> If the tangent complex of the **full**
> marked incidence problem has zero degree-one tangent space at every
> geometric exact-return point, then that locus is zero-dimensional and the
> finite-dessin field theorem applies.

On explicit charts this is a full-rank Jacobian calculation after removing
projective gauge directions.  On a log resolution where the complete
boundary/cluster marking is SNC, the familiar sufficient vanishing

\[
 H^1\!\left(S,T_S(-\log D')\right)=0
\]

may be used, provided the tangent conditions for the exceptional cluster,
marked sections, root local system, and Borcherds divisor are included as
well.  The universally correct object is the tangent complex/Jacobian of the
whole incidence scheme.  The unresolved surface marking can put three or
more curves through an ordinary face point, so a bare
\(T_S(-\log D)\) formula on that model is not generally valid.

The assertion must hold at **every** point of the exhaustive locus.  Showing
that one constructed surface is rigid does not rule out another component.
Likewise, a principal-level local system has only finitely many reductions,
but its integral lifts need not be finite.  They must be variables in the
bounded coefficient-order grammar, and the full incidence equations must cut
that lift space to dimension zero.  Prelisting finitely many
\(\mathcal O_k\)-lattices would make the field theorem circular.

## 4. Intrinsic divisor-lattice bonus, and its limitation

There is a separate unconditional surface-native bound when the arithmetic
lattice is genuinely algebraic.  If an \(\mathcal O_k\)-lattice of rank
\(n+1\) embeds faithfully in \(\operatorname{NS}(S)\), then

\[
  [k:\mathbb Q](n+1)\le \rho(S).
\]

For a complex K3 surface, \(\rho(S)\le20\).  Hence a degree-greater-than-one
intrinsic K3 divisor carrier is impossible for \(n\ge10\).  This is genuine
surface effectiveness and is independent of Prasad/ABSW.

It does **not** apply to the universal decorated Bertini surface, where the
quadratic root system is retained as an external local system rather than
embedded in \(\operatorname{NS}(S)\).  Treating an external marking as a
divisor-lattice embedding would manufacture a false dimension theorem.

### A conditional finite-Gram bridge from intrinsic surfaces

The intrinsic bound can be combined with two more genuinely effective
coordinates.  Put \(d_{\max}=\lfloor\rho(S)/(n+1)\rfloor\).  For a finite
Coxeter angle \(\pi/m\),

\[
 a_{ij}a_{ji}=4\cos^2(\pi/m),\qquad
 [\mathbb Q(4\cos^2(\pi/m)):\mathbb Q]=\varphi(m)/2
\]

for \(m\ge3\), with \(m=2\) separate.  Since this cyclic product lies in
\(k_V\subseteq k\), every finite dihedral label satisfies
\(\varphi(m)\le2d_{\max}\).  Hence the finite-angle alphabet

\[
 \mathcal M(d_{\max})=
 \{2\}\cup\{m\ge3:\varphi(m)\le2d_{\max}\}
\]

is finite and effectively enumerable.  This does not bound dotted
ultraparallel entries; it only bounds finite Coxeter labels.

Next suppose the facet sheets are separated and reduced on a projective
carrier with a fixed integral ample polarization \(H\).  If

\[
 \deg_H D_{\rm fac}
   =H^{s-1}\cdot\sum_{i=1}^{f}D_i\le\Delta ,
\]

then every nonempty distinct \(D_i\) has positive integral degree, so
\(f\le\Delta\).  On a stack or with a rational polarization one must instead
prove a uniform positive lower bound for the degree of one separated sheet;
without fixed normalization the claimed facet bound is invalid.

For \(n\ge3\), these two bounds give a finite unit-normal geometric Gram
image.  There are only
finitely many marked Coxeter matrices with at most \(\Delta\) nodes and
entries in \(\mathcal M(d_{\max})\cup\{\infty\}\).  A marked matrix has at
most one complete finite-volume hyperbolic reflection realization up to
isometry: apply Mostow--Prasad rigidity to the associated finite-covolume
reflection orbifolds (or their orientation-preserving index-two lattices).
This includes cusped \(n=3\); the cusp-completeness/finite-volume equations
must remain in the exact-return locus.  It fails in \(n=2\), where fixed
Coxeter data can have Teichmüller parameters.

Consequently, for an exhaustive intrinsic carrier class satisfying the
rank, normalized facet-degree, separation, and exact-return hypotheses, the
accepted unit-normal geometric Gram image and its
\(k_V,K_{\rm geom}\) field lists are finite.  Arbitrary primitive-root and
asymmetric Cartan coordinate fields are not asserted finite unless a
canonical bounded normalization is added.  The
census is effective in principle: enumerate the finite matrices and use
exact real algebraic elimination plus the literal Vinberg
finite-volume/cusp test.  Mostow gives uniqueness, not permission to omit
that existence test.

For a compact intrinsic K3 carrier in \(H^9\), \(k\ne\mathbb Q\) and
\[
  [k:\mathbb Q]\,10\le\rho\le20
\]
force \(d=2\), \(\rho=20\), and a zero-dimensional singular-K3 marked period
for each fixed primitive rank-twenty Neron--Severi lattice.  The union over
all such lattices and discriminants is not finite merely from \(\rho=20\).
Moreover
\[
 \mathcal M(2)=\{2,3,4,5,6,8,10,12\}.
\]
There is no such degree-greater-than-one intrinsic K3 carrier in \(H^{10}\)
or above.  These are carrier-relative conclusions, not global
nonexistence theorems for Coxeter polytopes.

## 5. Calibrations

### The `H11/26` line--quintic point

The raw chamber produces the exact incidence equation

\[
 a^2-b^2=0.
\]

The branch \(a=-b\) is nonreduced.  On \(a=b\ne0\), projective normalization
gives \(a=b=1\), the explicit line--quintic model over \(\mathbb Q\), and an
exact return of the same 26-wall chamber.  Thus the normalized accepted chart
contains a \(\mathbb Q\)-point and the recovered Vinberg field is
\(\mathbb Q\).  This input-specific computation does not by itself prove
trivial residual automorphisms or membership in the global rigid locus.

This does not yet prove that a variable-class stack
\(\mathscr H_{11,\mathcal D}\) is finite: the current calculation exhausts
the contraction endpoints of one input polytope, not every carrier chart for
every input in a bounded class.

### Bugaenko `H8/11`

The natural exact-Picard K3 carrier has a two-dimensional period family.  It
therefore fails the **strong full-carrier** zero-dimensional gate despite its
bounded Hilbert data.  It is not a counterexample to the weak field theorem:
the entire family can have constant Gram image and constant Vinberg field.
The separately chosen CM period and its coherent real multiplication pick an
algebraic carrier point only after extra Hodge/CM marking.  That extra
marking is irrelevant to the already fixed Bugaenko Gram field
\(\mathbb Q(\sqrt5)\), but it would have to be uniformly P-forced to support
the stronger carrier-field/descent theorem.

## 6. Exact logical ledger

Unconditional:

1. zero-dimensional finite type over \(\mathbb Q\) gives a finite Galois set;
2. orbit size bounds the field-of-moduli degree;
3. trivial automorphisms plus stack descent make the field of moduli a field
   of definition;
4. exact Gram/root markings recover the Vinberg field from cyclic products;
5. bounded Hilbert data gives finite type, and zero tangent everywhere gives
   zero dimension; and
6. a genuine \(\mathcal O_k\)-divisor-lattice embedding gives
   \([k:\mathbb Q](n+1)\le\rho(S)\).

Already sufficient for polytope fields, but not carrier fields: a
zero-dimensional scheme-theoretic image of the unit-normal geometric
Gram/cyclic-product
map.  This permits positive-dimensional surface fibres and requires no
carrier descent.

Still conditional for a general arithmetic Coxeter class:

1. a uniform, polytope-computable finite package \(\mathcal D\);
2. essential surjectivity: every target polytope has a point of the declared
   carrier stack;
3. exact return: every accepted point returns precisely its input chamber,
   including the literal Vinberg finite-volume test;
4. zero-dimensionality of at least the unit-normal geometric Gram image for
   \(k_V,K_{\rm geom}\)-effectivity, and of the entire exact-return locus only for
   carrier-field finiteness;
5. complete rigidification; and
6. explicit equations sufficient to run the enumeration.

## 7. Completion algorithm

1. From the P-only class, bound the Hilbert polynomials, number and degrees of
   marked curves/clusters, lattice types, root denominators, separating level,
   and Borcherds weight/support.
2. Build every Hilbert/flag-Hilbert and automorphic chart and impose the full
   incidence, divisor, reflection, cusp-cluster, and exact-return equations.
3. Quotient by projective gauge, or extract a unique gauge from forced marked
   data; then add full level/root ordering, exact Gram lift, and real
   signature/orientation labels, and prove essential surjectivity and exact
   return.
4. First compute the scheme-theoretic image of the unit-normal geometric Gram
   map (or the root-rescaling quotient).  A zero-dimensional image already
   gives finite \(k_V,K_{\rm geom}\) even if surface
   fibres move.  Compute the full Jacobian/tangent complex on every stratum
   only when finiteness/descent of the carrier fields \(M\) is also sought.
5. For the relevant zero-dimensional image or full locus, compute its reduced coordinate algebras,
   residue fields, and Galois orbits.  The algebra length is a certified upper
   bound for every field-of-moduli degree.
6. Recover \(k_V\) from cyclic Gram products, isolate all real embeddings, and
   run the exact Vinberg face/cusp/volume test at the distinguished one.
7. Output the finite list of fields and polytopes; emptiness is a
   nonexistence theorem, while surviving points are explicit examples.

This is the requested alternative to importing a finite field list from
Prasad's arithmetic-volume machinery or ABSW: **make the exact unit-normal geometric Gram image of the Coxeter--Borcherds dessin
finite, then read the polytope fields directly from its Galois orbits.**
Making the whole carrier stack finite is a stronger optional theorem.
