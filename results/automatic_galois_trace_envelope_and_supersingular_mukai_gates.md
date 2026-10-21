# Automatic Galois trace envelopes and the supersingular Mukai gates

**Status (2026-08-12).**  This note proves a new necessary theorem for any
K3/Mukai **action** carrier of an arithmetic Coxeter chamber.  One does not
have to put the full restriction-of-scalars lattice into the definition of
such a carrier: once rational cohomological isometries preserve a common
Coxeter slice, rationality of cohomology forces all of its Galois conjugates.
No separate algebraicity assumption on the chosen real basis is needed.  In
a mixed-characteristic K3 numerical Mukai lattice this
gives the sharp inequality

\[
                 [k:\mathbf Q](n+1)\leq 23.             \tag{0.1}
\]

Thus a proposal which stores only the distinguished real slice does **not**
evade the full-trace capacity gate.  This conclusion does not use a Vinberg,
Nikulin, Esselmann, or ABSW dimension/finiteness theorem.

The note then gives the exact lattice, crystalline-period, prime, kernel, and
hidden-wall gates that remain before (0.1) can be turned into a universal
carrier theorem.  Those gates are not consequences of reflection maximality.
In particular, reduction modulo a prime does not prove essential
surjectivity.  No reflection-maximal arithmetic chamber of trace rank at
least 24 is currently frozen in this repository, so the result is a sharp
necessary theorem rather than an explicit counterexample to the whole mixed
K3 species.

## 1. The action-carrier contract

Let \(P\subset\mathbf H^n\), with \(n\geq2\), be an irreducible
finite-volume arithmetic Coxeter chamber.  Put

\[
 m=n+1,\qquad d=[k:\mathbf Q],
\]

where \(k\) is its minimal Vinberg field.  Let

\[
 \rho_{\sigma _0}:W(P)\longrightarrow O(V_{\sigma _0},q_{\sigma _0})
\]

be the standard reflection representation at the distinguished embedding.
The arithmetic \(k\)-form has signature \((1,n)\) there and is definite of
the opposite sign at every other real embedding.

An **algebraic cohomological action carrier** for the purpose of this note
consists of

1. an integral lattice \((A,B)\) with nondegenerate rational form;
2. rational integral isometries \(g_i\in O(A)\), one for every marked side;
3. a common \(g_i\)-invariant nondegenerate real subspace
   \(E\subset A_{\mathbf R}\), of dimension \(m\), on which
   \(g_i\) is the corresponding side reflection and on which \(B\) is one
   positive multiple of \(q_{\sigma _0}\); and
4. no external real operators: the \(g_i\) above are the actual rational
   cohomological actions used for the return, rather than unrelated real
   matrices written only on \(E\).

There is no additional algebraicity hypothesis on a basis of \(E\).  The
intertwiner space

\[
 \operatorname {Hom}_{\langle s_i\rangle}
     (V_{\sigma _0,\mathbf C},A_{\mathbf C})
\]

is the solution space of linear equations with algebraic coefficients.  The
condition that the pulled-back form be nondegenerate is Zariski open, and
the required real signature is Euclidean open.  Since the given real
intertwiner lies in it, the same open set contains a real algebraic
intertwiner.  Thus the representation copy can be algebraized before taking
Galois conjugates.  By contrast, if the \(g_i\) themselves are allowed to be
external real matrices, the complete Gram can indeed be smuggled in as a
decoration and the theorem does not apply.

The general Galois-envelope statement below uses only this contract.  Its
signature sharpening, and hence its K3 capacity corollary, additionally
assume that the positive index of \((A,B)\) is at most two.  This is true for
the algebraic numerical Mukai lattice of a K3 surface or K3 category, but it
is not a formal property of an arbitrary integral action lattice.

The \(g_i\) need not define a strict functor-level action of \(W(P)\) on a
derived category.  It is enough that they preserve the common slice and have
the prescribed restrictions.  Formally one may regard them as a rational
representation of the free group on the sides.  Hence the theorem below also
applies to the weak, kernel-by-kernel version of a Fourier--Mukai carrier.

## 2. The character field is the Vinberg field

Choose side normals \(e_i\), and write

\[
 a_{ij}=\frac{2(e_i,e_j)}{(e_i,e_i)},\qquad
 P_i=e_i\otimes\left(x\mapsto\frac{2(x,e_i)}{(e_i,e_i)}\right).
\]

Then \(s_i=I-P_i\), and for every finite index sequence
\(i_1,\ldots,i_r\),

\[
 \operatorname {tr}(P_{i_1}\cdots P_{i_r})
   =a_{i_1i_2}a_{i_2i_3}\cdots a_{i_ri_1}.             \tag{2.1}
\]

Expanding a reflection word gives the exact identity

\[
\begin{split}
 \operatorname {tr}(s_{i_1}\cdots s_{i_r})
  =m+\sum_{\varnothing\ne\{j_1<\cdots<j_t\}\subseteq\{1,\ldots,r\}}
  (-1)^t
  a_{i_{j_1}i_{j_2}}\cdots a_{i_{j_t}i_{j_1}} .        \tag{2.2}
\end{split}
\]

The length-\(r\) cyclic product is the unique length-\(r\) term in (2.2).
Induction on \(r\) therefore recovers every directed Gram/Cartan cyclic
product from characters of reflection words.  Conversely, Vinberg's
construction supplies an \(m\)-dimensional model over the field generated by
these cycles.  Consequently

\[
 \boxed{
 k=\mathbf Q\bigl(\operatorname {tr}\rho(w):w\in W(P)\bigr).}    \tag{2.3}
\]

This is the character-theoretic form of the usual cycle-field definition;
it also shows why the abstract Coxeter matrix alone is insufficient when
dotted distances occur.

The standard representation is absolutely irreducible.  One way to see this
is that the orientation-preserving subgroup is a finite-covolume lattice in
\(SO(1,n)\), hence is Zariski dense, and the standard representation of the
complex orthogonal group is irreducible.  Equivalently, one can use the
connected nondegenerate geometric representation.

## 3. Automatic Galois-envelope theorem

### Theorem 3.1 (the full trace cost is forced)

Suppose the action-carrier contract of Section 1 holds.  Then \(A_{\mathbf Q}\)
contains a rational nondegenerate invariant subspace \(S_{\mathbf Q}\) with

\[
                 \dim_{\mathbf Q}S_{\mathbf Q}=dm.      \tag{3.1}
\]

Up to rational isometry its invariant form is a transfer

\[
             \operatorname {Tr}_{k/\mathbf Q}(\beta q)  \tag{3.2}
\]

for a scalar \(\beta\in k^*\) having the required distinguished sign.
Moreover \(S_{\mathbf Q}\cap A\) is stable under an order in \(k\).
Its distinguished real block has signature \((1,n)\), while every other
real block is definite.  If, in addition, the positive index of \(A\) is at
most two, then

\[
              \operatorname {sig}(B|_{S_{\mathbf R}})=(1,dm-1). \tag{3.3}
\]

#### Proof

By the intertwiner argument in Section 1, replace \(E\), if necessary, by a
real algebraic nondegenerate copy with the same distinguished signature.  Its
representation is absolutely irreducible.  Choose a finite Galois extension
\(L/\mathbf Q\) containing the normal closure of \(k\) and the coefficients
of this algebraic intertwiner.  The resulting simple \(L\)-submodule has
character in the distinguished copy of \(k\), even if the submodule itself
is defined only over the larger field \(L\).

Since the ambient matrices \(g_i\) are rational, applying
\(\operatorname {Gal}(L/\mathbf Q)\) produces a constituent for every real
embedding \(\sigma:k\hookrightarrow\mathbf R\).  Equation (2.3) makes the
constituents belonging to different embeddings pairwise nonisomorphic: an
isomorphism would identify their characters and hence the two embeddings of
\(k\).  Conjugates by elements having the same restriction to \(k\) have the
same character and can only increase the multiplicity of that simple type.
Thus coefficients of the intertwiner outside \(k\) create multiplicity, not
new character types.

The standard representation has Schur index one because Vinberg's
construction realizes it over \(k\).  Hence the simple rational module
belonging to this Galois orbit is
\(M=\operatorname {Res}_{k/\mathbf Q}V\).

For completeness, this descent does not assume that the whole ambient
free-group representation is semisimple.  Let \(\mathcal R\) be the
finite-dimensional unital \(\mathbf Q\)-algebra in
\(\operatorname {End}_{\mathbf Q}(A_{\mathbf Q})\) generated by the \(g_i\).
The socle of the \(\mathcal R\)-module \(A_{\mathbf Q}\) commutes with the
separable scalar extension \(L/\mathbf Q\).  The algebraic distinguished
simple submodule lies in the \(L\)-socle, so its full Galois orbit descends
to a nonzero \(M\)-isotypic part of the rational socle.  It has the form

\[
       \operatorname {Res}_{k/\mathbf Q}(V\otimes_k U)
\]

for a finite-dimensional multiplicity space \(U/k\).  On this part the
restriction of \(B\) is

\[
  (v\otimes u,w\otimes u')\longmapsto
  \operatorname {Tr}_{k/\mathbf Q}\bigl(q(v,w)h(u,u')\bigr)       \tag{3.2a}
\]

for a symmetric \(k\)-bilinear form \(h\) on \(U\).  This follows from
Schur's lemma; the adjoint involution is the identity on the scalar
commutant \(k\).  The original nondegenerate distinguished copy gives a real
multiplicity vector \(u_0\) with \(h_{\sigma _0}(u_0,u_0)\) nonzero and of
the required sign.  Weak approximation chooses \(u\in U(k)\) in that real
open set and off the quadric \(h(u,u)=0\).  The line \(ku\) therefore gives
one rational nondegenerate copy of \(M\).  This justifies the assertion that
the minimal envelope has rank exactly \(dm\), rather than only a multiple
of \(dm\), without a hidden complete-reducibility hypothesis.

Trace duality and Schur's lemma identify every rational invariant symmetric
form on \(M\) with (3.2).  Over the reals it decomposes orthogonally as

\[
 M_{\mathbf R}=\bigoplus_{\sigma:k\hookrightarrow\mathbf R}V_\sigma.
                                                               \tag{3.4}
\]

Indeed, a pairing between two different summands is a homomorphism
\(V_\sigma\to V_\tau^*\), hence vanishes because the standard
representations are self-dual and have different characters.  At
\(\sigma_0\) the chosen restriction has signature \((1,n)\).  At every
other embedding it is a nonzero real scalar multiple of the definite
arithmetic conjugate \(q_\sigma\).  This proves the unconditional block
statement.  Now assume that the positive index of \(A\) is at most two.  If
one of those scalars reversed the negative sign, its block would contribute
\(m=n+1\geq3\) new positive directions.  The distinguished block already
uses one positive direction, so this is impossible.  Every
nondistinguished block is therefore negative definite, proving (3.3).

Finally, the rational commutant contains \(k\);
the elements preserving the lattice \(S_{\mathbf Q}\cap A\) form an order in
\(k\).  \(\square\)

### Corollary 3.2 (sharp K3 capacities without assuming `TRACE`)

For a complex or finite-height projective K3, or a characteristic-zero K3
category, the algebraic Mukai lattice has negative index at most 20.  Thus

\[
                         d(n+1)\leq21.                  \tag{3.5}
\]

For a supersingular K3 numerical Mukai lattice the signature is
\((2,22)\), so

\[
                         \boxed{d(n+1)\leq23}.          \tag{3.6}
\]

Since \(n\ge2\), this would give the explicit carrier-relative bounds

\[
 d\le7,\qquad n\le\left\lfloor\frac{23}{d}\right\rfloor-1.       \tag{3.7}
\]

Thus a universal mixed K3 action compiler would, before any moduli argument,
prove a field-degree-seven and dimension-22 theorem for the raw chambers.
These are consequences of carrier existence, not facts assumed about all
arithmetic reflection groups.

In particular, a rank-24 Lorentzian trace block still cannot fit: it has
signature \((1,23)\), not merely rank 24.  Twisting, nonprimitive
embeddings, common rescaling, projected side vectors, and replacing a
surface by a bounded-rank genuine K3 category do not change this inertia
calculation.

The content of Theorem 3.1 is stronger than the elementary implication
"a marked full trace injection has rank \(d(n+1)\)": here the full Galois
envelope follows from a rational cohomological **action on one algebraic
real slice**.  Therefore a nonrational carrier cannot avoid the rank bill by
omitting a marked \(k\)-operator from its definition.

## 4. What a supersingular numerical embedding actually requires

Let \(N_{p,\sigma}\) be the Neron--Severi lattice of a supersingular K3
surface of Artin invariant \(\sigma\).  Rudakov--Shafarevich give

\[
 \operatorname {rk}N_{p,\sigma}=22,
 \quad \operatorname {sig}N_{p,\sigma}=(1,21),
 \quad \det N_{p,\sigma}=-p^{2\sigma},                 \tag{4.1}
\]

and its discriminant group is \(p\)-elementary (type I when \(p=2\)).
The numerical Mukai lattice

\[
             \widetilde N_{p,\sigma}=U\oplus N_{p,\sigma}       \tag{4.2}
\]

has rank 24, signature \((2,22)\), and the same \(p\)-primary
discriminant module.

Let \(T\hookrightarrow\widetilde N_{p,\sigma}\) be a primitive trace
envelope of rank \(R\), and put \(C=T^\perp\), of rank \(24-R\).  For every
prime \(\ell\ne p\), the ambient lattice is unimodular over
\(\mathbf Z_\ell\).  Primitive orthogonal gluing then gives

\[
 (A_{T,\ell},q_{T,\ell})\simeq
 (A_{C,\ell},-q_{C,\ell}),
 \qquad
 \boxed{\ell_\ell(A_T)\leq24-R}.                       \tag{4.3}
\]

Here \(\ell_\ell\) is the minimum number of generators of the
\(\ell\)-primary discriminant group.  At the rank-22 endpoint every
non-\(p\) part of the induced trace lattice must therefore be at most
two-generated.  Rank and real signature alone do not imply (4.3).

If \(H\subset A_T\oplus A_C\) is the isotropic gluing subgroup defining
the primitive extension, a side action \(s_i\in O(T)\) extends integrally
only if some \(c_i\in O(C)\) satisfies

\[
                         (\bar s_i,\bar c_i)H=H.         \tag{4.4}
\]

Thus an embedding of rational quadratic spaces does not construct the
integral side actions.  Passing to a finite-index lattice or multiplying the
form can move the obstruction, but does not make (4.3)--(4.4) automatic.

## 5. The crystalline period and Fourier--Mukai gates

The numerical lattice is not the supersingular K3 crystal.  For odd
characteristic, a marked supersingular K3 surface determines a strictly
characteristic isotropic subspace

\[
 K\subset N_0\otimes\overline{\mathbf F}_p,
 \qquad \dim K=\sigma,                                  \tag{5.1}
\]

in the \(2\sigma\)-dimensional discriminant quadratic space, with the
Frobenius incidence condition
\(\dim(K\cap F(K))=\sigma-1\).  Ogus's Torelli theorem identifies geometric
automorphisms with the appropriate lattice isometries that preserve this
period and the ample cone.  In particular, even at Artin invariant one the
period is not permission to use every abstract numerical isometry.

For twisted supersingular K3 surfaces, Bragg's derived crystalline Torelli
theorem similarly characterizes Fourier--Mukai equivalences by
orientation-preserving isomorphisms of **marked twisted Mukai crystals**.
It is a lifting theorem after the crystal isomorphism has been constructed;
it does not say that an arbitrary isometry of (4.2) preserves the crystal.
Consequently every proposed side kernel must pass

\[
 \boxed{\text{integral glue }(4.4)
   +\text{ characteristic-period preservation }
   +\text{ orientation}.}                              \tag{5.2}
\]

Brauer twists change the integral twisted Mukai lattice and the
characteristic subspace, so they may repair a particular gluing problem.
They do not change rank or real inertia, and they do not prove that one
period is simultaneously invariant under all side actions.

The general crystalline/derived Torelli results used above are stated in
odd characteristic (Bragg's paper assumes \(p\ne2\) unless otherwise
specified).  The Borcherds \(H^{21}\) control in characteristic two is valid
because its 168 long actions are explicit automorphisms and its 42 short
actions are explicit spherical twists; it is not evidence for a class-wide
characteristic-two lifting theorem.

Allowing **derived isogenies** instead of autoequivalences does not close the
route.  Bragg--Yang show that chains of twisted Fourier--Mukai equivalences
realize a broad class of rational cohomological isogenies in positive
characteristic.  This can relax integral gluing for an individual rational
isometry, but it starts only after compatible etale/crystalline data have
been supplied; the intermediate surfaces, Brauer orders, denominators, and
lengths of the chains are not uniformly bounded.  Moreover a correspondence
between different surfaces is not an action on one fixed cone.  If arbitrary
rational isogenies are admitted without bounds, the carrier returns to the
unrestricted rational-Chow species and loses the finite-type leverage.

Finally, a spherical twist realizes a root reflection only on numerical
cohomology.  If \(S\) is 2-spherical, then

\[
                  T_S(S)=S[-1],\qquad T_S^r(S)=S[-r],  \tag{5.3}
\]

so \(T_S\) is not an involution.  Individual numerical side lifts therefore
do not by themselves give a strict homomorphism
\(W(P)\to\operatorname {Aut}_{FM}(D^bX)\).  The automatic trace theorem does
not require this strict splitting, but a carrier contract which does require
it has an additional unsolved gate.

## 6. Why reduction modulo \(p\) is not a universal construction

Starting from a characteristic-zero K3 carrier and reducing it at a prime
does not supply the missing theorem.

1. **Supersingularity is not automatic.**  Specialization injects the old
   divisor lattice, but rank 22 occurs only at supersingular reductions.  No
   theorem chooses such a prime uniformly for an arbitrary K3 carrier with
   all its marked actions.
2. **Exact return is not stable under Picard jump.**  New divisor classes,
   especially new \((-2)\)-classes, may have mirrors cutting the selected
   slice.  The specialization can therefore return a strict subdivision of
   \(P\).
3. **Kernels have their own specialization condition.**  A numerical
   isometry on the special fibre must still pass (4.4) and (5.2).  Crystalline
   Torelli does not manufacture these properties from the generic fibre.
4. **The prime is a boundedness datum.**  The union

   \[
                 \coprod_p \mathcal M^{ss}_{p,\sigma}             \tag{6.1}
   \]

   over all characteristics is an infinite non-quasi-compact union.  Viewed
   inside mixed-characteristic K3 moduli, the supersingular-only locus is not
   a finite constructible stratification across all closed characteristics.
   A global bounded-carrier proof therefore needs a fixed finite set of
   primes, or a mixed-characteristic constructible extension on which exact
   extraction is controlled.

Thus "choose a convenient supersingular reduction" hides a new unbounded
choice and cannot be inserted into a quasi-compact AG finiteness argument.
Fixing characteristic two would remove this particular `PRIME` defect and
retain the \(H^{21}\) control, but it replaces it by two unproved class-wide
claims: that every forced trace envelope embeds into one of the ten
characteristic-two supersingular lattices with all side gluing, and that all
the resulting numerical actions have geometric or derived lifts.  Neither
is a consequence of the single Dolgachev--Kondo surface.

## 7. Boundedness still does not follow from K3 moduli

Even after `TRACE`, `EMBED`, `GLUE`, and `CRYSTAL` are solved, the decorated
carrier is not bounded merely because polarized K3 surfaces of fixed degree
form a finite-type moduli stack.  One must also bound

- the polarization degree and the number of side decorations;
- numerical classes and Hilbert polynomials of graphs/kernels;
- Tor amplitude or regularity of the Fourier--Mukai kernels;
- the finite set of allowed characteristics; and
- the exact no-hidden-wall predicate in families.

On one fixed K3, tensoring by \(L^r\) already gives diagonal kernels

\[
                 \mathcal O_\Delta\otimes p_2^*L^r     \tag{7.1}
\]

of unbounded numerical and Hilbert class.  Supersingular Picard jumps add
the further failure of local constancy noted in Section 6.  Hence the
decorated Fourier--Mukai stack is only a union of finite-type pieces until
these numerical data are bounded.

## 8. Exact verdict and controls

The automatic Galois-envelope theorem closes one ambiguity in the proposed
mixed K3 route:

\[
 \boxed{\text{rational cohomological side actions on a common exact slice}
 \Longrightarrow\text{ the full }d(n+1)\text{ trace cost}.}       \tag{8.1}
\]

It follows that universal essential surjectivity of this species would
already prove the new global arithmetic inequality (0.1).  That inequality
cannot be assumed in order to construct the carriers.

The two endpoint controls remain consistent with the theorem.

| chamber | \(d(n+1)\) | mixed K3 status |
|---|---:|---|
| Vinberg \(X_3\subset\mathbf H^{19}\) | 20 | complex K3; 24 curve twists and 12 further FM numerical side lifts |
| absolutely maximal Borcherds \(H^{21}\) | 22 | characteristic-two supersingular K3; 42 curve twists and 168 explicit automorphisms |

Neither has rank 24.  Neither proves universal embedding, simultaneous
crystalline preservation, bounded kernels, or exact hidden-wall exhaustion.
There is presently no frozen absolute-maximal arithmetic input of trace rank
at least 24 to refute the species by Corollary 3.2, and no theorem proving
that every absolute-maximal input satisfies the six remaining gates.

Accordingly the honest global contract is now

\[
\begin{array}{ll}
\texttt{AUTO-TRACE}&\text{proved by Theorem 3.1 after an algebraic action slice exists};\\
\texttt{CAPACITY}&d(n+1)\le23\text{ is then necessary};\\
\texttt{EMBED/GLUE}&\text{(4.3)--(4.4), not automatic};\\
\texttt{CRYSTAL/FM}&\text{(5.2), not automatic};\\
\texttt{EXACT}&\text{all new curve/derived mirrors are redundant};\\
\texttt{PRIME/BOUND}&\text{finite characteristics and bounded decorations};\\
\texttt{CONSTR}&\text{exact extraction is constant on finitely many strata}.
\end{array}                                                \tag{8.2}
\]

Only after all rows of (8.2) are proved uniformly would quasi-compact AG
boundedness finite-ize the returned Coxeter data.  The universal positive-
characteristic/mixed K3 theorem, and hence the requested global AG-native
finiteness theorem through this route, is not yet proved.

## 9. Primary sources and frozen dependencies

- Vinberg's cycle-field construction is reviewed with the explicit
  \(m\)-dimensional invariant \(k\)-space in Edoardo Dotti,
  [*On the commensurability of hyperbolic Coxeter groups*](https://arxiv.org/abs/2101.10024),
  especially the construction following the arithmeticity criterion.
- Kondo--Shimada,
  [*On a certain duality of Neron--Severi lattices of supersingular K3 surfaces*](https://www.algebraicgeometry.nl/2014-3/2014-3-016.pdf),
  Theorem 3.1 and Section 7: the Rudakov--Shafarevich lattice and Ogus
  characteristic-subspace/Torelli description.
- Daniel Bragg,
  [*Derived equivalences of twisted supersingular K3 surfaces*](https://arxiv.org/abs/1811.07379),
  Theorem 3.5.3: the marked, oriented crystalline lifting condition.
- Bragg--Yang,
  [*Twisted derived equivalences and isogenies between K3 surfaces in positive characteristic*](https://arxiv.org/abs/2102.01193):
  the rational derived-isogeny relaxation and its cohomological input.
- Dolgachev--Kondo,
  [*A supersingular K3 surface in characteristic 2 and the Leech lattice*](https://arxiv.org/abs/math/0112283):
  the explicit \(H^{21}\) control.
- [`mukai_crystalline_universal_carrier_capacity_and_lift_obstruction.md`](mukai_crystalline_universal_carrier_capacity_and_lift_obstruction.md)
  supplies the independent signature, discriminant, strict-action, and
  kernel-decoration audit which this note strengthens by `AUTO-TRACE`.
- [`rank4_k3_exact_maximal_eigen_interval_infinite_fields.md`](rank4_k3_exact_maximal_eigen_interval_infinite_fields.md)
  is a warning that even fixed polarization square, exact eigen-nef return,
  and relative reflection maximality do not finite-ize fields in hyperbolic
  dimension one; the present theorem is stated for \(n\ge2\).
