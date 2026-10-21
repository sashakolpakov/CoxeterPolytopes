# Hostile scope audit: automatic Galois trace envelopes and supersingular Mukai gates

**Verdict (2026-08-12): `AUTO-TRACE` PASS; universal mixed-K3 compiler
OPEN.**  This audit freezes and rechecks

```text
results/automatic_galois_trace_envelope_and_supersingular_mukai_gates.md
81cb7b966c31f1b1767dd46a0e3dddfa152689f2680cbaaf9e0bab98429c36ee
```

Two conclusions must be kept separate.

First, for an arbitrary nondegenerate rational action lattice, rational side
isometries preserving one exact arithmetic Coxeter slice force a rational
nondegenerate invariant envelope

\[
 S_{\mathbf Q}\simeq\operatorname {Res}_{k/\mathbf Q}V,
 \qquad \dim_{\mathbf Q}S=d(n+1),                     \tag{A.1}
\]

with transfer form \(\operatorname {Tr}_{k/\mathbf Q}(\beta q)\) and an
order in \(k\) acting on its integral lattice.  This general conclusion does
not determine the signs of all definite nondistinguished real blocks.

Second, if the ambient positive index is at most two and \(n\geq2\), every
nondistinguished block must be negative.  Only with this additional
K3/Mukai hypothesis does one obtain

\[
 \operatorname {sig}S=(1,d(n+1)-1),\qquad d(n+1)\le23.           \tag{A.2}
\]

The owner theorem does **not** establish `EMBED`, `GLUE`, `CRYSTAL`,
`EXACT`, `PRIME`, or boundedness.  No frozen absolutely reflection-maximal
arithmetic chamber of trace rank at least 24 was found.  Thus (A.2) is not
currently an explicit counterexample to the whole proposed species.

## 1. Exact hypothesis ledger

The general envelope theorem assumes:

1. \(P\subset\mathbf H^n\) is irreducible, finite-volume, arithmetic, and
   \(n\geq2\);
2. its standard reflection representation has minimal Vinberg field \(k\),
   degree \(d\), and dimension \(m=n+1\);
3. \((A,B)\) is an integral lattice with nondegenerate rational form;
4. the marked side operators \(g_i\in O(A)\) are rational integral
   isometries; and
5. they preserve a common nondegenerate real \(m\)-plane \(E\), on which
   they are the actual standard side reflections and \(B|_E\) is a positive
   multiple of the distinguished arithmetic form.

Positive index at most two is **not** needed for (A.1).  It is an additional
hypothesis used only for the signature conclusion (A.2).

The ambient \(g_i\) need not satisfy Coxeter relations.  Treating them as
images of free generators is enough, so `AUTO-TRACE` covers weak
kernel-by-kernel Fourier--Mukai side lifts at the numerical level.  It does
not cover unrelated external real wall matrices attached as decorations to
a rational object.  Without hypotheses 4--5, the Gram can be carried as
unbounded external data.

## 2. Character-field replay

For

\[
 P_i=e_i\otimes\frac{2(-,e_i)}{(e_i,e_i)},\qquad s_i=I-P_i,
\]

direct multiplication of rank-one operators gives, for every ordered index
sequence,

\[
 \operatorname {tr}(P_{i_1}\cdots P_{i_r})
 =a_{i_1i_2}\cdots a_{i_ri_1}.                        \tag{A.3}
\]

Expanding \(\prod_j(I-P_{i_j})\) expresses its trace as the sum over ordered
subsequences.  The full length-\(r\) cyclic monomial is the unique
length-\(r\) term.  Induction gives

\[
 \mathbf Q[\text{Cartan cycles}]
 \subseteq\mathbf Q[\operatorname {tr}\rho(w)].       \tag{A.4}
\]

Every term in the reverse expansion is a Cartan cycle, so the reverse
inclusion holds.  Vinberg's cycle-field construction identifies the
fraction field with the minimal Vinberg field:

\[
       k=\mathbf Q(\operatorname {tr}\rho(w):w\in W(P)).         \tag{A.5}
\]

This uses the cycle-field construction, not a Vinberg, Nikulin, Esselmann,
or ABSW reflection-finiteness theorem.

The standard representation is absolutely irreducible: the
orientation-preserving subgroup is a finite-covolume lattice in
\(SO^+(1,n)\), so Borel density reduces the assertion to irreducibility of
the standard complex orthogonal representation.  Borel density is a
representation-theoretic input, not a reflection-finiteness import.

## 3. The transcendental-slice attack fails only under rational actions

Choosing \(E\) with a transcendental basis does not evade the theorem.  The
intertwining equations

\[
                         g_i\phi=\phi s_i              \tag{A.6}
\]

have algebraic coefficients.  Their solution space is defined over the real
algebraic numbers.  Nondegeneracy of \(\phi^*B\) is Zariski open and the
distinguished signature is Euclidean open.  A real solution therefore gives
a real algebraic solution in the same open set.  Absolute irreducibility
makes its pulled-back invariant form a scalar multiple of \(q\).

This algebraizes a representation copy.  It does not algebraize an arbitrary
externally supplied real Gram.  Dropping rational ambient actions drops the
theorem.

## 4. Larger coefficient fields create multiplicity, not new trace types

The algebraic intertwiner need not be defined over \(k\).  Choose a finite
Galois extension \(L/\mathbf Q\) containing the normal closure of \(k\) and
all coefficients of the intertwiner.  Its image is a simple
\(L\)-submodule whose character is valued in the distinguished copy of
\(k\).

For \(\gamma\in\operatorname {Gal}(L/\mathbf Q)\), the conjugate submodule
has character \(\gamma\circ\chi\).  By (A.5), two such simple modules are
isomorphic exactly when the corresponding automorphisms have the same
restriction to \(k\).  Thus the distinct character types are indexed by the
\(d\) embeddings of \(k\); extra conjugates over the same embedding only
increase a multiplicity.

Vinberg's \(k\)-model gives Schur index one.  Equivalently, restriction of
scalars gives a rational module whose complexification is the multiplicity-
one sum of the \(d\) conjugate standard modules.  Galois transitivity and
absolute irreducibility make it the rational simple module

\[
                       M=\operatorname {Res}_{k/\mathbf Q}V,
 \qquad \dim_{\mathbf Q}M=dm.                         \tag{A.7}
\]

Hence a coefficient field larger than \(k\) cannot increase the *minimal*
forced rank and cannot suppress any of the \(d\) character types.

## 5. Rational descent has no hidden semisimplicity assumption

Let \(\mathcal R\subset\operatorname {End}_{\mathbf Q}(A_{\mathbf Q})\) be
the finite-dimensional unital algebra generated by the \(g_i\).  The full
ambient free-group representation need not be semisimple.  The proof uses
the socle instead.

The Jacobson radical of a finite-dimensional algebra commutes with the
separable extension \(L/\mathbf Q\), so

\[
 \operatorname {soc}_{\mathcal R}(A_{\mathbf Q})\otimes_{\mathbf Q}L
 =\operatorname {soc}_{\mathcal R_L}(A_L).             \tag{A.8}
\]

The algebraized distinguished simple submodule lies in the right side.  Its
full Galois orbit therefore descends to a nonzero \(M\)-isotypic part of the
rational socle.  Write that part as

\[
       \operatorname {Res}_{k/\mathbf Q}(V\otimes_kU).           \tag{A.9}
\]

Schur's lemma and trace duality put the restricted symmetric form in the
shape

\[
 B(v\otimes u,w\otimes u')
 =\operatorname {Tr}_{k/\mathbf Q}
       \bigl(q(v,w)h(u,u')\bigr),                     \tag{A.10}
\]

where \(h\) is symmetric and \(k\)-bilinear.  The adjoint involution on the
scalar commutant \(k\) is the identity, so no hidden quaternionic or unitary
multiplicity algebra occurs.

The original nondegenerate real copy supplies \(u_0\) with
\(h_{\sigma_0}(u_0,u_0)\) nonzero and of the required distinguished sign.
Weak approximation supplies \(u\in U(k)\) in that real open set and off the
quadric \(h(u,u)=0\).  The line \(ku\) gives a nondegenerate rational copy of
exact rank \(dm\).

This closes three candidate gaps: ambient nonsemisimplicity, a coefficient
field larger than \(k\), and multiplicity greater than one.

## 6. Signature replay and the role of positive index

Different real embeddings are orthogonal because a cross-pairing would give
an intertwiner \(V_\sigma\to V_\tau^*\), while these self-dual constituents
have different characters.  The distinguished block has signature
\((1,n)\); every other arithmetic conjugate is definite.  This much is
unconditional.

Now assume the ambient positive index is at most two.  A positive
nondistinguished block would contribute \(m=n+1\geq3\) positive directions,
in addition to the one already used by the distinguished block.  This is
impossible, so all nondistinguished blocks are negative definite and

\[
                         \operatorname {sig}M=(1,dm-1).          \tag{A.11}
\]

Without the positive-index hypothesis, (A.1), the transfer form, and block
definiteness survive, but (A.11) does not.  This split is explicit in the
frozen theorem.

Dimension one is also excluded for a reason.  When \(m=2\), one positive
definite conjugate uses exactly the one spare Mukai positive direction and
is not ruled out by index two.  The frozen dimension-one infinite-field
family is therefore not contradicted.

For a complex or finite-height algebraic K3 Mukai lattice, negative index at
most 20 gives \(dm\le21\).  For a supersingular numerical Mukai lattice of
signature \((2,22)\), (A.11) gives the sharp necessary bound \(dm\le23\).
Nonprimitive embeddings, common rescaling, Brauer twists, or projected side
vectors do not alter this inertia count once an exact rational action slice
exists.

## 7. `EMBED/GLUE` remains a genuine universal gap

For a supersingular K3 of characteristic \(p\) and Artin invariant
\(\sigma\),

\[
 \widetilde N_{p,\sigma}=U\oplus N_{p,\sigma},\qquad
 \operatorname {sig}=(2,22),\qquad |\det|=p^{2\sigma}.            \tag{A.12}
\]

If a rank-\(R\) trace lattice \(T\) embeds primitively and
\(C=T^\perp\), then for every \(\ell\ne p\) the ambient lattice is
unimodular over \(\mathbf Z_\ell\).  Discriminant gluing forces

\[
 (A_{T,\ell},q_{T,\ell})\simeq
 (A_{C,\ell},-q_{C,\ell}),\qquad
 \ell_\ell(A_T)\le24-R.                               \tag{A.13}
\]

At \(R=22\), every non-\(p\) discriminant part must be at most
two-generated.  Rank and real signature do not imply this.

Even when \(T\hookrightarrow\widetilde N_{p,\sigma}\) exists, a side action
extends only if it preserves the primitive gluing datum.  For the isotropic
gluing subgroup \(H\subset A_T\oplus A_C\), one needs

\[
                         (\bar s_i,\bar c_i)H=H         \tag{A.14}
\]

for some complement isometry \(c_i\).  Neither arithmeticity nor reflection
maximality supplies (A.13)--(A.14).  No universal `EMBED/GLUE` theorem
survived this audit.

## 8. `CRYSTAL/FM` and strict actions remain separate

A supersingular numerical lattice is not a marked K3 crystal.  A lift must
preserve the strictly characteristic period subspace, Frobenius data,
orientation, and the relevant cone.  Ogus and derived crystalline Torelli
are lifting criteria after compatible marked data are present; they do not
turn every numerical isometry into an automorphism or Fourier--Mukai
equivalence.

The general derived crystalline result cited in the owner report assumes
\(p\ne2\) unless separately specified.  The \(H^{21}\) characteristic-two
control is valid because its 168 long-side lifts are explicit automorphisms
and its 42 short-side lifts are explicit spherical twists.  This isolated
construction is not a characteristic-two universal lifting theorem.

A spherical twist has the required numerical reflection but

\[
                            T_S^r(S)=S[-r],             \tag{A.15}
\]

so it is not a functorial involution.  A strict homomorphism
\(W(P)\to\operatorname {Aut}_{FM}(D^bX)\) therefore has an additional
relation gate.  Derived-isogeny chains can relax an individual integral
extension, but do not give a simultaneous action on one cone or uniform
bounds on primes, Brauer orders, denominators, intermediate surfaces, or
chain length.  No universal `CRYSTAL/FM` bridge survived this audit.

## 9. `EXACT`, primes, and boundedness remain separate

Supersingular specialization can add divisor classes and new
\((-2)\)-walls.  Thus even a specialized action can return a strict
subdivision of the original chamber.  Torelli does not make the new walls
redundant.  There is no universal no-hidden-wall theorem here.

Allowing arbitrary characteristic also destroys the quasi-compactness
needed by the proposed finiteness argument: the disjoint union of
supersingular strata over all primes has infinitely many open-and-closed
components.  Fixing \(p=2\) removes that choice but leaves the unproved
universal embedding, gluing, and characteristic-two lifting claims.

On a fixed K3, the kernels

\[
                  \mathcal O_\Delta\otimes p_2^*L^r
\]

already have unbounded numerical and Hilbert class.  Finite-type K3 moduli
therefore do not make the decorated action stack finite type.  No universal
`EXACT/PRIME/BOUND` bridge survived this audit.

## 10. Search for an absolute-maximal rank-24 obstruction

The frozen endpoint controls are:

| raw chamber | trace rank | status relevant here |
|---|---:|---|
| Vinberg \(X_3\subset\mathbf H^{19}\) | 20 | complex K3 weak side lifts; not a rank obstruction |
| Borcherds \(H^{21}\) | 22 | absolutely reflection-maximal; characteristic-two weak side lifts |

The repository contains rank-24 trace calculations in proposed OG10 or
totally-real extension packages.  They are not a frozen raw absolutely
reflection-maximal arithmetic Coxeter chamber in the target class.  High-
trace finite-volume nonarithmetic Coxeter examples fail arithmeticity.

The honest search result is therefore

\[
 \boxed{
 \begin{array}{l}
 \text{no universal `EMBED/GLUE/CRYSTAL/EXACT` bridge is proved;}\\
 \text{no frozen absolute-maximal arithmetic chamber with }d(n+1)\ge24
 \text{ is available.}
 \end{array}}                                           \tag{A.16}
\]

Such a chamber would immediately refute the bounded-rank mixed-K3 action
species by (A.11), without a classical reflection-finiteness theorem.  Until
one is frozen or every remaining bridge is constructed, `AUTO-TRACE` is a
sharp necessary theorem, not an AG-native global finiteness theorem.

## 11. Final claim ledger

| claim | audit status |
|---|---|
| character field equals the cycle/Vinberg field | **proved** |
| transcendental basis evades Galois closure | **false under the rational-action contract** |
| a larger intertwiner field suppresses or enlarges the minimal trace type | **false; it creates multiplicity only** |
| rational action forces a rank-\(d(n+1)\) trace submodule | **proved** |
| arbitrary rational action lattice forces signature \((1,d(n+1)-1)\) | **false without a positive-index bound** |
| index at most two and \(n\ge2\) force that signature | **proved** |
| supersingular Mukai capacity \(d(n+1)\le23\) | **proved, carrier-relative** |
| all trace lattices embed into some supersingular Mukai lattice | **not proved** |
| all side actions extend through discriminant glue | **not proved** |
| numerical actions automatically preserve one crystalline period | **false as an automatic inference** |
| weak numerical reflections give a strict Coxeter action by kernels | **not proved** |
| supersingular reduction preserves the exact chamber | **not proved** |
| the decorated mixed-characteristic carrier stack is bounded | **not proved** |
| an absolute-maximal arithmetic rank-24 counterexample is frozen | **no** |
| global AG-native finiteness follows | **no** |

