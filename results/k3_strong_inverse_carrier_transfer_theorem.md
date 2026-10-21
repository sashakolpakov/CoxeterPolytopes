# Strong inverse-carrier theorem for complete K3 root chambers

## Result

There is one non-tautological class for which the proposed inverse-carrier
program is already a theorem rather than a heuristic: complete chambers of
primitive square-minus-two roots in a K3-embeddable even hyperbolic lattice.

Let $N$ be an even lattice of signature $(1,n)$, $2\leq n\leq19$, and
let ${\cal C}_N^+$ be one component of its positive cone.  Write

\[
 \Delta(N)=\{r\in N:r\text{ primitive},\ r^2=-2\},
 \qquad
 W_2(N)=\langle s_r:r\in\Delta(N)\rangle.
\]

Let $C\subset\overline{{\cal C}_N^+}$ be the closure of a connected
component of

\[
 {\cal C}_N^+\setminus\bigcup_{r\in\Delta(N)}r^\perp.
\]

Then the following are equivalent.

1. $N$ admits a primitive embedding in the K3 lattice
   
   \[
   \Lambda_{K3}=U^3\oplus E_8(-1)^2.
   \]
2. There is a projective complex K3 surface $X$ and a
   Néron--Severi marking
   
   \[
   \eta:\operatorname{NS}(X)\xrightarrow{\sim}N
   \]
   
   for which $\eta(\operatorname{Nef}(X))=C$.

Consequently, whenever $C$ has finitely many irredundant walls, its
unit-hyperboloid slice (equivalently, the projectivization of its positive
rays) is a genuine finite-sided hyperbolic Coxeter polyhedron in
$\mathbb H^n$, carried exactly by the full nef cone of $X$.  Conversely,
every full-nef K3 Coxeter polytope arises in this way.  The lattice types
admitting **finite-volume** strong carriers are precisely the primitively
K3-embeddable even hyperbolic lattices satisfying
$[O^+(N):W_2(N)]<\infty$; individual marked carriers additionally depend on
an embedding orbit, a chamber marking, and a point of the moduli locus below.

This theorem is about Néron--Severi markings.  It does not assert that root
reflections are induced by automorphisms of $X$.

## Moduli statement

For each primitive-embedding orbit $i:N\hookrightarrow\Lambda_{K3}$, put
$T=i(N)^\perp$.  If one retains a full $\Lambda_{K3}$-marking, the periods
with exact Picard lattice $i(N)$ lie in

\[
 \Omega_T^\circ=
 \{[\omega]\in\mathbf P(T_\mathbf C):
   \omega^2=0,
   (\omega,\bar\omega)>0\}
 \setminus\bigcup_{v\in\Lambda_{K3}\setminus i(N)}v^\perp.
\]

This exact-Picard locus is nonempty and has complex dimension
$20-\operatorname{rank}(N)$.  Choose one type-IV component
$\Omega_T^{\circ,+}$, and let $\Gamma_i$ be the arithmetic group induced by
K3-lattice isometries that fix $i(N)$ pointwise and preserve that component.
Impose the discrete chamber condition that the induced Néron--Severi
marking sends the ample chamber to $C$; this adds no period equation.  After
then forgetting the full K3 marking while retaining the pointwise
$N$-marking, the corresponding coarse exact-Picard carrier locus is
$\Gamma_i\backslash\Omega_T^{\circ,+}$ (with the usual finite stabilizers
and, if desired, stack structure).  Hence one complete integral root datum
gives an entire positive-dimensional family unless
$\operatorname{rank}(N)=20$.

The existence decision is finite and exact: $N$ embeds primitively if and
only if there is an even lattice $T$ of signature
$(2,20-\operatorname{rank}(N))$ with discriminant form
$q_T\cong-q_N$.  Stable-range length inequalities are sufficient shortcuts;
at the boundary one must apply the full local discriminant-form existence
test.  This is the exact discrete lattice gate; the period-domain
construction above is the soft, continuous part of the transfer theorem.

## Why the chamber is the nef cone

Period surjectivity produces a projective K3 surface with
$\operatorname{NS}(X)\cong N$.  The K3 Weyl-chamber theorem identifies its
ample cone with one component of the complement of all square-minus-two root
hyperplanes.  The group $W_2(N)$ acts transitively on these components.
Composing the Néron--Severi marking with the required Weyl element therefore
maps the closed nef cone to the prescribed $C$.

For an ample class $h\in C^\circ$, Riemann--Roch chooses the effective sign
of every square-minus-two root.  The simple roots bounding $C$ are exactly
the classes of irreducible smooth rational curves that give the nef walls.
Completeness of the root arrangement is essential: a finite list of
algebraic divisor classes, without the assertion that $C$ is a component
of the *full* root complement, does not identify a nef cone.

The converse follows from Hodge index and the standard description of the
ample cone: $\operatorname{NS}(X)$ is an even primitive hyperbolic sublattice
of the K3 lattice, and every codimension-one wall of the full nef cone is
orthogonal to an irreducible square-minus-two curve.

## The surface-algebra/volume dictionary

For any carrier in the theorem, the following conditions are equivalent:

\[
\begin{aligned}
 &\operatorname{vol}_{\rm hyp}(C/\mathbf R_{>0})<\infty;\\
 &[O^+(N):W_2(N)]<\infty;\\
 &\operatorname{Aut}(X)\text{ is finite};\\
 &\operatorname{Eff}(X)\text{ is rational polyhedral};\\
 &\operatorname{Cox}(X)\text{ is finitely generated}.
\end{aligned}
\]

Indeed $O^+(N)$ is an arithmetic lattice, $C$ is a fundamental chamber
for $W_2(N)$, and a subgroup of a lattice has finite covolume exactly when
it has finite index.  Global Torelli and Sterk identify the finite-index
condition with finiteness of the automorphism group.  The K3 cone and
Cox-ring theorems identify this in turn with rational polyhedrality of the
effective cone and finite generation of the Cox ring.

In particular, if a proposed finite wall list has separately been promoted
to the *complete* K3 nef cone, finite generation of the surface's
multigraded section algebra proves finite hyperbolic volume without a
diagram-by-diagram volume computation.  Vinberg's endpoint and full
parabolic-neighborhood tests certify finite volume of the polyhedron cut out
by the *displayed* walls; they do **not** prove that no additional
square-minus-two root hyperplane cuts its interior.  Full-nef promotion
therefore needs an independent complete-root enumeration with a stopping
criterion, an exact hidden-root separation argument, or a published complete
chamber theorem.

## Compactness transfer

Under the equivalent finite-volume conditions,

\[
 C/\mathbf R_{>0}\text{ is compact}
 \iff N_\mathbf Q\text{ is anisotropic}
 \iff X\text{ has no genus-one fibration}.
\]

A primitive nef isotropic class is the fibre class of a genus-one fibration
and gives an ideal ray of the chamber.  Conversely, if $N_{\mathbf Q}$ is
isotropic, clear denominators and divide by the content to obtain a primitive
integral isotropic class.  Riemann--Roch chooses an effective sign, and a
finite sequence of square-minus-two Weyl reflections moves it into the nef
closure; the resulting primitive nef isotropic class is the fibre class of a
genus-one fibration.  Equivalently, an ideal ray of a rational polyhedral K3
nef cone is rational and supplies the same class.

Meyer's theorem makes every indefinite rational quadratic space of dimension
at least five isotropic.  Therefore

\[
 \boxed{n\geq4\quad\Longrightarrow\quad
 \text{no compact full-nef K3 Coxeter polytope in }\mathbb H^n.}
\]

This is an abstract-polytope obstruction only for the stated full-nef K3
carrier class.  It does not exclude compact eigen-nef sections, quotient
sections, mixed divisor-reflection systems, rational-surface constructions,
or polytopes with no such carrier.

## Exact inverse-carrier formulation

Define the Coxeter/lattice class ${\cal K}_n$ without mentioning a surface:

1. $N$ is even of signature $(1,n)$ and embeds primitively in
   $\Lambda_{K3}$;
2. the permitted mirrors are all primitive square-minus-two roots of $N$;
3. a stated finite irredundant simple-root set cuts out one complete component
   of the full root-hyperplane complement; and
4. the normalized Gram data pass the Coxeter-angle conditions.

Let ${\cal M}_{N,C}$ be the marked K3 moduli problem above.  Then

```text
(N,C) belongs to K_n
    iff
M_(N,C) is nonempty, its marked nef cone is C, and C has the stated finite
irredundant wall set.
```

Let ${\cal K}_n^{\mathrm{fv}}$ be the finite-volume subclass, equivalently
the members also satisfying $[O^+(N):W_2(N)]<\infty$.  For this subclass,
classifying the integral **marked carrier data** reduces to the following
lattice problem:

```text
primitive-embedding orbits + 2-reflective Néron--Severi lattices
    + a complete chamber/simple-root system (up to W_2(N) if Weyl-unmarked)
    --> finite-volume marked full-nef K3 Coxeter polyhedra.
```

This arrow is not injective after forgetting the integral marking.  More
generally, the passage from carrier points or families to normalized abstract
Coxeter polyhedra is many-to-one: different period points,
primitive-embedding orbits, integral lattices or overlattices can yield the
same normalized abstract polyhedron.  Thus the theorem is an
existence-and-classification transfer for marked carrier data, followed by a
generally many-to-one forgetful map to abstract polytopes.

The continuous period space proves existence of the surfaces.  The discrete
classification has not disappeared: it is concentrated into discriminant
forms, 2-reflective lattices, and complete simple-root systems.  This is the
correct division of labour between algebraic geometry and Vinberg theory.

## Promotion boundary and next classes

The theorem deliberately excludes the mixed-root carriers already abundant
in the ledger.  A primitive square-minus-four reflection can be an integral
Picard reflection and an algebraic divisor class, but it cannot be a wall of
the full K3 nef cone.  Such examples remain weak or ample-intersection
carriers until a different geometric cone is identified.

The next non-tautological inverse-carrier classes are therefore:

1. complete invariant nef sections under a finite surface action;
2. complete eigenspace sections for a self-adjoint algebraic correspondence,
   with a finite visible effective cone;
3. finite Mori-chamber systems on bounded rational surfaces; and
4. complete MBM chambers on projective hyperkähler varieties.

Each extension needs its own analogue of the complete-root condition before
surface classification can imply a polytope theorem.

## Supporting results

- [`k3_weak_algebraic_carrier_theorem_audit.md`](k3_weak_algebraic_carrier_theorem_audit.md)
  proves the primitive-embedding/weak-carrier equivalence and records the
  exact nonstable discriminant-form decision.
- [`k3_nef_covolume_dictionary.md`](k3_nef_covolume_dictionary.md) proves the
  finite-volume/finite-automorphism/Cox-ring dictionary and identifies cusps
  with genus-one fibrations.
- [`vinberg_divisor_volume_gate.md`](vinberg_divisor_volume_gate.md) states
  the independent literal finite-volume tests used for every frozen example.
