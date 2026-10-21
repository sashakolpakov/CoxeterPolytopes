# A variable surface-free Enriques lift criterion

**Status.**  This is an exact existence criterion for complete finite-volume
nef chambers of complex Enriques surfaces.  It is a structural reformulation
of the Enriques period theorem, the K3 Torelli theorem, and Namikawa's
no-hidden-curve proposition.  It does not assert a new eighth Enriques
chamber: the finite-volume criterion below, together with the
Kondō--Martin classification, implies that every successful datum is one of
the seven already classified types.

The point of the result is different.  It identifies precisely what has to be
added to a surface-free Coxeter Gram matrix in order to make the proposed
``polytope if and only if surface'' principle true, and it separates the
finite lattice computation from the period-theoretic existence argument.

## 1. The fixed K3-cover lattice

Put

\[
 L:=U\oplus E_8(-1),\qquad \operatorname{sig}(L)=(1,9).
\]

This is the numerical lattice of every complex Enriques surface.  Fix the K3
lattice in the form

\[
 \Lambda=U\oplus L\oplus L
\]

and the involution

\[
 \tau(u,x,y)=(-u,y,x).
\]

Its invariant and anti-invariant lattices are

\[
 M=\Lambda^\tau=\{(0,x,x):x\in L\}\cong L(2)
   \cong U(2)\oplus E_8(-2),
\]

and

\[
 N=M^\perp=\{(u,x,-x):u\in U,\ x\in L\}
   \cong U\oplus U(2)\oplus E_8(-2).
\]

Thus `M` has signature `(1,9)`, `N` has signature `(2,10)`, and
`M+N` has index \(2^{10}\) in \(\Lambda\).  Equivalently, \(\Lambda\) is
the graph overlattice for the canonical anti-isometry between the
2-elementary discriminant forms of `M` and `N`.

For a downstairs primitive root \(r_i\in L\), \(r_i^2=-2\), write

\[
 u_i=(0,r_i,r_i)\in M,\qquad u_i^2=-4.
\]

A proposed lift is specified by \(v_i\in N\) with \(v_i^2=-4\) such
that

\[
 b_i^+=\frac{u_i+v_i}{2}\in\Lambda.
\]

Then also \(b_i^-=(u_i-v_i)/2\in\Lambda\), and automatically

\[
 (b_i^+)^2=(b_i^-)^2=-2,
 \qquad b_i^+\cdot b_i^-=0,
 \qquad \tau(b_i^+)=b_i^-.
\]

The integrality condition is exactly the discriminant-glue condition.  It
also implies that \(v_i\) is a primitive even-type vector of square `-4`
and divisibility two.  This is Namikawa's even `(-4)` vector, not an
optional decoration.

## 2. Why the downstairs Gram does not determine the lift

For two distinct walls put

\[
 k_{ij}=r_i\cdot r_j,
 \quad a_{ij}=b_i^+\cdot b_j^+=b_i^-\cdot b_j^-,
 \quad c_{ij}=b_i^+\cdot b_j^-=b_i^-\cdot b_j^+.
\]

If all four lift classes are represented by distinct irreducible curves,
then \(a_{ij},c_{ij}\) are nonnegative integers.  Since
\(u_i=b_i^++b_i^-\),

\[
 a_{ij}+c_{ij}=k_{ij},\qquad
 v_i\cdot v_j=2(a_{ij}-c_{ij})=4a_{ij}-2k_{ij}.
\]

Consequently, for a fixed downstairs wall Gram,

\[
 v_i\cdot v_j\in
 \{-2k_{ij},-2k_{ij}+4,\ldots,2k_{ij}-4,2k_{ij}\}.
\]

This is a finite alphabet.  In particular,

```text
k=0:  0
k=1: -2 or 2
k=2: -4, 0, or 4.
```

Changing the name of the two sheets over the `i`-th wall sends
\(v_i\mapsto-v_i\).  Hence anti-invariant Gram candidates should be
enumerated modulo simultaneous row-and-column sign changes.  The choices
around a cycle need not be removable by these sign changes; this is genuine
global lift data.

The tempting choice that places every lift on the same fixed sheet is
usually wrong.  For each of the seven complete Enriques wall systems in the
pinned Shimada data, that choice gives

\[
 \operatorname{Gram}(v_i)=2\operatorname{Gram}(r_i),
\]

of signature `(1,9)` on its rank-ten span.  Its orthogonal complement in
`N` has only one positive direction, whereas a K3 period needs two.  Thus
every individual wall lifts, but the canonical simultaneous lift is period
obstructed in all seven cases.  This exact obstruction is replayed by
[`certify_enriques_variable_lift_package.py`](../scripts/certify_enriques_variable_lift_package.py).

## 3. Surface-free Coxeter input

Let

\[
 R=(r_1,\ldots,r_m),\qquad r_i\in L,\quad r_i^2=-2,
\]

be a finite ordered set of distinct primitive roots.  Fix a component
\(\mathcal C^+\) of the positive cone and define

\[
 C_R=\{x\in L_{\mathbb R}:x\cdot r_i\geq0\text{ for all }i\}.
\]

The **finite-volume Coxeter gate** means all of the following.

1. `C_R` is pointed, full-dimensional, meets \(\mathcal C^+\), and every
   displayed inequality is irredundant.
2. Distinct roots have nonnegative integral intersection.  Thus the only
   finite Coxeter angles are `m=2` and `m=3`; intersection two is parallel
   and larger intersection is ultraparallel.
3. The entire homogeneous cone lies on the selected time sheet:

   \[
   C_R\subset\overline{\mathcal C^+}.
   \]

   For a rational polyhedral cone this is certified by exact double
   description: every extreme ray is finite or ideal, none is hyperideal,
   and all causal rays pair positively with one strict interior witness.
4. Independently, the literal Vinberg elliptic-edge and full affine
   zero-neighborhood tests certify finite volume.

The direct cone condition in item 3 is the hypothesis used in Namikawa's
no-hidden-curve theorem.  The Vinberg replay in item 4 prevents a finite
vertex sample or an incomplete cusp list from being mistaken for a volume
proof.

Choose once and for all an integral strict interior vector
\(h\in C_R\cap L\), and let \(H=(0,h,h)\in M\).  Such an `h` is obtained
by exact rational feasibility followed by clearing denominators.

## 4. A good lift package

A **good Enriques lift package** for `R` consists of vectors
\(v_1,\ldots,v_m\in N\) satisfying the following finite conditions.

### L1. Simultaneous glue and sheet intersections

For every `i`,

\[
 v_i^2=-4,\qquad b_i^\pm=(u_i\pm v_i)/2\in\Lambda.
\]

For every two distinct signed lifts, the pairing is nonnegative.  Equivalently,
the same-sheet/cross-sheet numbers above are nonnegative integers summing to
\(k_{ij}\).

### L2. A negative, rootless anti-invariant span

Let

\[
 A=\operatorname{Sat}_N\langle v_1,\ldots,v_m\rangle.
\]

Require `A` to be negative definite and to contain no vector of square
`-2`.  In particular \(\operatorname{rk}A\leq10\).

The use of the saturated span, rather than an arbitrary larger
anti-invariant Picard lattice, is sufficient for an existence theorem.  If
one starts from a carrier with additional anti-invariant algebraic classes,
one can deform its period to a very general point with exact algebraic
anti-invariant lattice `A`; Section 6 explains why the selected curves remain
irreducible.

### L3. The generic Picard lattice and a regular invariant chamber

Put

\[
 S_A=\Lambda\cap(M_{\mathbb Q}\oplus A_{\mathbb Q}).
\]

This is the primitive closure of the algebraic invariant and
anti-invariant data in the K3 lattice.  Require

\[
 H^2>0,\qquad H\cdot b_i^\pm>0,
\]

and require `H` to be root-regular in `S_A`:

\[
 \{e\in S_A:e^2=-2,\ H\cdot e=0\}=\varnothing.
\]

The last check is finite because \(H^\perp\cap S_A\) is negative definite.

### L4. Finite irreducibility, not merely effectivity

For every \(b=b_i^\pm\), require that there is no root \(e\in S_A\) with

\[
 e^2=-2,\qquad 0<H\cdot e<H\cdot b,
 \qquad b\cdot e<0. \tag{*}
\]

This is also a finite enumeration.  The degree \(H\cdot e\) belongs to a
fixed finite interval, and after projecting to the negative-definite lattice
\(H^\perp\) the equation \(e^2=-2\) gives only finitely many lattice points.

Condition `(*)` is the exact negative-component test for irreducibility; it
must not be replaced by a two-root decomposition test.  For example, a
reducible root can have the form \(b=e+F\), where \(e^2=-2\), \(F^2=0\),
and \(e\cdot F=0\).  Then \(b\cdot e=-2\), so `(*)` detects the rational
component even though the complementary class is isotropic rather than a
root.  More generally, if an effective divisor of square `-2` is reducible,
at least one irreducible component `E` has negative intersection with the
whole divisor: otherwise summing the nonnegative component intersections
would give \(b^2\geq0\).  Such an `E` is a smooth rational curve and obeys
`(*)`.  Conversely, an `H`-positive root `e` satisfying `(*)` is effective;
its negative intersection with an effective representative of `b` forces a
common component, and the strict degree inequality rules out irreducible
`b`.

This condition cannot be dropped: a primitive effective `(-2)` class can be
a reducible ADE tree, or can contain an isotropic residual divisor as above.
An even `(-4)` lift guarantees an effective root with a two-sheet
decomposition; it does not by itself say that either half is an irreducible
curve.

## 5. The exact criterion

### Theorem (variable Enriques lift criterion)

Let `R` pass the finite-volume Coxeter gate of Section 3.  Then the
following are equivalent.

1. There is a numerically marked complex Enriques surface `Y` whose complete
   set of smooth rational curves has numerical classes exactly
   \(r_1,\ldots,r_m\).  Consequently the projectivized full nef cone is the
   finite-volume Coxeter polytope defined by `R`.
2. `R` admits a good Enriques lift package satisfying `L1`--`L4`.

Everything in item 2 is stated in the fixed abstract lattices
\(L,M,N,\Lambda\); no surface occurs in the input.

### Proof: a good package produces the surface

Let \(a=\operatorname{rk}A\).  The orthogonal complement
\(A_N^\perp\) has signature `(2,10-a)`.  Its type-IV period domain

\[
 \Omega_A=
 \{[\omega]\in\mathbb P(A_N^\perp\otimes\mathbb C):
   \omega^2=0,\ \omega\cdot\bar\omega>0\}
\]

has complex dimension \(10-a\).

If \(a<10\), choose a point outside every rational hyperplane belonging to
an element of \(N\setminus A\).  These are countably many proper analytic
hyperplanes, so their complement is dense and nonempty.  If \(a=10\), the
domain is zero-dimensional; its oriented positive two-plane is already the
whole real complement of `A`, so an integral element orthogonal to the period
belongs to `A`.  In both cases one obtains

\[
 N\cap\omega^\perp=A,
 \qquad
 \operatorname{NS}(X)=
 \Lambda\cap(M_{\mathbb Q}\oplus A_{\mathbb Q})=S_A
\]

for a marked K3 surface supplied by period surjectivity.  Rootlessness of
`A` says precisely that no anti-invariant algebraic `(-2)` class is forced.

Because `H` is positive and root-regular, it lies in the interior of one
root chamber of the positive cone of `S_A`.  Use the refined marked K3
period theorem, equivalently re-mark by the K3 Weyl group, so that this is
the Kähler chamber and `H` is ample.  The lattice involution `tau` is a Hodge
isometry, fixes `H`, and therefore preserves that chamber.  Strong Torelli
realizes it by an automorphism \(\iota\) of `X`.

Ohashi's free-involution criterion now applies:

- the invariant lattice is the primitive lattice
  \(M\cong U(2)\oplus E_8(-2)\);
- no root of `NS(X)` is orthogonal to `M`, because
  \(\operatorname{NS}(X)\cap N=A\) and `A` is rootless; and
- `M` contains the ample class `H`.

Hence \(\iota\) is fixed-point-free and \(Y=X/\iota\) is an Enriques
surface.

Each \(b_i^\pm\) is a root of `NS(X)` with positive `H`-degree, so K3
Riemann--Roch makes it effective.  If an effective divisor of class
\(b_i^\pm\) were reducible, one of its rational components `E` would satisfy

\[
 0<H\cdot E<H\cdot b_i^\pm,
 \qquad b_i^\pm\cdot E<0,
\]

contradicting `L4`.  Conversely, a root satisfying `(*)` is effective and
an effective representative shares a component with the divisor of class
`b`; the strict `H`-degree inequality makes that impossible when `b` is an
irreducible curve.  Thus `L4` is exact rather than merely sufficient, and
the two roots are represented by disjoint smooth rational curves exchanged
by \(\iota\).  Their quotient is a smooth rational curve on `Y` of class
`r_i`.

It remains to exclude unlisted curves.  Namikawa's Proposition 6.9 applies
verbatim: distinct smooth rational curves with inequality cone contained in
\(\overline{\mathcal C^+}\) are all the smooth rational curves on the
surface.  Indeed, any additional nodal class would pair nonnegatively with
all selected curves and hence lie in `C_R`, contradicting its square `-2`.
Therefore the selected inequalities are the complete nef chamber.  The
independent cone and Vinberg gates give its finite hyperbolic volume.

### Proof: a surface produces a good package

Conversely, let `Y` have the stated complete finite-volume nef chamber, and
let \(X\to Y\) be its K3 cover.  A nodal curve \(C_i\) splits as two
disjoint smooth rational curves \(C_i^+,C_i^-\), exchanged by the covering
involution.  Under an Enriques marking put

\[
 v_i=[C_i^+]-[C_i^-]\in N.
\]

Then `L1` and the finite two-sheet formula follow immediately.  The
saturated span `A` is negative definite because it lies in the orthogonal
complement of the positive period plane; it is rootless because a free
involution admits no anti-invariant algebraic root.

Initially the anti-invariant algebraic lattice may be larger than `A`.
Move to a very general period in \(\Omega_A\), so the new Picard lattice is
exactly `S_A`.  It is a sublattice of the initial Picard lattice containing
all \(b_i^\pm\).  A new decomposition of one of these roots would already
have been a decomposition in the initial Picard lattice, where the class was
represented by an irreducible curve.  Hence `L4` survives the deformation.
An integral ample class in the strict interior supplies `H`, and `L3`
follows.  This gives a good package.

## 6. What is finite, what is continuous, and what is not yet an implementation

For a fixed wall Gram, the same-sheet/cross-sheet formula gives finitely many
anti-invariant Gram matrices.  For each one:

1. quotient by its rational radical and enumerate its finitely many even
   saturated overlattices;
2. test negative definiteness, rank at most ten, and rootlessness;
3. enumerate primitive embeddings in `N` with the prescribed mod-two glue
   classes \([v_i/2]\leftrightarrow[u_i/2]\); Nikulin's discriminant-form
   description and finiteness of lattice class numbers make this a finite
   arithmetic problem up to `O(N)`;
4. construct `S_A`, test `H`-regularity, and perform the bounded root
   enumeration in `L4`; and
5. replay exact cone and literal Vinberg volume tests.

Thus a *proposed* lift package has a finite exact verifier, and a fixed
surface-free wall datum has a finite theoretical compiler up to the standard
primitive-embedding/class-group enumeration.  That compiler has not yet been
implemented end to end in this repository; the present machine certificate
freezes the ambient lattice, the lift identities, the finite intersection
alphabet, and the canonical-lift obstruction.

The continuous part is softer.  For a surviving `A`, periods form a type-IV
domain of dimension \(10-\operatorname{rk}A\).  Requiring `A` to be contained
in the algebraic anti-invariant lattice gives an orthogonal modular special
cycle.  Requiring it to be *exactly* the algebraic anti-invariant lattice
removes a countable collection of Noether--Lefschetz hyperplanes.  This
very-general locus is dense and nonempty, but in a marked domain it is not in
general a single finite Zariski-open condition.  Effectivity and
irreducibility are restored by adjoining the Kähler-chamber witness `H` and
the finite root test `L4`.

Two endpoints deserve explicit mention.

- `rank(A)=0` is the ordinary very-general, unnodal Enriques locus.  It can
  carry no nonempty proposed nodal wall set and its full hyperbolic chamber
  has infinite volume.
- `rank(A)=10` gives a zero-dimensional period domain.  The Baire-category
  avoidance argument is replaced by the elementary positive-plane argument
  in the proof above; it is not an exception to the criterion.

## 7. Consequences and limitations

The criterion explains exactly why the one-wall slogan

```text
downstairs (-2) root + compatible anti-invariant (-4) vector
```

is necessary but not sufficient.  It omits simultaneous negative
definiteness, sheet-intersection data, period rootlessness, Kähler
orientation, irreducibility, and completeness.

For full-nef finite-volume Enriques chambers, the criterion does not open an
unclassified infinite source of examples.  Kondō's Torelli bridge and
Martin's classification imply

\[
 \operatorname{vol}(P(Y))<\infty
 \Longleftrightarrow
 [O^+(L):W_Y]<\infty
 \Longleftrightarrow
 |\operatorname{Aut}(Y)|<\infty,
\]

so every success is one of Types I--VII.  The variable lift theorem is
nevertheless useful in three ways:

1. it is an intrinsic, non-lookup explanation of what a surface carrier
   certificate contains;
2. failures of its finite stages give genuine carrier obstructions without
   confusing them with abstract-polytope nonexistence; and
3. the same package architecture can be adapted to K3, quotient-nef, or
   higher-dimensional hyperkähler sections where no seven-type exhaustion is
   available.

## 8. Primary-source anchors

- W. Barth and C. Peters, *Automorphisms of Enriques Surfaces*, Invent.
  Math. 73 (1983), 383--411,
  [author PDF](https://www-fourier.ujf-grenoble.fr/~peters/Articles/pub11.pdf):
  Section 1.2 gives the fixed `Lambda,M,N` eigenspaces and discriminant
  glue; Lemma 2.4 identifies invariant ample cones; Sections 2.4--2.5 give
  the period complement and exchanged nodal pairs.
- Yukihiko Namikawa, *Periods of Enriques Surfaces*, Math. Ann. 270
  (1985), 201--222,
  [doi:10.1007/BF01456182](https://doi.org/10.1007/BF01456182):
  Proposition 2.16 (even `(-4)` vectors), Proposition
  4.5 and Remark 4.6 (K3 lifts and the distinction between roots and
  irreducible curves), Proposition 6.2 (the nodal Heegner divisor), and
  Proposition 6.9 (the no-hidden-curve cone criterion).
- Valery Gritsenko and Klaus Hulek, *Moduli of polarized Enriques
  surfaces*, 2015, [arXiv:1502.02723](https://arxiv.org/abs/1502.02723),
  Sections 2--3: the lattices `M,N`, the period image
  \(\mathcal D_N\setminus\mathcal H_{-2}\), and the even `(-4)` nodal
  Heegner divisor.
- Hisanori Ohashi, *On the Number of Enriques Quotients of a K3 Surface*,
  Publ. RIMS 43 (2007), 181--200,
  [doi:10.2977/PRIMS/1199403814](https://doi.org/10.2977/PRIMS/1199403814),
  Proposition 2.2: the primitive
  `U(2)+E8(-2)` embedding, no orthogonal root, and invariant ample-class
  criterion for a free involution.
- Shigeyuki Kondō, *Enriques surfaces with finite automorphism groups*,
  Japan. J. Math. 12 (1986), 191--282, Section 1.3 and Theorem 1.5:
  the nodal reflection chamber and the Torelli finite-index bridge.
- Gebhard Martin, *Enriques surfaces with finite automorphism group in
  positive characteristic*, Algebraic Geometry 6 (2019), 592--649,
  introduction and Main Theorem: over the complex numbers, finite reflection
  index is equivalent to finite automorphism group, and the complete nodal
  graph is one of seven types.

## 9. Replayable finite artifacts

- [finite lift certifier](../scripts/certify_enriques_variable_lift_package.py)
- [exact lift certificate](../data/enumerations/enriques_variable_lift_package.json)
- [regression tests](../tests/test_enriques_variable_lift_package.py)
- [full-nef finite-volume classification](enriques_full_nef_finite_volume_classification.md)
- [seven-image exact theorem](enriques_finite_automorphism_nef_transfer_theorem.md)
