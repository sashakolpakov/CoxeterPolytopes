# Marked-surface compression for arithmetic Coxeter dessins

**Status (2026-08-09): repaired decorated-carrier theorem.** A hostile
audit found that the former single-point cusp statement was incompatible
with log resolution. The valid construction records a cusp by its complete
pre-resolution incidence, equivalently by the labeled exceptional boundary
cluster and its contraction to the Baily--Borel cusp. It also keeps the
original reflection stack and the marked root/descent data separate from the
neat surface. The construction is a noncanonical faithful association, not
a functor. It gives no classification pushback by itself.

## The repaired theorem

Let \(n\geq2\), and let

\[
 {\cal C}(P)=(\mathscr X,D,C,\mathbb V,\partial)
\]

be the orthogonal Coxeter dessin of a marked compact or finite-volume
arithmetic Coxeter polytope \(P\subset\mathbb H^n\). The input includes:

- the coefficient order and its fraction field \(k\), the quadratic module
  \((L,q)\), the marked root lines \(\{\pm r_i\}\), their distinguished
  reflections, and the chosen positive-cone/halfspace orientation;
- a finite, sheet-separated Heegner support and, after increasing the neat
  level if necessary, separated marked isotropic rays; and
- in the noncompact case, the completed rational-isotropic cusp list and
  the cusp-wall set \(I(c)\) at each ideal vertex \(c\).

Choose a neat finite cover, a projective toroidal compactification, and a
projective log resolution of the finite marked mirror-boundary arrangement.
After a finite extension \(K\) of the reflex field there is a smooth
projective surface \(S_P/K\), obtained by a complete intersection when
\(n>2\) and followed by finitely many blowups, with the following finite
decoration.

1. Each facet sheet gives a labeled curve \(C_i\subset S_P\).
2. Every face \(F\) whose relative interior lies in \(\mathbb H^n\) has a
   distinct marked point \(p_F\). Precisely the curves \(C_i\) with
   \(i\in I(F)\) pass through \(p_F\).
3. An ideal vertex is **not** represented by a purported common point after
   log resolution. It is represented by the labeled exceptional boundary
   cluster
   \[
      {\cal K}_c=
      (E_c,\Delta_c,\{\widetilde D_i:i\in I(c)\},
       \pi:E_c\longrightarrow c),
   \]
   where \(E_c\) is the reduced boundary/exceptional fiber over the
   Baily--Borel cusp, \(\Delta_c\) is the relevant labeled dual complex,
   and the attachment of every strict mirror transform is recorded. One
   chooses marked points in the open parts of the finitely many nonempty
   cluster strata needed to realize this incidence on \(S_P\). Equivalently,
   one may retain the pre-resolution record \((c,I(c))\) together with the
   resolution map. The contraction and component labels distinguish the
   different cusps.
4. Accidental curve intersections away from the marked finite-face points
   and marked cluster strata are removed by iterated blowups. The remaining
   marked incidence sets recover the ordinary face poset, while the cluster
   contractions identify exactly the ideal vertices and their incident
   faces.
5. The arithmetic decoration on the surface is the local system only on
   \[
      S_P^\circ=S_P\cap X,
   \]
   together with its \(\mathcal O_k\)-action, quadratic form, boundary
   monodromy/logarithmic extension data, and the marked root-line/reflection
   correspondence. No local-system extension across a nontrivial toroidal
   boundary is asserted.
6. The map to the original reflection stack is retained. Dihedral inertia
   is read there. It is not inertia of the neat surface: a neat cover has
   trivial finite inertia. Alternatively every finite label \(m_{ij}\) is
   computed directly from the two marked root reflections.
7. The oriented real chamber is retained as external ambient arithmetic
   data (equivalently, the positive component and oriented root
   halfspaces). It decides which complex mirror intersections are actual
   faces. It is not claimed to be an intrinsic real component of a generic
   surface section.

Thus, after making auxiliary choices, there is a faithful decorated
association

\[
 P\rightsquigarrow
 (S_P,\{C_i\},Z_P,\{{\cal K}_c\},
  \mathbb V|_{S_P^\circ},\text{root/descent marking},C)
 \longmapsto P.
\]

Faithfulness here is deliberately modest: the finite arithmetic root and
chamber marking is retained, rather than reconstructed from the ordinary
Neron--Severi lattice or the unmarked curve configuration. Since the
points, sections, field extension, and resolution are choices, the arrow is
not a functor without adding a category of choices.

## Retraction: why one boundary point is false

The earlier version asserted that, after resolving the mirror-boundary
arrangement, one could choose a boundary point incident to all cusp-wall
curves. This is false even for the arithmetic Coxeter triangle of type
\((2,3,\infty)\).

At its ideal vertex the two cusp-wall mirrors and the toroidal boundary give
three distinct local divisor branches. A log resolution makes the total
marked divisor simple normal crossings. On a smooth surface an SNC divisor
has at most two irreducible components through a point. Hence no point of
the resolved surface lies simultaneously on the boundary and on both strict
wall transforms. The resolution replaces the old triple incidence by an
exceptional cluster.

The same dimension count is universal: even a simple ideal vertex of an
\(n\)-polytope has \(n\) wall branches, and including the boundary gives
\(n+1\) branches in the \(n\)-dimensional Shimura variety. Log resolution
must separate that incidence. Reducible affine cusp links can have still
more wall branches. Therefore the cluster, its labels, and its contraction
map are essential; a single post-resolution cusp point cannot replace them.

## Proof of the repaired construction

### 1. The finite-volume boundary statement

For a finite-volume finite-sided hyperbolic polytope,

\[
 \overline P\cap\partial\mathbb H^n
 =\{\text{finitely many ideal vertices}\}.
\]

There is no positive-dimensional face contained in the sphere at infinity.
An edge or higher-dimensional face can, of course, have ideal vertices in
its closure while its relative interior lies in \(\mathbb H^n\). Those are
ordinary faces in the rest of the construction.

Every ideal vertex in the arithmetic quadratic realization gives a
rational isotropic ray \(\ell_c=\mathbb Q u_c\). Under the determinant-line
real-slice embedding, a ray tending to \([u_c]\) limits to the
**zero-dimensional Baily--Borel stratum indexed by the isotropic line**
\(\ell_c\), not to a positive-dimensional face at infinity. The ambient
type-IV variety can also have one-dimensional Baily--Borel strata indexed
by isotropic planes; those are different boundary strata and are not being
identified with a hyperbolic ideal vertex here.

All incident wall roots satisfy \((r_i,u_c)=0\), so the closures of their
special divisors contain the isotropic-line cusp before resolution. A
projective toroidal resolution turns this common limiting incidence into the
finite labeled cluster \({\cal K}_c\). Distinct marked cusp rays can be
separated on a sufficiently deep principal congruence cover: after choosing
primitive representatives, take a level at which no two are congruent up to
sign. If one does not make that cover, the ray labels must remain part of
the decoration even when several rays map to the same coarse cusp. The full
dual complex and the labeled contraction retain exactly the information
which strict-transform incidence loses.

### 2. Algebraic representatives of finite-face and cluster strata

For an ordinary face \(F\), let \(I(F)\) be its incident facet set. A point
in the relative interior of \(F\) on the real slice witnesses nonemptiness of
the corresponding complex Heegner intersection. On the sheet-separated
cover, remove the finitely many nonincident marked branches. The remaining
locally closed algebraic stratum is nonempty.

There is no hidden nonsimplicity at such a finite face. At
\(x\in\operatorname{relint}(F)\), a normal slice to \(F\) is a chamber for
the finite real reflection group generated by the incident facet
reflections. Its rank is \(\operatorname{codim}F\). Every chamber of a
finite real reflection group is simplicial: its simple roots form a basis,
also when the group is reducible. Consequently

\[
 |I(F)|=\operatorname{codim}F,
\]

and the incident facet roots are linearly independent. The corresponding
holomorphic Heegner equations have independent differentials at the marked
real-slice point, so the **selected facet branches** are simple normal
crossings there. The sheet-separated hypothesis is essential: an additional
Heegner branch through the same face could destroy SNC even though the
facet branches themselves are independent. With those unwanted branches
excluded, the log resolution can be chosen to be an isomorphism over this
open locus.

A nonempty finite-type stratum over a number field has a closed point, whose
residue field is finite over that field. The same statement applies to the
open strata of each finite exceptional cluster. Since there are finitely
many faces, cusps, components, and strata, one common finite extension \(K\)
makes a chosen finite set of geometric representatives rational. If one
wants to remain over the smaller field, take the corresponding finite
Galois-stable reduced zero-cycles instead. This is an algebraic-point
statement; it does not say that the chosen point itself lies on the original
real chamber.

### 3. Assigned-first-jet Bertini

Let \(Z\) be the resulting finite reduced set of geometric points, including
ordinary-face representatives and the representatives needed for the cusp
clusters, and write \(T=|Z|\). Embed the resolved compactification in
\(\mathbb P^N\). Forms of degree \(2T-1\) simultaneously separate all first
jets at \(Z\).

Indeed, for each ordered pair \(i\ne j\), choose a linear form
\(\ell_{ij}\) which vanishes at \(p_j\) and not at \(p_i\), and put

\[
 Q_i=\prod_{j\ne i}\ell_{ij}^{\,2}.
\]

Then \(Q_i\) has degree \(2(T-1)\), vanishes to second order at every
\(p_j\ne p_i\), and is a unit modulo \(\mathfrak m_{p_i}^2\). Multiplication
by \(Q_i\) therefore takes linear forms onto the full ambient first-jet space
at \(p_i\), without changing any other first jet. Summing \(Q_iL_i\) proves
the claim. Restriction surjects onto the tangent first jets of the embedded
smooth compactification. Surjectivity also descends for a Galois-stable
zero-cycle; alternatively it is immediate after the common splitting
extension \(K\).

At each selected point choose a two-dimensional tangent plane avoiding the
finitely many forbidden marked-divisor tangencies. Prescribe \(n-2\)
independent differentials with that common kernel. Bertini away from the
finite base set, together with the prescribed independent jets at the base
points, gives a smooth complete-intersection surface through every selected
point and not contained in a marked divisor. For \(n=2\), use the resolved
compactification itself.

Consequently, if all \(n-2\) hypersurfaces have degree at most \(2T-1\),

\[
 \deg S\leq \deg(\widetilde X)(2T-1)^{n-2}
\]

before surface cleanup blowups. Here \(T\) includes the chosen cluster
strata. The former formula in which \(T\) counted only one putative cusp
point is retracted. A general complementary linear cut of each chosen
stratum supplies an explicit Galois-stable zero-cycle bounded by the degree
of its projective closure, so one can replace \(T\) by the sum of those
stratum degrees after first choosing disjoint reduced cuts.

### 4. Cleanup and combinatorial recovery

Distinct sheet-separated marked divisors have no common component on a
general surface section. Their unmarked intersection is finite. Blow up
unmarked intersection points, and repeat at infinitely near unmarked
intersections until the strict transforms are separated there. Local
intersection multiplicity drops at every such blowup, so the sum of the
pairwise intersection numbers bounds the number of cleanup steps. No
marked finite-face or cluster point is blown up.

For ordinary faces, the map

\[
 F\longmapsto I(F)
\]

is injective and reverses inclusion: a face of a convex polytope is the
intersection with all facets containing it. Thus the marked incident sets
recover the ordinary face poset. For an ideal vertex, \(I(c)\) is recovered
from the wall attachments in \(\mathcal K_c\), and an ordinary face \(F\)
contains \(c\) in its closure exactly when \(I(F)\subseteq I(c)\). The
cluster contraction distinguishes cusps with different boundary points.
The marked quadratic module and oriented halfspaces then recover all finite,
parallel, and ultraparallel Gram data.

## What the surface does not intrinsically recover

Three distinctions are mandatory.

1. A neat cover erases reflection torsion. Dihedral inertia belongs to the
   retained original arithmetic stack, not to the neat surface.
2. A bare quadratic local system and a collection of labeled curves do not
   canonically identify facet root vectors or a minimal field of definition.
   The \(\mathcal O_k\)-action, descent datum, root lines, reflections, and
   curve-to-root correspondence are part of the decoration. The marked
   field \(k\) is recovered from that coefficient data; no claim is made that
   it is the intrinsic minimal field of the chosen surface over \(K\).
3. For \(n>4\), a generic complex surface section need not even meet the
   \(n\)-dimensional real hyperbolic slice by an expected-dimension count.
   The oriented chamber is therefore retained ambient data, not an intrinsic
   real component of \(S_P\).

In particular this theorem is a finite algebraic carrier/encoding theorem,
not a theorem that the ordinary divisor lattice of \(S_P\) determines \(P\).

## Borcherds/Belyi refinement

If the facet principal part passes the Bruinier cusp-form obstruction, the
Borcherds automorphic discriminant restricts to a meromorphic automorphic
section on \(S_P^\circ\), with its compactification orders recorded
separately. A ratio of two sections of the same automorphic line bundle and
character gives a rational function when the ratio extends meromorphically
to \(S_P\). This is a possible Belyi analogue, not automatically a
three-branch-value map. Compensator divisors, boundary orders, hidden
real-slice cuts, and extension across the cluster must all be audited.

## Remaining leverage gate

The repaired theorem supplies a universal **decorated** surface carrier at a
fixed separated Heegner support. It does not yet make these surfaces easier
to classify and does not give a polytope theorem beyond the input
arithmetic/finite-volume facts. Quantitative pushback now requires bounds
for the level, field extension, degrees and number of all finite-face and
cusp-cluster strata, curve degrees, cleanup intersections, and Borcherds
correction in terms of

\[
 n,\quad [k:\mathbb Q],\quad |\operatorname{disc}L|,\quad
 \operatorname{level}(L),\quad\#\operatorname{facets}(P).
\]

Only such bounds, followed by a Hilbert/Hurwitz/moduli classification that
returns a new example or obstruction, would provide the requested
classification leverage.
