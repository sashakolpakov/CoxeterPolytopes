# Absolute facet degree, Mostow rigidity, and an effective Coxeter-dessin census

**Status (2026-08-10): proved conditional leverage theorem for intrinsic
divisor-lattice carriers in dimensions at least three.** The coefficient
lattice must occur intrinsically in the Neron--Severi group, the facet sheets
must be distinct marked divisors, and their absolute polarized degree must be
bounded. Under those hypotheses the surface bounds the field degree, all
finite Coxeter labels, and the number of facets. Mostow--Prasad strong
rigidity then supplies the missing zero-dimensionality: a fixed Coxeter
matrix has at most one finite-covolume marked reflection realization in a
fixed dimension \(n\geq3\).

This does not assert that every arithmetic Coxeter polytope has such an
intrinsic carrier. The universal Bertini carrier retains the root lattice as
an external local system and does not satisfy the Neron--Severi hypothesis.
Dimension two is also excluded: fixed-angle Coxeter polygons can have
positive-dimensional length moduli.

The strong-rigidity input is the Mostow--Prasad theorem for finite-covolume
hyperbolic lattices. It is unrelated to Prasad's arithmetic covolume formula;
neither that formula nor ABSW/Borel--Prasad finiteness is used here.

## 1. Absolute FACET degree

Let \(X\) be a projective \(n\)-fold, let \(A\) be an ample Cartier divisor,
and suppose an exact-support neat carrier for
\(P\subset\mathbb H^n\) has pairwise distinct irreducible facet divisors
\(H_1,\ldots,H_f\). Put

\[
 Z_{\rm FAC}=\sum_{i=1}^f q_iH_i,\qquad
 q_i\in\mathbb Z_{>0},\qquad q_{\min}=\min_iq_i,
\tag{1.1}
\]

and define

\[
 \mathcal C_{\rm FAC}(X,A,Z_{\rm FAC})
 =\frac{Z_{\rm FAC}\cdot A^{n-1}}{q_{\min}}.
\tag{1.2}
\]

This is a Hilbert invariant rather than an inserted wall count:

\[
 \chi\!\left(Z_{\rm FAC},\mathcal O(mA)\right)
 =\frac{Z_{\rm FAC}\cdot A^{n-1}}{(n-1)!}m^{n-1}
  +O(m^{n-2}).
\tag{1.3}
\]

### Proposition 1.1

Every separated facet sheet consumes at least one positive integral unit:

\[
 \boxed{f\leq\mathcal C_{\rm FAC}.}
\tag{1.4}
\]

Indeed, ampleness and integrality give
\(H_i\cdot A^{n-1}\in\mathbb Z_{>0}\), so

\[
 Z_{\rm FAC}\cdot A^{n-1}
 =\sum_iq_i(H_i\cdot A^{n-1})\geq q_{\min}f.
\]

Separation is essential. On a coarse quotient, several real sheets can lie
in one Heegner component. At an exact-support level, equality of the
components attached to \(F_i\) and \(F_j\) would put both facets in the
closed real trace of one marked divisor, contradicting exact support.

For an intrinsic surface carrier \((S,A)\), this becomes

\[
 \mathcal C_{\rm FAC}(S,A,Z_{\rm FAC})
 =\frac{A\cdot\sum_iq_iC_i}{q_{\min}},
\tag{1.5}
\]

where the \(C_i\) are distinct marked curves. If the \(n\)-fold is cut by
\(n-2\) members of \(|mA|\), then

\[
 A|_S\cdot Z_{\rm FAC}|_S
 =m^{n-2}Z_{\rm FAC}\cdot A^{n-1}.
\tag{1.6}
\]

Higher-face multiple points need not be SNC. Their full incidence scheme, or
the labelled exceptional clusters of a log resolution, must be retained.
This does not change (1.6).

## 2. A Borcherds upper certificate

Let \(\bar X\) be a smooth projective compactification of a neat carrier.
Suppose a meromorphic Borcherds section \(s\) of \(\bar{\mathcal L}\) has

\[
 \operatorname{div}(s)
 =D^+_{\rm int}+B^+-D^-_{\rm int}-B^-,
\tag{2.1}
\]

where \(B^\pm\) are toroidal boundary divisors and
\(Z_{\rm FAC}\leq D^+_{\rm int}\) coefficientwise. Intersecting the divisor
identity with \(A^{n-1}\) gives

\[
 \boxed{
 q_{\min}f
 \leq Z_{\rm FAC}\cdot A^{n-1}
 \leq c_1(\bar{\mathcal L})\cdot A^{n-1}
      +D^-_{\rm int}\cdot A^{n-1}
      +B^-\cdot A^{n-1}.}
\tag{2.2}
\]

Boundary zeros only improve the inequality. Interior poles cannot be
omitted: positive and negative Heegner terms can cancel in the line class.

In the compact anisotropic lane there is no boundary. If the compensated
principal part is effective in the interior, take \(D^-_{\rm int}=0\).
After killing a finite multiplier, take \(A=\lambda\) and
\(\bar{\mathcal L}=\lambda^w\). Then

\[
 \boxed{
 f\leq
 \frac{Z_{\rm FAC}\cdot\lambda^{n-1}}{q_{\min}}
 \leq\frac{w\lambda^n}{q_{\min}}.}
\tag{2.3}
\]

Safe compensators do not occur in the definition (1.2). They only provide a
global effective section whose absolute line degree bounds the marked
facet degree. In the rational noncompact lane one instead uses a projective
toroidal compactification and an integral ample \(A\); every boundary pole
must appear in (2.2). Automorphic boundary divisors are not
positive-dimensional ideal faces of \(P\): a finite-volume polytope meets
infinity only in ideal vertices.

### Stack and cover normalization

On the original DM stack, intersection numbers can have stabilizer
denominators. Either pass to a neat cover or insert a stated common
Cartier/inertia denominator. An unbounded denominator yields no uniform
lower bound.

The degree must remain absolute on the sheet-separating cover. If
\(p:Y\to X\) is finite etale of degree \(d\), then

\[
 c_1(p^*\mathcal L)(p^*A)^{n-1}
 =d\,c_1(\mathcal L)A^{n-1}.
\tag{2.4}
\]

Dividing by \(d\) is fatal: a degree-\(d\) cover can split one divisor into
\(d\) sheets while \(d^{-1}\deg(p^*H)=\deg(H)\) stays fixed. Thus none of

\[
 \frac{w\lambda^n}{d q_{\min}},\qquad
 \frac{w}{q_{\min}},\qquad
 \text{the unseparated coarse Heegner degree}
\tag{2.5}
\]

bounds the separated sheet count without an independent bound on \(d\).
The valid quantity is \(w\lambda^n/q_{\min}\) on a fixed
canonical/minimal separating carrier, or the absolute FACET degree itself.

## 3. The intrinsic surface bounds all finite angles

Let \(k\) be the Vinberg ground field, \(d=[k:\mathbb Q]\), and suppose the
full \(\mathcal O_k\)-lattice of field rank \(n+1\) embeds faithfully in
\(\operatorname{NS}(S)\). If \(\rho(S)\leq R\), comparison of
\(\mathbb Z\)-ranks gives

\[
 \boxed{d(n+1)\leq\rho(S)\leq R.}
\tag{3.1}
\]

This is false for a carrier on which the root lattice is merely an external
local system.

For a finite Coxeter label \(m\geq3\), a length-two cyclic Gram product is

\[
 4\cos^2(\pi/m)=2+2\cos(2\pi/m)\in k.
\tag{3.2}
\]

It generates the maximal real cyclotomic field of conductor \(m\), of degree
\(\varphi(m)/2\). Consequently

\[
 \varphi(m)\leq2d\leq\frac{2R}{n+1}.
\tag{3.3}
\]

Moreover,

\[
 \frac{\varphi(m)^2}{m}
 =\prod_{p^a\parallel m}p^{a-2}(p-1)^2\geq\frac12,
\tag{3.4}
\]

so \(m\leq8d^2\). The exact finite alphabet is therefore

\[
 \mathcal A_d=\{2\}\cup
 \{m\geq3:\varphi(m)\leq2d\},\qquad m\leq8d^2.
\tag{3.5}
\]

The symbol \(\infty\) is added for disjoint facets. It does not prescribe the
dotted distance: \(g_{ij}=-1\) is parallel and \(g_{ij}<-1\) is
ultraparallel.

## 4. Strong rigidity fixes every dotted entry

### Theorem 4.1

Fix \(n\geq3\) and an ordered Coxeter matrix \(M=(m_{ij})\). There is at most
one full-dimensional finite-volume hyperbolic reflection realization of
\(M\), up to an isometry carrying every ordered facet to its counterpart.

#### Proof

Let \(W_M\) be the abstract Coxeter group and let

\[
 \rho_a:W_M\hookrightarrow\operatorname{Isom}(\mathbb H^n),
 \qquad a=1,2,
\]

be two realizations. The Coxeter fundamental-domain theorem makes both
representations faithful, including for nonsimple and cusped chambers. The
even subgroup \(W_M^+\), the kernel of Coxeter word parity, is canonical;
both \(\rho_a(W_M^+)\) are finite-covolume orientation-preserving lattices.

Mostow--Prasad strong rigidity supplies an isometry \(g\) satisfying

\[
 g\rho_1(h)g^{-1}=\rho_2(h)\qquad(h\in W_M^+).
\tag{4.1}
\]

For a standard reflection \(s_i\), put

\[
 c_i=\rho_2(s_i)^{-1}g\rho_1(s_i)g^{-1}.
\]

Conjugation by either reflection realizes the same automorphism
\(h\mapsto s_ihs_i\) of \(W_M^+\), so \(c_i\) centralizes
\(\rho_2(W_M^+)\). A finite-covolume hyperbolic lattice is Zariski dense,
and its centralizer in the projective isometry group is trivial. Thus
\(c_i=1\) for every \(i\). The isometry \(g\) conjugates all ordered
reflections, so it identifies the reflecting hyperplanes, oriented chamber,
and outward unit-normal geometric Gram matrix. \(\square\)

This is global rigidity of complete finite-volume lattices, not local
rigidity. Nonuniform three-dimensional representation varieties can have
incomplete deformation directions; the full cusp-completeness gate leaves
the unique Mostow point. The theorem fixes parallel/ultraparallel data even
though the abstract label \(\infty\) carries no distance.

The statement fails in dimension two: polygons with at least four sides can
have fixed angles and positive-dimensional length moduli.

## 5. Effective intrinsic-surface finiteness

### Theorem 5.1

Fix \(R,B\in\mathbb Z_{>0}\). Consider marked compact or finite-volume
arithmetic Coxeter polytopes \(P\subset\mathbb H^n\), \(n\geq3\), admitting
an exact intrinsic surface carrier such that:

1. the full \(\mathcal O_{k(P)}\)-lattice of rank \(n+1\) embeds faithfully
   in \(\operatorname{NS}(S)\), with \(\rho(S)\leq R\);
2. every facet is a distinct marked irreducible curve/divisor on a
   sheet-separated carrier; and
3. \(\mathcal C_{\rm FAC}\leq B\).

Then the class contains finitely many marked polytopes. Their exact outward
unit-normal geometric Gram matrices, compact/cuspidal alternatives, parallel
versus ultraparallel tags, and Vinberg ground fields are effectively
enumerable. This does not assert finiteness of arbitrary primitive-root
coordinates or asymmetric Cartan entries: diagonal root rescaling can vary
those without changing the reflection arrangement or cyclic products.

#### Proof

Equation (3.1) bounds \(n\leq R-1\) and
\(d\leq\lfloor R/(n+1)\rfloor\). Proposition 1.1 gives \(f\leq B\), and
(3.5) gives a finite computable angle alphabet. Hence there are finitely
many ordered Coxeter matrices.

For every matrix and dimension, introduce variables for its \(\infty\)-Gram
entries and impose:

- diagonal entries \(1\), the fixed values
  \(-\cos(\pi/m_{ij})\), and \(g_{ij}\leq-1\) for \(\infty\)-pairs;
- rank \(n+1\), Lorentzian inertia, and spanning of the normal space;
- a nonempty irredundant oriented chamber; and
- the literal compact or finite-volume condition, including every cusp's
  full parabolic neighbourhood.

These are finitely many algebraic equalities and real semialgebraic
inequalities. One may equivalently impose in a Klein chart that the
Euclidean facet-halfspace intersection lies in the closed unit ball. Gauge
fix the Lorentz action by enumerating a spanning \(n+1\)-subset of normals
and putting it in a triangular sign-normalized frame. Real quantifier
elimination gives a terminating exact decision and algebraic sample points.
Theorem 4.1 says that the accepted Gram set for a fixed matrix and dimension
has at most one point, up to the finite gauge choices.

Finally compute

\[
 k(P)=\mathbb Q\bigl(\operatorname{Cyc}(G(P))\bigr)
\tag{5.1}
\]

from cyclic Gram products, isolate its real embeddings, and retain the one
passing Vinberg's distinguished-signature and arithmetic-integrality tests.
This enumerates the fields from the finite Coxeter-dessin census itself.
\(\square\)

Different matrices can define isomorphic abstract groups or differ only by
a wall permutation. Exact Gram-isometry and permutation tests remove these
finite duplicates.

## 6. Sharp K3 endpoint

For a complex K3 surface, \(R=20\). An intrinsic full-trace carrier for a
compact arithmetic \(P\subset\mathbb H^9\) satisfies \(10d\leq20\).
A rational compact quadratic chamber in this rank is excluded by Meyer's
isotropy theorem. Hence

\[
 \boxed{d=2,\qquad\rho(S)=20.}
\tag{6.1}
\]

For each fixed primitive Picard-rank-twenty Neron--Severi lattice, the carrier
is a Picard-rank-twenty ("singular" in classical K3 terminology)
zero-dimensional marked period. The union over all such lattices and
discriminants is not finite merely from \(\rho=20\). Its finite angle
alphabet is exactly

\[
 \boxed{\{2,3,4,5,6,8,10,12\}.}
\tag{6.2}
\]

For \(n\geq10\), compactness forces \(d\geq2\), while
\(d(n+1)\geq22>20\). Thus

\[
 \boxed{\text{no compact }H^{10+}\text{ polytope has an intrinsic
 full-trace K3 divisor-lattice carrier}.}
\tag{6.3}
\]

This is a theorem about this carrier lane, not nonexistence of compact
Coxeter polytopes. It focuses the K3 search at \(H^9\): quadratic field,
Picard rank twenty, the eight finite labels in (6.2), and bounded absolute
facet degree.

## 7. Gram finiteness versus carrier finiteness

Theorem 5.1 proves that the returned unit-normal geometric Coxeter Gram image
is finite. It does
not make every surface or automorphic carrier above a Gram point rigid. A
fixed Picard lattice can support a positive-dimensional K3 period family,
and a bounded Hilbert polynomial gives only a finite-type Hilbert scheme,
which can contain algebraic points of unbounded degrees.

Fields of definition for the full carriers require the complete marked
Hilbert/incidence/Borcherds stack to be a zero-dimensional fine
Hurwitz-style locus. Full level, an ordered root frame, and either a
canonical projective gauge forced by the marked data or a projective-group
quotient must kill all automorphisms; a freely chosen frame would add a
positive-dimensional torsor. Its geometric points then form a finite
Galois set, and trivial automorphisms make the field of moduli a field of
definition. Without rigidification, a residual gerbe can obstruct descent.
The Coxeter ground fields need no extra hypothesis: (5.1) computes them from
the finite unit-normal geometric Gram image.

## 8. Completion path

1. Prove that a declared P-only class has intrinsic, rather than external,
   \(\mathcal O_k\)-lattice surface carriers and exact same-P return.
2. Derive \(R\) from its bounded surface/Hilbert class and bound the absolute
   FACET degree \(B\), directly or by the pole-audited certificate (2.2).
3. Enumerate (3.5), all matrices with \(f\leq B\), and their gauge-fixed
   semialgebraic Gram loci.
4. Apply the full finite-volume gate before Theorem 4.1; do not mistake an
   incomplete cusp deformation for another lattice.
5. Compute cyclic-product fields and replay Vinberg on every survivor.
6. If fields of definition of the surfaces are required, prove the separate
   zero-dimensional fine-carrier gate in Section 7.

Once intrinsic-carrier existence and a class-uniform absolute FACET degree
bound are proved, no separate dotted-height theorem or imported field
catalogue is needed in dimensions \(n\geq3\).

## 9. Exact arithmetic replay

The totient cutoff, the rational/quadratic angle alphabets, and the sharp K3
\(H^9/H^{10}\) rank consequences are replayed by
[the certificate script](../scripts/certify_intrinsic_surface_coxeter_dessin_bounds.py)
and
[its tests](../tests/test_intrinsic_surface_coxeter_dessin_bounds.py).
The Mostow--Prasad step is a proof theorem, not a numerical inference.

## 10. Primary inputs

- G. D. Mostow, *Strong Rigidity of Locally Symmetric Spaces*, Annals of
  Mathematics Studies **78**, Princeton University Press, 1973.
- G. Prasad, *Strong rigidity of Q-rank 1 lattices*, Invent. Math. **21**
  (1973), 255--286. This is the rigidity theorem, not an arithmetic volume
  formula.
- E. B. Vinberg, *Hyperbolic reflection groups*, Russian Math. Surveys
  **40** (1985), 31--75, <https://www.mathnet.ru/rm2140>.
- J. H. Bruinier,
  [*Regularized theta lifts for orthogonal groups over totally real fields*](https://arxiv.org/abs/0908.3076),
  J. Reine Angew. Math. **672** (2012), 177--222.
- [Facet-volume and Coxeter-degree bridge](facet_volume_coxeter_degree_bridge.md).
