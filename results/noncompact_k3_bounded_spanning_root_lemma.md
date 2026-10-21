# Finite volume forces a bounded spanning packet of physical K3 roots

**Status (2026-08-11): proved, classification-free, with the uniform
constant \(B=13\).** Let \(X\) be a complex projective K3 surface of Picard
rank

\[
                         3\le \rho\le20,
\]

and suppose its exact full-Picard root chamber is a finite-sided
finite-volume Coxeter polyhedron in \(\mathbf H^{\rho-1}\). Then there are
actual facet curves \(E_1,\ldots,E_\rho\) such that

\[
 \operatorname {span}_{\mathbf Q}\{E_1,\ldots,E_\rho\}
      =\operatorname {NS}(X)_{\mathbf Q},
 \qquad
                   0\le E_i.E_j\le13\quad(i\ne j).       \tag{0.1}
\]

The proof has two elementary pieces. A facet of minimal visual aperture has
a **closed star** with all normalized pairings strictly below \(14\). In the
noncompact case, that star still spans: take two distinct extreme points of
the chosen facet. The complete incident packet at a finite vertex spans its
orthogonal hyperplane, and boundedness of the horospherical link gives the
same statement at a nonsimple ideal vertex. The two vertices determine two
different hyperplanes, whose sum is the whole Lorentz space.

Thus the cusp-local divisibility escape is closed by a second vertex of one
global narrow facet. No affine-diagram classification, reflective-lattice
classification, facet-count bound, cusp simplicity, compactness, arithmetic
volume formula, or reflection maximality is used.

As an immediate consequence, only finitely many Picard hosts occur in this
exact finite-volume full-nef lane: the selected root lattice \(R\) has one of
finitely many nonsingular Gram matrices and

\[
                R\subseteq\operatorname {NS}(X)\subseteq R^*.
                                                               \tag{0.2}
\]

This proves the global repair lemma requested by the
[cusp-divisibility gate](noncompact_k3_cusp_divisibility_host_gate.md). It
does not prove that every raw maximal arithmetic chamber has an exact
full-nef K3 carrier; that remains the separate Property 4 problem.

## 1. Exact returned package and conventions

Put \(L=\operatorname {NS}(X)\), so \(L_{\mathbf R}\) has signature
\((1,\rho-1)\). Choose a positive-cone component. The hypothesis is that the
projectivized nef/root chamber

\[
                         P\subset\mathbf H^{\rho-1}
\]

is full-dimensional, finite-sided, and finite-volume, and that its listed
facets are exactly the walls of the full Picard root chamber. Each facet is
therefore represented by the class of an irreducible smooth rational curve

\[
                            E^2=-2.                       \tag{1.1}
\]

Distinct such curves have nonnegative integral intersection. Equality of
the returned chamber with the full nef chamber is important here: a divisor
shadow or a merely projected wall does not put the literal square-minus-two
normal in \(L\).

For two facet supports \(F,G\), with nonzero normals \(f,g\), write

\[
 b(F,G)=\frac{2\lvert(f,g)\rvert}
                   {\sqrt{\lvert f^2g^2\rvert}}.          \tag{1.2}
\]

For the K3 normalization (1.1), this is simply

\[
                            b(F,G)=F.G                    \tag{1.3}
\]

for distinct physical facets. Intersecting Coxeter walls have
\(b=2\cos(\pi/q)\le2\), parallel walls have \(b=2\), and ultraparallel walls
have \(b=2\cosh d>2\).

All face incidences below are taken in the projective closure of the Klein
model. Thus two facets which meet only at a cusp do meet for purposes of the
closed star.

## 2. The minimal-aperture closed star has pairing below 14

Fix \(O\in\operatorname {int}P\). The visual cone from \(O\) to a facet
support is a round cone of angular opening \(0<\theta_F<\pi\). There are
finitely many facets, so choose \(E\) with minimal opening

\[
                         \theta_E\le\theta_F
             \quad\hbox{for every facet }F.              \tag{2.1}
\]

Let \(F,G\) be facets whose closed faces both meet \(E\). First suppose
their supporting hyperplanes intersect in hyperbolic space or are
asymptotic. The ordinary hyperplane formula gives, respectively,

\[
             b(F,G)=2\lvert\cos\alpha\rvert<2
             \quad\hbox{or}\quad b(F,G)=2,               \tag{2.2}
\]

so these pairs already satisfy the desired estimate.

It remains to treat ultraparallel supports. Their visual cones are disjoint.
Let \(\theta_{FG}>0\) be the angular gap between those cones. Choose a ray
through a point of \(\overline F\cap\overline E\) and one through a point of
\(\overline G\cap\overline E\). Both rays lie in the closed visual cone of
\(E\), so the gap between the two cones is no greater than the angle between
these two rays. Consequently

\[
       0<\theta_{FG}\le\theta_E
       \le\min\{\theta_F,\theta_G\}<\pi.                 \tag{2.3}
\]

For this ultraparallel pair, the two-dimensional plane through the cone
axes gives the exact identity

\[
 b(F,G)=
 4\frac{
 \sin((\theta_F+\theta_{FG})/2)
 \sin((\theta_G+\theta_{FG})/2)}
 {\sin(\theta_F/2)\sin(\theta_G/2)}-2.                  \tag{2.4}
\]

This ordinary visual-cone calculation depends only on the three supporting
hyperplanes. It is not being asserted for overlapping visual cones; those
were disposed of by (2.2). For
\(0\le c\le a<\pi\),

\[
              \frac{\sin((a+c)/2)}{\sin(a/2)}<2.         \tag{2.5}
\]

Indeed, if \(a\le\pi/3\), the quotient is at most
\(\sin(a)/\sin(a/2)=2\cos(a/2)<2\); if \(a>\pi/3\), it is at most
\(1/\sin(a/2)<2\). Applying (2.5) to both factors in (2.4) yields

\[
                              b(F,G)<14.                  \tag{2.6}
\]

Pairs involving \(E\) itself satisfy the sharper Coxeter bound \(b\le2\).
Consequently every off-diagonal normalized pairing in the closed star

\[
 \mathcal S(E)=\{E\}\cup
 \{F:\overline F\cap\overline E\ne\varnothing\}          \tag{2.7}
\]

is strictly below \(14\). This is the same elementary narrow-facet lemma
used in the independently audited
[compact field theorem](k3_compact_narrow_facet_field_finiteness.md); no
classification is hidden in it. In particular, compactness is not used
anywhere in this section.

## 3. Every ideal boundary point has a complete incident span

The only genuinely noncompact input is the following standard cusp fact.

### Lemma 3.1 (finite-volume cusp span)

Let \([u]\in\overline P\cap\partial_\infty\mathbf H^{\rho-1}\), and let
\(I(u)\) be the set of all facets whose closures contain \([u]\). Then the
normals of those facets span \(u^\perp\), or equivalently

\[
                   \bigcap_{F\in I(u)}H_F=\mathbf Ru.     \tag{3.1}
\]

This includes nonsimple cusps.

### Proof

Move \([u]\) to infinity in upper half-space. The incident supports become
vertical Euclidean hyperplanes. Every nonincident support is a bounded
hemisphere and ceases to impose a condition above a sufficiently high
horosphere. Hence the chamber above that horosphere is

\[
                          C\times(T,\infty),              \tag{3.2}
\]

where \(C\subset\mathbf R^{\rho-2}\) is the full-dimensional Euclidean
polyhedron cut out by the incident walls. Its hyperbolic volume is

\[
 \operatorname {vol}_{\mathrm E}(C)
       \int_T^\infty y^{-(\rho-1)}\,dy.                   \tag{3.3}
\]

For \(\rho\ge3\) the integral is finite and positive. An unbounded
full-dimensional Euclidean polyhedron has infinite Euclidean volume, so
finite volume of \(P\) forces \(C\) to be bounded.

Choose a second null vector \(v\) with \((u,v)=1\), and put
\(A=(\mathbf Ru+\mathbf Rv)^\perp\). Every incident normal can be written

\[
                              r_i=a_i u+w_i,
                  \qquad w_i\in A.                       \tag{3.4}
\]

The \(w_i\) are the oriented Euclidean side normals of \(C\). Boundedness of
\(C\) says that they span \(A\) and admit a positive dependence whose
\(w_i\)-part vanishes. Evaluating the corresponding strict inequalities at
an interior point shows that the \(u\)-coefficient in that dependence is
nonzero. Thus the \(r_i\) span \(A+\mathbf Ru=u^\perp\), proving (3.1).
\(\square\)

This is exactly the bounded-horosphere argument used in the independently
audited
[noncompact ground-field collapse](exact_trace_noncompact_ground_field_collapse.md).
It applies to every ideal point of the finite-volume chamber, not only to a
chosen simple cusp.

## 4. Two vertices make the narrow star span

### Lemma 4.1 (two-vertex star completion)

For \(\dim P\ge2\), the closed star \(\mathcal S(E)\) spans the entire
Lorentz space \(L_{\mathbf R}\).

### Proof

In the Klein model, the closure

\[
                 Q=\overline E\subset\overline{\mathbb B}^{\rho-1}
\]

is a compact convex set of affine dimension \(\rho-2\ge1\). Every compact
convex set of positive dimension has at least two distinct extreme points;
choose \(x\ne y\) in \(Q\).

An extreme point inside the Klein ball is an actual finite vertex of the
finite-sided facet, hence of \(P\). The normals of all facets through it
span its orthogonal hyperplane. Otherwise their common orthogonal subspace
would contain a direction besides the vertex vector; a sufficiently small
displacement in that direction remains timelike, preserves all active
equalities, and preserves the strict inactive inequalities. This would give
a segment of \(Q\) through the alleged extreme point.

An extreme point on the ideal sphere belongs to the ideal boundary of \(P\).
Lemma 3.1 makes its horospherical section bounded, so it is an actual
finite-volume ideal vertex, possibly nonsimple, and its complete incident
packet spans its orthogonal hyperplane. Thus every extreme point of \(Q\) is
an actual finite or ideal vertex with the required incident-span conclusion.

Apply this to \(x\) and \(y\). Every facet incident to either point belongs
to \(\mathcal S(E)\), because it meets \(E\) at that projective point.
Therefore

\[
 \operatorname {span}_{\mathbf R}\mathcal S(E)
                    \supseteq x^\perp+y^\perp.            \tag{4.1}
\]

The points are distinct projectively, so \(\mathbf Rx\cap\mathbf Ry=0\).
Nondegeneracy of the Lorentz form gives

\[
 (x^\perp+y^\perp)^\perp
       =(x^\perp)^\perp\cap(y^\perp)^\perp
       =\mathbf Rx\cap\mathbf Ry=0.                      \tag{4.2}
\]

Hence

\[
                         x^\perp+y^\perp=L_{\mathbf R},   \tag{4.3}
\]

proving the lemma. \(\square\)

This is the missing compact-to-cusped bridge. A single cusp packet is one
dimension short, but the two endpoint packets of one narrow facet cannot
share that missing direction. Because every root in both packets lies in
the same closed star, Section 2 bounds all cross-pairings between the two
packets as well as all pairings inside them.

## 5. Extraction of the physical K3 basis

By Lemma 4.1, the physical facet roots in \(\mathcal S(E)\) span
\(L_{\mathbf R}\). They are integral vectors, so the same matrix-rank
calculation over \(\mathbf Q\) shows that they span \(L_{\mathbf Q}\).
Extract \(\rho\) linearly independent members

\[
                         E_1,\ldots,E_\rho.              \tag{5.1}
\]

They are actual irreducible smooth rational curves and hence satisfy
\(E_i^2=-2\) and \(E_i.E_j\in\mathbf Z_{\ge0}\) for \(i\ne j\). Equations
(1.3) and (2.6) give

\[
                         0\le E_i.E_j<14.                 \tag{5.2}
\]

The fixed physical square \(E_i^2=-2\) is what makes the normalized quantity
equal to the literal intersection. Integrality sharpens (5.2) to

\[
                         0\le E_i.E_j\le13.               \tag{5.3}
\]

This proves (0.1), with the rank-independent choice

\[
                              \boxed{B(\rho)=13}.          \tag{5.4}
\]

## 6. Finite Picard hosts and complete root chambers

Let \(R=\sum_i\mathbf ZE_i\). Full rational spanning makes its Gram matrix

\[
                         G=(E_i.E_j)                       \tag{6.1}
\]

nonsingular. Its diagonal is \(-2\), and each off-diagonal entry belongs to
\(\{0,1,\ldots,13\}\). For fixed \(\rho\) there are at most

\[
                         14^{\binom\rho2}                  \tag{6.2}
\]

such symmetric matrices before nonsingularity and Lorentz-signature
filters. Since \(R\subseteq L\) and \(L\) is integral,

\[
                           R\subseteq L\subseteq R^*.      \tag{6.3}
\]

The discriminant group \(R^*/R\) is finite. Every possible \(L/R\) is one
of its finitely many subgroups; evenness, signature, and primitive
K3-lattice embeddability only remove candidates. Taking the finite union
over \(3\le\rho\le20\) proves finite Picard-host finiteness.

For a fixed \(L\) and positive-cone component, the group generated by all
square-minus-two root reflections acts transitively on the chambers of the
complete root-mirror arrangement: crossing a wall applies its reflection,
and a finite generic path between chamber interiors crosses finitely many
locally finite mirrors. Thus a **literal complete K3 root chamber** has one
Weyl-chamber orbit per oriented lattice. No reflection-maximality hypothesis
is needed for this downstream conclusion. Same-host maximality is relevant
only in a weaker native-host formulation that has not already returned the
complete K3 root arrangement.

## 7. Hostile scope and failure modes

The proof deliberately retains the following boundaries.

1. **Finite-sided exact return.** Minimal visual aperture is chosen from the
   finite returned facet set. No claim is made here for an arbitrary
   infinitely-faced convex set without an attained minimum.
2. **Projective closed incidence.** Facets meeting \(E\) only at an ideal
   point belong to the star. Replacing (2.7) by affine-interior ridge
   adjacency discards exactly the cusp packets used in Lemma 4.1.
3. **Full incident cusp packet.** A displayed affine subdiagram need not
   span \(u^\perp\). Lemma 3.1 uses every wall through the cusp and the
   bounded horospherical section.
4. **Finite volume.** Without it, \(C\) in (3.2) can be unbounded and the
   incident span can be smaller than \(u^\perp\). The family
   \(U(d)\oplus\langle-2\rangle\) in the cusp-divisibility gate remains the
   sharp local warning.
5. **Two distinct extreme points.** One ideal packet spans only a
   codimension-one degenerate hyperplane. The second point is the precise
   transverse repair; it is not optional.
6. **Literal physical roots.** The conversion of normalized pairing to an
   integral bound uses the same square-\(-2\) vectors as the facet normals.
   Merely projected physical roots do not suffice; the compact
   [projected-root scaling family](k3_compact_projected_root_scaling_no_go.md)
   shows this sharply.
7. **Rank at least three.** The chosen facet must have positive dimension.
   This matches the target \(3\le\rho\le20\); no rank-two statement is being
   smuggled in.
8. **No Property 4 conclusion.** The theorem starts after an exhaustive
   exact full-Picard/full-nef carrier has been returned. It does not prove
   that every raw maximal arithmetic chamber receives such a carrier or
   canonically minimalizes under the MMP.

Nonsimple ideal links cause no loss: Lemma 3.1 uses their entire bounded
Euclidean link, never a simplex assumption. Compactness is also unnecessary;
if both chosen extreme points happen to be finite, the same proof simply
specializes to the compact star-spanning argument.

## 8. Reproducibility

The exact replay is
[certifier](../scripts/certify_noncompact_k3_bounded_spanning_root_lemma.py),
with
[archive](../data/enumerations/noncompact_k3_bounded_spanning_root_lemma.json)
and
[regressions](../tests/test_noncompact_k3_bounded_spanning_root_lemma.py).

The replay checks all finite/finite, finite/ideal, ideal/finite, and
ideal/ideal two-vertex linear-algebra cases in every rank \(3\) through
\(20\). It also freezes the hostile countercheck that duplicating one
vertex packet leaves rank \(\rho-1\), the integral \(13\) bound, the finite
Gram-box count, and the exact scope flags. The analytic visual-cone and
bounded-horosphere arguments are the proof above; the executable is a
regression certificate for their algebraic interfaces, not a numerical
substitute for them.
