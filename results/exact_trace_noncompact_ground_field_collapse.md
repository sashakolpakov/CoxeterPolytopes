# Every exact admissible cusped trace chamber has ground field Q

**Status (2026-08-11): proved in every hyperbolic dimension.**  Let a
finite projective module over an order in a totally real field carry an
admissible form: Lorentzian at one distinguished embedding and positive
definite at every other real embedding.  If integral field-linear
reflections return an exact finite-volume noncompact hyperbolic chamber,
then

\[
                              \boxed{k=\mathbf Q}.          \tag{0.1}
\]

The proof is one fixed-line argument.  A finite-volume cusp has a bounded
horospherical section, so the fixed hyperplanes of all incident reflections
meet in exactly its null line.  Because the reflections are \(k\)-linear,
that line is defined over \(k\).  Its self-pairing is zero in \(k\), which
is impossible after any nonidentity embedding makes the form positive
definite.

No primitive wall normals, free-module basis, simplicial cusp, reflection
maximality, saturation, K3 surface, Coxeter classification, arithmetic
volume formula, or imposed height is used.  For exact full-trace K3
carriers this collapses the finite-volume noncompact field set to the
singleton \(\{\mathbf Q\}\).  The K3 Picard ceiling then gives the separate
carrier-relative dimension bound \(n\le19\).

## 1. Exact returned package

Let \(k\) be a totally real number field, let \(\mathcal O\subset k\) be an
order, and let \(M\) be a finite projective \(\mathcal O\)-module such that

\[
                              V=M\otimes_{\mathcal O}k
\]

has dimension \(n+1\), where \(n\ge2\).  Let \(h\) be a nondegenerate symmetric \(k\)-form
on \(V\).  Assume that at one real embedding \(\sigma _0\) it has signature
\((n,1)\), while \(h^\sigma\) is positive definite for every
\(\sigma\ne\sigma _0\).

The case \(n=1\) is vacuous: every noncompact hyperbolic interval has
infinite length.  Thus \(n\ge2\) loses no finite-volume chamber.

Suppose a finite set of \(\mathcal O\)-linear \(h\)-reflections

\[
                 s_i\in\operatorname {Aut}_{\mathcal O}(M)
\]

returns exactly the facets of a finite-sided finite-volume hyperbolic
chamber \(P\subset\mathbf H^n\) at \(\sigma _0\).  Thus the \(+1\)
eigenspace of \(s_i\) at the distinguished place is the support of the
corresponding facet.  Assume \(P\) is noncompact.  No minimal-field or
cyclic-generation hypothesis is needed.

## 2. A noncompact finite-volume chamber supplies a full cusp line

### Lemma 2.1 (bounded horosphere forces a one-dimensional intersection)

Let \(P\subset\mathbf H^n\) be a finite-sided, full-dimensional,
finite-volume noncompact convex polyhedron.  Then there is an ideal point
\([u_0]\in\partial_\infty P\) such that, if \(I\) is the set of all facets
whose closures contain \([u_0]\),

\[
                   \bigcap_{i\in I}H_i=\mathbf R u_0       \tag{2.1}
\]

in the ambient Lorentz space.  This includes nonsimple cusps.

### Proof

Noncompactness gives an interior sequence with an accumulation point
\([u_0]\) in the ideal boundary.
Move it to infinity in the upper-half-space model.  Facet supports through
\([u_0]\) become vertical Euclidean hyperplanes.  Every other facet support
is a bounded hemisphere, and its chosen half-space ceases to impose a
condition above some height because \([u_0]\) lies in the closure of
\(P\).  Hence, above one horosphere, the chamber is

\[
                              C\times(T,\infty),            \tag{2.2}
\]

where \(C\subset\mathbf R^{n-1}\) is the full-dimensional Euclidean
polyhedron cut out by the incident vertical walls.  The hyperbolic volume
there is

\[
             \operatorname {vol}_{\mathrm E}(C)
             \int_T^\infty y^{-n}\,dy.                    \tag{2.3}
\]

An unbounded full-dimensional Euclidean polyhedron has infinite Euclidean
volume, so finite hyperbolic volume forces \(C\) to be bounded.

Here is the exact linear consequence.  Choose a second null vector \(v_0\)
with \(h(u_0,v_0)=-1\), and put

\[
 E=(\mathbf R u_0+\mathbf R v_0)^\perp.
\]

An incident wall normal has the form

\[
                              e_i=a_i u_0+w_i,
                              \qquad w_i\in E.             \tag{2.4}
\]

Normalize the horosphere so the corresponding oriented wall inequality of
\(C\) is \(\langle w_i,x\rangle\le a_i\).  Boundedness says that its
recession cone

\[
                  \{z:\langle w_i,z\rangle\le0\text{ for all }i\}
                  =\{0\}.
\]

Consequently the \(w_i\) span \(E\) and zero lies in the interior of their
convex hull.  Thus there are
\(\lambda_i\ge0\), not all zero, with

\[
                         \sum_i\lambda_iw_i=0.             \tag{2.5}
\]

Orient the inequalities and evaluate them at an interior point of \(C\).
Their positive combination is strict, so
\(\sum_i\lambda_i a_i\ne0\).  Equations (2.4)--(2.5) give

\[
              0\ne\sum_i\lambda_i e_i
                    =\left(\sum_i\lambda_i a_i\right)u_0. \tag{2.6}
\]

Therefore the incident normals span
\(E+\mathbf R u_0=u_0^\perp\), and their common orthogonal is exactly
\((u_0^\perp)^\perp=\mathbf R u_0\).  This proves (2.1).  Extra incident
facets do not change the conclusion. \(\square\)

The proof also shows why merely displaying an affine subdiagram is not
enough: the walls must be the full incident cusp packet whose horospherical
section is bounded.

## 3. The fixed line is defined over the returned field

Let \(I\) be the incident set from Lemma 2.1 and form the \(k\)-subspace

\[
                         W=\bigcap_{i\in I}\ker(s_i-1)
                         \subset V.                       \tag{3.1}
\]

This definition uses only the returned reflection operators.  In
particular, it does not choose wall normals in \(M\), primitive or
otherwise.

Package all maps \(s_i-1\) into one \(k\)-linear map

\[
                 A:V\longrightarrow\bigoplus_{i\in I}V.
\]

Extension of scalars by \(\sigma _0:k\hookrightarrow\mathbf R\) is exact,
so

\[
 W\otimes_{k,\sigma _0}\mathbf R
 =\bigcap_{i\in I}\ker(s_i^{\sigma _0}-1)
 =\bigcap_{i\in I}H_i
 =\mathbf R u_0.                                         \tag{3.2}
\]

Consequently \(W\) is a one-dimensional \(k\)-space.  Choose
\(0\ne u\in W\).  Its distinguished image spans the cusp line, hence

\[
                         \sigma _0\bigl(h(u,u)\bigr)=0.
\]

The embedding \(\sigma _0\) is injective, so already

\[
                              h(u,u)=0\quad\hbox{in }k.    \tag{3.3}
\]

## 4. Definite conjugates force degree one

Suppose \(\sigma:k\hookrightarrow\mathbf R\) is different from
\(\sigma _0\).  A nonzero vector over a field remains nonzero after scalar
extension, so \(u^\sigma\ne0\).  Applying \(\sigma\) to (3.3) gives

\[
                              h^\sigma(u^\sigma,u^\sigma)=0.
\]

This contradicts positive definiteness of \(h^\sigma\).  Therefore no
nonidentity real embedding exists.  Since \(k\) is totally real,
\([k:\mathbf Q]=1\), proving (0.1).

### Theorem 4.1 (noncompact ground-field collapse)

Every package satisfying Section 1 and returning an exact finite-volume
noncompact chamber has coefficient field \(k=\mathbf Q\).  The statement
holds for simple or nonsimple cusps and in every hyperbolic dimension.

### Corollary 4.2 (full-trace K3 field set)

Every exact primitive full-trace complex-K3 eigen-nef carrier whose returned
chamber is finite-volume and noncompact has the singleton field set
\(\{\mathbf Q\}\).  After this collapse, the trace rank is

\[
                         n+1=\rho(X)\le20,
                         \qquad n\le19.                   \tag{4.1}
\]

The K3 Picard-rank estimate is not needed for the field conclusion; it is
used only for (4.1).  This is a bound for exact full-trace K3 carriers.  It
does not supersede the distinct rational weak-log branch with bound
\(n\le25\), and it does not assert that every raw chamber has a K3 carrier.

## 5. Scope

The theorem requires the full exact cusp packet to act \(k\)-linearly on one
admissible returned form.  It does not apply to a partial wall subsystem, an
affine diagram not known to be a genuine finite-volume cusp, a form lacking
definite nonidentity conjugates, or a chamber not returned over the declared
field.  It says nothing about compact chambers, carrier construction,
Property 4, essential surjectivity, chamber finiteness, or a dimension bound
for raw chambers without a supplied carrier.  The dimension statement (4.1)
is only the subsequent K3 carrier corollary.
The theorem itself proves a field collapse, not a classification of the
rational chambers.
