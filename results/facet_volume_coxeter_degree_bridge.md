# Facets bound hyperbolic volume and Coxeter degree

## Theorem

Let \(n\ge 2\), and let \(P\subset \mathbb H^n\) be a full-dimensional
finite-volume convex polytope with exactly \(f\) facets.  Redundant defining
halfspaces may instead be counted; this only weakens the bound.  Put

\[
U_n(f)=
\begin{cases}
\displaystyle {f-m\choose m}+{f-m-1\choose m-1},&n=2m,\\[6pt]
\displaystyle 2{f-m-1\choose m},&n=2m+1,
\end{cases}
\]

and let \(v_n\) be the volume of the regular ideal \(n\)-simplex.  Then

\[
 \boxed{\operatorname{Vol}(P)
 \le {U_n(f)-1\choose n}v_n.}
 \tag{1}
\]

This applies without simplicity, and to both compact and noncompact
finite-volume polytopes.  In the noncompact case every point of the
projective closure on \(\partial_\infty\mathbb H^n\) is an ideal **vertex**;
there is no positive-dimensional ideal face.

If \(P\) is a Coxeter chamber, \(W(P)\le W_0\), and \(Q\) is a fundamental
chamber (or measurable fundamental domain) of the lattice \(W_0\), then

\[
 [W_0:W(P)]={\operatorname{Vol}(P)\over\operatorname{Vol}(Q)}
 \le { {U_n(f)-1\choose n}v_n\over\operatorname{Vol}(Q)}.
 \tag{2}
\]

Thus, for a fixed overgroup \(W_0\), bounded facet number gives bounded
Coxeter/Schreier-dessin degree.  Uniformity while \(W_0\) varies additionally
requires an effective positive lower bound for \(\operatorname{Vol}(Q)\);
(1) alone does not provide that denominator.

## Proof, including the cusp point

Use the Klein ball \(B\subset\mathbb R^n\).  Extend the \(f\) hyperbolic
halfspaces defining \(P\) to affine halfspaces, and let their Euclidean
intersection be \(R\), so \(P=B\cap R\).

Finite volume forces \(R\subset\overline B\).  Indeed, if some point of
\(R\) lay outside \(\overline B\), joining it to an interior point of \(P\)
and using a small ball in the interior of \(R\) would put a relatively open
patch of \(\partial B\) in \(R\).  The hyperbolic cone from the interior
point to that patch has infinite volume, a contradiction.  Consequently
\(R\) is bounded and \(\overline P=R\).  Since the Euclidean ball is strictly
convex, a point of \(R\cap\partial B\) cannot be a nontrivial convex
combination of vertices of \(R\); hence it is itself a vertex.  This proves
both that the closure is an honest Euclidean polytope and that its infinity
set consists only of ideal vertices.

Let \(N\) be the number of vertices of \(R\).  Polarity and McMullen's Upper
Bound Theorem give

\[
 N\le f_{n-1}(C(f,n))=U_n(f).
\]

This is valid for nonsimple polytopes.  (The elementary substitute
\(N\le {f\choose n}\) also suffices.)

Choose one closure vertex \(a\).  A pulling triangulation, with \(a\) pulled
first, triangulates \(R\) without new vertices and makes every top simplex
contain \(a\).  Hence it has at most \({N-1\choose n}\) top simplices.  Every
one is a hyperbolic simplex whose vertices are finite or ideal.  By the
Haagerup--Munkholm maximal-simplex theorem its volume is at most \(v_n\),
with equality only for a regular ideal simplex.  Additivity across the
triangulation proves (1).  Both groups in (2) have finite covolume, so the
usual fundamental-domain decomposition proves the exact volume/index ratio
and (2).

## Low-dimensional checks and sharper constants

- For \(n=2\), Gauss--Bonnet gives the sharp bound
  \[
  \operatorname{Area}(P)=(f-2)\pi-\sum_i\theta_i\le(f-2)\pi.
  \]
  It is attained by an ideal \(f\)-gon.  Thus one may replace (1) by
  \(V(2,f)=(f-2)\pi\).

- For \(n=3\), \(N\le2f-4\).  Pulling from a vertex and triangulating the
  boundary sphere shows that the number of tetrahedra is at most
  \(2N-7\le4f-15\): a triangulated \(N\)-vertex sphere has \(2N-4\)
  triangles and at least three are incident to the pulled vertex.  Therefore
  \[
  \operatorname{Vol}(P)\le(4f-15)v_3,
  \qquad v_3=3\Lambda(\pi/3)\approx1.0149416064.
  \]
  At \(f=4\) this is the sharp regular-ideal-tetrahedron bound.

- In every dimension, \(f=n+1\) gives \(U_n(f)=n+1\) and (1) reduces to
  \(\operatorname{Vol}(P)\le v_n\), the sharp simplex theorem.

## Scope audit

No Coxeter angles, arithmeticity, or simplicity are used in (1).  Finite
volume is essential: a projective halfspace intersection extending outside
the Klein ball can meet the sphere in an open patch and have infinite
hyperbolic volume.  Dimension \(n=1\) is excluded because an interval with an
ideal endpoint has infinite length.  The triangulation uses only genuine
finite and ideal closure vertices; it introduces neither hyperideal nor
Steiner vertices.

## Primary inputs

- P. McMullen, *The maximum numbers of faces of a convex polytope*,
  Mathematika **17** (1970), 179--184,
  <https://doi.org/10.1112/S0025579300002850>.
- U. Haagerup and H. J. Munkholm, *Simplices of maximal volume in hyperbolic
  n-space*, Acta Math. **147** (1981), 1--11,
  <https://doi.org/10.1007/BF02392865>.

