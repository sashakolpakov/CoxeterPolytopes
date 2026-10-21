# Ordinary physical chambers stop exactly at the anticanonical cusp

**Status (2026-08-11): proved classification-free surface-free dimension
bound and exact carrier theorem for the standard marked ordinary lane.**
Let

\[
 L=I_{1,n}=\mathbf ZH\oplus\bigoplus_{i=1}^n\mathbf ZE_i,
 \qquad K=-3H+E_1+\cdots+E_n,
 \qquad F=-K .                                         \tag{0.1}
\]

For every relatively reflection-maximal finite-volume integral Coxeter
chamber whose oriented facets have only the ordinary physical colors

\[
 (C^2,K.C)=(-1,-1)\quad\hbox{or}\quad(-2,0),          \tag{0.2}
\]

one has

\[
                              \boxed{n\leq9}.          \tag{0.3}
\]

For every surviving \(3\leq n\leq9\), the raw chamber has an explicit
smooth rational carrier \(S/\mathbf C\), with the prescribed canonical
marking and

\[
                         P=\mathbf P(\operatorname {Nef}S).             \tag{0.4}
\]

Thus the marked ordinary property-4 problem is closed in hyperbolic
dimensions at least three.  The proof uses a one-line Lorentzian escape for
\(n\geq10\), the already proved timelike rigidity for \(3\leq n\leq8\),
and the already proved characteristic-cusp saturation for \(n=9\).  It uses
no Nikulin--Vinberg reflection classification, ABSW theorem, chamber census,
or surface census.

Twice-primitive square-\(-4\) facets are not ordinary: they reverse the sign
of \(F\) on a wall.  They are exactly the remaining adapted-host boundary.

## 1. The finite-volume wall-cone lemma

Let \(V\) have signature \((1,N)\), let \(\mathcal C^+\) be a component of
its positive cone, and let

\[
 \mathcal W=\{x\in V:x.r_j\geq0\text{ for }1\leq j\leq f\}             \tag{1.1}
\]

be the homogeneous cone of a finite-sided hyperbolic polyhedron \(P\).
Choose \(A\in\mathcal C^+\cap\operatorname {int}\mathcal W\).

### Lemma 1.1 (spacelike escape and null extremality)

If \(P\) has finite hyperbolic volume, then:

1. \(\mathcal W\) contains no vector \(w\) with \(w^2<0\); and
2. every nonzero causal vector in \(\mathcal W\) lies in the future
   component; in particular every nonzero null vector in \(\mathcal W\)
   spans an extreme ray and hence is an ideal vertex of \(P\).

### Proof

Suppose first that \(w\in\mathcal W\) and \(w^2<0\).  For every \(t\geq0\),

\[
                  (A+tw).r_j=A.r_j+t(w.r_j)>0,          \tag{1.2}
\]

so \(A+tw\) stays strictly inside every wall inequality.  Its square is

\[
             (A+tw)^2=A^2+2tA.w+t^2w^2.               \tag{1.3}
\]

It is positive at \(t=0\) and tends to \(-\infty\).  At its first positive
zero \(t_0\), the vector \(z=A+t_0w\) is a nonzero future-null vector and
still satisfies every inequality strictly.  Since there are finitely many
walls, a whole open neighborhood of \([z]\) on the ideal sphere satisfies
all the inequalities.  A smaller ideal patch has a radial collar inside
\(P\).  In the ball model the volume density is a positive constant times

\[
                         (1-|x|^2)^{-N}\,dx,            \tag{1.4}
\]

whose radial integral diverges on every positive-measure ideal patch.
This contradicts finite volume.

Now let \(0\ne w\in\mathcal W\) be causal.  If \(w\) were past-directed and
not proportional to \(-A\), the Lorentzian plane spanned by \(A,w\) would
contain a spacelike positive combination \(A+tw\), still in the strict wall
cone by (1.2).  Proportionality to \(-A\) is also impossible: \(A.r_j>0\)
whereas \(w.r_j\geq0\).  Thus every causal vector in \(\mathcal W\) is
future-directed.

Finally let \(w\) be future-null and suppose its ray is not extreme.  Write

\[
                         w=u+v                           \tag{1.5}
\]

with nonzero, nonproportional \(u,v\in\mathcal W\).  The first paragraph
gives \(u^2,v^2\geq0\), and the second makes both future-directed.  Reverse
Lorentz Cauchy--Schwarz then says that a null sum is possible only when
\(u,v\) are proportional null vectors, a contradiction.  Hence
\(\mathbf R_{\geq0}w\) is extreme.  Extreme future-null rays of (1.1) are
precisely the ideal vertices of the projectivized chamber. \(\square\)

Only the first assertion is needed for the dimension bound.  The causal
orientation and extremality assertions remove the last extra premise in
the rank-ten cusp theorem: ordinary coloring itself forces \(F\) to be the
actual future cusp.

## 2. The raw ordinary dimension bound

### Theorem 2.1 (anticanonical escape)

Let \(n\geq3\), and let \(P\) be a finite-sided hyperbolic polyhedron in the
chosen positive component of \(I_{1,n}\otimes\mathbf R\).  Assume every
oriented facet class satisfies (0.2).  If \(P\) has finite volume, then

\[
                              n\leq9.                   \tag{2.1}
\]

This conclusion does not require reflection maximality, integrality of the
facet reflections, or existence of a surface.

### Proof

For every facet,

\[
 F.C=-K.C=
 \begin{cases}
  1,&C^2=-1,\\
  0,&C^2=-2.
 \end{cases}                                           \tag{2.2}
\]

Thus \(F\in\mathcal W\), the raw homogeneous wall cone.  But

\[
                              F^2=K^2=9-n.              \tag{2.3}
\]

For \(n\geq10\), \(F\) is spacelike.  Lemma 1.1 puts a strict open ideal
patch inside \(P\), contradicting finite volume.  Hence \(n\leq9\).
\(\square\)

The proof is deliberately surface-free.  It does not call \(F\) nef and
does not infer anything from a proposed Mori cone.  It uses only the
canonical wall colors and the finite-volume geometry of their homogeneous
halfspaces.

## 3. Exact carrier closure in every surviving rank

We now add the genuine group hypothesis:

> the side-reflection group of \(P\) has no proper reflection overgroup
> inside \(O^+(I_{1,n})\).

### Theorem 3.1 (ordinary P-first exact carrier)

Under (0.1)--(0.2), relative reflection maximality, and finite volume, for
every \(n\geq3\):

1. \(n\leq9\);
2. \(P\) is conjugate by an integral reflection isometry fixing \(K\) to
   the standard simplex with facets
   \[
    E_n,\quad E_i-E_{i+1}\ (1\leq i<n),\quad
    H-E_1-E_2-E_3;                                    \tag{3.1}
   \]
3. there is a smooth projective rational surface \(S/\mathbf C\) with a
   marking \(\operatorname {NS}(S)\cong I_{1,n}\) carrying \(K_S\) to
   \(K\), realizing all facets as distinct effective irreducible curves,
   and satisfying (0.4).

### Proof

Theorem 2.1 removes \(n\geq10\).
For every surviving \(n\leq9\), equation (2.2) puts the nonzero causal
vector \(F\) in \(\mathcal W\).  Lemma 1.1 makes it future-directed.
Consequently \(K=-F\) is past-directed and

\[
                              K\notin P.                \tag{3.0}
\]

Thus canonical escape is a conclusion of the raw finite-volume hypotheses,
not an extra premise.

For \(3\leq n\leq8\), the
[timelike anticanonical-vertex theorem](rational_pfirst_timelike_anticanonical_vertex_closure.md)
proves exactly assertions 2--3.  Its classification-free reflection descent
first carries the oriented chamber to the standard simplex.  The image of
\(F\) satisfies every standard homogeneous inequality; the
reduced-characteristic inequality then forces that image to equal \(F\).
Thus the conjugating reflection fixes \(K\), without an unproved initial
sheet choice.  Blowing up the first \(n\) points of the flex-cubic base
cluster gives a pencil in class \(-K\) with no fixed curve.  This
independently proves \(-K\) nef, and maximal-reflection return gives the
entire nef cone.

Let \(n=9\).  Equation (2.2) puts the primitive characteristic null vector

\[
                     F=(3,-1,\ldots,-1)               \tag{3.2}
\]

in the homogeneous wall cone.  Lemma 1.1 makes it future-directed and its
ray an ideal vertex of \(P\); no anticanonical-cusp or orientation premise
is being added.  The
[rank-ten characteristic-cusp theorem](rational_pfirst_anticanonical_cusp_saturation_closure.md)
then identifies \(P\), by a reflection fixing \(F\) and \(K=-F\), with the
same standard simplex (3.1).  Resolving all nine basepoints of the flex
cubic pencil realizes \(F\) as a genus-one fibre, the last exceptional as
the unique section, and the other nine walls as the components of its
affine fibre.  The cited theorem proves directly that these are all negative
curves and returns the full nef cone.

This covers every surviving \(n\). \(\square\)

## 4. What is now closed, and the exact remaining boundary

For the standard plane-type canonical marking, the entire ordinary
physical maximal lane in dimensions \(n\geq3\) is now P-first:

\[
\boxed{
 \begin{array}{c}
 \text{finite volume + relative reflection maximality}\\
 \text{+ ordinary canonical colors}
 \end{array}
 \Longrightarrow
 \begin{array}{c}
 n\leq9,\\
 \text{explicit rational carrier and exact same-}\!P\text{ return}.
 \end{array}}
\]

The trichotomy is intrinsic:

\[
\begin{array}{c|c|c}
3\leq n\leq8&F^2>0&\text{finite anticanonical vertex},\\
n=9&F^2=0&\text{characteristic anticanonical cusp},\\
n\geq10&F^2<0&\text{spacelike wall-cone escape; infinite volume}.
\end{array}                                             \tag{4.1}
\]

A twice-primitive physical boundary \(B^2=-4\), \(K.B=2\), has

\[
                              F.B=-2<0.                 \tag{4.2}
\]

It therefore ejects \(F\) from the wall cone and destroys the escape
argument at exactly one sign.  The unresolved property-4 problem is now
sharply the exponent-two adapted-host lane with such walls.  In that lane
one must still obtain either:

1. a polytope-intrinsic nef/semiample certificate for
   \(D=-(K+B/2)\); or
2. a genuinely reflection-maximal nonrealizable packet.

No claim is made here for arbitrary nonstandard characteristic orbits,
dimensions \(n<3\), twice-\((-4)\) facets, or other commensurable hosts.

## 5. Frozen dependencies and replay

The positive branches are pinned at:

    timelike report
      6450ab0629c539fcb3fff00e6aa2bdd0ddf194f17bdaaa2c1f2265d2b4aebb29
    timelike hostile audit
      2d7eddae058f81f06fded910b4d8ed45219b30e011af5be6a232f20655df0499
    cusp report
      9d4142c8c4647bfced451ac1ff29184a6a278d9a51cdb1a0638276f3eb813032
    cusp hostile audit
      1cc36249914c82ba42145e64cb0e1a43a888d4ac0c02d6b02084cade3ea498d8

The branch arithmetic and dependency pins are replayed by
[the exact certifier](../scripts/certify_rational_pfirst_ordinary_dimension_bound_and_exact_carrier.py).
