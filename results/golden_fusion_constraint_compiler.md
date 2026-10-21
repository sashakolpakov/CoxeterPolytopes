# Exact golden vertex-fusion constraint compiler

## Outcome

There is now a reusable exact constraint compiler for the global step that
the local golden Lannér census does not see. Its input is

- a simple-polytope nerve, given by its maximal facet sets;
- a finite domain in \(\{2,3,5\}\) for every intersecting facet pair;
- a finite domain in \(\mathbf Z[\phi]\) for every dotted pair; and
- optionally, selected labelled placements of the eight order-four/five
  golden Lannér seeds relevant to missing four-/five-faces (the exhaustive
  class also has four rank-three triangles, which cannot occur in this layer).

It produces exact, denominator-free rank-one equations and oriented sign
conditions at every resolved finite vertex. It has both a finite-domain
propagator and a symbolic two-rational-component form suitable for an SMT or
CSP back end. This is a bounded implementation result, not a new compact
\(\mathbb H^9\) example and not a nonexistence theorem.

The implementation is
`coxeter_polytopes.geometric.golden_fusion`, with calibrations in
`tests/test_golden_fusion.py`.

## 1. Denominator-free vertex fusion

Let \(I\) be the nine facets through a finite vertex of a proposed
14-facet polytope, and let \(J\) be its five-facet complement. Order the
Gram matrix as

\[
G=\begin{pmatrix}H&V\\V^t&C\end{pmatrix},\qquad H=G_I.
\]

The assigned finite labels make \(H\) constant. It must be positive definite
at both real places. Put

\[
\Delta=\det H,
\qquad A=\operatorname{adj}(H),
\]

and, for \(j,k\in J\), define the Schur numerators

\[
D_j=2\Delta-v_j^tAv_j,
\qquad
R_{jk}=\Delta h_{jk}-v_j^tAv_k.                 \tag{1}
\]

The actual Schur entries are \(D_j/\Delta\) and \(R_{jk}/\Delta\). Since
\(H\) is positive definite at both places, \(\Delta\) is positive at both
places. The complete fusion constraints are

\[
\begin{aligned}
D_j^{(0)}&<0, & D_j^{(1)}&>0,\\
R_{jk}^2&=D_jD_k, & R_{jk}^{(0)}&<0.
\end{aligned}                                                   \tag{2}
\]

Here \((0)\) is the hyperbolic embedding and \((1)\) the conjugate
embedding. There are only \(\binom52=10\) golden equations per vertex,
equivalently 20 rational polynomial equations after writing every element as
\(a+b\phi\) and reducing \(\phi^2=\phi+1\).

### Why the unsquared signs matter

The equations \(R_{jk}^2=D_jD_k\) alone only fix absolute values. At the
identity place, put \(a_j=\sqrt{-D_j}>0\). The last inequality in (2) then
forces

\[
R_{jk}=-a_ja_k.
\]

Consequently the entire numerator Schur matrix is \(-aa^t\), not merely a
matrix whose principal two-by-two minors vanish. In particular all its
minors vanish and it has the required negative rank-one inertia. Algebraic
conjugation preserves rank; the positive conjugate diagonals then make the
conjugate Schur matrix positive semidefinite of rank one.

This is exactly the vertex-anchor orientation condition already used in the
two-Gram presentation. Choose an outside facet \(j_0\) as anchor. The
coefficient of an outside root \(r_k\) along the anchor direction is

\[
x_{j_0k}
=\frac{(S_I)_{j_0k}}{(S_I)_{j_0j_0}}
=\frac{R_{j_0k}}{D_{j_0}}>0.                    \tag{3}
\]

Thus the unsquared negative \(R\)-sign is precisely the assertion that all
five outside facets have one coherent time orientation.

### The ten-wall determinant square class

There is a useful corollary that does not require the unknown pair product
\(h_{jk}\). Since

\[
D_j=\det G_{I\cup\{j\}},
\]

rank one gives

\[
\frac{D_j}{D_a}
=\left(\frac{R_{ja}}{D_a}\right)^2.             \tag{4}
\]

Thus the five nonzero ten-wall anchor determinants must all lie in one class
of \(F^*/F^{*2}\), where \(F=\mathbf Q(\sqrt5)\). The compiler has an exact
square-root decision in \(F\) and audits this condition as soon as two
outside incidence rows are separately resolved. It can therefore discard an
outside attachment before any complement-pair fusion value is selected.

## 2. CSP propagation

For each vertex with resolved internal labels, the compiler stores
\(H,\Delta,A\) exactly in \(\mathbf Q(\sqrt5)\). As its five incidence rows
\(v_j\) resolve, it performs the following finite propagation:

1. compute every available \(D_j\) and reject unless its two place signs are
   \((-,+)\);
2. reject any two resolved outside attachments whose \(D\)-ratio is not a
   square in \(F\), even if other incidence rows remain unresolved;
3. for each of the ten complement pairs, test every domain value of
   \(h_{jk}\) against both the exact equation and the unsquared identity
   sign in (2);
4. delete every failing value, or return a contradiction if the domain
   becomes empty; and
5. iterate, since a deletion can resolve rows for other vertices.

The useful feature is that, after the incidence rows are fixed, each
remaining \(h_{jk}\) occurs only in its own \(R_{jk}\). The nonlinear global
rank condition therefore becomes ten small table filters rather than a
14-by-14 determinant calculation.

If the rows are not yet resolved, `symbolic_systems()` emits expressions
\(a+b\phi\), allowed coefficient tables for every domain variable, ten
golden equations, and all 20 sign conditions. A solver may use Boolean
selectors for the coefficient tables and split each golden equality into
its two rational components.

This formulation complements the vertex-plus-anchor kernel compiler. The
kernel equations are conditionally linear and are likely preferable for a
single global SMT solve; the five-by-five vertex formulation is symmetric,
local to the nerve, and particularly effective as exact domain propagation
and as an independent verifier.

## 3. Finite Lannér fusion layer

The five order-four and three order-five unlabelled seed types expand to
small exact labelled tables:

\[
54\quad\text{placements on a fixed four-set},
\qquad
240\quad\text{placements on a fixed five-set}.
\]

Selecting a placement intersects all six or ten internal pair domains with
its required labels. Conflicts terminate immediately.

There is also an overlap test. For every pair of selected seeds whose union
Gram entries are resolved, the union is a principal submatrix of the target
Gram. Hence it may have at most one negative identity-place direction and no
conjugate negative direction. This is a finite exact elimination on seed
overlaps; unresolved unions are deferred rather than guessed.

Together, the currently implemented eliminations are:

- exact seed/domain conflict;
- resolved seed-union inertia at both real places;
- pairwise ten-wall determinant square classes;
- the five diagonal Schur sign pairs at every ready vertex;
- ten oriented two-by-two Schur equations at every ready vertex; and
- iterative domain deletion across overlapping vertices.

## 4. Exact calibration

### Bugaenko's compact \(\mathbb H^8\) chamber

The compiler is dimension-independent, so the known 11-wall Bugaenko Gram
matrix gives a stringent calibration. All 41 eight-facet vertices compile;
all 41 pass. Each three-facet complement produces three golden equations,
six rational equations, and nine sign conditions.

At the first vertex, every ratio \(R_{j_0k}/D_{j_0}\) agrees exactly with the
corresponding coordinate computed by the pre-existing vertex-anchor audit,
and is positive. This directly tests (3), including the unsquared sign.

The unique dotted product is \(-2\phi\). If its domain is artificially
widened to

\[
\{-2\phi,-\phi^2\},
\]

both candidates pass the elementary dotted-pair arithmetic gates, but
vertex fusion deletes \(-\phi^2\) and retains only \(-2\phi\). This is a
small concrete example of genuinely global propagation.

### The local 14-wall countercertificate

For

\[
\Delta_3\times\Delta_3\times\Delta_1^3,
\]

all 128 nine-facet vertices pass the previously listed local ellipticity,
ADE, pairwise curve-intersection-sign/effectivity-sign, and Lannér tests. The
new compiler rejects the first tested vertex before complement-pair fusion:
the ten-wall determinants for outside facets 3 and 9 have ratio

\[
\frac{10+3\phi}{11},
\]

which is not a square in \(F\). If this early corollary is disabled, the full
oriented rank-one equation also rejects the pair \((3,7)\).

There is an even cheaper seed-overlap certificate. The two disjoint
\(A_4(3,5,3)\) seeds on facets \(0,\ldots,3\) and \(4,\ldots,7\) have union
inertias

\[
(6,2,0)\quad\text{and}\quad(8,0,0)
\]

at the identity and conjugate places. The two identity negative directions
already contradict any Lorentzian ambient Gram matrix.

## 5. Scope and next use

The compiler does not yet enumerate all 9/14 nerves, derive a complete
a-priori bound on dotted domains, select one Lannér placement for every
missing face, prove K3 effectivity, or certify that a resulting reflection
chamber has finite index. Its exact role is narrower: it converts the first
nonlocal Gram obstruction into a solver-ready finite constraint layer.

The next bounded experiment is therefore clear. Feed each surviving 9/14
nerve into the 54/240 seed tables, propagate the vertex equations after each
choice, and use the kernel formulation only for branches that survive these
cheap exact filters.

## Reproduction

```sh
PYTHONPATH=src python3 -m pytest -q tests/test_golden_fusion.py
```

The calibration currently contains seven tests, all exact.
