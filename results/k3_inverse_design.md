# Inverse-designing the K3 surface from the Coxeter constraints

## Outcome

The surface should be constructed **after** the Coxeter combinatorics, not
searched independently. In the real-quadratic K3 model the chain of data is

\[
\text{nerve and labels}
\longrightarrow (M,h,\{r_i\})
\longrightarrow N=\operatorname{Tr}_{F/\mathbb Q}(\alpha h)
\longrightarrow T_{\mathrm{tr}}=N^\perp
\longrightarrow X_{T_{\mathrm{tr}}}.
\]

Here \(F=\mathbb Q(\sqrt5)\), \(\mathcal O_F=\mathbb Z[\phi]\), and
\(\alpha=-\phi/\sqrt5\). For an \(\mathbb H^9\) target, \(M\) has rank 10
over \(\mathcal O_F\), \(N\) has rank 20, and \(T_{\mathrm{tr}}\) has rank
two. Thus a successful design does not leave a positive-dimensional K3
search: it determines a singular K3 surface, up to the standard finite marked
ambiguity, through its positive-definite transcendental lattice
\(T_{\mathrm{tr}}\). An Inose model can then be recovered from
\(T_{\mathrm{tr}}\).

This makes the example track a finite lattice-and-diagram problem. It is not a
search over equations of surfaces.

## What the K3 gives for free

The phrase “free hyperbolicity” is accurate but incomplete.

1. The Hodge index theorem gives
   \(\operatorname{NS}(X)_\mathbb R\) signature \((1,\rho-1)\), hence a
   hyperbolic positive cone.
2. Riemann--Roch and the K3 Weyl theorem identify the ample cone as a chamber
   of the complete \((-2)\)-root arrangement.
3. The K3 period theorem realizes every admissible even lattice of signature
   \((1,\rho-1)\) that embeds primitively in the K3 lattice.

The K3 does **not** give compactness. Nor does it say that a prescribed finite
list of mirrors is complete. Those are the Vinberg parts of the problem.

Andreev's theorem is a complete realization theorem only in hyperbolic
dimension three. In dimension nine its role is replaced by Vinberg's exact
Gram conditions, elliptic vertex diagrams, Lannér missing faces, rank and
inertia, and a separate chamber-completeness certificate.

## The curve-lift rule

Normalize every field root by

\[
h(r_i,r_i)=2.
\]

Its trace-lattice square is

\[
q(r_i,r_i)=\operatorname{Tr}(2\alpha)=-2.
\]

For two finite-angle facets, the only equal-norm labels available over
\(\mathbb Q(\sqrt5)\) are \(2,3,5\). Their field and trace products are

| Coxeter label | \(h(r_i,r_j)\) | \(q(r_i,r_j)\) |
|---:|---:|---:|
| \(2\) | \(0\) | \(0\) |
| \(3\) | \(-1\) | \(1\) |
| \(5\) | \(-\phi\) | \(1\) |

Thus, after choosing one trace root above each restricted wall, a 3-edge and
a 5-edge have the same curve-level shadow: two smooth rational curves meeting
once.  This one-root shadow is necessary but has forgotten the golden
decoration.

There is a stronger complete statement.  In an even golden trace lattice,
every restricted wall has the two orthogonal trace roots

\[
r_i,\qquad (\phi-1)r_i.
\]

Once the norm-two wall list is proved to be a complete chamber, the K3 Weyl
framing makes both roots irreducible effective curve classes automatically.
The full paired curve lift distinguishes the two labels: a 3-edge lifts to
\(A_2\sqcup A_2\), while a 5-edge lifts to \(A_4\).  The universal shell and
Weyl-framing proof is in `results/k3_eigen_nef_factory_theorem.md`.

## The missing datum: the golden correspondence

The correct object to design is

\[
(X,\tau,\delta_1,\ldots,\delta_n),
\qquad \tau^2-\tau-I=0,
\]

where \(\tau\) is self-adjoint on \(\operatorname{NS}(X)\). Besides the
ordinary curve matrix

\[
Q_{ij}=q(\delta_i,\delta_j),
\]

record the symmetric integral **correspondence matrix**

\[
B_{ij}=q(\tau\delta_i,\delta_j).
\]

Here \(\tau\) is extra divisor-lattice structure, not necessarily an
automorphism of \(X\). Rationally it can be represented by an algebraic
correspondence supported on products of divisors in \(X\times X\); the
inverse-design calculation itself only needs its integral action on
\(\operatorname{NS}(X)\).

If \(h(r_i,r_j)=a+b\phi\), then

\[
Q_{ij}=-(a+b),\qquad B_{ij}=-(a+2b),
\]

and conversely

\[
h(r_i,r_j)
=(B_{ij}-2Q_{ij})+(Q_{ij}-B_{ij})\phi.
\]

For the complete paired root list
\((r_i,(\phi-1)r_i)\), the intersection matrix is

\[
\widetilde Q=
\begin{pmatrix}Q&B-Q\\B-Q&2Q-B\end{pmatrix}.
\]

Consequently a completed eigen-nef facet system must satisfy
\(Q_{ij}\le B_{ij}\le2Q_{ij}\) for distinct walls.  This is implemented by
`paired_trace_curve_gram`.

Thus the two integral matrices recover the full golden Gram matrix. On a
finite edge,

\[
\begin{array}{c|ccc}
m&2&3&5\\ \hline
(Q_{ij},B_{ij})&(0,0)&(1,1)&(1,2).
\end{array}
\]

For a dotted edge \(h(r_i,r_j)=a+b\phi\), the necessary curve condition is

\[
Q_{ij}=-(a+b)\in\mathbb Z_{\ge0}.
\]

This converts the inverse problem into coupled integral equations

\[
\tau^2=\tau+I,\qquad \tau^tQ=Q\tau,\qquad B=Q\tau
\]

in a divisor basis, together with prescribed entries of \(Q\) and \(B\).
The exact conversions are implemented by the functions
trace_curve_gram, correspondence_gram, and recover_field_gram in
coxeter_polytopes.geometric.k3_design.

On a redundant set of wall generators one should not impose a rank condition
on \(Q\) itself.  Instead, the two matrices canonically present the trace
root span through

\[
\mathcal D(Q,B)=\begin{pmatrix}Q&B\\B&Q+B\end{pmatrix},
\]

the Gram matrix on \((\delta_i,\tau\delta_i)\).  For 14 walls in
\(\mathbb H^9\), \(\mathcal D\) must have rank 20 and inertia
\((1,19,8)\).  Its saturated integral radical gives the relations, and its
quotient gives the canonical rank-20 root-span lattice with its descended
golden action.  All possible ambient lattices are then the finite collection
of \(\tau\)-stable even overlattices.  This stronger presentation gate is
implemented in `coxeter_polytopes.geometric.k3_presentation` and proved in
`results/k3_root_presentation_gate.md`.

### Two-Gram reconstruction lemma

Let \((N,q)\) be a nondegenerate integral symmetric lattice and let
\(\tau\in\operatorname{End}_{\mathbb Z}(N)\) be \(q\)-self-adjoint with
\(\tau^2=\tau+I\). Then \(\phi x:=\tau x\) makes \(N\) a torsion-free
\(\mathbb Z[\phi]\)-module, and there is a unique symmetric
\(\mathbb Z[\phi]\)-bilinear form satisfying

\[
q(x,y)=\operatorname{Tr}_{F/\mathbb Q}(\alpha h(x,y)).
\]

It is

\[
h(x,y)=
\bigl(q(\tau x,y)-2q(x,y)\bigr)
+\bigl(q(x,y)-q(\tau x,y)\bigr)\phi.
\]

Indeed, if \(h(x,y)=a+b\phi\), the two traces \(q(x,y)\) and
\(q(\tau x,y)\) are \(-(a+b)\) and \(-(a+2b)\), so they recover \(a,b\)
uniquely. Self-adjointness gives symmetry, while replacing \(x\) by
\(\tau x\) and using \(\tau^2=\tau+I\) proves
\(h(\tau x,y)=\phi h(x,y)\). The radicals of \(q\) and \(h\) agree. If
\(q\) is even, then \(h(x,x)\in2\mathbb Z[\phi]\); unimodularity does not
follow and requires a separate discriminant calculation.

Thus the two-Gram package is not merely necessary: together with the
self-adjoint golden action it is equivalent to the integral golden lattice
datum. What remains insufficient is the passage from that lattice datum to
effective irreducible curves, a common nef chamber, and a complete compact
Vinberg chamber.

## Vertex contraction constraint

Let \(I\) be the facets through a finite vertex. Their identity-place field
Gram matrix is positive definite because the Coxeter subdiagram is elliptic.
The conjugate-place matrix is also positive definite by arithmeticity. The
two embeddings of \(\alpha\) are both negative, so the trace intersection
matrix

\[
Q_I=\bigl(\operatorname{Tr}(\alpha h(r_i,r_j))\bigr)_{i,j\in I}
\]

is negative definite. Geometrically, these are precisely the rational curves
simultaneously contracted by the nef class represented by the vertex.

After forgetting whether a single curve-intersection edge came from a
Coxeter 3-edge or 5-edge, every finite vertex graph must therefore be a
disjoint union of ADE diagrams. This gives the useful dictionary

\[
\text{elliptic Coxeter vertex}
\quad\rightsquigarrow\quad
\text{contractible ADE curve configuration}.
\]

Once a vertex has already been certified as a labeled elliptic Coxeter
diagram with labels in \(\{2,3,5\}\), this ADE condition is automatic from
the finite Coxeter classification: \(H_3,H_4,I_2(5)\) forget respectively to
\(A_3,A_4,A_2\).  It remains useful as a prelabeling filter and as the
geometric contraction interpretation, but it is not an additional rejection
test after labeled ellipticity.

The exact audit is implemented in
coxeter_polytopes.geometric.k3_design. On Bugaenko's chamber it recovers the
trace intersections \(0,1,2\); all 41 eight-facet vertex matrices are
negative definite, and the field Gram matrix has inertias
\((8,1,2)\) and \((9,0,2)\) at the two real places. The two integral matrices
\(Q,B\) reconstruct the archived golden Gram matrix entry by entry.

## Forced finite cores in the 14-facet target

The 9/14 combinatorics force at least one missing face of order four or five.
In the norm-two golden setting its standard finite subdiagram cannot be
\(B_3\) or \(F_4\), because those require a label-four edge. Hence every
candidate contains an \(H_3\) or \(H_4\) core. Restriction of scalars gives
the exact integral lifts

\[
H_3\longmapsto D_6(-1),\qquad
H_4\longmapsto E_8(-1).
\]

The first is primitive in any even trace lattice; the second is unimodular
and splits off. Thus the \(H_4\) branch begins with
\(N=E_8(-1)\oplus N_{12}\), reducing the undetermined golden part from rank
ten to rank six. Exact bases, the three possible order-five Lanner
extensions, and their complement norms are recorded in
`results/golden_k3_forced_cores.md`.

## The finite seed alphabet and the global fusion equation

The core dichotomy can be completed exactly.  Enumerating every graph on four
and five nodes, every \(3/5\) decoration, and imposing Lanner inertia at the
identity place plus positive definiteness at the conjugate place leaves eight
blocks:

\[
\begin{gathered}
A_4(3,5,3),\ A_4(5,3,5),\ D_4(5,3,3),\\
C_4(5,3,3,3),\ C_4(5,3,5,3),\\
A_5(5,3,3,3),\ A_5(5,3,3,5),\ D_5(5,3,3,3).
\end{gathered}
\]

Their doubled trace discriminants are respectively
\(11,19,16,31,55,4,20,16\) in absolute value.  Thus the example search starts
from a finite list of Lorentzian divisor seeds, not from arbitrary surfaces.

There is also a sharp global equation.  Choose any nine-facet vertex \(I\)
and let \(J\) be the five complementary facets.  If the completed golden Gram
matrix is written

\[
h=\begin{pmatrix}H&C\\C^t&D\end{pmatrix},\qquad H=h_I,
\]

then

\[
S_I=D-C^tH^{-1}C
\]

must have rank one, with inertias \((0,1,4)\) and \((1,0,4)\) at the two real
places.  Equivalently, all two-by-two minors of \(S_I\) vanish.  These
equations couple the five outside walls and are stronger than auditing the
five ten-wall extensions separately.

This global step is essential.  The exact formal model
\(\Delta_3^2\times\Delta_1^3\) passes the missing-face, local Lanner,
elliptic-vertex, trace-ADE, nonnegative pairwise curve-intersection, and
conjugate-positivity tests, but has field inertia \((9,5,0)\): its five local
negative directions never fuse.  It is therefore not a K3 configuration or a
hyperbolic polytope.  The census, Schur proof, and countercertificate are in
`results/golden_lanner_global_constraints.md`.

The rank-one condition now has a denominator-free exact compiler.  If
\(\Delta=\det H\), \(A=\operatorname{adj}(H)\), and \(v_j\) is the incidence
column of an outside wall, put

\[
D_j=2\Delta-v_j^tAv_j,\qquad
R_{jk}=\Delta h_{jk}-v_j^tAv_k.
\]

Then \(D_j^{(0)}<0<D_j^{(1)}\),
\(R_{jk}^2=D_jD_k\), and the common-side orientation requires
\(R_{jk}^{(0)}<0\).  In particular all five \(D_j\), which are the
ten-wall anchor determinants, lie in one square class over
\(\mathbb Q(\sqrt5)\).  The finite-domain and symbolic implementations are
recorded in `results/golden_fusion_constraint_compiler.md`.

The eight seeds themselves all survive exact K3 lattice theory.  Each has a
unique primitive K3 embedding orbit and no proper \(\tau\)-stable even
self-overlattice.  More concretely, the discriminant-19 rank-20 carrier
contains primitive copies of all eight, with negative-definite
\(\tau\)-stable complements and explicit isotropic gluing; the
determinant-11 path has the analogous \(H_{11}\) completion.  Thus seed
embedding is not where a candidate dies.  The remaining constraints are the
simultaneous labeled fusion, common nef orientation, extra-root closure, and
compactness; see `results/golden_lanner_k3_seeds.md`.

## The inverse construction

For a proposed compact \(d\)-polytope with \(n\) facets:

1. Enumerate a simple nerve satisfying the known missing-face restrictions.
2. Color finite edges only by \(2,3,5\), prescribe their integral pairs
   \((Q_{ij},B_{ij})\), and attach integral pairs to dotted edges.
3. Form \(\mathcal D(Q,B)\), require rank \(2(d+1)\), take its saturated
   radical quotient \(N_0\), and enumerate the finite set of
   \(\tau\)-stable even overlattices.
4. Recover the golden Gram matrix and impose Vinberg's equations: diagonal
   2, rank \(d+1\), identity inertia \((d,1)\), positive-definite conjugate
   form, elliptic Gram matrices at every vertex, and the vertex-anchor Schur
   and same-side equations.
5. Apply the curve-lift filter: all off-diagonal trace products must be
   nonnegative integers, and every finite vertex trace matrix must be negative
   definite (equivalently ADE in the simply-laced case).
6. Recover the integral \(\mathcal O_F\)-lattice \(M\) from each surviving
   trace overlattice and check every root reflection is integral.
7. Compute the discriminant form of each \(N\) and test primitive embedding in
   \(\Lambda_{K3}=U^3\oplus E_8(-1)^2\).
8. In dimension nine, identify the rank-two complement
   \(T_{\mathrm{tr}}=N^\perp\) and construct the corresponding singular K3
   \(X_{T_{\mathrm{tr}}}\).
9. Finally run exact norm-two Vinberg completion and compactness verification.

Steps 1--5 are the combinatorial design stage. Step 7 is where most proposed
diagrams should die. The surface itself appears only at step 8, after the
discrete data have essentially forced it.

## What remains nonautomatic

Passing the two-Gram, presentation, and ADE filters still does not prove:

1. that one of the finitely enumerated trace overlattices has the required
   rank-two discriminant-form complement;
2. that this lattice embeds primitively in the K3 lattice;
3. that the listed roots are walls of one common restricted chamber;
4. that every norm-two mirror has been accounted for; or
5. that the chamber is compact.

The first two are finite lattice problems. The last three require the exact
restricted-root/Vinberg certificate.

## Scope and the higher-dimensional boundary

A quadratic eigenspace of hyperbolic dimension \(d\) consumes divisor rank
\(2(d+1)\). Since a complex K3 has Picard rank at most 20, this construction
has the sharp internal ceiling

\[
d\le 9.
\]

So this is the right inverse factory for the \(\mathbb H^9\) boundary, but a
K3 surface cannot carry the same construction in \(\mathbb H^{10}\). To move
past nine one must either use a different surface with unbounded Picard rank
and a complete cone theorem (toric or more generally Mori-dream surfaces are
the clean candidates), or use an OG10 hyperkähler variety, whose \(b_2=24\)
permits a quadratic eigenspace in \(\mathbb H^{10}\). In both extensions the
same principle applies: start with the Vinberg diagram and solve its
divisor-lift equations; do not scan varieties first.
