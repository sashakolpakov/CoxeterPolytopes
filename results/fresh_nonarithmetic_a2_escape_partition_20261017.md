# All actual S6-facet profiles preserve the 55-wall escape

The exact 55-wall chamber has a strict spacelike escape (E_{55}). Every additional Coxeter wall which induces an actual facet of its face

\[
S=\{0,1,2,3,5,6\}=H_4+A_2
\]

has strictly negative pairing with (E_{55}), provided its halfspace and the retained chamber have a common strict future interior. This includes all seven possible contact profiles with (S), every finite label permitted by those profiles, and arbitrary real coefficient fields.

This is a fixed-face obstruction. It does **not** exclude completing the ambient chamber by walls which do not induce facets of (S). It also does not exclude closing the entire three-dimensional face while separate ambient walls cut the escape. No finite-volume example or general noncompletion theorem is claimed.

The [primary checker](../scripts/fresh_nonarithmetic_a2_escape_partition_20261017.py) and [certificate](fresh_nonarithmetic_a2_escape_partition_20261017.json) passed generation and fresh saved replay. They reuse frozen physical algebraic normal recipes and outward rational interval arithmetic. Every decision is an exact sign certificate; the earlier numerical scouts do not supply an acceptance test. The narrow conceptual read by `frontier_literature` checked the maximum principle, the seven-profile reduction and the common-interior scope; it was not a separate arithmetic replay.

## Retained data and complete profile reduction

The retained geometry is the [55-wall exact cone](fresh_nonarithmetic_nonright_s6_20261016_cone.json), with [independent complete-cone review](fresh_nonarithmetic_nonright_s6_20261016_cone_review.md). The checker uses its literal escape, multiplied by a saved positive rational power of two. Each old normal is positively normalized to square two. Thus normalization changes neither its halfspace nor any sign.

If a new square-two normal (r) induces an actual facet of the hyperbolic (S)-face, the Gram matrix of (S\cup\{r\}) is positive definite. Equivalently its projection to (S^\perp) has positive norm. The (H_4) inverse is entrywise positive and has diagonal

\[
14+6\sqrt5,\quad21+9\sqrt5,\quad10+4\sqrt5,\quad3+\sqrt5.
\]

Every diagonal exceeds two. A nonzero Coxeter product has absolute value at least one, so any nonright (H_4) contact already gives projection cost greater than two. Therefore (r) is right to all four (H_4) roots. The (A_2) inverse is

\[
J_{A_2}=\frac13\begin{pmatrix}2&1\\1&2\end{pmatrix}.
\]

Two nonzero contacts have cost at least two. A single contact of label six or greater also has cost at least two. Parallel and dotted products are larger still. Hence the complete list is the zero profile, or a single label (3,4,5) at node5 or node6. This is the same complete projection lemma recorded in the [A2 profile package](fresh_nonarithmetic_h3_a2_profiles_20261016.md), here checked again against the literal support.

Write (\phi=(1+\sqrt5)/2), and normalize a nonzero contact by (c=2\cos(\pi/m)\in\{1,\sqrt2,\phi\}). Then

\[
r=c(p_i+z),\qquad z\in S^\perp,
\qquad
p_5=-\frac{2n_5+n_6}{3},\quad
p_6=-\frac{n_5+2n_6}{3},
\]

and the required residual square is

\[
z^2=\frac2{c^2}-\frac23\ge\frac2{\phi^2}-\frac23>\frac9{100}.
\]

All later identities are expressed in the four exact coordinates

\[
\xi_j=\langle z,\pi_S(n_j)\rangle,
\qquad j\in(7,10,30,47),
\]

where \(\pi_S\) is orthogonal projection to \(S^\perp\). The projected Gram is nonsingular, so these are genuine coordinates on the Lorentz four-space. Its inverse (J) satisfies (z^2=\xi^TJ\xi). Every inverse pivot is certified nonzero by outward rational intervals.

## Zero profile

Let (v) be the actual retained finite vertex with support (S\cup\{8,10,44\}), constructed as the projection of the old center onto the orthogonal complement of those nine normals. Its Gram is positive definite, its norm and future pairing are strictly negative, and all46 nonincident old-wall inequalities are checked strictly.

The checker establishes the exact identity

\[
\pi_S(E_{55})=\lambda_7\pi_S(n_7)+\lambda_{47}\pi_S(n_{47})
+\lambda_{51}\pi_S(n_{51})+\lambda_vv,
\]

with all four coefficients strictly positive. Their exact definitions are the solution of a nonsingular four-dimensional linear system; the certificate gives enclosing rational intervals.

For completeness, finite-vertex persistence does not assume that the old selected center survives. Let (X) be any common strict future interior point and write

\[
X=av+\sum_{i\in\operatorname{supp}(v)}b_i n_i.
\]

Here (a>0), and (b=G^{-1}(\langle X,n_i\rangle)\le0), because the positive-definite acute Gram has entrywise nonnegative inverse. For an added acute normal (r),

\[
0>\langle X,r\rangle
=a\langle v,r\rangle+\sum b_i\langle n_i,r\rangle
\ge a\langle v,r\rangle.
\]

Consequently (\langle v,r\rangle<0). When (r\perp S), the positive decomposition above gives (\langle E_{55},r\rangle<0). This is the only part of the certificate needing the common strict future interior assumption.

## Node6 profiles

The checker establishes

\[
\pi_S(E_{55})=\sum_{j\in\{7,30,47,52\}}\lambda_j\pi_S(n_j),
\qquad\lambda_j>0.
\]

Old-wall acuteness of (p_6+z) therefore yields

\[
\langle p_6+z,E_{55}\rangle
\le\langle p_6,E_{55}\rangle-
\sum_j\lambda_j\langle p_6,n_j\rangle<0.
\]

The last bound is approximately \(-0.0002055050052\) in the saved rational scaling of the escape. The exact upper endpoint is negative. This excludes all three node6 profiles, and in fact every positive contact scale in this direction, without using the norm equation or any additional angle test.

## Node5 profiles: a bounded Lorentz quadratic problem

Suppose (\langle p_5+z,E_{55}\rangle\ge0). Keep only acuteness to old walls (7,10,30,47,53), together with this escape-cutting inequality. These six necessary inequalities define an outer polytope (P\subset\mathbb R^4); no claim that every point of (P) is a Coxeter wall is needed.

The first four coefficient rows are the exact unit rows in the chosen coordinates. Write the last two rows as (a_{53}) and (-a_E). The checker certifies

\[
a_E-2^{-15}a_{53}>0
\]

entrywise. Thus the six rows have a strictly positive dependence, with coefficients

\[
(a_E-2^{-15}a_{53},\ 2^{-15},\ 1).
\]

Since the first four rows span the dual space, their common recession cone is zero, so (P) is bounded.

All15 choices of four active rows have certified nonsingular matrices. Exactly nine solutions are feasible; each has exactly four active rows and the remaining two inequalities strict. Their barycenter is strict for all six inequalities. This proves that (P) is a full-dimensional simple polytope and that the vertex census is complete. Two vertices form an edge exactly when they share three active rows, giving18 edges and degree four at every vertex.

Every vertex satisfies

\[
\xi^TJ\xi<\frac1{50}.
\]

Seventeen edge directions have strictly positive square, so the quadratic is convex on them and has its maximum at an endpoint. The remaining edge has negative direction square, but its stationary parameter is less than (-102), outside the segment ([0,1]); the quadratic is monotone along that edge. Thus every edge also has norm below (1/50).

There is no untested maximum in a higher-dimensional face. Every linear subspace of a Lorentz four-space of dimension at least two contains a positive-square direction. At a relative-interior maximum, sufficiently short displacements in both signs of that direction would remain in the face, while the average of their quadratic values exceeds the alleged maximum. Therefore a maximum on the compact polytope occurs on a vertex or edge. It follows that

\[
z^2<\frac1{50}<\frac9{100}\le\frac2{c^2}-\frac23,
\]

contradicting the required norm. This excludes all three node5 profiles, including the relaxed continuum (1\le c\le\phi).

## Consequence for construction

No coherent additional wall which is an actual facet of the retained (H_4+A_2) three-face can cut or touch the current (E_{55}). Any finite family consisting only of such walls leaves a strict spacelike escape and hence cannot complete the ambient chamber to finite volume. Successful ambient completion must include at least one wall outside this actual-face-facet domain. Coordinated closure of the three-face may still be useful in combination with those separate ambient cutters.

This result replaces a finite sample failure with an all-profile obstruction for this specific escape and retained face. It does not address arbitrary changes to retained walls or other source faces.
