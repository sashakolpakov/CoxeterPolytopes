# A fixed-realization metric gate for H9/14 records 43, 44, and 47

## Scope

This note supplies one exact homogeneous facet realization for each of the
three multiwedge polytopes and records inexpensive Lorentzian obstructions for
those particular projective realizations.  It does **not** exclude another
realization of the same combinatorial type.

## Integral facet covectors

All inequalities below have the form `ell_i >= 0`; their coefficient rows,
including the constant term, form a (14\times10) matrix (A_R).

For records 43 and 47 use affine coordinates

\[
(x,y,z,t_1,t_2,s_1,s_2,s_3,s_4).
\]

The common rows are

\[
\begin{array}{c|l}
i&\ell_i\\ \hline
0&x+y+z-1\\
1&y\\
2&z\\
3&s_1\\4&s_2\\5&s_3\\6&s_4\\
7&1-s_1-s_2-s_3-s_4\\
8&t_1\\9&t_2\\10&x-t_1-t_2\\
11&10-x-y-z.
\end{array}
\]

For record 43 append

\[
\ell_{12}=7-y,\qquad \ell_{13}=16-2y-z,
\]

and for record 47 append

\[
\ell_{12}=9-y,\qquad \ell_{13}=9-z.
\]

For record 44 use affine coordinates

\[
(x,y,t_1,t_2,t_3,s_1,s_2,s_3,s_4)
\]

and rows

\[
\begin{array}{c|l}
0&2-x\\1&3-x-y\\
2&s_1\\3&s_2\\4&s_3\\5&s_4\\
6&1-s_1-s_2-s_3-s_4\\
7&t_1\\8&t_2\\9&t_3\\10&y-t_1-t_2-t_3\\
11&2-y\\12&x\\13&x+y-1.
\end{array}
\]

Exact enumeration gives rank ten, a four-dimensional dependency kernel, and
exactly 90 simple vertices for each matrix.  The complements of the nine
active facets at those vertices agree exactly with the frozen 90-member
record-43, record-44, or record-47 blocker.

## Exact Coxeter metric equations

Let (K_R) be any (14\times4) basis of

\[
\ker(A_R^{\mathsf T}).
\]

Let (C) be a normalized Coxeter Gram matrix, with diagonal two.  For
intersecting walls,

\[
C_{ij}=-2\cos(\pi/m_{ij}),\qquad m_{ij}\in\{2,3,\ldots\},
\]

while a disjoint pair has (C_{ij}\leq-2), strictly below (-2) in the
compact case.  Positive rescaling of the displayed facet covectors is encoded
by (U=\operatorname{diag}(u_0,\ldots,u_{13})), (u_i>0).

There is a Lorentzian covector form (Q) for this fixed projective realization
with normalized wall Gram (C) if and only if

\[
C U K_R=0,\qquad
\operatorname{rank}C=10,\qquad
\operatorname{inertia}C=(9,1,4).
\tag{1}
\]

Indeed, after rescaling the facet rows, both (C) and
(A_RQA_R^{\mathsf T}) must have the same four-dimensional kernel.  Conversely,
equality of these kernels gives a unique nondegenerate form on the
ten-dimensional quotient.

For compactness, each of the 90 nine-wall vertex principal submatrices must be
positive definite.  Finite volume also permits the appropriate affine
parabolic rank-eight vertex matrices.  These vertex conditions additionally
place every vertex ray in the finite or ideal part of one Lorentzian sheet.

Thus a fixed label assignment gives a finite exact algebraic system: equation
(1), the six dotted inequalities, positivity of the (u_i), and exact
principal-minor/inertia tests.  Local Lanner and finite-Coxeter
classifications sharply restrict the discrete stage.  Before claiming an
unconditionally finite enumeration one must still bound any surviving
(I_2(m)) label parameter, or state an explicit finite label domain.

## Exact right-angle obstructions

Writing (Q) as a symmetric (10\times10) matrix gives 55 linear variables.
The condition (m_{ij}=2) is the linear equation

\[
a_iQa_j^{\mathsf T}=0.
\]

The exact ranks are:

| record | all 85 meeting pairs right: rank/nullity | right outside the four A5 blocks: rank/nullity | forced zero wall squares in the latter system |
|---:|---:|---:|---|
| 43 | 55 / 0 | 40 / 15 | 0, 1, 2, 11, 12, 13 |
| 44 | 55 / 0 | 39 / 16 | 0, 1, 11, 12, 13 |
| 47 | 55 / 0 | 42 / 13 | 0, 1, 2, 11, 12, 13 |

Consequently none of the displayed realizations is right-angled.  More
strongly, the natural sparse ansatz in which every non-right finite edge is
contained in one of the four order-five Lanner blocks is impossible: several
facet normals would be isotropic instead of spacelike.

There is also a construction-level obstruction.  The five simplex-factor
covectors span a five-space, the nine residual covectors span a six-space,
their intersection has dimension one, and together they span all ten
dimensions.  If all 45 cross-factor pairs were right-angled, the common line
would lie in the radical of (Q).  Hence any Coxeter metric on these fixed
product realizations needs non-right edges crossing the two factors.

For example, in every record the first dependency column can be chosen as

\[
a_0-a_1-a_2+(a_3+\cdots+a_7)-a_8-a_9-a_{10}=0
\]

for records 43 and 47, and

\[
a_0-a_1+(a_2+\cdots+a_6)-(a_7+\cdots+a_{10})=0
\]

for record 44.  If the ten displayed pairs incident to wall zero are all
right-angled, pairing the relation with (a_0) forces (a_0^2=0).  At least
one additional (m>2) edge incident to wall zero is therefore mandatory in
this realization.

## Reproduction and next step

Run

```sh
PYTHONPATH=src:. python3 scripts/explore_h9_14_multiwedge_metric.py
```

The next useful finite search is to enumerate compatible Lanner placements,
add the extra cross-factor non-right edges forced above, and apply the positive
diagonal Gale equations (CUK_R=0) before any expensive nonlinear Gram or
hyperbolic-volume computation.
