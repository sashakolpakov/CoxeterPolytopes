# A sharp finite-ended cusp-port bound for one- and two-ideal triangles

Let P be a finite-sided, finite-volume acute hyperbolic Coxeter n-polytope.
Write T_j for the number of triangular two-faces with exactly j ideal
vertices. At an ideal vertex v write its combinatorial link as

\[
 L_v=\prod_{i=1}^{s_v}\Delta^{r_{vi}},\qquad
 a_v=\#\{\text{edges from v to a finite vertex}\},\qquad
 R_5(v)=\sum_i\min(r_{vi},5).
\]

Let T_j(v) count those triangles through v. Then

\[
 \boxed{T_1(v)+T_2(v)\le a_vR_5(v)},\qquad
 \boxed{T_1+2T_2\le\sum_{v\ \mathrm{ideal}}a_vR_5(v)}.       \tag{1}
\]

In particular, if every ideal vertex is simple, then

\[
 \boxed{T_1+2T_2\le5E_{fi}}.                              \tag{2}
\]

The constant5 in (2), and the uniform per-factor constant used in (1),
are attained by the classical finite-volume H9 simplex L24. No new
dimension bound or new hyperbolic example is claimed. For nonsimple
cusps, the factor sum in (1) is essential to this proof and is not
replaced by5.

## Factor faces and the valid completion argument

A port at v is an old edge through v, hence a vertex of L_v. Color it
finite or ideal according to its other endpoint. An old two-face through
v gives an edge of L_v. Every edge of a product of simplices varies in
one factor, so it lies in a unique maximal simplex factor fiber F.
Fixing the other factor coordinates makes F a simplex Delta^{r_i}.
The face correspondence for vertex links gives an actual (r_i+1)-face
P_F of P with link F at v. Thus **v is simple in P_F**, even if it is
nonsimple in P. Every finite vertex of P_F is also simple.

These local facts for general finite-volume acute chambers are recorded
in [Potyagailo–Vinberg, printed p.64](https://ems.press/content/serial-article-files/43032).
They do not say that arbitrary induced faces have Coxeter angles. The
argument below is combinatorial; its later simplex-size restrictions
come from the original-wall-star lemmas, not from treating P_F as a
Coxeter polytope in its own right.

Fix a finite-ended port va within F. Suppose vab and vac are triangular
faces. The three ports generate a three-face G of P_F. The vertices v
and a are simple in G and each already has the three neighbors among
v,a,b,c. If G had another vertex, deleting b,c would disconnect {v,a}
from it. This contradicts
[Balinski's three-connectivity theorem, printed p.432](https://msp.org/pjm/1961/11-2/pjm-v11-n2-p03-s.pdf).
Thus G is a tetrahedron. In particular vbc is triangular, regardless of
whether b,c are finite or nonsimple ideal vertices elsewhere in P_F.

All triangular neighbors of the chosen port consequently form a clique.
Let that clique together with va have k ports and generate the k-face A.
Only v and a need be simple: both already have degree k among the k+1
displayed vertices. If A contained any further vertex, deletion of the
other k−1 displayed vertices would disconnect {v,a}. Balinski's
k-connectivity excludes this. Therefore A is the whole k-simplex.
This also handles k=2 directly as the initial triangle; k=1 has no
triangle contribution.

These maximal simplex blocks **within the fixed fiber** have disjoint
sets of finite-ended ports. Indeed, if two blocks shared such a port,
the preceding completion applied to all its triangular neighbors would
merge them. They may share ideal-ended ports. Every triangle through v
with at least one other finite vertex belongs to exactly one such block.
The two types counted are T_1(v) and T_2(v).

No maximality in the ambient face poset is required when applying the
frozen ideal-count classification: it bounds **every actual simplex
face**, including these fiber-maximal blocks. No degree assertion is
made about their undisplayed ideal neighbors.

## The coefficient calculation

In one block let f>=1 be its number of finite-ended ports and i>=0 its
number of ideal-ended ports. The block simplex has dimension k=f+i,
with f finite vertices and i+1 ideal vertices, because v itself is ideal.
Its contribution is exactly

\[
 T_1(v;A)+T_2(v;A)=\binom f2+fi
       =f\left({f-1\over2}+i\right).                     \tag{3}
\]

The [original-star theorem](fresh_multi_ideal_triangle_scope_20261017.md)
gives k<=9. The [seven-simplex theorem](fresh_finite_edge_cusp_tradeoff10_20261017.md)
allows at most three ideals or all eight vertices ideal; the latter
case is excluded here by f>=1. The
[eight- and nine-simplex bounds](fresh_finite_edge_cusp_tradeoff12_20261017.md)
give at most three ideals as well. Hence for k>=7, i<=2 and

\[
 {f-1\over2}+i={k+i-1\over2}\le5.
\]

For k<=6 the same expression is at most k−1<=5. Thus the new coefficient
does not require extending the six-simplex theorem to the unproved
one-finite/six-ideal case. That case, with f=1,i=5, is explicitly retained
in the exact checker. Also k<=r_i+1 inside this factor face, so the
expression is at most r_i. Consequently each block contributes at most
f*min(r_i,5).

Sum over the disjoint finite-port sets in a fiber. For any one factor i,
each finite-ended port of the full link belongs to exactly one fiber of
that factor. Summing over all fibers therefore contributes at most
a_v*min(r_i,5), and then summing over factors proves the local bound in
(1). A T_1 triangle is counted at one cusp and a T_2 triangle at two;
this gives the global identity and inequality. Finally
sum_v a_v=E_fi. A simple cusp has one factor, proving (2).

The [standalone Fraction checker](../scripts/fresh_cusp_finite_port_triangle_capacity_20261017.py)
enumerates all28 allowed (k,f,i) types under the inherited safe envelopes,
and verifies the coefficient formula. The checks through factor rank64
only confirm the explicit all-rank argument above; no rank cutoff is
used in its proof. Formal maximizing types are not asserted realizable.

## Sharpness and the previous counterexample

The [literal L24 certificate](fresh_finite_edge_cusp_tradeoff_28over3_20261017.md)
already proves that this rational Gram defines a finite-volume H9 Coxeter
simplex with seven finite and three ideal vertices. Its full geometry was
also checked by the
[independent Fraction review](review_finite_edge_28over3_20261017.md).
The current checker reconstructs its ten-by-ten Gram and inverse without
project arithmetic, verifies all three positive affine kernels, and
recounts the triangle incidences.

At each cusp, the link is Delta^8, with seven finite-ended and two
ideal-ended ports. Thus T_1(v)=21, T_2(v)=14, and

\[
 T_1(v)+T_2(v)=35=5a_v.
\]

Globally T_1=63, T_2=21 and E_fi=21, so
T_1+2T_2=105=5E_fi. This supplies actual sharpness, not only an abstract
port distribution.

The [finite-vertex T2 counterexample](fresh_finite_vertex_t2_port_counterexample_20261017.md)
does not contradict the proof. Its two triangles share the finite–ideal
edge A--Iinf, and Iinf has the square link Delta^1 x Delta^1. At Iinf,
the two triangular link edges {A,I1} and {A,I2} belong to different factor
fibers. They never enter one application of the simple-at-cusp completion
argument. The checker reconstructs these fibers from the saved full
edge/face incidences and verifies (1) at all three cusps. We have not
restored triangle-port transitivity at the finite vertex A.

## Exact residual budget and limits

Combine (1) with the sharp finite-edge inequality and the original-face
budget from the [finite-volume report](fresh_finite_volume_simplex_truncation_20261017.md).
With its notation b_n, alpha_n, v_f and E(L_v), one obtains

\[
 \boxed{\begin{aligned}
 2T_3+Q\ \ge\ &(b_n-14n/3)v_f
       +\sum_v a_v\bigl(14/3-R_5(v)\bigr)\\
       &+\alpha_n\sum_v E(L_v)+T_1.
 \end{aligned}}                                                   \tag{4}
\]

Indeed the previous lower demand on 2(T_2+T_3)+Q loses at most
sum_v a_v R_5(v)−T_1 when the T_2 contribution is subtracted. This is a
new necessary tradeoff, not an upper payment bound for T_3 or Q. It is
vacant for the all-ideal problem in the sense that a_v=T_1=T_2=0 there;
the remaining all-ideal triangles and quadrilaterals still need control.
For simple cusps R_5(v)<=5, while multiple small factors can make R_5(v)
as large as n−1. No multiplicity-free generalization of (2) is claimed.

The two primary geometric sources above were opened on2026-10-17. The
one-factor completion argument received an independent conceptual check
from construction_audit. The parent independently read the full proof and
source, including factor-face construction, finite-port disjointness,
ideal-count scope, sharpness and residual algebra, and found no gap.
Neither read is described as a second arithmetic implementation.

**FINAL/FROZEN.** Generation and a fresh saved-data replay pass. Syntax,
all18 direct source hashes, eight local report links and the two primary
source links were checked. This is an internally reviewed incidence
lemma with exact counting and inherited geometric sharpness, not an
ambient dimension improvement or external refereeing.
