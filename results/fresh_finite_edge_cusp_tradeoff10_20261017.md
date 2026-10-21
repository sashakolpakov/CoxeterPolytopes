# Seven-simplex ideal vertices and the finite-edge coefficient10

Let P be a finite-sided, finite-volume acute hyperbolic Coxeter polytope.
Every actual seven-dimensional simplex face has at most three ideal
vertices, unless all eight of its vertices are ideal. In particular, a
seven-simplex face with at least two finite vertices has at most three
ideal vertices.

Combining this theorem with the preceding eight- and nine-simplex results
gives the uniform finite-edge triangle bound

\[
 2(T_0+T_1)\le10E_{ff}
       =5n v_f-5E_{fi}\le5n v_f.                 \tag{1}
\]

Here T_j counts actual triangular two-faces with exactly j ideal vertices,
E_ff and E_fi count finite–finite and finite–ideal edges, and v_f counts
finite vertices. This is a complete internally reviewed lemma and exact
incidence inequality. It is not an ambient dimension bound. Triangles
with two or three ideal vertices and the quadrilateral budget remain
separate obligations.

The [checker](../scripts/fresh_finite_edge_cusp_tradeoff10_20261017.py) and
[saved certificate](fresh_finite_edge_cusp_tradeoff10_20261017.json) bind
the complete predecessor proofs and independently reconstruct the small
arithmetic in the four last cases using SymPy, without importing their
producer arithmetic. All predecessor files remain unchanged.

## Actual simplex stars and the exhaustive case partition

The [original-wall simplex-star proof](fresh_multi_ideal_triangle_scope_20261017.md)
supplies eight black roots B and a positive-definite supporting white
block W for an actual seven-simplex face. Every black deletion in B union
W is PSD; it is positive definite exactly when that actual vertex is
finite. The whole star is indefinite. Any inclusion-minimal indefinite
subdiagram L contains every black root and has order at most10.

This argument does not assume that the induced simplex is a Coxeter
simplex, nor that an ambient ideal vertex is simple. The star description
uses simplicity at the positive-dimensional faces involved. The order
ceiling and complete high-rank classification are explicit inherited
inputs, documented in the
[high-rank census](fresh_high_rank_quasilanner_census_20261017.md) and the
[earlier tradeoff proof](fresh_finite_edge_cusp_tradeoff_20261017.md).
The classical primary table is
[Chein1969, printed pp.15–16](https://www.numdam.org/item/M2AN_1969__3_3_3_0.pdf).
The present checker does not rerun the146262-attachment census.

There are four order-eight, four order-nine and three order-ten minimal
types. Choosing their internal whites gives

\[
 4+4\binom91+3\binom{10}2=175
\]

black/white cases. The complete partition is:

| Mechanism | Cases | Conclusion |
|---|---:|---|
| At least two affine black deletions | 77 | At most3 ideal vertices |
| Exactly one affine black deletion | 85 | At most3, except one case permitting all8 |
| Finite-black L25 reductions | 9 | At most3 |
| Affine B7, D7, E7 black blocks | 3 | At most1,2,3 respectively |
| Finite D8 black block | 1 | At most3 |

The [171-case package](fresh_finite_edge_seven_simplex_20261017.md)
proves the first three rows. The
[four-case package](fresh_seven_simplex_affine_and_d8_20261017.md)
proves the last two. The current certificate includes all175 choices,
their applicable mechanism, their ideal-count envelope and the envelope
remaining when at least two vertices are finite.

## The162 affine-deletion cases and the all-ideal qualification

If deleting a black b from L gives a connected affine diagram, its
strictly positive kernel and the PSD actual vertex Gram force every
external white root to be orthogonal to L minus b. Indeed the pairing
with the positive null combination is zero and every summand is
nonpositive. With two such black deletions their supports cover L, so
external whites are orthogonal to all of L and cannot change vertex
status.

With a unique affine black deletion b, eliminating the external white
block changes only diagonal b, by subtracting a scalar tau>=0. For each
other BLACK deletion j the vertex condition is precisely

\[
 \tau\le\frac1{(G_{L\setminus\{j\}}^{-1})_{bb}}. \tag{2}
\]

Internal-white deletions are not imposed as additional vertex tests.
All thresholds and all minimum ties are saved and exactly checked in
the171-case certificate. Among these162 cases,130 have ideal upper
bound2 and31 have upper bound3. The remaining case is L19 with internal
white8: below its common threshold it has one ideal vertex; at the
threshold seven more vertices become ideal simultaneously. Thus its
possible ideal counts are1 or8, not an intermediate number.

This is why the theorem retains the all-eight-ideal alternative. No
realization claim for that boundary is needed for (1), since its simplex
contains no finite edge.

## Arbitrary supporting whites in the nine finite-black cases

For L25 with internal whites a,9, where a=0,...,8, let F be its affine-E8
subdiagram on0,...,8. Then B=F minus a is positive definite, and the
existing white9 has a single label3 contact at F node8.

If any supporting white u makes B+u indefinite, that diagram is an
order-nine minimal indefinite diagram with positive-definite white
deletion. Its exact31 possible black/white choices are also reconstructed
by this checker; each has at most three ideal black vertices after the
affine-kernel or scalar-threshold reduction. In particular the all-ideal
exception above cannot occur with this positive-definite B.

Otherwise every B+u is PSD. The complete all-label attachment proof in
the171-case package uses two literal finite inverses of F: deleting8
has inverse diagonal(30,8,14,4,20,12,6,2), and deleting3 has inverse block
[[2,1],[1,1]] on7,8. Nonnegative inverse entries and the projection budget2
leave only three direct-white profiles: label3 at8, label4 at8, or label3
at7. The existing label3-at8 root excludes the other two by negative
two-root residual determinants at an actual black vertex.

Consequently all direct whites U have the same contact vector. Put Z
for all other external whites; Z is orthogonal to F. Eliminating Z from
the positive-definite supporting-white block gives an arbitrary matrix

\[
 K=G_U-R^tG_Z^{-1}R\succ0.
\]

At actual black deletion j the remaining vertex matrix is
K-c_j 11^t, where

\[
 (c_0,...,c_8)=
 (5/6,8/9,7/8,1,4/5,3/4,2/3,1/2,0).            \tag{3}
\]

Among j not equal to a there is a unique maximum. PSD at that maximum
forces positive definiteness at every smaller cost: a null vector there
would be both orthogonal to1 and null for K. Thus this branch has at most
one ideal vertex. The proof permits arbitrarily many direct whites and
arbitrarily long indirect tails; no finite white-rank cutoff is used.

## Independent reconstruction of the four remaining cases

For affine B7, D7 or E7 black block B, every white touching B makes B+u
an order-nine minimal indefinite diagram. The exact census, including its
literal isomorphism maps, gives respectively the following complete
direct-white contact supports:

| Black block | Supports; every contact is label3 |
|---|---|
| affine B7 | {0}, {1} |
| affine D7 | {0}, {1}, {6}, {7} |
| affine E7 | {4}, {7}, {4,7} |

Repeated profiles are impossible: at a maximizing deletion their
projection cost exceeds1, and two identical norm-two rows already have
negative residual determinant at native right contact. Nonright native
contacts worsen it. For distinct profiles the exact pair determinants
exclude all but a single B7 profile; a D7 cross-end pair with native
right contact; or the two single-leaf E7 profiles with native magnitude0
or1. Every set of three distinct profiles contains an excluded pair.
The recorded zero-determinant thresholds0 or1 exclude all larger Coxeter
labels at once.

Let U be the remaining direct whites and Z the arbitrary whites right
to B. Eliminating Z subtracts the same matrix

\[
 D=R^tG_Z^{-1}R\succeq0
\]

from every direct-white residual H_b. In addition D is entrywise
nonnegative, because G_Z is a positive-definite matrix with nonpositive
offdiagonals and R is nonpositive. For one direct white, a scalar
threshold and its ties give the stated bounds. For a D7 pair and the
E7 pair with native magnitude1, an existing residual has a strictly
positive null vector. PSD after subtracting D forces D=0; the ideal
counts remain2 and3 respectively.

The only remaining affine case is the E7 right pair. Write
D=[[x,z],[z,y]]. Deletions4 and7 have uncorrected matrices diag(2,1/2)
and diag(1/2,2), hence

\[
 0\le x,y\le1/2,\qquad0\le z\le\sqrt{xy}\le1/2.
\]

At deletions0,2,3,5,6 the corrected determinant lower bounds are,
respectively,5/16,11/36,1/4,11/36,1/4; their diagonals are positive.
At deletion1 the matrix before correction is

\[
 \begin{pmatrix}9/8&-1/8\\-1/8&9/8\end{pmatrix}.
\]

Its corrected diagonals are at least5/8 and its offdiagonal magnitude at
most5/8. Singularity would require x=y=z=1/2, which gives determinant
-1/4 at deletion4 and is impossible. Only deletions4 and7 can be ideal.
This establishes the affine E7 upper bound3, including its nonright pair.

For the finite-D8 case, first use the same31 finite-black rank-nine
alternatives if some B+u is indefinite. Otherwise existing white9
extends B to affine E8. Its standard-to-original node map is

\[
 (0,1,2,3,4,5,6,7,8)\mapsto(2,0,1,9,3,4,5,6,7).
\]

Internal white9 is standard node3, and the other existing white is the
label4-at8 profile. The same complete E8 attachment list applies. At
actual black deletion1 this existing profile, paired with any additional
profile A,B,C, has residual determinant -4/3,-28/9,-10/9 at native right
contact. Thus no second direct white exists. Eliminating every other
white leaves a scalar correction. The existing profile costs at black
nodes0,1,2,4,5,6,7,8 are

\[
 5/3,16/9,7/4,8/5,3/2,4/3,1,0.
\]

Only node1 maximizes this list, so this branch has at most one ideal
vertex. The alternative rank-nine branch has at most three.

The present checker reconstructs all these inverses, matrices,
provenance maps, pair exclusions, positive null vectors and31 rank-nine
threshold lists with a separate SymPy implementation. The arbitrary-Z
argument is the written all-size proof, not a finite experiment on tails.
Together the175 cases prove the theorem in the opening paragraph.

## From actual simplex blocks to coefficient10

The [finite-ended-port theorem](fresh_finite_ended_simplex_ports_20261017.md)
assigns every finite edge to a unique maximal actual simplex with at
least two finite vertices. Include isolated finite edges as dimension-one
blocks. For a block with f finite and i ideal vertices,

\[
 \sum\binom f2=E_{ff},\qquad
 T_0+T_1=\sum\left[\binom f3+i\binom f2\right]. \tag{4}
\]

Its twice-triangle contribution per finite edge is
2i+2(f-2)/3. The original-star theorem bounds block dimension by9.
The admissible domains and their maxima are:

| Block dimension | Required domain | Maximum coefficient |
|---|---|---:|
| At most6 | f>=2, f+i<=7 | 10 |
| 7 | f>=2, f+i=8, i<=3 | 8 |
| 8 | f>=2, f+i=9, i<=3 | 26/3 |
| 9 | (f,i)=(9,1),(8,2),(7,3) | 28/3 |

The dimension-eight restriction is the
[coefficient12 predecessor](fresh_finite_edge_cusp_tradeoff12_20261017.md);
the dimension-nine possibilities are the
[coefficient14 predecessor](fresh_finite_edge_cusp_tradeoff_20261017.md).
The checker enumerates every integer pair in these domains. Therefore
(4) gives the first inequality in (1). Finite vertices are simple, so
2E_ff+E_fi=n v_f gives the equality and final inequality there.

The only formal equality type in this relaxation is a six-simplex with
(f,i)=(2,5). This report neither realizes nor excludes that type. Thus
10 is the maximum of the proved count envelope, not a claim of geometric
sharpness.

## Exact remaining cusp and quadrilateral budget

Retain the notation of the
[finite-volume truncation report](fresh_finite_volume_simplex_truncation_20261017.md).
For n>=4 put alpha_n=5/eta_n-1 and

\[
 \eta_n=\begin{cases}4(n-1)/(n-2)&n\text{ even},\\
                     4n/(n-1)&n\text{ odd},\end{cases}
 \qquad
 b_n=\begin{cases}n(n-6)/8&n\text{ even},\\
                 (n-1)(n-5)/8&n\text{ odd}.\end{cases}
\]

Let Q count quadrilateral two-faces and E(L_v) the number of edges in
the Euclidean link of cusp v. That report's inequality and the face-size
identity give

\[
 b_n v_f+\alpha_n\sum_v E(L_v)
 \le D_5(P)
 \le 10E_{ff}+2(T_2+T_3)+Q.
\]

Consequently the exact remaining necessary budget is

\[
 2(T_2+T_3)+Q\ge
 (b_n-5n)v_f+5E_{fi}+\alpha_n\sum_v E(L_v).       \tag{5}
\]

No bound for the left side is supplied here. In particular (5) does
not imply an unconditional dimension exclusion. Ideal-vertex incidence
is retained explicitly; it is not replaced by a simple-vertex formula.

## Verification and freeze

**FINAL/FROZEN.** Generation and fresh saved-data replay PASS. Syntax,
recursive source hashes and local report links PASS. The root agent
independently read the full171-case proof; frontier_literature independently
reviewed the root's four-case proof and recomputed its critical arithmetic
using the separate implementation recorded here. Both internal proof
reviews found no gap. The finite-ended-port, original-star, high-rank
classification and earlier eight-/nine-simplex results remain explicit
inherited dependencies.

This is computer-assisted exact arithmetic together with internally
reviewed written geometric arguments. It is not external peer review or
proof-assistant verification, and makes no novelty claim or new ambient
dimension claim. The earlier reports' historical outstanding-case lines
describe their individual scopes; the current composition resolves the
four remaining seven-simplex cases without editing those frozen reports.
