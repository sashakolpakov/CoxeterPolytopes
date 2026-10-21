# Unrestricted symmetric-quad moves and their carrier error

18 October 2026. New written lemmas, independent review pending. These
statements hold at mixed finite corners and do not assume globally ideal
incident edges or a bounded three-face spectrum. No dimension improvement,
optimization, computation, Lean execution, or frozen-file edit is claimed.

The inputs are acuteness of outward wall normals, the distinct-normal
product gap (product0 or at most-1), the actual P3 triangle-clique lemma,
and the bad-triangle attachment lemma from the reviewed
[weighted](right_angled_mixed_weighted_triangles_research_20261018.md) and
[rigidity](right_angled_mixed_low_facet_rigidity_research_20261018.md) packets.
In particular, two bad triangular facets sharing their ideal/ideal edge
make their entire right-angled three-polyhedron P3.

## Actual cosh3 quad moves without a three-face assumption

At a finite vertex q, take orthonormal outward incident normals n_i.
For distinct i,j define

    R_i=2q-n_i-2n_j,  R_j=2q-2n_i-n_j,
    q_ij=3q-2n_i-2n_j.                            (1)

If R_i and R_j are actual outward facet normals, they and the unchanged
n_l (l outside {i,j}) are mutually orthonormal. The acute positive-Gram
projection lemma makes their future unit intersection q_ij an actual
finite vertex. The null endpoints

    u=q-n_i=(n_j+R_j)/2,
    v=q-n_j=(n_i+R_i)/2

are feasible by acuteness. In the two-face plane, the four selected
inequalities enclose exactly conv(q,u,q_ij,v). Convexity and the selected
walls therefore prove that the original two-face is this quadrilateral.
Its finite vertices are opposite and cosh(distance)=3.

Conversely a two-face with two opposite finite and two ideal vertices
whose finite diagonal has cosh3 has exactly the normals (1). These
quads give a defensible undirected move set in every right-angled
polytope, including corners with finite/finite edges and larger carriers.

## A common bad neighbor forces the quad and its carrier to be canonical

Suppose ij is a missing edge of G(q), the two-face ij is a quadrilateral
with two opposite finite vertices, and some h is a common bad neighbor
of i,j. This already forces the quad to have cosh3.

Indeed its ideal endpoints are q-n_i and q-n_j. Its opposite side
normals have forms

    r_j=lambda(q-n_i)-n_j,
    r_i=mu(q-n_j)-n_i.

Both coefficients exceed1: coefficient1 gives the bad-triangle wall
q-n_i-n_j and cuts off the alleged fourth vertex. Their products with
the known walls m_jh=q-n_j-n_h and m_ih=q-n_i-n_h are1-lambda and1-mu.
The normal-product gap forces lambda,mu>=2. The two new walls meet
orthogonally at the opposite finite vertex, so

    (lambda-1)(mu-1)=1.

Thus lambda=mu=2, giving (1).

The entire three-face on {h,i,j} is now the reflected double of P3,
denoted D3. The seven known normals are n_h,n_i,n_j,m_hi,m_hj,R_i,R_j.
Its two finite vertices are q and q_ij; its four ideal vertices are

    q-n_i, q-n_j, q-n_h, q_ij-n_h.

The first three are already actual. The last equals R_i+m_hi and is a
future-null sum of outward normals, hence feasible. The seven-wall model
is the double-P3 model proved in the preceding
[actual-move note](right_angled_mixed_actual_move_support_research_20261018.md).
Its six vertices are therefore all in the original three-face, forcing
equality by convexity. This argument assumes no facet count or nu bound.

Consequently an induced path (i,h,j) in G(q) is backed by a symmetric
quad if and only if its actual three-face is D3. If its missing-pair
two-face is any Q22 quadrilateral, this backing is automatic; a larger
carrier cannot retain that Q22 together with these two bad triangles.

## Invariant multiplicities and exact corner displacement

For the move (1), the incident normals R_i,R_j replace n_i,n_j. For every
l outside the pair,

    q_ij-R_i-n_l=q-n_i-n_l,
    q_ij-R_j-n_l=q-n_j-n_l.                        (2)

As with P3 oppositions, the normal formula for a bad triangle and its
null-endpoint feasibility show that (2) preserves the cross edges
exactly. The common bad-neighbor set of i,j is therefore identical at
both ends of the quad move.

Let d_ij be its size. Give the undirected quad edge weight d_ij/2.
Different D3 carriers with the same missing pair share one quad and one
opposite vertex; their multiplicity is deliberately included in this
single symmetric weight. At q let

    t = number of triangles of G(q),
    b_i = number of triangles containing i,
    P = number of induced two-edge paths,
    c_i = number of such paths with i an endpoint of the missing pair.

Let P_b,c_i^b restrict the last two counts to backed paths. Combining
all actual P3 oppositions of weight1 with these quad moves gives exactly

    D_actual(q)=(t+P_b)q-sum_i (b_i+c_i^b)n_i.      (3)

Its sum over all finite vertices is zero by undirected cancellation.
The hypothetical field backing every induced path is

    D_full(q)=(t+P)q-sum_i (b_i+c_i)n_i.            (4)

The missing-path error is explicit:

    D_full-D_actual=Mq-sum_i c_i^miss n_i,
    M=P-P_b, sum_i c_i^miss=2M,
    ||c^miss||_2<=sqrt(2)M.                        (5)

These are combinatorial corner formulas. Formula (4) is not automatically
causal for every graph. For example K8 disjoint union K2 disjoint union
K2 has t=56,P=0 and spatial square8*21 squared=3528>56 squared. No actual
polytope realization of that graph is asserted. It pinpoints a remaining
obstruction to a universal proof based solely on the field (4).

The field does repair the triangle concentration example K5 joined to
independent7: t=80, b=(34 five times,10 seven times), P=105 and
c=(0 five times,30 seven times). Its full time coefficient is185 and
spatial square5*34 squared+7*40 squared=16980<185 squared. The triangle-
only field had spatial square6480>80 squared. This is evidence for the
quad mechanism, not an unconditional statement that these105 paths are
backed.

## Universal charge for the missing paths

For a non-P3 three-face F, let Z(F) count its finite corners whose
three-direction graph is an induced path. Each such corner has three
distinct ideal neighboring vertices and two bad triangular facets.
Distinct bad triangles use distinct ideal/ideal edges; otherwise the
attachment lemma would make F=P3. With p finite vertices and c cusps,

    E_FI>=3Z, E_II>=2Z, 4c=E_FI+2E_II>=7Z,
    p>=Z, p+4c=2nu-4.

Elimination yields

    4Z<=nu-2.                                     (6)

For integer nu>=12 this implies Z<=2(nu-11): at nu12, Z<=2; for
nu>=13, (nu-2)/4<=2(nu-11).

The nu11 error is fully determined at any induced-path corner. Its
three ideal neighbors give c>=3, so F+c=11 and p=4F-26 leave F=7 or8.
If F7,c4,p2, the three ideal edges at the chosen finite vertex imply
there is no finite/finite edge; the other finite vertex is also entirely
ideal-edge. The prior rigidity proof gives D3, which is backed.

The only unbacked nu11 case is therefore F8,c3,p6. All three cusps are
neighbors of the chosen corner. Any edge between two of these cusps
lies in their common face with that corner, so its bad-triangle graph
equals the ideal/ideal-edge graph: the finite corner is simple in three
dimensions, so its two incident edges lie in one facet containing their
endpoints a,b and q; if ab is also an edge, that polygon is exactly the
triangle qab. This graph is a two-edge path: hence
E_II=2, E_FI=8, and E_FF=(3p-E_FI)/2=5. No second induced-path corner
exists, since its bad triangles would reuse these ideal/ideal edges.

Let E12=sum_(three-faces with nu>=12)(nu-11), and e_FF count the ambient
finite/finite edges. Every ambient finite edge belongs to binom(n-1,2)
three-faces. Summing the preceding bounds gives the universal estimate

    sum_q M(q) <= 2E12 + (1/5)binom(n-1,2)e_FF.     (7)

The sharper second term can instead retain the actual finite-edge
incidences of just the F8,c3 error carriers. Inequality (7) makes no
global ideal-edge assumption. It does not yet show that the right side
is too small to cancel the causal field.

A companion bound covers every all-ideal-edge corner of a three-face
that is neither P3 nor D3. For nu11 such a carrier is again F8,c3,p6.
If its ideal/ideal graph has j=0,1,2 edges, it has3+j finite/finite
edges. For j>=1 it has at most one all-ideal-edge corner by the same
attachment argument; for j=0 it has at most two, since three would give
a K3,3 subgraph of its planar edge graph. Thus its number is at most
(2/3)E_FF. For nu>=12, if A3 counts all-ideal-edge corners, then
3A3<=4c and p>=A3 give A3<=(nu-2)/2<=5(nu-11). Consequently the total
number of these noncanonical all-ideal triple incidences is at most

    5E12 + (2/3)binom(n-1,2)e_FF.                  (8)

## Strong restrictions at a maximum-energy finite vertex

Fix any future unit u and choose q maximizing H(q)=-u dot q among
finite vertices. Put p_i=-u dot n_i/H(q), so sum p_i squared<1. Every
actual move imposes a necessary inequality:

    P3 triple S:                sum_(i in S)p_i>=1;
    symmetric-quad pair {i,j}: p_i+p_j>=1;
    finite edge with cosh=sqrt2: p_i>=sqrt2-1.

These follow directly from H(q')<=H(q), and hold at mixed corners and
with arbitrary larger carriers elsewhere.

In particular, two disjoint symmetric-quad pairs are impossible: their
four coordinates would have square sum at least1. The backed missing
pairs form an intersecting family of pairs, hence a star or a triangle.
Other useful forbidden disjoint packings are a quad pair with two P3
triples (norm cost1/2+1/3+1/3>1), or a quad pair with a P3 triple and
a cosh-sqrt2 finite-edge direction (cost1/2+1/3+(sqrt2-1) squared>1).

A backed missing pair also cannot have K5 in its common bad-neighbor
graph: summing the ten P3 triple constraints there gives coordinate
square cost at least5/9, disjoint from the pair's cost at least1/2.

For example, suppose G(q) is K4 joined to independent8. Its P3
constraints alone are consistent with square sum<1 (choose4/9 on the
core and1/9 on the eight leaves, giving8/9). But any backed missing
pair would be two leaves. Its complement has two disjoint P3 triples,
each using two core vertices and one other leaf. The forbidden packing
therefore excludes every quad pair. All112 induced paths at q are
unbacked, and must be carried by the F8,c3 error type or nu>=12 faces.
This is a concrete local-defect theorem, not a construction or a
universal dimension exclusion.
