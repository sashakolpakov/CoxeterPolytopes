# Actual opposition components and the nine-part corner pattern

18 October 2026. Written structural deductions; independent review pending.
No new dimension bound or example is claimed. No optimization, numerical
search, Lean process, or frozen-file modification was used for this note.

The graph G(q) has one vertex for each of the twelve incident facet
directions at a finite vertex q. An edge ij means that the corresponding
two-face is a triangle with one finite and two ideal vertices. The
P3-opposition graph has the finite vertices as its vertices and one
undirected edge for each canonical three-face P3. Its degree at q is
t(G(q)), the number of triangles of G(q).

The geometric inputs are the reviewed P3 model and the low-three-face
classification in the [edge-transition packet](right_angled_mixed_nu10_edge_transition_research_20261018.md),
and, only for the quantitative boundary estimates, the row bound in the
[reciprocal-square packet](right_angled_mixed_opposition_power_two_research_20261018.md).

## Exact transfer through an actual P3

Normalize q squared=-1 and take orthonormal outward incident normals n_i.
If S={1,2,3} is a triangle in G(q), the opposite vertex and its new
incident normals are

    q'=2q-n_1-n_2-n_3,
    n'_i=q-sum_(h in S minus {i}) n_h   (i in S),
    n'_j=n_j                         (j outside S).

These are actual vertex/facet data, rather than a formal change of basis.
For i in S and j outside S,

    q'-n'_i-n'_j=q-n_i-n_j.                         (1)

For completeness, at any finite corner a wall with normal q-n_i-n_j is
equivalent to a bad triangular two-face in these directions. The forward
implication follows from the ideal endpoints

    q-n_i = n_j+(q-n_i-n_j),
    q-n_j = n_i+(q-n_i-n_j).

Every other wall has nonpositive product with these sums by acuteness.
The selected walls themselves also do, and the resulting three-wall
region in the two-face plane is exactly the ideal-endpoint triangle.
Convexity therefore gives the whole two-face. The converse is the known
triangle normal formula. Applying this at both corners in (1) proves:

**Every edge between S and its complement is preserved exactly by a P3
opposition.** Within S there is still a triangle. For u,v outside S,
the new edge uv is instead equivalent to the presence of the wall

    2q-sum_(i in S union {u,v}) n_i.                (2)

Equation (2) is not assumed to be present merely because the old five
directions form a clique.

Call a finite vertex small when every incident three-face has facet-plus-
cusp sum at most10. Its corner graph is K_k plus isolates, where k is its
number of ideal incident edges: the only small three-face with two ideal
directions at that corner is P3 or an R7 leaf, and both have the required
bad triangle between those directions.

Consequently, if q' has an incoming P3 opposition from a small corner q,
the three vertices S are true twins in G(q'): they form a triangle and
have identical neighbors outside S. Thus the number of incoming small
corners is at most

    sum_C binom(|C|,3),                            (3)

over true-twin classes C of G(q'). If the source is a complete K_12
corner, S consists of universal vertices of G(q').

## A K9-free corner is unavoidable in each opposition component

There is a stronger qualitative obstruction that does not use the
reciprocal-square constant. At any corner put t=t(G(q)), and let b_i
count the triangles of G(q) containing direction i. The exact sum of its
actual opposition displacements is

    D(q)=sum_(q' P3 neighbor) (q'-q)=t q-sum_i b_i n_i,
    D(q) squared=sum_i b_i squared-t squared.        (4a)

If t>0 and sum b_i squared<=t squared, D(q) is nonzero future causal:
its product with the future unit vector q is -t<0. In any finite
opposition component the sum of all D(q) is zero by undirected edge
cancellation. Thus each nontrivial component has a corner with

    sum_i b_i squared>t squared.                   (4b)

In particular, some direction belongs to more than t/3 of its triangles;
otherwise sum b_i squared<=(t/3)sum b_i=t squared.

**Graph lemma.** In a graph on at most twelve vertices containing K9,
sum b_i squared<=t squared. Equality is possible only when all its
triangles are contained in the chosen K9.

To prove it, the chosen K9 contributes84 triangles and28 to each of its
nine incidence counts. Let M be the number of additional triangles,
c_i their incidence counts at the old nine vertices, d_j their counts
at the at most three outside vertices, and O=sum d_j. Then

    O>=M,  sum c_i=3M-O,
    c_i<=3*8+binom(3,2)=27,
    d_j<=L=min(M,binom(11,2))=min(M,55).

The first bound on c_i counts triangles using one other old vertex and
one outside vertex, or two outside vertices. Expanding gives exactly

    sum b_i squared-t squared
      =sum c_i squared+sum d_j squared-M squared-56O
      <=27(3M-O)+LO-M squared-56O
      =81M-M squared+(L-83)O
      <=M(L-M-2)<=-2M.

The penultimate step uses L<=55<83 and O>=M. If M=0, the expression is
zero; if M>0, it is strictly negative. This proves the graph lemma.

**Consequently every nontrivial actual P3-opposition component in
dimensions9--12 contains a K9-free corner graph.** This requires neither
a minimum-facet hypothesis nor canonical five-faces. A path starting at
any K9 corner must eventually leave the K9-containing corner set.
It is an actual geometric exit statement, not a clique-moment assertion.

For the nine-part pattern below, b=(52,52,52,52,52,52,43,43,43,43,43,43),
so sum b_i squared=27318<190 squared. Complete K12 corners are also
strictly causal. More generally, any nonzero triangle pattern with the
same triangle-incidence count in all twelve directions is strictly
causal. Thus a component of balanced six doubled parts, whose degree
is160 (below the reciprocal cap), still cannot close. The obstruction
detects incidence balance that a degree cap alone misses.

## The Schur cap is componentwise, with a boundary term

This quantitative section and the following nine-part boundary estimates
are specifically in dimension12. The qualitative K9-free exit above and
the restricted quad theorem below separately cover dimensions9--12.

Write

    A=516992/3125=165.43744.

On each actual P3-opposition component, the invariant
h=3q-sum_i n_i is spacelike. Fix a future null z perpendicular to h and
put g(q)=(-z dot q)^(-2)>0. The reciprocal-square proof gives the
pointwise row bound

    sum_(q' actual P3 neighbor of q) g(q')/g(q) <= A.       (4)

For every actual edge its ratio lies in [1/4,4], since the corresponding
first-power energy ratio lies in [1/2,2]. Summing (4) over one component
and using a+b>=2 when ab=1 proves its average degree is at most A.
In particular, if f_P3 counts finite vertices incident to at least one
P3, then

    2B <= A f_P3.                                  (5)

Triangle-free corner graphs cannot dilute this cap. An even more useful
version holds for a vertex subset W. Let e(W) count internal opposition
edges and partial(W) its edge boundary. Summing (4) only over W gives

    2e(W)+(1/4)|partial(W)| <= A|W|,                (6)

because each internal edge contributes a reciprocal pair and each
boundary edge contributes at least1/4 at its endpoint in W. The choice
of z is made separately on each component before summing over W.

## Forced escape from the nine-part pattern

Let W consist of all corners whose graph is the complete multipartite
graph with six singleton parts and three doubled parts. Equivalently,
G(q) is K_12 minus three disjoint edges. Its triangle count is

    220-3*10=190:

a triple cannot contain two of the disjoint missing edges. All twelve
directions are nonisolated in G(q), hence all twelve incident edges are
ideal and there are no finite/finite moves at such a corner.

Every actual P3 opposite of a corner in W is again all-ideal: for a
triangle S, each outside vertex has at most one nonneighbor in S, so
at least two cross edges survive by (1). The three directions in S
are ideal by the P3 model. In particular, such an opposite cannot have
an empty graph.

At q in W the only true-twin class of size at least3 is its six universal
vertices. Formula (3) therefore allows at most binom(6,3)=20 incoming
small corners. Every one of these has k=12, since the common neighbors
outside S are all nine outside directions. Hence at most20 opposition
edges at q lead to complete K_12 corners.

The degree sum on W is190|W|=2e(W)+|partial(W)|. Combining with (6),

    |partial(W)| >= (4/3)(190-A)|W|
                 = (102344/3125)|W|.               (7)

At most20|W| of these edges lead to complete corners. It follows that

    at least (39844/3125)|W|                       (8)

opposition edges from W reach additional all-ideal corners whose graphs
are neither complete nor the nine-part pattern. Thus any nonempty W
forces another noncomplete corner type. In particular, actual graph
support restricted to {empty graph, K_12, nine-part pattern} is impossible
when W is nonempty. The weaker componentwise cap alone already excludes
support restricted to {empty graph, nine-part pattern}.

These are obstructions to actual corner-type support. The old77-type
multipartite relaxation preserved clique counts, not the activation
coordinate 1[t(G)>0] or actual connected components. Consequently (5)--(8)
do not by themselves exclude every graph distribution with the same old
clique moments, and no improved universal dimension bound follows here.

## Finite-edge support also cannot close on the two old types

Call a finite/finite edge clean when every three-face containing it has
sum at most10. Its three-face carriers are Q8 or R7. Their (3finite,
1ideal) quadrilaterals flip the ideal status of every other direction
between its endpoints, giving k'=n-1-k. Moreover, for any pair of new
ideal directions at the other endpoint, their three-face with the edge
is an R7 leaf. They therefore form a bad triangle. Both endpoint graphs
of a clean finite edge are consequently complete on their ideal
directions, with all other directions isolated.

Let v_0 count empty graphs, v_t count graphs K_t plus n-t isolates, and
Y(q) count incident three-faces with sum at least11. For n>=4,

    (n-1)v_0 <= v_(n-1)+2v_(n-2)+3 sum_(G(q) empty) Y(q).  (9)

To prove this, first take empty-graph vertices with k<=1. They have at
least n-1 finite edges, and all clean ones arrive at v_(n-1) or v_(n-2),
where the finite-edge degree is at most1 or2. The number of unclean
finite edges at q is at most3Y(q). At an empty-graph vertex with k>=2,
every triple containing at least two ideal directions is large, so

    Y(q)>=binom(k,3)+(n-k)binom(k,2)>=n-2.

The last expression is nondecreasing for k>=2: its successive difference
is k(n-k-1). Thus (n-1) such vertices is bounded by three times their
large-carrier incidence total. Adding the two cases proves (9).

For n=12 every q in W has Y(q)=30 exactly. The elementary large-three-
face estimate p(F)<=14(nu(F)-10) from the
[excess-fraction packet](right_angled_mixed_carrier_excess_fraction_research_20261018.md)
then gives the additional universal support budget

    42E >= 11v_0-v_11-2v_10+90|W|,                (10)

where E=sum_(nu>=11)(nu-10). This budget is subordinate to the actual
component escape (7)--(8); no new optimizer run is proposed.

## A restricted extension using actual symmetric-quad moves

The following statement has an additional global hypothesis. It must not
be applied merely because one corner or one opposition component has
only ideal incident edges:

**There is no finite-volume right-angled polytope in dimensions9--12
with at least one finite vertex, with no finite/finite edges anywhere,
and with every three-face of sum at most11.**

First classify a three-face with nu=F+c=11 when all its finite vertices
have only ideal incident edges. Write p for its finite vertices. Euler
gives p=4F-26, and its finite/ideal edges give3p<=4c. Thus F<=7, and
F>=7 from p>=0, so F=7,c=4,p=2. Its polygon-corner count is
3p+4c=22, hence it has six triangles and one quadrilateral. Every triangle
has at most one finite corner, because finite/finite edges are absent.
Each finite vertex is incident to at most two bad triangular facets:
three would give the already classified whole P3, contrary to F=7.
There are consequently at most four bad triangles. The finite-corner
total is six, so the quad contains both finite vertices, opposite each
other, and the other facets are four bad triangles and two all-ideal
triangles.

This combinatorics also fixes its geometry. At one finite vertex q take
incident normals n_i,n_j,n_k, where n_i is the quad facet and j,k is the
missing pair of the corner path graph. Its ideal quad endpoints are
u=q-n_j and v=q-n_k. If p is its opposite finite vertex, write
p=a q-c n_j-d n_k. The right angle at p and its unit timelike norm give
(a-c)(a-d)=1 and c squared+d squared=a squared-1. Since a>1 and c,d>=0,
these imply

    c=1+s, d=1+1/s, a=1+s+1/s, s>0.

The new incident quad-side normals at p are

    r_j=d(q-n_k)-n_j,  r_k=c(q-n_j)-n_k.

Its bad triangle in directions n_i,r_j has opposite wall

    p-n_i-r_j=s(q-n_j)-n_i.

The bad triangle at q in directions i,j has opposite wall
q-n_i-n_j. These two actual unit normals have product1. Acuteness
therefore forces them to be the same wall, so s=1. The three-face is
the reflected double of P3, with its finite vertices at distance
cosh(L)=3 and

    p=3q-2n_j-2n_k.                                (11)

The new incident normals are

    n'_j=2q-n_j-2n_k, n'_k=2q-2n_j-n_k;

all other ambient incident normals stay unchanged. Formula (11) is the
reflection in the absent bad-wall normal q-n_j-n_k. Distinct three-face
carriers sharing the missing pair share this single quad and this
single opposite vertex; they are not counted as separate moves.

Under the global hypotheses in bold, every three-face is P3 or this
double. The small Q8 and R7 types have finite/finite edges and are
excluded. At every finite corner every three-subset of G(q) is therefore
a triangle or a two-edge path. The complement of G(q) is a matching of
r edges, 0<=r<=floor(n/2). Every missing pair supplies precisely one
actual symmetric-quad move (11).

Form the sum of all actual P3 displacements, plus one half the sum of
these quad displacements. This is an undirected weighted move graph, so
the sum of the displacement vectors over all finite vertices is zero.
At a corner its time coefficient a and spatial coefficients are

    a=binom(n,3)-r(n-3),
    b_un=binom(n-1,2)-r                (n-2r directions),
    b_pair=binom(n-2,2)-r+2            (2r directions).

The causal margin a squared-sum b_i squared is

    n=9:   7r squared+6r,
    n=10:  1440-168r+15r squared,
    n=11:  4950-488r+25r squared,
    n=12:  12100-1008r+37r squared.

The first is nonnegative, with equality only at r=0. The other three
decrease throughout the allowed real r-interval and have positive
minimum values975,3135,7384 at r=5,5,6 respectively. The time
coefficient is positive in every case, so each displacement sum is
nonzero future causal. Their total cannot be zero. This proves the
restricted statement. The argument introduces actual bounded-distance
quad moves, but does not establish the corresponding classification or
causal inequality at arbitrary mixed finite corners.
