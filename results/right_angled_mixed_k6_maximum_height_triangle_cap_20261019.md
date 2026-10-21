# Clique-six triangle cap at a maximum-height mixed corner

19 October 2026. New paper proof only. No graph enumeration, optimizer,
scientific calculation, Lean or build was run. This is a local
maximum-height theorem, not a bound on the ordinary global average
number of canonical oppositions and not a dimension exclusion.

## Abstract theorem and geometric application

Let G be a graph on at most twelve vertices, with clique number exactly
six. Suppose there are real numbers p_i such that

    sum_i p_i^2 < 1,
    p_i+p_j+p_k >= 1 for every triangle {i,j,k} of G.

Then G has at most56 triangles.

The same statement applies to the actual bad-triangle graph at a
maximum-height finite vertex of a finite-volume right-angled polytope
in dimension at most twelve. Indeed, for a future unit timelike u,
maximize E(q)=-u.q over the finite vertices, put
p_i=-u.n_i/E(q), and use the actual canonical P3 opposite
q'=2q-n_i-n_j-n_k. Maximality gives the displayed triple inequality,
and the orthonormal finite-corner frame gives sum p_i^2<1. This uses
only actual finite P3 opposites as competitors. It does not maximize
over ideal vertices.

We may replace negative p_i by zero: the norm does not increase and
every triangle inequality is preserved. We may also add isolated
vertices of coordinate zero to reach twelve vertices. Hence all
coordinates below are nonnegative, on a twelve-vertex graph.

First note that G cannot have three vertex-disjoint triangles. Each
triangle costs at least1/3 of the squared norm by Cauchy; three would
contradict its strict upper bound1.

## Six-clique and outside-edge reduction

Fix S=K6 and write O for the other six vertices. For x in O put
N_x=N(x) intersect S and a_x=|N_x|. Since there is no K7, a_x<=5.
The graph on O is trianglefree: an outside triangle and two disjoint
triangles in S would be three disjoint triangles. Thus O has at most
nine edges. This is the elementary trianglefree bound e<=n^2/4:
for every edge xy, d(x)+d(y)<=n; summing gives sum d^2<=ne, while
Cauchy gives sum d^2>=4e^2/n.

For an outside edge xy put c_xy=|N_x intersect N_y|. The exact triangle
count is

    t(G)=20+sum_(x in O) C(a_x,2)+sum_(xy in E(O)) c_xy.

Suppose first some outside edge has positive c, and choose one with
maximum c.

* If c>=3, every other outside vertex z has a_z<=1. Otherwise choose
  i in N_x intersect N_y with at least two elements of N_z left after
  deleting i. This is possible because either a_z>=3, or a_z=2 and
  the common-neighbor set has at least three elements. The triangle
  xyi, a triangle consisting of z and two of those remaining elements,
  and the three unused vertices of S are disjoint triangles. Since
  x,y plus a common five-clique would give K7, c<=4. The outside
  pair-count contribution is at most20. Every outside edge other
  than xy has at most one common S-neighbor, so its total common-
  neighbor contribution is at most4+8=12. Hence t(G)<=52.

* If c=2, every other outside vertex has a_z<=2 by the same disjoint-
  triangle argument. Also a_x+a_y<=8 and a_x,a_y<=5, so
  C(a_x,2)+C(a_y,2)<=13. The other four vertices contribute at most
  four. At most nine outside edges each contribute at most two.
  Thus t(G)<=20+17+18=55.

* If c=1, the other four outside vertices again have a_z<=2. Here
  a_x+a_y<=7, so their pair contribution is at most11; the others
  contribute at most four, and the outside edges at most nine.
  Therefore t(G)<=44.

These bounds allow arbitrary noncontributing outside edges. No
assumption that O is independent has been made.

It remains to consider c_xy=0 for every outside edge. Then

    t(G)=20+sum C(a_x,2).

Assume for contradiction that t(G)>=57. The pair sum is at least37.
Some outside vertex x has a_x=5, since six vertices of degree at most
four would contribute at most36. Moreover at least two other outside
vertices have degree at least four: otherwise the pair sum is at most
10+10+4*3=32.

Let s be the unique vertex of S omitted by N_x, and let S'=S minus s.
Every outside vertex y with a_y>=2 must omit s as well. If not, choose
another outside vertex z, distinct from x,y, with a_z>=4. Such z exists
by the preceding paragraph. Take a pair {s,t} in N_y. Of the four
remaining vertices of S at least two lie in N_z; assign those to z.
Assign the two vertices left to x. They belong to N_x, since s was
already used. These are three vertex-disjoint triangles, impossible.

We may therefore regard the same triangle family as a five-clique S'
and seven outside vertices: the old vertex s and the old O. At least
two outside vertices, s and x, are adjacent to all five core vertices.
Every triangle is either in this five-clique or contains two core
vertices and one outside vertex. Original outside edges have no common
S-neighbor; an edge from s to an old outside vertex can only go to a
vertex of old S-degree one and likewise produces no such triangle.
Writing b_y for each of the seven new core degrees, we have exactly

    t(G)=10+sum_y C(b_y,2),
    sum_y C(b_y,2)>=47,   #{y:b_y=5}>=2.

## Five exact ordered-coordinate norm certificates

Order the core coordinates v1>=v2>=v3>=v4>=v5>=0. An outside vertex
with core degree j>=2 has coordinate at least

    ell_j(v)=[1-v_(j-1)-v_j]_+,    [a]_+=max(a,0).

Indeed, the two smallest coordinates in its j-element neighbor set
are no larger than v_(j-1),v_j. Their triangle with that outside
vertex is an actual graph triangle. Choosing particular outside
vertices consequently gives the lower bound

    F(v)=sum_(i=1)^5 v_i^2 + sum_(j=2)^5 n_j ell_j(v)^2

for the full squared norm. All omitted coordinates contribute
nonnegatively. Each F is convex on all of R5.

Let k be the number of universal outside vertices (core degree5).
The count at least47 forces one of the following selected patterns.
Counts refer only to the outside vertices selected for a norm lower
bound; additional vertices can be omitted.

| Case | (n5,n4,n3,n2) | Forced by the pair count |
|---|---|---|
| k>=5 | (5,0,0,0) | Select five universal vertices. |
| k=4 | (4,0,0,1) | Some other vertex has core degree at least2. |
| k=3, n4>=3 | (3,3,0,0) | Select three degree4 vertices. |
| k=3, n4<=2 | (3,2,2,0) | Reaching47 requires exactly two degree4 and two degree3 vertices. |
| k=2 | (2,4,1,0) | The five remaining vertices must include four of degree4 and another of degree at least3. |

For a vertex of degree exceeding the selected j, keep any j neighbors;
the same relaxed norm constraint holds. The elementary count in the
last two rows is: 30+6+3+3+3=45<47, whereas with two degree4 vertices
both remaining contributions must be3; and 20+3*6+2*3=44<47,
whereas four degree4 vertices require the last contribution at least3.

Here are exact convex supporting-plane certificates. P_i denotes the
sum of the first i components of grad F(v*), for i=1,2,3,4. In all
rows the full gradient sum is zero.

| (n5,n4,n3,n2) | v* | F(v*) | (P1,P2,P3,P4) |
|---|---|---|---|
| (5,0,0,0) | (2,2,2,2,2)/5 | 1 | (4,8,12,6)/5 |
| (4,0,0,1) | (2,2,2,2,2)/5 | 1 | (2,4,8,4)/5 |
| (3,3,0,0) | (12,12,12,12,12)/29 | 30/29 | (24,48,42,6)/29 |
| (3,2,2,0) | (14,14,14,14,14)/33 | 35/33 | (28,36,24,2)/33 |
| (2,4,1,0) | (16,16,16,16,14)/37 | 39/37 | (32,54,36,0)/37 |

The entries follow by differentiating the displayed quadratic F at
the indicated point, where each ell_j is positive. To see directly
why this verifies a global lower bound rather than merely a candidate
stationary point, convexity gives

    F(v)>=F(v*)+grad F(v*).(v-v*).

Summation by parts and the zero full gradient sum give

    grad F(v*).(v-v*)
      =sum_(i=1)^4 P_i*((v_i-v_(i+1))-(v_i*-v_(i+1)*)) >= 0.

Every P_i is nonnegative and v is ordered. In the first four rows v*
is constant; in the last row its sole nonzero successive difference
is at i=4, where P4=0. Thus the last inequality holds in every row
without assuming that an unknown minimizer has equal coordinates,
or that all positive-part terms remain active there.

Each forced pattern has full norm at least1. This contradicts
sum p_i^2<1 and closes the c=0 case. Together with the outside-edge
cases it proves t(G)<=56 whenever omega(G)=6.

## Frontier and limitations

Combined with the earlier maximum-height packet
`right_angled_mixed_k8_maximum_height_carrier_defect_20261019.md`, the
clique-eight, clique-seven and clique-six cases have triangle bounds
56,50 and56 respectively. Clique nine is excluded by its own norm
cost. The cases of clique number at most five are not proved here.

This proof does not establish the conjectured pure numerical statement
that at most56 arbitrary triples can have sum at least1 under squared
norm less than1. It uses graph triangles, the fixed six-clique, and
three-disjoint-triangle obstructions. It also supplies no ordinary
global average bound: only actual maximum-height corners have the
triple inequalities used above. A quantitative global carrier or
boundary estimate is still needed for the mixed dimension problem.
