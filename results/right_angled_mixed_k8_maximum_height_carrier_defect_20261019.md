# Mixed RA maximum-height corners: K8 rigidity and a universal finite-edge gap

19 October2026. **New written candidate; review pending.** No scientific
job, optimization, graph census, checker, Lean or build ran. All earlier
sources and prepared audit/checker bytes remain unchanged. These are
local actual-move obstructions, not a new dimension ceiling or a global
unweighted triangle average.

## 1. Every finite edge beside a bad triangle is long

At an actual finite corner q use q^2=-1 and orthonormal outward incident
normals ni. Suppose directions i,j give a bad triangle. Its actual
opposite facet has unit normal

    m_ij=q-ni-nj.

Let direction f be any finite--finite edge at q, of length L>0.
Necessarily f differs from i,j, since those directions end ideally.
All other incident facets stay along the edge. Consequently its other
endpoint and new endpoint wall are exactly

    q_f=cosh(L)*q-sinh(L)*nf,
    r_f=sinh(L)*q-cosh(L)*nf.

The orientation is fixed by q.r_f=-sinh(L)<0. These are actual
edge/endpoint facets, not formal candidate walls. The two distinct
unit wall normals have product

    r_f.m_ij=-sinh(L)<0.

The acute normal-product gap forces sinh(L)>=1, hence

    cosh(L)>=sqrt(2).                           (1)

No three-face size assumption is used. In particular every finite edge
shorter than arcosh(sqrt2) has empty bad-triangle graphs at BOTH endpoints.
The actual-normal gap, rather than a preferred native arithmetic alphabet,
is the reason for the exclusion.

## 2. Disjoint actual moves at a maximum-height vertex

Fix a future unit timelike u and choose q maximizing H(q)=-u.q among
all actual finite vertices. Put pi=-u.ni/H(q), so sum pi^2<1. Actual
P3 oppositions and symmetric quad moves give respectively

    sum_(i in a P3 triple)pi>=1,
    pi+pj>=1                 for a backed missing pair{i,j}.   (2)

These follow by comparing their actual opposite finite vertices with q.
Their reversals and displacement formulas are the reviewed ones. A
finite edge gives pf>=tanh(L/2), so(1) gives

    pf>=sqrt(2)-1                              (3)

whenever G(q) has any edge. Thus(3) concerns ALL finite incident edges
at such a corner, not only the previously selected cosh-sqrt2 edges.

Three pairwise disjoint graph triangles are impossible: each actual
P3 triple in(2) costs at least1/3 of the squared p norm. A backed pair
disjoint from a K5 is also impossible: its cost is at least1/2, while
the clique's summed triple constraints cost at least5/9. All these
coordinates are disjoint, so their total would exceed1.

## 3. A K8 corner has exactly56 triangles and all edges ideal

Suppose G(q) contains a clique S of size8. Its triple constraints give
sum_(i in S)pi>=8/3 and squared norm at least8/9.

Any graph triangle not contained in S uses at most two of its vertices.
At least six vertices of S remain; split six into two disjoint triples.
These two clique triangles and the new triangle would be three
pairwise disjoint actual P3 triples, impossible by Section2. Therefore

    every graph triangle lies inside S,
    t(G(q))=binom(8,3)=56.                      (4)

In particular every outside direction has at most ONE neighbor in S.
If there were a finite incident edge, (3) would add squared norm at
least(sqrt2-1)^2>1/9 outside S, contradicting sum pi^2<1. Hence ALL
incident edges at this corner end ideally.

Choose an outside direction h and any two distinct i,j in S. Their
triple has the bad edge ij, but cannot be a triangle by(4). If h
is adjacent to one of i,j it is an induced path; its missing pair
contains h and a vertex of S. A symmetric quad there is impossible
because the remaining seven S vertices contain a disjoint K5, as in
Section2. Thus this path is unbacked. If h has no neighbor in that
pair the triple has just one bad edge.

In either case the three-face has an all-ideal-edge finite corner
and is neither P3 nor D3. The reviewed local polygon-excess theorem
therefore forces nu>=12. These are distinct actual three-faces for
different direction triples. Consequently

    number of nu>=12 three-faces through q
       >=binom(8,2)*(n-8),

giving at least112 in dimension12 and84 in dimension11. This improves
the type of the local defect, not its conversion to an ordinary global
vertex fraction. Its existence here is conditional on a K8 at a
height-maximizing vertex; no such vertex is asserted universally.

## 4. The clique-number-seven maximum-height case has at most50 triangles

Suppose now omega(G(q))=7 and choose S=K7. Again no three disjoint graph
triangles exist. A triangle containing at most one S vertex would leave
six clique vertices for two more disjoint triangles. Hence every graph
triangle is either inside S or consists of one outside vertex h and
two S neighbors. Let a_h=|N(h) intersect S|<=6. Thus

    t(G(q))=35+sum_(h outside S)binom(a_h,2).     (5)

If some a_h>=4 and another outside vertex z has a_z>=2, take a
triangle on z and two of its S neighbors. Five S vertices remain,
and h retains at least a_h-2>=2 neighbors among them. Use two of
these with h as a second triangle; the other three S vertices form
a third. They are disjoint, a contradiction. Thus every other
outside vertex has a<=1, and the sum in(5) is at most binom(6,2)=15.

Otherwise all outside vertices have a<=3. In dimension12 there are
at most five, so that sum is again at most5*binom(3,2)=15. The same
bound holds in dimension11, with only four outside vertices. Therefore

    omega(G(q))=7 ==> t(G(q))<=50               for n11 or12.  (6)

This proof uses only actual P3 constraints and disjoint triangle
packing, with no finite graph classification.

## 5. Remaining obstruction

Equations(4) and(6) settle the K8 and K7-containing cases of the local
triangle-count56 problem at a height maximum. K9 is already impossible
there by its9/9 clique norm cost. Graphs of clique number at most6
remain unresolved; they may have many combinatorial triangles even
without three disjoint ones. For example K5 joined to an independent7
has80 triangles and no three vertex-disjoint triangles, although its
stronger triangle-cover norm constraints exclude it at a height maximum.
That single example is not a proof for the whole remaining class.

Most importantly, a maximum-height count is not an unweighted global
triangle cap. The present112/84 local defects do not imply a fixed
fraction E12/f. A bounded-ratio/component or ordinary incidence argument
is still needed before these constraints can exclude ambient dimension12
or11. The unrestricted dimension ceiling stays12, and no new example
is claimed.
