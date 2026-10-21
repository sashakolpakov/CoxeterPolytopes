# A maximum-energy vertex requires many large three-face carriers

Date: 2026-10-18. Written geometric argument with 120 exact integer cases.
No optimizer, Lean verification or new dimension ceiling.

## Theorem

Let P be a finite-volume right-angled polytope with at least one finite
vertex. Call a three-face large when its facet-plus-cusp sum nu3 is at
least eleven. Choose any future unit timelike vector u and maximize
E(q)=-u.q over the future unit representatives of the finite vertices.
Every such maximizing vertex lies in at least

    27 distinct large three-faces in dimension eleven,
    40 distinct large three-faces in dimension twelve.

In particular the unnormalized nonnegative excess
E3=sum_non-P3(nu3-10) is at least 27 or 40 respectively. This is a local
defect at an actual extreme-energy vertex. It does not imply a positive
fraction E3/f, and does not reach the reciprocal-square triangle cap.
The finite-vertex hypothesis is explicit; it is satisfied in the mixed
case of interest here.

## Energy constraints on actual moves

At a maximizing finite vertex q, let n_i be its orthonormal outward
normals and put p_i=-u.n_i/E(q). As u has square -1,

    sum p_i^2=1-1/E(q)^2<1.

Let k be the number of ideal incident edges and m=n-k the number of finite
incident edges. Let g count just those finite edges whose cosh length is
sqrt(2). Call them good finite edges. The R7 rigidity proof in the separate
small-three-face-dimension packet shows that every finite edge in R7 is
good. An edge need not lie in R7 to be called good here.

For a good edge the actual other endpoint is q_j=sqrt(2)q-n_j. The
maximality of E gives p_j>=t=sqrt(2)-1. Hence g*t^2<1 and in particular
g<=5, because 6t^2>1.

Every triangle in the bad-triangle corner graph G gives an actual P3
opposite q_S=2q-sum_S n_i. Maximality gives sum_S p_i>=1. If a clique
of size r>=3 is present, sum these inequalities over all its triples:
the p sum on the clique is at least r/3. Cauchy therefore costs at least
r/9 of the squared p norm. Its directions are ideal and are disjoint
from the g good finite directions, so necessarily

    r/9+g t^2<1.

This forbids cliques of respective sizes 9,8,6,5,3,3 for g=0,1,2,3,4,5.
Thus G, restricted to its k ideal directions, has clique number at most

    r(g)=8,7,5,4,2,2.

There is no use of an ideal vertex as an energy competitor. Only actual
finite neighbors and actual finite P3 opposites enter these constraints.

## Counting large carriers through this corner

Write e,tG for the edge and triangle counts of G on the k ideal directions.
Let Y(q) count the incident large three-faces. The classification of the
small nu3<=10 carriers is P3, Q8 and R7. Their finite corner ideal degrees
are respectively 3, 1 and either 0 or 2; an R7 degree-two corner has its
ideal pair forming a bad triangle. This gives three disjoint sources of
large carriers.

1. Among three ideal directions, exactly tG triples are P3. Every other
   triple must be large, contributing binomial(k,3)-tG.
2. Among two ideal directions and one finite direction, a finite edge
   whose cosh length is not sqrt(2) cannot be R7, so all its ideal pairs
   are large. For each good finite edge, at least the nonedges of G give
   large carriers. This contributes at least m*binomial(k,2)-g*e.
3. An all-finite triple meeting an edge whose cosh length is not sqrt(2)
   cannot be R7, and cannot be either of the other two small types. It
   contributes binomial(m,3)-binomial(g,3).

Triples with one ideal direction are ignored. Therefore

    Y(q)>=binomial(k,3)+m binomial(k,2)+binomial(m,3)-binomial(g,3)
           -[tG+g e].

This count uses actual ideal endpoints, actual finite edge lengths and
the carrier types. An arbitrary graph clique vector alone does not encode
it.

## The finite extremal check

For any graph on k vertices with clique number at most r, both the edge
and triangle counts are maximized by the balanced complete r-partite
graph. A self-contained proof suffices here. Clone a nonadjacent vertex
class toward the class with larger contribution to tG+g e; no clique
meets both classes, so the objective does not decrease. Cloning preserves
the clique bound. Repeating yields a complete multipartite graph with at
most r parts. Smoothing two part sizes differing by more than one
increases edges and increases triangles by the edge increase times the
total size of the other parts. Hence balanced r parts maximize both
counts. Empty parts are permitted when r>k.

Substitute these exact maxima with r=r(g) in the preceding lower bound.
There are only 57 cases (k,g) in dimension eleven and 63 in dimension
twelve, since 0<=k<=n and 0<=g<=min(n-k,5). The checker evaluates all
120 with integers. The minima are:

| Dimension | Minimum Y(q) | Minimizing (k,m,g) of this bound |
| --- | ---: | --- |
| 11 | 27 | (10,1,1), (11,0,0) |
| 12 | 40 | (11,1,1), (12,0,0) |

The balancing and integer calculation supply a lower bound; they do not
assert that the minimizing graph/energy configurations are realizable.

## Relation to the zero-excess theorem and remaining gap

The separate causal-move theorem proves dimension at most eight if all
three-faces have nu3<=10. The present theorem locates many distinct
departures from that carrier system at every maximum-energy vertex. It
also says that such a vertex's bad graph cannot contain K9, regardless
of any global carrier hypothesis.

An energy-extreme vertex can support an unbounded share of the geometric
displacement error relative to ordinary vertex counts. No estimate here
turns the local 27/40 defect into a fixed positive fraction at all finite
vertices. Such a conversion, or an overlap obstruction between the large
carriers, is still required before this can meet a global triangle cap.
No additional linear optimization is claimed or needed for this packet.

## Reproduction

Run `python3 scripts/check_right_angled_mixed_max_energy_defect_research_20261018.py`.
It checks the six Q(sqrt(2)) clique exclusions and all 120 integer cases,
then compares with the compact saved JSON. The Lorentz, move and carrier
arguments above are written proofs, using the previously established
small-carrier facts. No Lean process or optimizer is invoked.
