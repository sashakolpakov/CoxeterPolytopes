# P3 opposition invariants, clique vectors and a universal triangle cap

Date: 2026-10-18. Written universal geometry and exact rational calculations;
no Lean verification and no improved dimension ceiling.

## Opposition and its invariant

Let P be a finite-volume right-angled n-polytope. At a finite vertex q write
q.q=-1 and let n_1,...,n_n be the orthonormal outward normals. A triangle S
in its bad-triangle corner graph gives an actual P3 opposite vertex

    q_S = 2q-sum_{i in S} n_i.

Its three new incident normals are m_ij=q-n_i-n_j, i,j in S; the remaining
n_l are unchanged. Consequently

    h = 3q-sum_i n_i

is invariant across this opposition edge. Always h.h=n-9 and h.q=-3.
This applies to every actual P3 opposition edge, including incomplete corner
graphs. Form the finite graph whose vertices are finite vertices of P and
whose edges are P3 oppositions. Each P3 contributes one edge; distinct
triples at a vertex give distinct neighbors. Thus this graph has B edges
and its vertex-degree sum is 2B, where B counts canonical P3 faces.

If a complete K_n corner has every five-direction face canonical P5, all
its P3 opposites have complete corner graph as well. Their three internal
edges come from the P3 itself; mixed edges follow from the original bad
walls; and edges between two unchanged directions have opposite wall
2q-sum_{S} n_i-n_a-n_b, supplied by the corresponding canonical P5.

No finite opposition component in dimension n>=9 can consist entirely of
complete corners. Indeed each such corner would have degree binomial(n,3),
and its average neighbor is

    (2-9/n)q+(3/n)h.

Summing over a hypothetical closed component of N vertices gives
(9-n)sum q=3Nh. At n=9 this contradicts h!=0; at n>9 its left side has
negative square and its right side positive square. Equivalently, choose a
fixed future unit timelike u in the component and put E(q)=-u.q. Then
E(h)=3 for u equal to one component vertex. A complete corner's mean
neighbor energy is (2-9/n)E(q)+9/n>E(q), contradicting a maximum-energy
complete corner. Therefore propagation of canonical five-faces at complete
corners must encounter a defect. This conclusion alone gives no uniform
fraction of defects and no numerical dimension exclusion.

## A universal moment bound in dimension twelve

In dimension twelve, h.h=3. On each opposition component choose a fixed
future null vector z perpendicular to h; such a vector exists since the
orthogonal complement of a spacelike vector has Lorentz signature. Put

    E(q)=-z.q>0,       p_i=-z.n_i/E(q).

The equations z.h=0 and z.z=0 give sum p_i=3 and sum p_i^2=1. For every
formal triple S, whether or not it is an actual bad-graph triangle, q_S is
future unit timelike and

    r_S=E(q_S)/E(q)=2-sum_{i in S}p_i>0.

Let x_i=p_i-1/4. Then sum x_i=0, sum x_i^2=1/4 and
r_S=mu+d_S with mu=5/4 and d_S=-sum_{i in S}x_i. The zero-sum
Cauchy--Schwarz bound gives |d_S|<=3/4: the squared norm of the
three-subset indicator after subtracting its mean is 9/4. Hence

    1/2 <= r_S <= 2.

Average over all 220 triples. Direct expansion by repeated-index pattern
gives the exact identities

    average d = 0,
    average d^2 = 9/176,
    average d^3 = -(27/220) sum x_i^3,
    average d^4 = 3/440-(3/220) sum x_i^4 <= 3/440.

For clarity the subset-inclusion probabilities for one, two and three
specified indices are 1/4, 1/22, 1/220. Thus the third power-sum
coefficient is 1/4-3/22+2/220=27/220. The fourth coefficients are
1/4-7/22+12/220=-3/220 and 3/22-6/220=6/55, the latter multiplied
by (sum x_i^2)^2=1/16. Also
|sum x_i^3|<=sum |x_i|^3<=(sum x_i^2)^(3/2)=1/8.

Use the exact reciprocal expansion, whose remainder is nonnegative:

    1/(mu+d) = 1/mu-d/mu^2+d^2/mu^3-d^3/mu^4
               +d^4/(mu^4(mu+d)).

The lower bound r>=1/2 and the moment bounds yield

    average (1/r)
       <= 4/5 +(64/125)(9/176)
               +(256/625)(27/1760)+(512/625)(3/440)
        = 28808/34375.

Set g(q)=1/E(q). Actual opposition neighbors are a subset of the 220
formal triples, all terms positive. At every finite vertex,

    sum_actual_neighbors g(q_S)/g(q) <= 1267552/6875.

Sum over all finite vertices, using one choice of z for each component.
Each undirected actual edge contributes g(v)/g(u)+g(u)/g(v)>=2. Therefore
the following universal necessary condition holds:

    2B <= (1267552/6875) f < 184.372 f.

No graph realization or arithmetic assumption enters this proof. It uses
the actual Lorentz normal geometry and the P3 opposition involution.

## Why 77 multipartite vectors suffice for linear clique objectives

For graphs on twelve vertices consider the vector of clique counts
(e,t,k4,k5,...). Its convex hull equals the hull of the vectors of complete
multipartite graphs, one for each of the 77 integer partitions of twelve.
This statement concerns clique vectors; it does not assert that an
individual corner graph is multipartite.

To prove it, maximize any real linear combination of clique counts.
Partition vertices into maximal false-twin classes, each an independent
set with identical external neighborhoods. Between two such classes all
edges or no edges are present. If classes A,B have no edges between them,
no clique meets both, so with the rest fixed the objective is

    constant + |A| d_A + |B| d_B.

Replace the class with lower per-vertex coefficient by clones of the
higher one. The objective does not decrease, the two classes merge, and
all other false-twin classes remain intact. Repeating strictly decreases
the number of classes until they are pairwise adjacent. The resulting
graph is complete multipartite. Thus every linear support function of the
two finite hulls agrees, proving the claimed equality. Signed objective
coefficients cause no problem. Extra statistics that are not clique
counts require separate justification and are not included in this claim.

## A second exact aggregate repair

The triangle cap rules out both earlier normalized witness assignments,
which had f approximately 96.60042 and B=9000. It still does not make the
known scalar constraints inconsistent. Change only a6 to 4956/5 and D to
4000 in the original witness's free coordinates. Then

    c=1,
    (a6,a7,a8,a9,a10,a11)=(4956/5,100,5,1/20,1/100000,1/10^11),
    f=4890020999979/50000000000,
    (T,B,C,D,K,J)=(5590,9000,8000,4000,64000,0).

Keep the finite-corner graph weights 8000/99 on K12, 125/36 on K4,4,4,
5/2 on K6,6, and the remaining f minus their sum on the empty graph.
They realize T,2B,5C,K exactly. In particular K=16D=8C. The polygon
classes are the same four classes as before with their counts recomputed
from the new face vector.

The exact companion checker independently repeats all rows of the earlier
aggregate diagnostic: DS, all stated strict Nikulin--Khovanskii rows,
facet and facet-plus-cusp minima, weighted triangle/vertex densities,
canonical facet/sum gaps and incidence caps, polygon incidences, and the
integer-scaled truncated h/g and Macaulay growth conditions. It also checks
the twenty-facet refined gap with J=0, the sharp 11K packing inequality,
and the new triangle cap. Everything passes. In particular their slacks
are respectively

    (1267552/6875)f-2B = 30840620924379/976562500000 > 0,
    (792f+160D)/11-K   = 84112078997921/68750000000 > 0.

Thus no feasibility LP is needed to decide these rows: an exact feasible
assignment is already available. A linear optimization of their minimum
triangle density can still give a target for a sharper geometric cap.

The assignment is not a geometric example. In fact its canonical K5
saturation at K12 corners conflicts with the closed-cluster propagation
argument if one tries to realize all its graph data literally. That
conflict has not yet been converted to a quantitative scalar inequality
that excludes all nearby aggregate assignments. No gluing, cusp-width
compatibility or universal dimension reduction is asserted.

## Reproduction and scope

Run `python3 scripts/check_right_angled_mixed_opposition_moments_research_20261018.py`
for a fresh exact replay. It checks the moment coefficients, all 77
multipartite clique vectors, and every row of the surviving assignment.
The JSON records the exact values and smallest slacks. The checker uses
only standard-library rational/integer arithmetic and invokes no optimizer,
Lean process or background computation. The Lorentz, projection and
symmetrization arguments above are written proofs rather than claims of
universal machine verification.
