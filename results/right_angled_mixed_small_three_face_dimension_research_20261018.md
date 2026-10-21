# A structural dimension-eight theorem for the small three-face carriers

Date: 2026-10-18. Written universal Lorentz/adjacency proof with nine tiny
exact base checks. No Lean verification or unrestricted dimension reduction.

## Theorem and scope

**If every three-face of a finite-volume right-angled hyperbolic n-polytope
has facet-plus-cusp sum at most ten, then n<=8.**

Consequently any such polytope in dimension eleven or twelve has a
three-face with facet-plus-cusp sum at least eleven. Equivalently the
nonnegative excess E=sum_non-P3(nu3-10) is positive in these dimensions.
This is a structural contradiction for the entire zero-excess carrier
system, not another optimization of a triangle density. It does not give
a positive lower bound on E/f or reach the current universal triangle cap.

The proof uses the previously established three-dimensional classification:
nu3<=10 consists of P3, the eight-facet/two-cusp all-quadrilateral Q8 type,
and the seven-facet/three-cusp R7 type. Their finite ideal-edge degrees are
respectively {3,3}, eight copies of 1, and {0,2,2,2}. The detailed facet
and degree proof is in the frozen nu10 edge-transition packet. Isometric
uniqueness of Q8 is not needed here.

## R7 fixes its finite-edge lengths

Let q be the central finite vertex of R7 and let n_1,n_2,n_3 be its
orthonormal outward normals. Its finite edge endpoints are

    q_i=cosh(t_i)q-sinh(t_i)n_i,          t_i>0.

The bad triangular facet at q_i has outward unit normal

    r_i=sinh(t_i)q-cosh(t_i)n_i.

Indeed it is orthogonal to the other two old facet normals and passes
through q_i; orientation fixes the indicated sign. The three bad triangular
facets meet pairwise only at their shared ideal vertex on the all-ideal
triangle. They are an opposite pair in that cusp rectangle, so

    r_i.r_j=-1=-sinh(t_i)sinh(t_j).

All three positive hyperbolic sines must therefore equal one. Thus every
finite edge of R7 has cosh length sqrt(2). No general polyhedral rigidity
theorem is invoked.

## What the carrier hypothesis forces at every ambient vertex

Suppose n>=4 and fix a finite vertex q. Let k count its incident edges with
ideal other endpoint and m=n-k its finite/finite edges. A triple of ideal
directions has three ideal edges at q. Among the three small carrier types
only P3 has such a finite corner; hence every such triple is an actual P3
face, with opposite q_S=2q-sum_{i in S}n_i.

Every finite edge has cosh length sqrt(2). If k>=2, choose two ideal
directions together with this edge. The resulting three-face has degree
two at q, so is an R7 leaf. If k<=1, there are at least three finite
directions; choose two others with the edge. The resulting all-finite
corner is the center of R7. Both cases give the asserted length. Therefore
the finite endpoint in direction j is

    q_j=sqrt(2)q-n_j.

There is at least one finite vertex: every three-face under the hypothesis
is one of the three listed types, and each has finite vertices. A finite
vertex of a face is also a finite vertex of the ambient polytope.

## A causal displacement at each finite vertex

Put t=sqrt(2)-1. At every finite vertex define the vector

    D_q =sum_finite_neighbors(q_j-q)
          +(1/20)sum_P3_opposites(q_S-q).

Each sum is over actual vertices of the polytope. Both move relations are
undirected. The triple opposites are distinct, and every P3 supplies the
same relation in reverse. Consequently

    sum_q D_q=0.

On the other hand, in the orthonormal Lorentz frame at q,

    D_q=Aq-(1/20)binomial(k-1,2)sum_ideal n_i-sum_finite n_j,
    A=mt+(1/20)binomial(k,3)>0,

where the first spatial coefficient is zero for k<3. Write C=binomial(k,3)
and d=binomial(k-1,2), taking d=0 for k<3. Its squared time minus squared
spatial part is

    Delta=A^2-[m+k d^2/400]
         =m^2 t^2+m[(C/10)t-1]+d^2 k(k-9)/3600.

For every n>=9, Delta>=0, as verified below. Thus D_q is a nonzero future
causal vector. A finite sum of nonzero future causal vectors cannot be
zero, contradicting the cancellation identity. This proves the theorem.

## Only nine finite checks are needed

For k>=9 the last term in Delta is nonnegative. If m>0, then C>=84 and
(C/10)t-1>0, so Delta>0. If m=0 it is nonnegative, zero only for k=9;
A remains positive, so the vector is still nonzero future null there.

For 0<=k<=8, one has m>=9-k. The exact checker evaluates Delta at the
nine base values m=9-k and finds each strictly positive in Q(sqrt(2)).
It also checks that Delta(m+1)-Delta(m)>0 at these base values. Every
subsequent increment increases by 2t^2>0. Hence all larger m give a
strictly positive Delta as well. These finite checks are solely coefficient
arithmetic; the graph/Lorentz argument above supplies the universal proof.

## Limits of the obstruction

The theorem shows that the R7/Q8/P3 carrier system cannot close in
dimensions eleven or twelve, even with arbitrary distributions of actual
finite ideal-edge degrees. It improves on a full-K9-only propagation
argument by allowing mixed finite-edge and P3 moves at every corner.

For a general polytope, larger three-faces can remove the required P3
moves or alter finite-edge lengths. The proof has not bounded their
geometric displacement error by a dimension-independent multiple of their
unweighted count. Thus E>0 is established, but a quantitative positive
fraction E/f remains open. It would be incorrect to substitute this
conditional dimension-eight statement for the universal mixed upper bound
twelve, or to claim that it meets the reciprocal-square triangle cap.

Replay with
`python3 scripts/check_right_angled_mixed_small_three_face_dimension_research_20261018.py`.
The compact JSON records all nine exact base cases. No optimizer, Lean,
large numerical search or background worker is used.
