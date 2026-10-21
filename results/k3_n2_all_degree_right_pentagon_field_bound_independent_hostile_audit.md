# Independent hostile audit: all-degree compact right pentagons

**Verdict (2026-08-11): PASS under the literal full-trace and cyclic-ground-
field hypotheses.**  The recurrence uses the degree-dependent norm floor in
both variables; the general field argument uses the complete cyclic packet,
not one pair coordinate.

## 1. Degree and norm floor

The rank-three module has full trace rank (3d=\rho(X)\le20), so (d\le6).
For a nonadjacent wall pair put (x_i=A_i-1).  Integral reflections give
(4x_i\in\mathcal O), independently of a primitive choice of normal: the
pair cycle is the trace of a product of the integral rank-one operators
(T_j=1-s_j).  Its distinguished conjugate is positive, while all other
conjugates lie in ((-4,0)).  Therefore

\[
                    1\le |N(4x_i)|<4^d(x_i)_0,
                    \qquad (x_i)_0>4^{-d}.                \tag{1.1}
\]

The same floor applies to every one of the five coordinates.

## 2. Five Gram identities and the corrected upper bound

Deleting one normal leaves a four-vertex tree in the nonzero Cartan graph.
I independently expanded its Gram determinant:

\[
 1-A_{i-1}-A_i-A_{i+1}+A_{i-1}A_{i+1}=0.
\]

With (x_j=A_j-1), the five rotated identities are exactly

\[
                        x_{i-1}x_{i+1}=x_i+1.              \tag{2.1}
\]

Fix (i), write (z=x_{i+2}), and put (\varepsilon=4^{-d}).  Two
instances of (2.1) give

\[
 x_{i+3}={z+1\over x_i z-1}.
\]

Now both (x_{i+3}>\varepsilon) and (z>\varepsilon) are required:

\[
 x_i<\varepsilon^{-1}+{1+\varepsilon\over\varepsilon z}
     <\varepsilon^{-2}+2\varepsilon^{-1}
     =(4^d+1)^2-1.                                       \tag{2.2}
\]

At (d=6), this yields

\[
                 x_i<16,785,408,
                 \qquad A_i<4097^2=16,785,409.
\]

The smaller provisional value (69,648) mixed a quadratic (1/16) floor
for one variable with the degree-six floor for the other and is not used.

## 3. Complete cyclic packet and field finiteness

The nonzero graph is the five-cycle (C_5).  Decomposing closed walks at
repeated vertices shows that its simple nonconstant cycles are the five
pair cycles and the directed length-five cycle (c_5).  All are traces of
products of the integral operators (T_i=1-s_i), so they lie in the
returned order without assuming integral directed Cartan entries.

Symmetry equates the two orientations and gives

\[
                   c_5^2=4^5\prod_{i=0}^4 A_i.
\]

Consequently every conjugate of every packet element has magnitude at most

\[
 H=32\cdot4097^5=36,938,546,139,294,728,224.
\]

The theorem assumes that this complete packet generates (k).  A minimal
generating subpacket has at most two elements because (d\le6).  Among
(\theta_1+c\theta_2), (c=0,\ldots,15), one element separates all
embeddings: at most ({6\choose2}=15) coefficients are bad.  Its conjugates
have magnitude at most (16H), and Vandermonde gives

\[
 |\operatorname{Disc}(k)|
 \le
 (1,182,033,476,457,431,303,168)^{30}.
\]

Hermite--Minkowski proves an actual finite field set.  A single pair is not
claimed to generate when (d>2).

For (d=2), the recurrence gives (1/16<x_i<288).  A pair cycle has its
two conjugates in disjoint intervals ((4,1156)) and ((0,4)), so it does
generate the quadratic field and the sharper frozen bound
(|\operatorname{Disc}(k)|\le1,343,277) follows.

The result is restricted to compact right-angled pentagons and the cyclic
ground field.  It does not cover non-right/noncompact pentagons, gratuitous
operator extensions, arbitrary polygons, or construction of a K3 carrier;
declared-host maximality is retained but unused in the bound.

The independent executable audit is
[audit_k3_n2_all_degree_right_pentagon_field_bound_independent_hostile.py](../scripts/audit_k3_n2_all_degree_right_pentagon_field_bound_independent_hostile.py),
with separate tests in
[test_k3_n2_all_degree_right_pentagon_field_bound_independent_hostile.py](../tests/test_k3_n2_all_degree_right_pentagon_field_bound_independent_hostile.py).
